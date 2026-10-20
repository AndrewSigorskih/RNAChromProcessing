
"""
Симуляция входных данных для модуля поиска новых хроматин-ассоциированных РНК.
 
Выход (в --outdir, имена с --prefix):
  {prefix}.known_genes.gtf     разметка "известных" генов (только записи gene) -> вход модуля
  {prefix}.sampleNN.bam(.bai)  отсортированные и проиндексированные BAM с RNA-ридами -> вход модуля
  {prefix}.samples.tsv         манифест: какая бамка в какой ориентации (sense / antisense)
  {prefix}.unknown_genes.bed   BED6 с "неизвестными" генами -> только для проверки
  {prefix}.truth_counts.tsv    все гены + число ридов в каждой бамке -> для проверки
  {prefix}.genome.fa(.fai)     случайный референс (для IGV / samtools)
 
Модель:
  * Длины генов ~ логнормальное распределение с заданным средним.
  * Гены одного стренда не перекрываются (с учётом хвоста и --min-gap).
    "Неизвестные" гены не перекрываются с известными ни на каком стренде
    (если не задан --allow-antisense), чтобы фильтрация по GTF их не съела.
  * Покрытие гена: равномерно по телу + "полимеразный след" за 3'-концом:
    плотность стыкуется с телом и экспоненциально затухает
    (--readthrough-decay), обрезана на --readthrough-length.
  * Экспрессия: базовая логнормальная на ген + логнормальный шум на образец.
  * Синглтоны: равномерно по геному, случайный стренд.
  * Риды single-end, без сплайсинга и мисматчей; флаг стренда = стренд РНК
    (для "antisense"-бамок, т.е. протоколов, читающих реверс-комплемент, — инвертирован).
  * XS:A — стренд для StringTie (см. --xs-mode).
  * Тэги для отладки: XO:Z (known/unknown/singleton), XG:Z (gene_id или '.').
"""
from __future__ import annotations
 
import argparse
import math
import os
import sys
from typing import Literal, TypedDict
 
import numpy as np
import numpy.typing as npt
import pysam
 
Strand = Literal["+", "-"]
 
 
class Gene(TypedDict):
    id: str
    kind: Literal["known", "unknown"]
    chrom: str
    start: int  # 0-based
    end: int    # half-open
    strand: Strand
 
 
Reads = dict[str, np.ndarray]
 
 
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Simulate RNA-chromatin RNA-part BAMs with known/unknown genes and singletons.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("-p", "--prefix", default="sim", help="prefix for all output files")
    p.add_argument("-b", "--n-bams", type=int, default=3, help="number of BAM files")
    p.add_argument("-n", "--n-reads", type=int, default=100_000, help="reads per BAM")
    p.add_argument("--n-rc-bams", type=int, default=0,
                   help="how many of the BAMs (the last ones) mimic protocols that read the reverse "
                        "complement of RNA (strand flags inverted)")
    p.add_argument("--chrom-length", type=int, default=1_000_000, help="length of each chromosome")
    p.add_argument("--n-chroms", type=int, default=1, help="number of chromosomes")
    p.add_argument("--n-known", type=int, default=100, help="number of known (annotated) genes")
    p.add_argument("--n-unknown", type=int, default=30, help="number of unknown (unannotated) genes")
    p.add_argument("--known-frac", type=float, default=0.6, help="fraction of reads from known genes")
    p.add_argument("--singleton-frac", type=float, default=0.01,
                   help="fraction of singleton reads; the rest go to unknown genes")
    p.add_argument("--mean-gene-length", type=float, default=5_000,
                   help="mean gene length (lognormal). Human mean (~27 kb+) does not fit into the "
                        "default 1 Mb genome with 130 genes; scale --chrom-length/--n-chroms with it")
    p.add_argument("--gene-length-sigma", type=float, default=0.8, help="sigma of log gene length")
    p.add_argument("--min-gene-length", type=int, default=300)
    p.add_argument("--mean-read-length", type=float, default=30)
    p.add_argument("--read-length-sd", type=float, default=None, help="default: 10%% of mean")
    p.add_argument("--readthrough-length", type=int, default=1000,
                   help="max extension of coverage past the 3' end (0 disables)")
    p.add_argument("--readthrough-decay", type=float, default=300,
                   help="exponential decay scale of readthrough coverage, bp")
    p.add_argument("--min-gap", type=int, default=500, help="min distance between gene footprints")
    p.add_argument("--allow-antisense", action="store_true",
                   help="allow unknown genes to overlap known genes on the opposite strand "
                        "(tests strand-aware filtering)")
    p.add_argument("--expr-sigma", type=float, default=1.0, help="sigma of per-gene log expression")
    p.add_argument("--sample-noise-sigma", type=float, default=0.3,
                   help="sigma of per-sample log expression noise")
    p.add_argument("--xs-mode", choices=["read", "rna"], default="read",
                   help="XS strand tag: 'read' = alignment strand as written in the BAM (for antisense "
                        "BAMs it must be flipped together with the flag, like an aligner run without "
                        "strandness); 'rna' = true RNA strand (like HISAT2 with correct --rna-strandness)")
    p.add_argument("--seed", type=int, default=1337)
    a = p.parse_args()
 
    if not 0 <= a.known_frac <= 1 or not 0 <= a.singleton_frac <= 1:
        p.error("fractions must be in [0, 1]")
    if a.known_frac + a.singleton_frac > 1 + 1e-9:
        p.error("--known-frac + --singleton-frac must be <= 1")
    if not 0 <= a.n_rc_bams <= a.n_bams:
        p.error("--n-rc-bams must be between 0 and --n-bams")
    if a.read_length_sd is None:
        a.read_length_sd = max(1.0, 0.1 * a.mean_read_length)
    return a
 
 
# ---------------------------------------------------------------- genes
 
def footprint(start: int, end: int, strand: Strand, tail: int, gap: int) -> tuple[int, int]:
    """Interval occupied by a gene incl. readthrough tail and spacing."""
    if strand == "+":
        return start - gap, end + tail + gap
    return start - tail - gap, end + gap
 
 
def overlaps(iv: tuple[int, int], occupied: list[tuple[int, int]]) -> bool:
    s, e = iv
    return any(s < oe and os_ < e for os_, oe in occupied)
 
 
def place_genes(rng: np.random.Generator, a: argparse.Namespace, chroms: list[str]) -> list[Gene]:
    mu = math.log(a.mean_gene_length) - a.gene_length_sigma ** 2 / 2
    max_len = a.chrom_length // 4
    # footprints per (chrom, strand)
    occ: dict[tuple[str, str], list[tuple[int, int]]] = {(c, s): [] for c in chroms for s in "+-"}
    plan: tuple[tuple[Literal["known", "unknown"], int, str], ...] = (
        ("known", a.n_known, "KG"), ("unknown", a.n_unknown, "UG"))
    genes: list[Gene] = []
    for kind, n, pref in plan:
        for i in range(n):
            for _attempt in range(10_000):
                L = int(np.clip(rng.lognormal(mu, a.gene_length_sigma), a.min_gene_length, max_len))
                chrom = chroms[rng.integers(len(chroms))]
                strand: Strand = "+" if rng.random() < 0.5 else "-"
                start = int(rng.integers(0, a.chrom_length - L))
                end = start + L
                fp = footprint(start, end, strand, a.readthrough_length, a.min_gap)
                opp = "-" if strand == "+" else "+"
                if overlaps(fp, occ[(chrom, strand)]):
                    continue
                if kind == "unknown" and not a.allow_antisense and overlaps(fp, occ[(chrom, opp)]):
                    continue
                occ[(chrom, strand)].append(fp)
                genes.append(Gene(id=f"{pref}{i + 1:04d}", kind=kind, chrom=chrom,
                                  start=start, end=end, strand=strand))
                break
            else:
                sys.exit(f"ERROR: could not place {kind} gene #{i + 1} without overlaps. "
                         f"Genome too small for {a.n_known + a.n_unknown} genes of mean length "
                         f"{a.mean_gene_length:.0f}: increase --chrom-length/--n-chroms or "
                         f"decrease gene number/length.")
    genes.sort(key=lambda g: (chroms.index(g["chrom"]), g["start"]))
    return genes
 
 
# ---------------------------------------------------------------- reads
 
def read_lengths(rng: np.random.Generator, a: argparse.Namespace, n: int) -> npt.NDArray[np.int_]:
    lo = max(15, int(a.mean_read_length / 2))
    hi = max(lo + 1, int(a.mean_read_length * 2))
    return np.clip(np.rint(rng.normal(a.mean_read_length, a.read_length_sd, n)), lo, hi).astype(int)
 
 
def gene_read_centers(rng: np.random.Generator, g: Gene, n: int, tail: int,
                      decay: float) -> npt.NDArray[np.int_]:
    """Genomic centers of n reads from gene g (body uniform + decaying 3' tail)."""
    L = g["end"] - g["start"]
    if tail > 0 and decay > 0:
        tail_mass = decay * (1 - math.exp(-tail / decay))
        is_tail = rng.random(n) < tail_mass / (L + tail_mass)
    else:
        tail_mass, is_tail = 0.0, np.zeros(n, bool)
    u = rng.random(n) * L  # offset along transcript, 5'->3'
    nt = int(is_tail.sum())
    if nt:
        v = rng.random(nt)
        u[is_tail] = L - decay * np.log1p(-v * (1 - math.exp(-tail / decay)))
    ui = u.astype(np.int_)
    return g["start"] + ui if g["strand"] == "+" else g["end"] - 1 - ui
 
 
def simulate_bam_reads(rng: np.random.Generator, a: argparse.Namespace, genes: list[Gene],
                       chroms: list[str], base_expr: npt.NDArray[np.float64]
                       ) -> tuple[Reads, npt.NDArray[np.int_]]:
    n = a.n_reads
    n_known, n_single, n_unknown = rng.multinomial(
        n, [a.known_frac, a.singleton_frac, 1 - a.known_frac - a.singleton_frac])
 
    expr = base_expr * rng.lognormal(0, a.sample_noise_sigma, len(genes))
    counts = np.zeros(len(genes), int)
    for kind, total in (("known", n_known), ("unknown", n_unknown)):
        idx = np.array([i for i, g in enumerate(genes) if g["kind"] == kind], int)
        if total and len(idx) == 0:
            sys.exit(f"ERROR: reads requested for {kind} genes, but there are none")
        if len(idx):
            w = expr[idx] * np.array([genes[i]["end"] - genes[i]["start"] for i in idx])
            counts[idx] = rng.multinomial(total, w / w.sum())
 
    chrom_i: list[np.ndarray] = []
    center: list[np.ndarray] = []
    strand: list[np.ndarray] = []
    origin: list[np.ndarray] = []
    gid: list[np.ndarray] = []
    for i, (g, c) in enumerate(zip(genes, counts)):
        if not c:
            continue
        chrom_i.append(np.full(c, chroms.index(g["chrom"])))
        center.append(gene_read_centers(rng, g, c, a.readthrough_length, a.readthrough_decay))
        strand.append(np.full(c, g["strand"] == "-"))
        origin.append(np.full(c, g["kind"]))
        gid.append(np.full(c, g["id"]))
    if n_single:
        chrom_i.append(rng.integers(0, len(chroms), n_single))
        center.append(rng.integers(0, a.chrom_length, n_single))
        strand.append(rng.random(n_single) < 0.5)
        origin.append(np.full(n_single, "singleton"))
        gid.append(np.full(n_single, "."))
 
    chrom_a = np.concatenate(chrom_i)
    center_a = np.concatenate(center)
    is_rev = np.concatenate(strand)  # True = RNA on minus strand
    origin_a = np.concatenate(origin)
    gid_a = np.concatenate(gid)
    rl = read_lengths(rng, a, len(center_a))
    pos = np.clip(center_a - rl // 2, 0, a.chrom_length - rl)
    order = np.lexsort((pos, chrom_a))
    return dict(chrom=chrom_a[order], pos=pos[order], len=rl[order], rev=is_rev[order],
                origin=origin_a[order], gene=gid_a[order]), counts
 
 
def write_bam(path: str, reads: Reads, chroms: list[str], genome: dict[str, str], sample: str,
              antisense: bool, xs_mode: Literal["read", "rna"]) -> None:
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": c, "LN": len(genome[c])} for c in chroms],
              "RG": [{"ID": sample, "SM": sample}],
              "PG": [{"ID": "simulate_chromatin_rna", "PN": "simulate_chromatin_rna"}]}
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for k in range(len(reads["pos"])):
            c = chroms[reads["chrom"][k]]
            p, l = int(reads["pos"][k]), int(reads["len"][k])
            rna_rev = bool(reads["rev"][k])
            rev = rna_rev ^ antisense  # alignment strand
            xs = "-" if (rev if xs_mode == "read" else rna_rev) else "+"
            seg = pysam.AlignedSegment(out.header)
            seg.query_name = f"{sample}_{k:09d}"
            seg.flag = 16 if rev else 0
            seg.reference_id = int(reads["chrom"][k])
            seg.reference_start = p
            seg.mapping_quality = 60
            seg.cigartuples = [(0, l)]
            seg.query_sequence = genome[c][p:p + l]
            seg.query_qualities = pysam.qualitystring_to_array("I" * l)
            seg.set_tags([("RG", sample), ("NM", 0), ("XS", xs, "A"),
                          ("XO", str(reads["origin"][k])), ("XG", str(reads["gene"][k]))])
            out.write(seg)
    pysam.index(path)
 
 
# ---------------------------------------------------------------- main
 
def main() -> None:
    a = parse_args()
    rng = np.random.default_rng(a.seed)
    os.makedirs(a.outdir, exist_ok=True)
    pre = os.path.join(a.outdir, a.prefix)
    chroms = [f"chr{i + 1}" for i in range(a.n_chroms)]
 
    genome = {c: "".join(rng.choice(list("ACGT"), a.chrom_length)) for c in chroms}
    with open(f"{pre}.genome.fa", "w") as f:
        for c in chroms:
            f.write(f">{c}\n")
            s = genome[c]
            f.writelines(s[i:i + 60] + "\n" for i in range(0, len(s), 60))
    pysam.faidx(f"{pre}.genome.fa")
 
    genes = place_genes(rng, a, chroms)
    base_expr = rng.lognormal(0, a.expr_sigma, len(genes))
 
    with open(f"{pre}.known_genes.gtf", "w") as f:
        for g in genes:
            if g["kind"] == "known":
                f.write(f'{g["chrom"]}\tsim\tgene\t{g["start"] + 1}\t{g["end"]}\t.\t{g["strand"]}\t.\t'
                        f'gene_id "{g["id"]}"; gene_name "{g["id"]}"; gene_biotype "simulated";\n')
    with open(f"{pre}.unknown_genes.bed", "w") as f:
        for g in genes:
            if g["kind"] == "unknown":
                f.write(f'{g["chrom"]}\t{g["start"]}\t{g["end"]}\t{g["id"]}\t0\t{g["strand"]}\n')
 
    all_counts: list[npt.NDArray[np.int_]] = []
    manifest: list[tuple[str, str, str]] = []
    for b in range(a.n_bams):
        sample = f"sample{b + 1:02d}"
        antisense = b >= a.n_bams - a.n_rc_bams
        reads, counts = simulate_bam_reads(rng, a, genes, chroms, base_expr)
        path = f"{pre}.{sample}.bam"
        write_bam(path, reads, chroms, genome, sample, antisense, a.xs_mode)
        all_counts.append(counts)
        manifest.append((sample, os.path.basename(path), "antisense" if antisense else "sense"))
        print(f"[{sample}] {len(reads['pos'])} reads "
              f"({'antisense' if antisense else 'sense'}) -> {path}", file=sys.stderr)
 
    with open(f"{pre}.samples.tsv", "w") as f:
        f.write("sample\tbam\torientation\n")
        f.writelines("\t".join(m) + "\n" for m in manifest)
    with open(f"{pre}.truth_counts.tsv", "w") as f:
        f.write("gene_id\tkind\tchrom\tstart\tend\tstrand\t" +
                "\t".join(m[0] for m in manifest) + "\n")
        for i, g in enumerate(genes):
            f.write(f'{g["id"]}\t{g["kind"]}\t{g["chrom"]}\t{g["start"]}\t{g["end"]}\t{g["strand"]}\t' +
                    "\t".join(str(c[i]) for c in all_counts) + "\n")
 
    print(f"Done: {a.n_known} known + {a.n_unknown} unknown genes on {a.n_chroms} chrom(s), "
          f"{a.n_bams} BAM(s) in {a.outdir}", file=sys.stderr)
 
 
if __name__ == "__main__":
    main()
