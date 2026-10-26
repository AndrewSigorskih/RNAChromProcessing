"""
Strand-aware Jaccard index between a truth BED and a predicted BED.
 
Global metric (base-pair level, same as `bedtools jaccard -s`):
    intervals in each file are merged per (chrom, strand), then
    J = |truth ∩ pred| / |truth ∪ pred|   in bp.
Overlaps on opposite strands do not count.
Also reports bp-level recall (|∩|/|truth|) and precision (|∩|/|pred|):
recall < 1 means predictions miss parts of true genes, precision < 1 means
predictions spill outside them (e.g. into the 3' readthrough tail).
 
Optional --per-interval TSV: for every truth interval, the merged predicted
intervals overlapping it on the same strand, its own Jaccard, and how far the
predicted 5'/3' ends are from the true ones (positive = predicted end lies
downstream of the true end, in transcript orientation).
 
Both files need at least 6 columns (chrom, start, end, name, score, strand);
strand must be '+' or '-'.
"""
 
import argparse
import sys
from collections import defaultdict
from typing import NamedTuple, TextIO
 
Key = tuple[str, str]           # (chrom, strand)
Interval = tuple[int, int]      # 0-based, half-open
 
 
class Record(NamedTuple):
    chrom: str
    start: int
    end: int
    name: str
    strand: str
 
 
def read_bed(path: str) -> list[Record]:
    records: list[Record] = []
    with open(path) as f:
        for n, line in enumerate(f, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                sys.exit(f"ERROR: {path}:{n}: need >= 6 columns (BED6) to know the strand")
            strand = cols[5]
            if strand not in ("+", "-"):
                sys.exit(f"ERROR: {path}:{n}: strand must be '+' or '-', got {strand!r}")
            start, end = int(cols[1]), int(cols[2])
            if end <= start:
                sys.exit(f"ERROR: {path}:{n}: end <= start")
            records.append(Record(cols[0], start, end, cols[3], strand))
    return records
 
 
def merge(intervals: list[Interval]) -> list[Interval]:
    out: list[Interval] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out
 
 
def by_key(records: list[Record]) -> dict[Key, list[Interval]]:
    groups: dict[Key, list[Interval]] = defaultdict(list)
    for r in records:
        groups[(r.chrom, r.strand)].append((r.start, r.end))
    return {k: merge(v) for k, v in groups.items()}
 
 
def total_length(intervals: list[Interval]) -> int:
    return sum(e - s for s, e in intervals)
 
 
def intersection_length(a: list[Interval], b: list[Interval]) -> int:
    """Overlap in bp between two sorted, merged interval lists."""
    i = j = total = 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s < e:
            total += e - s
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return total
 
 
def overlapping(iv: Interval, merged: list[Interval]) -> list[Interval]:
    return [m for m in merged if m[0] < iv[1] and iv[0] < m[1]]
 
 
def write_per_interval(out: TextIO, truth: list[Record], pred: dict[Key, list[Interval]]) -> None:
    out.write("name\tchrom\tstart\tend\tstrand\tn_pred\tpred_start\tpred_end\t"
              "overlap_bp\tjaccard\tdelta_5p\tdelta_3p\n")
    for t in sorted(truth, key=lambda r: (r.chrom, r.start)):
        hits = overlapping((t.start, t.end), pred.get((t.chrom, t.strand), []))
        base = f"{t.name}\t{t.chrom}\t{t.start}\t{t.end}\t{t.strand}\t{len(hits)}"
        if not hits:
            out.write(f"{base}\t.\t.\t0\t0.0000\t.\t.\n")
            continue
        inter = intersection_length([(t.start, t.end)], hits)
        union = (t.end - t.start) + total_length(hits) - inter
        ps, pe = hits[0][0], hits[-1][1]
        # transcript orientation: + means the predicted end is downstream of the true one
        if t.strand == "+":
            d5, d3 = ps - t.start, pe - t.end
        else:
            d5, d3 = t.end - pe, t.start - ps
        out.write(f"{base}\t{ps}\t{pe}\t{inter}\t{inter / union:.4f}\t{d5}\t{d3}\n")
 
 
def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0].strip(),
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("truth", help="truth BED6, e.g. sim.unknown_genes.bed")
    p.add_argument("pred", help="predicted BED6 from the module")
    p.add_argument("--per-interval", metavar="TSV",
                   help="also write per-truth-interval stats to this file ('-' for stdout)")
    a = p.parse_args()
 
    truth_records = read_bed(a.truth)
    truth = by_key(truth_records)
    pred = by_key(read_bed(a.pred))
 
    t_bp = sum(total_length(v) for v in truth.values())
    p_bp = sum(total_length(v) for v in pred.values())
    i_bp = sum(intersection_length(v, pred.get(k, [])) for k, v in truth.items())
    u_bp = t_bp + p_bp - i_bp
 
    def ratio(x: int, y: int) -> str:
        return f"{x / y:.4f}" if y else "NA"
 
    print("truth_bp\tpred_bp\tintersection_bp\tunion_bp\tjaccard\trecall_bp\tprecision_bp")
    print(f"{t_bp}\t{p_bp}\t{i_bp}\t{u_bp}\t{ratio(i_bp, u_bp)}\t{ratio(i_bp, t_bp)}\t{ratio(i_bp, p_bp)}")
 
    if a.per_interval:
        if a.per_interval == "-":
            write_per_interval(sys.stdout, truth_records, pred)
        else:
            with open(a.per_interval, "w") as f:
                write_per_interval(f, truth_records, pred)
 
 
if __name__ == "__main__":
    main()
