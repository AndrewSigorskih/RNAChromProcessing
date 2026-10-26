import logging
import shutil
import re
from csv import QUOTE_NONE
from functools import partial
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
from typing_extensions import Annotated

import pandas as pd
from pydantic import BaseModel, Field, field_validator

from .AnnotInfo import AnnotInfo
from .Labeller import Labeller
from ..utils.PoolExecutor import PoolExecutor
from ..plots import (
    set_style_white, plot_distance_to_closest, plot_length_distribution,
    plot_tpm_expressions
)
from ..utils import (
    exit_with_error, move_exist_ok, run_command, run_get_stdout, validate_tool_path
)

logger = logging.getLogger()
STRINGTIE_STAGES = (
    'stringtie_raw', 'stringtie_merge', 'bed_transforms', 'stringtie_cov',
    'xrna', 'plots'
)
TPM_PAT = re.compile(r'(?<=TPM \")([0-9]*[.])?[0-9]+(?=\";)')
X_ID_PAT = re.compile(r'(?<=gene_id \")[0-9a-zA-Z_]+(?=\";)')


class StringtieTool(BaseModel):
    tool_path: Annotated[Optional[Path], Field(validate_default=True)] = None
    stringtie_threads: int = 1
    cpus: Optional[int] = None
    
    @field_validator('tool_path')
    @classmethod
    def validate_stringtie_path(cls, val: Optional[Path]) -> Path:
        val_fn = partial(validate_tool_path, tool_name='stringtie')
        return val_fn(val)
    

    def run_stringtie(self,
                      inputs: Tuple[Path, Path],
                      out_file: Path) -> int:
        in_bam, gtf_annot = inputs
        cmd = (
            f'{self.tool_path} -o {out_file} -G {gtf_annot} '
            f'-p {self.stringtie_threads} {in_bam}'
        )
        return_code = run_command(cmd, shell=True)
        return return_code


    def run_stringtie_merge(self,
                            assembly_list: Path,
                            output_file: Path) -> int:
        logger.debug('Started stringtie merge.')
        cmd = (
            f'{self.tool_path} --merge -p {self.stringtie_threads} '
            f'-o {output_file} {assembly_list}'
        )
        return_code = run_command(cmd, shell=True)
        return return_code
    
    def run_stringtie_cov(self,
                          inputs: Tuple[Path, Path],
                          out_file: Path) -> int:
        inp_bam, inp_gtf = inputs
        cmd = (
            f'{self.tool_path} -e -B -p {self.stringtie_threads} '
            f'-G {inp_gtf} -o {out_file} {inp_bam}'
        )
        return_code = run_command(cmd, shell=True)
        return return_code


class StringtiePipeline:
    def __init__(self,
                 work_pth: Path,
                 executor: PoolExecutor,
                 stringtie: StringtieTool):
        self.executor = executor
        self.stringtie_tool = stringtie
        for subdir_name in STRINGTIE_STAGES:
            subdir = work_pth / subdir_name
            setattr(self, subdir_name, subdir)
            subdir.mkdir()
    
    def run_stringtie(self,
                      input_bams: List[Path],
                      gtf_annot: Path) -> List[Path]:
        inputs = list(
            zip(
                input_bams,
                (gtf_annot for _ in range(len(input_bams)))
            )
        )
        outputs = [
            self.stringtie_raw / f'{file.stem}.gtf' for file in input_bams
        ]
        self.executor.run_function(
            self.stringtie_tool.run_stringtie,
            inputs, outputs, override_cpus=self.stringtie_tool.cpus
        )
        return outputs
    
    def run_stringtie_merge(self,
                            raw_gtfs: List[Path]) -> None:
        assembly_lst: Path = self.stringtie_merge / 'assembly.lst'
        with open(assembly_lst, 'w') as f:
            print(*raw_gtfs, sep='\n', file=f)
        merged_file: Path = self.stringtie_merge / 'merged.gtf'
        self.stringtie_tool.run_stringtie_merge(assembly_lst, merged_file)
        assembly_lst.unlink()

    def handle_intervals(self,
                         bed_annot: Path) -> None:
        logger.debug('Started processing raw stringtie output.')
        raw_gtf: Path = self.stringtie_merge / 'merged.gtf'
        raw_bed: Path = self.bed_transforms / 'raw.bed'
        nonoverlap_bed: Path = self.bed_transforms / 'non-overlap.bed'
        counts_bed: Path = self.bed_transforms / 'counts.bed'
        xrnas: Path = self.xrna / 'xrna.bed'
        # gtf -> sorted bed
        gtf_header = ['chr', 'type', 'start', 'end', 'score', 'strand', 'misc']
        tab = pd.read_csv(
            raw_gtf, sep='\t', header=None, skiprows=2, usecols=[0,2,3,4,5,6,7],
            names=gtf_header
        )
        tab = tab[tab['type'] == 'transcript']
        tab = tab[['chr', 'start', 'end', 'misc', 'score', 'strand']]
        tab = tab.sort_values(by=['chr','start'], inplace=False)
        tab.to_csv(raw_bed, sep='\t', index=False, header=False)
        # sorted bed -> merge overlaps
        cmd = f'bedtools merge -s -c 6 -o distinct -i {raw_bed} > {nonoverlap_bed}'
        ret_1 = run_command(cmd, shell=True)
        # merge overlaps -> intersections w/ annotation
        cmd = (
            f'bedtools coverage -a {nonoverlap_bed} -b {bed_annot} '
            f'-s -counts > {counts_bed}'
        )
        ret_2 = run_command(cmd, shell=True)
        # safety
        if ret_1 or ret_2:
            exit_with_error('Bedtools operations failed!')
        # save results
        tab = pd.read_csv(
            counts_bed, sep='\t', header=None,
            names=['chr', 'start', 'end', 'strand', 'counts']
        )
        tab = tab[tab['counts'] == 0]
        tab = tab.drop(['counts'], axis=1)
        tab.to_csv(xrnas, sep='\t', index=False, header=False)
        # count intermediate results
        raw_counts = run_get_stdout(f'wc -l < {raw_bed}', shell=True)
        merged_counts = run_get_stdout(f'wc -l < {nonoverlap_bed}', shell=True)
        final_counts = tab.shape[0]
        logger.debug(
            f'{raw_counts.strip()} raw stringtie transcripts, {merged_counts.strip()} '
            f'merged intervals, {final_counts} dont overlap with annotation.'
        )
        # cleanup
        for tmp_file in (raw_bed, nonoverlap_bed, counts_bed):
            tmp_file.unlink()

    def assign_names(self) -> None:
        labeller = Labeller()
        xrnas: Path = self.xrna / 'xrna.bed'
        logger.debug('Started XRNA labelling.')
        tab = pd.read_csv(xrnas, sep='\t', header=None, names=['chr', 'start', 'end', 'strand'])
        tab = tab.sort_values(by=['chr', 'start'])  # do we actually need it here?
        tab['name'] = '.'
        tab['score'] = 100
        for i in tab.index:
            tab.at[i, 'name'] = labeller.next_label(tab.at[i, 'chr'], tab.at[i, 'start'])
        tab = tab[['chr', 'start', 'end', 'name', 'score', 'strand']]
        tab.to_csv(xrnas, sep='\t', index=False, header=False)

    def closest_gene(self,
                     bed_annot: Path) -> None:
        xrnas_bed: Path = self.xrna / 'xrna.bed'
        xrnas_tab: Path = self.xrna / 'xrna.tab'
        closest_res: Path = self.bed_transforms / 'closest.bed'
        logger.debug('Started determining closest gene.')
        # ignore overlaps (just in case, there wont be any),
        # require same strand, report signed distance, report first tie
        cmd = (
            f'bedtools closest -io -s -D a -t first -a {xrnas_bed} '
            f'-b {bed_annot} > {closest_res}'
        )
        ret = run_command(cmd, shell=True)
        if ret: exit_with_error('Running bedtools closest failed!')
        # process results
        closest_res_header = [
            'chr', 'start', 'end', 'name', 'score', 'strand', 'closest_gene_chr',
            'closest_gene_start', 'closest_gene_end', 'closest_gene_name',
            'closest_gene_score', 'closest_gene_strand', 'closest_gene_dist'
        ]
        tab = pd.read_csv(closest_res, sep='\t', header=None, names=closest_res_header)
        tab = tab.drop(columns=['score', 'closest_gene_chr', 'closest_gene_strand', 'closest_gene_score'])
        tab['closest_gene_side'] = '5\''
        tab.loc[tab['closest_gene_dist'] < 0, 'closest_gene_side'] = '3\''
        tab = tab.set_index('name', drop=True)
        tab.to_csv(xrnas_tab, sep='\t', index=True, header=True)
        # cleanup
        closest_res.unlink()

    def run_stringtie_cov(self, input_bams: List[Path]) -> None:
        xrnas_bed: Path = self.xrna / 'xrna.bed'
        xrnas_gtf: Path = self.xrna / 'xrna.gtf'
        xrnas_tab: Path = self.xrna / 'xrna.tab'
        # create X-rna GTF file
        tab = pd.read_csv(xrnas_bed, sep='\t', header=None,
                          names=['chr', 'start', 'end', 'name', 'score', 'strand'])
        tab['feature'] = 'transcript'
        tab['frame'] = '.'
        tab['attr'] = tab['name'].apply(
            lambda x: f'gene_id "{x}"; transcript_id "{x}"'
        )
        gtf_header = ['chr', 'name', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attr']
        tab = tab[gtf_header]
        tab.to_csv(xrnas_gtf, sep='\t', index=False, header=False, quoting=QUOTE_NONE)
        # run stringtie-cov
        inputs = list(
            zip(
                input_bams,
                (xrnas_gtf for _ in range(len(input_bams)))
            )
        )
        outputs: List[Path] = [
            self.stringtie_cov / f'{file.stem}.gtf' for file in input_bams
        ]
        self.executor.run_function(
            self.stringtie_tool.run_stringtie_cov,
            inputs, outputs, override_cpus=self.stringtie_tool.cpus
        )
        # add fields to xrna table
        tab = pd.read_csv(xrnas_tab, sep='\t', index_col=0)
        for output in outputs:
            name = output.stem
            cov_tbl = pd.read_csv(output, sep='\t', header=None,
                                  skiprows=2, names=gtf_header)
            cov_tbl = cov_tbl[cov_tbl['feature'] == 'transcript']
            cov_tbl['TPM'] = cov_tbl['attr'].apply(
                lambda x: re.search(TPM_PAT, x).group(0)
            ).astype(float)
            cov_tbl['X_id'] = cov_tbl['attr'].apply(
                lambda x: re.search(X_ID_PAT, x).group(0)
            )
            cov_tbl = cov_tbl.set_index('X_id')
            tab[f'{name}_TPM'] = cov_tbl['TPM']
        tab.to_csv(xrnas_tab, sep='\t', index=True, header=True)

    def make_plots(self, prefix: str) -> None:
        xrnas_tab: Path = self.xrna / 'xrna.tab'
        tab = pd.read_csv(xrnas_tab, sep='\t', index_col=0)
        # plot data
        set_style_white()
        plot_length_distribution(tab, self.plots, prefix)
        plot_distance_to_closest(tab, self.plots, prefix)
        plot_tpm_expressions(tab, self.plots, prefix)

    def run(self,
            input_bams: List[Path],
            annot_info: AnnotInfo,
            prefix: str) -> None:
        raw_gtfs = self.run_stringtie(input_bams, annot_info.gtf_annotation)
        self.run_stringtie_merge(raw_gtfs)
        self.handle_intervals(annot_info.annot_bed)
        self.assign_names()
        self.closest_gene(annot_info.annot_bed)
        self.run_stringtie_cov(input_bams)
        self.make_plots(prefix)

    def save_outputs(self,
                     prefix: str,
                     output_dir: Path,
                     keep_extras: Set[str]) -> None:
        KEEP_FULL: Dict[str, Path]  = {
            'stringtie_raw': self.stringtie_raw,
            'stringtie_merge': self.stringtie_merge,
            'stringtie_cov': self.stringtie_cov,
            'plots': self.plots
        }
        keep_extras = keep_extras & set(KEEP_FULL)
        keep_extras.add('plots') # save plots anyway
        for name in keep_extras:
            src = KEEP_FULL[name]
            dst = output_dir / name
            dst.mkdir(exist_ok=True)
            move_exist_ok(src, dst)
        main_outputs_ext = ('bed', 'tab')
        for ext in main_outputs_ext:
            src = self.xrna / f'xrna.{ext}'
            dst = output_dir / f'{prefix}.{ext}'
            shutil.move(src, dst)
