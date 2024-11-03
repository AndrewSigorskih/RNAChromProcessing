import concurrent.futures
import logging
import re
from os import chdir, listdir, symlink
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Dict, List, Literal, Optional, Tuple

import pandas as pd
from pydantic import BaseModel, PositiveInt, field_validator

from ..utils import (
    check_file_exists, exit_with_error, find_in_list, run_command
)
from ..utils.parsing import fetch_genes_from_bedrc, fetch_genes_from_gtf
from ..utils.run_utils import VERBOSE
from ..plots import rna_strand_barplot, rna_strand_boxplot, set_style_white

CONTACTS_COLS = ('rna_chr', 'rna_start', 'rna_end', 'rna_strand')
CONTACTS_BED_COLS = ('rna_chr', 'rna_start', 'rna_end', 'name', 'score', 'rna_strand')
BED_COLS = ('chr', 'start', 'end', 'name', 'score', 'strand')

CHUNKSIZE = 10_000_000


logger = logging.getLogger('strand')
logging.getLogger("matplotlib").setLevel(logging.WARNING)


def _contacts_to_bed(inp_file: Path, prep_file: Path) -> None:
    # will throw exception if file is not a valid contacts table
    dat_iter = pd.read_csv(inp_file, sep='\t', usecols=CONTACTS_COLS,
                           chunksize=CHUNKSIZE)
    for chunk in dat_iter:
        chunk['name'] = '.'
        chunk['score'] = 0
        chunk.loc[:,CONTACTS_BED_COLS].to_csv(
            prep_file, sep='\t', index=False, header=False, mode='a'
        )


def _prepare_input(inp_file: Path, prep_file: Path) -> Optional[Path]:
    if inp_file.suffix == '.bed':
        symlink(inp_file, prep_file)
    elif inp_file.suffix == '.tab':
        try:
            _contacts_to_bed(inp_file, prep_file)
        except Exception:  # actually ValueError, maybe smth else?
            logger.warning(f'{inp_file} is not a contacts table!')
            return None
    else:
        logger.warning(f'Unknown file extension: {inp_file}')
        return None
    return prep_file


class DetectStrand(BaseModel):
    prefix: str = 'strand'
    input_dir: Path
    output_dir: Path
    base_dir: Path = Path('.').resolve()
    cpus: PositiveInt = 1
    plots_format: Literal['png', 'svg'] = 'png'

    gene_annotation: Path
    genes_list: Path

    def __init__(self,
                 exp_groups: Dict[str, List[str]],
                 **kwargs):
        super().__init__(**kwargs)
        # directories
        self.__setup_dirs()
        chdir(self.base_dir)
        # tmp directory
        self._work_dir = TemporaryDirectory(dir=self.base_dir)
        self._work_pth = Path(self._work_dir.name)
        # genes
        self.__load_genes()
        # input files
        self.__read_inputs(exp_groups)
        
    @field_validator('base_dir', 'input_dir', 'output_dir')
    @classmethod
    def resolve_path(cls, pth: str) -> Path:
        return Path(pth).resolve()
    
    @field_validator('gene_annotation', 'genes_list')
    @classmethod
    def check_files(cls, pth: str) -> Path:
        check_file_exists(pth)
        return Path(pth).resolve()

    def __setup_dirs(self) -> None:
        if not self.input_dir.exists():
            exit_with_error('Input directory does not exist!')
        if not self.base_dir.exists():
            self.base_dir.mkdir(parents=True)
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def __load_genes(self) -> Path:
        # fetch genes from annotation file
        file_ext = self.gene_annotation.suffix
        if file_ext == '.gtf':
            gene_annot = fetch_genes_from_gtf(self.gene_annotation, self.genes_list)
        elif file_ext == '.bedrc':
            gene_annot = fetch_genes_from_bedrc(self.gene_annotation, self.genes_list)
        else:
            exit_with_error(
                f"Unknown gene annotation format: {file_ext}! "
                "Supported formats: gtf, bedrc."
            )
        logger.info(f'{gene_annot.shape[0]} genes selected from annotation file.')

        # save annotation as tmp bed file
        gene_annot['score'] = 100
        self._gene_names = gene_annot['name'].values
        self._bed_annot = self._work_pth / 'annotation.bed'
        gene_annot.loc[:,BED_COLS].to_csv(self._bed_annot, sep='\t', index=False, header=False)

    def __read_inputs(self, exp_groups: Dict[str, List[str]]) -> None:
        self._files_map: Dict[Tuple[str, str], Path] = {}
        files_list = listdir(self.input_dir)
        for group, id_list in exp_groups.items():
            for file_id in id_list:
                filename = find_in_list(file_id, files_list)
                if not filename:
                    logger.warning(f'{file_id} not found in input directory, skipping..')
                    continue
                inp_file = self.input_dir / filename
                prep_file = (self._work_pth / filename).with_suffix('.bed')
                prep_file = _prepare_input(inp_file, prep_file)
                if not prep_file:
                    logger.warning(f'Could not read {prep_file}, skipping..')
                    continue
                self._files_map[(group, file_id)] = prep_file
        if not self._files_map:
            exit_with_error('Could not find any if the listed files in input directory!')
        logger.info(f'{len(self._files_map)} files found in input directory')

    def _get_coverage(self, key: Tuple[str, str]) -> int:
        input_bed = self._files_map[key]
        name = input_bed.stem
        res_same = self._work_pth / f'{name}_same.bed'
        res_anti = self._work_pth / f'{name}_anti.bed'
        # run coverage both ror same and anti strand
        cmd_1 = (
            f'bedtools coverage -a {self._bed_annot} -b {input_bed} '
            f'-s -counts > {res_same}'
        )
        cmd_2 = (
            f'bedtools coverage -a {self._bed_annot} -b {input_bed} '
            f'-S -counts > {res_anti}'
        )
        ret_1, ret_2 = run_command(cmd_1, shell=True), run_command(cmd_2, shell=True)
        if ret := (ret_1 or ret_2):
            return ret
        # collect results
        same_dat = pd.read_csv(
            res_same, sep='\t', header=None, index_col='index',
            usecols=[3, 6], names=['index', 'counts']
        )['counts']
        anti_dat = pd.read_csv(
            res_anti, sep='\t', header=None, index_col='index',
            usecols=[3, 6], names=['index', 'counts']
        )['counts']
        for gene_name in self._gene_names:
            same_val = same_dat[gene_name]
            anti_val = anti_dat[gene_name]
            logger.log(VERBOSE, f'Gene {gene_name}: {same_val} reads on same strand, {anti_val} on antisence.')
            self._raw_result.at[key, gene_name] = (same_val, anti_val)
        # rm tmp files
        input_bed.unlink()
        res_same.unlink()
        res_anti.unlink()
        return 0

    def calculate(self):
        # do preps
        self._raw_result = pd.DataFrame(
            data=None,
            columns=self._gene_names,
            index=pd.MultiIndex.from_tuples(
                self._files_map.keys(),
                names=['group', 'id']
            )
        )
        # run calcultaions in parallel
        logger.debug(f'Estimating selected genes coverage with {self.cpus} parallel jobs..')
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.cpus) as executor:
            futures = [
                executor.submit(self._get_coverage, inp)
                for inp in self._files_map.keys()
            ]
            results = [
                future.result() for future in concurrent.futures.as_completed(futures)
            ]
        if any(x != 0 for x in results):
            msg = f'One or several calls of bedtools coverage returned non-zero process exit code!'
            exit_with_error(msg)
    
    def counts_to_values(self) -> None:
        same_wins = lambda x: x[0] > x[1]
        anti_wins = lambda x: x[1] > x[0]
        self._result = pd.DataFrame(data=None, index=self._raw_result.index)
        self._result['same'] = self._raw_result.apply(lambda row: row.apply(same_wins).sum(),
                                                    axis=1)
        self._result['anti'] = self._raw_result.apply(lambda row: row.apply(anti_wins).sum(),
                                                    axis=1)
        
        mask = (self._result['same'] - self._result['anti'])/(self._result['same'] + self._result['anti'])
        self._result['strand'] = 'UNKNOWN'
        self._result.loc[mask > 0.75, 'strand'] = 'SAME'
        self._result.loc[mask < -0.75, 'strand'] = 'ANTI'

    def run(self) -> None:
        self.calculate()
        self.counts_to_values()
        # save outputs
        self._result.to_csv(f'{self.output_dir}/{self.prefix}_wins.tsv', sep='\t')
        self._raw_result.to_csv(f'{self.output_dir}/{self.prefix}_raw_counts.tsv', sep='\t')
        # make plots
        set_style_white()
        logger.debug('Finished calculations, plotting results..')
        rna_strand_barplot(self._result, len(self._gene_names),
                           self.output_dir, self.prefix, self.plots_format)
        rna_strand_boxplot(self._raw_result, len(self._gene_names),
                           self.output_dir, self.prefix, self.plots_format)
        logger.info('Done.')
