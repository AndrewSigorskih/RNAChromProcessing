import logging
import shutil
from collections import defaultdict
from os import listdir
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set

import pandas as pd

from .AnnotInfo import AnnotInfo, GroupInfo, SampleInfo
from ..utils.PoolExecutor import PoolExecutor
from ..utils import (
    exit_with_error, find_in_list,  move_exist_ok, run_command, run_get_stdout, validate_tool_path
)

CHUNKSIZE = 10_000_000
STAGES = (
    'bam_to_bed', 'filter_bed', 'filter_ids', 'filter_bam', 'merge_bam', 'sort_bam',
)
INTERMEDIATE_FILE_KEYS = (
    'bed', 'bed_filtered_ids', 'filtered_ids', 'filtered_bam', 'merged_bam'
)
logger = logging.getLogger()


def find_samples_files(
        samples: List[SampleInfo],
        file_dir: Path,
        field_name: str,
        file_ext: Optional[str] = None
    ) -> None:
    files = listdir(file_dir)
    if file_ext:
        files = [x for x in files if x.endswith(file_ext)]
    for sample in samples:
        file_name = find_in_list(sample.sample_id, files)
        if not file_name:
            logger.warning(f'Could not find file {sample.sample_id} in {file_dir}!')
            continue
        sample.files_map[field_name] = file_dir / file_name


def _bam_to_bed(sample: SampleInfo, output: Path) -> int:
    cmd = f"bedtoold bamtobed -i {sample.files_map['input_bam']} > {output}"
    sample.files_map['bed'] = output
    return_code = run_command(cmd, shell=True)
    return return_code


def _filter_bed(sample: SampleInfo, annot_bed: Path) -> int:
    # require same or opposite strandness
    strand_flag = '-s' if sample.true_strand else '-S'
    counts_file = sample.files_map['bed_filtered_ids'].with_suffix('.counts')
    input_bed = sample.files_map['bed']
    coverage_cmd = (
        f'bedtools coverage -a {input_bed} -b {annot_bed} '
        f'{strand_flag} -counts > {counts_file}'
    )
    return_code = run_command(coverage_cmd, shell=True)
    if return_code != 0:
        return return_code # error will be managed by PoolExecutor
    # process resulting table by chunks
    for chunk in pd.read_csv(
        counts_file, sep='\t', header=None, usecols=[3, 7], chunksize=CHUNKSIZE
    ):
        chunk = chunk[chunk[7] == 0][3]
        #chunk = '@' + chunk.apply(str)
        chunk.to_csv(sample.files_map['bed_filtered_ids'], header=False, index=False, mode='a')
    #counts_file.unlink() # TODO remove this tmp file
    return 0


def _filter_ids(sample: SampleInfo, output: Path) -> int:
    sample.files_map['filtered_ids'] = output
    if 'input_filter_lst' in sample.files_map:
        cmd = (
            r"awk 'FNR==NR{a[$1];next}($1 in a){print}' "
            f"{sample.files_map['bed_filtered_ids']} {sample.files_map['input_filter_lst']} > {output}"
        )
        return_code = run_command(cmd, shell=True)
        return return_code
    else:
        try:
            shutil.copy(sample.files_map['bed_filtered_ids'], output)
        except Exception as _:
            import traceback
            logger.critical(traceback.format_exc())
            return 1
        return 0
    

def _filter_bam(sample: SampleInfo, output: Path) -> int:
    sample.files_map['filtered_bam'] = output
    cmd = f"samtools view -Sbh -N {sample.files_map['filtered_ids']} {sample.files_map['input_bam']} > {output}"
    return_code = run_command(cmd, shell=True)
    return return_code


def _merge_bam(in_files: Iterable[str], out_file: Path) -> int:
    cmd = [
        'samtools', 'merge', '-o', str(out_file), *in_files
    ]
    return_code = run_command(cmd)
    return return_code


def _sort_bam(in_file: Path, out_file: Path) -> int:
    cmd = ['samtools', 'sort', str(in_file), '-o', str(out_file)]
    return_code = run_command(cmd)
    return return_code


class XBamFilter:
    def __init__(self,
                work_pth: Path,
                executor: PoolExecutor):
        self.executor = executor
        for subdir_name in STAGES:
            subdir = work_pth / subdir_name
            setattr(self, subdir_name, subdir)
            subdir.mkdir()


    def run_bamtobed(self, samples_list: List[SampleInfo]) -> None:
        logger.debug("Started bam to bed conversion.")
        outputs = [self.bam_to_bed / f"{sample.sample_id}.bed" for sample in samples_list]
        self.executor.run_function(_bam_to_bed, samples_list, outputs)


    def run_filter_bed(self, samples_list: List[SampleInfo], annot_bed: Path) -> None:
        logger.debug("Started filtering bed files by genomic coordinates.")
        annot_inputs = [annot_bed for _ in range(len(samples_list))]
        for sample in samples_list:
            sample.files_map['bed_filtered_ids'] = self.filter_bed / f"{sample.sample_id}.lst"
        self.executor.run_function(_filter_bed, annot_inputs, samples_list)


    def run_filter_ids(self, samples_list: List[SampleInfo]) -> None:
        logger.debug("Started filtering id files by lists of provided ids, if any.")
        outputs = [self.filter_ids / f"{sample.sample_id}.lst" for sample in samples_list]
        self.executor.run_function(_filter_ids, samples_list, outputs)

    
    def run_filter_bam(self, samples_list: List[SampleInfo]) -> None:
        logger.debug("Started filtering bam files.")
        outputs = [self.filter_bam / f"{sample.sample_id}.lst" for sample in samples_list]
        self.executor.run_function(_filter_bam, samples_list, outputs)

    
    def run_merge_bam(self, replics_dct: Dict[str, List[SampleInfo]]) -> None:
        logger.debug("Started merging bam files.")
        inputs, outputs = [], []
        for group, samples in replics_dct.items():
            outputs.append(self.merge_bam / f'{group}.bam')
            inputs.append(
                str(sample.files_map['filtered_bam']) for sample in samples
            )  # generator of all bam files to be merged
        self.executor.run_function(_merge_bam, inputs, outputs)


    def run_sort_bam(self, replics_dct: Dict[str, List[SampleInfo]]) -> None:
        logger.debug("Started sorting bam files.")
        inputs = [self.merge_bam / f'{group}.bam' for group in replics_dct]
        outputs = [self.sort_bam / f'{group}.bam' for group in replics_dct]
        self.executor.run_function(_sort_bam, inputs, outputs)
        

    def run(self,
            file_ids: Set[str],
            annot_info: AnnotInfo,
            bam_input_dir: Path,
            filter_ids_input_dir: Optional[Path]) -> List[GroupInfo]:
        samples_list = [
            x for x in annot_info.samples if x.sample_id in file_ids          
        ]
        # manage inputs
        find_samples_files(samples_list, bam_input_dir, 'input_bam', '.bam')
        if filter_ids_input_dir:
            find_samples_files(samples_list, filter_ids_input_dir, 'input_filter_lst')
        samples_list = [sample for sample in samples_list if 'input_bam' in sample.files_map]
        if not samples_list:
            exit_with_error('No bam files found for provided ids!')
        else:
            logger.debug(f"Found {len(samples_list)} bam files, started processing.")

        # run proprocessing pipeline
        self.run_bamtobed(samples_list)
        self.run_filter_bed(samples_list, annot_info.annot_bed)
        self.run_filter_ids(samples_list)
        self.run_filter_bam(samples_list)
        # time to merge biological replicas
        replics_dct = defaultdict(list)
        for sample in samples_list:
            replics_dct[sample.sample_group].append(sample)
        self.run_merge_bam(replics_dct)
        self.run_sort_bam(replics_dct)

        # Create outputs to move further
        groups = []
        for group_id, samples in replics_dct.items():
            group = GroupInfo(group_id, samples, samples[0].true_strand)
            group.files_map['bam'] = self.sort_bam / f'{group_id}.bam'
            groups.append(group)

        return groups
