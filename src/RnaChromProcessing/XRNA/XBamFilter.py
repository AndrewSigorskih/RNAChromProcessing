import logging
from os import listdir
from pathlib import Path
from typing import List, Optional, Set

import pandas as pd

from .AnnotInfo import AnnotInfo, SampleInfo
from ..utils.PoolExecutor import PoolExecutor
from ..utils import (
    exit_with_error, find_in_list,  move_exist_ok, run_command, run_get_stdout, validate_tool_path
)

CHUNKSIZE = 10_000_000
STAGES = (
    'bam_to_bed', 'filter_bed', 'filter_ids', 'filter_bam', 'merge_bam', 'sort_bam',
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
    for chunk in pd.read_csv(counts_file, sep='\t', header=None,
                             usecols=[3, 7], chunksize=CHUNKSIZE):
        chunk = chunk[chunk[7] == 0][3]
        #chunk = '@' + chunk.apply(str)
        chunk.to_csv(sample.files_map['bed_filtered_ids'], header=False, index=False, mode='a')
    #counts_file.unlink() # TODO remove this tmp file
    return 0



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
        annot_inputs = [annot_bed for _ in range(len(samples_list))]
        for sample in samples_list:
            sample.files_map['bed_filtered_ids'] = self.filter_bed / f"{sample.sample_id}.lst"
        self.executor.run_function(_filter_bed, annot_inputs, samples_list)

    

    def run(self,
            file_ids: Set[str],
            annot_info: AnnotInfo,
            bam_input_dir: Path,
            filter_ids_input_dir: Optional[Path]):
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

