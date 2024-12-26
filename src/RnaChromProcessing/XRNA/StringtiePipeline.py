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

from .AnnotInfo import AnnotInfo, GroupInfo
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



def _convert_strand(in_gtf: Path, out_gtf: Path) -> None:
    REVERSE = {
        '+': '-',
        '-': '+',
        '.': '.'
    }
    with open(in_gtf, 'r') as in_f, open(out_gtf, 'w') as out_f:
        for line in in_f:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                fields[6] = REVERSE[fields[6]]
                line = '\t'.join(fields) + '\n'
            out_f.write(line)


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
                      group: GroupInfo,
                      gtf_annot: Path) -> int:
        in_bam = group.files_map['bam']
        out_file = group.files_map['gtf']
        raw_gtf = out_file.with_suffix('.tmp')        
        cmd = (
            f'{self.tool_path} -o {raw_gtf} -G {gtf_annot} '
            f'-p {self.stringtie_threads} {in_bam}'
        )
        return_code = run_command(cmd, shell=True)
        if return_code:
            return return_code
        # finally manage strandness HERE
        try:
            if group.true_strand:
                shutil.copy(raw_gtf, out_file)
            else:
                _convert_strand(raw_gtf, out_file)
        except Exception as _:
            import traceback
            logger.critical(traceback.format_exc())
            return 1
        return 0

    def run_stringtie_merge(self,
                            assembly_list: Path,
                            output_file: Path) -> int:
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
                      input_groups: List[GroupInfo],
                      gtf_annot: Path) -> None:
        logger.debug('Started stringtie inference.')
        gtf_inputs = [gtf_annot for _ in range(len(input_groups))]
        for group in input_groups:
            group.files_map['gtf'] = self.stringtie_raw / f'{group.files_map["bam"].stem}.gtf'
        self.executor.run_function(
            self.stringtie_tool.run_stringtie,
            input_groups, gtf_inputs, override_cpus=self.stringtie_tool.cpus
        )

    def run_stringtie_merge(self, input_groups: List[GroupInfo]) -> None:
        logger.debug('Started stringtie merge.')
        assembly_lst: Path = self.stringtie_merge / 'assembly.lst'
        with open(assembly_lst, 'w') as f:
            print(
                *(group.files_map['gtf'] for group in input_groups),
                sep='\n', file=f
            )
        merged_file: Path = self.stringtie_merge / 'merged.gtf'
        self.stringtie_tool.run_stringtie_merge(assembly_lst, merged_file)
        assembly_lst.unlink()


    def run(self,
            input_groups: List[GroupInfo],
            annot_info: AnnotInfo,
            prefix: str) -> None:
        self.run_stringtie(input_groups, annot_info.gtf_annotation)
        self.run_stringtie_merge(input_groups)