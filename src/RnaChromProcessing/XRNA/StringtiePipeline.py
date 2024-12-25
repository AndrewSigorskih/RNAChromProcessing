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


    def run(self,
            input_bams: List[GroupInfo],
            annot_info: AnnotInfo,
            prefix: str) -> None:
        pass