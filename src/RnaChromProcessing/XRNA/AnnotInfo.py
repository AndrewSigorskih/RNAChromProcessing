import logging
import os
from dataclasses import dataclass, Field
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
from pydantic import BaseModel, field_validator

from ..utils import exit_with_error, run_command
from ..utils.parsing import read_bedrc, read_gtf, write_bed

logger = logging.getLogger()


@dataclass
class SampleInfo:
    sample_id: str
    sample_group: str
    true_strand: bool
    files_map: Dict[str, Path] = Field(default_factory=dict)

    # TODO remove it ALL
    fq_file: Optional[Path] = None
    bed_file: Optional[Path] = None
    lst_file: Optional[Path] = None
    bam_file: Optional[Path] = None

    def update_field(self,
                     field: str,
                     val: Path) -> None:
        setattr(self, field, val)
    
    def inputs_ok(self) -> bool:
        return (self.bed_file is not None) and (self.fq_file is not None)
    

class GroupInfo:
    group_id: str
    samples: List[SampleInfo]
    true_strand: bool
    files_map: Dict[str, Path] = Field(default_factory=dict)


class AnnotInfo(BaseModel):
    gene_annotation: Path
    strand_info: Path
    samples: List[SampleInfo] = []

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._read_strand_info()

    @field_validator('gene_annotation', 'strand_info')
    @classmethod
    def check_file(cls, pth: str) -> Path:
        assert \
            (os.path.exists(pth) and (os.path.getsize(pth) > 0)),\
            f'File {pth} does not exist or is empty!'
        return Path(pth).resolve()

    def _read_strand_info(self) -> None:
        logger.debug("Started processing strand information file.")
        strand_info = pd.read_table(self.strand_info, sep='\t', index_col=[0,1]).sort_index()

        # identify and drop any ids for which strand is unknown
        unk_strand = strand_info[strand_info.strand == 'UNKNOWN'].index
        if (unk_num := len(unk_strand)) > 0:
            logger.warning(
                f'Will ignore {unk_num} files that failed strand detection: '
                f'{", ".join(x for x in unk_strand.get_level_values(1))}'
            )
            strand_info = strand_info.drop(index=unk_strand)

        # identify and drop any group with inconsistent strand status
        to_drop = []
        for group in strand_info.index.unique(level=0):
            if strand_info.loc[group,].strand.nunique() > 1:
                to_drop.append(group)
        if to_drop:
            logger.warning(
                "These groups have inconsistent strand status and will be ignored: "
                f"{', '.join(to_drop)}"
            )
            strand_info = strand_info.drop(index=to_drop)

        if not strand_info.shape[0]:
            exit_with_error("No suitable entries to process!")

        self.samples = [
            *(
                SampleInfo(sample_id, sample_group, True)
                for sample_group, sample_id in
                strand_info[strand_info.strand == 'SAME'].index.to_flat_index()
            ),
            *(
                SampleInfo(sample_id, sample_group, False)
                for sample_group, sample_id in
                strand_info[strand_info.strand == 'ANTI'].index.to_flat_index()
            )
        ]

    def prepare_annotation(self,
                           work_pth: Path) -> None:
        self._annot_bed: Path = work_pth / 'gene_annotation.bed'
        tmp_bed: Path = work_pth / 'tmp_annotation.bed'
        
        # convert annotation to bed
        file_ext = self.gene_annotation.suffix
        if file_ext == '.gtf':
            genes = read_gtf(self.gene_annotation)
        elif file_ext == '.bedrc':
            genes = read_bedrc(self.gene_annotation)
        else:
            exit_with_error(
                f"Unknown gene annotation format: {file_ext}! "
                "Supported formats: gtf, bedrc."
            )
        genes['score'] = 1
        write_bed(genes, tmp_bed)

        # sort annot bed file
        cmd = f'bedtools sort -i {tmp_bed} > {self._annot_bed}'
        return_code = run_command(cmd, shell=True)
        if return_code:
            exit_with_error(f'bedtools sort returned non-zero exit code: {return_code}')
        tmp_bed.unlink()
    
    @property
    def annot_bed(self) -> Path:
        if not hasattr(self, '_annot_bed'):
            msg = (
                f'{self.__class__.__name__} was not properly instantiated! '
                'Run prepare_annotation() method before requesting bed annotation.'
            )
            exit_with_error(msg)
        return self._annot_bed
