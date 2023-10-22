import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd
from pydantic import BaseModel

from ..utils import check_file_exists, exit_with_error, run_command

logger = logging.getLogger()


@dataclass
class SampleInfo:
    sample_id: str
    sample_group: str
    true_strand: bool
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


class AnnotInfo(BaseModel):
    gtf_annotation: Path
    strand_info: Path

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        check_file_exists(self.gtf_annotation)
        check_file_exists(self.strand_info)
        self._read_strand_info()


    def _read_strand_info(self) -> None:
        strand_info = pd.read_table(self.strand_info, sep='\t',
                                    index_col=[0,1])
        unk_strand = strand_info[strand_info.strand == 'UNKNOWN'].index.get_level_values(1)
        if (unk_num := len(unk_strand)) > 0:
            logger.debug(f'Will ignore {unk_num} files that failed strand detection:')
            logger.debug(f'{", ".join(x for x in unk_strand)}')

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
        # gtf to bed conversion
        genes = pd.read_csv(self.gtf_annotation, skiprows=5, sep='\t', 
                            header=None, usecols=[0,2,3,4,6,8])
        genes = genes[genes[2] == 'gene']
        genes[8] = genes[8].apply(lambda x: x.split(";")[0].split(" ")[-1].strip('""'))
        genes[5] = 1
        genes[[0, 3, 4, 8, 5, 6]].to_csv(tmp_bed, sep='\t', index=False, header=False)

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
