from shutil import rmtree
from logging import getLogger
from os import chdir, listdir
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import List

from pydantic import BaseModel, PositiveInt

from .stages import SamplePair
from ..utils import exit_with_error, find_in_list

logger = getLogger()


class BaseProcessor(BaseModel):
    input_dir: Path
    output_dir: Path
    base_dir: Path = Path('.').resolve()

    cpus: PositiveInt = 1
    
    rna_ids: List[str]
    dna_ids: List[str]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # check for errors
        self.__validate_inputs()
        # work dirs
        self.base_dir.mkdir(parents=True, exist_ok=True)
        chdir(self.base_dir)
        self._work_dir = TemporaryDirectory(dir=self.base_dir)
        self._work_pth = Path(self._work_dir.name).resolve()

    def __del__(self):
        rmtree(self._work_pth)

    def __validate_inputs(self):  # TODO probably rewrite validations in a more pydantic-ish way
        self.input_dir = self.input_dir.resolve()
        self.output_dir = self.output_dir.resolve()
        if not self.input_dir.exists():
            exit_with_error(f'Input directory doesnt exist: {self.input_dir}')
        if not self.output_dir.exists():
            self.output_dir.mkdir(parents=True)
        if len(self.rna_ids) != len(self.dna_ids):
            exit_with_error('DNA and RNA id arrays have different length!')

    def gather_inputs(self) -> List[SamplePair]:
        samples = []
        filenames = listdir(self.input_dir)
        for rna_id, dna_id in zip(self.rna_ids, self.dna_ids):
            rna_file = find_in_list(rna_id, filenames)
            dna_file = find_in_list(dna_id, filenames)
            if not rna_file or not dna_file:
                logger.warning(
                    f'Could not find one or both files {rna_id} {dna_id} '
                    'in input directory, skipping pair!'
                )
                continue
            samples.append(
                SamplePair(
                    self.input_dir / dna_file,
                    self.input_dir / rna_file
                )
            )
        if not samples:
            exit_with_error('Could not find any files in input directory, aborting!')
        return samples
