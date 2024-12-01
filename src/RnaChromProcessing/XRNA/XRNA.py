import logging
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Optional, Set

from pydantic import BaseModel, PositiveInt, field_validator

from .AnnotInfo import AnnotInfo
from .DataPreprocessing import HisatTool, PreprocessingPipeline
from .StringtiePipeline import StringtieTool, StringtiePipeline
from ..utils.PoolExecutor import PoolExecutor

logger = logging.getLogger()
logging.getLogger("matplotlib").setLevel(logging.WARNING)


class XRNAProcessor(BaseModel):
    file_ids: Set[str]
    bam_input_dir: Path
    filter_ids_input_dir: Optional[Path] = None
    
    output_dir: Path
    base_dir: Path = Path('.').resolve()
    cpus: PositiveInt = 1

    annotation: AnnotInfo
    ouputs_prefix: str = 'xrna'
    keep_extras: Set[str] = set() # TODO maybe remove

    hisat: HisatTool # TODO remove
    stringtie: StringtieTool

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # working dir
        if not self.base_dir.exists():
            self.base_dir.mkdir(parents=True)
        os.chdir(self.base_dir)
        self._work_dir = TemporaryDirectory(dir=self.base_dir)
        # outout dir
        if not self.output_dir.exists():
            self.output_dir.mkdir(parents=True)
        # other members
        work_pth = Path(self._work_dir.name)
        executor = PoolExecutor(self.cpus)
        self.annotation.prepare_annotation(work_pth)
        self._preprocessing = PreprocessingPipeline(
            work_pth, executor, self.hisat
        )
        self._pipeline = StringtiePipeline(
            work_pth, executor, self.stringtie
        )

    @field_validator('file_ids', mode='after')
    @classmethod
    def check_length(cls, ids: Set[str]):
        assert len(ids) > 0, "No input file ids were provided!"
        return ids
    
    @field_validator('bam_input_dir', 'filter_ids_input_dir', mode='after')
    @classmethod
    def check_inp_dir_exists(cls, path: Optional[Path]):
        if path:
            assert path.exists() and path.is_dir(), f"{path} should be an existing directory!"
        return path

    def run(self):
        # run pipeline
        prepared_bams = self._preprocessing.run(
            self.rna_ids, self.bed_input_dir, self.fq_input_dir, self.annotation
        )
        self._pipeline.run(
            prepared_bams, self.annotation, self.ouputs_prefix
        )
        # save outputs
        self._preprocessing.save_outputs(self.output_dir, self.keep_extras)
        self._pipeline.save_outputs(self.ouputs_prefix, self.output_dir, self.keep_extras)
        logger.info('Done.')
