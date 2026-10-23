import logging
import os
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Optional, Set

from pydantic import BaseModel, PositiveInt, field_validator

from .AnnotInfo import AnnotInfo
from .StringtiePipeline import StringtieTool, StringtiePipeline
from .XBamFilter import XBamFilter
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
    keep_intermediate: bool = False

    stringtie: StringtieTool

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # working dir
        if not self.base_dir.exists():
            self.base_dir.mkdir(parents=True)
        os.chdir(self.base_dir)
        self._work_dir = TemporaryDirectory(dir=self.base_dir, delete=True)
        # other members
        work_pth = Path(self._work_dir.name)
        executor = PoolExecutor(self.cpus)
        self.annotation.prepare_annotation(work_pth)
        self._preprocessing = XBamFilter(work_pth, executor)
        self._pipeline = StringtiePipeline(work_pth, executor, self.stringtie)

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
        groups = self._preprocessing.run(
            self.file_ids,
            self.annotation,
            self.bam_input_dir, 
            self.filter_ids_input_dir
        )
        self._pipeline.run(
            groups,
            self.annotation,
            self.ouputs_prefix
        )
        # save outputs
        collect_outputs(
            Path(self._work_dir.name),
            self.output_dir,
            self.ouputs_prefix,
            self.keep_intermediate,
        )
        logger.info('Done.')


def _merge_tree(src: Path, dst: Path) -> None:
    """Copy src into dst, overwriting only files whose relative paths match."""
    dst.mkdir(parents=True, exist_ok=True)
    for item in src.rglob("*"):
        target = dst / item.relative_to(src)
        if item.is_dir():
            target.mkdir(parents=True, exist_ok=True)
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(item, target)


def collect_outputs(
        workdir: Path,
        output_dir: Path,
        output_prefix: str,
        keep_intermediate: bool = False,
    ) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    for suffix in ("bed", "tab"):
        destination = output_dir / f"{output_prefix}.{suffix}"
        shutil.copyfile(workdir / "xrna" / f"xrna.{suffix}", destination)

    for subdir in workdir.iterdir():
        if not subdir.is_dir() or subdir.name == "xrna":
            continue
        if keep_intermediate or subdir.name == "plots":
            _merge_tree(subdir, output_dir / subdir.name)
