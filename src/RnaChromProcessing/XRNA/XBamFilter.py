from pathlib import Path

from .AnnotInfo import AnnotInfo, SampleInfo
from ..utils.PoolExecutor import PoolExecutor



class XBamFilter:
    def __init__(self,
                work_pth: Path,
                executor: PoolExecutor):
        pass