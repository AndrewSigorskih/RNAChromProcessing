import logging
import subprocess
from pathlib import Path

from typing import Any, Iterable, List, Union

DATE_FORMAT = '%Y-%m-%d %H:%M:%S'
LOGGING_FORMAT = '%(asctime)s | %(levelname)-8s | %(module)s | %(message)s'
VERBOSE = 5

logger = logging.getLogger()
PathLike = Union[str, Path]


def configure_logger(logger: logging.Logger, verbose: int = 0) -> None:
    level: int = logging.INFO if not verbose else logging.DEBUG if verbose == 1 else VERBOSE
    logging.addLevelName(VERBOSE, "VERBOSE")
    handler = logging.StreamHandler()
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter(fmt=LOGGING_FORMAT,
                                           datefmt=DATE_FORMAT))
    logger.setLevel(level)
    logger.addHandler(handler)


def find_in_list(id: str, lst: List[str]):
    """A simple helper to find file from
    list of files regardless of extension"""
    try:
        return next(x for x in lst if x.startswith(id))
    except StopIteration:
        return None


def run_command(cmd: Union[List[str], str],
                **subprocess_args: Any) -> int:
    log_cmd(cmd)
    return_code: int = subprocess.run(cmd, **subprocess_args).returncode
    return return_code


def run_get_stdout(cmd: Union[List[str], str],
                   **subprocess_args: Any) -> str:
    log_cmd(cmd)
    result = subprocess.run(cmd, capture_output=True,
                            text=True, **subprocess_args)
    return result.stdout


def log_cmd(cmd: Union[Iterable[PathLike], str]) -> None:
    if isinstance(cmd, list):
        cmd_str = " ".join(str(x) for x in cmd)
    else:
        cmd_str = cmd
    logger.log(VERBOSE, f'Running command: {cmd_str}')
