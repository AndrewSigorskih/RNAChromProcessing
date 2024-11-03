import gzip
import shutil
import os
from pathlib import Path
from typing import Optional

from .errors import exit_with_error


def check_file_exists(path: str) -> None:
    if not (os.path.exists(path) and (os.path.getsize(path) > 0)):
        message = f'File {path} does not exist or is empty!'
        exit_with_error(message)


def gzip_file(src: str, dst: str, remove_src=True) -> None:
    """Gzip file to specified name"""
    with open(src, 'rb') as infile, gzip.open(dst, 'wb') as outfile:
        shutil.copyfileobj(infile, outfile)
    if remove_src:
        os.remove(src)


def make_directory(path: str, exist_ok: bool = True) -> None:
    try:
        os.makedirs(path, exist_ok=exist_ok)
    except OSError as e:
        exit_with_error(f'Could not create directory: {e}')


def move_exist_ok(src_dir: str, dest_dir: str):
    """Move directory src to dest even if dest_dir
    and some of its contents exist"""
    fileList = os.listdir(src_dir)
    for i in fileList:
        src = os.path.join(src_dir, i)
        dest = os.path.join(dest_dir, i)
        if os.path.exists(dest):
            if os.path.isdir(dest):
                move_exist_ok(src, dest)
                continue
            else:
                os.remove(dest)
        shutil.move(src, dest_dir)
    

def remove_suffixes(pth: Path) -> str:
    '''remove all suffixes from Path object'''
    return str(pth).removesuffix(''.join(pth.suffixes))


def validate_tool_path(path: Optional[str], tool_name: str) -> Path:
    """A helper function used in custom validators of tool configs"""
    path = path or shutil.which(tool_name)
    if not path:
        exit_with_error(f'Cannot deduce path to {tool_name} executable!')
    check_file_exists(path)
    return Path(path)
