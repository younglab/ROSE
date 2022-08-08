#!/usr/bin/env python

import errno
import pathlib
import os

from pathlib import Path


def check_file(
    filename: str
) -> str:
    """Checks if file exists

    Args:
        filename (str): Complete path to file
    Raises:
        FileNotFoundError: If file does not exist
    Returns:
        str: Original file name
    """

    if not os.path.isfile(filename):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filename)

    return filename


def check_path(
    filename: pathlib.PosixPath
) -> str:
    """Checks the path of a given file and creates new directories if necessary

    Args:
        filename (pathlib.PosixPath): Complete path to file
    Returns:
        str: Original file name
    """

    if not filename.is_dir():
        filename.mkdir()

    return str(filename)


def get_path() -> str:
    """Get project path

    Returns:
        str: Full path to the project
    """

    idx = str(Path.cwd()).split("/").index("ROSE")
    path = "/".join(str(Path.cwd()).split("/")[:idx+1])

    return path