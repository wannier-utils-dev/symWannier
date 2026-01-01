"""Shared I/O utilities for symwannier."""

import os
import gzip
from typing import Tuple, IO


def open_text_or_gz(path: str, desc: str = "file") -> Tuple[IO[str], str]:
    """Open a text file, falling back to gzip if ``.gz`` exists.

    Parameters
    ----------
    path : str
        Base path (without .gz).
    desc : str, optional
        Description used in error message.

    Returns
    -------
    fp : IO[str]
        Opened file-like object in text mode.
    used_path : str
        Actual path used (plain or .gz).

    Raises
    ------
    FileNotFoundError
        If neither the plain file nor the gzipped file exists.
    """
    if os.path.exists(path):
        return open(path, "r"), path
    gz_path = path + ".gz"
    if os.path.exists(gz_path):
        return gzip.open(gz_path, "rt"), gz_path
    raise FileNotFoundError(f"{desc} not found: {path}(.gz)")
