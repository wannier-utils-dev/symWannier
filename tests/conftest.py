import os
import shutil
from pathlib import Path

import pytest

from symwannier.wannierize import Wannierize


@pytest.fixture
def test_data_dir():
    """Path to bundled test input files."""
    return Path(__file__).parent / "inputs"


@pytest.fixture
def copy_inputs(test_data_dir):
    def _copy(material: str, dest: Path) -> None:
        for ext in ["nnkp", "isym", "iamn", "immn", "ieig", "win"]:
            plain = test_data_dir / f"{material}.{ext}"
            gz = test_data_dir / f"{material}.{ext}.gz"
            if plain.exists():
                shutil.copy(plain, dest / plain.name)
            elif gz.exists():
                shutil.copy(gz, dest / gz.name)
            else:
                raise FileNotFoundError(f"missing input: {plain}(.gz)")

    return _copy


@pytest.fixture
def run_wannier(copy_inputs, tmp_path):
    def _run(
        material: str, lsym: bool = True, num_iter=None, win_overrides=None, **kwargs
    ):
        copy_inputs(material, tmp_path)
        cwd = os.getcwd()
        try:
            os.chdir(tmp_path)
            wann = Wannierize(prefix=material, lsym=lsym, **kwargs)
            if num_iter is not None:
                wann.win.num_iter = num_iter
            for key, value in (win_overrides or {}).items():
                setattr(wann.win, key, value)
            wann.run()
        finally:
            os.chdir(cwd)
        return wann, tmp_path

    return _run
