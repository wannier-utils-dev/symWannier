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
            source = test_data_dir / f"{material}.{ext}"
            destination = dest / f"{material}.{ext}"
            if source.exists():
                shutil.copy(source, destination)
            else:
                shutil.copy(Path(f"{source}.gz"), Path(f"{destination}.gz"))
    return _copy


@pytest.fixture
def run_wannier(copy_inputs, tmp_path):
    def _run(
        material: str,
        lsym: bool = True,
        num_iter=None,
        projectability_disentangle: bool = False,
    ):
        copy_inputs(material, tmp_path)
        cwd = os.getcwd()
        try:
            os.chdir(tmp_path)
            wann = Wannierize(
                prefix=material,
                lsym=lsym,
                projectability_disentangle=projectability_disentangle,
            )
            if num_iter is not None:
                wann.win.num_iter = num_iter
            wann.run()
        finally:
            os.chdir(cwd)
        return wann, tmp_path
    return _run
