import os
import shutil
import tempfile
from pathlib import Path

import pytest

from symwannier.wannierize import Wannierize


def _preferred_tmp_root() -> Path | None:
    configured = os.environ.get("SYMWANNIER_PYTEST_BASETEMP")
    if configured:
        return Path(configured).expanduser()

    candidate = Path("/home2") / os.environ.get("USER", "") / "tmp" / "pytest"
    if candidate.parent.exists() and os.access(candidate.parent, os.W_OK):
        return candidate

    return None


@pytest.fixture
def test_data_dir():
    """Path to bundled test input files."""
    return Path(__file__).parent / "inputs"


@pytest.fixture
def tmp_path(tmp_path_factory):
    """Create per-test work directories with a local fast-storage override."""
    preferred_root = _preferred_tmp_root()
    if preferred_root is None:
        yield tmp_path_factory.mktemp("test")
        return

    preferred_root.mkdir(parents=True, exist_ok=True)
    path = Path(tempfile.mkdtemp(prefix="pytest-", dir=preferred_root))
    try:
        yield path
    finally:
        shutil.rmtree(path, ignore_errors=True)


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
