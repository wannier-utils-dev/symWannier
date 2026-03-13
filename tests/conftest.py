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


def pytest_collection_modifyitems(items):
    """Mark all non-slow tests as fast for simpler selection."""
    fast = pytest.mark.fast
    for item in items:
        if item.get_closest_marker("slow") is None:
            item.add_marker(fast)


@pytest.fixture
def test_data_dir():
    """Path to bundled test input files."""
    return Path(__file__).parent / "inputs"


@pytest.fixture
def work_dir(tmp_path_factory):
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
    def _copy(case_name: str, dest: Path) -> None:
        for ext in ["nnkp", "isym", "iamn", "immn", "ieig", "win"]:
            plain = test_data_dir / f"{case_name}.{ext}"
            gz = test_data_dir / f"{case_name}.{ext}.gz"
            if plain.exists():
                shutil.copy(plain, dest / plain.name)
            elif gz.exists():
                shutil.copy(gz, dest / gz.name)
            else:
                raise FileNotFoundError(f"missing input: {plain}(.gz)")

    return _copy


@pytest.fixture
def run_wannier(copy_inputs, work_dir):
    def _run(
        case_name: str, lsym: bool = True, num_iter=None, win_overrides=None, **kwargs
    ):
        copy_inputs(case_name, work_dir)
        cwd = os.getcwd()
        try:
            os.chdir(work_dir)
            wann = Wannierize(prefix=case_name, lsym=lsym, **kwargs)
            if num_iter is not None:
                wann.win.num_iter = num_iter
            for key, value in (win_overrides or {}).items():
                setattr(wann.win, key, value)
            wann.run()
        finally:
            os.chdir(cwd)
        return wann, work_dir

    return _run
