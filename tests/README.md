# Running tests

This repository uses a `src/` layout. Running `python3 -m pytest` directly may fail with
`ModuleNotFoundError: symwannier` unless the package is installed or `src` is added to
`PYTHONPATH`.

## Quick run

From the repository root:

```bash
PYTHONPATH=src python3 -m pytest -q -m fast
```

This is the simplest way to run the default regression suite without creating a virtual
environment. It skips the large `Fe SW+PD` end-to-end test marked as `slow`.

## Slow regression

The repository also includes a larger static regression case,
`tests/test_fe_sw_pd.py`, which exercises the `atom_proj + irr_bz + -P -S` path.

Run it explicitly when needed:

```bash
PYTHONPATH=src python3 -m pytest -q tests/test_fe_sw_pd.py -m slow
```

This test takes several minutes because it runs the full disentanglement workflow.

## Editable install

If you want to run the tests in an isolated environment:

```bash
python3 -m venv .venv
.venv/bin/pip install -e '.[dev]'
.venv/bin/python -m pytest -q -m fast
```

## Run a single test file

```bash
PYTHONPATH=src python3 -m pytest -q tests/test_graphene.py
```
