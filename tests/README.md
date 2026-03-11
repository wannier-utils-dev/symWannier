# Running tests

This repository uses a `src/` layout. Running `python3 -m pytest` directly may fail with
`ModuleNotFoundError: symwannier` unless the package is installed or `src` is added to
`PYTHONPATH`.

## Quick run

From the repository root:

```bash
PYTHONPATH=src python3 -m pytest -q
```

This is the simplest way to run the full test suite without creating a virtual environment.

## Editable install

If you want to run the tests in an isolated environment:

```bash
python3 -m venv .venv
.venv/bin/pip install -e '.[dev]'
.venv/bin/python -m pytest -q
```

## Run a single test file

```bash
PYTHONPATH=src python3 -m pytest -q tests/test_graphene.py
```
