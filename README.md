# meadowlark_lut — Meadowlark SLM LUT calibration

Calibration workflows for Meadowlark SLM look-up tables, pinned to ModeLab hardware packages (`slm`, `cameras`).

Python package name: **`meadowlark_lut`** (GitHub: [MarKo7s/meadowlark-lut](https://github.com/MarKo7s/meadowlark-lut)).

---

## Installation

### Conda environment (recommended)

```bash
conda env create -f environment.yml
conda activate meadowlarkCalLUT
```

This creates **`meadowlarkCalLUT`** (Python 3.11) with numpy, scipy, matplotlib from conda, and installs:

- `slm[notebooks] @ git+https://github.com/MarKo7s/slm.git@v0.1.0`
- `cameras[gui,notebooks] @ git+https://github.com/MarKo7s/cameras.git@v0.2.1`

Allied Vision Vimba X / `vmbpy` is host-side (see the [cameras](https://github.com/MarKo7s/cameras) README); it is not a pip extra.

### From GitHub (tagged release)

After a release tag exists:

```bash
pip install "meadowlark_lut @ git+https://github.com/MarKo7s/meadowlark-lut.git@v0.1.0"
```

That install pulls the same pinned `slm` / `cameras` GitHub tags via `pyproject.toml` dependencies.

### Local development (editable install)

```bash
git clone git@github.com:MarKo7s/meadowlark-lut.git
cd meadowlark-lut
conda activate meadowlarkCalLUT
pip install -e .
```

`requirements.txt` remains available as a legacy mirror of the hardware GitHub pins; prefer `environment.yml` or `pip install -e .`.

---

## Versioning

The package version is defined in **one place only**: `pyproject.toml` → `[project].version`.

Do **not** edit `__init__.py` on each release. `meadowlark_lut.__version__` is read from pip metadata after install (`importlib.metadata`).

```bash
pip show meadowlark_lut
python -c "import meadowlark_lut; print(meadowlark_lut.__version__)"
```

Use [semantic versioning](https://semver.org/): `MAJOR.MINOR.PATCH`.

Hardware compatibility is frozen by pinning `slm` and `cameras` to git tags in `pyproject.toml` / `environment.yml`. Bump those pins deliberately when you intentionally adopt newer hardware APIs.

## Releasing a new version

1. Add an entry for the new version at the top of `CHANGELOG.md`.
2. Bump `version` in `pyproject.toml`.
3. Commit all changes (including the changelog).
4. Run the release script from the repo root:

```bash
python scripts/release.py --from-changelog
```

The script reads the version from `pyproject.toml`, pushes `main`, creates annotated git tag `vX.Y.Z`, and pushes the tag. With `--from-changelog`, the tag message is taken from the matching `CHANGELOG.md` section.

Install a released tag:

```bash
pip install "meadowlark_lut @ git+https://github.com/MarKo7s/meadowlark-lut.git@vX.Y.Z"
```

Dry run:

```bash
python scripts/release.py --from-changelog --dry-run
```

Optional GitHub Release:

```bash
gh release create vX.Y.Z --title "meadowlark_lut X.Y.Z" --notes-file CHANGELOG.md
```

**Requirements before release:** clean working tree; tag `vX.Y.Z` must not already exist on GitHub.
