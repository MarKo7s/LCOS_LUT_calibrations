# Changelog

All notable changes to meadowlark_lut are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.1.0] - 2026-08-06

### Added

- ModeLab pip/git-tag packaging (`pyproject.toml`, `CHANGELOG.md`, `scripts/release.py`).
- Conda environment `meadowlarkCalLUT` (`environment.yml`): conda provides Python + pip only; numpy/scipy/matplotlib and hardware pins install via pip (avoids conda Qt vs PySide6 on Windows).
- Pinned hardware deps: `slm[notebooks]@v0.1.0`, `cameras[gui,notebooks]@v0.2.1`.
