# meadowlark_lut — Meadowlark SLM LUT calibration

Calibration workflows for Meadowlark SLM look-up tables, pinned to ModeLab hardware packages (`slm`, `cameras`).

Python package name: **`meadowlark_lut`** (GitHub: [MarKo7s/meadowlark-lut](https://github.com/MarKo7s/meadowlark-lut)).

---

## Camera polarisation-interference LUT

Acquisition (`acquisition/E01_measure_LCOS_response_cam_pol_interference.ipynb`) pistons the LCOS through linear driver levels (typically 256 samples from 0–2047) and records mean camera power in two viewport ROIs (`Ih`, `Iv`; left/right SLM halves). Each ROI sits on a spatial interference fringe. Intensity vs voltage is a cosine of the SLM phase \(\phi(V)\):

\[
I(V)=I_{\mathrm{dc}}+A\cos\bigl(\phi(V)+\psi\bigr)
\]

where \(\psi\) is the local fringe phase at that ROI. Analysis notebooks turn \(I(V)\) into \(\phi(V)\), average repeats, invert to \(V(\phi)\) on a uniform \(0\ldots 2\pi\) grid, and write a 256-entry `.lut`. Shared helpers live in `analysis/utilities.py`.

### 1. Single-channel interference (current default)

Notebook: `analysis/E01_get_LUT_cam_pol_interference.ipynb`.

Uses **one** trace per measurement (typically `Ih` or `Iv`). Assumes \(\phi(V)\) is **monotonic increasing** over the sweep.

1. Optional Savitzky–Golay smooth (`window_length=11`, `polyorder=5`) to knock down camera noise without shifting fringe peaks much.
2. Min–max scale onto \([-1,1]\):

\[
\tilde I = 2\,\frac{I-I_{\min}}{I_{\max}-I_{\min}}-1
\]

3. Folded phase \(\alpha=\arccos(\tilde I)\in[0,\pi]\). As \(\phi\) climbs through a full cycle, \(\alpha\) is a triangle wave \(0\to\pi\to 0\), not a \(2\pi\) wrap.
4. Monotonic unwrap (`unwrap_phase`): integrate **absolute** steps so the downward half of the triangle still counts as increasing \(\phi\):

\[
\phi[0]=0,\qquad
\phi[n]=\sum_{k=1}^{n}\lvert\alpha[k]-\alpha[k-1]\rvert
\]

Over one cosine fringe that sums to \(\approx 2\pi\). `numpy.unwrap` is the wrong tool here: \(\alpha\) never jumps by \(> \pi\).

5. Repeat for every pickle / left–right / H–V trace, then **average** the unwrapped curves and take that mean as \(\phi(V)\). Convert to cycles (\(/\ 2\pi\)), spline-smooth, invert on the \(2\pi\) interval of interest (+ buffer), quantise to integer driver levels → LUT.

**Prior:** SLM phase only increases with voltage, so small noise dips are treated as extra positive phase (slight stroke overestimate) rather than reversals. A constant piston offset is absorbed later via `start_phase_cycle`.

### 2. Quadrature detection

Notebook: `analysis/E01_get_LUT_cam_pol_interference_quadrature_detection.ipynb`.

Uses **both** ROIs on the same fringe field as in-phase / quadrature channels. If the boxes are separated by a quarter of the fringe period \(\Lambda\),

\[
\Delta x=\Lambda/4 \quad\Rightarrow\quad \delta\approx\pi/2
\]

and \((I_x,I_y)\) trace a circle (after scaling). Separation \(\Lambda/2\) is a \(\pi\) shift (\(I\) and \(1-I\)) — a line in the Lissajous plane, not I/Q. Each ROI must be \(\ll\Lambda/4\) or spatial averaging mixes the quadratures.

Model (Heydemann):

\[
\begin{aligned}
I_x &= x_0 + a\cos\phi \\
I_y &= y_0 + b\cos(\phi-\delta),\qquad \varepsilon=\cos\delta .
\end{aligned}
\]

Pipeline (`quadrature_detection`):

1. Fit the implicit ellipse (Fitzgibbon / Halír, constraint \(4AC-B^2=1\))

\[
A I_x^2 + B I_x I_y + C I_y^2 + D I_x + E I_y + F = 0 .
\]

2. Convert coefficients to \((x_0,y_0,a,b,\varepsilon)\). Geometric tilt \(\theta\) and major/minor axes are only for the Lissajous plot; they are not \(\varepsilon\).
3. Centre, scale, and invert:

\[
X=\frac{I_x-x_0}{a},\quad
Y=\frac{I_y-y_0}{b},\quad
\sin\phi=\frac{Y-\varepsilon X}{\sqrt{1-\varepsilon^2}} .
\]

Then \(\phi=\mathrm{atan2}(\sin\phi,\,X)\in(-\pi,\pi]\). If the unwrapped sweep decreases, the \(\sin\phi\) sign is flipped (SLM \(\phi(V)\) should increase).
4. `numpy.unwrap` on this 4-quadrant \(\phi\) (real \(\pm 2\pi\) jumps). Subtract \(\phi[0]\). Average left/right (or repeats) → same LUT inversion as method 1.

When \(\lvert\varepsilon\rvert\to 1\) (needle Lissajous) the \(\sin\phi\) formula divides by a small number and `atan2` is noisy; the recovered \(\phi(V)\) collapses to the single-channel cosine result. Open the ellipse optically (\(\Lambda/4\) ROI offset, or a QWP) so \(\lvert\varepsilon\rvert\) is small; then quadrature is the better estimator (all four quadrants, no \(\arccos\) fold).

### 3. Hilbert transform

Notebook: `analysis/E01_get_LUT_cam_pol_interference_hilbert.ipynb` (placeholder).

Not implemented yet. Intended path: analytic signal of a **single** cosine interferogram to synthesise a quadrature, then `atan2` + unwrap — no second ROI and no assumption of a spatial \(\Lambda/4\).

---

## Installation


### Conda environment (recommended)

```bash
conda env create -f environment.yml
conda activate meadowlarkCalLUT
```

**Protocol:** conda installs **Python + pip only**. Scientific packages and ModeLab hardware pins are installed with **pip** so PySide6 (from `slm` / `cameras`) does not clash with conda Qt DLLs.

This creates **`meadowlarkCalLUT`** (Python 3.11) and pip-installs numpy, scipy, matplotlib plus:

- `slm[notebooks] @ git+https://github.com/MarKo7s/slm.git@v0.1.0`
- `cameras[gui,notebooks] @ git+https://github.com/MarKo7s/cameras.git@v0.2.1`

Allied Vision Vimba X / `vmbpy` is host-side (not a pip extra). After creating the env:

```bash
pip install "C:/Program Files/Allied Vision/Vimba X/api/python/vmbpy-1.2.1-py3-none-win_amd64.whl"
```

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
