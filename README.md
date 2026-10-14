# NTN-Nulling

This repository contains simulation tools and notebooks for evaluating interference nulling from terrestrial base stations to Non-Terrestrial Network User Equipments (NTN-UEs). The focus is on analyzing nulling performance under different lambda values, antenna configurations, and detection strategies.

## 📁 Project Structure

### 📊 Nulling CDF Evaluation

- **Lambda_CDF.ipynb**  
  Computes the CDF of nulling gain under different lambda values using perfect channel knowledge (threshold set to the true channel gain).

- **Lambda_CDF_det.ipynb**  
  Computes the CDF of nulling gain using estimated channels (`h_hat`) and thresholds derived from a false alarm probability.

### 📡 Antenna Configuration

- **Large_Antenna_Array.ipynb**  
  Demonstrates a 32×32 antenna array setup for nulling and beamforming evaluation.

### 🛰️ Satellite Geometry and Scene Setup

- **satellite_projection.py**  
  Computes the satellite "look" direction and position for NTN-UEs.

- **SceneConfigSionna.py**  
  Used by `Lambda_CDF.ipynb` to configure the ray tracing scene for nulling evaluations.

- **SceneConfigSionna2.py**  
  Used by `Pos_for_RM.ipynb` to configure scenarios with selected user positions.

- **sionnautils/**  
  Utility functions used by `SceneConfigSionna.py` for generating grids and locating TN/NTN UEs.
  A small acknowledgment: some utilities are adapted from material by Sundeep Rangan.

### 🗺️ Mapping and Visualization

- **map_and_pics.ipynb**  
  Plots the radio map, ray tracing map, and satellite visibility map for the simulated area. 

- **Pos_for_RM.ipynb**  
  Computes user positions and nulling vectors, serving as input to `map_and_pics.ipynb`.(In `SceneConfigSionna2.py` manually picked three TN-UEs for three sectors ), saved TN-UE postions, NTN-UE postions and Beamforming Vectors.

### 🔍 Signal Detection

- **detection_use_time_channel.ipynb**  
  Implements signal detection using time-domain CIR (`cir_to_time`) and analyzes paths and tap delay profiles.

## 📌 Notes

- All simulations assume a known or estimated channel model.
- Nulling vectors are computed using SVD-based techniques.
- Visualizations require matplotlib and ray-tracing output from the Sionna simulator.

## 🔧 Sionna 2.2 environment

The local Conda environment `sionna22` is an independent copy of `sionna20`,
with Sionna and its required dependencies upgraded. Activate it with:

```bash
conda activate sionna22
```

In Jupyter or VS Code, select **Python (sionna22)** and restart the kernel
when switching environments.

| Component | Version selected for `sionna22` |
|---|---|
| Python | 3.11.15 |
| Sionna / Sionna RT | 2.2.0 / 2.2.0 |
| PyTorch | 2.14.1+cu130 |
| Mitsuba / Dr.Jit | 3.9.1 / 1.5.0 |

Sionna 2.2 requires Python >=3.11 and PyTorch >=2.9.1. Python 3.13/3.14
have official wheels for the selected PyTorch, Mitsuba and Dr.Jit versions,
but changing the Python minor version requires reinstalling the original
environment's binary extensions. Python 3.11.15 is retained to preserve the
complete existing toolset. PyTorch 2.14.1 uses CUDA 13.0; the local RTX 5080
Laptop GPU and driver 595.84 support the tested configuration. A newer Python
version alone is not evidence of faster ray tracing.

Sources: [Sionna requirements](https://github.com/NVlabs/sionna/blob/v2.2.0/pyproject.toml),
[RT requirements](https://github.com/NVlabs/sionna-rt/blob/v2.2.0/pyproject.toml),
[PyTorch 2.14.1](https://pypi.org/project/torch/2.14.1/).

To reproduce the environment under a new, unused name, clone the original
environment first, then upgrade the selected packages:

```bash
conda create -n sionna22 --clone sionna20 --copy -y
conda activate sionna22
python -m pip install 'sionna==2.2.0' 'sionna-rt==2.2.0' 'torch==2.14.1'
python scripts/patch_sionna22_inr.py
python -m ipykernel install --user --name sionna22 --display-name 'Python (sionna22)'
python -m pip check
```

## 🔁 Preserve the custom INR metric

The original `radio_map.py` is kept as a record of the project's earlier
customization. **Do not copy that entire file over Sionna 2.2.** It would
replace newer upstream fixes and uses an older internal layout.

Instead, `scripts/patch_sionna22_inr.py` adds the custom metric to the official
2.2.0 files in `sionna/rt/radio_map_solvers/`. It verifies upstream checksums,
backs up each original file as `*.py.pre-ntn-inr`, and refuses to overwrite
unrecognized edits or a different Sionna RT version. Applying it again is safe.

The original formula is preserved:

```python
# All terrestrial base stations contribute interference at the NTN receiver.
# Linear, unitless output with shape [1, ...].
inr_ntn = sum(received_power_from_each_base_station) / thermal_noise_power

values = rm.inr_ntn.numpy()          # rm.inr is an alias
rm.show(metric="inr_ntn")           # planar maps, displayed in dB
rm.cdf(metric="inr_ntn")
scene.preview(radio_map=rm, rm_metric="inr_ntn")
# scene.render(..., radio_map=rm, rm_metric="inr_ntn") also works.
```

Both `"inr_ntn"` and `"inr"` work as metric names. INR is a sum across all
transmitters, so selecting a valid `tx` for display still returns that same
aggregate map. INR has no serving-transmitter association: use `"rss"` for
strongest-source association, or `tx_association=False` when sampling cells
using INR thresholds. INR uses `10*log10(INR)` in dB, without the RSS-to-dBm offset.

After reinstalling Sionna RT 2.2.0, reapply and verify the patch, then restart
the notebook kernel:

```bash
python scripts/patch_sionna22_inr.py
python scripts/patch_sionna22_inr.py --check
MPLBACKEND=Agg python -m unittest discover -s tests -p test_sionna22_inr.py -v
```

The tests cover the aggregate formula, planar/mesh shapes, aliases, dB plots,
CDF edge cases, cell sampling, scene rendering/preview, PyTorch CUDA operations,
and a small GPU path/radio-map simulation. They do not run the full research
notebooks. Local package snapshots, original sources, and validation output
are stored in `.local/sionna22-migration/` (excluded from Git).

