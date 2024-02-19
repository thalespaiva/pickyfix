# WeightedFix

This is an attempt to create an improved decoder for BIKE, called WeightedFix, based on PickyFix.

## Commands

Install Python dependencies and launch virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
pip install statsmodels ipython pandas tqdm seaborn
```

Run mini DFR experiment:

```bash
cd analysis
python dfr_experiment.py ../bike ../data/setup/dfr_experiment_weightedfix_simple.csv reproduced/tmp_dfr || rm -rf reproduced/tmp_dfr
```

Cleanup:

```bash
cd ../bike && make clean && cd -
rm -rf reproduced/tmp_dir
```