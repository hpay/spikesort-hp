"""Run Kilosort 4 on one raw binary. Called from runKilosort4.m (MATLAB) with the kilosort4
conda environment's python.

    python runKilosort4.py <settings.json>

settings.json: data_file, probe_file (Kilosort2 .mat channel map), results_dir,
bad_channels (raw CH#, 1-indexed), and ks4_settings (passed to run_kilosort as settings;
must include n_chan_bin and fs).
"""
import json
import sys
from pathlib import Path

from kilosort import io, run_kilosort

with open(sys.argv[1]) as f:
    cfg = json.load(f)

probe = io.load_probe(cfg['probe_file'])
bad = cfg.get('bad_channels', [])
bad_rows = [c - 1 for c in (bad if isinstance(bad, list) else [bad])]   # rows in the binary file, 0-indexed
Path(cfg['results_dir']).mkdir(parents=True, exist_ok=True)

run_kilosort(settings=cfg['ks4_settings'], probe=probe, filename=cfg['data_file'],
             results_dir=cfg['results_dir'], bad_channels=bad_rows or None)
