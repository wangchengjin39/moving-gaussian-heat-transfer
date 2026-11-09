"""Assistant-authored comparison of archived idealized-model probe exports.

No fitting or experimental validation. Raw exports are kept unchanged.
"""
from pathlib import Path
import hashlib
import json
import re
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
V = ROOT / "validation"
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)


def read_probe(path):
    rows = []
    for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 3:
            try:
                rows.append((int(fields[0]), float(fields[1]), float(fields[2])))
            except ValueError:
                continue
    a = np.asarray(rows, dtype=float)
    assert a.ndim == 2 and a.shape[1] == 3 and len(a) > 1
    assert np.isfinite(a).all() and np.all(np.diff(a[:, 1]) > 0)
    assert np.array_equal(a[:, 0], np.arange(1, len(a) + 1))
    assert abs(a[-1, 1] - 40.0) < 1e-10
    return a[:, 1], a[:, 2]


def mesh_blocks(text):
    lines = text.splitlines()
    blocks = []
    for token in ("nblock,", "eblock,"):
        start = next(i for i, s in enumerate(lines) if s.lower().startswith(token))
        end = next(i for i in range(start + 2, len(lines)) if lines[i].strip() == "-1")
        blocks.append("\n".join(lines[start:end + 1]))
    return blocks


arch_fast = V / "probe_1600W_4mms_40s_h125_two_steps.txt"
slow = V / "probe_1600W_2mms_40s_h125_tableS0p625.txt"
input_slow = V / "moving_gaussian_40s_1600W_h125_tableS0p625_solver_input.dat"
input_fast = V / "moving_1600W_4mms_two_steps_2026-11-05_solver_input.dat"
texts = [p.read_text(encoding="utf-8", errors="replace") for p in (input_slow, input_fast)]
assert mesh_blocks(texts[0]) == mesh_blocks(texts[1])
for text in texts:
    assert re.search(r"^local,12,0,5.e-002,1.e-002,2.8e-002,", text, re.M | re.I)
    for field, value in (("KXX", 20), ("DENS", 7800), ("C", 500)):
        assert re.search(rf"^MP,{field},1,{value},", text, re.M | re.I)
    assert re.search(r"^tunif,20\.", text, re.M | re.I)
for name in ("moving_gaussian_40s_1600W_h125_tableS0p625_probe_CAERep.xml",
             "moving_1600W_4mms_two_steps_2026-11-05_CAERep.xml"):
    root = ET.parse(V / name).getroot()
    probes = [e for e in root.iter() if e.attrib.get("Type") == "ResultProbe"]
    assert any(e.findtext("LocationCoordSysId", "").strip() == "53" for e in probes)

records = []
curves = []
fig, ax = plt.subplots(figsize=(9, 5.4), layout="constrained")
for speed, path, color in ((2, slow, "#2563eb"), (4, arch_fast, "#d97706")):
    t, temp = read_probe(path)
    peak = float(temp.max())
    peak_times = t[temp == peak]
    pass_time = 40 / speed
    heating_end = 80 / speed
    records.append({
        "speed_mm_s": speed, "power_W": 1600, "record_count": len(t),
        "first_time_s": float(t[0]), "heating_end_s": heating_end,
        "nominal_input_J": 1600 * heating_end, "source_pass_time_s": pass_time,
        "peak_temperature_C": peak, "peak_rise_K": peak - 20,
        "first_highest_sample_time_s": float(peak_times[0]),
        "last_highest_sample_time_s": float(peak_times[-1]),
        "highest_sample_count": len(peak_times),
        "first_sample_peak_lag_s": float(peak_times[0] - pass_time),
        "temperature_40s_C": float(temp[-1]),
        "raw_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    })
    curves.append(pd.DataFrame({"speed_mm_s": speed, "time_s": t, "temperature_C": temp}))
    ax.plot(t, temp, color=color, lw=1.8, label=f"{speed} mm/s")
    ax.scatter(peak_times, temp[temp == peak], s=22, color=color, zorder=4)
    ax.axvline(pass_time, color=color, ls=":", lw=1, alpha=.7)
ax.axvline(20, color="0.35", ls="--", lw=1, label="4 mm/s heating ends (20 s)")
ax.set(xlabel="Time (s)", ylabel="Probe temperature (deg C)", xlim=(0, 40),
       title="Same probe, 1600 W: speed comparison\nConstant properties; insulated plate; numerical model")
ax.grid(alpha=.25)
ax.legend()
fig.savefig(OUT / "probe_1600W_speed_comparison.png", dpi=180)
plt.close(fig)
pd.DataFrame(records).to_csv(OUT / "probe_1600W_speed_comparison_summary.csv", index=False)
pd.concat(curves, ignore_index=True).to_csv(OUT / "probe_1600W_speed_comparison_curves.csv", index=False)
a, b = records
report = {
    "cases": records, "same_mesh_and_constant_materials_checked": True,
    "same_probe_m": [.05, .01, .028],
    "fast_minus_slow_peak_K": b["peak_temperature_C"] - a["peak_temperature_C"],
    "fast_peak_rise_reduction_percent": 100 * (1 - b["peak_rise_K"] / a["peak_rise_K"]),
    "fast_minus_slow_first_peak_sample_time_s": b["first_highest_sample_time_s"] - a["first_highest_sample_time_s"],
    "fast_minus_slow_final_temperature_K": b["temperature_40s_C"] - a["temperature_40s_C"],
    "scope": "Two idealized numerical cases, not experimental validation or proof of full-field convergence.",
    "peak_time_caution": "Highest sampled times refer to rounded exports, not exact continuous peak times or guaranteed brackets.",
    "energy_caution": "Nominal P*time; continuous finite-plate Gaussian tails cause a small input truncation. Body energy is checked separately by validate_cooling_energy.py.",
    "authorship": "Comparison script by assistant; not an independent student programming assignment.",
}
(OUT / "speed_comparison_check.json").write_text(
    json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(report, ensure_ascii=False, indent=2))
