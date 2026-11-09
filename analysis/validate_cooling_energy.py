"""Assistant-authored body-energy check on archived affine SOLID279 exports.

Reads relative archived inputs; no solver invocation or source-file modification.
The integration method follows validate_moving_gaussian_h125_energy.py, with
per-value rounding bounds for Mechanical's significant-digit text exports.
"""
from pathlib import Path
import hashlib
import json
import math
import re

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
V = ROOT / 'validation'
MESH = V / 'moving_1600W_4mms_two_steps_2026-11-05_solver_input.dat'
LOG = V / 'moving_1600W_4mms_two_steps_2026-11-05_file.nlh'
SOURCES = [(20., V / 'temperature_1600W_4mms_20s_h125.txt'),
           (40., V / 'temperature_1600W_4mms_40s_h125.txt')]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rounding_half_unit(token):
    # Half the last printed decimal unit; assumes nearest rounding.
    mantissa, *exponent = token.lower().split('e')
    exponent = int(exponent[0]) if exponent else 0
    decimals = len(mantissa.split('.')[1]) if '.' in mantissa else 0
    return .5 * 10. ** (exponent - decimals)


lines = MESH.read_text(encoding='utf-8', errors='replace').splitlines()
coords = {}
i = next(k for k, s in enumerate(lines) if s.lower().startswith('nblock,')) + 2
while not lines[i].strip().startswith('-1'):
    fields = lines[i].split()
    coords[int(fields[0])] = list(map(float, fields[1:4]))
    i += 1
elements = []
i = next(k for k, s in enumerate(lines) if s.lower().startswith('eblock,')) + 2
while not lines[i].strip().startswith('-1'):
    fields = list(map(int, lines[i].split()))
    assert len(fields) == 21
    elements.append(fields[1:])
    i += 1
assert (len(coords), len(elements)) == (92529, 20480)
assert any(s.lower() == 'et,1,279' for s in lines)
text = '\n'.join(lines)
for field, value in [('KXX', 20), ('DENS', 7800), ('C', 500)]:
    assert re.search(rf'^MP,{field},1,{value},', text, re.M | re.I)
assert re.search(r'^tunif,20\.', text, re.M | re.I)

node_ids = np.array(sorted(coords), dtype=int)
node_xyz = np.array([coords[n] for n in node_ids])
lookup = {n: k for k, n in enumerate(node_ids)}
conn = np.array([[lookup[n] for n in e] for e in elements])
xyz = node_xyz[conn]
corners = xyz[:, :8]
origin = corners[:, 0]
u = corners[:, 1] - origin
v = corners[:, 3] - origin
w = corners[:, 4] - origin
expected = np.stack([origin, origin+u, origin+u+v, origin+v,
                     origin+w, origin+u+w, origin+u+v+w, origin+v+w], axis=1)
assert np.allclose(corners, expected, rtol=0, atol=1e-10)
edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),
         (0,4),(1,5),(2,6),(3,7)]
for k, (a, b) in enumerate(edges):
    assert np.allclose(xyz[:, 8+k], (corners[:, a]+corners[:, b])/2,
                       rtol=0, atol=1e-10)
volumes = np.abs(np.einsum('ij,ij->i', u, np.cross(v, w)))
assert np.all(volumes > 0)
volume = float(volumes.sum())
assert np.isclose(volume, 4e-5, rtol=1e-10, atol=1e-15)

# Exact mean weights for affine twenty-node serendipity bricks.
weights = np.array([-1/8]*8 + [1/6]*12)
signs = np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
                  [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]])
points, gw = np.polynomial.legendre.leggauss(3)
integrated = np.zeros(20)
for a, wa in zip(points, gw):
    for b, wb in zip(points, gw):
        for c, wc in zip(points, gw):
            p = np.array([a, b, c])
            vals = list(np.prod(1+signs*p, axis=1) * (signs@p-2)/8)
            for first, second in edges:
                s1, s2 = signs[first], signs[second]
                varying = np.flatnonzero(s1 != s2)[0]
                fixed = [j for j in range(3) if j != varying]
                vals.append((1-p[varying]**2) * np.prod(1+s1[fixed]*p[fixed])/4)
            integrated += np.array(vals)*wa*wb*wc/8
assert np.allclose(integrated, weights, rtol=0, atol=1e-12)
nodal_volume_weights = np.zeros(len(node_ids))
np.add.at(nodal_volume_weights, conn.ravel(),
          (volumes[:, None]*weights).ravel())
assert np.isclose(nodal_volume_weights.sum(), volume, rtol=1e-10)
assert np.allclose(node_xyz.T@nodal_volume_weights/volume, [.05, .005, .02])

log_rows = []
for line in LOG.read_text().splitlines():
    fields = line.split()
    if len(fields) == 3:
        try:
            log_rows.append(list(map(float, fields)))
        except ValueError:
            pass
log = np.array(log_rows)
rho, cp, initial = 7800., 500., 20.
nominal = 1600.*20.


def finite_plate_power(t):
    xc = .010+.004*t
    return 1600.*.5*(math.erf((.1-xc)/.005)+math.erf(xc/.005))*math.erf(4)


reference, reference_quad_error = quad(finite_plate_power, 0, 20, epsabs=1e-7)
records = []
for time, path in SOURCES:
    rows = [s.split() for s in path.read_text(encoding='latin-1').splitlines()[1:] if s.strip()]
    a = np.array([[float(s) for s in row] for row in rows])
    assert a.shape == (len(coords), 5) and np.isfinite(a).all()
    assert np.array_equal(a[:, 0], node_ids)
    assert np.allclose(a[:, 1:4], node_xyz, rtol=0, atol=1e-10)
    t = a[:, 4]
    eps = np.array([rounding_half_unit(row[4]) for row in rows])
    lr = log[np.isclose(log[:, 0], time, rtol=0, atol=1e-7)]
    assert len(lr) == 1, (time, lr)
    for idx, column in [(np.argmax(t), 1), (np.argmin(t), 2)]:
        assert abs(t[idx]-lr[0, column]) <= eps[idx]+.00051, (time, idx)
    integral = float(nodal_volume_weights@(t-initial))
    independently_summed = float(volumes@((t[conn]-initial)@weights))
    assert abs(integral-independently_summed) < 1e-12
    energy = rho*cp*integral
    bound = rho*cp*float(np.abs(nodal_volume_weights)@eps)
    records.append({
        'time_s': time, 'temperature_export': path.name,
        'temperature_sha256': sha(path), 'nodes': len(coords),
        'temperature_min_C': float(t.min()), 'temperature_max_C': float(t.max()),
        'time_check': 'Export extrema match archived global result extrema at requested time; header omits time.',
        'archived_time_min_C': float(lr[0, 2]), 'archived_time_max_C': float(lr[0, 1]),
        'unweighted_node_mean_C': float(t.mean()),
        'volume_mean_temperature_C': initial+integral/volume,
        'stored_heat_J': energy,
        'export_rounding_energy_bound_J': bound,
        'export_rounding_mean_temperature_bound_K': bound/(rho*cp*volume),
        'difference_from_continuous_finite_plate_reference_J': energy-reference,
        'difference_from_nominal_J': energy-nominal,
    })

drift = records[1]['stored_heat_J']-records[0]['stored_heat_J']
drift_bound = sum(r['export_rounding_energy_bound_J'] for r in records)
report = {
    'authorship': 'Assistant; repetitive integration, not independent user programming work.',
    'input_sha256': sha(MESH), 'global_history_sha256': sha(LOG),
    'nodes': len(coords), 'elements': len(elements), 'element': 'SOLID279 affine twenty-node brick',
    'volume_m3': volume, 'mass_kg': rho*volume, 'heat_capacity_J_K': rho*cp*volume,
    'method': 'Integrate actual twenty-node shape functions over each verified affine element; sum volumes.',
    'nominal_energy_J': nominal, 'nominal_mean_temperature_C': initial+nominal/(rho*cp*volume),
    'continuous_finite_plate_input_reference_J': reference,
    'continuous_reference_mean_temperature_C': initial+reference/(rho*cp*volume),
    'reference_quadrature_error_estimate_J': reference_quad_error,
    'cases': records, 'stored_heat_change_20_to_40s_J': drift,
    'relative_stored_heat_change_percent': 100*drift/records[0]['stored_heat_J'],
    'export_rounding_bound_for_change_J': drift_bound,
    'change_within_export_rounding_bound': abs(drift) <= drift_bound,
    'limitations': [
        'Temperature export lacks time metadata: requested times cross-checked by matching extrema against archived time history.',
        'Rounding bounds assume nearest rounding to the last printed decimal digit; not bounds on solver discretization error.',
        'Finite-plate reference integrates the continuous Gaussian; actual TABLE load interpolation/time integration can differ.',
        'Post-stop energy agreement at these two times does not prove local temperature accuracy or full-field convergence.',
        'Constant properties, no latent heat, and insulated model; no experimental welding-quality validation.',
    ],
    'method_sources': [
        'https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_elem/Hlp_E_SOLID279.html',
        'https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_thry/thy_shp8.html',
    ],
}
out = ROOT / 'outputs' / 'cooling_energy_check.json'
out.parent.mkdir(exist_ok=True)
out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(report, ensure_ascii=False, indent=2))
