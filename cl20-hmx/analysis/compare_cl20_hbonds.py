#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 6: initial-origin intermolecular C-H...O contact comparison.

Requires Python >= 3.8, NumPy and Matplotlib.
Usage:
  python3 compare_cl20_hbonds.py --root /path/to/trajectories --check-only
  python3 -u compare_cl20_hbonds.py --root /path/to/trajectories
  python3 compare_cl20_hbonds.py --plot-only /path/to/comparison_mean_sd.csv

Layout under --root (default: current directory):
  pcl20/{1,2,3}-100ps-1500k/
  ycl20/{1,2,3}-100ps-1500k/
  md/1-100-ps-1500k/, md/5-100ps-1500k/, md/6-100ps-1500k/
Each run contains input.lammps and dump.lammpstrj. The input's read_data
path must resolve to the intact initial structure; --data GROUP=PATH can
provide an override. Numeric atom types default to 1:C,2:H,3:N,4:O.

METHOD (unchanged from the original analysis)
  C...O <= 3.5 A; full C-H...O angle >= 140 deg; original C-H <= 1.30 A.
  Initial molecules are reconstructed by covalent-distance connectivity and
  checked against the CL-20/HMX formulas. Molecular origins remain fixed.
  Count unique C-H...O triplets between different initial molecules;
  include CL20->CL20, CL20->HMX and HMX->CL20, excluding HMX->HMX.
  Normalize by the initial CL-20 count (16 for the manuscript systems).
  Average within 1 ps bins per run, then calculate the mean and sample SD
  (ddof=1) across three runs. Bins are [left,right), with the final right
  endpoint included. Time zero is dump step 0 unless --start-step is supplied.
  Periodic boundaries, numerical integrity and duplicate-run checks are kept.

OUTPUT
  Per-frame and per-bin CSVs, comparison_mean_sd.csv, analysis_details.json,
  and comparison_1500K.png/.pdf/.svg. Inputs are read only.
  --plot-only reads saved means/SDs without recounting, renormalizing,
  smoothing or modifying the source CSV, and writes into a new folder.
  All three curves are SOLID. Legend order: cocrystal (black), beta (green),
  gamma (light blue). Shading is mean +/- sample SD; the y-axis starts at 0.
  Font: DejaVu Sans, 10 pt; size: 6.5 x 4.5 in; lines: 1.7 pt; PNG: 600 dpi.

LEGACY INPUT COMPATIBILITY
  Use --restore-sed-numbers ONLY for trajectories known to have undergone
  the old whole-line replacements 1->C, 2->H, 3->N, 4->O in numeric fields.
  Restoration is opt-in, in memory, and logged; input files are never changed.
  Clean numeric trajectories do not need this option.
"""

import argparse
import csv
import hashlib
import io
import itertools
import json
import math
import re
import shlex
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path('.')
RUNS = {
    'pcl20': ['1-100ps-1500k', '2-100ps-1500k', '3-100ps-1500k'],
    'ycl20': ['1-100ps-1500k', '2-100ps-1500k', '3-100ps-1500k'],
    'md': ['1-100-ps-1500k', '5-100ps-1500k', '6-100ps-1500k'],
}
LABELS = {'pcl20': '\u03b2-CL-20', 'ycl20': '\u03b3-CL-20',
          'md': 'CL-20/HMX co-crystal'}
LEGEND_ORDER = ('md', 'pcl20', 'ycl20')
PLOT_FONT = 'DejaVu Sans'
PLOT_FONT_SIZE = 10  # pt: axis labels, tick labels and combined legend
CUTOFFS = {('C', 'H'): 1.25, ('C', 'C'): 1.85, ('C', 'N'): 1.80,
           ('N', 'N'): 1.75, ('N', 'O'): 1.65}
MASSES = {'C': 12.011, 'H': 1.008, 'N': 14.007, 'O': 15.999}
FORMULAS = {'CL20': Counter(C=6, H=6, N=12, O=12),
            'HMX': Counter(C=4, H=8, N=8, O=8)}
FIELDS = ['CL20_to_CL20', 'CL20_to_HMX', 'HMX_to_CL20',
          'total_excluding_HMX_HMX', 'contacts_per_initial_CL20',
          'attached_original_CH_CL20', 'attached_original_CH_HMX']
SED_NUMBER_INVERSE = {'C': '1', 'H': '2', 'N': '3', 'O': '4'}


class NumberRestoration:
    """Invert the demonstrated sed substitutions in explicitly numeric fields."""
    def __init__(self, enabled=False, csv_path=None):
        self.enabled = enabled
        self.csv_path = csv_path
        self.handle = None
        self.writer = None
        self.state = {'enabled': enabled, 'mapping': dict(SED_NUMBER_INVERSE), 'total': 0,
                      'counts_by_field': {}, 'first_examples': [], 'audit_csv': None,
                      'scope': 'Numeric fields used in this analysis; in-memory recovery only'}

    def __enter__(self):
        return self

    def __exit__(self, *args):
        if self.handle is not None:
            self.handle.close()

    def restore(self, token, step, line_number, field, raw_line):
        if not self.enabled or token not in SED_NUMBER_INVERSE:
            return token
        restored = SED_NUMBER_INVERSE[token]
        self.state['total'] += 1
        counts = self.state['counts_by_field']
        counts[field] = counts.get(field, 0) + 1
        example = {'step': step, 'line_in_stream': line_number, 'field': field,
                   'token_before': token, 'token_restored': restored, 'raw_line': raw_line.strip()}
        if len(self.state['first_examples']) < 12:
            self.state['first_examples'].append(example)
        if self.csv_path is not None:
            if self.handle is None:
                self.handle = self.csv_path.open('x', newline='', encoding='utf-8-sig')
                self.writer = csv.DictWriter(self.handle, fieldnames=list(example))
                self.writer.writeheader()
                self.state['audit_csv'] = self.csv_path.name
            self.writer.writerow(example)
        return restored


class DumpReadError(ValueError):
    def __init__(self, path, step, line_number, field, token, raw_line, header=''):
        self.detail = {'file': str(path), 'step': step, 'line_in_stream': line_number,
                       'field': field, 'invalid_token': token, 'raw_line': raw_line.strip(),
                       'atoms_header': header}
        message = ('Invalid numeric field in %s\n  step=%s; line_in_stream=%s; field=%s; token=%r'
                   '\n  RAW: %s' % (path, step, line_number, field, token, raw_line.strip()))
        if header:
            message += '\n  HEADER: ' + header
        message += '\nNo frame was skipped; the source dump was not modified.'
        super().__init__(message)


class Box:
    def __init__(self, origin, matrix):
        self.origin = np.asarray(origin, dtype=float)
        self.matrix = np.asarray(matrix, dtype=float)  # lattice vectors as columns
        if not np.isfinite(self.matrix).all() or np.linalg.det(self.matrix) <= 1e-9:
            raise ValueError('Invalid or left-handed simulation cell')
        self.inverse = np.linalg.inv(self.matrix)
        self.recip_norm = np.linalg.norm(self.inverse, axis=1)

    def nearby_image(self, delta, radius):
        """Exact nearest image if an image is within radius; no fixed 27-image assumption.

        For any vector v inside the radius, |(inverse @ v)[i]| <= radius *
        norm(inverse[i]). This bounds every lattice shift worth enumerating.
        Pairs with no image within radius need only be rejected by the caller.
        """
        frac = np.asarray(delta) @ self.inverse.T
        frac -= np.rint(frac)
        base = frac @ self.matrix.T
        bounds = np.floor(radius * self.recip_norm + 0.5 + 1e-12).astype(int)
        if np.prod(2 * bounds + 1) > 729:
            raise ValueError('Cell is too skewed/small for this image search; inspect the cell')
        best = base.copy()
        best2 = np.sum(best * best, axis=-1)
        for shift in itertools.product(*(range(-int(b), int(b) + 1) for b in bounds)):
            if shift == (0, 0, 0):
                continue
            candidate = base - self.matrix @ np.asarray(shift)
            d2 = np.sum(candidate * candidate, axis=-1)
            take = d2 < best2
            best = np.where(take[..., None], candidate, best)
            best2 = np.minimum(best2, d2)
        return best


def parse_assignments(items, allowed):
    result = {}
    for item in items or []:
        key, sep, value = item.partition('=')
        if not sep or key not in allowed or not value:
            raise ValueError('Expected GROUP=value; GROUP must be ' + ', '.join(allowed))
        if key in result:
            raise ValueError('Duplicate override for ' + key)
        result[key] = value
    return result


def read_input(path, data_override, dt_override):
    """Read literal/simple-variable inputs without executing any LAMMPS or shell code."""
    variables, commands, pending = {}, [], ''
    for line in path.read_text(encoding='utf-8').splitlines():
        line = line.split('#', 1)[0].strip()
        if not line:
            continue
        pending += line.rstrip('&').strip() + ' '
        if line.endswith('&'):
            continue
        line, pending = pending.strip(), ''
        line = re.sub(r'\$\{([^}]+)\}|\$([A-Za-z])',
                      lambda m: variables.get(m.group(1) or m.group(2), m.group(0)), line)
        words = shlex.split(line)
        if words and words[0] == 'variable' and len(words) == 4:
            if words[2] in ('string', 'index', 'equal'):
                variables[words[1]] = words[3]
        if words:
            commands.append(words)
    if pending:
        raise ValueError('Unfinished continuation in ' + str(path))
    for words in commands:
        if words[0] in ('include', 'jump', 'if', 'replicate', 'create_atoms', 'delete_atoms'):
            raise ValueError('Automatic input reader does not support ' + words[0] + ': ' + str(path))
        if words[:2] == ['shell', 'cd']:
            raise ValueError('shell cd changes the input working directory; inspect ' + str(path))
    units = {w[1] for w in commands if w[0] == 'units'}
    if len(units) != 1 or next(iter(units)) not in ('metal', 'real'):
        raise ValueError('Need an explicit, constant units metal or units real in ' + str(path))
    units = next(iter(units))
    if dt_override is None:
        try:
            dts = {float(w[1]) for w in commands if w[0] == 'timestep'}
        except ValueError:
            raise ValueError('Cannot resolve timestep; supply --dt-fs after checking ' + str(path))
        if len(dts) != 1:
            raise ValueError('Need one explicit timestep; otherwise supply --dt-fs after checking ' + str(path))
        dt_fs = next(iter(dts)) * (1000.0 if units == 'metal' else 1.0)
    else:
        dt_fs = dt_override
    if not math.isfinite(dt_fs) or dt_fs <= 0:
        raise ValueError('Invalid timestep in ' + str(path))
    styles = {w[1] for w in commands if w[0] == 'atom_style'}
    style = next(iter(styles)) if len(styles) == 1 else None
    reads = [w for w in commands if w[0] == 'read_data']
    if data_override is not None:
        data = data_override
    else:
        if len(reads) != 1 or len(reads[0]) != 2 or '$' in reads[0][1]:
            raise ValueError('Cannot resolve read_data; supply --data GROUP=/absolute/initial.lmp for ' + str(path))
        data = (path.parent / reads[0][1]).resolve()
    if not data.is_file():
        raise ValueError('Initial data file not found: ' + str(data) + '; use --data GROUP=/absolute/path')
    targets, seeds = [], []
    run_seen = False
    for w in commands:
        if w[0] == 'run':
            run_seen = True
        if w[0] in ('timestep', 'reset_timestep') and run_seen:
            raise ValueError('A timestep/reset_timestep after run makes automatic time mapping ambiguous: ' + str(path))
        if w[0] == 'fix' and len(w) > 6 and w[3] in ('nvt', 'npt', 'nph') and 'temp' in w:
            pos = w.index('temp')
            targets.append([float(w[pos + 1]), float(w[pos + 2])])
        if w[0] == 'velocity' and len(w) > 4 and w[2] == 'create':
            seeds.append(w[4])
    if targets and any(not np.allclose(t, [1500.0, 1500.0]) for t in targets):
        raise ValueError('Input contains a thermostat target other than constant 1500 K: ' + str(path))
    return data, dt_fs, style, {'units': units, 'dt_fs': dt_fs, 'thermostat_targets_K': targets,
                               'velocity_create_seeds': seeds}


def read_data(path, style_hint, type_map):
    lines = path.read_text(encoding='utf-8').splitlines()
    header, sections, current, atom_style = [], {}, None, None
    for raw in lines[1:]:  # LAMMPS always consumes the first line as the title
        line, _, comment = raw.partition('#')
        line = line.strip()
        if not line:
            continue
        if re.match(r'^[A-Za-z]', line):
            current = line
            sections.setdefault(current, [])
            if current == 'Atoms' and comment.strip():
                atom_style = comment.strip().split()[0]
        elif current is None:
            header.append(line.split())
        else:
            sections[current].append(line.split())
    count = next((int(w[0]) for w in header if len(w) == 2 and w[1] == 'atoms'), None)
    if count is None or count <= 0:
        raise ValueError('Missing/invalid atom count: ' + str(path))
    bounds = {}
    tilt = [0.0, 0.0, 0.0]
    for w in header:
        if w[-2:] in (['xlo', 'xhi'], ['ylo', 'yhi'], ['zlo', 'zhi']):
            bounds[w[-2][0]] = list(map(float, w[:2]))
        if w[-3:] == ['xy', 'xz', 'yz']:
            tilt = list(map(float, w[:3]))
    if set(bounds) != set('xyz'):
        raise ValueError('Need orthogonal or restricted triclinic data box: ' + str(path))
    origin = [bounds[k][0] for k in 'xyz']
    lengths = [bounds[k][1] - bounds[k][0] for k in 'xyz']
    xy, xz, yz = tilt
    box = Box(origin, [[lengths[0], xy, xz], [0, lengths[1], yz], [0, 0, lengths[2]]])
    rows = sections.get('Atoms', [])
    if len(rows) != count:
        raise ValueError('Atoms section count mismatch: ' + str(path))
    style = atom_style or style_hint
    layouts = {'atomic': (1, 2, 5), 'charge': (1, 3, 6), 'full': (2, 4, 7),
               'molecular': (2, 3, 6), 'bond': (2, 3, 6), 'angle': (2, 3, 6)}
    if style not in layouts:
        raise ValueError('Unsupported/unknown atom_style ' + str(style) + ' in ' + str(path))
    ti, xi, ncols = layouts[style]
    records = []
    for w in rows:
        if len(w) not in (ncols, ncols + 3):
            raise ValueError('Unexpected Atoms columns for ' + style + ': ' + str(path))
        records.append((int(w[0]), int(w[ti]), *map(float, w[xi:xi + 3])))
    records.sort()
    ids = np.array([w[0] for w in records], dtype=np.int64)
    types = np.array([w[1] for w in records], dtype=np.int64)
    xyz = np.array([w[2:] for w in records], dtype=float)
    if len(set(ids)) != count or ids[0] < 1 or not np.isfinite(xyz).all():
        raise ValueError('Duplicate/invalid atom IDs or coordinates: ' + str(path))
    elements = np.array([type_map[int(t)] for t in types])
    for w in sections.get('Masses', []):
        atom_type, mass = int(w[0]), float(w[1])
        if atom_type in set(types) and abs(mass - MASSES[type_map[atom_type]]) > 0.35:
            raise ValueError('Masses disagree with --type-map at atom type ' + str(atom_type) + ': ' + str(path))
    return ids, types, elements, xyz, box


def make_system(path, style, type_map):
    ids, types, elements, xyz, box = read_data(path, style, type_map)
    n = len(ids)
    delta = box.nearby_image(xyz[None, :, :] - xyz[:, None, :], max(CUTOFFS.values()))
    d2 = np.sum(delta * delta, axis=2)
    cut = np.zeros((n, n))
    for (a, b), radius in CUTOFFS.items():
        mask = ((elements[:, None] == a) & (elements[None, :] == b))
        mask |= ((elements[:, None] == b) & (elements[None, :] == a))
        cut[mask] = radius
    edges = (cut > 0) & (d2 <= cut * cut) & (d2 > 1e-12)
    np.fill_diagonal(edges, False)
    mol = np.full(n, -1, dtype=int)
    kinds = []
    for start in range(n):
        if mol[start] >= 0:
            continue
        index, component, stack = len(kinds), [], [start]
        mol[start] = index
        while stack:
            atom = stack.pop()
            component.append(atom)
            for neighbor in np.flatnonzero(edges[atom]):
                if mol[neighbor] < 0:
                    mol[neighbor] = index
                    stack.append(int(neighbor))
        formula = Counter(elements[component])
        matched = [name for name, expected in FORMULAS.items() if formula == expected]
        if not matched:
            raise ValueError('Initial connectivity gives ' + str(dict(formula)) + ' in ' + str(path)
                             + '. Need the intact initial structure and correct type map; no molecule IDs were guessed.')
        kinds.append(matched[0])
    kind = np.array(kinds)[mol]
    H, O = np.flatnonzero(elements == 'H'), np.flatnonzero(elements == 'O')
    carbons = []
    for h in H:
        neighbors = np.flatnonzero(edges[h] & (elements == 'C'))
        if len(neighbors) != 1:
            raise ValueError('Expected exactly one original C bonded to H id ' + str(ids[h]))
        carbons.append(int(neighbors[0]))
    source_h, source_o = kind[H], kind[O]
    inter = mol[H, None] != mol[O][None, :]
    masks = [inter & (source_h[:, None] == a) & (source_o[None, :] == b)
             for a, b in [('CL20', 'CL20'), ('CL20', 'HMX'), ('HMX', 'CL20')]]
    return {'ids': ids, 'types': types, 'elements': elements, 'xyz': xyz, 'box': box,
            'H': H, 'O': O, 'C': np.array(carbons), 'source_h': source_h,
            'masks': masks, 'counts': dict(Counter(kinds)), 'mol': mol, 'kind': kind}


def read_box(header, rows):
    words = header.split()
    flags = [w for w in words[3:] if re.fullmatch('[pfsm]{2}', w)]
    if flags != ['pp', 'pp', 'pp']:
        raise ValueError('This crystal-analysis script requires periodic boundaries in all three directions')
    if 'abc' in words:
        if 'origin' not in words or any(len(r) != 4 for r in rows):
            raise ValueError('Malformed general triclinic dump box')
        a = np.array(rows)
        return Box(a[:, 3], a[:, :3].T)
    if 'xy' in words:
        if any(len(r) != 3 for r in rows):
            raise ValueError('Malformed restricted triclinic dump box')
        xlb, xhb, xy = rows[0]
        ylb, yhb, xz = rows[1]
        zl, zh, yz = rows[2]
        xl, xh = xlb - min(0, xy, xz, xy + xz), xhb - max(0, xy, xz, xy + xz)
        yl, yh = ylb - min(0, yz), yhb - max(0, yz)
    else:
        if any(len(r) != 2 for r in rows):
            raise ValueError('Malformed orthogonal dump box')
        (xl, xh), (yl, yh), (zl, zh) = rows
        xy = xz = yz = 0.0
    return Box([xl, yl, zl], [[xh - xl, xy, xz], [0, yh - yl, yz], [0, 0, zh - zl]])


def frames(handle, system, restoration=None):
    """Stream strict LAMMPS custom/atom text frames; atom rows may be unsorted."""
    line_number = 0
    source = getattr(handle, 'name', '<stream>')
    restoration = restoration if restoration is not None else NumberRestoration()

    def read_line():
        nonlocal line_number
        line = handle.readline()
        if line:
            line_number += 1
        return line

    def read_integer(field, step):
        line = read_line()
        token = restoration.restore(line.strip(), step, line_number, field, line)
        try:
            return int(token)
        except ValueError:
            raise DumpReadError(source, step, line_number, field, line.strip(), line) from None

    while True:
        line = read_line()
        if not line:
            return
        if not line.strip():
            continue
        if line.strip() != 'ITEM: TIMESTEP':
            raise ValueError('Expected ITEM: TIMESTEP, got ' + line[:100])
        step = read_integer('TIMESTEP', None)
        if read_line().strip() != 'ITEM: NUMBER OF ATOMS':
            raise ValueError('Missing NUMBER OF ATOMS')
        n = read_integer('NUMBER OF ATOMS', step)
        if n != len(system['ids']):
            raise ValueError('Dump/data atom count mismatch at step ' + str(step))
        header = read_line().strip()
        if not header.startswith('ITEM: BOX BOUNDS'):
            raise ValueError('Missing BOX BOUNDS at step ' + str(step))
        box_rows = []
        for axis in range(3):
            raw = read_line()
            values = []
            for j, token in enumerate(raw.split()):
                numeric = restoration.restore(token, step, line_number,
                                              'BOX BOUNDS row%d col%d' % (axis + 1, j + 1), raw)
                try:
                    value = float(numeric)
                    if not math.isfinite(value):
                        raise ValueError('non-finite')
                except ValueError:
                    raise DumpReadError(source, step, line_number,
                                        'BOX BOUNDS row%d col%d' % (axis + 1, j + 1), token, raw) from None
                values.append(value)
            box_rows.append(values)
        box = read_box(header, box_rows)
        header = read_line().strip()
        if not header.startswith('ITEM: ATOMS '):
            raise ValueError('Missing ATOMS at step ' + str(step))
        cols = header.split()[2:]
        if len(cols) != len(set(cols)) or 'id' not in cols:
            raise ValueError('Need unique ATOMS column names including id')
        coord_sets = [('x', 'y', 'z'), ('xu', 'yu', 'zu'), ('xs', 'ys', 'zs'), ('xsu', 'ysu', 'zsu')]
        chosen = next((c for c in coord_sets if set(c).issubset(cols)), None)
        if chosen is None:
            raise ValueError('Need x/y/z, xu/yu/zu, xs/ys/zs or xsu/ysu/zsu')
        first_atom_line = line_number + 1
        raw_rows = [read_line() for _ in range(n)]
        rows = [line.split() for line in raw_rows]
        if any(len(w) != len(cols) for w in rows):
            raise ValueError('Truncated or malformed atom records at step ' + str(step))
        arr = np.asarray(rows)
        try:
            ids = arr[:, cols.index('id')].astype(np.int64)
        except (ValueError, OverflowError):
            ids = np.empty(n, dtype=np.int64)
            for row_index, w in enumerate(rows):
                token = restoration.restore(w[cols.index('id')], step, first_atom_line + row_index,
                                            'id', raw_rows[row_index])
                try:
                    ids[row_index] = int(token)
                except (ValueError, OverflowError):
                    raise DumpReadError(source, step, first_atom_line + row_index, 'id',
                                        w[cols.index('id')], raw_rows[row_index], header) from None
        order = np.argsort(ids)
        if not np.array_equal(ids[order], system['ids']):
            raise ValueError('Missing/duplicate/unknown atom ID at step ' + str(step))
        arr = arr[order]
        for field in ('type', 'element'):
            if field in cols:
                actual = arr[:, cols.index(field)]
                if field == 'type':
                    # sed can leave numeric and symbolic types in the same frame.
                    ok = np.all((actual == system['types'].astype(str)) | (actual == system['elements']))
                else:
                    ok = np.array_equal(actual, system['elements'])
                if not ok:
                    raise ValueError(field + ' mapping disagrees with initial data at step ' + str(step))
        try:
            xyz = arr[:, [cols.index(c) for c in chosen]].astype(float)
        except ValueError:
            recovered = np.empty((n, 3), dtype=float)
            for row_index, w in enumerate(rows):
                for j, field in enumerate(chosen):
                    token = w[cols.index(field)]
                    numeric = restoration.restore(token, step, first_atom_line + row_index,
                                                  field, raw_rows[row_index])
                    try:
                        recovered[row_index, j] = float(numeric)
                    except ValueError:
                        raise DumpReadError(source, step, first_atom_line + row_index, field,
                                            token, raw_rows[row_index], header) from None
            xyz = recovered[order]
        if chosen[0] in ('xs', 'xsu'):
            xyz = box.origin + xyz @ box.matrix.T
        if not np.isfinite(xyz).all():
            row_index, column_index = np.argwhere(~np.isfinite(xyz))[0]
            original_row = int(order[row_index])
            field = chosen[column_index]
            raise DumpReadError(source, step, first_atom_line + original_row, field,
                                rows[original_row][cols.index(field)], raw_rows[original_row], header)
        yield step, xyz, box, cols


def peek_ends(path, system, restore_sed_numbers=False):
    with path.open('r', encoding='utf-8') as handle:
        reader = frames(handle, system, NumberRestoration(restore_sed_numbers))
        first = next(reader)
        second = next(reader)
    with path.open('rb') as handle:
        size = path.stat().st_size
        handle.seek(max(0, size - 4 * 1024 * 1024))
        tail = handle.read().decode('utf-8')
    markers = list(re.finditer(r'(?m)^ITEM: TIMESTEP\r?$', tail))
    if not markers:
        raise ValueError('Cannot locate final frame in ' + str(path))
    last = next(frames(io.StringIO(tail[markers[-1].start():]), system,
                       NumberRestoration(restore_sed_numbers)))
    return first, second, last


def count_contacts(xyz, box, system, args):
    h, c, o = system['H'], system['C'], system['O']
    ch = box.nearby_image(xyz[h] - xyz[c], args.ch_max)
    ch2 = np.sum(ch * ch, axis=1)
    attached = (ch2 > 1e-12) & (ch2 <= args.ch_max ** 2)
    co = box.nearby_image(xyz[o][None, :, :] - xyz[c][:, None, :], args.distance)
    co2 = np.sum(co * co, axis=2)
    # Use one consistent image for the C-H-O triangle; independently wrapping
    # H...O could produce an angle inconsistent with the chosen donor/acceptor.
    ho = co - ch[:, None, :]
    ho2 = np.sum(ho * ho, axis=2)
    denom = np.sqrt(ch2[:, None] * ho2)
    cosine = np.divide(np.einsum('ij,ikj->ik', -ch, ho), denom,
                       out=np.ones_like(denom), where=denom > 1e-12)
    valid = (attached[:, None] & (co2 <= args.distance ** 2) & (denom > 1e-12)
             & (cosine <= math.cos(math.radians(args.angle_min))))
    counts = [int(np.count_nonzero(valid & m)) for m in system['masks']]
    total = sum(counts)
    active = [int(np.count_nonzero(attached & (system['source_h'] == k))) for k in ('CL20', 'HMX')]
    return counts + [total, total / system['counts']['CL20']] + active


def preflight(args, type_map, data_overrides):
    jobs, cache, errors, signatures = [], {}, [], {}
    for group, runs in RUNS.items():
        for replicate, name in enumerate(runs, start=1):
            directory = args.root / group / name
            try:
                input_path, traj = directory / 'input.lammps', directory / 'dump.lammpstrj'
                if not input_path.is_file() or not traj.is_file():
                    raise ValueError('Need both input.lammps and dump.lammpstrj in ' + str(directory))
                override = data_overrides.get(group)
                override = (args.root / override).resolve() if override else None
                data, dt, style, info = read_input(input_path, override, args.dt_fs)
                key = (str(data), style)
                if key not in cache:
                    cache[key] = make_system(data, style, type_map)
                system = cache[key]
                counts = system['counts']
                if counts.get('CL20', 0) == 0 or ((group == 'md') != (counts.get('HMX', 0) > 0)):
                    raise ValueError('Initial composition does not match group ' + group + ': ' + str(counts))
                first, second, last = peek_ends(traj, system, args.restore_sed_numbers)
                stride = second[0] - first[0]
                t0 = (first[0] - args.start_step) * dt / 1000.0
                t1 = (last[0] - args.start_step) * dt / 1000.0
                if stride <= 0 or last[0] <= second[0]:
                    raise ValueError('Need increasing timesteps and at least three frames')
                if t0 > stride * dt / 1000.0 + 1e-8 or t1 < args.duration_ps - 1e-8:
                    raise ValueError('Time coverage %.6g to %.6g ps does not cover 0 to %.6g ps. '
                                     'Inspect timestep/time origin; --start-step sets the zero step.'
                                     % (t0, t1, args.duration_ps))
                if stride * dt / 1000.0 > args.bin_ps + 1e-8:
                    raise ValueError('Dump spacing exceeds --bin-ps; use a wider bin')
                endpoint_hash = hashlib.sha256()
                for frame in (first, last):
                    endpoint_hash.update(frame[1].tobytes())
                    endpoint_hash.update(frame[2].matrix.tobytes())
                signature = endpoint_hash.hexdigest()
                if signature in signatures:
                    print('  WARNING: identical first/last coordinates in %s and %s; '
                          'independence needs checking.' % (signatures[signature], directory), flush=True)
                signatures[signature] = str(directory)
                info.update(group=group, replicate=replicate, run=name, trajectory=str(traj),
                            initial_data=str(data), initial_counts=counts, atom_count=len(system['ids']),
                            first_step=first[0], last_step=last[0], first_time_ps=t0, last_time_ps=t1,
                            dump_interval_ps=stride * dt / 1000.0, dump_columns=first[3],
                            initial_data_sha256=hashlib.sha256(data.read_bytes()).hexdigest())
                jobs.append({'group': group, 'replicate': replicate, 'name': name, 'path': traj,
                             'system': system, 'dt': dt, 'stride': stride, 'info': info})
                print('[OK] %s/%s | %s | dt=%.6g fs | %.6g--%.6g ps | seed=%s'
                      % (group, name, counts, dt, t0, t1, info['velocity_create_seeds'] or 'not found'), flush=True)
                print('     initial data: ' + str(data), flush=True)
                if not info['thermostat_targets_K']:
                    print('     NOTE: 1500 K target was not verified from an nvt/npt temp command.', flush=True)
            except (OSError, ValueError, KeyError, IndexError, StopIteration) as exc:
                errors.append(str(directory) + ': ' + (str(exc) or 'empty/short input'))
                print('[ERROR] ' + errors[-1], flush=True)
    if errors:
        raise ValueError('%d preflight check(s) failed. No trajectories were analyzed. '
                         'Inspect the [ERROR] lines and the corresponding input/data files.' % len(errors))
    for group in RUNS:
        group_jobs = [job for job in jobs if job['group'] == group]
        if len({json.dumps(j['system']['counts'], sort_keys=True) for j in group_jobs}) != 1:
            raise ValueError('Replicates in ' + group + ' have different initial molecule counts')
        seeds = [j['info']['velocity_create_seeds'] for j in group_jobs]
        if seeds[0] and all(s == seeds[0] for s in seeds):
            print('NOTE: all %s runs show the same velocity-create seed. '
                  'Check whether distinct initial states/other random seeds establish independent runs.' % group,
                  flush=True)
    print('Preflight passed for all nine trajectories. Full-frame integrity is checked during analysis.', flush=True)
    return jobs


def update_coordinate_hash(digest, step, xyz, box):
    # Hash the coordinate sequence to detect duplicate trajectories.
    digest.update(str(step).encode('ascii'))
    digest.update(xyz.tobytes())
    digest.update(box.matrix.tobytes())


def analyze(job, args, out):
    group, rep, path = job['group'], job['replicate'], job['path']
    nbin = int(round(args.duration_ps / args.bin_ps))
    sums = np.zeros((nbin, len(FIELDS)))
    numbers = np.zeros(nbin, dtype=np.int64)
    raw_path = out / ('%s_run%d_frames.csv' % (group, rep))
    n_frames, n_used, previous = 0, 0, None
    digest = hashlib.sha256()
    before = path.stat()
    restoration = NumberRestoration(args.restore_sed_numbers, out / ('%s_run%d_restored_numbers.csv' % (group, rep)))
    job['info']['numeric_restoration'] = restoration.state
    with restoration, raw_path.open('w', newline='', encoding='utf-8-sig') as output, path.open('r', encoding='utf-8') as handle:
        writer = csv.writer(output)
        writer.writerow(['step', 'time_ps'] + FIELDS)
        for step, xyz, box, columns in frames(handle, job['system'], restoration):
            if previous is not None and step - previous != job['stride']:
                raise ValueError('Missing/irregular/duplicate timesteps in %s at %s; no interpolation was performed'
                                 % (path, step))
            previous = step
            t = (step - args.start_step) * job['dt'] / 1000.0
            n_frames += 1
            update_coordinate_hash(digest, step, xyz, box)
            if t < -1e-8 or t > args.duration_ps + 1e-8:
                continue
            values = count_contacts(xyz, box, job['system'], args)
            writer.writerow([step, '%.10g' % t] + values)
            bin_index = min(nbin - 1, max(0, int(math.floor((t + 1e-10) / args.bin_ps))))
            sums[bin_index] += values
            numbers[bin_index] += 1
            n_used += 1
            if n_used == 1 or n_used % 1000 == 0:
                print('  %s run%d: %d frames counted, t=%.3f ps' % (group, rep, n_used, t), flush=True)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError('Trajectory changed while reading: ' + str(path))
    if np.any(numbers == 0):
        raise ValueError('Empty time bins in ' + str(path) + '; no endpoint extrapolation was performed')
    means = sums / numbers[:, None]
    centers = (np.arange(nbin) + 0.5) * args.bin_ps
    with (out / ('%s_run%d_bins.csv' % (group, rep))).open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.writer(handle)
        writer.writerow(['time_ps', 'bin_left_ps', 'bin_right_ps', 'frame_count'] + FIELDS)
        for i, t in enumerate(centers):
            writer.writerow([t, i * args.bin_ps, (i + 1) * args.bin_ps, numbers[i]] + means[i].tolist())
    job['info'].update(frames_read=n_frames, frames_used=n_used, trajectory_coordinate_sha256=digest.hexdigest(),
                       frames_csv_sha256=hashlib.sha256(raw_path.read_bytes()).hexdigest(),
                       bins_csv_sha256=hashlib.sha256((out / ('%s_run%d_bins.csv' % (group, rep))).read_bytes()).hexdigest())
    if restoration.state['total']:
        print('[RESTORED] %s run%d: %d numeric tokens; %s'
              % (group, rep, restoration.state['total'], restoration.state['counts_by_field']), flush=True)
    return means


def write_summary(out, times, results):
    with (out / 'comparison_mean_sd.csv').open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.writer(handle)
        columns = ['time_ps']
        for group in RUNS:
            columns += [group + '_' + field + '_' + stat for field in FIELDS for stat in ('mean', 'sd')]
        writer.writerow(columns)
        stats = {g: (arr.mean(axis=0), arr.std(axis=0, ddof=1)) for g, arr in results.items()}
        for i, t in enumerate(times):
            row = [t]
            for group in RUNS:
                mean, sd = stats[group]
                row.extend(x for j in range(len(FIELDS)) for x in (mean[i, j], sd[i, j]))
            writer.writerow(row)


def plot_mean_sd(out, times, stats, labels, args):
    """Shared plotting code for a fresh analysis and saved-CSV redraws."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    colors = {'pcl20': '#008a5c', 'ycl20': '#81cfff', 'md': '#000000'}
    styles = {group: '-' for group in RUNS}  # Figure 6: all curves solid
    plt.rcParams.update({'font.family': PLOT_FONT, 'font.size': PLOT_FONT_SIZE,
                         'axes.labelsize': PLOT_FONT_SIZE, 'legend.fontsize': PLOT_FONT_SIZE,
                         'xtick.labelsize': PLOT_FONT_SIZE, 'ytick.labelsize': PLOT_FONT_SIZE,
                         'axes.linewidth': 0.8, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                         'xtick.direction': 'in', 'ytick.direction': 'in'})
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    lines = {}
    for group in RUNS:
        mean, sd = stats[group]
        ax.fill_between(times, mean - sd, mean + sd, color=colors[group], alpha=0.16, linewidth=0)
        lines[group], = ax.plot(times, mean, color=colors[group], linestyle=styles[group],
                               linewidth=1.7, label=labels[group])
    ax.set(xlim=(0, args.duration_ps), xlabel='Time (ps)',
           ylabel='C–H···O contacts per initial CL-20 molecule')
    ax.legend([lines[group] for group in LEGEND_ORDER],
              [labels[group] for group in LEGEND_ORDER], frameon=False, loc='upper right')
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.tick_params(top=True, right=True)
    fig.tight_layout()
    for extension in ('png', 'pdf', 'svg'):
        fig.savefig(out / ('comparison_1500K.' + extension), dpi=600, bbox_inches='tight')
    plt.close(fig)
    print('Plot font: %s; axis labels, tick labels and legend: %g pt; PNG: 600 dpi.'
          % (PLOT_FONT, PLOT_FONT_SIZE), flush=True)


def plot_comparison(out, times, results, labels, args):
    stats = {}
    for group in RUNS:
        values = results[group][:, :, 4]
        stats[group] = (values.mean(axis=0), values.std(axis=0, ddof=1))
    plot_mean_sd(out, times, stats, labels, args)


def replot_summary(args, labels):
    """Redraw the combined figure using only the original saved mean/SD CSV."""
    source = args.plot_only.expanduser().resolve()
    if source.is_dir():
        source = source / 'comparison_mean_sd.csv'
    required = ['time_ps']
    for group in RUNS:
        required += [group + '_contacts_per_initial_CL20_' + stat
                     for stat in ('mean', 'sd')]
    with source.open('r', newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        missing = [name for name in required if name not in (reader.fieldnames or [])]
        if missing:
            raise ValueError('Missing columns in %s: %s' % (source, ', '.join(missing)))
        rows = []
        for line, row in enumerate(reader, start=2):
            try:
                rows.append([float(row[name]) for name in required])
            except (TypeError, ValueError) as error:
                raise ValueError('Invalid numeric value in %s at CSV line %d'
                                 % (source, line)) from error
    if not rows:
        raise ValueError('No data rows in ' + str(source))
    data = np.asarray(rows, dtype=float)
    if not np.isfinite(data).all() or np.any(data < 0):
        raise ValueError('Saved times, means and SDs must be finite and nonnegative: ' + str(source))
    times = data[:, 0]
    width = float(times[1] - times[0]) if len(times) > 1 else float(2 * times[0])
    if width <= 0 or not np.allclose(times, (np.arange(len(times)) + 0.5) * width,
                                     rtol=0, atol=1e-8):
        raise ValueError('Expected increasing, equally spaced bin centers starting at half a bin: '
                         + str(source))
    args.bin_ps, args.duration_ps = width, float(len(times) * width)
    stats = {group: (data[:, 1 + 2 * i], data[:, 2 + 2 * i])
             for i, group in enumerate(RUNS)}
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    out = (args.outdir or source.parent / ('replot_' + stamp)).expanduser().resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError('Output directory must be new or empty: ' + str(out))
    out.mkdir(parents=True, exist_ok=True)
    print('Redrawing saved means and SDs from: ' + str(source), flush=True)
    plot_mean_sd(out, times, stats, labels, args)
    print('Redraw finished. PNG/PDF/SVG saved in: ' + str(out), flush=True)
    print('Combined figure: ' + str(out / 'comparison_1500K.png'), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='compare_cl20_hbonds 2.0 (Figure 6, solid curves)')
    parser.add_argument('--root', type=Path, default=ROOT, help='Parent of pcl20, ycl20 and md')
    parser.add_argument('--check-only', action='store_true', help='Check inputs and first/last frames, then exit')
    parser.add_argument('--plot-only', type=Path, nargs='?', const=Path('.'), metavar='CSV_OR_FOLDER',
                        help='Redraw saved comparison_mean_sd.csv only; defaults to the current folder')
    parser.add_argument('--restore-sed-numbers', action='store_true',
                        help='Invert confirmed 1->C,2->H,3->N,4->O sed substitutions in numeric fields while reading')
    parser.add_argument('--data', action='append', metavar='GROUP=PATH',
                        help='Override initial data for pcl20, ycl20 or md; repeat as needed')
    parser.add_argument('--label', action='append', metavar='GROUP=TEXT', help='Override a legend label')
    parser.add_argument('--dt-fs', type=float, help='Explicit timestep in fs; default reads units/timestep from each input')
    parser.add_argument('--start-step', type=int, default=0, help='Dump step corresponding to time zero (default: 0)')
    parser.add_argument('--duration-ps', type=float, default=100.0)
    parser.add_argument('--bin-ps', type=float, default=1.0, help='Within-run averaging interval (default: 1 ps)')
    parser.add_argument('--distance', type=float, default=3.5, help='Maximum C...O distance, A (default: 3.5)')
    parser.add_argument('--angle-min', type=float, default=140.0, help='Minimum FULL C-H...O angle (default: 140 degrees)')
    parser.add_argument('--ch-max', type=float, default=1.30, help='Maximum current original C-H distance, A')
    parser.add_argument('--type-map', default='1:C,2:H,3:N,4:O', help='Numeric LAMMPS type to element')
    parser.add_argument('--outdir', type=Path, help='New output directory; default is timestamped under --root')
    args = parser.parse_args()
    labels = dict(LABELS, **parse_assignments(args.label, RUNS))
    if args.plot_only is not None:
        if args.check_only:
            parser.error('--plot-only cannot be combined with --check-only')
        replot_summary(args, labels)
        return
    for value in (args.duration_ps, args.bin_ps, args.distance, args.ch_max):
        if not math.isfinite(value) or value <= 0:
            parser.error('Time/distance parameters must be positive and finite')
    if not 0 < args.angle_min <= 180:
        parser.error('--angle-min must be in (0,180]')
    nbin = args.duration_ps / args.bin_ps
    if nbin < 1 or not math.isclose(nbin, round(nbin), abs_tol=1e-8):
        parser.error('--duration-ps must be an integer multiple of --bin-ps')
    args.root = args.root.expanduser().resolve()
    type_map = {int(k.strip()): v.strip() for k, v in (item.split(':') for item in args.type_map.split(','))}
    if not type_map or any(v not in MASSES for v in type_map.values()):
        parser.error('--type-map may contain only C, H, N, O')
    data_overrides = parse_assignments(args.data, RUNS)
    print('Definition: C...O <= %.3g A; C-H...O >= %.3g deg; original C-H <= %.3g A'
          % (args.distance, args.angle_min, args.ch_max), flush=True)
    print('Count initial-origin intermolecular contacts; HMX->HMX excluded. Time zero: step %d.'
          % args.start_step, flush=True)
    if args.restore_sed_numbers:
        print('Confirmed sed-number recovery enabled: C->1, H->2, N->3, O->4 in numeric fields. '
              'Type/element labels retained; original dumps are read only.', flush=True)
    jobs = preflight(args, type_map, data_overrides)
    if args.check_only:
        print('Check completed. No full trajectory analysis or figure was performed.', flush=True)
        return
    # Fail early if plotting dependencies are absent, before reading gigabytes.
    import matplotlib
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    out = (args.outdir or args.root / ('hbond_compare_1500k_' + stamp)).expanduser().resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError('Output directory is not empty; use a new --outdir to preserve previous results: ' + str(out))
    out.mkdir(parents=True, exist_ok=True)
    print('Output directory: ' + str(out), flush=True)
    config = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}
    method = {'configuration': config, 'labels': labels, 'type_map': type_map,
              'connectivity_cutoffs_A': {'-'.join(k): v for k, v in CUTOFFS.items()},
              'method': __doc__, 'status': 'running', 'numpy_version': np.__version__,
              'matplotlib_version': matplotlib.__version__,
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    report = out / 'analysis_details.json'
    report.write_text(json.dumps(dict(method, trajectories=[j['info'] for j in jobs]), indent=2), encoding='utf-8')
    results = {group: [] for group in RUNS}
    seen = {}
    try:
        for job in jobs:
            values = analyze(job, args, out)
            digest = job['info']['trajectory_coordinate_sha256']
            if digest in seen:
                raise ValueError('Identical full coordinate/step sequences: %s and %s. '
                                 'Cannot count the same trajectory as an independent replicate.'
                                 % (seen[digest], job['path']))
            seen[digest] = str(job['path'])
            results[job['group']].append(values)
        results = {g: np.stack(v) for g, v in results.items()}
        times = (np.arange(int(round(nbin))) + 0.5) * args.bin_ps
        write_summary(out, times, results)
        plot_comparison(out, times, results, labels, args)
        method['status'] = 'completed'
    except Exception as exc:
        method['status'], method['error'] = 'failed', str(exc)
        raise
    finally:
        report.write_text(json.dumps(dict(method, trajectories=[j['info'] for j in jobs]), indent=2), encoding='utf-8')
    print('Finished all nine trajectories. Combined figure: ' + str(out / 'comparison_1500K.png'), flush=True)
    print('PDF/SVG, per-frame/per-bin CSVs and comparison_mean_sd.csv were also saved.', flush=True)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError, IndexError, StopIteration, ImportError) as error:
        print('\nERROR: ' + (str(error) or error.__class__.__name__), file=sys.stderr)
        sys.exit(1)
