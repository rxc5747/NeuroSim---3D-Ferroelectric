#!/usr/bin/env python3
"""
build_rw_xlsx.py -- turn the read/write-operator-mode sweep CSV into a workbook
with every latency and energy auto-scaled to its natural unit.

    python3 build_rw_xlsx.py rw_sq_5_2.csv [-o 3D_FeRAM_rw_square.xlsx] [--grid-scheme vanilla]

Works on the output of sweep.sh (mixed aspect ratios) AND sweep_square.sh
(rows == cols, swept over subarray size S and plane count n) -- the CSV header
is the same.  The Grid_* sheets only make sense for the square sweep: they
pivot every headline metric into an S (down) x n (across) table, one table per
metric, for ONE write scheme (--grid-scheme, default vanilla).  Those tables
are the numbers behind the paper figures; plot_rw_sweep.py draws them straight
from the CSV.

Sheets
    README                 what each sheet is, units, how the CSV maps in
    Grid_Energy            S x n pivots: read fJ/bit, write fJ/bit, true write fJ/bit
    Grid_Latency           S x n pivots: read / write macro latency (ns)
    Grid_Area              S x n pivots: macro area, density, array efficiency, leakage
    Grid_Sense             S x n pivots: dV (mV), sense verdict, V_half/V_coercive
    Summary                one row per run, headline numbers + the 3-way area split
    Physics                cell/stack physics actually used by that run
    Read_Breakdown         read energy split: switch / linear / transistor DC / lines / restore
    Write_Breakdown        write energy split: cells / inhibit / lines, + scheme separators
    Write_Schemes          the three write schemes side by side per configuration
    Latency_Breakdown      select / charge / settle / sense / latch / restore / selectWrite
    Block_Read_Energy      15 peripheral blocks x runs
    Block_Write_Energy
    Block_Read_Latency
    Block_Write_Latency
    Block_Leakage
    Block_Area
    Area                   core / periphery / macro area, density, overhead
    Scaling                read & write fJ/bit vs n, one table per geometry

CHANGES in THIS revision (square sweep)
    * Grid_* sheets: S x n pivot tables for the square sweep (see above).
      Non-square runs (EXTRA=...) are listed under each grid, not dropped.
    * The Scaling caption no longer quotes numbers from one particular run;
      it states the mechanism and lets the table carry the values.
    * Scaling tables are one per geometry as before, but a square sweep also
      gets one table per n (size sweep at fixed n), since that is the other
      axis of the paper figures.

EARLIER (post-audit)
    * FIX 1  sweep.sh's write-group header now matches rw_main's emit order, so
      the last six Write_Breakdown columns stop reporting each other's values.
      nHalfSel is exposed as its own column.
    * FIX 2  columns belonging to one breakdown SHARE a unit (unit_groups), so
      a total and its components can be added up on the page.  This is what hid
      the V/2 inhibit (85-96% of the core write) behind a unit change, and what
      made the restore look as though it had been dropped from the totals.
    * FIX 4  planeDecoder joins BLOCKS, immediately after sslSwitchMatrix.
    * FIX 5  'PASS (clipped)' is recognised and shaded amber.
    * FIX 7  half-select bias, field and V_half/V_coercive are surfaced.
    * FIX 10 the Block_* LATENCY total is labelled as a sum, not a critical path.
    * The Scaling caption said read fJ/bit FALLS with n.  It rises.  Fixed.

EARLIER changes
    * Summary now carries CORE ARRAY / PERIPHERY / TOTAL area plus array
      efficiency.  Periphery is DERIVED here as macro - core, so no C++ change
      and no sweep.sh field renumbering is needed.
    * areaFlag catches the NeuroSim area-model failure at small row counts
      (the 8x8x8 run reports a negative macro area because four row-direction
      decoders each return ~-9.8e8 um^2).  Those rows are flagged, not hidden.
    * Physics picks up senseMargin / senseSigma / nSeg / nSegments
      automatically IF the CSV has them.  Runs made before the statistical
      sense model was added simply do not show those columns.
    * Latency_Breakdown note corrected: the write pulse is
      param->writePulseWidth, not the NLS tau.  tauSwitch is a printed
      cross-check only.
"""
import argparse, csv, math, os, sys
from collections import OrderedDict

import openpyxl
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter

# --------------------------------------------------------------------------- #
#  unit scaling                                                               #
# --------------------------------------------------------------------------- #
LADDERS = {
    's': [(1.0, 's'), (1e-3, 'ms'), (1e-6, 'us'), (1e-9, 'ns'), (1e-12, 'ps')],
    'J': [(1.0, 'J'), (1e-3, 'mJ'), (1e-6, 'uJ'), (1e-9, 'nJ'), (1e-12, 'pJ'),
          (1e-15, 'fJ'), (1e-18, 'aJ')],
    'W': [(1.0, 'W'), (1e-3, 'mW'), (1e-6, 'uW'), (1e-9, 'nW'), (1e-12, 'pW')],
    'F': [(1.0, 'F'), (1e-6, 'uF'), (1e-9, 'nF'), (1e-12, 'pF'), (1e-15, 'fF'),
          (1e-18, 'aF')],
    'C': [(1.0, 'C'), (1e-9, 'nC'), (1e-12, 'pC'), (1e-15, 'fC'), (1e-18, 'aC')],
    'A': [(1.0, 'A'), (1e-3, 'mA'), (1e-6, 'uA'), (1e-9, 'nA'), (1e-12, 'pA')],
    'V': [(1.0, 'V'), (1e-3, 'mV'), (1e-6, 'uV')],
    'm2': [(1.0, 'm^2'), (1e-12, 'um^2'), (1e-18, 'nm^2')],
}

def pick_unit(values, kind, ref_mode='median'):
    """Choose the ladder rung so the reference magnitude lands in [1, 1000).

    ref_mode='median'  standalone column -- robust to one huge/tiny run.
    ref_mode='max'     SHARED unit across a group of columns: the largest member
                       lands in [1,1000) and every smaller member prints below
                       it, so the reader can add the row up by eye.

    FIX 2: the per-column median is what put the write cells term in fJ and the
    V/2 inhibit in pJ in the same row, which made an 85%-of-the-total term look
    like 1%.  It also put the restore energy in aJ next to fJ neighbours and the
    restore latency in ps next to ns neighbours.  Breakdown columns must share
    one unit; see unit_groups in table().
    """
    vals = [abs(v) for v in values if isinstance(v, float) and v != 0 and math.isfinite(v)]
    if not vals:
        return LADDERS[kind][0]
    vals.sort()
    ref = vals[-1] if ref_mode == 'max' else vals[len(vals) // 2]
    best = LADDERS[kind][0]
    for scale, name in LADDERS[kind]:
        if ref / scale >= 1.0:
            best = (scale, name)
            break
    else:
        best = LADDERS[kind][-1]
    return best

# --------------------------------------------------------------------------- #
#  styling                                                                    #
# --------------------------------------------------------------------------- #
HDR_FILL = PatternFill('solid', fgColor='1F3864')
HDR_FONT = Font(bold=True, color='FFFFFF', size=10)
SUB_FILL = PatternFill('solid', fgColor='D9E1F2')
SUB_FONT = Font(bold=True, size=9, color='1F3864')
KEY_FILL = PatternFill('solid', fgColor='FFF2CC')
FAIL_FILL = PatternFill('solid', fgColor='FFC7CE')
PASS_FILL = PatternFill('solid', fgColor='C6EFCE')
WARN_FILL = PatternFill('solid', fgColor='FFE699')
THIN = Side(style='thin', color='BFBFBF')
BOX = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)

BAD_TOKENS = ('FAIL', 'NO SIGNAL', 'NEGATIVE')
# FIX 5: rw_main now reports 'PASS (clipped)' when Qsw/capSense exceeds the
# rail headroom -- the signal is readable but the raw number was unphysical.
WARN_TOKENS = ('PASS (clipped)',)
OK_TOKENS = ('PASS', 'OK')

def write_header(ws, row, labels, units=None, widths=None, freeze=None):
    for j, lab in enumerate(labels, start=1):
        c = ws.cell(row=row, column=j, value=lab)
        c.fill, c.font, c.border = HDR_FILL, HDR_FONT, BOX
        c.alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)
    if units is not None:
        for j, u in enumerate(units, start=1):
            c = ws.cell(row=row + 1, column=j, value=u)
            c.fill, c.font, c.border = SUB_FILL, SUB_FONT, BOX
            c.alignment = Alignment(horizontal='center')
    if widths:
        for j, w in enumerate(widths, start=1):
            ws.column_dimensions[get_column_letter(j)].width = w
    if freeze:
        ws.freeze_panes = freeze

def put(ws, r, c, v, fmt=None, fill=None, bold=False):
    cell = ws.cell(row=r, column=c, value=v)
    cell.border = BOX
    if fmt:
        cell.number_format = fmt
    if fill:
        cell.fill = fill
    if bold:
        cell.font = Font(bold=True, size=10)
    return cell

# --------------------------------------------------------------------------- #
#  CSV                                                                        #
# --------------------------------------------------------------------------- #
# FIX 4: planeDecoder is the LAST entry of rw_main's P[] array, so it must sit
# immediately after sslSwitchMatrix here and in sweep.sh's BLOCKS.  Get the
# position wrong and every per-block column after it is misattributed.
BLOCKS = ['rblDecoder', 'rblSwitchMatrix', 'rslDecoder', 'rslSwitchMatrix',
          'wblDecoder', 'wblSwitchMatrix', 'wwlDecoder', 'wwlSwitchMatrix',
          'wplDecoder', 'wplSwitchMatrix', 'sslDecoder', 'sslSwitchMatrix',
          'planeDecoder',
          'currentSenseAmp', 'outputLatch', 'cellArray']

TEXT_KEYS = ('cell', 'rdo', 'scheme', 'senseVerdict')

def fnum(s):
    try:
        v = float(s)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None

def has(rows, key):
    """True if at least one run actually carries this column."""
    return any(r.get(key) is not None for r in rows)

def load(path):
    with open(path, newline='') as fh:
        rdr = csv.DictReader(fh)
        rows = [r for r in rdr if r.get('cell') and r.get('n')]
    if not rows:
        sys.exit(f"{path}: no data rows -- did the sweep actually run?")
    for r in rows:
        for k, v in list(r.items()):
            if k in TEXT_KEYS:
                r[k] = (v or '').strip()
            else:
                r[k] = fnum(v)

    # Guard: a verdict is only meaningful if the signal it judged is a real
    # number. capNode == 0 (an uncomputed bit-line capacitance) makes
    # dVsense = Qsw/0 = inf, and `inf > vmin` reports PASS on nothing at all.
    for r in rows:
        if r.get('dVsense_V') is None or r.get('capNode_F') in (None, 0.0):
            r['senseVerdict'] = 'NO SIGNAL'

    # ------------------------------------------------------------------ area
    # PERIPHERY IS DERIVED, not emitted.  areaCore_um2 and areaMacro_um2 are
    # already in the CSV, so computing the difference here keeps the CSV column
    # order -- and therefore sweep.sh's awk field numbers -- completely
    # untouched.
    #
    # areaFlag catches the NeuroSim area-model failure that shows up at small
    # row counts: below roughly 16 rows the array is narrower than one
    # peripheral cell row and four row-direction decoders each return about
    # -9.8e8 um^2, so the 8x8x8 run reports a macro area of -6.45e9 um^2.
    # Flag it rather than letting a negative periphery look like a result.
    for r in rows:
        core, macro = r.get('areaCore_um2'), r.get('areaMacro_um2')
        if core is None or macro is None:
            r['areaPeriph_um2'] = None
            r['arrayEff_pct'] = None
            r['areaFlag'] = 'n/a'
            continue
        periph = macro - core
        r['areaPeriph_um2'] = periph
        r['arrayEff_pct'] = (100.0 * core / macro) if macro > 0 else None
        r['areaFlag'] = 'OK' if (core > 0 and periph > 0 and macro > 0) else 'NEGATIVE'

    rows.sort(key=lambda r: (r['n'], r['rows'], r['cols'], r['scheme'] or ''))
    return rows

def cfg(r):
    return f"n={int(r['n'])}  {int(r['rows'])}x{int(r['cols'])}"


def rows_by_cfg(rows):
    """Tag each row with a printable configuration name, preserving sort order."""
    for r in rows:
        r['_cfg'] = f"n={int(r['n'])} {int(r['rows'])}x{int(r['cols'])}"
    return rows

# --------------------------------------------------------------------------- #
#  generic table writer: index columns + scaled value columns                 #
# --------------------------------------------------------------------------- #
IDX = [('cell', 'cell', 10), ('rdo', 'mode', 8), ('n', 'planes', 8),
       ('rows', 'rows', 8), ('cols', 'cols', 8), ('scheme', 'write scheme', 16)]

def table(wb, title, rows, cols, note=None, extra_idx=(), unit_groups=()):
    """cols: list of (csv_key, label, kind) where kind is a LADDERS key or 'raw'.

    unit_groups: lists of csv_keys that MUST share one unit.  A breakdown whose
    components print in fJ while its total prints in pJ reads as though the
    components are negligible -- which is exactly how the V/2 inhibit (85-96%
    of the core write energy) came to look like 1% of the cells term, and how
    the restore energy (aJ) and restore latency (ps) came to look like they had
    been dropped from totals they were actually inside.
    """
    ws = wb.create_sheet(title)
    idx = list(IDX) + list(extra_idx)
    r0 = 1
    if note:
        c = ws.cell(row=1, column=1, value=note)
        c.font = Font(italic=True, size=9, color='555555')
        ws.merge_cells(start_row=1, start_column=1,
                       end_row=1, end_column=max(4, len(idx) + len(cols)))
        r0 = 3

    labels = [l for _, l, _ in idx] + [l for _, l, _ in cols]

    # FIX 2 -- resolve the shared units first, so a breakdown and its total
    # always print in the same unit and can be added up on the page.
    kind_of = {k: kd for k, _l, kd in cols}
    shared = {}
    for g in unit_groups:
        kd = kind_of.get(g[0])
        if kd in (None, 'raw'):
            continue
        pooled = [r.get(k) for r in rows for k in g if kind_of.get(k) == kd]
        u = pick_unit(pooled, kd, ref_mode='max')
        for k in g:
            if kind_of.get(k) == kd:
                shared[k] = u

    units, scales = [''] * len(idx), [None] * len(idx)
    for key, _lab, kind in cols:
        if kind == 'raw':
            units.append('')
            scales.append(None)
        else:
            sc, un = shared.get(key) or pick_unit([r.get(key) for r in rows], kind)
            units.append(un)
            scales.append(sc)
    widths = [w for _, _, w in idx] + [max(11, len(l) + 2) for _, l, _ in cols]
    write_header(ws, r0, labels, units, widths,
                 freeze=ws.cell(row=r0 + 2, column=len(idx) + 1).coordinate)

    for i, r in enumerate(rows):
        rr = r0 + 2 + i
        for j, (key, _lab, _w) in enumerate(idx, start=1):
            v = r.get(key)
            put(ws, rr, j, int(v) if isinstance(v, float) and float(v).is_integer() else v,
                fill=SUB_FILL if j <= 6 else None)
        for j, ((key, _lab, kind), sc) in enumerate(zip(cols, scales[len(idx):]),
                                                    start=len(idx) + 1):
            v = r.get(key)
            if kind == 'raw':
                fill = None
                if isinstance(v, str):
                    if v in OK_TOKENS:
                        fill = PASS_FILL
                    elif v in WARN_TOKENS:
                        fill = WARN_FILL
                    elif v in BAD_TOKENS:
                        fill = FAIL_FILL
                # a negative area is a broken model, not a small number
                elif isinstance(v, float) and key.startswith('area') and v < 0:
                    fill = FAIL_FILL
                put(ws, rr, j, v, fill=fill)
            else:
                put(ws, rr, j, None if v is None else v / sc, fmt='0.00000')
    return ws

def block_table(wb, title, field, kind, note):
    """One sheet with runs down and the 15 peripheral blocks across."""
    def key(b):
        return f'{b}_{field}'
    ws = wb.create_sheet(title)
    c = ws.cell(row=1, column=1, value=note)
    c.font = Font(italic=True, size=9, color='555555')
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=18)

    allv = [r.get(key(b)) for r in ROWS for b in BLOCKS]
    sc, un = pick_unit(allv, kind)
    # FIX 10: for LATENCY the TOTAL column is a naive sum over blocks, not the
    # macro critical path -- 22546 ps vs the Summary's 19170 ps at 256x512x512,
    # because the decoders run in parallel and the Summary takes their max.
    # Energy totals DO match the Summary and need no caveat.
    tot_label = 'TOTAL (sum, NOT critical path)' if kind == 's' else 'TOTAL'
    labels = [l for _, l, _ in IDX] + BLOCKS + [tot_label]
    units = [''] * len(IDX) + [un] * (len(BLOCKS) + 1)
    widths = [w for _, _, w in IDX] + [16] * len(BLOCKS) + [14]
    write_header(ws, 3, labels, units, widths, freeze='G5')

    for i, r in enumerate(ROWS):
        rr = 5 + i
        for j, (k, _l, _w) in enumerate(IDX, start=1):
            v = r.get(k)
            put(ws, rr, j, int(v) if isinstance(v, float) and float(v).is_integer() else v,
                fill=SUB_FILL)
        tot = 0.0
        for j, b in enumerate(BLOCKS, start=len(IDX) + 1):
            v = r.get(key(b))
            if v is not None:
                tot += v
            cell = put(ws, rr, j, None if v is None else v / sc, fmt='0.0000')
            if isinstance(v, float) and v < 0:
                cell.fill = FAIL_FILL
        put(ws, rr, len(IDX) + len(BLOCKS) + 1, tot / sc, fmt='0.0000',
            fill=KEY_FILL, bold=True)

    # share-of-total block, so the dominant periphery is obvious at a glance
    base = 5 + len(ROWS) + 2
    ws.cell(row=base, column=1, value='share of total (%)').font = Font(bold=True, size=10)
    write_header(ws, base + 1, labels[:len(IDX)] + BLOCKS, None,
                 None, freeze=None)
    for i, r in enumerate(ROWS):
        rr = base + 2 + i
        for j, (k, _l, _w) in enumerate(IDX, start=1):
            v = r.get(k)
            put(ws, rr, j, int(v) if isinstance(v, float) and float(v).is_integer() else v,
                fill=SUB_FILL)
        tot = sum(r.get(key(b)) or 0.0 for b in BLOCKS)
        for j, b in enumerate(BLOCKS, start=len(IDX) + 1):
            v = r.get(key(b)) or 0.0
            put(ws, rr, j, (100.0 * v / tot) if tot else None, fmt='0.00')
    return ws

# --------------------------------------------------------------------------- #
#  S x n pivot sheets (square sweep)                                          #
# --------------------------------------------------------------------------- #
def derive(r):
    """Per-run derived quantities shared by the Grid sheets and the Area sheet."""
    bits = r['rows'] * r['cols'] * r['n']
    r['_bits'] = bits
    macro, core = r.get('areaMacro_um2'), r.get('areaCore_um2')
    ok = r.get('areaFlag') == 'OK'
    r['_density'] = (bits / macro) if (ok and macro) else None
    r['_overhead'] = (100.0 * (macro - core) / core) if (ok and core) else None
    t = r.get('eWriteTrue_J')
    r['_trueFJb'] = (t / r['cols'] * 1e15) if (t and r.get('cols')) else None
    return r

def grid_sheet(wb, title, rows, scheme, metrics, note, p1=0.5):
    """One S (down) x n (across) table per metric, for one write scheme.

    metrics: list of (key, label, scale, number_format, is_text).  A cell is
    left blank when that (S, n) was not run.  Non-square runs are listed after
    the tables so an EXTRA= run is not silently dropped.
    """
    ws = wb.create_sheet(title)
    ws.column_dimensions['A'].width = 16
    sq = [r for r in rows if r['rows'] == r['cols'] and r['scheme'] == scheme
          and (r.get('p1') is None or r['p1'] == p1)]
    other = [r for r in rows if r['rows'] != r['cols'] and r['scheme'] == scheme
             and (r.get('p1') is None or r['p1'] == p1)]
    sizes = sorted({int(r['rows']) for r in sq})
    ns = sorted({int(r['n']) for r in sq})
    for j in range(2, len(ns) + 3):
        ws.column_dimensions[get_column_letter(j)].width = 13
    c = ws.cell(row=1, column=1, value=f'{note}  Write scheme: {scheme}; p1 = q1 = {p1}. '
                'Rows = subarray size S (rows = cols = S); columns = plane count n. '
                'Blank = not run.')
    c.font = Font(italic=True, size=9, color='555555')
    c.alignment = Alignment(wrap_text=True, vertical='top')
    ws.merge_cells(start_row=1, start_column=1, end_row=2, end_column=max(6, len(ns) + 1))
    if not sq:
        ws.cell(row=4, column=1, value=f'no square runs for scheme "{scheme}"')
        return ws
    at = {(int(r['rows']), int(r['n'])): r for r in sq}

    rr = 4
    for key, label, scale, fmt, is_text in metrics:
        ws.cell(row=rr, column=1, value=label).font = Font(bold=True, size=11, color='1F3864')
        rr += 1
        write_header(ws, rr, ['S \\ n'] + [str(n) for n in ns])
        rr += 1
        for s in sizes:
            put(ws, rr, 1, s, fill=SUB_FILL, bold=True)
            for j, n in enumerate(ns, start=2):
                r = at.get((s, n))
                v = None if r is None else r.get(key)
                if is_text:
                    fill = None
                    if isinstance(v, str):
                        fill = (PASS_FILL if v in OK_TOKENS else WARN_FILL if v in WARN_TOKENS
                                else FAIL_FILL if v in BAD_TOKENS else None)
                    put(ws, rr, j, v, fill=fill)
                else:
                    if isinstance(v, float):
                        v = v * scale
                    cell = put(ws, rr, j, v, fmt=fmt)
                    if r is not None and key.startswith('area') and isinstance(v, float) and v < 0:
                        cell.fill = FAIL_FILL
                    if r is not None and key in ('_density', 'arrayEff_pct', '_overhead') \
                            and r.get('areaFlag') != 'OK':
                        cell.fill = FAIL_FILL
            rr += 1
        rr += 2

    if other:
        ws.cell(row=rr, column=1,
                value='non-square runs in this CSV (not in the grids above):').font = Font(
                    bold=True, size=10)
        rr += 1
        write_header(ws, rr, ['n', 'rows', 'cols'] + [m[1] for m in metrics])
        rr += 1
        for r in sorted(other, key=lambda x: (x['n'], x['rows'], x['cols'])):
            put(ws, rr, 1, int(r['n'])); put(ws, rr, 2, int(r['rows'])); put(ws, rr, 3, int(r['cols']))
            for j, (key, _l, scale, fmt, is_text) in enumerate(metrics, start=4):
                v = r.get(key)
                put(ws, rr, j, v if (is_text or not isinstance(v, float)) else v * scale,
                    fmt=None if is_text else fmt)
            rr += 1
    return ws

# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('csv')
    ap.add_argument('-o', '--out', default=None)
    ap.add_argument('--grid-scheme', default='vanilla',
                    choices=['vanilla', 'erase-program', 'differential'],
                    help='write scheme shown in the Grid_* pivot sheets (default vanilla)')
    a = ap.parse_args()
    out = a.out or os.path.splitext(os.path.basename(a.csv))[0] + '.xlsx'

    global ROWS
    ROWS = load(a.csv)
    for r in ROWS:
        derive(r)
    schemes_present = sorted({r['scheme'] for r in ROWS})
    grid_scheme = a.grid_scheme if a.grid_scheme in schemes_present else schemes_present[0]
    n_square = sum(1 for r in ROWS if r['rows'] == r['cols'])
    ncol = len(open(a.csv).readline().split(','))
    wb = openpyxl.Workbook()
    wb.remove(wb.active)

    nbad = sum(1 for r in ROWS if r.get('areaFlag') == 'NEGATIVE')

    # ---------------------------------------------------------------- README
    ws = wb.create_sheet('README')
    ws.column_dimensions['A'].width = 26
    ws.column_dimensions['B'].width = 108
    lines = [
        ('3D FeRAM -- read/write operator mode sweep', ''),
        ('', ''),
        ('source CSV', os.path.basename(a.csv)),
        ('runs', str(len(ROWS))),
        ('CSV columns', str(ncol)),
        ('operation', 'SINGLE-ROW read and SINGLE-ROW write. activityRowRead = activityRowWrite '
                      '= 1/numRow, so every peripheral block charges exactly one activation.'),
        ('periphery included', 'row/plane decoders, RBL/RSL/WBL/WWL/WPL/SSL switch matrices, '
                               'current sense amplifier, output latch, cell array.'),
        ('periphery EXCLUDED', 'all CIM hardware -- no multilevel WL drivers, no ADC/accumulation, '
                               'no adder tree, no shift-add, no CIM input DFF scan chain.'),
        ('', ''),
        ('tiers', ''),
        ('  core', 'the array plus the wordline/plateline/bitline drive that reaches it.'),
        ('  macro', 'core + the non-CIM readout periphery (sense amp, latch, decoders).'),
        ('', ''),
        ('area split', 'core array area = areaArray + areaStaircase. periphery area = macro - core '
                       '= decoders + switch matrices + sense amplifiers + output latch. '
                       'total area = macro. array efficiency = core / total.'),
        ('area caveat', 'NeuroSim has no guard for an array narrower than one peripheral cell row. '
                        'Below about 16 rows, four row-direction decoders each return roughly '
                        '-9.8e8 um^2 AND two switch matrices about -1.27e9 um^2 (six negative '
                        'blocks, not four) and the macro area goes negative. Those runs are marked '
                        'NEGATIVE in the areaFlag column and their area, periphery, efficiency '
                        'and density numbers are meaningless. Fix with a MIN_BLK clamp in '
                        'SubArray::CalculateArea before trusting them.'
                        + (f'  >>> {nbad} run(s) in THIS workbook are affected. <<<' if nbad else
                           '  No run in this workbook is affected.')),
        ('', ''),
        ('units', 'READ THE GREY SUB-HEADER. Every energy, latency, power, capacitance, charge '
                  'and current column is auto-scaled to its own unit. Standalone columns use the '
                  'median magnitude of that column; columns that belong to one breakdown now SHARE '
                  'a single unit chosen from the largest member, so a breakdown and its total can '
                  'be added up on the page. Before that fix the write cells term printed in fJ '
                  'next to a V/2 inhibit in pJ -- making a term worth 85-96% of the total look '
                  'like 1% -- and the restore printed in aJ next to fJ neighbours. Area columns '
                  'are already in um^2 and are left unscaled.'),
        ('', ''),
        ('sheets', ''),
        ('  Grid_*', f'S x n pivot tables for the SQUARE sweep (rows = cols = S), one table per '
                     f'metric, write scheme "{grid_scheme}" (change with --grid-scheme). '
                     f'{n_square} of {len(ROWS)} runs are square. These are the tables the '
                     'paper figures are drawn from; plot_rw_sweep.py reads the same CSV.'),
        ('  Summary', 'headline read/write energy and latency, the three-way area split '
                      '(core array / periphery / total), array efficiency, leakage and the '
                      'sense verdict.'),
        ('  Physics', 'C_FE, C_node, Qsw, bias, unselected-plate drive, the two read currents, '
                      'the NLS switching time and the pulse count that each run actually used. '
                      'Sense margin and yield sigma appear here if the run emitted them.'),
        ('  Read_Breakdown', 'switch / linear / transistor-DC / lines / restore.'),
        ('  Write_Breakdown', 'cells / inhibit / lines, plus the flip and drive fractions, '
                             'switches per write, writes to failure, the true cost, the '
                             'half-selected cell count and the V_half/V_coercive disturb ratio.'),
        ('  Write_Schemes', 'the three write schemes side by side per configuration, best of '
                            'three highlighted.'),
        ('  Latency_Breakdown', 'select / charge / settle / sense amp / latch / restore / write select.'),
        ('  Block_*', 'the 16 peripheral blocks, one column each, plus share-of-total. '
                      'The TOTAL column is a plain sum: for ENERGY it equals the Summary macro '
                      'number, for LATENCY it does NOT, because the macro latency is a critical '
                      'path and the decoders run in parallel.'),
        ('  Area', 'core, periphery and macro area, overhead, efficiency, bit density.'),
        ('  Scaling', 'A: headline metrics against plane count n, one table per geometry. '
                      'B: against subarray size S at fixed n (square runs).'),
        ('', ''),
        ('write schemes', ''),
        ('  erase-program', 'two unconditional phases with a global plate polarity in each. '
                            'switch = p1 + q1, drive = 2, blind. A cell that should stay 1 '
                            'switches TWICE.'),
        ('  vanilla', 'THE STANDARD 1T1C FeRAM WRITE. One plate pulse 0 -> Vw -> 0 with the bit '
                      'line held at data: the rising edge writes the cells whose BL is low, the '
                      'falling edge those whose BL is high. switch = Hamming, drive = 1, blind -- '
                      'it gets the Hamming switching count for free, with no read. '
                      '[NVDRAM IEDM 2023 Sec. IV]'),
        ('  differential', 'read-modify-write: drive only the bits that change. '
                           'switch = drive = Hamming, but owes a read. Buys disturb immunity, '
                           'not switching count.'),
        ('  all three', 'take 2 plate phases -- a shared plate still has to apply both '
                        'polarities, so a scheme skips CELLS, not PHASES.'),
        ('', ''),
        ('sense criterion', 'dVsense = Qsw/C_node > senseAmpResolution, against a reference column. '
                            'NOT Vint > Vth -- no ferroelectric memory is read single-ended '
                            'against an absolute threshold. dVsense is CLAMPED to the rail '
                            'headroom (readVoltage - pillar bias); a run whose unloaded swing '
                            'would have exceeded that reports PASS (clipped) in amber rather than '
                            'a physically impossible signal.'),
        ('disturb criterion', 'V_half = writeVoltage/inhibitDivider across tFE, against '
                              'V_coercive = ecFE*tFE. A ratio above 1 means the V/2 scheme is '
                              'WRITING the half-selected cells, not inhibiting them, and the '
                              '"inhibit" energy is then switching energy under a wrong name. '
                              'See the Write_Breakdown sheet and the disturb table sweep.sh '
                              'prints at the end of a run.'),
        ('plane decoder', 'a log2(n)-bit RowDecoder is instantiated in rw_main and appended to '
                          'the block list. Without it nothing in the model selects among the n '
                          'planes, and write latency, t_select and the periphery area were all '
                          'insensitive to n -- the last of those in the wrong direction.'),
        ('sense is worst case', 'C_node = n*C_FE + C_gate assumes every plate is DC-driven, so '
                                'every one of the n capacitors is an AC short to ground and loads '
                                'the sense node with a full C_FE. dV therefore falls as 1/n. '
                                'A segmented pillar (a select device every nSeg planes) makes the '
                                'load nearly independent of total n. If the run emitted nSeg, the '
                                'Physics sheet shows which model produced these numbers.'),
    ]
    for i, (k, v) in enumerate(lines, start=1):
        c = ws.cell(row=i, column=1, value=k)
        c.font = Font(bold=True, size=11 if i == 1 else 10,
                      color='1F3864' if i == 1 else '000000')
        c.alignment = Alignment(vertical='top')
        c2 = ws.cell(row=i, column=2, value=v)
        c2.alignment = Alignment(wrap_text=True, vertical='top')

    # ----------------------------------------------------------------- Grids
    grid_sheet(wb, 'Grid_Energy', ROWS, grid_scheme, [
        ('eReadMacro_fJb', 'read energy, macro (fJ/bit)', 1.0, '0.00', False),
        ('eWriteMacro_fJb', 'write energy, macro (fJ/bit)', 1.0, '0.00', False),
        ('_trueFJb', 'true write energy incl. owed read (fJ/bit)', 1.0, '0.00', False),
        ('eReadCore_J', 'read energy, core (pJ per row access)', 1e12, '0.0000', False),
        ('eWriteCore_J', 'write energy, core (pJ per row access)', 1e12, '0.0000', False),
    ], note='Energy per bit of one single-row access.')
    grid_sheet(wb, 'Grid_Latency', ROWS, grid_scheme, [
        ('tReadMacro_s', 'read latency, macro (ns)', 1e9, '0.000', False),
        ('tWriteMacro_s', 'write latency, macro (ns)', 1e9, '0.000', False),
        ('tWriteTrue_s', 'true write latency incl. owed read (ns)', 1e9, '0.000', False),
        ('t_select_s', 'row/plane select (ns)', 1e9, '0.000', False),
        ('t_settle_s', 'line settle (ns)', 1e9, '0.000', False),
    ], note='Macro critical-path latency of one single-row access.')
    grid_sheet(wb, 'Grid_Area', ROWS, grid_scheme, [
        ('areaMacro_um2', 'total (macro) area (um^2)', 1.0, '0.00', False),
        ('areaCore_um2', 'core array area (um^2)', 1.0, '0.00', False),
        ('_density', 'bit density (Mbit/mm^2, all n planes)', 1.0, '0.0', False),
        ('arrayEff_pct', 'array efficiency (%)', 1.0, '0.00', False),
        ('_overhead', 'periphery overhead (%)', 1.0, '0.00', False),
        ('leakMacro_W', 'macro leakage (nW)', 1e9, '0.00', False),
        ('areaFlag', 'area ok?', 1.0, None, True),
    ], note='Area, density and leakage per macro. Red = NeuroSim area model failed (see README).')
    grid_sheet(wb, 'Grid_Sense', ROWS, grid_scheme, [
        ('dVsense_V', 'sense signal dV (mV)', 1e3, '0.00', False),
        ('senseVerdict', 'sense verdict', 1.0, None, True),
        ('disturbRatio', 'V_half / V_coercive', 1.0, '0.000', False),
        ('nHalfSel', 'half-selected cells per write', 1.0, '0', False),
        ('switchesPerWrite', 'reversals per cell per write', 1.0, '0.000', False),
        ('writesToFailure', 'writes to failure', 1.0, '0.00E+00', False),
    ] + [m for m in [
        ('degenLoss_frac', "far-row '1' current loss (%)", 100.0, '0.00', False),
        ('rblDropRatio', 'RBL IR drop / (V_RBL - V_dsat)  (>1 = page read fails)', 1.0, '0.000', False),
        ('t_senseAmp_s', 'sense time develop + regen (ns)', 1e9, '0.000', False),
    ] if has(ROWS, m[0])],
       note='Read margin and write disturb. dV depends on n only (flat pillar), not on S. '
            'The far-row current loss and the RBL page-current drop are the column-path limits '
            '(2026-09-25).')

    # --------------------------------------------------------------- Summary
    table(wb, 'Summary', ROWS, [
        ('flipFrac', 'flip fraction', 'raw'),
        ('eReadCore_J', 'read energy, core', 'J'),
        ('eReadMacro_J', 'read energy, macro', 'J'),
        ('eReadMacro_fJb', 'read fJ/bit', 'raw'),
        ('eWriteCore_J', 'write energy, core', 'J'),
        ('eWriteMacro_J', 'write energy, macro', 'J'),
        ('eWriteMacro_fJb', 'write fJ/bit', 'raw'),
        ('tReadCore_s', 'read latency, core', 's'),
        ('tReadMacro_s', 'read latency, macro', 's'),
        ('tWriteCore_s', 'write latency, core', 's'),
        ('tWriteMacro_s', 'write latency, macro', 's'),
        ('leakMacro_W', 'macro leakage', 'W'),
        # ---- the three-way area split -------------------------------------
        ('areaCore_um2', 'core array area (um^2)', 'raw'),
        ('areaPeriph_um2', 'periphery area (um^2)', 'raw'),
        ('areaMacro_um2', 'total area (um^2)', 'raw'),
        ('arrayEff_pct', 'array efficiency (%)', 'raw'),
        ('areaFlag', 'area ok?', 'raw'),
        # -------------------------------------------------------------------
        ('dVsense_V', 'sense signal', 'V'),
        ('senseVerdict', 'sense', 'raw'),
    ], unit_groups=[['eReadCore_J', 'eReadMacro_J', 'eWriteCore_J', 'eWriteMacro_J'],
                    ['tReadCore_s', 'tReadMacro_s', 'tWriteCore_s', 'tWriteMacro_s']],
       note='One row per run. p1 = fraction of 1s already stored, q1 = fraction of 1s '
            'being written. fJ/bit is per column of the single activated row. '
            'AREA: core array = cell array + staircase; periphery = total - core = decoders, '
            'switch matrices, sense amplifiers and output latch; total = macro. '
            'A NEGATIVE areaFlag means NeuroSim s area model failed for that geometry '
            '(array narrower than one peripheral cell row) -- ignore its area, efficiency '
            'and density.')

    # --------------------------------------------------------------- Physics
    phys_cols = [
        ('capFE_F', 'C_FE per cap', 'F'),
        ('capNode_F', 'C_node (pillar)', 'F'),
        ('Qsw_C', 'Qsw = 2Pr*A', 'C'),
        ('dVsense_V', 'dV = Qsw/C_node', 'V'),
    ]
    # These four only exist once the statistical / segmented sense model is in
    # rw_main.cpp.  Older CSVs simply do not get the columns.
    phys_cols += [c for c in [
        ('senseMargin', 'margin  dV/Vmin', 'raw'),
        ('senseSigma', 'yield sigma', 'raw'),
        ('nSeg', 'planes per segment', 'raw'),
        ('nSegments', 'segments in pillar', 'raw'),
    ] if has(ROWS, c[0])]
    phys_cols += [
        ('vIntBias_V', 'pillar bias', 'V'),
        ('vUnsel_V', 'unselected plate drive', 'V'),
        ('iRead1_A', 'I(read 1)', 'A'),
        ('iRead0_A', 'I(read 0)', 'A'),
        ('tauSwitch_s', 'NLS switching time', 's'),
        ('numWritePulses', 'write pulses', 'raw'),
        ('senseVerdict', 'sense', 'raw'),
    ]
    # SENSE/LINE FIX diagnostics (only in CSVs made after 2026-09-25)
    phys_cols += [c for c in [
        ('capJunctionTr_F', 'S/D cap per transistor', 'F'),
        ('capWBLpar_F', 'plate strip parasitic (no cells)', 'F'),
        ('capWBLwire_F', 'plate strip wire-to-ground', 'F'),
        ('rPlateDrv_ohm', 'R_plate (resTg)', 'raw'),
        ('rPlateEff_ohm', 'R_plate NeuroSim timed', 'raw'),
        ('tauLatch_s', 'latch tau = C/gm', 's'),
        ('senseRefFrac', 'reference fraction', 'raw'),
        ('senseMode', 'sense mode (0 clamp, 1 integrate)', 'raw'),
        ('gClamp_S', 'clamp gm (S)', 'raw'),
        ('rRslPath_ohm', 'R of the column, far row (ohm)', 'raw'),
        ('vSourceIR_V', 'far-row source IR drop', 'V'),
        ('degenLoss_frac', "far-row '1' current loss (frac)", 'raw'),
        ('vDropRBL_V', 'RBL IR drop, page current', 'V'),
        ('rblDropRatio', 'RBL drop / (V_RBL - V_dsat)', 'raw'),
    ] if has(ROWS, c[0])]
    table(wb, 'Physics', ROWS, phys_cols,
          note='C_FE = eps0*epsFE*A/tFE, derived -- not hardcoded. C_node = n*C_FE + C_gate '
               '(flat pillar) or nSeg*C_FE + C_gate + segment load (segmented pillar). '
               'The unselected plates are DRIVEN to vUnsel (self-boosted), so the pillar sits at '
               'the bias point at any n and the unselected capacitors see ~0 V -- but a DC-driven '
               'plate is an AC short, so each unselected capacitor still loads the node with a '
               'full C_FE. That is what makes the flat-pillar dV fall as 1/n, and it is the '
               'worst case. '
               'I(read 1)/I(read 0) is capped by the subthreshold slope: 10^(dV/SS), and clamps '
               'at ION_sat once 0.5*dV exceeds SS*log10(ION/Ith). '
               'FIX 2: the two current columns now SHARE one unit -- they used to print in uA '
               'and nA side by side, so 10 next to 148.339 read backwards.',
          unit_groups=[['iRead1_A', 'iRead0_A'], ['capFE_F', 'capNode_F']])

    # ------------------------------------------------------------ breakdowns
    rd_cols = [
        ('eRd_cellSwitch_J', 'cell switching', 'J'),
        ('eRd_cellLinear_J', 'cell linear C*V^2', 'J'),
        ('eRd_transistorDC_J', 'read transistor DC', 'J'),
        ('eRd_lines_J', 'line charging (incl. unselected strips)', 'J'),
        ('eRd_unselStrips_J', '  of which (n-1) unselected strips to vUnsel', 'J'),
        ('eRd_restore_J', 'restore / writeback (amortised 1/refresh)', 'J'),
        ('eReadCore_J', 'core total', 'J'),
        ('eSenseAmp_J', 'sense amplifier', 'J'),
        ('eLatch_J', 'output latch', 'J'),
        ('eDecRd_J', 'decoders', 'J'),
        ('ePrecharge_J', 'precharge', 'J'),
        ('eReadMacro_J', 'macro total', 'J'),
    ]
    rd_cols = [c for c in rd_cols if c[0] != 'eRd_unselStrips_J' or has(ROWS, c[0])]
    table(wb, 'Read_Breakdown', ROWS, rd_cols,
          unit_groups=[[k for k, _l, _kd in rd_cols]],
       note='Only the switching term and the transistor DC term depend on the data. The linear '
            'C*V^2 term is paid by every cell in the row whatever it holds. Precharge is 0 for '
            '2T-nC: the column precharge is inside line charging and overlaps decode. '
            'LINE CHARGING (2026-09-25) = selected strip parasitics (its MFM caps are the cell '
            'terms, no longer counted twice) + the (n-1) unselected strips of the accessed row '
            'driven to vUnsel to bias the pillar node (wire-to-ground only: node and plates move '
            'together) + RBL bias + RSL precharge.  The unselected-strip term grows linearly with '
            'n and is the read analogue of the write s V/2 inhibit.  The transistor DC term is '
            'V_RBL * I * (sense window) and the window is now the physical develop + regenerate '
            'time, not the 1 ns NeuroSim floor. '
            'RESTORE: this is a FULL writeback amortised over param->qndroRefreshInterval reads '
            '(1e6), so it is ~5e-18 J against a ~2.5e-13 J read -- five orders down, and it IS '
            'inside the core total. The UNAMORTISED cost of one restore is this number times the '
            'refresh interval; quote both in the paper, because the amortised figure on its own '
            'makes a destructive read look free.')

    wr_extra = [('eWr_unselStrips_J', '  of which (n-1) unselected strips at Vw/2', 'J')] \
               if has(ROWS, 'eWr_unselStrips_J') else []
    table(wb, 'Write_Breakdown', ROWS, [
        ('flipFrac', 'flip fraction', 'raw'),
        ('eWr_cells_J', 'cells (Qsw*V + C*V^2)', 'J'),
        ('eWr_inhibit_J', 'V/2 inhibit on unselected', 'J'),
        ('eWr_lines_J', 'line charging (incl. unselected strips)', 'J')] + wr_extra + [
        ('eWriteCore_J', 'core total', 'J'),
        ('eDecWr_J', 'decoders', 'J'),
        ('eWriteMacro_J', 'macro total', 'J'),
        ('numWritePulses', 'plate phases', 'raw'),
        ('driveFrac', 'field applications per cell', 'raw'),
        ('switchesPerWrite', 'reversals per cell per write', 'raw'),
        ('writesToFailure', 'writes to failure', 'raw'),
        ('eWriteTrue_J', 'true cost (write + any owed read)', 'J'),
        ('tWriteTrue_s', 'true latency (write + any owed read)', 's'),
        # FIX 1 / FIX 7 -- these four were either missing or mislabelled.
        ('nHalfSel', 'half-selected cells', 'raw'),
        ('vHalfSel_V', 'half-select bias', 'V'),
        ('eHalfSel_MVcm', 'half-select field (MV/cm)', 'raw'),
        ('disturbRatio', 'V_half / V_coercive', 'raw'),
    ], unit_groups=[['eWr_cells_J', 'eWr_inhibit_J', 'eWr_lines_J', 'eWriteCore_J',
                     'eDecWr_J', 'eWriteMacro_J', 'eWriteTrue_J'] + [c[0] for c in wr_extra]],
       note='THREE schemes. erase-program: switch = p1+q1, drive = 2 (every cell, both phases), '
            'blind. vanilla (the standard 1T1C single plate cycle): switch = Hamming, drive = 1, '
            'blind -- it gets the Hamming switching count for free. differential RMW: '
            'switch = drive = Hamming, but owes a read. All three take 2 plate phases, because a '
            'shared plate still has to apply both polarities. '
            'INHIBIT: V/2 self-boosted on nHalfSel = numCol*(n-1) capacitors -- row-independent, '
            'because a deselected pillar\'s select transistor is off. It is 85-96% of the core '
            'write energy, so read it against the cells term in the SAME unit (they now share one). '
            'DISTURB: V_half/V_coercive > 1 means the scheme is not inhibiting, it is writing '
            'those cells. V/3 was measured to fail on the vertical stack, so the fix is a lower '
            'write voltage or a segmented pillar, not a bigger divider.')

    # SENSE/LINE FIX: the new diagnostics only exist in CSVs made after the
    # 2026-09-25 rw_main; older CSVs simply do not get those columns.
    lat_cols = [
        ('t_select_s', 'row/plane select', 's'),
        ('t_charge_s', 'plate charge (FE floor + strip RC)', 's'),
        ('t_plateSupply_s', 'plate-driver supply bound', 's'),
        ('t_settle_s', 'line settle (RBL far end + column)', 's'),
        ('t_settleCol_s', 'column settle, worst-case row', 's'),
        ('capColSense_F', 'C_col (source line)', 'F'),
        ('t_senseAmp_s', 'sense amplifier (develop + regen)', 's'),
        ('t_senseDev_s', 'develop 60 mV: C_int*Vmin/(f*dI)', 's'),
        ('t_senseRegen_s', 'latch regeneration', 's'),
        ('t_senseNeuroSim_s', 'NeuroSim CSA floor (not used)', 's'),
        ('capSenseNode_F', 'C_int (column + amp input)', 'F'),
        ('iDiff_A', 'I1 - I0 - Ileak', 'A'),
        ('iLeakRows_A', 'row leakage', 'A'),
        ('readOv_V', 'read overdrive', 'V'),
        ('t_latch_s', 'output latch', 's'),
        ('t_restore_s', 'restore (amortised 1/refresh)', 's'),
        ('tReadMacro_s', 'read macro total', 's'),
        ('t_selectWrite_s', 'write select', 's'),
        ('tauSwitch_s', 'NLS tau (cross-check only)', 's'),
        ('tWriteMacro_s', 'write macro total', 's'),
    ]
    lat_cols = [c for c in lat_cols if c[0] in ('capSenseNode_F', 'iDiff_A', 'iLeakRows_A', 'readOv_V')
                or has(ROWS, c[0])]
    lat_units = [[k for k, _l, kd in lat_cols if kd == 's'],
                 # row leakage is SUBTRACTED from the usable difference, so the
                 # two must be comparable without a unit conversion in the head
                 ['iDiff_A', 'iLeakRows_A']]
    if has(ROWS, 'capColSense_F'):
        lat_units.append(['capSenseNode_F', 'capColSense_F'])
    table(wb, 'Latency_Breakdown', ROWS, lat_cols,
          unit_groups=lat_units,
       note='The write pulse is param->writePulseWidth (a flat 10 ns), NOT the NLS tau -- '
            'tauSwitch is printed as a cross-check and is not used by any latency here. '
            'Read macro = select -> charge -> settle -> sense amp -> latch -> restore, in series, '
            'and every term above prints in ONE shared unit so the chain adds up on the page. '
            'row/plane select is max(decoder delays), so a 512-row array and a 512-column array '
            'give the same number: same 9-bit decoder, either orientation. '
            'plate charge = MAX(param->chargeDelay, plate-driver supply bound) + 2.2*resRow*capWBL: '
            'the 15 ns ferroelectric floor plus the strip RC (kLineRc*R*C, distributed 90%).  The old '
            '2.2*R_PLATE*capSense term (the whole n dependence, up to 9.4 ns at n = 256) is GONE -- '
            'the pillar node is not charged through the plate driver on a read.  '
            'column settle (senseMode 0, clamped source line) = 3*(C_col/g_clamp + 0.5*R*C) for the '
            'far row; (senseMode 1, column integrates) = 0.5*R*C.  sense amplifier = time for the '
            'usable current 0.5*(I1-I0) (mid-point reference, far-row IR degeneration included in '
            'mode 0) to develop senseAmpVmin (60 mV) on the integrating node, plus the latch '
            'regeneration tau*ln(Vdd/2/Vmin).  Both the column settle and the development are '
            'LINEAR in the row count; NeuroSim s flat 1 ns is printed for reference only.')

    # ---------------------------------------------------------------- blocks
    block_table(wb, 'Block_Read_Energy', 'eRead_J', 'J',
                'Per-block read energy for ONE single-row read. Each block charges exactly one '
                'activation because activityRowRead = 1/numRow.')
    block_table(wb, 'Block_Write_Energy', 'eWrite_J', 'J',
                'Per-block write energy for ONE single-row write.')
    block_table(wb, 'Block_Read_Latency', 'tRead_s', 's',
                'Per-block read latency. These are the block critical paths, not a sum -- the '
                'macro latency is the serial chain select -> charge -> settle -> sense -> latch.')
    block_table(wb, 'Block_Write_Latency', 'tWrite_s', 's',
                'Per-block write latency.')
    block_table(wb, 'Block_Leakage', 'leak_W', 'W',
                'Per-block static leakage power.')
    block_table(wb, 'Block_Area', 'area_m2', 'm2',
                'Per-block layout area. Cells shaded red are NEGATIVE -- NeuroSim s area model '
                'has no guard when the array is narrower than one peripheral cell row, which is '
                'what breaks the small-row configurations.')

    # --------------------------------------------------------- Write_Schemes
    ws = wb.create_sheet('Write_Schemes')
    ws.column_dimensions['A'].width = 12
    for j in range(2, 14):
        ws.column_dimensions[get_column_letter(j)].width = 14
    c = ws.cell(row=1, column=1,
                value='The three write schemes side by side. ENERGY barely separates them -- the '
                      'V/2 half-select dominates the write and all three take 2 plate phases. '
                      'ENDURANCE is where they differ: erase-then-program switches 1.0 per write '
                      'however little the data changes, while vanilla gets the Hamming count for '
                      'free (no read) and differential pays a read for the same count.')
    c.font = Font(italic=True, size=9, color='555555')
    c.alignment = Alignment(wrap_text=True, vertical='top')
    ws.merge_cells(start_row=1, start_column=1, end_row=2, end_column=13)

    ORDER = ['erase-program', 'vanilla', 'differential']
    SHORT = ['erase-prog', 'vanilla', 'differential']
    by_cfg = OrderedDict()
    for r in rows_by_cfg(ROWS):
        by_cfg.setdefault(r['_cfg'], {})[r['scheme']] = r

    groups = [('write fJ/bit', 'eWriteMacro_fJb', '0.00'),
              ('switches per write', 'switchesPerWrite', '0.000'),
              ('writes to failure', 'writesToFailure', '0.00E+00'),
              ('true write fJ/bit (incl. owed read)', '_trueFJb', '0.00')]

    r0 = 4
    ws.cell(row=r0, column=1, value='configuration').fill = HDR_FILL
    ws.cell(row=r0, column=1).font = HDR_FONT
    col = 2
    for gname, _key, _fmt in groups:
        ws.merge_cells(start_row=r0, start_column=col, end_row=r0, end_column=col + 2)
        cc = ws.cell(row=r0, column=col, value=gname)
        cc.fill, cc.font = HDR_FILL, HDR_FONT
        cc.alignment = Alignment(horizontal='center')
        for k, sname in enumerate(SHORT):
            sc = ws.cell(row=r0 + 1, column=col + k, value=sname)
            sc.fill, sc.font, sc.border = SUB_FILL, SUB_FONT, BOX
            sc.alignment = Alignment(horizontal='center')
        col += 3
    ws.freeze_panes = 'B6'

    rr = r0 + 2
    for cfg_name, per in by_cfg.items():
        put(ws, rr, 1, cfg_name, fill=SUB_FILL, bold=True)
        col = 2
        for _g, key, fmt in groups:
            best, vals = None, []
            for sch in ORDER:
                r = per.get(sch)
                v = None
                if r is not None:
                    v = r.get(key)          # _trueFJb is filled by derive()
                vals.append(v)
            good = [v for v in vals if isinstance(v, float)]
            if good:
                best = max(good) if key == 'writesToFailure' else min(good)
            for k, v in enumerate(vals):
                cell = put(ws, rr, col + k, v, fmt=fmt)
                if best is not None and isinstance(v, float) and abs(v - best) < 1e-12:
                    cell.fill = PASS_FILL
            col += 3
        rr += 1
    ws.cell(row=rr + 1, column=1,
            value='green = best of the three for that metric').font = Font(
                italic=True, size=9, color='555555')

    # ------------------------------------------------------------------ Area
    # (_bits, _density, _overhead were filled by derive() in main)
    table(wb, 'Area', ROWS, [
        ('_bits', 'bits in macro', 'raw'),
        ('areaCore_um2', 'core array area (um^2)', 'raw'),
        ('areaPeriph_um2', 'periphery area (um^2)', 'raw'),
        ('areaMacro_um2', 'total area (um^2)', 'raw'),
        ('arrayEff_pct', 'array efficiency (%)', 'raw'),
        ('_overhead', 'periphery overhead (%)', 'raw'),
        ('_density', 'density (bit/um^2 = Mbit/mm^2)', 'raw'),
        ('leakMacro_W', 'macro leakage', 'W'),
        ('areaFlag', 'area ok?', 'raw'),
    ], note='bit/um^2 and Mbit/mm^2 are the same number. Density counts all n planes. '
            'array efficiency = core / total; periphery overhead = periphery / core. '
            'Density, efficiency and overhead are left BLANK for any run whose areaFlag is '
            'NEGATIVE, because NeuroSim s area model has no guard for an array narrower than '
            'one peripheral cell row and returns large negative decoder areas there.')

    # --------------------------------------------------------------- Scaling
    ws = wb.create_sheet('Scaling')
    ws.column_dimensions['A'].width = 18
    for j in range(2, 12):
        ws.column_dimensions[get_column_letter(j)].width = 15
    c = ws.cell(row=1, column=1,
                value='Part A: energy per bit against plane count n, one table per geometry. '
                      'Part B: against subarray size S at fixed n (square runs only). '
                      'Mechanisms to read the tables against: the read cost per bit rises with n '
                      'because the read transistor is biased so the "1" state sits at I_ON and '
                      'I(read 0) climbs as the sense signal shrinks (1/n on a flat pillar), so '
                      'the read-transistor DC term grows and outweighs amortising the line and '
                      'decoder cost over more stored bits. The write cost per bit rises faster '
                      'with n, for a different reason: the V/2 inhibit charges numCol*(n-1) '
                      'capacitors and dominates the core write. Against S, everything a '
                      'single-row access pays scales with the column count and divides back out '
                      'of fJ/bit; what survives is the half-select population and the periphery '
                      'amortisation, so density improves with S while fJ/bit moves only weakly.')
    c.font = Font(italic=True, size=9, color='555555')
    c.alignment = Alignment(wrap_text=True, vertical='top')
    ws.merge_cells(start_row=1, start_column=1, end_row=3, end_column=11)

    SCAL_HDR = ['write scheme', 'read fJ/bit', 'write fJ/bit', 'read latency (ns)',
                'write latency (ns)', 'dV (mV)', 'sense', 'total area (um^2)',
                'array eff (%)', 'density (Mbit/mm^2)']

    def scal_row(ws, rr, r, first):
        put(ws, rr, 1, first, fill=SUB_FILL)
        put(ws, rr, 2, r['scheme'], fill=SUB_FILL)
        put(ws, rr, 3, r['eReadMacro_fJb'], fmt='0.000')
        put(ws, rr, 4, r['eWriteMacro_fJb'], fmt='0.000')
        put(ws, rr, 5, (r['tReadMacro_s'] or 0) * 1e9, fmt='0.000')
        put(ws, rr, 6, (r['tWriteMacro_s'] or 0) * 1e9, fmt='0.000')
        put(ws, rr, 7, (r['dVsense_V'] or 0) * 1e3, fmt='0.00')
        sv = r['senseVerdict']
        put(ws, rr, 8, sv, fill=PASS_FILL if sv in OK_TOKENS else
            WARN_FILL if sv in WARN_TOKENS else FAIL_FILL)
        ok = r.get('areaFlag') == 'OK'
        put(ws, rr, 9, r.get('areaMacro_um2'), fmt='0.00', fill=None if ok else FAIL_FILL)
        put(ws, rr, 10, r.get('arrayEff_pct') if ok else None, fmt='0.00')
        put(ws, rr, 11, r.get('_density') if ok else None, fmt='0.0')

    rr = 5
    ws.cell(row=rr, column=1, value='A. plane-count scaling, one table per geometry').font = Font(
        bold=True, size=12, color='1F3864')
    rr += 2
    geoms = OrderedDict()
    for r in ROWS:
        geoms.setdefault((int(r['rows']), int(r['cols'])), []).append(r)
    for (rows_, cols_), grp in geoms.items():
        ws.cell(row=rr, column=1, value=f'{rows_} rows x {cols_} cols').font = Font(
            bold=True, size=11, color='1F3864')
        rr += 1
        write_header(ws, rr, ['n (planes)'] + SCAL_HDR)
        rr += 1
        for r in sorted(grp, key=lambda x: (x['scheme'], x['n'])):
            scal_row(ws, rr, r, int(r['n']))
            rr += 1
        rr += 2

    by_n = OrderedDict()
    for r in ROWS:
        if r['rows'] == r['cols']:
            by_n.setdefault(int(r['n']), []).append(r)
    if by_n:
        ws.cell(row=rr, column=1, value='B. subarray-size scaling (square runs), one table per n').font = Font(
            bold=True, size=12, color='1F3864')
        rr += 2
        for n, grp in by_n.items():
            ws.cell(row=rr, column=1, value=f'n = {n} planes').font = Font(
                bold=True, size=11, color='1F3864')
            rr += 1
            write_header(ws, rr, ['S (rows = cols)'] + SCAL_HDR)
            rr += 1
            for r in sorted(grp, key=lambda x: (x['scheme'], x['rows'])):
                scal_row(ws, rr, r, int(r['rows']))
                rr += 1
            rr += 2

    wb.save(out)
    print(f"{out}   {len(wb.sheetnames)} sheets, {len(ROWS)} runs, {ncol} CSV columns")
    if nbad:
        print(f"    WARNING: {nbad} run(s) have a NEGATIVE macro area -- NeuroSim's area model "
              f"failed for that geometry. See the README sheet.")
    for s in wb.sheetnames:
        print(f"    {s}")

if __name__ == '__main__':
    main()
