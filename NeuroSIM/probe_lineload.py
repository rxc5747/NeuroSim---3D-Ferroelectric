#!/usr/bin/env python3
"""
probe_lineload.py -- read-only. Prints, for SubArray.cpp:
  (a) the whole line-capacitance / geometry block, and
  (b) every Initialize / CalculateLatency / CalculatePower call on the six line
      drivers and their decoders, WITH the full argument list (calls that wrap
      over several lines are joined), so the LOAD each driver is given is visible.

    python3 probe_lineload.py SubArray.cpp
"""
import re, sys

SRC = sys.argv[1] if len(sys.argv) > 1 else 'SubArray.cpp'
L = open(SRC, errors='replace').read().split('\n')

BLOCKS = ['wblSwitchMatrix', 'wwlSwitchMatrix', 'wplSwitchMatrix', 'sslSwitchMatrix',
          'rblSwitchMatrix', 'rslSwitchMatrix',
          'wblDecoder', 'wwlDecoder', 'wplDecoder', 'sslDecoder',
          'rblDecoder', 'rslDecoder']
CAPS = ['capWBL', 'capRBL', 'capRSL', 'capWWL', 'capWPL', 'capSSL',
        'capWBLTotal', 'totalCapCrossWBL', 'totalCapCrossRBL',
        'resWBL', 'resRBL', 'resRSL', 'resWWL', 'resWPL', 'resSSL',
        'resRow', 'resCol', 'lengthRow', 'lengthCol',
        'sourceCap', 'gateCap', 'DrainCap', 'drainCap', 'capCol', 'capRow']

def uncommented(i):
    return not L[i].lstrip().startswith('//')

def join_call(i):
    """Join from line i until parentheses balance, so multi-line calls print whole."""
    s, depth, started = '', 0, False
    for j in range(i, min(i + 25, len(L))):
        t = L[j].split('//')[0]
        s += (' ' if s else '') + t.strip()
        for ch in t:
            if ch == '(':
                depth += 1; started = True
            elif ch == ')':
                depth -= 1
        if started and depth <= 0:
            return re.sub(r'\s+', ' ', s), j
    return re.sub(r'\s+', ' ', s), i

print('=' * 100)
print(f'  {SRC}')
print('=' * 100)

print('\n########## A. every assignment to a line capacitance / resistance / geometry\n')
for i, ln in enumerate(L):
    if not uncommented(i):
        continue
    for c in CAPS:
        if re.search(r'\b' + c + r'\s*(=|\+=)[^=]', ln):
            print(f'{i+1:>6}: {ln.rstrip()}')
            break

print('\n########## B. the neighbourhood of the capacitance block (context)\n')
hits = [i for i, ln in enumerate(L) if uncommented(i)
        and re.search(r'\bcapW(BL|WL|PL)\s*=', ln)]
if hits:
    a, b = max(0, min(hits) - 30), min(len(L), max(hits) + 15)
    for i in range(a, b):
        print(f'{i+1:>6}: {L[i].rstrip()}')
else:
    print('  (no capWBL/capWWL/capWPL assignment found)')

print('\n########## C. every driver / decoder call, with its full argument list\n')
seen = set()
for i, ln in enumerate(L):
    if not uncommented(i):
        continue
    for b in BLOCKS:
        m = re.search(r'\b' + b + r'\s*\.\s*(Initialize|CalculateArea|CalculateLatency'
                      r'|CalculatePower)\s*\(', ln)
        if m and (i, b) not in seen:
            seen.add((i, b))
            call, end = join_call(i)
            tag = f'{b}.{m.group(1)}'
            print(f'{i+1:>6}: [{tag}]')
            print(f'        {call[:400]}')
print('\n########## D. what to send back')
print('  Paste sections A, B and C above. That is everything needed to see')
print('  which load each of the six line drivers is actually handed.')
