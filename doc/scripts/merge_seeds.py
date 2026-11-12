"""Merge runs of one setting with different seeds (equal numbers of events) into one output directory.

  python3 merge_seeds.py <outdir> <rundir1> <rundir2> ...

Every histogram is per event, so the merged value is the mean over the runs and the
error is sqrt(sum of squared errors)/N.  Histogram files have rows "x value error"
(or "x value"); pairmap_* files have rows "|p1+p2| k value".  summary.dat rows
"key mean error" are merged the same way and the event count in its header is the sum.
Other lines (comments, Sherpa's histogram header) are taken from the first run.  Files whose
rows differ between the runs (multiplicity_*, filled per multiplicity class) are skipped."""
import os, re, sys

def numbers(line):
    try: return [float(v) for v in line.split()]
    except ValueError: return None

out, runs = sys.argv[1], sys.argv[2:]
os.makedirs(out, exist_ok=True)
n = len(runs)
for name in sorted(os.listdir(runs[0])):
    if not name.endswith('.dat'): continue
    files = [open(os.path.join(r, name)).read().splitlines() for r in runs]
    if name == 'summary.dat':
        # rows "key mean error"; a key missing in a run (a rare count that stayed 0) counts as 0
        head = [l for l in files[0] if l.startswith('#')]
        total = sum(float(re.search(r'events = ([0-9.e+]+)', f[0]).group(1)) for f in files)
        head[0] = re.sub(r'events = [0-9.e+]+', f'events = {total:g}', head[0])
        keys, vals = [], {}
        for f in files:
            for l in f:
                if l.startswith('#') or not l.strip(): continue
                k, m, e = l.split()
                if k not in vals: keys.append(k); vals[k] = [0., 0.]
                vals[k][0] += float(m)/n; vals[k][1] += float(e)**2
        rows = [f'{k}  {vals[k][0]:g}  {vals[k][1]**0.5/n:g}' for k in keys]
        open(os.path.join(out, name), 'w').write('\n'.join(head + rows) + '\n')
        continue
    if any(len(f) != len(files[0]) for f in files):
        print('skipped (rows differ between runs):', name); continue
    merged = []
    for i, line in enumerate(files[0]):
        rows = [numbers(f[i]) for f in files]
        if line.startswith('#') or rows[0] is None or len(rows[0]) not in (2, 3):
            merged.append(line); continue
        if name.startswith('pairmap_'):
            merged.append(f'{rows[0][0]:g} {rows[0][1]:g} {sum(r[2] for r in rows)/n:g}')
        elif len(rows[0]) == 2:
            merged.append(f'{rows[0][0]:g}  {sum(r[1] for r in rows)/n:g}')
        else:
            merged.append(f'{rows[0][0]:g}  {sum(r[1] for r in rows)/n:g}  {sum(r[2]**2 for r in rows)**0.5/n:g}')
    open(os.path.join(out, name), 'w').write('\n'.join(merged) + '\n')
print(f'merged {n} runs into {out}')
