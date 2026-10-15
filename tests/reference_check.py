"""Independent SciPy comparison to JS vectors; not a runtime dependency.

Run the Node core suite first, then: python tests/reference_check.py
Uses the repo's existing Python paired/unpaired implementations.
"""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'python'))
from sensei import paired, unpaired

vectors = json.loads((ROOT / '.artifacts/reference-vectors.json').read_text())
largest = (0.0, None)
for i, case in enumerate(vectors):
    p = case['p']
    def raw(offset, sign, test_type):
        args = ([p['N0'], p['N1']], [p['a0'], p['a1']], [p['b0'], p['b1']])
        kwargs = dict(alpha=p['alpha'], test_type=test_type.replace('_', '-'), offset=offset, sign=sign)
        if p['paired']:
            return float(paired.calc_fn_rate_beta(case['m0'], *args, rho=p['rho'], **kwargs))
        return float(unpaired.calc_fn_rate_beta([case['m0'], case['m1']], *args, **kwargs))
    if p['test'] == 'tost_test':
        expected = min(raw(p['upper'], -1, 'one_sided') + raw(p['lower'], 1, 'one_sided'), 1)
    elif p['test'] == 'sesoi_test':
        sign = (1 if p['offset'] > 0 else -1) * (1 if p['verify'] else -1)
        expected = raw(p['offset'], sign, 'one_sided')
    else:
        expected = raw(0, 0, p['test'])
    error = abs(expected - case['actual'])
    assert math.isfinite(expected) and error < 2e-7, (i, error, case, expected)
    if error > largest[0]:
        largest = (error, i)
print(f'PASS: {len(vectors)} independent SciPy comparisons; maximum absolute error {largest[0]:.3g} (vector {largest[1]}).')
