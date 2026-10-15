const { test } = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const { createHash } = require('node:crypto');
const core = require('../assets/sensei-core.js');
const jStat = require('../assets/vendor/jstat.min.js').jStat;
const root = path.resolve(__dirname, '..');
const original = fs.readFileSync(path.join(root, 'table_beta_legacy.html'), 'utf8');
const source = original.match(/<script type="text\/javascript">([\s\S]*?)<\/script>/)[1];
const legacy = vm.createContext({ jStat, console: { log() {} } });
vm.runInContext(source, legacy);
test('legacy calculation and event scripts are unchanged from the original checkout', () => {
  const scripts = Array.from(original.replace(/\r\n/g, '\n').matchAll(/<script(?: type="text\/javascript")?>([\s\S]*?)<\/script>/g), m => m[1]);
  assert.deepEqual(scripts.map(s => createHash('sha256').update(s).digest('hex')), [
    'ea481396cff54b3333fbc58fa92c1dfc7dbc3d5109519f59c612ed0d12b564e2',
    'c6583509315c28810b32828eecf80d2c57ca15d7cc2128cfa1cf95b46d2aae06'
  ]);
});
let seed = 85491;
const random = () => { seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0; return seed / 2 ** 32; };
const vectors = [];
function expected(p, m0, m1) {
  const raw = (offset, sign, type) => p.paired
    ? legacy.calcBetaRawPaired(m0, p.N0, p.N1, p.a0, p.a1, p.b0, p.b1, p.alpha, p.rho, type, offset, sign)
    : legacy.calcBetaRaw(m0, m1, p.N0, p.N1, p.a0, p.a1, p.b0, p.b1, p.alpha, type, offset, sign);
  if (p.test === 'tost_test') return Math.min(raw(p.upper, -1, 'one_sided') + raw(p.lower, 1, 'one_sided'), 1);
  if (p.test === 'sesoi_test') return raw(p.offset, (p.verify ? 1 : -1) * Math.sign(p.offset), 'one_sided');
  return raw(0, 0, p.test);
}
test('10,000 full-precision calculations match the original implementation exactly', () => {
  for (let i = 0; i < 10000; i++) {
    const p = { paired: i % 2 === 0, N0: [1, 100, 384, 1000, 100000][i % 5], N1: 1 + Math.floor(random() * 10000), a0: 1.01 + random() * 100, a1: 1.01 + random() * 100, b0: 1.01 + random() * 300, b1: 1.01 + random() * 300, alpha: [0.05, 0.01, 0.1][i % 3], rho: [-1, -0.5, 0, 0.5, 1][i % 5], test: ['one_sided', 'two_sided', 'sesoi_test', 'tost_test'][Math.floor(i / 2) % 4], offset: i % 3 ? 0.01 : -0.05, verify: i % 3 !== 0, lower: -0.1, upper: 0.1 };
    const m0 = 2 + Math.floor(random() * 199), m1 = 2 + Math.floor(random() * 199);
    const actual = core.falseNegative(p, m0, m1), old = expected(p, m0, m1);
    assert.equal(actual, old, JSON.stringify({ i, p, m0, m1 }));
    assert.ok(Number.isFinite(actual) && actual >= 0 && actual <= 1);
    if (i % 50 === 0) vectors.push({ p, m0, m1, actual });
  }
});
test('mean/SD and beta conversions retain original precision and rounding', () => {
  for (const [mean, sd] of [[0.03, 0.015], [0.05, 0.01], [0.4, 0.05], [0.001, 0.0002], [0.99, 0.001], [0.5, 0.25]]) {
    assert.deepEqual(core.betaStatParam(mean, sd), Array.from(legacy.betaStatParam(mean, sd)));
    const [a, b] = core.betaStatParam(mean, sd).map(v => core.round(v, 3));
    assert.deepEqual(core.betaParamStat(a, b), Array.from(legacy.betaParamStat(a, b)));
    assert.equal(core.round(mean, 4), legacy.round(mean, 4));
  }
});
test('null effect, rare cells, SESOI directions, TOST cap, and correlation boundaries', () => {
  const base = { paired: false, N0: 1000, N1: 1000, a0: 7, b0: 10, a1: 10, b1: 7, alpha: 0.05, rho: 0.5, test: 'one_sided', offset: 0.01, verify: true, lower: -0.1, upper: 0.1 };
  const scenarios = [base, { ...base, a1: 7, b1: 10 }, { ...base, a0: 1.1, b0: 9999, a1: 2.1, b1: 9999 }, { ...base, a0: 0.5, b1: 0.5 }, { ...base, test: 'tost_test', lower: -0.0001, upper: 0.0001 }];
  for (const p of scenarios) for (const paired of [false, true]) for (const testType of ['one_sided', 'two_sided', 'sesoi_test', 'tost_test']) for (const verify of [false, true]) for (const rho of [-1, 0, 0.5, 1]) {
    const q = { ...p, paired, test: testType, verify, rho };
    for (const m0 of [2, 5, 12, 50, 200]) {
      const actual = core.falseNegative(q, m0, 13);
      assert.equal(actual, expected(q, m0, 13));
      vectors.push({ p: q, m0, m1: 13, actual });
    }
  }
  assert.ok(Math.abs(core.falseNegative({ ...base, a1: 7, b1: 10 }, 10, 10) - 0.95) < 1e-7);
  assert.equal(core.falseNegative(scenarios[4], 5, 5), 1);
  assert.ok(core.falseNegative({ ...base, paired: true, rho: 0.9 }, 10, 10) < core.falseNegative({ ...base, paired: true, rho: 0 }, 10, 10));
});
test('write vectors for independent SciPy validation and manual checks', () => {
  fs.mkdirSync(path.join(root, '.artifacts'), { recursive: true });
  fs.writeFileSync(path.join(root, '.artifacts/reference-vectors.json'), JSON.stringify(vectors));
  console.log(`${vectors.length} vectors available for tests/reference_check.py`);
  const p = { paired: false, N0: 1000, N1: 1000, a0: 7, b0: 10, a1: 10, b1: 7, alpha: 0.05, rho: 0.5, test: 'one_sided', offset: 0.01, verify: true, lower: -0.1, upper: 0.1 };
  console.log('Manual check values:', JSON.stringify([
    ['unpaired default M0=5,M1=5', core.round(core.falseNegative(p, 5, 5), 3)],
    ['unpaired default M0=12,M1=12', core.round(core.falseNegative(p, 12, 12), 3)],
    ['paired default M=5', core.round(core.falseNegative({ ...p, paired: true }, 5, 5), 3)],
    ['unpaired two-sided M0=5,M1=5', core.round(core.falseNegative({ ...p, test: 'two_sided' }, 5, 5), 3)],
    ['unpaired SESOI verify M0=5,M1=5', core.round(core.falseNegative({ ...p, test: 'sesoi_test' }, 5, 5), 3)],
    ['unpaired SESOI falsify M0=5,M1=5', core.round(core.falseNegative({ ...p, test: 'sesoi_test', verify: false }, 5, 5), 3)],
    ['unpaired TOST M0=5,M1=5', core.round(core.falseNegative({ ...p, test: 'tost_test' }, 5, 5), 3)]
  ]));
});
