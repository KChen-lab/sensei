/* Original Sensei web arithmetic; preserve operation order and parameter rounding.
 * See table_beta_legacy.html for the independent, original implementation. */
(function (root) {
  'use strict';
  const stats = typeof module === 'object' && module.exports ? require('./vendor/jstat.min.js').jStat : root.jStat;
  function round(x, p) { return Math.round(x * Math.pow(10, p)) / Math.pow(10, p); }
  function betaStatParam(m, std) {
    m = parseFloat(m);
    var v = parseFloat(std);
    v = v * v;
    var a = ((1 - m) / v - 1 / m) * m * m;
    var b = a * (1 / m - 1);
    return [a, b];
  }
  function betaParamStat(a, b) {
    a = parseFloat(a); b = parseFloat(b);
    var m = a / (a + b);
    var v = a * b / Math.pow(a + b, 2) / (a + b + 1);
    return [m, Math.sqrt(v)];
  }
  function calcBetaRaw(M0, M1, N0, N1, a0, a1, b0, b1, alpha, test_type, offset, sign) {
    var Ep0 = a0 / (a0 + b0), Ep1 = a1 / (a1 + b1);
    var Vp0 = a0 * b0 * (a0 + b0 + N0) / ((N0 * (a0 + b0) * (a0 + b0)) * (a0 + b0 + 1));
    var Vp1 = a1 * b1 * (a1 + b1 + N1) / ((N1 * (a1 + b1) * (a1 + b1)) * (a1 + b1 + 1));
    var Et = (Ep1 - (Ep0 + offset)) / Math.sqrt(Vp0 / M0 + Vp1 / M1);
    if (sign == 0) Et = Math.abs(Et); else Et = sign * Et;
    var nu = (Vp0 / M0 + Vp1 / M1) ** 2 / ((Vp0 / M0) ** 2 / (M0 - 1) + (Vp1 / M1) ** 2 / (M1 - 1));
    var t_star = stats.studentt.inv(1 - (test_type === 'one_sided' ? alpha : alpha / 2), nu);
    return stats.studentt.cdf(t_star - Et, nu);
  }
  function calcBetaRawPaired(M, N0, N1, a0, a1, b0, b1, alpha, rho, test_type, offset, sign) {
    var Ep0 = a0 / (a0 + b0), Ep1 = a1 / (a1 + b1);
    var Vp0 = a0 * b0 * (a0 + b0 + N0) / ((N0 * (a0 + b0) * (a0 + b0)) * (a0 + b0 + 1));
    var Vp1 = a1 * b1 * (a1 + b1 + N1) / ((N1 * (a1 + b1) * (a1 + b1)) * (a1 + b1 + 1));
    var cov = rho * Math.sqrt(a0 * b0 / ((a0 + b0) * (a0 + b0) * (a0 + b0 + 1))) * Math.sqrt(a1 * b1 / ((a1 + b1) * (a1 + b1) * (a1 + b1 + 1)));
    var Et = (Ep1 - (Ep0 + offset)) * Math.sqrt(M) / Math.sqrt(Vp0 + Vp1 - 2 * cov);
    if (sign == 0) Et = Math.abs(Et); else Et = sign * Et;
    var nu = M - 1;
    var t_star = stats.studentt.inv(1 - (test_type === 'one_sided' ? alpha : alpha / 2), nu);
    return stats.studentt.cdf(t_star - Et, nu);
  }
  function falseNegative(p, m0, m1) {
    const raw = (offset, sign, type) => p.paired
      ? calcBetaRawPaired(m0, p.N0, p.N1, p.a0, p.a1, p.b0, p.b1, p.alpha, p.rho, type, offset, sign)
      : calcBetaRaw(m0, m1, p.N0, p.N1, p.a0, p.a1, p.b0, p.b1, p.alpha, type, offset, sign);
    if (p.test === 'tost_test') return Math.min(raw(p.upper, -1, 'one_sided') + raw(p.lower, 1, 'one_sided'), 1);
    if (p.test === 'sesoi_test') return raw(p.offset, (p.verify ? 1 : -1) * Math.sign(p.offset), 'one_sided');
    return raw(0, 0, p.test);
  }
  const api = { round, betaStatParam, betaParamStat, calcBetaRaw, calcBetaRawPaired, falseNegative };
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.Sensei = api;
})(typeof globalThis === 'object' ? globalThis : this);
