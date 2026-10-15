/* Build-free, browser-only UI. Scientific arithmetic lives in sensei-core.js. */
(() => {
  'use strict';
  const $ = id => document.getElementById(id);
  const states = new Map();
  const tooltip = $('result-tooltip');
  let tooltipCell = null, tooltipHideTimer;
  const chartObserver = new ResizeObserver(entries => {
    for (const { target, contentRect } of entries) {
      if (target.data && contentRect.width > 0 && contentRect.height > 0) {
        Plotly.relayout(target, { width: Math.floor(contentRect.width), height: Math.floor(contentRect.height) });
      }
    }
  });
  const suffix = paired => paired ? '_paired' : '';
  function field(id, label, value, hint = '', attrs = '') {
    const step = /^N[01](?:_paired)?$/.test(id) ? '1' : '0.01';
    return `<div class="field"><label for="${id}"${hint ? ` class="help-label" title="${hint}"` : ''}>${label}</label><input id="${id}" type="number" step="${step}" value="${value}" ${attrs}></div>`;
  }
  function range(min, max, label) {
    return `<div class="field"><span class="field-label">${label}</span><div class="range"><input type="number" step="1" min="2" id="${min}" value="5" aria-label="${label}, minimum"><span>to</span><input type="number" step="1" min="2" id="${max}" value="12" aria-label="${label}, maximum"></div></div>`;
  }
  function group(g, s, paired) {
    const name = g ? 'Experimental group' : 'Control group';
    return `<section class="control-section" data-group="${g}" data-source="shape"><h2><span class="group-dot ${g ? 'exp' : ''}"></span>${name}</h2>
      ${paired ? '' : range(`M${g}_min`, `M${g}_max`, `Samples M${g}`)}
      ${field(`N${g}${s}`, 'Cells per sample N' + g, 1000, 'Total cells profiled per biological sample, not just the cell type of interest.', 'min="1"')}
      <h3>Expected cell abundance</h3>
      ${field(`mean${g}${s}`, 'Mean proportion', '', 'Mean cell type proportion across biological samples. 0.03 means 3%.')}
      ${field(`std${g}${s}`, 'Standard deviation', '', 'Standard deviation of cell type proportion across biological samples. Enter a proportion, not a percentage.')}
      <details class="beta-params"><summary>Or set beta distribution parameters</summary><div class="shape-fields">
        ${field(`a${g}${s}`, 'a', g ? 10 : 7, '', 'aria-label="' + name + ', beta a"')}
        ${field(`b${g}${s}`, 'b', g ? 7 : 10, '', 'aria-label="' + name + ', beta b"')}
      </div><p class="hint">Mean/SD and a/b update each other. Calculations use a/b.</p></details>
    </section>`;
  }
  function setup(paired) {
    const s = suffix(paired), container = $(paired ? 'paired' : 'main');
    const previousPlot = $('plot' + s);
    if (previousPlot) {
      chartObserver.unobserve(previousPlot);
      if (typeof Plotly !== 'undefined' && previousPlot.data) Plotly.purge(previousPlot);
    }
    const radios = [['one_sided', 'One-sided t-test'], ['two_sided', 'Two-sided t-test'], ['sesoi_test', 't-test with SESOI'], ['tost_test', 'TOST for equivalence']]
      .map(([id, text], i) => `<label for="${id}${s}"><input type="radio" name="test_type${s}" id="${id}${s}" value="${id}" ${i === 0 ? 'checked' : ''}>${text}</label>`).join('');
    container.innerHTML = `<aside class="controls" aria-label="Study parameters">
      <section class="control-section"><div class="section-heading"><h2>${paired ? 'Paired study' : 'Unpaired study'}</h2><button class="reset" type="button" id="reset${s}">Reset defaults</button></div>
      <p class="hint">${paired ? 'Matched samples from the same participants. M is the number of pairs.' : 'Independent control and experimental groups. Sample sizes may differ.'}</p>
      ${paired ? range('M_min_paired', 'M_max_paired', 'Pairs M') + field('rho_paired', 'Correlation ρ', 0.5, 'Correlation of true cell proportions between matched samples, from −1 to 1.', 'min="-1" max="1"') : ''}
      <p class="hint group-note">Enter proportions as decimals: <strong>0.03 = 3%</strong>.</p></section>
      ${group(0, s, paired)}${group(1, s, paired)}
      <section class="control-section"><h2>Statistical test</h2>
        ${field('alpha' + s, 'Significance level α', 0.05, 'False positive rate, commonly 0.05 or 0.01.')}
        ${field('beta' + s, 'False negative threshold βth', 0.1, 'Highlight results strictly below this false negative rate. Power = 1 − β.')}
        <div class="radio-list" role="group" aria-label="Statistical test type">${radios}</div>
        <div id="sesoi_param${s}" class="test-params" hidden><div class="radio-list" role="group" aria-label="SESOI direction">
          <label for="sesoi_verify${s}"><input type="radio" name="sesoi_type${s}" id="sesoi_verify${s}" checked>Verify</label>
          <label for="sesoi_falsify${s}"><input type="radio" name="sesoi_type${s}" id="sesoi_falsify${s}">Falsify</label></div>
          ${field('sesoi_bound' + s, 'SESOI offset (≠ 0)', 0.01)}<p class="hint">0.01 = one percentage point. Negative offsets reverse the test direction.</p></div>
        <div id="tost_param${s}" class="test-params" hidden>
          ${field('tost_lower_bound' + s, 'Lower offset', -0.1)}${field('tost_upper_bound' + s, 'Upper offset', 0.1)}<p class="hint">Bounds are offsets from control, in proportions.</p></div>
      </section></aside>
      <div class="results-column">
        <div id="error${s}" class="status-box" role="status" hidden></div>
        <div id="warning${s}" class="status-box warning" role="status" hidden></div>
        <section class="panel chart-panel"><div class="panel-head"><div><h2>Expected cell abundance</h2><p class="hint">Beta distributions across biological samples</p></div></div><div id="plot${s}" class="chart" role="img" aria-label="Control and experimental beta distributions"></div><p class="hint chart-note">Proportion of the cell type of interest; density indicates how likely each proportion is.</p></section>
        <section class="panel result-panel"><div class="panel-head"><div><h2>False negative rate <span class="muted">β</span></h2><p class="hint">Lower is better. Statistical power = 1 − β.</p></div><button type="button" class="button" id="export${s}">Download CSV</button></div>
          <div class="result-meta"><span id="legend${s}"></span><span id="count${s}"></span></div>
          <p class="table-axes">${paired ? 'Columns: number of matched pairs M' : 'Rows: control samples M0 · Columns: experimental samples M1'}</p>
          <div class="table-scroll" tabindex="0" role="region" aria-label="False negative results table"><div id="result${s}"></div></div>
          <div id="summary${s}" class="result-summary" role="status"></div>
          <p class="footnote">Values rounded to 3 decimal places. Highlighting uses unrounded β &lt; βth.</p>
        </section>
        <details class="panel method-note"><summary>How to read these results</summary><p>Each value is the estimated probability of failing to reach statistical significance under the chosen test and expected cell abundances. A value of 0.1 means 10% false negatives, or 90% statistical power.</p><p>Sample sizes count biological samples${paired ? ' (matched pairs)' : ' in each group'}, not cells. More cells per sample and more biological samples address different sources of variability. <a href="#instructions">Read the guide &amp; methods</a>.</p></details>
      </div>`;
    states.set(paired, { paired, s, container, rows: [], params: null });
    chartObserver.observe($('plot' + s));
    for (let g = 0; g < 2; g++) syncGroup(paired, g, 'shape');
    container.oninput = event => {
      const input = event.target;
      if (!input.matches('input')) return;
      const match = /^(mean|std|a|b)([01])(?:_paired)?$/.exec(input.id);
      if (match) syncGroup(paired, Number(match[2]), ['a', 'b'].includes(match[1]) ? 'shape' : 'moments');
      update(paired);
    };
    $('reset' + s).addEventListener('click', () => { setup(paired); update(paired); });
    $('export' + s).addEventListener('click', () => download(paired));
  }
  function syncGroup(paired, g, source) {
    const s = suffix(paired), section = $(paired ? 'paired' : 'main').querySelector(`[data-group="${g}"]`);
    section.dataset.source = source;
    const get = key => Number($(key + g + s).value);
    const filled = keys => keys.every(key => $(key + g + s).value !== '');
    if (source === 'moments') {
      const m = get('mean'), sd = get('std');
      if (!filled(['mean', 'std']) || !(m > 0 && m < 1 && sd > 0 && sd * sd < m * (1 - m))) return;
      const [a, b] = Sensei.betaStatParam(m, sd);
      $('a' + g + s).value = Sensei.round(a, 3); $('b' + g + s).value = Sensei.round(b, 3);
    } else {
      if (!filled(['a', 'b']) || !(get('a') > 0 && get('b') > 0)) return;
      const [m, sd] = Sensei.betaParamStat(get('a'), get('b'));
      $('mean' + g + s).value = Sensei.round(m, 4); $('std' + g + s).value = Sensei.round(sd, 4);
    }
  }
  function read(paired) {
    const state = states.get(paired), { s, container } = state;
    const errors = [];
    container.querySelectorAll('[aria-invalid]').forEach(input => input.removeAttribute('aria-invalid'));
    function invalid(ids, message) { errors.push(message); ids.forEach(id => $(id).setAttribute('aria-invalid', 'true')); }
    function number(id) {
      const input = $(id), value = input.value === '' ? NaN : Number(input.value);
      if (!Number.isFinite(value)) invalid([id], 'Complete all active fields with finite numbers.');
      return value;
    }
    const p = { paired, test: container.querySelector(`input[name="test_type${s}"]:checked`).value };
    for (const key of ['N0', 'N1', 'a0', 'a1', 'b0', 'b1', 'alpha', 'beta']) p[key] = number(key + s);
    for (const key of ['N0', 'N1']) if (!(Number.isSafeInteger(p[key]) && p[key] > 0)) invalid([key + s], 'Cells per sample must be positive whole numbers.');
    for (const key of ['a0', 'a1', 'b0', 'b1']) if (!(p[key] > 0)) invalid([key + s], 'Beta parameters a and b must be greater than zero. Very small values may round to zero.');
    for (const key of ['alpha', 'beta']) if (!(p[key] > 0 && p[key] < 1)) invalid([key + s], 'Significance level and false negative threshold must be between 0 and 1.');
    for (let g = 0; g < 2; g++) {
      if (container.querySelector(`[data-group="${g}"]`).dataset.source !== 'moments') continue;
      const m = number('mean' + g + s), sd = number('std' + g + s);
      if (!(m > 0 && m < 1 && sd > 0 && sd * sd < m * (1 - m))) invalid(['mean' + g + s, 'std' + g + s], 'Mean must be between 0 and 1; SD must be positive with SD² < mean × (1 − mean).');
    }
    const ranges = paired ? [['M_min_paired', 'M_max_paired']] : [['M0_min', 'M0_max'], ['M1_min', 'M1_max']];
    p.ranges = ranges.map(ids => {
      const [lo, hi] = ids.map(number);
      if (!(Number.isSafeInteger(lo) && Number.isSafeInteger(hi) && lo >= 2 && hi >= lo)) invalid(ids, 'Sample ranges must use whole numbers ≥ 2, with maximum ≥ minimum.');
      return [lo, hi];
    });
    const size = p.ranges.reduce((total, [lo, hi]) => total * (hi - lo + 1), 1);
    if (size > 10000) invalid(ranges.flat(), 'Choose a smaller range: at most 10,000 result cells can be displayed at once.');
    if (paired) { p.rho = number('rho_paired'); if (!(p.rho >= -1 && p.rho <= 1)) invalid(['rho_paired'], 'Paired correlation must be between −1 and 1.'); }
    if (p.test === 'sesoi_test') {
      p.offset = number('sesoi_bound' + s); p.verify = $('sesoi_verify' + s).checked;
      if (p.offset === 0) invalid(['sesoi_bound' + s], 'Use a nonzero SESOI offset, or choose an ordinary t-test.');
    }
    if (p.test === 'tost_test') {
      p.lower = number('tost_lower_bound' + s); p.upper = number('tost_upper_bound' + s);
      if (!(p.lower <= 0 && p.upper >= 0 && p.lower < p.upper)) invalid(['tost_lower_bound' + s, 'tost_upper_bound' + s], 'TOST requires lower ≤ 0 ≤ upper and lower < upper.');
    }
    return { p, errors: [...new Set(errors)] };
  }
  function update(paired) {
    hideTooltip();
    const state = states.get(paired), { s, container } = state;
    const test = container.querySelector(`input[name="test_type${s}"]:checked`).value;
    $('sesoi_param' + s).hidden = test !== 'sesoi_test'; $('tost_param' + s).hidden = test !== 'tost_test';
    const { p, errors } = read(paired);
    state.rows = []; state.params = null;
    $('error' + s).hidden = !errors.length; $('error' + s).textContent = errors.join(' ');
    $('export' + s).disabled = !!errors.length; $('warning' + s).hidden = true;
    if (errors.length) {
      $('result' + s).innerHTML = '<p class="empty">Update the highlighted inputs to calculate results.</p>';
      $('count' + s).textContent = ''; $('summary' + s).textContent = ''; $('legend' + s).textContent = '';
      if (typeof Plotly !== 'undefined') Plotly.purge('plot' + s);
      return;
    }
    const sequence = ([lo, hi]) => Array.from({ length: hi - lo + 1 }, (_, i) => lo + i);
    const rowSizes = paired ? [null] : sequence(p.ranges[0]);
    const colSizes = sequence(p.ranges[paired ? 0 : 1]);
    const rows = rowSizes.map(m0 => colSizes.map(m1 => ({ m0: paired ? m1 : m0, m1, beta: Sensei.falseNegative(p, paired ? m1 : m0, m1) })));
    if (rows.flat().some(cell => !Number.isFinite(cell.beta) || cell.beta < 0 || cell.beta > 1)) {
      $('error' + s).hidden = false; $('error' + s).textContent = 'These parameters exceed the numerical limits of the model. Adjust beta parameters or sample sizes.';
      $('result' + s).innerHTML = '<p class="empty">No finite result is available for these inputs.</p>';
      $('count' + s).textContent = ''; $('summary' + s).textContent = ''; $('legend' + s).textContent = '';
      $('export' + s).disabled = true;
      if (typeof Plotly !== 'undefined') Plotly.purge('plot' + s);
      return;
    }
    state.rows = rows; state.params = p;
    const passing = rows.flat().filter(c => c.beta < p.beta);
    $('legend' + s).innerHTML = `<span class="legend-chip"></span>Highlighted: β &lt; ${p.beta} (power &gt; ${Sensei.round((1 - p.beta) * 100, 6)}%)`;
    $('count' + s).textContent = `${passing.length} of ${rows.flat().length} designs meet threshold`;
    $('result' + s).innerHTML = `<table class="result-table"><caption>False negative rates for ${paired ? 'paired' : 'unpaired'} study</caption><thead><tr><th scope="col">${paired ? 'Pairs M' : 'Control / Exp.'}</th>${colSizes.map(m => `<th scope="col">${m}</th>`).join('')}</tr></thead><tbody>${rows.map((row, i) => `<tr><th scope="row">${paired ? 'False negative β' : rowSizes[i]}</th>${row.map((c, j) => `<td class="${c.beta < p.beta ? 'meets' : ''}" tabindex="${i === 0 && j === 0 ? 0 : -1}" data-m0="${c.m0}" data-m1="${c.m1}" data-beta="${c.beta}">${Sensei.round(c.beta, 3)}</td>`).join('')}</tr>`).join('')}</tbody></table>`;
    const smallest = passing.slice().sort((a, b) => (a.m0 + a.m1) - (b.m0 + b.m1) || a.m0 - b.m0)[0];
    $('summary' + s).textContent = smallest ? `Smallest ${paired ? 'pair count' : 'total sample count'} meeting threshold in this range: ${paired ? smallest.m0 + ' matched pairs' : smallest.m0 + ' control + ' + smallest.m1 + ' experimental samples'} · power ${((1 - smallest.beta) * 100).toFixed(1)}%.` : 'No design in this range meets the threshold. Try a larger sample range or more cells per sample.';
    if ([p.a0, p.a1, p.b0, p.b1].some(v => v <= 1)) {
      $('warning' + s).hidden = false; $('warning' + s).textContent = 'At least one beta parameter is ≤ 1. The paper recommends a and b > 1 for unimodal distributions; inspect the plot and model assumptions.';
    }
    plot(paired);
  }
  // Display only: cap decimal places at six before using scientific notation.
  function formatBeta(value) {
    if (value > 0 && value < 0.001) {
      const decimals = -Math.floor(Math.log10(value));
      return decimals > 6 ? value.toExponential(0) : value.toFixed(decimals).replace(/0+$/, '');
    }
    return value.toFixed(3);
  }
  function hideTooltip() {
    clearTimeout(tooltipHideTimer);
    if (tooltipCell) tooltipCell.removeAttribute('aria-describedby');
    tooltipCell = null;
    tooltip.hidden = true;
  }
  function showTooltip(cell) {
    hideTooltip();
    const paired = cell.closest('.workspace').id === 'paired';
    const p = states.get(paired).params;
    if (!p) return;
    const beta = Number(cell.dataset.beta), meets = beta < p.beta;
    const items = paired ? [['Matched pairs', cell.dataset.m0]] : [['Control samples', cell.dataset.m0], ['Experimental samples', cell.dataset.m1]];
    items.push(['False negative rate β', formatBeta(beta)], ['Statistical power', `${((1 - beta) * 100).toFixed(2)}%`]);
    tooltip.innerHTML = `<p class="tooltip-heading">${paired ? 'Paired' : 'Unpaired'} study design</p><dl>${items.map(([label, value]) => `<div><dt>${label}</dt><dd>${value}</dd></div>`).join('')}</dl><p class="tooltip-status${meets ? ' meets' : ''}">${meets ? 'Meets threshold' : 'Does not meet threshold'}<span>Target: β &lt; ${p.beta}</span></p>`;
    tooltipCell = cell;
    cell.setAttribute('aria-describedby', tooltip.id);
    tooltip.hidden = false;
    const anchor = cell.getBoundingClientRect(), box = tooltip.getBoundingClientRect(), gap = 10;
    const left = Math.max(gap, Math.min(anchor.left + (anchor.width - box.width) / 2, innerWidth - box.width - gap));
    const preferredTop = anchor.top >= box.height + 2 * gap ? anchor.top - box.height - gap : anchor.bottom + gap;
    const top = Math.max(gap, Math.min(preferredTop, innerHeight - box.height - gap));
    tooltip.style.left = `${Math.round(left)}px`;
    tooltip.style.top = `${Math.round(top)}px`;
  }
  // One shared popup stays outside the scrolling table and fits the viewport.
  const resultCell = target => target.closest('.result-table td[data-beta]');
  document.addEventListener('pointerover', event => {
    const cell = resultCell(event.target);
    if (cell) showTooltip(cell);
    else if (tooltip.contains(event.target)) clearTimeout(tooltipHideTimer);
  });
  document.addEventListener('pointermove', event => {
    const cell = resultCell(event.target);
    if (cell && cell !== tooltipCell) showTooltip(cell);
  });
  document.addEventListener('pointerout', event => {
    if (!tooltipCell || (event.relatedTarget && (tooltipCell.contains(event.relatedTarget) || tooltip.contains(event.relatedTarget)))) return;
    if (document.activeElement !== tooltipCell) tooltipHideTimer = setTimeout(hideTooltip, 120);
  });
  document.addEventListener('focusin', event => {
    const cell = resultCell(event.target);
    if (!cell) { hideTooltip(); return; }
    cell.closest('table').querySelectorAll('td[tabindex="0"]').forEach(previous => { previous.tabIndex = -1; });
    cell.tabIndex = 0;
    showTooltip(cell);
  });
  document.addEventListener('click', event => {
    const cell = resultCell(event.target);
    if (cell) showTooltip(cell);
    else if (!tooltip.contains(event.target)) hideTooltip();
  });
  document.addEventListener('keydown', event => {
    if (event.key === 'Escape') { hideTooltip(); return; }
    const cell = resultCell(event.target);
    if (!cell) return;
    const row = cell.parentElement, table = cell.closest('table');
    const next = event.key === 'ArrowLeft' ? row.cells[cell.cellIndex - 1]
      : event.key === 'ArrowRight' ? row.cells[cell.cellIndex + 1]
      : event.key === 'ArrowUp' ? table.rows[row.rowIndex - 1]?.cells[cell.cellIndex]
      : event.key === 'ArrowDown' ? table.rows[row.rowIndex + 1]?.cells[cell.cellIndex]
      : null;
    if (next?.matches('td[data-beta]')) { event.preventDefault(); next.focus({ preventScroll: true }); next.scrollIntoView({ block: 'nearest', inline: 'nearest' }); }
  });
  window.addEventListener('scroll', () => {
    if (!tooltipCell) return;
    const anchor = tooltipCell.getBoundingClientRect();
    const table = tooltipCell.closest('.table-scroll').getBoundingClientRect();
    const visible = anchor.bottom > Math.max(0, table.top) && anchor.top < Math.min(innerHeight, table.bottom)
      && anchor.right > Math.max(0, table.left) && anchor.left < Math.min(innerWidth, table.right);
    if (visible) showTooltip(tooltipCell); else hideTooltip();
  }, true);
  window.addEventListener('resize', hideTooltip);
  function plot(paired) {
    const state = states.get(paired), { s, params: p, container } = state;
    if (!p || container.hidden) return;
    if (typeof Plotly === 'undefined') { $('plot' + s).textContent = 'Distribution chart could not load. The result table is available.'; return; }
    const x = []; for (let v = 0; v <= 1; v += 0.005) x.push(v);
    const density = (a, b) => x.map(v => { const y = jStat.beta.pdf(v, a, b); return Number.isFinite(y) ? y : null; });
    const y0 = density(p.a0, p.b0);
    const traces = [{ x, y: y0, name: 'Control', type: 'scatter', line: { color: '#2c5fa8', width: 2 }, fill: 'tozeroy', fillcolor: 'rgba(44,95,168,.055)' }, { x, y: density(p.a1, p.b1), name: 'Experimental', type: 'scatter', line: { color: '#c46b32', width: 2 }, fill: 'tozeroy', fillcolor: 'rgba(196,107,50,.045)' }];
    const shifted = (offset, name, dash) => traces.push({ x: x.map(v => v + offset), y: y0, name, type: 'scatter', line: { color: '#2c5fa8', dash, width: 1.5 } });
    if (p.test === 'sesoi_test') shifted(p.offset, 'Control + SESOI', 'dot');
    if (p.test === 'tost_test') { shifted(p.lower, 'Lower bound', 'dot'); shifted(p.upper, 'Upper bound', 'dash'); }
    const chart = $('plot' + s);
    Plotly.react(chart, traces, { autosize: true, width: chart.clientWidth, height: chart.clientHeight, margin: { l: 45, r: 15, b: 42, t: 32 }, font: { family: 'system-ui, sans-serif', size: 11, color: '#5d6677' }, paper_bgcolor: '#fff', plot_bgcolor: '#fff', legend: { orientation: 'h', x: 0, y: 1.22 }, xaxis: { title: { text: 'Cell type proportion', font: { size: 11 } }, gridcolor: '#edf0f4', zerolinecolor: '#dde1e8', tickformat: '.0%' }, yaxis: { title: { text: 'Density', font: { size: 11 } }, gridcolor: '#edf0f4', rangemode: 'tozero', zerolinecolor: '#dde1e8' }, hovermode: 'closest' }, { displayModeBar: false });
  }
  function download(paired) {
    const { rows, params: p } = states.get(paired); if (!p) return;
    const lines = [['study', 'test', 'control_samples', 'experimental_samples', 'cells_control', 'cells_experimental', 'a_control', 'b_control', 'a_experimental', 'b_experimental', 'alpha', 'beta_threshold', 'rho', 'sesoi_offset', 'sesoi_direction', 'tost_lower', 'tost_upper', 'false_negative_rate', 'power', 'meets_threshold']];
    for (const c of rows.flat()) lines.push([paired ? 'paired' : 'unpaired', p.test, c.m0, c.m1, p.N0, p.N1, p.a0, p.b0, p.a1, p.b1, p.alpha, p.beta, paired ? p.rho : '', p.offset ?? '', p.test === 'sesoi_test' ? (p.verify ? 'verify' : 'falsify') : '', p.lower ?? '', p.upper ?? '', c.beta, 1 - c.beta, c.beta < p.beta]);
    const url = URL.createObjectURL(new Blob([lines.map(row => row.join(',')).join('\r\n') + '\r\n'], { type: 'text/csv;charset=utf-8' }));
    const a = document.createElement('a'); a.href = url; a.download = `sensei-${paired ? 'paired' : 'unpaired'}.csv`; a.click(); setTimeout(() => URL.revokeObjectURL(url), 1000);
  }
  function navigate() {
    hideTooltip();
    const page = location.hash === '#paired' ? 'paired' : location.hash === '#instructions' ? 'instructions' : 'main';
    for (const id of ['main', 'paired', 'instructions']) $(id).hidden = id !== page;
    document.querySelectorAll('[data-page]').forEach(link => { if (link.dataset.page === page) link.setAttribute('aria-current', 'page'); else link.removeAttribute('aria-current'); });
    if (page !== 'instructions') plot(page === 'paired');
  }
  setup(false); setup(true); update(false); update(true);
  window.addEventListener('hashchange', navigate); navigate();
})();
