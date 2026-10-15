/* Real browser integration tests against both independently implemented pages.
 * All non-loopback requests are blocked to prove the pages are self-contained.
 */
const { chromium } = require('playwright');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const { pathToFileURL } = require('node:url');
const server = require('../scripts/serve.cjs');
const output = path.resolve(__dirname, '../.artifacts');
fs.mkdirSync(output, { recursive: true });
const cases = [];
for (const paired of [false, true]) {
  for (const test of ['one_sided', 'two_sided', 'sesoi_test', 'tost_test']) {
    cases.push({ paired, test, name: `${paired ? 'paired' : 'unpaired'} ${test}` });
    if (test === 'sesoi_test') {
      cases.push({ paired, test, verify: false, name: `${paired ? 'paired' : 'unpaired'} SESOI falsify` });
      cases.push({ paired, test, offset: -0.05, name: `${paired ? 'paired' : 'unpaired'} SESOI negative verify` });
      cases.push({ paired, test, offset: -0.05, verify: false, name: `${paired ? 'paired' : 'unpaired'} SESOI negative falsify` });
    }
  }
  cases.push({ paired, test: 'one_sided', fields: { N0: 100, N1: 384, mean0: 0.03, std0: 0.015, mean1: 0.05, std1: 0.01, alpha: 0.01, beta: 0.2 }, name: `${paired ? 'paired' : 'unpaired'} rare-cell moments, unequal cells` });
  cases.push({ paired, test: 'two_sided', fields: { a0: 7, b0: 10, a1: 7, b1: 10 }, name: `${paired ? 'paired' : 'unpaired'} equal means` });
  cases.push({ paired, test: 'tost_test', fields: { tost_lower_bound: -0.0001, tost_upper_bound: 0.0001 }, name: `${paired ? 'paired' : 'unpaired'} TOST capped at one` });
}
for (const rho of [-1, 0, 1]) cases.push({ paired: true, test: 'one_sided', fields: { rho }, name: `paired correlation ${rho}` });
cases.push({ paired: false, test: 'two_sided', ranges: { M0_min: 2, M0_max: 7, M1_min: 10, M1_max: 14 }, name: 'unpaired unequal sample ranges' });
cases.push({ paired: true, test: 'one_sided', ranges: { M_min_paired: 2, M_max_paired: 25 }, name: 'paired range starts at two' });
let checks = 0, cells = 0;
async function apply(page, scenario, modern) {
  if (scenario.paired) {
    if (modern) await page.getByRole('link', { name: 'Paired t-test', exact: true }).click();
    else await page.getByRole('button', { name: 'Paired t-Test', exact: true }).click();
    await page.locator('#paired').waitFor({ state: 'visible' });
  }
  // Set all related values before firing the documented edit event; prevents
  // transient invalid moment pairs from changing the final conversion order.
  await page.evaluate(({ scenario, modern }) => {
    const s = scenario.paired ? '_paired' : '';
    const changed = [];
    for (const [key, value] of Object.entries(scenario.fields || {})) {
      const id = key === 'rho' ? 'rho_paired' : key + s;
      document.getElementById(id).value = value; changed.push(id);
    }
    for (const [id, value] of Object.entries(scenario.ranges || {})) { document.getElementById(id).value = value; changed.push(id); }
    if (scenario.offset !== undefined) { document.getElementById('sesoi_bound' + s).value = scenario.offset; changed.push('sesoi_bound' + s); }
    document.getElementById(scenario.test + s).checked = true;
    if (scenario.verify === false) document.getElementById('sesoi_falsify' + s).checked = true;
    changed.push(scenario.test + s);
    for (const id of changed) document.getElementById(id).dispatchEvent(new Event(modern ? 'input' : 'change', { bubbles: true }));
  }, { scenario, modern });
}
async function results(page, paired, modern) {
  return page.evaluate(({ paired, modern }) => {
    const s = paired ? '_paired' : '';
    const rows = Array.from(document.querySelectorAll('#result' + s + ' tr'));
    return rows.slice(1).map(row => Array.from(row.children).slice(1).map(cell => ({ value: Number(cell.textContent), pass: modern ? cell.classList.contains('meets') : cell.hasAttribute('bgcolor') })));
  }, { paired, modern });
}
(async () => {
  await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
  const base = `http://127.0.0.1:${server.address().port}`;
  const browser = await chromium.launch({ channel: process.env.PLAYWRIGHT_CHANNEL || (process.platform === 'win32' ? 'msedge' : undefined), headless: true });
  const context = await browser.newContext({ viewport: { width: 1440, height: 1100 }, acceptDownloads: true });
  context.setDefaultTimeout(10000);
  const errors = [], blocked = [], environmentRequests = [];
  await context.route('**/*', route => {
    if (route.request().url().startsWith(base + '/')) return route.continue();
    // The host's system AdGuard can inject these scripts into isolated Edge.
    // Still block them; distinguish them from requests initiated by this app.
    const url = route.request().url();
    if (new URL(url).hostname === 'local.adguard.org') environmentRequests.push(url);
    else blocked.push(url);
    return route.abort();
  });
  const modern = await context.newPage(), old = await context.newPage();
  for (const page of [modern, old]) page.on('pageerror', error => errors.push(error.message));
  try {
    for (const scenario of cases) {
      await Promise.all([modern.goto(base + '/table_beta.html'), old.goto(base + '/table_beta_legacy.html')]);
      await Promise.all([apply(modern, scenario, true), apply(old, scenario, false)]);
      const [actual, expected] = await Promise.all([results(modern, scenario.paired, true), results(old, scenario.paired, false)]);
      assert.ok(actual.length > 0 && actual[0].length > 0, scenario.name);
      assert.deepEqual(actual, expected, scenario.name);
      cells += actual.flat().length; checks++;
      const s = scenario.paired ? '_paired' : '';
      for (const page of [modern, old]) {
        const available = await page.evaluate(id => !!document.getElementById(id).data, 'plot' + s);
        assert.ok(available, `${scenario.name}: missing plot at ${page.url()}; browser errors: ${errors.join('; ')}`);
      }
      const curves = await Promise.all([modern, old].map(page => page.evaluate(id => document.getElementById(id).data.map(trace => ({ x: trace.x, y: trace.y })), 'plot' + s)));
      assert.deepEqual(curves[0], curves[1], scenario.name + ' distribution curves');
      const params = await Promise.all([modern, old].map(page => page.evaluate(s => Object.fromEntries(['a0', 'a1', 'b0', 'b1', 'mean0', 'mean1', 'std0', 'std1'].map(key => [key, document.getElementById(key + s).value])), s)));
      assert.deepEqual(params[0], params[1], scenario.name + ' conversion fields');
    }
    console.log(`Passed ${cases.length} comparison scenarios (${cells} table cells). Checking interactions.`);
    await modern.goto(base + '/table_beta.html');
    for (const paired of [false, true]) {
      const s = paired ? '_paired' : '';
      await modern.getByRole('link', { name: paired ? 'Paired t-test' : 'Unpaired t-test', exact: true }).click();
      await modern.locator(paired ? '#paired' : '#main').waitFor({ state: 'visible' });
      const steps = await modern.locator((paired ? '#paired' : '#main') + ' input[type=number]').evaluateAll(inputs => inputs.map(input => ({ id: input.id, step: input.step })));
      for (const { id, step } of steps) assert.equal(step, /^[MN]/.test(id) ? '1' : '0.01', id + ' spinner step');
      for (const [id, initial, expected] of [[`alpha${s}`, '0.05', '0.06'], [`beta${s}`, '0.1', '0.11'], [`mean0${s}`, '0.41', '0.42'], [`std0${s}`, '0.11', '0.12'], [`N0${s}`, '1000', '1001'], ...(paired ? [['rho_paired', '0.5', '0.51']] : [['M0_min', '5', '6']])]) {
        await modern.locator('#' + id).fill(initial);
        await modern.locator('#' + id).press('ArrowUp');
        assert.equal(await modern.locator('#' + id).inputValue(), expected, id + ' ArrowUp');
      }
      await modern.locator('#reset' + s).click();
      // Click the actual upper/lower native spinner buttons as well as using keys.
      const alphaInput = modern.locator('#alpha' + s);
      await alphaInput.hover();
      const inputBox = await alphaInput.boundingBox();
      await alphaInput.click({ position: { x: inputBox.width - 16, y: inputBox.height / 2 - 4 } });
      assert.equal(await alphaInput.inputValue(), '0.06');
      await alphaInput.click({ position: { x: inputBox.width - 16, y: inputBox.height / 2 + 4 } });
      assert.equal(await alphaInput.inputValue(), '0.05');
      await modern.locator('#reset' + s).click(); checks++;
      const cell = modern.locator('#result' + s + ' td[data-beta]').first(), tooltip = modern.getByRole('tooltip');
      assert.equal(await cell.getAttribute('title'), null, 'no browser tooltip competes with the card');
      await cell.hover(); await tooltip.waitFor({ state: 'visible' });
      assert.deepEqual(await tooltip.locator('dt').allTextContents(), paired ? ['Matched pairs', 'False negative rate β', 'Statistical power'] : ['Control samples', 'Experimental samples', 'False negative rate β', 'Statistical power']);
      const rate = Number(await cell.getAttribute('data-beta'));
      assert.ok((await tooltip.locator('dd').allTextContents()).includes(rate.toFixed(3)));
      assert.ok((await tooltip.locator('dd').allTextContents()).includes(((1 - rate) * 100).toFixed(2) + '%'));
      assert.ok((await tooltip.locator('.tooltip-status').textContent()).includes('Does not meet threshold'));
      assert.equal(await cell.getAttribute('aria-describedby'), 'result-tooltip');
      if (!paired) await tooltip.screenshot({ path: path.join(output, 'sensei-table-tooltip.png') });
      await modern.locator('h1').hover(); await tooltip.waitFor({ state: 'hidden' });
      await cell.focus(); await tooltip.waitFor({ state: 'visible' });
      await cell.press('ArrowRight');
      assert.equal(await modern.locator('#result' + s + ' td[tabindex="0"]').count(), 1);
      assert.equal(await modern.locator('#result' + s + ' td').nth(1).evaluate(el => el === document.activeElement), true);
      await modern.keyboard.press('Escape'); await tooltip.waitFor({ state: 'hidden' });
      const passingCell = modern.locator('#result' + s + ' td.meets').last();
      await passingCell.hover(); await tooltip.waitFor({ state: 'visible' });
      assert.ok((await tooltip.locator('.tooltip-status.meets').textContent()).startsWith('Meets threshold'));
      await alphaInput.fill('0.06'); await tooltip.waitFor({ state: 'hidden' });
      await modern.locator('#reset' + s).click(); checks++;
    }
    await modern.getByRole('link', { name: 'Unpaired t-test', exact: true }).click();
    await modern.locator('#main').waitFor({ state: 'visible' });
    await modern.locator('h1').hover();
    await modern.evaluate(() => window.scrollTo(0, 0));
    console.log('Passed decimal/integer spinner buttons and keyboard checks, plus structured tooltips in both study modes.');
    await modern.screenshot({ path: path.join(output, 'sensei-desktop.png'), fullPage: true });
    // Exercise actual user edits and recovery, separately from batch scenarios.
    await modern.locator('#N0').fill('384');
    assert.notEqual((await results(modern, false, true))[0][0].value, 0.307);
    await modern.getByRole('button', { name: 'Reset defaults', exact: true }).first().click();
    assert.equal((await results(modern, false, true))[0][0].value, 0.307); checks++;
    const invalid = [['N0', '0'], ['M0_min', '1'], ['M0_min', '5.5'], ['M0_min', '13'], ['alpha', '1'], ['beta', '-0.1'], ['std0', '1'], ['mean0', '0'], ['N0', ''], ['M0_max', '100000']];
    for (const [id, value] of invalid) {
      await modern.locator('#' + id).fill(value);
      assert.equal(await modern.locator('#error').isVisible(), true, id + '=' + value);
      assert.equal(await modern.locator('#export').isDisabled(), true);
      assert.equal(await modern.locator('#result td').count(), 0);
      await modern.locator('#reset').click();
      assert.equal((await results(modern, false, true))[0][0].value, 0.307); checks++;
    }
    await modern.locator('#sesoi_test').check(); await modern.locator('#sesoi_bound').fill('0');
    assert.equal(await modern.locator('#error').isVisible(), true); checks++;
    await modern.locator('#reset').click(); await modern.locator('#tost_test').check(); await modern.locator('#tost_upper_bound').fill('-0.2');
    assert.equal(await modern.locator('#error').isVisible(), true); checks++;
    await modern.locator('#reset').click();
    console.log('Passed invalid input recovery. Checking CSV downloads.');
    for (const paired of [false, true]) {
      if (paired) await modern.getByRole('link', { name: 'Paired t-test', exact: true }).click();
      const s = paired ? '_paired' : '';
      const waiting = modern.waitForEvent('download'); await modern.locator('#export' + s).click(); const download = await waiting;
      const csv = await fs.promises.readFile(await download.path(), 'utf8');
      const lines = csv.trim().split(/\r?\n/), headers = lines[0].split(',');
      assert.equal(lines.length, paired ? 9 : 65);
      const betaIndex = headers.indexOf('false_negative_rate');
      const firstValue = await modern.locator('#result' + s + ' td[data-beta]').first().getAttribute('data-beta');
      assert.equal(Number(lines[1].split(',')[betaIndex]), Number(firstValue)); checks++;
    }
    await modern.locator('#rho_paired').fill('1.1'); assert.equal(await modern.locator('#error_paired').isVisible(), true); checks++;
    await modern.locator('#reset_paired').click();
    await modern.getByRole('link', { name: 'Unpaired t-test', exact: true }).click();
    // Strict threshold comparison is performed before display rounding.
    const exact = await modern.locator('#result td[data-beta]').first().getAttribute('data-beta');
    await modern.locator('#beta').fill(exact); assert.equal(await modern.locator('#result td').first().getAttribute('class'), ''); checks++;
    await modern.locator('#reset').click();
    // Narrow mobile, tablet, desktop: overflow belongs to the table, not body.
    for (const width of [320, 390, 768, 1440]) {
      await modern.setViewportSize({ width, height: 1000 });
      await modern.waitForFunction(() => Math.abs(document.querySelector('#plot .svg-container').getBoundingClientRect().width - document.getElementById('plot').clientWidth) < 2);
      const dimensions = await modern.evaluate(() => ({ scroll: document.documentElement.scrollWidth, client: document.documentElement.clientWidth }));
      assert.ok(dimensions.scroll <= dimensions.client + 1, 'page overflow at ' + width + ': ' + JSON.stringify(dimensions)); checks++;
      await modern.screenshot({ path: path.join(output, `sensei-${width}.png`), fullPage: true });
      await modern.locator('#result td').last().hover();
      await modern.getByRole('tooltip').waitFor({ state: 'visible' }).catch(error => { throw new Error(`Tooltip failed at viewport width ${width}: ${error.message}`); });
      const popupBounds = await modern.getByRole('tooltip').boundingBox();
      assert.ok(popupBounds.x >= 0 && popupBounds.x + popupBounds.width <= width, 'tooltip fits width ' + width);
      assert.ok(popupBounds.y >= 0 && popupBounds.y + popupBounds.height <= 1000, 'tooltip fits viewport height');
      if (width === 390) await modern.getByRole('tooltip').screenshot({ path: path.join(output, 'sensei-tooltip-mobile.png') });
      await modern.keyboard.press('Escape');
      checks++;
    }
    await modern.goto(base + '/table_beta.html#paired'); await modern.locator('#paired').waitFor({ state: 'visible' }); checks++;
    await modern.getByRole('link', { name: 'Guide & methods', exact: true }).click(); await modern.locator('#instructions').waitFor({ state: 'visible' }); checks++;
    await modern.goBack(); await modern.locator('#paired').waitFor({ state: 'visible' }); checks++;
    await modern.getByRole('link', { name: 'Original website', exact: false }).first().click();
    assert.ok(modern.url().endsWith('table_beta_legacy.html'));
    await modern.getByRole('link', { name: 'Go to the modern website' }).click();
    assert.ok(modern.url().endsWith('table_beta.html')); checks++;
    const ids = await modern.evaluate(() => Array.from(document.querySelectorAll('[id]')).map(el => el.id)); assert.equal(new Set(ids).size, ids.length); checks++;
    await modern.goto(base + '/index.html'); await modern.waitForURL(base + '/table_beta.html'); checks++;
    await modern.keyboard.press('Tab');
    assert.equal(await modern.evaluate(() => document.activeElement.textContent), 'Skip to calculator'); checks++;
    await modern.keyboard.press('Enter'); assert.equal(await modern.locator('#workspace').evaluate(el => el === document.activeElement), true); checks++;
    // File URLs exercise the no-server, double-click workflow too.
    const fileContext = await browser.newContext();
    const filePage = await fileContext.newPage();
    filePage.on('pageerror', error => errors.push(error.message));
    for (const filename of ['table_beta.html', 'table_beta_legacy.html']) {
      await filePage.goto(pathToFileURL(path.resolve(__dirname, '..', filename)).href);
      assert.equal((await results(filePage, false, filename === 'table_beta.html'))[0][0].value, 0.307); checks++;
    }
    await fileContext.close();
    assert.deepEqual(errors, [], 'uncaught browser errors'); assert.deepEqual(blocked, [], 'external network dependencies');
    console.log(`PASS: ${cases.length} old/new browser scenarios; ${cells} table values and highlight decisions; all distribution curves and conversions match. ${checks} integration checks passed. No uncaught errors or application external requests. ${environmentRequests.length} host AdGuard requests blocked separately. Browser: ${browser.version()}.`);
  } finally { await browser.close(); server.close(); }
})().catch(error => { console.error(error); server.close(); process.exitCode = 1; });
