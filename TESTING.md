# Review the website locally

Nothing has been pushed or deployed. The published website will stay unchanged until you choose to publish this checkout.

## Open the calculator

The fastest option is to double-click `table_beta.html` in this folder. Open `table_beta_legacy.html` in a second tab to compare. Both pages include their dependencies locally.

For a local HTTP preview, run this from `C:\Files\projects\sensei` in PowerShell using the Python installation already on this computer:

```powershell
& 'C:\Users\lshh\miniconda3\envs\py313_torch291_cuda130\python.exe' -m http.server 8000 --bind 127.0.0.1
```

Then open:

- Modern: http://127.0.0.1:8000/table_beta.html
- Original: http://127.0.0.1:8000/table_beta_legacy.html

Keep that terminal open while testing. Press Ctrl+C afterward. If port 8000 is occupied, use 8001 in both the command and URLs. On other computers, `python -m http.server 8000 --bind 127.0.0.1` or `node scripts/serve.cjs` works.

The “Original website” link opens the preserved calculator. Its banner links back to the modern one. Inputs are independent between the two pages; the links do not transfer your parameters.

## Manual numerical checks

Start each case from defaults. In the modern page, click **Reset defaults** in the active tab; in the original, reload the page, then choose the relevant tab. Unless specified otherwise, retain 1,000 cells per sample, α = 0.05, βth = 0.1, control a = 7 / b = 10, experimental a = 10 / b = 7. Paired correlation defaults to 0.5.

| Case | What to change | Cell to check | Expected β on both pages |
| --- | --- | --- | --- |
| Default unpaired | Nothing | Control 5, experimental 5 | **0.307**, unhighlighted |
| Larger unpaired | Nothing | Control 12, experimental 12 | **0.03**, highlighted |
| Default paired | Choose Paired t-test | 5 matched pairs | **0.146**, unhighlighted |
| Two-sided | Unpaired, select Two-sided t-test | Control 5, experimental 5 | **0.469** |
| SESOI verify | Unpaired, select SESOI; offset 0.01; Verify | Control 5, experimental 5 | **0.353** |
| SESOI falsify | Same as above, select Falsify | Control 5, experimental 5 | **0.998** |
| TOST | Unpaired, select TOST; lower −0.1, upper 0.1 | Control 5, experimental 5 | **1** |
| Rare-cell example | Unpaired one-sided; control mean 0.03 / SD 0.015; experimental mean 0.05 / SD 0.01 | Control 5, experimental 5; then 12 and 12 | **0.366**; **0.049** |
| No expected difference | Unpaired one-sided; set experimental a = 7, b = 10 in beta parameters | Control 5, experimental 5 | **0.95** |

For the rare-cell case, changing mean/SD must yield control a = 3.85, b = 124.483 and experimental a = 23.7, b = 450.3. In the original page, commit each edit with Tab or click outside the field. New-page inputs update immediately. Ignore intermediate results while editing a related pair of fields.

**Rounding matters:** defaults have exact a/b values. The displayed mean (0.4118 or 0.5882) and SD (0.116) are rounded. Re-entering these displayed numbers converts them back to slightly different a/b, just as in the original. For exact comparisons, reset defaults or enter identical a/b on both pages. Mean/SD → a/b rounds to three decimal places; a/b → mean/SD rounds to four.

## Interaction checks

1. Compare every value and highlighted cell in the default 8 × 8 unpaired table. **26 of 64 designs** should be highlighted. The modern smallest-total-sample summary should show **8 control + 9 experimental samples** in this range.
2. Change only βth from 0.1 to 0.2. Table values must stay the same; more cells should be highlighted. Change α to 0.01; false negative rates should rise for otherwise identical parameters.
3. In Paired t-test, change correlation from 0.5 to 0, then to 0.9. With the default unequal means, higher positive correlation should reduce false negatives. Change cell counts independently in both groups and compare the two pages.
4. Switch tabs and return. Each tab keeps its own inputs. Use browser Back/Forward and direct `#paired` / `#instructions` links. Reset changes only the active calculator.
5. Download CSV in both study modes. The default unpaired file has 64 data rows; paired has 8. It includes exact a/b, cell counts, sample counts, α, βth, relevant offsets/correlation, unrounded false negative rate, power, and threshold status. The file retains more precision than the table.
6. On the modern page, try zero cells, sample size 1, an inverted range, mean 0, an impossible SD, α = 1, or correlation 1.1. Results must clear with an error and CSV must be disabled. Correct the input or reset; results must recover. Blank active fields behave the same way.
7. Reduce the browser width to about 390 pixels. Inputs stack above results, the plot resizes, and the table scrolls horizontally within its panel. The entire page should fit the screen. Use Tab to check visible focus indicators and the “Skip to calculator” link.
8. Use the spinner buttons or Up/Down keys in numeric inputs: decimal parameters move in steps of **0.01**, while sample sizes and cell counts move by **1**. You can still type finer precision, including the rare-cell SDs above.
9. Hover or tap a table value to see a styled card with separate lines for sample counts, false negative rate, power, and threshold status. The card should fit the screen in both study modes. Use Tab to focus a value, arrow keys to move between cells, and Escape to dismiss the card. Changing parameters or leaving the table should clear stale details.

The popup shows β to three decimal places, or through the first nonzero digit for small values (up to six decimal places, then scientific notation). For example: `0.307`, `0.00004`, `3e-9`. This only changes display; calculations, highlighting, and CSV precision remain unchanged. The paper and original-site links are grouped below the title bar with a brief note describing the visual update.

The modern page adds validation for invalid distributions, noninteger sample/cell counts, ranges below 2, α/βth outside (0,1), correlation outside [−1,1], zero SESOI, and invalid TOST bounds. It caps the displayed table at 10,000 result cells. These input-handling differences are intentional; valid supported study designs use the original calculations. Positive a/b values ≤ 1 are accepted with a model-assumption note; infinite endpoint densities are omitted from the chart. The legacy UI and its validation behavior remain as before.

## Automated regression tests

Node.js 20+ runs the numerical tests without installing dependencies:

```powershell
node --test tests/core.test.cjs
```

For this computer's bundled Node installation, if `node` is not on PATH:

```powershell
& 'C:\Users\lshh\.cache\codex-runtimes\codex-primary-runtime\dependencies\node\bin\node.exe' --test tests/core.test.cjs
```

The test loads the original functions independently from the legacy HTML, checks their frozen source hashes, compares 10,000 deterministic generated designs and 1,600 boundary/mode cases, and checks parameter conversions. It writes `.artifacts/reference-vectors.json`. Run the independent reference comparison afterward with Python and SciPy (the environment above already has them):

```powershell
& 'C:\Users\lshh\miniconda3\envs\py313_torch291_cuda130\python.exe' tests/reference_check.py
```

Browser tests use Playwright only as a development dependency; it is never loaded by the website:

```powershell
npm install
npm run test:browser
```

On Windows these tests use the installed Edge in an isolated headless session. On other systems, install Playwright Chromium with `npx playwright install chromium` first. To choose Chrome on Windows, set `$env:PLAYWRIGHT_CHANNEL = 'chrome'` before running the suite.

This computer already has Playwright in the bundled runtime, so it can also run without installation:

```powershell
$env:NODE_PATH = 'C:\Users\lshh\.cache\codex-runtimes\codex-primary-runtime\dependencies\node\node_modules'
& 'C:\Users\lshh\.cache\codex-runtimes\codex-primary-runtime\dependencies\node\bin\node.exe' tests/browser.test.cjs
```

The suite starts and stops its own loopback server, compares 25 old/new browser scenarios, exercises real inputs, validation/recovery, exact threshold behavior, downloads, navigation, direct-file loading, and 320/390/768/1440-pixel widths. It saves screenshots in `.artifacts/`. All external requests are blocked during HTTP tests. Requests injected by this computer's system AdGuard are still blocked and are reported separately from application requests; no security settings are changed.

## Compatibility baseline and scope

`table_beta_legacy.html` is the website previously served at `table_beta.html`, not the older experimental file already named `table_beta_old.html`. Its calculation and interaction scripts were checked against both the initial Git checkout and the live GitHub Pages HTML on 2026-10-14: both are identical. Only three dependency URLs and a return-navigation banner changed. Libraries are bundled at the versions resolved by the original URLs; see `assets/vendor/README.md` for versions and hashes.

The scientific approximation, SESOI directions, TOST sum capped at 1, defaults, conversion rounding, three-decimal table display, and strict unrounded β < βth highlighting are preserved. This is a compatibility check, not a new validation of the underlying scientific model. The Python and MATLAB research code has not been modified. Earlier experimental HTML pages retain their older approximation and are not expected to match this calculator.

The design follows the navy tab bar, light background, white panels, system typography, and blue controls of `financial_widgets`. No framework, build pipeline, server computation, or remote deployment was added.

## Recorded verification — 2026-10-14

- Numerical suite: **5 tests passed**. All **11,600** generated and boundary/mode comparisons matched the original implementation exactly at full JavaScript precision. Conversion and frozen-source checks passed.
- Independent Python/SciPy comparison: **1,800 vectors passed**; maximum absolute difference **3.71 × 10⁻⁸**, below the 2 × 10⁻⁷ test tolerance.
- Browser suite: **25 comparison scenarios**, **798 table values and highlight decisions**, matching distribution curves and conversion fields, and **64 integration checks passed** in headless Microsoft Edge **155.0.4283.45**. This includes actual spinner-button clicks, decimal/integer keyboard steps, multiline tooltips in both study modes, and popup bounds at all four screen widths. No uncaught application errors or application external requests occurred. The host injected 110 AdGuard requests, all blocked and reported separately.
- Direct-file opening passed for both pages. Layout and chart resizing passed at **320, 390, 768, and 1440 pixels**. Keyboard navigation, reset/recovery, CSV contents, tab history, root redirect, and old/new navigation passed.
- Both legacy scripts matched the live published calculator and the initial checkout. Local asset/page links resolved. `git diff --check` passed.
- Screenshots: `.artifacts/sensei-desktop.png`, `.artifacts/sensei-320.png`, `.artifacts/sensei-390.png`, `.artifacts/sensei-768.png`, `.artifacts/sensei-1440.png`, `.artifacts/sensei-table-tooltip.png`, and `.artifacts/sensei-tooltip-mobile.png`. This ignored folder contains local verification outputs, not website dependencies.

Browser automation was performed in Edge; Safari and Firefox were not run. There has been no remote push or deployment.
