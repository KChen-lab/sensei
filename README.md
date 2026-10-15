# Sensei

Sample size estimation for single-cell studies. The website is static HTML, CSS, and JavaScript: no build step, backend, or upload of study parameters.

- [Calculator](https://kchen-lab.github.io/sensei/table_beta.html)
- [Methods paper](https://link.springer.com/article/10.1186/s12859-021-04526-5)
- [Local testing and manual comparison cases](TESTING.md)

The modern calculator stays at `table_beta.html`. Its “Original website” link opens `table_beta_legacy.html`, which preserves the previous calculator's scientific and interaction code. Both pages use the same bundled statistical and plotting libraries. Earlier experimental pages (`single_beta.html`, `table_beta_old.html`, and `table_beta_unpaired_only.html`) remain available and unchanged; some use an older normal approximation and are not the compatibility baseline.

To preview locally with Node.js, run `node scripts/serve.cjs`, then open http://127.0.0.1:8000/table_beta.html. Alternatively use `python -m http.server 8000 --bind 127.0.0.1` from this directory. Directly opening the HTML files also works. Installation is needed only for automated browser tests; see `TESTING.md`.

The `python/` and `matlab/` implementations and research data are unchanged.
