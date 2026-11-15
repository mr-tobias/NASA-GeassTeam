# GEASS Team — Earth Trend Detective

Static project website for the NASA Space Apps Challenge project.

<<<<<<< Updated upstream
This website is a documentation/showcase page. It does **not** run the Matplotlib desktop GUI in the browser.

## Important scientific and technical notes
- No analysis results or trend values are fabricated or included.
- NDVI is a greenness indicator and does not, by itself, prove vegetation health or explain causes.
- Verify the raster CRS and coordinate transforms before interpreting map locations. The supplied code's map labels assume geographic coordinates, which may not match a MODIS sinusoidal grid.
- The OLS autocorrelation adjustment is approximate; assess its assumptions before making inferential claims.
- The regional Kendall test is not adjusted for serial correlation.
- The supplied source was syntax-checked before packaging, but its runtime behavior has not been validated against the user's complete dataset.

=======
## Files
- `index.html` — website
- `style.css` — dark charcoal/crimson styling
- `assets/ndvi-trend-map.png` — supplied NDVI trend map, used as the interface image
- `Earth-Trend-Detective-Project-Report.pdf` — supplied 8-page project report, linked by download buttons

## Publish with GitHub Pages
1. Upload the *contents* of this folder to the root of your repository.
2. Confirm `index.html` is at the repository root and keep the `assets` folder and PDF beside it.
3. In GitHub, open **Settings → Pages**.
4. Choose **Deploy from a branch**, select `main` and `/(root)`, then save.

The report PDF is included unchanged from the supplied file. The screenshot is also used as supplied. The screenshot displays a 2000–2012 window, while the report describes the broader 2001–2025 project period; the page explicitly distinguishes them.
>>>>>>> Stashed changes
