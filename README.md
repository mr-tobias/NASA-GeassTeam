# GEASS Team — MODIS Vegetation Trend Project Website

A responsive, static project website documenting the supplied Python NDVI analysis application.

## Files
- `index.html` — website content and layout
- `style.css` — responsive visual design
- `script.js` — mobile navigation
- `src/ndvi_trend_analyzer.py` — supplied Python application source

## Preview locally
Open `index.html` in a modern browser. No build step is required.

For a local server, from this folder run:
```bash
python -m http.server 8000
```
Then visit `http://localhost:8000`.

## Deploy
Upload the contents of this folder to any static host (GitHub Pages, Netlify, or compatible static hosting). Keep the `src` folder in place so the source-code download link works.

## Python application requirements
The desktop application imports:
- numpy
- rasterio
- scipy
- matplotlib

Install them in a Python environment with:
```bash
python -m pip install numpy rasterio scipy matplotlib
```

Run:
```bash
python src/ndvi_trend_analyzer.py
```

The script expects NDVI and Pixel Reliability GeoTIFF files in the configured `DATA_DIR`. Before running, review the configuration at the top of the script:
- `DATA_DIR` currently points to `data/modis_fake`; replace it with the directory containing your actual files.
- `NDVI_TAG`, `QA_TAG`, and `DATE_RE` must match the filenames of your downloaded data.
- The script expects matching NDVI and QA files with the same `doyYYYYDDD` date token.

This website is a documentation/showcase page. It does **not** run the Matplotlib desktop GUI in the browser.

## Important scientific and technical notes
- No analysis results or trend values are fabricated or included.
- NDVI is a greenness indicator and does not, by itself, prove vegetation health or explain causes.
- Verify the raster CRS and coordinate transforms before interpreting map locations. The supplied code's map labels assume geographic coordinates, which may not match a MODIS sinusoidal grid.
- The OLS autocorrelation adjustment is approximate; assess its assumptions before making inferential claims.
- The regional Kendall test is not adjusted for serial correlation.
- The supplied source was syntax-checked before packaging, but its runtime behavior has not been validated against the user's complete dataset.

## Credits
Project content and Python source supplied by the user. Visual layout prepared as a static project website.
\n\nThe footer-area attribution section includes the NASA Space Apps Challenge and NASA logos supplied for this project. They identify the event context and do not imply endorsement.\n