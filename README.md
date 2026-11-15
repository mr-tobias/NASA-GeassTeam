# GEASS Team — MODIS Vegetation Trend Project Website

A responsive, static project website documenting the supplied Python NDVI analysis application.

This website is a documentation/showcase page. It does **not** run the Matplotlib desktop GUI in the browser.

## Important scientific and technical notes
- No analysis results or trend values are fabricated or included.
- NDVI is a greenness indicator and does not, by itself, prove vegetation health or explain causes.
- Verify the raster CRS and coordinate transforms before interpreting map locations. The supplied code's map labels assume geographic coordinates, which may not match a MODIS sinusoidal grid.
- The OLS autocorrelation adjustment is approximate; assess its assumptions before making inferential claims.
- The regional Kendall test is not adjusted for serial correlation.
- The supplied source was syntax-checked before packaging, but its runtime behavior has not been validated against the user's complete dataset.

