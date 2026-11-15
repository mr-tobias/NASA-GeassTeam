import csv
import glob
import os
import re
import textwrap
import warnings
from collections import defaultdict
from datetime import datetime

import numpy as np
import rasterio
from scipy import stats


# Configuration


DATA_DIR = "data/modis_fake"                
NDVI_TAG = "250m_16_days_NDVI"          
QA_TAG = "250m_16_days_pixel_reliability"
DATE_RE = re.compile(r"doy(\d{4})(\d{3})")   
SCALE = 0.0001                          # MODIS NDVI scale factor
VALID_RANGE = (-0.2, 1.0)               # physically sensible NDVI after scaling
GOOD_QA = (0,)                          
MIN_COMPOSITES = 12                     
MAX_SIDE = 500                          
ALPHA = 0.05
MIN_VALID_FRAC = 0.6
BACKEND = "TkAgg"
AGGREGATIONS = {"Annual mean": np.nanmean, "Annual peak": np.nanmax}
OUTPUT_DIR = "outputs"                  
CITIES = [("San Juan", -31.5375, -68.5364)]   


# Data: MODIS GeoTIFF time series


def block_mean(a, k):
    
    if k == 1:
        return a
    h, w = (a.shape[0] // k) * k, (a.shape[1] // k) * k
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmean(a[:h, :w].reshape(h // k, k, w // k, k), axis=(1, 3))


def load_modis_pair(ndvi_path, qa_path, block=1):
    with rasterio.open(ndvi_path) as n_src, rasterio.open(qa_path) as q_src:
        if (n_src.width, n_src.height, n_src.crs, n_src.transform) != \
           (q_src.width, q_src.height, q_src.crs, q_src.transform):
            raise ValueError(f"NDVI and QA not aligned: {ndvi_path}")
        ndvi = n_src.read(1, masked=True).astype("float32") * SCALE
        qa = q_src.read(1, masked=True)
        bad = (np.ma.getmaskarray(ndvi) | np.ma.getmaskarray(qa)
               | ~np.isin(qa.data, GOOD_QA)
               | (ndvi.data < VALID_RANGE[0]) | (ndvi.data > VALID_RANGE[1]))
        out = np.where(bad, np.nan, ndvi.data).astype("float32")
        meta = {"crs": n_src.crs, "transform": n_src.transform * n_src.transform.scale(block)}
    return block_mean(out, block), meta


def find_pairs(folder):
    
    qa_by_date = {}
    for p in glob.glob(os.path.join(folder, f"*{QA_TAG}*.tif")):
        m = DATE_RE.search(os.path.basename(p))
        if m:
            qa_by_date[m.groups()] = p
    by_year = defaultdict(list)
    for p in sorted(glob.glob(os.path.join(folder, f"*{NDVI_TAG}*.tif"))):
        m = DATE_RE.search(os.path.basename(p))
        if m and m.groups() in qa_by_date:
            by_year[int(m.group(1))].append((p, qa_by_date[m.groups()]))
    return by_year


def load_annual_ndvi(folder=DATA_DIR):
    
    by_year = find_pairs(folder)
    if len(by_year) < 12:
        raise SystemExit(f"Found {len(by_year)} usable years in '{folder}'. The trend tool needs a "
                         "time series (>= 12 years), not a single composite. Download the full "
                         "NDVI + pixel_reliability series from AppEEARS for your San Juan area.")
    with rasterio.open(next(iter(by_year.values()))[0][0]) as s:
        block = max(1, int(np.ceil(max(s.width, s.height) / MAX_SIDE)))
    years = np.array(sorted(by_year))
    annual = {name: [] for name in AGGREGATIONS}
    meta = None
    for y in years:
        comps = []
        for ndvi_p, qa_p in by_year[y]:
            arr, meta = load_modis_pair(ndvi_p, qa_p, block)
            comps.append(arr)
        cube = np.stack(comps)                      # (composites, H, W)
        count = np.isfinite(cube).sum(axis=0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            for name, fn in AGGREGATIONS.items():
                a = fn(cube, axis=0)
                annual[name].append(np.where(count >= MIN_COMPOSITES, a, np.nan).astype("float32"))
        del cube, comps
        print(f"  {y}: {len(by_year[y])} composites")
    return years, {k: np.stack(v) for k, v in annual.items()}, meta



# Statistics (unchanged logic from the temperature version)


def trend_stats(arr, x, min_frac=MIN_VALID_FRAC):
    arr = np.asarray(arr, dtype=float)
    x = np.asarray(x, dtype=float)
    ok = np.isfinite(arr)
    n = ok.sum(axis=0)
    with warnings.catch_warnings(), np.errstate(invalid="ignore", divide="ignore"):
        warnings.simplefilter("ignore", RuntimeWarning)
        xx = np.where(ok, x[:, None, None], np.nan)
        xm = np.nanmean(xx, axis=0)
        ym = np.nanmean(arr, axis=0)
        dx = xx - xm
        sxx = np.nansum(dx ** 2, axis=0)
        b = np.nansum(dx * (arr - ym), axis=0) / sxx

        resid = arr - (ym + b * dx)
        ss = np.nansum(resid ** 2, axis=0)
        r1 = np.nansum(resid[1:] * resid[:-1], axis=0) / ss
        r1 = np.clip(np.nan_to_num(r1), 0.0, 0.99)
        n_eff = np.minimum(np.maximum(n * (1 - r1) / (1 + r1), 4.0), n)

        se = np.sqrt(ss / (n_eff - 2) / sxx)
        p = 2 * stats.t.sf(np.abs(b / se), n_eff - 2)

    bad = (n < max(5, int(min_frac * len(x)))) | ~np.isfinite(b) | ~np.isfinite(p)
    return np.where(bad, np.nan, b), np.where(bad, np.nan, p), n


def fdr_mask(p, alpha=ALPHA):
    sig = np.zeros(p.shape, dtype=bool)
    valid = np.isfinite(p)
    pv = p[valid]
    m = pv.size
    if m == 0:
        return sig
    ps = np.sort(pv)
    passed = np.where(ps <= alpha * np.arange(1, m + 1) / m)[0]
    if passed.size:
        sig[valid] = pv <= ps[passed.max()]
    return sig


def region_series(annual, row, col, half_px):
   
    r0, r1 = max(0, row - half_px), min(annual.shape[1], row + half_px + 1)
    c0, c1 = max(0, col - half_px), min(annual.shape[2], col + half_px + 1)
    if r0 >= r1 or c0 >= c1:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmean(annual[:, r0:r1, c0:c1], axis=(1, 2)), (r0, r1, c0, c1)


def analyze_series(x, y):
    ok = np.isfinite(y)
    xs, ys = np.asarray(x, float)[ok], np.asarray(y, float)[ok]
    if len(xs) < 8:
        return None
    b, p, _ = trend_stats(ys[:, None, None], xs, min_frac=0.0)
    ts = stats.theilslopes(ys, xs, 0.95)
    tau, p_mk = stats.kendalltau(xs, ys)
    return {
        "n": len(xs),
        "ols": float(b[0, 0]), "p_ols": float(p[0, 0]),
        "ols_x": float(xs.mean()), "ols_y": float(ys.mean()),
        "ts": float(ts[0]), "ts_icpt": float(ts[1]),
        "ts_lo": float(ts[2]), "ts_hi": float(ts[3]),
        "tau": float(tau), "p_mk": float(p_mk),
    }


def plain_summary(res, start, end, alpha=ALPHA):
    
    if res is None:
        return "Not enough valid data in this region and window.", "0.4"
    direction = "greening" if res["ts"] > 0 else "browning"
    amt = abs(res["ts"] * 10)
    s_ols, s_mk = res["p_ols"] < alpha, res["p_mk"] < alpha
    colour = "darkgreen" if res["ts"] > 0 else "saddlebrown"
    if s_ols and s_mk:
        return (f"This region is {direction}: NDVI changed by about {amt:.3f} per decade "
                f"between {start} and {end} (statistically significant).", colour)
    if s_ols or s_mk:
        return (f"Possible {direction} (about {amt:.3f} NDVI per decade), but the two tests "
                "disagree. Treat with caution.", "darkorange")
    return ("No clear change in NDVI in this window. A longer record might change this.", "0.4")


def nice_limit(v):
    return max(0.02, float(np.ceil(v * 100) / 100)) if np.isfinite(v) else 0.1


def pixel_size_m(meta, shape):
    t = meta["transform"]
    if meta["crs"] is not None and meta["crs"].is_geographic:
        lat_mid = t.f + t.e * shape[0] / 2
        return abs(t.a) * 111320 * np.cos(np.deg2rad(lat_mid))
    return abs(t.a)



# Interactive app


def main():
    import matplotlib
    matplotlib.use(BACKEND)
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.widgets import Button, CheckButtons, RadioButtons, Slider

    print("Loading MODIS NDVI series (one year at a time)...")
    years, cubes, meta = load_annual_ndvi()
    agg_names = list(AGGREGATIONS)
    H, W = next(iter(cubes.values())).shape[1:]
    t = meta["transform"]
    extent = (t.c, t.c + t.a * W, t.f + t.e * H, t.f)       
    px_m = pixel_size_m(meta, (H, W))

    state = {"pt": (H // 2, W // 2)}                         

    fig = plt.figure(figsize=(14, 7.5))
    ax = fig.add_axes([0.04, 0.34, 0.58, 0.60])
    cax = fig.add_axes([0.12, 0.235, 0.40, 0.02])
    ax_ts = fig.add_axes([0.70, 0.52, 0.28, 0.40])
    ax_info = fig.add_axes([0.70, 0.17, 0.28, 0.28])
    ax_info.axis("off")

    blank = np.full((H, W), np.nan)
    mesh = ax.imshow(blank, extent=extent, cmap="RdYlGn", vmin=-0.1, vmax=0.1,
                     origin="upper", interpolation="nearest")
    cbar = fig.colorbar(mesh, cax=cax, orientation="horizontal")
    overlay = ax.imshow(blank, extent=extent, cmap=ListedColormap(["white"]), vmin=0, vmax=1,
                        alpha=0.65, origin="upper", interpolation="nearest")
    (box_line,) = ax.plot([], [], color="k", lw=1.6)
    title = ax.set_title("")
    if meta["crs"] is not None and meta["crs"].is_geographic:
        ax.set_aspect(1 / np.cos(np.deg2rad(t.f + t.e * H / 2)))   # true-shape map
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    for label, clat, clon in CITIES:                          # reference points
        if extent[0] <= clon <= extent[1] and extent[2] <= clat <= extent[3]:
            ax.plot(clon, clat, marker="*", ms=13, color="k", mec="white", zorder=5)
            ax.annotate(label, (clon, clat), xytext=(6, 6), textcoords="offset points",
                        fontsize=9, fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, lw=0))

    s_end = Slider(fig.add_axes([0.12, 0.14, 0.42, 0.025]), "End year",
                   int(years[0]) + 9, int(years[-1]), valinit=int(years[-1]), valstep=1)
    s_win = Slider(fig.add_axes([0.12, 0.09, 0.42, 0.025]), "Window (yrs)",
                   10, len(years), valinit=min(len(years), 20), valstep=1)
    s_half = Slider(fig.add_axes([0.12, 0.04, 0.42, 0.025]), "Region ± (km)",
                    1, 25, valinit=5, valstep=1)
    radio = RadioButtons(fig.add_axes([0.60, 0.03, 0.09, 0.14]), agg_names, active=0)
    radio.ax.set_title("NDVI metric", fontsize=9)
    chk = CheckButtons(fig.add_axes([0.70, 0.04, 0.27, 0.10]),
                       ["Grey out non-significant cells",
                        "FDR correction (multiple tests)"], [True, True])
    btn_save = Button(fig.add_axes([0.60, 0.20, 0.09, 0.04]), "Save PNG + CSV")
    status = fig.text(0.70, 0.015, "", fontsize=8, color="0.3")
    fig.text(0.04, 0.005, "Click the map to choose a region  ·  sliders change the time window "
             "and region size  ·  Save exports the figure and the region's yearly values",
             fontsize=8, color="0.4")

    def refresh(_=None):
        end, win, half_km = int(s_end.val), int(s_win.val), float(s_half.val)
        metric = radio.value_selected
        vals = cubes[metric]
        start = max(int(years[0]), end - win + 1)
        sel = (years >= start) & (years <= end)
        x = years[sel].astype(float)

        # ---- map ----
        b, p, _ = trend_stats(vals[sel], x)
        valid = np.isfinite(b)
        show_grey, use_fdr = chk.get_status()
        with np.errstate(invalid="ignore"):
            sig = fdr_mask(p, ALPHA) if use_fdr else (p < ALPHA)
        data = b * 10                                         # NDVI change per decade
        lim = nice_limit(np.nanpercentile(np.abs(data[valid]), 98)) if valid.any() else 0.1
        mesh.set_data(np.ma.masked_invalid(data))
        mesh.set_clim(-lim, lim)
        overlay.set_data(np.ma.masked_where(~(valid & ~sig), np.ones_like(data)))
        overlay.set_visible(show_grey)
        cbar.set_label(f"Change in NDVI ({metric.lower()}) per decade  [greening > 0]")
        frac = 100 * (valid & sig).sum() / max(valid.sum(), 1)
        warn = "   ⚠ short window: mostly noise" if win < 15 else ""
        title.set_text(f"NDVI trend {start}–{end}  |  {frac:.0f}% of cells significant"
                       f" at {int((1 - ALPHA) * 100)}%{warn}")

        # ---- region ----
        ax_ts.clear()
        ax_info.clear()
        ax_info.axis("off")
        row, col = state["pt"]
        half_px = max(0, int(round(half_km * 1000 / px_m)))
        out = region_series(vals, row, col, half_px)
        if out is None:
            ax_ts.text(0.5, 0.5, "Click the map to analyse a region",
                       ha="center", va="center", transform=ax_ts.transAxes)
            box_line.set_data([], [])
        else:
            y, (r0, r1, c0, c1) = out
            xl, xr = t.c + t.a * c0, t.c + t.a * c1
            yt, yb = t.f + t.e * r0, t.f + t.e * r1
            box_line.set_data([xl, xr, xr, xl, xl], [yb, yb, yt, yt, yb])
            res = analyze_series(x, y[sel])
            ax_ts.plot(years, y, color="0.7", lw=1, marker="o", ms=2.5, label="all years")
            ax_ts.plot(years[sel], y[sel], color="tab:green", lw=1.3, marker="o", ms=3,
                       label=f"{start}–{end}")
            if res:
                xs = np.array([start, end], float)
                ax_ts.plot(xs, res["ols_y"] + res["ols"] * (xs - res["ols_x"]),
                           "k-", lw=2, label="OLS trend")
                ax_ts.plot(xs, res["ts_icpt"] + res["ts"] * xs, "b--", lw=1.6, label="Theil–Sen")
            lon0, lat0 = t.c + t.a * (col + 0.5), t.f + t.e * (row + 0.5)
            ax_ts.set_title(f"Region {lat0:.2f}°, {lon0:.2f}°  (±{half_km:.0f} km)", fontsize=10)
            state["last"] = dict(metric=metric, lat=lat0, lon=lon0, half_km=half_km,
                                 start=start, end=end, y=y.copy(), res=res)
            ax_ts.set_ylabel(f"NDVI ({metric.lower()})")
            ax_ts.grid(alpha=0.3)
            ax_ts.legend(fontsize=8, loc="upper left")
            summary, colour = plain_summary(res, start, end)
            ax_info.text(0, 1, textwrap.fill(summary, 44), va="top", ha="left", fontsize=10,
                         fontweight="bold", color=colour, transform=ax_info.transAxes)
            if res:
                ax_info.text(0, 0.58, (
                    f"Window {start}–{end} ({res['n']} valid years)\n"
                    f"OLS (autocorr-adj.): {res['ols'] * 10:+.3f} NDVI/dec, p={res['p_ols']:.3g}\n"
                    f"Theil–Sen: {res['ts'] * 10:+.3f} NDVI/dec\n"
                    f"   95% CI [{res['ts_lo'] * 10:+.3f}, {res['ts_hi'] * 10:+.3f}]\n"
                    f"Mann–Kendall: tau={res['tau']:+.2f}, p={res['p_mk']:.3g}"),
                    va="top", ha="left", fontsize=9, family="monospace",
                    transform=ax_info.transAxes)
        fig.canvas.draw_idle()

    def on_click(ev):
        if ev.inaxes is not ax or ev.button != 1:
            return
        tb = getattr(fig.canvas, "toolbar", None)
        if tb is not None and getattr(tb, "mode", "") != "":
            return
        col = int((ev.xdata - t.c) / t.a)
        row = int((ev.ydata - t.f) / t.e)
        if 0 <= row < H and 0 <= col < W:
            state["pt"] = (row, col)
            refresh()

    def save_current(_):
        last = state.get("last")
        if last is None:
            status.set_text("Click the map to choose a region first.")
            fig.canvas.draw_idle()
            return
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        base = os.path.join(OUTPUT_DIR, "ndvi_region_" + datetime.now().strftime("%Y%m%d_%H%M%S"))
        fig.savefig(base + ".png", dpi=150)
        res = last["res"]
        with open(base + ".csv", "w", newline="", encoding="utf-8") as f:
            f.write(f"# NDVI region analysis, metric: {last['metric']}\n")
            f.write(f"# centre lat/lon: {last['lat']:.4f}, {last['lon']:.4f}; "
                    f"half-width: {last['half_km']:.0f} km\n")
            f.write(f"# window: {last['start']}-{last['end']}\n")
            if res:
                f.write(f"# OLS slope: {res['ols'] * 10:+.4f} NDVI/decade (p={res['p_ols']:.3g}); "
                        f"Theil-Sen: {res['ts'] * 10:+.4f} NDVI/decade; "
                        f"Mann-Kendall p={res['p_mk']:.3g}\n")
            w = csv.writer(f)
            w.writerow(["year", "region_mean_ndvi", "in_window"])
            for yr, v in zip(years, last["y"]):
                w.writerow([int(yr), f"{v:.5f}" if np.isfinite(v) else "",
                            int(last["start"] <= yr <= last["end"])])
        status.set_text(f"Saved {base}.png and .csv")
        fig.canvas.draw_idle()

    btn_save.on_clicked(save_current)
    for s in (s_end, s_win, s_half):
        s.on_changed(refresh)
    radio.on_clicked(refresh)
    chk.on_clicked(refresh)
    fig.canvas.mpl_connect("button_press_event", on_click)

    refresh()
    plt.show()


if __name__ == "__main__":
    main()
