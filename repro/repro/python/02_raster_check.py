"""Raster vs station check. Run from the repository root after 01_descriptives_hotspots_grid.py. Set RASTER and CLIMATEPLUS_KIT if your files are elsewhere."""
import numpy as np, pandas as pd, rasterio
from scipy.stats import spearmanr, pearsonr, rankdata
import os
RASTER = os.environ.get('RASTER','data/pollution_int.tif')
GRID = os.path.join(os.environ.get('CLIMATEPLUS_KIT','data/kit'),'research/gis/kaggle-v2/Pollution_Intensity_Coordinates_Updated.csv')
r = rasterio.open(RASTER); a = r.read(1); a = np.where(a < -1e30, np.nan, a)
g = pd.read_csv(GRID, encoding='utf-8-sig'); st = pd.read_csv('results/stations.csv')
v = np.array([x[0] for x in r.sample(zip(g.Longitude, g.Latitude))]); v = np.where(v < -1e30, np.nan, v); ok = np.isfinite(v)
d = g.Pollution_Intensity_Value.values[ok] - v[ok]
print('raster', r.width, 'x', r.height, 'px, res deg', r.res[0], '| valid px', int(np.isfinite(a).sum()))
print('grid vs raster: equal share', round(float((abs(d) < 1e-4).mean()), 4), '| pearson', round(pearsonr(g.Pollution_Intensity_Value.values[ok], v[ok])[0], 3))
def samp(lon, lat, k):
    i, j = r.index(lon, lat); w = a[max(i-k, 0):i+k+1, max(j-k, 0):j+k+1]; return np.nanmean(w)
for km, k in ((0, 0), (1, 2), (2, 5)):
    st['rv'] = [samp(x, y, max(int(round(km / 0.43)), 0)) for x, y in zip(st.lon, st.lat)]
    print('buffer', km, 'km', {m: tuple(round(float(q), 3) for q in spearmanr(st.rv, st[m])) for m in ('p2', 'DOmed', 'WQI')})
st['rv'] = [samp(x, y, 5) for x, y in zip(st.lon, st.lat)]
print('LORO range', round(min(spearmanr(st[st.River != q].rv, st[st.River != q].p2)[0] for q in st.River.unique()), 3),
      round(max(spearmanr(st[st.River != q].rv, st[st.River != q].p2)[0] for q in st.River.unique()), 3))
st['dc'] = np.hypot((st.lat-23.78)*111, (st.lon-90.40)*111*np.cos(np.radians(23.8)))
def partial(x, y, z):
    rx, ry, rz = [rankdata(q) for q in (x, y, z)]
    return np.corrcoef(rx-np.polyval(np.polyfit(rz, rx, 1), rz), ry-np.polyval(np.polyfit(rz, ry, 1), rz))[0, 1]
print('partial rho given distance to centre: hypoxia', round(partial(st.rv, st.p2, st.dc), 3), '| DO', round(partial(st.rv, st.DOmed, st.dc), 3))
rng = np.random.default_rng(3); obs = spearmanr(st.rv, st.p2)[0]
c = sum(spearmanr(st.rv, st.assign(p2=st.groupby('River').p2.transform(lambda s: rng.permutation(s.values))).p2)[0] >= obs for _ in range(2000))
print('within-river permutation p', (c+1)/2001)
