// =====================================================================
// Sentinel-2 dry-season composite sampled at the 60 stations (to see what the
// pollution raster follows: vegetation, built-up land or water)
// 1. Upload stations_for_gee.csv as a table asset (X = lon, Y = lat).
// 2. Paste in the Code Editor, edit the asset ID, press Run, then run the export task.
// =====================================================================
var STATIONS_RAW = ee.FeatureCollection('users/YOUR_NAME/stations_for_gee');  // <-- EDIT
// Rebuild each point from its lon/lat columns, so the geometry is right even if the
// upload's X/Y setting was wrong. The 'lon' and 'lat' columns must be numbers.
function withPoint(f) {
  return ee.Feature(ee.Geometry.Point([ee.Number(f.get('lon')), ee.Number(f.get('lat'))]),
                    f.toDictionary());
}
var STATIONS = STATIONS_RAW.map(withPoint);
print('first station point', STATIONS.first());   // coordinates should be near [90.4, 23.7]
var START = '2023-12-01', END = '2024-03-01';   // dry season next to the 2024 Landsat mosaic

var S2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterDate(START, END)
  .filterBounds(STATIONS.geometry().bounds().buffer(5000))
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
  .map(function (img) {
    var scl = img.select('SCL');
    var clear = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
    return img.updateMask(clear);
  });
print('images in composite', S2.size());

var m = S2.median().select(['B3', 'B4', 'B8', 'B11']).divide(10000);
var green = m.select('B3'), red = m.select('B4'), nir = m.select('B8'), swir = m.select('B11');
var ndvi  = nir.subtract(red).divide(nir.add(red)).rename('ndvi');          // vegetation
var ndbi  = swir.subtract(nir).divide(swir.add(nir)).rename('ndbi');        // built-up / bare
var mndwi = green.subtract(swir).divide(green.add(swir)).rename('mndwi');   // water
var ndti  = red.subtract(green).divide(red.add(green)).rename('ndti');      // turbidity-ish
var stack = ee.Image.cat([ndvi, ndbi, mndwi, ndti]);

var radii = [100, 500, 2000];     // metres around each station
var bands = ['ndvi', 'ndbi', 'mndwi', 'ndti'];

var out = STATIONS.map(function (f) {
  var g = f.geometry();
  var o = ee.Feature(g).set({station: f.get('station'), river: f.get('river'),
                             p_hyp: f.get('p_hyp')});
  radii.forEach(function (r) {
    var d = stack.reduceRegion({reducer: ee.Reducer.mean(), geometry: g.buffer(r),
                                scale: 20, maxPixels: 1e8, tileScale: 4});
    bands.forEach(function (b) { o = o.set(b + '_r' + r, d.get(b)); });
  });
  return o;
});
print('first rows', out.limit(3));

var cols = ['station', 'river', 'p_hyp'];
radii.forEach(function (r) { bands.forEach(function (b) { cols.push(b + '_r' + r); }); });

Export.table.toDrive({collection: out, description: 's2_indices_per_station',
                      fileFormat: 'CSV', selectors: cols});
