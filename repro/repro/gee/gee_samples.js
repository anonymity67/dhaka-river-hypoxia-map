// =====================================================================
// Sentinel-2 water indices for every water sample (Earth Engine Code Editor)
// 1. Upload samples_for_gee.csv as a table asset (Assets > NEW > Table upload,
//    X column = lon, Y column = lat). Copy the asset ID below.
// 2. Paste this whole script in code.earthengine.google.com and press Run.
// 3. In the Tasks tab, click RUN next to the export, then pick Drive.
//
// Source rule: surface-reflectance images (SR) are used when any exist within the
// window; otherwise top-of-atmosphere images (TOA) are used. The column 'source'
// says which. Early samples (2017 to 2018) will often be TOA.
// =====================================================================
var SAMPLES_RAW = ee.FeatureCollection('users/YOUR_NAME/samples_for_gee');  // <-- EDIT
// Rebuild each point from its lon/lat columns, so the geometry is right even if the
// upload's X/Y setting was wrong. The 'lon' and 'lat' columns must be numbers.
function withPoint(f) {
  return ee.Feature(ee.Geometry.Point([ee.Number(f.get('lon')), ee.Number(f.get('lat'))]),
                    f.toDictionary());
}
var SAMPLES = SAMPLES_RAW.map(withPoint);
print('first sample point', SAMPLES.first());   // coordinates should be near [90.4, 23.7]
var WINDOW_DAYS = 5;      // use images within +/- this many days of the sample
var BUFFER_M    = 100;    // radius around the station, in metres

var S2_SR  = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var S2_TOA = ee.ImageCollection('COPERNICUS/S2_HARMONIZED');

// Same indices for both sources; 'clear' = pixel not cloud or shadow
function indices(img, clear) {
  var s     = img.select(['B3', 'B4', 'B5', 'B11']).divide(10000);
  var green = s.select('B3'), red = s.select('B4'), re1 = s.select('B5'), swir = s.select('B11');
  var mndwi = green.subtract(swir).divide(green.add(swir)).rename('mndwi');   // water > 0
  var ndci  = re1.subtract(red).divide(re1.add(red)).rename('ndci');          // chlorophyll proxy
  var ndti  = red.subtract(green).divide(red.add(green)).rename('ndti');      // turbidity proxy
  var water = mndwi.gt(0).and(clear);
  return ee.Image.cat([red.rename('red'), green.rename('green'), mndwi, ndci, ndti])
           .toFloat().updateMask(water);
}
function prepSR(img) {
  var scl = img.select('SCL');
  var clear = scl.neq(0).and(scl.neq(1)).and(scl.neq(3))
                 .and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  return indices(img, clear);
}
function prepTOA(img) {
  var qa = img.select('QA60');   // cloud flags (valid for 2017 to 2021)
  var clear = qa.bitwiseAnd(1 << 10).eq(0).and(qa.bitwiseAnd(1 << 11).eq(0));
  return indices(img, clear);
}

// An all-masked image, so the median always has the same bands even when no
// image exists within the window (those samples then get empty values).
var BLANK = ee.Image.constant([0, 0, 0, 0, 0]).toFloat()
  .rename(['red', 'green', 'mndwi', 'ndci', 'ndti'])
  .updateMask(ee.Image.constant(0));

var out = SAMPLES.map(function (f) {
  var date  = ee.Date(f.get('date'));
  var start = date.advance(-WINDOW_DAYS, 'day');
  var end   = date.advance(WINDOW_DAYS + 1, 'day');
  var geom  = f.geometry().buffer(BUFFER_M);
  var colSR  = S2_SR.filterBounds(f.geometry()).filterDate(start, end);
  var colTOA = S2_TOA.filterBounds(f.geometry()).filterDate(start, end);
  var useSR  = colSR.size().gt(0);
  var col    = ee.ImageCollection(ee.Algorithms.If(useSR, colSR.map(prepSR), colTOA.map(prepTOA)));
  var med    = col.merge(ee.ImageCollection([BLANK])).median();
  var means  = med.reduceRegion({reducer: ee.Reducer.mean(), geometry: geom,
                                 scale: 10, maxPixels: 1e8, tileScale: 4});
  var nwat   = med.select('mndwi').reduceRegion({reducer: ee.Reducer.count(), geometry: geom,
                                 scale: 10, maxPixels: 1e8, tileScale: 4});
  return f.set(means).set({
    n_img: col.size(),
    source: ee.Algorithms.If(useSR, 'SR', 'TOA'),
    n_water_px: nwat.get('mndwi')
  });
});

print('samples', SAMPLES.size());
print('first rows', out.limit(5));

Export.table.toDrive({
  collection: out,
  description: 's2_indices_per_sample',
  fileFormat: 'CSV',
  selectors: ['sample_id', 'station', 'river', 'date', 'season', 'do', 'hyp',
              'source', 'n_img', 'n_water_px', 'red', 'green', 'mndwi', 'ndci', 'ndti']
});
