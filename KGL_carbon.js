/*
===============================================================================
Forest Degradation and Carbon Sink Accounting Using Google Earth Engine
Case Study: Mount Kigali Forest, Nyarugenege District, Rwanda
===============================================================================
*/

// ============================================================================
// STUDY AREA DEFINITION - MOUNT KIGALI FOREST
// ============================================================================

// Define Mount Kigali Forest boundary (Nyarugenege District, Rwanda)
var mountKigali = ee.Geometry.Rectangle([
  30.0, -1.95,   // Southwest corner
  30.15, -1.85   // Northeast corner
]);

// Center the map on Mount Kigali
Map.centerObject(mountKigali, 12);
Map.addLayer(mountKigali, {color: 'red'}, 'Study Area - Mount Kigali Forest', false);

// ============================================================================
// TIME PERIODS FOR ANALYSIS
// ============================================================================

var timePeriods = {
  'Period_1_2015_2016': {start: '2015-01-01', end: '2016-12-31', label: '2015-2016'},
  'Period_2_2017_2018': {start: '2017-01-01', end: '2018-12-31', label: '2017-2018'},
  'Period_3_2019_2020': {start: '2019-01-01', end: '2020-12-31', label: '2019-2020'},
  'Period_4_2021_2022': {start: '2021-01-01', end: '2022-12-31', label: '2021-2022'},
  'Period_5_2023_2024': {start: '2023-01-01', end: '2024-12-31', label: '2023-2024'}
};

// ============================================================================
// LANDSAT DATA PROCESSING FUNCTIONS
// ============================================================================

// Function to mask clouds and apply scaling factors
function maskL8sr(image) {
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);
  
  // Apply scaling factors
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true)
              .updateMask(qaMask)
              .updateMask(saturationMask);
}

// Function to calculate vegetation indices
function addVegetationIndices(image) {
  // NDVI
  var ndvi = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');
  
  // EVI
  var evi = image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': image.select('SR_B5'),
      'RED': image.select('SR_B4'),
      'BLUE': image.select('SR_B2')
    }).rename('EVI');
  
  // SAVI
  var savi = image.expression(
    '((NIR - RED) / (NIR + RED + 0.5)) * (1 + 0.5)', {
      'NIR': image.select('SR_B5'),
      'RED': image.select('SR_B4')
    }).rename('SAVI');
  
  // MSAVI
  var msavi = image.expression(
    '(2 * NIR + 1 - sqrt(pow((2 * NIR + 1), 2) - 8 * (NIR - RED))) / 2', {
      'NIR': image.select('SR_B5'),
      'RED': image.select('SR_B4')
    }).rename('MSAVI');
  
  return image.addBands([ndvi, evi, savi, msavi]);
}

// Function to get Landsat composite for a specific time period
function getLandsatComposite(startDate, endDate) {
  return ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate(startDate, endDate)
    .filterBounds(mountKigali)
    .filter(ee.Filter.lt('CLOUD_COVER', 20))
    .map(maskL8sr)
    .map(addVegetationIndices)
    .median()
    .clip(mountKigali);
}

// ============================================================================
// FOREST COVER MAPPING
// ============================================================================

// Function to create forest cover map
function createForestCoverMap(image, threshold) {
  threshold = threshold || 0.4;
  var forestCover = image.select('NDVI').gt(threshold).rename('Forest_Cover');
  return forestCover;
}

// Function to estimate above-ground biomass
function estimateAGB(image) {
  var ndvi = image.select('NDVI');
  var evi = image.select('EVI');
  
  // Biomass estimation using multiple regression approach
  // Based on tropical forest allometric equations
  var agb = image.expression(
    '(NDVI > 0.2) ? (150 * pow(NDVI, 2) + 50 * EVI + 25) : 0', {
      'NDVI': ndvi,
      'EVI': evi
    }).rename('AGB');
  
  return agb.updateMask(agb.gt(0));
}

// Function to convert biomass to carbon
function biomassToCarbon(agbImage) {
  var carbonFraction = 0.47; // IPCC default for tropical forests
  var rootShootRatio = 0.24; // IPCC default for tropical forests
  
  var agc = agbImage.multiply(carbonFraction).rename('AGC');
  var bgc = agc.multiply(rootShootRatio).rename('BGC');
  var totalCarbon = agc.add(bgc).rename('Total_Carbon');
  
  return ee.Image([agc, bgc, totalCarbon]);
}

// ============================================================================
// COLOR PALETTES AND VISUALIZATION PARAMETERS
// ============================================================================

// Beautiful color palettes for different map types
var palettes = {
  forest: ['#ffffe5', '#d9f0a3', '#addd8e', '#78c679', '#41ab5d', '#238443', '#005a32'],
  ndvi: ['#d73027', '#f46d43', '#fdae61', '#fee08b', '#d9ef8b', '#a6d96a', '#66bd63', '#1a9850'],
  biomass: ['#fff7f3', '#fde0dd', '#fcc5c0', '#fa9fb5', '#f768a1', '#dd3497', '#ae017e', '#7a0177'],
  carbon: ['#f7fcf0', '#e0f3db', '#ccebc5', '#a8ddb5', '#7bccc4', '#4eb3d3', '#2b8cbe', '#08589e'],
  change: ['#d73027', '#f46d43', '#fdae61', '#fee090', '#ffffbf', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4'],
  degradation: ['#006837', '#31a354', '#78c679', '#c2e699', '#ffffcc', '#fec44f', '#fe9929', '#ec7014', '#cc4c02']
};

// Visualization parameters for different map types
var visParams = {
  forest: {min: 0, max: 1, palette: palettes.forest},
  ndvi: {min: -0.2, max: 0.8, palette: palettes.ndvi},
  biomass: {min: 0, max: 300, palette: palettes.biomass},
  carbon: {min: 0, max: 150, palette: palettes.carbon},
  change: {min: -0.5, max: 0.5, palette: palettes.change},
  degradation: {min: 0, max: 1, palette: palettes.degradation}
};

// ============================================================================
// MAIN ANALYSIS WORKFLOW
// ============================================================================

// Store results for each time period
var forestCoverResults = {};
var carbonStockResults = {};
var compositeResults = {};

print('Starting Forest Degradation and Carbon Analysis for Mount Kigali...');

// Process each time period
Object.keys(timePeriods).forEach(function(periodKey) {
  var period = timePeriods[periodKey];
  
  print('Processing ' + period.label + '...');
  
  // Get Landsat composite for the period
  var composite = getLandsatComposite(period.start, period.end);
  compositeResults[periodKey] = composite;
  
  // Create forest cover map
  var forestCover = createForestCoverMap(composite);
  forestCoverResults[periodKey] = forestCover;
  
  // Estimate biomass and carbon
  var agb = estimateAGB(composite);
  var carbon = biomassToCarbon(agb);
  carbonStockResults[periodKey] = carbon;
  
  // Add layers to map for visualization
  Map.addLayer(composite, {bands: ['SR_B4', 'SR_B3', 'SR_B2'], min: 0, max: 0.3}, 
               'RGB ' + period.label, false);
  Map.addLayer(composite.select('NDVI'), visParams.ndvi, 
               'NDVI ' + period.label, false);
  Map.addLayer(forestCover, visParams.forest, 
               'Forest Cover ' + period.label, false);
  Map.addLayer(carbon.select('Total_Carbon'), visParams.carbon, 
               'Carbon Stock ' + period.label, false);
});

// ============================================================================
// CHANGE DETECTION ANALYSIS
// ============================================================================

print('Calculating change detection...');

// Get first and last periods for change analysis
var firstPeriodKey = Object.keys(timePeriods)[0];
var lastPeriodKey = Object.keys(timePeriods)[Object.keys(timePeriods).length - 1];

var firstForest = forestCoverResults[firstPeriodKey];
var lastForest = forestCoverResults[lastPeriodKey];
var firstCarbon = carbonStockResults[firstPeriodKey].select('Total_Carbon');
var lastCarbon = carbonStockResults[lastPeriodKey].select('Total_Carbon');

// Calculate annual change rates
var yearDiff = 9; // 2015 to 2024 = 9 years
var forestChangeRate = lastForest.subtract(firstForest).divide(yearDiff).rename('Forest_Change_Rate');
var carbonChangeRate = lastCarbon.subtract(firstCarbon).divide(yearDiff).rename('Carbon_Change_Rate');

// Identify degradation hotspots (areas with significant forest loss)
var degradationHotspots = forestChangeRate.lt(-0.1).rename('Degradation_Hotspots');

// Add change layers to map
Map.addLayer(forestChangeRate, visParams.change, 'Forest Change Rate (%/year)', false);
Map.addLayer(carbonChangeRate, {min: -10, max: 10, palette: palettes.change}, 
             'Carbon Change Rate (Mg C/ha/year)', false);
Map.addLayer(degradationHotspots, {min: 0, max: 1, palette: ['green', 'red']}, 
             'Degradation Hotspots', false);

// ============================================================================
// FRAGMENTATION ANALYSIS
// ============================================================================

print('Calculating fragmentation indices...');

// Function to calculate fragmentation metrics
function calculateFragmentation(forestImage) {
  // Convert to binary
  var forest = forestImage.gt(0);
  
  // Calculate patch density using focal statistics
  var patchDensity = forest.focal_mean(1, 'square', 'pixels', 3).rename('Patch_Density');
  
  // Calculate edge density
  var edges = forest.subtract(forest.focal_min(1, 'square', 'pixels', 3)).rename('Edges');
  var edgeDensity = edges.focal_sum(1, 'square', 'pixels', 5).rename('Edge_Density');
  
  // Calculate connectivity
  var connectivity = forest.focal_mean(1, 'square', 'pixels', 7).rename('Connectivity');
  
  return ee.Image([patchDensity, edgeDensity, connectivity]);
}

// Calculate fragmentation for the latest period
var latestFragmentation = calculateFragmentation(forestCoverResults[lastPeriodKey]);

Map.addLayer(latestFragmentation.select('Patch_Density'), visParams.degradation, 
             'Patch Density', false);
Map.addLayer(latestFragmentation.select('Connectivity'), visParams.degradation, 
             'Forest Connectivity', false);

// ============================================================================
// BEAUTIFUL MAP CREATION AND EXPORT FUNCTIONS
// ============================================================================

// Function to create a beautiful styled map
function createStyledMap(image, title, palette, min, max, units) {
  units = units || '';
  
  // Create a panel for the map
  var panel = ui.Panel({
    style: {
      width: '100%',
      height: '100%',
      padding: '0px'
    }
  });
  
  // Add title
  var titleLabel = ui.Label({
    value: title,
    style: {
      fontSize: '24px',
      fontWeight: 'bold',
      textAlign: 'center',
      margin: '10px',
      color: '#2c3e50'
    }
  });
  
  // Add subtitle with study area info
  var subtitleLabel = ui.Label({
    value: 'Mount Kigali Forest, Nyarugenege District, Rwanda',
    style: {
      fontSize: '16px',
      textAlign: 'center',
      margin: '5px',
      color: '#7f8c8d'
    }
  });
  
  return {title: titleLabel, subtitle: subtitleLabel, image: image, 
          palette: palette, min: min, max: max, units: units};
}

// Function to export images to Google Drive
function exportImageToDrive(image, description, folder) {
  folder = folder || 'GEE_Mount_Kigali_Analysis';
  
  Export.image.toDrive({
    image: image,
    description: description,
    folder: folder,
    scale: 30,
    region: mountKigali,
    maxPixels: 1e9,
    crs: 'EPSG:4326',
    fileFormat: 'GeoTIFF'
  });
  
  print('Export task created: ' + description);
}

// Function to export maps as high-resolution images
function exportMapToDrive(image, title, palette, min, max, description, folder) {
  folder = folder || 'GEE_Mount_Kigali_Maps';
  
  // Create visualization
  var visualized = image.visualize({
    min: min,
    max: max,
    palette: palette
  });
  
  Export.image.toDrive({
    image: visualized,
    description: description + '_Map',
    folder: folder,
    scale: 30,
    region: mountKigali,
    dimensions: '2048x2048',
    maxPixels: 1e9,
    crs: 'EPSG:4326',
    fileFormat: 'PNG'
  });
  
  print('Map export task created: ' + description + '_Map');
}

// ============================================================================
// EXPORT ALL RESULTS
// ============================================================================

print('Creating export tasks...');

// Export forest cover maps for each period
Object.keys(timePeriods).forEach(function(periodKey) {
  var period = timePeriods[periodKey];
  var forestCover = forestCoverResults[periodKey];
  var carbon = carbonStockResults[periodKey];
  var composite = compositeResults[periodKey];
  
  // Export GeoTIFF files
  exportImageToDrive(forestCover, 'Forest_Cover_' + periodKey, 'Mount_Kigali_TIF');
  exportImageToDrive(carbon.select('Total_Carbon'), 'Carbon_Stock_' + periodKey, 'Mount_Kigali_TIF');
  exportImageToDrive(composite.select('NDVI'), 'NDVI_' + periodKey, 'Mount_Kigali_TIF');
  exportImageToDrive(composite.select(['SR_B4', 'SR_B3', 'SR_B2']), 'RGB_' + periodKey, 'Mount_Kigali_TIF');
  
  // Export beautiful maps
  exportMapToDrive(forestCover, 'Forest Cover ' + period.label, palettes.forest, 0, 1, 
                   'Forest_Cover_Map_' + periodKey, 'Mount_Kigali_Maps');
  exportMapToDrive(carbon.select('Total_Carbon'), 'Carbon Stock ' + period.label + ' (Mg C/ha)', 
                   palettes.carbon, 0, 150, 'Carbon_Stock_Map_' + periodKey, 'Mount_Kigali_Maps');
  exportMapToDrive(composite.select('NDVI'), 'NDVI ' + period.label, palettes.ndvi, -0.2, 0.8, 
                   'NDVI_Map_' + periodKey, 'Mount_Kigali_Maps');
});

// Export change detection results
exportImageToDrive(forestChangeRate, 'Forest_Change_Rate', 'Mount_Kigali_TIF');
exportImageToDrive(carbonChangeRate, 'Carbon_Change_Rate', 'Mount_Kigali_TIF');
exportImageToDrive(degradationHotspots, 'Degradation_Hotspots', 'Mount_Kigali_TIF');

// Export change maps
exportMapToDrive(forestChangeRate, 'Annual Forest Change Rate (%/year)', palettes.change, -0.5, 0.5, 
                 'Forest_Change_Rate_Map', 'Mount_Kigali_Maps');
exportMapToDrive(carbonChangeRate, 'Annual Carbon Change Rate (Mg C/ha/year)', palettes.change, -10, 10, 
                 'Carbon_Change_Rate_Map', 'Mount_Kigali_Maps');
exportMapToDrive(degradationHotspots, 'Forest Degradation Hotspots', ['green', 'red'], 0, 1, 
                 'Degradation_Hotspots_Map', 'Mount_Kigali_Maps');

// Export fragmentation analysis
exportImageToDrive(latestFragmentation, 'Fragmentation_Indices_2023_2024', 'Mount_Kigali_TIF');
exportMapToDrive(latestFragmentation.select('Connectivity'), 'Forest Connectivity Index', 
                 palettes.degradation, 0, 1, 'Forest_Connectivity_Map', 'Mount_Kigali_Maps');

// ============================================================================
// COMPREHENSIVE STATISTICS CALCULATION
// ============================================================================

print('Calculating comprehensive statistics...');

// Function to calculate area statistics
function calculateAreaStats(image, scale, region) {
  scale = scale || 30;
  region = region || mountKigali;
  
  var areaImage = ee.Image.pixelArea().divide(10000); // Convert to hectares
  var stats = areaImage.updateMask(image.gt(0)).reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
  });
  
  return stats;
}

// Calculate forest area for each period
print('\n=== FOREST AREA STATISTICS ===');
Object.keys(timePeriods).forEach(function(periodKey) {
  var period = timePeriods[periodKey];
  var forestArea = calculateAreaStats(forestCoverResults[periodKey]);
  print(period.label + ' Forest Area: ' + forestArea.get('area').getInfo().toFixed(2) + ' hectares');
});

// Calculate total carbon stock for each period
print('\n=== CARBON STOCK STATISTICS ===');
Object.keys(timePeriods).forEach(function(periodKey) {
  var period = timePeriods[periodKey];
  var carbonStats = carbonStockResults[periodKey].select('Total_Carbon').reduceRegion({
    reducer: ee.Reducer.mean().combine({
      reducer2: ee.Reducer.sum(),
      sharedInputs: true
    }),
    geometry: mountKigali,
    scale: 30,
    maxPixels: 1e9
  });
  
  var meanCarbon = carbonStats.get('Total_Carbon_mean');
  var totalCarbon = carbonStats.get('Total_Carbon_sum');
  
  print(period.label + ' - Mean Carbon: ' + ee.Number(meanCarbon).format('%.2f').getInfo() + 
        ' Mg C/ha, Total Carbon: ' + ee.Number(totalCarbon).divide(1000).format('%.2f').getInfo() + ' Gg C');
});

// Calculate change statistics
print('\n=== CHANGE DETECTION STATISTICS ===');
var forestChangeStats = forestChangeRate.reduceRegion({
  reducer: ee.Reducer.mean().combine({
    reducer2: ee.Reducer.minMax(),
    sharedInputs: true
  }),
  geometry: mountKigali,
  scale: 30,
  maxPixels: 1e9
});

print('Annual Forest Change Rate:');
print('  Mean: ' + ee.Number(forestChangeStats.get('Forest_Change_Rate_mean')).format('%.4f').getInfo() + ' %/year');
print('  Min: ' + ee.Number(forestChangeStats.get('Forest_Change_Rate_min')).format('%.4f').getInfo() + ' %/year');
print('  Max: ' + ee.Number(forestChangeStats.get('Forest_Change_Rate_max')).format('%.4f').getInfo() + ' %/year');

var carbonChangeStats = carbonChangeRate.reduceRegion({
  reducer: ee.Reducer.mean().combine({
    reducer2: ee.Reducer.minMax(),
    sharedInputs: true
  }),
  geometry: mountKigali,
  scale: 30,
  maxPixels: 1e9
});

print('Annual Carbon Change Rate:');
print('  Mean: ' + ee.Number(carbonChangeStats.get('Carbon_Change_Rate_mean')).format('%.2f').getInfo() + ' Mg C/ha/year');
print('  Min: ' + ee.Number(carbonChangeStats.get('Carbon_Change_Rate_min')).format('%.2f').getInfo() + ' Mg C/ha/year');
print('  Max: ' + ee.Number(carbonChangeStats.get('Carbon_Change_Rate_max')).format('%.2f').getInfo() + ' Mg C/ha/year');

// Calculate degradation hotspot area
var hotspotArea = calculateAreaStats(degradationHotspots);
print('Degradation Hotspot Area: ' + hotspotArea.get('area').getInfo().toFixed(2) + ' hectares');

// ============================================================================
// CORRELATION ANALYSIS
// ============================================================================

print('\n=== CORRELATION ANALYSIS ===');

// Sample points for correlation analysis
var samplePoints = ee.FeatureCollection.randomPoints({
  region: mountKigali,
  points: 1000,
  seed: 42
});

// Extract values for correlation analysis (using latest period)
var latestComposite = compositeResults[lastPeriodKey];
var latestForest = forestCoverResults[lastPeriodKey];
var latestCarbon = carbonStockResults[lastPeriodKey].select('Total_Carbon');
var latestFragmentation = calculateFragmentation(latestForest);

var correlationData = latestComposite.select(['NDVI', 'EVI'])
  .addBands(latestForest)
  .addBands(latestCarbon)
  .addBands(latestFragmentation)
  .sampleRegions({
    collection: samplePoints,
    scale: 30,
    projection: 'EPSG:4326'
  });

print('Correlation analysis data sampled at 1000 random points');
print('Use the sampled data for statistical analysis in your preferred software');

// ============================================================================
// FINAL INSTRUCTIONS
// ============================================================================

print('\n=== ANALYSIS COMPLETE ===');
print('🎉 Forest Degradation and Carbon Analysis for Mount Kigali Forest completed!');
print('\n📊 RESULTS SUMMARY:');
print('• Processed 5 time periods (2015-2024)');
print('• Generated forest cover maps, carbon stock maps, and change detection');
print('• Calculated fragmentation indices and degradation hotspots');
print('• Created comprehensive statistics and correlation data');

print('\n📁 EXPORT TASKS CREATED:');
print('Check the Tasks tab to run all export jobs');
print('Files will be saved to your Google Drive in:');
print('  • Mount_Kigali_TIF/ - GeoTIFF files for analysis');
print('  • Mount_Kigali_Maps/ - Beautiful PNG maps for presentation');

print('\n🗺️ MAP LAYERS:');
print('Toggle layer visibility in the Layers panel to explore results');
print('All layers are loaded but hidden by default for performance');

print('\n📈 NEXT STEPS:');
print('1. Run export tasks in the Tasks tab');
print('2. Download files from Google Drive');
print('3. Use statistics for your research paper');
print('4. Create presentations with the beautiful maps');

print('\n🔬 STUDY AREA INFO:');
print('Location: Mount Kigali Forest, Nyarugenege District, Rwanda');
print('Coordinates: 30.0°E to 30.15°E, -1.95°S to -1.85°S');
print('Analysis Resolution: 30m (Landsat 8)');
print('Temporal Coverage: 2015-2024 (9 years)');

// Display the study area boundary prominently
Map.addLayer(mountKigali, {color: 'red', fillColor: '00000000'}, 'Mount Kigali Forest Boundary', true);

// Set a good default view
Map.centerObject(mountKigali, 13);
Map.setOptions('SATELLITE');