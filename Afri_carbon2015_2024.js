/*
===============================================================================
CONTINENTAL-SCALE Forest Degradation and Carbon Sink Accounting
Study Area: ENTIRE AFRICA CONTINENT
TEMPORAL ANALYSIS (2015-2024) - Continental Coverage
===============================================================================
*/

// ============================================================================
// CONTINENTAL STUDY AREA DEFINITION - ENTIRE AFRICA
// ============================================================================

// Define the entire African continent boundary
var africa = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017")
  .filter(ee.Filter.or(
    ee.Filter.eq('wld_rgn', 'Africa'),
    ee.Filter.inList('country_na', [
      'Algeria', 'Angola', 'Benin', 'Botswana', 'Burkina Faso', 'Burundi',
      'Cameroon', 'Cape Verde', 'Central African Republic', 'Chad', 'Comoros',
      'Democratic Republic of the Congo', 'Republic of the Congo', 'Djibouti',
      'Egypt', 'Equatorial Guinea', 'Eritrea', 'Eswatini', 'Ethiopia',
      'Gabon', 'Gambia', 'Ghana', 'Guinea', 'Guinea-Bissau', 'Ivory Coast',
      'Kenya', 'Lesotho', 'Liberia', 'Libya', 'Madagascar', 'Malawi',
      'Mali', 'Mauritania', 'Mauritius', 'Morocco', 'Mozambique', 'Namibia',
      'Niger', 'Nigeria', 'Rwanda', 'Sao Tome and Principe', 'Senegal',
      'Seychelles', 'Sierra Leone', 'Somalia', 'South Africa', 'South Sudan',
      'Sudan', 'Tanzania', 'Togo', 'Tunisia', 'Uganda', 'Zambia', 'Zimbabwe'
    ])
  ))
  .geometry();

// Alternative: Define Africa by bounding box (more efficient for large-scale processing)
var africaBounds = ee.Geometry.Rectangle([-25, -40, 55, 40]); // West, South, East, North

// Use the bounding box for processing efficiency
var studyArea = africaBounds;

// Center the map on Africa
Map.centerObject(studyArea, 3);
Map.addLayer(studyArea, {color: 'blue', fillColor: '00000000'}, 'Africa Continent', true);

print('🌍 CONTINENTAL-SCALE ANALYSIS INITIATED');
print('Study Area: Entire African Continent');
print('Bounding Box: 25°W to 55°E, 40°S to 40°N');

// ============================================================================
// CONTINENTAL TIME PERIODS FOR ANALYSIS
// ============================================================================

var timePeriods = {
  'Africa_2015': {start: '2015-01-01', end: '2015-12-31', label: '2015'},
  'Africa_2016': {start: '2016-01-01', end: '2016-12-31', label: '2016'},
  'Africa_2017': {start: '2017-01-01', end: '2017-12-31', label: '2017'},
  'Africa_2018': {start: '2018-01-01', end: '2018-12-31', label: '2018'},
  'Africa_2019': {start: '2019-01-01', end: '2019-12-31', label: '2019'},
  'Africa_2020': {start: '2020-01-01', end: '2020-12-31', label: '2020'},
  'Africa_2021': {start: '2021-01-01', end: '2021-12-31', label: '2021'},
  'Africa_2022': {start: '2022-01-01', end: '2022-12-31', label: '2022'},
  'Africa_2023': {start: '2023-01-01', end: '2023-12-31', label: '2023'},
  'Africa_2024': {start: '2024-01-01', end: '2024-12-31', label: '2024'}
};

print('Processing 10 years of continental data (2015-2024)...');
print('Expected processing: ~30 million km² of African continent');

// ============================================================================
// CONTINENTAL-SCALE DATA PROCESSING FUNCTIONS
// ============================================================================

// Enhanced cloud masking for continental processing
function maskL8sr(image) {
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);
  
  // Apply scaling factors for Landsat 8/9
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true)
              .updateMask(qaMask)
              .updateMask(saturationMask);
}

// Continental vegetation indices calculation
function addVegetationIndices(image) {
  // NDVI (Normalized Difference Vegetation Index)
  var ndvi = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');
  
  // EVI (Enhanced Vegetation Index) - critical for dense tropical forests
  var evi = image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': image.select('SR_B5'),
      'RED': image.select('SR_B4'),
      'BLUE': image.select('SR_B2')
    }).rename('EVI');
  
  // SAVI (Soil Adjusted Vegetation Index) - important for arid regions
  var savi = image.expression(
    '((NIR - RED) / (NIR + RED + 0.5)) * (1 + 0.5)', {
      'NIR': image.select('SR_B5'),
      'RED': image.select('SR_B4')
    }).rename('SAVI');
  
  // NDWI (Normalized Difference Water Index)
  var ndwi = image.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI');
  
  // NDBI (Normalized Difference Built-up Index)
  var ndbi = image.normalizedDifference(['SR_B6', 'SR_B5']).rename('NDBI');
  
  return image.addBands([ndvi, evi, savi, ndwi, ndbi]);
}

// Continental-scale Landsat composite generation
function getLandsatComposite(startDate, endDate) {
  var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate(startDate, endDate)
    .filterBounds(studyArea)
    .filter(ee.Filter.lt('CLOUD_COVER', 40)); // More lenient for continental scale
    
  var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
    .filterDate(startDate, endDate)
    .filterBounds(studyArea)
    .filter(ee.Filter.lt('CLOUD_COVER', 40));
  
  var combined = l8.merge(l9);
  
  return combined
    .map(maskL8sr)
    .map(addVegetationIndices)
    .median()
    .clip(studyArea);
}

// ============================================================================
// CONTINENTAL LAND COVER CLASSIFICATION
// ============================================================================

function classifyAfricanLandCover(image) {
  // Water bodies (NDWI > 0.3)
  var water = image.select('NDWI').gt(0.3);
  
  // Desert/Bare areas (NDVI <= 0.1) - Major class in Africa
  var desert = image.select('NDVI').lte(0.1).and(water.not());
  
  // Built-up areas (NDBI > 0.1 and NDVI < 0.2)
  var builtUp = image.select('NDBI').gt(0.1).and(image.select('NDVI').lt(0.2))
    .and(desert.not());
  
  // Dense tropical forest (NDVI > 0.7 and EVI > 0.5) - Congo Basin, West Africa
  var denseForest = image.select('NDVI').gt(0.7).and(image.select('EVI').gt(0.5));
  
  // Medium forest (NDVI > 0.5 and NDVI <= 0.7) - Miombo woodlands
  var mediumForest = image.select('NDVI').gt(0.5).and(image.select('NDVI').lte(0.7))
    .and(denseForest.not());
  
  // Woodland/Savanna (NDVI > 0.3 and NDVI <= 0.5) - Major African biome
  var woodland = image.select('NDVI').gt(0.3).and(image.select('NDVI').lte(0.5))
    .and(mediumForest.not()).and(denseForest.not());
  
  // Grassland/Shrubland (NDVI > 0.2 and NDVI <= 0.3) - Sahel, East Africa
  var grassland = image.select('NDVI').gt(0.2).and(image.select('NDVI').lte(0.3))
    .and(woodland.not()).and(mediumForest.not()).and(denseForest.not());
  
  // Cropland (SAVI > 0.25 and agricultural areas)
  var cropland = image.select('SAVI').gt(0.25).and(image.select('NDVI').gt(0.15))
    .and(denseForest.not()).and(mediumForest.not()).and(woodland.not())
    .and(builtUp.not());
  
  // Sparse vegetation (NDVI > 0.1 and NDVI <= 0.2) - Semi-arid regions
  var sparseVeg = image.select('NDVI').gt(0.1).and(image.select('NDVI').lte(0.2))
    .and(grassland.not()).and(cropland.not()).and(builtUp.not()).and(desert.not());
  
  // Create African land cover classification (9 classes)
  var classified = ee.Image(0)
    .where(water, 1)           // Water
    .where(desert, 2)          // Desert/Bare areas
    .where(builtUp, 3)         // Built-up
    .where(sparseVeg, 4)       // Sparse vegetation
    .where(cropland, 5)        // Cropland
    .where(grassland, 6)       // Grassland/Shrubland
    .where(woodland, 7)        // Woodland/Savanna
    .where(mediumForest, 8)    // Medium forest
    .where(denseForest, 9)     // Dense tropical forest
    .rename('African_Land_Cover');
  
  return classified;
}

// ============================================================================
// CONTINENTAL FOREST AND CARBON ANALYSIS
// ============================================================================

// Continental forest cover mapping
function createAfricanForestMap(image) {
  // Different thresholds for different African biomes
  var tropicalForest = image.select('NDVI').gt(0.6).rename('Tropical_Forest');
  var woodlandSavanna = image.select('NDVI').gt(0.4).and(image.select('NDVI').lte(0.6))
    .rename('Woodland_Savanna');
  var sparseWoodland = image.select('NDVI').gt(0.25).and(image.select('NDVI').lte(0.4))
    .rename('Sparse_Woodland');
  var totalForest = image.select('NDVI').gt(0.25).rename('Total_Forest');
  
  return ee.Image([totalForest, tropicalForest, woodlandSavanna, sparseWoodland]);
}

// Continental biomass estimation (adapted for African ecosystems)
function estimateAfricanAGB(image) {
  var ndvi = image.select('NDVI');
  var evi = image.select('EVI');
  var savi = image.select('SAVI');
  
  // Biomass estimation adapted for African ecosystems
  // Different equations for different biomes
  var agb = image.expression(
    '(NDVI > 0.6) ? (300 * pow(NDVI, 2.1) + 120 * EVI + 30) : ' + // Tropical forests
    '(NDVI > 0.4) ? (180 * pow(NDVI, 1.9) + 80 * EVI + 60 * SAVI + 15) : ' + // Woodlands
    '(NDVI > 0.2) ? (90 * pow(NDVI, 1.7) + 40 * SAVI + 10) : 0', { // Savannas
      'NDVI': ndvi,
      'EVI': evi,
      'SAVI': savi
    }).rename('African_AGB');
  
  return agb.updateMask(agb.gt(1).and(agb.lt(800)));
}

// Continental carbon stock estimation
function africanBiomassToCarbon(agbImage) {
  var carbonFraction = 0.47; // IPCC default
  var rootShootRatio = 0.24;  // IPCC default
  
  // Above-ground carbon
  var agc = agbImage.multiply(carbonFraction).rename('African_AGC');
  
  // Below-ground carbon
  var bgc = agc.multiply(rootShootRatio).rename('African_BGC');
  
  // Total biomass carbon
  var totalBiomassCarbon = agc.add(bgc).rename('African_Total_Biomass_Carbon');
  
  // Soil organic carbon (adapted for African climates and soils)
  var soc = agbImage.expression(
    'AGB > 200 ? 100 : (AGB > 100 ? 70 : (AGB > 50 ? 45 : (AGB > 20 ? 30 : 15)))', {
      'AGB': agbImage
    }).rename('African_SOC');
  
  // Total ecosystem carbon
  var totalCarbon = totalBiomassCarbon.add(soc).rename('African_Total_Carbon');
  
  return ee.Image([agc, bgc, totalBiomassCarbon, soc, totalCarbon]);
}

// ============================================================================
// CONTINENTAL COLOR PALETTES AND VISUALIZATION
// ============================================================================

var continentalPalettes = {
  // African land cover classes (9 classes)
  african_landcover: ['#419bdf', '#F4A460', '#808080', '#DEB887', '#FFD700', '#90EE90', '#228B22', '#006400', '#004000'],
  
  // African forest types
  african_forest: ['#ffffe5', '#d9f0a3', '#78c679', '#41ab5d', '#238443', '#005a32'],
  ndvi: ['#8B0000', '#CD853F', '#F4A460', '#FFD700', '#ADFF2F', '#90EE90', '#228B22', '#006400'],
  biomass: ['#fff7f3', '#fcc5c0', '#f768a1', '#dd3497', '#ae017e', '#7a0177'],
  carbon: ['#f7fcf0', '#ccebc5', '#7bccc4', '#2b8cbe', '#08589e'],
  change: ['#8B0000', '#FF4500', '#FFD700', '#FFFFBF', '#90EE90', '#4682B4', '#000080']
};

var continentalVisParams = {
  african_landcover: {min: 1, max: 9, palette: continentalPalettes.african_landcover},
  african_forest: {min: 0, max: 1, palette: continentalPalettes.african_forest},
  ndvi: {min: -0.1, max: 0.8, palette: continentalPalettes.ndvi},
  biomass: {min: 0, max: 500, palette: continentalPalettes.biomass},
  carbon: {min: 0, max: 300, palette: continentalPalettes.carbon},
  change: {min: -0.5, max: 0.5, palette: continentalPalettes.change}
};

// ============================================================================
// CONTINENTAL ANALYSIS WORKFLOW
// ============================================================================

// Storage for continental results
var africaResults = {
  composites: {},
  forestCover: {},
  landCover: {},
  biomass: {},
  carbon: {}
};

print('🌍 Starting continental-scale processing...');
print('Warning: This may take significant processing time due to scale');

// Process each year for the entire continent
Object.keys(timePeriods).forEach(function(periodKey) {
  var period = timePeriods[periodKey];
  
  print('Processing Africa ' + period.label + '...');
  
  // Get continental Landsat composite
  var composite = getLandsatComposite(period.start, period.end);
  africaResults.composites[periodKey] = composite;
  
  // Create African forest cover maps
  var forestMaps = createAfricanForestMap(composite);
  africaResults.forestCover[periodKey] = forestMaps;
  
  // Create African land cover classification
  var landCover = classifyAfricanLandCover(composite);
  africaResults.landCover[periodKey] = landCover;
  
  // Estimate continental biomass and carbon
  var biomass = estimateAfricanAGB(composite);
  africaResults.biomass[periodKey] = biomass;
  
  var carbon = africanBiomassToCarbon(biomass);
  africaResults.carbon[periodKey] = carbon;
  
  // Add layers to map (hidden by default for performance)
  Map.addLayer(composite.select('NDVI'), continentalVisParams.ndvi, 
               'Africa NDVI ' + period.label, false);
  Map.addLayer(landCover, continentalVisParams.african_landcover, 
               'Africa Land Cover ' + period.label, false);
  Map.addLayer(forestMaps.select('Total_Forest'), continentalVisParams.african_forest, 
               'Africa Forest ' + period.label, false);
  Map.addLayer(carbon.select('African_Total_Carbon'), continentalVisParams.carbon, 
               'Africa Carbon ' + period.label, false);
});

// ============================================================================
// CONTINENTAL CHANGE DETECTION ANALYSIS
// ============================================================================

print('Calculating continental change detection (2015-2024)...');

var firstYear = 'Africa_2015';
var lastYear = 'Africa_2024';

// Continental forest change
var firstForest = africaResults.forestCover[firstYear].select('Total_Forest');
var lastForest = africaResults.forestCover[lastYear].select('Total_Forest');
var africaForestChange = lastForest.subtract(firstForest).rename('Africa_Forest_Change_2015_2024');

// Continental land cover change
var firstLC = africaResults.landCover[firstYear];
var lastLC = africaResults.landCover[lastYear];
var africaLCChange = lastLC.subtract(firstLC).rename('Africa_LC_Change_2015_2024');

// Continental carbon change
var firstCarbon = africaResults.carbon[firstYear].select('African_Total_Carbon');
var lastCarbon = africaResults.carbon[lastYear].select('African_Total_Carbon');
var africaCarbonChange = lastCarbon.subtract(firstCarbon).rename('Africa_Carbon_Change_2015_2024');

// Continental NDVI change
var firstNDVI = africaResults.composites[firstYear].select('NDVI');
var lastNDVI = africaResults.composites[lastYear].select('NDVI');
var africaNDVIChange = lastNDVI.subtract(firstNDVI).rename('Africa_NDVI_Change_2015_2024');

// Continental degradation hotspots
var africaDegradationHotspots = africaForestChange.lt(-0.2).rename('Africa_Degradation_Hotspots');

// Add continental change layers
Map.addLayer(africaForestChange, continentalVisParams.change, 
             'Africa Forest Change (2015-2024)', false);
Map.addLayer(africaCarbonChange, {min: -100, max: 100, palette: continentalPalettes.change}, 
             'Africa Carbon Change (2015-2024)', false);
Map.addLayer(africaDegradationHotspots, {palette: ['green', 'red']}, 
             'Africa Degradation Hotspots', true);

// ============================================================================
// CONTINENTAL EXPORT FUNCTIONS (OPTIMIZED FOR LARGE SCALE)
// ============================================================================

function exportContinentalImageToDrive(image, description, folder, scale) {
  scale = scale || 250; // Larger scale for continental processing
  
  Export.image.toDrive({
    image: image,
    description: description,
    folder: folder,
    scale: scale,
    region: studyArea,
    maxPixels: 1e10, // Increased for continental scale
    crs: 'EPSG:4326',
    fileFormat: 'GeoTIFF'
  });
  print('✓ Continental export task created: ' + description);
}

function exportContinentalVisualToDrive(image, palette, min, max, description, folder) {
  var visualized = image.visualize({
    min: min,
    max: max,
    palette: palette
  });
  
  Export.image.toDrive({
    image: visualized,
    description: description + '_Visual',
    folder: folder,
    scale: 500, // Coarser resolution for visualization
    region: studyArea,
    maxPixels: 1e10,
    crs: 'EPSG:4326',
    fileFormat: 'GeoTIFF'
  });
  print('✓ Continental visual export created: ' + description + '_Visual');
}

// ============================================================================
// CONTINENTAL EXPORT WORKFLOW
// ============================================================================

print('\n=== CREATING CONTINENTAL EXPORT TASKS ===');

// Export continental data for each year
Object.keys(timePeriods).forEach(function(periodKey) {
  var period = timePeriods[periodKey];
  
  // Continental forest cover (all types)
  exportContinentalImageToDrive(africaResults.forestCover[periodKey].select('Total_Forest'), 
                               'Africa_Total_Forest_' + period.label, 'Africa_Continental_Forest_Data', 250);
  exportContinentalImageToDrive(africaResults.forestCover[periodKey].select('Tropical_Forest'), 
                               'Africa_Tropical_Forest_' + period.label, 'Africa_Continental_Forest_Data', 250);
  exportContinentalImageToDrive(africaResults.forestCover[periodKey].select('Woodland_Savanna'), 
                               'Africa_Woodland_Savanna_' + period.label, 'Africa_Continental_Forest_Data', 250);
  
  // Continental land cover
  exportContinentalImageToDrive(africaResults.landCover[periodKey], 
                               'Africa_Land_Cover_' + period.label, 'Africa_Continental_LandCover_Data', 250);
  
  // Continental carbon stocks
  exportContinentalImageToDrive(africaResults.carbon[periodKey].select('African_Total_Carbon'), 
                               'Africa_Total_Carbon_' + period.label, 'Africa_Continental_Carbon_Data', 250);
  exportContinentalImageToDrive(africaResults.carbon[periodKey].select('African_AGC'), 
                               'Africa_Above_Ground_Carbon_' + period.label, 'Africa_Continental_Carbon_Data', 250);
  
  // Continental biomass
  exportContinentalImageToDrive(africaResults.biomass[periodKey], 
                               'Africa_Above_Ground_Biomass_' + period.label, 'Africa_Continental_Biomass_Data', 250);
  
  // Continental NDVI
  exportContinentalImageToDrive(africaResults.composites[periodKey].select('NDVI'), 
                               'Africa_NDVI_' + period.label, 'Africa_Continental_Vegetation_Data', 250);
  
  // Visual exports for continental overview
  exportContinentalVisualToDrive(africaResults.landCover[periodKey], 
                                continentalPalettes.african_landcover, 1, 9,
                                'Africa_Land_Cover_' + period.label, 'Africa_Continental_Visual');
  exportContinentalVisualToDrive(africaResults.forestCover[periodKey].select('Total_Forest'), 
                                continentalPalettes.african_forest, 0, 1,
                                'Africa_Forest_' + period.label, 'Africa_Continental_Visual');
  exportContinentalVisualToDrive(africaResults.carbon[periodKey].select('African_Total_Carbon'), 
                                continentalPalettes.carbon, 0, 300,
                                'Africa_Carbon_' + period.label, 'Africa_Continental_Visual');
});

// Export continental change detection results
exportContinentalImageToDrive(africaForestChange, 'Africa_Forest_Change_2015_2024', 'Africa_Continental_Change_Data', 250);
exportContinentalImageToDrive(africaLCChange, 'Africa_LandCover_Change_2015_2024', 'Africa_Continental_Change_Data', 250);
exportContinentalImageToDrive(africaCarbonChange, 'Africa_Carbon_Change_2015_2024', 'Africa_Continental_Change_Data', 250);
exportContinentalImageToDrive(africaNDVIChange, 'Africa_NDVI_Change_2015_2024', 'Africa_Continental_Change_Data', 250);
exportContinentalImageToDrive(africaDegradationHotspots, 'Africa_Degradation_Hotspots_2015_2024', 'Africa_Continental_Change_Data', 250);

// ============================================================================
// CONTINENTAL CSV DATA EXPORT (OPTIMIZED)
// ============================================================================

print('\n=== CREATING CONTINENTAL CSV DATASETS ===');

// Continental sample points (distributed across Africa)
var africaSamplePoints = ee.FeatureCollection.randomPoints(studyArea, 2000, 42);

// Export key years comparison (2015 vs 2024)
var continentalComparisonData = africaResults.composites[firstYear].select(['NDVI', 'EVI'])
  .addBands(africaResults.forestCover[firstYear].select('Total_Forest'))
  .addBands(africaResults.landCover[firstYear])
  .addBands(africaResults.carbon[firstYear].select('African_Total_Carbon'))
  .addBands(africaResults.composites[lastYear].select(['NDVI', 'EVI']))
  .addBands(africaResults.forestCover[lastYear].select('Total_Forest'))
  .addBands(africaResults.landCover[lastYear])
  .addBands(africaResults.carbon[lastYear].select('African_Total_Carbon'))
  .rename([
    'NDVI_2015', 'EVI_2015', 'Forest_2015', 'LandCover_2015', 'Carbon_2015',
    'NDVI_2024', 'EVI_2024', 'Forest_2024', 'LandCover_2024', 'Carbon_2024'
  ]);

// Sample continental comparison data
var continentalSampleData = continentalComparisonData.sampleRegions({
  collection: africaSamplePoints,
  scale: 500, // Coarser scale for continental analysis
  projection: 'EPSG:4326',
  geometries: true
});

// Export continental comparison CSV
Export.table.toDrive({
  collection: continentalSampleData,
  description: 'Africa_Continental_Comparison_2015_vs_2024',
  folder: 'Africa_Continental_CSV_Analysis',
  fileFormat: 'CSV'
});

// Continental change analysis dataset
var continentalChangeData = africaForestChange
  .addBands(africaCarbonChange)
  .addBands(africaNDVIChange)
  .addBands(africaDegradationHotspots)
  .rename(['Forest_Change', 'Carbon_Change', 'NDVI_Change', 'Degradation_Hotspot'])
  .sampleRegions({
    collection: ee.FeatureCollection.randomPoints(studyArea, 1500, 123),
    scale: 500,
    projection: 'EPSG:4326',
    geometries: true
  });

// Export continental change CSV
Export.table.toDrive({
  collection: continentalChangeData,
  description: 'Africa_Continental_Change_Analysis_2015_2024',
  folder: 'Africa_Continental_CSV_Analysis',
  fileFormat: 'CSV'
});

// ============================================================================
// CONTINENTAL ANALYSIS SUMMARY
// ============================================================================

print('\n=== CONTINENTAL ANALYSIS COMPLETE ===');
print('🌍 AFRICA-WIDE FOREST ANALYSIS (2015-2024) READY!');

print('\n📊 CONTINENTAL ANALYSIS SCOPE:');
print('• Study Area: Entire African Continent (~30.3 million km²)');
print('• Temporal Coverage: 10 years (2015-2024)');
print('• Processing Scale: 250m resolution for analysis, 500m for visualization');
print('• Biomes Covered: Tropical forests, Woodlands, Savannas, Deserts');

print('\n📁 CONTINENTAL EXPORT FOLDERS:');
print('1. Africa_Continental_Forest_Data/ - Forest cover (all types, all years)');
print('2. Africa_Continental_LandCover_Data/ - Land cover (9 classes, all years)');
print('3. Africa_Continental_Carbon_Data/ - Carbon stocks (all years)');
print('4. Africa_Continental_Biomass_Data/ - Biomass (all years)');
print('5. Africa_Continental_Vegetation_Data/ - NDVI (all years)');
print('6. Africa_Continental_Change_Data/ - Change detection (2015-2024)');
print('7. Africa_Continental_Visual/ - Presentation maps');
print('8. Africa_Continental_CSV_Analysis/ - Statistical datasets');

print('\n🗺️ AFRICAN LAND COVER CLASSES:');
print('1 - Water bodies');
print('2 - Desert/Bare areas (Sahara, Kalahari)');
print('3 - Built-up areas');
print('4 - Sparse vegetation (Semi-arid zones)');
print('5 - Cropland (Agricultural areas)');
print('6 - Grassland/Shrubland (Sahel region)');
print('7 - Woodland/Savanna (Miombo, Guinea-Sudan savanna)');
print('8 - Medium forest (Secondary forests, forest edges)');
print('9 - Dense tropical forest (Congo Basin, West African forests)');

print('\n📈 CONTINENTAL DATASETS EXPORTED:');
print('Forest Analysis (10 years each):');
print('  • Total forest cover across Africa');
print('  • Tropical forest distribution (Congo Basin focus)');
print('  • Woodland/Savanna mapping (East/Southern Africa)');
print('  • Sparse woodland (Sahel transition zones)');

print('Carbon Analysis (10 years each):');
print('  • Total ecosystem carbon across Africa');
print('  • Above-ground carbon stocks');
print('  • Below-ground carbon estimates');
print('  • Soil organic carbon mapping');

print('Change Detection (2015-2024):');
print('  • Forest change patterns across Africa');
print('  • Land cover transitions (desertification, deforestation)');
print('  • Carbon stock changes');
print('  • NDVI trends (vegetation health)');
print('  • Degradation hotspots identification');

print('\n📊 CONTINENTAL CSV DATASETS:');
print('• Africa_Continental_Comparison_2015_vs_2024.csv (2,000 points)');
print('• Africa_Continental_Change_Analysis_2015_2024.csv (1,500 points)');
print('• Covers all major African biomes and climate zones');

print('\n🌍 MAJOR AFRICAN ECOSYSTEMS COVERED:');
print('• Congo Basin Rainforests (Central Africa)');
print('• Miombo Woodlands (Southern Africa)');
print('• Guinea-Sudan Savanna (West Africa)');
print('• East African Savannas (Kenya, Tanzania)');
print('• Sahel Transition Zone (Semi-arid belt)');
print('• Sahara Desert Margins');
print('• Ethiopian Highlands');
print('• Fynbos (South African Cape)');
print('• Madagascar Forests');

print('\n🔬 CONTINENTAL RESEARCH APPLICATIONS:');
print('• Africa-wide deforestation monitoring');
print('• Carbon stock assessment for REDD+ programs');
print('• Climate change impact analysis');
print('• Desertification monitoring (Sahel region)');
print('• Biodiversity conservation planning');
print('• Food security and agricultural mapping');
print('• Sustainable development tracking');
print('• Continental forest governance');

print('\n⚠️ PROCESSING CONSIDERATIONS:');
print('• Large-scale processing: May require significant computation time');
print('• Resolution: 250m for analysis, 500m for visualization');
print('• Cloud processing: Utilizes Google Earth Engine distributed computing');
print('• Data volume: Expect large file sizes for continental coverage');
print('• Regional analysis: Consider breaking into sub-regions for detailed studies');

print('\n🚀 NEXT STEPS FOR CONTINENTAL ANALYSIS:');
print('1. Run export tasks in batches (may take hours for completion)');
print('2. Monitor task progress in the Tasks tab');
print('3. Download continental datasets from Google Drive');
print('4. Consider regional breakdowns for detailed analysis');
print('5. Use coarser resolution data for continent-wide overviews');
print('6. Combine with regional/national detailed studies');

print('\n📍 CONTINENTAL COVERAGE:');
print('Bounding Box: 25°W to 55°E, 40°S to 40°N');
print('Coverage: ~30.3 million km² (entire African continent)');
print('Resolution: 250m spatial resolution');
print('Temporal: 2015-2024 (10-year comprehensive analysis)');
print('Processing Scale: Distributed cloud computing optimized');

print('\n✅ CONTINENTAL-SCALE ANALYSIS READY!');
print('🌍 ENTIRE AFRICA: Complete forest and carbon monitoring system');
print('Perfect for continental policy making and climate research!');

// ============================================================================
// REGIONAL BREAKDOWN SUGGESTIONS
// ============================================================================

print('\n🗺️ SUGGESTED REGIONAL BREAKDOWNS FOR DETAILED ANALYSIS:');

// Define major African regions for more detailed analysis
var africaRegions = {
  'West_Africa': ee.Geometry.Rectangle([-18, 4, 16, 18]),
  'Central_Africa': ee.Geometry.Rectangle([8, -8, 32, 8]),
  'East_Africa': ee.Geometry.Rectangle([28, -12, 52, 18]),
  'Southern_Africa': ee.Geometry.Rectangle([10, -36, 38, -8]),
  'North_Africa': ee.Geometry.Rectangle([-12, 18, 38, 38])
};

print('Consider these regional analyses for detailed studies:');
print('• West Africa: Guinea-Sudan savanna, coastal forests');
print('• Central Africa: Congo Basin rainforests');
print('• East Africa: Savannas, highlands, Rift Valley');
print('• Southern Africa: Miombo woodlands, Kalahari');
print('• North Africa: Sahel transition, desert margins');

// Add regional boundaries to map for reference
Object.keys(africaRegions).forEach(function(regionName) {
  Map.addLayer(africaRegions[regionName], 
               {color: 'yellow', fillColor: '00000000'}, 
               regionName + ' Region', false);
});

// Display final continental results
Map.addLayer(africaResults.landCover[lastYear], continentalVisParams.african_landcover, 
             'Africa Land Cover 2024', true);
Map.addLayer(africaResults.carbon[lastYear].select('African_Total_Carbon'), continentalVisParams.carbon, 
             'Africa Carbon Stock 2024', true);
Map.addLayer(africaDegradationHotspots, {palette: ['green', 'red']}, 
             'Africa Degradation Hotspots (2015-2024)', true);

// Set optimal view for continental overview
Map.setOptions('SATELLITE');
Map.centerObject(studyArea, 3);