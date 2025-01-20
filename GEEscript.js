var ROI = ROI6

/* ---- I. Functions ---- */
print( 'SECTION 1 - FUNCTIONS')
// SOME FUNCTIONS ARE NOT CURRENTLY USED

function maskLandsatSR(image) {
  var cloudsBitMask = (1 << 3);
  var cloudShadowBitMask = (1 << 4);
  var snowBitMask = (1 << 5);
  var waterBitMask = (1 << 7);
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(cloudsBitMask).eq(0)
                .and(qa.bitwiseAnd(cloudShadowBitMask).eq(0))
                .and(qa.bitwiseAnd(snowBitMask).eq(0))
                .and(qa.bitwiseAnd(waterBitMask).eq(0))
                //.and(dem.gt(demthr))
  return image.updateMask(mask);
}
 
function applyScaleFactors(image) {
  var opticalBands = image.select(['NIR','RED']).multiply(0.0000275).add(-0.2);
  return image.addBands(opticalBands, null, true);
}

function addNDVI(image){
  var ndvi= image.normalizedDifference(['NIR','RED']).multiply(10000).toInt().rename('NDVI');
  return image.addBands(ndvi);  
}

function addNIRv(image){
  var nirv = image.expression('NIR*((NIR - RED)/(NIR + RED))', {
    RED: image.select('RED'),
    NIR: image.select('NIR')
  })
  return image.addBands(nirv.multiply(10000).toInt().rename('NDVI'))
}

function addKNDVI_tanh(image) {
  // Compute (nir-red)^2
  var red = image.select('RED')
  var nir = image.select('NIR')
  var D2 = nir.subtract(red).pow(2);
  var sigma = nir.add(red).multiply(0.5);
  //var sigma = ee.Number(0.15)
  var kndvi = D2.divide(sigma.multiply(2.0).pow(2)).tanh();
  return image.addBands(kndvi.select([0], ['NDVI']).multiply(10000).toInt().rename('NDVI'));
}

function reclassify(image){
  return image.where(image.gt(-10000).and(image.lte(10000)), 1);
}

function brdfCorrect(image) {
      var inputBandNames = image.bandNames() 
      var constants = {
        pi: Math.PI
      }

      var coefficientsByBand = {
        //'blue': {fiso: 0.0774, fgeo: 0.0079, fvol: 0.0372},
        //'green': {fiso: 0.1306, fgeo: 0.0178, fvol: 0.0580},
        'RED': {fiso: 0.1690, fgeo: 0.0227, fvol: 0.0574},
        'NIR': {fiso: 0.3093, fgeo: 0.0330, fvol: 0.1535},
        //'swir1': {fiso: 0.3430, fgeo: 0.0453, fvol: 0.1154},
        //'swir2': {fiso: 0.2658, fgeo: 0.0387, fvol: 0.0639}
      }
      
      var corners = findCorners()
      
      viewAngles()
      solarPosition()
      sunZenOut()
      set('relativeSunViewAz', 'i.sunAz - i.viewAz')
      rossThick('kvol', 'i.sunZen', 'i.viewZen', 'i.relativeSunViewAz')
      rossThick('kvol0', 'i.sunZenOut', 0, 0)
      liThin('kgeo', 'i.sunZen', 'i.viewZen', 'i.relativeSunViewAz')
      liThin('kgeo0', 'i.sunZenOut', 0, 0)
      adjustBands()
      return image.select(inputBandNames)
      

      function viewAngles() {
        var maxDistanceToSceneEdge = 1000000
        var maxSatelliteZenith = 7.5
        var upperCenter = pointBetween(corners.upperLeft, corners.upperRight)
        var lowerCenter = pointBetween(corners.lowerLeft, corners.lowerRight)
        var slope = slopeBetween(lowerCenter, upperCenter)
        var slopePerp = ee.Number(-1).divide(slope)
        set('viewAz',
          ee.Image(ee.Number(Math.PI / 2).subtract((slopePerp).atan())))
    
        var leftLine = toLine(corners.upperLeft, corners.lowerLeft)
        var rightLine = toLine(corners.upperRight, corners.lowerRight)
        var leftDistance = ee.FeatureCollection(leftLine).distance(maxDistanceToSceneEdge)
        var rightDistance = ee.FeatureCollection(rightLine).distance(maxDistanceToSceneEdge)
        var viewZenith = rightDistance.multiply(maxSatelliteZenith * 2)
          .divide(rightDistance.add(leftDistance))
          .subtract(maxSatelliteZenith)
        set('viewZen',
          viewZenith.multiply(Math.PI).divide(180))
      }
    
      function solarPosition() {
        // Ported from http://pythonfmask.org/en/latest/_modules/fmask/landsatangles.html
        var date = ee.Date(ee.Number(image.get('system:time_start')))
        var secondsInHour = 3600
        set('longDeg',
          ee.Image.pixelLonLat().select('longitude'))
        set('latRad',
          ee.Image.pixelLonLat().select('latitude')
            .multiply(Math.PI).divide(180))
        set('hourGMT',
          ee.Number(date.getRelative('second', 'day')).divide(secondsInHour))
        set('jdp', // Julian Date Proportion
          date.getFraction('year'))
        set('jdpr', // Julian Date Proportion in Radians
          'i.jdp * 2 * {pi}')
        set('meanSolarTime',
          'i.hourGMT + i.longDeg / 15')
        set('localSolarDiff',
          '(0.000075 + 0.001868 * cos(i.jdpr) - 0.032077 * sin(i.jdpr)' +
          '- 0.014615 * cos(2 * i.jdpr) - 0.040849 * sin(2 * i.jdpr))' +
          '* 12 * 60 / {pi}')
        set('trueSolarTime',
          'i.meanSolarTime + i.localSolarDiff / 60 - 12')
        set('angleHour',
          'i.trueSolarTime * 15 * {pi} / 180')
        set('delta',
          '0.006918 - 0.399912 * cos(i.jdpr) + 0.070257 * sin(i.jdpr) - 0.006758 * cos(2 * i.jdpr)' +
          '+ 0.000907 * sin(2 * i.jdpr) - 0.002697 * cos(3 * i.jdpr) + 0.001480 * sin(3 * i.jdpr)')
        set('cosSunZen',
          'sin(i.latRad) * sin(i.delta) ' +
          '+ cos(i.latRad) * cos(i.delta) * cos(i.angleHour)')
        set('sunZen',
          'acos(i.cosSunZen)')
        set('sinSunAzSW',
          toImage('cos(i.delta) * sin(i.angleHour) / sin(i.sunZen)')
            .clamp(-1, 1))
        set('cosSunAzSW',
          '(-cos(i.latRad) * sin(i.delta)' +
          '+ sin(i.latRad) * cos(i.delta) * cos(i.angleHour)) / sin(i.sunZen)')
        set('sunAzSW',
          'asin(i.sinSunAzSW)')
        setIf('sunAzSW',
          'i.cosSunAzSW <= 0',
          '{pi} - i.sunAzSW',
          'sunAzSW')
        setIf('sunAzSW',
          'i.cosSunAzSW > 0 and i.sinSunAzSW <= 0',
          '2 * {pi} + i.sunAzSW',
          'sunAzSW')
        set('sunAz',
          'i.sunAzSW + {pi}')
        setIf('sunAz',
          'i.sunAz > 2 * {pi}',
          'i.sunAz - 2 * {pi}',
          'sunAz')
      }
    
      function sunZenOut() {
        // https://nex.nasa.gov/nex/static/media/publication/HLS.v1.0.UserGuide.pdf
        set('centerLat',
          ee.Number(
            ee.Geometry(image.get('system:footprint'))
              .bounds().centroid(30).coordinates().get(0))
            .multiply(Math.PI).divide(180))
        set('sunZenOut',
          '(31.0076' +
          '- 0.1272 * i.centerLat' +
          '+ 0.01187 * pow(i.centerLat, 2)' +
          '+ 2.40E-05 * pow(i.centerLat, 3)' +
          '- 9.48E-07 * pow(i.centerLat, 4)' +
          '- 1.95E-09 * pow(i.centerLat, 5)' +
          '+ 6.15E-11 * pow(i.centerLat, 6)) * {pi}/180')
      }
    
      function rossThick(bandName, sunZen, viewZen, relativeSunViewAz) {
        var args = {sunZen: sunZen, viewZen: viewZen, relativeSunViewAz: relativeSunViewAz}
        cosPhaseAngle('cosPhaseAngle', sunZen, viewZen, relativeSunViewAz)
        set('phaseAngle',
          'acos(i.cosPhaseAngle)')
        set(bandName,
          '(({pi}/2 - i.phaseAngle) * i.cosPhaseAngle + sin(i.phaseAngle)) ' +
          '/ (cos({sunZen}) + cos({viewZen})) - {pi}/4', args)
      }
    
      function liThin(bandName, sunZen, viewZen, relativeSunViewAz) {
        // From https://modis.gsfc.nasa.gov/data/atbd/atbd_mod09.pdf
        var args = {
          sunZen: sunZen,
          viewZen: viewZen,
          relativeSunViewAz: relativeSunViewAz,
          'h/b': 2,
        }
    
        anglePrime('sunZenPrime', sunZen)
        anglePrime('viewZenPrime', viewZen)
        cosPhaseAngle('cosPhaseAnglePrime', 'i.sunZenPrime', 'i.viewZenPrime', relativeSunViewAz)
        set('distance',
          'sqrt(pow(tan(i.sunZenPrime), 2) + pow(tan(i.viewZenPrime), 2)' +
          '- 2 * tan(i.sunZenPrime) * tan(i.viewZenPrime) * cos({relativeSunViewAz}))', args)
        set('temp',
          '1/cos(i.sunZenPrime) + 1/cos(i.viewZenPrime)')
        set('cosT',
          toImage('{h/b} * sqrt(pow(i.distance, 2) + pow(tan(i.sunZenPrime) * tan(i.viewZenPrime) * sin({relativeSunViewAz}), 2))' +
            '/ i.temp', args)
            .clamp(-1, 1))
        set('t', 'acos(i.cosT)')
        set('overlap',
          '(1/{pi}) * (i.t - sin(i.t) * i.cosT) * (i.temp)')
        setIf('overlap', 'i.overlap > 0', 0)
        set(bandName,
          'i.overlap - i.temp' +
          '+ (1/2) * (1 + i.cosPhaseAnglePrime) * (1/cos(i.sunZenPrime)) * (1/cos(i.viewZenPrime))')
      }
    
      function anglePrime(name, angle) {
        var args = {'b/r': 1, angle: angle}
        set('tanAnglePrime',
          '{b/r} * tan({angle})', args)
        setIf('tanAnglePrime', 'i.tanAnglePrime < 0', 0)
        set(name,
          'atan(i.tanAnglePrime)')
      }
    
      function cosPhaseAngle(name, sunZen, viewZen, relativeSunViewAz) {
        var args = {
          sunZen: sunZen,
          viewZen: viewZen,
          relativeSunViewAz: relativeSunViewAz
        }
        set(name,
          toImage('cos({sunZen}) * cos({viewZen})' +
            '+ sin({sunZen}) * sin({viewZen}) * cos({relativeSunViewAz})', args)
            .clamp(-1, 1))
      }
    
      function adjustBands() {
        for (var bandName in coefficientsByBand)
          applyCFactor(bandName, coefficientsByBand[bandName])
      }
    
      function applyCFactor(bandName, coefficients) {
        brdf('brdf', 'kvol', 'kgeo', coefficients)
        brdf('brdf0', 'kvol0', 'kgeo0', coefficients)
        set('cFactor',
          'i.brdf0 / i.brdf', coefficients)
        set(bandName,
          '{bandName} * i.cFactor', {bandName: 'i.' + bandName})
      }
    
      function brdf(bandName, kvolBand, kgeoBand, coefficients) {
        var args = merge(coefficients, {
          // kvol: 'i.' + kvolBand,
          kvol: '3 * i.' + kvolBand,     // check this multiplication factor.  Is there an 'optimal' value?  Without a factor here, there is not enough correction.
          kgeo: 'i.' + kgeoBand
        })
        return set(bandName,
          '{fiso} + {fvol} * {kvol} + {fgeo} * {kvol}', args)
      }
    
      function findCorners() {
        var footprint = ee.Geometry(image.get('system:footprint'))
        var bounds = ee.List(footprint.bounds().coordinates().get(0))
        var coords = footprint.coordinates()
    
        var xs = coords.map(function (item) {
          return x(item)
        })
        var ys = coords.map(function (item) {
          return y(item)
        })
    
        function findCorner(targetValue, values) {
          var diff = values.map(function (value) {
            return ee.Number(value).subtract(targetValue).abs()
          })
          var minValue = diff.reduce(ee.Reducer.min())
          var idx = diff.indexOf(minValue)
          return coords.get(idx)
        }
    
        var lowerLeft = findCorner(x(bounds.get(0)), xs)
        var lowerRight = findCorner(y(bounds.get(1)), ys)
        var upperRight = findCorner(x(bounds.get(2)), xs)
        var upperLeft = findCorner(y(bounds.get(3)), ys)
        return {
          upperLeft: upperLeft,
          upperRight: upperRight,
          lowerRight: lowerRight,
          lowerLeft: lowerLeft
        }
      }
  
      function x(point) {
        return ee.Number(ee.List(point).get(0))
      }
    
      function y(point) {
        return ee.Number(ee.List(point).get(1))
      }
    
      function pointBetween(pointA, pointB) {
        return ee.Geometry.LineString([pointA, pointB]).centroid().coordinates()
      }
    
      function slopeBetween(pointA, pointB) {
        return ((y(pointA)).subtract(y(pointB))).divide((x(pointA)).subtract(x(pointB)))
      }
    
      function toLine(pointA, pointB) {
        return ee.Geometry.LineString([pointA, pointB])
      }

// ************** COMMON HELPERS **************

      function set(name, toAdd, args) {
        toAdd = toImage(toAdd, args)
        image = image.addBands(toAdd.rename(name), null, true)
      }
    
      function setIf(name, condition, trueValue, falseValue) {
        condition = toImage(condition)
        var trueMasked = toImage(trueValue).mask(toImage(condition))
        var falseMasked = toImage(falseValue).mask(invertMask(condition))
        var value = trueMasked.unmask(falseMasked)
        set(name, value)
        
    
        function invertMask(mask) {
          return mask.multiply(-1).add(1)
        }
      }
    
      function toImage(band, args) {
        if ((typeof band) === 'string') {
          if (band.indexOf('.') > -1 || band.indexOf(' ') > -1 || band.indexOf('{') > -1) {
            band = image.expression(format(band, args), {i: image})
          } else
            band = image.select(band)
        }
        return ee.Image(band)
      }
    
      function format(s, args) {
        if (!args) args = {}
        var allArgs = merge(constants, args)
        var result = s.replace(/{([^{}]*)}/g,
          function (a, b) {
            var replacement = allArgs[b]
            if (replacement == null) {
              print('Undeclared argument: ' + b, 's: ' + s, args)
              return null
            }
            return allArgs[b]
          }
        )
        if (result.indexOf('{') > -1)
          return format(result, args)
        return result
      }
      
      function merge(o1, o2) {
        function addAll(target, toAdd) {
          for (var key in toAdd) target[key] = toAdd[key]
        }
    
        var result = {}
        addAll(result, o1)
        addAll(result, o2)
        return result
      }
    
      function show(band, min, max) {
        Map.addLayer(toImage(band), {min: min ? min : -1, max: max ? max : 1}, band)
      }
    }
    
function getCal(factors, bandName, factorNumber) {
  return ee.List(factors.get(bandName)).getNumber(factorNumber);
}

function applyCal(image, factors, bandName) {
  return ee.Image(getCal(factors, bandName, 0))
    .add(image.select(bandName).multiply(getCal(factors, bandName, 1)))
    .add(image.select(bandName).pow(2).multiply(getCal(factors, bandName, 2)))
    .add(image.select(bandName).pow(3).multiply(getCal(factors, bandName, 3)))
    .rename(bandName);
}

function applyCrossCalibrationLT5(image) {
  
  // Temperate mountains
  //var factors = ee.Dictionary({
    //RED: [-0.0042,	1.0049, 0.0000, 0.0000],
    // NIR: [0.0070,	0.8904,	0.4760,	-0.5883], 
    // not correcting because it degrades cross-sensor calibration
  //});
  
  // Arctic biome
    var factors = ee.Dictionary({
    RED: [-0.0068,	1.0075, 0.0000, 0.0000],
    //GREEN : [-0.0011, 0.8733, 0.3266, 0.0000],
    NIR: [0.0057,	0.9686,	0.0000,	0.0000], 
    //SWIR: [0.0010, 0.9791, 0.0000, 0.0000]
  });

  var cal = applyCal(image, factors, 'RED')
    .addBands(applyCal(image, factors, 'NIR'))
    //.addBands(applyCal(image, factors, 'GREEN'))
    //.addBands(applyCal(image, factors, 'SWIR'))

  return image.addBands(cal, null, true)
}

function applyCrossCalibrationLC8(image) {
  
  // Temperate mountains
  //var factors = ee.Dictionary({
  //  RED:   [0.0027,	1.0517,	-0.1658, 0.0000],
  //  NIR:   [0.0140,	0.8557, 0.0000, 0.0000],
  //});
  
  // Artic biome
  var factors = ee.Dictionary({
    RED: [-0.0012, 1.1592, -0.9166, 1.8795],
    //GREEN : [-0.0005, 1.0412, 0.0000, 0.0000],
    NIR: [0.0221, 0.8442, 0.1811, 0.0000], 
    //SWIR: [0.0071, 0.9698, 0.0685, 0.0000]
  });

  var cal = applyCal(image, factors, 'RED')
    .addBands(applyCal(image, factors, 'NIR'))
    //.addBands(applyCal(image, factors, 'GREEN'))
    //.addBands(applyCal(image, factors, 'SWIR'))

  return image.addBands(cal, null, true)
}

function addVariables(image) {
  // Compute time in fractional years since the epoch.
  var date = image.date();
  var years = date.difference(ee.Date('1970-01-01'), 'year');
  // Return the image with the added bands.
  return image
  // Add a time band.
  .addBands(ee.Image(years).rename('t')).float()
  // Add a constant band.
  .addBands(ee.Image.constant(1));
};

function fitHarmonic(LX,year){
  var window = windLength;
  var y = ee.Number(year);
  // Filter a time window of X years around the target year
  var yearWindowCol = LX.filter(
    ee.Filter.calendarRange(y.subtract(window),y.add(window),'year'))
  .map(addVariables);
  // First, add the harmonic variables (the third and fourth terms of equation 2)
  // to the image collection.
  // Use these independent variables in the harmonic regression.
  var harmonicIndependents = ee.List(['constant', 't', 'cos', 'sin']);
  // Add harmonic terms as new image bands.
  var harmonicLandsat = yearWindowCol.map(function(image) {
    var timeRadians = image.select('t').multiply(2 * Math.PI);
      return image
        .addBands(timeRadians.cos().rename('cos'))
        .addBands(timeRadians.sin().rename('sin'));
    });
  // Fit the model with a linear trend, using the linearRegression() reducer:
  // Name of the dependent variable.
  var dependent = ee.String('NDVI');
  var harmonicTrend = harmonicLandsat
    .select(harmonicIndependents.add(dependent))
    // The output of this reducer is a 4x1 array image.
    .reduce(ee.Reducer.linearRegression({
     numX: harmonicIndependents.length(),
     numY: 1
    }));
    
  // Turn the array image into a multi-band image of coefficients.
  var harmonicTrendCoefficients = harmonicTrend.select('coefficients')
    .arrayProject([0])
    .arrayFlatten([harmonicIndependents]);
  // Compute fitted values.
  var fittedHarmonic = harmonicLandsat.map(function(image) {
    return image.addBands(
      image.select(harmonicIndependents)
        .multiply(harmonicTrendCoefficients)
        .reduce('sum')
        .rename('fitted'));
  });
  return fittedHarmonic;
};

function adjustFit(yearCol,maxfit){
  // "adjustments factor" (difference btw fit and maxfit at every obs date)
  var adjustedFit = yearCol.map(function(image) {
    image = image.addBands(
      maxfit.subtract(image.select('fitted'))
      .rename('adjustment'));
    image = image.addBands(
      image.expression('adjusted = b("NDVI")+b("adjustment")'));
      return image;
  });
  return adjustedFit;
};

function maskCor(yearCol,themask){
  var maskYearCol = yearCol.map(function(image) {
    return image.updateMask(themask)
  });
  return maskYearCol
}

/* ---- II. Landsat ---- */

print( 'SECTION 2 - LANDSAT TIME SERIES')
var startYear =1984; // start of the time series
var endYear = 2024; // end of the time series
var startMonth = 6; // start of month to use
var endMonth = 9; // end of month to use
var CloudLT = 80; // don't use image with average cloud cover above ...
var demthr = 0; // mask pixels under ... m
var windLength = 8; // length of moving window (8 + 1 + 8 = 17 years)
var limObs = 100; // above ..., do not compute corrections (NOT USED)
var GrowSeasTh = 0.75; // remove observations too far from maximum NDVI
var minVI = 1500; // becarful as we multiply NDVI per 10000 (1500 = 0.15)

// LOAD LANDSAT COLLECTIONS
var L5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
.filter(ee.Filter.calendarRange(startMonth,endMonth,'month'))
.filter(ee.Filter.calendarRange(startYear,endYear,'year'))
.filterBounds(ROI)
.filterMetadata('CLOUD_COVER', 'less_than', CloudLT)
.select(['SR_B4', 'SR_B3','QA_PIXEL'], ['NIR', 'RED','QA_PIXEL']);

var L7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
.filter(ee.Filter.calendarRange(startMonth,endMonth,'month'))
.filter(ee.Filter.calendarRange(startYear,endYear,'year'))
.filterBounds(ROI)
.filterMetadata('CLOUD_COVER', 'less_than', CloudLT)
.select(['SR_B4', 'SR_B3','QA_PIXEL'], ['NIR', 'RED','QA_PIXEL']);

var L8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
.filter(ee.Filter.calendarRange(startMonth,endMonth,'month'))
.filter(ee.Filter.calendarRange(startYear,endYear,'year'))
.filterBounds(ROI)
.filterMetadata('CLOUD_COVER', 'less_than', CloudLT)
.select(['SR_B5', 'SR_B4', 'QA_PIXEL'], ['NIR', 'RED','QA_PIXEL']);

// APPLY CORRECTIONS AND COMPUTE NDVI
var L5c = L5
.map(maskLandsatSR)
.map(applyScaleFactors)
//.map(brdfCorrect)
.map(applyCrossCalibrationLT5)
.map(addNDVI)
//.map(addKNDVI_tanh)
//.map(addNIRv)
.select(['NDVI']);

var L7c = L7
.map(maskLandsatSR)
.map(applyScaleFactors)
//.map(brdfCorrect)
.map(addNDVI)
//.map(addKNDVI_tanh)
//.map(addNIRv)
.select(['NDVI']);

var L8c = L8
.map(maskLandsatSR)
.map(applyScaleFactors)
//.map(brdfCorrect)
.map(applyCrossCalibrationLC8)
.map(addNDVI)
//.map(addKNDVI_tanh)
//.map(addNIRv)
.select(['NDVI']);

var LXm = L5c.merge(L7c).merge(L8c)

print(LXm)

// prepare start and end of time series
var startDate = ee.Date('1984-01-01');
var endDate = ee.Date('2025-01-01');

// convert to millis
var start = startDate.millis()
var end = endDate.millis()

// how much is 8 days in millis ?
var millisIn8Days = 8*24*60*60*1000;

// create sequence of date advanced every 8 days
var DATEList = ee.List.sequence(start, end , millisIn8Days).map(function(dateMillis){
  return ee.Date(dateMillis).millis();
});

// create empty image based on dates
var LXtmp = DATEList.map(function(date){
  return ee.Image().set('system:time_start',date).rename('NDVI')
});

// convert list to image collection
var LXempty = ee.ImageCollection(LXtmp)

var LX = LXm.merge(LXempty).sort('system:time_start')
.filter(ee.Filter.calendarRange(startMonth,endMonth,'month'))

var LXDates = LX.map(function(img){
  var datenum = ee.Number.parse(img.date().format('YDH')); // .millis() too strict
  var out = img.multiply(0).add(datenum); // ee.Image.constant(datenum).toUint16().copyProperties(img, ['system:time_start'])
  return ee.Image(out).toInt().copyProperties(img, ['system:time_start']); 
});

/* ---- II. Computing NDVImax ---- */
print( 'SECTION 4 - COMPUTING NDVIMAX')
/*
var years = ee.List.sequence(startYear, endYear);

// compute max NDVI and DOY for every year 
var annualList = years.map(function(year) {
  
  // filter an annual subset of the collection 
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate = ee.Date.fromYMD(year, 12, 31);
  var yearCol = LX.filterDate(startDate, endDate);
  
  // reduce this collection to an image where pixels represent the annual max NDVI
  var ndviMax = yearCol.select('NDVI')
  .max()
  .rename('NDVI_max')
  .addBands(ee.Image.constant(ee.Number(year)).toInt().rename('year'))
  .set('year', year)

  return ndviMax;
});

var annualRaw = ee.ImageCollection.fromImages(annualList);


print(annualRaw)
var visS = {
  min: 0,
  max: 10000,
  palette: ['white','green']};
Map.addLayer(annualRaw.select('NDVI_max'),visS,'ndvimax')

var fit_raw = annualRaw.sort('year').select(['year','NDVI_max'])
    .reduce(ee.Reducer.sensSlope());
var slope_raw = fit_raw.select('slope').clip(ROI)

var visSlope = {
  min: -100,
  max: 100,
  palette: ['red','white','green']};
Map.addLayer(slope_raw,visSlope,'slope ')
*/
/*
Export.image.toDrive({
        image: slope_raw.clip(ROI),
        description: 'SLOPE_UMIUJAQ',
        folder:'SLOPE_UMIUJAQ',
        scale: 30,
        maxPixels:1e11,
        region: ROI
      });
*/


/* ---- III. Computing NDVImax from Harmonic model ---- */
print( 'SECTION 4bis - COMPUTING NDVIMAX FROM HARMONIC MODEL')

var years = ee.List.sequence(startYear, endYear);

var annualListHarmo = years.map(function(year) {
  
  ////////////////////////////////
  // ---- PREPARE DATASETS ---- //
  ////////////////////////////////
  
  // fit harmonic model to multiyear NDVI time series 
  var fittedHarmonic = fitHarmonic(LX,year);
  
  // extract raw ndvi and ndvi fitted by harmonic model for year of interest
  var yearCol = fittedHarmonic.select(['fitted', 'NDVI']).filter(
    ee.Filter.calendarRange(year,year,'year'));

  // count clean pixels to split data with cor / without cor
  //var countCol = yearCol.select('NDVI').count().rename('count')
  var avgVI = yearCol.select('NDVI').mean().rename('avgVI')
  
  // make mask to split between pixels on which to apply corrections or not
  var maskDonotCorrect = avgVI.lte(minVI)
  var maskDoCorrect = avgVI.gt(minVI)
  
  //////////////////////////////////////
  // ---- CORRECTION OF NDVI MAX ---- //
  //////////////////////////////////////
    
  // remove pixels with more than a given number of observations from correction loop [NOT USED !!!]
  var yearColdoCor = maskCor(yearCol,maskDoCorrect)
    
  // compute annual max from the fitted model to estimate "growing season"
  var maxfit = yearColdoCor.select('fitted').max();
  
  // compute the 'growing season' period based on maxfit (defined as the period that will be used to compute max)
  var quantfit = maxfit.multiply(GrowSeasTh)
  
  // keep only pixels that are within the 'growing season' to apply correction
  var yearColdoCor = yearColdoCor.map(function(image){
    var mask = image.gte(quantfit)
    return image.updateMask(mask)
  })
  
  var maxadj = adjustFit(yearColdoCor,maxfit).select('adjusted').max().rename('NDVI_max');

  /////////////////////////////////////////
  // ---- NO CORRECTION OF NDVI MAX ---- //
  /////////////////////////////////////////
  
  // remove pixels with less than 8 observations
  var yearColdonotCor = maskCor(yearCol,maskDonotCorrect)
  
  // keep max of unmasked
  var maxraw = yearColdonotCor.select('NDVI').reduce(ee.Reducer.percentile([90])).rename('NDVI_max');

  //////////////////////////////////////////////////////////
  // ---- MERGE CORRECTED AND NOT CORRECTED NDVI MAX ---- //
  //////////////////////////////////////////////////////////
  
  // merge both products
  var maxfinal = ee.ImageCollection([maxadj.select('NDVI_max'),maxraw.select('NDVI_max')])
  var maxfinal = maxfinal.median()
  
  /////////////////////////////////////////
  // ---- FINALIZE THE MAXIMUM NDVI ---- //
  /////////////////////////////////////////
  
  // reduce this collection to an image where pixels represent the annual max NDVI 
  // (or x-th quantile)
  var ndvimax = maxfinal
    .rename('NDVI_max')
    .toInt()
    // add year band (cast as same type as NDVI for linearFit)
    .addBands(ee.Image.constant(ee.Number(year)).rename('year')).toInt()
    // set the year as a property for further filtering
    .set('year', year)
    
  return ndvimax;
});

var annualHarmo = ee.ImageCollection.fromImages(annualListHarmo);
print(annualHarmo)

/*
var visS = {
  min: 0,
  max: 10000,
  palette: ['white','green']};
var toMap = annualHarmo.select('NDVI_max').filterMetadata('year', 'equals', 2024)
print(toMap)
Map.addLayer(toMap,visS,'ndvimax')


var fit_harmo = annualHarmo.sort('year').select(['year','NDVI_max'])
    .reduce(ee.Reducer.sensSlope());
var slope_harmo = fit_harmo.select('slope').clip(ROI)

var visS = {
  min: -100,
  max: 100,
  palette: ['red','white','green']};
Map.addLayer(slope_harmo,visS,'slope')
*/

/* ---- IV. Export results ---- */
print( 'SECTION 6 - EXPORT RESULTS')

/*
Export.image.toDrive({
image:slope_harmo,
description: 'slope_greening',
folder: 'POUR_MARILIE',
fileNamePrefix: 'slope_greening',
region: ROI,
scale: 30,
maxPixels: 1e10});

*/

var years = [1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024]

for(var y in years){ 
var year = ee.Number(years[y])
var exportYear = ee.String(year)
var im = annualHarmo
          .filter(ee.Filter.rangeContains('year',year, year))
          .first()
Export.image.toDrive({
        image: im.select('NDVI_max').toInt().clip(ROI),
        description: 'NDVImax_GRADIENT_NORDIQUE_ROI6_'+exportYear.getInfo(),
        folder:'NDVImax_GRADIENT_NORDIQUE_ROI6',
        scale: 30,
        maxPixels:1e11,
        region: ROI
      });
}
