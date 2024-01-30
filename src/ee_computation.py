import shapely
from shapely.geometry import Polygon
import ee
import geemap

ee.Initialize


def NDWIS2(image):
    ndwi = image.normalizedDifference(['B8','B12']).rename('NDWI')
    return image.addBands(ndwi)


def NDVIS2(image):
    ndvi = image.normalizedDifference(['B8','B4']).rename('NDVI')
    return image.addBands(ndvi)

def sclCloudMask_ndvi(image):
    ndvi = image.select('NDVI')
    crs = (ndvi.projection()).crs()
    scl = image.select('SCL')
    reScl = scl.resample('bilinear').reproject(crs = crs, scale = 10)
    mask = reScl.gt(3) and reScl.lte(7)
    maskedNdvi = ndvi.updateMask(mask)
    return ee.Image(maskedNdvi)


def sclCloudMask_ndwi(image):
    ndvi = image.select('NDWI')
    crs = (ndvi.projection()).crs()
    scl = image.select('SCL')
    reScl = scl.resample('bilinear').reproject(crs = crs, scale = 10)
    mask = reScl.gt(3) and reScl.lte(7)
    maskedNdvi = ndvi.updateMask(mask)
    return ee.Image(maskedNdvi)


def saveEEImage(img, aoi, path):
    geemap.ee_export_image(img,
                           path,
                           scale=10,
                           region=aoi)
    return print('Saved Successfully')

##### Landsat Functions #####
def landsat_calculate_ndvi(image):
    ndvi = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI')
    return image.addBands(ndvi)

# # Function to calculate NDWI
def landsat_calculate_ndwi(image):
    ndwi = image.normalizedDifference(['SR_B5', 'SR_B6']).rename('NDWI')
    return image.addBands(ndwi)

### MODIS ###
def calculate_ndvi_modis(image):
    ndvi = image.normalizedDifference(['sur_refl_b02', 'sur_refl_b01']).rename('ndvi')
    return image.addBands(ndvi)

def calculate_ndwi_modis(image):
    ndwi = image.normalizedDifference(['sur_refl_b02', 'sur_refl_b03']).rename('ndwi')
    return image.addBands(ndwi)