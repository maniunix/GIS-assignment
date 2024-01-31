import ee
import geemap
import geopandas as gpd
import os
import glob
import shapely
from shapely.geometry import Polygon
from ee_computation import NDWIS2, NDVIS2, sclCloudMask_ndvi, sclCloudMask_ndwi, landsat_calculate_ndwi, landsat_calculate_ndvi, saveEEImage, calculate_ndvi_modis, calculate_ndwi_modis, ExportEEtoDrive


# Authenticate and initialize Earth Engine
# ee.Authenticate() 
ee.Initialize()

# Create a map for visualization using geemap
Map = geemap.Map()

# Function to flatten 3D geometry to 2D
def flatten_geom(geometry: shapely.geometry.Polygon):
    '''
    Takes the Geometry in 3D and Returns 2D Geometry

    Input: Geometry
    Output: list of Coordinates
    '''
    poly_wkt = geometry.exterior
    flat_geom = list(Polygon([xy[0:2] for xy in list(poly_wkt.coords)]).exterior.coords)
    return flat_geom

# Function to clip an image using a specified area of interest (aoi)
def clipImage(image):
    '''
    Clip the image using aoi
    Input: Image
    Output: Clipped Image
    '''
    return image.clip(aoi)

# Function to save an Earth Engine image to a local path
def saveEEImage(img: ee.Image, aoi: ee.Geometry.Polygon, path: str):
    '''
    Save the Image in a local path 
    Input: Image, Area of Interest, path
    '''
    geemap.ee_export_image(img, path, scale=10, region=aoi)
    return print('Saved Successfully')

# Function to get NDVI or NDWI index raster from Sentinel-2 imagery
def getIndexRaster(aoi: ee.Geometry.Polygon, sdate: str, edate: str, indice: str):
    '''
    Returns the Index Image i.e. NDVI, NDWI

    Input: Area of Interest, Start Date, End Date, Index [NDVI/NDWI]
    '''
    if indice == "NDVI":
        dataset = ee.ImageCollection('COPERNICUS/S2_SR').filterDate(sdate, edate).filterBounds(aoi)
        data = dataset.map(NDVIS2).map(sclCloudMask_ndvi).select('NDVI').map(clipImage).mean()
        saveEEImage(data, aoi, f'{indice}_image.tif')
    else:
        dataset = ee.ImageCollection('COPERNICUS/S2_SR').filterDate(sdate, edate).filterBounds(aoi)
        data = dataset.map(NDWIS2).map(sclCloudMask_ndwi).select('NDWI').map(clipImage).mean()
        saveEEImage(data, aoi, f'{indice}_image.tif')
    return data

# Read shapefile and convert to Earth Engine Geometry
shapefile_path = "AOI/Kalahandi_AOI.shp"
gdf = gpd.read_file(shapefile_path)
gdf = gdf.to_crs(4326)
flat_geometry = flatten_geom(gdf.geometry[0])
aoi = ee.Geometry.Polygon(flat_geometry)

# Set start and end dates for the analysis
sdate = "2022-01-01"
edate = "2022-12-31"

# Get NDVI and NDWI images from Sentinel-2
ndvi_s2_img = getIndexRaster(aoi, sdate, edate, "NDVI")
ndwi_s2_img = getIndexRaster(aoi, sdate, edate, "NDWI")

# Landsat 9 analysis
landsat9 = (
    ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
    .filterDate(sdate, edate)
    .filterBounds(aoi)
).map(landsat_calculate_ndvi).map(landsat_calculate_ndwi)

ndvi_image_l9 = landsat9.select('NDVI').median().clip(aoi)
ndwi_image_l9 = landsat9.select('NDWI').median().clip(aoi)

# Export Landsat NDVI image to Google Drive
ExportEEtoDrive(ndvi_image_l9, "NDVI_landsat",aoi)
ExportEEtoDrive(ndwi_image_l9, "NDWI_landsat",aoi)



# MODIS analysis
modis = ee.ImageCollection('MODIS/006/MOD09GA') \
    .filterDate(sdate, edate) \
    .filterBounds(aoi) \
    .select(['sur_refl_b01', 'sur_refl_b02'])

modis_scaled = modis.map(lambda img: img.multiply(0.0001))

modis_ndvi = modis_scaled.map(calculate_ndvi_modis).select('ndvi').mean().clip(aoi)
modis_ndwi = modis_scaled.map(calculate_ndwi_modis).select('ndwi').mean().clip(aoi)


ExportEEtoDrive(modis_ndvi,"modis_NDVI_image", aoi)
ExportEEtoDrive(modis_ndwi,"modis_NDWI_image", aoi)