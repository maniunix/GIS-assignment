# For Earth Engine
import ee
import geemap

# For Vector Operations
import geopandas as gpd

# For File Operations
import os
import glob

import shapely
from shapely.geometry import Polygon
from ee_computation import NDWIS2, NDVIS2, sclCloudMask_ndvi, sclCloudMask_ndwi, landsat_calculate_ndwi, landsat_calculate_ndvi, saveEEImage, calculate_ndvi_modis,calculate_ndwi_modis

ee.Initialize()
Map = geemap.Map()


def flatten_geom(geometry) -> shapely.geometry.Polygon:
    '''
    Takes the Geometry in 3D and Returns 2D Geometry

    Input: Shapely.Geometry.Polygon
    Output: list of Coordinates
    '''
    poly_wkt = geometry.exterior
    flat_geom = list(Polygon([xy[0:2]
                     for xy in list(poly_wkt.coords)]).exterior.coords)
    return flat_geom


def clipImage(image):
    return image.clip(aoi)


def saveEEImage(img, aoi, path):
    geemap.ee_export_image(img,
                           path,
                           scale=10,
                           region=aoi)
    return print('Saved Successfully')


def getIndexRaster(aoi, sdate, edate, indice):
    if indice == "NDVI":
        dataset = ee.ImageCollection(
            'COPERNICUS/S2_SR').filterDate(sdate, edate).filterBounds(aoi)
        data = dataset.map(NDVIS2).map(sclCloudMask_ndvi).select(
            'NDVI').map(clipImage).mean()
        saveEEImage(data, aoi, f'{indice}_image.tif')
        # geemap.ee_export_image(data, f'{indice}_image.tif', scale = 10, region = aoi)
    else:
        dataset = ee.ImageCollection(
            'COPERNICUS/S2_SR').filterDate(sdate, edate).filterBounds(aoi)
        data = dataset.map(NDWIS2).map(sclCloudMask_ndwi).select(
            'NDWI').map(clipImage).mean()
        saveEEImage(data, aoi, f'{indice}_image.tif')
        # geemap.ee_export_image(data, f'{indice}_image.tif', scale = 10, region = aoi)
    return data


shapefile_path = "AOI/Kalahandi_AOI.shp"

gdf = gpd.read_file(shapefile_path)
gdf = gdf.to_crs(4326)

flat_geometry = flatten_geom(gdf.geometry[0])

aoi = ee.Geometry.Polygon(flat_geometry)

sdate = "2022-01-01"
edate = "2022-12-01"

ndvi_s2_img = getIndexRaster(aoi, sdate, edate, "NDVI")
ndwi_s2_img = getIndexRaster(aoi, sdate, edate, "NDWI")

###### Landsat ######

landsat9 = (
    # ee.ImageCollection('LANDSAT/LC09/C02/T1_TOA')
    ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
    .filterDate(sdate, edate)
    .filterBounds(aoi)
).map(landsat_calculate_ndvi).map(landsat_calculate_ndwi)

ndvi_image = landsat9.select('NDVI').median().clip(aoi)
ndwi_image = landsat9.select('NDWI').median().clip(aoi)

task = ee.batch.Export.image.toDrive(
    image=ndvi_image,
    description='ndwi_export',
    folder='ee_demos',
    region=aoi,
    scale=30,
    crs='EPSG:4326'
)

task.start()


#### MODIS ####
modis = ee.ImageCollection('MODIS/006/MOD09GA') \
    .filterDate('2022-01-01', '2022-01-30') \
    .filterBounds(aoi) \
    .select(['sur_refl_b01', 'sur_refl_b02']) 

modis_scaled = modis.map(lambda img: img.multiply(0.0001))

modis_ndvi = modis_scaled.map(calculate_ndvi_modis).select('ndvi').mean().clip(aoi)

task = ee.batch.Export.image.toDrive(
    image=modis_ndvi,
    description='modis_ndvi',
    folder='ee_demos',
    region=aoi,
    scale=250,
    crs='EPSG:4326'
)

task.start()