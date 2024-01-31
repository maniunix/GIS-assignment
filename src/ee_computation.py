import shapely
from shapely.geometry import Polygon
import ee
import geemap

ee.Initialize()

import ee
import geemap

# Function to calculate NDWI using Sentinel-2 bands
def NDWIS2(image):
    """
    Calculates Normalized Difference Water Index (NDWI) using 
    Sentinel-2 bands B8 (Near-Infrared) and B12 (Shortwave Infrared).

    Parameters:
        image (ee.Image): Sentinel-2 image.

    Returns:
        ee.Image: Image with an added band representing NDWI.
    """
    ndwi = image.normalizedDifference(['B8', 'B12']).rename('NDWI')
    return image.addBands(ndwi)

# Function to calculate NDVI using Sentinel-2 bands
def NDVIS2(image):
    """
    Calculates Normalized Difference Vegetation Index (NDVI) using 
    Sentinel-2 bands B8 (Near-Infrared) and B4 (Red).

    Parameters:
        image (ee.Image): Sentinel-2 image.

    Returns:
        ee.Image: Image with an added band representing NDVI.
    """
    ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI')
    return image.addBands(ndvi)

# Function to mask clouds in NDVI image based on Sentinel-2 Scene Classification Map (SCL)
def sclCloudMask_ndvi(image):
    """
    Masks clouds in the NDVI image using 
    Sentinel-2 Scene Classification(SCL).

    Parameters:
        image (ee.Image): Sentinel-2 image with NDVI band.

    Returns:
        ee.Image: NDVI image with clouds masked.
    """
    ndvi = image.select('NDVI')
    crs = (ndvi.projection()).crs()
    scl = image.select('SCL')
    reScl = scl.resample('bilinear').reproject(crs=crs, scale=10)
    mask = reScl.gt(3) and reScl.lte(7)
    maskedNdvi = ndvi.updateMask(mask)
    return ee.Image(maskedNdvi)

# Function to mask clouds in NDWI image based on Sentinel-2 Scene Classification Map (SCL)
def sclCloudMask_ndwi(image):
    """
    Masks clouds in the NDWI image using
    Sentinel-2 Scene Classification(SCL).

    Parameters:
        image (ee.Image): Sentinel-2 image with NDWI band.

    Returns:
        ee.Image: NDWI image with clouds masked.
    """
    ndwi = image.select('NDWI')
    crs = (ndwi.projection()).crs()
    scl = image.select('SCL')
    reScl = scl.resample('bilinear').reproject(crs=crs, scale=10)
    mask = reScl.gt(3) and reScl.lte(7)
    maskedNdwi = ndwi.updateMask(mask)
    return ee.Image(maskedNdwi)

# Function to export an Earth Engine image to Google Drive
def saveEEImage(img, aoi, path):
    """
    Exports an Earth Engine image to Local System.

    Parameters:
        img (ee.Image): Earth Engine image.
        aoi (ee.Geometry): Area of interest.
        path (str): Google Drive path for export.

    Returns:
        None
    """
    geemap.ee_export_image(img, path, scale=10, region=aoi)
    return print('Saved Successfully')

# Landsat Functions

# Function to calculate NDVI using Landsat bands
def landsat_calculate_ndvi(image):
    """
    Calculates Normalized Difference Vegetation Index (NDVI) using 
    Landsat bands SR_B5 (Near-Infrared) and SR_B4 (Red).

    Parameters:
        image (ee.Image): Landsat image.

    Returns:
        ee.Image: Image with an added band representing NDVI.
    """
    ndvi = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI')
    return image.addBands(ndvi)

# Function to calculate NDWI using Landsat bands
def landsat_calculate_ndwi(image):
    """
    Calculates Normalized Difference Water Index (NDWI) using
    Landsat bands SR_B5 (Near-Infrared) and SR_B6 (Shortwave Infrared).

    Parameters:
        image (ee.Image): Landsat image.

    Returns:
        ee.Image: Image with an added band representing NDWI.
    """
    ndwi = image.normalizedDifference(['SR_B5', 'SR_B6']).rename('NDWI')
    return image.addBands(ndwi)


# MODIS Functions
# Function to calculate NDVI using MODIS bands
def calculate_ndvi_modis(image):
    """
    Calculates Normalized Difference Vegetation Index (NDVI) using 
    MODIS bands sur_refl_b02 (Near-Infrared) and sur_refl_b01 (Red).

    Parameters:
        image (ee.Image): MODIS image.

    Returns:
        ee.Image: Image with an added band representing NDVI.
    """
    ndvi = image.normalizedDifference(['sur_refl_b02', 'sur_refl_b01']).rename('ndvi')
    return image.addBands(ndvi)


# Function to calculate NDWI using MODIS bands
def calculate_ndwi_modis(image):
    """
    Calculates Normalized Difference Water Index (NDWI) using 
    MODIS bands sur_refl_b02 (Near-Infrared) and sur_refl_b03 (Shortwave Infrared).

    Parameters:
        image (ee.Image): MODIS image.

    Returns:
        ee.Image: Image with an added band representing NDWI.
    """
    ndwi = image.normalizedDifference(['sur_refl_b02', 'sur_refl_b03']).rename('ndwi')
    return image.addBands(ndwi)



def ExportEEtoDrive(eeimage,img_name, aoi):
        """
        Takes the Image as an Input and Save it to Drive

        Parameters:
            eeimage: ee.Image
            img_name: Output Name of the Image 
            aoi: ee.Geometry.Polygon
        """
        task = ee.batch.Export.image.toDrive(
        image=eeimage,
        description=img_name,
        folder='ee_demos',
        region=aoi,
        scale=250,
        crs='EPSG:4326'
    )
        task.start()
        return print("Image saved to Drive")
