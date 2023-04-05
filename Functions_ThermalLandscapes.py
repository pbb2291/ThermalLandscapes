import rasterio as rio
import numpy as np
import geopandas as gpd
from pathlib import Path
import xarray as xr
import rioxarray as rioxr
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import pandas as pd
import warnings
import seaborn as sns
from scipy.stats import ks_2samp
import time
from osgeo import gdal


# Define a simple evenness function
# compares the difference in coverage 
# where x is a pandas df
def evenness(x):
    # Filter for non-zero values
    X = x.loc[x>=0.01]
    # if more than 2:
    if len(X)>=2:
        # test difference between the values
        # This is the mean difference between the values
        # Divided by the sum of the values
        # 1 - is to make it an descending scale (higher values are better)
        # even = 1 - (np.mean(np.abs(np.diff((X)))) / np.sum(X))
        
        # Actually, screw that! Use coeff. of variation
         # NOTE: -1* need to rescale so scale is descending
        even = -1*(np.std(X)/np.mean(X))
    else:
        even = np.nan
    return even


def classifyVegType(h):
    
    h = np.round(h, 5)
    
    # "bare ground" pixels
    if ((h > -0.03) & (h <= 0.05)):
        veg = 'ground'
        
    # grass
    elif ((h > 0.05) & (h <= 0.5)):
        veg = 'grass'
        
    # shrub
    elif ((h > 0.5) & (h <= 3)):
        veg = 'shrub'
        
    # tree (just using the top bracket to filter out noise points)
    elif ((h > 3) & (h <= 20)):
        veg = 'tree'
        
    else:
        veg = np.nan
        
    return veg

def computeProximityRaster(inf, of):
# Distance to nearest ground from veg type (Ground, Grass, Woody) file
# Copied and pasted from: 
# https://github.com/marcelovilla/gdal-py-snippets/blob/master/scripts/compute_proximity.py
# https://gdal.org/programs/gdal_proximity.html

    ds = gdal.Open(inf, 0)
    
    try:
        
        band = ds.GetRasterBand(1)
        gt = ds.GetGeoTransform()
        sr = ds.GetProjection()
        cols = ds.RasterXSize
        rows = ds.RasterYSize

        # create empty proximity raster
        driver = gdal.GetDriverByName('GTiff')
        out_ds = driver.Create(of, cols, rows, 1, gdal.GDT_Float32)
        out_ds.SetGeoTransform(gt)
        out_ds.SetProjection(sr)
        out_band = out_ds.GetRasterBand(1)

        # compute proximity
        gdal.ComputeProximity(band, out_band, ['VALUES=1', 'DISTUNITS=PIXEL'])

        print(f'Done with {Path(inf).name}.\n Output {of}.\n')
        
    except Exception as e:
        
        print(f'Issue with {Path(inf).name}.\n\t{e}\n')
        



def readThermalImg(imgf):
        
    # Load row image
    img_xr = rioxr.open_rasterio(imgf,
                              parse_coordinates=True,
                              masked=True)

    # extract the temp band and convert to celcius
    img_c = (img_xr.sel(band=4)*0.04) - 273.15

    # filter 0 degree C and below as nans
    img_c = img_c.where(img_c > 0)

    return img_c

# for plotting difference in temp values between overlapping images
def plotDiffImages(pix1, pix2, diff, projstr, outstr1, outstr2):
    fig, ax = plt.subplots()
    pix1.plot.pcolormesh(ax=ax, alpha=0.5, add_colorbar=False)
    pix2.plot.pcolormesh(ax=ax, alpha=0.5, add_colorbar=False)
    diff.plot.pcolormesh(ax=ax, alpha=1, cmap='PiYG', add_colorbar=True)
    ax.axis('equal')
    fig.savefig(f'./figs/ImageDifferenceFigures/{projstr}/TempDiff_{projstr}_{outstr1}_{outstr2}.png', dpi=300)
    plt.close()


# Define a function to perform kstest
# NOTE: This isn't really comparing the exact overlap between the 2 images,
# Just the values in the overlap zone...
# For that reason, probably shouldn't trust it!
# 7/21/22
def KSTest_OverlappingPixels(polygon, pix1, pix2):
    
    # Read in both images
    img1 = readThermalImg(imgf1)
    img2 = readThermalImg(imgf2)

    # Clip the pixels of each using the intersection polygon
    pix1 = img1.rio.clip(polygon, drop=True)
    pix2 = img2.rio.clip(polygon, drop=True)
    
    # Flatten data and remove nas
    vals1 = pix1.data.flatten()[np.isnan(pix1.data.flatten(), where=False)]
    vals2 = pix2.data.flatten()[np.isnan(pix1.data.flatten(), where=False)]
    
    if ((vals1.size > 0) & (vals2.size > 0)):
        
        # Run a ks test to compare the values
        results = ks_2samp(vals1, vals2)
    
    else:
        
        results = np.nan
    
    # Return the ks results
    return results

# Function for saving square matrix results
def saveMatrixResults(results, projname, outstr, rowcolnames, m, n):
    # Combine results, and save
    arr = np.array(results)
    # Reshape into a matrix
    arr = arr.reshape(m, n)
    # make into a dataframe, labelling each row,col combination with the image names
    df = pd.DataFrame(arr, index=rowcolnames, columns=rowcolnames)
    # save
    df.to_csv(f'./data/out/{projstr}_{outstr}.csv')
    return df


# Cohen's D - A Measure of effect size 
# aka: a measure of standardized mean difference
# mean1 - mean2 / pooled sample Std
# https://stackoverflow.com/questions/21532471/how-to-calculate-cohens-d-in-python/33002123#33002123
# Other useful links:
# http://ethen8181.github.io/machine-learning/ab_tests/causal_inference/matching.html
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html
# https://en.wikipedia.org/wiki/Effect_size#Difference_family:_Effect_sizes_based_on_differences_between_means
def cohen_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx-1)*np.std(x, ddof=1) ** 2 + (ny-1)*np.std(y, ddof=1) ** 2) / dof)


def meanDiff_ImagesandExclosures(pix1_overlap, pix2_overlap, i1, i2, shp, printvalues=False, includepartial=False):
    
    # This Function clips the overlap zones of 2 images with exclosure polygons
    # It returns:
    # the pixel values (OPdiff_img1, OPdiff_img2)
    # The mean differences between images in the overlap zone (imgdiff) 
    # The mean differences between exclosure treatments (df_img12)
    # NOTE: This may not work well if the area overlaps more than 2 exclosure treatments (Open & Partial & Full all at once)
    
    # Inputs pix1_overlap and pix2_overlap are the overlapping pixels of 2 image files
    # i1 and i2 are the full paths to both img files
    # shp is a geodataframe of the shapefile polygons
    # print values will spit out a mean absolute difference if True
    
    # 7/25 - Filter shapefile for only the full and open exclosure types
    # Andrew doesn't want partial in the equation
    # if includepartial:
    #     shp_filter = shp.copy()
    # else:
    #     shp_filter = shp[shp.Exclosure != 'Partial']
    
    shp_filter=shp.copy()
    
    # Initialize df list for loop
    df_list = []

    # For each polygon treatment in the shapefile
    for ID, g in zip(shp_filter.index, shp_filter.geometry):

        # For each image with overlapping pixels
        for pix, i in zip([pix1_overlap, pix2_overlap], [i1, i2]):

            # Try clipping the image with the feature in the shapefile
            try:
                
                p = pix.rio.clip([g], shp.crs, drop=True)

                # Initialize an empty df
                df = pd.DataFrame()

                # img name
                iname = Path(i).name

                # Make an exploded dataframe, with each pixel noted by it's treatment and imgfile
                df = pd.DataFrame({'Temperature':p.data.flatten(),
                                   'Exclosure':shp_filter.Exclosure.iloc[ID],
                                   'imgf':iname})

                # Drop na rows
                df.dropna(axis=0, how='any', inplace=True)

                # store DataFrame in list
                df_list.append(df.copy(deep=True))

            # if it fails, no pixels, move on
            except:
                continue

    # Concat all the dfs to make a full df 
    df_img12 = pd.concat(df_list, ignore_index=True, sort=False)

    # set nodata values again (rioxarray changes it to 3.4028234663852886e+38)
    df_img12[df_img12 == 3.4028234663852886e+38] = np.nan

    # Group by Exclosure and Images
    df_img12_g = df_img12.groupby(['Exclosure','imgf'])
    df_img12_img_g = df_img12.groupby(['imgf'])
    
    # If there are 2 images and 2 exclosure types to work with:
    if (df_img12_g.mean().shape[0] >= 4):
    
        # Take the absolute mean difference in exclosures per image
        OPdiff_img1 = df_img12_g.mean().iloc[0] - df_img12_g.mean().iloc[2]
        OPdiff_img2 = df_img12_g.mean().iloc[1] - df_img12_g.mean().iloc[3]

        # Just return single value
        OPdiff_img1 = OPdiff_img1.values[0]
        OPdiff_img2 = OPdiff_img2.values[0]

        # Mean difference between images
        imgdiff = df_img12_img_g.mean().iloc[0] - df_img12_img_g.mean().iloc[1]
        imgdiff = imgdiff.values[0]
    
    else:
        
        OPdiff_img1 = np.nan
        OPdiff_img2 = np.nan
        imgdiff = np.nan
        
    # Print the means
    # What we want to see here is that the differences in exclosure temperatures
    # are consistent between images    
    if printvalues:
        print(f'Overlap Zones \n Exclosure Temp Difference: \n \t {Path(i1).name}: {OPdiff_img1} \n \t {Path(i2).name}: {OPdiff_img2}')
        print(f' Difference Between Images: \n \t {imgdiff}')
    
    # Return df of pixels, and mean difference values
    return OPdiff_img1, OPdiff_img2, imgdiff, df_img12