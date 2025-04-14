import glob
import os
import rasterio as rio
import rasterio.merge as merge
import math
import argparse

# Define function to actually perform merge on raster files by resolution:
def mergeRasterFiles(indir, outstr, outpath, epsg):

    # String Names of resolution in Files
    resstr = ['010m', '025m', '050m', '1m']

    for res in resstr:
        # Catch for Noise, because only 1m res for this one
        # Also, if files are missing, let's you know.
        try:
            ## Make list of files in current directory by Resolution
            filelist = []

            for file in glob.glob(indir + '/*' + res + '*.tif'):

                filelist.append(file)

            ### Merge Files Together
            # Set CRS in advance since it wasn't read in
            # EPSG 32736 corresponds to WGS84 UTM36S
            krugerCRS = rio.crs.CRS.from_epsg(epsg)

            # Make list of Raster Files to Merge
            filestomerge = []

            for i in filelist:
                currtif = rio.open(i, mode='r+')
                currtif.crs = krugerCRS
                # if its a "Count" data product set nodata to 0
                if "Count" in i:
                    currtif.nodata = 0
                else:
                    currtif.nodata = -9999
                
#                 # Note: added for the special case of Catherine's Karingani data 08/05/21
#                 # Which needs to be divided by 1000
#                 if "SouthSouth" in i:
#                     currtif.calc
                
                filestomerge.append(currtif)
                
            # Merge
            mosaic, out_trans = merge.merge(filestomerge, nodata=currtif.nodata)

            # Write the output raster
            with rio.open(outpath + '/' + outstr + res + '.tif',
                          mode='w+',
                          driver='Gtiff',
                          width=mosaic.shape[2],
                          height=mosaic.shape[1],
                          crs=krugerCRS,
                          count=1,
                          dtype='float32',
                          transform=out_trans,
                          nodata=-9999) as dest:
                dest.write(mosaic)
                
        except:    
#         except Exception as e:
#             print(e.with_traceback, e.args)
            print(f"Unable to merge {outstr} files at {res} resolution.")

            continue

# Define function to Create folders and Directories
def createMergedFilesandFolder(indir, outdir):
    # Get absolute path from OS, in case it's a relative path:
    indir = os.path.abspath(indir) + '/'
    outdir = os.path.abspath(outdir) + '/'
    # print(f'Indir: {indir}')

    # Make Indir String Names for merge function

    # Aka: Gather all folder names that contain tif files in current directory
    indirlist = []
    # at the same time, make out strings for all files
    outstrlist = []
    for dirName in os.walk(indir):
        if len(dirName[2]) != 0:
            if '.tif' in dirName[2][0]:
                indirlist.append(dirName[0] + '/')
                tempoutname = dirName[2][0].split('_')
                outstrlist.append(tempoutname[0] + '_' + tempoutname[1] + '_')

    # Make output directories for merged files
    pcountnames = ['all', 'ground', 'noise', 'pulse']
    foldernames = ['CHM', 'DSM', 'DTM', 'Intensity', 'PCount']

    # Make an output directories list
    outdirlist = []
    for i in foldernames:
        if ((i != 'PCount') | (i != 'Pcount')):
            outdirlist.append(os.path.abspath(outdir + i + '/') + '/')

        if ((i == 'PCount') | (i == 'Pcount')):
            for j in pcountnames:
                outdirlist.append(os.path.abspath(outdir + i + '/' + j + '/') + '/')

    # Make folders, if they don't exist:
    # if not os.path.exists(outdir + 'Merged/'):
    #     os.mkdir(os.path.abspath(outdir + 'Merged/') + '/')

    for outname in outdirlist:
        if os.path.exists(outname) == False:
            os.mkdir(os.path.abspath(outname))
        else:
            print(f'{os.path.abspath(outname)} already exists.')
            continue

    # Remove the extra foldername (/PCount/) in outdirlist
    outputdirs = outdirlist[0:4] + outdirlist[5:]

    # Print a text file reporting what went where:
    with open(outdir + "/MergeReport.txt", "w") as text_file:
        print(f"Indir:\n {indirlist},\n Outstr=\n {outstrlist},\n Outpath=\n {outputdirs}", file=text_file)

    # iterate over string resolution names
    # uses zip to iterate over multiple lists:
    # https://www.geeksforgeeks.org/python-iterate-multiple-lists-simultaneously/
    for infolder, outstr, outpath in zip(indirlist, outstrlist, outputdirs):
        try:
            mergeRasterFiles(infolder, outstr, outpath)
        except:
            print(f"Cannot merge files in {infolder}.")

# Define function to Create folders and Directories
def mergefiles(indir, outdir):

    # Get absolute path from OS, in case it's a relative path:
    indir = os.path.abspath(indir) + '/'
    outdir = os.path.abspath(outdir) + '/'

    # Get all tifs in the current directory:
    flist = []
    for f in glob.glob(indir + '*.tif'):
        flist.append(f)

    # String Names of resolution in Files
    resstr = ['010m', '025m', '050m', '1m']

    for res in resstr:
        # Catch for Noise, because only 1m res for this one
        # Also, if files are missing, let's you know.
        try:
            ## Make list of files in current directory by Resolution
            filelist = []

            for file in glob.glob(indir + '*' + res + '*.tif'):
                filelist.append(file)

            ### Merge Files Together
            # Set CRS in advance since it wasn't read in
            # EPSG 32736 corresponds to WGS84 UTM36S
            krugerCRS = rio.crs.CRS.from_epsg(32736)

            # Make list of Raster Files to Merge
            filestomerge = []

            for i in filelist:
                currtif = rio.open(i, mode='r+')
                currtif.crs = krugerCRS
                # if its a "Count" data product set nodata to 0
                if "Count" in i:
                    currtif.nodata = 0
                else:
                    currtif.nodata = -9999

                filestomerge.append(currtif)

            # Merge
            mosaic, out_trans = merge.merge(filestomerge, nodata=-9999)

            # Write the output raster
            with rio.open(outpath + outstr + res + '.tif',
                          mode='w+',
                          driver='Gtiff',
                          width=mosaic.shape[2],
                          height=mosaic.shape[1],
                          crs=krugerCRS,
                          count=1,
                          dtype='float32',
                          transform=out_trans,
                          nodata=-9999) as dest:
                dest.write(mosaic)
        except:
            print(f"Unable to merge files in {indir} at {res} resolution. \n")
            continue

# UTM Zone finder function
# see https://stackoverflow.com/questions/40132542/get-a-cartesian-projection-accurate-around-a-lat-lng-pair
def getUTMcodefromWGS84(lon, lat):
    utm_band = str((math.floor((lon + 180) / 6 ) % 60) + 1)
    if len(utm_band) == 1:
        utm_band = '0'+utm_band
    if lat >= 0:
        epsg_code = '326' + utm_band
    else:
        epsg_code = '327' + utm_band
    return epsg_code

# Resample to 20 m pixels (downsampling)
# https://gis.stackexchange.com/questions/329945/should-resampling-downsampling-a-raster-using-rasterio-cause-the-coordinates-t
# https://gis.stackexchange.com/questions/329434/creating-an-in-memory-rasterio-dataset-from-numpy-array/329439#329439
def resample_raster(raster, scale=10, outpath='./test_DTM.tif'):
    t = raster.transform

    # rescale the metadata
    transform = Affine(t.a * scale, t.b, t.c, t.d, t.e * scale, t.f)
    height = raster.height // scale
    width = raster.width // scale

    profile = raster.profile
    profile.update(transform=transform, driver='GTiff', height=height, width=width)

    # Note changed order of indexes, arrays are band, row, col order not row, col, band
    data = raster.read(
            out_shape=(raster.count, height, width),
            resampling=Resampling.bilinear,
        )

    # Write the output resampled raster
    with rio.open(outpath, 'w+', **profile) as dest:
        dest.write(data)
