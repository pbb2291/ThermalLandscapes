# Import Libraries
import os
import glob
import datetime
from osgeo import gdal, osr
import json
from matplotlib import pyplot as plt
import numpy as np
import subprocess
# import pdal
import laspy
import PyPDF2
import pandas as pd
import geopandas as gpd
import zipfile
import fiona
import georinex
import math

import Classes
import Functions
import struct
from datetime import datetime, timedelta
import geojson
import rasterio as rio
from rasterio.mask import mask as riomask
from rasterio.features import shapes
from rasterio.enums import Resampling
import seaborn as sns
import shutil

class Trajectory:

    def __init__(self, projdir=None, loadfromjson=False):

        if loadfromjson == False:

            self.projdir = projdir
            self.alltrajs = None
            self.alltrajs_utm = None

            # Field names for all sbet files
            # Fields for 17 columns of data
            # Note: lat and lon are in radians, make a new field for deg and convert with math.deg()
            self.field_names = ('time', 'latitude', 'longitude', 'altitude',
                           'x_vel', 'y_vel', 'z_vel',
                           'roll', 'pitch', 'platform_heading', 'wander_angle',
                           'x_acceleration', 'y_acceleration', 'z_acceleration',
                           'x_angular_rate', 'y_angular_rate', 'z_angular')

            # find all sbet files
            self.findsbetfiles()

            # Get starting positions and time durations of trajectories
            self.getinitialpositions()
            self.getflightduration()

            # Load in all the full trajectories (disabled, but it's fairly fast to load if you need them)
            # self.getalltrajs()  # in WGS84 - fast
            # self.getUTMtraj()  # in UTM - slower

    def findsbetfiles(self):

        # find all sbet files:
        self.sbet_files = []
        for f in glob.glob(self.projdir + "/**/*SBET_*.OUT", recursive=True):
            self.sbet_files.append(f)

    def getinitialpositions(self):

        # initiate dict of initial values of each sbet
        sbet_dict = {k: [] for k in self.field_names}

        # for each sbet file
        for f in glob.glob(self.projdir + "/**/*SBET*.OUT", recursive=True):

            # Open up the current sbet to get the initial x,y,z location
            sbet_file = open(f, 'rb')
            sbet_data = sbet_file.read()

            # Unpack the first packet of data in the file
            # Note: 17d means it's looking for 17 fields of double precision numbers, 8 bytes each
            values = struct.unpack('17d', sbet_data[0:8 * 17])

            # Add the values to the fields in a dict
            for key, value in zip(self.field_names, values):
                # Turn lat long radians to degrees (more understandable)
                if key == 'latitude':
                    sbet_dict[key].append(math.degrees(value))
                if key == 'longitude':
                    sbet_dict[key].append(math.degrees(value))
                elif (key != 'latitude') and (key != 'longitude') and (key != 'filepath'):
                    sbet_dict[key].append(value)

        self.startpositions = sbet_dict

    def getflightduration(self):

        # get all trajectories (if not already made)
        if not self.alltrajs:
            self.getalltrajs()

        # Get the FIRST and LAST position from the Sbet dict:
        idx_first = np.where(self.alltrajs['time'] == np.min(self.alltrajs['time']))[0][0]
        idx_last = np.where(self.alltrajs['time'] == np.max(self.alltrajs['time']))[0][0]

        self.duration_sec = self.alltrajs['time'][idx_last] - self.alltrajs['time'][idx_first]
        self.duration_hrs = self.duration_sec/3600

    def getUTMtraj(self):
        
        if not hasattr(self, 'alltrajs'):
            self.getalltrajs()

        # Convert to geodf
        sbet_geodf = gpd.GeoDataFrame(self.alltrajs['time'],
                                      geometry=gpd.points_from_xy(self.alltrajs['longitude'],
                                                                  self.alltrajs['latitude'],
                                                                  self.alltrajs['altitude'],
                                                                  crs="EPSG:4326")
                                      )

        # determine utmzone of project (if not already known)
        utmzone = Functions.getUTMcodefromWGS84(lon=self.startpositions['longitude'][0],
                                                lat=self.startpositions['latitude'][0])

        # Reproject to UTM
        sbet_geodf = sbet_geodf.to_crs(f"EPSG:{utmzone}")

        # Steal the coordinates for a dict, and stick them in:
        self.alltrajs_utm = {'easting': sbet_geodf.geometry.x,
                             'northing': sbet_geodf.geometry.y}

    # Add a function for loading all trajectories and plotting
    def getalltrajs(self):

        # loop through data packets and store info
        # Initilize Empty Dict (note: DON"T use .fromkeys() function!)
        # Or with a lambda function:
        sbet_dict = {k: [] for k in self.field_names}

        # Great walk-through on how to read sbet file with struct.unpack:
        # http: // vislab - ccom.unh.edu / ~schwehr / rt / python - binary - files.html
        for f in self.sbet_files:

            sbet_file = open(f, 'rb')
            sbet_data = sbet_file.read()

            # Get number of data packets (ie: rows of data, each with 17 fields)
            # Each field stores a value of 8 bytes, a double precision value
            ndatapackets = len(sbet_data) / (8 * 17)

            # for every 10th row of data
            for i in range(0, int(ndatapackets), 10):
                # offset tells it where to start
                offset = int((i) * (8 * 17))
                # get current row of values
                values = list(struct.unpack('17d', sbet_data[(offset + 0):(offset + (8 * 17))]))
                # stick in a dictionary
                for key, value in zip(self.field_names, values):
                    # Turn lat long radians to degrees (cause more understandable)
                    if key == 'latitude':
                        sbet_dict[key].append(math.degrees(value))
                    if key == 'longitude':
                        sbet_dict[key].append(math.degrees(value))
                    elif (key != 'latitude') and (key != 'longitude'):
                        sbet_dict[key].append(value)

        self.alltrajs = sbet_dict

    def plottraj(self, utm=False, fig=None, ax=None):

        if not hasattr(self, 'alltrajs'):
            self.getalltrajs()

        if utm == True:

            if not self.alltrajs_utm:
                self.getUTMtraj()

            # plot WGS84 lat long
            if not ax:
                fig, ax = plt.subplots()

            a = ax.scatter(x=self.alltrajs_utm['easting'],
                           y=self.alltrajs_utm['northing'],
                           c=self.alltrajs['time'],
                           s=1.5,
                           alpha=0.6)

            ax.set_xlabel('Easting')
            ax.set_ylabel('Northing')

            ax.set_aspect('equal', 'box')

            fig.colorbar(a, label='GPS Time')

        if utm == False:

            # plot WGS84 lat long
            if not ax:
                fig, ax = plt.subplots()

            a = ax.scatter(x=self.alltrajs['longitude'],
                           y=self.alltrajs['latitude'],
                           c=self.alltrajs['time'],
                           s=1.5,
                           alpha=0.6)

            ax.set_xlabel('longitude')
            ax.set_ylabel('latitude')

            ax.set_aspect('equal', 'box')

            fig.colorbar(a, label='GPS Time')

    def flatten(self):

        # initialize empty dict to fill
        dict = {}

        # for all attributes
        for attr in self.__dict__:
            if ((attr != 'alltrajs') & (attr != 'allUTMtrajs')):
                # stick them in the dictionary
                dict[attr] = self.__getattribute__(attr)

        # Return the dictionary object of itself
        return dict


class BaseStation:

    def __init__(self, projdir=None, ausposfile=None, loadfromjson=False):

        if loadfromjson == False:

            self.projdir = projdir

            # Get your location
            self.getlocation()

            if ausposfile == None:

                # Find the auspos pdf
                self.ausposfile = glob.glob(self.projdir + "/**/**stonex/*.pdf", recursive=True)[0]

                if not self.ausposfile:

                    print('No base station auspos pdf report file found.')

            else:

                self.ausposfile = ausposfile

            try:

                # Get your positional error
                self.getposerror()

            except:

                print(f'Unable to extract information from base station pdf file: {self.ausposfile} \n')

    def getlocation(self):

        # Read Rinex File
        self.rinexfile = glob.glob(self.projdir + '/**/*.*O', recursive=True)[0]
        rinex_head = georinex.rinexheader(self.rinexfile)

        # First row of the geodataframe is the BaseStation
        # Other rows are the rest of the trajectory (can store it if you want it, but doesn't have roll pitch yaw or accuracy, unlike sbet)
        self.lon = rinex_head['position'][0]
        self.lat = rinex_head['position'][1]
        self.alt = rinex_head['position'][2]
        
        # set Epsg for ITRF2014
        self.epsg = "EPSG:7789"

        # Get time information for base station
        # Note:times are datetime objects, use strftime function to convert to a string
        self.startdatetime = rinex_head['t0'].strftime("%m/%d/%Y %H:%M:%S")
        
        try:
            
            self.enddatetime = rinex_head['t1'].strftime("%m/%d/%Y %H:%M:%S")
            # Calc duration, note: make sure to do end - start time here, not the other way around
            duration = rinex_head['t1'] - rinex_head['t0']
            self.duration_sec = duration.seconds
            self.duration_hrs = duration.seconds/3600
            
        except:
            
            print(f'Missing endtime information in basestation file {self.rinexfile} \n')
            
        self.antenna = rinex_head['ANT # / TYPE']

    def getposerror(self):

        # Read Pdf
        # Note - this is not a sureproof way, it's really difficult to read in pdfs with python
        # but as long as report is auspos, and formats stay roughly the same, then we are good.
        pdf = PyPDF2.PdfFileReader(self.ausposfile)

        # Get Uncertainty of position
        p4 = pdf.getPage(3)
        text4 = p4.extractText()

        # Split text from third page to find station location
        text4 = text4.split('\n') # split text into list by line
        for txt, idx in zip( text4, range(len(text4))):
            if txt == "StationLongitude(East)(m)Latitude(North)(m)EllipsoidalHeight(Up)(m)":
                self.poserrstr = text4[idx + 1] # get the next row, with the error info
                self.xerr = float(self.poserrstr[4:9])
                self.yerr = float(self.poserrstr[9:15])
                self.zerr = float(self.poserrstr[15:21])
                self.radialerr = round((self.xerr**2 + self.yerr**2 + self.xerr**2)**0.5, 3) # the max error you could be from your position

    def flatten(self):

        # initialize empty dict to fill
        dict = {}

        for attr in self.__dict__:

                # stick it in the dictionary
                dict[attr] = self.__getattribute__(attr)

        # Return the dictionary object of itself
        return dict


class MergedTile:

    def __init__(self, filepath=None, epsg="32736"):

        self.filepath = filepath

        # NOTE: Need to adjust this so that you can set the epsg for each project
        # Needs to come from the project, not the merged tile class
        self.epsg = epsg

        # Make path string consistent with abspath()
        if filepath != None:

            self.filepath = os.path.abspath(filepath)

            # Grab information from filename
            self.CollectFileInfo()

            # Calculate raster information
            self.CalcRastInfo()

    def CollectFileInfo(self):

        # Get file name from filepath for consistency
        self.filename = self.filepath.split('/')[-1]

        # Set Processing Date
        currdate = datetime.now()  # curr date
        self.dateprocess = currdate.strftime("%Y-%m-%d")  # format current time

        # Get info from filename
        self.res = self.filename.split('_')[-1].split('.')[0]
        self.site = self.filename.split('_')[0]
        self.product = self.filename.split('_')[1]

        # Set sensor
        lidarstrlist = ['CHM', 'DSM', 'DTM', 'Intensity', "count"]
        if any(lidarstr in self.filename for lidarstr in lidarstrlist):
            self.sensor = "Lidar - Riegl Vux1LR"

        if ( ("RGB" in self.filename) | ("rgb" in self.filename) ) :
            self.sensor = "RGB Cam - Sony A6000"

        if ( ("thermal" in self.filename) | ("Thermal" in self.filename) ):
            self.sensor = "Thermal Cam - FLIR Tau 2 640"

    def CalcRastInfo(self):

        # Open with Gdal
        tif = gdal.Open(self.filepath)
        band = tif.GetRasterBand(1)
        # self.stats = band.GetMetadata_List()  # gets list of approx stats
        stats = band.ComputeStatistics(True)  # Calcs true stats (turn off true for approx stats)
        self.min = round(stats[0], 3)
        self.max = round(stats[1], 3)
        self.avg = round(stats[2], 3)
        self.std = round(stats[3], 3)

        # Get CRS and Spatial Ref
        # self.epsg = epsg # NOTE: commented out because we do this up front now
        spatref = osr.SpatialReference()
        spatref.ImportFromEPSG(int(self.epsg))
        self.crs = spatref.GetName()

        # Get Bounding Box and Center Coords
        # Note: This will probably only work as long as your in UTM, but may not in lat long... see:
        # https://gis.stackexchange.com/questions/57834/how-to-get-raster-corner-coordinates-using-python-gdal-bindings
        ulx, xres, xskew, uly, yskew, yres = tif.GetGeoTransform()
        lrx = ulx + (tif.RasterXSize * xres)  # number of cols * res
        lry = uly + (tif.RasterYSize * yres)  # note, yres is negative for some reason...
        self.ulcoords = [ulx, uly]  # upper left coords
        self.lrcoords = [lrx, lry]  # lower right coords
        self.centercoords = [np.mean([ulx, lrx]), np.mean([uly, lry])]  # center coords

        # Extents for matplotlib - plt.imshow(CHM, extent=self.extent)
        # REORDERED 9/27/2020
        # note: correct order is: left, right, top , bottom
        self.extent = [self.ulcoords[1],
                       self.lrcoords[1],
                       self.ulcoords[0],
                       self.lrcoords[0]]

        # Get Number of Hectares covered:
        bandarray = tif.ReadAsArray()
        nodatavalue = band.GetNoDataValue()
        self.datafootprint_msq = (abs(xres) * abs(yres)) * len(bandarray[bandarray != nodatavalue]) # data footprint in m2 (no data rmvd)
        self.datafootprint_msq_sum = (abs(xres) * abs(yres)) * np.sum(bandarray != nodatavalue)
        self.datafootprint_ha = self.datafootprint_msq / 10000  # convert to hectares

    # Method for Getting Band data
    def GetTifandBand(self):

        # Open with Gdal
        self.tif = gdal.Open(self.filepath)
        self.band = self.tif.GetRasterBand(1)
        self.bandarray = self.band.ReadAsArray()

    # A quick plot function for mergedtiles, v1
    def plot(self, cmin=None, cmax=None, cmap='magma', interp=None, returnaxis=False, fig=None, ax=None):

        # If no color scale provided, make one
        if cmin == None:
            cmin = self.min
            cmax = self.max

        # Get tif and band info from gdal
        self.GetTifandBand()

        # read band data as array
        data = self.band.ReadAsArray()

        # Convert to float if not already
        data = data.astype('float32')

        # set nodata to nan
        data[data == self.band.GetNoDataValue()] = np.nan

        # plot
        # if there's already a fig and axis to plot on
        if ((fig == None) & (ax == None)):

            # make them
            fig, ax = plt.subplots()

        # else, use the fig and axis inputs
        if ((fig != None) & (ax != None)):

            a = ax.imshow(data, cmap=cmap,
                          vmin=cmin, vmax=cmax,
                          interpolation=interp,
                          extent=self.extent)

            fig.colorbar(a, ax=ax, label=self.product)

        # if you want to return the fig and axis for the plot
        if returnaxis==True:

            return fig, ax

    # Flatten method for saving with json
    def flatten(self):
        # Return the dictionary object of itself
        return self.__dict__


class Project:

    def __init__(self, projdir, projstr=None,
                 projtags=[], loadfromjson=True,
                 region=[], date=[],time=[],
                 BuildingSize=[],TerrainAngle=[],IterationAngle=[],
                 epsg=None, productdir='/n/davies_lab/Lab/data/products_merged/'):

        self.BuildingSize = BuildingSize
        self.TerrainAngle = TerrainAngle
        self.IterationAngle = IterationAngle

        # if not loading from json, make metadata and merged files
        if loadfromjson == False:

            # Note - the idea would be that you only need projtags and region to grab project files for a collection
            self.projdir = os.path.abspath(projdir)
            self.projtags = projtags
            self.productdir = os.path.abspath(productdir)

            # Create project string from directory path if none given:
            self.projstr = projstr
            if projstr is None:
                self.projstr = self.projdir.split('_')[-1]

            # Eventually, add a way to get the region from the location of plot, using a map of kruger
            self.region = region

            # Get base station data
            self.BaseStationData = BaseStation(projdir=self.projdir)

            # Get Trajectory
            self.Trajectory = Trajectory(projdir=self.projdir)

            # Calc distance to basestation
            # Note: this also sets self.epsg!!!
            self.getdisttostart()

            # Get the start date and time from the trajectory
            self.getdatetime()

            # Get SE version from date
            self.getSEversion()

            # Get Ground Parameters from macro file
            self.getgroundparams()

            # Merge Tiles and get Tile Data (if not already present)
            self.mergetiles()

            # Get Info about the Merged Rasters
            self.gettiles()
            
            # Read accuracy files
            self.readaccuracyfiles()

            # Make extent shape and kml files:
            try:
                self.makeextentshpandkml()
            except:
                print('Unable to make shapefile of extent. \n')

            # Find las tiles
            self.findlasfiles()

            # Save the metadata as a json for another day...
            self.metadatajson()

        # Else, just load from the json file
        else:

            # Get input proj dir
            self.projdir = os.path.abspath(projdir)

            # Get metadatapath
            self.metadatapath = glob.glob(self.projdir + '/docs/*_meta.json')[0]

            # Load attributes up
            self.load()

    def findmacros(self, kind='ground'):

        macrofilelist = []

        folder = self.projdir

        if kind=='ground':
            for name in glob.glob(folder + '/**/*ClassifyGround*.mac', recursive=True):
                if ('byline' in name) == False:
                    macrofilelist.append(name)

        # Add functionality for Noise later (since it's more complicated)
        # if type == 'noise':
        #     for name in glob.glob(indir + '/**/*Noise*.mac', recursive=True):
        #             macrofilelist.append(name)

        self.groundmacrofiles = macrofilelist

    def getgroundparams(self, macrofilepath=None):

        # if not inputting a ground macro path, go find ground macros
        if macrofilepath==None:
            self.findmacros(kind='ground')

        # For loop in case of multiple ground files:
        for macfile in self.groundmacrofiles:
            idx = open(macfile, 'r').read().find('FnScanClassifyGround(')
            paramlist = open(macfile, 'r').read()[(idx + 21):(idx + 21 + 45)].split(',')
            self.BuildingSize.append(float(paramlist[4]))
            self.TerrainAngle.append(float(paramlist[5]))
            self.IterationAngle.append(float(paramlist[6]))

    def mergetiles(self):

        # Check if there is a merged lattice folder already
        if not hasattr(self, 'mergedfolder'):
            # If not, set and make a merged folder for the project in the product directory
            self.mergedfolder = self.productdir + '/' + str(self.date) + '_' + self.projstr
            # Check the path doesn't exist before remaking merged folder
            if not os.path.exists(self.mergedfolder):
                os.mkdir(self.mergedfolder)

        # Start filling merged folder with tiles
        print("Merging Tiles... \n")

        # Find the lattice folder directory in the projdir
        latticefolder = glob.glob(self.projdir + '/**/*attices/', recursive=True)[0]

        # Walk the lattice folder
        latfolderlist = []
        for root, dirs, files in os.walk(latticefolder, topdown=False):
            # For each folder of rasters
            for name in dirs:
                # get current lattice folder
                latfolder = os.path.join(root, name)
                # As long as there's not a merged folder floating around-
                if not (('Merged' in latfolder) | ('merged' in latfolder)):
                    # And as long as it has tif files in it:
                    if glob.glob(latfolder + '/*.tif'):
                        # stick the folder in the list
                        latfolderlist.append(os.path.abspath(latfolder))

        # For each folder of raster products
        for folder in latfolderlist:

            # initialize f for each folder
            f = None

            try:

                # Get the name of the first tif file
                f = os.path.abspath(glob.glob(folder + '/*.tif')[0])

                # Use the name to make the new name of a directory in the merged folder
                prodname = f.split('/')[-1].split('_')[1]
                
            except:
                
                print(f'Issue with tif files in {folder} \n')
            
            try:
                # Make the new directory for this product
                outdir = self.mergedfolder + '/' + prodname
                os.mkdir(self.mergedfolder + '/' + prodname)

            except OSError as exc:
                # The below doesn't seem to work, but it's supposed to
                # properly deal with exception of folder existing
                # https://stackoverflow.com/questions/18973418/os-mkdirpath-returns-oserror-when-directory-does-not-exist
#                 if exc.errno != errno.EEXIST:
#                     raise
                pass

            # Merge the files in each folder by resolution:
            Functions.mergeRasterFiles(indir=folder,
                                       outpath=outdir,
                                       outstr=self.projstr + '_' + prodname + '_',
                                       epsg=self.epsg)

            print(f"{prodname} tiles merged in: {outdir} \n")

        # else:
        #
        #     print(f"Tiles already merged in: {self.mergedfolder}")

    def addprojecttag(self, tag, save=True):

        #Add Tag
        self.projtags.append(tag)

        # Re-dump json file, if save == True
        if save == True:
            self.metadatajson()
        # Note: This might be inefficient, or it could overwrite certain changes
        # May be better to open metadatafile and add to it

    # loading in tiles by filename, v2
    def gettiles(self):

        # Make Product Labels by looking for tifs in Merged Folder
        # initialize list of merged files
        mergedfilepaths = []  
        # initialize list of product names for each file, to be iterated through below
        prodstrs = []  
        
        # for each merged tif file
        for tif in glob.glob(self.mergedfolder + '/**/*.tif', recursive=True):
            # standardize the path str with abspath()
            tif = os.path.abspath(tif)
            # add the file path to the mergedfilepath attribute
            mergedfilepaths.append(tif)
            # get the label of the current product from the file path string, and add it to a list of product names
            prodstrs.append(tif.split('/')[-1].split('_')[-2])

        # Initialize a temporary dict
        tempobjdict = {}

        # Get unique labels of each available merged product
        self.productlabels = list(set(prodstrs))

        # for each product (CHM, DSM, DTM, etc.)
        for prodstr in self.productlabels:

            # initialize a list within the dict, with the product name as the key
            tempobjdict[prodstr] = {}

        # for each file and for the product name of each file
        for f, prodstr in zip(mergedfilepaths, prodstrs):

            # Make a merged tile and pass the project epsg to each tile
            temptile = MergedTile(f, epsg=self.epsg)
            tempobjdict[prodstr][temptile.res] = temptile

        # for each key (product) in the temp dictionary
        for key in tempobjdict:
            # set a new attribute
            self.__setattr__(key, tempobjdict[key])

    # Flatten method for saving with json
    # v3 03/01/2021
    def flatten(self):

        # FLATTEN:
        output_dict = {}

        # for every attribute in the project object
        for att in self.__dict__:

            # if the current attribute contains MergedTile objects
            if att in self.productlabels:

                # initialize a dict within this dict named by the current product (CHM, DSM, etc.)
                output_dict[att] = {}

                # for each resolution of MergedTile object in the dictionary
                for res in getattr(self, att):

                    # print(f"{att, res} is in here.")

                    # unpack the content of the class object and flatten into a dictionary
                    tiledict = self.__getattribute__(att)[res].flatten()

                    # delete unwanted keys
                    # (those made by gdal that shouldn't be saved in json)
                    for i in ['band', 'tif', 'bandarray']:
                        if i in tiledict.keys():
                            tiledict.__delitem__(i)

                    # Store the flat tile in the output dict
                    # with product name then tile resolution as the keys
                    output_dict[att][res] = tiledict

            if att == 'BaseStationData':

                output_dict[att] = self.BaseStationData.flatten()

            if att == 'Trajectory':

                output_dict[att] = self.Trajectory.flatten()

            # else, if the attribute is not an object
            # else: For some reason, else statement does not work on its own! ends up overwriting products
            if not((att in self.productlabels) | (att == 'BaseStationData') | (att == 'LasTile') | (att == 'Trajectory')):

                # Save the value of the attribute in the output dictionary with the attribute name as the key
                output_dict[att] = getattr(self, att)

        # Return the output dictionary
        return output_dict

    # Save as json in docs folder of projdir
    def metadatajson(self):

        # Make a \docs\ sub-directory within the project directory. (if not already made)
        if not os.path.exists(self.projdir + '/docs'):
            os.mkdir(self.projdir + '/docs')

        with open(self.projdir + '/docs/' + self.projstr + "_meta.json", "w", encoding="utf-8") as of:
            # Save dictionary as json, indent=4 adds indentation
            json.dump(self.flatten(), of, indent=4)
            
        # Make a copy of the docs folder with the json in the products folder
        shutil.copytree(self.projdir + '/docs', self.mergedfolder + '/docs', dirs_exist_ok=True, copy_function=shutil.copy2)

    def load(self):

        # Json loading:
        with open(self.metadatapath, "r", encoding="utf-8") as f:

            input_dict = json.load(f)

            # Load from json
            # For each key in the input dictionary
            for key in input_dict:

                # If the key is a product, and therefore a dictionary
                if key in input_dict['productlabels']:

                    # Initialize empty dict
                    rdict = {}

                    # For each resolution of the product
                    for res in input_dict[key]:

                        # initialize a new tile object
                        tile = MergedTile()

                        # fill attributes of that tile object with dict from the loaded json
                        tile.__dict__ = input_dict[key][res]

                        #  fill rdict
                        rdict[res] = tile

                    # Save the dictionary of MergedTile objects back in the Project Object
                    setattr(self, key, rdict)

                if key == 'BaseStationData':

                    # initialize a new BaseStation object
                    # use loadfromjson flag to tell it not to search for functions
                    temp = BaseStation(loadfromjson=True)

                    # fill attributes of that object with dict from the loaded json
                    temp.__dict__ = input_dict[key]
                    
                    # # for each attribute in the BaseStation dict
                    # for att in input_dict[key]:
                    #
                    #     # fill attributes of that tile object with dict from the loaded json
                    #     tile.__dict__ = input_dict[key][res]
                    #
                    #     # fill the station dictionary with current attribute value
                    #     stationdict[att] = inputdict[key][att]

                    # then add the dictionary back into the project object!
                    setattr(self, key, temp)

                if key == 'Trajectory':

                    # initialize a new Trajectory object
                    temp = Trajectory(loadfromjson=True)

                    # fill attributes of that object with dict from the loaded json
                    temp.__dict__ = input_dict[key]

                    # then add the dictionary back into the project object!
                    setattr(self, key, temp)

                # Otherwise if it not is a product, simply load in the attribute
                # Again, strangely, else: doesn't work! Need to explicitly call out scenario
                # else:
                if not((key in input_dict['productlabels']) | (key == 'BaseStationData') | (key == 'Trajectory')):

                    # Transfer value to the attribute
                    setattr(self, key, input_dict[key])

    # updates the project with new MergedTiles and the metadata json with new project info
    def addrastertometa(self, rasterfilepath):

        # Get prodstr from raster filename
        prodstr = rasterfilepath.split('/')[-1].split('_')[-2]

        # Get res from raster filename
        filename = rasterfilepath.split('_')[-1].split('.')[0]
        res = filename.split('_')[-1]

        if prodstr in self.__dict__.keys():

            # deposit current dictionary of product in temp
            temp = self.__getattribute__(prodstr)

            # Add new product string to existing dict
            temp[res] = MergedTile(rasterfilepath)

            # Overwrite old dictionary with new updated dictionary
            self.__setattr__(prodstr, temp)

        # Else, if there isn't an existing dictionary for this additional product
        else:

            # Generate new product string and a new dict, if it doesn't already exist
            self.__setattr__(prodstr, {res:MergedTile(rasterfilepath)})

            # add prodstr to prodlabels
            self.productlabels.append(prodstr)

        # Save the metadata as a json for another day...
        self.metadatajson()

    # method for Calculating Hillshade from DEM, and depositing in merged raster folder
    def hillshade(self, res=['010m', '025m', '050m', '1m'], update=True):
        # update means it adds it to the project and save it to json
        # else, it will just make an instance of it and not save it to json
        # either way, it adds a new tif to the merged folder

        # Commented out, just moved to function call above
        # if res is None:
        #     res = ['010m', '025m', '050m', '1m']

        # If res is just 1 str, make it a list so that the loop works below
        if type(res) is str:
            res = list(res)

        # make an output directory
        outdir = self.mergedfolder + "/Hillshade"

        if not os.path.exists(outdir):
            os.mkdir(outdir)

        for resstr in res:

            # make target outfile
            of = outdir + '/' + self.projstr + '_Hillshade_' + resstr + '.tif'

            # Make Hillshade
            # Use subprocess to call gdal command line:
            # uses default values for the moment, but can choose resolution of product
            # Note: may not work on the server
            # https://joeyklee.github.io/broc-cli-geo/guide/XX_digital_elevation_models.html
            args = ['gdaldem', 'hillshade', self.DTM[resstr].filepath, of]
            subprocess.call(args=args)

            if update==True:

                # Add hillshade to raster tiles list & metadata json
                self.addrastertometa(rasterfilepath=of)

    # New function For project class
    # also calculates the EPSG of the project (so that's automatic now!)
    def getdisttostart(self):

        sbet_geodf = gpd.GeoDataFrame(self.Trajectory.startpositions,
                                      geometry=gpd.points_from_xy(self.Trajectory.startpositions['longitude'],
                                                                  self.Trajectory.startpositions['latitude'],
                                                                  self.Trajectory.startpositions['altitude'],
                                                                  crs="EPSG:4326"))
        # determine utm zone from WGS84 to with function (from Functions script)
        # NOTE: may want to use this in your project class!
        self.epsg = Functions.getUTMcodefromWGS84(lon=self.Trajectory.startpositions['longitude'][0],
                                                 lat=self.Trajectory.startpositions['latitude'][0])

        # Reproject to UTM
        sbet_geodf_utm = sbet_geodf.to_crs(f"EPSG:{self.epsg}")

        # now, get the base station coords
        bs_geodf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(x=[self.BaseStationData.lon],
                                                                y=[self.BaseStationData.lat],
                                                                z=[self.BaseStationData.alt],
                                                                crs=self.BaseStationData.epsg))

        # Reproject to UTM
        bs_geodf_utm = bs_geodf.to_crs(f"EPSG:{self.epsg}")

        # Now that they're both in UTM, calc the distance from the BS to the traj start
        self.DistancetoBaseStat = np.sqrt((sbet_geodf_utm.geometry.x[0] - bs_geodf_utm.geometry.x[0]) ** 2 +
                                          (sbet_geodf_utm.geometry.y[0] - bs_geodf_utm.geometry.y[0]) ** 2)

    def getdatetime(self):

        # Get the first time position in GPS TOW from the trajectory
        start_GPSTOW = np.min(self.Trajectory.startpositions['time'])

        # Get the current date from the base station
        curr_date = self.BaseStationData.startdatetime[0:10]

        # Make a basestation datetime, and use it to count back to the beginning of the week (don't actually need the time here)
        BS_datetime = datetime(year=int(curr_date[-4:]),
                               day=int(curr_date[3:5]),
                               month=int(curr_date[0:2]))

            # Start the time of that week (since GPS TOW counts seconds since the start of the week)
            # Subtract number of days  isoweekday() of base station date to count back to the beg. of week
        start_of_week = datetime(year=int(curr_date[-4:]),
                                 day=int(curr_date[3:5]),
                                 month=int(curr_date[0:2])) - timedelta(days=BS_datetime.isoweekday())

        # Use datetime's timedelta to add the GPS TOW seconds to the start of the week datetime
        curr_datetime = start_of_week + timedelta(seconds=start_GPSTOW)

        # And also, get the time of the end of the flight, using the duration from the trajactory
        end_datetime = curr_datetime + timedelta(seconds=self.Trajectory.duration_sec)

        # Set the current date and time! Easy.
        self.date = curr_datetime.strftime("%Y-%m-%d")
        self.time = curr_datetime.strftime("%H:%M:%S")
        self.enddate = end_datetime.strftime("%Y-%m-%d")
        self.endtime = end_datetime.strftime("%H:%M:%S")

    def getSEversion(self):

        timeofswitchtoSEv6 = '2020-03-26'

        if self.date >= timeofswitchtoSEv6:
            self.SEversion = "6.0.3"
        if self.date < timeofswitchtoSEv6:
            self.SEversion = "5.0.4"

    def makeextentshpandkml(self):

        # read in DTM for full extent of data
        DTM = rio.open(self.DTM['1m'].filepath)

        # Resample DTM to 10 m pixels and save it (to speed up below process)
        scale = 10
        outf = self.projdir + f'/docs/{self.projstr}_DTM_{scale}m_rescaled.tif'
        Functions.resample_raster(DTM, scale=scale, outpath=outf)

        # Open the resampled DTM
        DTM_scaled = rio.open(outf)

        # Extract Shapes from Raster features
        # https://rasterio.readthedocs.io/en/latest/topics/features.html?highlight=polygon#extracting-shapes-of-raster-features
        # Note: this shape is a Generator object: https://realpython.com/introduction-to-python-generators/
        # https://gis.stackexchange.com/questions/187877/how-to-polygonize-raster-to-shapely-polygons
        # note: != True flips the True values to False (it masks with True, but rio's "shapes" expects False)
        # This is a generator object that you can iterate through
        DTM_mask_shapes = shapes(DTM_scaled.read(1),
                                 mask=DTM_scaled.read(1, masked=True).mask != True,
                                 transform=DTM_scaled.transform)

        # Use Geojson to make each pixel into a polygon (polygonize)
        results_poly = []
        for i in DTM_mask_shapes.__iter__():
            results_poly.append(geojson.Polygon(coordinates=i[0]['coordinates']))

        # Use geopandas to dissolve them all into one polygon
        # Group is just a unique index for each pixel
        DTM_poly_gdf = gpd.GeoDataFrame({'group': np.repeat(1, len(results_poly))},
                                        geometry=results_poly)
        DTM_poly_gdf_dissolve = DTM_poly_gdf.dissolve(by='group')

        # Plot it to check it:
        # DTM_poly_gdf_dissolve.plot()

        # Set the CRS
        DTM_poly_gdf_dissolve = DTM_poly_gdf_dissolve.set_crs(f'EPSG:{self.epsg}')

        # Now export it to a shp file
        DTM_poly_gdf_dissolve.to_file(self.projdir + f'/docs/{self.projstr}_extent')

        # Export as KML
        # make kmls read/writable in geopandas
        gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
        gpd.io.file.fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'
        try:
            DTM_poly_gdf_dissolve.to_file(self.projdir + f'/docs/{self.projstr}_extent.kml', driver='KML')
        except:
            pass

    def findlasfiles(self):

        try:
            # find output las file directory
            self.lasdir = glob.glob(self.projdir + '/**/*loud*/', recursive=True)[0]

            # list alll lasfiles within it
            self.lasfiles = glob.glob(self.lasdir + '/*.las')
        
        except:
            
            print('Issue finding las directory or files \n')
                  
    def readaccuracyfiles(self):
        
        try: 
            
            # get paths of all report files
            accfiles = glob.glob(self.projdir + '/**/*eport*/*.txt', recursive=True)
            
            # for each accuracy file
            for accf in accfiles:
                
                # grab the name for below
                fname = accf.split('/')[-1]
                
                # open it and read the lines
                with open(accf, 'r') as f:
                    lines = f.readlines()
                    # if it's the initial accuracy, record the mean
                    if ('initial' in fname) | ('00' in fname):
                        self.initialaccuracy = float(lines[3][-8:-1])
                    
                    # if it's the last, record the mean
                    if ('fluct' in fname) | ('inal' in fname):
                        self.finalaccuracy = float(lines[3][-8:-1])
                
        except:
            
            print('Issue finding or reading accuracy reports. \n')
            
    

# Collection Class
class Collection:

    def __init__(self,
                 cdir=[],
                 searchby=None,
                 pdirs=[],
                 ctags=[],
                 reinitialize=False,
                 shapes=[]):

        # Add option to searchby ='tags'

        # self.indir = os.path.abspath(indir)
        self.cdir = os.path.abspath(cdir)
        self.pdirs = pdirs

        # if there's only 1 ctag as a str,
        # then make sure to turn it into a list
        if type(ctags) is str:
            ctags = [ctags]

        self.ctags = ctags
        self.searchby = searchby
        self.shapes = shapes

        # Flag to force it to initialize projects or not
        self.reinitialize = reinitialize

        # Initiate Dict of Projects
        self.projects = {}
        # initiate list of project names
        self.projstrs = []
        # Get projects
        self.initprojects()
        # Make accuracy table
        self.makeaccuracytable()

    def initprojects(self):

        # if there was a directory passed for the collection
        if self.cdir:

            # Get the project folders of all projs in cdirs path
            # works by searching for the Terrasolid dgn file
            # Note: Need to make this more sophisticated in the future
            # for now, it just grabs all folders
            pdirs = glob.glob(self.cdir + '/**/*.dgn', recursive=True)

            # For each project in the directory of projects
            for p in pdirs:

                # Get the path of the project
                projdir = '/'.join(p.split('/')[0:-2])

                # Commented out os.remove - because you don't actually need it!
                # And it messes up cases where you have some projects have jsons, while others do not.
                # If reinitialization == True
                # (ie: need to force a new metadata file creation with new paths)
                # delete the existing metadata file
                if self.reinitialize:
                    metafilepath = glob.glob(projdir + '/**/*_meta.json', recursive=True)
                    try:
                        os.remove(metafilepath[0])
                        print(f'Overwriting {metafilepath[0]}... \n')
                    except:
                        print('No json file deleted. \n')

                    # One way to do it: ask for input, but this would be annoying when collecting on a large scale
                    # input = Input(f'Overwrite json in {projdir}? Y or N')
                    # if (overwritejson == 'Y') | (overwritejson == 'y'):
                    #     os.remove(metafilepath[0])

                # if a metadata file already exists
                if glob.glob(projdir + '/**/*_meta.json', recursive=True):

                    print(f'Loading project in: {projdir} \n')

                    # then, load the project from the file
                    proj = Project(projdir, loadfromjson=True)

                    # And append the collection tag to the project (if it doesn't already exist)
                    for t in self.ctags:
                        if not t in proj.projtags:
                            proj.addprojecttag(t)

                # if the Project is new (no metadata file)
                if not glob.glob(projdir + '/**/*_meta.json', recursive=True):

                    print(f'Initiating project in: {projdir} \n')

                    # then, create a project and include the collection tags
                    proj = Project(projdir, projtags=self.ctags, loadfromjson=False)

                # Finally, add the project dir to the pdirs att.
                # and add the project to the projects dictionary
                self.projects[proj.projstr] = proj
                self.pdirs.append(projdir)
                self.projstrs.append(proj.projstr)

        # TBD
        # If no cdirs were passed to the argument
        if self.searchby == 'tags':

            # Search tags and populate the project dictionaries
            self.searchtags()

            
    def searchtags(self):

        # TBD
        print(self.searchtags)

        
    def statstable(self, shape=None,
                   res=['1m', '010m'],
                   productlist=['CHM', 'DTM', 'DSM', 'AvgInt1Ret', 'PCount', 'GPCount','PulseCount'],
                   saveclippedband=False,
                   writeoutclippedrasters=False,
                   odir = '/n/davies_lab/Users/pbb/collections/test/rasters/'):
        
        # if input is a str, make it a list
        if type(res) is str: 
            res = list(res)
        
        # If no input shape to clip rasters:
        if shape==None:
        
            # Initialize Dict
            stats_dict = {'site':[],
                          'product':[],
                          'res':[],
                          'min':[],
                          'max':[],
                          'avg':[],
                          'std':[],
                          'datafootprint_msq':[],
                          'datafootprint_ha':[]}

            # for each project in the collection
            # Note: try except loops protect against missing attributes at each res
            for key in self.projects.keys():

                    proj = self.projects[key]

                    # for each resolution chosen
                    for r in res:

                        # for each of these products
                        for prod in productlist:
                            
                            try:
                                
                                # get the product for the current project
                                projprod = proj.__getattribute__(prod)
                                
                            except:
                                
                                print(f'Issue loading {proj.projstr} {prod} \n')

                            # Then, loop through the stats and record them
                            # Note: outer try loop protects against missing products at a given resolution
                            # Inner try prevents against missing attributes
                            try:
                                for attr in stats_dict:
#                                 for attr in projprod[r].__dict__:
                                    if attr != 'sensor':
                                        try:
                                            stats_dict[attr].append(projprod[r].__getattribute__(attr))
                                        except:
                                            stats_dict[attr].append(np.nan)
                            except:
                                pass
       
            # Save dataframe to collection object  
            self.statsbysite = pd.DataFrame(stats_dict)

            # delete some unecessary columns
#             self.statsbysite.drop(columns=['epsg','filepath','dateprocess','datafootprint_msq_sum'], inplace=True)
        
        # Else, Use a shape to clip all rasters and re-calculate stats
        # note: shape is assumed to be a geojson polygon file
        else:
            
            # if you want to save the clipped band pixel data below
            if saveclippedband == True:
                
                # initialize another dict
                self.clippedband_dict = {'site':[],
                                         'product':[],
                                         'res':[],
                                         'data':[],
                                         'transform':[],
                                         'meta':[]}
                
                # if you are writing out the rasters as saved tifs
                if writeoutclippedrasters==True:
                    
                    # add filepath attribute
                    self.clippedband_dict['filepath'] = []
                                        
                
            # Add the shape to the collection class for safekeeping
            shape = os.path.abspath(shape)
            self.shapes.append(shape)
            
            print(f'Loading {shape}... \n')

            with open(shape, 'r') as f:
                poly = geojson.load(f)
            
            # Now Clip raster and Calculate Stats
            
            # Initialize Dict
            stats_dict = {'site':[],
                          'product':[],
                          'res':[],
                          'min':[],
                          'max':[],
                          'avg':[],
                          'std':[],
                          'datafootprint_msq':[],
                          'datafootprint_ha':[]}
            
            # for each project in the collection
            # Note: try except loops protect against missing attributes at each res
            for key in self.projects.keys():

                    proj = self.projects[key]      
                    print(f'Calculating stats for {proj.projstr} \n')

                    # for each resolution chosen
                    for r in res:

                        # for each of these products
                        for prod in productlist:

                            try:
                                # get the product for the current project
                                projprod = proj.__getattribute__(prod)
                                
                            except:
                                
                                print(f'Issue loading {proj.projstr} {prod} \n')
                                
                            try:

                                # Load it with rio
                                tif = rio.open(projprod[r].filepath)

                                # clip raster with polygon
                                # Written with help from: https://rasterio.readthedocs.io/en/latest/topics/masking-by-shapefile.html
                                tif_clip, tif_transform = riomask(dataset=tif,
                                                                  shapes=[poly],
                                                                  crop=True)

                                # Re-mask array
                                bandarray = np.ma.masked_where(tif_clip < 0, tif_clip)

                                # Fill out stats_dict table
                                # Could fill with try except loops...
                                stats_dict['site'].append(proj.projstr)
                                stats_dict['product'].append(prod)
                                stats_dict['res'].append(r)
                                stats_dict['avg'].append(np.ma.mean(bandarray))
                                stats_dict['std'].append(np.ma.std(bandarray))
                                stats_dict['min'].append(np.ma.min(bandarray))
                                stats_dict['max'].append(np.ma.max(bandarray))
                                
                                # Calc data footprint (number of data pixels * area of pixel)
                                # Note: Only works with UTM
                                pixelsize = np.abs(tif.bounds.right - tif.bounds.left) / tif.width
                                stats_dict['datafootprint_msq'].append(np.ma.count(bandarray)*(pixelsize**2))
                                stats_dict['datafootprint_ha'].append( (np.ma.count(bandarray)*(pixelsize**2)) / 10000)
                                
#                               Add option to save band array data
                                if saveclippedband == True:
        
                                    self.clippedband_dict['site'].append(proj.projstr)
                                    self.clippedband_dict['product'].append(prod)
                                    self.clippedband_dict['res'].append(r)
                                    self.clippedband_dict['data'].append(tif_clip)
                                    self.clippedband_dict['transform'].append(tif_transform)
                            
                                    # Update meta info before saving
                                    meta = tif.meta
                                    
                                    meta.update({"height": bandarray.shape[1],
                                                 "width": bandarray.shape[2],
                                                 "transform": tif_transform})
                                    
                                    self.clippedband_dict['meta'].append(meta)
                                
                                    if writeoutclippedrasters == True:
                                        
                                        # Write out clipped rasters to your lab user folder:
                                        # Note: odir is by default the directory in the function
                                        # some day, put it all in an organized collection class
                                        self.clippedband_dict['filepath'].append(odir + f'{proj.projstr}_{prod}_{r}_clip.tif')

                                        with rio.open(odir + f'{proj.projstr}_{prod}_{r}_clip.tif', "w", **meta) as dest:
                                            dest.write(tif_clip)

                            except Exception as e:
                                print(e.with_traceback, e.args)
                                print(f'Issue with {proj.projstr} {prod} at {r} \n')
                                
            
            # Save dataframe to collection
#             self.stats_dict = stats_dict
            self.statsbysite = pd.DataFrame(stats_dict)
            print('Done calculating stats.')
            
            
    def densityplot(self, res=['1m'], kind='CHM',
                    subsample=True, sample='Block', save=False, returnfig=False):

        if type(res) is str:
            res =[res]

        # for each resolution selected
        # (1 plot per resolution)
        for r in res:

            # initiate arrays for matplotlib and seaborne
            arraylist = []
            labels = []

            # for each project in the collection
            for key in self.projects.keys():

                proj = self.projects[key]

                if kind == 'CHM':
                    # open up the raster with Rio
                    tif = rio.open(proj.CHM[r].filepath)
                    band = tif.read(1)
                    # Re-mask np array:
                    bandarray = np.ma.masked_where(band < 0, band)

                elif kind == 'PCount':
                    # open up the raster with Rio
                    tif = rio.open(proj.PCount[r].filepath)
                    band = tif.read(1)
                    # Re-mask np array:
                    bandarray = np.ma.masked_where(band <= 0, band)

                elif kind == 'GPCount':
                    # open up the raster with Rio
                    tif = rio.open(proj.GPCount[r].filepath)
                    band = tif.read(1)
                    # Re-mask np array:
                    bandarray = np.ma.masked_where(band <= 0, band)

                else:

                    print(f'Unrecognized kind of raster: {kind}. \n Select from: \"CHM\",\"PCount\",\"GPCount\" \n')

                    continue

                # If you want a subsample to plot, and not the whole raster
                if subsample == True:

                    # if you want a random subsample (like 50% of pixels) to plot instead
                    if sample == 'Random':
                        bandarray = bandarray.flatten()[bandarray.flatten() > 0]
                        bandarray = np.random.choice(bandarray,
                                                    int(0.50 * np.size(bandarray)),
                                                    replace=False)  # grab random pixels from the array (50% of size)
                        # Re-mask np array:
                        # bandarray = np.ma.masked_where(bandarray < 0, bandarray)

                    # if you want to plot a square section (centered on the raster)
                    if sample == 'Block':

                        # This method seems to be closest to center (even though it's not quite...)
                        bounds = list(tif.bounds._asdict().values())
                        x = np.mean([bounds[0], bounds[2]])
                        y = np.mean([bounds[1], bounds[3]])

                        # Make a geojson polygon of the 300x300m block in the center of the raster
                        # Geojson formatting: https://python-geojson.readthedocs.io/en/latest/#coords
                        # NOTE: This assumes a UTM projection
                        d = 300
                        # d = block dim size

                        block_poly = geojson.Polygon([[(x - d/2, y + d/2),
                                                      (x + d/2, y + d/2),
                                                      (x + d/2, y - d/2),
                                                      (x - d/2, y - d/2),
                                                      (x - d/2, y + d/2)]])

                        # NOTE: you probably want to use the same polygon for BuffCamp Project
                        # so maybe save this and use it for later:
                        # add a "shape" feature in the collection class
                        # could also add a clipped raster folder
                        # with open('./test_outputs/geojson_block.json', 'w') as of:
                        #     geojson.dump(block_poly, of)

                        # clip raster with block polygon
                        tif_clip = riomask(dataset=tif,
                                             shapes=[block_poly],
                                             crop=True)

                        # re-open up the raster with Rio
                        band_clip = tif_clip[0]
                        
                        if 'Count' in kind:
                            # Re-mask array
                            bandarray = np.ma.masked_where(band_clip <= 0, band_clip)
                            # Multiply by 1000 for pcount since values are in thousands
                            bandarray = bandarray*1000
                            
                        else:
                            bandarray = np.ma.masked_where(band_clip < 0, band_clip)

                    else:

                        print(f'Unrecognized subsample: {sample}. \n Select from: \"Random\" or \"Block\" \n')

                # then, add these values to the arraylist
                arraylist.append(bandarray.flatten())

                # And append the label and resolution
                labels.append(proj.projstr + '_' + r)

        # Make a seaborn density plot
        # Nice tutorial on how to do this: https://towardsdatascience.com/histograms-and-density-plots-in-python-f6bda88f5ac0
        # Iterate through the arraylist
        fig1, ax1 = plt.subplots()
        for array, i in zip(arraylist, range(len(arraylist))):

            # Draw the density plot
            sns.kdeplot(array, fill=True, ax=ax1)

        # Plot formatting
        ax1.legend(labels)
        ax1.set_ylabel('Density')
        ax1.grid()
        plt.title(kind + ' ' + r)
        
        if 'Count' in kind:
            ax1.set_xlabel('Point Count per Pixel')
        else:
            ax1.set_xlabel('Height Above Ground [m]')
        
        if returnfig == True:
            return fig1, ax1
        
        if save == True:

            outdir = input('Output directory for saving figs: ')

            try:
                outdir = os.path.abspath(outdir)
            except:
                print(f'{outdir} is not a valid path. \n')

            fig1.savefig(outdir + '/' + kind + '_DensityPlot.png', dpi=300)
            
    
    def makeaccuracytable(self):
        
        self.accuracytable = {'projstr':[],
                              'initial':[],
                              'final':[],
                              'diff':[]}
        
        for p in self.projects:
            
            self.accuracytable['projstr'].append(self.projects[p].projstr)
            
            try:
                self.accuracytable['final'].append(self.projects[p].finalaccuracy)
            except:
                self.accuracytable['final'].append(np.nan)
                
            try:
                self.accuracytable['initial'].append(self.projects[p].initialaccuracy)
            except:
                self.accuracytable['initial'].append(np.nan)
            
            try:
                self.accuracytable['diff'].append(self.projects[p].initialaccuracy - self.projects[p].finalaccuracy)
            except:
                self.accuracytable['diff'].append(np.nan)
            
            
            
class XSection:

    def __init__(self,
                 proj=None,
                 centerxy=None,
                 xsize=0.5,
                 ysize=30,
                 res=0.5,
                 points=[],
                 quantiles=[0.25, 0.5, 0.75, 1],
                 outdir=None):

        # Initialize attributes
        self.proj = proj
        self.centerxy = centerxy
        self.xsize = xsize
        self.ysize = ysize
        self.points = points
        self.res = res
        self.quantiles = quantiles
        self.outdir = outdir

        # Set outdirectory for figures (if you use save below)
        # Note: the plot function asks for this if it's empty and you try to save)
        self.outdir = outdir

        # Check if it's a project object being fed in
        if isinstance(self.proj, Classes.Project):
            # if so, grab some useful attributes
            self.projstr = self.proj.projstr
            self.lasdir = self.proj.lasdir

            # If no center for the XS was defined
            if centerxy is None:
                # set it as the center of the given project
                self.centerxy = proj.DTM['010m'].centercoords

        # if the folder of las directories exists
        if not self.lasdir is None:
            # make it an absolute path
            self.lasdir = os.path.abspath(self.lasdir)
            # and go fetch the las tiles you need
            self.getlastiles()
            # and go fetch the points you need from each tile
            self.getpoints()
            # Compute height above ground
            self.normheights()
            # compute percentile heights
            self.heightpercentiles()

    def getlastiles(self):

        # initialize list of lasfiles
        self.lasfiles = []
        self.headers = []
        # for each las file in the directory
        for f in glob.glob(self.lasdir + '/*.las'):
            # open the las file for reading
            with laspy.open(f) as l:
                
                # if the XS is long in the x direction
                if self.xsize >= self.ysize:
                    
                    # use the header and the min/max of the Xsection to check if the cross section falls in that tile
                    if ((((self.centerxy[0] + self.xsize/2) < l.header.maxs[0]) & ((self.centerxy[0] + self.xsize/2) > l.header.mins[0]) &
                          (self.centerxy[1] < l.header.maxs[1]) & (self.centerxy[1] > l.header.mins[1]) ) |
                        (((self.centerxy[0] - self.xsize/2) < l.header.maxs[0]) & ((self.centerxy[0] - self.xsize/2) > l.header.mins[0]) &
                            ((self.centerxy[1]) < l.header.maxs[1]) & (self.centerxy[1] > l.header.mins[1]))):
                        self.lasfiles.append(f)
                        self.headers.append(l.header)
                        
                else:
                    
                    # use the header and the min/max of the Xsection to check if the cross section falls in that tile
                    if ((((self.centerxy[1] + self.ysize/2) < l.header.maxs[1]) & ((self.centerxy[1] + self.ysize/2) > l.header.mins[1]) &
                          (self.centerxy[0] < l.header.maxs[0]) & (self.centerxy[0] > l.header.mins[0]) ) |
                        (((self.centerxy[1] - self.ysize/2) < l.header.maxs[1]) & ((self.centerxy[1] - self.ysize/2) > l.header.mins[1]) &
                            ((self.centerxy[0]) < l.header.maxs[0]) & (self.centerxy[0] > l.header.mins[0]))):
                        self.lasfiles.append(f)
                        self.headers.append(l.header)
                        

    def getpoints(self):
        
        self.las_x_scaled = []
        self.las_y_scaled = []
        self.las_z_scaled = []
        self.classification = []
        
        # open up the las files
        for f in self.lasfiles:

            las = laspy.read(f)

            # subset the points to only default and ground points within the cross section
            # NOTE: 
            self.points = las.points[(las.x <= self.centerxy[0] + self.xsize/2) &
                                     (las.x >= self.centerxy[0] - self.xsize/2) &
                                     (las.y <= self.centerxy[1] + self.ysize/2) &
                                     (las.y >= self.centerxy[1] - self.ysize/2) &
                                     (las.classification != 7)]

             # Scale x, y, z and ground points
            self.las_x_scaled.append(self.points.array['X'] * self.headers[0].x_scale + self.headers[0].x_offset)
            self.las_y_scaled.append(self.points.array['Y'] * self.headers[0].y_scale + self.headers[0].y_offset)
            self.las_z_scaled.append(self.points.array['Z'] * self.headers[0].z_scale + self.headers[0].z_offset)

            self.classification.append(self.points.classification)
            
        # Concatenate (in the case of multiple files being selected)
        if len(self.classification) > 1:
            self.classification = np.concatenate(self.classification)
            self.las_x_scaled = np.concatenate(self.las_x_scaled)
            self.las_y_scaled = np.concatenate(self.las_y_scaled)
            self.las_z_scaled = np.concatenate(self.las_z_scaled)
        else:
            #otherwise, just pop them out of the list (they are numpy arrays)
            self.classification = self.classification[0]
            self.las_x_scaled = self.las_x_scaled[0]
            self.las_y_scaled = self.las_y_scaled[0]
            self.las_z_scaled = self.las_z_scaled[0]
            
        
        # Make a set of ground points
        self.ground_x_scaled = self.las_x_scaled[self.classification==2]
        self.ground_y_scaled = self.las_y_scaled[self.classification==2]
        self.ground_z_scaled = self.las_z_scaled[self.classification==2]

    def plotpoints(self, norm=False, returnfigandax=True, topdown=False, color=None, clabel='Class', fig=None, ax=None, colorbar=True):
        
        # if no figure is provided, make one
        if not fig:
            fig, ax = plt.subplots()

        if color is None:
            color = self.classification

        if topdown==True:
            a = ax.scatter(x=self.las_x_scaled, y=self.las_y_scaled,
                           c=color,
                           s=2, cmap='viridis', alpha=0.6)
            ax.axis('scaled')
            if colorbar==True:
                fig.colorbar(a, ax=ax, label=clabel)
            ax.set_xlabel('Easting [m]')
            ax.set_ylabel('Northing [m]')

        else:
            if norm==True:
                try:
                    z = self.las_z_norm
                except:
                    self.normheights()
                    z = self.las_z_norm

            else:
                z = self.las_z_scaled

            if self.xsize >= self.ysize:
                a = ax.scatter(self.las_x_scaled, z, s=2, c=color, cmap='viridis', alpha=0.6)
                if colorbar==True:
                    fig.colorbar(a, ax=ax, label=clabel)
                ax.set_xlabel('Easting [m]')
                ax.axis('scaled')
                ax.set(ylim=(-0.2, np.max(z)+1))

            else:
                a = ax.scatter(self.las_y_scaled, z, s=2, c=color, cmap='viridis', alpha=0.6)
                if colorbar==True:
                    fig.colorbar(a, ax=ax, label=clabel)
                ax.set_xlabel('Northing [m]')
                ax.axis('equal')
                ax.set(ylim=(-0.2, np.max(z)+1))

        if returnfigandax==True:
            return fig, ax

        # plt.tight_layout()

    def makeshape(self):
        print("TBD \n")

        # Make a geojson polygon to clip a raster
        # Geojson formatting: https://python-geojson.readthedocs.io/en/latest/#coords
        # NOTE: This assumes a UTM projection
        # self.shape_json = geojson.Polygon([[(self.centerxy[0] + self.xsize/2, self.centerxy[1] + self.ysize/2),
        #                                    self.centerxy[0] - self.xsize/2,
        #                                    self.centerxy[1] + self.ysize/2,
        #                                    self.centerxy[1] - self.ysize/2]])

    def normheights(self, method='DTM'):
        # print("TBD")

        # using nearest:
        if method=='nearest':

            # use nearest ground point to norm heights
            self.las_z_norm = []
            # for each point
            for x, y, z in zip(self.las_x_scaled,
                               self.las_y_scaled,
                               self.las_z_scaled):

                # find the nearest ground point
                self.diff = (x - self.ground_x_scaled)**2 + \
                            (y - self.ground_y_scaled)**2 + \
                            (z - self.ground_z_scaled)**2

                idxmin = np.argmin(self.diff)

                # subtract the current ground z value
                self.las_z_norm.append(z - self.las_z_scaled[idxmin])

        # using DTM
        if method == 'DTM':

            # load the 10 cm DTM (replace with project points)
            DTM = rio.open(self.proj.DTM['010m'].filepath,
                           masked=True)

            # make a coordinate list of tuples from x and y
            # https: // geopandas.readthedocs.io / en / latest / gallery / geopandas_rasterio_sample.html
            coord_list = [(x, y) for x, y in zip(self.las_x_scaled, self.las_y_scaled)]

            # NOTE: this is a generator object, which prevents you from loading everything into memory
            # if you really want to use this right, iterate through it in the loop below instead
            # https://realpython.com/introduction-to-python-generators/
            DTM_elev = rio.sample.sample_gen(DTM,
                                             xy=coord_list,
                                             indexes=1,
                                             masked=False)

            # initialize list of heights
            self.las_z_norm = []

            # for each terrain elevation & point elevation
            for e, z, c in zip(DTM_elev, self.las_z_scaled, self.classification):
                # if it's a ground point
                if c == 2:
                    # set it to height 0
                    self.las_z_norm.append(0)

                else:
                    # subtract the point's elevation from the terrain to get height above ground
                    self.las_z_norm.append(z - e[0])

    def heightpercentiles(self):
        # print('TBD')

        # normalize z to height if not already done
        if not hasattr(self, 'las_z_norm'):
            self.normheights()

        # make a dataframe
        df = pd.DataFrame(data={'x': self.las_x_scaled,
                                'y': self.las_y_scaled,
                                'z': self.las_z_norm})

        # if the transect goes along the x (Easting) direction
        if self.xsize >= self.ysize:
            # use x to bin
            var = df.x

        else:
            # use y
            var = df.y

        # Bin by x coordinate
        # https://stackoverflow.com/questions/16947336/binning-a-dataframe-in-pandas-in-python/
        # generate correct number of bins for given res
        self.nbins = int(np.ceil((var.max() - var.min()) / self.res) + 1)

        # make bin edges
        self.binedges = np.linspace(var.min(), var.max(), num=self.nbins)

        # also, make bin centers (for plotting)
        self.bincenters = self.binedges + self.res/2

        # group by bin
        self.df_group = df.groupby(np.digitize(var, self.binedges))

        # Calc percentiles for each bin, and use "unstack" to flip output
        self.percentiles = pd.DataFrame(self.df_group.quantile(self.quantiles).z).unstack().to_dict('list')

        # Check for empty bins (empty areas without points in them)
        # will show if the total number of groups (number of group keys)
        # is less than the max value of the group ids
        if len(self.df_group.groups.keys()) != list(self.df_group.groups.keys())[-1]:

            print('Found empty bins. Filling... \n')

            # Make group id col for percentile dict
            bin_id = list(self.df_group.groups.keys())

            # for each potential group number
            for g in list(np.arange(1, self.nbins)):

                # check if the group is included, and if not
                if g not in bin_id:

                    # insert a group with 0s for each col in percentile dict
                    for k in self.percentiles.keys():
                        # insert a 0 value at the index of the bin
                        # Note: the bin index is the bin number - 1 (since python is 0 indexed)
                        self.percentiles[k].insert(g - 1, 0)


    def plotpercentiles(self, plotpoints=True, returnfigandax=True, fill=False, savefig=False, fig=None, ax=None, colorbar=True, fillcolor=None):

        if plotpoints == True:
            # if no figure provided
            if not fig:
                fig, ax = self.plotpoints(norm=True,
                                          topdown=False,
                                          color=self.las_z_norm,
                                          clabel='Height [m]',
                                          colorbar=colorbar)
            else:
                self.plotpoints(norm=True,
                                topdown=False,
                                color=self.las_z_norm,
                                clabel='Height [m]',
                                fig=fig, ax=ax,
                                colorbar=colorbar)
        else:
            # else, if there's not a fig provided, make one
            if not fig:
                fig, ax = plt.subplots()

        if fill == False:

            # loop through percentiles and plot each column
            for col in self.percentiles.keys():
                b = ax.plot(self.bincenters,
                            list(self.percentiles[col]),
                            '--',
                            label=f'RH{str(int(col[1]*100))}',
                            alpha=0.8)
                ax.legend()
        else:

            # convert into an array
            self.perc_array = np.array(list(self.percentiles.values()))

            # Make square bin edges (binedges at both sides of the res) for plotting
            self.binedges_square = []
            for left, right in zip(self.binedges, self.binedges + self.res):
                self.binedges_square.append(left)
                self.binedges_square.append(right)

            # loop through percentiles and plot each column
            for col in np.arange(0, self.perc_array.shape[0]):
                # If it's the first column of RHs
                if col == 0:
                    # Set the bottom of the plot fill to be 0 m (ground)
                    minvals = np.zeros(self.perc_array.shape[1])
                else:
                    # Else, use the next lowest RH as the min fill values
                    minvals = self.perc_array[col - 1]

                maxvals = list(self.perc_array[col])

                # duplicate values to match square bins
                minvals_square = []
                maxvals_square = []
                for minv, maxv in zip(minvals, maxvals):
                    # append each value twice to duplicate
                    minvals_square.append(minv)
                    minvals_square.append(minv)
                    maxvals_square.append(maxv)
                    maxvals_square.append(maxv)

                # Line plot (innacurate, but looks good)
                # b = ax.fill_between(self.bincenters,
                #                     minvals,
                #                     maxvals,
                #                     '--',
                #                     label=f'RH{str(int(list(self.percentiles.keys())[col][1] * 100))}',
                #                     alpha=0.3)
                
                # if there's no fillcolor, make it blue
                if not fillcolor:
                    fillcolor='C0'
                
                # Bar plot (accuracte), but looks less good)
                b = ax.fill_between(x=np.array(self.binedges_square),
                                    y2=np.array(minvals_square),
                                    y1=np.array(maxvals_square),
                                    label=f'RH{str(int(list(self.percentiles.keys())[col][1] * 100))}',
                                    alpha=0.25, color=fillcolor)
                ax.legend()

        # Print to check sizes
         # print(f'maxvals: {len(maxvals)} \n minvals: {len(minvals)} \n bincenters {len(self.bincenters)}')

        if savefig==True:
            # Check if there's an outdirectory, and ask user for input if not
            if not self.outdir:
                self.outdir = input('Specify directory path for saving figures: ')
                self.outdir = os.path.abspath(self.outdir)

            fig.savefig(self.outdir + f'/{self.proj.projstr}_{self.xsize}x{self.ysize}m_{self.res}res.png', dpi=300)

        if returnfigandax==True:
            return fig, ax