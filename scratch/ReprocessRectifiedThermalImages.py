import numpy as np
import xarray as xr
import sys
sys.path.append("C://Users//peter//Documents//python//Metadata_DaviesLab_server//")
from Classes import Trajectory, Project
from matplotlib import pyplot as plt
import rioxarray as rioxr
from rasterio.enums import Resampling
from pathlib import Path
import time

# # # START USER INPUTS

# set project directory
projdir = Path("C:/Users/peter/DaviesLab/data/UAV/2020/Kruger_Jan/20200124-081112_MakhohlolaEP")

# rectified image directory (as a Path object)
imgdir_rect = Path(f'{projdir}/TerraSolid/Terraphoto/thermal/rect')

# Set orientation file location (created in terrasolid, with tphoto export orientation)
imgorientf = Path(f"{projdir}/TerraSolid/Terraphoto/thermal/MakhEP_thermal_ImageOrientation.txt")

# Set air temp of site (in celcius)
curr_airtemp = 22
# Set relative humidity (from 0-100%)
curr_relhumid = 0

# # # END USER INPUTS

# # # DEFINE FUNCTIONS

# Function for calculating temp from signal (note in kelvin)
def tempcalc(signal, R=366545, B=1428, F=1, O=342):
    temp = (B / np.log((R / (signal - O)) + F))
    return temp

# Calculating signal (phi or energy)
def phicalc(temp_kelvin, R=366545, B=1428, F=1, O=342):
    phi = (R / (np.exp(B / temp_kelvin) - F)) + O
    return phi

#Function to calculate distance raster
def dist_calc(elev_xr, xcoord, ycoord, Xcam, Ycam, Zcam):
    dist = np.sqrt(((Xcam - xcoord)**2 + (Ycam - ycoord)**2 + (Zcam - elev_xr)**2))
    return dist

# # # END FUNCTIONS

# # # LOAD INPUTS

# Open image orientation file, where each line is:
# time x y z H R P image
with open(imgorientf) as f:
    imgorient_str = f.readlines()

# load in the orientation of the images (note only for rectified images)
time_tsolid = []
camx_l = []
camy_l = []
camz_l = []
camH_l = []
camR_l = []
camP_l = []
imgname_tsolid = []

for i in imgorient_str:
     time_tsolid.append(float(i.split(' ')[0]))
     camx_l.append(float(i.split(' ')[1]))
     camy_l.append(float(i.split(' ')[2]))
     camz_l.append(float(i.split(' ')[3]))
     camH_l.append(float(i.split(' ')[4]))
     camR_l.append(float(i.split(' ')[5]))
     camP_l.append(float(i.split(' ')[6]))
     imgname_tsolid.append(i.split(' ')[7].replace('\n', ''))

# Load merged DSM file 
MakhEP = Project(projdir=projdir)
DSM_xr = rioxr.open_rasterio(MakhEP.DSM['010m'].filepath, parse_coordinates=True, masked=True)
DSM_xr = DSM_xr.where(DSM_xr > 0)

# # # PROCESSING

# Make a list of paths to all the rectified multi-band thermal images
# to loop over below
imgpaths = list(imgdir_rect.glob('*.tif'))

# For each rectified image
for imgpath in imgpaths:

    # start time
    starttime = time.time()

    # get imgname from path
    imgname = imgpath.parts[-1].split('.')[0]

    # get index of image in orientation file
    idx_tsolid = imgname_tsolid.index(imgname)

    # Use index to grab camx,y, and z for later
    camx = camx_l[idx_tsolid]
    camy = camy_l[idx_tsolid]
    camz = camz_l[idx_tsolid]

    # #  Open img
    img = rioxr.open_rasterio(filename=imgpath, parse_coordinates=True, masked=True)
    # extract the temp band and convert to kelvin
    img_kelvin = (img.sel(band=4)*0.04)
    # filter 0s as nans
    img_kelvin = img_kelvin.where(img_kelvin > 0)
    # Make signal raster
    img_signal = xr.apply_ufunc(phicalc, img_kelvin)

    # # Clip DSM to image extent and resample to image resolution
    # Get bounds of img
    bounds = img_signal.rio.bounds()
    # Clip using boundary
    # https://corteva.github.io/rioxarray/stable/examples/clip_box.html
    DSM_clip = DSM_xr.rio.clip_box(
        minx=bounds[0],
        miny=bounds[1],
        maxx=bounds[2],
        maxy=bounds[3],
    )
    # Resample DSM_clip to match img resolution
    DSM_resamp = DSM_clip.rio.reproject(
        img_signal.rio.crs,
        shape=(img_signal.rio.height, img_signal.rio.width),
        resampling=Resampling.bilinear
    )
    # Make ndarrays for x and y variables
    X_array, Y_array = np.meshgrid(DSM_resamp.x, DSM_resamp.y)
    
    # Make a distance xarray raster of the same shape and size
    dist_xr = np.sqrt(((camx - X_array)**2 + (camy- Y_array)**2 + (camz - DSM_resamp)**2))

    # # Calculating Tau (transmissivity of air)
    # MJ - Note that some of these atmos. Constants may be different for your camera:
    # Note: Temp in C
    def taucalc(distance, airtemp=curr_airtemp, RH=curr_relhumid,
                ATA1=0.006569, ATA2=0.01262, ATB1=0.002276,
                ATB2=0.00667, ATX=1.9):
        ch2o = (RH/100) * np.exp((1.5587 + 0.06939*(airtemp) - 0.00027816*(airtemp)**2 + 0.00000068455*(airtemp)**3))
        tau = ATX * np.exp(-np.sqrt(distance) * (ATA1 + ATB1 * np.sqrt(ch2o))) + (1 - ATX) * np.exp(-np.sqrt(distance) * (ATA2 + ATB2*np.sqrt(ch2o)))
        return tau

    # Make a tau (air transmissivity) raster
    tau_xr = xr.apply_ufunc(taucalc, dist_xr, curr_airtemp)

    # Plot real quick
    # fig, ax = plt.subplots()
    # tau_xr.plot(ax=ax)

    # # Make an air signal (air energy) raster
    def airsignalcalc(airtau, airtemp_kelvin=curr_airtemp):
        # Calc phi (energy)of air at given temp
        airphi = phicalc(airtemp_kelvin)
        # Tau of window (assuming 100% transmissivity thru the window)
        tauwin = 1
        # Calc airsignal
        airsignal = airphi * (1 - airtau) / tauwin
        return airsignal

    airsignal_xr = xr.apply_ufunc(airsignalcalc, tau_xr, curr_airtemp)

    # # Subtract airsignal from signal raster to get temp of targets
    img_signal_target = img_signal - airsignal_xr.sel(band=1).data

    # Convert img_signal_target to temp
    # Target temp (in celcius)
    target_temp_xr = xr.apply_ufunc(tempcalc, img_signal_target) - 273.15
    target_temp_xr = target_temp_xr.where(target_temp_xr > 0)

    # Export target temp to a geotif
    target_temp_xr.rio.to_raster(f"./data/out/processedthermalimages/{imgname}.tif")

    # end time
    endtime = time.time()

    # print(f'{endtime - starttime} seconds elapsed for {imgname} \n')
    # time.sleep(5)