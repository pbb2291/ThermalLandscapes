from dotenv import load_dotenv
import os
import numpy as np
import pandas as pd

# Function for calculating temp from signal (note in kelvin)
def tempcalc(signal, R=366545, B=1428, F=1, O=342):
    temp = (B / np.log((R / (signal - O)) + F))
    return temp

# Calculating signal (phi or energy)
# note temp is in kelvin
def phicalc(temp_kelvin, R=366545, B=1428, F=1, O=342):
    phi = (R / (np.exp(B / temp_kelvin) - F)) + O
    return phi

# # Calculating Tau (transmissivity of air)
# MJ - "Note that some of these atmos. Constants may be different for your camera"
# Note: Temp in C
def taucalc(distance, airtemp, RH,
            ATA1=0.006569,
            ATA2=0.01262,
            ATB1=0.002276,
            ATB2=0.00667,
            ATX=1.9):
    ch2o = ((RH/100) * np.exp((1.5587 + 0.06939*(airtemp)) - (0.00027816*(airtemp)**2) + (0.00000068455*(airtemp)**3)))
    tau = ATX * np.exp(-np.sqrt(distance) * (ATA1 + ATB1 * np.sqrt(ch2o))) + (1 - ATX) * np.exp(-np.sqrt(distance) * (ATA2 + ATB2*np.sqrt(ch2o)))
    # if np.isfinite(tau):
    #     tau = 1
    # if not np.isnan(tau):
    #     tau = 1
    return tau

# # Make an air signal (air energy) raster
def airsignalcalc(airtau, airtemp_kelvin):
    # Calc phi (energy)of air at given temp
    airphi = phicalc(airtemp_kelvin)
    # Tau of window (assuming 100% transmissivity thru the window)
    tauwin = 1
    # Calc airsignal
    airsignal = airphi * (1 - airtau) / tauwin
    return airsignal


def calcTargetTemp(phitotal, 
                   etarget,
                   phiair,
                   phisky,
                   tauair=1,
                   esky=1,
                   phiwin=1,
                   reflwin=1,
                   tauwin=1):
    
    # p = (phitotal / (etarget*tauair*tauwin)) -
    #     ((phisky * esky * (1 - etarget)) / etarget) -
    #     ((phiair * (1 - tauair))/(etarget * tauair)) -
    #     ((phiwin * (1 - reflwin - tauwin)) / (etarget * tauair * tauwin))
    
    # Sky temperature is the temperature you'd record with an infrared thermometer if you pointed it at the sky -- 
    # typically considerably colder than ambient air temperature. Clear skies are colder than cloudy skies; 12 degrees C is ~clouds.
    
    # Phirefl is reflected energy off the target. In reality, it's coming from everywhere, 
    # but we simplify and say it's from the sky. 
    
    # Phiair is energy the sensor sees in the parcel of air it's looking through to the target.
    
    # Updated 1/28/24 -
    # This version adds in the sky term,
    # which was previously left out.
    p = (phitotal / (etarget*tauair*tauwin)) - ((phisky * esky * (1 - etarget)) / etarget) - ((phiair * (1 - tauair))/(etarget * tauair)) 
    
    # convert to temp
    t = tempcalc(p)
        
    return t

# Assign emissivities to veg types,
# and see if it makes a difference
# Mainly using emissivity values from Griend and Owe, 1993, Table 2
# Also, this paper for more context
# https://www.researchgate.net/publication/341253627_Applications_of_UAV_Thermal_Imagery_in_Precision_Agriculture_State_of_the_Art_and_Future_Research_Outlook/figures?lo=1
def assign_emissivity(vegtype):
    
    e = 1.0 
    
    if ((vegtype == 'tree') |( vegtype == 'shrub')):
        e = 0.97 # from Griend and Owe 1993 - Table 2 - average of all shrub values
        
    elif vegtype == 'grass':
        e= 0.96 # from Griend and Owe 1993 - Table 2 - Long Grass
        
    elif vegtype == 'ground':
        e = 0.91 # from Griend and Owe 1993 - Table 2 - Bare Soil
        
    return e


def merge_imgtraj_pixel_dfs(imgtraj_df, pixel_df):

    # Pair down and clean image traj df
    # so that you can merge it with pixel_df
    imgtraj_df = imgtraj_df[['Imagef', 'Altitude', 'AirTemp', 'RH', 'STRD', 'DateTime']]
    
    # make an image file column for joining dfs 
    imgtraj_df['imgf'] = imgtraj_df.Imagef.apply(lambda i: i + '.tif')
    
    # drop column without a use
    imgtraj_df = imgtraj_df.drop('Imagef', axis=1)
    
    # now, make a new pixel_df that's got weather data in it
    pixel_df = pixel_df.merge(imgtraj_df)
    
    return pixel_df


def correct_temps(imgtraj_df, pixel_df):

    # load relevant file paths from .env file
    load_dotenv()

    # set all your file paths
    DATADIR = os.getenv("DATADIR")
    FIGDIR = os.getenv("FIGDIR")
    STATSDIR = os.getenv("STATSDIR")
    OUTDIR = os.getenv("OUTDIR")
    TRAJDIR = os.getenv("TRAJDIR")

    # merge the two inputs processing
    df = merge_imgtraj_pixel_dfs(imgtraj_df=imgtraj_df, pixel_df=pixel_df)

    # This is the air temp that FLIR used to calculate temperature in the images
    # likewise, the other parameters below are the FLIR camera parameters
    # (or our best guess at them) for converting signal -> temp
    default_airtemp = 22

    # # # Convert temp -> signal -> temp

    # Add 'distance' column by subtracting elevation and altitude
    # Note: Had to add an else here to fix distance calculation... cause it comes up negative for RoanCamp
    pixel_df['Distance'] = [a - e if ((a-e)>0) else 100.0 for a, e in zip(pixel_df['Altitude'], pixel_df['Elevation'])]

    # Create a signal (phi) column
    df['phi'] = df.Temperature.apply(lambda t: phicalc(temp_kelvin = t + 273.15))

    # Create a tau (air transmissivity) column
    df['tau'] = df.apply(lambda x: taucalc(distance = x['Distance'],
                                                       airtemp = x['AirTemp'],
                                                       RH = x['RH']/100),
                                     axis=1)

    
    # Assign emissivities to veg types (using function defined above)
    pixel_df['eVeg'] = pixel_df.VegType.apply(lambda v: assign_emissivity(v))
    
    # # # Compute Emissivity by vegtype
    pixel_df[f'Temperature_eVeg'] = pixel_df.apply(lambda x: calcTargetTemp(phitotal=x['phi'],
                                                                            etarget=x['eVeg'],
                                                                            phiair=airsignalcalc(airtemp_kelvin=default_airtemp + 273.15,
                                                                                                 airtau=x['tau']),
                                                                            phisky=(x['STRD']/3600), # new addition of thermal radiation downwarss, converted to Watts (/3600 seconds)
                                                                            tauair=x['tau']) - 273.15,
                                                   axis=1)
    

    return pixel_df


if __name__ == "__main__":

    import matplotlib.pyplot as plt

    # Buffalo camp for testing
    s = ['BuffaloCamp']
    s_it = ['BuffaloCamp100m']
    
    # # # Make a new pixel dataframe with weather data, everything you'd need for corrections
    # read in image trajectories
    imgtraj_df = pd.read_csv(f'{DATADIR}/in/ImageTrajectories/{s_it}_ImageTrajectory_withWeather.csv')
    
    # read in image pickle file 
    pixel_df = pd.read_pickle(f'{DATADIR}/out/pixelpickles/{s}_pixels.obj')

    # # # Run correction
    pixel_df_correct = correct_temps(imgtraj_df=imgtraj_df, pixel_df=pixel_df)
    
    # # # Plot and save figure
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
    
    sns.kdeplot(data=pixel_df, x='Temperature', ax=ax1, hue='Exclosure',
                palette='plasma_r',
                fill=True,
                hue_order=exclosure_hues,
                legend=True,
                common_norm=False)

    sns.kdeplot(data=pixel_df, x='Temperature_eVeg', ax=ax2, hue='Exclosure',
                palette='plasma_r',
                fill=True,
                hue_order=exclosure_hues,
                legend=True,
                common_norm=False)

    ax3.scatter(pixel_df_correct.Temperature,
                pixel_df_correct.Temperature_eVeg,
                c = pixel_df_correct.Temperature - pixel_df_correct.Temperature_eVeg)

    ax3.plot([0, pixel_df_correct.Temperature.max()],
             [0, pixel_df_correct.Temperature_eVeg.max()],
             'k')

    ax3.set_xlabel('Temp Before')
    ax3.set_xlabel('Temp After Correction')

    fig.tight_layout()

    plt.show()

    # fig.savefig(f'{FIGDIR}/KDEPlot_SensitivityTests_{s}.png', dpi=300)