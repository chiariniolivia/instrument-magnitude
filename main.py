#dependancies
import numpy as np
import os.path
import csv


from astropy.io import fits
from matplotlib import pyplot as plt
from matplotlib.patches import Circle

from glob import glob
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit
from astropy.nddata import CCDData
from sys import argv 
from photutils.detection import DAOStarFinder,IRAFStarFinder


#functions
def loadFits(inputpath):
    return fits.open(f'{inputpath}')[0]


def plotImg(image, numStd, titleStr):
  mean, std, median = np.mean(image), np.std(image), np.median(image)

  plt.figure()
  plt.imshow(image, cmap='gray', vmin=mean-numStd*std, vmax=mean+numStd*std)
  plt.title(titleStr + f"\nMean:{mean:0.2f}; Median: {median:0.2f}; STD: {std:0.2f}")
  plt.tight_layout()
  plt.colorbar()

def extract_radial_data(data, xC, yC):
    #Get matrix of integer indices associated with subFrame
    y, x = np.indices((data.shape))

    #Generate matrix of radius values
    r = np.sqrt((x - xC)**2 + (y - yC)**2)

    #Force integer values (np.sqrt gives floats)
    r = r.astype(int)

    #Generate a histogram of radius bin, each weighed by corresponding counts
    weightedRadiusHistogram = np.bincount(r.ravel(), weights=data.ravel())
    unweightedRadiusHistogram = np.bincount(r.ravel())

    #Get average for each radius bin
    averageCountsPerRadiusBin = weightedRadiusHistogram / unweightedRadiusHistogram
    return averageCountsPerRadiusBin

def gaussian_1d(x, mu, sigma, amplitude, offset):
    #Model function as gaussian with amplitude A and offset G
    return amplitude * np.exp( -((x-mu)/sigma)**2/2 ) + offset


def fit_gaussian_1d(x_data, y_data, p0):
    # p0 behaves by taking a best guess at params (mu, sigma, amplitude, offset)
    params, _ = curve_fit(gaussian_1d, x_data, y_data, p0)

    #Calculate coefficient of determination
    res = y_data - gaussian_1d(x_data, *params)
    sumSqrs_res = np.sum(res*res)
    totSumSqrs = np.sum((y_data-np.mean(y_data))**2)
    R2 = 1.0 - (sumSqrs_res / totSumSqrs)

    return params, R2



#main script

inputpath = None

"""
ASSUMPTIONS:
    -all fits from within input path are same source
    -all darks, flats, and biases are provided within directory given as input
"""
if __name__ == '__main__':

    print("rough instrument magnitude calculator now running")
    try:
        inputPath=argv[1]
    except IndexError:
        print("No file path provided as argument to python script.\n Exiting")
        exit(0)



    #dictonaries of all .fits files
    darkDicionary=dict()
    flatDicionary=dict()
    lightDicionary=dict()
    
    biasFits=[]



    """

    add thing to exlucde ultraviolet or ingrared frames if flats not taken for lgiht expsoure times

    """

    for root, dirs, files in os.walk(inputPath):
        for file in files:
            if file.endswith(".fits"):
                hdul=loadFits(os.path.join(root, file))
                time=hdul.header["EXPTIME"]

                #find the type of image, then find its exposure time. 
                #see if any other of that type of image have been saved for that time. if not, make a list for that time.abs
                #append list for that exposure time with current file
                if(hdul.header["FRAMETYP"]=="Light"):
                    if hdul.header['FILTER'] =='U':
                        continue
                    if   f"{hdul.header['FILTER']}-{float (time):0.2f}" not  in lightDicionary:
                        lightDicionary[f"{hdul.header['FILTER']}-{float (time):0.2f}"]=[]
                    lightDicionary[f"{hdul.header['FILTER']}-{float (time):0.2f}"].append(hdul.data)

                elif(hdul.header["FRAMETYP"]=="Dark"):
                    if  f'{float (time):0.2f}' not in darkDicionary:
                        darkDicionary[f'{float (time):0.2f}']=[]
                    darkDicionary[f'{float (time):0.2f}'].append(hdul.data)

                elif(hdul.header["FRAMETYP"]=="Flat"):
                    #if hdul.header['FILTER'] =='U':
                      # continue
                    if f"{hdul.header['FILTER']}-{float (time):0.2f}" not in flatDicionary:
                        flatDicionary[f"{hdul.header['FILTER']}-{float (time):0.2f}"]=[]
                    flatDicionary[f"{hdul.header['FILTER']}-{float (time):0.2f}"].append(hdul.data)         
                
                elif(hdul.header["FRAMETYP"]=="Bias"):
                    biasFits.append(hdul.data)

    timesNotInDark = []
    darkFlatDict = dict()
    for key in flatDicionary:
        throwaway, time = key.split('-')
        time = f'{float (time):0.2f}'
        if time not in darkDicionary:
            timesNotInDark.append(time)
        else:
            darkFlatDict[key]=darkDicionary[time]
    darkLightDict = dict()
    for key in lightDicionary:
        throwaway, time = key.split('-')
        time = f'{float (time):0.2f}'
        if time not in darkDicionary:
            timesNotInDark.append(time)
        else:
            darkLightDict[key]=darkDicionary[time]
    print(f'{timesNotInDark}')

    """
    something here to check if all needed redux files are there, i.e. if 5s exposure, 5s dark,
    need to make sure all flats have a dark AND all lights have a dark, if not then make a guess
    ensure darkFlatDict and flatDictionary have the same keys, needs to add filter
    ensure darkLightDict and lightDictionary have the same keys, needs to add filter
    """

    exit(0)

    #master darks for flat integration times
    masterDarkFlatdict=dict()
    for key in darkFlatDict.keys(): 
        masterDarkFlatdict[key]=None
        darkFlatStack = darkFlatDict[key]
        masterDarkFlat = np.median(darkFlatStack, axis=0)
        masterDarkFlatdict[key]= masterDarkFlat
        plotImg(masterDarkFlatdict[key], 2, "Master Dark for Flats in Filter --")

    
    #master darks for light integration times
    masterDarkLightdict=dict()
    for key in darkLightDict.keys():
        masterDarkLightdict[key]=None
        darkLightStack = darkLightDict[key]
        masterDarkLight = np.median(darkLightStack, axis=0)
        masterDarkLightdict[key]=masterDarkLight
        plotImg(masterDarkLightdict[key], 2, "Master Dark for Light Image in time ---")

    #master flats for filter-time pairs
    masterFlatdict=dict()
    for key in flatDicionary.keys():
        masterFlatdict[key]=None
        flatStack=flatDicionary[key]
        masterFlat= np.median(darkFlatStack, axis=0) - masterDarkFlatdict[key]
        C = np.median(masterFlat)
        masterFlatdict[key]=masterFlat/C
        plotImg(masterDarkFlatdict[key], 2, "Master Flat in Filter --")

    #master lights
    scienceFramelist=[]
    for key in lightDicionary.keys():
        lightStack=lightDicionary[key]
        masterLight = (np.median(lightStack, axis=0) - masterDarkLightDict[key])/masterFlatdict[key]
        scienceFramelist.append(masterLight)
        plotImg(scienceFramelist[-1], 2, "Top of the Telescope Light Frame for ---")
    
    #calculate statistics for science frames
    
    #consts.
    dw = 25
    radii = np.linspace(0, 2*dw, 2*dw)
    Y, X = np.ogrid[:dw*2, :dw*2]
    dist = np.sqrt((X-dw)**2 + (Y-dw)**2)
    ones = np.ones((dw*2, dw*2))
    radius = 15

    #loop through all science frames
    
    """ 
    generate a CSV file, each science frame will have its own file. each source within the frame its own row. coresponding images with numbered stars as well.


    """


    for scienceFrame in scienceFramelist:
        mean, median, std, max = np.mean(scienceFrame), np.median(scienceFrame), np.std(scienceFrame), np.max(scienceFrame)
        sourceList = DAOStarFinder( scienceFrame,threshold=median, fwhm=20.0, sky=mean, exclude_border=True, brightest=10, peakmax=max)

        

    
    