#dependancies
import numpy as np
import os.path

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
    darkDicionary=dict(None)
    flatDicionary=dict(None)
    lightDicionary=dict(None)
    
    biasFits=[]

    for root, dirs, files in os.walk(inputPath):
        for file in files:
            if file.endswith(".fits"):
                hdul=loadFits(os.path.join(root, file))
                time=hdul.header["EXPTIME"]

                #find the type of image, then find its exposure time. 
                #see if any other of that type of image have been saved for that time. if not, make a list for that time.abs
                #append list for that exposure time with current file

                if(hdul.header["FRAMETYP"]=="Light"):
                    if lightDicionary.get(f"{hdul.header['FILTER']}-{time}")==None:
                        lightDicionary[f"{hdul.header['FILTER']}-{time}"]=[]
                    lightDicionary[f"{hdul.header['FILTER']}-{time}"].append(hdul.data)

                elif(hdul.header["FRAMETYP"]=="Dark"):
                    if darkicionary.get(time)==None:
                        lightDicionary[time]=[]
                    darkDicionary[time].append(hdul.data)

                elif(hdul.header["FRAMETYP"]=="Flat"):
                    if flatDicionary.get(f"{hdul.header['FILTER']}-{time}")==None:
                        flatDicionary[f"{hdul.header['FILTER']}-{time}"]=[]
                    flatDicionary[f"{hdul.header['FILTER']}-{time}"].append(hdul.data)         
                
                elif(hdul.header["FRAMETYP"]=="Bias"):
                    biasFits.append(hdul.data)
        



    darkTimes = ''
    for key in darkDicionary.keys():
        darkTimes=f'{darkTime},{key}' 
    
    flatTimes = ''
    for key in flatDicionary.keys():
        darkTimes=f'{flatTimes},{key}' 

    lightTimes = ''
    for key in lightDicionary.keys():
        lightTimes=f'{lightTimes},{key}'

    print(f'All dark time: \n {darkTimes}')
    print(f'All flat time: \n {flatTimes}')
    print(f'All light time: \n {lightTimes}')


