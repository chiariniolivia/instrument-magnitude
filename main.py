#dependancies
import numpy as np
import os.path
import csv
import datetime 

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
ignore = None
"""
ASSUMPTIONS:
    -all fits from within input path are same source
    -all darks, flats, and biases are provided within directory given as input
"""

"""
ARGUMENTS:
    [1] path to all fits files
    [2] filters to ignore {UVBRI}
"""

if __name__ == '__main__':

    print("rough instrument magnitude calculator now running")
    try:
        inputpath=argv[1]
    except IndexError:
        print("No file path provided as argument to python script.\n Exiting")
        exit(0)

    try:
        ignore=argv[2]
    except:
        pass
    

    #dictonaries of all .fits files
    darkDicionary=dict()
    flatDicionary=dict()
    lightDicionary=dict()
    
    biasFits=[]



    dTime=f'{datetime.datetime.now()}'
    year=f'{dTime[:4]}'
    month=f'{dTime[5:7]}'
    day=f'{dTime[8:10]}'
    hour=f'{dTime[11:13]}'
    minute=f'{dTime[14:16]}'
    second=f'{dTime[17:19]}'
    runTime =f'{year}{month}{day}T{hour}{minute}{second}'
    savePath=f'{inputpath}/inst-mag-script-{runTime}/'
    os.mkdir(savePath)
    for root, dirs, files in os.walk(inputpath):
        for file in files:
            if file.endswith(".fits"):
                hdul=loadFits(os.path.join(root, file))
                time=hdul.header["EXPTIME"]

                #find the type of image, then find its exposure time. 
                #see if any other of that type of image have been saved for that time. if not, make a list for that time.abs
                #append list for that exposure time with current file
                if(hdul.header["FRAMETYP"]=="Light"):
                    if hdul.header['FILTER'] in ignore:
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
                    if hdul.header['FILTER'] in ignore:
                        continue
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
    print(f'Times not represented in darks: {timesNotInDark}')

    #master darks for flat integration times
    masterDarkFlatdict=dict()
    for key in darkFlatDict: 
        masterDarkFlatdict[key]=None
        darkFlatStack = darkFlatDict[key]
        masterDarkFlat = np.median(darkFlatStack, axis=0)
        masterDarkFlatdict[key]= masterDarkFlat
        plotImg(masterDarkFlatdict[key], 2, f"Master Dark for Flats in {key}")
        plt.savefig(f'{savePath}/masterDarkFlat-{key}.png')
        plt.close()
    
    #master darks for light integration times
    masterDarkLightDict=dict()
    for key in darkLightDict:
        masterDarkLightDict[key]=None
        darkLightStack = darkLightDict[key]
        masterDarkLight = np.median(darkLightStack, axis=0)
        masterDarkLightDict[key]=masterDarkLight
        plotImg(masterDarkLightDict[key], 2, f"Master Dark for Light Image in {key}")
        plt.savefig(f'{savePath}/masterDarkLight-{key}.png')
        plt.close()

    #master flats for filter-time pairs
    masterFlatdict=dict()
    for key in flatDicionary:
        lightfilter, time = key.split('-')
        masterFlatdict[lightfilter]=None
        flatStack=flatDicionary[key]
        masterFlat= np.median(flatStack, axis=0) - masterDarkFlatdict[key]
        C = np.median(masterFlat)
        masterFlatdict[lightfilter]=masterFlat/C
        plotImg(masterFlatdict[lightfilter], 2, f"Master Flat in Filter {key}")
        plt.savefig(f'{savePath}/masterFlat-{key}.png')
        plt.close()


    #master lights
    scienceFramelist=[]
    for key in lightDicionary.keys():
        lightfilter, time = key.split('-')
        lightStack=lightDicionary[key]
        masterLight = (np.median(lightStack, axis=0) - masterDarkLightDict[key])/masterFlatdict[lightfilter]
        scienceFramelist.append((masterLight,float (time), lightfilter))

    #calculate statistics for science frames
    #consts.
    dw = 25
    radii = np.linspace(0, 2*dw, 2*dw)
    Y, X = np.ogrid[:dw*2, :dw*2]
    dist = np.sqrt((X-dw)**2 + (Y-dw)**2)
    ones = np.ones((dw*2, dw*2))
    radius = 15
    fields= ['source id','instrument magnitude', 'R^2']
    #loop through all science frames
    
    #show top of the telescope light frame
    for scienceFrame in scienceFramelist:
        mean, median, std, maximum = np.mean(scienceFrame[0]), np.median(scienceFrame[0]), np.std(scienceFrame[0]), np.max(scienceFrame[0])
        starFind = DAOStarFinder(threshold=median, fwhm=20.0, sky=mean, exclude_border=True, brightest=10, peakmax=maximum)
        sourceList = starFind(scienceFrame[0])

        os.mkdir(f'{savePath}/scienceFrame-{scienceFrame[2]}/')
        print(f'working on {savePath}/scienceFrame-{scienceFrame[2]}/')

        with open(f'{savePath}/scienceFrame-{scienceFrame[2]}/output.csv', 'a') as f:
            writer = csv.writer(f)
            writer.writerow(fields)


        plt.figure()
        plt.imshow(scienceFrame[0], cmap='gray', vmin=mean-0.5*std, vmax=mean+0.5*std)
        plt.title("Top-of-Telescope Light Frame [Counts]" + f"\nMean:{mean:0.2f}; Median: {median:0.2f}; STD: {std:0.2f}")
        plt.tight_layout()

        for source in sourceList:
            xc, yc = source[2], source[1]
            loc = (xc, yc)
            plt.gca().add_patch(Circle(loc[::-1] ,radius=radius,fill=False, edgecolor='m', alpha=0.5, zorder=100, lw=2.0, linestyle="-"))
            plt.text(yc+5, xc+5, f"{source[0]}")
        plt.savefig(f'{savePath}/scienceFrame-{scienceFrame[2]}/topOfTelescope.png')
        plt.close()
        
        for source in sourceList:
            sourceID = source[0]
            xc,yc = source[2], source[1]
            loc = (xc,yc)
            subFrame = scienceFrame[0][int(xc-dw):int(xc+dw), int(yc-dw):int(yc+dw)]
            radialData_raw=extract_radial_data(subFrame, xC=dw, yC=dw)[:dw]
            radialData = np.concatenate((radialData_raw[::-1], radialData_raw))
            p0 = [dw, 1, maximum, mean]
            params, R2 = fit_gaussian_1d(radii, radialData, p0)
            counts = np.sum(subFrame, where=dist<radius)
            nPix = np.sum(ones, where=dist<radius)
            countFlux = counts/nPix/scienceFrame[1]

            #Need to account for atmospheric extinction
            instMag = -2.5*np.log10(countFlux)   
            outputs = [sourceID, instMag, R2]
            with open(f'{savePath}/scienceFrame-{scienceFrame[2]}/output.csv', 'a') as f:
                writer = csv.writer(f) 
                writer.writerow(outputs)


            plt.figure()
            plt.subplot(1,2,1)
            plt.imshow(subFrame - params[3], cmap='gray')
            plt.gca().add_patch(Circle((dw, dw),radius=radius,fill=False, edgecolor='m', alpha=0.5, zorder=100, lw=2.0, linestyle="--"))
            xLabels = np.concatenate((np.linspace(0,dw,5)[::-1], np.linspace(0,dw,5)[1:]))

            plt.xticks(np.linspace(0,dw*2,len(xLabels)), xLabels, rotation=45)
            plt.yticks(np.linspace(0,dw*2,len(xLabels)), xLabels)



            plt.subplot(1,2,2)

            plt.plot(radii, radialData-params[3], 'b.')
            plt.plot(radii, gaussian_1d(radii, *params[:-1], 0), 'r')
            plt.grid(1)

            plt.axvline(x = params[0]-radius, color = 'm', linestyle="--")
            plt.axvline(x = params[0]+radius, color = 'm', linestyle="--")
            plt.xticks(np.linspace(0, dw*2, len(xLabels)), xLabels, rotation=45)
            plt.suptitle(f"Filter V PhotUtils Source ID {sourceID} w/o Background\nFit $R^2=${R2:0.4f}; InstMag = {instMag:0.3f}")
            plt.savefig(f'{savePath}/scienceFrame-{scienceFrame[2]}/source-{sourceID}.png')
            plt.close()