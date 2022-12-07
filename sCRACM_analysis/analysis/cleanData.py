"""
To extract raw data from XSG files, process, and save in .csv file format. 

Functions: 
Helpers:
    - _check_keys: checks if entries in dic are mat-objects, if yes, todict is called to change to nested dicts 
    - _todict: recursive, contructs from nested dicts 
    - loadmat: use function instead of scipy.io.loadmat, cures problem of not properly recovering python dicts from mat and .m files. 
    - cart2pol: cartesian to polar coordinate transformer 
    - pol2cart: polar to cartesian coordinate transformer 
    - somaPositionTransformer: takes raw soma coordinates from sCRACM map (designated during experiment) and applies rotation and pattern offset 
    - analysisPlots: plots sCRACM maps (mean, min, onset, charge) for a single cell 
    - mousePoints: called by measure_points, stores left mouse click coordinates when clicking on .tif image. Measures real world distance between points 
    - measure_points: driver function for mousePoints, used to measure distance between microscope .tif image (4x LSPS rig calibration only)
    - plotMappingArray: plots mapping array over acute slice image. Used to determine layer1Row. 
    - loadXSG: helper function called by averageMaps, takes amplifer number. Will open a matlab XSG file from a single mapping sweep, and saves salient data in python dictionary.
    - averageMaps: helper function used to average multiple sweeps from a single experiment into one map.
Execution functions: 
    -     - analyzeSCRACMmap: used to average sweeps from a single experiment into one map, optioal save flag to save .csv file for cell, used by all later functions
"""

#imports 
import numpy as np
import pandas as pd
import scipy as sc
from scipy import io
import math
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
import sys
from tkinter import filedialog
from tkinter import Tk
import cv2
import csv
from sCRACM_analysis.sCRACM_global import path_database,filename_cell_db,path_ephys_data


def _check_keys( dict):
    """
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    """
    for key in dict:
        if isinstance(dict[key], io.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict

def _todict(matobj):
    """
    A recursive function which constructs from matobjects nested dictionaries
    """
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, io.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

def loadmat(filename):
    """
    this function should be called instead of direct scipy.io .loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    """
    data = io.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)


#cartesian to polar coordinates transformer 
def cart2pol(x, y):
    '''
    Parameters:
    - x: float, x coord. of vector end
    - y: float, y coord. of vector end
    Returns:
    - r: float, vector amplitude
    - theta: float, vector angle
    '''

    z = x + y * 1j
    r,theta = np.abs(z), np.angle(z)

    return r,theta

def pol2cart(r,theta):
    '''
    Parameters:
    - r: float, vector amplitude
    - theta: float, vector angle
    Returns:
    - x: float, x coord. of vector end
    - y: float, y coord. of vector end
    '''

    z = r * np.exp(1j * theta)
    x, y = z.real, z.imag
    return x,y


def somaPositionTransformer(somaX,somaY,spatialRotation,xPatternOffset,yPatternOffset):
    """"
    Helper function, will take original soma coordinates (x,y) aquired during experiment
    apply rotation and pattern offset to allow for accurate plotting on top of map.
    """
    somaXoffset = somaX - xPatternOffset
    somaYoffset = somaY - yPatternOffset
    rotationAngleRadians = (-1)*spatialRotation*(math.pi/180)

    rho,theta = cart2pol(somaXoffset,somaYoffset)
    somaXnew,somaYnew = pol2cart(rho,theta+rotationAngleRadians)
    return somaXnew,somaYnew

def analysisPlots(map_dict,ampNum):
    """"
    Helper function, plots analysis plots (mean, min,onset,charge)
    after averaging on single cell 
    """
    #to plot 
    mean= map_dict['mean']
    minimum = map_dict['minimum']
    onset = map_dict['onset']
    integral = map_dict['integral']

    to_plot = [mean,minimum,onset,integral]

    mapSpacing = 50
    rcAddOn =0
    somax = map_dict['soma1Coordinates'][0]
    somay = -map_dict['soma1Coordinates'][1]
    #print('somaY cord that plotting function is getting: '+str(somay))
    fig,ax = plt.subplots(2,2,sharex = 'col',sharey='row',figsize=(16,12))
    ax = ax.flatten()

    [r,c] = np.shape(mean)
    dx = mapSpacing #50
    dy = mapSpacing #50 
    xdata = np.arange(0,(c*dx),dx) - int((c-1) * dx/2)
    ydata = np.arange(0,(r*dx),dx) - int((r-1) * dx/2)
    xmax = (c+rcAddOn)*dx/2   #xmax = (map size in x direction * map spacing) - (mapsize in x)*50/2 * note: I removed the c-1 and r-1 to fix the problem with the soma plotting
    xmin = -xmax # x data is centered on 0, so min and max are the same here, with signs reversed
    ymax = (r+rcAddOn)*dx/2 #ydata also centered on zero
    ymin = -ymax

    for n in range(np.shape(to_plot)[0]):
        img = ax[n].imshow(to_plot[n],aspect =
                      'auto', cmap = 'magma_r',extent = [xmin,xmax,ymax,ymin]);
        img.axes.tick_params(axis = 'both',which = 'both',bottom=False,top=False,labelbottom=False,labelleft=False)
        #img.axes.tick_params(axis = 'y',which = 'both',left=False,right=False,labelright=False)
        ax[n].plot(somax,somay,'^',markersize=10, color='b');
        cbaxes = inset_axes(ax[n],width = "3%", height = "75%",loc = 'lower right')
        color_bar = plt.colorbar(img,cax=cbaxes);
        color_bar.ax.yaxis.set_tick_params(color = 'white')
        plt.setp(plt.getp(color_bar.ax.axes,'yticklabels'),color='white')
    ax[0].set_title('Mean EPSC',color = 'white')
    ax[1].set_title('Minimum',color = 'white')
    ax[2].set_title('Onset',color = 'white')
    ax[3].set_title('Charge',color = 'white')


def mousePoints(event,x,y,flags,img):
    """helper function, stores mouse click coordinates 
    and converts to real world distances between points. Calibrated to 4x objective on
    LSPS rig. """
    global previous_point
    #left button click 
    if event == cv2.EVENT_LBUTTONDOWN:
        if previous_point!= None:
            x2,y2 = previous_point
            dist = ((x-x2)**2 + (y-y2)**2 )**0.5
            real_dist =(dist/1600)*2657 #note: the horizontal width is for the 4x objective 
            #on the rig. image width defaul is 1600 (rig camara default is 1200by1600)
            print('distance:')
            print(f'{real_dist} um')
            cv2.line(img,(x,y),(x2,y2),(0,0,0),2)
            #redraw point to hide begining of line 
            cv2.circle(img,(x2,y2),3,(0,0,255),-1)
        
        cv2.circle(img,(x,y),3,(0,0,255),-1)
        cv2.imshow('image',img)
        previous_point = (x,y)

def measure_points():
    """Driver function used to measure micron distance between 
    two mouse clicks on a .tif image"""
    global previous_point
    
    previous_point = None #initialize the point_counter, should always start at None
    print('Please select cell location .tiff image captured during experiment')
    Tk().withdraw() #hide Tk root folder 
    img_path = filedialog.askopenfilename();
    img=cv2.imread(img_path)
    

    
    cv2.imshow('image', img)
    cv2.setMouseCallback('image',mousePoints,img)
    cv2.waitKey(0)
    cv2.destroyAllWindows()
    cv2.waitKey(1);
    return img_path


def plotMappingArray(img_path,map_analysis_dict):
    """
    Will plot mapping array over acute slice image. Use to determine layer1Row number! 
    """
    
    im = plt.imread(img_path)
    R = im/(np.max(im[:])/110)
   # variables coded from the workspace variable: 
    XRange =map_analysis_dict['horizontalVideoDistance']
    YRange = map_analysis_dict['verticalVideoDistance']
    xSpacing = map_analysis_dict['xSpacing']
    ySpacing = map_analysis_dict['ySpacing']
    theta = map_analysis_dict['spatialRotation']
    pattern = map_analysis_dict['pattern'][0]
    xPatternOffset = map_analysis_dict['xPatternOffset'] #check the affind tranform lines 
    yPatternOffset = map_analysis_dict['yPatternOffset']
    rotation = (theta/180*np.pi) #for rotating the tracing (user rotation plus experiment rotation)
    
    #initilize figure 
    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(111)
    extent = (-.5*XRange,.5*XRange, -.5*YRange,.5*YRange)
    #plot acture slice image 
    im = plt.imshow(R, extent = extent,cmap = 'gray'); #extent key argument most similar to Xdata and Ydata
    ax.set_title(map_analysis_dict['experimentNumber'], color = 'w');
    ax.set_aspect(1)
    ax.set_xlim([-XRange/2,XRange/2])
    ax.set_ylim([-YRange/2,YRange/2])

    #create coordinats for map array from map pattern size! 
    arrayX= np.arange(pattern.shape[1])*xSpacing - 0.5 * (pattern.shape[1]-1)*xSpacing #create 1 row of x cords 
    xorig= np.reshape(np.tile(arrayX,(pattern.shape[0],1)),(pattern.size),order='F') #repeat coords for all rows, and reshape into 1D array, need fortran like indexing for ordered coords
    arrayY= np.arange(pattern.shape[0])*ySpacing - 0.5 * (pattern.shape[0]-1)*ySpacing
    yorig= np.reshape(np.tile(arrayY,(pattern.shape[1],1)),(pattern.size)) #don't need fortran like indexing here! 

    #preform rotation
    [rho,theta2] = cart2pol(xorig,yorig)
    [xpoints,ypoints] = pol2cart(rho,theta2)
    xpoints = xpoints+xPatternOffset
    ypoints = ypoints+yPatternOffset
    #plot mapping array 
    ax.plot(xpoints,ypoints,'o',markersize = 4);
    plt.axis('off')

def loadXSG(ampNumber):
    """ This helper function is called by averageMaps() and takes the amplifier number from an experiment, 
     and opens a pop up dialogue for you to select a specific .XSG map file (Matlab data file). It will parse the 
     selected file, and return a dictionary of parameters used for sCRACM and LSPS map averaging.
    """
    #1: open dialogue box
    print('Please select XSG file for single map sweep')
    Tk().withdraw() #hide Tk root folder
    xsg_path = filedialog.askopenfilename()
    #load using loadmat function 
    xsg_file = loadmat(xsg_path) 
    
    #make new dict, keep only important info
    #assing to dict using class object, allows for easy assignment
    class NewClass(object): pass
    map_dict = NewClass()

    #info needed: 

    #sample rate 
    map_dict.sampleRate = xsg_file['header']['ephys']['ephys']['sampleRate'] #sample rate 
    sr = xsg_file['header']['ephys']['ephys']['sampleRate'] 
    # trace data 
    ampNumber = 1 #note, hard coded for now. Fix when turning into a function! 
    if ampNumber == 1:
        traceData = xsg_file['data']['ephys']['trace_1']
    elif ampNumber == 2: 
        traceData = xsg_file['data']['ephys']['trace_2']
    traceLength = xsg_file['header']['mapper']['mapper']['isi'] * sr
    #reshape trace, make each laser flash data a individual row!
    # (768000,1) reshaped to (4000,192)
    #note: you must use Fortran-like index ordering to get correct ordering
    #becuase we want "unstack" the data columnwise. 
    map_dict.traces= np.reshape(traceData,(int(traceLength),int(len(traceData)/traceLength)),order = 'F')
    
    #possibly add step to create a new dictionary with just the relevant information
    # possibly add step to reshape the traceData from [4000,] to [4000,1] array
    traceLength = xsg_file['header']['mapper']['mapper']['isi'] * sr
    # pattern
    map_dict.pattern = xsg_file['header']['mapper']['mapper']['mapPatternArray']
    # spacing (x and y)
    map_dict.xSpacing = xsg_file['header']['mapper']['mapper']['xSpacing']
    map_dict.ySpacing = xsg_file['header']['mapper']['mapper']['ySpacing']
    # laser power
    map_dict.laserPower = xsg_file['header']['mapper']['mapper']['specimenPlanePower']
    #sampling rate: defined earlier: include in the data structure as well
    # soma coordinates 
    map_dict.soma1Coordinates = xsg_file['header']['mapper']['mapper']['soma1Coordinates']
    #print(f"from raw xsg: "+str(xsg_file['header']['mapper']['mapper']['soma1Coordinates']))
    map_dict.soma2Coordinates = xsg_file['header']['mapper']['mapper']['soma2Coordinates']
    # experiment number 
    map_dict.experimentNumber = xsg_file['header']['xsg']['xsg']['initials'] + xsg_file['header']['xsg']['xsg']['experimentNumber']
    # pattern rotation
    map_dict.patternRotation = xsg_file['header']['mapper']['mapper']['patternRotation']
    # pattern flip bool
    map_dict.patternFlip = xsg_file['header']['mapper']['mapper']['patternFlip']
    # spatial rotation
    map_dict.spatialRotation = xsg_file['header']['mapper']['mapper']['spatialRotation']
    # pattern offset 
    map_dict.xPatternOffset = xsg_file['header']['mapper']['mapper']['xPatternOffset']
    map_dict.yPatternOffset = xsg_file['header']['mapper']['mapper']['yPatternOffset']
    # voltage step
    map_dict.voltageStep = xsg_file['header']['ephys']['ephys']['pulseParameters']['amplitude']
    # voltage step start time and width 
    map_dict.VStepStartT= xsg_file['header']['ephys']['ephys']['pulseParameters']['squarePulseTrainDelay']
    map_dict.VstepWidth = xsg_file['header']['ephys']['ephys']['pulseParameters']['squarePulseTrainWidth']
    # horizontal and vertical video distance 
    try: 
        map_dict.horizontalVideoDistance = xsg_file['header']['imagingSys']['imagingSys']['xMicrons']
        map_dict.verticalVideoDistance = xsg_file['header']['imagingSys']['imagingSys']['yMicrons']
    except:
        map_dict.horizontalVideoDistance = 2657 #for 4x objective when the header is missing 
        map_dict.verticalVideoDistance = 1993 
    return map_dict.__dict__    

def averageMaps(numMap):
    """
    takes a user defined number of map sweeps from a single experiment and averages map responses together. 
    Assumptions: 1. stimulation for a single point occures after 100ms 
                 2.  max response of single stimulation even takes place within 50ms of stimulation, or 1500 from start
                 3. stimulation window is 1000 ms 
    Other important details: 
                1. A significant response must be > 6 std above baseline 
                2. A significant reponse for a single point stimulation must be found in more than half of the maps, or it will be considered a spontaneous even and will be wiped out. 
    """
    ampNumber = 1 #default for now. If duel patch experiments, pull this to 
    loadedData = {}
    plotFlag = 1 #when != 0, plot individual maps 

    #load in the xsg files cooresponding to the maps you want to investigate 
    for i in range(numMap):
        loadedData[i] = loadXSG(ampNumber)

    # get small voltage injection to measure Rs for quality control 
    RsVoltageStep = loadedData[0]['voltageStep']
    VStepStartT = loadedData[0]['VStepStartT']
    VStepWidth = loadedData[0]['VstepWidth'] #fix type, capitalieze Step here and in previous function
    samplingRate = loadedData[0]['sampleRate']
    VStepStartInd = VStepStartT*samplingRate
    VStepEndInd = VStepStartInd + VStepWidth*samplingRate

    noFlipTracker = 0


    #Preallocate variables: 
    plotFlag = 1
    _,c=loadedData[0]['traces'].shape
    laserPower = []
    patternFlip = np.empty((1,numMap))
    pattern = []
    dataArray = []
    Racc = np.empty((numMap,c))
    Rtotal = np.empty((numMap,c))
    Rinput = np.empty((numMap,c))
    Rinput_ratio = np.empty((numMap,c))

    for XSG_num in range(numMap):
        laserPower.append(loadedData[XSG_num]['laserPower'])
        patternFlip[0,XSG_num] = loadedData[XSG_num]['patternFlip']
        pattern.append(loadedData[XSG_num]['pattern'])
        _,c = loadedData[XSG_num]['traces'].shape
        dataArray.append(loadedData[XSG_num]['traces'])
        # use 5 ms before the end of the small votlage step to calculate steady state current and Rtotal
        I_Rtotal = np.mean(dataArray[XSG_num][int(VStepEndInd)-50:int(VStepEndInd)-1],axis = 0)
        #use immediate current after voltage step to calculate access resistance
        I_Racc = np.min(dataArray[XSG_num][int(VStepStartInd):int(VStepStartInd)+50],axis = 0)
        I_baseline = np.mean(dataArray[XSG_num][int(VStepStartInd)-50:int(VStepStartInd)-1],axis = 0)
        RsVoltageStep_arr = np.full((1,c),RsVoltageStep)
        #multiple by 1,000 to convert to mega Ohms
        Racc[XSG_num,:] = RsVoltageStep_arr/(I_Racc-I_baseline)*1000
        Rtotal[XSG_num,:] = RsVoltageStep_arr/(I_Rtotal-I_baseline)*1000
        Rinput[XSG_num,:] = Rtotal[XSG_num,:]-Racc[XSG_num,:]
        Rinput_ratio[XSG_num,:] = Rinput[XSG_num,:]/Rtotal[XSG_num,:]
        rows,NumTraces = dataArray[XSG_num].shape
        if noFlipTracker==0:#if non-flipped pattern has not been found, keep searching
            if patternFlip[0,XSG_num]==0: #if current pattern is not flipped 
                standardPattern = pattern[XSG_num] #use it as the standard pattern
                noFlipTracker = 1 # update that non-flipped pattern has been found (no search in later loops)
        if plotFlag: #if plotflag is something other than zero
            test = 0 #put in code to plot traces as one map, one trial! 
    if noFlipTracker==0: #if still no non-flipped pattern is found (all maps where flipped)
        standardPattern = np.flip(pattern[0],axis = 1) #flip first pattern to get standard pattern


    baselineStartIndex = 900 
    baselineEndIndex = 999
    stimOnInd = 1000
    peakSearchEnd = 1500 #min val should come within 50ms of stimulation! 

    # subtract mean baseline from dataArray 
    #preallocate some variables: 

    baselineMedians = np.empty((numMap,c))
    SD = []
    bsArray = []

    #note: consider when baseline is influced by previous large response: may need to adjust this 
    # next part 
    for XSG_num in range(numMap):
        baselineMedians[XSG_num,:]= np.median(dataArray[XSG_num][baselineStartIndex:baselineEndIndex],axis = 0)
        baselineSDs = np.std(dataArray[XSG_num][baselineStartIndex:baselineEndIndex],axis = 0)
        SD.append(np.mean(baselineSDs))
        #baseline subtraced array
        baselineArray = np.full((rows,c),baselineMedians[XSG_num])
        bsArray.append(dataArray[XSG_num]-baselineArray)

    eventDetectFactor = 6 # 6 SD, note, this may not match BCJ's analysis, check

    if numMap ==1:
        averageMap = bsArray[0]
        averageBaseline = baselineMedians[0]
    else: 
        #preallocate variables 
        traceStack = np.empty((bsArray[0].shape[0],c,numMap)) #specify size explicitly here (4000,numtraces,numMaps) 
        baselineStack = np.empty((c,numMap)) #specify size so you can index into it to assign values
        eventDetectFlag = np.empty((c,numMap)) #explicite size preallocation 

        for j in range(numMap):
            for i in range(NumTraces):
                a,b = np.nonzero(standardPattern==i+1) #find 2d index of trace i in the standard map 
                #pattern. Find the number (acquisition sequence) at position (a,b) in map pattern j 
                #and extract cooresponding trace. From IG: note: the standard pattern goes from 1-192, 
                # since it is indexed using Matlab 1 indexing, so add +1 to the iterator! 
                trace = bsArray[j][:,pattern[j][a,b]-1] #IG: note, -1 to index back into bsArray, which is 
                #indexed using 0th indexing! 
                trace = np.reshape(trace,(bsArray[0].shape[0])) #reshape from (4000,1) to (4000,) array 
                #to be able to broadcast into trace stack variable 
                traceStack[:,i,j]=trace #made a list of all the traces. Possible make into 3d array instead if needed
                baseline = baselineMedians[j][i]
                baselineStack[i,j] = baseline
                eventArr = trace[stimOnInd:peakSearchEnd]
                # valley is defined by min within seach window - mean of 5ms baseline before stim
                valleyofEvent = np.min(eventArr)-np.mean(trace[stimOnInd-50:stimOnInd])
                eventDetectFlag[i,j] = np.abs(valleyofEvent)>np.abs(SD[j]*eventDetectFactor)

    #wipe out spontaneous events: detected in less than half of iterations 
    #sum the event-detection logic numbers of each tarce by row (aka across maps)
    #possibly wiping out way more responses than necessary, IG: look into this 
    
    sumofEventDetectFlag = np.sum(eventDetectFlag,axis=1)

    for i in range(NumTraces):
        # if sum of event detection flags across all maps is less than half number of maps, 
        #and is not 0, there is random events 
        if sumofEventDetectFlag[i] <= numMap/2 and sumofEventDetectFlag[i] != 0:
            #for each trace, find indices of map iterations that have spontaneous events.
            #mapIDsofSpontaneousEvent is an array of interations numbers of maps with a spontaneous
            #event during that trace (stim position)
            mapIDofSpontaneousEvent = np.nonzero(eventDetectFlag[i,:]==1)
            NoSponEvents = np.shape(mapIDofSpontaneousEvent)[0] #check to make sure this is correct! 
            for x in range(NoSponEvents): #x is index in mapIDsofSpontaneousEvent
                mapID = mapIDofSpontaneousEvent[x] #mapID is one of the map iDs that have spontaneous eent.
                #replace spon event in the detection window with zero! 
                traceStack[stimOnInd:peakSearchEnd,i,mapID] = 0

    #plot traces as maps again after wiping out spontaneous events 

    plotFlag2 = 1 
    if plotFlag2 ==1:
        for j in range(numMap):
            #plotting code, colors one trial here 
            something = 1

    #avereging: 
    #note: traceStack alignes with Matlab. 
    averageMap = np.mean(traceStack,axis = 2)
    averageBaseline = np.mean(baselineStack,axis=1)

    #dict to save for the rest of the maps 
    class NewClass(object): pass
    map_dict = NewClass()

    map_dict.experimentNumber = loadedData[0]['experimentNumber']
    map_dict.pattern = pattern
    map_dict.standardPattern = standardPattern
    map_dict.Racc = Racc
    map_dict.Rtotal = Rtotal
    map_dict.Rinput_ratio = Rinput_ratio
    map_dict.soma1Coordinates = loadedData[0]['soma1Coordinates']
    map_dict.ampNumber = 1
    map_dict.xSpacing = loadedData[0]['xSpacing']
    map_dict.ySpacing = loadedData[0]['ySpacing']
    map_dict.xPatternOffset = loadedData[0]['xPatternOffset']
    map_dict.yPatternOffset = loadedData[0]['yPatternOffset']
    map_dict.laserPower = laserPower
    map_dict.traces = dataArray
    map_dict.averageMap = averageMap
    map_dict.baseline = averageBaseline
    map_dict.samplingRate = loadedData[0]['sampleRate']
    #map_dict.experimentNumber = loadedData[0]['experimentNumber']
    map_dict.spatialRotation = loadedData[0]['spatialRotation']
    map_dict.horizontalVideoDistance = loadedData[0]['horizontalVideoDistance']
    map_dict.verticalVideoDistance = loadedData[0]['verticalVideoDistance']
    map_dict = map_dict.__dict__
    #print('soma coords after average: '+str(map_dict['soma1Coordinates']))
    #QC PLOTTING:
    #R access 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(map_dict['Racc'],'o');
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    fig.suptitle('Racc',color = 'white')
    
    #Input ratio
    fig2 = plt.figure()
    ax = fig2.add_subplot(111)
    ax.plot(map_dict['Rinput_ratio']);
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    fig2.suptitle('Rinput_ratio',color = 'white')
    return map_dict    

def analyzeSCRACMmap(numMaps,save=False): 
    """
    Functioned used to import, average, clean, and store data from a single sCRACM experiment (one cell).
    Remember, when saving the analysis.csv files, be sure that your naming convention matches the CellID in 
    your database, as these files need to be callable for future analysis. Store all analysis files in single directory! 
    """
    map_dict =averageMaps(numMaps)

    rows,numTraces = map_dict['averageMap'].shape
    sr = map_dict['samplingRate']
    responsOnsetEnd = 1501 #changed from 1250, some distance response may be that slow (YC, LSPS)
    baselineStartIndex = 901 #IG: changed from 950 to 900 to increase baseline, same as BJC. #IG (note, changed in python indexing))
    baselineEndIndex = 1001
    stimOnInd = 1001
    synEndInd = 1751 #integrate responses 75ms after stim
    chargeEndInd = 1751
    baselineMedians = np.median(map_dict['averageMap'][baselineStartIndex:baselineEndIndex,:],axis=0)
    baselineSDs = np.std(map_dict['averageMap'][baselineStartIndex:baselineEndIndex,:],axis=0)
    SD = np.mean(baselineSDs)
    print(f'Baseline SD: {SD} ')

    #Thresholds, these are vectors remember 
    dirLevel = 6 #note, 15 has also been used when analyzing sCRACM data, see BJC's BCS
    dirNegThresh = np.reshape(baselineMedians - dirLevel*SD, (1,numTraces))
    meanThresh = np.mean(dirNegThresh)
    print(f'Average response threshold is {meanThresh}')
    threshold = np.repeat(dirNegThresh,rows,axis=0)
    responseCounter = 0
    noResponseCounter = 0
    avgTraceArr = map_dict['averageMap'] #copy of map matrix, should have size [rows,cols]

    #preallocate variables
    r,c = map_dict['pattern'][0].shape
    integral = np.empty((1,numTraces))
    minimum = np.empty((1,numTraces))
    minOnset = np.empty((1,numTraces))
    onset = np.empty((1,numTraces))
    onset90 = np.empty((1,numTraces))
    onset10 = np.empty((1,numTraces))
    riseTime1090 = np.empty((1,numTraces))
    traceMean = np.empty((1,numTraces))
    mapOnset = np.empty((r,c))
    mapMinOnset = np.empty((r,c))
    mapMin = np.empty((r,c))
    mapMean = np.empty((r,c))
    mapIntegral = np.empty((r,c))

    #note: average maps are different between Matlab and Python: difference is 
    # small, rounding point error? present in the averaging somewhere. Look into it later 
    for i in range(numTraces):
        integral[:,i] = np.trapz(map_dict['averageMap'][stimOnInd:chargeEndInd,i])/sr
        logical = map_dict['averageMap'][stimOnInd:responsOnsetEnd,i]<threshold[stimOnInd:responsOnsetEnd,i]
        # logical is vector of 1,0, 1 when comparison holds true 
        if not any(logical): #when logical is all false (aka, no response for that trace), this statement returns True
            noResponseCounter = noResponseCounter +1
            onset[:,i] = float('nan')
            minimum[:,i]=float('nan')
            minOnset[:,i]=float('nan')
            onset90[:,i]=float('nan')
            onset10[:,i]=float('nan')
            riseTime1090[:,i]=float('nan')
            traceMean[:,i]=0
        else:
            responseCounter = responseCounter+1
            onset[:,i] = np.argmax(logical>0)/sr #find indecies where value exceeds 0 (first true), divide by sample rate 
            #find min and min index for trace 
            Minimum,MinOnset = np.min(map_dict['averageMap'][stimOnInd:synEndInd,i]), np.argmin(map_dict['averageMap'][stimOnInd:synEndInd,i])
            # need try, except statments for onset 90 and 10. Not worrying about this right now
            #
            #
            #
            #
            #

            minimum[:,i] = Minimum
            minOnset[:,i] = MinOnset/sr
            traceMean[:,i] = np.mean(map_dict['averageMap'][stimOnInd:synEndInd,i])

        #put values in the position of LSPS/sCRACM map
        a,b = np.nonzero(map_dict['standardPattern']==i+1)
        mapOnset[a[0],b[0]] = onset[:,i] 
        mapMinOnset[a[0],b[0]] = minOnset[:,i]
        mapMin[a[0],b[0]] = minimum[:,i]
        mapMean[a[0],b[0]] = traceMean[:,i]
        mapIntegral[a[0],b[0]] = integral[:,i]

    print(f'{responseCounter} sites have response, {noResponseCounter} sites do not')

    # transform soma location from coordinates saved from image capture to relative postion on map
    #check to catch for empty soma designation 
    if np.size(map_dict['soma1Coordinates'])>0: #if you forgot to mark soma, then size will be 0
        print('soma coords designated online')
    else:
        print('Soma position was not selected during experiment, please enter manual points and fix post-hoc')
        x=input("Enter temporary soma X coord (if unsure use 0)")
        y=input('Enter temporary soma y coord (if unsure use 0)')
        map_dict['soma1Coordinates'] = [int(x),int(y)]
    newX,newY= somaPositionTransformer(map_dict['soma1Coordinates'][0],map_dict['soma1Coordinates'][1],map_dict['spatialRotation'],map_dict['xPatternOffset'],map_dict['yPatternOffset'])
    newSomaCoords = [newX,newY]
    #save data to new_dict: 
    class NewClass(object): pass
    map_analysis_dict = NewClass()
    
    map_analysis_dict.experimentNumber = map_dict['experimentNumber']
    map_analysis_dict.amplifierNumber = map_dict['ampNumber']
    map_analysis_dict.negativeControl = float('nan')
    #map_analysis_dict.soma1Coordinates = map_dict['soma1Coordinates']
    map_analysis_dict.soma1Coordinates = newSomaCoords
    map_analysis_dict.horizontalVideoDistance = map_dict['horizontalVideoDistance']
    map_analysis_dict.verticalVideoDistance = map_dict['verticalVideoDistance']
    map_analysis_dict.distanceToPia = float('nan') #user defined after measuring
    map_analysis_dict.cortexThickness = float('nan') #user defined after measuring 
    map_analysis_dict.layer1Row = float('nan') #user defined after measuring 
    map_analysis_dict.laserPower = map_dict['laserPower']
    map_analysis_dict.pulseWidth = float('nan') #not applicable for now
    map_analysis_dict.xSpacing = 50
    map_analysis_dict.ySpacing = 50
    map_analysis_dict.numberOfMaps = numMaps
    map_analysis_dict.xPatternOffset = map_dict['xPatternOffset']
    map_analysis_dict.yPatternOffset = map_dict['yPatternOffset']
    map_analysis_dict.spatialRotation = map_dict['spatialRotation']
    map_analysis_dict.pattern = map_dict['pattern']
    map_analysis_dict.onset = mapOnset
    map_analysis_dict.minOnset = mapMinOnset
    map_analysis_dict.minimum = mapMin
    map_analysis_dict.mean = mapMean
    map_analysis_dict.integral = mapIntegral
    
    map_analysis_dict = map_analysis_dict.__dict__
    
   # plot out maps (not working right now when loading the function!)
    analysisPlots(map_analysis_dict,map_analysis_dict['amplifierNumber'])
    
    # meausure soma depth and cortical thickness.  
    previous_point = None
    img_path=measure_points()
    plotMappingArray(img_path,map_analysis_dict)
    
    # save data in .csv file. Be consistent about naming! Should be: CEllID_analysis.csv
    if save ==True: 
        fileName=input("File name (include file extension ie .csv): ")
        pathToSave = path_ephys_data+fileName 
        field_names=map_analysis_dict.keys() #variable names 
        with open(pathToSave,'w') as csvfile:
            writer = csv.DictWriter(csvfile,fieldnames = field_names)
            writer.writeheader()
            writer.writerows([map_analysis_dict])
    return map_analysis_dict