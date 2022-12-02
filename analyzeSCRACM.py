"""
Averaging and plotting sCRACM data 

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
    - loadXSG: helper function called by averageMaps, takes amplifer number. Will open a matlab XSG file from a single mapping sweep, and saves salient data in python dictionary.
    - averageMaps: helper function used to average multiple sweeps from a single experiment into one map. 
    - load_sCRACMdatabase: loads databse with meta data, includes dates and cell IDs
    - annotate_layers: gives layer assignments to cells based on cortical depth
    - sCRACM_cellShifter: shifts sCRACM maps up or down (adds rows/cols) to alignm by soma
    - sCRACM_cellShifterpia: shifts sCRACM maps up or down (adds rows/cols) to alignm by pia
    - mapAverage: takes stack of maps and averages into single averaged map
    - get_cells_to_average: makes list of file locations for cells 
    - average_map_stack: returns stack of aligned maps, based on input dependent filters
Execution/plotting functions: 
    - analyzeSCRACMmap: used to average sweeps from a single experiment into one map, optioal save flag to save .csv file for cell, used by all later functions. 
    - average_map: plots average maps (soma and pia aligned) based in input dependent filters
    - page_through_maps: returns interactive plot, allowing user to scroll through all maps 
    making up a particular average map 
    - collapse_map: collapses all cells from an input source on the x-axis and plots 
    heatmaps of each cell, arranged from superfical to deep layer cells 
    - bead_comparison: plots paired charts showing total integration of bead positive cells 
    against negative controls per input source, per layer. Returns ttest stats for each layer

Citations and such: 
Much of the sCRACM map averaing has been adapted, with permission, from Petreanu et al., 2009: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2745650/  (original code writen in Matlab)

"""

#imports 
import numpy as np
import pandas as pd
import scipy as sc
from scipy import optimize, stats,io
import math
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
import sys
from tqdm.notebook import tqdm
from plotly import graph_objects as go
from tkinter import filedialog
import cv2
import csv


#define path to analysis files. Hard coded for now 
analysis_path = '/Users/iangingerich/Dropbox (OHSU)/Ian_Tianyi_shared/Databases/Analyzed_cells/'

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
    #print('new coords inside transform function:'+str(somaXnew)+','+str(somaYnew))
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
        img.axes.tick_params(color = 'white',labelcolor = 'white')
        ax[n].plot(somax,somay,'^',markersize=10, color='b');
        cbaxes = inset_axes(ax[n],width = "3%", height = "75%",loc = 'lower right')
        color_bar = plt.colorbar(img,cax=cbaxes);
        color_bar.ax.yaxis.set_tick_params(color = 'white')
        plt.setp(plt.getp(color_bar.ax.axes,'yticklabels'),color='white')
    ax[0].set_title('Mean EPSC',color = 'white')
    ax[1].set_title('Minimum',color = 'white')
    ax[2].set_title('Onset',color = 'white')
    ax[3].set_title('Charge',color = 'white')    

#function to store pixel values when we do a left mouse click 
def mousePoints(event,x,y,flags,params):
    """function to capture coordinates of mouse clicks on a .tif image. Needs a driver function
    to be called."""
    global previous_point
    #left button click 
    if event == cv2.EVENT_LBUTTONDOWN:
        if previous_point:
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
    img_path = filedialog.askopenfilename();
    img=cv2.imread(img_path)
    

    
    cv2.imshow('image', img)
    cv2.setMouseCallback('image',mousePoints,img)
    cv2.waitKey(0)
    cv2.destroyAllWindows()
    cv2.waitKey(1);

def loadXSG(ampNumber):
    """ This helper function is called by averageMaps() and takes the amplifier number from an experiment, 
     and opens a pop up dialogue for you to select a specific .XSG map file (Matlab data file). It will parse the 
     selected file, and return a dictionary of parameters used for sCRACM and LSPS map averaging.
    """
    #1: open dialogue box
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
    map_analysis_dict.onset = mapOnset
    map_analysis_dict.minOnset = mapMinOnset
    map_analysis_dict.minimum = mapMin
    map_analysis_dict.mean = mapMean
    map_analysis_dict.integral = mapIntegral
    
    map_analysis_dict = map_analysis_dict.__dict__
    
   # plot out maps 
    analysisPlots(map_analysis_dict,map_analysis_dict['amplifierNumber'])
    
    # meausure soma depth and cortical thickness.  
    previous_point = None
    measure_points()
    
    # save data in .csv file. Be consistent about naming! Should be: CEllID_analysis.csv
    if save ==True: 
        root_dir = '/Users/iangingerich/Desktop/sCRACM_analysis/sample_data/' #hard coded for now: change to proper root dir! 
        fileName=input("File name (include file extension ie .csv): ")
        pathToSave = root_dir+fileName
        field_names=map_analysis_dict.keys() #variable names 
        with open(pathToSave,'w') as csvfile:
            writer = csv.DictWriter(csvfile,fieldnames = field_names)
            writer.writeheader()
            writer.writerows([map_analysis_dict])
    return map_analysis_dict


def load_sCRACMdatabase():
    """"
    Helper function, will load database will all meta data. File path is hard coded, change to 
    appropriate path, where database file lives. 
    """
    sCRACMdf = pd.read_csv('/Users/iangingerich/Dropbox (OHSU)/Ian_Tianyi_shared/Databases/sCRACM_database1.csv')
    return(sCRACMdf)

def annotate_layers(df):
    """
    Will annotate loaded DataFrame with layer identies of each cell. Layers are based on normalized 
    soma depth. Database must contain soma depth and cortical thickness information to be annotated. Use load_sCRACMdatabase function to load in df.
    """
    df['ratio_thickness'] = df['piadistance']/df['cortical thickness']  #find the normalized cortical depth of the soma. #soma depth and cortical thickness are defined during data processing 
    labels = ['L1','L2','L3','L5','L6','Claustrum'] #label. Note, this is for the agranular insula (AI), where L4 is absent 
    bins = [0,0.13,0.25,0.43,0.68,0.87,1] #normalized thickness bins are for AI only, based on 3 nissl brain layer annotations. 
    #note: layer annotations are variable depending on dorsal/ventral axis. For more acurate layer identities, map cells to insular flatmap first, and query thickness based on coordinates! 

    df['layer_bin']= pd.cut(df['ratio_thickness'],bins=bins,labels = labels) #annotate the df using the labels and bins variables 
    
# Helper functions listed first 
def  sCRACM_cellShifter(map_n,xy,spacing,trueRowNumber):
    """
    Helper function used when averaging multiple sCRACM maps together. used when averaging maps of different 
    sizes. Will shift maps up or down to align individual maps to each other. Used when aligning individual maps
    to the soma of each recording. 
    inputs: 
    map_n: single average map from a LSPS or sCRACM experiment 
    xy: soma coordinates to base shifting on 
    trueRowNumber: true number of rows the map had when recording the cell 
    """
    
    [r,c] = np.shape(map_n) #row and columns for the map you are shifting 
    newmap = np.zeros([2*r,2*c]) # new map size, double the size of the original map
    #go through maps one at a time
    #calculate how many rows and columns to shift by. x for rows, y for columns 

    shiftRows = np.round(xy[1]/spacing)+(r-trueRowNumber)/2 #calulate the number of rows to shift by, x-coordinate of the soma, divided by the map spacing (usually 50um) added to the difference 
    # between map row and the true row number divided by 2. 
    shiftCols = -np.round(xy[0]/spacing) # -y coord divided by the map spacing (50um usually)
   
   #place the original map onto the new map, with shift applied  
    newmap[int(1 + np.round(r/2)+shiftRows):int(r+np.round(r/2)+shiftRows+1),int(1+np.round(c/2)+shiftCols):int(c+np.round(c/2)+shiftCols)+1] = map_n 
    return newmap

# sCRACM cell shifter for pia alignment
def  sCRACM_cellShifterPia(map_n,xy,spacing,trueRowNumber):
    """
    Helper function used when averaging multiple sCRACM maps together. Used when aligning to pia 
    """

    [r,c] = np.shape(map_n)
    newmap = np.zeros([2*r,2*c])
    #go through maps one at a time
    #calculate how man rows and columns to shift by. x for rows, y for columns 

    shiftRows = np.round(xy[1]/spacing)-(trueRowNumber/2) #note: this line is the differnce between the above function. The true row number is subraced from the coordinate instead of diff and added.
    shiftCols = -np.round(xy[0]/spacing)
    #place the original map onto the new map, with shift 
    newmap[int(1 + np.round(r/2)+shiftRows):int(r+np.round(r/2)+shiftRows+1),int(1+np.round(c/2)+shiftCols):int(c+np.round(c/2)+shiftCols+1)] = map_n
    return newmap


# mapAverage helper function: takes stack of maps to be averaged, after shifting, normalization, and col/row adding, and averages together into a single 24 by 12 array
def mapAverage(mapStack):
    """
    Helper function called when plotting average sCRACM maps. Will take (x,y,n) array of maps, average by last dim,
    returning single (x,y) array, which represents the average map. 
    inputs: 
    mapStack: a stack of maps to average, will be (n,24,12) in size with n being the number of maps you are averaging 
    """
    #mapStack = maps
    [planes,rows,cols] = np.shape(mapStack) #take the dimensions of the stack
    mapAvg = np.zeros([rows,cols]) #make a new array for the average map
    for r in range(0,rows): #iterate through the rows and columns of the mapStack, and take the mean at each point, return as mapAvg. 
        for c in range(0,cols):
            mapAvg[r,c] = np.mean(mapStack[:,r,c])
    return mapAvg 
# note: this above function could be imporved. Intead of looping, just use a vectorized approach such as: np.mean(mapStack, axis = 0). Taking the mean of the first axis 


def get_cells_to_average(input_source, layer,
                         by_layer_bin = True,
                         by_layer_online = False, 
                         by_Mtype = False,
                         bead_positive = False):
    """
    Helper function: contructs file paths for cells to be analyzed.
    
    :param str input_source: input source for map str(contraAI,Ent,Orbital,S1,PL,Piri,etc.)
    :para str layer: the layer you want to average across (L2,L3,L4,L5,L6)
    :param bool by_layer_bin: layer defined by relative cortical depth, as calculated by BCJ, default = True
    :param bool by_layer_online: layer defined by experimenter upon data collection, default = False
    :param bool by_Mtype: average based on cell type classification, if present, default = False
    :param bool bead_positive: will average only bead positive cells sitting in designated layer and receiving input from designated source, default = False
    """
    #load database 
    sCRACMdf = load_sCRACMdatabase()
    #calculate relative cortical thickenss
    sCRACMdf['ratio_thickness'] = sCRACMdf['piadistance']/sCRACMdf['cortical thickness']
    #annotate cells based on relative thickness. Note: see caviates with the layer designations in the documentation for annotate_layers()  
    annotate_layers(sCRACMdf) 
    
    # average based on relative cortical depth assignment, input source, make sure there is a map to analyze, and the map has responses. 
    if by_layer_bin == True:
        printout = 'by_layer_bin'
        cells_to_average_df = sCRACMdf[(sCRACMdf['InjectionSite']==input_source) & (sCRACMdf['MapForAnalysis']==True) & (sCRACMdf['Response']==True)
                                   &(sCRACMdf['layer_bin']== layer)]# & sCRACMdf['Bead Positive?']== bead_positive] 
        cell_IDs = analysis_path + cells_to_average_df['CellID'] + '_analysis.csv' #append together file names to load 

        #Nested IF statments: by layer bin first, and then add additional filters
        # average based on cellular M_type 
        if by_Mtype != False: #if by_Mtype is something other than false (default), then filter df my m_type as well. 
            printout = 'by_Mtype: '+str(by_Mtype)
            cells_to_average_df = cells_to_average_df[(cells_to_average_df['M_type_final']== by_Mtype)] #filter based on newest M_type assignments 
            cell_IDs = analysis_path + cells_to_average_df['CellID'] + '_analysis.csv'


      # if you want to get the file names of bead positive cells, bead positve= True during input

        elif bead_positive == True:
            printout = 'bead_positive'
            cells_to_average_df = cells_to_average_df[(cells_to_average_df['Injection Type']=='virus+bead')& (cells_to_average_df['Bead Positive?']== True)]

            cell_IDs = analysis_path + cells_to_average_df['CellID'] + '_analysis.csv'
            
        elif bead_positive == 'negative_control': #note: right now this only works for cells that have a negative control designated in their .m file. For full python code, the inport settings 
            # also needs to have a negative control variable to filtering 
            printout = 'negative_control'
            cells_to_average_df = cells_to_average_df[(cells_to_average_df['Injection Type']=='virus+bead')& (cells_to_average_df['Bead Positive?']== False)]
            cell_IDs = analysis_path + cells_to_average_df['CellID'] + '_analysis.csv'
            
    

    # average based on online designation of layer, not the average layer consensus. This should be the option if using this code to average cells recorded in places outside the AI 
    elif by_layer_online == True:
        printout = 'by_layer_online'
        cells_to_average_df = sCRACMdf[(sCRACMdf['InjectionSite']==input_source) & (sCRACMdf['MapForAnalysis']==True) & (sCRACMdf['Response']==True)
                                   &(sCRACMdf['Layer estimation']== layer)]#& sCRACMdf['Bead Positive?']==bead_positive] #filter based on mapfor analysis and input source
        cell_IDs = analysis_path + cells_to_average_df['CellID'] + '_analysis.csv'
    
        #nested IF statement: by online layer designation first, and then other filters!
        # average based on cellular M_type 
        if by_Mtype != False:
            printout = 'by_Mtype: '+str(by_Mtype)
            cells_to_average_df = cells_to_average_df[(cells_to_average_df['M_type_final']== by_Mtype)] #filter based on mapfor analysis and input source
            cell_IDs = analysis_path + cells_to_average_df['CellID'] + '_analysis.csv'


        
        elif bead_positive == True:
            printout = 'bead_positive'
            cells_to_average_df = cells_to_average_df[(cells_to_average_df['Injection Type']=='virus+bead')& (cells_to_average_df['Bead Positive?']== True)]

            cell_IDs = analysis_path + cells_to_average_df['CellID'] + '_analysis.csv'
            
        elif bead_positive == 'negative_control':
            printout = 'negative_control'
            cells_to_average_df = cells_to_average_df[(cells_to_average_df['Injection Type']=='virus+bead')& (cells_to_average_df['Bead Positive?']== False)]
            cell_IDs = analysis_path + cells_to_average_df['CellID'] + '_analysis.csv'
    
    #must select some kind of averaging method here! 
    else: 
        printout = 'please select variable to average on'
    print('Input Source: ' +input_source +' '+ layer+ ', Average method:' + printout)
    
    return cells_to_average_df, cell_IDs,printout


# Helper function that returns all map arrays of a particular input_source, layer, but does not avearge them. Used for the page through function! 
def average_map_stack(input_source, layer,
                         by_layer_bin = True,
                         by_layer_online = False, 
                         by_Mtype = False,
                         bead_positive = False):
    """
    apapted from Petreanu et al., 2009: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2745650/
    Loads cellular data into Python dicts (one per cell), extracts map and soma data, returns 
    a (x,y,n) array of all maps to be averaged. 
    
    :param str input_source: input source for map str(contraAI,Ent,Orbital,S1,PL,Piri,etc.)
    :para str layer: the layer you want to average across (L2,L3,L4,L5,L6)
    :param bool by_layer_bin: layer defined by relative cortical depth, as calculated by BCJ, default = True
    :param bool by_layer_online: layer defined by experimenter upon data collection, default = False
    :param bool by_Mtype: average based on cell type classification, if present, default = False
    :param bool bead_positive: will average only bead positive cells sitting in designated layer and receiving input from designated source, default = False
    """
    
    #first run function to get all cell files you want to use. Inputs are the same input to average_map_stack function 
    cells_to_average_df,cell_IDs,printout=get_cells_to_average(input_source=input_source, 
                         layer = layer,
                         by_layer_bin = by_layer_bin,
                         by_layer_online = by_layer_online, 
                         by_Mtype = by_Mtype,
                         bead_positive = bead_positive)
    
    numCells = np.shape(cell_IDs)[0] #total number of cells to grab 
   
    # initiate variables 
    mapSpacing = np.zeros((1,numCells)) #spacing
    layer1Row = np.zeros((1,numCells))
    XY = np.zeros([numCells,2]) #initializing array for soma coordinates
    piaXY = np.zeros([numCells,2]) #initializing array for soma coordinates, where y will be shifted for pia alignment 
    trueRowNumber_arr = []
    M = np.zeros([numCells,24,20]) #array where all the maps will be held during for loop 
    max_row_n = 24
    distanceToPia = []
    # start for loop to load in cells now
    for i,n in enumerate(cell_IDs):
        #load cells 
        cellA = m2p.mstruct2pydict(n) #use parcer to take .m files that have already been processed. Note: for full python, upade to load in datafiles of individual maps! 
        cellA_mean = cellA['mean'] #define the map 
        XY[i] = cellA['soma1Coordinates'] #save coordinates 
        [a,b] = cellA_mean.shape
        trueRowNumber = cellA_mean.shape[0]
        layer1Row[0,i] = cellA['layer1Row'][0]
        # add zeros to cells that are less that 24 rows! 
        if cellA_mean.shape[0] < max_row_n:
                add_arr = np.zeros(((24-cellA_mean.shape[0]),b)) #number of rows of zeros to add, note, change ending from 12 to 1, as we are adding after collapsing
                average_map_new = np.concatenate((cellA_mean,add_arr),axis = 0) # add rows of zeros here
                shift_unit = max_row_n - trueRowNumber
                trueRowNumber_arr += [trueRowNumber]
        
        # keep cells that are 24 rows as is 
        else: 
            average_map_new = cellA_mean
            trueRowNumber_arr += [trueRowNumber]
        addColumnFlag =1 # add columns to smaller maps, or to put averaged maps more in the center of image 
        if addColumnFlag: 
            maxColumnNumber = 20
            [row, col] = np.shape(average_map_new)
            if col < maxColumnNumber:
                placeHolderMatrix = np.zeros([row,maxColumnNumber])
                half_col_to_add = (maxColumnNumber-col)/2
                placeHolderMatrix[:,int(half_col_to_add):int(maxColumnNumber-half_col_to_add)] = average_map_new
                average_map_new = placeHolderMatrix
        #end if/else statement 
        #remove responses outside of tissue here 
        average_map_new[0:int((cellA['layer1Row'][0])-1),:] = 0 
        mapSpacing = cellA['xSpacing'] #assumes dx and dy are the same, which they are, both 50.  
        #next step normalizes maps to min value 
        minVal = np.min(average_map_new) #find min for each map
        M[i,:,:] = average_map_new/minVal 
        distanceToPia.append(cellA['distanceToPia'])
    
    rcAddOn = 0 #must be even interger 
    [r,c]=np.shape(M[1])
    X = XY[:,0]
    Y = -XY[:,1]
    dx = mapSpacing #50
    # shift maps for alignmnent 
    # preallocate newmap variable 
    ###############################################################
    # take out alignment options when making a function! 
    newmap = np.zeros([numCells,2*r,2*c]) 
    X_soma = np.zeros(np.shape(X))
    Y_soma = np.zeros(np.shape(Y))
    newmap_soma = np.zeros(([numCells,2*r,2*c]))

    for n in range(0,numCells):
        newmap_soma[n,:,:] = sCRACM_cellShifter(M[n,:,:],XY[n,:],dx,trueRowNumber_arr[n]) #use the cell shifter to align all maps by their soma coordinates (center on soma)
        X_soma[n] = X[n]-XY[n,0]
        Y_soma[n]= Y[n]+XY[n,1]
    maps_soma = newmap_soma[:,int((1+r/2)-rcAddOn/2):int((r+r/2)+rcAddOn/2+1),int((1+c/2)-rcAddOn/2):int((c+c/2)+rcAddOn/2+1)]
    
    
    return maps_soma, cells_to_average_df['CellID']

########## average sCRACM maps here, will call above helper function to average based on a specfic parameter! ###############

# fixed to pull .csv analysis files instead of .m files 
def average_map(input_source, layer,
                         by_layer_bin = True,
                         by_layer_online = False, 
                         by_Mtype = False,
                         bead_positive = False,
                         interp = False):
    
    """
    Plotting function, returns averaged heatmaps based on input source, and layer. Can average using 
    other optional catagorical parameters such as cell type classification, retrograde bead information, 
    and return an interpolated or pixelated heatmap. 
    
    pull cell IDs and file paths to average sCRACM maps
    :param str input_source: input source for map str(contraAI,Ent,Orbital,S1,PL,Piri,etc.)
    :para str layer: the layer you want to average across (L2,L3,L4,L5,L6)
    :param bool by_layer_bin: layer defined by relative cortical depth, as calculated by BCJ, default = True
    :param bool by_layer_online: layer defined by experimenter upon data collection, default = False
    :param bool by_Mtype: average based on cell type classification, if present, default = False
    :param bool bead_positive: will average only bead positive cells sitting in designated layer and receiving input from designated source, default = False
    """
    
    
    cells_to_average_df,cell_IDs,printout =get_cells_to_average(input_source = input_source,
                          layer=layer,
                          by_layer_bin = by_layer_bin,
                          by_layer_online=by_layer_online,
                          by_Mtype=by_Mtype,
                          bead_positive=bead_positive)
    numCells = np.shape(cell_IDs)[0]
    
    if interp == True:
        interpolation = 'gaussian'
    else:
        interpolation = None
    # initiate variables 
    mapSpacing = np.zeros((1,numCells)) #spacing
    layer1Row = np.zeros((1,numCells))
    XY = np.zeros([numCells,2]) #initializing array for soma coordinates
    piaXY = np.zeros([numCells,2]) #initializing array for soma coordinates, where y will be shifted for pia alignment 
    trueRowNumber_arr = []
    M = np.zeros([numCells,24,20]) #array where all the maps will be held during for loop 
    max_row_n = 24
    distanceToPia = []
    # start for loop to load in cells now 
    for i,n in enumerate(cell_IDs):
        #load cells 
        cellA =pd.read_csv(n) # read individual cell into df with nested cells 
        #format map data correctly 
        cellA_mean = np.fromstring(cellA['mean'].values[0].replace('\n','').replace('] [','\n').replace('[[','').replace(']]',''),sep=' ') # parse map array into correct datatype
        cols = 12 #hard code number of columns for now 
        rows = int(cellA_mean[0]/cols) #number of rows for the map 
        cellA_mean = np.reshape(cellA_mean,(rows,cols)) #reshape map to be (row,cols)
        #get soma coords in correct formate 
        soma1Coords=np.fromstring(cellA['soma1Coordinates'][0].replace('[','').replace(']',''),sep = ',')
        XY[i] = soma1Coords #save soma coordinates  
        [a,b] = cellA_mean.shape
        trueRowNumber = cellA_mean.shape[0]
        layer1Row[0,i] = cellA['layer1Row'][0]
        # add zeros to cells that are less that 24 rows! 
        if cellA_mean.shape[0] < max_row_n:
                add_arr = np.zeros(((24-cellA_mean.shape[0]),b)) #number of rows of zeros to add, note, change ending from 12 to 1, as we are adding after collapsing
                average_map_new = np.concatenate((cellA_mean,add_arr),axis = 0) # add rows of zeros here
                shift_unit = max_row_n - trueRowNumber
                trueRowNumber_arr += [trueRowNumber]
        
        # keep cells that are 24 rows as is 
        else: 
            average_map_new = cellA_mean
            trueRowNumber_arr += [trueRowNumber]
        addColumnFlag =1 # add columns to smaller maps, or to put averaged maps more in the center of image 
        if addColumnFlag: 
            maxColumnNumber = 20
            [row, col] = np.shape(average_map_new)
            if col < maxColumnNumber:
                placeHolderMatrix = np.zeros([row,maxColumnNumber])
                half_col_to_add = (maxColumnNumber-col)/2
                placeHolderMatrix[:,int(half_col_to_add):int(maxColumnNumber-half_col_to_add)] = average_map_new
                average_map_new = placeHolderMatrix    
        #end if/else statement 
        #remove responses outside of tissue here 
        average_map_new[0:int((cellA['layer1Row'][0])-1),:] = 0 
        #note: in original code, the mapSpacing variable is designated (n), need to save somewhere else? 
        mapSpacing = cellA['xSpacing'] #assumes dx and dy are the same, which they are, both 50.  
        #next step normalizes maps to min value 
        minVal = np.min(average_map_new) #find min for each map
        M[i,:,:] = average_map_new/minVal 
        distanceToPia.append(cellA['distanceToPia'])

        # end for loop now 
    rcAddOn = 0 #must be even interger 
    [r,c]=np.shape(M[1])
    X = XY[:,0]
    Y = -XY[:,1]
    dx = mapSpacing #50
    # shift maps for alignmnent 
    # preallocate newmap variable 
    ###############################################################
    newmap = np.zeros([numCells,2*r,2*c]) 
    #initialize figure 
    fig  = plt.figure(figsize=(30,30))
    fig.suptitle('Input: '+input_source+', Averaged '+ printout + ' (n='+str(numCells)+')',color = 'black',fontsize = 20)
    ax1 = plt.subplot2grid((30,30),(2,0),colspan =13,rowspan =13) #average map subplot
    ax2 = plt.subplot2grid((30,30),(0,0),colspan = 13,rowspan = 2) #subplot for intensity by column
    ax3 = plt.subplot2grid((30,30),(2,13),colspan =2,rowspan = 13) #subplot for intensity by row 
    # for pia alignment 
    ax4 = plt.subplot2grid((30,30),(2,15),colspan =13,rowspan =13) #average map subplot
    ax5 = plt.subplot2grid((30,30),(0,15),colspan = 13,rowspan = 2) #subplot for intensity by column
    ax6 = plt.subplot2grid((30,30),(2,28),colspan = 2,rowspan = 13) #subplot for intensity by row 
    
    
    # preallocated size of maps for soma alignment
    X_soma = np.zeros(np.shape(X))
    Y_soma = np.zeros(np.shape(Y))
    newmap_soma = np.zeros(([numCells,2*r,2*c]))
    # align by soma 
    for n in range(0,numCells):
        newmap_soma[n,:,:] = sCRACM_cellShifter(M[n,:,:],XY[n,:],dx,trueRowNumber_arr[n])
        X_soma[n] = X[n]-XY[n,0]
        Y_soma[n]= Y[n]+XY[n,1]
    maps_soma = newmap_soma[:,int((1+r/2)-rcAddOn/2):int((r+r/2)+rcAddOn/2+1),int((1+c/2)-rcAddOn/2):int((c+c/2)+rcAddOn/2+1)]
    
    rowSumByCell = np.sum(maps_soma, axis =2)[...,None] # summation by row of each cell, gives shape (numCells,24,1)
    colSumByCell = np.transpose(np.sum(maps_soma, axis =1)[...,None],(0,2,1)) #summation by column, gives shape (numcell, 1, 12)
    # calculate standard error of row summations 
    rowSumStdByCell = np.std(rowSumByCell, axis =0)
    colSumStdByCell = np.std(colSumByCell, axis = 0)
    NofCells,_,_ = np.shape(maps_soma)
    rowSumSEByCell = rowSumStdByCell/np.sqrt(NofCells)
    colSumSEByCell = colSumStdByCell/np.sqrt(NofCells)
    
    avgMap_soma = mapAverage(maps_soma)
    rowSum = np.abs(np.sum(avgMap_soma, axis = 1))[...,None]
    rowSum_1 = np.abs(np.sum(avgMap_soma, axis = 1))
    colSum = np.abs(np.sum(avgMap_soma, axis = 0))[...,None].T
    colSum_1 = np.abs(np.sum(avgMap_soma, axis = 0))
    #setting up figure limits and such 
    [r,c] = np.shape(avgMap_soma)
    dx = mapSpacing #50
    dy = mapSpacing #50 
    xdata = np.arange(0,(c*dx),dx) - int((c-1) * dx/2)
    ydata = np.arange(0,(r*dx),dx) - int((r-1) * dx/2)
    xmax = (c+rcAddOn)*dx/2   #xmax = (map size in x direction * map spacing) - (mapsize in x)*50/2 * note: I removed the c-1 and r-1 to fix the problem with the soma plotting
    xmin = -xmax # x data is centered on 0, so min and max are the same here, with signs reversed
    ymax = (r+rcAddOn)*dx/2 #ydata also centered on zero
    ymin = -ymax
    img = ax1.imshow(avgMap_soma,aspect =
                  'auto', cmap = 'magma',extent = [xmin,xmax,ymax,ymin],interpolation = interpolation);

    #colorbar code 
    cbaxes = inset_axes(ax1,width = "3%", height = "75%",loc = 'lower left')
    color_bar = plt.colorbar(img,cax=cbaxes);
    #color_bar = plt.colorbar(img,orientation = 'vertical')
    cbytick_obj = plt.getp(color_bar.ax.axes, 'yticklabels');     #getting color bar tick object          
    plt.setp(cbytick_obj, color='w'); #setting colorbar ticks to white

    # plot soma position for soma alignment
    for n in range(0,np.shape(X_soma)[0]):
        ax1.plot(X_soma[n],Y_soma[n], 'k^',markersize=8, color='b');
    # plot for input intensity by row 
    rowNum,_ = np.shape(rowSum)
    yarr = np.linspace(0.5, (rowNum-0.5), num = rowNum)
    yarr = yarr[::-1] #flip 
    
    ax3.errorbar(rowSum.reshape((-1)), yarr, xerr= rowSumSEByCell.reshape((-1)), fmt = '-o');
    ax3.set_title('Relative intensity', color = 'black');
   
    #plot for input intenstiy by column 
    _,colNum = np.shape(colSum)
    yarr = np.linspace(0.5,(colNum-0.5),num = colNum)
    ax2.errorbar(yarr,colSum.T.reshape((-1)),yerr=colSumSEByCell.T.reshape((-1)),fmt = '-o');
    ax2.set_title('Soma alignment', color = 'black', fontsize = 20)
    
  #####################################################  
# now is the plotting code for pia alignment! 
    # initialize X and Y coordinates for pia 
    X_pia = np.zeros(np.shape(X)[0])
    Y_pia = np.zeros(np.shape(Y)[0])
    newmap_pia = np.zeros(([numCells,2*r,2*c]))
    for n in range(0,numCells):
        piaXY[n,:] = [XY[n,0],distanceToPia[n]+ XY[n,1]]
        newmap_pia[n,:,:] = sCRACM_cellShifterPia(M[n,:,:],piaXY[n,:],dx,trueRowNumber_arr[n])
        X_pia[n] = X[n]-piaXY[n,0]
        Y_pia[n] = Y[n]+piaXY[n,1]-(r/2)*dx
    maps_pia = newmap_pia[:,int((1+r/2)-rcAddOn/2):int((r+r/2)+rcAddOn/2+1),int((1+c/2)-rcAddOn/2):int((c+c/2)+rcAddOn/2+1)]
    
    #############
    # average map code here 
    rowSumByCell = np.sum(maps_pia, axis =2)[...,None] # summation by row of each cell, gives shape (numCells,24,1)
    colSumByCell = np.transpose(np.sum(maps_pia, axis =1)[...,None],(0,2,1)) #summation by column, gives shape (numcell, 1, 12)
    # calculate standard error of row summations 
    rowSumStdByCell = np.std(rowSumByCell, axis =0)
    colSumStdByCell = np.std(colSumByCell, axis = 0)
    NofCells,_,_ = np.shape(maps_pia)
    rowSumSEByCell = rowSumStdByCell/np.sqrt(NofCells)
    colSumSEByCell = colSumStdByCell/np.sqrt(NofCells)

    #averaging step using helper function mapAverage 
    avgMap_pia = mapAverage(maps_pia)
    rowSum = np.abs(np.sum(avgMap_pia, axis = 1))[...,None]
    colSum = np.abs(np.sum(avgMap_pia, axis = 0))[...,None].T
    #setting up figure limits and such 
    [r,c] = np.shape(avgMap_pia)
    dx = mapSpacing #50
    dy = mapSpacing #50 
    xdata = np.arange(0,(c*dx),dx) - int((c-1) * dx/2)
    ydata = np.arange(0,(r*dx),dx) - int((r-1) * dx/2)
    xmax = (c+rcAddOn)*dx/2   #xmax = (map size in x direction * map spacing) - (mapsize in x)*50/2 * note: I removed the c-1 and r-1 to fix the problem with the soma plotting
    xmin = -xmax # x data is centered on 0, so min and max are the same here, with signs reversed
    ymax = (r+rcAddOn)*dx/2 #ydata also centered on zero
    ymin = -ymax

    img = ax4.imshow(avgMap_pia,aspect =
                  'auto', cmap = 'magma',extent = [xmin,xmax,ymax,ymin],interpolation = interpolation);

    #colorbar code 
    cbaxes = inset_axes(ax4,width = "3%", height = "75%",loc = 'lower left')
    color_bar = plt.colorbar(img,cax=cbaxes);
    cbytick_obj = plt.getp(color_bar.ax.axes, 'yticklabels');     #getting color bar tick object          
    plt.setp(cbytick_obj, color='w'); #setting colorbar ticks to white

    # plot soma positions 
    for n in range(0,np.shape(X_pia)[0]):
        ax4.plot(X_pia[n],Y_pia[n], 'k^',markersize=8, color = 'b');
    
    # plot for input intensity by row 
    rowNum,_ = np.shape(rowSum)
    yarr = np.linspace(0.5, (rowNum-0.5), num = rowNum)
    yarr = yarr[::-1] #flip 
    ax6.errorbar(rowSum.reshape((-1)), yarr, xerr= rowSumSEByCell.reshape((-1)), fmt = '-o');
    ax6.set_title('Relative intensity', color = 'black');

    #plot for input intenstiy by column 
    _,colNum = np.shape(colSum)
    yarr = np.linspace(0.5,(colNum-0.5),num = colNum)
    ax5.errorbar(yarr,colSum.T.reshape((-1)),yerr=colSumSEByCell.T.reshape((-1)),fmt = '-o');
    ax5.set_title('Pia alignment', color = 'black', fontsize = 20);
    #fig.savefig('/Users/iangingerich/Desktop/Poster_materials/'+input_source+layer+'.svg',format = 'svg',dpi = 300)

def page_through_maps(input_source, layer,
                         by_layer_bin = True,
                         by_layer_online = False, 
                         by_Mtype = False,
                         bead_positive = False):
    
    """"
    Will plot all cell input maps for a certain group based on categorical variables beloning 
    to the cell type. Outputs a scrollable figure utilizing Plotly. Note: this is compatable with 
    Jupyder Notebook only, for plotting on other platforms, consult Plotly's API. 
    
    :param str input_source: source of excitation (contraAI, Ent, Orbital, S1, PL, Piri)
    :param str layer: layer based in relative cortical depth of recorded cell 
    :param bool by_layer_bin: layer defined by relative cortical depth, as calculated by BCJ, default = True
    :param bool by_layer_online: layer defined by experimenter upon data collection, default = False
    :param bool by_Mtype: average based on cell type classification, if present, default = False
    :param bool bead_positive: will average only bead positive cells sitting in designated layer and receiving input from designated source, default = False
    """
    #get map data: 
    
    
    M,cell_id =average_map_stack(input_source,layer,
                                 by_layer_bin = by_layer_bin,
                                 by_layer_online = by_layer_online, 
                                 by_Mtype = by_Mtype,
                                 bead_positive = bead_positive)
    
    #initailize size 
    [e,r,c] = np.shape(M)
    dx = 50 #50
    dy = 50 #50 
    xdata = np.arange(0,(c*dx),dx) - int((c-1) * dx/2)
    ydata = np.arange(0,(r*dx),dx) - int((r-1) * dx/2)

    #flip horizontally, since y axis is inverted 
    M2 = np.zeros([e,24,20]) 
    for i in range(e):
        M2[i,:,:] = np.flip(M[i,:,:],axis = 0)

    #make frames 
    test_frames = [go.Frame(data = go.Heatmap(z = M2[i,:,:],x= xdata,y=ydata), name = cell) for i,cell in zip(range(e),cell_id.values)]


    #make sliders
    sliders =[{"steps": [{"args": [[cell],{"frame": {"duration": 0, "redraw": True},
                                                "mode": "immediate",},{'title':'CellID: '+str(cell)}],
                             "label": cell, "method": "animate",}
                            for frame,cell in zip(test_frames,cell_id.values)],},
             {"currentvalue": {'prefix':'CellID: '}}]



    fig = go.Figure(data = test_frames[0].data,frames = test_frames)
    fig.update_layout(sliders=sliders,yaxis = dict(scaleanchor = 'x'))
    fig.add_trace(go.Scatter(x = [0,0],y =[0,0],mode = 'markers', marker = dict(color=['#ff0000'])))
    fig.layout.title = 'Input: '+ input_source + ', Layer: '+layer+', n = ' + str(e)
    fig.show()

####################################################################################################################################################################
def collapse_map(input_source):
    
    """
    Plotting function. Will take all input maps from a particular input source, and collapse over 
    the x-axis. Returns a heatmap of all cells, sorted by layer and cortical depth. 
    """
   # Read in sCRACM database. Make sure it is up to date! 
    sCRACMdf = load_sCRACMdatabase() 
    contraidf = sCRACMdf[(sCRACMdf['InjectionSite']==input_source) & (sCRACMdf['MapForAnalysis']==True) & sCRACMdf['Response']==True]
    #Build analysis files to load from analyzed cells folder
    #note: make sure analysis_path is set in the ini file, so the path is correct!  
    IDs = analysis_path + contraidf['CellID'] + '_analysis.csv' #should be reading in the analysis.csv files 
    cell_ID = contraidf['CellID']
    # convert cell IDs from DF to np.array to use for labels in the figure
    cell_ID_array = cell_ID.to_numpy()[...,None] 
    # initialize variable to store nested dict. struct. Could maybe change to dataframe in the future 
    max_row_n = 24 #used to add extra rows to smaller maps. Update to size of map instead of 24. 
    collapsed_arr = [[] for n in range(0,24)]
    soma_y_arr = []
    trueRowNumber_arr = []
    # change back to IDs to get full list of paths to files! 
    for n in IDs: 
        loaded_cell = pd.read_csv(n) #note: n is the path name, so passing it n will result in the file being read. 
        # parse map array into correct datatype
        average_map = np.fromstring(loaded_cell['mean'].values[0].replace('\n','').replace('] [','\n').replace('[[','').replace(']]',''),sep=' ') 
        #DANGER: cols variable is hard coded, if averaging larger maps, this need to be flexible! 
        cols = 12 #hard code number of columns for now 
        rows = int(average_map[0]/cols) #number of rows for the map 
        average_map = np.reshape(average_map,(rows,cols)) #reshape map to be (row,cols)
        # get soma coordinates and parse into correct format
        soma1Coords=np.fromstring(loaded_cell['soma1Coordinates'][0].replace('[','').replace(']',''),sep = ',')
        [a,b] = average_map.shape
        trueRowNumber = average_map.shape[0]
        trueRowNumber_arr = [trueRowNumber_arr,trueRowNumber]
        # if cell is less than 24 rows, add zeros, and adjust soma position, else, keep as is 
        if average_map.shape[0] < max_row_n:
            add_arr = np.zeros(((24-average_map.shape[0]),b)) #number of rows of zeros to add, note, change ending from 12 to 1, as we are adding after collapsing
            average_map_new = np.concatenate((average_map,add_arr),axis = 0) # add rows of zeros here
            shift_unit = max_row_n - trueRowNumber
            shift_val = (shift_unit*0.5)*50 # value to subtract from Y coordinate 
            soma_Y = (-1*soma1Coords[1])-shift_val #transformed Y-coordinate to plot on top of the heatmap
        else: 
            soma_Y = (-1*soma1Coords[1]) 
            average_map_new = average_map
        layer1Row = loaded_cell['layer1Row'][0]
        average_map_new[0:int((loaded_cell['layer1Row'][0])-1),:] = 0 #remove data over Pia by setting sites to 0, note, instead of 3 it should be the layer 1 number in the original dictinary, minus 1! 
        collapse_map = np.sum(average_map_new,axis = 1)[...,None] #last bit here is transposing 1D sum array into column. Rework code so collapse happens before adding zeros so larger maps can be averaged! 
        soma_y_arr = np.append(soma_y_arr,soma_Y) # read in soma coordinates after each round        
        collapsed_arr = np.concatenate((collapsed_arr,collapse_map), axis=1)
    # end for loop here 
    # sort soma by depth
    soma_y_sorted = np.sort(soma_y_arr,axis = -1)
    soma_order = np.argsort(soma_y_arr,axis = -1)
    sorted_cell_ID = cell_ID_array[soma_order,:]
    sorted_collapsed_arr = collapsed_arr[:,soma_order] #sort by column depth
    # plotting code below: 
    [r,c] = sorted_collapsed_arr.shape # find out the shape, aka map size you used (example 12 by 12)
    dx = 50 #both dx and dy should be 50 for sCRACM maps! 
    dy = 50
    xdata = np.arange(0,(c*dx),dx) - int((c-1) * dx/2)
    #extent tuples 
    xmax = (c*dx)-int((c) * dx/2)   #xmax = (map size in x direction * map spacing) - (mapsize in x)*50/2 * note: I removed the c-1 and r-1 to fix the problem with the soma plotting
    xmin = -xmax # x data is centered on 0, so min and max are the same here, with signs reversed
    ymax = (r*dy)-int((r) * dy/2) #ydata also centered on zero
    ymin = -ymax
    _,ax = plt.subplots(figsize = (24,6)) #create plot, use underscore to exclude getting first output from subplots
    img = ax.imshow(sorted_collapsed_arr,aspect =
              'auto', cmap = 'plasma',extent = [xmin,xmax,ymax,ymin]);
    
    ax.plot(xdata,soma_y_sorted,'k^',markersize=8);
    ax.set_xticks(xdata)
    ax.set_xticklabels(sorted_cell_ID,rotation =90, color='b')
    ytick = [ymin,ymax]
    ax.set_yticks(ytick)
    yticklab = ['pia','CC']
    ax.set_yticklabels(yticklab,color='b')
    color_bar = plt.colorbar(img, ax=ax)
    cbytick_obj = plt.getp(color_bar.ax.axes, 'yticklabels')     #getting color bar tick object          
    plt.setp(cbytick_obj, color='b') #setting colorbar ticks to white
    plt.title(input_source+': Total current collapsed on x axis', color = 'b')
    
    # below for normalized plots 
    _,ax = plt.subplots(figsize = (24,6)) #create plot, use underscore to exclude getting first output from subplots
    img = ax.imshow(sorted_collapsed_arr/sorted_collapsed_arr.min(axis=0),aspect =
              'auto', cmap = 'plasma',extent = [xmin,xmax,ymax,ymin]);

    ax.plot(xdata,soma_y_sorted,'k^',markersize=8);
    ax.set_xticks(xdata)
    ax.set_xticklabels(sorted_cell_ID,rotation =90, color='w')
    ytick = [ymin,ymax]
    ax.set_yticks(ytick)
    yticklab = ['pia','CC']
    ax.set_yticklabels(yticklab,color='w')
    color_bar = plt.colorbar(img, ax=ax)
    cbytick_obj = plt.getp(color_bar.ax.axes, 'yticklabels')     #getting color bar tick object          
    plt.setp(cbytick_obj, color='b') #setting colorbar ticks to white
    plt.title(input_source+': Normalized current collapsed on x axis',color='w');



#note: for bead comparison code: need to fix dataframe slicing: Use explicit indexing! 
#updated so function works by calling analysis.csv files instead of .m files. 
def bead_comparison(input_source,save = False):
    
    """
    This function will take cells of a particular input source, and compare between 
    bead negative and bead postive cells. Note: for bead positive cells, cell ID for bead negative controls 
    must be annotated manually in .m analysis file for function to work properly. 
    Function will return paired charts and ttest stats on a layer by layer basis. 
    """
    
    SCdf = load_sCRACMdatabase()
    #filter for bead postive only 
    bead_df = SCdf[(SCdf['Bead Positive?']==True)& (SCdf['MapForAnalysis']==True)]
    annotate_layers(bead_df) #annotate layer identities for each cell. 


    #get total integration for bead postive and bead 
    #note: need to change analysis_path to be set in the ini file, as a global!  
    ids = analysis_path + bead_df['CellID'] + '_analysis.csv'
    #take sum of the map
    #take sum of the bead negative control cell's map as well 
    #and add the value to the bead df 
    mapSum = []
    mapSumControl = []
    for n in ids:
        loaded_cell = pd.read_csv(n)
        mapSum.append(np.sum(np.fromstring(loadedCell['mean'].values[0].replace('\n','').replace('] [','\n').replace('[[','').replace(']]',''),sep=' ')))
        negID = analysis_path+str(loaded_cell['negativeControl'][0])+'_analysis.csv'
        mapSumControl.append(np.sum(np.fromstring(loadedControl['mean'].values[0].replace('\n','').replace('] [','\n').replace('[[','').replace(']]',''),sep=' ')))
    bead_df['bead_positive'] = mapSum
    bead_df['bead_negative']
    



    #slice by input source and normalize to negative control cell
    #reshape dataframe for seaborn plotting 
    Cdf=bead_df[(bead_df['InjectionSite']==input_source)]
    Cdf['bead_positive'] = Cdf['bead_positive']/Cdf['bead_negative']
    Cdf['bead_negative'] = Cdf['bead_negative']/Cdf['bead_negative']
    #Cdf.groupby('layer_bin').apply(lambda df: stats.ttest_rel(df['bead_negative'],df['bead_positive']))
    ContraMelt=Cdf.melt(id_vars = ['InjectionSite','layer_bin','AnimalID','CellID'], value_vars = ['bead_negative','bead_positive'])

    #plotting code: 

    unique= ContraMelt['CellID'].unique()
    colors = {}
    for i in unique:
        color = {
            i : 'tab:grey'
      }
        colors.update(color)
    def annotate(data,**kws):
        n = len(data)/2
        ax = plt.gca()
        ax.text(.1,.9,f"N = {n} pairs",transform = ax.transAxes)

    g = sns.FacetGrid(ContraMelt,col = 'layer_bin', height = 5,aspect=.6,sharey=False)
    g.map_dataframe(sns.lineplot,'variable','value',hue = 'CellID', alpha = 0.75,palette = colors,marker='o')
    g.map_dataframe(sns.lineplot,x='variable',y='value',linewidth = 2,ci = False,marker = 'o')
    g.map_dataframe(annotate)
    g.set_axis_labels(None,'absolute total normalized integration')
    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle(input_source+" bead comparison by layer ");


    #absolute total integration 
    Cdf=bead_df[(bead_df['InjectionSite']==input_source)]
    Cdf['bead_positive'],Cdf['bead_negative'] = Cdf['bead_positive'].abs(),Cdf['bead_negative'].abs()
    ContraMelt=Cdf.melt(id_vars = ['InjectionSite','layer_bin','AnimalID','CellID'], value_vars = ['bead_negative','bead_positive'])
    unique= ContraMelt['CellID'].unique()
    colors = {}
    for i in unique:
        color = {
            i : 'tab:grey'
      }
        colors.update(color)
    def annotate(data,**kws):
        n = len(data)/2
        ax = plt.gca()
        ax.text(.1,.9,f"N = {n} pairs",transform = ax.transAxes)

    g = sns.FacetGrid(ContraMelt,col = 'layer_bin', height = 5,aspect=.6,sharey=False)
    g.map_dataframe(sns.lineplot,'variable','value',hue = 'CellID', alpha = 0.75,palette = colors,marker='o')
    g.map_dataframe(sns.lineplot,x='variable',y='value',linewidth = 2,ci = False,marker = 'o')
    g.map_dataframe(annotate)
    g.set_axis_labels(None,'absolute total normalized integration')
    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle(input_source + " bead comparison by layer ");
    stat=Cdf.groupby('layer_bin').apply(lambda df: stats.ttest_rel(df['bead_negative'],df['bead_positive']))
    #stat2=Cdf.groupby('layer_bin').apply(lambda df: stats.wilcoxon(df['bead_positive']-1,zero_method='wilcox'))
    if save == True:
        g.savefig('/Users/iangingerich/Desktop/Poster_materials/bead_comparison_by_layer'+'_'+input_source+'.pdf',format = 'pdf',dpi = 300)
    return stat
