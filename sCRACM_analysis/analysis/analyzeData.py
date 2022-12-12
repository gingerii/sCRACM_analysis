"""
Averaging and plotting of cleaned sCRACM data.

Functions: 
Helpers:
    - forwardAffineTransform: preform affine transform on column-shaped vectors 
    - checkMetaData: checks that manual annotations are present in .csv analysis files once loaded.
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
    - plotSWC: Will plot annotations and traced cell from an SWC file. 
    - overlaySWC: Overlays sCRACM map with SWC tracing and .tif actute slice image 

Citations and such: 
Much of the sCRACM map averaing has been adapted, with permission, from Petreanu et al., 2009: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2745650/  (original code writen in Matlab)

"""

#imports 
import numpy as np
import pandas as pd
import scipy as sc
from scipy import optimize, stats,io
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
import sys
from plotly import graph_objects as go
from tkinter import filedialog
from tkinter import Tk
import csv
from sCRACM_analysis.sCRACM_global import path_database,filename_cell_db,path_ephys_data
from sCRACM_analysis.analysis.cleanData import cart2pol,pol2cart



def forwardAffineTransform(T,v1,v2):
    """
    Preform an affine tranformation of column shaped vectors. 
    """
    if v1.shape[1] != 1 or v2.shape[1] != 1:
        print('Vectors must be column-shaped!')
        return
    elif v1.shape[0] != v2.shape[0]:
        print('Vectors must be of equal length!')
        return

    vecSize = v1.shape[0]

    concVec = np.concatenate((v1,v2),axis=1)
    onesVec = np.ones((vecSize,1))

    U = np.concatenate((concVec,onesVec),axis=1)

    retMat = np.dot(U,T[:,0:2])

    return (retMat[:,0].reshape((vecSize,1)), retMat[:,1].reshape((vecSize,1)))


def checkMetaData(cell):
    """
    Function checks meta data of .csv files for missing values. 
    """
    ID = cell['experimentNumber'][0]
    if cell['layer1Row'].isnull().values.any() == True:
        print(f'Cell {ID} is missing layer1Row annotation, fix in .csv')
    if cell['distanceToPia'].isnull().values.any() == True:
        print(f'Cell {ID} is missing distance to pia annotation, fix in .csv')
    if cell['cortexThickness'].isnull().values.any() == True:
        print(f'Cell {ID} is missing cortical thickness annotation, fix in .csv')


def load_sCRACMdatabase():
    """"
    Helper function, will load database will all meta data. File path is hard coded, change to 
    appropriate path, where database file lives. 
    """

    sCRACMdf = pd.read_csv(path_database+filename_cell_db)
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
        cell_IDs = path_ephys_data + cells_to_average_df['CellID'] + '_analysis.csv' #append together file names to load 

        #Nested IF statments: by layer bin first, and then add additional filters
        # average based on cellular M_type 
        if by_Mtype != False: #if by_Mtype is something other than false (default), then filter df my m_type as well. 
            printout = 'by_Mtype: '+str(by_Mtype)
            cells_to_average_df = cells_to_average_df[(cells_to_average_df['M_type_final']== by_Mtype)] #filter based on newest M_type assignments 
            cell_IDs = path_ephys_data + cells_to_average_df['CellID'] + '_analysis.csv'


      # if you want to get the file names of bead positive cells, bead positve= True during input

        elif bead_positive == True:
            printout = 'bead_positive'
            cells_to_average_df = cells_to_average_df[(cells_to_average_df['Injection Type']=='virus+bead')& (cells_to_average_df['Bead Positive?']== True)]

            cell_IDs = path_ephys_data + cells_to_average_df['CellID'] + '_analysis.csv'
            
        elif bead_positive == 'negative_control': #note: right now this only works for cells that have a negative control designated in their .m file. For full python code, the inport settings 
            # also needs to have a negative control variable to filtering 
            printout = 'negative_control'
            cells_to_average_df = cells_to_average_df[(cells_to_average_df['Injection Type']=='virus+bead')& (cells_to_average_df['Bead Positive?']== False)]
            cell_IDs = path_ephys_data + cells_to_average_df['CellID'] + '_analysis.csv'
            
    

    # average based on online designation of layer, not the average layer consensus. This should be the option if using this code to average cells recorded in places outside the AI 
    elif by_layer_online == True:
        printout = 'by_layer_online'
        cells_to_average_df = sCRACMdf[(sCRACMdf['InjectionSite']==input_source) & (sCRACMdf['MapForAnalysis']==True) & (sCRACMdf['Response']==True)
                                   &(sCRACMdf['Layer estimation']== layer)]#& sCRACMdf['Bead Positive?']==bead_positive] #filter based on mapfor analysis and input source
        cell_IDs = path_ephys_data + cells_to_average_df['CellID'] + '_analysis.csv'
    
        #nested IF statement: by online layer designation first, and then other filters!
        # average based on cellular M_type 
        if by_Mtype != False:
            printout = 'by_Mtype: '+str(by_Mtype)
            cells_to_average_df = cells_to_average_df[(cells_to_average_df['M_type_final']== by_Mtype)] #filter based on mapfor analysis and input source
            cell_IDs = path_ephys_data + cells_to_average_df['CellID'] + '_analysis.csv'


        
        elif bead_positive == True:
            printout = 'bead_positive'
            cells_to_average_df = cells_to_average_df[(cells_to_average_df['Injection Type']=='virus+bead')& (cells_to_average_df['Bead Positive?']== True)]

            cell_IDs = path_ephys_data + cells_to_average_df['CellID'] + '_analysis.csv'
            
        elif bead_positive == 'negative_control':
            printout = 'negative_control'
            cells_to_average_df = cells_to_average_df[(cells_to_average_df['Injection Type']=='virus+bead')& (cells_to_average_df['Bead Positive?']== False)]
            cell_IDs = path_ephys_data + cells_to_average_df['CellID'] + '_analysis.csv'
    
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
        cellA = pd.read_csv(n)
        checkMetaData(cellA) #check to make sure manual annotations are present!
        #formate map data 
        cellA_mean = np.fromstring(cellA['mean'].values[0].replace('\n','').replace('] [','\n').replace('[[','').replace(']]',''),sep=' ') # parse map array into correct datatype
        cols = 12 #hard code number of columns for now 
        rows = int(cellA_mean.shape[0]/cols) #number of rows for the map 
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
        checkMetaData(cellA) #check meta data to make sure manual annotations are present!
        #format map data correctly 
        cellA_mean = np.fromstring(cellA['mean'].values[0].replace('\n','').replace('] [','\n').replace('[[','').replace(']]',''),sep=' ') # parse map array into correct datatype
        cols = 12 #hard code number of columns for now 
        rows = int(cellA_mean.shape[0]/cols) #number of rows for the map 
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
        mapSpacing = cellA['xSpacing'][0] #assumes dx and dy are the same, which they are, both 50.  
        #next step normalizes maps to min value 
        minVal = np.min(average_map_new) #find min for each map
        M[i,:,:] = average_map_new/minVal 
        distanceToPia.append(cellA['distanceToPia'])

        # end for loop now 
    rcAddOn = 0 #must be even interger 
    [r,c]=np.shape(M[1])
    X = XY[:,0]
    Y = -XY[:,1]
    dx = mapSpacing #should always be 50um on LSPS rig 
    # shift maps for alignmnent 
    # preallocate newmap variable 
    ###############################################################
    newmap = np.zeros([numCells,2*r,2*c]) 
    #initialize figure 
    fig  = plt.figure(figsize=(30,30))
    fig.suptitle('Input: '+input_source+', Averaged '+ printout + ' (n='+str(numCells)+')',color = 'w',fontsize = 30)
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
        ax1.plot(X_soma[n],Y_soma[n],"^",markersize=8, color='b');
    # plot for input intensity by row 
    rowNum,_ = np.shape(rowSum)
    yarr = np.linspace(0.5, (rowNum-0.5), num = rowNum)
    yarr = yarr[::-1] #flip 
    
    ax3.errorbar(rowSum.reshape((-1)), yarr, xerr= rowSumSEByCell.reshape((-1)), fmt = '-o');
    ax3.set_title('Relative intensity', color = 'w');
   
    #plot for input intenstiy by column 
    _,colNum = np.shape(colSum)
    yarr = np.linspace(0.5,(colNum-0.5),num = colNum)
    ax2.errorbar(yarr,colSum.T.reshape((-1)),yerr=colSumSEByCell.T.reshape((-1)),fmt = '-o');
    ax2.set_title('Soma alignment', color = 'w', fontsize = 20)
    
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
        ax4.plot(X_pia[n],Y_pia[n],"^",markersize=8, color = 'b');
    
    # plot for input intensity by row 
    rowNum,_ = np.shape(rowSum)
    yarr = np.linspace(0.5, (rowNum-0.5), num = rowNum)
    yarr = yarr[::-1] #flip 
    ax6.errorbar(rowSum.reshape((-1)), yarr, xerr= rowSumSEByCell.reshape((-1)), fmt = '-o');
    ax6.set_title('Relative intensity', color = 'w');

    #plot for input intenstiy by column 
    _,colNum = np.shape(colSum)
    yarr = np.linspace(0.5,(colNum-0.5),num = colNum)
    ax5.errorbar(yarr,colSum.T.reshape((-1)),yerr=colSumSEByCell.T.reshape((-1)),fmt = '-o');
    ax5.set_title('Pia alignment', color = 'w', fontsize = 20);
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
    test_frames = [go.Frame(data = go.Heatmap(z = M2[i,:,:],x= xdata,y=ydata,colorscale='magma'), name = cell) for i,cell in zip(range(e),cell_id.values)]


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
    IDs = path_ephys_data + contraidf['CellID'] + '_analysis.csv' #should be reading in the analysis.csv files 
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
        # check that necessary data has been manually annotated! 
        checkMetaData(loaded_cell) #check meta data to make sure manual annotations are present! 
        # parse map array into correct datatype
        average_map = np.fromstring(loaded_cell['mean'].values[0].replace('\n','').replace('] [','\n').replace('[[','').replace(']]',''),sep=' ') 
        #DANGER: cols variable is hard coded, if averaging larger maps, this need to be flexible! 
        cols = 12 #hard code number of columns for now 
        rows = int(average_map.shape[0]/cols) #number of rows for the map 
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
              'auto', cmap = 'magma_r',extent = [xmin,xmax,ymax,ymin]);
    
    ax.plot(xdata,soma_y_sorted,'^',markersize=8,color = 'b');
    ax.set_xticks(xdata)
    ax.set_xticklabels(sorted_cell_ID,rotation =90, color='w')
    ytick = [ymin,ymax]
    ax.set_yticks(ytick)
    yticklab = ['pia','CC']
    ax.set_yticklabels(yticklab,color='w')
    color_bar = plt.colorbar(img, ax=ax)
    cbytick_obj = plt.getp(color_bar.ax.axes, 'yticklabels')     #getting color bar tick object          
    plt.setp(cbytick_obj, color='w') #setting colorbar ticks to white
    plt.title(input_source+': Total current collapsed on x axis', color = 'w')
    
    # below for normalized plots 
    _,ax = plt.subplots(figsize = (24,6)) #create plot, use underscore to exclude getting first output from subplots
    img = ax.imshow(sorted_collapsed_arr/sorted_collapsed_arr.min(axis=0),aspect =
              'auto', cmap = 'magma',extent = [xmin,xmax,ymax,ymin]);

    ax.plot(xdata,soma_y_sorted,'^',markersize=8,color = 'w');
    ax.set_xticks(xdata)
    ax.set_xticklabels(sorted_cell_ID,rotation =90, color='w')
    ytick = [ymin,ymax]
    ax.set_yticks(ytick)
    yticklab = ['pia','CC']
    ax.set_yticklabels(yticklab,color='w')
    color_bar = plt.colorbar(img, ax=ax)
    cbytick_obj = plt.getp(color_bar.ax.axes, 'yticklabels')     #getting color bar tick object          
    plt.setp(cbytick_obj, color='w') #setting colorbar ticks to white
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
    annotate_layers(SCdf) #annotate layer identities for each cell. 
    #filter for bead postive only 
    bead_df = SCdf.loc[(SCdf['Bead Positive?']==True)& (SCdf['MapForAnalysis']==True)].copy()
    


    #get total integration for bead postive and bead 
    #note: need to change analysis_path to be set in the ini file, as a global!  
    ids = path_ephys_data + bead_df['CellID'] + '_analysis.csv'
    #take sum of the map
    #take sum of the bead negative control cell's map as well 
    #and add the value to the bead df 
    mapSum = []
    mapSumControl = []
    for n in ids:
        loaded_cell = pd.read_csv(n)
        checkMetaData(loaded_cell) #check metadata to make sure manual annotations are present! 
        mapSum.append(np.sum(np.fromstring(loaded_cell['mean'].values[0].replace('\n','').replace('] [','\n').replace('[[','').replace(']]',''),sep=' ')))
        negID = path_ephys_data+str(loaded_cell['negativeControl'][0])+'_analysis.csv'
        loadedControl = pd.read_csv(negID)
        mapSumControl.append(np.sum(np.fromstring(loadedControl['mean'].values[0].replace('\n','').replace('] [','\n').replace('[[','').replace(']]',''),sep=' ')))
    bead_df['bead_positive'] = mapSum
    bead_df['bead_negative'] = mapSumControl
    



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
    g.map_dataframe(sns.lineplot,x='variable',y='value',linewidth = 2,errorbar = ('ci',False),marker = 'o')
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
    g.map_dataframe(sns.lineplot,x='variable',y='value',linewidth = 2,errorbar = ('ci',False),marker = 'o')
    g.map_dataframe(annotate)
    g.set_axis_labels(None,'absolute total normalized integration')
    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle(input_source + " bead comparison by layer ");
    stat=Cdf.groupby('layer_bin').apply(lambda df: stats.ttest_rel(df['bead_negative'],df['bead_positive']))
    #stat2=Cdf.groupby('layer_bin').apply(lambda df: stats.wilcoxon(df['bead_positive']-1,zero_method='wilcox'))
    if save == True:
        g.savefig('/Users/iangingerich/Desktop/Poster_materials/bead_comparison_by_layer'+'_'+input_source+'.pdf',format = 'pdf',dpi = 300)
    return stat


def plotSWC(axis=False,
        flip=False,
        spatial_rotation = 0,
        input_rotation = 0, 
        shrinkX = 1,
        shrinkY = 1,
        tweekSomaX = 0,
        tweekSomaY = 0,
        soma1Coordinates = [0,0],
        file_path = False):  
    """
    Plots swc file, used to check if tracing needs to be flipped 
    over x-axis before overlaying with .tif image and sCRACM heatmap.
    """
    
    if file_path == False:
        print('Please select traced neuron SWC file')
        Tk().withdraw()
        file_path = filedialog.askopenfilename()
    #read file 
    file_data = np.loadtxt(file_path)
    #make all dendrites have same annotation! 
    logic = file_data[:,1] ==4
    file_data[logic,1] = 3
    array = file_data
    if flip ==True:
        #this line will flip over the x-axis: 
        array[:,2] = array[:,2]*-1
    if axis == False: 
        fig = plt.figure()
        axis = fig.add_subplot(111)
    axis.set_facecolor('black')

    #total rotation should be your tweek factor, plus the spatialRotation stored in the cell's
    #.csv analysis file. Default is no rotation. 
    rotation = input_rotation + (spatial_rotation/180*np.pi)
    #using above functions: 
    [rho,theta]=cart2pol(array[:,2], array[:,3])
    #apply rotation, and convert back 
    [newX,newY]=pol2cart(rho,theta+rotation)
    #assign rotation adjusted coordinates back to swc array
    array[:,2] = newX
    array[:,3]= newY

    # find center of mass 
    contourArray = []
    [r,c] = np.shape(array)
    A = np.arange(r) #range to iterate over (rows)
    # find where there is a soma annotation, pull out the x,y,z coordinates, and average
    arraySoma=array[array[:,1] == 1,:] #subset data using indexing instead of for loop (faster)
    somaX = arraySoma[:,2]
    somaY=arraySoma[:,3]
    #find center 
    centerX = np.mean(somaX)
    centerY = np.mean(somaY)
    # apply shrinkage 
    array[:,2] = array[:,2]/shrinkX
    array[:,3] = array[:,3]/shrinkY
    # center soma position respect to the map
    array[:,2] = array[:,2] - centerX + soma1Coordinates[0] + tweekSomaX 
    array[:,3] = array[:,3] - centerY + soma1Coordinates[1] + tweekSomaY 
   #plotting  
    contourArray = []
    for n in A:
        if array[n,1] == 1: #ie, its a soma
            somaX[n] = array[n,2]
            somaY[n] = array[n,3]
        elif array[n,1] ==2 or array[n,1] == 3:
            parent = int(array[n,6])
            if parent != -1: 
                x1 = array[parent-1,2]
                y1 = array[parent-1,3]
                x2 = array[n,2]
                y2 = array[n,3]
                hline, = axis.plot([x1,x2],[y1,y2],color = 'w')
                if array[n,1] == 2:
                    hline.set_color('w') # for coloring axon red
                elif array[n,1] == 3: 
                    hline.set_color('w')
        elif array[n,1] > 3: # its some contour
            contourArray.append([array[n,1],array[n,2],array[n,3],array[n,4]])
    contourArray = np.array(contourArray) # convert from list to np.array 
    [r,c] = np.shape(contourArray)  
    if r ==0 or c == 0: 
        print('no contours annotated')# must not be any contours, skip to next part 
    else: 
        B = np.unique(contourArray[:,0])
        for n in B:
            inds = np.where(contourArray[:,0]==n)
            #inds = contourArray[:,0].index(n)
            startRow = np.min(inds)
            endRow = np.max(inds)
            Xdata = contourArray[startRow:endRow,1]
            Ydata = contourArray[startRow:endRow,2]
            #Zdata = contourArray[startRow:endRow,3]
            hline, = axis.plot(Xdata,Ydata,color = 'w')
    return file_path

def overlaySWC(CellID,
               input_rotation = 0,
               shrinkX = 1,
               shrinkY = 1,
               tweekSomaX = 0,
               tweekSomaY = 0):
    
    print('check orientation of tracing')
    #plot out the SWC and see if it needs to be flipped 
    file_path=plotSWC()
    plt.show()
    #capture flip after checking plotting 
    user_input=input('Does the SWC need to be flipped (yes/no)')
    if user_input.lower()=='yes':
        print('SWC will be flipped over x-axis')
        flip = True
    elif user_input.lower() == 'no':
        print('SWC is correct as is')
        flip = False
        
    # Picture 
    print('Please select the acute slice .tif image')
    Tk().withdraw()
    picture_path = filedialog.askopenfilename()
    im = plt.imread(picture_path)
    R = im/(np.max(im[:])/110)
    
    #load analysis .csv file for single cell 
    cell_path = os.path.join(path_ephys_data,CellID+'_analysis.csv')
    loaded_cell = pd.read_csv(cell_path)
    #define some variables 
    XRange=loaded_cell.horizontalVideoDistance[0]
    YRange = loaded_cell.verticalVideoDistance[0]
    mapSpacing = loaded_cell.xSpacing[0]
    spatialRotation = loaded_cell.spatialRotation[0]
    xPatternOffset = loaded_cell.xPatternOffset[0]
    yPatternOffset = loaded_cell.yPatternOffset[0]
    soma1Coordinates = np.fromstring(loaded_cell.soma1Coordinates[0].replace('[','').replace(']',''),sep = ',')
    #sCRACM map 
    M=np.fromstring(loaded_cell['mean'].values[0].replace('\n','').replace('] [','\n').replace('[[','').replace(']]',''),sep=' ') 
    cols = 12 #hard code number of columns for now 
    rows = int(M.shape[0]/cols) #number of rows for the map 
    M = np.reshape(M,(rows,cols)) #reshape map to be (row,cols)
    
    #initialize figure 
    [Ypix,Xpix] = np.shape(R)
    xVideoScaleFactor = XRange/Xpix
    yVideoScaleFactor = YRange/Ypix

    %matplotlib inline
    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(111)

    extent = (-.5*XRange,.5*XRange, -.5*YRange,.5*YRange)
    im = plt.imshow(R, extent = extent,cmap = 'gray') #extent key argument most similar to Xdata and Ydata
    im.axes.tick_params(axis = 'both',which='both',bottom=False,top=False,labelbottom=False,labelleft=False)
    ax.set_title(loaded_cell['experimentNumber'][0], color = 'w');
    ax.set_aspect(1)
    ax.set_xlim([-XRange/2,XRange/2])
    ax.set_ylim([-YRange/2,YRange/2])

    ############################# MAP ##########################################
    #get averaged map
    #M = test_cell['mean']
    [yn,xn] = np.shape(M)
    # initalize arrays for patch values 
    xgrid = np.empty((4,(xn*yn)))
    xgrid[:] = np.nan
    ygrid = np.empty((4,(xn*yn)))
    ygrid[:] = np.nan
    cgrid = np.empty((1,xn*yn))
    cgrid[:] = np.nan
    # iterate through all x/y combos of the stimulation grid 

    tmp = []
    for x in range(xn):
        for y in range(yn):
            idx = (x) * yn + y #index for storing in x/y grids (note differnce in indexing)
            tmp = np.append(tmp,idx)
            xgrid[:,idx] = np.squeeze(np.array([(np.array([-1,1,1,-1]) * mapSpacing/2) + ((x+1)*mapSpacing)]).T) #Coords for 4 corners of the patch
            ygrid[:,idx] = np.squeeze(np.array([(np.array([-1,-1,1,1]) * mapSpacing/2) + ((y+1)*mapSpacing)]).T) #Cords for 4 corners of the patch
            cgrid[:,idx] = M[y,x] #test to make sure this indexing is correct! 


    #translate grid according to metadata. NOTE inversion of y-axis! 
    xgrid1 = xgrid - (xn*mapSpacing/2) - mapSpacing/2
    ygrid1 = -(ygrid -(yn*mapSpacing/2)-mapSpacing/2)
    # make sure all values at or above zero are transparent (replace zeros with Nan)
    cgrid[cgrid == 0] = np.nan
    # rotate grid according to metadata (note, this data in only contained in the workspace variable!)
    tform = np.array([[np.cos(np.deg2rad(spatialRotation)),np.sin(np.deg2rad(spatialRotation)),0],[-np.sin(np.deg2rad(spatialRotation)),np.cos(np.deg2rad(spatialRotation)),0],[xPatternOffset,yPatternOffset,1]])
    ygrid_vector=ygrid1[:].T.reshape(-1,1)
    xgrid_vector = xgrid1[:].T.reshape(-1,1)
    #transform coordinates 
    [xt,yt]=forwardAffineTransform(tform,xgrid_vector,ygrid_vector)
    #reshape 
    xgrid2 = xt.reshape(4,((xn*yn)),order='F')
    ygrid2 = yt.reshape(4,((xn*yn)),order='F')
    #color valus: 
    Z = np.vstack((cgrid,cgrid,cgrid,cgrid))


    # plot the grid using pcolormesh.
    #note: pcolormesh update #16258 depreciates shading = 'flat' when data and coords have the same shape 
    # and introduces shading ='nearest' which interpolates the coordinates. However, the logic for 
    # the interpolation breaks down if they are not monotonically incerasing, which is the case here. 
    # to fix, drop the last column of Z and use shading='flat'
    cm =ax.pcolormesh(xgrid2,ygrid2,Z[:-1,:-1],cmap = 'magma_r', shading = 'flat')
    cbaxes = inset_axes(ax,width = '3%',height = '50%', loc = 1)
    color_bar = plt.colorbar(cm,cax=cbaxes)
    color_bar.ax.yaxis.set_tick_params(color = 'white')
    plt.setp(plt.getp(color_bar.ax.axes,'yticklabels'),color='white')
    #fig.colorbar(cm,cax = cbaxes)
    ax.set_aspect('equal')
    #hgrid = 



    ######################## Neurolucida Tracing ########################
    plotSWC(axis = ax,
            flip = flip,
            spatial_rotation = spatialRotation,
            input_rotation = input_rotation,
            shrinkX = 1,
            shrinkY=1,
            tweekSomaX=0,
            tweekSomaY=0,
            soma1Coordinates=soma1Coordinates,
            file_path=file_path)

   

