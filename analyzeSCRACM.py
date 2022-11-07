"""
Averaging and plotting sCRACM data 

Functions: 
Helpers:
    - load_sCRACMdatabase: loads databse with meta data, includes dates and cell IDs
    - annotate_layers: gives layer assignments to cells based on cortical depth
    - sCRACM_cellShifter: shifts sCRACM maps up or down (adds rows/cols) to alignm by soma
    - sCRACM_cellShifterpia: shifts sCRACM maps up or down (adds rows/cols) to alignm by pia
    - mapAverage: takes stack of maps and averages into single averaged map
    - get_cells_to_average: makes list of file locations for cells 
    - average_map_stack: returns stack of aligned maps, based on input dependent filters
Execution/plotting functions: 
    - average_map: plots average maps (soma and pia aligned) based in input dependent filters
    - page_through_maps: returns interactive plot, allowing user to scroll through all maps 
    making up a particular average map 
    - collapse_map: collapses all cells from an input source on the x-axis and plots 
    heatmaps of each cell, arranged from superfical to deep layer cells 
    - bead_comparison: plots paired charts showing total integration of bead positive cells 
    against negative controls per input source, per layer. Returns ttest stats for each layer
"""

from IPython.core.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))
import numpy as np
from scipy import optimize, stats
import math
import seaborn as sns
import matplotlib.pylab as plt
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
import os
import sys
from tqdm.notebook import tqdm
from IPython.core.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))
np.set_printoptions(threshold=sys.maxsize)
sys.path.insert(0, '/Users/iangingerich/Desktop/morph_analysis/') #to direct python to location of interpreter! 
import mstruct2pydict as m2p 
from plotly import graph_objects as go
analysis_path = '/Users/iangingerich/Dropbox (OHSU)/Ian_Tianyi_shared/Databases/Analyzed_cells/'

# 

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
    soma depth. Database must contain soma depth and cortical thickness information to be annotated. 
    """
    df['ratio_thickness'] = df['piadistance']/df['cortical thickness']
    labels = ['L2','L3','L5','L6','Claustrum']
    bins = [0,0.231,0.4080,0.65275,0.848617,1]
    df['layer_bin']= pd.cut(df['ratio_thickness'],bins=bins,labels = labels)
    
# Helper functions listed first 
def  sCRACM_cellShifter(map_n,xy,spacing,trueRowNumber):
    """
    Helper function used when averaging multiple sCRACM maps together. used when averaging maps of different 
    sizes. Will shift maps up or down to align individual maps to each other. Used when aligning individual maps
    to the soma of each recording. 
    """
    
    [r,c] = np.shape(map_n)
    newmap = np.zeros([2*r,2*c])
    #go through maps one at a time
    #calculate how many rows and columns to shift by. x for rows, y for columns 

    shiftRows = np.round(xy[1]/spacing)+(r-trueRowNumber)/2
    shiftCols = -np.round(xy[0]/spacing)
    #place the original map onto the new map, with shift 
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

    shiftRows = np.round(xy[1]/spacing)-(trueRowNumber/2)
    shiftCols = -np.round(xy[0]/spacing)
    #place the original map onto the new map, with shift 
    newmap[int(1 + np.round(r/2)+shiftRows):int(r+np.round(r/2)+shiftRows+1),int(1+np.round(c/2)+shiftCols):int(c+np.round(c/2)+shiftCols+1)] = map_n
    return newmap


# mapAverage helper function: takes stack of maps to be averaged, after shifting, normalization, and col/row adding, and averages together into a single 24 by 12 array
def mapAverage(mapStack):
    """
    Helper function called when plotting average sCRACM maps. Will take (x,y,n) array of maps, average by last dim,
    returning single (x,y) array, which represents the average map. 
    """
    #mapStack = maps
    [planes,rows,cols] = np.shape(mapStack)
    mapAvg = np.zeros([rows,cols])
    for r in range(0,rows): # maybe at +1
        for c in range(0,cols):
            mapAvg[r,c] = np.mean(mapStack[:,r,c])
    return mapAvg 




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
    #annotate cells based on relative thickness. Note: layer 1 cells not excluded in this code!  
    labels = ['L2','L3','L5','L6','Claustrum']
    bins = [0,0.231,0.4080,0.65275,0.848617,1]
    sCRACMdf['layer_bin']= pd.cut(sCRACMdf['ratio_thickness'],bins=bins,labels = labels)
    
    # average based on relative cortical depth assignment
    
   
    
    if by_layer_bin == True:
        printout = 'by_layer_bin'
        cells_to_average_df = sCRACMdf[(sCRACMdf['InjectionSite']==input_source) & (sCRACMdf['MapForAnalysis']==True) & (sCRACMdf['Response']==True)
                                   &(sCRACMdf['layer_bin']== layer)]# & sCRACMdf['Bead Positive?']== bead_positive] 
        cell_IDs = analysis_path + cells_to_average_df['CellID'] + '_analysis.m'

        #Nested IF statments: by layer bin first, and then add additional filters
        # average based on cellular M_type 
        if by_Mtype != False:
            printout = 'by_Mtype: '+str(by_Mtype)
            cells_to_average_df = cells_to_average_df[(cells_to_average_df['M_type_new']== by_Mtype)] #filter based on mapfor analysis and input source
            cell_IDs = analysis_path + cells_to_average_df['CellID'] + '_analysis.m'


      

        elif bead_positive == True:
            printout = 'bead_positive'
            cells_to_average_df = cells_to_average_df[(cells_to_average_df['Injection Type']=='virus+bead')& (cells_to_average_df['Bead Positive?']== True)]

            cell_IDs = analysis_path + cells_to_average_df['CellID'] + '_analysis.m'
            
        elif bead_positive == 'negative_control':
            printout = 'negative_control'
            cells_to_average_df = cells_to_average_df[(cells_to_average_df['Injection Type']=='virus+bead')& (cells_to_average_df['Bead Positive?']== False)]
            cell_IDs = analysis_path + cells_to_average_df['CellID'] + '_analysis.m'
            
    

    # average based on online designation of layer 
    elif by_layer_online == True:
        printout = 'by_layer_online'
        cells_to_average_df = sCRACMdf[(sCRACMdf['InjectionSite']==input_source) & (sCRACMdf['MapForAnalysis']==True) & (sCRACMdf['Response']==True)
                                   &(sCRACMdf['Layer estimation']== layer)]#& sCRACMdf['Bead Positive?']==bead_positive] #filter based on mapfor analysis and input source
        cell_IDs = analysis_path + cells_to_average_df['CellID'] + '_analysis.m'
    
        #nested IF statement: by online layer designation first, and then other filters!
        # average based on cellular M_type 
        if by_Mtype != False:
            printout = 'by_Mtype: '+str(by_Mtype)
            cells_to_average_df = cells_to_average_df[(cells_to_average_df['M_type_new']== by_Mtype)] #filter based on mapfor analysis and input source
            cell_IDs = analysis_path + cells_to_average_df['CellID'] + '_analysis.m'


        
        elif bead_positive == True:
            printout = 'bead_positive'
            cells_to_average_df = cells_to_average_df[(cells_to_average_df['Injection Type']=='virus+bead')& (cells_to_average_df['Bead Positive?']== True)]

            cell_IDs = analysis_path + cells_to_average_df['CellID'] + '_analysis.m'
            
        elif bead_positive == 'negative_control':
            printout = 'negative_control'
            cells_to_average_df = cells_to_average_df[(cells_to_average_df['Injection Type']=='virus+bead')& (cells_to_average_df['Bead Positive?']== False)]
            cell_IDs = analysis_path + cells_to_average_df['CellID'] + '_analysis.m'
    
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
    
    cells_to_average_df,cell_IDs,printout=get_cells_to_average(input_source=input_source, 
                         layer = layer,
                         by_layer_bin = by_layer_bin,
                         by_layer_online = by_layer_online, 
                         by_Mtype = by_Mtype,
                         bead_positive = bead_positive)
    
    numCells = np.shape(cell_IDs)[0]
   
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
        cellA = m2p.mstruct2pydict(n) #use Mike's parcer to take .m files that have already been processed 
        cellA_mean = cellA['mean'] #define the map 
        XY[i] = cellA['soma1Coordinates'] #save coordinates 
        [a,b] = cellA_mean.shape
        trueRowNumber = cellA_mean.shape[0]
        layer1Row[0,i] = cellA['layer1Row']
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
        average_map_new[0:int((cellA['layer1Row'])-1),:] = 0 
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
        newmap_soma[n,:,:] = sCRACM_cellShifter(M[n,:,:],XY[n,:],dx,trueRowNumber_arr[n])
        X_soma[n] = X[n]-XY[n,0]
        Y_soma[n]= Y[n]+XY[n,1]
    maps_soma = newmap_soma[:,int((1+r/2)-rcAddOn/2):int((r+r/2)+rcAddOn/2+1),int((1+c/2)-rcAddOn/2):int((c+c/2)+rcAddOn/2+1)]
    
    
    return maps_soma, cells_to_average_df['CellID']

########## average sCRACM maps here, will call above helper function to average based on a specfic parameter! ###############

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
        cellA = m2p.mstruct2pydict(n) #use Mike's parcer to take .m files that have already been processed 
        cellA_mean = cellA['mean'] #define the map 
        XY[i] = cellA['soma1Coordinates'] #save coordinates 
        [a,b] = cellA_mean.shape
        trueRowNumber = cellA_mean.shape[0]
        layer1Row[0,i] = cellA['layer1Row']
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
        average_map_new[0:int((cellA['layer1Row'])-1),:] = 0 
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
    IDs = analysis_path + contraidf['CellID'] + '_analysis.m'
    cell_ID = contraidf['CellID']
    # convert cell IDs from DF to np.array to use for labels in the figure
    cell_ID_array = cell_ID.to_numpy()[...,None] 
    # initialize variable to store nested dict. struct. Could maybe change to dataframe in the future 
    max_row_n = 24 #used to add extra rows to smaller maps. Update to size of map instead of 24. 
    collapsed_arr = [[] for n in range(0,24)]
    soma_y_arr = []
    trueRowNumber_arr = []
    # change back to IDs to get full list of paths to files! 
    for n in tqdm(IDs): 
        loaded_cell = m2p.mstruct2pydict(n) #note: n is the path name, so passing it n will result in the .m file being read. 
        average_map = loaded_cell['mean'] #take just the average and soma coordinates!
        [a,b] = average_map.shape
        trueRowNumber = average_map.shape[0]
        trueRowNumber_arr = [trueRowNumber_arr,trueRowNumber]
        # if cell is less than 24 rows, add zeros, and adjust soma position, else, keep as is 
        if average_map.shape[0] < max_row_n:
            add_arr = np.zeros(((24-average_map.shape[0]),b)) #number of rows of zeros to add, note, change ending from 12 to 1, as we are adding after collapsing
            average_map_new = np.concatenate((average_map,add_arr),axis = 0) # add rows of zeros here
            shift_unit = max_row_n - trueRowNumber
            shift_val = (shift_unit*0.5)*50 # value to subtract from Y coordinate 
            soma_Y = (-1*(loaded_cell['soma1Coordinates'][1]))-shift_val #transformed Y-coordinate to plot on top of the heatmap
        else: 
            soma_Y = (-1*(loaded_cell['soma1Coordinates'][1])) 
            average_map_new = average_map
        layer1Row = loaded_cell['layer1Row']
        average_map_new[0:int((loaded_cell['layer1Row'])-1),:] = 0 #remove data over Pia by setting sites to 0, note, instead of 3 it should be the layer 1 number in the original dictinary, minus 1! 
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
    #calculate ratio thickness and put into a column 
    bead_df['ratio_thickness'] = bead_df['piadistance']/bead_df['cortical thickness']
    # annotate to give each cell a layer identity based on bins 
    labels = ['L2','L3','L5','L6','Claustrum']
    bins = [0,0.231,0.4080,0.65275,0.848617,1]
    bead_df['layer_bin']= pd.cut(bead_df['ratio_thickness'],bins=bins,labels = labels)


    #get total integration for bead postive and bead negative 
    ids = analysis_path + bead_df['CellID'] + '_analysis.m'
    bead_df['bead_positive'] = [np.sum(m2p.mstruct2pydict(n)['mean']) for n in ids]
    bead_df['bead_negative']= [np.sum(m2p.mstruct2pydict((analysis_path + m2p.mstruct2pydict(n)['negativeControl'] + '_analysis.m'))['mean']) for n in ids]

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
