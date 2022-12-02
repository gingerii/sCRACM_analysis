
def load_sCRACMdatabase():
    """"
    Helper function, will load database will all meta data. File path is hard coded, change to 
    appropriate path, where database file lives. 
    """
    sCRACMdf = pd.read_csv('/Users/iangingerich/Desktop/sCRACM_analysis/sample_data/test_database1.csv')
    return(sCRACMdf)





# Helper function that returns all map arrays of a particular input_source, layer, but does not avearge them. Used for the page through function! 
# Updated, to pull from .csv analysis files instead of .m files 
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
        cellA = pd.read_csv(n) # read individual cell into df with nested cells 
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
        newmap_soma[n,:,:] = sCRACM_cellShifter(M[n,:,:],XY[n,:],dx,trueRowNumber_arr[n]) #use the cell shifter to align all maps by their soma coordinates (center on soma)
        X_soma[n] = X[n]-XY[n,0]
        Y_soma[n]= Y[n]+XY[n,1]
    maps_soma = newmap_soma[:,int((1+r/2)-rcAddOn/2):int((r+r/2)+rcAddOn/2+1),int((1+c/2)-rcAddOn/2):int((c+c/2)+rcAddOn/2+1)]
    
    
    return maps_soma, cells_to_average_df['CellID']



