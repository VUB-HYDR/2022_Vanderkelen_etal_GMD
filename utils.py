"""
Utils and functions to for MizuRoute postprocessing on Cheyenne

Inne Vanderkelen -  March 2021
"""

import numpy as np

def set_plot_param():
    """Set my own customized plotting parameters"""
    
    import matplotlib as mpl
    mpl.rc('xtick',labelsize=12)
    mpl.rc('ytick',labelsize=12)
    mpl.rc('axes',titlesize=16)
    mpl.rc('axes',labelsize=12)
    mpl.rc('axes',edgecolor='grey')
    mpl.rc('grid', color='lightgray')
    mpl.rc('axes',labelcolor='dimgrey')
    mpl.rc('xtick',color='dimgrey')
    mpl.rc('ytick',color='dimgrey')
    mpl.rc('text',color='dimgrey')
    mpl.rc('legend',fontsize=12, frameon=False)

        
def calc_nse(obs,mod): 
    """Nash-Sutcliffe Efficiency""" 
    
    return (1 - np.sum((mod-obs)**2)/np.sum((obs-np.mean(obs))**2))

def calc_rmse(obs,mod): 
    """Root Mean Square Error"""
    return np.sqrt(np.mean((mod-obs)**2))

def calc_kge(obs,mod): 
    """Kling-Gupta efficiency https://agrimetsoft.com/calculators/Kling-Gupta%20efficiency"""
    term1 =   (np.corrcoef(obs, mod)[0,1]-1)**2
    term2 = (np.std(mod)/np.std(obs)-1)**2
    term3 = (np.mean(mod)/np.mean(obs)-1)**2
    return 1-np.sqrt(term1 + term2 + term3)

def calc_bias(obs,mod): 
    """Mean bias""" 
    
    return np.mean(mod)-np.mean(obs)

def get_rescontrolled_pfafs(ntopo,pfaf_reservoirs,threshold): 
    """Get pfafstetter codes of main stream until which the reservoir has influence
        based on: 1. river mouth is reached (there is no downstream pfaf code)
                  2. Next reservoir on stream network is reached 
                  3. length threshold is exceeded (only inlcude segments within threshold)
                  
        input: da of river topo, list of pfaf codes of reservoirs and length threshold (in m)
        output: list of outlets corresponding to reservoir list"""
    import math
    count = 1
    print('')
    print('-------- Finding reservoir influenced streams --------')
    
    # initialise list with reservoir dependend stream segements
    controlled_pfaf = []

    for pfaf_res in pfaf_reservoirs: 
        print('processing '+str(count)+ ' of '+str(len(pfaf_reservoirs)),end='\r')

        # transform ntopo to river topo dataframe
        df_ntopo = ntopo[['PFAF','seg_id','Tosegment','Length']].to_dataframe()
        df_ntopo['PFAF'] = np.char.strip(df_ntopo['PFAF'].values.astype(str))

        # get downstream lookup table 
        downstream_lookup = get_downstream_PFAF(df_ntopo)

        # initialise
        pfaf_current = pfaf_res
        total_length = 0
        outlet_found = False


        # travel downstream from reservoir pfaf to identify outlet (end of reservoir influence)
        while not outlet_found: 

            # add current pfaf to list of controlled pfaf
            controlled_pfaf.append(pfaf_current)
            
            # next downstream segment
            pfaf_down = downstream_lookup[downstream_lookup['PFAF']==pfaf_current]['PFAF_downstream'].values[0]
            
            # res pfaf has no downstream
            if math.isnan(float(pfaf_down)): 
                pfaf_outlet = pfaf_current
                outlet_found = True

            else: 
                # add length of downstream segment to total length 
                total_length = total_length + df_ntopo[df_ntopo['PFAF']==pfaf_down]['Length'].values[0]

                # check if pfaf_current is outlet: 

                # 1. pfaf is river mouth (no outlet downstream)
                if math.isnan(float(downstream_lookup[downstream_lookup['PFAF']==pfaf_down]['PFAF_downstream'].values[0])): 

                    pfaf_outlet = pfaf_down
                    outlet_found = True
                    controlled_pfaf.append(pfaf_down)

                # 2. pfaf is other reservoir
                elif pfaf_down in pfaf_reservoirs:
                    pfaf_outlet = pfaf_down 
                    controlled_pfaf.append(pfaf_down)
                    outlet_found = True

                # 3. length threshold is exceeded for downstream segment (so include current)
                elif total_length > threshold:

                    pfaf_outlet = pfaf_current
                    outlet_found = True

                # move downstream
                else: 
                    pfaf_current = pfaf_down
                    
        count=count+1

    return controlled_pfaf

        
def get_downstream_PFAF(df_ntopo):
    """Get PFAF code of directly downstream segment
    input df of river segments with seg_id and Tosegment """

    to_segment = df_ntopo[['PFAF','seg_id']].rename(columns={'seg_id':'Tosegment', "PFAF":"PFAF_downstream"})
    df_PFAF_downstream = df_ntopo.merge(to_segment[["PFAF_downstream","Tosegment"]], on='Tosegment',how='left')
    return df_PFAF_downstream[["PFAF","PFAF_downstream"]]




### CLM FUNCTIONS


# function to cut out analysis period out of data-array (1900-2015)
def extract_anaperiod(da, stream, nspinupyears): 
    
    if nspinupyears == 0 :
        # no spin up 
        da = da[:-1,:,:]
        
    elif stream == 'h1' : # this option still to test 
        # daily timesteps
        # last day of previous year is also saved in variable therefore add one
        nspinupdays = (nspinupyears * 365) + 1
    
        # exclude spin up year and last timestep ()
        da = da[nspinupdays:-1,:,:]
        
    else: 
        # spin up with monthly timestep
        # first month of first year is not saved in variable therefore substract one
        nspinupmonths = (nspinupyears * 12) - 1
    
        # exclude spin up year and last timestep ()
        da = da[nspinupmonths:-1,:,:]

    return da