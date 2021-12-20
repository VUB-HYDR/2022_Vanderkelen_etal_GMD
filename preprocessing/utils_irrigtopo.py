"""
Utils and functions to determine irrigation topology 

Inne Vanderkelen -  March 2021
"""

# import modules
import pfaf.pfafstetter as pfaf # decode package from Naoki, see https://github.com/nmizukami/pfaf_decode 
import pandas as pd
import numpy as np
import geopandas as gpd


###################
# 1. Helper functions

def get_mainstream(pfafs, pfaf_outlet):
    """get list of mainstream pfaf codes"""
    mainstream = []
    for p in pfafs: 
        #print(pfaf.get_tributary(p, pfaf_outlet))
        if pfaf.get_tributary(p, pfaf_outlet)=='-999': 
                mainstream.append(p)
    return mainstream

def get_downstream_PFAF(river_shp):
    """Get PFAF code of directly downstream segment
    input df of river segments with seg_id and Tosegment """
    
    to_segment = river_shp[['PFAF','seg_id']].rename(columns={'seg_id':'Tosegment', "PFAF":"PFAF_downstream"})
    df_PFAF_downstream = river_shp.merge(to_segment[["PFAF_downstream","Tosegment"]], on='Tosegment',how='left')
    return df_PFAF_downstream[["PFAF","PFAF_downstream"]]


def get_pfafs_start2end(pfaf_start,pfaf_end,river_shp, include_end=True):
    """Get pfafs of stream network of pfaf_start to pfaf_end"""

    # get PFAF codes and downstream PFAF codes
    df_PFAF_downsteam = get_downstream_PFAF(river_shp)

    # get list of pfafs between start and end
    pfaf_start2end = [pfaf_start]

    # initialise
    pfaf_current = pfaf_start 
    pfaf_down = pfaf_start
    # go downstream
    while pfaf_down != pfaf_end:
        if (df_PFAF_downsteam['PFAF']==pfaf_current).sum() > 0: 
            pfaf_down =  df_PFAF_downsteam.loc[df_PFAF_downsteam['PFAF']==pfaf_current,'PFAF_downstream'].values[0]
            pfaf_start2end.append(pfaf_down)
            pfaf_current = pfaf_down
            
    if not include_end: pfaf_start2end.remove(pfaf_start)
        
    return pfaf_start2end

def get_streamlength_total(pfaf_start,pfaf_end,river_shp): 
    """Get river stream length between start pfaf and end pfaf, including length of pfaf end and start
    input: start and end pfaf code, df with river network (pfaf, to_segment, length)"""
    
    # get PFAF codes between start and end
    pfafs_start2end = get_pfafs_start2end(pfaf_start,pfaf_end,river_shp)
    
    # calculate total legth between 2 segments (including the lengths of the segments themselves)
    return river_shp.loc[river_shp['PFAF'].isin(pfafs_start2end),'Length'].sum()



###################
# 2. Functions to select reservoir pfafs and corresponding outlets

def get_pfafs_res(river_shp): 
    """return list of reservoir pfaf codes based on islake variable (to be replaced!!)"""
# Select reservoirs  based on islake. 
    return list(river_shp.loc[river_shp['islake']==1,'PFAF'].values)

def get_outlets(river_shp,pfaf_reservoirs,threshold): 
    """Get pfafstetter codes of outlets (end of influence where reservoir serves)
        based on: 1. river mouth is reached (there is no downstream pfaf code)
                  2. Next reservoir on stream network is reached 
                  3. length threshold is exceeded (only inlcude segments within threshold)
                  
        input: gdf of river shps, list of pfaf codes of reservoirs and length threshold (in m)
        output: list of outlets corresponding to reservoir list"""
    import math
    pfaf_outlets = []
    count = 1
    print('-------- Searching for outlets --------')
    print('')
    for pfaf_res in pfaf_reservoirs: 
        print('processing '+str(count)+ ' of '+str(len(pfaf_reservoirs)))
        print('reservoir: '+pfaf_res)
        # get downstream lookup table 
        downstream_lookup = get_downstream_PFAF(river_shp)

        # initialise
        pfaf_current = pfaf_res
        total_length = 0
        outlet_found = False


        # travel downstream from reservoir pfaf to identify outlet (end of reservoir influence)
        while not outlet_found: 

            # next downstream segment
            pfaf_down = downstream_lookup[downstream_lookup['PFAF']==pfaf_current]['PFAF_downstream'].values[0]
            
            # res pfaf has no downstream
            if math.isnan(float(pfaf_down)): 
                pfaf_outlet = pfaf_current
                outlet_found = True
                print('reservoir has no downstream')
            else: 
                # add length of downstream segment to total length 
                total_length = total_length + river_shp[river_shp['PFAF']==pfaf_down]['Length'].values[0]

                # check if pfaf_current is outlet: 

                # 1. pfaf is river mouth (no outlet downstream)
                if math.isnan(float(downstream_lookup[downstream_lookup['PFAF']==pfaf_down]['PFAF_downstream'].values[0])): 
                    print('river mouth reached')
                    print('')
                    pfaf_outlet = pfaf_down
                    outlet_found = True

                # 2. pfaf is other reservoir
                elif pfaf_down in pfaf_reservoirs:
                    print('reservoir reached')
                    print('')
                    pfaf_outlet = pfaf_down 
                    outlet_found = True

                # 3. length threshold is exceeded for downstream segment (so include current)
                elif total_length > threshold:
                    print('length threshold exceeded')
                    print('')
                    pfaf_outlet = pfaf_current
                    outlet_found = True

                # move downstream
                else: 
                    pfaf_current = pfaf_down
        pfaf_outlets.append(pfaf_outlet)
        count=count+1
    print('---------- All outlets found! ----------')

    return pfaf_outlets



###################
# 3. Functions to determine conditions for topology selection

def get_pfaf_betweenres(pfaf_res,pfaf_outlet,all_pfafs, include_end=True): 
    """
    return list of pfaf codes from all_pfafs that are located between pfaf_res (upper) and pfaf_outlet (lower)
    including upper reach (providing reservoir), but excluding outlet

    """
    # find all pfaf codes of segments upstream of Palisades
    pfaf_subbasin = pfaf.get_subbasin(all_pfafs, pfaf_outlet, include_closed=False)

    # identify upstream basin of reservoir and exclude from pfafs_subbasin
    pfaf_betweenres= []
    for i in pfaf_subbasin: 
        if not pfaf.check_upstream(i,pfaf_res) and i!=pfaf_res: # only include reaches downstream of reservoir 
            pfaf_betweenres.append(i)
    if include_end: 
        pfaf_betweenres.append(pfaf_res)
    return pfaf_betweenres

def exclude_higher_botelev(pfafs, pfaf_res, seg_topo): 
    """ Exclude segments that have higher BottomElev than reservoir BotElev
    input: pfafs: list with all segments
            pfaf_res: reservoir pfaf code
            seg_topo: (geo)pandas dataframe with BotElev and PFAF per segment
    output: pfaf_belowres: list with pfafs where higher botelev are excluded. 
    """
    botelev_res = seg_topo.loc[seg_topo['PFAF']==pfaf_res,'BotElev'].values[0]
    pfaf_belowres = seg_topo.loc[seg_topo.PFAF.isin(pfafs) & (seg_topo.BotElev < botelev_res), 'PFAF']
    
    return pfaf_belowres


def exclude_upstream_other_res(pfafs, pfaf_reservoirs):
    """ 
    Exlcude segments upstream of other reservoirs and segments of other reservoirs themselves as well
    input: pfafs: list of current selection of pfafs, pfaf_reservoir (list of all reservoir pfafs). 
    output:  pfafs_without_upstream_other_res
    """
    # identify pfafs of reservoirs which are in currently identified network (to reduce length of loop below)
    pfaf_reservoirs_network = [] # pfaf_reservoirs in pfaf_belowres 
    for p in pfaf_reservoirs: 
        if p in pfafs: 
            pfaf_reservoirs_network.append(p) 

            # identify pfafs upstream of other reservoirs (-> to exclude)
    pfafs_upstream_other_res = []
    for i in pfafs: 
        for p in pfaf_reservoirs_network: 
            if pfaf.check_upstream(i,p) or (i==p): # exclude upstream of other reservoir and other reservoir itself
                pfafs_upstream_other_res.append(i)

    # exclude those values
    pfafs_without_upstream_other_res = []
    for i in pfafs: 
        if not i in pfafs_upstream_other_res: 
            pfafs_without_upstream_other_res.append(i)

    return pfafs_without_upstream_other_res

def godown_tributaries(tributary,tributaries_res,river_shp,length_trib,threshold):
    """go down on tributaries and apply length condition (recursive function)
    output: list of dependent tributaries"""

    # get downstream values for all pfafs
    downstream_lookup = get_downstream_PFAF(river_shp)

    # get all subtributaries of segment that flows in resstream
    all_subtributaries = pfaf.get_tributaries(tributary, 1)

    for i in range(1,len(all_subtributaries.keys())): 
        pfaf_inresstream = list(all_subtributaries.keys())[i]

        # add segment to reservoir tributaries (serverd by reservoir) # HERE THE CONDITION of first tributary to main stem can be adjusted. 
        tributaries_res.append(pfaf_inresstream)

        # go down in subtributaries and apply threshold length condition
        subtributaries_considered = []

        for key in all_subtributaries[pfaf_inresstream].keys():
            # save considered tributaries to not loop over them twice
            subtributary = all_subtributaries[pfaf_inresstream][key]

            # check if substream tributary flows into tributary
            downstream_subtributary = downstream_lookup[downstream_lookup['PFAF'].isin(subtributary)]['PFAF_downstream'].values.tolist() # get all substream pfafs for subtributaries
            intributary =  [True for trib in downstream_subtributary if trib in tributaries_res]
            #if not intributary: print(subtributary)
            if subtributary == '31249843': print('subtribfound')    
            if (subtributary not in subtributaries_considered) and (intributary): 

                if len(subtributary) == 1: 
                    length_trib_withsubtrib1 = length_trib + river_shp.loc[river_shp['PFAF'] == subtributary[0],'Length'].values[0]
                    if length_trib_withsubtrib1 < threshold: 
                        tributaries_res.append(subtributary[0])
                        # also include other subtributary (part of main stream) if exists   
                        tributary_upstream = str(int(subtributary[0])+1) 
                        if (tributary_upstream in river_shp['PFAF'].values): 
                            length_trib_withsubtrib2 = length_trib + river_shp.loc[river_shp['PFAF'] == tributary_upstream,'Length'].values[0]
                            if (length_trib_withsubtrib2 < threshold): 
                                tributaries_res.append(tributary_upstream)

                    return tributaries_res
                else: 
                    # find tributaries again
                    tributaries_res = godown_tributaries(subtributary,tributaries_res,river_shp,length_trib,threshold)            

            subtributaries_considered.append(subtributary)
    return tributaries_res

        
  

 # this function replaces the get_pfaf_between_res 
def get_pfaf_resstream_and_in_tributary_threshold(pfaf_res, pfaf_outlet, river_shp, threshold):
    """ get reservoir stream and all tributaries on threshold distance from reservoir main stream
    input:   pfaf_res, pfaf_outlet, river_shp, threshold
    output: tributaries_res: list with all pfafs that are on mainstream or on tributary wihtin treshold distance from reservoir main stream
    """  
    resstream = get_pfafs_start2end(pfaf_res,pfaf_outlet,river_shp, include_end=False)

    # get lookup list of downstream PFAFS
    downstream_lookup = get_downstream_PFAF(river_shp)

    tributaries_considered = [] # all considered tributaries
    tributaries_res = [] # all tributaries dependend on reservoir
    count = 1
    for pfaf_resstream in resstream: 
        print('processing tributaries of reach '+str(count)+' of '+str(len(resstream)), end='\r')

        if pfaf_resstream != pfaf_res: 
            # determine tributaries
            subbasin = pfaf.get_subbasin(river_shp['PFAF'], pfaf_resstream ,include_closed = False)
            all_tributaries = pfaf.get_tributaries(subbasin, 1)

            # check if dictionary is not empty
            if all_tributaries:
                for key in all_tributaries[pfaf_resstream].keys():
                    # get list of pfafs in one tributary
                    tributary = all_tributaries[pfaf_resstream][key]

                    # get list of downstream pfafs of tributary pfafs
                    downstream_tributary = downstream_lookup[downstream_lookup['PFAF'].isin(tributary)]['PFAF_downstream'].values.tolist()
                    # check if one or more segments tributary flow directly into resstream
                    inresstream =  [True for trib in downstream_tributary if trib in resstream]

                    # get only new tributaries and tributaries that directly flow into resstream
                    if (tributary not in tributaries_considered) and inresstream:
                        # NEXT TO DO only include if one pfaf of tributary list flows into resstream

                        if len(tributary)==1: 
                            
                            # save length of tributary river
                            trib_length = river_shp.loc[river_shp['PFAF'] ==tributary[0],'Length']

                            # add segment to reservoir demand (possibly here we can put a condition on only adding segment if it is below certain length)
                            tributaries_res.append(tributary[0])

                            # POSSIBLY: extend condition based on length of segment? TO be determined. 
                        else: 

                            # get all subtributaries of segment that flows in resstream
                            all_subtributaries = pfaf.get_tributaries(tributary, 1)
                            for i in range(len(all_subtributaries.keys())): 
                             
                                pfaf_inresstream = list(all_subtributaries.keys())[i]

                                # initialise length_trib with first segment length  
                                length_trib = river_shp.loc[river_shp['PFAF'] == pfaf_inresstream,'Length'].values[0]

                                # add segment to reservoir tributaries (serverd by reservoir) # HERE THE CONDITION of first tributary to main stem can be adjusted. 
                                tributaries_res.append(pfaf_inresstream)

                                # go down in subtributaries and apply threshold length condition
                                subtributaries_considered = []

                                for key in all_subtributaries[pfaf_inresstream].keys():
                                    # get list of pfafs in one tributary
                                    subtributary = all_subtributaries[pfaf_inresstream][key]

                                   # check if substream tributary flows into tributary
                                    downstream_subtributary = downstream_lookup[downstream_lookup['PFAF'].isin(subtributary)]['PFAF_downstream'].values.tolist() # get all substream pfafs for subtributaries
                                    # concatenate resstream list and tributaries_res (if filled)
                                    considered = resstream+tributaries_res if tributaries_res else resstream
                                    intributary =  [True for trib in downstream_subtributary if trib in considered]

                                    if (subtributary not in subtributaries_considered) and intributary:
                                        subtributaries_considered.append(subtributary)

                                        if len(subtributary) == 1: 
                                            length_trib = length_trib + river_shp.loc[river_shp['PFAF'] == subtributary[0],'Length'].values[0]
                                            if length_trib < threshold: 
                                                tributaries_res.append(subtributary[0])
                                        else: 
                                            # find tributaries again
                                            tributaries_res = godown_tributaries(subtributary,tributaries_res,river_shp,length_trib,threshold)  
                        # save considered tributaries to not loop over them twice
                        tributaries_considered.append(tributary)
                    
                    # make sure to not include pfaf_res in tributaries_res
                    if pfaf_res in tributaries_res: tributaries_res.remove(pfaf_res)
        count = count+1
    return tributaries_res + resstream
           




###################
# 4. Functions to determine topology based on conditions


def get_seg_dependency(pfaf_reservoirs, pfaf_outlets, seg_topo, threshold): 
    
    """Create dependency dictionary with for each river segment pfaf corresponding reservoir pfafs 
        based on the following rules
        1. Segments must be downstream of reservoir
        2. Bottom elevation of segment cannot be higher than bottom elevation of reservoir segment
        
        seg_topo is pandas df with information on river topology. 
        """
    # initialise dictionary with reservoir dependency per river
    dependency_dict = {}
    count = 1
    # loop over reservoir and their respective outlets
    for pfaf_res,pfaf_outlet in zip(pfaf_reservoirs, pfaf_outlets):

        print('processing reservoir '+str(count)+ ' of '+str(len(pfaf_reservoirs)))
        ### CONDITIONS TO SELECT SEGMENTS PER RESERVOIR based on topology
        
        # find all reaches between (res) and (outlet)
        #pfaf_res2outlet = get_pfaf_betweenres(pfaf_res,pfaf_outlet,seg_topo['PFAF'],include_end=False)
        
        # find all reaches between res and outlet with distance threshold for second order tributaries. 
        pfaf_in_threshold = get_pfaf_resstream_and_in_tributary_threshold(pfaf_res, pfaf_outlet, seg_topo, threshold)
        
        # Exclude segments that have higher BottomElev than reservoir BotElev
        pfaf_belowres = exclude_higher_botelev(pfaf_in_threshold, pfaf_res, seg_topo)

        # Exlcude segments upstream of other reservoirs and segments of other reservoirs themselves as well
        #pfafs_dependent = exclude_upstream_other_res(pfaf_belowres.values,pfaf_reservoirs)
        
        # define here name of final selection
        pfafs_selection = pfaf_belowres

        #### END CONDITIONS
        
        # get pfaf codes of all segments depenent on reservoir
        pfaf_segment_res = seg_topo.loc[seg_topo.PFAF.isin(pfafs_selection),'PFAF'].values

        for pfaf_seg in pfaf_segment_res: 
            # check if river pfaf has already dependend res
            if pfaf_seg in dependency_dict: 
                # append to already existent list
                dependency_dict[pfaf_seg] =  dependency_dict[pfaf_seg] + [pfaf_res]
            else: 
                dependency_dict[pfaf_seg] = [pfaf_res]
        count = count+1
    print('')
    return dependency_dict

def get_res_dependency(pfaf_reservoirs, seg_dependency_dict): 
    
    """Get list of segments dependend per reservoir
    input:  1. list with pfaf codes of reservoirs, 
            2. dictionary with per segment dependend reservoirs 
            3. dictionary with weights per segment dependend reservoir
    output: res_dependency_dict[pfaf_res] = [pfaf_to_sum]
    dictionary with first level keys: reservoir pfafs
                            values: contributing river segments pfafs
                            """


    # dictionary to store dependency per reservoir
    res_dependency_dict = {}

    for pfaf_res in pfaf_reservoirs: 

        pfafs_to_sum = []    # initialise list of segments to sum by
        for pfaf_seg, lookup_pfaf_res in seg_dependency_dict.items():

            if pfaf_res in lookup_pfaf_res: 
                pfafs_to_sum.append(pfaf_seg)

        # save dependend segments per reservoir in dict 
        res_dependency_dict[pfaf_res] = pfafs_to_sum
    return res_dependency_dict



# function contained within function
def get_nhrus_per_res(pfaf_reservoirs, seg_dependency_dict):
    """Calculcate number of segments per reservoir 
        output: dictionary with keys: reservoir pfaf, value: number of segments dependend on it""" 

    # get hrus per reservoir ordered
    res_dependency_dict = get_res_dependency(pfaf_reservoirs, seg_dependency_dict)

    nseg_res_dict = {}
    nseg_res = []
    for pfaf_res in pfaf_reservoirs: 

        nseg_res.append(len(res_dependency_dict[pfaf_res]))
        nseg_res_dict[pfaf_res] = len(res_dependency_dict[pfaf_res])

    return  nseg_res_dict



def get_weights_per_seg(dependency_dict, seg_topo, pfaf_reservoirs, weigh_smax_with_nseg = True):
    """Assign weights to dependend reservoirs for each river segment
    returns dict with river segments and respective res weights. 
    !!! look up df has to be replaced with GRanD maximal storages
    """
    
    if weigh_smax_with_nseg: # if weighting, get number of segments per reservoir
        nseg_res_dict = get_nhrus_per_res(pfaf_reservoirs, dependency_dict)
    
    weights_dict = {}

    # identify segments with multiple reservoirs
    for pfaf_seg in dependency_dict: 
        print('processing segment '+str(pfaf_seg)+ ' of '+str(len(dependency_dict)),end='\r' )
        if len(dependency_dict[pfaf_seg])>1: 
            max_storages = [] # initialise list wich is going to be divided by list of total storages
            for pfaf_res in dependency_dict[pfaf_seg]: 
                # create list of maximum reservoirs storages based on PFAF codes
                
                # weigh max storage with number of segments the reservoir is providing to
                if weigh_smax_with_nseg: 
                    nseg =  nseg_res_dict[pfaf_res]
                    max_storage = seg_topo.loc[seg_topo['PFAF'] == pfaf_res,'lake_Vol'].values[0] / nseg

                else: # just take max storages 
                    max_storage = seg_topo.loc[seg_topo['PFAF'] == pfaf_res,'lake_Vol'].values[0]
                max_storages.append(max_storage)
            # calculate weights and save in new dictionary
            weights_dict[pfaf_seg] = (max_storages/sum(max_storages)).tolist()
            
    return weights_dict
  
    
def get_weights_per_seg(dependency_dict, seg_topo, pfaf_reservoirs, weigh_smax_with_nseg = True):
    """Assign weights to dependend reservoirs for each river segment
    returns dict with river segments and respective res weights. 
    !!! look up df has to be replaced with GRanD maximal storages
    """
    
    if weigh_smax_with_nseg: # if weighting, get number of segments per reservoir
        nseg_res_dict = get_nhrus_per_res(pfaf_reservoirs, dependency_dict)
    
    weights_dict = {}
    count = 1
    # identify segments with multiple reservoirs
    for pfaf_seg in dependency_dict: 
        
        print('processing segment '+str(count)+ ' of '+str(len(dependency_dict.keys())),end='\r' )
        if len(dependency_dict[pfaf_seg])>1: 
            max_storages = [] # initialise list wich is going to be divided by list of total storages
            for pfaf_res in dependency_dict[pfaf_seg]: 
                # create list of maximum reservoirs storages based on PFAF codes
                
                # weigh max storage with number of segments the reservoir is providing to
                if weigh_smax_with_nseg: 
                    nseg =  nseg_res_dict[pfaf_res]
                    max_storage = seg_topo.loc[seg_topo['PFAF'] == pfaf_res,'lake_Vol'].values[0] / nseg

                else: # just take max storages 
                    max_storage = seg_topo.loc[seg_topo['PFAF'] == pfaf_res,'lake_Vol'].values[0]
                max_storages.append(max_storage)
            # calculate weights and save in new dictionary
            weights_dict[pfaf_seg] = (max_storages/sum(max_storages)).tolist()
            count = count+1
    print('')
    return weights_dict
  

def get_res_dependency_and_weights(pfaf_reservoirs, seg_dependency_dict, weights_dict): 
    
    """Get list of segments dependend per reservoir, both in pfaf as in weights
    input:  1. list with pfaf codes of reservoirs, 
            2. dictionary with per segment dependend reservoirs 
            3. dictionary with weights per segment dependend reservoir
    output: res_dependency_dict[pfaf_res][pfaf_to_sum] = weight
    dictionary with first level keys: reservoir pfafs
                            second level keys: contributing river segments pfafs
                            values: correspoding weights 
                            """
    

    # dictionary to store dependency per reservoir
    res_dependency_dict = {}

    for pfaf_res in pfaf_reservoirs: 

        pfafs_to_sum = []    # initialise list of segments to sum by
        weights_to_sum  = [] # initialise list of segments to sum by

        for pfaf_seg, lookup_pfaf_res in seg_dependency_dict.items():

            if pfaf_res in lookup_pfaf_res: 
                pfafs_to_sum.append(pfaf_seg)

                # get weights (if only reservoir weight = 1 )
                if len(lookup_pfaf_res) == 1: 
                    weight=1
                else: 
                    weight = weights_dict[pfaf_seg][lookup_pfaf_res.index(pfaf_res)]

                weights_to_sum.append(weight)
        # save dependend segments per reservoir in dict 
        res_dependency_dict[pfaf_res] = dict(zip(pfafs_to_sum, weights_to_sum)) 

    return res_dependency_dict

    
    

###################
# 5. Functions to calculate irrigation demand
 
def calc_demand_per_res(hrus_shp, res_topo_dict): 
    """Calculate irrigation demand per reservoir
    input: dataframe with hru pfaf codes and irrigation demand ('QIRRIGmean')
            dictionary with reservoir dependencies and weights
    output: dictionary with reservoirs and corresponding water demands
    """
    demand_dict = {}
    count = 1
    for reservoir in res_topo_dict: 
        print('calculating demand for reservoir '+str(count)+ ' of '+str(len(res_topo_dict)),end='\r' )
        demand_tot = 0
        # calculate sum over all river segments with weights applied
        for seg_tosum in res_topo_dict[reservoir]:

            # get weight to apply to segment
            seg_weight = res_topo_dict[reservoir][seg_tosum]
            if seg_tosum in hrus_shp['PFAF'].values: 
                demand_tosum = hrus_shp.loc[hrus_shp['PFAF']==seg_tosum, 'basinIrrig'].values[0] * seg_weight
                demand_tot = demand_tot + demand_tosum
        # store in dict per reservoir
        demand_dict[reservoir] = demand_tot
        count = count+1
    return demand_dict
