"""
Utils for plotting with geopandas 

Inne Vanderkelen  - March 2021
"""

# import modules
import matplotlib.pyplot as plt
import geopandas as gpd
import cartopy.crs as ccrs

# Functions
def set_plot_param():
    """Set my own customized plotting parameters"""
    
    import matplotlib as mpl
    mpl.rc('axes',edgecolor='grey')
    mpl.rc('axes',labelcolor='dimgrey')
    mpl.rc('xtick',color='dimgrey')
    mpl.rc('xtick',labelsize=12)
    mpl.rc('ytick',color='dimgrey')
    mpl.rc('ytick',labelsize=12)
    mpl.rc('axes',titlesize=14)
    mpl.rc('axes',labelsize=12)
    mpl.rc('legend',fontsize='large')
    mpl.rc('text',color='dimgrey')


def scale_bar(ax, length=None, location=(0.5, 0.05), linewidth=2):
    """
    ax is the axes to draw the scalebar on.
    length is the length of the scalebar in km.
    location is center of the scalebar in axis coordinates.
    (ie. 0.5 is the middle of the plot)
    linewidth is the thickness of the scalebar.
    """
    #Get the limits of the axis in lat long
    llx0, llx1, lly0, lly1 = ax.get_extent(ccrs.PlateCarree())
    #Make tmc horizontally centred on the middle of the map,
    #vertically at scale bar location
    sbllx = (llx1 + llx0) / 2
    sblly = lly0 + (lly1 - lly0) * location[1]
    tmc = ccrs.TransverseMercator(sbllx, sblly)
    #Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(tmc)
    #Turn the specified scalebar location into coordinates in metres
    sbx = x0 + (x1 - x0) * location[0]
    sby = y0 + (y1 - y0) * location[1]

    #Calculate a scale bar length if none has been given
    #(Theres probably a more pythonic way of rounding the number but this works)
    if not length: 
        length = (x1 - x0) / 5000 #in km
        ndim = int(np.floor(np.log10(length))) #number of digits in number
        length = round(length, -ndim) #round to 1sf
        #Returns numbers starting with the list
        def scale_number(x):
            if str(x)[0] in ['1', '2', '5']: return int(x)        
            else: return scale_number(x - 10 ** ndim)
        length = scale_number(length) 

    #Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbx - length * 500, sbx + length * 500]
    #Plot the scalebar
    ax.plot(bar_xs, [sby, sby], transform=tmc, color='dimgray', linewidth=linewidth)
    #Plot the scalebar label
    ax.text(sbx, sby, str(length) + ' km', transform=tmc,
            horizontalalignment='center', verticalalignment='bottom', color='dimgray',weight='bold')
    
def north_arrow(ax, location=(0.5,0.5), arrow_length=0.1): 
    x = location[0]
    y = location[1]
    ax.annotate('N', xy=(x, y), xytext=(x, y-(arrow_length+arrow_length*0.2)),
                arrowprops=dict(facecolor='dimgray', edgecolor='dimgray', width=5, headwidth=15),
                ha='center', va='center', fontsize=20, color='dimgray',
                xycoords=ax.transAxes)