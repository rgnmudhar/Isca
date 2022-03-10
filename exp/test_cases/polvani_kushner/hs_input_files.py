# -*- coding: utf-8 -*-s
import numpy as np
#from stephen import create_timeseries as cts
import xarray as xr
#from data_handling_updates import model_constants as mc
import matplotlib.pyplot as plt
import sh
from pylab import rcParams

"""
def hs_default_forcing(t_zero=315., t_strat=200., delh=60., delv=10., eps=0., sigma_b=0.7, pref=1.e5, write_netcdf=True):
    
    ozone_file = '/Isca/input/rrtm_input_files/ozone_1990_notime.nc'
    data = xr.open_dataset( ozone_file, decode_times=False)
    
    template = data.ozone_1990 * 0. + 1.
    
    sin_lat = np.sin(data.lat * np.pi/180.)
    cos_lat = np.cos(data.lat * np.pi/180.)
    sin_lat_2 = np.sin(data.lat * np.pi/180.) ** 2.
    cos_lat_2 = np.cos(data.lat * np.pi/180.) ** 2.
    
    t_star = (t_zero - delh*sin_lat_2 - eps*sin_lat) * template
    
    tstr = (t_strat - eps*sin_lat) * template
    
    p_norm = 100.*data.pfull/pref
    
    the = t_star - delv*cos_lat_2 * np.log(p_norm)
    teq = the * p_norm ** mc.kappa
    
    teq = teq.where(teq > tstr, tstr)
    #teq.mean('lon').plot.contourf(x='lat', y='pfull', yincrease=False, levels=np.arange(200.,316.,5.,))
    #plt.show()
    
    # NB filename should be 32 characters or less
    
    coord_list = ["pfull", "lat", "lon"]
    teq = xr.Dataset(
         data_vars=dict(
             teq = (coord_list, teq.transpose('pfull','lat','lon').values)
         ),
         coords=data.coords
    )
    
    if write_netcdf:
        teq.to_netcdf('/disco/share/rg419/ArctiCONNECT_data/teq.nc', format="NETCDF3_CLASSIC",
             encoding = {"teq": {"dtype": 'float32', '_FillValue': None},
                        "lat": {'_FillValue': None}, "lon": {'_FillValue': None},
                        "latb": {'_FillValue': None}, "lonb": {'_FillValue': None},
                        "pfull": {'_FillValue': None}, "phalf": {'_FillValue': None}}
                    )
    
    return teq
    
    

def st_forcing(t_zero=315., t_strat=200., delht=60., delhs=40., delv=10., eps=0., sigma_b=0.7, pref=1.e5, del_p=800., p_th=50., write_netcdf=True):
    
    ozone_file = 'Isca/input/rrtm_input_files/ozone_1990_notime.nc'
    data = xr.open_dataset( ozone_file, decode_times=False)
    
    template = data.ozone_1990 * 0. + 1.
    
    sin_lat = np.sin(data.lat * np.pi/180.)
    cos_lat = np.cos(data.lat * np.pi/180.)
    sin_lat_2 = np.sin(data.lat * np.pi/180.) ** 2.
    cos_lat_2 = np.cos(data.lat * np.pi/180.) ** 2.
    
    delh_p = (delht + 0.5 * (delhs - delht) * (1. + np.tanh((data.pfull - del_p)/p_th)))*template
    delh_p.mean(('lat','lon')).plot()
    plt.show()
    
    t_star = (t_zero - delh_p*sin_lat_2 - eps*sin_lat) * template
    
    tstr = (t_strat - eps*sin_lat) * template
    
    p_norm = 100.*data.pfull/pref
    
    the = t_star - delv*cos_lat_2 * np.log(p_norm)
    teq = the * p_norm ** mc.kappa
    
    teq = teq.where(teq > tstr, tstr)
    
    #teq.mean('lon').plot.contourf(x='lat', y='pfull', yincrease=False, levels=np.arange(200.,316.,5.,))
    #plt.show()
    
    
    coord_list = ["pfull", "lat", "lon"]
    teq = xr.Dataset(
         data_vars=dict(
             teq = (coord_list, teq.transpose('pfull','lat','lon').values)
         ),
         coords=data.coords
    )
    
    
    # NB filename should be 32 characters or less
    filename = 'dht' + str(delht) + 'dhs' + str(delhs) + 'p' + str(del_p) + 'g' + str(p_th)
    print(len(filename))
    #filename='heating_test'
    teq = teq.rename({"teq" : filename})
    
    
    if write_netcdf:
        teq.to_netcdf('/disco/share/rg419/ArctiCONNECT_data/'+ filename + '.nc', format="NETCDF3_CLASSIC",
             encoding = {filename: {"dtype": 'float32', '_FillValue': None},
                        "lat": {'_FillValue': None}, "lon": {'_FillValue': None},
                        "latb": {'_FillValue': None}, "lonb": {'_FillValue': None},
                        "pfull": {'_FillValue': None}, "phalf": {'_FillValue': None}}
                    )
    
    return teq
"""    
    



def polar_heating(y_wid=15., th_mag=4., p_top = 800., p_th = 50., p_ref=800., save_output=True):
    
    # Parameter sweep
    # 1. Vary p_top - depth of forcing: 0:200:800 (1000 would be no forcing!)
    # 2. Vary th_mag - magnitude of forcing in K/s when centred on 800hPa  (0.5,1.,1.5,2.) Could try negative values too
    # 3. Vary y_wid - decay of forcing away from pole (10., 15., 20.)
    # 4. Vary p_th - sets vertical gradient of forcing at cap (25.,50.,75.) Sensitivity check that steepness of transition doesn't cause unexpected behaviour
    
    ozone_file = '/home/links/rm811/Isca/input/rrtm_input_files/ozone_1990_notime.nc'
    data = xr.open_dataset(ozone_file, decode_times=False)
    
    template = data.ozone_1990 * 0. + 1.
    
    # Vary with latitude in similar way to Orlanski and Solman 2010
    #heat_lat = np.exp(-((data.lat - 90.)/y_wid)**2.) + np.exp(-((data.lat + 90.)/y_wid)**2.) * template
    heat_lat = np.exp(-((data.lat - 90.)/y_wid)**2.) * template # for north pole only
        
    # fix so that the function has magnitude 1 when p_top = p_ref, and otherwise scales to give constant net energy input
    # p_ref can be varied to alter the total input, which I think will be (1000-p_ref) * cp/g * th_mag
    
    if p_top==0.:  # If the heating is going right to the model top then make heating uniform in height, scaled to fit p_ref level
        polar_heating = th_mag * (1000. - p_ref)/1000. * heat_lat /86400.
    else:
        polar_heating = 0.5 * th_mag * (1000. - p_ref)/(1000. - p_top) * heat_lat * (1. + np.tanh((data.pfull - p_top)/p_th)) /86400.
    
    coord_list = ["pfull", "lat", "lon"]
    polar_heating = xr.Dataset(
         data_vars=dict(
             polar_heating = (coord_list, polar_heating.transpose('pfull','lat','lon').values)
         ),
         coords=data.coords
    )
    
    if save_output:
        # NB filename should be 32 characters or less
        filename = 'w' + str(int(y_wid)) + 'a' + str(int(th_mag)) + 'p' + str(int(p_top)) + 'f' + str(int(p_ref)) + 'g' + str(int(p_th))
        print(len(filename))
        #filename='heating_test'
        polar_heating = polar_heating.rename({"polar_heating" : filename})
        
        polar_heating.to_netcdf('/home/links/rm811/Isca/input/polar_heating/' + filename + '.nc', format="NETCDF3_CLASSIC",
             encoding = {filename: {"dtype": 'float32', '_FillValue': None},
                    "lat": {'_FillValue': None}, "lon": {'_FillValue': None},
                    "latb": {'_FillValue': None}, "lonb": {'_FillValue': None},
                    "pfull": {'_FillValue': None}, "phalf": {'_FillValue': None}}
                )
    
    return filename

def altitude(p):
    """Finds altitude from pressure using z = -H*log10(p/p0) """
        
    z = np.empty_like(p)
    
    for i in range(p.shape[0]):
        z[i] = -H*np.log((p[i])/p0)
        
    # Make into an xarray DataArray
    z_xr = xr.DataArray(z, coords=[z], dims=['pfull'])
    z_xr.attrs['units'] = 'km'
    
    #below is the inverse of the calculation
    #p[i] = p0*np.exp((-1)*z[i]*(10**3)/((R*T/g)))
    
    return z_xr

def plot_polar_heating(filename):
    
    name = filename
    file = '/home/links/rm811/Isca/input/polar_heating/' + name + '.nc'
    ds = xr.open_dataset(file)

    lat = ds.coords['lat'].data
    lon = ds.coords['lon'].data
    p = ds.coords['pfull'].data


    upper_p = ds.coords['pfull'].sel(pfull=1, method='nearest') # in order to cap plots at pressure = 1hPa
    z = altitude(p)
    upper_z = -H*np.log(upper_p/p0)

    heat = ds.variables[name]
    heatz = heat.mean(dim='lon').data   
    
    # Plot for comparison with Screen & Simmonds (2010)
    plt.figure(figsize=(10,8))
    plt.contourf(lat, p, heatz, cmap='Reds', levels=np.arange(0, 5e-5, 2e-6))
    plt.xlim(0, 90)
    plt.ylim(max(p), 100)
    #plt.yscale('log')
    plt.colorbar(label="Heating (K/s)")
    plt.xlabel(r'Latitude ($\degree$N)', fontsize='x-large')
    plt.ylabel("Pressure (hPa)", fontsize='x-large')
    plt.title(name, fontsize='x-large')
    plt.tick_params(axis='both', labelsize = 'x-large', which='both', direction='in')

    return plt.show()
    
"""
def plot_polar_heating(th_mag=1., del_ps = [0.,400.,800.]):
    
    plot_dir = '/scratch/rg419/plots/ArctiCONNECT/dry_experiments/hs_input_files/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    rcParams['figure.figsize'] = 10, 4
    rcParams['font.size'] = 11
    
    fig, ((ax1, ax2, ax3)) = plt.subplots(1, 3, sharey='row')
    axes = [ax1, ax2, ax3]
    
    for i in range(3):
        polar_ht = polar_heating(th_mag=th_mag, del_p = del_ps[i], save_output=False)
        polar_ht = polar_ht.polar_heating.isel(lon=0) * 86400.
        
        f1 = polar_ht.plot.contourf(x='lat', y='pfull', ax=axes[i], yincrease=False, levels=np.arange(0.,1.1,0.1), add_labels=False, add_colorbar=False)
        axes[i].set_xlim(0.,90.)
        axes[i].set_title(str(int(del_ps[i])) + ' hPa')
        axes[i].set_xlabel('Latitude')
        
    ax1.set_ylabel('Pressure, hPa')
    plt.subplots_adjust(left=0.07, right=0.97, top=0.9, bottom=0.1, hspace=0.3, wspace=0.15)
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.1, pad=0.2, aspect=30, shrink=0.5)
    cb1.set_label('Heating, K/day')
    
    filename = 'w' + str(15.) + 'a' + str(th_mag) + 'g' + str(50.)
    
    plt.savefig(plot_dir + filename + '.pdf', format='pdf')
    plt.close()
    
    
    
plot_polar_heating()
"""

H = 8
p0 = 1000

filename = polar_heating()
plot_polar_heating(filename)