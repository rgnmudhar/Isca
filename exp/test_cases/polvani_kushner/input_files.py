"""
Script that creates a file for prescribed polar heating based on Orlanski and Solman (2010) and Ruth Geen's H-S work
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

def plot_polar_heating(filename):
    
    name = filename
    file = '/home/links/rm811/Isca/input/polar_heating/' + name + '.nc'
    ds = xr.open_dataset(file)

    lat = ds.coords['lat'].data
    lon = ds.coords['lon'].data
    p = ds.coords['pfull'].data

    heat = ds.variables[name]
    
    # Plot for comparison with Screen & Simmonds (2010)
    plt.figure(figsize=(10,8))
    plt.contourf(lat, p, heat, cmap='Reds', levels=np.arange(0, 5e-5, 2e-6))
    plt.xlim(0, 90)
    plt.ylim(max(p), 100)
    #plt.yscale('log')
    plt.colorbar(label="Heating (K/s)")
    plt.xlabel(r'Latitude ($\degree$N)', fontsize='x-large')
    plt.ylabel("Pressure (hPa)", fontsize='x-large')
    plt.title(name, fontsize='x-large')
    plt.tick_params(axis='both', labelsize = 'x-large', which='both', direction='in')

    return plt.show()

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

def heat_perturb(q_0=6, m=2, y_cen=45, p_0=800, p_t=200, save_output=True):
    """
    Tropospheric diabatic heating perturbation to induce NH winter-like wave activity.
    Based on Lindgren et al. (2018).
    """
    # Parameters
    # 1. q_0 - magnitude of forcing in K/day
    # 2. m - longitudinal wave number (1 or 2)
    # 3. y_cen - center of heating in latitude
    # 4. p_0 and p_t - lower and upper bounds in pressure (hPa) respectively
    y_wid = 0.175*360/(2*np.pi)

    file = '/home/links/rm811/Isca/input/dimensions.nc'
    ds = xr.open_dataset(file, decode_times=False)
    
    template = ds.data * 0. + 1.
    heat_lat = np.exp(-0.5 * ((ds.lat - y_cen)/y_wid)**2.) * template
    heat_lon = np.cos(m * np.deg2rad(ds.lon)) * template
    heat_p = np.sin(np.pi * np.log(ds.pfull/p_0)/np.log(p_t/p_0)) * template

    heat = q_0 * heat_lat * heat_lon *  heat_p /86400.  # convert to K/s
    heat = heat.where(heat['pfull']<p_0, 0)
    heat = heat.where(heat['pfull']>p_t, 0)

    coord_list = ["pfull", "lat", "lon"]
    heat = xr.Dataset(
         data_vars=dict(
             heat = (coord_list, heat.transpose('pfull','lat','lon').values)
         ),
         coords=ds.coords
    )
    
    if save_output:
        # NB filename should be 32 characters or less
        filename = 'q' + str(int(q_0)) + 'm' + str(int(m)) + 'y' + str(int(y_cen)) + 'l' + str(int(p_0)) + 'u' + str(int(p_t))
        print(len(filename))
        heat = heat.rename({"heat" : filename})
        
        heat.to_netcdf('/home/links/rm811/Isca/input/asymmetry/' + filename + '.nc', format="NETCDF3_CLASSIC",
             encoding = {filename: {"dtype": 'float32', '_FillValue': None},
                    "lat": {'_FillValue': None}, "lon": {'_FillValue': None},
                    "latb": {'_FillValue': None}, "lonb": {'_FillValue': None},
                    "pfull": {'_FillValue': None}, "phalf": {'_FillValue': None}}
                )
    
    return filename

def plot_vertical(filename):
    
    name = filename
    file = '/home/links/rm811/Isca/input/asymmetry/' + name + '.nc'
    ds = xr.open_dataset(file)
    q_0 = int(name.partition("q")[2][0])/86400
    inc = 0.5e-5

    lat = ds.coords['lat'].data
    p = ds.coords['pfull'].data
    heat = ds.sel(lon=180, method='nearest').variables[name]
    
    # Plot
    plt.figure(figsize=(10,8))
    plt.contourf(lat, p, heat, cmap='Reds', levels=np.arange(0, q_0+inc, inc))
    plt.ylim(max(p), 100)
    plt.yscale('log')
    plt.colorbar(label="Heating (K/s)")
    plt.xlabel(r'Latitude ($\degree$)', fontsize='x-large')
    plt.ylabel('Pressure (hPa)', fontsize='x-large')
    plt.title(name, fontsize='x-large')
    plt.tick_params(axis='both', labelsize = 'x-large', which='both', direction='in')

    return plt.show()

def plot_horizontal(filename):
    
    name = filename
    file = '/home/links/rm811/Isca/input/asymmetry/' + name + '.nc'
    ds = xr.open_dataset(file)
    q_0 = int(name.partition("q")[2][0])/86400

    lat = ds.coords['lat'].data
    lon = ds.coords['lon'].data
    heat = ds.sel(pfull=500, method='nearest').variables[name]

    #Plot
    plt.figure(figsize=(10,8))
    plt.contourf(lon, lat, heat, cmap='RdBu_r', levels=21)
    plt.colorbar(label="Heating (K/s)")
    plt.xlabel(r'Longitude ($\degree$)', fontsize='x-large')
    plt.ylabel(r'Latitude ($\degree$)', fontsize='x-large')
    plt.title(name, fontsize='x-large')
    plt.tick_params(axis='both', labelsize = 'x-large', which='both', direction='in')

    return plt.show()

if __name__ == '__main__': 
    option = input("a) polar heating, b) heating perturbation or c) idealised topography?")

    H = 8
    p0 = 1000

    if option =='a':
        filename = polar_heating()
        plot_polar_heating(filename)
    elif option =='b':
        filename = heat_perturb()
        plot_horizontal(filename)
        plot_vertical(filename)
    elif option =='c':
        filename = ideal_topo()
        plot(filename)
