"""
Script that creates a file for prescribed polar heating based on Orlanski and Solman (2010) and Ruth Geen's H-S work
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

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

H = 8
p0 = 1000

filename = polar_heating()
plot_polar_heating(filename)
