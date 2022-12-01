"""
Script that checks column integrated heat input of prescribed polar heating
"""
from glob import glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from scipy import integrate
from aostools.constants import *

def open_heat_static(filename):
    file = '/home/links/rm811/Isca/input/' + folder + '/' + filename + '.nc'
    ds = xr.open_dataset(file)
    lon = ds.lon
    heat = ds.variables[filename].transpose('lat','lon','pfull')
    return lon, heat

def open_heat_evolution(filename):
    folder = '/disco/share/rm811/isca_data/PK_e0v4z13_' + filename
    files = sorted(glob(folder+'/run*/*.nc'))
    files_subset = files[:12]
    ds = xr.open_mfdataset(files_subset)
    t = ds.time.data
    heat = ds.local_heating
    return t, heat

def diff(heat, var):
    return heat.differentiate(var, edge_order=2)

def integral(heat, ticks, ax):
    return integrate.trapezoid(y = heat, x = ticks, axis = ax)

def time_integrate(heat, t):
    dTdt = diff(heat, 'time')
    return integral(dTdt, t, 0)

def surf_integrate(heat, lon, p):
    n = len(lon)
    if n == 1:
        heat = heat.mean('lon')
        s = integral(heat, p, 1)
        surf = s
        for j in range(127):
            surf = np.vstack((surf, s))
        surf = surf.transpose() 
    elif n > 1:
        surf = integral(heat, p, 2)
    return surf

def ratio(x0, x1):
    return x1/x0

def integration_levels(filenames, level):
    # set-up the 3D grid first
    file = '/home/links/rm811/Isca/input/dimensions_3d.nc'
    ds = xr.open_dataset(file, decode_times=False)
    lat = ds.lat
    coslat = np.cos(np.deg2rad(lat.data))
    lon = ds.lon
    p = ds.pfull

    n = len(filenames)
    int1 = []
    for i in range(n):
       heat_lon, heat = open_heat_static(filenames[i])
       int1.append(surf_integrate(heat, heat_lon, p))

    int2 = []
    for i in range(n):
        int2.append(integral(int1[i], lon, 1))

    if level == 'lat': 
        finished = int2
    
    elif level == 'full':
        int3 = []
        for i in range(n):
            int3_sub = []
            for j in range(len(int2[i])):
                int3_sub.append((int2[i][j]*coslat[j]))
            int3.append(np.sum(int3_sub)/np.sum(coslat))
        finished = int3
        print(finished)

    return n, lat, finished

def plot_vslat(filenames, labels, colors):
    n, lat, int = integration_levels(filenames, 'lat')
    fig, ax = plt.subplots(figsize=(6,6))
    for i in range(n):
        ax.plot(lat, int[i], color=colors[i], label=labels[i])
    ax.set_xlim(0, max(lat))
    ax.set_xlabel(r'latitude ($\degree$N)', fontsize='x-large')
    ax.set_ylim(bottom=0)
    ax.set_ylabel(r'integrated heat (W m$^{-2}$)', fontsize='x-large')
    plt.legend(fontsize='x-large')
    ax.tick_params(axis='both', labelsize = 'x-large', which='both', direction='in')
    plt.savefig('heatcheck1_polar.pdf', bbox_inches = 'tight')
    return plt.close()

def plot_vsdefault(filenames, labels, colors):
    n, lat, int = integration_levels(filenames, 'full')
    ratios = []
    for i in range(n):
        ratios.append(ratio(int[0], int[i]))
    print(ratios)
    fig, ax = plt.subplots(figsize=(8,6))
    ax.bar(labels[1:], ratios[1:], color=colors[1:])
    ax.set_xticklabels(labels[1:], rotation=30, ha='right')
    ax.axhline(ratios[0], color=colors[0], linewidth=1, linestyle='--')
    ax.text(0, ratios[0]+0.02, labels[0], color=colors[0], fontsize='x-large')
    ax.set_ylim(bottom=0)
    ax.set_ylabel(r'integrated heat ratio', fontsize='x-large')
    ax.tick_params(axis='both', labelsize = 'x-large', which='both', direction='in')
    plt.savefig('heatcheck2_polar.pdf', bbox_inches = 'tight')
    return plt.close()

def plot_multiheat():
    basis = 'PK_e0v4z13_'
    perturb = 'q6m2y45l800u200'
    T_ctrl = xr.open_dataset('/disco/share/rm811/processed/'+basis+perturb+'_Ttz.nc', decode_times=False).temp[0]
    heat_names = ['w15a4p800f800g50', 'w15a4p400f800g50', 'w30a4p800f800g50']
    heats = []
    for i in range(len(heat_names)):
        lat, lon, p, heat = open_heat_static(heat_names[i])
        heats.append(heat.mean('lon').transpose())

    h = 4/86400
    inc = 0.25e-5
    h_range = np.arange(0, h+inc, inc)

    # Plot
    plt.figure(figsize=(6.5,6))
    plt.contour(lat, T_ctrl.pfull, T_ctrl, colors='gainsboro', levels=np.arange(160, 330, 10), linewidths=1)
    cs = plt.contourf(lat, p, heats[0], cmap='Reds', levels=h_range)
    cb = plt.colorbar(cs, location='right')
    cb.set_label(label=r'Heating (K s$^{-1}$)', size='xx-large')
    cb.ax.tick_params(labelsize='xx-large')
    plt.contour(lat, p, heats[1], colors='#4D0099', levels=h_range, linewidths=1, alpha=0.25)
    plt.contour(lat, p, heats[2], colors='#4D0099', levels=h_range, linewidths=1, alpha=0.25)
    plt.xlim(0, 90)
    plt.xticks([10, 30, 50, 70, 90], ['10', '30', '50', '70', '90'])
    plt.xlabel(r'Latitude ($\degree$N)', fontsize='xx-large')
    plt.ylim(max(p), 100)
    plt.yscale('log')
    plt.ylabel('Pressure (hPa)', fontsize='xx-large')
    plt.tick_params(axis='both', labelsize = 'xx-large', which='both', direction='in')
    plt.savefig('heatings.pdf', bbox_inches = 'tight')
    return plt.close()

folder = 'polar_heating'
filenames = ['w15a4p800f800g50', 'w10a4p800f800g50_scaled', 'w20a4p800f800g50_scaled', 'w30a4p800f800g50_scaled', 'w15a4p600f800g50', 'w15a4p400f800g50']
labels = ['default', r'$\phi_{w} = 10\degree$', r'$\phi_{w} = 20\degree$', r'$\phi_{w} = 30\degree$', r'$p_{top} = 600$ hPa', r'$p_{top} = 400$ hPa']
colors = ['k', '#B30000', '#FF9900', '#FFCC00', '#00B300', '#0099CC', '#4D0099', '#CC0080']

plot_vslat(filenames, labels, colors)
plot_vsdefault(filenames, labels, colors)
