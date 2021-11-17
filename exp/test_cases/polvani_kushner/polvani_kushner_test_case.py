import numpy as np
import os
from isca import DryCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 16
RESOLUTION = 'T42', 40  # T42 horizontal resolution, 25 levels in pressure - NOTE: changed from 25 to 40 pressure levels

# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = DryCodeBase.from_directory(GFDL_BASE)

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = DryCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics

exp_name = 'Polvani_Kushner_topo' # updated experiment name
exp = Experiment(exp_name, codebase=cb)

exp.inputfiles = [os.path.join(GFDL_BASE,'input/land_masks/era_land_t42.nc')]

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True) # added diagnostic field for height

diag.add_field('hs_forcing', 'teq', time_avg=True) # added diagnostic field as sanity check on Teq

exp.diag_table = diag

# define namelist values as python dictionary
# wrapped as a namelist object.
namelist = Namelist({
    'main_nml': {
        'dt_atmos': 600,
        'days': 30,
        'calendar': 'thirty_day',
        'current_date': [2000,1,1,0,0,0]
    },

    'atmosphere_nml': {
        'idealized_moist_model': False  # False for Newtonian Cooling.  True for Isca/Frierson
    },

    'spectral_dynamics_nml': {
        'damping_order'           : 4,                      # default: 2
        #'water_correction_limit'  : 200.e2,                 # default: 0 - NOTE: needs to be commented out???
        'do_water_correction': False,                       # NOTE: added but is it necessary for P-K???
        'reference_sea_level_press': 1.0e5,                  # default: 101325
        'valid_range_t'           : [50., 800.],           # default: (100, 500) - NOTE: was changed to go from 50
        'initial_sphum'           : 0.0,                  # default: 0
        'vert_coord_option'       : 'uneven_sigma',         # default: 'even_sigma'
        'scale_heights': 11.0,                               # NOTE: was changed to be 11.0 from 6.0
        'exponent': 3.0,                                    # NOTE: was changed to be 3.0 from 7.5
        'surf_res': 0.5
    },

    'spectral_init_cond_nml': { # namelist required for topography added
        'topog_file_name': 'era_land_t42.nc', #input file name
        'topography_option': 'input' # take topography from input file
    },

    # configure the relaxation profile
    'hs_forcing_nml': {
        't_zero': 315.,    # temperature at reference pressure at equator (default 315K)
        't_strat': 200.,   # stratosphere temperature (default 200K)
        'delh': 60.,       # equator-pole temp gradient (default 60K)
        'delv': 10.,       # lapse rate (default 10K)
        'eps': 10.,         # stratospheric latitudinal variation (default 0K) - NOTE: was changed to 10
        'sigma_b': 0.7,    # boundary layer friction height (default p/ps = sigma = 0.7)
        'equilibrium_t_option': 'Polvani_Kushner', # add new option for polvani_kushner relaxation

        # negative sign is a flag indicating that the units are days
        'ka':   -40.,      # Constant Newtonian cooling timescale (default 40 days)
        'ks':    -4.,      # Boundary layer dependent cooling timescale (default 4 days)
        'kf':   -1.,       # BL momentum frictional timescale (default 1 days)

        'z_ozone': 15.,     # added height of stratospheric heating source

        'do_conserve_energy':   True,  # convert dissipated momentum into heat (default True)

        'sponge_flag': True     # added sponge layer for simple damping in upper levels
    },

    'diag_manager_nml': {
        'mix_snapshot_average_fields': False
    },

    'fms_nml': {
        'domains_stack_size': 600000                        # default: 0
    },

    'fms_io_nml': {
        'threading_write': 'single',                         # default: multi
        'fileset_write': 'single',                           # default: multi
    }
})

exp.namelist = namelist
exp.set_resolution(*RESOLUTION)

#Lets do a run!
if __name__ == '__main__':
    exp.run(1, num_cores=NCORES, use_restart=False)
    for i in range(2, 13): #~1y worth
        exp.run(i, num_cores=NCORES)  # use the restart i-1 by default

