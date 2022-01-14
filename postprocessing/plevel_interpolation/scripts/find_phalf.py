"""
# The following finds phalfs from prescribed list of pressure levels
plevs = "3 16 51 138 324 676 1000 1266 2162 3407 5014 6957 9185 10000 11627 14210 16864 19534 20000 22181 24783 27331 29830 32290 34731 37173 39637 42147 44725 47391 50164 53061 56100 59295 62661 66211 70000 73915 78095 82510 85000 87175 92104 97312"
plevs_list = list(plevs.split())
plevs_list = list(map(int, plevs_list))

phalfs_list = []
phalfs_list.append(0)
for i in range(len(plevs_list)-1):
    step = (plevs_list[i+1] - plevs_list[i])/2
    phalf = int(plevs_list[i] + step)
    phalfs_list.append(phalf)
phalfs_list.append(100000)

phalfs = ""

for i in phalfs_list:
    phalfs += str(i) + " "
"""

# The following converts Isca-output sigma levels to pressure levels to interpolate to
import xarray as xr
ds = xr.open_dataset('/home/links/rm811/scratch/isca_data/PK_eps0_vtx2_zoz18_1y_heat_test/run0012/atmos_monthly.nc', decode_times=False)

p = ds.pfull.data
ph = ds.phalf.data
pfull = ""
phalf = ""

for i in p:
    pfull += str(int(i * 100)) + " "

print(pfull)

for i in ph:
    phalf += str(int(i * 100)) + " "

print(phalf)