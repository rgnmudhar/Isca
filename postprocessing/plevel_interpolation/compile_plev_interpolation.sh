#!/usr/bin/env bash
#Script that compiles the plev.x executable

cd ./exec

source $GFDL_BASE/src/extra/env/$GFDL_ENV

#compiler=${GFDL_MKMF_TEMPLATE:-ia64} #from stephen's change to .sh file
#template=mkmf.template.${compiler} #from stephen's change to .sh file  
#../bin/mkmf -p plev.x -t $GFDL_BASE/src/extra/python/isca/templates/$template -c "-Duse_netCDF" -a ../src ../src/path_names ../src/shared/mpp/include ../src/shared/include #from stephen's change to .sh file  
#../bin/mkmf -p plev.x -t ../bin/mkmf.template.ia64 -c "-Duse_netCDF" -a ../src ../src/path_names ../src/shared/mpp/include ../src/shared/include #from old version of .sh file

../bin/mkmf -p plev.x -t $GFDL_BASE/src/extra/python/isca/templates/mkmf.template.gfort -c "-Duse_netCDF" -a ../src ../src/path_names ../src/shared/mpp/include ../src/shared/include #from emily's change to .sh file plus point to fortran rather than intel compiler

make -f Makefile
