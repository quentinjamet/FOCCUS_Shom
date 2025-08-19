#!/usr/bin/env python

import numpy as np
import xarray as xr
import datetime
import os
import glob
from scipy.interpolate import griddata
import json


#########################
def read_params(filename):
    '''
    Decode input parameters from json file.
    '''
    with open(filename, "r") as jsonfile:
        params = json.load(jsonfile)
    var_info  = params["var_info"]
    date_info = params["date_info"]
    ens_info  = params["ens_info"]
    dir_info  = params["dir_info"]  

    return var_info, date_info, ens_info, dir_info

################################
def date_range_list(date_info):
    '''
    Return list of datetime.date objects (inclusive) between start_date and end_date (inclusive).
    '''
    
    start_date = datetime.date(date_info['yr_ini'], date_info['mth_ini'], date_info['day_ini'])
    end_date   = datetime.date(date_info['yr_end'], date_info['mth_end'], date_info['day_end'])
    step       = date_info['step']
    
    date_list = []
    curr_date = start_date
    while curr_date <= end_date:
        date_list.append(curr_date)
        curr_date += datetime.timedelta(days=step)
    return date_list


########################################
def make_regular_mesh(da, var_info):
    '''
    Generate a regulat lon/lon mesh and associated mask based on temperature field.
    '''

    print("-- Define regular lat/lon grid --")
    dlon = eval("(da.%s[:, 1:]-da.%s[:, :-1]).mean().data" % \
                (var_info['longitude'], var_info['longitude']) )
    dlat = eval("(da.%s[1:, :]-da.%s[:-1, :]).mean().data" % \
                (var_info['latitude'], var_info['latitude']) )
    lon = eval("np.arange(da.%s.min(), da.%s.max()+dlon, dlon)" % \
                (var_info['longitude'], var_info['longitude']) )
    lat = eval("np.arange(da.%s.min(), da.%s.max()+dlat, dlat)" % \
                (var_info['latitude'], var_info['latitude']) )
    depth = eval("da.%s.data" % (var_info['z']) )
    #
    [nz, ny, nx] = [da.dims[var_info['z']], len(lat), len(lon)]

    #-- construct associated land/ocean mask --
    print("-- Generate 3D land/ocean mask on regular grid --")
    msk = eval("xr.where(da.%s.isel(%s=0) > 0, 1., 0.)" % (var_info['thetao'], var_info['time'] ))
    msk = eval("msk.stack(yx=('%s', '%s'))" % (var_info['y'], var_info['x']) )
    msk_reg = np.zeros([nz, ny, nx])
    for kkk in range(nz):
        msk_reg[kkk, ...] = eval("griddata((msk.%s, msk.%s), msk[kkk, :], (lon[None, :], lat[:, None]), method='nearest')" %
                                 (var_info['longitude'], var_info['latitude']) )
    msk_reg = np.where(msk_reg > 0.999, 1., np.nan)

    return lat, lon, depth, msk_reg


#############################################################
def create_dataset(dir_info, var_info, ens_info, nt=1):
    '''
    Create a CMEMS-like dataset on a regular lat/lon grid.
    '''

    #-- temporary temperature dataarray used to construct mask --
    filelist = sorted( glob.glob("%s/*%s*" % (dir_info["dir_in"], var_info['thetao']) ) )
    tmp = xr.open_dataset(filelist[0])

    #-- extract attributrs history --
    attrs_his = tmp.attrs['history']

    #-- make regular lat/lon mesh and mask --
    [lat, lon, depth, mask] = make_regular_mesh(tmp, var_info)

    #-- get dimensions and mesh --
    [nz, ny, nx] = [len(depth), len(lat), len(lon)]
    [dlat, dlon]   = [lat[1]-lat[0], lon[0]-lon[1]]
   
    if ens_info["ens_in"] and ens_info["ens_out"]:
        nmem = len(ens_info["number"])

    if ens_info["ens_out"]:
        variables=dict(
                zos=(["number", "time", "latitude", "longitude"], np.full((nmem, nt, ny, nx), np.nan) ),
                uo=(["number", "time", "depth", "latitude", "longitude"], np.full((nmem, nt, nz, ny, nx), np.nan) ),
                vo=(["number", "time", "depth", "latitude", "longitude"], np.full((nmem, nt, nz, ny, nx), np.nan) ),
                thetao=(["number", "time", "depth", "latitude", "longitude"], np.full((nmem, nt, nz, ny, nx), np.nan) ),
                so=(["number", "time", "depth", "latitude", "longitude"], np.full((nmem, nt, nz, ny, nx), np.nan) ),
            )
        coordinates=dict(
             depth=depth.astype('float'),
             latitude=lat.astype('float'),
             longitude=lon.astype('float'),
             time=np.zeros(nt).astype('float'),
             number=ens_info["number"].astype('int16'),
         )
    else:
        variables=dict(
                zos=(["time", "latitude", "longitude"], np.full((nt, ny, nx), np.nan) ),
                uo=(["time", "depth", "latitude", "longitude"], np.full((nt, nz, ny, nx), np.nan) ),
                vo=(["time", "depth", "latitude", "longitude"], np.full((nt, nz, ny, nx), np.nan) ),
                thetao=(["time", "depth", "latitude", "longitude"], np.full((nt, nz, ny, nx), np.nan) ),
                so=(["time", "depth", "latitude", "longitude"], np.full((nt, nz, ny, nx), np.nan) ),
            )
        coordinates=dict(
             depth=depth.astype('float'),
             latitude=lat.astype('float'),
             longitude=lon.astype('float'),
             time=np.zeros(nt).astype('float'),
         )

    #-- create dataset --
    ds = xr.Dataset(
        coords=coordinates,
        data_vars=variables,
        attrs=dict(title='Hourly mean fields from Ensemble Global Ocean Physics Analysis ',
                   comments='Retrieve from Mercator Ocean International ftp',
                   history=[attrs_his, \
                            "%s -- convertion to CMEMS product type by Quentin Jamet (quentin.jamet@shom.fr)" % \
                            datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S') ],
                   easting='longitude',
                   field_type='mean',
                   northing='latitude',
                   Conventions='CF-1.6',
               )
    )
   
    #-- add land/ocean mask --
    ds["mask"] = (["depth", "latitude", "longitude"], mask)
    
    #-- variables attributes --
    if ens_info["ens_out"]:
        ds.number.attrs      = {'long_name': 'ensemble member numerical id',
                            'units': '1',
                            'standard_name': 'realization'}
    ds.depth.attrs       = {'axis': 'Z',
                         'long_name': 'Depth',
                         'positive': 'down',
                         'standard_name': 'depth',
                         'unit_long': 'Meters',
                         'units': 'm'}
    ds.latitude.attrs    = {'axis': 'Y',
                        'long_name': 'Latitude',
                        'standard_name': 'latitude',
                        'unit_long': 'Degrees North',
                        'units': 'degrees_north',
                        'step': dlat,
                        'valid_max': lat.max(),
                        'valid_min': lat.min()}
    ds.longitude.attrs   = {'axis': 'X',
                        'long_name': 'Longitude',
                        'standard_name': 'longitude',
                        'unit_long': 'Degrees East',
                        'units': 'degrees_east',
                        'step': dlon,
                        'valid_max': lon.max(),
                        'valid_min': lon.min()}
    ds.zos.attrs         = {'cell_methods': 'area: mean',
                       'long_name': 'Sea surface height',
                       'standard_name': 'sea_surface_height_above_geoid',
                       'unit_long': 'Meters',
                       'units': 'm'}
    ds.uo.attrs          = {'cell_methods': 'area: mean',
                        'long_name': 'Eastward velocity',
                        'standard_name': 'eastward_sea_water_velocity',
                        'unit_long': 'Meters per second',
                        'units': 'm s-1'}
    ds.vo.attrs          = {'cell_methods': 'area: mean',
                        'long_name': 'Northward velocity',
                        'standard_name': 'northward_sea_water_velocity',
                        'unit_long': 'Meters per second',
                        'units': 'm s-1'}
    ds.thetao.attrs      = {'cell_methods': 'area: mean',
                        'long_name': 'Temperature',
                        'standard_name': 'sea_water_potential_temperature',
                        'unit_long': 'Degrees Celsius',
                        'units': 'degrees_C'}
    ds.so.attrs          = {'cell_methods': 'area: mean',
                        'long_name': 'Salinity',
                        'standard_name': 'sea_water_salinity',
                        'unit_long': 'Practical Salinity Unit',
                        'units': '1e-3'}
    
    
    #-- The end! --
    return ds



################################
def convert_to_cmems_type(filename):
    '''
    Reorganize and interpolate input dataset to a CMEMS-like dataset, i.e. all variables are interpolated on a regular, non-staggered latitude/longitude grid with conventional variable names.
    This script assumes an input dataset with one file per variables (i.e. ssh, temperature, salinity, zonal and meridional velocity), per model output dump (i.e. no time dimension in the original dataset), and in the case of ensemble simulations one file per ensemble member.
    Input file names should be in form 'XXXmyname*VARIABLE*_YYYYMMDD-YYYYMMDD*.nc', with XXX the ensemble member number, VARIABLE the name of the variable (see var_info), YYYYMMDD the date.
    In case of ensemble, the resulting CMEMS-like output dataset can either be saved with one netcdf file per ensemble member (i.e. ens_info["ens_out"] == False, default), 
    or one netcdf file with all ensemble members (i.e. ens_info["ens_out"] == True) with 'number' as the ensemble dimension.
    
    Input: 
    - filename: J[A]SON file containing the following parameters: 
        - dir_info: Input, Output and scratch directories informations, as well as output filename.
		To be provided in a dictionary format: dir_info = dict({"dir_in": '/my/directory/with/original/dataset/',
                                                                        "dir_out": 'my/directory/where/to/store/data/',
                                                                        "fileN": filename.nc).
		Assumes all data are directly accessible from "dir_in" by their file name, no sub-directories.
        - date_info: information on initial and final dates (i.e. [year, month, day]) to consider, and the time stepping increment (in days).
                To be provided in a dictionary format: date_info = dict({"yr_ini": YYYY, "mth_ini": MM, "day_ini": DD, "yr_end": YYYY, "mth_end": MM, "day_end": DD, "step": X})
        - ens_info: information on the type of input and output dataset (ensemble or single run), how to handle it in the output netcdf file, and the number of ensemble members in the input dataset, if necessary.
		To be provided in a dictionary format: ens_info = dict({"ens_in": True/False, "ens_out": True/False, "number": np.arange(nmem)}), with 'nmem' the number of ensemble members.
		If ens_info["ens_in"] ==  True and ens_info["ens_out"] == False, one netcdf file per ensemble member will be written out.
		If end_info["ens_in"] ==  True and ens_info["ens_out"] == True , ensemble memebers are agredated within the same netcdf file.
                'number' refers to the list of ensemble members IN THE INPUT DATASET to consider (can either be continuous, e.g. number=np.arange(10), or discontinuous, e.g. number=[2, 5, 10])
        - var_info: List of variable names associated with temperature ('thetao'), salinity ('so'), 
                    zonal ('uo') and meridional ('vo') velocity, sea surface elevation ('zos'), 
                    time, longitude, latitude, x and y grid index.
                Should be in the format: 
                var_info = dict({"thetao": "",
                                "so": "", 
                                "uo": "", 
                                "vo": "", 
                                "zos": "", 
                                "time": "", 
                                "longitude": "", 
                                "latitude": "",
                                "x": "", 
                                "y": "",
                                "z": "",
                                "x_u": "",
                                "y_u": "",
                                "z_u": "",
                                "x_v": "",
                                "y_v": "",
                                "z_v": "",})
    '''
 
    #-- Read parameters from json file --
    [var_info, date_info, ens_info, dir_info] = read_params(filename)   

 
    #-- check the type of data and associated list of variables --
    if var_info is None:
        print("STOP -- Please provide the name of variables associated to the type of data to load.")
        print("STOP -- Should be in a dictionary format (see description).")
        return

    
    #-- check time period to consider and generate list of dates to consider --
    if date_info is None: 
        print("STOP -- Please provide dates in a dictionary.")
        return
    if (datetime.datetime(date_info['yr_end'], date_info['mth_end'], date_info['day_end']).toordinal() - datetime.datetime(date_info['yr_ini'], date_info['mth_ini'], date_info['day_ini']).toordinal()) < 0:
        print("-- End date (%04.i/%02.i/02.i) is prior to start date (%04.i/%02.i/%02.i) -- Please provide time increasing dates." % (date_info['yr_end'], date_info['mth_end'], date_info['day_end'], date_info['yr_ini'], date_info['mth_ini'], date_info['day_in']))
        return
    #
    date_list = date_range_list(date_info)    


    #-- check output file name --
    if dir_info is None:
        print("STOP -- Please provide necessary informations on directories and output netcdf file (see description for details).")
        return
   
    #-- check and extract ensemble informations --
    if ens_info["ens_in"]:
        if ens_info["number"] is None:
            print("STOP -- Please provide ensemble members numbers")
            return
        else:
            nmem = len(ens_info["number"])
    else:
        ens_info["ens_out"] = False     # to insure consistancy between input and outputs 

    #-- initiate xarray dataset for output --
    ds = create_dataset(dir_info, var_info, ens_info, nt=1)

    
    #-- define CMEMS-type variable names --
    varList = ["zos", "thetao", "so", "uo", "vo"]

    #--------------
    # Loop on dates
    #--------------
    for idate in date_list:
        idate = idate.strftime("%Y%m%d")
        print("-- Convert date: %s" % idate )
        
        #-------------------------
        # Loop on ensemble members (break if not an ensemble)
        #-------------------------
        if ens_info["ens_out"]:
            #- cleanup dataset -
            for ivar in varList:
                ds[ivar] = eval("(ds.%s.dims, np.full(ds.%s.shape, np.nan) )" % (ivar, ivar) )
                           
        for imem in range(nmem):
            print('---- Ensemble member: %03.i' % (ens_info["number"][imem]) )
            if not ens_info["ens_out"]:
                #- cleanup dataset -
                for ivar in varList:
                    ds[ivar] = eval("(ds.%s.dims, np.full(ds.%s.shape, np.nan) )" % (ivar, ivar) )

            #------------------
            # Loop on variables 
            #------------------
            for ivar in varList:

                print('------ Interpolate variable: %s' % ivar)

                #- file name -
                if ens_info["ens_in"]:
                    tmpfile = str("%s/%03.i*%s*%s-%s*.nc" % (dir_info["dir_in"], ens_info["number"][imem], var_info[ivar], idate, idate) )
                else:
                    tmpfile = str("%s/*%s*%s-%s*.nc"      % (dir_info["dir_in"], var_info[ivar], idate, idate) )
    
                #------------------------
                # Loop on vertical levels (break for ssh)
                #------------------------
                for kkk in range(ds.dims['depth']):
                    
                    #- open file and select depth if not SSH -
                    if ivar == "zos":
                        tmpda = eval("xr.open_mfdataset('%s').%s" % (tmpfile, ivar) )
                    elif ivar == "uo":
                        tmpda = eval("xr.open_mfdataset('%s').%s.isel(%s=%i)" % (tmpfile, ivar, var_info['z_u'], kkk) )
                    elif ivar == "vo":      
                        tmpda = eval("xr.open_mfdataset('%s').%s.isel(%s=%i)" % (tmpfile, ivar, var_info['z_v'], kkk) )
                    else:
                        tmpda = eval("xr.open_mfdataset('%s').%s.isel(%s=%i)" % (tmpfile, ivar, var_info['z'], kkk) )

                    #- extract time dimension from SSH -
                    if ivar == "zos": 
                        ds["time"] = eval("tmpda.%s.data" % (var_info['time']) )
                        ds.time.attrs = dict(time_origin=eval("tmpda.%s.attrs['time_origin']" % var_info['time']))

                    #- prepare for interpolation (remove time dimension and stack on the horizontal) -
                    if ivar == "uo":
                        tmpda = eval("tmpda.isel(%s=0).stack(yx=('%s', '%s')).fillna(0.)" % \
                                     (var_info['time'], var_info['y_u'], var_info['x_u']) )
                    elif ivar == "vo":
                        tmpda = eval("tmpda.isel(%s=0).stack(yx=('%s', '%s')).fillna(0.)" % \
                                     (var_info['time'], var_info['y_v'], var_info['x_v']) )
                    else:
                        tmpda = eval("tmpda.isel(%s=0).stack(yx=('%s', '%s')).dropna(dim='yx', how='any')" % \
                                     (var_info['time'], var_info['y'], var_info['x']) )

                    #- break kkk loop if deeper than depest level (for tracers only) -
                    if tmpda.size == 0.:
                        break
                    
                    #- interpolate -
                    if ivar == 'uo' or ivar == 'vo':
                       tmpvar = eval("griddata((tmpda.%s, tmpda.%s), tmpda, (ds.longitude.data[None, :], ds.latitude.data[:, None]), method='linear')" % \
                                     (var_info['longitude'], var_info['latitude']))
                    else:
                       tmpvar = eval("griddata((tmpda.%s, tmpda.%s), tmpda, (ds.longitude.data[None, :], ds.latitude.data[:, None]), method='nearest')" % \
                                     (var_info['longitude'], var_info['latitude']))

                    #- apply land/ocean mask -
                    tmpvar *= ds.mask[kkk, ...]
 
                    #- write into CMEMS-like dataset -
                    if ens_info['ens_out']:
                        if ivar == "zos": 
                            ds[ivar][imem, 0,  ...]      = tmpvar
                        else:
                            ds[ivar][imem, 0, kkk,  ...] = tmpvar
                    else:
                        if ivar == "zos":
                            ds[ivar][0, ...]       = tmpvar 
                        else:
                            ds[ivar][0, kkk,  ...] = tmpvar 

                    #- break loop if SSH -
                    if ivar == "zos":
                        break

                    #-- END LOOP ON VERTICAL LEVELS --

                #-- END LOOP ON VARIABLES --

            #-- save one netcdf file per ensemble member --
            if not ens_info["ens_out"]:
                if ens_info["ens_in"]:
                     tmpout = str("%s/%03.i%s_%s-%s.nc" % (dir_info["dir_out"], ens_info["number"][imem], dir_info["fileN"], idate, idate) )
                else:
                     tmpout = str("%s/%s_%s-%s.nc"      % (dir_info["dir_out"], dir_info["fileN"], idate, idate) )
                #-- write to netcdf --
                print("(%s) -- write to netcdf: %s" % (datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'), tmpout) )
                ds.to_netcdf(path=tmpout, engine='netcdf4', compute=False)

            #-- break loop if not an ensemble --
            if not ens_info["ens_in"]:
                break

            #-- END LOOP ON ENSEMBLE MEMBERS --
           
        if ens_info["ens_out"]:
            tmpout = str("%s/%s_%s-%s.nc" % (dir_in, dir_info["fileN"], idate, idate) )
            #-- write to netcdf --
            print("(%s) -- write to netcdf: %s" % (datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'), tmpout) )
            ds.to_netcdf(path=tmpout, engine='netcdf4', compute=False)
 
    return
