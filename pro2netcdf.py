#!/usr/bin/env python
'''
transform all *.pro from one experiment into a more nice *.netcdf file

Folder structure:
    .../Snowpack-3.5.0/
        data/
            output/
                ArithMean/
                    DepthHoar_R/
                        WETMOD2_dWL1.pro
                        WETMOD2_dWL2.pro
                        ...
                        WETMOD2_hWL1.pro
                        ...
                    MeltCrust_R/
                        ...
                    SurfHoar_R/
                        ...
                    LayerCH_R/
                GeoMean/
                    ...
Author: @Hannes Hohenwarter
'''

import sys
import getopt
import xarray as xr
import datatable as dt
import numpy as np
import os
import fnmatch
import re
import pandas as pd


def toNetcdf(filepath, flist):
    '''
    :param filepath: path to *.pro files
    :param flist: list of *.pro files with the same type
    :return: ds (xarray of data)
    '''

    # sort filenames in flist according do run number e.g. dWL1, dWL2,... etc.
    regex = re.compile(r'\d+')  # set re.compiler to find \d+ matches 1 or more unicode decimal digit
    run = list()
    for f in flist:  # get run number out of file name
        run.append(int(regex.findall(f)[1]))  # two numbers in filename, second important
    flist = [x for _, x in sorted(zip(run, flist))]  # sort flist by number of run
    run.sort()  # sort list of runs

    # variable dictionary of key number in *.pro
    var_dict = {'height': 501, 'density': 502, 'temp': 503, 'eliD': 504, 'lwc': 506, 'bond_size': 511,
                'grain_size': 512, 'grain_type': 513,  'stress': 517, 'shear_strength': 601, 'Sn38': 532,
                'Sk38': 533, 'dz': None}
    variables = [None] * len(var_dict)
    h_variables = [None] * len(var_dict)

    i = 0
    for file in flist:
        # load *.pro as datatable; skip Header of file!; fill missing values (NaN)
        datatable = dt.fread(filepath + file, skip_to_line=48, fill=True)
        # get number of Snowpack elements (whithout ground)
        index = np.max(np.array(datatable[dt.f.C0 == 508, 1].to_list()).astype(int))

        # ---- extract data from datatable

        # set up coordinates and arrays at first loop
        if i == 0:
            #  create coordinates (dates and nElems) only at first loop
            dates = datatable[dt.f.C0 == 500, 1].to_list()
            date = [pd.to_datetime(date, format='%d.%m.%Y %H:%M:%S') for date in dates[0]]

            j = 0
            # read data from first *.pro file and save into np.array
            for value in var_dict.values():
                ind = get_ind(value, index)  # index depending on variable
                variables[j] = np.array(datatable[dt.f.C0 == value, ind[0]:ind[1]])
                j += 1
            variables[-1] = np.c_[variables[0][:, 0], np.diff(variables[0])]

        else:
            j = 0
            # read data from all other *.pro files and save into np.array
            for value in var_dict.values():
                ind = get_ind(value, index)
                h_variables[j] = np.array(datatable[dt.f.C0 == value, ind[0]:ind[1]])
                j += 1
            h_variables[-1] = np.c_[h_variables[0][:, 0], np.diff(h_variables[0])]

            # compare dimensions of arrays: (different size due to different Ne in *.pro)
            dim_new = h_variables[0].shape
            dim_old = variables[0].shape
            dim_diff = dim_new[1] - dim_old[1]

            # adjust array dimensions to put both together as one multidim array, fill empty spots with NaN
            if dim_diff > 0:  # add NaN's to previous matrix
                if len(dim_new) == len(dim_old):
                    add = np.empty((dim_new[0], np.abs(dim_diff)))
                else:
                    add = np.empty((dim_new[0], np.abs(dim_diff), dim_old[2]))  # matrix to add
                add.fill(np.nan)

                for j in np.arange(len(variables)):
                    variables[j] = np.append(variables[j], add, axis=1)

            elif dim_diff < 0:  # add NaN's to new matrix
                add = np.empty((dim_new[0], np.abs(dim_diff)))
                add.fill(np.nan)

                for j in np.arange(len(h_variables)):
                    h_variables[j] = np.append(h_variables[j], add, axis=1)

            else:
                pass

            # stack data
            for j in np.arange(len(variables)):
                variables[j] = np.dstack((variables[j], h_variables[j]))

        i += 1  # count loop index

    # create Dataset (xArray)
    dim = variables[0].shape
    elem = np.arange(dim[1])  # number of elements = second axis of variables

    ds = xr.Dataset()
    ds.coords['date'] = date
    ds.coords['elem'] = elem
    ds.coords['run'] = run

    # Units
    unit_dict = {'height': '[cm]', 'dz': '[cm]', 'density': '[kg m-3]', 'temp': '[degC]', 'lwc': '[%]',
                 'bond_size': '[mm]', 'grain_size': '[mm]', 'grain_type': 'Swiss Code F1F2F3',
                 'stress': '[kPa]', 'shear_strength': '[kPa]', 'Sn38': '', 'Sk38': ''}
    ds.attrs = unit_dict

    j = 0
    # write multidim variable array into xArray
    for key in var_dict.keys():
        ds[key] = (('date', 'elem', 'run'), variables[j])
        j += 1

    # ds = ds.transpose('elem', 'date', 'run')  # transpose ds for date beeing x axis or not -\o/-

    # ---- save xarray as netcdf
    return ds


def get_ind(value, ind):
    '''https://drive.google.com/drive/folders/1WIj9hyrWnedMdjRfPQARzW6obvJjRv7x?usp=sharing
    get line index where snow data starts and ends to neglect ground data
    '''

    ground = [502, 503, 504, 506, 515, 516, 517, 518, 519, 520, 521]  # variables containing ground info

    if value == 501:
        inds = int(-ind)
        inde = None
    elif value in ground:
        inds = int(-(ind+1))
        inde = int(-1)
    else:
        inds = int(2)
        inde = int(2+ind)

    return [inds, inde]


def main(argv):

    # default: run all Variations of Profiles and Parameters if none specific is chosen
    # via command line input
    profiles = ['MeltCrust_R', 'SurfHoar_R', 'DepthHoar_R', 'LayerCh_R']
    types = ['kWL', 'hWL', 'dWL', 'roh']
    mean_typ = 'ArithMean'
    profile_in = False
    mode = 'RE'

    # get input from command line. default: above if none chosen
    try:
        opts, args = getopt.getopt(argv, "hs:m:", ["sfile=" "mfile="])
    except getopt.GetoptError:
        print('pro2netcdf.py -s <Name_of_snowprofile> -m <Mode ("RE" or "B">')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('pro2netcdf.py -s <Name_of_snowprofile> -m <Mode ("RE" or "B">')
            sys.exit()
        elif opt in ("-s", "--sfile"):
            profiles = [arg]
            profile_in = True
        elif opt in ("-m", "--mfile"):
            mode = arg

    if not profile_in and mode == 'B':  # use Bucket Files if no Profile as input. else use input
        profiles = ['MeltCrust', 'SurfHoar', 'DepthHoar', 'LayerCh']

    # define Paths: folders of *.pro files and output folder
    path = '/Users/fuzi/Documents/Simulations/Snowpack-3.5.0/data/output/'
    outpath = '/Users/fuzi/Documents/Simulations/Snowpack-3.5.0/data/output/NetCDF/'

    if not os.path.exists(outpath):  # create output folder
        os.makedirs(outpath)

    for profile in profiles:
        filepath = path + mean_typ + '/' + profile + '/'  # path for *.pro files
        if not os.path.exists(filepath):  # check if path exists
            print(filepath + ' file does not exist')
            continue
        else:
            # get all files in 'path' including 'kWL' and '.pro'
            for type in types:
                # list of *.pro files in folder depending on type ('hWL',...)
                flist = fnmatch.filter(os.listdir(filepath), 'WETMOD2*' + type + '*.pro')

                if not flist:
                    continue
                else:
                    # safe all *.pro files whith same type into an *.nc file
                    ds = toNetcdf(filepath, flist)
                    outfile = outpath + profile + '_' + type + '.nc'
                    ds.to_netcdf(path=outfile)  # save as NETCDF4 file


if __name__ == "__main__":
    main(sys.argv[1:])
