'''
Tools for running MySnowpack.py
Author: @Hannes Hohenwarter
'''

import configparser
import numpy as np
import shutil


def config_ini(inifile, meteopath, station, snowpath, snowfile, outpath):
    '''
    function to adjust *.ini file for running SNOWPACK
    :param inifile: Path to *.ini file
    :param meteopath: Path to meteo file
    :param station: station name
    :param snowpath: Path to *.sno files
    :param snowfile: *.sno filename
    :param outpath: Path to save SNOWPACK output
    :return: adjust Paths ect. in *.ini file and safe *.ini
    '''

    # load *.ini via ConfigParser
    config = configparser.ConfigParser()
    config.read(inifile)

    # chanche path names in *.ini file
    config['INPUT']['METEOPATH'] = meteopath
    config['INPUT']['STATION1'] = station
    config['INPUT']['SNOWPATH'] = snowpath
    config['INPUT']['SNOWFILE1'] = snowfile
    config['OUTPUT']['METEOPATH'] = outpath

    # safe *.ini file
    write_ini(inifile, config)


def write_ini(path_ini, config):
    '''
    Safe *.ini file
    :param path_ini: Path to write *.ini file
    :param config: data from ConfigParser
    :return: save *.ini file
    '''
    with open(path_ini, 'w') as configfile:
        config.write(configfile)


def run_snowpack(inifile, enddate):
    '''
    :param inifile: Path to *.ini file
    :param enddate: End date for Simulation
    :return Running shell script of SNOWPACK line by line
    '''
    import os

    os.system('TOOL="valgrind --tool=callgrind --simulate-cache=yes"')
    os.system('TOOL="/software/bin/valgrind --tool=memcheck \
        --leak-check=full --show-reachable=yes --track-origins=yes --log-fd=2 "')
    os.system('TOOL="time"')
    os.system('TOOL=""')
    os.system('${TOOL} snowpack -c ' + inifile + ' -e ' + enddate)


def out_extension(inifile, para, index):
    '''
    :param inifile: Path to *.ini file
    :param para: name of Parameter
    :param index: Index of model run
    :return: change output extension of *.ini depending on Parameter and model run
    '''
    config = configparser.ConfigParser()
    config.read(inifile)
    config['OUTPUT']['experiment'] = para + str(index)

    write_ini(inifile, config)


def get_enddate(file):
    '''
    :param file: Path to *.smet file
    :return: enddate from Meteo data
    '''
    with open(file, 'r') as f:
        lines = f.read().splitlines()
        last_line = lines[-1]
        enddate = last_line.split("\t")[0]

    return enddate


def reset_sno(snowfile, snowpath):
    ''' reset *.sno file by copying original one into working folder '''
    orig = snowpath + "Orig/" + snowfile
    shutil.copy(orig, snowpath + snowfile)


def Sn_wet(ds, kind='auto'):
    ''' Calculate Sn_wet based on Yamanoi and Endo 2002
    :param ds:  xArray Dataset
    :param kind: =auto: get sigma & tau from ds;
                 =manu: calc sigma via Jamison & Johnson 2001
    :return: Sn_wet
    '''

    if kind == 'auto':
        Sn = ((-1) * ds.shear_strength * np.exp(-0.235 * ds.lwc)) / ds.stress

    elif kind == 'manu':
        Sn = ((-1) * shear_str(ds) * np.exp(-0.235 * ds.lwc)) / ds.stress

    return Sn


def Sn_dry(ds, kind='auto'):
    ''' Calculate Sn_dry
       :param ds:  xArray Dataset
       :param kind: =auto: get sigma & tau from ds;
                    =manu: calc sigma via Jamison & Johnson 2001
       :return: Sn_dry
    '''
    if kind == 'auto':
        Sn = (-1) * ds.shear_strength / ds.stress
    elif kind == 'manu':
        Sn = (-1) * shear_str(ds) / ds.stress

    return Sn


def Sn38(ds, kind='auto'):
    ''' Calculate Sn38 by hand
    :param ds: xArray Dataset
    :param kind: =auto: get sigma & tau from ds;
                 =manu: calc sigma via Jamison & Johnson 2001
    :return: Sn38
    '''

    if kind == 'auto':
        Sn38 = (-1) * ds.shear_strength / ds.stress * np.cos(38)
    elif kind == 'manu':
        Sn38 = (-1) * shear_str(ds) / ds.stress * np.cos(38)

    return Sn38


def shear_str(ds):
    ''' Calculate shear strength after Jamison & Johnson 2001
    :param ds: xArray Dataset
    :return: shear_strength
    '''

    sigma = 22 * (ds.density / 917) ** 2.2  # equ. from Jamison and Johnson 2001

    return sigma


def orig_run(profile, typ):
    ''' Get Model run equal to original snowpack setting.
    :param profile: name of snowprofile
    :param typ: changed parameter
    :return: orignal model run (int)
    '''

    # Original model run for different types and profiles
    orig_dWL = {'DepthHoar_R': 3, 'SurfHoar_R': 2, 'MeltCrust_R': 2, 'LayerCh_R': None,
                'MeltCrust_hr_R': 2, 'SurfHoar_hr_R': 2}
    orig_kWL = {'DepthHoar_R': 6, 'SurfHoar_R': 5, 'MeltCrust_R': 3, 'LayerCh_R': 3}

    if typ == 'hWL':
        return 6
    elif typ == 'rho':
        return 5
    elif typ == 'dWL':
        return orig_dWL[profile]
    elif typ == 'kWL':
        return orig_kWL[profile]


def get_ext_time(ds, var, extr):
    '''
    get timespan when wet snow indices exceed critical values
    :param ds: xarray dataset loaded from netcdf
    :param var: lwc, Sn, ect.
    :param extr: max, min, mean
    :return: xarray dataset dim= date, run
    '''

    import xarray as xr

    if var == 'lwc':
        if extr == 'max':
            ds_out = ds[var].max(dim='elem')
            ds_out.values = xr.where(ds_out.values >= 6, 1, np.nan)
            return ds_out
        elif extr == 'mean':
            ds_out = ds[var].mean(dim='elem')
            ds_out.values = xr.where(ds_out.values >= 3, 1, np.nan)
            return ds_out
    elif var == 'Sn_wet':
        ds[var] = Sn_wet(ds)
        ds.attrs[var] = ''
        if extr == 'min':
            ds_out = ds[var].min(dim='elem')
            ds_out.values = xr.where(ds_out.values <= 1, 1, np.nan)
            return ds_out


# Parameter Ranges for sensitivity analysis:
# Density ranges ± 20%; 5% steps
rho = {'DepthHoar': np.linspace(0.8, 1.2, 9), 'SurfHoar': np.linspace(0.8, 1.2, 9),
       'MeltCrust': np.linspace(0.8, 1.2, 9), 'LayerCh': np.linspace(0.8, 1.2, 9)}

# height of Slab 5-30cm
hWL = {'DepthHoar': np.arange(5, 31, 5), 'SurfHoar': np.arange(5, 31, 5),
       'MeltCrust': np.arange(5, 31, 5), 'LayerCh': np.arange(5, 31, 5)}  # shifted * 1e-2 into main program

# height of WL: max corn size - 3cm; increments depending on corn size
# MeltCrust: height of WL > 3mm to be distinct layer
dWL = {'DepthHoar': np.arange(0.4, 3.3, 0.4), 'SurfHoar': np.arange(0.25, 3.1, 0.25),
       'MeltCrust': np.arange(0.3, 3.1, 0.3), 'LayerCh': None, 'MeltCrust_hr': np.arange(0.3, 3.1, 0.3),
       'SurfHoar_hr': np.arange(0.25, 3.1, 0.25)}

# corn size (diameter) dependent on corn form 0.25 mm increments
kWL = {'DepthHoar': np.arange(1.5, 5.1, 0.5), 'SurfHoar': np.arange(1.5, 2.6, 0.25),
       'MeltCrust': np.arange(0.5, 1.3, 0.25), 'LayerCh': np.arange(0.5, 1.3, 0.25)}

# density relation WL/Slab for original setting
density_rel = {'DepthHoar_R': 1.03, 'SurfHoar_R': 0.38, 'MeltCrust_R': 0.80, 'LayerCh_R': 1.19}
