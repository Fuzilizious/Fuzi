#!/usr/bin/env python

'''
Main Programm for Master Thesis
no beautiful coding and pretty much hard coded for my purpose ;)

Author: @Hannes Hohenwarter

Folder Structure:
.../Snowpack-3.5.0/
        data/
            cfgfiles/
                Orig/
                    *.ini files (originals) -> backup
                *.ini files (working)
            input/
                Orig/
                    *.sno (originals) -> starting snowpack stratification
                *.sno (working) -> created at first run
                *.smet -> meteo data
            output/ -> outpath
                Folders are created for each experiment
'''

import sys
import getopt
import SnowTools as st
import sno
import os
import numpy as np


def main(argv):

    # default: run all Variations of Profiles and Parameters if none specific is chosen
    # via command line input

    # Name of the original *.sno files: '_R' for using Richards Equation
    snowprofile = ['MeltCrust_R', 'SurfHoar_R', 'DepthHoar_R', 'LayerCH_R']

    # parameter to change
    # kWL: grain Size of the WL
    # hWL: height of the slab
    # dWL: thickness of the WL
    # density of the WL
    paras = ['kWL', 'hWL', 'dWL', 'roh']

    # get input from command line. default: above if none chosen
    try:
        opts, args = getopt.getopt(argv, "hs:p:", ["sfile=", "pfile="])
    except getopt.GetoptError:
        print('MySnowpack.py -s <Name_of_snowprofile> -p <Parameter to change>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('MySnowpack.py -s <Name_of_snowprofile> -p <Parameter to change (roh, hWL, dWL, kWL)>')
            sys.exit()
        elif opt in ("-s", "--sfile"):
            snowprofile = [arg]
        elif opt in ("-p", "--pfile"):
            paras = [arg]

    # Paths for input data
    meteopath = "/Users/fuzi/Documents/Simulations/Snowpack-3.5.0/data/input/"
    station = "WETMOD2.smet"
    snowpath = "/Users/fuzi/Documents/Simulations/Snowpack-3.5.0/data/input/"

    for profile in snowprofile:  # loop over chosen profiles
        # Paths for *.ini and *.sno files
        inifile = "/Users/fuzi/Documents/Simulations/Snowpack-3.5.0/data/cfgfiles/" + profile + ".ini"
        snowfile = profile + ".sno"
        outpath = "/Users/fuzi/Documents/Simulations/Snowpack-3.5.0/data/output/" + profile + "/"

        mode = 'Bucket'
        if profile[-1] == 'R':  # check last letter of profile is 'R' for Richards Mode
            mode = 'Richards'   # select water transport mode by filname ending 'R'
            profile = profile[:-2]  # rename snowprofile to get parameters from lists in SnowTools.py

        for para in paras:  # loop over chosen parameters
            if not os.path.exists(outpath):  # create output directory
                os.makedirs(outpath)

            # adjust *.ini file
            st.config_ini(inifile, meteopath, station, snowpath, snowfile, outpath)

            # run Snowpack in loop
            index = 1
            parameter = getattr(st, para)  # get values of parameters

            if parameter[profile] is None:  # don't loop for empty parameter
                continue

            for p in parameter[profile]:  # loop over parameters to change
                # read *.sno file
                [head, data, timestamp] = sno.read_sno(snowpath + 'Orig/' + snowfile, mode)
                # change parameter
                if para == 'roh':
                    data = sno.adjustDensity(data, p)
                elif para == 'dWL':
                    [head, data] = sno.adjustLayerThickness(head, data, np.round(p * 1e-2, 3), para)
                elif para == 'hWL':
                    [head, data] = sno.adjustLayerThickness(head, data, np.round(p * 1e-2, 3), para, 0)  # 0 -> Slab
                elif para == 'kWL':
                    data = sno.adjustCornSize(data, p)
                # save adjusted *.sno file
                sno.write_sno(snowpath + snowfile, head, data, timestamp)

                st.out_extension(inifile, para, index)  # set filename extension according to loop run
                st.run_snowpack(inifile, st.get_enddate(meteopath + station))  # run SNOWPACK

                index += 1


if __name__ == "__main__":
    main(sys.argv[1:])
