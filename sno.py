'''
Functions to adjust *.sno file
Author: @Hannes Hohenwarter
'''


import numpy as np


def read_sno(sno_file, mode):
    ''' Read *.sno file
    input: path to *.sno file
    output:
        - head: list with strings of header
        - timestamp: list (str)
        - data: np.array(float) of values in file
    '''

    #  line where data starts depending on mode (Richards/ Bucket)
    split = 26  # row split for Bucket
    if mode == 'Richards':
        split = 31  # row split for Richards (soil Layers)

    # open data and read header
    f = open(sno_file, 'r')
    lines = f.readlines()
    f.close()
    head = lines[0:split]                                  # get header
    # names = lines[24].strip().split()[2:]             # get var names

    # read data
    data = lines[split].strip().split()
    for line in lines[split+1:]:
        data = np.vstack((data, line.strip().split()))
    timestamp = data[:, 0]  # get timestamp
    data = data[:, 1:].astype(float)  # save data values as floating point numbers
    data = np.flipud(data)  # make top layer of Snowpack top layer of file

    return head, data, timestamp


def write_sno(outpath, head, data, timestamp):
    ''' save *.sno file
    :param outpath: Path to safe file
    :param head: list with strings of header
    :param data: np.array(float) of values in file
    :param timestamp: list (str)
    :return: save *.sno data to file
    '''

    data = np.flipud(data)  # make top Snowpack Layer bottom layer again for saving
    data = np.hstack((timestamp.reshape(len(timestamp), 1), data.round(5)))  # combine values with timestamp

    f = open(outpath, 'w')     # 'wt' write text
    f.writelines(head)  # write header
    np.savetxt(f, data, fmt='%s')  # write data
    f.close()


def adjustLayerThickness(head, data, thick, para, *args):
    ''' Change thickness of layer
    :param head: *.sno header list (str)
    :param data: np.array(float) of values in *.sno file
    :param thick: new layer thickness
    :param para: parameter to change (hWL or dWL)
    :param args: Index of layer to change
    :return: - head: header of *.sno file list(str)
             - data: adjusted data of *.sno file np.array(float)
    '''

    if len(args) == 0:  # if no layer defined -> change WL
        layer = findWL(data)
    elif len(args) == 1:
        layer = int(args[0])
    else:
        print('too many input Variables (data, thick, *Layer_index)')

    dz = data[layer, 0] - thick
    nelem = np.sum(data[:, 15])  # get current number of elements of layer to change

    data[layer, 0] = thick     # layer thickness = data[:,0]

    # adjust bottom layer to maintain overall height with dWL except if DepthHoar
    if get_SPname(head) != 'DepthHoar' and para == 'dWL':
        data[-1, 0] += dz
        data = adjustNe(data, -1)

    # calculate overall height and change in header
    height = np.sum(data[:, 0])
    head[11] = head[11].replace(head[11][-9:-1], str('%.6f' % height))

    # adjust number of elements according to layer thickness
    data = adjustNe(data, layer)

    # correct number of elements in bottom layer if neccesary to maintain overall number of elements
    if nelem != np.sum(data[:, 15]):
        data[-1, 15] = int(nelem - np.sum(data[0:-1, 15]))

    return head, data


def adjustDensity(data, p, *args):
    ''' Change density of layer
    :param data: np.array(float) of values in *.sno file
    :param p: density values (list like)
    :param args: Index of layer to change
    :return: adjusted data of *.sno file (np.array(float))
    '''

    if len(args) == 0:  # if no layer defined -> change WL
        layer = findWL(data)
    elif len(args) == 1:
        layer = int(args[0])
    else:
        print('too many input Variables (data, density, *Layer_index)')

    data[layer, 2] = np.round(data[layer, 2] * p, 4)   # multiply Vol_Frac_I with percentage (0.8-1.2)
    data[layer, 4] = 1-data[layer, 2]       # calculate Vol_Frac_V

    return data


def adjustCornSize(data, d, *args):
    ''' Change corn size of layer
    :param data: np.array(float) of values in *.sno file
    :param d: corn size values (list like)
    :param args: Index of layer to change
    :return: adjusted data of *.sno file (np.array(float))
    '''

    if len(args) == 0:  # if no layer defined -> change WL
        layer = findWL(data)
    elif len(args) == 1:
        layer = int(args[0])
    else:
        print('too many input Variables (data, d, *Layer_index)')

    rg = d/2    # calculate Grain Radius from diameter
    cornform = data[layer, 13]  # get grain type of layer

    data[layer, 9] = rg
    data[layer, 10] = get_rb(rg, cornform)  # change rb as well

    return data


def adjustNe(data, layer):
    ''' Adjust number of elements of a layer
    :param data: np.array(float) of values in *.sno file
    :param layer: layer index to change
    :return: adjusted data of *.sno file (np.array(float))
    '''

    z = data[layer, 0] * 100  # layer thickness in [cm]
    # use 4 Layers per cm
    data[layer, 15] = np.round(z * 4, 0)         # Ne = data[:,15]

    return data


def findWL(data):
    '''
    :param data: np.array(float) of values in *.sno file
    :return: Boolean Array if weak layer
    '''

    # define WL by the grain type
    ind1 = data[:, 13] == 0  # FC
    ind2 = data[:, 13] == 1  # DH
    ind3 = data[:, 13] == 3  # SH

    return ind1 + ind2 + ind3


def get_rb(rg, cornform):
    '''
    :param rg: grain radius
    :param cornform: grain type
    :return: rb: bond radius
    '''
    # marker according to snowpack metamorphism.cc
    # PP=DF=0, FC=DH=1, RG = 2, SH = 3
    rb_rg = {'0': 0.25, '1': 0.45, '2': 0.45, '3': 0.4}   # rb/rg relation from Walter

    rb = rb_rg[str(int(cornform))] * rg

    return rb


def get_SPname(head):
    '''
    :param head: list with strings of header
    :return: name of Snowprofile
    '''
    begin = head[2].find('=') + 2
    end = head[2].find('\n')
    SPname = head[2][begin: end]

    return SPname


