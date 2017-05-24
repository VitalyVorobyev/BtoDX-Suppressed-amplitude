#! /usr/bin/python3
""" dt plot """

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pylab as plt

DATA_PATH = '/home/vitaly/B0toD0pipi/btoucbard/data/'

def fitfcn(x, dm, c, s):
    """ Time-dependent A_CP """
    return c * np.cos(x * dm) + s * np.sin(x * dm)

def make_fit(dt, acp):
    """ Run fit """
    popt, pcov = curve_fit(fitfcn, dt, acp, [0.5, 0., 0.7])
    for idx, par in enumerate(popt):
        print(par, '+-', np.sqrt(pcov[idx, idx]))
    return popt

def read_data(infile_name):
    """ Read dt and tag from text file """
    data = []
    with open(infile_name) as infile:
        for line in infile:
            time, tag = line.split()
            data.append([float(time), int(tag)])
    return np.array(data, dtype=((np.float, np.float)))

def dtplot(data):
    """ dt histograms for two tags """
    posi_data, bins = np.histogram(data[data[:, 1] == 1][:, 0], bins=40)
    nega_data, bins = np.histogram(data[data[:, 1] == -1][:, 0], bins=40)
    xvals = 0.5 * (bins[:-1] + bins[1:])
    plt.figure(num=1, figsize=(6, 4))
    plt.plot(xvals, posi_data, 'bo')
    plt.plot(xvals, nega_data, 'rp')
    plt.tight_layout(pad=1.08)
    plt.grid()

    plt.figure(num=2, figsize=(6, 4))
    diff_data = (posi_data - nega_data) / (posi_data + nega_data)
    dm, c, s = make_fit(xvals, diff_data)
    plt.plot(xvals, diff_data, 'ko')
    xfcn = np.linspace(xvals[0], xvals[-1]+0.01, 200)
    plt.plot(xfcn, fitfcn(xfcn, dm, c, s), 'b-')
    plt.tight_layout(pad=1.08)
    plt.grid()

    plt.show()

DATA = read_data(DATA_PATH + 'data.txt')
dtplot(DATA)
