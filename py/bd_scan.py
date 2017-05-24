#! /usr/bin/python3
""" deld_B phase scan """

import re
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pylab as plt

redeld = re.compile(r'db [0-9]+')
betafitre = re.compile(r'beta = [0-9]+.[0-9]+ \+- [0-9]+.[0-9]+')

DATA_PATH  = '/home/vitaly/B0toD0pipi/btoucbard/data/'

def fitfcn(delb, amp, phase, offset):
    """ Bias fit function """
    return amp * np.cos(delb/180.*np.pi - phase) + offset

def fitfcn2(delb, amp_cos, amp_sin, offset):
    """ Bias fit function """
    return amp_cos * np.cos(delb/180.*np.pi) +\
           amp_sin * np.cos(delb/180.*np.pi) + offset

def get_data(fname):
    """ Read fit results """
    delb, betav, betae = [], [], []
    with open(fname) as infile:
        for line in infile:
            if redeld.match(line):
                delb.append(int(line[2:]))
            elif betafitre.match(line):
                _, _, val, _, err = line.split()
                betav.append(float(val))
                betae.append(float(err))
    return [np.array(delb), np.array(betav), np.array(betae)]

def db_scan_plot(posi_data, nega_data):
    """ Show fit bias """
    x = np.arange(0., 360., 1.)

    delb, betav, betae = posi_data
    betav = preprocess_beta(betav)
    plt.errorbar(delb, betav, yerr=betae, fmt='bo')
    amp, phase, offset = make_fit(delb, betav)
    plt.plot(x, fitfcn(x, amp, phase, offset), 'k-')

    delb, betav, betae = nega_data
    betav = preprocess_beta(betav)
    plt.errorbar(delb, betav, yerr=betae, fmt='ro')
    amp, phase, offset = make_fit(delb, betav)
    plt.plot(x, fitfcn(x, amp, phase, offset), 'k-')

    plt.xticks(np.arange(0., 361., 45))

    plt.grid()
    plt.tight_layout()
    plt.show()

def preprocess_beta(betav, trueval=23.):
    """ Eliminate umbiguity and buntract true value """
    betav[betav > 90] -= 90
    betav[betav > 50] = 90 - betav[betav > 50]
    betav -= trueval
    return betav

def make_fit(delb, betav):
    """ Run fit """
    popt, pcov = curve_fit(fitfcn, delb, betav, [2.5, 0., 0.])
    for idx, par in enumerate(popt):
        print(par, '+-', np.sqrt(pcov[idx, idx]))
    return popt

POSI_DATA = get_data(DATA_PATH + 'db_scan_posicp.txt')
NEGA_DATA = get_data(DATA_PATH + 'db_scan_negacp.txt')
db_scan_plot(POSI_DATA, NEGA_DATA)
