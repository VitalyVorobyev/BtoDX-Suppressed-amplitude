#! /usr/bin/python3
""" deld_B phase scan """

import sys
import re
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pylab as plt

REDELD = re.compile(r'db [0-9]+')
BETAFITRE = re.compile(r'beta = [0-9]+.[0-9]+ \+- [0-9]+.[0-9]+')
FITRESRE1 = re.compile(r'beta = -?\d*\.\d* \+- \d*\.\d*')
NEWDELB = re.compile(r'New delb \d*')

DATA_PATH = '/home/vitaly/B0toD0pipi/btoucbard/data/'
LOGS_PATH = '/home/vitaly/B0toD0pipi/btoucbard/logs/'

def fitfcn(delb, amp, phase, offset):
    """ Bias fit function """
    return amp * np.cos(delb/180.*np.pi - phase) + offset

def fitfcn2(delb, amp_cos, amp_sin, offset):
    """ Bias fit function """
    return amp_cos * np.cos(delb/180.*np.pi) +\
           amp_sin * np.sin(delb/180.*np.pi) + offset

def get_data(fname):
    """ Read fit results """
    delb, betav, betae = [], [], []
    with open(fname) as infile:
        for line in infile:
            if REDELD.match(line):
                delb.append(int(line[2:]))
            elif BETAFITRE.match(line):
                _, _, val, _, err = line.split()
                betav.append(float(val))
                betae.append(float(err))
    return [np.array(delb), np.array(betav), np.array(betae)]

def get_bias(data):
    data = np.array(data)
    return data[1::2] - data[::2]

def get_data2(fname):
    """ Read fit results """
    delb, betav, betae = [], [], []
    with open(fname) as infile:
        for line in infile:
            if NEWDELB.match(line):
                delb.append(int(line.split()[2]))
            elif FITRESRE1.match(line):
                _, _, val, _, err = line.split()
                betav.append(float(val))
                betae.append(float(err))
    return [np.array(delb), get_bias(betav), get_bias(betae)]

def db_scan_plot(posi_data, nega_data=None, prep=False):
    """ Show fit bias """
    x = np.arange(0., 360., 1.)

    delb, betav, betae = posi_data
    if prep:
        betav = preprocess_beta(betav)
    plt.errorbar(delb, betav, yerr=betae, fmt='bo')
    amp, phase, offset = make_fit(delb, betav)
    plt.plot(x, fitfcn(x, amp, phase, offset), 'k-')

    if nega_data is not None:
        delb, betav, betae = nega_data
        if prep:
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

def cp_scan():
    posi_data = get_data(DATA_PATH + 'cpp_delb_scan.txt')
    nega_data = get_data(DATA_PATH + 'cpn_delb_scan.txt')
    db_scan_plot(posi_data, nega_data)

def kspp_scan():
    # data = get_data2(LOGS_PATH + 'fit_dh_wf_rb_02_nocp_delb_scan.txt')
    datap = get_data2(LOGS_PATH + 'fit_dh_wf_rb_02_cpp_delb_scan.txt')
    datan = get_data2(LOGS_PATH + 'fit_dh_wf_rb_02_cpn_delb_scan.txt')
    db_scan_plot(datap, datan)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1] == 'cp':
            cp_scan()
        else:
            kspp_scan()
    else:
        cp_scan()
