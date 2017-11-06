#! /usr/bin/python3
""" Analyse results of a toy MC simulation """

import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

FITRESREX = re.compile(r'.*: -?\d{1}\.\d* -> -?\d{1}\.\d* \+- \d{1}\.\d*')
FITRESRE = re.compile(r'beta: -?\d*\.\d* -> -?\d*\.\d* \+- 0')
FITRESRE1 = re.compile(r'beta = -?\d*\.\d* \+- \d*\.\d*')
FITRESRE2 = re.compile(r'.* = -?\d*\.\d* \+- \d*\.\d*')
NEWDELB = re.compile(r'New delb \d*')

def fill(fname):
    """ Collect beta fit results """
    return [[float(line.split()[-3]), float(line.split()[-1])]\
            for line in open(fname).readlines()\
            if FITRESRE1.match(line)]

def full_fill(fname):
    """ Collect fit resuts of the full fit procedure """
    rsmpl, rcorr = {}, {}
    smplfl = True
    for line in open(fname).readlines():
        if 'Corrected' in line:
            smplfl = False
        if 'Simple' in line:
            smplfl = True
        if FITRESRE2.match(line):
            sline = line.split()
            varname, value, error = sline[0], float(sline[2]), float(sline[4])
            if smplfl:
                if varname in rsmpl:
                    rsmpl[varname].append([value, error])
                else:
                    rsmpl[varname] = [[value, error]]
            else:
                if varname in rcorr:
                    rcorr[varname].append([value, error])
                else:
                    rcorr[varname] = [[value, error]]
    return (rsmpl, rcorr)

def get_diff(var, rsmpl, rcorr, refval=None):
    """ Calculate average difference and find max difference """
    if (var not in rsmpl) or (var not in rcorr):
        print('get_diff: wrong variable {}'.format(var))
        return
    smpl = np.array(rsmpl[var])[:, 0]
    corr = np.array(rcorr[var])[:, 0]
    diff = corr - smpl
    if refval is not None:
        corr -= refval
        print('Corrected sample offset: {:.3f} +- {:.3f}'.format(\
            np.mean(corr), np.std(corr)/np.sqrt(len(corr)))\
        )
    return (diff, np.mean(diff), np.max(abs(diff)))

def plot_fitr(fitr):
    """ Plot fit results """
    smpl, full = fitr[::2], fitr[1::2]
    print(np.mean(smpl), '+-', np.std(smpl) / np.sqrt(len(smpl)))
    print(np.mean(full), '+-', np.std(full) / np.sqrt(len(full)))

    plt.figure(num=2)
    plt.hist(smpl, bins=40, color='blue', histtype='step')
    plt.hist(full, bins=40, color='red', histtype='step')
    plt.show()

def plot_bin_scan(fitr):
    """ Bin scan """
    smpl, full = fitr[::2, 0], fitr[1::2, 0]
    esmpl, efull = fitr[::2, 1], fitr[1::2, 1]
    bins = np.arange(len(smpl))
    plt.figure(num=1)
    plt.errorbar(bins-0.05, smpl, xerr=0, yerr=esmpl, ecolor='crimson', markerfacecolor='crimson',\
                 markeredgecolor='crimson', linestyle=' ', marker='o', label='simple')
    plt.errorbar(bins+0.05, full, xerr=0, yerr=efull, ecolor='steelblue', markerfacecolor='steelblue',\
                 markeredgecolor='steelblue', linestyle=' ', marker='o', label='corrected')
    plt.legend(loc='best')
    plt.grid()
    plt.tight_layout()

def plot_fit_diff(fitr):
    """ Fit difference for each Dalitz plot bin """
    x = fitr[::2, 0] - fitr[1::2, 0]
    ex = 0.5 * (fitr[::2, 1] + fitr[1::2, 1])
    bins = np.arange(len(x))
    plt.figure()
    plt.errorbar(bins, x, xerr=0, yerr=ex, ecolor='steelblue',\
                 markerfacecolor='steelblue',\
                 markeredgecolor='steelblue', linestyle=' ', marker='o')
    plt.ylim(ymax=np.max(1.5*abs(x)), ymin=-np.max(1.5*abs(x)))
    plt.grid()
    plt.tight_layout()
    plt.show()

# FITR = np.array(fill(sys.argv[1]))
# plot_bin_scan(FITR)
# plot_fit_diff(FITR)
# # print(len(FITR), "events")
# # plot_fitr(FITR)
# plt.show()

def read_dt_dist(fname):
    """ Read toy MC events """
    dt, tag, dbin, bbin = [[] for _ in range(4)]
    for line in open(fname).readlines():
        vals = line.split()
        dt.append(float(vals[0]))
        tag.append(int(vals[1]))
        dbin.append(int(vals[2]))
        bbin.append(int(vals[3]))
    return np.array([dt, tag, dbin, bbin])

def pdf(dt, sindbeta, tau, dm, tag):
    """ Signal PDF """
    return (1. + tag * sindbeta * np.sin(dm * dt)) * np.exp(-abs(dt) / tau) / (2. * tau)

def pdfp(dt, sindbeta, offset):
    """ Signal PDF """
    tau, dm = 1.520, 0.505
    ddt = dt - offset
    return (1. + sindbeta * np.sin(dm * ddt)) * np.exp(-abs(ddt) / tau) / (2. * tau)

def pdfn(dt, sindbeta, offset):
    """ Signal PDF """
    tau, dm = 1.520, 0.505
    ddt = dt - offset
    return (1. - sindbeta * np.sin(dm * ddt)) * np.exp(-abs(ddt) / tau) / (2. * tau)

def make_fit(dt, tag):
    """ Fit PDF to data """
    nbins = 200
    edge = [-10., 10.]
    raw_hist, _ = np.histogram(dt[tag==1], bins=nbins, range=edge, normed=False)
    raw_hist[raw_hist < 1] = 1
    hist, bins = np.histogram(dt[tag==1], bins=nbins, range=edge, normed=True)
    scale = raw_hist[100] / hist[100]
    yerr = np.sqrt(raw_hist) / scale
    x = 0.5 * (bins[:-1] + bins[1:])
    popt, pcov = curve_fit(pdfp, x, hist, [0.72, 0.])
    print(popt[0], '+-', np.sqrt(pcov[0, 0]))
    plt.figure(num=1)
    plt.errorbar(x, hist, xerr=0, yerr=yerr, fmt='k.')
    plt.plot(x, pdfp(x, *popt), color='orchid')

    plt.figure(num=2)
    pull = (hist - pdfp(x, *popt)) / yerr
    plt.plot(x, pull, 'k.')
    plt.grid()
    plt.tight_layout()

    raw_hist, _ = np.histogram(dt[tag==1], bins=nbins, range=edge, normed=False)
    raw_hist[raw_hist < 1] = 1
    hist, bins = np.histogram(dt[tag==-1], bins=nbins, range=edge, normed=True)
    scale = raw_hist[100] / hist[100]
    yerr = np.sqrt(raw_hist) / scale
    x = 0.5 * (bins[:-1] + bins[1:])
    popt, pcov = curve_fit(pdfn, x, hist, [0.72, 0.])
    print(popt[0], '+-', np.sqrt(pcov[0, 0]))
    plt.figure(num=1)
    plt.errorbar(x, hist, xerr=0, yerr=yerr, ecolor='silver', markerfacecolor='steelblue',\
                 markeredgecolor='steelblue', linestyle=' ', marker='.')
    plt.plot(x, pdfn(x, *popt))
    plt.grid()
    plt.tight_layout()
    plt.semilogy()

    plt.figure(num=3)
    pull = (hist - pdfn(x, *popt)) / yerr
    plt.plot(x, pull, marker='.', markerfacecolor='steelblue',\
                 markeredgecolor='steelblue', linestyle=' ')
    plt.grid()
    plt.tight_layout()
    print(np.sin(2. * 23. * np.pi / 180.))

    plt.show()

def wf_full_test(fname, rb):
    """ Analyse result of the full fit with 17 free parameters """
    print('Reading {}'.format(fname))
    print('  rB = {:.3f}'.format(rb))
    rsmpl, rcorr = full_fill(fname)
    diff, mean, maxd = get_diff('beta', rsmpl, rcorr, 23.)
    print('Mean bias: {:.3f}'.format(mean))
    print('Max  bias: {:.3f}'.format(maxd))

    dhist, bins = np.histogram(diff, bins=20)
    dwidth = 0.75 * (bins[1] - bins[0])
    delx = 0.5 * (bins[1:] + bins[:-1])

    plt.figure()
    plt.title(r'$\beta$ fit bias (deg), $r_B = {:.2f}$'.format(rb))
    plt.bar(delx, dhist, dwidth)
    plt.grid()
    plt.tight_layout()
    plt.show()

def wf_test(fname):
    """ Systematics because of WF """
    rb = fname.split('_')[2]
    rb = int(rb) * 10**(-len(rb))
    print('rb = {}'.format(rb))
    data = np.array(fill(fname))
    fsmpl, fcorr = data[:100, 0], data[100:200, 0]
    fsmple, fcorre = data[:100, 1], data[100:200, 1]

    mask = abs(fcorr - fsmpl) < .5
    fsmpl = fsmpl[mask]
    fcorr = fcorr[mask]
    fsmple = fsmple[mask]
    fcorre = fcorre[mask]

    plt.rc('font', size=16)

    diff = fsmpl - fcorr
    maxdev = np.argmax(abs(diff))
    print('max deviation: {:.3f}, sample {}'.format(diff[maxdev], maxdev))
    print('mean deviation: {:.3f}'.format(np.mean(abs(diff))))
    print('  smpl: {:.3f} +- {:.3f}'.format(fsmpl[maxdev], fsmple[maxdev]))
    print('  corr: {:.3f} +- {:.3f}'.format(fcorr[maxdev], fcorre[maxdev]))

    dhist, bins = np.histogram(diff, bins=20)
    dwidth = 0.75 * (bins[1] - bins[0])
    delx = 0.5 * (bins[1:] + bins[:-1])

    plt.figure(num=1)
    plt.title(r'Fit difference (deg), $r_B = {:.4f}$'.format(rb))
    plt.bar(delx, dhist, dwidth)
    plt.grid()
    plt.tight_layout()

    # pull = diff / fsmple
    # phist, bins = np.histogram(pull, bins=20)
    # pwidth = 0.75 * (bins[1] - bins[0])
    # px = 0.5 * (bins[1:] + bins[:-1])

    # plt.figure(num=2)
    # plt.title(r'Pull distibution, $r_B = {:.3f}$'.format(rb))
    # plt.bar(px, phist, pwidth)
    # plt.grid()
    # plt.tight_layout()

    plt.show()

def plot_mean_and_max():
    """ Max and mean deviations """
    plt.rc('font', size=16)
    data = np.array(
        [
            [0.0000, 0.000, 0.002],
            [0.0001, 0.000, 0.000],
            [0.0005, 0.001, 0.008],
            [0.0008, 0.001, 0.003],
            [0.0012, 0.002, 0.005],
            [0.0100, 0.016, 0.041],
            [0.0150, 0.025, 0.068],
            [0.0200, 0.039, 0.096],
            [0.0300, 0.053, 0.140],
            [0.0400, 0.080, 0.242],
            [0.0500, 0.098, 0.202]
        ])
    rb, mean, maxd = data[:, 0], data[:, 1], data[:, 2]
    plt.plot(rb, mean, marker='.', markerfacecolor='steelblue',\
                 markeredgecolor='steelblue', linestyle=' ', markersize=15)
    plt.plot(rb, maxd, marker='.', markerfacecolor='orchid',\
                 markeredgecolor='orchid', linestyle=' ', markersize=15)
    plt.grid()
    plt.xlabel(r'$r_B$')
    plt.ylabel(r'$\delta\beta$ (deg)')
    # plt.ylim([0., 0.15])
    plt.tight_layout()
    # plt.semilogx()
    # plt.semilogy()

    plt.show()

def plot_full_mean_and_max():
    """ Max and mean deviations """
    plt.rc('font', size=16)
    data = np.array(
        [
            [0.00, 0.000, 0.000],
            [0.02, 0.064, 0.168],
            [0.05, 0.161, 0.527],
            [0.10, 0.363, 1.189]
        ])
    rb, mean, maxd = data[:, 0], data[:, 1], data[:, 2]
    plt.plot(rb, mean, marker='.', markerfacecolor='steelblue',\
                 markeredgecolor='steelblue', linestyle=' ', markersize=15,
                 label='fullMean')
    plt.plot(rb, maxd, marker='.', markerfacecolor='orchid',\
                 markeredgecolor='orchid', linestyle=' ', markersize=15,
                 label='fullMax')
    plt.plot([0.02, 0.02], [0.036, 0.143], marker='.', markerfacecolor='k',\
                 markeredgecolor='k', linestyle=' ', markersize=15,
                 label='posiCP')
    plt.plot([0.02, 0.02], [0.052, 0.149], marker='.', markerfacecolor='r',\
                 markeredgecolor='r', linestyle=' ', markersize=15,
                 label='noCP')
    plt.plot([0.02, 0.02], [0.069, 0.268], marker='.', markerfacecolor='b',\
                 markeredgecolor='b', linestyle=' ', markersize=15,
                 label='negaCP')
    plt.grid()
    plt.xlabel(r'$r_B$')
    plt.ylabel(r'$\delta\beta$ (deg)')
    plt.tight_layout()
    plt.show()

# plot_mean_and_max()

# wf_test(sys.argv[1])

wf_full_test(sys.argv[1], float(sys.argv[2]))
# plot_full_mean_and_max()

# DPATH = '/home/vitaly/B0toD0pipi/btoucbard/data/'
# FNAME = 'cm_x_0.000000_y_0.000000.txt'
# FNAME1 = 'ddist_rb0_02_seed_8648_idx_0.txt'
# FNAME2 = 'cm_x_0.000000_y_0.000000.txt'
# DATA = read_dt_dist(DPATH + FNAME)
# make_fit(DATA[0], DATA[1])
# plot_dt(DATA[0], DATA[1])
