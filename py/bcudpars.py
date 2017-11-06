#! /usr/bin/python3
""" Analysis of the binned Dalitz plot parameters """

import numpy as np
import matplotlib.pyplot as plt

def style(font_size):
    """ Plot style """
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rc('font', size=font_size)
    plt.style.use('seaborn-white')

def cfg_dict():
    """ Title dict """
    return {
        'K+rf' : [r'$k_{j}$',               'kjp',    [0., 0.2]],
        'K-rf' : [r'$k_{-j}$',              'kjn',    [0., 0.2]],
        'K+wf' : [r'$\overline{k}{}_{j}$',  'kbarjp', [0., 0.2]],
        'K-wf' : [r'$\overline{k}{}_{-j}$', 'kbarjn', [0., 0.2]],
        'Cwf' : [r'$\overline{c}{}_{j}$',   'cbarj',  [-1., 1.]],
        'Swf' : [r'$\overline{s}{}_{j}$',   'sbarj',  [-1., 1.]],
        'Crf' : [r'$c_j$',                  'cj',     [-1., 1.]],
        'Srf' : [r'$s_j$',                  'sj',     [-1., 1.]],
        'Ctp' : [r'$\widetilde{c}{}_j$',    'ctjp',   [-1., 1.]],
        'Ctn' : [r'$\widetilde{c}{}_{-j}$', 'ctjn',   [-1., 1.]],
        'Stp' : [r'$\widetilde{s}{}_j$',    'stjp',   [-1., 1.]],
        'Stn' : [r'$\widetilde{s}{}_{-j}$', 'stjn',   [-1., 1.]],
        'Cpp' : [r'$c^{\prime}_j$',         'cpjp',   [-1., 1.]],
        'Cpn' : [r'$c^{\prime}_{-j}$',      'cpjn',   [-1., 1.]],
        'Spp' : [r'$s^{\prime}_j$',         'spjp',   [-1., 1.]],
        'Spn' : [r'$s^{\prime}_{-j}$',      'spjn',   [-1., 1.]]
    }

def parnames():
    """ List of parameter names """
    return cfg_dict().keys()

def pars_path():
    """ Path to directory with text files """
    return '/home/vitaly/B0toD0pipi/B0toD0pipiFeas/params/'

def pars_fname(seed, idx):
    """ File name """
    return 'wfpars_wf_tblSymABBC_bdpp_wf_seed_' + str(seed)\
           + '_idx_' + str(idx) + '.txt'

def read_table(seed, idx):
    """ Read a parameters set from text file """
    pars = {}
    for line in open(pars_path() + pars_fname(seed, idx)).readlines():
        content = line.split()
        idx = 1 if len(content) == 13 else 0
        for _ in range(4):
            val = float(content[idx+2][:-1])
            if content[idx] in pars:
                pars[content[idx]].append(val)
            else:
                pars[content[idx]] = [val]
            idx += 3
    return pars

def read_data(seed, idxlo, idxhi):
    """ Read data """
    data = {}
    for key in cfg_dict().keys():
        data[key] = [[] for _ in range(8)]
    for idx in range(idxlo, idxhi):
        table = read_table(seed, idx)
        for key, arr in table.items():
            for jdx in range(8):
                data[key][jdx].append(arr[jdx])
    return data

def make_hist(name, binn, data, num):
    """ Distribution of parameter value """
    array = np.array(data[name][binn - 1])
    print('bin {}, par {}: mean {:.3f}, stddev {:.3f}'.format(binn, name,\
          np.mean(array), np.std(array)))
    hist, bins = np.histogram(array, bins=20)
    width = 0.75 * (bins[1] - bins[0])
    x = 0.5 * (bins[1:] + bins[:-1])
    plt.figure(num=num)
    plt.bar(x, hist, width=width, color='steelblue')
    plt.title(cfg_dict()[name][0] + ' bin ' + str(binn))
    plt.tight_layout()
    plt.gca().yaxis.grid()

def save_plot(pref, fig):
    """ Save pdf, eps and png """
    pdi = 150
    fig.savefig('pics/' + pref + '.eps', format='eps', dpi=pdi)
    fig.savefig('pics/pdf/' + pref + '.pdf', format='pdf', dpi=pdi)
    fig.savefig('pics/png/' + pref + '.png', format='png', dpi=pdi)

def plot_mean_std(name, data):
    """ Plot mean and std dev for each bin """
    mean, stdd = [], []
    maxv, minv = [], []
    for binn in range(8):
        array = np.array(data[name][binn])
        mean.append(np.mean(array))
        stdd.append(np.std(array))
        minv.append(min(array))
        maxv.append(max(array))
    fig = plt.figure(num=1, figsize=(8, 5))
    col = 'black'
    plt.errorbar(range(1, 9), mean, xerr=0, yerr=stdd, ecolor=col,\
                 markerfacecolor=col, markersize=10,\
                 markeredgecolor=col, linestyle=' ', marker='o')
    if (max(stdd) > 0.001):
        colmax = 'royalblue'
        plt.plot(range(1, 9), maxv, marker='v', markerfacecolor=colmax, markersize=10,\
                    markeredgecolor=colmax, linestyle=' ')
        colmin = 'firebrick'
        plt.plot(range(1, 9), minv, marker='^', markerfacecolor=colmin, markersize=10,\
                    markeredgecolor=colmin, linestyle=' ')
    plt.title(cfg_dict()[name][0])
    ymin, ymax = cfg_dict()[name][2]
    plt.ylim(ymin=ymin, ymax=ymax)
    plt.xticks(range(1, 9))
    plt.gca().set_xlabel(r'$\textrm{Dalitz plot bin index}$')
    delta = (ymax - ymin) / 8
    plt.yticks(np.arange(ymin, ymax + 0.01, delta))
    plt.tight_layout(pad=0.15)
    plt.grid()
    save_plot('mean_std_' + cfg_dict()[name][1], fig)

def make_plots(data):
    """ Plots for each variable """
    for key in parnames():
        plot_mean_std(key, data)
        plt.show()

style(24)
DATA = read_data(8648, 0, 100)
make_plots(DATA)
# plot_mean_std('Cwf', DATA)
# for idx in range(8):
#     make_hist('Spp', idx + 1, DATA, idx + 1)
# plt.show()
