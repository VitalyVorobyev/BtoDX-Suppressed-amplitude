#! /usr/bin/python3
""" Explore local minima in the full fit procedure """

import re
import numpy as np
import matplotlib.pyplot as plt

DATA_PATH  = '/home/vitaly/B0toD0pipi/btoucbard/data/'

btfitre = re.compile(r'beta = [0-9]+.[0-9]+ \+- [0-9]+.[0-9]+')
# rbfitre = re.compile(r'  rb = [0-9]+.[0-9]+ \+- [0-9]+.[0-9]+')
# zerorbfitre = re.compile(r'  rb = 0 \+- [0-9]+.[0-9]+')
rbfitre = re.compile(r'  rb = [0-9,.,e,-]* \+- [0-9,.,e,-]*')
dbfitre = re.compile(r'delb = [-]?[0-9]+.[0-9]+ \+- [0-9]+.[0-9]+')

def parse_fit_result(line):
    vals = line.split()
    return [float(vals[2]), float(vals[4])]

def get_data(fname):
    """ Read fit results """
    beta, rb, db = [], [], []
    for line in open(fname):
        if btfitre.match(line):
            beta.append(parse_fit_result(line))
        elif rbfitre.match(line):
            rb.append(parse_fit_result(line))
        elif dbfitre.match(line):
            db.append(parse_fit_result(line))
    return [np.array(beta), np.array(rb), np.array(db)]

def plot_hists(data, num, pars):
    """ Fit histograms for beta, rb and db """
    title, lims, step = pars
    plt.figure(num=num)
    plt.title(title)
    plt.xticks(np.arange(lims[0], lims[1]+0.01, step))
    plt.hist(data, bins=200)
    plt.grid(axis='x')
    plt.tight_layout()

def plot_scatter(data1, data2, num, parx, pary):
    titlex, limsx, stepx = parx
    titley, limsy, stepy = pary
    plt.figure(num=num)
    plt.scatter(data1, data2)
    plt.xticks(np.arange(limsx[0], limsx[1]+0.01, stepx))
    plt.yticks(np.arange(limsy[0], limsy[1]+0.01, stepy))
    plt.title(titlex + ' vs. ' + titley)
    plt.grid()
    plt.tight_layout()

DATA = get_data(DATA_PATH + 'full_fit_500_times.txt')
# beta, rb, db = DATA[0][:, 0], DATA[1][:, 0], DATA[2][:, 0]

# beta[beta > 45] = 90. - beta[beta > 45]

beta = DATA[0][:, 0]
rb = DATA[1][:, 0]
db = DATA[2][:, 0]

beta = beta[rb>0.0001]
db = db[rb>0.0001]
rb = rb[rb>0.0001]
print(len(db))
while len(db[db>180]) > 0:
    db[db>180] = db[db>180] - 360
while len(db[db<-170]) > 0:
    db[db<-170] = db[db<-170] + 360
# db[db>100] = db[db>100] - 180

beta_pars = [r'$\beta$', [0, 90], 5]
rb_pars = [r'$r_{B}$', [0, 0.4], 0.05]
db_pars = [r'$\delta_{B}$', [-30, 210], 30]

plot_hists(beta, 1, beta_pars)
plot_hists(rb, 2, rb_pars)
plot_hists(db, 3, db_pars)
print(len(beta), len(rb), len(db))
plot_scatter(rb, db, 4, rb_pars, db_pars)

plot_scatter(rb, beta, 5, rb_pars, beta_pars)

plt.show()
