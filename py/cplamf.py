#! /usr/bin/python3
""" Corrections for CP specific decays """

import sys
import numpy as np
import matplotlib.pyplot as plt

from lamf import ccoef, scoef, lamdaf
from lamf import ccoef_cp, scoef_cp

DEG_TO_RAD = np.pi / 180.
PIC_PATH = '/home/vitaly/B0toD0pipi/btoucbard/pics/'

def params(delb, xid, zero=False):
    """ Dictionary with parameters """
    return {
        'deld' : 0. if xid == 1 else np.pi,
        'rd' : 1,
        'rb' : 0.05 if not zero else 0.,
        'delb' : delb,
        'gamma' : 71. * DEG_TO_RAD,
        'beta' : 23.  * DEG_TO_RAD
    }

def cos_cp_params(delb, xid):
    """ Dictionary with parameters """
    return {
        'xid' : xid,
        'rb' : 0.05,
        'delb' : delb,
        'gamma' : 71. * DEG_TO_RAD
    }

def sin_cp_params(delb, xid):
    """ Dictionary with parameters """
    return {
        'xid' : xid,
        'rb' : 0.05,
        'delb' : delb,
        'gamma' : 71. * DEG_TO_RAD,
        'beta' : 23.  * DEG_TO_RAD
    }

plt.rc('font', size=14)

def save_plt(fig, pref, xid):
    """ Save plot """
    cptype = 'posi' if xid == 1 else 'nega'
    fig.savefig(PIC_PATH + '_'.join(['cp', cptype, pref]) + '.eps', format='eps', dpi=150)
    fig.savefig(PIC_PATH + 'pdf/' + '_'.join(['cp', cptype, pref]) + '.pdf', format='pdf', dpi=150)
    fig.savefig(PIC_PATH + 'png/' + '_'.join(['cp', cptype, pref]) + '.png', format='png', dpi=150)

def cos_plot(xid):
    xaxis = np.linspace(0., 361., 400)
    delb = DEG_TO_RAD * xaxis
    p, p0 = params(delb, xid), params(delb, xid, True)
    pcp = cos_cp_params(delb, xid)
    lamb = lamdaf(**p)
    lamb0 = lamdaf(**p0)
    exact_coef = ccoef(lamb)
    first_order_coef = ccoef_cp(**pcp)
    fig = plt.figure(num=1)
    plt.plot(xaxis, exact_coef, 'b-', label='Exact')
    plt.plot(xaxis, ccoef(lamb0), 'k-', label='Zero rb')
    plt.plot(xaxis, first_order_coef, 'r-', label='First order')
    plt.plot(xaxis, ccoef(lamb0)+exact_coef-first_order_coef, 'y-', label='First order bias')
    plt.xticks(np.arange(0, 361, 45))
    plt.gca().legend(loc='best', shadow=True)
    plt.xlabel(r'$\delta_B\ (\mathrm{deg})$')
    plt.title('cos coefficient CP ' + ('+1' if xid == 1 else '-1'))
    plt.grid()
    plt.tight_layout()
    save_plt(fig, 'cos', xid)

def sin_plot(xid):
    xaxis = np.linspace(0., 361., 400)
    delb = DEG_TO_RAD * xaxis
    p, p0 = params(delb, xid), params(delb, xid, True)
    pcp = sin_cp_params(delb, xid)
    lamb = lamdaf(**p)
    lamb0 = lamdaf(**p0)
    exact_coef = scoef(lamb)
    first_order_coef = scoef_cp(**pcp)
    fig = plt.figure(num=2)
    plt.plot(xaxis, exact_coef, 'b-', label='Exact')
    plt.plot(xaxis, scoef(lamb0), 'k-', label='Zero rb')
    # plt.plot(xaxis, exact_scoef(rd, rb, deld, delb, gamma, beta), 'g:', label='High level')
    plt.plot(xaxis, first_order_coef, 'r-', label='First order')
    plt.plot(xaxis, scoef(lamb0)+exact_coef-first_order_coef, 'y-', label='First order bias')
    plt.xticks(np.arange(0, 361, 45))
    plt.gca().legend(loc='best', shadow=True)
    plt.xlabel(r'$\delta_B\ (\mathrm{deg})$')
    plt.title('sin coefficient CP ' + ('+1' if xid == 1 else '-1'))
    plt.grid()
    plt.tight_layout()
    save_plt(fig, 'sin', xid)

XID = int(sys.argv[1])
cos_plot(XID)
sin_plot(XID)
plt.show()
