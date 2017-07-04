#! /usr/bin/python3
""" Corrections for D0 -> K- pi+ """

import numpy as np
import matplotlib.pyplot as plt

from lamf import ccoef, scoef, lamdaf
from lamf import scoef_kpi, scoef_pik

DEG_TO_RAD = np.pi / 180.
PIC_PATH = '/home/vitaly/B0toD0pipi/btoucbard/pics/'

def params(delb, dcpv, zero=False):
    """ Dictionary with parameters """
    if not dcpv:
        rd, deld = 0.063, 10
    else:
        rd, deld = 1. / 0.063, -10
    return {
        'rd' : rd,
        'deld' : deld * DEG_TO_RAD,
        'rb' : 0.05 if not zero else 0.,
        'delb' : delb, #30. * DEG_TO_RAD,
        'gamma' : 71. * DEG_TO_RAD,
        'beta' : 23.  * DEG_TO_RAD
    }

plt.rc('font', size=14)

def save_plt(fig, pref):
    """ Save plot """
    fig.savefig(PIC_PATH + '_'.join(['kpi', pref]) + '.eps', format='eps', dpi=150)
    fig.savefig(PIC_PATH + 'pdf/' + '_'.join(['kpi', pref]) + '.pdf', format='pdf', dpi=150)
    fig.savefig(PIC_PATH + 'png/' + '_'.join(['kpi', pref]) + '.png', format='png', dpi=150)

def cos_plot(kpi, num):
    xaxis = np.linspace(0., 361., 400)
    delb = DEG_TO_RAD * xaxis
    p = params(delb, not kpi)
    p0 = params(delb, not kpi, True)
    lamb = lamdaf(**p)
    lamb0 = lamdaf(**p0)
    exact_coef = ccoef(lamb)
    fig = plt.figure(num=num)
    plt.plot(xaxis, exact_coef, 'b-', label=r'Exact')
    plt.plot(xaxis, ccoef(lamb0), 'k-', label=r'Zero rb')
    plt.xticks(np.arange(0, 361, 45))
    plt.gca().legend(loc='best', shadow=True)
    plt.xlabel(r'$\delta_B\ (\mathrm{deg})$')
    plt.title('cos coefficient ' + (r'$K\pi$' if kpi else r'$\pi K$'))
    plt.grid()
    plt.tight_layout()
    save_plt(fig, 'cos_' + ('kpi' if kpi else 'pik'))

def sin_plot():
    xaxis = np.linspace(0., 361., 400)
    delb = DEG_TO_RAD * xaxis
    p = params(delb, False)
    p0 = params(delb, False, True)
    lamb = lamdaf(**p)
    lamb0 = lamdaf(**p0)
    exact_coef = scoef(lamb)
    first_order_coef = scoef_kpi(**p)

    fig = plt.figure(num=3, figsize=(10, 7))
    plt.plot(xaxis, exact_coef, 'b-', label='Exact Kpi')
    plt.plot(xaxis, scoef(lamb0), 'k-', label='Zero rb Kpi')
    plt.plot(xaxis, first_order_coef, 'b:', label='First order Kpi')
    plt.plot(xaxis, scoef(lamb0)+exact_coef-first_order_coef, 'y-', label='First order bias Kpi')

    q = params(delb, True)
    q0 = params(delb, True, True)
    qlamb = lamdaf(**q)
    qlamb0 = lamdaf(**q0)
    qexact_coef = scoef(qlamb)
    qfirst_order_coef = scoef_pik(**q)

    plt.plot(xaxis, qexact_coef, 'r-', label='Exact piK')
    plt.plot(xaxis, scoef(qlamb0), 'm-', label='Zero rb piK')
    plt.plot(xaxis, qfirst_order_coef, 'r:', label='First order piK')
    plt.plot(xaxis, scoef(qlamb0)+qexact_coef-qfirst_order_coef, 'g-', label='First order bias piK')

    plt.xticks(np.arange(0, 361, 45))
    plt.gca().legend(loc='best', shadow=True)
    plt.xlabel(r'$\delta_B\ (\mathrm{deg})$')
    plt.title(r'sin coefficient $K\pi$')
    plt.grid()
    plt.tight_layout()
    save_plt(fig, 'sin')

# cos_plot(True, 1)
# cos_plot(False, 2)
sin_plot()
plt.show()
