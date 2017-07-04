#! /usr/bin/python3
""" Corrections for rd ~ 1 """

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from lamf import ccoef, scoef, lamdaf
from lamf import ccoef_fo, scoef_fo

plt.rc('font', size=14)

DEG_TO_RAD = np.pi / 180.
PIC_PATH = '/home/vitaly/B0toD0pipi/btoucbard/pics/'

def params(delb, deld, zero=False):
    """ Dictionary with parameters """
    return {
        'deld' : deld * DEG_TO_RAD,
        'rd' : 0.5,
        'rb' : 0.05 if not zero else 0.,
        'delb' : delb,
        'gamma' : 71. * DEG_TO_RAD,
        'beta' : 23.  * DEG_TO_RAD
    }

def cos_params(delb, deld, zero=False):
    """ Dictionary with parameters """
    return {
        'deld' : deld * DEG_TO_RAD,
        'rd' : 0.5,
        'rb' : 0.05 if not zero else 0.,
        'delb' : delb,
        'gamma' : 71. * DEG_TO_RAD,
    }

XAXIS = np.linspace(0., 361., 400)
DELB = DEG_TO_RAD * XAXIS

def cos_calc(deld):
    p, p0 = params(DELB, deld), params(DELB, deld, True)
    pcp = cos_params(DELB, deld)
    lamb = lamdaf(**p)
    lamb0 = lamdaf(**p0)
    exact_coef = ccoef(lamb)
    first_order_coef = ccoef_fo(**pcp)
    zero_rb = ccoef(lamb0)
    first_order_bias = zero_rb + exact_coef - first_order_coef
    return [zero_rb, exact_coef, first_order_coef, first_order_bias]

def sin_calc(deld):
    p, p0 = params(DELB, deld), params(DELB, deld, True)
    lamb = lamdaf(**p)
    lamb0 = lamdaf(**p0)
    exact_coef = scoef(lamb)
    first_order_coef = scoef_fo(**p)
    zero_rb = scoef(lamb0)
    first_order_bias = zero_rb + exact_coef - first_order_coef
    return [zero_rb, exact_coef, first_order_coef, first_order_bias]

def save_plt(fig, pref, deld):
    """ Save plot """
    fig.savefig(PIC_PATH+'_'.join(['rd_deld', str(deld), pref])+'.eps', format='eps', dpi=150)
    fig.savefig(PIC_PATH+'pdf/'+'_'.join(['rd_deld', str(deld), pref])+'.pdf', format='pdf', dpi=150)
    fig.savefig(PIC_PATH+'png/'+'_'.join(['rd_deld', str(deld), pref])+'.png', format='png', dpi=150)

def cos_plot(deld):
    zero_rb, exact_coef, first_order_coef, first_order_bias = cos_calc(deld)
    fig = plt.figure(num=1)
    plt.plot(XAXIS, exact_coef, 'b-', label='Exact')
    plt.plot(XAXIS, zero_rb, 'k-', label='Zero rb')
    plt.plot(XAXIS, first_order_coef, 'r-', label='First order')
    plt.plot(XAXIS, first_order_bias, 'y-', label='First order bias')
    plt.xticks(np.arange(0, 361, 45))
    plt.gca().legend(loc='best', shadow=True)
    plt.xlabel(r'$\delta_B\ (\mathrm{deg})$')
    plt.title(r'cos coefficient $\delta_{D}=' + str(deld) + r'$')
    plt.grid()
    plt.tight_layout()
    save_plt(fig, 'cos', deld)

def sin_plot(deld):
    zero_rb, exact_coef, first_order_coef, first_order_bias = sin_calc(deld)
    fig = plt.figure(num=2)
    plt.plot(XAXIS, exact_coef, 'b-', label='Exact')
    plt.plot(XAXIS, zero_rb, 'k-', label='Zero rb')
    plt.plot(XAXIS, first_order_coef, 'r-', label='First order')
    plt.plot(XAXIS, first_order_bias, 'y-', label='First order bias')
    plt.xticks(np.arange(0, 361, 45))
    plt.gca().legend(loc='best', shadow=True)
    plt.xlabel(r'$\delta_B\ (\mathrm{deg})$')
    plt.title(r'sin coefficient $\delta_{D}=' + str(deld) + r'$')
    plt.grid()
    plt.tight_layout()
    save_plt(fig, 'sin', deld)

FIG = plt.figure()
AX = plt.axes(xlim=(0, 2), ylim=(0, 100))
LINES = [plt.plot([], [])[0] for _ in range(4)]

def init():    
    for line in LINES:
        line.set_data([], [])
    return LINES

def cos_animate(i):
    data = cos_calc(i / 360.)
    for vals, line in zip(data, LINES):
        line.set_data(XAXIS, vals)
    return LINES

def cos_animat():
    """ Make animation """

# anim = animation.FuncAnimation(fig, animate, init_func=init,
#                                frames=100, interval=20, blit=True)

# plt.show()

DELD = int(sys.argv[1])
# for deld in range(0, 361, 10):
cos_plot(DELD)
sin_plot(DELD)
plt.show()
