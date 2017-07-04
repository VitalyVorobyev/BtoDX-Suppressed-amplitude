#! /usr/bin/python3
""" Numerical tests of analytic tramsformations """

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

DEG_TO_RAD = np.pi / 180.

def ccoef_cp_fcn(delb, C, S):
    """ Bias fit function """
    return C * np.cos(delb) + S * np.sin(delb)

def ccoef_fit(delb, ccoef, rb, gamma):
    """ Run fit """
    popt, pcov = curve_fit(ccoef_cp_fcn, delb, ccoef, [0., -2.])
    scale = rb * np.sin(gamma)
    for idx, par in enumerate(popt):
        print(par / scale, '+-', np.sqrt(pcov[idx, idx]) / scale)
    return popt

# """ Low level computations """
def acubard(rd, rb, deld, delb, gamma):
    """ B0 -> D0bar h0 decay amplitude """
    return 1. + rb * rd * np.exp(1j*(-deld + delb + gamma))

def aucbard(rd, rb, deld, delb, gamma):
    """ B0 -> D0 h0 decay amplitude """
    return rd * np.exp(-1j*deld) + rb * np.exp(1j*(delb - gamma))

def lamdaf(rd, rb, deld, delb, gamma, beta):
    """ lambda_f """
    return np.exp(-2.*1j*beta) *\
           aucbard(rd, rb, deld, delb, gamma) /\
           acubard(rd, rb, deld, delb, gamma)

def ccoef(lamf):
    """ C * cos(dt) """
    abssq = np.abs(lamf)**2
    return (1. - abssq) / (1. + abssq)

def scoef(lamf):
    """ S * sin(dt) """
    return 2. * np.imag(lamf) / (1. + np.abs(lamf)**2)
# """ end of low level computations """

# def cp_denominator(rb, delb, deld, gamma):
#     """ An exact expression """
#     return 1. + 2.*rb*np.cos(delb) * np.cos(gamma) * np.cos(deld) + rb**2

# def exact_ccoef_cp(rd, rb, deld, delb, gamma, beta):
#     """ An exact expression """
#     return -2. * rb * np.sin(gamma) * np.sin(delb) * np.cos(deld) /\
#            cp_denominator(rb, delb, deld, gamma)

# def exact_scoef_cp(rd, rb, deld, delb, gamma, beta):
#     """ An exact expression """
#     return -(np.sin(2.*beta+deld)+2.*rb*np.cos(delb)*np.sin(gamma+2.*beta)+\
#             rb**2*np.sin(2.*beta+2.*gamma)) / cp_denominator(rb, delb, deld, gamma)

# """ The first order expressions for CP specific decays """
def ccoef_cp(rb, xid, delb, gamma):
    """ Analytical result up to the first order of rb """
    return -2. * xid * rb * np.sin(delb) * np.sin(gamma)# +\
    # rb**2 * np.sin(2.*delb) * np.sin(2.*gamma)

def scoef_cp(rb, xid, delb, gamma, beta):
    """ Analytical result up to the first order of rb """
    return -xid*np.sin(2.*beta) - 2.*rb*np.cos(delb)*np.sin(gamma)*np.cos(2.*beta)
# """ end of approximations for CP """

# """ The first order expressions for D0 -> K pi """
def ccoef_kpi():
    """ Analytical result up to the first order of rb """
    return 1.

def scoef_kpi(rd, rb, deld, delb, gamma, beta):
    """ Analytical result up to the first order of rb """
    return -2. * rd * np.sin(deld + 2.*beta) - 2.*rb*np.sin(gamma + 2.*beta - delb)

def scoef_pik(rd, rb, deld, delb, gamma, beta):
    """ Analytical result up to the first order of rb """
    return -2. / rd * np.sin(deld + 2.*beta) - 2.*rb*np.sin(gamma + 2.*beta + delb)
# """ end of the first order expressions for D0 -> K pi """

# """ The first order expressions for rd ~ 1 """
def ccoef_fo(rd, rb, deld, delb, gamma):
    """ Analytical result up to the first order of rb """
    rdsq = rd**2
    return (1. - rdsq**2 - 4.*rb*rd*(np.cos(gamma - deld - delb) -\
    rdsq * np.cos(gamma - deld + delb))) / (1 + rdsq)**2

def scoef_fo(rd, rb, deld, delb, gamma, beta):
    """ Analytical result up to the first order of rb """
    rdsq = rd**2
    return -2.*rd*np.sin(deld+2.*beta) / (1. + rdsq) -\
    2.*rb / (1 + rdsq)**2 * (np.sin(gamma + 2.*beta - delb) -\
    2.*rdsq*np.cos(delb)*np.sin(2.*deld+2.*beta-gamma)+\
    rdsq**2*np.sin(gamma+2.*beta+delb))
# """ end of the first order expressions for rd ~ 1 """

# """ Exact high level expressions """
def one_plus_lamfsq(rd, rb, deld, delb, gamma):
    """ An exact expression for numerator """
    return 1. + rd**2 + rb**2 + rd**2 * rb**2 +\
           4.*rd*rb*np.cos(delb)*np.cos(gamma - deld)

def one_minus_lamfsq(rd, rb, deld, delb, gamma):
    """ An exact expression for numerator """
    return 1. - rd**2 - rb**2 + rd**2 * rb**2 -\
           4.*rd*rb*np.sin(delb)*np.sin(gamma - deld)

def im_lamf(rd, rb, deld, delb, gamma, beta):
    """ An exact expression for numerator """
    return -rd*np.sin(deld + 2.*beta) +\
            rb*np.sin(delb - gamma - 2.*beta) -\
            rb*rd**2*np.sin(delb + gamma + 2.*beta) -\
            rb**2*rd*np.sin(2.*beta + 2.*gamma - deld)

def exact_ccoef(rd, rb, deld, delb, gamma):
    """ An exact expression """
    return one_minus_lamfsq(rd, rb, deld, delb, gamma) /\
           one_plus_lamfsq(rd, rb, deld, delb, gamma)

def exact_scoef(rd, rb, deld, delb, gamma, beta):
    """ An exact expression """
    return 2.*im_lamf(rd, rb, deld, delb, gamma, beta) /\
           one_plus_lamfsq(rd, rb, deld, delb, gamma)
# """ end of exact high level expression """

def plot_cs(rd, rb, deld, gamma, beta):
    """ Plot coefficients C and S as a function of deld """
    gamma, beta = DEG_TO_RAD * gamma, DEG_TO_RAD * beta
    deld = DEG_TO_RAD * deld
    xaxis = np.linspace(0., 361., 400)
    delb = DEG_TO_RAD * xaxis
    lamb = lamdaf(rd, rb, deld, delb, gamma, beta)
    lamb0 = lamdaf(rd, 0, deld, delb, gamma, beta)
    plt.figure(num=1)
    plt.plot(xaxis, ccoef(lamb), 'b-')
    plt.plot(xaxis, ccoef(lamb0), 'k-')
    plt.plot(xaxis, exact_ccoef(rd, rb, deld, delb, gamma), 'r:')
    # plt.plot(xaxis, ccoef_posi_cp(rd, rb, deld, delb, gamma, beta), 'r:')
    # plt.plot(xaxis, ccoef(lamb0)+ccoef(lamb)-ccoef_posi_cp(rd, rb, deld, delb, gamma, beta), 'y-')
    plt.title('cos coefficient')
    plt.grid()
    plt.tight_layout()
    plt.figure(num=2)
    plt.plot(xaxis, scoef(lamb), 'b-')
    plt.plot(xaxis, scoef(lamb0), 'k-')
    plt.plot(xaxis, exact_scoef(rd, rb, deld, delb, gamma, beta), 'r:')
    # plt.plot(xaxis, scoef_posi_cp(rd, rb, deld, delb, gamma, beta), 'r:')
    # plt.plot(xaxis, scoef(lamb0)+scoef(lamb)-scoef_posi_cp(rd, rb, deld, delb, gamma, beta), 'y-')
    plt.title('sin coefficient')
    plt.grid()
    plt.tight_layout()

    print(ccoef(lamb0[0]), scoef(lamb0[0]))
    plt.show()

if __name__ == '__main__':
    # plot_cs(0.063, 0.05, 10., 71., 23.)
    plot_cs(1./0.063, 0.05, 10., 71., 23.)
    # plot_cs(1, 0.05, 0., 71., 23.)
    # plot_cs(rd=1, rb=0.05, deld=180., gamma=71., beta=23.)
