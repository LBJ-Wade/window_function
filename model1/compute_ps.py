import sys
import math
import json
import numpy as np
from spherical_bessel_integral import SphericalBesselIntegral
from scipy.interpolate import interp1d
from scipy.special import legendre
from scipy.integrate import quad
from tqdm import tqdm

def compute_power_spectrum(xis, *, kind='linear', k_unit=0, k_min=0.5, k_max=256, dk=2, boxsize=3400.0):
    if k_unit == 0:
        # k is in units of fundamental frequency
        kf = 2.0*math.pi/boxsize
        dk = dk*kf
        k_min = dk*k_min
        k_max = dk*k_max
        
    nk = round((k_max - k_min)/dk)
    ks = k_min + np.arange(nk)*dk
    
    r = xis['r']
    xi0 = xis['xi0']
    xi2 = xis['xi2']
    xi4 = xis['xi4']
    
    P0 = np.zeros_like(ks)
    P2 = np.zeros_like(ks)
    P4 = np.zeros_like(ks)

    for i, k in enumerate(tqdm(ks)):
        fac = 4.0*math.pi/k**3
            
        kr = k*r
        s = SphericalBesselIntegral(kr)

        # integrate(f, n, l) = \int f(x) x^n j_l(x) dx
        P0[i] = fac*s.integrate(xi0, 2, 0, kind=kind)
        P2[i] = fac*s.integrate(xi2, 2, 2, kind=kind)
        P4[i] = fac*s.integrate(xi4, 2, 4, kind=kind)

    P2 = -P2  # (-i)^l
    ps = {'k':ks, 'P0':P0, 'P2':P2, 'P4':P4}

    return ps

def apply_window(filename, xis):
    """
    Apply window function multipoles
    """

    a = np.loadtxt(filename)

    
    r_w = 0.5*(a[:, 0] + a[:, 1])
    r_xi  = xis['r']

    # window function is zero beyond the last point -- cut xi there
    idx = r_xi < r_w[-1]
    r_xi = r_xi[idx]
    xi0 = xis['xi0'][idx]
    xi2 = xis['xi2'][idx]
    xi4 = xis['xi4'][idx]

    # interpolate W on xi points
    W0 = interp1d(r_w, a[:, 2], bounds_error=False, fill_value=1.0)(r_xi)
    W2 = interp1d(r_w, a[:, 3], bounds_error=False, fill_value=0.0)(r_xi)
    W4 = interp1d(r_w, a[:, 4], bounds_error=False, fill_value=0.0)(r_xi)
    W6 = interp1d(r_w, a[:, 5], bounds_error=False, fill_value=0.0)(r_xi)
    W8 = interp1d(r_w, a[:, 6], bounds_error=False, fill_value=0.0)(r_xi)


    # DEBUG!!!
    xi0_w = xi0*W0 + (1.0/5.0)*xi2*W2 + (1.0/9.0)*xi4*W4
    xi2_w = (xi2*W2
             + xi2*(W0 + 2.0/7.0*W2 + 2.0/7.0*W4)
             + xi4*(2.0/7.0*W2 + 100.0/693.0*W4 + 25.0/143.0*W6))
    xi4_w = (xi0*W4
             + xi2*(18.0/35.0*W2 + 20.0/77.0*W4 + 45.0/143*W6)
             + xi4*(W0 + 20.0/77.7*W2 + 162.0/1001.0*W4
                    + 20.0/143.0*W6 + 490.0/2431.0*W8))

    return {'r':r_xi, 'xi0':xi0_w, 'xi2':xi2_w, 'xi4':xi4_w}

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='corr npz file name')
    parser.add_argument('--window', help='window function multipole file name')
    parser.add_argument('--print-xi', action='store_true', help='print xi')

    parser.add_argument('--interpolation', default='linear', help='interpolation method of spherical bessel integral')
    
    arg = parser.parse_args()
    xis = np.load(arg.filename)
    
    if arg.window:
        print('# window ', arg.window)
        xis = apply_window(arg.window, xis)

    if arg.print_xi:
        rs = xis['r']
        xi0 = xis['xi0']
        xi2 = xis['xi2']
        xi4 = xis['xi4']
        for i, r in enumerate(rs):
            print('%e %e %e %e' % (r, xi0[i], xi2[i], xi4[i]))
        sys.exit()
    
    ps = compute_power_spectrum(xis, kind=arg.interpolation, dk=2)

    ks = ps['k']
    P0 = ps['P0']
    P2 = ps['P2']
    P4 = ps['P4']
    for i, k in enumerate(ks):
        print('%e %e %e %e' % (k, P0[i], P2[i], P4[i]))

if __name__ == "__main__":
    main()

