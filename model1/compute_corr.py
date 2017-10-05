"""
Compute correlation function multipoles from power spectrum
"""

import sys
import math
import json
import numpy as np
from spherical_bessel_integral import SphericalBesselIntegral
from scipy.interpolate import interp1d
from scipy.special import legendre
from scipy.integrate import quad
from tqdm import tqdm

def compute_corr(ps, kind='linear'):
    """
    kind (str): interplation scheme for power spectrum
                'linear' or 'cubic'
    """
    k = ps['k']
    P0 = ps['P0']
    P2 = ps['P2']
    P4 = ps['P4']
    
    rmin = 0.01
    rmax = 1.0e5
    nbin = 1001
    rs = np.logspace(math.log10(rmin), math.log10(rmax), nbin)

    xi0 = np.zeros_like(rs)
    xi2 = np.zeros_like(rs)
    xi4 = np.zeros_like(rs)

    for i, r in enumerate(tqdm(rs)):
        fac = 1.0/(2.0*math.pi**2*r**3)
            
        kr = k*r
        s = SphericalBesselIntegral(kr)

        # integrate(f, n, l) = \int f(x) x^n j_l(x) dx
        xi0[i] = fac*s.integrate(P0, 2, 0, kind=kind)
        xi2[i] = fac*s.integrate(P2, 2, 2, kind=kind)
        xi4[i] = fac*s.integrate(P4, 2, 4, kind=kind)

    xi2 = -xi2  # i^l
    xis = {'r':rs, 'xi0':xi0, 'xi2':xi2, 'xi4':xi4}

    return xis

    
def compute_power_spectrum(filename, growth, b, f):
    """
    Compute power spectrum multipoles
    """

    a = np.loadtxt(filename)

    fac = {}
    
    for l in (0, 2, 4):
        pl = legendre(l)
        integ = quad(lambda mu: (b + f*mu**2)**2*pl(mu), 0, 1)[0]
        fac[l] = (2.0*l + 1.0)*integ

    ps = {'k': a[:, 0]}
    
    ps['P0'] = (fac[0]*growth**2)*a[:, 1]
    ps['P2'] = (fac[2]*growth**2)*a[:, 1]
    ps['P4'] = (fac[4]*growth**2)*a[:, 1]

    return ps
    
        

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('zbin', help='zbin name in parameter file')
    parser.add_argument('--param', default='../data/param.json', help='parameter file')
    parser.add_argument('--print-ps', action='store_true', help='print power spectrum multipoles and exit')
    parser.add_argument('--print-xi', action='store_true', help='print xi multipoles and exit')
    parser.add_argument('--ofilename', '-o', default='corr.npz', help='print xi multipoles and exit')

    arg = parser.parse_args()

    with open(arg.param) as f:
        d = json.load(f)

    filename = d['linear']
    param = d[arg.zbin] 


    ps = compute_power_spectrum(filename,
                                param['D'], param['b'], param['f'])

    if arg.print_ps:
        ks = ps['k']
        P0 = ps['P0']
        P2 = ps['P2']
        P4 = ps['P4']
        for i, k in enumerate(ks):
            print('%e %e %e %e' % (k, P0[i], P2[i], P4[i]))
        sys.exit()

    xis = compute_corr(ps)

    if arg.print_xi:
        rs = xis['r']
        xi0 = xis['xi0']
        xi2 = xis['xi2']
        xi4 = xis['xi2']

        for i, r in enumerate(rs):
            print('%e %e %e %e' % (r, xi0[i], xi2[i], xi4[i]))

        sys.exit()

    np.savez(arg.ofilename, r=xis['r'],
             xi0=xis['xi0'], xi2=xis['xi2'], xi4=xis['xi4'])
                                                                    
    print(arg.ofilename, 'written')

if __name__ == "__main__":
    main()

