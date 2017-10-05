import math
import numpy as np
from scipy.special import sici, spherical_jn
from scipy.integrate import quad
from scipy.interpolate import CubicSpline

# Reference
# arxiv:1703.06428
# Indefinite Integrals of Spherical Bessel Functions
# Bloomfield, Face and Moss (2017)

class SphericalBesselIntegral:
    """
    Arg:
        x (np.array): array of x, sorted in increasing order

    Methods:
        
    """
    def __init__(self, x):
        assert(np.all(x == np.sort(x)))
        assert(x[0] >= 0.0)
        self.x = np.copy(x)
        
        self._x1 = x[x < 1.0]
        self._x2 = x[x >= 1.0]
        self._idx2 = x >= 1.0
        
        self.si, self.ci = sici(x)
        self.cache = {}
        self.sinx = None
        self.cosx = None
        self.sin_integ_cache = {}
        self.cos_integ_cache = {}
        self.jl_cache = {}

    def sin(self):
        if self.sinx is None:
            self.sinx = np.sin(self.x)
        return self.sinx

    def cos(self):
        if self.cosx is None:
            self.cosx = np.cos(self.x)
        return self.cosx
    
    def sin_integ(self, n):
        """
        \int^x[i] x^n sin(x)
        n (int): n >= -1
        x (array)
        """

        assert(type(n) == int)
        
        if n == 0:
            return -self.cos()
        elif n == -1:
            return self.si
        elif n in self.sin_integ_cache:
            return self.sin_integ_cache[n]
        
        if n > 0:
            integ = n*self.cos_integ(n - 1) - self.x**n*self.cos()
        else:
            integ = 1.0/(n + 1)*(self.x**(n + 1)*self.sin() - self.cos_integ(n + 1))

        self.sin_integ_cache[n] = integ
        
        return integ

    def cos_integ(self, n):
        """
        \int^x[i] x^n cos(x)
        """

        assert(type(n) == int)
        #assert(n >= -1)

        if n == 0:
            return self.sin()
        elif n == -1:
            return self.ci
        elif n in self.cos_integ_cache:
            return self.cos_integ_cache[n]

        if n > 0:
            integ = self.x**n*self.sin() - n*self.sin_integ(n - 1)
        else:
            integ = 1.0/(n + 1)*(self.x**(n + 1)*self.cos() + self.sin_integ(n + 1))

        self.cos_integ_cache[n] = integ
        
        return integ

    def spherical_bessel(self, l):
        assert(l >= 0)
        if l in self.jl_cache:
            jl = self.jl_cache[l]
        else:
            jl = spherical_jn(l, self.x)
            self.jl_cache[l] = jl
            
        return jl

    def _integ_xn_recursive(self, n, l):
        if l == 0:
            result = self.sin_integ(n - 1)
            self.cache[(n, l)] = result
            return result

        result = ((l + n - 1)*self._integ_xn_recursive(n - 1, l - 1) -
                  self.x**n*self.spherical_bessel(l - 1))
        return result
        
    def jl_integ(self, n, l):
        """
        \int_^x[i] x^n j_l(x)  dx 
        
        n (int): n >= 0 
        l (int): l >= 0
        x (array):

        Returns 
        array of integrals to x[i]; lowerbound is arbitaray
        """

        if (n, l) in self.cache:
            return self.cache[(n, l)]

        assert(l >= 0)
        
        if n < l:
            # x < 1.0
            # recursive formula is not good because some sin cos terms
            # are divergent at r = 0

            result = np.zeros_like(self.x)
            
            for i, xi in enumerate(self._x1):
                result[i] = quad(lambda x: x**n*spherical_jn(l, x), 0, xi)[0]

            # x > 1.0
            result2 = self._integ_xn_recursive(n, l)
            
            # match the integration constant\
            n1= len(self._x1)
            x0 = self.x[n1]
            y0 = quad(lambda x: x**n*spherical_jn(l, x), 0, x0)[0]
            result2 += y0 - result2[n1]

            result[n1:] = result2[n1:]        
        else:
            ## recursive method is also not accurate for very small x
            n1 = np.argmin(self.x**(n + l) < 1.0e-2)

            if n1 > 0:
                #print('small x', self.x[:n1])
                #print('small x', n1)
                result = np.zeros(len(self.x))
                result2 = self._integ_xn_recursive(n, l)

                xsmall = self.x[:n1]
                #for i in range(n1):
                #    result[i] = quad(lambda x: x**n*spherical_jn(l, x), 0, self.x[i])[0]
                # j_l(x) ~ fac1*x**l + fac2*x**(l+2)
                fac1 = 2**l*(math.factorial(l)/math.factorial(2*l + 1))
                fac2 = -2**l*(math.factorial(l + 1)/math.factorial(2*l + 3))
                #fac3 = 2**l*(math.factorial(l + 2)/math.factorial(2*l + 5))
                
                result[:n1] = (fac1/(n + l + 1)*xsmall**(n + l + 1)
                               + fac2/(n + l + 3)*xsmall**(n + l + 3))


                # match integration constant
                x0 = self.x[n1]
                y0 = quad(lambda x: x**n*spherical_jn(l, x), 0, x0)[0]
                result2 += y0 - result2[n1]
                result[n1:] = result2[n1:]
            else:
                result = self._integ_xn_recursive(n,l)
            

        self.cache[(n, l)] = result
        
        return result

    def integrate(self, f, n, l, *, kind='linear'):
        if kind == 'linear':
            return self.integrate_linear(f, n, l)
        elif kind == 'cubic':
            return self.integrate_spline(f, n, l)
        else:
            raise TypeError('Unknown kind: ' + kind)

        
    def integrate_linear(self, f, n, l):
        """
        \integ x^n f(x) j_l(x) dx

        f(x) is a piecewise-linear function interpolating the input array f

        f (array)
        n (int): n >= 0
        l (int): l >= 0


        """

        assert(len(f) == len(self.x))

        # linear interpolation
        # piece-wise f(x) = a*x + b
        dx = self.x[1:] - self.x[:-1]
        a = (f[1:] - f[:-1])/dx
        b = (f[:-1]*self.x[1:] - f[1:]*self.x[:-1])/dx

        # x^n f(x) = a x^n+1  + b x^n
        a_integ = self.jl_integ(n + 1, l)
        b_integ = self.jl_integ(n, l)

        return np.sum(a*(a_integ[1:] - a_integ[:-1]) + b*(b_integ[1:] - b_integ[:-1]))


    def integrate_spline(self, f, n, l):
        """
        \integ x^n f(x) j_l(alpha*x) dx

        f(x) is a cubic polynomial interpolating the input array f

        f (array)
        n (int): n >= 0
        l (int): l >= 0


        """

        assert(len(f) == len(self.x))

        #
        # Cubic spline interpolation
        #
        spline = CubicSpline(self.x, f)
        
        # convert the coeffecient of c_n (x - x_0)^n to a_n x^n
        c = spline.c
        x0 = self.x[:-1]
        x0_2 = x0**2

        a = {}
        a[0] = c[3, :] - c[2, :]*x0 + c[1, :]*x0_2 - c[0, :]*x0**3
        a[1] = c[2, :] - 2.0*c[1, :]*x0 + 3.0*c[0, :]*x0_2
        a[2] = c[1, :] - 3.0*c[0, :]*x0
        a[3] = c[0, :]

        total = 0.0

        # \int^x_i x^(n+i) j_l(x) dx
        for k in range(4):
            integ = self.jl_integ(n + k, l)
            integ = integ[1:] - integ[:-1]

            #print(a[k][:20])
            #print(integ[:20])
            # remove possible numerical precision limit
            # ???
            #idx = np.absolute(integ) < 1.0e-13
            #integ[idx] = 0.0
            total += np.sum(a[k]*integ)
            #print('total', k, total)
#            a1_integ = self.jl_integ(n + 1, l)
#            a2_integ = self.jl_integ(n + 2, l)
#            a3_integ = self.jl_integ(n + 3, l)

        # integ[i] = sum_i \int_x(i)^x(i+1) a_k(i) x^(n+k) j_l(x)
        #i#nteg =  a0*(a0_integ[1:] - a0_integ[:-1])
        #print('integ', np.sum(integ))

        #integ += a1*(a1_integ[1:] - a1_integ[:-1])
        #print('integ', np.sum(integ))

        #integ += a2*(a2_integ[1:] - a2_integ[:-1])
        #print('integ', np.sum(integ))
                
        #integ += a3*(a3_integ[1:] - a3_integ[:-1])
        #print((a3[:20]))
        #print((a3_integ[1:] - a3_integ[:-1])[:20])
        #print((a3*(a3_integ[1:] - a3_integ[:-1]))[:20])
        #print(np.max(a3*(a3_integ[1:] - a3_integ[:-1])))
        #print(np.argmax(a3*(a3_integ[1:] - a3_integ[:-1])))
        #print('integ', np.sum(integ))
                

        return total #np.sum(integ)


def test_sin_integ(sb, x):
    import math


    for n in [-1, 0, 1, 2]:
        y = sb.sin_integ(n)
        y -= y[0]
        errs = []
        for i, xi in enumerate(x):
            if xi == 0.0 and n < 0:
                continue
            
            y_check = quad(lambda x: x**n*math.sin(x), 0, xi)[0]
            e = y[i] - y_check
            errs.append(e)
            assert(abs(e) < 1.0e-8)

        print("sin_integ n=%2d OK: %e" % (n, np.std(errs)))

def test_cos_integ(sb, x):
    import math
    from scipy.integrate import quad

    for n in [0, 1, 2]:
        y = sb.cos_integ(n)
        y -= y[0]
        errs = []
        for i, xi in enumerate(x):
            y_check = quad(lambda x: x**n*math.cos(x), 0, xi)[0]
            e = y[i] - y_check
            errs.append(e)
            assert(abs(e) < 1.0e-8)

        print("cos_integ n=%2d OK: %e" % (n, np.std(errs)))

def test_bessel_integ(sb, x):
    import math
    from scipy.integrate import quad

    l=0

    for l in [0, 2, 4]:
        for n in [0, 1, 2]:
            y = sb.jl_integ(n, l)
            y -= y[0]
            errs = []
            for i, xi in enumerate(x):
                #print('xi', xi)
                y_check = quad(lambda x: x**n*spherical_jn(l, x), 0, xi)[0]
                e = y[i] - y_check
                errs.append(e)
                if abs(e) >= 1.0e-8:
                    print('Error', e)
                    print('x**%d j_%d(x=%e) = %e, %e' %
                          (n, l, xi, y[i], y_check))

                assert(abs(e) < 1.0e-8)

            print("jnl_integ n=%d l=%d OK: %e" % (n, l, np.std(errs)))

def test_integrate(s):
    a = 2.0
    b = 1.0
    x = np.linspace(0.01, 10, 100)
    sb = SphericalBesselIntegral(x)
    f = a*x + b

    errs = []

    for l in [0, 2, 4]:
        for n in [0, 1, 2]:
            y = sb.integrate(f, n, l)
            y_check = quad(lambda x: (a*x + b)*x**n*spherical_jn(l, x), x[0], x[-1])[0]
            e = y - y_check
            errs.append(e)
            if abs(e) >= 1.0e-8:
                print('Error', e)
                print('f x**n*j_l(x); n=%d l=%d' % (n, l))

            assert(abs(e) < 1.0e-8)
    print("test_integrate OK: %e" % np.std(errs))

def test_integrate_spline(s):
    a0 = 1.0
    a1 = 2.0
    a2 = 3.0
    a3 = 4
    x = np.linspace(0.01, 10, 100)
    sb = SphericalBesselIntegral(x)
    f = a0 + a1*x + a2*x**2 + a3*x**3

    errs = []

    for l in [0, 2, 4]:
        for n in [0, 1, 2]:
            y = sb.integrate_spline(f, n, l)
            y_check = quad(lambda x: (a0 + a1*x + a2*x**2 + a3*x**3)*x**n*spherical_jn(l, x), x[0], x[-1])[0]
            e = y - y_check
            errs.append(e)
            if abs(e) >= 1.0e-8:
                print('Error', e)
                print('f x**n*j_l(x); n=%d l=%d' % (n, l))

            assert(abs(e) < 1.0e-8)
    print("test_integrate OK: %e" % np.std(errs))

    
def main():
    """Test SphericalBesselIntegral"""
    x = np.linspace(0.0, 10, 101)
    sb = SphericalBesselIntegral(x)

    #test_sin_integ(sb, x)
    #test_cos_integ(sb, x)
    #test_bessel_integ(sb, x)
    #test_integrate(sb)
    test_integrate_spline(sb)
    
if __name__ == "__main__":
    main()

