#ifndef HIST_H
#define HIST_H

#include <vector>
#include <cassert>
#include <cmath>

struct Multipole2 {
  double r, y00, y22, y44, y02, y04, y24;
};

struct Linear {
  double operator()(const double x) const { return x; }
  double inv(const double x) const { return x; }
};

struct Log {
  double operator()(const double x) const { assert(x>0); return log(x); }
  double inv(const double x) const { return exp(x); }
};

template <class F> class Histogram2 {
 public:
  Histogram2(const double rmin, const double rmax,
		     const int nbin, const F f);

  void add(const double r, const double mu, const double w) {
    int i= (int)((f(r) - fxmin)/(fxmax - fxmin)*nbin);
    if(0 <= i && i < nbin) {
      Multipole2& m= h[i];

      double mu2= mu*mu;
      double l_2= 7.5*mu2 - 2.5;
      double mu4= mu2*mu2;
      double l_4= 1.125*(35.0*mu2*mu2 - 30.0*mu2 + 3.0);
      m.r   += w*r;
      m.y00 += w;
      m.y22 += w*l_2*l_2;
      m.y44 += w*l_4*l_4;
      m.y02 += w*l_2;
      m.y04 += w*l_4;
      m.y24 += w*l_2*l_4;
    }
  }

  Multipole2 operator[](const int i) const {
    return h[i];
  }
  double x_left(const int i) const {
    return f.inv(fxmin + (fxmax - fxmin)*i/nbin);
  }
  double x_right(const int i) const {
    return f.inv(fxmin + (fxmax - fxmin)*(i+1)/nbin);
  }
  int size() const {
    return nbin;
  }
  double x_min() const { return f.inv(fxmin); }
  double x_max() const { return f.inv(fxmax); }

  double w_sum, w2_sum, nw2_sum;
  
 private:
  double fxmin, fxmax;
  int nbin;
  std::vector<Multipole2> h;
  
  F f;
};


template <class F>
Histogram2<F>::Histogram2(const double xmin, const double xmax,
					 const int nbin_, const F func):
 fxmin(func(xmin)), fxmax(func(xmax)), nbin(nbin_), f(func)
{
  Multipole2 m;
  m.r= m.y00= m.y22= m.y44= m.y02= m.y04= m.y24= 0.0;
    
  h.assign(nbin, m);
}

void write_hist(const char filename[], const Histogram2<Linear>& hist);

#endif
