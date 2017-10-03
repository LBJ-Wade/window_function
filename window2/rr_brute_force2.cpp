//
// Computes window function in Beutler et al (2014) arXiv:1312.4611v2
//

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <boost/program_options.hpp>


#include "read.h"
#include "histogram2.h"

static double rmax2= 0.0;
#include "inline_functions.h"

using namespace std;
using namespace std::chrono;
using namespace boost::program_options;
typedef time_point<high_resolution_clock> t;


template <class F>
void count_pairs_multipole(const vector<Particle>& v,
			   Histogram2<F>& hist)
{  
  double r, mu;
  for(vector<Particle>::const_iterator p= v.begin(); p != v.end(); ++p) {
   for(vector<Particle>::const_iterator q= p+1; q != v.end(); ++q) {
     if(r_mu(p->x, q->x, r, mu))
       hist.add(r, mu, p->w*q->w);
   }
  }
}


int main(int argc, char* argv[])
{
  //
  // command-line options (Boost program_options)
  //
  options_description opt("corr_brute_force [options] filename");
  opt.add_options()
    ("help,h", "display this help")
    ("filename,f", value<string>(), "filename")
    ("rand-factor", value<double>()->default_value(50.0),
     "number of randoms / data")
    ("rmin", value<double>()->default_value(0.0f), "r minimum")
    ("rmax", value<double>()->default_value(100.0f), "r maximum")
    ("nbin", value<int>()->default_value(10), "nbin")
    ("logbin", "use logarithmic binning")
    ("subsample", value<double>()->default_value(1.0), "subsample rate")
    ("seed", value<unsigned long>()->default_value(1), "random subsample seed")
    ("Pest", value<double>()->default_value(20000.0), "Pest for FKP weight")
    ("ofilename,o", value<string>()->default_value("rr.txt"), "output RR file name")
    ;
  
  positional_options_description p;
  p.add("filename", -1);
  
  variables_map vm;
  store(command_line_parser(argc, argv).options(opt).positional(p).run(), vm);
  notify(vm);

  if(vm.count("help") || ! vm.count("filename")) {
    cout << opt; 
    return 0;
  }

  const string filename= vm["filename"].as<string>();
  //
  // Read particle/halo file
  //
  
  vector<Particle> v;
  const unsigned long seed= vm["seed"].as<unsigned long>();
  const double subsample_rate= vm["subsample"].as<double>();
  const double Pest= vm["Pest"].as<double>();

  t ts= std::chrono::high_resolution_clock::now();

  read_ascii_fkp(filename.c_str(), Pest, subsample_rate, seed, v);

  t te= std::chrono::high_resolution_clock::now();
  double dt= std::chrono::duration<double>(te - ts).count();
  cout << "Time reading: " << dt << endl;
    
  cout << v.size() << " particles read.\n";
  assert(v.size() > 0);

  // one-point statistices
  double num_factor= vm["rand-factor"].as<double>()*subsample_rate;
  // correct the missmatch of the number density of the points / number density in p.nbar
  long double w_sum= 0.0;
  long double w2_sum= 0.0;
  long double nw2_sum= 0.0;
  
  for(vector<Particle>::const_iterator p= v.begin(); p != v.end(); ++p) {
    w_sum += p->w;
    w2_sum += p->w*p->w;
    nw2_sum += p->nbar*p->w*p->w;
  }

  printf("w_sum %Le\n", w_sum);
  printf("nw2_sum %Le\n", nw2_sum);

  //
  
  // pair counts
  const double rmin= vm["rmin"].as<double>();
  const double rmax= vm["rmax"].as<double>();
  rmax2= rmax*rmax;
  const int nbin= vm["nbin"].as<int>(); assert(nbin > 0);

  string ofilename= vm["ofilename"].as<string>();

  Histogram2<Linear> hist(rmin, rmax, nbin, Linear());
  hist.w_sum= w_sum;
  hist.w2_sum= w2_sum;
  hist.nw2_sum= num_factor*nw2_sum;

  

  ts= std::chrono::high_resolution_clock::now();

  count_pairs_multipole(v, hist);
  te= std::chrono::high_resolution_clock::now();

  write_hist(ofilename.c_str(), hist);

  dt= std::chrono::duration<double>(te - ts).count();
  cout << "Time pair count: " << dt << endl;
  
  return 0;
}
