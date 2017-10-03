#include <iostream>
#include <cstdio>
#include <cassert>
#include <gsl/gsl_rng.h>
#include "read.h"


using namespace std;

void read_pos_ascii_file(const char filename[], vector<Particle>& v)
{
  // Read ASCII File x y z
  
  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Error: file not found: " << filename << endl;
    throw filename;
  }

  Particle p;
  p.w= 1.0;
  char buf[256];
  while(fgets(buf, 255, fp)) {
    if(buf[0] == '#')
      continue;

    int ret= sscanf(buf, "%le %le %le", p.x, p.x+1, p.x+2);
    assert(ret == 3);

    v.push_back(p);
  }
  
  fclose(fp);
}

void read_ascii_fkp(const char filename[], const double Pest,
		    const double sample_rate, const unsigned long seed,
		    vector<Particle>& v)
{
  // Read ASCII File x y z nbar

  gsl_rng* r = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(r, seed);
  
  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Error: file not found: " << filename << endl;
    throw filename;
  }

  Particle p;
  char buf[256];
  while(fgets(buf, 255, fp)) {
    if(buf[0] == '#')
      continue;

    int ret= sscanf(buf, "%le %le %le %le", p.x, p.x+1, p.x+2, &p.nbar);
    assert(ret == 4);

    p.w = 1.0/(1.0 + p.nbar*Pest);

    if(gsl_rng_uniform(r) < sample_rate)
      v.push_back(p);
  }

  fclose(fp);
  gsl_rng_free(r);
    
}


