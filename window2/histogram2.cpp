#include <iostream>
#include <cstdio>
#include "histogram2.h"

using namespace std;

void write_hist(const char filename[], const Histogram2<Linear>& hist)
{
  FILE* const fp= fopen(filename, "w");
  if(fp == 0) {
    cerr << "Unable to write to: " << filename << endl;
    throw filename;
  }
  
  for(int i=0; i<hist.size(); ++i) {
    double xleft= hist.x_left(i);
    double xright= hist.x_right(i);
    double vol= 4.0*M_PI/3.0*(xright*xright*xright - xleft*xleft*xleft);
    double fac= 2.0/(vol*hist.nw2_sum);
    // 0.5 for 1/2*N*(N - 1) pairs

    const Multipole2& m= hist[i];
    double r= m.r / m.y00;
    fprintf(fp, "%le %le %le %le %le %le %le\n", r,
	    fac*m.y00, fac*m.y22, fac*m.y44, fac*m.y02, fac*m.y04, fac*m.y24);
    
    // Column 1: r mean
    // Column 2: RR00
    // Column 3: RR22
    // Columm 4: RR44
    // Column 5: RR02
    // Column 6: RR04
    // Column 7: RR24
  }  
  fclose(fp);
}
