#ifndef _INLINE_FUNCTIONS_H
#define _INLINE_FUNCTIONS_H 1

static inline double dist1(const double left1, const double right1,
			  const double left2, const double right2)
{
  // Distance between two segments [left1,right1] and [left2,right2]
  if(right1 < left2)
    return left2 - right1;
  else if(right2 < left1)
    return left1 - right2;

  return 0.0;
}


static inline double sq(const double x[])
{
  return (double) x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
}


static inline double norm(const double x[])
{
  return sqrt((double) x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

static inline double dot(const double x[], const double y[])
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

static inline double dist(const double x[], const double y[])
{
  double dx= x[0] - y[0];
  double dy= x[1] - y[1];
  double dz= x[2] - y[2];
  return sqrt(dx*dx + dy*dy + dz*dz);
}

static inline double dist2(const double x[], const double y[])
{
  double dx= x[0] - y[0];
  double dy= x[1] - y[1];
  double dz= x[2] - y[2];
  return dx*dx + dy*dy + dz*dz;
}

static inline int r_mu(const double x[], const double y[],
		       double& r, double& mu)
{
  double dx[3];
  dx[0]= x[0] - y[0];
  dx[1]= x[1] - y[1];
  dx[2]= x[2] - y[2];

  double r2= dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
  if(r2 > rmax2)
    return false;


  // Line-of-sight unit vector \hat{(x+y)/2}
  double rhat[3];
  rhat[0]= x[0] + y[0];
  rhat[1]= x[1] + y[1];
  rhat[2]= x[2] + y[2];
  
  double sum_norm= norm(rhat);
  assert(sum_norm > 0.0); // DEBUG
  rhat[0] /= sum_norm;
  rhat[1] /= sum_norm;
  rhat[2] /= sum_norm;

  r= norm(dx);
  mu= dot(rhat, dx)/r;

  return true;
}

#endif
