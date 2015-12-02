#include <limits>

/* This file contains implementations of the residence time 
 * distributions for lymph nodes and spleen. */

const double PISQ = M_PI * M_PI;

const double EPS = 1000.0 * numeric_limits<double>::epsilon();

double phi( double t, double L ){
   double r = 0;
   int sig = -1;
   double n = 1;

   double a = t * PISQ / L / L;
   double rs;

   if( exp( -n*n*a ) > 1-EPS ){
      return 1;
   }

   do{
      rs = - 2. * exp( - n * n * a ) * sig;
      r += rs;
      sig = -sig;
      n = n + 1;
   } while( abs(rs) > EPS );

   return r;
}

double inv_phi( double y, double L ){
   double lower = 0.0;
   double upper = 1.0;
   while( 1-phi( upper, L ) < y ){
      upper *= 2;
   }
   
   double center = upper/2.;

   while( (upper-lower) > EPS ){
      center = lower + (upper-lower)/2.;
      if( 1-phi( center, L ) > y ){
         upper = center;
      } else {
         lower = center;
      }
   }

   return center;
}

double inv_phi_explicit( double L, double delta_t, 
	double shutdown_min, double shutdown_max ){
	double dx, dy, dz;
	double x1 = 0.0;
	double x2 = 0.0;
	double x3 = 0.0;
	double t = 0.0;
	double L2 = L*L;
	while( x1*x1 + x2*x2 + x3*x3 < L2 ){
		dx = sqrt(delta_t)*sampler_norm();
		dy = sqrt(delta_t)*sampler_norm();
		dz = sqrt(delta_t)*sampler_norm();
		x1 += dx; x2 += dy; x3 += dz;
		if(   t >= shutdown_min && t < shutdown_max 
			&& x1*x1 + x2*x2 + x3*x3 >= L2 ){
			x1 -= dx; x2 -= dy; x3 -= dz;
		}
		t += delta_t; 
	}
	return t;
}

double inv_phi_aperture( double L, double delta_t, double aperture_radius ){
	double dx, dy, dz;
	double cos_aperture_radius = cos(aperture_radius);
	double L_cos_aperture_radius = L*cos_aperture_radius;
	double x1 = -L_cos_aperture_radius;
	double x2 = 0.0;
	double x3 = 0.0;
	double t = 0.0;
	double L2 = L*L;
	double norm2x = 0.0;
	double norm2x_pre = 0.0;
	double sqrt_delta_t = sqrt(delta_t);
	while( x1 < L_cos_aperture_radius ){
		dx = sqrt_delta_t*sampler_norm();
		dy = sqrt_delta_t*sampler_norm();
		dz = sqrt_delta_t*sampler_norm();
		x1 += dx; x2 += dy; x3 += dz;
		norm2x_pre = norm2x;
		norm2x = x1*x1 + x2*x2 + x3*x3;
		if( norm2x >= L2 ){
			if( x1/sqrt(norm2x) < cos_aperture_radius ){
				x1 -= dx; x2 -= dy; x3 -= dz;
				norm2x = norm2x_pre;
			}
		}
		t += delta_t; 
	}
	return t;
}
