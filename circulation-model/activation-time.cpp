#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include "unistd.h"
#include "getopt.h"
#include "SimpleRNG.cpp"
#include "INIReader.h"

using namespace std;

SimpleRNG rng;

inline double sampler_01(){
	return rng.GetUniform();
}

inline double sampler_norm(){
	return rng.GetNormal( 0.0, sqrt(2.0) );
}

#include "jacobi.h"
#include "params.h"

const int IN_BLOOD = 0, IN_SPLEEN = 1, IN_LN = 2, IN_DLN = 3;

double mod (double a, int b)
{
   int ret = fmod(a,b);
   if(ret < 0)
     ret+=b;
   return ret;
}

double blood_adjust( double t ){
	double x = mod(t,24);
	if( x < 7 || x > 22. ){
		return 2.0;
	} else {
		return 0.5;
	}
	//return 1.0; //mod(t,24)*(1.5+sin(2*PI*t/24));
}

double nu_t( double t ){
   double ymin = RHO_INITIAL;
   double ymax = RHO_FINAL;
   double tmin = 0.;
   double tmax = 4.5*24.;
   if( t <= tmin ){
      return ymin;
   } 
   if( t >= tmax ){
      return ymax;
   }
   double y = (t-tmin) / (tmax-tmin) * (ymax-ymin);
   return ymin + y;
}

inline double s01(){
   return rng.GetUniform();
}

inline double activation_sampler(){
    return rng.GetGamma( GAMMA_K, 1.0 );
}

void track_position( int state, int ** cell_distribution, int & t_discrete, double t ){
	if( OUTPUT_DISTRIBUTION ){
		while( t_discrete < t ){
			if( t_discrete >= 0 && t_discrete < T_MAX ){
				cell_distribution[state][t_discrete]++;
			}
			t_discrete ++;
		}
	}
}

typedef double (*smp_t)(void);

int main( int argc, char ** argv ){

	INIReader reader("param.ini");
    if (reader.ParseError() < 0) {
        std::cout << "Can't load 'test.ini'\n";
        return 1;
    }
    N = reader.GetInteger("","N",N);
    CUTOFF = reader.GetInteger("","CUTOFF",CUTOFF);
	T_MAX = reader.GetInteger("","T_MAX",T_MAX);
    OUTPUT_DISTRIBUTION = reader.GetBoolean("","OUTPUT_DISTRIBUTION",OUTPUT_DISTRIBUTION);
    OUTPUT_ACTIVATION = reader.GetBoolean("","OUTPUT_ACTIVATION",OUTPUT_ACTIVATION);
    ALPHA = reader.GetReal("","ALPHA",ALPHA);
    USE_BLOOD_ADJUST = reader.GetBoolean("","USE_BLOOD_ADJUST",USE_BLOOD_ADJUST);


	if( SEED < 0 ){
		rng.SetState( time( NULL ), getpid() );
	} else {
		rng.SetState( SEED, 12345 );
	}

   int ** cell_distribution;
	
	if( OUTPUT_DISTRIBUTION ){
		cell_distribution = new int*[5];
		for( int i = 0 ; i < 5 ; i ++ ){
			int tmax_min = (int)ceil( T_MAX )+1;
			cell_distribution[i] = new int[ tmax_min ];
			memset( cell_distribution[i], 0, tmax_min * sizeof(int) );
		}
	}

   for( int i = 0 ; i < N ; i ++ ){

      int state = IN_BLOOD, pre_state = IN_BLOOD, activated = 0;

      double t = (int) -(T_PRE*(1+s01()));

		if( START_IN_LN ){
			if( sampler_01() < nu_t( 0 ) / RHO_LYMPH ){
				state = IN_DLN; 
			} else {
				state = IN_LN;
			}
			pre_state = IN_BLOOD; t = 0.0;
		}
		
      int t_discrete = (int) t;

      double pre_t = 0.0, ta;

      while( !CUTOFF || t <= T_MAX ){

         pre_state = state;
         pre_t = t;
         double delta_t=0.0,ts;

         double rho_n, rho_dln, rho_ndln; 

         switch( state ){
         case IN_BLOOD :

               // perform kinetic monte carlo
               rho_dln = nu_t( USE_SWELLING ? t : 0 );
               rho_ndln = RHO_LYMPH - nu_t(0);

               // determine next reaction proportional to current rates
               ts = sampler_01()*(RHO_SPLEEN + rho_dln + rho_ndln);
               if( ts < RHO_SPLEEN ){
                  state = IN_SPLEEN;
               } else if( ts < RHO_SPLEEN+rho_ndln ){
                  state = IN_LN;
               } else {
                  state = IN_DLN;
               }

               // determine time step
               delta_t = log(1./s01())/(RHO_SPLEEN + rho_dln + rho_ndln);
               if( USE_BLOOD_ADJUST ) delta_t /= blood_adjust(t); 

               break;

         case IN_LN :
               // nothing ever happens in non-dLNs
              delta_t = CONSTANT_LN_TRANSIT_TIME ? R_LYMPH : inv_phi( sampler_01(), L_LYMPH );
              state = IN_BLOOD;

               break;


         case IN_SPLEEN :

               //delta_t = inv_phi( sampler_01(), L_SPLEEN );
               delta_t = inv_phi_aperture( L_SPLEEN/5.167, 0.01, M_PI/9. );

               // move to blood by default ...

               // ... but check for activation
               if( SYSTEMIC && ALPHA >= 0 && t + delta_t > 0 ){
                  ta = activation_sampler() * ALPHA;
                  if( ta < min( delta_t, t + delta_t ) ){
                     if( OUTPUT_ACTIVATION ) cout << (max(t,0.)+ta)/24. << endl;
                     activated = 1;
                  }
               }

               state = IN_BLOOD;

               break;

         case IN_DLN :

               if( CONSTANT_LN_TRANSIT_TIME ){
                  delta_t = R_LYMPH;
               }
               else{
						if( DLN_SHUTDOWN ){
							delta_t = inv_phi_explicit( L_LYMPH, RW_DELTA_T, DLN_SHUTDOWN_FROM-t, DLN_SHUTDOWN_TO-t  );
						} else {
							delta_t = inv_phi( sampler_01(), L_LYMPH );
						}
               }

               // ... but check for activation
               if( !SYSTEMIC && ALPHA >= 0 && t + delta_t > 0 ){
                  ta = activation_sampler() * ALPHA;
                  if( ta < min( delta_t, t + delta_t ) ){
                     if( OUTPUT_ACTIVATION ) cout << (max(t,0.)+ta)/24. << endl;
                     activated = 1;
                  }
               }

               state = IN_BLOOD;

               break;
         }
			
			t = t + delta_t;
			track_position( pre_state, cell_distribution, t_discrete, t );

			if( OUTPUT_TRAJECTORY ){
				if( !activated ){
					cout << i << " " << t << " " << state << endl;
				} else {
					cout << i << " " << max(t-delta_t,0.)+ta << " A" << endl;
				}
			}

         if( activated ){
            break;
         }
      }
   }
	if( OUTPUT_DISTRIBUTION ){
		for( int i = 0 ; i < T_MAX ; i ++ ){
			cout << i/24. << " " << cell_distribution[IN_BLOOD][i] << " " 
				<< cell_distribution[IN_SPLEEN][i] << " " 
				<< cell_distribution[IN_LN][i] << " " 
				<< cell_distribution[IN_DLN][i] << endl;
		}
	}
}
