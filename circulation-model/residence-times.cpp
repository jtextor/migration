#include <iostream>
#include <cstdio>
#include <cstring>
#include <ctime>
#include "boost/random.hpp"
#include "boost/random/exponential_distribution.hpp"

using namespace std;

boost::mt19937 rng( time(NULL) + getpid() );
boost::variate_generator<boost::mt19937&, 
boost::normal_distribution< double > > 
sampler_norm(rng, boost::normal_distribution<double>(0.0,sqrt(2.0)));

#include "params.h"
#include "jacobi.h"

int main( int argc, char ** argv ){
	double delta_t = 0.01;
	for( double t = 0.0 ; t < 48.0 ; t = t + delta_t ){
		cout << t << " " << (phi( t, L_SPLEEN )-phi(t+delta_t,L_SPLEEN))/delta_t << 
			" " << " " << (phi( t, L_LYMPH )-phi(t+delta_t,L_LYMPH))/delta_t << endl;
	}
}
