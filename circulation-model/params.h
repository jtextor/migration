double RHO_INITIAL = 1./39.;

double RHO_FINAL = 16./39.;

const int IN_BLOOD = 0, IN_SPLEEN = 1, IN_LN = 2, IN_DLN = 3;

double R_SPLEEN = 6; // hours
double L_SPLEEN = sqrt( 6 * R_SPLEEN );

// multiply these L-values by sqrt( 6000 ) to get the radius in 
// microns for a motility coeff. of 6000 um²/h = 100 um²/min 

double R_LYMPH = 13.5; // hours
double L_LYMPH = sqrt( 6 * R_LYMPH );

double RHO_SPLEEN = 1.0, RHO_LYMPH = 1.5;

double T_PRE = 7 * 24, T_MAX = 7.5 * 24;

int N = 10;
int N_MULT = 10;

// only relevant for shutdown (explicit transit simulation)

double RW_DELTA_T = .01;

int SYSTEMIC = 0;

int CUTOFF = 0;

//double ALPHA = 0.125;
//double GAMMA_K = 1.0;

double ALPHA=0.125;
double GAMMA_K=1.;

int USE_SWELLING = 0;

int OUTPUT_ACTIVATION = 1;
int OUTPUT_DISTRIBUTION = 0;
int OUTPUT_TRAJECTORY = 0;

int CONSTANT_LN_TRANSIT_TIME = 0;

int DLN_SHUTDOWN = 0;

double DLN_SHUTDOWN_FROM = 0.0;

double DLN_SHUTDOWN_TO = 18.0;


// do not simulate pre-circulation to attain
// equilibrium. Instead, start in a LN in all
// cases.
int START_IN_LN = 0;

// set the seed for the PRNG if necessary to
// keep illustration for figures constant
long SEED = -1;
