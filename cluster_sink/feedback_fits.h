/* Coefficients for the analytical fits of different feedback mechanisms rates, masses, energies and yields
 * as described in https://ui.adsabs.harvard.edu/abs/2022arXiv220300040H/abstract
 * This file was written by Marta Reina-Campos (reinacampos@mcmaster.ca) for GIZMO.
 */

#ifdef CLUSTER_SINK

// core-collapse SNe FB
#ifdef CLUSTER_SINK_SNII
// coefficients for the rate of core-collapse SNe
// Coefficients as,1 , as,2  and as,3 - in Gyr^-1 MSun^-1
static double SNII_coeff_asj[3] =  {0.39, 0.51, 0.18}; 
// Timescales ts,1 , ts,2  and ts,3 - in Myr
static double SNII_tsj[3] =  {3.7, 7.0, 44}; 
#endif

// SNIa FB
#ifdef CLUSTER_SINK_SNIa
// coefficients for the rate of core-collapse SNIa
// Coefficient aIa,1 - in Gyr^-1 MSun^-1
static double SNIa_coeff_aj[1] =  {0.0083}; 
// Timescale tIa,1 - in Myr
static double SNIa_tj[1] =  {44}; 
// Slope lambdaIa,1 - dimensionless
static double SNIa_slope[1] =  {-1.1}; 
#endif

// AGB&OB winds FB
#ifdef CLUSTER_SINK_WINDS
// non-z dependent coefficients for the mass-loss from winds
// Coefficients aa,1 and aa,2 - in Gyr^-1, except aa,2
static double WINDS_coeff_aaj[2] =  {0.01, 0.01}; 
// Timescales tw,1 , tw,2, tw,3 and ta - in Myr
static double WINDS_twj[4] =  {1.7, 4.0, 20, 1000}; 
// timescales for the velocity of injection
// Timescales tv,1 , and tv,2 - in Myr
static double WINDS_tvj[2] =  {3.0, 50.0}; 
#endif

#endif // CLUSTER_SINK