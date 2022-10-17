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
// yields ejected: He, C, N, O, Ne, Mg, Si, S, Ca, Fe
// coefficients acc,j - Table 1 in Hopkins+22
static double SNII_yields_accj[NUM_METAL_SPECIES-1][5] = { {0.461, 0.330, 0.358, 0.365, 0.359},
                                                    {0.237, 8.57e-3, 1.69e-2, 9.33e-3, 4.47e-3},
                                                    {1.07e-2, 3.48e-3, 3.44e-4, 3.72e-3, 3.50e-3},
                                                    {9.53e-2, 0.102, 9.85e-2, 1.73e-2, 8.2e-3},
                                                    {2.60e-2, 2.20e-2, 1.93e-2, 2.70e-3, 2.75e-3},
                                                    {2.89e-2, 1.25e-2, 5.77e-3, 1.03e-3, 1.03e-3},
                                                    {4.12e-4, 7.69e-3, 8.73e-3, 2.23e-3, 1.18e-3},
                                                    {3.63e-4, 5.61e-3, 5.49e-3, 1.26e-3, 5.75e-4},
                                                    {4.28e-5, 3.21e-4, 6.00e-4, 1.84e-4, 9.64e-5},
                                                    {5.46e-4, 2.18e-3, 1.08e-2, 4.57e-3, 1.83e-3} };
// timescales tcc,j - in Myr - Table 1 in Hopkins+22
static double SNII_yields_tccj[5] = {3.7, 8, 18, 30, 44};
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
// yields ejected: Z, He, C, N, O, Ne, Mg, Si, S, Ca, Fe
static double SNIa_yields[NUM_METAL_SPECIES] = {1, 0, 1.76e-2, 2.1e-6, 7.36e-2, 2.02e-3, 6.21e-3, 0.146, 7.62e-2, 1.29e-2, 0.558};
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