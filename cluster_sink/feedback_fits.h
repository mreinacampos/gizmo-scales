/* Coefficients for the analytical fits of different feedback mechanisms rates, masses, energies and yields
 * as described in https://ui.adsabs.harvard.edu/abs/2022arXiv220300040H/abstract
 * This file was written by Marta Reina-Campos (reinacampos@mcmaster.ca) for GIZMO.
 */

#ifdef CLUSTER_SINK

#ifdef CLUSTER_SINK_SNII
// Coefficients as,1 , as,2  and as,3 - in Gyr^-1 MSun^-1
static double SNII_coeff_asj[3] =  {0.39, 0.51, 0.18}; 
// Timescales ts,1 , ts,2  and ts,3 - in Myr
static double SNII_tsj[3] =  {3.7, 7.0, 44}; 
#endif

#endif // CLUSTER_SINK