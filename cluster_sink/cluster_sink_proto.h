
/*
 * Declarations of routines within the CLUSTER_SINK module
 * This file was written by Marta Reina-Campos (reinacampos@mcmaster.ca) for GIZMO.
 */

#ifdef CLUSTER_SINK
double determine_sne_rate(int i, double dt);

#ifdef CLUSTER_SINK_SNII
double determine_corecollapse_sne_rate(double age);
double determine_corecollapse_sne_total_ejected_mass(double age);
double determine_snii_yields(int k, double age);
#endif

#ifdef CLUSTER_SINK_SNIa
double determine_snia_rate(double age);
double determine_snia_yields(int k);
#endif

#ifdef CLUSTER_SINK_WINDS
double determine_winds_mass_loss_rate(double age, double zh);
double determine_winds_velocity_injection(double age, double zh);
double determine_winds_yields(int i, int k);
double determine_winds_HHe_production(double age, double z_CNO);
double determine_winds_CNO_production(double age, double z_CNO);
double determine_winds_HC_production(double age, double z_CNO);
#endif

#endif