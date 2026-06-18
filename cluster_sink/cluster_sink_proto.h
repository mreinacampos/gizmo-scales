
/*
 * Declarations of routines within the CLUSTER_SINK module
 * This file was written by Marta Reina-Campos (reinacampos@mcmaster.ca) for GIZMO.
 */

#ifdef CLUSTER_SINK
void continuous_star_formation_in_sinks(void);
double determine_mass_gas_reservoir(int i);
double determine_sne_rates(int i, double dt);
double evaluate_stellar_age_Gyr_for_msp(long i, int j);
double evaluate_initial_stellar_age_Gyr_for_msp(long i, int j);

// structure to contain the mass loss of a given MSP
struct fb_massloss_for_msp
{
    double mass_snii, mass_snia, mass_winds;
} *fb_dm;

void reduce_mass_from_msps(void);
void calculate_fb_mass_ejected_for_msps(struct fb_massloss_for_msp *fb_dm, int i, int j);

#ifndef CLUSTER_SINK_AVOID_MERGERS
void cluster_sink_allocate_merger_loop(void);
#endif


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
double determine_winds_yields(int i, double age, int k);
double determine_winds_HHe_production(double age, double z_CNO);
double determine_winds_CNO_production(double age, double z_CNO);
double determine_winds_HC_production(double age, double z_CNO);
#endif

#ifdef CLUSTER_SINK_RADIATION
double calculate_relative_light_to_mass_ratio(double age_in_gyr, int i, int j);
double determine_ionizing_flux_fraction(double age_in_gyr, int i);
#endif

#endif