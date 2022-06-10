#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/* Routines to calculate different feedback mechanisms rates, masses, energies and yields
 * based on the IMF-integrated analytical fits from https://ui.adsabs.harvard.edu/abs/2022arXiv220300040H/abstract
 * This file was written by Marta Reina-Campos (reinacampos@mcmaster.ca) for GIZMO.
 */

#ifdef CLUSTER_SINK

// load the header with all the coefficients
#include "./feedback_fits.h"


/** \brief Return the rate of SNe for a given star particle 
 *
 * \param i       index of the particle
 * \param dt       timestep of the particle in physical units
 * \return        total rate of SNe (SNII and Ia) in SNe / Myr / MSun
 */
double determine_sne_rate(int i, double dt) 
{
    double RSNe = 0., RSNII = 0., RSNIa = 0., age, p, n_sn_0 = 0;
    // determine the age in Myr
    age = evaluate_stellar_age_Gyr(P[i].StellarAge)*1e3;

#ifdef CLUSTER_SINK_SNII
    RSNII = determine_corecollapse_sne_rate(age); // determine rate from analytical fit - in SNe / Myr / MSun
#endif
#ifdef CLUSTER_SINK_SNIa
    RSNIa = determine_snia_rate(age); // determine rate from analytical fit - in SNe / Myr / MSun
#endif

    // total rate of SNe
    RSNe = RSNII + RSNIa;

    // determine number of SNII
    p = RSNII * (P[i].Mass*UNIT_MASS_IN_SOLAR) * (dt*UNIT_TIME_IN_MYR); // unit conversion factor
    n_sn_0=(float)floor(p);
    p-=n_sn_0;
    if(get_random_number(P[i].ID+6) < p) {n_sn_0++;} // determine if SNe occurs
    P[i].SNII_ThisTimeStep = n_sn_0; // assign to particle

    // determine number of SNIa
    p = RSNIa * (P[i].Mass*UNIT_MASS_IN_SOLAR) * (dt*UNIT_TIME_IN_MYR); // unit conversion factor
    n_sn_0=(float)floor(p);
    p-=n_sn_0;
    if(get_random_number(P[i].ID+6) < p) {n_sn_0++;} // determine if SNe occurs
    P[i].SNIa_ThisTimeStep = n_sn_0; // assign to particle


    // total number of SNe per timestep
    P[i].SNe_ThisTimeStep = P[i].SNII_ThisTimeStep + P[i].SNIa_ThisTimeStep; 
#ifdef CLUSTER_SINK_DEBUG
    // debug: only one SNe
    if (P[i].CumNumSNe > 0){ P[i].SNe_ThisTimeStep = 0; RSNe = 0; }
#endif
    return RSNe;
}


#ifdef CLUSTER_SINK_SNII
/** \brief Return the rate of core-collapse SNe for a given star particle 
 * using the tables in feedback_fits.h
 *
 * \param age       age of the stellar population in Myr
 * \return          rate of core-collapse SNe in SNe / Myr / MSun
 */
double determine_corecollapse_sne_rate(double age) 
{
    double rate = 0., slope;
    // power-law analytical fit
    if((age < SNII_tsj[0])|| (age > SNII_tsj[2])){ rate = 0.;}
    else if ((age >= SNII_tsj[0]) && (age <= SNII_tsj[1])) {
        slope = log(SNII_coeff_asj[1]/SNII_coeff_asj[0])/log(SNII_tsj[1]/SNII_tsj[0]);
        rate = SNII_coeff_asj[0]*pow( age / SNII_tsj[0], slope);
    } else if ((age >= SNII_tsj[1]) && (age <= SNII_tsj[2])) {
        slope = log(SNII_coeff_asj[2]/SNII_coeff_asj[1])/log(SNII_tsj[2]/SNII_tsj[1]);
        rate = SNII_coeff_asj[1]*pow( age / SNII_tsj[1], slope);
    }

    // return rate in SNe / Myr / MSun
    return rate*1e-3;
}

/** \brief Return the total mass ejected per core-collapse SNe
 *
 * \param age       age of the stellar population in Myr
 * \return          total mass ejected per SNII in MSun
 */
double determine_corecollapse_sne_total_ejected_mass(double age) 
{
    double mass, slope;
    // power-law analytical fit - eq. 2 in Hopkins+22
    if(age <= 6.5){ slope = 2.22;}
    else if (age > 6.5){ slope = 0.267;}
    mass = 10 * pow(age / 6.5, -slope);
    // return ejected mass in MSun
    return mass;
}
#endif // CLUSTER_SINK_SNII

#ifdef CLUSTER_SINK_SNIa
/** \brief Return the rate of SNIa for a given star particle 
 * using the tables in feedback_fits.h
 *
 * \param age       age of the stellar population in Myr
 * \return          rate of SNIa in SNe / Myr / MSun
 */
double determine_snia_rate(double age) 
{
    double rate;
    // power-law analytical fit - eq. 3 in Hopkins+22
    if(age < SNIa_tj[0]){ rate = 0.;}
    else if (age >= SNIa_tj[0]) {
        rate = SNIa_coeff_aj[0]*pow( age / SNIa_tj[0], SNIa_slope[0]);
    }

    // return rate in SNe / Myr / MSun
    return rate*1e-3;
}
#endif // CLUSTER_SINK_SNIa


#ifdef CLUSTER_SINK_WINDS
/** \brief Return the mass loss from AGB&OB winds for a given star particle 
 * using the tables in feedback_fits.h
 *
 * \param age       age of the stellar population in Myr
 * \param zh        metallicity of the stellar population - zh = 10^{[Fe/H]}
 * \return          mass loss in Myr^-1
 */
double determine_winds_mass_loss_rate(double age, double zh) 
{

    // metallicity-dependent coefficients aw,1, aw,2 and aw,3 - in Gyr^-1
    double WINDS_coeff_awj[3] = {3*pow(zh, 0.87), 20*pow(zh, 0.45), 0.6*zh}; 

    double mass_loss = 0., slope = 0.;
    // power-law analytical fit - eq. 4 in Hopkins+22
    if(age < WINDS_twj[0]){ mass_loss = WINDS_coeff_awj[0]; }
    else if((age > WINDS_twj[0]) && (age <= WINDS_twj[1])){
        slope = log(WINDS_coeff_awj[1]/WINDS_coeff_awj[0])/log(WINDS_twj[1]/WINDS_twj[0]);
        mass_loss = WINDS_coeff_awj[0] * pow(age / WINDS_twj[0], slope);
    } else if((age > WINDS_twj[1]) && (age <= WINDS_twj[2])){
        slope = log(WINDS_coeff_awj[2]/WINDS_coeff_awj[1])/log(WINDS_twj[2]/WINDS_twj[1]);
        mass_loss = WINDS_coeff_awj[1] * pow(age / WINDS_twj[1], slope);
    } else if(age > WINDS_twj[2]){
        slope = -3.1;
        mass_loss = WINDS_coeff_awj[2] * pow(age / WINDS_twj[2], slope);
    }

    mass_loss += WINDS_coeff_aaj[0] / ((1 + pow(age/WINDS_twj[3] , 1.1) * (1 + WINDS_coeff_aaj[1] * pow(age/WINDS_twj[3] , -1))));

    // return mass_loss in Myr^-1
    return mass_loss*1e-3;
}

/** \brief Return the injection velocity from AGB&OB winds for a given star particle 
 * using the tables in feedback_fits.h
 *
 * \param age       age of the stellar population in Myr
 * \param zh        metallicity of the stellar population - zh = 10^{[Fe/H]}
 * \return          velocity in km/s
 */
double determine_wind_velocity_injection(double age, double zh) 
{

    // analytical fit - eq. 5 in Hopkins+22
    double velocity = pow(zh, 0.12) * ( 3000 / ( 1 + pow(age/WINDS_tvj[0], 2.5) ) + 600 / (1+pow(zh, 3)*pow(age/WINDS_tvj[1], 6) + 11.2*pow(zh, 1.5)) + 30); 

    // return mass_loss in km/s
    return velocity;
}
#endif // CLUSTER_SINK_WINDS

#endif // CLUSTER_SINK
