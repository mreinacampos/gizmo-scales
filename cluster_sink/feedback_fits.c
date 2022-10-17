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
#ifdef CLUSTER_SINK_DEBUG_ONESNE
    // debug: only one SNe
    if (P[i].CumNumSNe > 0){ P[i].SNe_ThisTimeStep = 0; RSNe = 0; }
#endif

#ifdef CLUSTER_SINK_WINDS 
    // do wind FB when there are no SNe
    if (P[i].SNe_ThisTimeStep == 0){ RSNe = 1; P[i].SNe_ThisTimeStep += 0.2;}
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
    double rate = 0., slope = 0.;
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
    double mass = 0., slope = 0.;
    // power-law analytical fit - eq. 2 in Hopkins+22
    if(age <= 6.5){ slope = 2.22;}
    else if (age > 6.5){ slope = 0.267;}
    mass = 10 * pow(age / 6.5, -slope);
    // return ejected mass in MSun
    return mass;
}
/** \brief Return the yield k ejected by SNII
 *
 * \param k       index of the yield to return
 * \param age     age of the stellar population in Myr
 * \return        fraction of ejecta mass in species k
 */
double determine_snii_yields(int k, double age) 
{
    double yield = 0., slope = 0.;

    // power-law analytical fit - eq. 7 in Hopkins+22
    if((age <= SNII_yields_tccj[0]) || (age >= SNII_yields_tccj[4])){ yield = 0.;}
    else if ((age > SNII_yields_tccj[0]) && (age <= SNII_yields_tccj[1])) {
        slope = log(SNII_yields_accj[k-1][1]/SNII_yields_accj[k-1][0])/log(SNII_yields_tccj[1]/SNII_yields_tccj[0]);
        yield = SNII_yields_accj[k-1][0]*pow( age / SNII_yields_tccj[0], slope);
    } else if ((age > SNII_yields_tccj[1]) && (age <= SNII_yields_tccj[2])) {
        slope = log(SNII_yields_accj[k-1][2]/SNII_yields_accj[k-1][1])/log(SNII_yields_tccj[2]/SNII_yields_tccj[1]);
        yield = SNII_yields_accj[k-1][1]*pow( age / SNII_yields_tccj[1], slope);
    } else if ((age > SNII_yields_tccj[2]) && (age <= SNII_yields_tccj[3])) {
        slope = log(SNII_yields_accj[k-1][3]/SNII_yields_accj[k-1][2])/log(SNII_yields_tccj[3]/SNII_yields_tccj[2]);
        yield = SNII_yields_accj[k-1][2]*pow( age / SNII_yields_tccj[2], slope);
    } else if ((age > SNII_yields_tccj[3]) && (age <= SNII_yields_tccj[4])) {
        slope = log(SNII_yields_accj[k-1][4]/SNII_yields_accj[k-1][3])/log(SNII_yields_tccj[4]/SNII_yields_tccj[3]);
        yield = SNII_yields_accj[k-1][3]*pow( age / SNII_yields_tccj[3], slope);
    }

    return yield;
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
    double rate = 0.;
    // power-law analytical fit - eq. 3 in Hopkins+22
    if(age < SNIa_tj[0]){ rate = 0.;}
    else if (age >= SNIa_tj[0]) {
        rate = SNIa_coeff_aj[0]*pow( age / SNIa_tj[0], SNIa_slope[0]);
    }

    // return rate in SNe / Myr / MSun
    return rate*1e-3;
}

/** \brief Return the yield k ejected by SNIa
 *
 * \param k       index of the yield to return
 * \return        fraction of ejecta mass in species k
 */
double determine_snia_yields(int k) 
{
    return SNIa_yields[k];
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
double determine_winds_velocity_injection(double age, double zh) 
{

    // analytical fit - eq. 5 in Hopkins+22
    double velocity = pow(zh, 0.12) * ( 3000 / ( 1 + pow(age/WINDS_tvj[0], 2.5) ) + 600 / (1+pow(zh, 3)*pow(age/WINDS_tvj[1], 6) + 11.2*pow(zh, 1.5)) + 30); 

    // return mass_loss in km/s
    return velocity;
}


/** \brief Return the production of He from H by ABG/OB winds
 *
 * \param age     age of the star in Gyr
 * \param z_CNO   CNO-based metallicity
 * \return        production of He from H by ABG/OB winds in mass fraction
 */
double determine_winds_HHe_production(double age, double z_CNO) 
{
    // Table 2 in Hopkins+22
    // slope of the first piece of the piece-wise analytical fit
    double HHe_slope0 = 3; 
    // timescales of the piece-wise analytical fit - in Gyr
    double HHe_timescales[5] = {0.0028, 0.01, 2.3, 3.0, 100}; 
    // metallicity-dependent coefficients of the piece-wise analytical fit - in Gyr
    double HHe_coeff[5] = {0.4*DMIN( pow(z_CNO+0.001, 0.6), 2), 0.08, 0.07, 0.042, 0.042}; 

    double yield = 0, slope = 0;
    if (age <= HHe_timescales[0]){
        slope = HHe_slope0;
        yield = HHe_coeff[0] * pow(age/HHe_timescales[0], slope);
    } else if ( (age > HHe_timescales[0]) && (age <= HHe_timescales[1]) ){
        slope = log(HHe_coeff[1]/HHe_coeff[0])/log(HHe_timescales[1]/HHe_timescales[0]);
        yield = HHe_coeff[0] * pow(age/HHe_timescales[0], slope);
    } else if ( (age > HHe_timescales[1]) && (age <= HHe_timescales[2]) ){
        slope = log(HHe_coeff[2]/HHe_coeff[1])/log(HHe_timescales[2]/HHe_timescales[1]);
        yield = HHe_coeff[1] * pow(age/HHe_timescales[1], slope);
    } else if ( (age > HHe_timescales[2]) && (age <= HHe_timescales[3]) ){
        slope = log(HHe_coeff[3]/HHe_coeff[2])/log(HHe_timescales[3]/HHe_timescales[2]);
        yield = HHe_coeff[2] * pow(age/HHe_timescales[2], slope);
    } else if ( (age > HHe_timescales[3]) && (age <= HHe_timescales[4]) ){
        slope = log(HHe_coeff[4]/HHe_coeff[3])/log(HHe_timescales[4]/HHe_timescales[3]);
        yield = HHe_coeff[3] * pow(age/HHe_timescales[3], slope);
    }

    return yield;
}

/** \brief Return the production from CNO cycle by ABG/OB winds
 *
 * \param age     age of the star in Gyr
 * \param z_CNO   CNO-based metallicity
 * \return        production from CNO cycle by ABG/OB winds in mass fraction
 */
double determine_winds_CNO_production(double age, double z_CNO) 
{
    // Table 2 in Hopkins+22
    // slope of the first piece of the piece-wise analytical fit
    double CNO_slope0 = 3.5; 
    // timescales of the piece-wise analytical fit - in Gyr
    double CNO_timescales[6] = {0.001, 0.0028, 0.05, 1.9, 14, 100}; 
    // metallicity-dependent coefficients of the piece-wise analytical fit - in Gyr
    double CNO_coeff[6] = {0.2*DMIN( pow(z_CNO, 2)+1e-4, 0.9 ), 0.68*DMIN( pow(z_CNO+0.001, 0.1), 0.9 ), 0.4, 0.23, 0.065, 0.065}; 

    double yield = 0, slope = 0;
    if (age <= CNO_timescales[0]){
        slope = CNO_slope0;
        yield = CNO_coeff[0] * pow(age/CNO_timescales[0], slope);
    } else if ( (age > CNO_timescales[0]) && (age <= CNO_timescales[1]) ){
        slope = log(CNO_coeff[1]/CNO_coeff[0])/log(CNO_timescales[1]/CNO_timescales[0]);
        yield = CNO_coeff[0] * pow(age/CNO_timescales[0], slope);
    } else if ( (age > CNO_timescales[1]) && (age <= CNO_timescales[2]) ){
        slope = log(CNO_coeff[2]/CNO_coeff[1])/log(CNO_timescales[2]/CNO_timescales[1]);
        yield = CNO_coeff[1] * pow(age/CNO_timescales[1], slope);
    } else if ( (age > CNO_timescales[2]) && (age <= CNO_timescales[3]) ){
        slope = log(CNO_coeff[3]/CNO_coeff[2])/log(CNO_timescales[3]/CNO_timescales[2]);
        yield = CNO_coeff[2] * pow(age/CNO_timescales[2], slope);
    } else if ( (age > CNO_timescales[3]) && (age <= CNO_timescales[4]) ){
        slope = log(CNO_coeff[4]/CNO_coeff[3])/log(CNO_timescales[4]/CNO_timescales[3]);
        yield = CNO_coeff[3] * pow(age/CNO_timescales[3], slope);
    } else if ( (age > CNO_timescales[4]) && (age <= CNO_timescales[5]) ){
        slope = log(CNO_coeff[5]/CNO_coeff[4])/log(CNO_timescales[5]/CNO_timescales[4]);
        yield = CNO_coeff[4] * pow(age/CNO_timescales[4], slope);
    }

    return yield;
}

/** \brief Return the production of C from H by ABG/OB winds
 *
 * \param age     age of the star in Gyr
 * \param z_CNO   CNO-based metallicity
 * \return        production of C from H by ABG/OB winds in mass fraction
 */
double determine_winds_HC_production(double age, double z_CNO) 
{
    // Table 2 in Hopkins+22
    // slope of the first piece of the piece-wise analytical fit
    double HC_slope0 = 3; 
    // timescales of the piece-wise analytical fit - in Gyr
    double HC_timescales[4] = {0.005, 0.04, 10, 100}; 
    // metallicity-dependent coefficients of the piece-wise analytical fit - in Gyr
    double HC_coeff[4] = {1e-6, 0.001, 0.005, 0.005}; 

    double yield = 0, slope = 0;
    if (age <= HC_timescales[0]){
        slope = HC_slope0;
        yield = HC_coeff[0] * pow(age/HC_timescales[0], slope);
    } else if ( (age > HC_timescales[0]) && (age <= HC_timescales[1]) ){
        slope = log(HC_coeff[1]/HC_coeff[0])/log(HC_timescales[1]/HC_timescales[0]);
        yield = HC_coeff[0] * pow(age/HC_timescales[0], slope);
    } else if ( (age > HC_timescales[1]) && (age <= HC_timescales[2]) ){
        slope = log(HC_coeff[2]/HC_coeff[1])/log(HC_timescales[2]/HC_timescales[1]);
        yield = HC_coeff[1] * pow(age/HC_timescales[1], slope);
    } else if ( (age > HC_timescales[2]) && (age <= HC_timescales[3]) ){
        slope = log(HC_coeff[3]/HC_coeff[2])/log(HC_timescales[3]/HC_timescales[2]);
        yield = HC_coeff[2] * pow(age/HC_timescales[2], slope);
    } 

    return yield;
}

/** \brief Return the yield k ejected by ABG/OB winds
 *
 * \param i       index of the particle
 * \param k       index of the yield to return
 * \return        fraction of ejecta mass in species k
 */
double determine_winds_yields(int i, int k) 
{

    // determine the age in Gyr
    double age = evaluate_stellar_age_Gyr(P[i].StellarAge);
    double z_CNO = (P[i].Metallicity[2] + P[i].Metallicity[3] + P[i].Metallicity[4])/(All.SolarAbundances[2]+All.SolarAbundances[3]+All.SolarAbundances[4]);

    double yield = 0.;

    double y_HHe, y_HeC, y_HC, y_CN, y_ON, y_CNO, f_h0;
    y_HHe = determine_winds_HHe_production(age, z_CNO);
    y_CNO = determine_winds_CNO_production(age, z_CNO);
    y_HC = determine_winds_HC_production(age, z_CNO);

    // ratio of the initial O to C abundances
    double x_OC = (P[i].Metallicity[5] / P[i].Metallicity[3]);
    // secondary production of N from C and O
    y_CN = DMIN(1, 0.5 * y_CNO * (1 + x_OC));
    y_ON = y_CNO + (y_CNO - y_CN)/x_OC;
    // production of C from He and H
    y_HeC = y_HC;
    // initial hydrogen abundance: 1 - f_He,0 - f_Z,0
    f_h0 = 1 - P[i].Metallicity[1] - P[i].Metallicity[0];

    // assume initial surface abundances for total metallicity and heavy elements
    if ((k == 0) || (k > 4)) { yield = P[i].Metallicity[k]; }
    else if (k == 1) { // He
        yield = P[i].Metallicity[1] * (1 - y_HeC) + y_HHe * f_h0;
    } else if (k == 2) { // C
        yield = P[i].Metallicity[2] * (1 - y_CN) + y_HeC * P[i].Metallicity[1] + y_HC * f_h0 * (1 - y_HHe);
    } else if (k == 3) { // N
        yield = P[i].Metallicity[3] + y_CN * P[i].Metallicity[2] + y_ON * P[i].Metallicity[4];
    } else if (k == 4) { // O
        yield = P[i].Metallicity[4] * (1 - y_ON);
    }
    return yield;
}
#endif // CLUSTER_SINK_WINDS

#endif // CLUSTER_SINK
