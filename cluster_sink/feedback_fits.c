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

/** \brief Set the amount of mass and velocity to be injected from the multiple stellar populations
 *
 * \struct in     structure holding the input properties for the neighbour loop
 * \param i       index of the particle
 * \return        void
 */
void set_fb_input_quantities_from_msps(struct addFB_evaluate_data_in_ *in, int i, int fb_loop_iteration)
{

    int k, j;
    // set all bookkeeping quantities to zero - summ over all MSPs
    double total_mass_snii = 0, total_energy_snii = 0;
    double total_mass_snia = 0, total_energy_snia = 0;
    double total_mass_winds = 0, total_energy_winds = 0;
    double yields_snii[NUM_METAL_SPECIES], yields_snia[NUM_METAL_SPECIES], yields_winds[NUM_METAL_SPECIES]; 
    for(k=0;k<NUM_METAL_SPECIES;k++) { yields_snii[k] = 0.;  yields_snia[k] = 0.; yields_winds[k] = 0.;}

    // declare bookkeeping quantities to be used for each MSP
    double velocity_winds;

    //loop over every stellar population within the sink
    for (j = 0; j<CLUSTER_SINK_NUMMSP; j++){

        // zero out quantities
        velocity_winds = 0;

        if (P[i].MSP[j].Mass == 0) continue; // this MSP has no FB to produce

        // determine the age in Myr
        double age = evaluate_stellar_age_Gyr_for_msp(i, j)*1e3, zh;
        // metallicity of the stellar population - zh = 10^[Fe/H] = (N_Fe/N_H)_star / (N_Fe/N_H)_solar
        if (NUM_METAL_SPECIES > 1){
            zh = (P[i].MSP[j].Metallicity[NUM_METAL_SPECIES-1]/(1 - P[i].MSP[j].Metallicity[1] - P[i].MSP[j].Metallicity[0]))/(All.SolarAbundances[NUM_METAL_SPECIES-1]/(1 - All.SolarAbundances[0] - All.SolarAbundances[1]));
        } else {
            zh = (P[i].MSP[j].Metallicity[0]/All.SolarAbundances[0]);
        }

        // calculate the mass loss associated to each mechanism for this MSP
        struct fb_massloss_for_msp fb_dm;
        calculate_fb_mass_ejected_for_msps(&fb_dm, i, j);

#ifdef CLUSTER_SINK_SNII
        // total mass ejected by core-collapse SNe in code units
        total_mass_snii += fb_dm.mass_snii; // increase the total budget
        // total energy ejected per core-collapse SNe in code units
        double eps_z = DMAX(pow(zh + 1e-4, -0.12), 1);
        total_energy_snii += P[i].SNII_ThisTimeStep[j] * (eps_z*(1.e51/UNIT_ENERGY_IN_CGS));
#ifdef METALS
        // yields ejected: as ejecta masses
        double total_z = 0.;
        for(k=1;k<NUM_METAL_SPECIES;k++) { yields_snii[k] += fb_dm.mass_snii*determine_snii_yields(k, age); total_z += determine_snii_yields(k, age); }
        yields_snii[0] = fb_dm.mass_snii*1.02*total_z; // from App. A in Hopkins+18
        //yields_snii[0] = DMAX(1.02*total_z, P[i].Metallicity[0]*P[i].MSP_Mass[j]/mass_snii); // from App. A in Hopkins+18
#endif
#endif

#ifdef CLUSTER_SINK_SNIa
        // total mass ejected by SNIa in code units
        total_mass_snia += fb_dm.mass_snia; // increase the total budget
        // total energy ejected by SNIa in code units - assume 1e51ergs per SNIa
        total_energy_snia += P[i].SNIa_ThisTimeStep[j] * (1.e51/UNIT_ENERGY_IN_CGS);
#ifdef METALS
        // yields ejected: as ejecta masses
        for(k=0;k<NUM_METAL_SPECIES;k++) { yields_snia[k] += fb_dm.mass_snia*determine_snia_yields(k); }
#endif
#endif

#ifdef CLUSTER_SINK_WINDS
        // total mass ejected by winds in code units 
        total_mass_winds += fb_dm.mass_winds; // increase the total budget
        // wind injection velocity in code units 
        velocity_winds = determine_winds_velocity_injection(age, zh)/UNIT_VEL_IN_KMS;

        // total energy ejected by winds in code units -- assumed kinetic
        total_energy_winds += 0.5 * fb_dm.mass_winds * velocity_winds * velocity_winds;
#ifdef METALS
        // yields ejected: as ejecta masses
        for(k=0;k<NUM_METAL_SPECIES;k++) { yields_winds[k] += fb_dm.mass_winds*determine_winds_yields(i, age, k); }
#endif
#endif
        //printf("[MRC - set_fb_input_quantities_from_msps] - ThisTask %d, i %d, P[i].ID %d, P.Mass %g, P.Age %g - MSP j %d - Age %g, MSP_Mass %g - mass_snii %g, mass_snia %g, mass_winds %g\n", ThisTask,
        // i, P[i].ID, P[i].Mass, evaluate_stellar_age_Gyr(i), j, age, P[i].MSP_Mass[j], fb_dm.mass_snii, fb_dm.mass_snia, fb_dm.mass_winds);
    }

    // yields ejected: convert into dimensionless ejecta mass fractions - avoid NaNs
    if (total_mass_snii) for(k=0;k<NUM_METAL_SPECIES;k++) { yields_snii[k] /= total_mass_snii;}
    if (total_mass_snia) for(k=0;k<NUM_METAL_SPECIES;k++) { yields_snia[k] /= total_mass_snia;}
    if (total_mass_winds) for(k=0;k<NUM_METAL_SPECIES;k++) { yields_winds[k] /= total_mass_winds;}

    // total mass being ejected
    in->Msne = total_mass_snii + total_mass_snia + total_mass_winds;
    // velocity of the combined ejecta in code units
    in->SNe_v_ejecta = sqrt(2.*(total_energy_snii + total_energy_snia + total_energy_winds)/in->Msne); 
#ifdef METALS
    for(k=0;k<NUM_METAL_SPECIES;k++) {
        in->yields[k] = (yields_snii[k] * total_mass_snii + yields_snia[k] * total_mass_snia + yields_winds[k] * total_mass_winds)/in->Msne;
        //printf("[feedback_fits.c] - k %d, total_mass_snii %g, total_mass_snia %g, total_mass_winds %g - yields [%g, %g, %g] = [%g]\n", k, 
        //    total_mass_snii, total_mass_snia, total_mass_winds, yields_snii[k], yields_snia[k], yields_winds[k], in->yields[k]);
        assert((in->yields[k] >= 0.)&&(in->yields[k] < 1.));
    }
#endif
}

/** \brief Calculates the mass lost due to each feedbach mechanism from each multiple stellar populations
 *
 * \struct in     structure holding the mass loss associated to each mechanism
 * \param i       index of the particle
 * \param j       index of the multiple stellar population to look at
 * \return        void
 */
void calculate_fb_mass_ejected_for_msps(struct fb_massloss_for_msp *fb_dm, int i, int j)
{

    // determine the age in Myr
    double age = evaluate_stellar_age_Gyr_for_msp(i, j)*1e3, zh;
    // metallicity of the stellar population - zh = 10^[Fe/H] = (N_Fe/N_H)_star / (N_Fe/N_H)_solar
    if (NUM_METAL_SPECIES > 1){
        zh = (P[i].MSP[j].Metallicity[NUM_METAL_SPECIES-1]/(1 - P[i].MSP[j].Metallicity[1] - P[i].Metallicity[0]))/(All.SolarAbundances[NUM_METAL_SPECIES-1]/(1 - All.SolarAbundances[0] - All.SolarAbundances[1]));
    } else {
        zh = (P[i].MSP[j].Metallicity[0]/All.SolarAbundances[0]);
    }

    // zero out quantities
    fb_dm->mass_snii = 0; fb_dm->mass_snia = 0; fb_dm->mass_winds = 0;

#ifdef CLUSTER_SINK_SNII
    // total mass ejected by core-collapse SNe in code units
    fb_dm->mass_snii = P[i].SNII_ThisTimeStep[j] * (determine_corecollapse_sne_total_ejected_mass(age)/UNIT_MASS_IN_SOLAR); 
#endif

#ifdef CLUSTER_SINK_SNIa
    // total mass ejected by SNIa in code units - assume 1.4MSun per SNIa
    fb_dm->mass_snia = P[i].SNIa_ThisTimeStep[j] * (1.4/UNIT_MASS_IN_SOLAR); 
#endif

#ifdef CLUSTER_SINK_WINDS
    // particle timestep in Myr
    double dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i) * UNIT_TIME_IN_MYR;
    // total mass ejected by winds in code units 
    fb_dm->mass_winds = determine_winds_mass_loss_rate(age, zh) * P[i].MSP[j].Mass * dt; 
#endif
}

/** \brief Reduce the mass of the stellar particle that has done feedback and of its MSPs
 */
void reduce_mass_from_msps(void) 
{
    int i, j;
    struct fb_massloss_for_msp fb_dm;
    // loop over all active particles //
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if((P[i].Type!=4)&&(P[i].Type!=5)) {continue;} // has this sink output feedback in this timestep?
#ifdef BH_INTERACT_ON_GAS_TIMESTEP
        if(P[i].Type == 5 && !P[i].do_gas_search_this_timestep) {continue;}
#endif
        if(P[i].SNe_ThisTimeStep == 0) {continue;} // has this sink output feedback in this timestep?

        //loop over every stellar population within the sink
        for (j = 0; j<CLUSTER_SINK_NUMMSP; j++){

            if (P[i].MSP[j].Mass == 0) continue; // this MSP has no FB to produce

            // calculate the mass loss associated to each mechanism for this MSP
            calculate_fb_mass_ejected_for_msps(&fb_dm, i, j);

            // now remove the ejected mass from each MSP 
            P[i].MSP[j].Mass -= (fb_dm.mass_snii + fb_dm.mass_snia + fb_dm.mass_winds);
            assert(P[i].MSP[j].Mass >= 0); 

            // remove mass from the total mass too - already done in out2particle_addFB
            //P[i].Mass -= (fb_dm.mass_snii + fb_dm.mass_snia + fb_dm.mass_winds);
            //if((P[i].Mass<0)||(isnan(P[i].Mass))) {P[i].Mass=0;}
            
            //printf("[MRC - reduce_mass_from_msps] - ThisTask %d, i %d, P[i].ID %d - P.Mass %g, P.BH_Mass %g - MSP j %d - MSP_Mass %g - mass_snii %g, mass_snia %g, mass_winds %g\n", ThisTask, i, P[i].ID, P[i].Mass, P[i].BH_Mass, j, 
            //    P[i].MSP_Mass[j], fb_dm.mass_snii, fb_dm.mass_snia, fb_dm.mass_winds);
        }
    }
}


/** \brief Return the rate of SNe for a given star particle 
 *
 * \param i       index of the particle
 * \param dt       timestep of the particle in physical units
 * \return        total rate of SNe (SNII and Ia) in SNe / Myr / MSun
 */
double determine_sne_rates(int i, double dt) 
{
    int j;
    double RSNe = 0.;

    // loop over all multiple stellar populations
    for(j = 0; j<CLUSTER_SINK_NUMMSP; j++){
        
        if (P[i].MSP[j].Mass == 0) continue; // this MSP has no FB to produce

        // declare and define variables here
        double RSNII = 0., RSNIa = 0., age, p, n_sn_0 = 0;

        // determine the age in Myr
        age = evaluate_stellar_age_Gyr_for_msp(i, j)*1e3;

#ifdef CLUSTER_SINK_SNII
        RSNII = determine_corecollapse_sne_rate(age); // determine rate from analytical fit - in SNe / Myr / MSun
#endif
#ifdef CLUSTER_SINK_SNIa
        RSNIa = determine_snia_rate(age); // determine rate from analytical fit - in SNe / Myr / MSun
#endif

        // total rate of SNe
        RSNe += RSNII + RSNIa;

        // determine number of SNII
        p = RSNII * (P[i].MSP[j].Mass*UNIT_MASS_IN_SOLAR) * (dt*UNIT_TIME_IN_MYR); // unit conversion factor
        n_sn_0=(float)floor(p);
        p-=n_sn_0;
        if(get_random_number(P[i].ID+6) < p) {n_sn_0++;} // determine if SNe occurs
        P[i].SNII_ThisTimeStep[j] = n_sn_0; // assign to particle

        // determine number of SNIa
        p = RSNIa * (P[i].MSP[j].Mass*UNIT_MASS_IN_SOLAR) * (dt*UNIT_TIME_IN_MYR); // unit conversion factor
        n_sn_0=(float)floor(p);
        p-=n_sn_0;
        if(get_random_number(P[i].ID+6) < p) {n_sn_0++;} // determine if SNe occurs
        P[i].SNIa_ThisTimeStep[j] = n_sn_0; // assign to particle

        // total number of SNe per timestep over all MSPs
        P[i].SNe_ThisTimeStep += P[i].SNII_ThisTimeStep[j] + P[i].SNIa_ThisTimeStep[j]; 
#ifdef CLUSTER_SINK_DEBUG_ONESNE
        // debug: only one SNe
        if (P[i].CumNumSNe[j] > 0){ P[i].SNe_ThisTimeStep[j] = 0; RSNe = 0; }
#endif

#ifdef CLUSTER_SINK_WINDS 
        // do wind FB even when there are no SNe
        if (P[i].SNe_ThisTimeStep == 0){ RSNe += 1; P[i].SNe_ThisTimeStep += 0.2;}
#endif 
    }

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
double determine_winds_yields(int i, double age, int k) 
{

    // determine the CNO metallicity
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


# ifdef CLUSTER_SINK_RADIATION
/** \brief Calculate the bolometric luminosity for a SSP of a given age
 *
 * \param age_in_gyr        age of the SSP
 * \param i                 index of the particle
 * \return                  light-to-mass ratio in Lsun/Msun
 */
double calculate_relative_light_to_mass_ratio(double age_in_gyr, int i, int j){

    double zh, age_in_myr = age_in_gyr*1e3, light_to_mass = 0, slope = 0;
    // metallicity of the stellar population - zh = 10^[Fe/H] = (N_Fe/N_H)_star / (N_Fe/N_H)_solar
    if (NUM_METAL_SPECIES > 1){
        zh = (P[i].Metallicity[NUM_METAL_SPECIES-1]/(1 - P[i].Metallicity[1] - P[i].Metallicity[0]))/(All.SolarAbundances[NUM_METAL_SPECIES-1]/(1 - All.SolarAbundances[0] - All.SolarAbundances[1]));
    } else { // asume primordial composition
        zh = (P[i].Metallicity[0]/All.SolarAbundances[0]);
    }

    // metallicity-dependent coefficients aL,1, aL,2 and aL,3 - in LSun/MSun
    double RAD_coeff_alj[3] = {800, 1100*pow(zh, -0.1), 0.163}; 

    // bolometric luminosities per unit stellar mass - eq. 1 in Hopkins+ 22
    if (age_in_myr <= RAD_tlj[0]){ light_to_mass = RAD_coeff_alj[0]; }
    else if ((age_in_myr > RAD_tlj[0]) && (age_in_myr <= RAD_tlj[1])){
        slope = log(RAD_coeff_alj[1]/RAD_coeff_alj[0])/log(RAD_tlj[1]/RAD_tlj[0]);
        light_to_mass = RAD_coeff_alj[0]*pow(age_in_myr/RAD_tlj[0], slope);
    } else if (age_in_myr > RAD_tlj[1]){
        slope = -1.82*(1 - 0.1*(1 - 0.073*log(age_in_myr/RAD_tlj[1]))*log(age_in_myr/RAD_tlj[1]));
        double coeff_fl = 1 + 1.2*exp(- pow(log(age_in_myr/RAD_tlj[2])/RAD_coeff_alj[2], 2));
        light_to_mass = RAD_coeff_alj[1]*pow(age_in_myr/RAD_tlj[1], slope)*coeff_fl;
    }

    return light_to_mass;
}

/** \brief Determine the fraction of bolometric flux that corresponds to ionizing radiation
 *
 * \param age_in_gyr        age of the SSP
 * \param i                 index of the particle
 * \return                  fraction of ionizing flux
 */
double determine_ionizing_flux_fraction(double age_in_gyr, int i){

    double fraction = 0, age_in_myr = age_in_gyr*1e3;

    // Hopkins+ 22 - paragraph after eq. 1
    if (age_in_myr < RAD_ION_tion[0]){ fraction = RAD_ION_coeff_fion;}
    else if ((age_in_myr >= RAD_ION_tion[0]) && (age_in_myr <= RAD_ION_tion[1])){
        fraction = RAD_ION_coeff_fion * pow(age_in_myr/RAD_ION_tion[0], -2.9);
    } else if (age_in_myr > RAD_ION_tion[1]){ fraction = 0;}

    return fraction;
}
#endif // CLUSTER_SINK_RADIATION


#endif // CLUSTER_SINK
