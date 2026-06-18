#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/* Routines to calculate the formation of several stellar populations within sink particles
 * This file was written by Marta Reina-Campos (reinacampos@mcmaster.ca) for GIZMO.
 */

#ifdef CLUSTER_SINK

/** \brief Top-level routine to spawn the continuous stellar populations 
 *
 * \param i       index of the particle
 * \param dt       timestep of the particle in physical units
 * \return        total rate of SNe (SNII and Ia) in SNe / Myr / MSun
 */
void continuous_star_formation_in_sinks(void) 
{

    // loop over all active sinks
    // if there's at least XYZ of gas mass, convert to a stellar population
    // record: initial and current mass of stellar population, age and metallicity
    // think: which SFE to use?

    double sp_mass = 0;
    int i, j;
    // loop over particles //
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type != 5) {continue;} // only sinks can form multiple stellar populations

        // assume all accreted gas mass will form stars (SFE = 100%) 
        // calculate the amount of gas available within the sink
        sp_mass = determine_mass_gas_reservoir(i); 
        
        // if there's not enough mass, keep accreting
        if (sp_mass < All.ClusterSink_MinGasMass) {continue;}

        // loop over all the MSPs
        for (j = 0; j<CLUSTER_SINK_NUMMSP; j++){

            if (P[i].MSP[j].InitialMass > 0){ // for every existing MSP - checkrelative its initial age and metallicity
                if (evaluate_initial_stellar_age_Gyr_for_msp(i, j)*1e3 < All.ClusterSink_Delta_AgeInMyr){ // if the MSP is younger than set in the param file, add the mass here as a mass-weight
                    // only combine MSPs if the metallicity difference \Delta [Z/ZSun] is less than set in the param file
                    if (fabs(log10(P[i].Metallicity[0]/All.SolarAbundances[0]) - log10(P[i].MSP[j].InitialMetallicity_Z/All.SolarAbundances[0])) < All.ClusterSink_Delta_ZZSun){
                        double m0 = P[i].MSP[j].Mass, mf = P[i].MSP[j].Mass + sp_mass;
                        P[i].MSP[j].Age = (m0/mf)*P[i].MSP[j].Age + (sp_mass/mf)*All.Time;
                        for(int k=0;k<NUM_METAL_SPECIES;k++) {P[i].MSP[j].Metallicity[k] = (m0/mf)*P[i].MSP[j].Metallicity[k] + (sp_mass/mf)*P[i].Metallicity[k];}
                        P[i].MSP[j].Mass += sp_mass;
                        P[i].MSP[j].InitialMass += sp_mass;
                        assert(P[i].MSP[j].Age <= All.Time); // check that no MSP ends up with spurious ages
                        break;
                    } else { continue; } 
                } else { continue; } }

            // if no existing MSP is younger than 0.5Myr, create a new one
            P[i].MSP[j].InitialMass = sp_mass;
            P[i].MSP[j].Mass = sp_mass;
            P[i].MSP[j].Age = All.Time; // scale factor or time - needs to be evaluated with evaluate_stellar_age_Gyr_for_msp(i, j)
            P[i].MSP[j].InitialAge = All.Time; // scale factor or time - needs to be evaluated with evaluate_stellar_age_Gyr_for_msp(i, j)
            for(int k=0;k<NUM_METAL_SPECIES;k++) {P[i].MSP[j].Metallicity[k] = P[i].Metallicity[k];} // collecting the mass-weighted metallicity of accreted gas
            P[i].MSP[j].InitialMetallicity_Z = P[i].Metallicity[0]; // we only need total Z to determine whether to spawn other MSPs
            break;
        }

        if (j >= CLUSTER_SINK_NUMMSP){printf("[WARNING - formation_stellarpops.c] ThisTask %d, P.ID %d - MISSING MSP Mass %g because array is full\n", ThisTask, P[i].ID, sp_mass);}

    } // loop over active particles
}




/** \brief Calculate mass in the gas reservoir for a given sink
 *
 * \param i       index of the particle
 * \return        total mass available in gas
 */
double determine_mass_gas_reservoir(int i) 
{
    int k;
    double mass_gas, mass_msp; mass_gas = 0; mass_msp = 0;
    // calculate total mass currently in stars within the sink
    for (k=0; k<CLUSTER_SINK_NUMMSP; k++){ mass_msp += P[i].MSP[k].Mass; }
    // calculate the gas mass
    mass_gas = P[i].Mass - mass_msp;
    // avoid negative values from precision round-offs
    if (mass_gas < 0) {mass_gas = 0;} // printf("[WARNING - formation_stellarpops.c - mgas] mass_gas %g, P[i].Mass %g, mass_msp %g\n", mass_gas, P[i].Mass, mass_msp);}
    assert(mass_gas >= 0.);
    return mass_gas;
}

/** \brief Return the mass-weighted stellar age in Gyr for a given MSP, needed throughout for stellar feedback
 *
 * \param i       index of the particle
 * \param j       index of the MSP
 * \return        age of the MSP in Gyr
 */
double evaluate_stellar_age_Gyr_for_msp(long i, int j)
{
    double age = evaluate_time_since_t_initial_in_Gyr(P[i].MSP[j].Age);
    age = DMAX(age, 1.e-5); // set a floor for some routines
    return age;
}

/** \brief Return the initial stellar age in Gyr for a given MSP
 *
 * \param i       index of the particle
 * \param j       index of the MSP
 * \return        initial age of the MSP in Gyr
 */
double evaluate_initial_stellar_age_Gyr_for_msp(long i, int j)
{
    double age = evaluate_time_since_t_initial_in_Gyr(P[i].MSP[j].InitialAge);
    age = DMAX(age, 1.e-5); // set a floor for some routines
    return age;
}

#endif // CLUSTER_SINK
