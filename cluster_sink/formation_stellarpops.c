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

    double sp_mass;
    int i, j;
    // loop over particles //
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type != 5) {continue;} // only sinks can form multiple stellar populations

        // assume all accreted gas mass will form stars (SFE = 100%) 
        // calculate the amount of gas available within the sink
        sp_mass = determine_mass_gas_reservoir(i); 

        // if there's not enough mass, keep accreting
        if (sp_mass < All.ClusterSink_MinGasMass) continue;

        // loop until there's an empty entry
        for (j = 0; j<CLUSTER_SINK_NUMMSP; j++){
            if (P[i].MSP_InitialMass[j] > 0) continue;
            printf("[formation_stellarpops.c] - i %d j %d, sp_mass %g, P[i].MSP_InitialMass[j] %g\n", i, j, sp_mass, P[i].MSP_InitialMass[j]);
            // assign properties
            P[i].MSP_InitialMass[j] = sp_mass;
            P[i].MSP_Mass[j] = sp_mass;
            P[i].MSP_Age[j] = All.Time; // scale factor or time - needs to be evaluated with evaluate_stellar_age_Gyr_for_msp(i, j)
            P[i].MSP_Metallicity[j] = P[i].Metallicity[0]; // MRC - need to collect the mass-weighted metallicity of accreted gas
            printf("[formation_stellarpops.c] - i %d j %d, MSP_InitialMass[j] %g, MSP_Mass %g, MSP_Age %g, MSP_Metallicity %g\n",
             i, j, P[i].MSP_InitialMass[j], P[i].MSP_Mass[j], P[i].MSP_Age[j], P[i].MSP_Metallicity[j]);
            break;
        }

    // MRC - assert that mass budget is correct
    // assert that we can keep on going even if end of array is reached

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
    for (k=0; k<CLUSTER_SINK_NUMMSP; k++){ mass_msp += P[i].MSP_Mass[k]; }
    // calculate the gas mass
    mass_gas = P[i].Mass - mass_msp;
    if (mass_gas < 0) printf("[formation_stellarpops.c - mgas] mass_gas %g, P[i].Mass %g, mass_msp %g\n", mass_gas, P[i].Mass, mass_msp);
    if (abs(mass_gas) < 1e-8) {mass_gas = 0;} // avoid negative values from precision round-offs
    assert(mass_gas >= 0.);
    return mass_gas;
}

/** \brief Return the stellar age in Gyr for a given MSP, needed throughout for stellar feedback
 *
 * \param i       index of the particle
 * \param j       index of the MSP
 * \return        age of the MSP in Gyr
 */
double evaluate_stellar_age_Gyr_for_msp(long i, int j)
{
    double age = evaluate_time_since_t_initial_in_Gyr(P[i].MSP_Age[j]);
    age = DMAX(age, 1.e-5); // set a floor for some routines
    return age;
}

#endif // CLUSTER_SINK
