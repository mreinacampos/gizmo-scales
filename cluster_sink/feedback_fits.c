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


/** \brief Return the rate of core-collapse SNe for a given star particle 
 * using the tables in feedback_fits.h
 *
 * \param age       age of the stellar population in Myr
 * \return          rate of core-collapse SNe in SNe / Myr / MSun
 */
double determine_corecollapse_sne_rate(double age) 
{
    double rate, slope;
    // power-law analytical fit
    if((age < SNII_tsj[0])|| (age > SNII_tsj[2])){
        rate = 0.;
    } else if ((age >= SNII_tsj[0]) && (age <= SNII_tsj[1])) {
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
    // power-law analytical fit
    if(age <= 6.5){ slope = 2.22;}
    else if (age > 6.5){ slope = 0.267;}
    mass = 10 * pow(age / 6.5, -slope);
    // return ejected mass in MSun
    return mass;
}

#endif // CLUSTER_SINK
