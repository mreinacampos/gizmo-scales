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

#if !defined(CLUSTER_SINK_AVOID_MERGERS)
// function to allocate the arrays needed to collect the information of the MSPs that have to be appended
void cluster_sink_allocate_merger_loop(void) 
{
    for(int i=0; i<N_active_loc_BHs; i++)
    {
        if (BlackholeTempInfo[i].flag_SinkMerger_withMSP > 0){
            BlackholeTempInfo[i].append_MSP = (struct cluster_sink_multiple_stellar_population *) mymalloc("BHTempInfo.append_MSP", NTask * CLUSTER_SINK_NUMMSP_ACCRETE * sizeof(struct cluster_sink_multiple_stellar_population));
            memset( &BlackholeTempInfo[i].append_MSP[0], 0, NTask * CLUSTER_SINK_NUMMSP_ACCRETE * sizeof(struct cluster_sink_multiple_stellar_population) );
            if(BlackholeTempInfo[i].append_MSP == NULL) { terminate("Failed to allocate memory for BHTempInfo[Nbh].append_MSP");}
        }
    }// for(i=0; i<N_active_loc_BHs; i++)
}
#endif

#endif // CLUSTER_SINK
