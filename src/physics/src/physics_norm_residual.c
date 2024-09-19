/**
 * \file   physics_norm_residual.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief Residual discrete vector norm functions.
 * 
 * Created on August 29, 2018, 9:20 AM
 */

/* header files */
#include "physics_norm_residual.h"


double physics_norm_L2_residual(int dim,int nfields,int nelem_subgrid,int ncell,int* cell_info,double *residual){

    int i,j;
    double resnorm = 0.0;


    if(dim==2){

        for(i=0;i<ncell;++i){

            const int* cell_info_loc = &cell_info[d_physics_interface.extern_cell_size*i+0];
            //const int cell_type = cell_info[0];
            const int ind = cell_info_loc[2];

            const int nftm = nfields*nelem_subgrid*nelem_subgrid;
            double* res = &residual[ind];

            for(j=0;j<nftm;j++){
                resnorm += res[j]*res[j];
            }
        }

    }else{

        for(i=0;i<ncell;++i){

            const int* cell_info_loc = &cell_info[d_physics_interface.extern_cell_size*i+0];
            //const int cell_type = cell_info[0];
            const int ind = cell_info_loc[2];

            const int nftm = nfields*nelem_subgrid*nelem_subgrid*nelem_subgrid;
            double* res = &residual[ind];

            for(j=0;j<nftm;j++){
                resnorm += res[j]*res[j];
            }
        }

    }
    return resnorm;

}

