/**
 * \file   var_quad.h
 * \ingroup amr_group
 * \author akirby
 *
 * \brief   Quadrant solution data storage. 
 *          This is the data that is stored, communicated, adapted, visualized, 
 *          and checkpointed by p4est. NOTE: Keep minimum amount of data here.
 * 
 * Created on August 15, 2018, 12:40 PM
 */

#ifndef VAR_QUAD_H
#define VAR_QUAD_H


/* header files */
#include "var_defines.h"


/**
 * quad_data_t contains per qaudrant p4est related data.
 */
typedef struct {
    
    int tag;        /**< Grid refinement tag type indicator:
                     *   * Required for refinement and coarsening. 
                     */
    
    int type;       /**< Cell element type:
                     *  * Set by user as needed, e.g. cut-cell type 
                     */
    
    
#ifndef P4_TO_P8
    double soln[NFIELDS*NXPATCH*NXPATCH]; /**< Quadrant solution data:
                                           *    * Solution size:  
                                           * soln[(# of fields)*(# of 1D subgrid elements)^(spatial dimension)] 
                                           *    * NOTE: Size is required to be static constant  by p4est
                                           */
#else
    double soln[NFIELDS*NXPATCH*NXPATCH*NXPATCH]; /**< Quadrant solution data:
                                           *    * Solution size:  
                                           * soln[(# of fields)*(# of 1D subgrid elements)^(spatial dimension)] 
                                           *    * NOTE: Size is required to be static constant  by p4est
                                           */
#endif
    
}
quad_data_t; /**< Internal data type for per quadratant p4est data */


#endif /* VAR_QUAD_H */

