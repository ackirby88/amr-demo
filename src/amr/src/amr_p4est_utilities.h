/**
 * \file   amr_p4est_utilities.h
 * \author akirby
 *
 * Created on August 17, 2018, 2:44 PM
 */

#ifndef AMR_P4EST_UTILITIES_H
#define AMR_P4EST_UTILITIES_H


/* header files */
#include "main.h"
#include "var_ctx.h"
#include "var_quad.h"

#include <mpi.h>

#ifndef P4_TO_P8
#include <p4est_extended.h>
#else
#include <p8est_extended.h>
#endif

/** Calculates amr quadrant/octant geometry coordinates
 * 
 * @param [in]  p4est           p4est tree structure
 * @param [in]  which_tree      p4est tree in forest of trees
 * @param [in]  q               p4est cell quadrant
 * @param [out] xyz             cell geometry coordinates
 * @param [in]  display_quad    flag to display the cell geometry coordinates
 */
void p4est_utilities_quad_coordinates(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *q, double *xyz, int display_quad);

/** Calculates amr quadrant/octant lower corner geometry coordinates
 * 
 * @param [in]  p4est           p4est tree structure
 * @param [in]  which_tree      p4est tree in forest of trees
 * @param [in]  q               p4est cell quadrant
 * @param [out] xyz             cell lower geometry coordinates
 */
void p4est_utilities_get_xlo(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *q, double xyz[3]);

/** Provides local id of cell volume
 * 
 * @param [in] info     p4est cell volume information
 * @return local id of quadrant
 */
p4est_locidx_t p4est_utilities_get_local_id_volume(p4est_iter_volume_info_t *info);

/** Provides local id of quad from hanging face
 * 
 * @param [in] info         p4est face information
 * @param [in] which_side   p4est which side face is on
 * @param [in] which_hang   p4est which face is hanging
 * @return local id of quadrant
 */
p4est_locidx_t p4est_utilities_get_local_id_face_hang(p4est_iter_face_info_t *info, int which_side, int which_hang);

/** Provides local id of quad from full face
 * 
 * @param [in] info         p4est face information
 * @param [in] which_side   p4est which side face is on
 * @return local id of quadrant
 */
p4est_locidx_t p4est_utilities_get_local_id_face_full(p4est_iter_face_info_t *info, int which_side);

/** Displays the mesh statistics to the terminal
 * 
 * @param [in] p4est        p4est tree structure
 */
void p4est_utilities_mesh_stats(p4est_t *p4est);


#endif /* AMR_P4EST_UTILITIES_H */

