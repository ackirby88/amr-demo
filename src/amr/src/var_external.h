/**
 * \file   var_external.h
 * \ingroup amr_group
 * \author akirby
 * 
 * \brief   External solver data related to mpi and amr information: 
 *          real and ghost cell counts, mpi communication schedule data, 
 *          cell/face/edge/node AMR relation information.
 *
 * Created on August 15, 2018, 10:42 AM
 */

#ifndef VAR_EXTERNAL_H
#define VAR_EXTERNAL_H


/**
 * external_t contains external communication related data.
 */
typedef struct {
    
    int ncell;          /**< Number of real cells in array */
    int nface;          /**< Number of faces in array*/
    int nghost;         /**< Number of ghost cells in array */
    int nfringe;        /**< Ghost fringe depth */
    int dof;            /**< Number of real degrees of freedom */
    
    int nsend;          /**< MPI send count */
    int nrecv;          /**< MPI receive count */
    int *send;          /**< MPI send buffer: [0] = send cell id */
    int *recv;          /**< MPI receive buffer */
    int *sizes;         /**< Number of fields in a cell sizes[cell_ind] */
    
    double *soln;       /**< External solution buffer */
    double *cell_geom;  /**< Grid element coordinates */
    
    int *cell_info;     /**< Cell information on element (3 fields per element)    
                         * + cell_info[0] = cell type                       
                         * + cell_info[1] = grid level                      
                         * + cell_info[2] = solution index                      
                         */
    
    int *face_info;     /**< Face information on element (11 fields per face)
                         * + face_info[0] = face type
                         *      * type == 1: boundary
                         *      * type == 2: full-full
                         *      * type == 3: 2D hanging: check side (face_info[1])
                         *      * type == 5: 3D hanging: check side (face_info[1])
                         *                   
                         * + face_info[1] = side of element
                         *      * side == 0: xlo face
                         *      | Left  | Right |
                         *      | ----: | :---- |
                         *      | hang  | full  |
                         * 
                         *      * side == 1: xhi face
                         *      | Left  | Right |
                         *      | ----: | :---- |
                         *      | full  | hang  |
                         * 
                         *      * side == 2: xlo face
                         *      | Left  | Right |
                         *      | ----: | :---- |
                         *      | hang  | full  |
                         * 
                         *      * side == 3: xhi face
                         *      | Left  | Right |
                         *      | ----: | :---- |
                         *      | full  | hang  |
                         * 
                         *      * side == 4: xlo face
                         *      | Left  | Right |
                         *      | ----: | :---- |
                         *      | hang  | full  |
                         * 
                         *      * side == 5: xhi face
                         *      | Left  | Right |
                         *      | ----: | :---- |
                         *      | full  | hang  |
                         *    
                         * 
                         * + face_info[2] = this cell index in cell_info
                         * + face_info[3] = neighbor 1 cell index in cell_info (if present)
                         * + face_info[4] = neighbor 2 cell index in cell_info (if present)
                         * + face_info[5] = neighbor 3 cell index in cell_info (if present)
                         * + face_info[6] = neighbor 4 cell index in cell_info (if present)
                         * + face_info[7] = neighbor 1 side (if present)
                         * + face_info[8] = neighbor 2 side (if present)
                         * + face_info[9] = neighbor 3 side (if present)
                         * + face_info[10]= neighbor 4 side (if present)
                         */
    
    int *edge_info;     /**< Edge information on element ({NOT SET} fields per face) 
                         * + edge_info[0] = {NOT SET}
                         */
    
    int *corner_info;   /**< Corner information on element ({NOT SET} fields per face) 
                         * + corner_info[0] = {NOT SET}
                         */
    
}
external_t; /**< Internal data type for external communication data */


#endif /* VAR_EXTERNAL_H */

