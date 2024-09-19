/**
 * \file   physics_var_physics.h
 * \ingroup physics_group
 * \author akirby
 * 
 * \brief   Physics-based data: 
 *          variables read from user input related to physics constants 
 *          are placed here.
 *
 * Created on August 15, 2018, 10:15 AM
 */

#ifndef VAR_PHYSICS_H
#define VAR_PHYSICS_H


/**
 * physics_fluid_t contains data related to fluid simulation (flow constants) 
 * and user-based amr tagging function flags.
 */
typedef struct {

    int initial_condition;  /**< Initial condition type */

    int bc[6];          /**< Boundary condition types:
                         *   * bc[0] = bc_xlo
                         *   * bc[1] = bc_xhi
                         *   * bc[2] = bc_ylo
                         *   * bc[3] = bc_yhi
                         *   * bc[4] = bc_zlo
                         *   * bc[5] = bc_zhi
                         */

    int amr_ntag_methods;           /**< Number of AMR tagging methods */
    int *amr_tag_methods;           /**< AMR tagging methods */
    double *amr_tag_tolerances ;    /**< AMR tagging tolerances */

    double cfl;         /**< CFL condition */
    double mach;        /**< Freestream Mach number */
    double alpha;       /**< Angle of attack */
    double beta;        /**< Angle of yaw */
    double gamma;       /**< Ratio of specific heats */
    double density;     /**< Freestream density */
    double pressure;    /**< Freestream pressure */

}
physics_fluid_t; /**< data type for the fluid physics data */

physics_fluid_t d_physics_fluid; /**< Primary physics fluid storage */


#endif /* VAR_PHYSICS_H */

