/**
 * \file   physics_var_simulation.h
 * \ingroup physics_group
 * \author akirby
 *
 * \brief  User-specified simulation related data: 
 *         checkpoint/visualization/regrid flags, 
 *         time stepping flags.
 * 
 * Created on August 24, 2018, 2:55 PM
 */

#ifndef PHYSICS_VAR_SIMULATION_H
#define PHYSICS_VAR_SIMULATION_H

/**
 * physics_simulation_t contains user-specified 
 * simulation data and time marching flags
 */
typedef struct {
    
    int checkpoint_interval;    /**< User-specified solution checkpoint interval (read from input) */
    int visualization_interval; /**< User-specified solution visualization interval (read from input) */
    int regrid_interval;        /**< User-specified regrid interval (read from input) */
    
    int time_scheme;        /**< Time marching scheme */
    int time_steps;         /**< Number time steps in simulation */
    double dt;              /**< Time step size */
    
}
physics_simulation_t; /**< data type for the simulation data */


/**
 *  static primary storage unit for physics code base
 */
physics_simulation_t d_physics_simulation;

#endif /* PHYSICS_VAR_SIMULATION_H */

