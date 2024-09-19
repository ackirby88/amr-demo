/**
 * \file   var_initialize.h
 * \ingroup amr_group
 * \author akirby
 *
 * \brief   Solver initialization data: 
 *          restart simulation flag & file, amr input file.
 * 
 * Created on August 15, 2018, 2:32 PM
 */

#ifndef VAR_INITIALIZE_H
#define VAR_INITIALIZE_H


/**
 * initialize_t contains user input related data.
 */
typedef struct {
    
    int  restart_flag;                  /**< Restart simulation indicator */
    char restart_file[BUFF_SIZE];       /**< Restart file name */
    char input_file[BUFF_SIZE];         /**< Input file name */
    
}
initialize_t; /**< Internal data type for user input data */


#endif /* VAR_INITIALIZE_H */

