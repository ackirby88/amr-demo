/**
 * \file   var_solution.h
 * \ingroup amr_group
 * \author akirby
 *
 * \brief   Solution checkpoint data information:
 *          checkpoint_interval.
 * 
 * Created on August 16, 2018, 5:05 PM
 */

#ifndef VAR_SOLUTION_H
#define VAR_SOLUTION_H


/**
 * solution_t contains simulation checkpoint related data.
 */
typedef struct {
    
    int checkpoint_interval;    /**< Interval to checkpoint solution */
    
}
solution_t; /**< Internal data type for solution data */


#endif /* VAR_SOLUTION_H */

