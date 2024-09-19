/**
 * \file   physics_utilities.h
 * \author akirby
 *
 * Created on August 21, 2018, 2:08 PM
 */

#ifndef PHYSICS_UTILITIES_H
#define PHYSICS_UTILITIES_H

/* header files */
#ifndef DOXYGEN_IGNORE
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#endif


#define BUFFER_SIZE 256 /**< Size of file name */


/** Display PHYSICS module message to terminal
 * 
 * @param [in] mpi_rank          mpi rank
 */
void utilities_write_physics_message(int mpi_rank);

/** Input file read line helper function: one integer
 * 
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] integer     integer variable to read
 * @return                  returns 0 if successfully read variable
 */
int utilities_find_keyword_integer(char* filename,char *keyword,
                                   int *integer);

/** Input file read line helper function: two integers
 * 
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] integer1    integer variable #1 to read
 * @param [out] integer2    integer variable #2 to read
 * @return                  returns 0 if successfully read variable
 */
int utilities_find_keyword_two_integers(char* filename,char *keyword, 
                                        int *integer1,
                                        int *integer2);

/** Input file read line helper function: three integers
 * 
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] integer1    integer variable #1 to read
 * @param [out] integer2    integer variable #2 to read
 * @param [out] integer3    integer variable #3 to read
 * @return                  returns 0 if successfully read variable
 */
int utilities_find_keyword_three_integers(char* filename,char *keyword, 
                                          int *integer1,
                                          int *integer2, 
                                          int *integer3);

/** Input file read line helper function: four integers
 * 
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] integer1    integer variable #1 to read
 * @param [out] integer2    integer variable #2 to read
 * @param [out] integer3    integer variable #3 to read
 * @param [out] integer4    integer variable #4 to read
 * @return                  returns 0 if successfully read variable
 */
int utilities_find_keyword_four_integers(char* filename,char *keyword, 
                                         int *integer1,
                                         int *integer2, 
                                         int *integer3,
                                         int *integer4);

/** Input file read line helper function: five integers
 * 
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] integer1    integer variable #1 to read
 * @param [out] integer2    integer variable #2 to read
 * @param [out] integer3    integer variable #3 to read
 * @param [out] integer4    integer variable #4 to read
 * @param [out] integer5    integer variable #5 to read
 * @return                  returns 0 if successfully read variable
 */
int utilities_find_keyword_five_integers(char* filename,char *keyword, 
                                         int *integer1,
                                         int *integer2, 
                                         int *integer3,
                                         int *integer4,
                                         int *integer5);

/** Input file read line helper function: one double
 * 
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] dbl         double precision variable to read
 * @return                  returns 0 if successfully read variable
 */
int utilities_find_keyword_double(char* filename, char *keyword,
                                  double *dbl);

/** Input file read line helper function: two doubles
 * 
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] dbl1        double precision variable #1 to read
 * @param [out] dbl2        double precision variable #2 to read
 * @return                  returns 0 if successfully read variable
 */
int utilities_find_keyword_two_doubles(char* filename,char *keyword,
                                       double *dbl1, 
                                       double *dbl2);

/** Input file read line helper function: three doubles
 * 
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] dbl1        double precision variable #1 to read
 * @param [out] dbl2        double precision variable #2 to read
 * @param [out] dbl3        double precision variable #3 to read
 * @return                  returns 0 if successfully read variable
 */
int utilities_find_keyword_three_doubles(char* filename,char *keyword,
                                         double *dbl1, 
                                         double *dbl2, 
                                         double *dbl3);

/** Input file read line helper function: four doubles
 * 
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] dbl1        double precision variable #1 to read
 * @param [out] dbl2        double precision variable #2 to read
 * @param [out] dbl3        double precision variable #3 to read
 * @param [out] dbl4        double precision variable #4 to read
 * @return                  returns 0 if successfully read variable
 */
int utilities_find_keyword_four_doubles(char* filename,char *keyword,
                                        double *dbl1, 
                                        double *dbl2, 
                                        double *dbl3,
                                        double *dbl4);

/** Input file read line helper function: five doubles
 * 
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] dbl1        double precision variable #1 to read
 * @param [out] dbl2        double precision variable #2 to read
 * @param [out] dbl3        double precision variable #3 to read
 * @param [out] dbl4        double precision variable #4 to read
 * @param [out] dbl5        double precision variable #5 to read
 * @return                  returns 0 if successfully read variable
 */
int utilities_find_keyword_five_doubles(char* filename,char *keyword,
                                        double *dbl1, 
                                        double *dbl2, 
                                        double *dbl3,
                                        double *dbl4,
                                        double *dbl5);

/** Input file read line helper function: string
 * 
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] string      string to read
 * @return                  returns 0 if successfully read variable
 */
int utilities_find_keyword_string(  char *filename, char *keyword, 
                                    char *string);


#endif /* PHYSICS_UTILITIES_H */

