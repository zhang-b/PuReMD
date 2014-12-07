/*
 * =====================================================================================
 *
 *       Filename:  ndx.h
 *
 *    Description:  ndx.h file
 *
 *        Version:  1.0
 *        Created:  07/29/2014 11:49:48 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef __NDX_H_
#define __NDX_H_

#include "mytypes.h"

char Read_Ndx_File( FILE*, reax_groups*, control_params*, output_controls* );
char Make_Default_Groups( reax_groups*, control_params*, reax_system*,
                         output_controls* );

#endif
