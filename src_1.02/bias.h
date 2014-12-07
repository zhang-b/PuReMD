/*----------------------------------------------------------------------
  SeriallReax - Reax Force Field Simulator
      
  Copyright (2010) Purdue University
  Hasan Metin Aktulga, haktulga@cs.purdue.edu
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu
 
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of 
  the License, or (at your option) any later version.
               
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#ifndef __BIAS_H_
#define __BIAS_H_

#include "mytypes.h"

void Bias_Foo();
void Bias_COn_Combine(reax_system *, control_params *, simulation_data *, 
                 static_storage *, list **, output_controls *);
void Bias_COn_Decompose(reax_system *, control_params *, simulation_data *, 
                 static_storage *, list **, output_controls *);
void Bias_LJ_126(reax_system *, control_params *, simulation_data *, 
                 static_storage *, list **, output_controls *);
void Bias_Charge(reax_system *, control_params *, simulation_data *,
                 static_storage *, list **, output_controls *);
            
void Compute_AMD_Force(reax_system *, control_params *, simulation_data *,
                 static_storage *, list **);
void Compute_Bond_Boost_Force(reax_system *, control_params *, simulation_data *,
                 static_storage *, list **, output_controls *);
void Compute_Bond_Boost_Force_All(reax_system *, control_params *, simulation_data *,
                 static_storage *, list **, output_controls *);
void Compute_Bond_Boost_Force_All_Couple(reax_system *, control_params *, simulation_data *,
                 static_storage *, list **, output_controls *);
void Compute_Bond_Boost_Force_All_Couple_New(reax_system *, control_params *, simulation_data *,
                 static_storage *, list **, output_controls *);
#endif