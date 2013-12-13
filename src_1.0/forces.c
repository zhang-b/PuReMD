/*----------------------------------------------------------------------
 SerialReax - Reax Force Field Simulator

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

#include "forces.h"
#include "box.h"
#include "bond_orders.h"
#include "single_body_interactions.h"
#include "two_body_interactions.h"
#include "three_body_interactions.h"
#include "four_body_interactions.h"
#include "bias.h"
#include "list.h"
#include "print_utils.h"
#include "system_props.h"
#include "QEq.h"
#include "vector.h"
#include "math.h"

void Dummy_Interaction(reax_system *system, control_params *control,
		simulation_data *data, static_storage *workspace, list **lists,
		output_controls *out_control) {
}

void Init_Bonded_Force_Functions(control_params *control) {
	Interaction_Functions[0] = Calculate_Bond_Orders;
	Interaction_Functions[1] = Bond_Energy; //*/Dummy_Interaction;
	Interaction_Functions[2] = LonePair_OverUnder_Coordination_Energy;
	//*/Dummy_Interaction;
	Interaction_Functions[3] = Three_Body_Interactions; //*/Dummy_Interaction;
	Interaction_Functions[4] = Four_Body_Interactions; //*/Dummy_Interaction;
	if (control->hb_cut > 0)
		Interaction_Functions[5] = Hydrogen_Bonds; //*/Dummy_Interaction;
	else
		Interaction_Functions[5] = Dummy_Interaction;
	Interaction_Functions[6] = Dummy_Interaction; //empty
	Interaction_Functions[7] = Dummy_Interaction; //empty
	Interaction_Functions[8] = Dummy_Interaction; //empty
	Interaction_Functions[9] = Dummy_Interaction; //empty
}


void Compute_Bonded_Forces(reax_system *system, control_params *control,
		simulation_data *data, static_storage *workspace, list **lists,
		output_controls *out_control) {

	int i;
	// real t_start, t_end, t_elapsed;

#ifdef TEST_ENERGY
	/* Mark beginning of a new timestep in each energy file */
	fprintf( out_control->ebond, "step: %d\n%6s%6s%12s%12s%12s\n",
			data->step, "atom1", "atom2", "bo", "ebond", "total" );
	fprintf( out_control->elp, "step: %d\n%6s%12s%12s%12s\n",
			data->step, "atom", "nlp", "elp", "total" );
	fprintf( out_control->eov, "step: %d\n%6s%12s%12s\n",
			data->step, "atom", "eov", "total" );
	fprintf( out_control->eun, "step: %d\n%6s%12s%12s\n",
			data->step, "atom", "eun", "total" );
	fprintf( out_control->eval, "step: %d\n%6s%6s%6s%12s%12s%12s%12s%12s%12s\n",
			data->step, "atom1", "atom2", "atom3",
			"angle", "bo(12)", "bo(23)", "eval", "epen", "total" );
	fprintf( out_control->epen, "step: %d\n%6s%6s%6s%12s%12s%12s%12s%12s\n",
			data->step, "atom1", "atom2", "atom3",
			"angle", "bo(12)", "bo(23)", "epen", "total" );
	fprintf( out_control->ecoa, "step: %d\n%6s%6s%6s%12s%12s%12s%12s%12s\n",
			data->step, "atom1", "atom2", "atom3",
			"angle", "bo(12)", "bo(23)", "ecoa", "total" );
	fprintf( out_control->ehb, "step: %d\n%6s%6s%6s%12s%12s%12s%12s%12s\n",
			data->step, "atom1", "atom2", "atom3",
			"r(23)", "angle", "bo(12)", "ehb", "total" );
	fprintf( out_control->etor, "step: %d\n%6s%6s%6s%6s%12s%12s%12s%12s\n",
			data->step, "atom1", "atom2", "atom3", "atom4",
			"phi", "bo(23)", "etor", "total" );
	fprintf( out_control->econ, "step:%d\n%6s%6s%6s%6s%12s%12s%12s%12s%12s%12s\n",
			data->step, "atom1", "atom2", "atom3", "atom4",
			"phi", "bo(12)", "bo(23)", "bo(34)", "econ", "total" );
#endif 

	/* Implement all the function calls as function pointers */
	for (i = 0; i < NO_OF_INTERACTIONS; i++) {
		(Interaction_Functions[i])(system, control, data, workspace, lists,
				out_control);
#if defined(DEBUG_FOCUS)
		fprintf( stderr, "f%d-", i );
#endif
#ifdef TEST_FORCES
		(Print_Interactions[i])(system, control, data, workspace,
				lists, out_control);
#endif
	}
}

void Compute_NonBonded_Forces(reax_system *system, control_params *control,
		simulation_data *data, static_storage *workspace, list** lists,
		output_controls *out_control) {
	real t_start, t_elapsed;
#ifdef TEST_ENERGY
	fprintf( out_control->evdw, "step: %d\n%6s%6s%12s%12s%12s\n",
			data->step, "atom1", "atom2", "r12", "evdw", "total" );
	fprintf( out_control->ecou, "step: %d\n%6s%6s%12s%12s%12s%12s%12s\n",
			data->step, "atom1", "atom2", "r12", "q1", "q2", "ecou", "total" );
#endif

	t_start = Get_Time();
        if (control->qeq)
	    QEq(system, control, data, workspace, lists[FAR_NBRS], out_control);
	t_elapsed = Get_Timing_Info(t_start);
	data->timing.QEq += t_elapsed;
#if defined(DEBUG_FOCUS)
	fprintf( stderr, "qeq - " );
#endif

	if (control->tabulate == 0)
		vdW_Coulomb_Energy(system, control, data, workspace, lists, out_control);
	else
		Tabulated_vdW_Coulomb_Energy(system, control, data, workspace, lists,
				out_control);
#if defined(DEBUG_FOCUS)
	fprintf( stderr, "nonb forces - " );
#endif

#ifdef TEST_FORCES
	Print_vdW_Coulomb_Forces( system, control, data, workspace,
			lists, out_control );
#endif
}

/* This version of Compute_Total_Force computes forces from coefficients 
 accumulated by all interaction functions. Saves enormous time & space! */
void Compute_Total_Force(reax_system *system, control_params *control,
		simulation_data *data, static_storage *workspace, list **lists) {
	int i, pj;
	list *bonds = (*lists) + BONDS;

	for (i = 0; i < system->N; ++i)
		for (pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj)
			if (i < bonds->select.bond_list[pj].nbr) {
				if (control->ensemble == NVE || control->ensemble == NVT)
					Add_dBond_to_Forces(i, pj, system, data, workspace, lists);
				else
					Add_dBond_to_Forces_NPT(i, pj, system, data, workspace,
							lists);
			}
}

void Compute_AMD_Force(reax_system *system, control_params *control,
		simulation_data *data, static_storage *workspace, list **lists) {
	/**
	 * compute the AMD energy and force.
	 */
	int i, type_i;
	real q;
	real E_Pol, E_Pot;
	real delta_Pot, fscale;
	real s1, s2, s3, s4;
    // initiate parameters
    delta_Pot = 0.0;
    fscale = 0.0;
	/* Compute Potential Energy */
	E_Pol = 0.0;
	for (i = 0; i < system->N; i++) {
		q = system->atoms[i].q;
		type_i = system->atoms[i].type;

		E_Pol += (system->reaxprm.sbp[type_i].chi * q
				+ (system->reaxprm.sbp[type_i].eta / 2.0) * SQR(q))
				* KCALpMOL_to_EV;
	}

	E_Pot = data->E_BE + data->E_Ov + data->E_Un + data->E_Lp + data->E_Ang
			+ data->E_Pen + data->E_Coa + data->E_HB + data->E_Tor
			+ data->E_Con + data->E_vdW + data->E_Ele + data->E_Pol;
	/* compute delta V and Force scaling factor
	 * reference: JCTC 2011, 7, 890–897
	 */
	data->E_amd = control->amd_energy - data->Fragment_wat*120.0;
	if (data->E_amd > E_Pot) {
		if (control->amd_func == 1){
		/* s1 = Eb -V(r) */
		s1 = data->E_amd - E_Pot;
		/* s2 = (Eb -V(r))^p + a */
		s2 = POW(s1, control->amd_power) + control->amd_alpha;
		s3 = s1/s2/s2;
		s4 = control->amd_power*POW(s1, control->amd_power)-2*s2;
		delta_Pot = s1 * s1 / s2;
		fscale = 1 + s3*s4;
		}
		else if (control->amd_func == 2){
			s1 = data->E_amd - E_Pot;
			s2 = s1*s1;
			s3 = s2 + 2*s1;
			s4 = control->amd_alpha*(s3+1);
			delta_Pot = s2*(s1+1)/s4;
			fscale = 1-(s2+2*s1)/s4;
		}

		data->E_amd_delta = delta_Pot;
		data->F_amd_scale = fscale;
		/* scale the force */
		for (i = 0; i < system->N; ++i)
			rvec_Scale(system->atoms[i].f, fscale, system->atoms[i].f);
	}
	else
	{
		data->E_amd_delta = 0;
		data->F_amd_scale = 1;
	}
	//printf("in AMD no of water is %12.4f\n", data->F_amd_scale);
	//printf("in AMD no of water is %d\n", data->Fragment_wat);
}

int Find_Radicals(reax_system *system, control_params *control,
		simulation_data *data, static_storage *workspace, list **lists, 
        output_controls *out_control) {
  int i, no, nu;
  real ov, un;
  char *atp;
  nu = 0;
  no = 0;
  for( i=0; i < system->N; ++i ) {
      atp = system->reaxprm.sbp[ system->atoms[i].type ].name;
      ov = system->atoms[i].ov;
      un = system->atoms[i].un;
      if (strcmp(atp, "H") == 0) {
          if (un < -0.8)
              nu += 1;
      }
      else if (strcmp(atp, "O") == 0) {
          if (un < -0.8)
              nu += 1;
          if (ov > 0.3)
              no += 1;
      }
      else if (strcmp(atp, "C") == 0) {
          if (un < -0.8)
              nu += 1;
          if (ov > 0.4)
              no += 1;
      }
  }
  return nu + no;
}

void Compute_Bond_Boost_Force(reax_system *system, control_params *control,
		simulation_data *data, static_storage *workspace, list **lists, 
        output_controls *out_control) {
  int i, j, pj;
  int type_i, type_j;
  int adatom, adatom2; // label the boost atom
  int nbond, nrad; // Nb
  int start_i, end_i;
  real e, emax, r, re, r_max; // eta, eta_max, r, r_e
  real bo;
  real A, dA, V; // A(\eta^max), and \Delta A(\eta^max)
  real *rv, *rv_max;
  real q, P1, vmax; // q, P1 in equation 13
  real S1, S2, f1, f2, C1, C2;
  real bf; // boost force scale
  rvec df; // boost force
  reax_atom *atom1, *atom2;
  bond_order_data *bo_ij;
  two_body_parameters *twbp;
  list *bonds;

  nrad = Find_Radicals(system, control, data, workspace, lists, out_control);
  //printf("-------------------------step %d  -----------------\n", data->step);
  bonds = (*lists) + BONDS;

  // bond boost parameters
  q = control->bboost_q; // should read this from control file
  P1 = control->bboost_P1; 
  vmax = control->bboost_Vmax;

  // initiate parameters
  e = 0;
  adatom = 0;
  adatom2 = 0;
  nbond = 0;
  emax = 0.0;
  bo = 0.0;
  re = 0.001;

  A = 0.0;
  dA = 0.0;
  r_max = 0.0;

  // first get the max bond order
  for( i=0; i < system->N; ++i ) {
    re = 0.001;
    start_i = Start_Index(i, bonds);
    end_i = End_Index(i, bonds);
    for( pj = start_i; pj < end_i; ++pj ){
      if( i < bonds->select.bond_list[pj].nbr ) {
        j = bonds->select.bond_list[pj].nbr;
        type_i = system->atoms[i].type;
        type_j = system->atoms[j].type;
        twbp = &( system->reaxprm.tbp[type_i][type_j] );
        bo_ij = &( bonds->select.bond_list[pj].bo_data );
        re = twbp->r_e; // get r_e from ffield.ext
        r = bonds->select.bond_list[pj].d;
        e = (r - re)/re; // eta
        //printf("r = %f, re = %f, e = %f\n", r, re, e);
        bo = bo_ij->BO;
        if ( (fabs (emax) < fabs(e)) & (bo > 0.3)) {
            emax = e;
            adatom = i;
            adatom2 = j;
        }
      }
    }
  }
  /*
  printf("atom1 = %s ", system->reaxprm.sbp[ system->atoms[adatom].type ].name);
  printf("atom2 = %s ", system->reaxprm.sbp[ system->atoms[adatom2].type ].name);
  printf("bo = %f\n", bo);
  */

  V = 0.0; // boost energy

  if (fabs(emax) < q && nrad == 0) {
      i = adatom;
      vmax = control->bboost_Vmax;

      // calculate A, and dA
      S1 = emax/q;
      S2 = S1 * S1;
      C1 = 1 - (emax/q)*(emax/q);
      C2 = 1 - P1*P1*S2;
      A = C1 * C1 / C2;

      /*
      printf("S1 = %f, S2 = %f ", S1, S2);
      printf("C1 = %f, C2 = %f ", C1, C2);
      printf("A = %f, emax = %f\n", A, emax);
      */

      start_i = Start_Index(i, bonds);
      end_i = End_Index(i, bonds);
      for( pj = start_i; pj < end_i; ++pj ) {
        if( i < bonds->select.bond_list[pj].nbr ) {
            bo_ij = &( bonds->select.bond_list[pj].bo_data );
            bo = bo_ij->BO;
            if (bo > 0.3)
                nbond += 1;
        }
      }

      for( pj = start_i; pj < end_i; ++pj ) {
        if( i < bonds->select.bond_list[pj].nbr ) {
          bo_ij = &( bonds->select.bond_list[pj].bo_data );
          bo = bo_ij->BO;
          if (bo > 0.3) {
            j = bonds->select.bond_list[pj].nbr;
            r = bonds->select.bond_list[pj].d;
            rv = bonds->select.bond_list[pj].dvec;
            type_i = system->atoms[i].type;
            type_j = system->atoms[j].type;
            twbp = &( system->reaxprm.tbp[type_i][type_j] );
            vmax = twbp->v_max;
            re = twbp->r_e; // get r_e from ffield.ext

            e = (r - re)/re; // eta
            C1 = 1 - e*e/(q*q);
            V += vmax / nbond * C1;
            bf = 2*A*vmax*e/(nbond*q*re);

            /*
            printf("e = %.3f, emax = %.3f, C1 = %.3f, re = %.3f, V = %.3f, A = %.3f bf = %.3f, vmax = %.3f, \n",\
                    e, emax, C1, re, V, A, bf, vmax);
            */

            if (fabs(e - emax) < 0.00001) {
                C2 = 1 - e*e*P1*P1/(q*q);
                f1 = 2/(q*re);
                f2 = f1/2*P1*P1;
                dA = 2*e*C1/C2*(f1 - f2*C1/C2);
                rv_max = bonds->select.bond_list[pj].dvec;
                r_max = bonds->select.bond_list[pj].d;;
            }
            rvec_Scale(df, bf, rv);
            //rvec_Scale(df, 1/r, df);
            // assign the boost force
            atom1 = &( system->atoms[i] );
            atom2 = &( system->atoms[j] );
            /* printf("ofx = %8.3f ofy = %8.3f ofz = %8.3f \n", system->atoms[i].f[0],\
                    system->atoms[i].f[1], system->atoms[i].f[2]);
            */
            rvec_Add(system->atoms[i].f, df);
            //printf(" fx = %8.3f  fy = %8.3f  fz = %8.3f \n", df[0], df[1], df[2]);
            /*
            rvec_Scale(df, -1, df);
            rvec_Add(system->atoms[j].f, df);
            */
          }
        }
      }
      // add the contribution of evalope function A
      bf = dA * V;
      rvec_Scale(df, bf, rv_max);
      //rvec_Scale(df, 1/r_max, df);
      //printf("Max bf = %8.3f\n", bf);
      /*
      printf("Mfx = %8.3f Mfy = %8.3f Mfz = %8.3f dA = %8.3f V = %8.3f emax = %8.3f \n",\
             df[0], df[1], df[2], dA, V, emax);
      */
      atom1 = &( system->atoms[i] );
      rvec_Add(system->atoms[i].f, df);
  }

  fprintf( out_control->bboost, "%-10d%6d%6d%10.4f%10.4f%10.4f", \
  data->step, adatom + 1, adatom2 + 1, bo, emax, A*V );

  fprintf( out_control->bboost, " %4s %4s\n", \
  system->reaxprm.sbp[ system->atoms[adatom].type ].name, 
  system->reaxprm.sbp[ system->atoms[adatom2].type ].name);
  fflush( out_control->bboost);
}

void Compute_Bond_Boost_Force_All(reax_system *system, control_params *control,
		simulation_data *data, static_storage *workspace, list **lists, 
        output_controls *out_control) {
  int i, j, pj;
  int type_i, type_j;
  int adatom, adatom2; // label the boost atom
  int nbond, nrad; // Nb
  int start_i, end_i;
  real e, emax, r, re, r_max; // eta, eta_max, r, r_e
  real bo;
  real A, dA, V; // A(\eta^max), and \Delta A(\eta^max)
  real *rv, *rv_max;
  real q, P1, vmax; // q, P1 in equation 13
  real S1, S2, f1, f2, C1, C2;
  real T;
  real bf, bfactor; // boost force scale
  rvec df; // boost force
  reax_atom *atom1, *atom2;
  bond_order_data *bo_ij;
  two_body_parameters *twbp;
  list *bonds;

  nrad = Find_Radicals(system, control, data, workspace, lists, out_control);
  //printf("-------------------------step %d  -----------------\n", data->step);
  bonds = (*lists) + BONDS;

  // bond boost parameters
  q = control->bboost_q; // should read this from control file
  P1 = control->bboost_P1; 
  vmax = control->bboost_Vmax;

  // initiate parameters
  e = 0;
  adatom = 0;
  adatom2 = 0;
  bo = 0.0;
  re = 0.001;

  A = 0.0;
  dA = 0.0;
  r_max = 0.0;
  bfactor = 0.0;
  T = control->T_final;

  // first get the max bond order
  for( i=0; i < system->N; ++i ) {
    re = 0.001;
    A = 0.0;
    V = 0.0; // bost energy
    nbond = 0;
    emax = 0.0;
    start_i = Start_Index(i, bonds);
    end_i = End_Index(i, bonds);
    for( pj = start_i; pj < end_i; ++pj ){
      if( i < bonds->select.bond_list[pj].nbr ) {
        j = bonds->select.bond_list[pj].nbr;
        type_i = system->atoms[i].type;
        type_j = system->atoms[j].type;
        twbp = &( system->reaxprm.tbp[type_i][type_j] );
        bo_ij = &( bonds->select.bond_list[pj].bo_data );
        re = twbp->r_e; // get r_e from ffield.ext
        r = bonds->select.bond_list[pj].d;
        e = (r - re)/re; // eta
        bo = bo_ij->BO;
        if (bo > 0.3) {
            nbond += 1;
            if ( fabs(emax) < fabs(e) ) {
                emax = e;
            }
        }
      }
    }

    if (fabs(emax) < q && nbond > 0 && nrad == 0) {
      vmax = control->bboost_Vmax;

      // calculate A, and dA
      S1 = emax/q;
      S2 = S1 * S1;
      C1 = 1 - (emax/q)*(emax/q);
      C2 = 1 - P1*P1*S2;
      A = C1 * C1 / C2;

      for( pj = start_i; pj < end_i; ++pj ) {
        if( i < bonds->select.bond_list[pj].nbr ) {
          bo_ij = &( bonds->select.bond_list[pj].bo_data );
          bo = bo_ij->BO;
          if (bo > 0.3) {
            j = bonds->select.bond_list[pj].nbr;
            r = bonds->select.bond_list[pj].d;
            rv = bonds->select.bond_list[pj].dvec;
            type_i = system->atoms[i].type;
            type_j = system->atoms[j].type;
            twbp = &( system->reaxprm.tbp[type_i][type_j] );
            vmax = twbp->v_max;
            re = twbp->r_e; // get r_e from ffield.ext

            e = (r - re)/re; // eta
            C1 = 1 - e*e/(q*q);
            V += vmax / nbond * C1;
            bf = 2*A*vmax*e/(nbond*q*re);

            if (fabs(e - emax) < 0.00001) {
                C2 = 1 - e*e*P1*P1/(q*q);
                f1 = 2/(q*re);
                f2 = f1/2*P1*P1;
                dA = 2*e*C1/C2*(f1 - f2*C1/C2);
                rv_max = bonds->select.bond_list[pj].dvec;
                r_max = bonds->select.bond_list[pj].d;;
            }
            rvec_Scale(df, bf, rv);
            atom1 = &( system->atoms[i] );
            atom2 = &( system->atoms[j] );
            rvec_Add(system->atoms[i].f, df);
          }
        }
      }
      // add the contribution of evalope function A
      bf = dA * V;
      rvec_Scale(df, bf, rv_max);
      atom1 = &( system->atoms[i] );
      rvec_Add(system->atoms[i].f, df);
    }
    bfactor += exp(A*V/(T * 8.314 / 4184));
  }
  fprintf( out_control->bboost, "%-10d%6d%6d%10.4f%10.4f%10.4f", \
  data->step, adatom + 1, adatom2 + 1, bo, emax, bfactor );

  fprintf( out_control->bboost, " %4s %4s\n", \
  system->reaxprm.sbp[ system->atoms[adatom].type ].name, 
  system->reaxprm.sbp[ system->atoms[adatom2].type ].name);
  fflush( out_control->bboost);
}

void Compute_Bond_Boost_Force_All_Couple(reax_system *system, control_params *control,
		simulation_data *data, static_storage *workspace, list **lists, 
        output_controls *out_control) {
  int i, j, pj;
  int type_i, type_j;
  int adatom, adatom2; // label the boost atom
  int nbond, nrad; // Nb
  int start_i, end_i;
  real e, emax, r, re, r_max; // eta, eta_max, r, r_e
  real bo;
  real A, dA, V; // A(\eta^max), and \Delta A(\eta^max)
  real *rv, *rv_max;
  real q, P1, vmax; // q, P1 in equation 13
  real S1, S2, f1, f2, C1, C2;
  real T;
  real bf, bfactor; // boost force scale
  rvec df; // boost force
  reax_atom *atom1, *atom2;
  bond_order_data *bo_ij;
  two_body_parameters *twbp;
  list *bonds;

  nrad = Find_Radicals(system, control, data, workspace, lists, out_control);
  //nrad = 0;
  //printf("-------------------------step %d  -----------------\n", data->step);
  bonds = (*lists) + BONDS;

  // bond boost parameters
  q = control->bboost_q; // should read this from control file
  P1 = control->bboost_P1; 
  vmax = control->bboost_Vmax;

  // initiate parameters
  e = 0;
  adatom = 0;
  adatom2 = 0;
  bo = 0.0;
  re = 0.001;

  A = 0.0;
  dA = 0.0;
  r_max = 0.0;
  bfactor = 1.0;
  T = control->T_final;
  V = 0.0; // bost energy
  nbond = 0;
  emax = 0.0;

  // first get the max bond order
  for( i=0; i < system->N; ++i ) {
    start_i = Start_Index(i, bonds);
    end_i = End_Index(i, bonds);
    for( pj = start_i; pj < end_i; ++pj ){
      if( i < bonds->select.bond_list[pj].nbr ) {
        j = bonds->select.bond_list[pj].nbr;
        type_i = system->atoms[i].type;
        type_j = system->atoms[j].type;
        twbp = &( system->reaxprm.tbp[type_i][type_j] );
        vmax = control->bboost_Vmax;
        vmax = twbp->v_max;
        bo_ij = &( bonds->select.bond_list[pj].bo_data );
        re = twbp->r_e; // get r_e from ffield.ext
        r = bonds->select.bond_list[pj].d;
        e = (r - re)/re; // eta
        bo = bo_ij->BO;
        if (bo > 0.3 && vmax > 0) {
            nbond += 1;
            if ( fabs(emax) < fabs(e) ) {
                emax = e;
                adatom = i;
                adatom2 = j;
                /* for debug
                printf("i = %d, j = %d, re = %.2f, r = %.2f ", i, j, re, r);
                printf("x1 = %.2f, y1 = %.2f, z1 = %.2f, ", \
                system->atoms[i].x[0], system->atoms[i].x[1], system->atoms[i].x[2]);
                printf("x2 = %.2f, y2 = %.2f, z2 = %.2f\n", \
                system->atoms[j].x[0], system->atoms[j].x[1], system->atoms[j].x[2]);
                printf( " %4s %4s\n", \
                system->reaxprm.sbp[ system->atoms[i].type ].name, 
                system->reaxprm.sbp[ system->atoms[j].type ].name);
                */
                rv_max = bonds->select.bond_list[pj].dvec;
                r_max = bonds->select.bond_list[pj].d;;
            }
        }
      }
    }
  }
  int ntmp; 
  if (fabs(emax) < q && nbond > 0 && nrad == 0) {
    data->boost ++ ;
    // calculate A, and dA
    S1 = emax/q;
    S2 = S1 * S1;
    C1 = 1 - (emax/q)*(emax/q);
    C2 = 1 - P1*P1*S2;
    A = C1 * C1 / C2;

    ntmp = 0;
    for( i=0; i < system->N; ++i ) {
      start_i = Start_Index(i, bonds);
      end_i = End_Index(i, bonds);
      for( pj = start_i; pj < end_i; ++pj ) {
        if( i < bonds->select.bond_list[pj].nbr ) {
          j = bonds->select.bond_list[pj].nbr;
          type_i = system->atoms[i].type;
          type_j = system->atoms[j].type;
          twbp = &( system->reaxprm.tbp[type_i][type_j] );
          vmax = control->bboost_Vmax;
          vmax = twbp->v_max;
          bo_ij = &( bonds->select.bond_list[pj].bo_data );
          r = bonds->select.bond_list[pj].d;
          re = twbp->r_e; // get r_e from ffield.ext
          rv = bonds->select.bond_list[pj].dvec;
          bo = bo_ij->BO;
          if (bo > 0.3 && vmax > 0) {

            e = (r - re)/re; // eta
            C1 = 1 - e*e/(q*q);
            V += vmax / nbond * C1;
            bf = 2*A*vmax*e/(nbond*q*re);

            if (fabs(e - emax) < 0.00001) {
                C2 = 1 - e*e*P1*P1/(q*q);
                f1 = 2/(q*re);
                f2 = f1/2*P1*P1;
                dA = 2*e*C1/C2*(f1 - f2*C1/C2);
            }
            rvec_Scale(df, bf, rv);
            atom1 = &( system->atoms[i] );
            atom2 = &( system->atoms[j] );
            rvec_Add(system->atoms[i].f, df);
            ntmp++;
          }
        }
      }
    }
    // add the contribution of evalope function A
    bf = dA * V;
    rvec_Scale(df, bf, rv_max);
    atom1 = &( system->atoms[adatom]);
    rvec_Add(system->atoms[adatom].f, df);
  }
  else {
    data->boost = 0;
  }
  bfactor = exp(4184 * A * V/(T * 8.314));
  //bfactor = V;
  fprintf( out_control->bboost, "%-10d%6d%6d%6d%3d%8.4f%8.4f %24.4f", \
  data->step, nbond, adatom, adatom2, nrad, r_max, emax, bfactor );

  fprintf( out_control->bboost, " %4s %4s\n", \
  system->reaxprm.sbp[ system->atoms[adatom].type ].name, 
  system->reaxprm.sbp[ system->atoms[adatom2].type ].name);
  fflush( out_control->bboost);
}

void Validate_Lists(static_storage *workspace, list **lists, int step, int n,
		int Hmax, int Htop, int num_bonds, int num_hbonds) {
	int i, flag;
	list *bonds, *hbonds;

	bonds = *lists + BONDS;
	hbonds = *lists + HBONDS;

	/* far neighbors */
	if (Htop > Hmax * DANGER_ZONE) {
		workspace->realloc.Htop = Htop;
		if (Htop > Hmax) {
			fprintf(stderr,
					"step%d - ran out of space on H matrix: Htop=%d, max = %d",
					step, Htop, Hmax);
			exit(INSUFFICIENT_SPACE);
		}
	}

	/* bond list */
	flag = -1;
	workspace->realloc.num_bonds = num_bonds;
	for (i = 0; i < n - 1; ++i)
		if (End_Index(i, bonds) >= Start_Index(i + 1, bonds) - 2) {
			workspace->realloc.bonds = 1;
			if (End_Index(i, bonds) > Start_Index(i + 1, bonds))
				flag = i;
		}

	if (flag > -1) {
		fprintf(stderr, "step%d-bondchk failed: i=%d end(i)=%d str(i+1)=%d\n",
				step, flag, End_Index(flag, bonds),
				Start_Index(flag + 1, bonds));
		exit(INSUFFICIENT_SPACE);
	}

	if (End_Index(i, bonds) >= bonds->num_intrs - 2) {
		workspace->realloc.bonds = 1;

		if (End_Index(i, bonds) > bonds->num_intrs) {
			fprintf(stderr,
					"step%d-bondchk failed: i=%d end(i)=%d bond_end=%d\n",
					step, flag, End_Index(i, bonds), bonds->num_intrs);
			exit(INSUFFICIENT_SPACE);
		}
	}

	/* hbonds list */
	if (workspace->num_H > 0) {
		flag = -1;
		workspace->realloc.num_hbonds = num_hbonds;
		for (i = 0; i < workspace->num_H - 1; ++i)
			if (Num_Entries(i, hbonds) >= (Start_Index(i + 1, hbonds)
					- Start_Index(i, hbonds)) * DANGER_ZONE) {
				workspace->realloc.hbonds = 1;
				if (End_Index(i, hbonds) > Start_Index(i + 1, hbonds))
					flag = i;
			}

		if (flag > -1) {
			fprintf(stderr,
					"step%d-hbondchk failed: i=%d end(i)=%d str(i+1)=%d\n",
					step, flag, End_Index(flag, hbonds), Start_Index(flag + 1,
							hbonds));
			exit(INSUFFICIENT_SPACE);
		}

		if (Num_Entries(i, hbonds) >= (hbonds->num_intrs - Start_Index(i,
				hbonds)) * DANGER_ZONE) {
			workspace->realloc.hbonds = 1;

			if (End_Index(i, hbonds) > hbonds->num_intrs) {
				fprintf(stderr,
						"step%d-hbondchk failed: i=%d end(i)=%d hbondend=%d\n",
						step, flag, End_Index(i, hbonds), hbonds->num_intrs);
				exit(INSUFFICIENT_SPACE);
			}
		}
	}
}

void Init_Forces(reax_system *system, control_params *control,
		simulation_data *data, static_storage *workspace, list **lists,
		output_controls *out_control) {
	int i, j, pj;
	int start_i, end_i;
	int type_i, type_j;
	int Htop, btop_i, btop_j, num_bonds, num_hbonds;
	int ihb, jhb, ihb_top, jhb_top;
	int flag;
	real r_ij, r2, self_coef;
	real dr3gamij_1, dr3gamij_3, Tap;
	//real val, dif, base;
	real C12, C34, C56;
	real Cln_BOp_s, Cln_BOp_pi, Cln_BOp_pi2;
	real BO, BO_s, BO_pi, BO_pi2;
	real p_boc1, p_boc2;
	sparse_matrix *H;
	list *far_nbrs, *bonds, *hbonds;
	single_body_parameters *sbp_i, *sbp_j;
	two_body_parameters *twbp;
	far_neighbor_data *nbr_pj;
	//LR_lookup_table *t;
	reax_atom *atom_i, *atom_j;
	bond_data *ibond, *jbond;
	bond_order_data *bo_ij, *bo_ji;

	far_nbrs = *lists + FAR_NBRS;
	bonds = *lists + BONDS;
	hbonds = *lists + HBONDS;

	H = workspace->H;
	Htop = 0;
	num_bonds = 0;
	num_hbonds = 0;
	btop_i = btop_j = 0;
	p_boc1 = system->reaxprm.gp.l[0];
	p_boc2 = system->reaxprm.gp.l[1];

	for (i = 0; i < system->N; ++i) {
		atom_i = &(system->atoms[i]);
		type_i = atom_i->type;
		start_i = Start_Index(i, far_nbrs);
		end_i = End_Index(i, far_nbrs);
		H->start[i] = Htop;
		btop_i = End_Index(i, bonds);
		sbp_i = &(system->reaxprm.sbp[type_i]);
		ihb = ihb_top = -1;
		if (control->hb_cut > 0 && (ihb = sbp_i->p_hbond) == 1)
			ihb_top = End_Index(workspace->hbond_index[i], hbonds);

		for (pj = start_i; pj < end_i; ++pj) {
			nbr_pj = &(far_nbrs->select.far_nbr_list[pj]);
			j = nbr_pj->nbr;
			atom_j = &(system->atoms[j]);

			flag = 0;
			if ((data->step - data->prev_steps) % control->reneighbor == 0) {
				if (nbr_pj->d <= control->r_cut)
					flag = 1;
				else
					flag = 0;
			} else if ((nbr_pj->d = Sq_Distance_on_T3(atom_i->x, atom_j->x,
					&(system->box), nbr_pj->dvec)) <= SQR(control->r_cut)) {
				nbr_pj->d = sqrt(nbr_pj->d);
				flag = 1;
			}

			if (flag) {
				type_j = system->atoms[j].type;
				r_ij = nbr_pj->d;
				sbp_j = &(system->reaxprm.sbp[type_j]);
				twbp = &(system->reaxprm.tbp[type_i][type_j]);
				self_coef = (i == j) ? 0.5 : 1.0;

				/* H matrix entry */
				Tap = control->Tap7 * r_ij + control->Tap6;
				Tap = Tap * r_ij + control->Tap5;
				Tap = Tap * r_ij + control->Tap4;
				Tap = Tap * r_ij + control->Tap3;
				Tap = Tap * r_ij + control->Tap2;
				Tap = Tap * r_ij + control->Tap1;
				Tap = Tap * r_ij + control->Tap0;

				dr3gamij_1 = (r_ij * r_ij * r_ij + twbp->gamma);
				dr3gamij_3 = POW(dr3gamij_1, 0.33333333333333);

				H->entries[Htop].j = j;
				H->entries[Htop].val = self_coef * Tap * EV_to_KCALpMOL
						/ dr3gamij_3;
				++Htop;

				/* hydrogen bond lists */
				if (control->hb_cut > 0 && (ihb == 1 || ihb == 2) && nbr_pj->d
						<= control->hb_cut) {
					// fprintf( stderr, "%d %d\n", atom1, atom2 );
					jhb = sbp_j->p_hbond;
					if (ihb == 1 && jhb == 2) {
						hbonds->select.hbond_list[ihb_top].nbr = j;
						hbonds->select.hbond_list[ihb_top].scl = 1;
						hbonds->select.hbond_list[ihb_top].ptr = nbr_pj;
						++ihb_top;
						++num_hbonds;
					} else if (ihb == 2 && jhb == 1) {
						jhb_top = End_Index(workspace->hbond_index[j], hbonds);
						hbonds->select.hbond_list[jhb_top].nbr = i;
						hbonds->select.hbond_list[jhb_top].scl = -1;
						hbonds->select.hbond_list[jhb_top].ptr = nbr_pj;
						Set_End_Index(workspace->hbond_index[j], jhb_top + 1,
								hbonds);
						++num_hbonds;
					}
				}

				/* uncorrected bond orders */
				if (far_nbrs->select.far_nbr_list[pj].d <= control->nbr_cut) {
					r2 = SQR(r_ij);

					if (sbp_i->r_s > 0.0 && sbp_j->r_s > 0.0) {
						C12 = twbp->p_bo1 * POW(r_ij / twbp->r_s, twbp->p_bo2);
						BO_s = (1.0 + control->bo_cut) * EXP(C12);
					} else
						BO_s = C12 = 0.0;

					if (sbp_i->r_pi > 0.0 && sbp_j->r_pi > 0.0) {
						C34 = twbp->p_bo3 * POW(r_ij / twbp->r_p, twbp->p_bo4);
						BO_pi = EXP(C34);
					} else
						BO_pi = C34 = 0.0;

					if (sbp_i->r_pi_pi > 0.0 && sbp_j->r_pi_pi > 0.0) {
						C56 = twbp->p_bo5 * POW(r_ij / twbp->r_pp, twbp->p_bo6);
						BO_pi2 = EXP(C56);
					} else
						BO_pi2 = C56 = 0.0;

					/* Initially BO values are the uncorrected ones, page 1 */
					BO = BO_s + BO_pi + BO_pi2;

					if (BO >= control->bo_cut) {
						num_bonds += 2;
						/****** bonds i-j and j-i ******/
						ibond = &(bonds->select.bond_list[btop_i]);
						btop_j = End_Index(j, bonds);
						jbond = &(bonds->select.bond_list[btop_j]);

						ibond->nbr = j;
						jbond->nbr = i;
						ibond->d = r_ij;
						jbond->d = r_ij;
						rvec_Copy(ibond->dvec, nbr_pj->dvec);
						rvec_Scale(jbond->dvec, -1, nbr_pj->dvec);
						ivec_Copy(ibond->rel_box, nbr_pj->rel_box);
						ivec_Scale(jbond->rel_box, -1, nbr_pj->rel_box);
						ibond->dbond_index = btop_i;
						jbond->dbond_index = btop_i;
						ibond->sym_index = btop_j;
						jbond->sym_index = btop_i;
						++btop_i;
						Set_End_Index(j, btop_j + 1, bonds);

						bo_ij = &(ibond->bo_data);
						bo_ji = &(jbond->bo_data);
						bo_ji->BO = bo_ij->BO = BO;
						bo_ji->BO_s = bo_ij->BO_s = BO_s;
						bo_ji->BO_pi = bo_ij->BO_pi = BO_pi;
						bo_ji->BO_pi2 = bo_ij->BO_pi2 = BO_pi2;

						/* Bond Order page2-3, derivative of total bond order prime */
						Cln_BOp_s = twbp->p_bo2 * C12 / r2;
						Cln_BOp_pi = twbp->p_bo4 * C34 / r2;
						Cln_BOp_pi2 = twbp->p_bo6 * C56 / r2;

						/* Only dln_BOp_xx wrt. dr_i is stored here, note that
						 dln_BOp_xx/dr_i = -dln_BOp_xx/dr_j and all others are 0 */
						rvec_Scale(bo_ij->dln_BOp_s, -bo_ij->BO_s * Cln_BOp_s,
								ibond->dvec);
						rvec_Scale(bo_ij->dln_BOp_pi, -bo_ij->BO_pi
								* Cln_BOp_pi, ibond->dvec);
						rvec_Scale(bo_ij->dln_BOp_pi2, -bo_ij->BO_pi2
								* Cln_BOp_pi2, ibond->dvec);
						rvec_Scale(bo_ji->dln_BOp_s, -1., bo_ij->dln_BOp_s);
						rvec_Scale(bo_ji->dln_BOp_pi, -1., bo_ij->dln_BOp_pi);
						rvec_Scale(bo_ji->dln_BOp_pi2, -1., bo_ij->dln_BOp_pi2);

						/* Only dBOp wrt. dr_i is stored here, note that
						 dBOp/dr_i = -dBOp/dr_j and all others are 0 */
						rvec_Scale(bo_ij->dBOp, -(bo_ij->BO_s * Cln_BOp_s
								+ bo_ij->BO_pi * Cln_BOp_pi + bo_ij->BO_pi2
								* Cln_BOp_pi2), ibond->dvec);
						rvec_Scale(bo_ji->dBOp, -1., bo_ij->dBOp);

						rvec_Add(workspace->dDeltap_self[i], bo_ij->dBOp);
						rvec_Add(workspace->dDeltap_self[j], bo_ji->dBOp);

						bo_ij->BO_s -= control->bo_cut;
						bo_ij->BO -= control->bo_cut;
						bo_ji->BO_s -= control->bo_cut;
						bo_ji->BO -= control->bo_cut;
						workspace->total_bond_order[i] += bo_ij->BO; //currently total_BOp
						workspace->total_bond_order[j] += bo_ji->BO; //currently total_BOp
						bo_ij->Cdbo = bo_ij->Cdbopi = bo_ij->Cdbopi2 = 0.0;
						bo_ji->Cdbo = bo_ji->Cdbopi = bo_ji->Cdbopi2 = 0.0;

						/*fprintf( stderr, "%d %d %g %g %g\n",
						 i+1, j+1, bo_ij->BO, bo_ij->BO_pi, bo_ij->BO_pi2 );*/

						/*fprintf( stderr, "Cln_BOp_s: %f, pbo2: %f, C12:%f\n",
						 Cln_BOp_s, twbp->p_bo2, C12 );
						 fprintf( stderr, "Cln_BOp_pi: %f, pbo4: %f, C34:%f\n",
						 Cln_BOp_pi, twbp->p_bo4, C34 );
						 fprintf( stderr, "Cln_BOp_pi2: %f, pbo6: %f, C56:%f\n",
						 Cln_BOp_pi2, twbp->p_bo6, C56 );*/
						/*fprintf(stderr, "pbo1: %f, pbo2:%f\n", twbp->p_bo1, twbp->p_bo2);
						 fprintf(stderr, "pbo3: %f, pbo4:%f\n", twbp->p_bo3, twbp->p_bo4);
						 fprintf(stderr, "pbo5: %f, pbo6:%f\n", twbp->p_bo5, twbp->p_bo6);
						 fprintf( stderr, "r_s: %f, r_p: %f, r_pp: %f\n",
						 twbp->r_s, twbp->r_p, twbp->r_pp );
						 fprintf( stderr, "C12: %g, C34:%g, C56:%g\n", C12, C34, C56 );*/

						/*fprintf( stderr, "\tfactors: %g %g %g\n",
						 -(bo_ij->BO_s * Cln_BOp_s + bo_ij->BO_pi * Cln_BOp_pi +
						 bo_ij->BO_pi2 * Cln_BOp_pp),
						 -bo_ij->BO_pi * Cln_BOp_pi, -bo_ij->BO_pi2 * Cln_BOp_pi2 );*/
						/*fprintf( stderr, "dBOpi:\t[%g, %g, %g]\n",
						 bo_ij->dBOp[0], bo_ij->dBOp[1], bo_ij->dBOp[2] );
						 fprintf( stderr, "dBOpi:\t[%g, %g, %g]\n",
						 bo_ij->dln_BOp_pi[0], bo_ij->dln_BOp_pi[1],
						 bo_ij->dln_BOp_pi[2] );
						 fprintf( stderr, "dBOpi2:\t[%g, %g, %g]\n\n",
						 bo_ij->dln_BOp_pi2[0], bo_ij->dln_BOp_pi2[1],
						 bo_ij->dln_BOp_pi2[2] );*/

						Set_End_Index(j, btop_j + 1, bonds);
					}
				}
			}
		}

		H->entries[Htop].j = i;
		H->entries[Htop].val = system->reaxprm.sbp[type_i].eta;
		++Htop;

		Set_End_Index(i, btop_i, bonds);
		if (ihb == 1)
			Set_End_Index(workspace->hbond_index[i], ihb_top, hbonds);
		//fprintf( stderr, "%d bonds start: %d, end: %d\n",
		//     i, Start_Index( i, bonds ), End_Index( i, bonds ) );
	}

	// mark the end of j list
	H->start[i] = Htop;
	/* validate lists - decide if reallocation is required! */
	Validate_Lists(workspace, lists, data->step, system->N, H->m, Htop,
			num_bonds, num_hbonds);

#if defined(DEBUG_FOCUS)
	fprintf( stderr, "step%d: Htop = %d, num_bonds = %d, num_hbonds = %d\n",
			data->step, Htop, num_bonds, num_hbonds );

#endif
}

void Init_Forces_Tab(reax_system *system, control_params *control,
		simulation_data *data, static_storage *workspace, list **lists,
		output_controls *out_control) {
	int i, j, pj;
	int start_i, end_i;
	int type_i, type_j;
	int Htop, btop_i, btop_j, num_bonds, num_hbonds;
	int tmin, tmax, r;
	int ihb, jhb, ihb_top, jhb_top;
	int flag;
	real r_ij, r2, self_coef;
	real val, dif, base;
	real C12, C34, C56;
	real Cln_BOp_s, Cln_BOp_pi, Cln_BOp_pi2;
	real BO, BO_s, BO_pi, BO_pi2;
	real p_boc1, p_boc2;
	sparse_matrix *H;
	list *far_nbrs, *bonds, *hbonds;
	single_body_parameters *sbp_i, *sbp_j;
	two_body_parameters *twbp;
	far_neighbor_data *nbr_pj;
	LR_lookup_table *t;
	reax_atom *atom_i, *atom_j;
	bond_data *ibond, *jbond;
	bond_order_data *bo_ij, *bo_ji;

	far_nbrs = *lists + FAR_NBRS;
	bonds = *lists + BONDS;
	hbonds = *lists + HBONDS;

	H = workspace->H;
	Htop = 0;
	num_bonds = 0;
	num_hbonds = 0;
	btop_i = btop_j = 0;
	p_boc1 = system->reaxprm.gp.l[0];
	p_boc2 = system->reaxprm.gp.l[1];

	for (i = 0; i < system->N; ++i) {
		atom_i = &(system->atoms[i]);
		type_i = atom_i->type;
		start_i = Start_Index(i, far_nbrs);
		end_i = End_Index(i, far_nbrs);
		H->start[i] = Htop;
		btop_i = End_Index(i, bonds);
		sbp_i = &(system->reaxprm.sbp[type_i]);
		ihb = ihb_top = -1;
		if (control->hb_cut > 0 && (ihb = sbp_i->p_hbond) == 1)
			ihb_top = End_Index(workspace->hbond_index[i], hbonds);

		for (pj = start_i; pj < end_i; ++pj) {
			nbr_pj = &(far_nbrs->select.far_nbr_list[pj]);
			j = nbr_pj->nbr;
			atom_j = &(system->atoms[j]);

			flag = 0;
			if ((data->step - data->prev_steps) % control->reneighbor == 0) {
				if (nbr_pj->d <= control->r_cut)
					flag = 1;
				else
					flag = 0;
			} else if ((nbr_pj->d = Sq_Distance_on_T3(atom_i->x, atom_j->x,
					&(system->box), nbr_pj->dvec)) <= SQR(control->r_cut)) {
				nbr_pj->d = sqrt(nbr_pj->d);
				flag = 1;
			}

			if (flag) {
				type_j = system->atoms[j].type;
				r_ij = nbr_pj->d;
				sbp_j = &(system->reaxprm.sbp[type_j]);
				twbp = &(system->reaxprm.tbp[type_i][type_j]);
				self_coef = (i == j) ? 0.5 : 1.0;
				tmin = MIN( type_i, type_j );
				tmax = MAX( type_i, type_j );
				t = &(LR[tmin][tmax]);

				/* cubic spline interpolation */
				r = (int) (r_ij * t->inv_dx);
				if (r == 0)
					++r;
				base = (real) (r + 1) * t->dx;
				dif = r_ij - base;
				val = ((t->ele[r].d * dif + t->ele[r].c) * dif + t->ele[r].b)
						* dif + t->ele[r].a;
				val *= EV_to_KCALpMOL / C_ele;

				H->entries[Htop].j = j;
				H->entries[Htop].val = self_coef * val;
				++Htop;

				/* hydrogen bond lists */
				if (control->hb_cut > 0 && (ihb == 1 || ihb == 2) && nbr_pj->d
						<= control->hb_cut) {
					// fprintf( stderr, "%d %d\n", atom1, atom2 );
					jhb = sbp_j->p_hbond;
					if (ihb == 1 && jhb == 2) {
						hbonds->select.hbond_list[ihb_top].nbr = j;
						hbonds->select.hbond_list[ihb_top].scl = 1;
						hbonds->select.hbond_list[ihb_top].ptr = nbr_pj;
						++ihb_top;
						++num_hbonds;
					} else if (ihb == 2 && jhb == 1) {
						jhb_top = End_Index(workspace->hbond_index[j], hbonds);
						hbonds->select.hbond_list[jhb_top].nbr = i;
						hbonds->select.hbond_list[jhb_top].scl = -1;
						hbonds->select.hbond_list[jhb_top].ptr = nbr_pj;
						Set_End_Index(workspace->hbond_index[j], jhb_top + 1,
								hbonds);
						++num_hbonds;
					}
				}

				/* uncorrected bond orders */
				if (far_nbrs->select.far_nbr_list[pj].d <= control->nbr_cut) {
					r2 = SQR(r_ij);

					if (sbp_i->r_s > 0.0 && sbp_j->r_s > 0.0) {
						C12 = twbp->p_bo1 * POW(r_ij / twbp->r_s, twbp->p_bo2);
						BO_s = (1.0 + control->bo_cut) * EXP(C12);
					} else
						BO_s = C12 = 0.0;

					if (sbp_i->r_pi > 0.0 && sbp_j->r_pi > 0.0) {
						C34 = twbp->p_bo3 * POW(r_ij / twbp->r_p, twbp->p_bo4);
						BO_pi = EXP(C34);
					} else
						BO_pi = C34 = 0.0;

					if (sbp_i->r_pi_pi > 0.0 && sbp_j->r_pi_pi > 0.0) {
						C56 = twbp->p_bo5 * POW(r_ij / twbp->r_pp, twbp->p_bo6);
						BO_pi2 = EXP(C56);
					} else
						BO_pi2 = C56 = 0.0;

					/* Initially BO values are the uncorrected ones, page 1 */
					BO = BO_s + BO_pi + BO_pi2;

					if (BO >= control->bo_cut) {
						num_bonds += 2;
						/****** bonds i-j and j-i ******/
						ibond = &(bonds->select.bond_list[btop_i]);
						btop_j = End_Index(j, bonds);
						jbond = &(bonds->select.bond_list[btop_j]);

						ibond->nbr = j;
						jbond->nbr = i;
						ibond->d = r_ij;
						jbond->d = r_ij;
						rvec_Copy(ibond->dvec, nbr_pj->dvec);
						rvec_Scale(jbond->dvec, -1, nbr_pj->dvec);
						ivec_Copy(ibond->rel_box, nbr_pj->rel_box);
						ivec_Scale(jbond->rel_box, -1, nbr_pj->rel_box);
						ibond->dbond_index = btop_i;
						jbond->dbond_index = btop_i;
						ibond->sym_index = btop_j;
						jbond->sym_index = btop_i;
						++btop_i;
						Set_End_Index(j, btop_j + 1, bonds);

						bo_ij = &(ibond->bo_data);
						bo_ji = &(jbond->bo_data);
						bo_ji->BO = bo_ij->BO = BO;
						bo_ji->BO_s = bo_ij->BO_s = BO_s;
						bo_ji->BO_pi = bo_ij->BO_pi = BO_pi;
						bo_ji->BO_pi2 = bo_ij->BO_pi2 = BO_pi2;

						/* Bond Order page2-3, derivative of total bond order prime */
						Cln_BOp_s = twbp->p_bo2 * C12 / r2;
						Cln_BOp_pi = twbp->p_bo4 * C34 / r2;
						Cln_BOp_pi2 = twbp->p_bo6 * C56 / r2;

						/* Only dln_BOp_xx wrt. dr_i is stored here, note that
						 dln_BOp_xx/dr_i = -dln_BOp_xx/dr_j and all others are 0 */
						rvec_Scale(bo_ij->dln_BOp_s, -bo_ij->BO_s * Cln_BOp_s,
								ibond->dvec);
						rvec_Scale(bo_ij->dln_BOp_pi, -bo_ij->BO_pi
								* Cln_BOp_pi, ibond->dvec);
						rvec_Scale(bo_ij->dln_BOp_pi2, -bo_ij->BO_pi2
								* Cln_BOp_pi2, ibond->dvec);
						rvec_Scale(bo_ji->dln_BOp_s, -1., bo_ij->dln_BOp_s);
						rvec_Scale(bo_ji->dln_BOp_pi, -1., bo_ij->dln_BOp_pi);
						rvec_Scale(bo_ji->dln_BOp_pi2, -1., bo_ij->dln_BOp_pi2);

						/* Only dBOp wrt. dr_i is stored here, note that
						 dBOp/dr_i = -dBOp/dr_j and all others are 0 */
						rvec_Scale(bo_ij->dBOp, -(bo_ij->BO_s * Cln_BOp_s
								+ bo_ij->BO_pi * Cln_BOp_pi + bo_ij->BO_pi2
								* Cln_BOp_pi2), ibond->dvec);
						rvec_Scale(bo_ji->dBOp, -1., bo_ij->dBOp);

						rvec_Add(workspace->dDeltap_self[i], bo_ij->dBOp);
						rvec_Add(workspace->dDeltap_self[j], bo_ji->dBOp);

						bo_ij->BO_s -= control->bo_cut;
						bo_ij->BO -= control->bo_cut;
						bo_ji->BO_s -= control->bo_cut;
						bo_ji->BO -= control->bo_cut;
						workspace->total_bond_order[i] += bo_ij->BO; //currently total_BOp
						workspace->total_bond_order[j] += bo_ji->BO; //currently total_BOp
						bo_ij->Cdbo = bo_ij->Cdbopi = bo_ij->Cdbopi2 = 0.0;
						bo_ji->Cdbo = bo_ji->Cdbopi = bo_ji->Cdbopi2 = 0.0;

						Set_End_Index(j, btop_j + 1, bonds);
					}
				}
			}
		}

		H->entries[Htop].j = i;
		H->entries[Htop].val = system->reaxprm.sbp[type_i].eta;
		++Htop;

		Set_End_Index(i, btop_i, bonds);
		if (ihb == 1)
			Set_End_Index(workspace->hbond_index[i], ihb_top, hbonds);
	}

	// mark the end of j list
	H->start[i] = Htop;
	/* validate lists - decide if reallocation is required! */
	Validate_Lists(workspace, lists, data->step, system->N, H->m, Htop,
			num_bonds, num_hbonds);

#if defined(DEBUG_FOCUS)
	fprintf( stderr, "step%d: Htop = %d, num_bonds = %d, num_hbonds = %d\n",
			data->step, Htop, num_bonds, num_hbonds );
	//Print_Bonds( system, bonds, "sbonds.out" );
	//Print_Bond_List2( system, bonds, "sbonds.out" );
	//Print_Sparse_Matrix2( H, "H.out" );
#endif
}

void Estimate_Storage_Sizes(reax_system *system, control_params *control,
		list **lists, int *Htop, int *hb_top, int *bond_top, int *num_3body) {
	int i, j, pj;
	int start_i, end_i;
	int type_i, type_j;
	int ihb, jhb;
	real r_ij, r2;
	real C12, C34, C56;
	real BO, BO_s, BO_pi, BO_pi2;
	real p_boc1, p_boc2;
	list *far_nbrs;
	single_body_parameters *sbp_i, *sbp_j;
	two_body_parameters *twbp;
	far_neighbor_data *nbr_pj;
	reax_atom *atom_i, *atom_j;

	far_nbrs = *lists + FAR_NBRS;
	p_boc1 = system->reaxprm.gp.l[0];
	p_boc2 = system->reaxprm.gp.l[1];

	for (i = 0; i < system->N; ++i) {
		atom_i = &(system->atoms[i]);
		type_i = atom_i->type;
		start_i = Start_Index(i, far_nbrs);
		end_i = End_Index(i, far_nbrs);
		sbp_i = &(system->reaxprm.sbp[type_i]);
		ihb = sbp_i->p_hbond;

		for (pj = start_i; pj < end_i; ++pj) {
			nbr_pj = &(far_nbrs->select.far_nbr_list[pj]);
			j = nbr_pj->nbr;
			atom_j = &(system->atoms[j]);
			type_j = atom_j->type;
			sbp_j = &(system->reaxprm.sbp[type_j]);
			twbp = &(system->reaxprm.tbp[type_i][type_j]);

			if (nbr_pj->d <= control->r_cut) {
				++(*Htop);

				/* hydrogen bond lists */
				if (control->hb_cut > 0.1 && (ihb == 1 || ihb == 2)
						&& nbr_pj->d <= control->hb_cut) {
					jhb = sbp_j->p_hbond;
					if (ihb == 1 && jhb == 2)
						++hb_top[i];
					else if (ihb == 2 && jhb == 1)
						++hb_top[j];
				}

				/* uncorrected bond orders */
				if (nbr_pj->d <= control->nbr_cut) {
					r_ij = nbr_pj->d;
					r2 = SQR(r_ij);

					if (sbp_i->r_s > 0.0 && sbp_j->r_s > 0.0) {
						C12 = twbp->p_bo1 * POW(r_ij / twbp->r_s, twbp->p_bo2);
						BO_s = (1.0 + control->bo_cut) * EXP(C12);
					} else
						BO_s = C12 = 0.0;

					if (sbp_i->r_pi > 0.0 && sbp_j->r_pi > 0.0) {
						C34 = twbp->p_bo3 * POW(r_ij / twbp->r_p, twbp->p_bo4);
						BO_pi = EXP(C34);
					} else
						BO_pi = C34 = 0.0;

					if (sbp_i->r_pi_pi > 0.0 && sbp_j->r_pi_pi > 0.0) {
						C56 = twbp->p_bo5 * POW(r_ij / twbp->r_pp, twbp->p_bo6);
						BO_pi2 = EXP(C56);
					} else
						BO_pi2 = C56 = 0.0;

					/* Initially BO values are the uncorrected ones, page 1 */
					BO = BO_s + BO_pi + BO_pi2;

					if (BO >= control->bo_cut) {
						++bond_top[i];
						++bond_top[j];
					}
				}
			}
		}
	}

	*Htop += system->N;
	*Htop *= SAFE_ZONE;
	for (i = 0; i < system->N; ++i) {
		hb_top[i] = MAX( hb_top[i] * SAFE_HBONDS, MIN_HBONDS );
		*num_3body += SQR(bond_top[i]);
		bond_top[i] = MAX( bond_top[i] * 2, MIN_BONDS );
	}
	*num_3body *= SAFE_ZONE;
}

void Compute_Forces(reax_system *system, control_params *control,
		simulation_data *data, static_storage *workspace, list** lists,
		output_controls *out_control) {
	real t_start, t_elapsed;

	t_start = Get_Time();
	if (!control->tabulate)
		Init_Forces(system, control, data, workspace, lists, out_control);
	else
		Init_Forces_Tab(system, control, data, workspace, lists, out_control);
	t_elapsed = Get_Timing_Info(t_start);
	data->timing.init_forces += t_elapsed;
#if defined(DEBUG_FOCUS)
	fprintf( stderr, "init_forces - ");
#endif

	t_start = Get_Time();
	Compute_Bonded_Forces(system, control, data, workspace, lists, out_control);
	t_elapsed = Get_Timing_Info(t_start);
	data->timing.bonded += t_elapsed;
#if defined(DEBUG_FOCUS)  
	fprintf( stderr, "bonded_forces - ");
#endif

	t_start = Get_Time();
	Compute_NonBonded_Forces(system, control, data, workspace, lists,
			out_control);
	t_elapsed = Get_Timing_Info(t_start);
	data->timing.nonb += t_elapsed;
#if defined(DEBUG_FOCUS)
	fprintf( stderr, "nonbondeds - ");
#endif

	Compute_Total_Force(system, control, data, workspace, lists);
	// implement amd simulation
	if (control->amd)
	{
		Compute_AMD_Force(system, control, data, workspace, lists);
	}
    if (control->bboost == 1)
        Compute_Bond_Boost_Force(system, control, data, workspace, lists, out_control);
    else if (control->bboost == 2)
        Compute_Bond_Boost_Force_All(system, control, data, workspace, lists, out_control);
    else if (control->bboost == 3)
        Compute_Bond_Boost_Force_All_Couple(system, control, data, workspace, lists, out_control);
	//Print_Total_Force( system, control, data, workspace, lists, out_control );
        //Bias Potential
        //Bias_Spring(system, control, data, workspace, lists, out_control);
        if (control->bias_con_de)
            Bias_COn_Decompose(system, control, data, workspace, lists, out_control);
        else if (control->bias_con_com)
            Bias_COn_Combine(system, control, data, workspace, lists, out_control);
   if (control->bias_lj126 == 1)
       Bias_LJ_126(system, control, data, workspace, lists, out_control);
   if (control->bias_charge == 1)
       Bias_Charge(system, control, data, workspace, lists, out_control);
#if defined(DEBUG_FOCUS)
	fprintf( stderr, "totalforces - ");
	//Print_Total_Force( system, control, data, workspace, lists, out_control );
#endif

#ifdef TEST_FORCES
	Print_Total_Force( system, control, data, workspace, lists, out_control );
	Compare_Total_Forces( system, control, data, workspace, lists, out_control );
#endif
#if defined(DEBUG_FOCUS)  
	fprintf( stderr, "forces - ");
#endif
}
