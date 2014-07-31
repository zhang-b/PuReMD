/*----------------------------------------------------------------------
 * A set of bias potentials to accelate the simulation
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *----------------------------------------------------------------------*/
#include "bias.h"
#include "list.h"
#include "vector.h"

#define BIAS_V_MAX 40.0          /* the maximum velocity allowed */  
#define BIAS_ATOM_MAX 4000        /* maximum atom number */
#define BIAS_R_MIN 40.0          /* abitrary initial value for r_min */
#define BIAS_R1_CO_MAX  3.0      /* cut-off for non-bond distance */
#define BIAS_R2_CO_MAX  5.0      /* cut-off2 for non-bond distance */
#define BIAS_RE_CO_MAX  1.40     /* re for C-O distance */    
#define BIAS_LAST_CO 500         /* number of stablized steps */
#define BIAS_BOND_CUTOFF   0.6   /* cut-off for bond */
#define BIAS_BOND_CUTOFF2   0.8   /* cut-off for bond (more rigorous)*/

void Bias_Foo()
{
    /* A test function of the bias algorithm */
    printf("I am in bias.cpp !\n");
    return;
}

int Slow_Down_Atom(reax_atom *atom)
{
    /* Slow down the atoms if the atoms is over accelerated */
    real scale;
    real v;
    v = rvec_Norm(atom->v);
    if ( v > BIAS_V_MAX ){
        scale = BIAS_V_MAX/v;
        rvec_Scale(atom->v, scale, atom->v);
        return 1;
    }
    return 0;
}

int Reactive_O_Atom(int grp[], int n, int atom)
{
    /* Determine if an O atom is reactive */
    int i;
    int flag;
    flag = 0;
    for ( i = 0; i < n; i++) {
        if (atom == grp[i])
            flag++;
    }
    if (flag > 0)
    	return 0;
    else
        return 1;
}

void Bias_COn_Combine(reax_system *system, control_params *control,
                simulation_data *data, static_storage *workspace, list **lists,
        output_controls *out_control)
{
    /* Combine CO2 and O to generate a CO3 */

    int i, j;
    int start_i, end_i; 
    int pj;
    int step, interval;
    int n, n_o, n_c;
    int grp_o[BIAS_ATOM_MAX], grp_c[BIAS_ATOM_MAX];
    int flag;
    real bo;
    real r_ij, r_min, scale;
    rvec rv;
    char *i_elem, *j_elem;
    reax_atom *atom1, *atom2;
    list *bonds;
    bond_order_data *bo_ij;
    bonds = (*lists) + BONDS;

    step = data->step - data->prev_steps;
    interval = control->bias_con_com_interval;

    if( step%interval == 0 ){
        fprintf(out_control->bias, "Step %d COn Combination\n", step);
        data->bias_search = 1;
    }

    if ( data->bias_search ){
        data->bias_counter = 0;
        data->bias_success = 0;
        n_o = 0;
        n_c = 0;
        r_min = BIAS_R_MIN;

        /* First, build the group index */
        for (i = 0; i < BIAS_ATOM_MAX; i++){
            grp_o[i] = -1;
            grp_c[i] = -1;
        }

        for ( i=0; i < system->N; i++ ){
            n = 0;
            i_elem = system->reaxprm.sbp[system->atoms[i].type].name;
            if (strcmp(i_elem, "C") == 0){
                start_i = Start_Index(i, bonds);
                end_i = End_Index(i, bonds);
                for( pj = start_i; pj < end_i; ++pj ) {
                    j = bonds->select.bond_list[pj].nbr;
                    j_elem = system->reaxprm.sbp[system->atoms[j].type].name;
                    bo_ij = &( bonds->select.bond_list[pj].bo_data );
                    bo = bo_ij->BO;  
                    if ( bo > BIAS_BOND_CUTOFF2 && strcmp(j_elem, "O") == 0){
                        grp_o[n_o] = j;
                        n_o++;
                        n++;
                    }
                }
                if (n < control->bias_con_com_n){
                    grp_c[n_c] = i;
                    n_c++;
                }
            }
        }
    
        fprintf(out_control->bias, "Totally %d C in CO2 group, and %d O in CO2 group\n", n_c, n_o);
    
        /* Find the atoms (atom1 and atom2) to apply bias potential */
        for ( i = 0; i < n_c; i++){
            for ( j = 0; j < system->N; j++) {
                j_elem = system->reaxprm.sbp[system->atoms[j].type].name;
                flag = Reactive_O_Atom( grp_o, n_o, j);
                if (strcmp(j_elem, "O") == 0 && flag){
                    atom1 = &( system->atoms[grp_c[i]] );
                    atom2 = &( system->atoms[j] );
                    rvec_ScaledSum(rv, 1, atom1->x, -1, atom2->x); 
                    r_ij = rvec_Norm(rv);
                    if ( r_min > r_ij && r_ij < BIAS_R1_CO_MAX ){
                        r_min = r_ij;
                        data->bias_r = r_min;
                        data->bias_atom1 = grp_c[i];
                        data->bias_atom2 = j;
                        data->bias_success = 1;
                        data->bias_search = 0;
                    }
                } 
            }
        }
    }
    
    if (data->bias_success){
        atom1 = &( system->atoms[data->bias_atom1] );
        atom2 = &( system->atoms[data->bias_atom2] );
        rvec_ScaledSum(rv, 1, atom1->x, -1, atom2->x); 
        r_ij = rvec_Norm(rv);
        /* bias potential to drag the C (atom1) and O(atom2) together */
        if ( r_ij > BIAS_R2_CO_MAX )
            scale = 0.0;
        else if ( r_ij > BIAS_R1_CO_MAX && r_ij <= BIAS_R2_CO_MAX )
            scale = control->bias_con_com_vmax * 2 * 1.65;
        else if ( r_ij <= BIAS_R1_CO_MAX && r_ij > BIAS_RE_CO_MAX)
            scale = control->bias_con_com_vmax * 2 * ( r_ij - BIAS_RE_CO_MAX);
        else {
            scale = 0;
            data->bias_counter ++;
            if (data->bias_counter > BIAS_LAST_CO)
                data->bias_success = 0;
        }
        atom1->f[0] += scale * rv[0];
        atom1->f[1] += scale * rv[1];
        atom1->f[2] += scale * rv[2];
        Slow_Down_Atom(atom1);
        Slow_Down_Atom(atom2);

        fprintf(out_control->bias, "atom1 = %d atom2 = %d r_ij = %.2f scale = %.2f sucess = %d\n", 
                data->bias_atom1, data->bias_atom2, r_ij, scale, data->bias_counter);
        /*
        printf("v1x = %.2f v1y = %.2f v1z = %.2f ", atom1->v[0], atom1->v[1], atom1->v[2]);
        printf("v2x = %.2f v2y = %.2f v2z = %.2f\n", atom1->v[0], atom1->v[1], atom1->v[2]);
        */
    }
    fflush(out_control->bias);
    return;
}

void Bias_COn_Decompose(reax_system *system, control_params *control,
                simulation_data *data, static_storage *workspace, list **lists,
        output_controls *out_control)
{
    /* Apply a force to drag one C-O bond from CO3 to get a CO2 */
    int i, j, adatom1, adatom2, adatom_mol; 
    int pj;
    int start_i, end_i; 
    int n, step, interval;
    int flag;
    char *i_elem, *j_elem;
    char *atp1, *atp2;
    real bo, bo_min, bo_min_mol;
    real r, scale;
    rvec rv;
    reax_atom *atom1, *atom2;
    list *bonds;
    bond_order_data *bo_ij;
    bonds = (*lists) + BONDS;

    atp1 = system->reaxprm.sbp[control->bias_con_de_atom1 - 1].name;
    atp2 = system->reaxprm.sbp[control->bias_con_de_atom2 - 1].name;

    step = data->step - data->prev_steps;
    interval = control->bias_con_de_interval;
    
    if( step == 0 || data->bias_success)
        fprintf(out_control->bias, "Step %d COn decompostion\n", step);

    /* determin the drag atoms */
    if ( step%interval == 0){
        data->bias_counter = 0;
        data->bias_success = 0;
        bo_min = 2.0;
        flag = 0;
        for ( i=0; i<system->N; ++i){
            n = 0;
            bo_min_mol = 2.0;
            i_elem = system->reaxprm.sbp[system->atoms[i].type].name;
            if (strcmp(i_elem, atp1) == 0){
                start_i = Start_Index(i, bonds);
                end_i = End_Index(i, bonds);
                for( pj = start_i; pj < end_i; ++pj ){
                    j = bonds->select.bond_list[pj].nbr;
                    j_elem = system->reaxprm.sbp[system->atoms[j].type].name;
                    bo_ij = &( bonds->select.bond_list[pj].bo_data );
                    bo = bo_ij->BO;
                    if ( bo > control->bias_con_de_cutoff&& strcmp(j_elem, atp2) == 0){
                        n += 1;
                        if (bo_min_mol > bo) {
                            bo_min_mol = bo;
                            adatom_mol = j;
                        }
                    }
                }
            }
            if ( n > control->bias_con_de_n){
                adatom1 = i;
                adatom2 = adatom_mol;
                flag = 1;
                data->bias_success = 1;
                data->bias_n = n;
            }
        }

        /* Update the drag atoms C (adatom1) and O (adatom2) */
        if (flag){
            data->bias_atom1 = adatom1;
            data->bias_atom2 = adatom2;
            fprintf(out_control->bias, "Update atoms C = %d and atom O = %d\n", 
                    adatom1, adatom2);
        }
    }

    /* apply bias potential */
    if (data->bias_success) {
        n = 0;
        start_i = Start_Index(data->bias_atom1, bonds);
        end_i = End_Index(data->bias_atom1, bonds);
        for ( pj = start_i; pj < end_i; pj++ ){
            j = bonds->select.bond_list[pj].nbr;
            j_elem = system->reaxprm.sbp[system->atoms[j].type].name;
            atom1 = &( system->atoms[data->bias_atom1] );
            atom2 = &( system->atoms[j] );
            rvec_ScaledSum(rv, 1, atom1->x, -1, atom2->x); 
            r = rvec_Norm(rv);
            bo_ij = &( bonds->select.bond_list[pj].bo_data );
            bo = bo_ij->BO;
            scale = 0.0;
            if (strcmp(j_elem, atp2) == 0 && bo > control->bias_con_de_cutoff){
                n++;
                /* bias potential */
                if ( j == data->bias_atom2) {
                    if ( r < 1.8) 
                        scale = -control->bias_con_de_vmax * 2 * (r - 1.8);
                    else {
                        scale = 0;
                        data->bias_counter++;
                    }
                    atom2->f[0] += scale * rv[0];
                    atom2->f[1] += scale * rv[1];
                    atom2->f[2] += scale * rv[2];
                    atom1->f[0] += -1 * scale * rv[0];
                    atom1->f[1] += -1 * scale * rv[1];
                    atom1->f[2] += -1 * scale * rv[2];
                    Slow_Down_Atom(atom1);
                    Slow_Down_Atom(atom2);
                    fprintf(out_control->bias, "scale = %.2f bo = %.2f r = %.2f sucess %d\n",
                    scale, bo, r, data->bias_counter);
                }
            }
        }
        data->bias_n = n;
        if (data->bias_n <= control->bias_con_de_n)
            data->bias_success = 0; 
        if (data->bias_counter > BIAS_LAST_CO)
            data->bias_success = 0; 
        if (data->step % interval > interval * 0.9)
            data->bias_success = 0; 
    }
    fflush(out_control->bias);

    return;
}

void Bias_LJ_126(reax_system *system, control_params *control,
                simulation_data *data, static_storage *workspace, list **lists,
        output_controls *out_control)
{
    /* Combine CO2 and O to generate a CO3 */
    int i, j;
    int start_i, end_i; 
    int pj;
    real r_ij, r3, r6, r12;
    real e_lj126, f_lj126; 
    real epslon, sigma;
    far_neighbor_data *nbr_pj;
    list *far_nbrs;

    sigma = control->bias_lj126_sigma;
    epslon = control->bias_lj126_epsilon;

    far_nbrs = (*lists) + FAR_NBRS;
    for( i = 0; i < system->N; ++i ) {
        start_i = Start_Index(i, far_nbrs);
        end_i   = End_Index(i, far_nbrs);
        for( pj = start_i; pj < end_i; ++pj )
            if( far_nbrs->select.far_nbr_list[pj].d <= control->r_cut ) {
                nbr_pj = &( far_nbrs->select.far_nbr_list[pj] );
                j = nbr_pj->nbr;
                r_ij = nbr_pj->d / sigma;
                r3 = r_ij * r_ij * r_ij;
                r6 = r3 * r3;
                r12 = r6 * r6;
                e_lj126 = 4 * epslon * ( 1/r12 - 1/r6);
                f_lj126 = -4 * epslon * (12/r12 - 6/r6) * (1/r_ij);
                if (abs(f_lj126) > 20.0)
                    f_lj126 = f_lj126/abs(f_lj126) * 20.0;
                rvec_ScaledAdd( system->atoms[i].f, -f_lj126, nbr_pj->dvec );
                rvec_ScaledAdd( system->atoms[j].f,  f_lj126, nbr_pj->dvec );
            }
    }

    return;
}

void Bias_Charge(reax_system *system, control_params *control,
                simulation_data *data, static_storage *workspace, list **lists,
        output_controls *out_control)
{
    /* Combine CO2 and O to generate a CO3 */
    int i, j;
    int start_i, end_i; 
    int pj;
    real r_ij, r2; 
    real e_charge, f_charge; 
    real dfactor; // damping factor
    far_neighbor_data *nbr_pj;
    list *far_nbrs;

    dfactor = control->bias_charge_dfactor;

    far_nbrs = (*lists) + FAR_NBRS;
    for( i = 0; i < system->N; ++i ) {
        start_i = Start_Index(i, far_nbrs);
        end_i   = End_Index(i, far_nbrs);
        for( pj = start_i; pj < end_i; ++pj )
            if( far_nbrs->select.far_nbr_list[pj].d <= control->r_cut ) {
                nbr_pj = &( far_nbrs->select.far_nbr_list[pj] );
                j = nbr_pj->nbr;
                r_ij = nbr_pj->d;
                r2 = r_ij * r_ij;
                e_charge = dfactor / r_ij;
                f_charge = - dfactor / r2;
                if (abs(f_charge) > 20.0)
                    f_charge = f_charge/abs(f_charge) * 20.0;
                rvec_ScaledAdd( system->atoms[i].f, -f_charge, nbr_pj->dvec );
                rvec_ScaledAdd( system->atoms[j].f,  f_charge, nbr_pj->dvec );
            }
    }

    return;
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
        if ( (fabs (emax) < fabs(e)) & (bo > q)) {
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
    printf("max allowd %d radicals in simulation.\n");

    if (fabs(emax) < q && nbond > 0 && nrad <= control->bboost_nrad) {
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
        if (bo > q && vmax > 0) {
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
  if (fabs(emax) < q && nbond > 0 && nrad <= control->bboost_nrad) {
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
          if (bo > q && vmax > 0) {

            e = (r - re)/re; // eta
            C1 = 1 - e*e/(q*q);
            V += vmax / nbond * C1;
            bf = 2*A*vmax*e/(nbond*q*q*re);

            if (fabs(e - emax) < 0.00001) {
                C2 = 1 - e*e*P1*P1/(q*q);
                f1 = 2/(q*q*re);
                f2 = f1/2*P1*P1;
                dA = 2*e*C1/C2*(f1 - f2*C1/C2);
            }
            rvec_Scale(df, bf, rv);
            atom1 = &( system->atoms[i] );
            atom2 = &( system->atoms[j] );
            rvec_Add(system->atoms[i].f, df);
            rvec_ScaledAdd(system->atoms[j].f, -1, df);
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
    rvec_ScaledAdd(system->atoms[adatom2].f, -1, df);
  }
  else {
    data->boost = 0;
  }
  bfactor = exp(4184 * A * V/(T * 8.314));
  //bfactor = V;
  fprintf( out_control->bboost, "%-10d%6d%6d%6d%3d%8.4f%8.4f %44.4f", \
  data->step, nbond, adatom, adatom2, nrad, r_max, emax, bfactor );

  fprintf( out_control->bboost, " %4s %4s\n", \
  system->reaxprm.sbp[ system->atoms[adatom].type ].name, 
  system->reaxprm.sbp[ system->atoms[adatom2].type ].name);
  fflush( out_control->bboost);
}

