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
#define BIAS_BOND_CUTOFF   0.3   /* cut-off for bond */
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
                    if ( bo > BIAS_BOND_CUTOFF && strcmp(j_elem, atp2) == 0){
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
            if (strcmp(j_elem, "O") == 0 && bo > BIAS_BOND_CUTOFF){
                n++;
                /* bias potential */
                if ( j == data->bias_atom2) {
                    if ( r < 1.9) 
                        scale = -control->bias_con_de_vmax * 2 * (r - 2.0);
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
