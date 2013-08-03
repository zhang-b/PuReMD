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

void Bias_Foo()
{
    printf("I am in bias.cpp !\n");
    return;
}

int Reactive_O_Atom(int grp[], int n, int atom)
{
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
    int i, j;
    int start_i, end_i; 
    int pj;
    int step, interval;
    int n, n_o, n_c;
    int MAX = 400;
    int grp_o[MAX], grp_c[MAX];
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
    interval = 5000;

    if ( step%interval == 0){
        data->bias_counter = 0;
        data->bias_success = 0;
        n_o = 0;
        n_c = 0;
        r_min = 40.0;

        for (i=0; i<MAX; i++){
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
                    //printf("atom %s with bo %.2f in pj %d\n", j_elem, bo, pj);
                    if ( bo>0.80 && strcmp(j_elem, "O") == 0){
                        grp_o[n_o] = j;
                        n_o++;
                        n++;
                    }
                }
                if (n < 3){
                    grp_c[n_c] = i;
                    n_c++;
                }
            }
        }
    
        //printf("totally %d C in CO2 group, and %d O in CO2 group\n", n_c, n_o);
    
        for ( i = 0; i < n_c; i++){
            for ( j = 0; j < system->N; j++) {
                j_elem = system->reaxprm.sbp[system->atoms[j].type].name;
                flag = Reactive_O_Atom( grp_o, n_o, j);
                if (strcmp(j_elem, "O") == 0 && flag){
                    atom1 = &( system->atoms[grp_c[i]] );
                    atom2 = &( system->atoms[j] );
                    rvec_ScaledSum(rv, 1, atom1->x, -1, atom2->x); 
                    r_ij = rvec_Norm(rv);
                    if ( r_min > r_ij && r_ij < 3.0 ){
                        r_min = r_ij;
                        data->bias_r = r_min;
                        data->bias_atom1 = grp_c[i];
                        data->bias_atom2 = j;
                        data->bias_success = 1;
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
        if ( r_ij > 5.0 )
            scale = 0.0;
        else if ( r_ij > 3.0 && r_ij <= 5.0 )
            scale = control->bias_V * 2 * 1.65;
        else if ( r_ij <= 3.0 && r_ij > 1.35)
            scale = control->bias_V * 2 * ( r_ij - 1.35);
        else {
            scale = 0;
            data->bias_counter ++;
            if (data->bias_counter > 500)
                data->bias_success = 0;
        }
        //scale = 0.0;
        atom1->f[0] += scale * rv[0];
        atom1->f[1] += scale * rv[1];
        atom1->f[2] += scale * rv[2];
        //atom2->f[0] += 1 * scale * rv[0];
        //atom2->f[1] += 1 * scale * rv[1];
        //atom2->f[2] += 1 * scale * rv[2];

        /*
        if (abs(scale) > 0)
            printf("step %d n = %d atom1 = %d atom2 = %d r_ij = %.2f scale = %.2f\n", 
                   data->step, data->bias_counter, data->bias_atom1, 
                   data->bias_atom2, r_ij, scale);
        printf("v1x = %.2f v1y = %.2f v1z = %.2f ", atom1->v[0], atom1->v[1], atom1->v[2]);
        printf("v2x = %.2f v2y = %.2f v2z = %.2f\n", atom1->v[0], atom1->v[1], atom1->v[2]);
        */
    }
    return;
}

void Bias_COn_Decompose(reax_system *system, control_params *control,
                simulation_data *data, static_storage *workspace, list **lists,
        output_controls *out_control)
{
    /* 
     * 
     * */
    int i, j, adatom1, adatom2, adatom_mol; 
    int pj;
    int start_i, end_i; 
    int n, step, interval;
    int flag;
    char *i_elem, *j_elem;
    real bo, bo_min, bo_min_mol;
    real r, scale;
    rvec rv, df;
    reax_atom *atom1, *atom2;
    list *bonds;
    bond_order_data *bo_ij;
    bonds = (*lists) + BONDS;

    step = data->step - data->prev_steps;
    interval = 5000;

    if ( step%interval == 0){
        data->bias_success = 0;
        bo_min = 2.0;
        flag = 0;
        for ( i=0; i<system->N; ++i){
            n = 0;
            bo_min_mol = 2.0;
            i_elem = system->reaxprm.sbp[system->atoms[i].type].name;
            if (strcmp(i_elem, "C") == 0){
                start_i = Start_Index(i, bonds);
                end_i = End_Index(i, bonds);
                for( pj = start_i; pj < end_i; ++pj ){
                    j = bonds->select.bond_list[pj].nbr;
                    j_elem = system->reaxprm.sbp[system->atoms[j].type].name;
                    bo_ij = &( bonds->select.bond_list[pj].bo_data );
                    bo = bo_ij->BO;
                    if ( bo>0.30 && strcmp(j_elem, "O") == 0){
                        n += 1;
                        if (bo_min_mol > bo) {
                            bo_min_mol = bo;
                            adatom_mol = j;
                        }
                    }
                }
            }
            if ( n > control->bias_con ){
                adatom1 = i;
                adatom2 = adatom_mol;
                flag = 1;
                data->bias_success = 1;
                data->bias_n = n;
            }
        }

        if (flag){
            data->bias_atom1 = adatom1;
            data->bias_atom2 = adatom2;
            /*
            printf("%s ",system->reaxprm.sbp[system->atoms[adatom1].type].name);
            printf("%s\n",system->reaxprm.sbp[system->atoms[adatom2].type].name);
            */
        }
    }

    if (data->bias_success && data->bias_n > control->bias_con) {
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
            if (strcmp(j_elem, "O") == 0 && bo > 0.3){
                n++;
                if ( j == data->bias_atom2) {
                    if ( r < 1.9) 
                        scale = -control->bias_V * 2 * (r - 1.9) / r;
                    else
                        scale = 0;
                    atom2->f[0] += scale * rv[0];
                    atom2->f[1] += scale * rv[1];
                    atom2->f[2] += scale * rv[2];
                    atom1->f[0] += -1 * scale * rv[0];
                    atom1->f[1] += -1 * scale * rv[1];
                    atom1->f[2] += -1 * scale * rv[2];
                    //printf("step = %d scale = %.2f atomj = %d bo = %.2f r = %.2f adatom %d\n",data->step, scale, j, bo, r, data->bias_atom2 );
                    //printf("vx = %.2f vy = %.2f vz = %.2f\n",atom2->v[0], atom2->v[1], atom2->v[2]);
                }
                /*
                else {
                    scale = -control->bias_V * 2 * (r - 1.35) / r;
                    atom2->f[0] += scale * rv[0];
                    atom2->f[0] += scale * rv[1];
                    atom2->f[0] += scale * rv[2];
                }
                */
            }
        }
        data->bias_n = n;
    }
    return;
}

void Bias_Spring(reax_system *system, control_params *control,
                simulation_data *data, static_storage *workspace, list **lists,
        output_controls *out_control)
{
    int i, j;
    int pj;
    int n1, n2, n3;
    int start_i, end_i; 
    char *i_elem, *j_elem;
    real bo;
    list *bonds;
    bond_order_data *bo_ij;
    bonds = (*lists) + BONDS;
    
    for ( i=0; i<system->N; ++i){
        i_elem = system->reaxprm.sbp[system->atoms[i].type].name;
        n1 = 0;
        n2 = 0;
        n3 = 0;
        if (strcmp(i_elem, "O") == 0){
            start_i = Start_Index(i, bonds);
            end_i = End_Index(i, bonds);
            for( pj = start_i; pj < end_i; ++pj ){
                 j = bonds->select.bond_list[pj].nbr;
                 j_elem = system->reaxprm.sbp[system->atoms[j].type].name;
                 bo_ij = &( bonds->select.bond_list[pj].bo_data );
                 bo = bo_ij->BO;
                 /*
                 if (strcmp(j_elem, "CA") == 0){
                     if ( bo>=0.0 && bo<0.3)
                         nf++;
                     else if ( bo>=0.3 && bo<2)
                         nn++;
                 }
                 */
                 if ( bo>0.60){
                     if (strcmp(j_elem, "CA") == 0)
                         n1 += 1;
                     else if (strcmp(j_elem, "C") == 0)
                         n2 += 1;
                     else if (strcmp(j_elem, "O") == 0)
                         n3 += 1;
                 }
            }
        //printf("atom %5d with Ca = %d, C = %d, O = %d\n",i , n1, n2, n3);
        }
    }

    return;
}

