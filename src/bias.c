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

void Bias_COn(reax_system *system, control_params *control,
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
    interval = 500;
    if ( step%interval == 0){

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
            }
        }

        if (flag){
            data->bias_atom1 = adatom1;
            data->bias_atom2 = adatom2;
            /*
            printf("atom1 = %5d atom2 = %d, force = %.2f, rmax = %.2f ",adatom1 ,adatom2 , control->bias_V, r_max);
            printf("%s ",system->reaxprm.sbp[system->atoms[adatom1].type].name);
            printf("%s\n",system->reaxprm.sbp[system->atoms[adatom2].type].name);
            */
        }
    }
    atom1 = &( system->atoms[data->bias_atom1] );
    atom2 = &( system->atoms[data->bias_atom2] );
    rvec_ScaledSum(rv, 1, atom1->x, -1, atom2->x); 
    r = rvec_Norm(rv);
    if ( r < 2.5)
        scale = control->bias_V * (2.5 - r);
    else
        scale = 0;
    rvec_Scale(df, scale, rv);
    rvec_Add(atom2->f, df);
    printf("rx = %.2f ry = %.2f rz = %.2f r = %.2f\n", rv[0], rv[1], rv[2], r);
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
                 if ( bo>0.10){
                     if (strcmp(j_elem, "CA") == 0)
                         n1 += 1;
                     else if (strcmp(j_elem, "C") == 0)
                         n2 += 1;
                     else if (strcmp(j_elem, "O") == 0)
                         n3 += 1;
                 }
            }
        printf("atom %5d with Ca = %d, C = %d, O = %d\n",i , n1, n2, n3);
        }
    }

    return;
}

