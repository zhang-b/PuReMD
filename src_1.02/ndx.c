/*
 * =====================================================================================
 *
 *       Filename:  ndx.c
 *
 *    Description:  index for groups
 *
 *        Version:  1.0
 *        Created:  07/29/2014 05:20:44 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#include "ndx.h"

#define MAX_LINE 1024
#define MAX_GROUPS 10000
#define MAX_ATOMS   10000
#define IN 1
#define OUT 0

char Read_Ndx_File(FILE *fp, reax_groups *grps, control_params *control, 
                  output_controls *out_control) {
    char *s;
    char *ret1, *ret2;
    char *sep = "\t \n!=";
    char *tmp, *token;
    int n, natom, counter;
    int flag;

    natom = 0;
    counter = 0;
    flag = 0;
    
    s = (char*) calloc(MAX_LINE, sizeof(char));
    tmp = (char*) calloc(MAX_LINE, sizeof(char));
    token = (char*) calloc(MAX_LINE, sizeof(char));

    if ((fp = fopen("index.ndx", "r")) == NULL) {
        fprintf(stderr, "Error opening index file!\n");
        exit(FILE_NOT_FOUND_ERR);
    }
    
    while(!feof(fp)) {
        fgets(s, MAX_LINE, fp);
        ret1 = strchr(s, '[');
        ret2 = strchr(s, ']');
        if ((ret1 == NULL) || ret2 == NULL) {
            token = strtok(s, sep);
            while(token != NULL) {
                grps->atoms[grps->ngrps][counter] = atoi(token);
                token = strtok(NULL, sep);
                counter++;
            }
        }
        else {
            n = ret2 - ret1;
            strncpy(grps->names[grps->ngrps], ret1 + 2, n - 2);
            if (flag) {
                grps->natoms[grps->ngrps] = counter;
                grps->ngrps++;
                counter = 0;
            }
            flag = IN;
        }
    }
    grps->natoms[grps->ngrps] = counter - 1;
    //nlist[grps->ngrps] = counter;
    
    // Output
    
    /*
    for (i =0 ; i <= ng; i++) {
        printf("ng %d, %d\n", i, grps->natoms[i]);
        for (int j=0; j < grps->natoms[i]; j++)
            printf("%d", grps->atoms[i][j]);
        printf("\n");
    }
    */

    fclose(fp);

    free(s);
/*  
    free(grps->natoms);
    for (i = 0; i < MAX_GROUPS; i++)
        free(grps->names[i]);
        free(grps->atoms[i]);
    free(grps->names);
    free(grps->atoms);
    free(grps);
    return 0;
*/
    return 0;

}

char Make_Default_Groups(reax_groups *grps, control_params *control, 
                  reax_system* system, output_controls *out_control) {
    int i;

    grps = (reax_groups*) calloc(1, sizeof(reax_groups)); 
    grps->natoms = (int*) calloc(MAX_GROUPS, sizeof(int));
    grps->names = (char**) calloc(MAX_GROUPS, sizeof(char *));
    grps->atoms = (int**) calloc(MAX_GROUPS, sizeof(int*));
    grps->ngrps = 0;

    for (i = 0; i < MAX_GROUPS; i++) {
        grps->names[i] = (char*) calloc(MAX_LINE, sizeof(char));
        grps->atoms[i] = (int*) calloc(MAX_ATOMS, sizeof(int*));
    }
    grps->ngrps = 1;
    strcpy(grps->names[grps->ngrps], "system");
    for (i = 0; i<system->N; i++) {
        grps->atoms[grps->ngrps][i] = i;
    }
    grps->ngrps++;

    return 0;
}
