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
    int n, i, ng, natom, counter;
    int flag;

    ng = 0;
    natom = 0;
    counter = 0;
    flag = 0;
    
    s = (char*) calloc(MAX_LINE, sizeof(char));
    tmp = (char*) calloc(MAX_LINE, sizeof(char));
    token = (char*) calloc(MAX_LINE, sizeof(char));
    grps = (reax_groups*) calloc(1, sizeof(reax_groups)); 
    grps->natoms = (int*) calloc(MAX_GROUPS, sizeof(int));
    grps->names = (char**) calloc(MAX_GROUPS, sizeof(char *));
    grps->atoms = (int**) calloc(MAX_GROUPS, sizeof(int*));
    for (i = 0; i < MAX_GROUPS; i++) {
        grps->names[i] = (char*) calloc(MAX_LINE, sizeof(char));
        grps->atoms[i] = (int*) calloc(MAX_ATOMS, sizeof(int*));
    }

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
                grps->atoms[ng][counter] = atoi(token);
                token = strtok(NULL, sep);
                counter++;
            }
        }
        else {
            n = ret2 - ret1;
            strncpy(grps->names[ng], ret1 +2, n -2);
            if (flag) {
                grps->natoms[ng] = counter;
                ng++;
                counter = 0;
            }
            flag = IN;
        }
    }
    grps->natoms[ng] = counter - 1;
    //nlist[ng] = counter;
    
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

