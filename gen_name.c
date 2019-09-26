#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define header_lines 11
#define title 2
#define trailor 5
#define PI 3.1415
#define eps_die 78.54
#define R 8.314
#define kbol 1.38
#define e_charge 1.6
#define eps_zero 8.85
#define T 300.0
#define radius_ion 10
#define initial_energy 30345.39054428783736704468
#define initial_energy_bending 0.0
#define b_factor 2.494353
#define traj_number 1604  

main (int argc, char *argv[])

{
        int i, j, k, grid_points_ion, grid_points_surf, nx_ion, ny_ion, nz_ion, nx_surf, ny_surf, nz_surf;
        int replace_x, replace_y, replace_z, n_dir, *n_lines_in_dir;
        double origin_x_ion, origin_y_ion, origin_z_ion, origin_x_surf, origin_y_surf, origin_z_surf, dist_surf_ion, energy_new;
        double adsorption_energy, total_energy, bend_energy;
        double *data_ion, *data_surf, *data_output, counter, check_dist, radius, hx, hy, hz, z_low, z_high, z_surf, el_energy, fourth_term, mixing_term;
        char headers[11][150], tails[5][150];
        FILE *in1ptr, *in2ptr, *ionparaptr, *surfparaptr, *outptr, *paraptr;

        if ((paraptr = fopen(argv[1], "r")) == NULL)
                printf ("\nERROR - Cannnot open the trajectory file\n");

	                fscanf (paraptr, "%d", &n_dir);
	n_lines_in_dir = (int *) malloc(n_dir*(sizeof(int)));
        for (i=0;i<n_dir;i=i+1) {
		fscanf (paraptr, "%d", &n_lines_in_dir[i]);
	}

        printf ("./sasa_residue parameters_lines.dat sasa_average_prot.dat sasa_protA.dat sasa_protB.dat ");

	for (j=0;j<n_dir;j=j+1) {
	        for (i=0;i<n_lines_in_dir[j];i=i+1) {
        	                printf ("%d/prot.%d.rsa ", j+1, i);
	        }
	}

	printf ("\n");
        printf ("./sasa_residue parameters_lines.dat sasa_average_memb_prot.dat sasa_memb_protA.dat sasa_memb_protB.dat ");

        for (j=0;j<n_dir;j=j+1) {
                for (i=0;i<n_lines_in_dir[j];i=i+1) {
                                printf ("%d/prot_memb.%d.rsa ", j+1, i);
                }
        }







}


