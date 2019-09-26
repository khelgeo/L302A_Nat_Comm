#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>

#define n_popc 0 
#define n_popcat 46

#define n_sdpc 207
#define n_sdpcat 134

#define n_chol 40
#define n_cholat 29

#define n_water 0
#define n_waterat 3

#define n_eau 144

#define n_sod 2
#define n_sodat 1

#define bin_size_angle 1
#define n_grid_dist 801
#define n_grid_phi 360
#define n_grid_theta 180
#define PI 3.1415
#define NX 101
#define NY 101
#define NZ 101
#define HX 0.1
#define HY 0.1
#define HZ 0.1
#define Origin_x -5.0
#define Origin_y -5.0
#define Origin_z -5.0
#define esmall 0.05
#define n_replace 25
#define window 10
//#define window 100
#define jdist 16
#define jlip 15
#define cut_off_for_dist 3.5
#define cut_off_for_lip 2.0 


void Randomize();
double Random(double);
double GRandom(double, double);

main (int argc, char *argv[]) {

	int **angle_counter, i, j, k, t, n_atoms, bin_dist, bin_angle_phi, n, ***hist_dist_angle, counter_line[9], item, nk[3], n_of_param;
	double aa, putnumber, popc_data[400][50][3], sdpc_data[207][134][3], **bin_angle_theta, water_data[12820][3][3], sod_data[2][1][3], eau_data, chol_data[40][29][3];
	double popc_vel, chol_vel, water_vel, sod_vel, protein_vel, eau_vel, new_popc_data[207][52][3], new_chol_data[50][49][3], min_z, max_z, junk;
	double *bin_size_dist, *final_count_water, *el_density_chol, **count_water, **count_popc, *final_count_popc, *final_count_chol, **count_chol;
	double protein_com[3], angle_phi, angle_theta, x_coor, y_coor, theta, *data, PNz_SDPC, PNz_POPC, PNz_CHOL, *el_density_water, *el_density_popc;
	double box_x, box_y, box_z, dist_x, dist_y, dist_z, *r, distance, const_fact, rlower, rupper, nideal, ***grdf, rho, average_dist, phi, *dz, delta_z;
	double *grid_x, *grid_y, *grid_rdf, *grid_z, temp_element, min_x, min_y, max_x, max_y, **grdf_mat, tempX, tempY, tempZ, distance2D, z_coor, ala_in_window2, ala_in_window3, ala_in_window4;
	int CHOL_O, POPC_P, SDPE_C2, SDPC_P, POPC_N, SDPC_N, POPC_C34, SDPC_C34, CHOL_C27, *mark_lip, n_of_files, n_slice, CHOL_C3, CHOL_C13;
	double vect1_x, vect1_y, vect1_z, arccos, distance1, *round_max, *round_min, angle_final, dist_res, counter_ala, counter_leu, *dist_res_ala, *dist_res_ala2, *dist_res_ala3, *dist_res_leu; 
	double *min, *max, *dist_res_ala4, *dist_res_ala5, *dist_res_ala6, *dist_res_ala7, *dist_res_ala8, *dist_res_ala9, ala_in_window, leu_in_window;
	double ala_in_window5, ala_in_window6, ala_in_window7, ala_in_window8, ala_in_window9, ala_in_window10, *dist_res_ala10, *data_in_window, **data_in_file;
	char *popc_mol_name, *popc_atom_name[50];
        int *popc_mol_number, popc_atom_number;

        char *sdpc_mol_name, *sdpc_atom_name[134];
        int *sdpc_mol_number, sdpc_atom_number;

	char *chol_mol_name, *chol_atom_name[49];
	int *chol_mol_number, chol_atom_number;

	char *water_mol_name, *water_atom_name[3];
	int *water_mol_number, water_atom_number;

        char *eau_mol_name, *eau_atom_name[3];
        int *eau_mol_number, eau_atom_number;

        char *sod_mol_name, *sod_atom_name[1];
        int *sod_mol_number, sod_atom_number, frame, *n_bins, count_lip_numbers, count_gate_breaking;

	char first[80];

	FILE *paramptr[16], *inptr, *incholptr, *eldppcptr, *elwatptr, *paraptr, *elcholptr, *outptr, *leuptr;

	Randomize();

        if ((paraptr = fopen(argv[1], "r")) == NULL)
                printf ("\nERROR - Cannnot open the trajectory file\n");

        fscanf (paraptr, "%d\n", &n_of_files);
        fscanf (paraptr, "%d\n", &n_of_param);
	n_of_param = n_of_param + 1;

	POPC_P = 7;
	POPC_C34 = 49;
	CHOL_O = 0;
	CHOL_C27 = 48; 
        CHOL_C3 = 2;
        CHOL_C13 = 16;

	mark_lip = (int *) malloc(n_popc*(sizeof(int)));
	popc_mol_number = (int *) malloc(n_popc*(sizeof(int)));
	sdpc_mol_number = (int *) malloc(n_sdpc*(sizeof(int)));
	chol_mol_number = (int *) malloc(n_chol*(sizeof(int)));
	water_mol_number = (int *) malloc(n_water*(sizeof(int)));
	eau_mol_number = (int *) malloc(n_eau*(sizeof(int)));
	sod_mol_number = (int *) malloc(n_sod*(sizeof(int)));

        data_in_file = (double **) malloc(n_of_files * sizeof(double *));
        data_in_file[0] = malloc(n_of_files * n_of_param * sizeof(double));
        for(i = 1; i < n_of_files; ++i)
                data_in_file[i] = data_in_file[0] + i * n_of_param;

	data_in_window = (double *) malloc(n_of_param*(sizeof(double)));
	min = (double *) malloc(n_of_param*(sizeof(double)));
	max = (double *) malloc(n_of_param*(sizeof(double)));
	round_max = (double *) malloc(n_of_param*(sizeof(double)));
	round_min = (double *) malloc(n_of_param*(sizeof(double)));
	n_bins = (int *) malloc(n_of_param*(sizeof(int)));
	bin_size_dist = (double *) malloc(n_of_param*(sizeof(double)));

	popc_mol_name = (char *) malloc(5*(sizeof(char)));
	sdpc_mol_name = (char *) malloc(5*(sizeof(char)));
	chol_mol_name = (char *) malloc(5*(sizeof(char)));
	water_mol_name = (char *) malloc(5*(sizeof(char)));
	eau_mol_name = (char *) malloc(5*(sizeof(char)));
	sod_mol_name = (char *) malloc(5*(sizeof(char)));

	for (i=0; i<n_popcat; ++i) 
		popc_atom_name[i] = (char *) malloc(5*(sizeof(char)));
	for (i=0; i<n_sdpcat; ++i) 
		sdpc_atom_name[i] = (char *) malloc(5*(sizeof(char)));
	for (i=0; i<n_cholat; ++i) 
		chol_atom_name[i] = (char *) malloc(5*(sizeof(char)));
	for (i=0; i<n_waterat; ++i) 
		water_atom_name[i] = (char *) malloc(5*(sizeof(char)));
        for (i=0; i<n_waterat; ++i)
                eau_atom_name[i] = (char *) malloc(5*(sizeof(char)));
        for (i=0; i<n_sodat; ++i)
                sod_atom_name[i] = (char *) malloc(5*(sizeof(char)));


	if ((inptr = fopen(argv[2], "r")) == NULL)
                printf ("\nERROR - Cannnot open the trajectory2 file\n");

	for (i=0; i<n_of_files; ++i) {
		for (j=0; j<n_of_param; ++j) {
			fscanf (inptr, "%lf ", &data_in_file[i][j]);
		}
		fscanf (inptr, "\n");
	}
		
	count_gate_breaking = 0;
	count_lip_numbers = 0;
	
	for (i=0; i<n_of_files; ++i) {
		if (data_in_file[i][jdist] > cut_off_for_dist) count_gate_breaking = count_gate_breaking + 1;
		if (data_in_file[i][jlip] > cut_off_for_lip) count_lip_numbers = count_lip_numbers + 1;
	}
	
	printf ("%lf  %lf\n", 1.0*count_gate_breaking/(1.0*n_of_files), 1.0*count_lip_numbers/(1.0*n_of_files));



}

/************************************************************************************************************/

double Random(double range)
{
  return (random()*range)/2147483647L;
}

/************************************************************************************************************/

void Randomize()
{
  srandom(time(NULL));
}

/************************************************************************************************************/

double GRandom(double mean, double sigma)
{
  double v1 = Random(2.0) - 1.0;
  double v2 = Random(2.0) - 1.0;
  double rsq = v1*v1 + v2*v2;

  while (rsq >= 1.0 || rsq == 0.0)
    {
      v1 = Random(2.0) - 1.0;
      v2 = Random(2.0) - 1.0;
      rsq = v1*v1 + v2*v2;
    }

  return mean+v1*sigma*sqrt(-2.0*log(rsq)/rsq);
}

