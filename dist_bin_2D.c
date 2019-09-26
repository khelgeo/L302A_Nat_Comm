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
#define n_of_param 1


void Randomize();
double Random(double);
double GRandom(double, double);

main (int argc, char *argv[]) {

	int ***angle_counter, i, j, k, t, n_atoms, bin_dist, bin_angle_phi, n, ***hist_dist_angle, counter_line[9], item, nk[3];
	double aa, putnumber, popc_data[400][50][3], sdpc_data[207][134][3], **bin_angle_theta_TM4, **bin_angle_theta_TM3, water_data[12820][3][3], sod_data[2][1][3], eau_data, chol_data[40][29][3];
	double popc_vel, chol_vel, water_vel, sod_vel, protein_vel, eau_vel, new_popc_data[207][52][3], new_chol_data[50][49][3], min_TM4_z, max_TM4_z, junk;
	double *bin_size_dist, *final_count_water, *el_density_chol, **count_water, **count_popc, *final_count_popc, *final_count_chol, **count_chol;
	double protein_com[3], angle_phi, angle_theta, x_coor, y_coor, theta, *data, PNz_SDPC, PNz_POPC, PNz_CHOL, *el_density_water, *el_density_popc, *max_TM3, *min_TM3;
	double box_x, box_y, box_z, dist_x, dist_y, dist_z, *r, distance, const_fact, rlower, rupper, nideal, ***grdf, rho, *average_val, phi, *dz, delta_z;
	double *grid_x, *grid_y, *grid_rdf, *grid_z, temp_element, min_TM4_x, *round_max_TM3, *round_min_TM3, min_TM4_y, max_TM4_x, max_TM4_y, **grdf_mat, tempX, tempY, tempZ, distance2D, z_coor, ala_in_window2, ala_in_window3, ala_in_window4;
	int CHOL_O, POPC_P, SDPE_C2, SDPC_P, POPC_N, SDPC_N, POPC_C34, SDPC_C34, CHOL_C27, *mark_lip, n_of_files, n_slice, CHOL_C3, CHOL_C13;
	double vect1_x, vect1_y, vect1_z, arccos, distance1, *round_max_TM4, *round_min_TM4, angle_final, dist_res, counter_ala, counter_leu, *dist_res_ala, *dist_res_ala2, *dist_res_ala3, *dist_res_leu; 
	double *min_TM4, *max_TM4, *dist_res_ala4, *dist_res_ala5, *dist_res_ala6, *dist_res_ala7, *dist_res_ala8, *dist_res_ala9, ala_in_window, leu_in_window, **data_in_file_TM3;
	double ala_in_window5, ala_in_window6, ala_in_window7, ala_in_window8, ala_in_window9, ala_in_window10, *dist_res_ala10, *data_in_window, **data_in_file_TM4;
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
        int *sod_mol_number, sod_atom_number, frame, *n_bins_TM4, *n_bins_TM3;

	char first[80];

	FILE *paramptr[24], *inptr, *incholptr, *eldppcptr, *elwatptr, *paraptr, *elcholptr, *outptr, *leuptr;

	Randomize();

        if ((paraptr = fopen(argv[1], "r")) == NULL)
                printf ("\nERROR - Cannnot open the trajectory file\n");

        fscanf (paraptr, "%d\n", &n_of_files);

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

        data_in_file_TM4 = (double **) malloc(n_of_files * sizeof(double *));
        data_in_file_TM4[0] = malloc(n_of_files * n_of_param * sizeof(double));
        for(i = 1; i < n_of_files; ++i)
                data_in_file_TM4[i] = data_in_file_TM4[0] + i * n_of_param;


        data_in_file_TM3 = (double **) malloc(n_of_files * sizeof(double *));
        data_in_file_TM3[0] = malloc(n_of_files * n_of_param * sizeof(double));
        for(i = 1; i < n_of_files; ++i)
                data_in_file_TM3[i] = data_in_file_TM3[0] + i * n_of_param;



	data_in_window = (double *) malloc(n_of_param*(sizeof(double)));
	average_val = (double *) malloc(n_of_param*(sizeof(double)));
	min_TM4 = (double *) malloc(n_of_param*(sizeof(double)));
	min_TM3 = (double *) malloc(n_of_param*(sizeof(double)));
	max_TM4 = (double *) malloc(n_of_param*(sizeof(double)));
	max_TM3 = (double *) malloc(n_of_param*(sizeof(double)));
	round_max_TM4 = (double *) malloc(n_of_param*(sizeof(double)));
	round_max_TM3 = (double *) malloc(n_of_param*(sizeof(double)));
	round_min_TM4 = (double *) malloc(n_of_param*(sizeof(double)));
	round_min_TM3 = (double *) malloc(n_of_param*(sizeof(double)));
	n_bins_TM4 = (int *) malloc(n_of_param*(sizeof(int)));
	n_bins_TM3 = (int *) malloc(n_of_param*(sizeof(int)));
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
			fscanf (inptr, "%lf %lf %lf %lf ", &junk, &data_in_file_TM4[i][j], &junk, &data_in_file_TM3[i][j]);
		}
		fscanf (inptr, "\n");
	}
		

		for (j=0; j<n_of_param; ++j) {
        		max_TM4[j] = -5000;
        		max_TM3[j] = -5000;
        		min_TM4[j] = 5000;
        		min_TM3[j] = 5000;
			bin_size_dist[j]=2.0;
			average_val[j]=0;
		}

	for (j=0; j<n_of_param; ++j) {
        	for (i=0; i<n_of_files; ++i) {
                	if (data_in_file_TM4[i][j] > max_TM4[j]) max_TM4[j] = data_in_file_TM4[i][j];
                	if (data_in_file_TM4[i][j] < min_TM4[j]) min_TM4[j] = data_in_file_TM4[i][j];
                	if (data_in_file_TM3[i][j] > max_TM3[j]) max_TM3[j] = data_in_file_TM3[i][j];
                	if (data_in_file_TM3[i][j] < min_TM3[j]) min_TM3[j] = data_in_file_TM3[i][j];
        	}
	}


	for (j=0; j<n_of_param; ++j) {
        	round_min_TM4[j] = floor(min_TM4[j]) - 5;
        	round_max_TM4[j] = ceil(max_TM4[j]) + 5;
        	round_min_TM3[j] = floor(min_TM3[j]) - 5;
        	round_max_TM3[j] = ceil(max_TM3[j]) + 5;
        	n_bins_TM4[j] = (round_max_TM4[j]-round_min_TM4[j])/bin_size_dist[j];
        	n_bins_TM3[j] = (round_max_TM3[j]-round_min_TM3[j])/bin_size_dist[j];
//		printf ("%d %d %lf %lf %lf %lf\n", n_bins_TM3[j], n_bins_TM4[j], round_min_TM3[j], round_max_TM3[j], round_min_TM4[j], round_max_TM4[j]);
	}

        bin_angle_theta_TM4 = (double **) malloc(500 * sizeof(double *));
        bin_angle_theta_TM4[0] = malloc(500 * n_of_param * sizeof(double));
        for(i = 1; i < 500; ++i)
                bin_angle_theta_TM4[i] = bin_angle_theta_TM4[0] + i * n_of_param;

        bin_angle_theta_TM3 = (double **) malloc(500 * sizeof(double *));
        bin_angle_theta_TM3[0] = malloc(500 * n_of_param * sizeof(double));
        for(i = 1; i < 500; ++i)
                bin_angle_theta_TM3[i] = bin_angle_theta_TM3[0] + i * n_of_param;

        angle_counter = (int ***) malloc( 50 * sizeof(int **) );
        angle_counter[0] = (int **) malloc( 50 * 50 * sizeof( int *) );
        angle_counter[0][0] = (int *) malloc( 50 * 50 * 50 * sizeof(int) );
        for ( j = 1  ;  j < 50  ;  ++j ) {
            angle_counter[0][j] = angle_counter[0][j-1] + 50;
        }
        for ( i = 1  ;  i < 50  ;  ++i ) {
            angle_counter[i] = angle_counter[i-1] + 50;
            angle_counter[i][0] = angle_counter[i-1][50-1] + 50;
            for ( j = 1  ;  j < 50  ;  ++j ) {
                angle_counter[i][j] = angle_counter[i][j-1] + 50;
            }
        }

/*        cof = (double ***) malloc( nk[2] * sizeof(double **) );
        cof[0] = (double **) malloc( nk[2] * nk[1] * sizeof( double *) );
        cof[0][0] = (double *) malloc( nk[2] * nk[1] * nk[0] * sizeof(double) );
        for ( j = 1  ;  j < nk[1]  ;  ++j ) {
            cof[0][j] = cof[0][j-1] + nk[0];
        }
        for ( i = 1  ;  i < nk[2]  ;  ++i ) {
            cof[i] = cof[i-1] + nk[1];
            cof[i][0] = cof[i-1][nk[1]-1] + nk[0];
            for ( j = 1  ;  j < nk[1]  ;  ++j ) {
                cof[i][j] = cof[i][j-1] + nk[0];
            }
        }*/


	for (k=0; k<n_of_param; ++k) {
		for (j=0; j<50; ++j) {
			for (i=0; i<50; ++i) {
				angle_counter[k][i][j]=0;
			}
		}
	}

	for (j=0; j<n_of_param; ++j) {
        	for (i=0; i<n_bins_TM4[j]; ++i) {
                	bin_angle_theta_TM4[i][j] = round_min_TM4[j] + i*bin_size_dist[j];
        	}
	}



        for (j=0; j<n_of_param; ++j) {
                for (i=0; i<n_bins_TM3[j]; ++i) {
                        bin_angle_theta_TM3[i][j] = round_min_TM3[j] + i*bin_size_dist[j];
                }
        }


	for (k=0; k<n_of_param; ++k) {
        	for (i=0; i<n_of_files; ++i) {
                	for (j=0; j<n_bins_TM4[k]; ++j) {
                		for (t=0; t<n_bins_TM3[k]; ++t) {

                        		if (bin_angle_theta_TM4[j][k] <= data_in_file_TM4[i][k] && data_in_file_TM4[i][k] < bin_angle_theta_TM4[j+1][k] && bin_angle_theta_TM3[t][k] <= data_in_file_TM3[i][k] && data_in_file_TM3[i][k] < bin_angle_theta_TM3[t+1][k]) {
	//				printf ("%d %lf\n", j, bin_angle_theta_TM4[j][k]);	
					angle_counter[k][t][j] = angle_counter[k][t][j] + 1;
					}
				}

                	}
		}
        }

	for (k=0; k<n_of_param; ++k) {
       		if ((paramptr[k] = fopen(argv[3+k], "w")) == NULL)
                printf ("\nERROR - Cannnot open the trajectory file %d\n", k);
	}
	
	for (k=0; k<n_of_param; ++k) {
        	for (i=0; i<n_bins_TM3[k]; ++i) {
        		for (j=0; j<n_bins_TM4[k]; ++j) {
                		aa = 1.0*angle_counter[k][i][j]/n_of_files;
			//	if (angle_counter[k][i][j] > 0) printf ("%d %lf %d\n", j, bin_angle_theta_TM4[j][k], angle_counter[k][i][j]);
                		fprintf (paramptr[k], "%lf %lf %lf %d\n", bin_angle_theta_TM3[i][k], bin_angle_theta_TM4[j][k], aa, angle_counter[k][i][j]);
        		}
		}
	}


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

