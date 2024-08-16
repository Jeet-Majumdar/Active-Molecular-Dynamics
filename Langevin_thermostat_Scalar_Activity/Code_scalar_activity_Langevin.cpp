#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <random>
#include <ctime>
#include <stdio.h>
#include <omp.h>
#include <cstring>

#define seed 23410981

using namespace std;

// Compile using "code.cpp -lm -fopenmp -o parallel.o"

/*
Declaration of all Glabal Variables
*/

//==============================================================================
//                           Global Physical Constants
//==============================================================================


const double Avogadro = 1.0; 
const double Boltzmann = 1.0; 

//************************************************ Functions Begin Here *****************************************************

double smallest(double x, double y, double z){
  return x < y ? (x < z ? x : z) : (y < z ? y : z);
}

//
// Generate a random number between 0 and 1
// return a uniform number in [0,1].
double unifRand()
{
    return rand() / double(RAND_MAX);
}
//
// Generate a random number in a real interval.
// param a one end point of the interval
// param b the other end of the interval
// return a inform rand numberin [a,b].
double unifRand(double a, double b)
{
    return (b-a)*unifRand() + a;
}
//
// Generate a random integer between 1 and a given value.
// param n the largest value
// return a uniform random value in [1,...,n]
long unifRand(long n)
{

    if (n < 0) n = -n;
    if (n==0) return 0;
    long guard = (long) (unifRand() * n) +1;
    return (guard > n)? n : guard;
}
//
// Reset the random number generator with the system clock.
void Seed()
{
    srand(time(0));
}


/* Prints usage information */
void usage ( void ) {
  fprintf(stdout,"Scalar_Activity_Berendsen usage:\n");
  fprintf(stdout,"Scalar_Activity_Berendsen.o [options]\n\n");
  fprintf(stdout,"Options:\n");
  fprintf(stdout,"\t -rho             [real]\t\tNumber density [Default: 0.8]\n");
  fprintf(stdout,"\t -relax           [real]\t\tRelax constant for Langevin Thermostat [Default: 10]\n");
  fprintf(stdout,"\t -dt              [real]\t\tTime step [Default: 0.001]\n");
  fprintf(stdout,"\t -ns              [integer]\t\tNumber of integration steps [Default: 5000]\n");
  fprintf(stdout,"\t -p               [real]\t\tPercentage of Hot particles [Default: 100]\n");
  fprintf(stdout,"\t -Th              [real]\t\tTemperature of Hot particles [Default: 5.0]\n");
  fprintf(stdout,"\t -Tc              [real]\t\tTemperature of Cold particles [Default: 5.0]\n");
  fprintf(stdout,"\t -dumpFreq        [integer]\t\tDump frequency [Default: 10]\n");
  fprintf(stdout,"\t -dumpFilename    [string]\t\tDump File name [Default: dump.lammpstrj] \n");
  fprintf(stdout,"\t -r               [integer]\t\tRestart Flag: 1 = Start from an old configuration (provide -restartFilename) [Default: 0]\n");
  fprintf(stdout,"\t -restartFilename [string]\t\tRestart File name if -r is 1 [Default: restart.jeet] \n");
}


double** create_2D_Array(int m, int n){
	double** array = (double**)malloc(m*sizeof(double*));
	for(int i = 0; i<m; i++){
		array[i] = (double*)malloc(n*sizeof(double));
	}
	return array;
	//access with: double **force_F = create2D_Array(rows, column)
}


void delete_2D_Array(double **array, int m, int n){
	for(int i = 0; i<m; i++){
		free(array[i]);
	}
	free(array);
}

void readInputData(FILE *&filename_F, int &natoms_F, int dim_F, int *&group_F, double **&box_F, double **&pos_F, double **&vels_F, double *&mass_F){
	char it1[10],time1[10], it2[20],num[10],of[10],atom1[10],it3[10],box[10],bound[10],pp1[10],pp2[10],pp3[10],it4[10],atom2[10],id[10],type[10],masses[10],XU[10],YU[10],ZU[10],VX[10],VY[10],VZ[10],act_x[10],act_y[10],act_z[10];
	int s1,s2, id_temp, type_temp;

	fscanf(filename_F, "%s %s", it1, time1);
	fscanf(filename_F,"%d",&s1);
	fscanf(filename_F,"%s %s %s %s",it2,num,of,atom1);
	fscanf(filename_F,"%d",&s2);
	fscanf(filename_F,"%s %s %s %s %s %s",it3,box,bound,pp1,pp2,pp3);
	fscanf(filename_F,"%lf %lf",&box_F[0][0], &box_F[0][1]);
	fscanf(filename_F,"%lf %lf",&box_F[1][0], &box_F[1][1]);
	fscanf(filename_F,"%lf %lf",&box_F[2][0], &box_F[2][1]);
	fscanf(filename_F,"%s %s %s %s%s %s %s %s %s %s",it4,atom2,id,type,XU,YU,ZU,VX,VY,VZ);
	//fscanf(filename_F,"%s %s %s %s %s %s %s %s %s %s %s",it4,atom2,id,type,masses,XU,YU,ZU,VX,VY,VZ); // In case you want to print mass too.
	
	if(natoms_F != s2) perror("Wrong Restart Fiile!!");
	natoms_F = s2;
	for(int j=0;j<natoms_F;j++)
	{
	fscanf(filename_F,"%d %d %lf %lf %lf %lf %lf %lf",&id_temp, &group_F[j], &pos_F[j][0], &pos_F[j][1], &pos_F[j][2], &vels_F[j][0], &vels_F[j][1], &vels_F[j][2]);
	mass_F[j] = 1.0;
	//printf("%d %d %lf %lf %lf %lf %lf %lf\n",id_temp, group_F[j], pos_F[j][0], pos_F[j][1], pos_F[j][2], vels_F[j][0], vels_F[j][1], vels_F[j][2]);
	}
}

void writeOutput(FILE *&filename_F, int step_F, int natoms_F, int dim_F, double **&box_F, double **&pos_F, double **&vels_F, int *&group_F){


	fprintf(filename_F, "ITEM: TIMESTEP\n");
	fprintf(filename_F, "%4d\n", step_F);

	fprintf(filename_F, "ITEM: NUMBER OF ATOMS\n");
	fprintf(filename_F, "%4d\n", natoms_F);

	fprintf(filename_F, "ITEM: BOX BOUNDS  pp pp pp \n");
	for(int i = 0; i<dim_F; i++){
		fprintf(filename_F, "%0.5f %0.5f\n", box_F[i][0], box_F[i][1]);
	}
	if (dim_F==2)  fprintf(filename_F,"0.0 0.0 \n");

	fprintf(filename_F,"ITEM: ATOMS ID type x y z vx vy vz\n");

	for(int i=0; i<natoms_F; i++){
		fprintf(filename_F, " %5d %2d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n", (i+1), group_F[i], pos_F[i][0], pos_F[i][1], pos_F[i][2], vels_F[i][0], vels_F[i][1], vels_F[i][2]);
        }
}

void variable_Initialization(double **&pos_F, double **&vels_F, int *&group_F, double **&forces_F, double **&box_F,\
	double *&mass_F, double *&radius_F, int natoms_F, int dim_F, double mass_atomic_F, double rcut_off_F, \
	double radii_F, double &ke_Energy_hot_F, double temp_hot_F, double relax_F, double dt_F, double epsilon_F, \
	double sigma_F, double temp_cold_F, double &ke_Energy_cold_F, int n_hot_atoms_F){

	std::default_random_engine generator (seed);
	std::normal_distribution<double> distribution (0.0,1.0);
	
	int n_temp=2;
	while ((n_temp*n_temp*n_temp) < natoms_F) n_temp++;
	
	int ii_dim[3]= {0, 0, 0};	
	/* Assign particle positions */
	for (int i=0;i<natoms_F;i++) {
		for(int j=0; j<dim_F; j++)   pos_F[i][j] = ((double)ii_dim[j]+0.5)*(box_F[j][1] - box_F[j][0])/n_temp;
		ii_dim[0]++;
		if (ii_dim[0]==n_temp) {
			ii_dim[0]=0;
		        ii_dim[1]++;
			if (ii_dim[1]==n_temp) {
				ii_dim[1]=0;
				ii_dim[2]++;
		        }
	      }
    	}


	// Initialize Velocity with random variable from a gaussian distribution and Fill up Group array simultaneously.

		for(int i=0; i<natoms_F; i++){
			for(int j=0; j<dim_F; j++){
				vels_F[i][j] = distribution(generator);
			}
			if(i < n_hot_atoms_F) 
				group_F[i] = 1;		// First n_hot_atoms_F as HOT
			else	
				group_F[i] = 2;		// Remaining atoms as COLD
		}

	// Subtract Velocity of COM from all
		double avg_vels[3] = {0.0, 0.0, 0.0};
		for(int j=0; j<dim_F; j++){
    		for(int i=0; i<natoms_F; i++){
    			avg_vels[j] = avg_vels[j] + vels_F[i][j];
			}
		}

		ke_Energy_hot_F = 0;
		ke_Energy_cold_F = 0;
		for(int j=0; j<dim_F; j++){
    		for(int i=0; i<natoms_F; i++){
    		vels_F[i][j] -= avg_vels[j]/natoms_F;
    		if(i < n_hot_atoms_F) 
    			ke_Energy_hot_F += vels_F[i][j]*vels_F[i][j];
    		else
    			ke_Energy_cold_F += vels_F[i][j]*vels_F[i][j];
		}
		}
		
		ke_Energy_hot_F *= 0.5; 
		ke_Energy_cold_F *= 0.5;

	// Scale velocity to a given starting temperature by default
		double scale_T_hot = 2.0*(ke_Energy_hot_F)/n_hot_atoms_F/((double)dim_F);
		double factor_scale_vel_hot = sqrt(temp_hot_F/scale_T_hot);
		double scale_T_cold = (natoms_F - n_hot_atoms_F)>0?2.0*(ke_Energy_hot_F)/(natoms_F - n_hot_atoms_F)/((double)dim_F):0.0;
		double factor_scale_vel_cold = (natoms_F - n_hot_atoms_F)>0?sqrt(temp_cold_F/scale_T_cold):0.0;
		ke_Energy_hot_F = 0;	
		ke_Energy_cold_F = 0;
		for(int j=0; j<dim_F; j++){
    		for(int i=0; i<natoms_F; i++){
	    		if(i < n_hot_atoms_F){
	    			vels_F[i][j] *= factor_scale_vel_hot;
	    			ke_Energy_hot_F += vels_F[i][j]*vels_F[i][j];
	    		}else{
	    			vels_F[i][j] *= factor_scale_vel_cold;
	    			ke_Energy_cold_F += vels_F[i][j]*vels_F[i][j];
	    		}
		}
		}
		ke_Energy_hot_F *= 0.5;
		ke_Energy_cold_F *= 0.5;

	// Initialize Mass

		for(int i=0; i<natoms_F; i++){
			mass_F[i] = mass_atomic_F;
		}

	// Initialize Radius
		for(int i=0; i<natoms_F; i++){
			radius_F[i] = 1 * radii_F;
		}
}



void compute_Force(double **&box_F, double **&pos_F, double **&vels_F, double *&mass_F, double **&force_F, double rcut_off_F, double epsilon_F, double sigma_F, int natoms_F, double &pot_Energy_F, double &pressure_virial_F, double energy_correction_F,  double energy_at_cutoff_F){
		
	
	pot_Energy_F = 0.0;
	pressure_virial_F = 0.0;
	double rij[3];
	
	

	#pragma omp parallel default(shared) private(rij)
	{
	#pragma omp for
	for (int i=0; i<(natoms_F-1); i++) {
		double tempo_virial = 0.0;
		double tempo_pot_E = 0.0;
	    for (int j=i+1; j<(natoms_F); j++) {
		rij[0] = pos_F[i][0] - pos_F[j][0];
		rij[1] = pos_F[i][1] - pos_F[j][1];
		rij[2] = pos_F[i][2] - pos_F[j][2];

		if(rij[0] > (box_F[0][1] - box_F[0][0])*0.5)	rij[0] -= (box_F[0][1] - box_F[0][0]);
		else if(rij[0] < -1.0*(box_F[0][1] - box_F[0][0])*0.5)	rij[0] += (box_F[0][1] - box_F[0][0]);

		if(rij[1] > (box_F[1][1] - box_F[1][0])*0.5)	rij[1] -= (box_F[1][1] - box_F[1][0]);
		else if(rij[1] < -1.0*(box_F[1][1] - box_F[1][0])*0.5)	rij[1] += (box_F[1][1] - box_F[1][0]);

		if(rij[2] > (box_F[2][1] - box_F[2][0])*0.5)	rij[2] -= (box_F[2][1] - box_F[2][0]);
		else if(rij[2] < -1.0*(box_F[2][1] - box_F[2][0])*0.5)	rij[2] += (box_F[2][1] - box_F[2][0]);

		double r2 = (rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);

		if(r2 < rcut_off_F * rcut_off_F){	
		  double r6i   = 1.0/(r2*r2*r2);
		  double r_by_sigma2 = r2/sigma_F/sigma_F;
		  double r_by_sigma6i   = 1.0/(r_by_sigma2*r_by_sigma2*r_by_sigma2);

		  tempo_pot_E    += 4.0*epsilon_F*(r_by_sigma6i*r_by_sigma6i - r_by_sigma6i) - energy_at_cutoff_F;
		  double factor_F = 48.0*epsilon_F*(pow(sigma_F,12)*r6i*r6i - 0.5*pow(sigma_F,6)*r6i);

		  force_F[i][0] += rij[0]*factor_F/r2;
		  force_F[j][0] -= rij[0]*factor_F/r2;
		  force_F[i][1] += rij[1]*factor_F/r2;
		  force_F[j][1] -= rij[1]*factor_F/r2;
		  force_F[i][2] += rij[2]*factor_F/r2;
		  force_F[j][2] -= rij[2]*factor_F/r2;	
		  tempo_virial +=  factor_F;
		}
	    }

	    #pragma omp critical
		{
			pot_Energy_F += tempo_pot_E;
			pressure_virial_F += tempo_virial;
		}

		}

		
	}
	
	pot_Energy_F = 	pot_Energy_F + natoms_F*energy_correction_F;
}


void periodic_check(double **&pos_F, double **&box_F, int natoms_F){
        for(int i=0; i<natoms_F; i++){

		
                if(pos_F[i][0] > box_F[0][1])
                        pos_F[i][0] = fmod ((pos_F[i][0] - box_F[0][0]),(box_F[0][1] - box_F[0][0])) + box_F[0][0];                    
		else if(pos_F[i][0] < box_F[0][0])
                        pos_F[i][0] = fmod ((pos_F[i][0] - box_F[0][0]),(box_F[0][1] - box_F[0][0])) + box_F[0][1];
                else    pos_F[i][0] = pos_F[i][0];

                if(pos_F[i][1] > box_F[1][1])
                        pos_F[i][1] = fmod ((pos_F[i][1] - box_F[1][0]),(box_F[1][1] - box_F[1][0])) + box_F[1][0];
                else if(pos_F[i][1] < box_F[1][0])
                        pos_F[i][1] = fmod ((pos_F[i][1] - box_F[1][0]),(box_F[1][1] - box_F[1][0])) + box_F[1][1];
                else    pos_F[i][1] = pos_F[i][1];

                if(pos_F[i][2] > box_F[2][1])
                        pos_F[i][2] = fmod ((pos_F[i][2] - box_F[2][0]),(box_F[2][1] - box_F[2][0])) + box_F[2][0];
                else if(pos_F[i][2] < box_F[2][0])
                        pos_F[i][2] = fmod ((pos_F[i][2] - box_F[2][0]),(box_F[2][1] - box_F[2][0])) + box_F[2][1];
                else    pos_F[i][2] = pos_F[i][2];



        }

}


void integrate(int natoms_F, double **&pos_F, double **&vels_F, double *&mass_F, double **&force_F, \
	int *&group_F, double relax_F, double dt_F, double rcut_off_F, \
	double epsilon_F, double sigma_F, double &pot_Energy_F, double **&box_F, double gfric_F, double noise_hot_F, \
	double noise_cold_F, double &pressure_virial_F, double energy_correction_F, double energy_at_cutoff_F){

	std::default_random_engine generator (seed);
        std::normal_distribution<double> distribution (0.0,1.0);

	for(int i=0; i<natoms_F; i++){

		pos_F[i][0] += vels_F[i][0]*dt_F + 0.5*dt_F*dt_F*force_F[i][0]/mass_F[i];
		pos_F[i][1] += vels_F[i][1]*dt_F + 0.5*dt_F*dt_F*force_F[i][1]/mass_F[i];
		pos_F[i][2] += vels_F[i][2]*dt_F + 0.5*dt_F*dt_F*force_F[i][2]/mass_F[i];

		vels_F[i][0] = vels_F[i][0] - vels_F[i][0]*relax_F*0.5*dt_F/mass_F[i] + 0.5*dt_F*force_F[i][0]/mass_F[i];
		vels_F[i][1] = vels_F[i][1] - vels_F[i][1]*relax_F*0.5*dt_F/mass_F[i] + 0.5*dt_F*force_F[i][1]/mass_F[i];
		vels_F[i][2] = vels_F[i][2] - vels_F[i][2]*relax_F*0.5*dt_F/mass_F[i] + 0.5*dt_F*force_F[i][2]/mass_F[i];
        }

	void periodic_check(double**&, double**&, int);
    	periodic_check(pos_F, box_F, natoms_F);	


	for(int i = 0; i<natoms_F; i++){
		if(group_F[i] == 1){
			force_F[i][0] = 2.0*noise_hot_F*(unifRand() - 0.5);
			force_F[i][1] = 2.0*noise_hot_F*(unifRand() - 0.5);
			force_F[i][2] = 2.0*noise_hot_F*(unifRand() - 0.5);
		}
		else{
			force_F[i][0] = 2.0*noise_cold_F*(unifRand() - 0.5);
			force_F[i][1] = 2.0*noise_cold_F*(unifRand() - 0.5);
			force_F[i][2] = 2.0*noise_cold_F*(unifRand() - 0.5);
		}
	}

     
	void compute_Force(double**&, double**&, double**&, double*&, double**&, double, double, double, int, double &, double &, double, double);
	compute_Force(box_F, pos_F, vels_F, mass_F, force_F, rcut_off_F, epsilon_F, sigma_F, natoms_F, pot_Energy_F, pressure_virial_F, energy_correction_F, energy_at_cutoff_F);
	
	for(int i = 0; i<natoms_F; i++){
		vels_F[i][0] = vels_F[i][0] - vels_F[i][0]*relax_F*0.5*dt_F/mass_F[i] + 0.5*dt_F*force_F[i][0]/mass_F[i];
        vels_F[i][1] = vels_F[i][1] - vels_F[i][1]*relax_F*0.5*dt_F/mass_F[i] + 0.5*dt_F*force_F[i][1]/mass_F[i];
		vels_F[i][2] = vels_F[i][2] - vels_F[i][2]*relax_F*0.5*dt_F/mass_F[i] + 0.5*dt_F*force_F[i][2]/mass_F[i];
	}

}



double compute_Temp(double **&vels_F, double *&mass_F, int dim_F, int natoms_F, int *&group_F , \
	double &ke_Energy_hot_F, double &ke_Energy_cold_F, double &temp_hot_F, double &temp_cold_F){

	int count_hot = 0;
	double avg_vels_hot[3] = {0.0, 0.0, 0.0};
	double avg_vels_cold[3] = {0.0, 0.0, 0.0};
	for(int j=0; j<dim_F; j++){
		count_hot = 0;
    	for(int i=0; i<natoms_F; i++){
    		if(group_F[i] == 1){
    			count_hot = count_hot + 1;
    			avg_vels_hot[j] = avg_vels_hot[j] + vels_F[i][j];
    		}
    		else{
    			avg_vels_cold[j] = avg_vels_cold[j] + vels_F[i][j];
    		}
		}
	if(count_hot!=0) avg_vels_hot[j] = avg_vels_hot[j]/count_hot; else avg_vels_hot[j] = 0.0;
	if((natoms_F - count_hot)!=0) avg_vels_cold[j] = avg_vels_cold[j]/(natoms_F - count_hot); else avg_vels_cold[j] = 0.0;
	}

	double ke_hot = 0.0;
	double ke_cold = 0.0;
	for(int j=0; j<dim_F; j++){
    	for(int i=0; i<natoms_F; i++){
    		if(group_F[i] == 1) ke_hot += mass_F[i]*(vels_F[i][j] - avg_vels_hot[j])*(vels_F[i][j] - avg_vels_hot[j]);
    		else ke_cold += mass_F[i]*(vels_F[i][j] - avg_vels_cold[j])*(vels_F[i][j] - avg_vels_cold[j]);
    	}
	}


	double ke_system = ke_hot + ke_cold;
	if(count_hot!=0) temp_hot_F =  ke_hot / (dim_F * Boltzmann * count_hot); else temp_hot_F = 0.0;
	if((natoms_F - count_hot)!=0) temp_cold_F =  ke_cold / (dim_F * Boltzmann * (natoms_F - count_hot)); else temp_cold_F = 0.0;
	double ins_temp_system =  ke_system / (dim_F * Boltzmann * natoms_F);
        ke_Energy_hot_F = ke_hot * 0.5;
        ke_Energy_cold_F = ke_cold * 0.5;
	return ins_temp_system;
}


//************************************************ Functions End Here *******************************************************




//==============================================================================
//                           Parameters
//==============================================================================

// ======= User Parameters =======
int dim = 3;
int natoms = 864;		// 864;
double rho = 0.8;
double sigma = 1;
double epsilon = 1;
double mass_atomic = 1.0;
double percentage_hot = 100;								// Percentage of Hot particles. Remaining will be cold
double temp_hot = 5.0;
double temp_cold = 5.0;
double relax = 10;
double dt = 0.001;
double rcut_off = 4.0; //1.58113883;   // 2^(1/6)*sigma
int nsteps = 1000;
int output_freq = 10;
const char*  dump_File = "dump.lammpstrj";
const char*  input_File = "data.jeet";					// Input DATA File
const char*  output_File = "restart.jeet";					// Output DATA File.
int restart_flag = 0;									// Flag whether to start from restart file (1 = Yes; 0 = No)
double box_x_low = 0;
double box_x_high = pow((natoms/rho),0.3333333);
double box_y_low = 0;
double box_y_high = pow((natoms/rho),0.3333333);
double box_z_low = 0;
double box_z_high = pow((natoms/rho),0.3333333);
double radii = 1.0;

int main(int argc, char **argv){

/* Here we parse the command line arguments;  If
   you add an option, document it in the usage() function! */
  for (int i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-rho")) rho=atof(argv[++i]);
    else if (!strcmp(argv[i],"-relax")) relax=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dt")) dt=atof(argv[++i]);
    else if (!strcmp(argv[i],"-ns")) nsteps=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-p")) percentage_hot=atof(argv[++i]);
    else if (!strcmp(argv[i],"-Th")) temp_hot = atof(argv[++i]);
    else if (!strcmp(argv[i],"-Tc")) temp_cold=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dumpFreq")) output_freq=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-dumpFilename")) dump_File=(argv[++i]);
    else if (!strcmp(argv[i],"-r")) restart_flag=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-restartFilename")) input_File=(argv[++i]);
    else if (!strcmp(argv[i],"-h")) {
      usage(); exit(0);
    }
    else {
      fprintf(stderr,"Error: Command-line argument '%s' not recognized.\n",
	      argv[i]);
      exit(-1);
    }
  }

        

	box_x_high = pow((natoms/rho),0.3333333);
	box_y_high = pow((natoms/rho),0.3333333);
	box_z_high = pow((natoms/rho),0.3333333);
	
	double **pos;
	double **vels;
	double *mass;
	double **forces;
	double **box;
	int *group;
	double *radius = (double*)malloc(natoms*sizeof(double));

	int n_hot_atoms = (int)percentage_hot*natoms/100;		// natoms - n_hot_atoms = number of cold atoms

    double pot_Energy = 0.0;
    double ke_Energy_hot = 0.0;
    double ke_Energy_cold = 0.0;
	double pressure_virial = 0.0;
	double temp_system = 0.0;
	
	double volume = (box_x_high - box_x_low)*(box_y_high - box_y_low)*(box_z_high - box_z_low);
	double gfric = 1.0-relax*dt*0.5/mass_atomic;
	double noise_hot = sqrt(6.0*relax*temp_hot*Boltzmann/dt);
	double noise_cold = sqrt(6.0*relax*temp_cold*Boltzmann/dt);
	
	double rr3 = 1.0/pow(rcut_off,3);
 	double energy_correction = 8.0*M_PI*rho*pow(sigma,3)*epsilon*(pow(sigma,9)*rr3*rr3*rr3/9.0 - pow(sigma,3)*rr3/3.0);
 	double pressure_correction = 16.0/3.0*M_PI*rho*rho*pow(sigma,3)*epsilon*(2./3.*rr3*rr3*rr3*pow(sigma,9) - pow(sigma,3)*rr3);	
 	double energy_at_cutoff = 4*epsilon*(pow(sigma,12)*rr3*rr3*rr3*rr3 - pow(sigma,6)*rr3*rr3);

	FILE *outfile_Dump;
        outfile_Dump = fopen(dump_File, "w+");

	FILE *outfile;
        outfile = fopen(output_File, "w+");

	if(outfile_Dump==NULL || outfile==NULL) perror("error");

	// Memory allocation
	pos = create_2D_Array(natoms, dim);
	vels = create_2D_Array(natoms, dim);
	box =  create_2D_Array(dim, 2);
   	mass = (double*)malloc(natoms*sizeof(double));
   	group = (int*)malloc(natoms*sizeof(int));			// 1 for HOT, 2 for COLD

	FILE *infile;
     if(restart_flag == 1){ 
     	if(fopen(input_File, "r") == NULL)  perror("error");
     	infile = fopen(input_File, "r");
     	readInputData(infile, natoms, dim, group, box, pos, vels, mass);
     	fclose(infile);
     }else{

		// Initialize Box
		box[0][0] = box_x_low; box[0][1] = box_x_high;
		box[1][0] = box_y_low; box[1][1] = box_y_high;
		if(dim == 3){
			box[2][0] = box_z_low, box[2][1] = box_z_high;
			}
		variable_Initialization(pos, vels, group, forces, box, mass, radius, natoms, dim, mass_atomic, rcut_off, radii, ke_Energy_hot, temp_hot, relax, dt, epsilon, sigma, temp_cold, ke_Energy_cold, n_hot_atoms);
    }
    
    forces = create_2D_Array(natoms,dim);

    double start = omp_get_wtime();

    int step = 0;
    //writeOutput(outfile, step, natoms, dim, box, pos, vels, group);
	
	void compute_Force(double**&, double**&, double**&, double*&, double**&, double, double, double, int, double &, double &, double, double);
	compute_Force(box, pos, vels, mass, forces, rcut_off, epsilon, sigma, natoms, pot_Energy, pressure_virial, energy_correction, energy_at_cutoff);	

	printf("Step    T_HOT    T_COLD    T_SYS     PE    KE    P \n");
	while (step <= nsteps){
		//Compute output (temperature)
	    temp_system = compute_Temp(vels, mass, dim, natoms, group, ke_Energy_hot, ke_Energy_cold, temp_hot, temp_cold );
	    	
		printf("%5d    %2.6lf    %2.6lf    %2.6lf    %2.6lf    %2.6lf    %lf\n", step, temp_hot, temp_cold, temp_system, pot_Energy/natoms, (ke_Energy_hot+ke_Energy_cold)/natoms, temp_system*Boltzmann*rho + pressure_virial/volume/dim);
		
	    if(step % output_freq == 0) writeOutput(outfile_Dump, step, natoms, dim, box, pos, vels, group);

		//Move the system in time

		integrate(natoms, pos, vels, mass, forces, group, relax, dt, rcut_off, epsilon, sigma, pot_Energy, box, gfric, noise_hot, noise_cold, pressure_virial, energy_correction, energy_at_cutoff);
		
		step += 1;
		
	}

	double end =  omp_get_wtime();

	cout<< "\n\nWork took "<< (end - start) << "\n";

writeOutput(outfile, step, natoms, dim, box, pos, vels, group);
fclose(outfile);
fclose(outfile_Dump);
return 0;
}
