#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <array>
#include <sstream>
#include "dcd.h"
#include "gaussianNoise.hpp"
#include "forces.hpp"

using namespace std;
std::string line;
ifstream infile;

void get_configuration(double *x_coor, double *y_coor, double *z_coor, int bead_number){
    infile.open("coor_column.dat");// file containing numbers in 3 columns 
    int ctr = 0;
    while (std::getline(infile, line))
        {
            std::stringstream ss(line);
            double a, b, c;
            ss >> a;
            ss >> b;
            ss >> c;
            x_coor[ctr] = a;
            y_coor[ctr] = b;
            z_coor[ctr] = c;
            ctr+=1;
        }
}

int main(void){
    long long int seed = 10000;
    srand(seed);
    int bead_number = 300;
    double kbT = 4.1418;   // pN*nm
    double epsilon_lj = 0.10*4.1418; // pN*nm
    double sigma_lj = 27.0; //nm
    double rcut = 3.0*sigma_lj;
    double separation = 20.0; // nm initial separation between particles
    double lk = 2.0*sigma_lj;
    double d0=1.01*sigma_lj; // nm bond length
    double harmonic_constant = 3000.0*epsilon_lj/(pow(d0, 2)); // spring constant in the harmonic force 
    double dt=0.17;   // [dt] = fs
    double gamma=0.002924; // 
    double angle_constant=33.89454997; // 
    double m=1.0;  // Daltons, mass of 10 base pairs that are approx in the 3.18 nm of the DNA
    double theta0=3.14159265358979323846;
    double angleenergy=0;
    double box_size[3];
    long long int N_steps=500; // number of steps in the simulation
    int factor = 10;
    void *dcd = open_dcd_write ( "traj.dcd", "dcd", 1, dt, bead_number);
    double rer=0;
    double relax = 33.2*dt; 
    double k1 = 1.0/(relax*0.01);
    double k2 = (1.0/50)*k1;
    double p_um = 0;
    double p_mu = 0;
    double lambda = k1 + k2;
    p_um = 1.0 - (k2+k1*exp(-lambda*relax))/(lambda);
    p_mu = 1.0 - (k1+k2*exp(-lambda*relax))/(lambda);
    double p_um_1D;
    double p_mu_1D;
    p_um_1D = 1.0 - (k2+k1*exp(-lambda*dt))/(lambda);
    p_mu_1D = 1.0 - (k1+k2*exp(-lambda*dt))/(lambda);
    /// cooperative rates
    double alpha = 2.00*k2/k1;
    double k1c = alpha*k1;
    double k2c = k2;
    double lambdac = k1c+k2c;
    double p_um_cooperative = 1.0 - (k2c+k1c*exp(-lambdac*relax))/(lambdac);
    double p_mu_cooperative = 1.0 - (k1c+k2c*exp(-lambdac*relax))/(lambdac);
    double p_um_1D_cooperative = 1.0 - (k2c+k1c*exp(-lambdac*dt))/(lambdac);
    double p_mu_1D_cooperative = 1.0 - (k1c+k2c*exp(-lambdac*dt))/(lambdac);
    double dist;
    double cut = 1.122; 
    int epigenetic_identity[300] = {3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,   1 ,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3};
    int tmp_epigenetic_identity[300] = {3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,   1 ,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3};
    // create coordinate arrays
    double x_coor[bead_number];
    double y_coor[bead_number];
    double z_coor[bead_number];
    double vx[bead_number];
    double vy[bead_number];
    double vz[bead_number];
    
    // create list of directional distances
    double bond_x[bead_number-1];
    double bond_y[bead_number-1];
    double bond_z[bead_number-1];
    double correl_angle[bead_number-1];
    
    // create list of bond vectors
    double bond_r[bead_number-1];
    
    // create harmonic terms
    double harter_x[bead_number];
    double harter_y[bead_number];
    double harter_z[bead_number];
    // create and initialize lennard-jones terms
    double vdW_x[bead_number];
    double vdW_y[bead_number];
    double vdW_z[bead_number];
    // create and initialize angle terms
    double angle_x[bead_number];
    double angle_y[bead_number];
    double angle_z[bead_number];
    // create new coordinate arrays, for time evolution
    double new_x_coor[bead_number];
    double new_y_coor[bead_number];
    double new_z_coor[bead_number];
    double new_vx[bead_number];
    double new_vy[bead_number];
    double new_vz[bead_number];
    //double ranr_x=0, ranr_y=0, ranr_z=0; 
    double ranr_x[bead_number];
    double ranr_y[bead_number];
    double ranr_z[bead_number];
    int DNAReplicate = 0; 
    int DNAReplicateTime = 40000; 
    double rerDNA; 
    string epi_str = " ";
    for (int i=0; i<bead_number; ++i){
        new_x_coor[i]=0;
        new_y_coor[i]=0;
        new_z_coor[i]=0;
        new_vx[i]=0;
        new_vy[i]=0;
        new_vz[i]=0;
        angle_x[i]=0;
        angle_y[i]=0;
        angle_z[i]=0;
        vx[i]=generateGaussianNoise(0.0, 1.0, 1)*sqrt(kbT*0.6022/m);
        vy[i]=generateGaussianNoise(0.0, 1.0, 1)*sqrt(kbT*0.6022/m);
        vz[i]=generateGaussianNoise(0.0, 1.0, 1)*sqrt(kbT*0.6022/m);
    }
    
    for (int i=0; i<bead_number-1; ++i){
        correl_angle[i]=0;
    }

    get_configuration(x_coor, y_coor, z_coor, bead_number);
    
    ofstream myfile3;
    stringstream ss;
    ss << "n" << bead_number << "_end_to_end.dat";
    myfile3.open (ss.str());
    
    ofstream myfile6;
    myfile6.open ("epigenetic_identity.dat");

    double v1x, v1y, v1z, modv1, v2x, v2y, v2z, modv2, v1dotv2, theta;
    for (int k=1; k<N_steps; ++k){
          for(int p=0; p<bead_number-1; ++p){
              bond_x[p]=x_coor[p+1]-x_coor[p];
              bond_y[p]=y_coor[p+1]-y_coor[p];
              bond_z[p]=z_coor[p+1]-z_coor[p];
          }
 
          /////////////////////
          // DNA REPLICATION //
          /////////////////////
          if (DNAReplicate){
              if (k!=0 and k%DNAReplicateTime == 0){
                  cout << "replicating DNA at step: " << k << endl;
                  for(int p=0; p<bead_number; ++p){
                      rerDNA =  ((double)rand()/(double)RAND_MAX);     
                      if (rerDNA <= 0.5 and epigenetic_identity[p]!=1){
                          epigenetic_identity[p] =  3;
                          tmp_epigenetic_identity[p] =  3;
                      }
                  }
              }
          }
          
	  for(int p=0; p<bead_number-1; ++p){
              bond_x[p]=x_coor[p+1]-x_coor[p];
              bond_y[p]=y_coor[p+1]-y_coor[p];
              bond_z[p]=z_coor[p+1]-z_coor[p];
          }

          /////////////////////////////////
          // WRITE DATA                  //
          ///////////////////////////////// 
          if (k%factor==0) {
	      epi_str = "";
              for (int w=0; w<bead_number; ++w){
                  epi_str.append(std::to_string(epigenetic_identity[w]));
                  epi_str.append(1u,' ');
              }
              myfile6 << epi_str << endl; 
              myfile3 << sqrt(pow((x_coor[bead_number-1]-x_coor[0]), 2)+pow((y_coor[bead_number-1]-y_coor[0]), 2)+pow((z_coor[bead_number-1]-z_coor[0]), 2)) << endl;
              write_timestep (dcd, box_size, x_coor, y_coor, z_coor);
          }
    
          if ( k%1==0 ) {
              cout << "step: " << k << endl;
              myfile3 << sqrt(pow((x_coor[bead_number-1]-x_coor[0]), 2)+pow((y_coor[bead_number-1]-y_coor[0]), 2)+pow((z_coor[bead_number-1]-z_coor[0]), 2)) << endl;
              /////////////////////////////////
              // FORWARD TRANSITION          //
              ///////////////////////////////// 
              for ( int j=0; j<bead_number; ++j ){
                  if (epigenetic_identity[j] == 1){
                      for ( int m=0; m<bead_number; ++m){
                          dist = sqrt(pow((x_coor[j]-x_coor[m]), 2)+pow((y_coor[j]-y_coor[m]), 2)+pow((z_coor[j]-z_coor[m]), 2) );
                          if (dist < cut*sigma_lj and (  epigenetic_identity[m] == 3 )){
                              rer =  ((double)rand()/(double)RAND_MAX);
                              if (abs(j-m)>=2 and rer<p_um){
                                  tmp_epigenetic_identity[m]=2;
                              }
                              else if (abs(j-m) == 1 and rer < p_um_1D){ 
                                  tmp_epigenetic_identity[m]=2;
                              }
                          }
                      }
                  }
                  if (epigenetic_identity[j] == 2){
                      for ( int m=0; m<bead_number; ++m){
                          dist = sqrt(pow((x_coor[j]-x_coor[m]), 2)+pow((y_coor[j]-y_coor[m]), 2)+pow((z_coor[j]-z_coor[m]), 2) );
                          if (dist < cut*sigma_lj and (  epigenetic_identity[m] == 3 )){
                              rer =  ((double)rand()/(double)RAND_MAX);
                              if (abs(j-m)>=2 and rer<p_um_cooperative){
                                  tmp_epigenetic_identity[m]=2;
                              }
                              else if (abs(j-m) == 1 and rer < p_um_1D_cooperative){ 
                                  tmp_epigenetic_identity[m]=2;
                              }
                          }
                      }
                  }
              }
              for ( int j=0; j<bead_number; ++j ){
                   epigenetic_identity[j] = tmp_epigenetic_identity[j];
              }
              /////////////////////////////////
              // BACKWARD TRANSITION         //
              ///////////////////////////////// 
              for ( int j=0; j<bead_number; ++j ){
                  rer =  ((double)rand()/(double)RAND_MAX); 
                  if (rer <= p_mu and epigenetic_identity[j] != 1){
                      epigenetic_identity[j]=3; 
                      tmp_epigenetic_identity[j]=3; 
                  }
              } 
              for ( int j=0; j<bead_number; ++j ){
                  if (epigenetic_identity[j] == 3){
                      for ( int m=0; m<bead_number; ++m){
                          dist = sqrt(pow((x_coor[j]-x_coor[m]), 2)+pow((y_coor[j]-y_coor[m]), 2)+pow((z_coor[j]-z_coor[m]), 2) );
                          if (dist < cut*sigma_lj and (  epigenetic_identity[m] == 2 )){
                              rer =  ((double)rand()/(double)RAND_MAX);
                              if (abs(j-m)>=2 and rer<p_mu_cooperative){
                                  tmp_epigenetic_identity[m]=3;
                              }
                              else if (abs(j-m) == 1 and rer < p_mu_1D_cooperative){ 
                                  tmp_epigenetic_identity[m]=3;
                              }
                          }
                      }
                  }

              }
              for ( int j=0; j<bead_number; ++j ){
                   epigenetic_identity[j] = tmp_epigenetic_identity[j];
              }

          }

          /////////////////////////////////
          // INTEGRATE                   //
          ///////////////////////////////// 

          harmonic_force(x_coor, y_coor, z_coor, bond_x, bond_y, bond_z, bond_r, bead_number, harmonic_constant, harter_x, harter_y, harter_z, d0);
          kratky_porod(x_coor, y_coor, z_coor, bead_number, angle_x, angle_y, angle_z, kbT, lk, sigma_lj);
          lj_force(x_coor, y_coor, z_coor, vdW_x, vdW_y, vdW_z, bead_number, epsilon_lj, sigma_lj, rcut);
          
          for (int i=0; i<bead_number; ++i){
              ranr_x[i]=generateGaussianNoise(0.0, 1.0, 1);
              ranr_y[i]=generateGaussianNoise(0.0, 1.0, 1);
              ranr_z[i]=generateGaussianNoise(0.0, 1.0, 1);
          }

          for (int i=0; i<bead_number; ++i){
              new_vx[i] = vx[i] - gamma*dt/2.0*vx[i] + dt/(2.0*m)*0.6022*(harter_x[i]+angle_x[i]+vdW_x[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma*m)/dt)*ranr_x[i];
              new_x_coor[i]=x_coor[i]+new_vx[i]*dt;
              x_coor[i]=new_x_coor[i];
              vx[i]=new_vx[i];
              
              new_vy[i] = vy[i] - gamma*dt/2.0*vy[i] + dt/(2.0*m)*0.6022*(harter_y[i]+angle_y[i]+vdW_y[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma*m)/dt)*ranr_y[i];
              new_y_coor[i]=y_coor[i]+new_vy[i]*dt;
              y_coor[i]=new_y_coor[i];
              vy[i]=new_vy[i];
              
              new_vz[i] = vz[i] - gamma*dt/2.0*vz[i] + dt/(2.0*m)*0.6022*(harter_z[i]+angle_z[i]+vdW_z[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma*m)/dt)*ranr_z[i];
              new_z_coor[i]=z_coor[i]+new_vz[i]*dt;
              z_coor[i]=new_z_coor[i];
              vz[i]=new_vz[i];
          }

          harmonic_force(x_coor, y_coor, z_coor, bond_x, bond_y, bond_z, bond_r, bead_number, harmonic_constant, harter_x, harter_y, harter_z, d0);
          kratky_porod(x_coor, y_coor, z_coor, bead_number, angle_x, angle_y, angle_z, kbT, lk, sigma_lj);
          lj_force(x_coor, y_coor, z_coor, vdW_x, vdW_y, vdW_z, bead_number, epsilon_lj, sigma_lj, rcut);

          for (int i=0; i<bead_number; ++i){
              new_vx[i] = vx[i] - gamma*dt/2.0*vx[i] + dt/(2.0*m)*0.6022*(harter_x[i]+angle_x[i]+vdW_x[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma*m)/dt)*ranr_x[i];
              vx[i]=new_vx[i];          

              new_vy[i] = vy[i] - gamma*dt/2.0*vy[i] + dt/(2.0*m)*0.6022*(harter_y[i]+angle_y[i]+vdW_y[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma*m)/dt)*ranr_y[i];
              vy[i]=new_vy[i];          

              new_vz[i] = vz[i] - gamma*dt/2.0*vz[i] + dt/(2.0*m)*0.6022*(harter_z[i]+angle_z[i]+vdW_z[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma*m)/dt)*ranr_z[i];
              vz[i]=new_vz[i];          
          }

    }
myfile3.close();
myfile6.close();
return 0;
}
