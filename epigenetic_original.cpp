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


void harmonic_force(double *x_coor, double *y_coor, double *z_coor, double *bond_x, double *bond_y, double *bond_z, double *bond_r, int bead_number, double harmonic_constant, double *harter_x, double *harter_y, double *harter_z, double d0){
    for(int i=0; i<bead_number; ++i){
        harter_x[i]=0;
        harter_y[i]=0;
        harter_z[i]=0;
        }
    for(int i=0; i<bead_number-1; ++i){
        bond_x[i]=0;
        bond_y[i]=0;
        bond_z[i]=0;
        }
    //create bond distance and vectors
    for(int i=0; i<bead_number-1; ++i){
        bond_x[i]=x_coor[i+1]-x_coor[i];
        bond_y[i]=y_coor[i+1]-y_coor[i];
        bond_z[i]=z_coor[i+1]-z_coor[i];
        bond_r[i]=sqrt(pow(bond_x[i], 2)+pow(bond_y[i], 2)+pow(bond_z[i], 2));
    }
    // create list of basic harmonic terms and fill the array with zeros
    double bht_x[bead_number-1];
    double bht_y[bead_number-1];
    double bht_z[bead_number-1];
    for(int i=0; i<bead_number-1; ++i){
        bht_x[i]=0;
        bht_y[i]=0;
        bht_z[i]=0;
    }
    // calculate basic harmonic terms
    for(int i=0; i<bead_number-1; ++i){
        bht_x[i]=(-1)*(harmonic_constant)*bond_x[i]*(bond_r[i]-d0)/bond_r[i];
        bht_y[i]=(-1)*(harmonic_constant)*bond_y[i]*(bond_r[i]-d0)/bond_r[i];
        bht_z[i]=(-1)*(harmonic_constant)*bond_z[i]*(bond_r[i]-d0)/bond_r[i];
    }
    
    for(int i=0; i<bead_number; ++i){
        if (i==0) {
            harter_x[0]=(-1)*bht_x[0];
            harter_y[0]=(-1)*bht_y[0];
            harter_z[0]=(-1)*bht_z[0];
        }
        else if (i==bead_number-1) {
            harter_x[i]=bht_x[i-1];
            harter_y[i]=bht_y[i-1];
            harter_z[i]=bht_z[i-1];
        }
        else {
            harter_x[i]=bht_x[i-1]+(-1)*bht_x[i];
            harter_y[i]=bht_y[i-1]+(-1)*bht_y[i];
            harter_z[i]=bht_z[i-1]+(-1)*bht_z[i];
        }
    }
}

double harmonic_energy(double *x_coor, double *y_coor, double *z_coor, double *bond_x, double *bond_y, double *bond_z, double *bond_r, int bead_number, double harmonic_constant, double *harter_x, double *harter_y, double *harter_z, double d0){
    double total_harmonic=0;
    //create bond distance and vectors
    for(int i=0; i<bead_number-1; ++i){
        bond_x[i]=x_coor[i+1]-x_coor[i];
        bond_y[i]=y_coor[i+1]-y_coor[i];
        bond_z[i]=z_coor[i+1]-z_coor[i];
        bond_r[i]=sqrt(pow(bond_x[i], 2)+pow(bond_y[i], 2)+pow(bond_z[i], 2));
        total_harmonic+=0.5*harmonic_constant*pow((bond_r[i]-d0), 2);
    }
    // create list of basic harmonic terms and fill the array with zeros
return total_harmonic;    
}

void kratky_porod(double *x_coor, double *y_coor, double *z_coor, int bead_number, double *angle_x, double *angle_y, double *angle_z, double kbT, double lk, double sigma_lj){
    double v1x=0;
    double v1y=0;
    double v1z=0;
    double v2x=0;
    double v2y=0;
    double v2z=0;
    double modv1=0;
    double modv2=0;
    double v1dotv2=0;
 
    for(int i=0; i<bead_number; ++i){
        angle_x[i]=0;
        angle_y[i]=0;
        angle_z[i]=0;
    }

    for(int i=0; i<bead_number-2; ++i){
        v1x=x_coor[i]-x_coor[i+1];     
        v1y=y_coor[i]-y_coor[i+1];     
        v1z=z_coor[i]-z_coor[i+1];
        modv1=sqrt(pow(v1x, 2)+pow(v1y, 2)+pow(v1z, 2));
        v2x=x_coor[i+2]-x_coor[i+1];     
        v2y=y_coor[i+2]-y_coor[i+1];     
        v2z=z_coor[i+2]-z_coor[i+1];    
        modv2=sqrt(pow(v2x, 2)+pow(v2y, 2)+pow(v2z, 2));
        v1dotv2=v1x*v2x+v1y*v2y+v1z*v2z;

        angle_x[i]+=(-1)*kbT*lk/(2.0*sigma_lj) * ( v2x/(modv1*modv2)  - ((v1dotv2/(pow((modv1*modv2) , 2))) * (v1x/modv1)*modv2) ); 
        angle_y[i]+=(-1)*kbT*lk/(2.0*sigma_lj) * ( v2y/(modv1*modv2)  - ((v1dotv2/(pow((modv1*modv2) , 2))) * (v1y/modv1)*modv2) );
        angle_z[i]+=(-1)*kbT*lk/(2.0*sigma_lj) * ( v2z/(modv1*modv2)  - ((v1dotv2/(pow((modv1*modv2) , 2))) * (v1z/modv1)*modv2) );
    
        angle_x[i+1]+= (-1) * ( (-1)*kbT*lk/(2.0*sigma_lj) * ( v2x/(modv1*modv2)  - ((v1dotv2/(pow((modv1*modv2) , 2))) * (v1x/modv1)*modv2) ) + (-1)*kbT*lk/(2.0*sigma_lj)* ( v1x/(modv1*modv2) - ((v1dotv2/(pow((modv1*modv2) , 2))) * (v2x/modv2)*modv1) )      );
        angle_y[i+1]+= (-1) * ( (-1)*kbT*lk/(2.0*sigma_lj) * ( v2y/(modv1*modv2)  - ((v1dotv2/(pow((modv1*modv2) , 2))) * (v1y/modv1)*modv2) ) + (-1)*kbT*lk/(2.0*sigma_lj)* ( v1y/(modv1*modv2) - ((v1dotv2/(pow((modv1*modv2) , 2))) * (v2y/modv2)*modv1) )      );
        angle_z[i+1]+= (-1) * ( (-1)*kbT*lk/(2.0*sigma_lj) * ( v2z/(modv1*modv2)  - ((v1dotv2/(pow((modv1*modv2) , 2))) * (v1z/modv1)*modv2) ) + (-1)*kbT*lk/(2.0*sigma_lj)* ( v1z/(modv1*modv2) - ((v1dotv2/(pow((modv1*modv2) , 2))) * (v2z/modv2)*modv1) )   );
 
        angle_x[i+2]+=(-1)*kbT*lk/(2.0*sigma_lj)* ( v1x/(modv1*modv2) - ((v1dotv2/(pow((modv1*modv2) , 2))) * (v2x/modv2)*modv1) );  
        angle_y[i+2]+=(-1)*kbT*lk/(2.0*sigma_lj)* ( v1y/(modv1*modv2) - ((v1dotv2/(pow((modv1*modv2) , 2))) * (v2y/modv2)*modv1) );  
        angle_z[i+2]+=(-1)*kbT*lk/(2.0*sigma_lj)* ( v1z/(modv1*modv2) - ((v1dotv2/(pow((modv1*modv2) , 2))) * (v2z/modv2)*modv1) );  
    }
}

void angle(double *x_coor, double *y_coor, double *z_coor, int bead_number, double angle_constant, double *angle_x, double *angle_y, double *angle_z, double theta0){
    double v1x=0;
    double v1y=0;
    double v1z=0;
    double v2x=0;
    double v2y=0;
    double v2z=0;
    double modv1=0;
    double modv2=0;
    double v1dotv2=0;
    double theta=0;
    double v1crossv2_x=0;
    double v1crossv2_y=0;
    double v1crossv2_z=0;
    double v1crossv2_mod=0;
 
    for(int i=0; i<bead_number; ++i){
        angle_x[i]=0;
        angle_y[i]=0;
        angle_z[i]=0;
    }

    for(int i=0; i<bead_number-2; ++i){
        v1x=x_coor[i]-x_coor[i+1];     
        v1y=y_coor[i]-y_coor[i+1];     
        v1z=z_coor[i]-z_coor[i+1];
        modv1=sqrt(pow(v1x, 2)+pow(v1y, 2)+pow(v1z, 2));
        v2x=x_coor[i+2]-x_coor[i+1];     
        v2y=y_coor[i+2]-y_coor[i+1];     
        v2z=z_coor[i+2]-z_coor[i+1];    
        modv2=sqrt(pow(v2x, 2)+pow(v2y, 2)+pow(v2z, 2));
        v1dotv2=v1x*v2x+v1y*v2y+v1z*v2z;
        if (v1dotv2/(modv1*modv2) > 1.0 || v1dotv2/(modv1*modv2) < -1.0){
            //cout <<  "original: " << v1dotv2/(modv1*modv2) << endl;
            theta=acos(round(v1dotv2/(modv1*modv2)));      //theta is in radians, so should theta0 be
            //cout <<  "rounded: "<< round(v1dotv2/(modv1*modv2))  << endl;
        }
        else{
            theta=acos(v1dotv2/(modv1*modv2));      //theta is in radians, so should theta0 be
        }
        v1crossv2_x=v1y*v2z-v1z*v2y;
        v1crossv2_y=v1z*v2x-v1x*v2z;
        v1crossv2_z=v1x*v2y-v1y*v2x;
        v1crossv2_mod=sqrt(pow(v1crossv2_x, 2)+pow(v1crossv2_y, 2)+pow(v1crossv2_z, 2));

        angle_x[i]+=2.0*angle_constant*(theta-theta0)*1.0/v1crossv2_mod*(v2x-v1dotv2/pow(modv1, 2)*v1x)*v1x/modv1;
        angle_y[i]+=2.0*angle_constant*(theta-theta0)*1.0/v1crossv2_mod*(v2y-v1dotv2/pow(modv1, 2)*v1y)*v1y/modv1;
        angle_z[i]+=2.0*angle_constant*(theta-theta0)*1.0/v1crossv2_mod*(v2z-v1dotv2/pow(modv1, 2)*v1z)*v1z/modv1;
       
        angle_x[i+1]+=(-1)*( 2.0*angle_constant*(theta-theta0)*1.0/v1crossv2_mod*(v2x-v1dotv2/pow(modv1, 2)*v1x)*v1x/modv1 + 2.0*angle_constant*(theta-theta0)*1.0/v1crossv2_mod*(v1x-v1dotv2/pow(modv2, 2)*v2x)*v2x/modv2 );
        angle_y[i+1]+=(-1)*( 2.0*angle_constant*(theta-theta0)*1.0/v1crossv2_mod*(v2y-v1dotv2/pow(modv1, 2)*v1y)*v1y/modv1 + 2.0*angle_constant*(theta-theta0)*1.0/v1crossv2_mod*(v1y-v1dotv2/pow(modv2, 2)*v2y)*v2y/modv2  );
        angle_z[i+1]+=(-1)*( 2.0*angle_constant*(theta-theta0)*1.0/v1crossv2_mod*(v2z-v1dotv2/pow(modv1, 2)*v1z)*v1z/modv1 + 2.0*angle_constant*(theta-theta0)*1.0/v1crossv2_mod*(v1z-v1dotv2/pow(modv2, 2)*v2z)*v2z/modv2  );
   
        angle_x[i+2]+=2.0*angle_constant*(theta-theta0)*1.0/v1crossv2_mod*(v1x-v1dotv2/pow(modv2, 2)*v2x)*v2x/modv2;
        angle_y[i+2]+=2.0*angle_constant*(theta-theta0)*1.0/v1crossv2_mod*(v1y-v1dotv2/pow(modv2, 2)*v2y)*v2y/modv2;
        angle_z[i+2]+=2.0*angle_constant*(theta-theta0)*1.0/v1crossv2_mod*(v1z-v1dotv2/pow(modv2, 2)*v2z)*v2z/modv2;

        //cout << "theta: " << theta << " v1dotv2: " << v1dotv2 << " modv1: " << modv1 << " modv2: " << modv2 << " ratio: " << v1dotv2/(modv1*modv2)  << endl;
        //cout << "theta: " << theta << " v1dotv2: " << v1dotv2 << " modv1: " << modv1 << " modv2: " << modv2 << " ratio: " << v1dotv2/(modv1*modv2)  << endl;
    }
}

double angle_energy(double *x_coor, double *y_coor, double *z_coor, int bead_number, double angle_constant, double theta0, double angleenergy){
    double v1x=0;
    double v1y=0;
    double v1z=0;
    double v2x=0;
    double v2y=0;
    double v2z=0;
    double modv1=0;
    double modv2=0;
    double v1dotv2=0;
    double theta=0;
    angleenergy=0;

    for(int i=0; i<bead_number-2; ++i){
        v1x=x_coor[i]-x_coor[i+1];     
        v1y=y_coor[i]-y_coor[i+1];     
        v1z=z_coor[i]-z_coor[i+1];
        modv1=sqrt(pow(v1x, 2)+pow(v1y, 2)+pow(v1z, 2));
        v2x=x_coor[i+2]-x_coor[i+1];     
        v2y=y_coor[i+2]-y_coor[i+1];     
        v2z=z_coor[i+2]-z_coor[i+1];    
        modv2=sqrt(pow(v2x, 2)+pow(v2y, 2)+pow(v2z, 2));
        v1dotv2=v1x*v2x+v1y*v2y+v1z*v2z;
        if (v1dotv2/(modv1*modv2) > 1.0 || v1dotv2/(modv1*modv2) < -1.0){
            theta=acos(round(v1dotv2/(modv1*modv2)));      //theta is in radians, so should theta0 be
        }
        else{
            theta=acos(v1dotv2/(modv1*modv2));      //theta is in radians, so should theta0 be
        }
        angleenergy+=1.0/2.0*angle_constant*pow((theta-theta0), 2);
        //cout << "angle energy: " << 1.0/2.0*pow((theta-theta0), 2)  << endl;
    }
    return angleenergy;
}

void epigenetic_force(double *x_coor, double *y_coor, double *z_coor, double *vdW_x, double *vdW_y, double *vdW_z, int bead_number, double epsilon_lj, double sigma_lj, int *epigenetic_identity, double kbT, double active_e, double inactive_e){

double rvector=0;
    for (int i=0; i<bead_number; ++i){
        vdW_x[i]=0;
        vdW_y[i]=0;
        vdW_z[i]=0;
        for (int j=0; j<bead_number; j++){
            if (j != i-1 && j != i+1 && j != i){
                rvector=sqrt(pow((x_coor[i]-x_coor[j]), 2)+pow((y_coor[i]-y_coor[j]), 2)+pow((z_coor[i]-z_coor[j]), 2));
                    if (rvector>pow(2, 1.0/6.0)*sigma_lj and ( (epigenetic_identity[i] != epigenetic_identity[j]))) { 
                        vdW_x[i]+=0;
                        vdW_y[i]+=0;
                        vdW_z[i]+=0;
                    }
                    else if (rvector<=pow(2, 1.0/6.0)*sigma_lj and ((epigenetic_identity[i]!=epigenetic_identity[j]))){
                        vdW_x[i]+=24.0*inactive_e*(x_coor[i]-x_coor[j])*1/pow(rvector, 2)*(2*pow(sigma_lj, 12)*pow(rvector, -12)-pow(sigma_lj, 6)*pow(rvector, -6));
                        vdW_y[i]+=24.0*inactive_e*(y_coor[i]-y_coor[j])*1/pow(rvector, 2)*(2*pow(sigma_lj, 12)*pow(rvector, -12)-pow(sigma_lj, 6)*pow(rvector, -6));
                        vdW_z[i]+=24.0*inactive_e*(z_coor[i]-z_coor[j])*1/pow(rvector, 2)*(2*pow(sigma_lj, 12)*pow(rvector, -12)-pow(sigma_lj, 6)*pow(rvector, -6));
                    }
                    else if (rvector<=1.8*sigma_lj and ((epigenetic_identity[i]==3 and epigenetic_identity[j]==3)  )){
                        vdW_x[i]+=24.0/0.8858*active_e*(x_coor[i]-x_coor[j])*1/pow(rvector, 2)*(2*pow(sigma_lj, 12)*pow(rvector, -12)-pow(sigma_lj, 6)*pow(rvector, -6));
                        vdW_y[i]+=24.0/0.8858*active_e*(y_coor[i]-y_coor[j])*1/pow(rvector, 2)*(2*pow(sigma_lj, 12)*pow(rvector, -12)-pow(sigma_lj, 6)*pow(rvector, -6));
                        vdW_z[i]+=24.0/0.8858*active_e*(z_coor[i]-z_coor[j])*1/pow(rvector, 2)*(2*pow(sigma_lj, 12)*pow(rvector, -12)-pow(sigma_lj, 6)*pow(rvector, -6));
                    }
                    else if (rvector<=1.8*sigma_lj and ( (epigenetic_identity[i]==2 and epigenetic_identity[j]==2) )){
                        vdW_x[i]+=24.0/0.8858*inactive_e*(x_coor[i]-x_coor[j])*1/pow(rvector, 2)*(2*pow(sigma_lj, 12)*pow(rvector, -12)-pow(sigma_lj, 6)*pow(rvector, -6));
                        vdW_y[i]+=24.0/0.8858*inactive_e*(y_coor[i]-y_coor[j])*1/pow(rvector, 2)*(2*pow(sigma_lj, 12)*pow(rvector, -12)-pow(sigma_lj, 6)*pow(rvector, -6));
                        vdW_z[i]+=24.0/0.8858*inactive_e*(z_coor[i]-z_coor[j])*1/pow(rvector, 2)*(2*pow(sigma_lj, 12)*pow(rvector, -12)-pow(sigma_lj, 6)*pow(rvector, -6));
                    }
                    else if (rvector>1.8*sigma_lj and ((epigenetic_identity[i]==3 and epigenetic_identity[j]==3) or (epigenetic_identity[i]==2 and epigenetic_identity[j]==2) )){
                        vdW_x[i]+=0;
                        vdW_y[i]+=0;
                        vdW_z[i]+=0;
                    }
            }
        }
    }
}

void lj_force(double *x_coor, double *y_coor, double *z_coor, double *vdW_x, double *vdW_y, double *vdW_z, int bead_number, double epsilon_lj, double sigma_lj, double rcut){

double rvector=0;

    for (int i=0; i<bead_number; ++i){
        vdW_x[i]=0;
        vdW_y[i]=0;
        vdW_z[i]=0;
        for (int j=0; j<bead_number; j++){
            if (j != i-1 && j != i+1 && j != i){
                rvector=sqrt(pow((x_coor[i]-x_coor[j]), 2)+pow((y_coor[i]-y_coor[j]), 2)+pow((z_coor[i]-z_coor[j]), 2));
                    if (rvector>rcut) { 
                        vdW_x[i]+=0;
                        vdW_y[i]+=0;
                        vdW_z[i]+=0;
                    }
                    else{
                        vdW_x[i]+=24.0*epsilon_lj*(x_coor[i]-x_coor[j])*1/pow(rvector, 2)*(2*pow(sigma_lj, 12)*pow(rvector, -12)-pow(sigma_lj, 6)*pow(rvector, -6));
                        vdW_y[i]+=24.0*epsilon_lj*(y_coor[i]-y_coor[j])*1/pow(rvector, 2)*(2*pow(sigma_lj, 12)*pow(rvector, -12)-pow(sigma_lj, 6)*pow(rvector, -6));
                        vdW_z[i]+=24.0*epsilon_lj*(z_coor[i]-z_coor[j])*1/pow(rvector, 2)*(2*pow(sigma_lj, 12)*pow(rvector, -12)-pow(sigma_lj, 6)*pow(rvector, -6));
                    }
            }
        }
    }
}

double lj_energy(double *x_coor, double *y_coor, double *z_coor, double *vdW_x, double *vdW_y, double *vdW_z, int bead_number, double epsilon_lj, double sigma_lj){
    double total_lj=0;
    double rvector=0;
    for (int i=0; i<bead_number; ++i){
        for (int j=0; j<i; j++){
        //for (int j=0; j<bead_number; j++){
            if (j != i-1 && j != i+1 && j != i){
                //cout << "lj interaction between particle i: " << i << " and particle j:  " << j << endl; 
                rvector=sqrt(pow((x_coor[i]-x_coor[j]), 2)+pow((y_coor[i]-y_coor[j]), 2)+pow((z_coor[i]-z_coor[j]), 2));
                if (rvector>pow(2, 1.0/6.0)*sigma_lj) {
                    total_lj+=0.0; 
                }
                else{
                    total_lj+=4.0*epsilon_lj*( pow((sigma_lj/rvector), 12) - pow((sigma_lj/rvector), 6)  );
                }
            }
        }
    }
return total_lj;    
}

int main(void){
for (int g=300; g<301; g+=40){
long long int seed = 10000;
    srand(seed);
    int bead_number=g; //number of beads
    double kbT = 4.1418;   // pN*nm
    double epsilon_lj = 0.10*4.1418; // pN*nm
    double sigma_lj = 27.0; //nm
    double rcut = 3.0*sigma_lj;
    double separation = 20.0; // nm initial separation between particles
    double lk = 2.0*sigma_lj;
    double d0=1.01*sigma_lj; // nm bond length
    double harmonic_constant = 3000.0*epsilon_lj/(pow(d0, 2)); // spring constant in the harmonic force // [harmonic_constant]= kcal/mol/Angstrom^2
    double dt=0.17;   // [dt] = fs
    double gamma=0.002924; // [epsilon] = N*s/m, here it is expressed in mol*Angstrom^2/kcal/fs
    double angle_constant=33.89454997; // pN*nm
    double m=1.0;  // Daltons, mass of 10 base pairs that are approx in the 3.18 nm of the DNA
    double theta0=3.14159265358979323846;
    double angleenergy=0;
    double box_size[3];
    long long int N_steps=100000000; // number of steps in the simulation
    int factor = 10000;
    void *dcd = open_dcd_write ( "traj.dcd", "dcd", 1, dt, bead_number);
    int epigenetic_identity[300] = {3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,   1 ,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3};
   int tmp_epigenetic_identity[300] = {3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,   1 ,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3};
    double rer=0;
    double relax = 33.2*dt; //formally, in timestep units, the relaxation time for lp=1sigma is 20700 timesteps!
    double k1 = 1.0/(relax*100.0);
    double k2 = (1.0/10000.0)*k1;
    double p_um = 0;
    double p_mu = 0;
    double lambda = k1 + k2;
    p_um = 1.0 - (k2+k1*exp(-lambda*relax))/(lambda);
    p_mu = 1.0 - (k1+k2*exp(-lambda*relax))/(lambda);
    cout << "k1= "<< k1 << endl;
    cout << "k2= "<< k2 << endl;
    cout << "p_um= "<< p_um << endl;
    cout << "p_mu= "<< p_mu << endl;
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
    box_size[0]=sigma_lj+sigma_lj/2.0;
    box_size[1]=sigma_lj+sigma_lj/2.0;
    box_size[2]=sigma_lj+sigma_lj/2.0;
    int DNAReplicate = 0; 
    int DNAReplicateTime = 40000; 
    double rerDNA; 
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

    //create_configuration(x_coor, y_coor, z_coor, bead_number, separation);
    get_configuration(x_coor, y_coor, z_coor, bead_number);
    //for (int i=0; i<bead_number; ++i){
    //    cout << x_coor[i] << " " << y_coor[i] << " " << z_coor[i] << endl;  
    //}
 
    //ofstream myfile1;
    //stringstream aa;
    //aa << "n" << bead_number << "_angle_energy.dat";
    //myfile1.open (aa.str());
    
    ofstream myfile3;
    stringstream ss;
    ss << "n" << bead_number << "_end_to_end.dat";
    myfile3.open (ss.str());
    
    //ofstream myfile4;
    //stringstream cc;
    //cc << "traj.xyz";
    //myfile4.open (cc.str());
    
    //ofstream myfile5;
    //myfile5.open ("positions.dat");
    
    ofstream myfile6;
    myfile6.open ("epigenetic_identity.dat");

    double v1x, v1y, v1z, modv1, v2x, v2y, v2z, modv2, v1dotv2, theta;
    for (int k=1; k<N_steps; ++k){
        //cout << "step: " << k << endl;
        for(int s=0; s<bead_number; ++s) {
            //cout << " coordinates particle: " << s << " " << x_coor[s] << " " << y_coor[s] << " " << z_coor[s] << endl;
        }
    //for (int k=1; k<5; ++k){
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

          if (k%factor==0) {
              myfile6 << epigenetic_identity[0] << " "     << epigenetic_identity[1] << " "     << epigenetic_identity[2] << " "     << epigenetic_identity[3] << " "     << epigenetic_identity[4] << " "     << epigenetic_identity[5] << " "     << epigenetic_identity[6] << " "     << epigenetic_identity[7] << " "     << epigenetic_identity[8] << " "     << epigenetic_identity[9] << " "     << epigenetic_identity[10] << " "     << epigenetic_identity[11] << " "     << epigenetic_identity[12] << " "     << epigenetic_identity[13] << " "     << epigenetic_identity[14] << " "     << epigenetic_identity[15] << " "     << epigenetic_identity[16] << " "     << epigenetic_identity[17] << " "     << epigenetic_identity[18] << " "     << epigenetic_identity[19] << " "     << epigenetic_identity[20] << " "     << epigenetic_identity[21] << " "     << epigenetic_identity[22] << " "     << epigenetic_identity[23] << " "     << epigenetic_identity[24] << " "     << epigenetic_identity[25] << " "     << epigenetic_identity[26] << " "     << epigenetic_identity[27] << " "     << epigenetic_identity[28] << " "     << epigenetic_identity[29] << " "     << epigenetic_identity[30] << " "     << epigenetic_identity[31] << " "     << epigenetic_identity[32] << " "     << epigenetic_identity[33] << " "     << epigenetic_identity[34] << " "     << epigenetic_identity[35] << " "     << epigenetic_identity[36] << " "     << epigenetic_identity[37] << " "     << epigenetic_identity[38] << " "     << epigenetic_identity[39] << " "     << epigenetic_identity[40] << " "     << epigenetic_identity[41] << " "     << epigenetic_identity[42] << " "     << epigenetic_identity[43] << " "     << epigenetic_identity[44] << " "     << epigenetic_identity[45] << " "     << epigenetic_identity[46] << " "     << epigenetic_identity[47] << " "     << epigenetic_identity[48] << " "     << epigenetic_identity[49] << " "     << epigenetic_identity[50] << " "     << epigenetic_identity[51] << " "     << epigenetic_identity[52] << " "     << epigenetic_identity[53] << " "     << epigenetic_identity[54] << " "     << epigenetic_identity[55] << " "     << epigenetic_identity[56] << " "     << epigenetic_identity[57] << " "     << epigenetic_identity[58] << " "     << epigenetic_identity[59] << " "     << epigenetic_identity[60] << " "     << epigenetic_identity[61] << " "     << epigenetic_identity[62] << " "     << epigenetic_identity[63] << " "     << epigenetic_identity[64] << " "     << epigenetic_identity[65] << " "     << epigenetic_identity[66] << " "     << epigenetic_identity[67] << " "     << epigenetic_identity[68] << " "     << epigenetic_identity[69] << " "     << epigenetic_identity[70] << " "     << epigenetic_identity[71] << " "     << epigenetic_identity[72] << " "     << epigenetic_identity[73] << " "     << epigenetic_identity[74] << " "     << epigenetic_identity[75] << " "     << epigenetic_identity[76] << " "     << epigenetic_identity[77] << " "     << epigenetic_identity[78] << " "     << epigenetic_identity[79] << " "     << epigenetic_identity[80] << " "     << epigenetic_identity[81] << " "     << epigenetic_identity[82] << " "     << epigenetic_identity[83] << " "     << epigenetic_identity[84] << " "     << epigenetic_identity[85] << " "     << epigenetic_identity[86] << " "     << epigenetic_identity[87] << " "     << epigenetic_identity[88] << " "     << epigenetic_identity[89] << " "     << epigenetic_identity[90] << " "     << epigenetic_identity[91] << " "     << epigenetic_identity[92] << " "     << epigenetic_identity[93] << " "     << epigenetic_identity[94] << " "     << epigenetic_identity[95] << " "     << epigenetic_identity[96] << " "     << epigenetic_identity[97] << " "     << epigenetic_identity[98] << " "     << epigenetic_identity[99] << " "     << epigenetic_identity[100] << " "     << epigenetic_identity[101] << " "     << epigenetic_identity[102] << " "     << epigenetic_identity[103] << " "     << epigenetic_identity[104] << " "     << epigenetic_identity[105] << " "     << epigenetic_identity[106] << " "     << epigenetic_identity[107] << " "     << epigenetic_identity[108] << " "     << epigenetic_identity[109] << " "     << epigenetic_identity[110] << " "     << epigenetic_identity[111] << " "     << epigenetic_identity[112] << " "     << epigenetic_identity[113] << " "     << epigenetic_identity[114] << " "     << epigenetic_identity[115] << " "     << epigenetic_identity[116] << " "     << epigenetic_identity[117] << " "     << epigenetic_identity[118] << " "     << epigenetic_identity[119] << " "     << epigenetic_identity[120] << " "     << epigenetic_identity[121] << " "     << epigenetic_identity[122] << " "     << epigenetic_identity[123] << " "     << epigenetic_identity[124] << " "     << epigenetic_identity[125] << " "     << epigenetic_identity[126] << " "     << epigenetic_identity[127] << " "     << epigenetic_identity[128] << " "     << epigenetic_identity[129] << " "     << epigenetic_identity[130] << " "     << epigenetic_identity[131] << " "     << epigenetic_identity[132] << " "     << epigenetic_identity[133] << " "     << epigenetic_identity[134] << " "     << epigenetic_identity[135] << " "     << epigenetic_identity[136] << " "     << epigenetic_identity[137] << " "     << epigenetic_identity[138] << " "     << epigenetic_identity[139] << " "     << epigenetic_identity[140] << " "     << epigenetic_identity[141] << " "     << epigenetic_identity[142] << " "     << epigenetic_identity[143] << " "     << epigenetic_identity[144] << " "     << epigenetic_identity[145] << " "     << epigenetic_identity[146] << " "     << epigenetic_identity[147] << " "     << epigenetic_identity[148] << " "     << epigenetic_identity[149] << " "     << epigenetic_identity[150] << " "     << epigenetic_identity[151] << " "     << epigenetic_identity[152] << " "     << epigenetic_identity[153] << " "     << epigenetic_identity[154] << " "     << epigenetic_identity[155] << " "     << epigenetic_identity[156] << " "     << epigenetic_identity[157] << " "     << epigenetic_identity[158] << " "     << epigenetic_identity[159] << " "     << epigenetic_identity[160] << " "     << epigenetic_identity[161] << " "     << epigenetic_identity[162] << " "     << epigenetic_identity[163] << " "     << epigenetic_identity[164] << " "     << epigenetic_identity[165] << " "     << epigenetic_identity[166] << " "     << epigenetic_identity[167] << " "     << epigenetic_identity[168] << " "     << epigenetic_identity[169] << " "     << epigenetic_identity[170] << " "     << epigenetic_identity[171] << " "     << epigenetic_identity[172] << " "     << epigenetic_identity[173] << " "     << epigenetic_identity[174] << " "     << epigenetic_identity[175] << " "     << epigenetic_identity[176] << " "     << epigenetic_identity[177] << " "     << epigenetic_identity[178] << " "     << epigenetic_identity[179] << " "     << epigenetic_identity[180] << " "     << epigenetic_identity[181] << " "     << epigenetic_identity[182] << " "     << epigenetic_identity[183] << " "     << epigenetic_identity[184] << " "     << epigenetic_identity[185] << " "     << epigenetic_identity[186] << " "     << epigenetic_identity[187] << " "     << epigenetic_identity[188] << " "     << epigenetic_identity[189] << " "     << epigenetic_identity[190] << " "     << epigenetic_identity[191] << " "     << epigenetic_identity[192] << " "     << epigenetic_identity[193] << " "     << epigenetic_identity[194] << " "     << epigenetic_identity[195] << " "     << epigenetic_identity[196] << " "     << epigenetic_identity[197] << " "     << epigenetic_identity[198] << " "     << epigenetic_identity[199] << " "     << epigenetic_identity[200] << " "     << epigenetic_identity[201] << " "     << epigenetic_identity[202] << " "     << epigenetic_identity[203] << " "     << epigenetic_identity[204] << " "     << epigenetic_identity[205] << " "     << epigenetic_identity[206] << " "     << epigenetic_identity[207] << " "     << epigenetic_identity[208] << " "     << epigenetic_identity[209] << " "     << epigenetic_identity[210] << " "     << epigenetic_identity[211] << " "     << epigenetic_identity[212] << " "     << epigenetic_identity[213] << " "     << epigenetic_identity[214] << " "     << epigenetic_identity[215] << " "     << epigenetic_identity[216] << " "     << epigenetic_identity[217] << " "     << epigenetic_identity[218] << " "     << epigenetic_identity[219] << " "     << epigenetic_identity[220] << " "     << epigenetic_identity[221] << " "     << epigenetic_identity[222] << " "     << epigenetic_identity[223] << " "     << epigenetic_identity[224] << " "     << epigenetic_identity[225] << " "     << epigenetic_identity[226] << " "     << epigenetic_identity[227] << " "     << epigenetic_identity[228] << " "     << epigenetic_identity[229] << " "     << epigenetic_identity[230] << " "     << epigenetic_identity[231] << " "     << epigenetic_identity[232] << " "     << epigenetic_identity[233] << " "     << epigenetic_identity[234] << " "     << epigenetic_identity[235] << " "     << epigenetic_identity[236] << " "     << epigenetic_identity[237] << " "     << epigenetic_identity[238] << " "     << epigenetic_identity[239] << " "     << epigenetic_identity[240] << " "     << epigenetic_identity[241] << " "     << epigenetic_identity[242] << " "     << epigenetic_identity[243] << " "     << epigenetic_identity[244] << " "     << epigenetic_identity[245] << " "     << epigenetic_identity[246] << " "     << epigenetic_identity[247] << " "     << epigenetic_identity[248] << " "     << epigenetic_identity[249] << " "     << epigenetic_identity[250] << " "     << epigenetic_identity[251] << " "     << epigenetic_identity[252] << " "     << epigenetic_identity[253] << " "     << epigenetic_identity[254] << " "     << epigenetic_identity[255] << " "     << epigenetic_identity[256] << " "     << epigenetic_identity[257] << " "     << epigenetic_identity[258] << " "     << epigenetic_identity[259] << " "     << epigenetic_identity[260] << " "     << epigenetic_identity[261] << " "     << epigenetic_identity[262] << " "     << epigenetic_identity[263] << " "     << epigenetic_identity[264] << " "     << epigenetic_identity[265] << " "     << epigenetic_identity[266] << " "     << epigenetic_identity[267] << " "     << epigenetic_identity[268] << " "     << epigenetic_identity[269] << " "     << epigenetic_identity[270] << " "     << epigenetic_identity[271] << " "     << epigenetic_identity[272] << " "     << epigenetic_identity[273] << " "     << epigenetic_identity[274] << " "     << epigenetic_identity[275] << " "     << epigenetic_identity[276] << " "     << epigenetic_identity[277] << " "     << epigenetic_identity[278] << " "     << epigenetic_identity[279] << " "     << epigenetic_identity[280] << " "     << epigenetic_identity[281] << " "     << epigenetic_identity[282] << " "     << epigenetic_identity[283] << " "     << epigenetic_identity[284] << " "     << epigenetic_identity[285] << " "     << epigenetic_identity[286] << " "     << epigenetic_identity[287] << " "     << epigenetic_identity[288] << " "     << epigenetic_identity[289] << " "     << epigenetic_identity[290] << " "     << epigenetic_identity[291] << " "     << epigenetic_identity[292] << " "     << epigenetic_identity[293] << " "     << epigenetic_identity[294] << " "     << epigenetic_identity[295] << " "     << epigenetic_identity[296] << " "     << epigenetic_identity[297] << " "     << epigenetic_identity[298] << " "     << epigenetic_identity[299] << endl;

              myfile3 << sqrt(pow((x_coor[bead_number-1]-x_coor[0]), 2)+pow((y_coor[bead_number-1]-y_coor[0]), 2)+pow((z_coor[bead_number-1]-z_coor[0]), 2)) << endl;
              write_timestep (dcd, box_size, x_coor, y_coor, z_coor);
              //myfile4 << bead_number << endl;
              //myfile4 << "frame " << k/1000<< endl;
              //for (int j=0; j<bead_number; ++j){
              //    myfile4 << "   " << j+1 << "   " << x_coor[j] << "    " << y_coor[j] << "    " << z_coor[j] << endl;
              //}
          }
    
          if ( k%1==0 ) {
              cout << "step: " << k << endl;
              myfile3 << sqrt(pow((x_coor[bead_number-1]-x_coor[0]), 2)+pow((y_coor[bead_number-1]-y_coor[0]), 2)+pow((z_coor[bead_number-1]-z_coor[0]), 2)) << endl;
              /////////////////////
              // 1D modification //
              ///////////////////// 
              //for ( int j=0; j<bead_number; ++j ){
              //    rer =  ((double)rand()/(double)RAND_MAX);
              //    if (epigenetic_identity[j] == 2  or epigenetic_identity[j] == 1){
              //        if ((j>=1 and j<=(bead_number-2)) and epigenetic_identity[j-1]==3  ){
              //            if (rer < p_um_1D){tmp_epigenetic_identity[j-1]=2;}     
              //        }
              //        if ((j>=1 and j<=(bead_number-2)) and epigenetic_identity[j+1]==3){                                                                          
              //            if (rer < p_um_1D){tmp_epigenetic_identity[j+1]=2;}     
              //        } 
              //        if (j==0 and epigenetic_identity[j+1]==3  ){                                          
              //            if (rer < p_um_1D){tmp_epigenetic_identity[j+1]=2;}     
              //        }
              //        if (j==(bead_number-1) and epigenetic_identity[j-1]==3  ){                                          
              //            if (rer < p_um_1D) {tmp_epigenetic_identity[j-1]=2;}     
              //        }
              //    }
              //}
              /////////////////////
              // 3D modification //
              ///////////////////// 
              //for ( int j=0; j<bead_number; ++j ){
              //    if (epigenetic_identity[j] == 2  or epigenetic_identity[j] == 1){
              //        for ( int m=0; m<bead_number; ++m){
              //            dist = sqrt(pow((x_coor[j]-x_coor[m]), 2)+pow((y_coor[j]-y_coor[m]), 2)+pow((z_coor[j]-z_coor[m]), 2) );
              //            if (dist < cut*sigma_lj and (  epigenetic_identity[m] == 3 )){
              //                rer =  ((double)rand()/(double)RAND_MAX);
              //                if (abs(j-m)>=2 and rer<p_um){
              //                //    cout << "this is 3D spreading between particles: " << j << " " << m << "ranr: " << rer << "p_um: " << p_um  << "dist: " << dist << endl;
              //                    tmp_epigenetic_identity[m]=2;
              //                }
              //                //cout << "random number value: " << rer << endl; 
              //                else if (abs(j-m) == 1 and rer < p_um_1D){ 
              //                    tmp_epigenetic_identity[m]=2;
              //                    //cout << "particle changing identity: " << m << endl;
              //                }
              //                //cout << "new epigenetic identity: " << epigenetic_identity[m] << endl;
              //            }
              //        }
              //    }
              //}
              /////////////////////////////////
              // 3D modification cooperative //
              ///////////////////////////////// 
              for ( int j=0; j<bead_number; ++j ){
                  if (epigenetic_identity[j] == 1){
                      for ( int m=0; m<bead_number; ++m){
                          dist = sqrt(pow((x_coor[j]-x_coor[m]), 2)+pow((y_coor[j]-y_coor[m]), 2)+pow((z_coor[j]-z_coor[m]), 2) );
                          if (dist < cut*sigma_lj and (  epigenetic_identity[m] == 3 )){
                              rer =  ((double)rand()/(double)RAND_MAX);
                              if (abs(j-m)>=2 and rer<p_um){
                              //    cout << "this is 3D spreading between particles: " << j << " " << m << "ranr: " << rer << "p_um: " << p_um  << "dist: " << dist << endl;
                                  tmp_epigenetic_identity[m]=2;
                              }
                              //cout << "random number value: " << rer << endl; 
                              else if (abs(j-m) == 1 and rer < p_um_1D){ 
                                  tmp_epigenetic_identity[m]=2;
                                  //cout << "particle changing identity: " << m << endl;
                              }
                              //cout << "new epigenetic identity: " << epigenetic_identity[m] << endl;
                          }
                      }
                  }
                  if (epigenetic_identity[j] == 2){
                      for ( int m=0; m<bead_number; ++m){
                          dist = sqrt(pow((x_coor[j]-x_coor[m]), 2)+pow((y_coor[j]-y_coor[m]), 2)+pow((z_coor[j]-z_coor[m]), 2) );
                          if (dist < cut*sigma_lj and (  epigenetic_identity[m] == 3 )){
                              rer =  ((double)rand()/(double)RAND_MAX);
                              if (abs(j-m)>=2 and rer<p_um_cooperative){
                              //    cout << "this is 3D spreading between particles: " << j << " " << m << "ranr: " << rer << "p_um: " << p_um  << "dist: " << dist << endl;
                                  tmp_epigenetic_identity[m]=2;
                              }
                              //cout << "random number value: " << rer << endl; 
                              else if (abs(j-m) == 1 and rer < p_um_1D_cooperative){ 
                                  tmp_epigenetic_identity[m]=2;
                                  //cout << "particle changing identity: " << m << endl;
                              }
                              //cout << "new epigenetic identity: " << epigenetic_identity[m] << endl;
                          }
                      }
                  }
              }
              for ( int j=0; j<bead_number; ++j ){
                   epigenetic_identity[j] = tmp_epigenetic_identity[j];
              }
              //attempt to unmodify particle
              for ( int j=0; j<bead_number; ++j ){
                  rer =  ((double)rand()/(double)RAND_MAX); 
                  if (rer <= p_mu and epigenetic_identity[j] != 1){
                      //cout << "particle has been unsilenced: " << j << endl;
                      epigenetic_identity[j]=3; 
                      tmp_epigenetic_identity[j]=3; 
                  }
              } 
              /////////////////////////////////////
              // unmodify particle cooperatively //              
              /////////////////////////////////////
              for ( int j=0; j<bead_number; ++j ){
                  if (epigenetic_identity[j] == 3){
                      for ( int m=0; m<bead_number; ++m){
                          dist = sqrt(pow((x_coor[j]-x_coor[m]), 2)+pow((y_coor[j]-y_coor[m]), 2)+pow((z_coor[j]-z_coor[m]), 2) );
                          if (dist < cut*sigma_lj and (  epigenetic_identity[m] == 2 )){
                              rer =  ((double)rand()/(double)RAND_MAX);
                              if (abs(j-m)>=2 and rer<p_mu_cooperative){
                              //    cout << "this is 3D spreading between particles: " << j << " " << m << "ranr: " << rer << "p_um: " << p_um  << "dist: " << dist << endl;
                                  tmp_epigenetic_identity[m]=3;
                              }
                              //cout << "random number value: " << rer << endl; 
                              else if (abs(j-m) == 1 and rer < p_mu_1D_cooperative){ 
                                  tmp_epigenetic_identity[m]=3;
                                  //cout << "particle changing identity: " << m << endl;
                              }
                              //cout << "new epigenetic identity: " << epigenetic_identity[m] << endl;
                          }
                      }
                  }

              }
              for ( int j=0; j<bead_number; ++j ){
                   epigenetic_identity[j] = tmp_epigenetic_identity[j];
              }

          }

          harmonic_force(x_coor, y_coor, z_coor, bond_x, bond_y, bond_z, bond_r, bead_number, harmonic_constant, harter_x, harter_y, harter_z, d0);
          kratky_porod(x_coor, y_coor, z_coor, bead_number, angle_x, angle_y, angle_z, kbT, lk, sigma_lj);
          //angle(x_coor, y_coor, z_coor, bead_number, angle_constant, angle_x, angle_y, angle_z, theta0);
          lj_force(x_coor, y_coor, z_coor, vdW_x, vdW_y, vdW_z, bead_number, epsilon_lj, sigma_lj, rcut);
          //epigenetic_force(x_coor, y_coor, z_coor, vdW_x, vdW_y, vdW_z, bead_number, epsilon_lj, sigma_lj, epigenetic_identity, kbT, active_e, inactive_e);
          
          for (int i=0; i<bead_number; ++i){
              ranr_x[i]=generateGaussianNoise(0.0, 1.0, 1);
              ranr_y[i]=generateGaussianNoise(0.0, 1.0, 1);
              ranr_z[i]=generateGaussianNoise(0.0, 1.0, 1);
          }

          for (int i=0; i<bead_number; ++i){
              //myfile5 << fixed << setprecision(7) << k*dt << "   "  << x_coor[0] << "\n";
              //ranr_x=generateGaussianNoise(0.0, 1.0, 1);
              new_vx[i] = vx[i] - gamma*dt/2.0*vx[i] + dt/(2.0*m)*0.6022*(harter_x[i]+angle_x[i]+vdW_x[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma*m)/dt)*ranr_x[i];
              new_x_coor[i]=x_coor[i]+new_vx[i]*dt;
              x_coor[i]=new_x_coor[i];
              vx[i]=new_vx[i];
              
              //ranr_y=generateGaussianNoise(0.0, 1.0, 1);
              new_vy[i] = vy[i] - gamma*dt/2.0*vy[i] + dt/(2.0*m)*0.6022*(harter_y[i]+angle_y[i]+vdW_y[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma*m)/dt)*ranr_y[i];
              new_y_coor[i]=y_coor[i]+new_vy[i]*dt;
              y_coor[i]=new_y_coor[i];
              vy[i]=new_vy[i];
              
              //ranr_z=generateGaussianNoise(0.0, 1.0, 1);
              new_vz[i] = vz[i] - gamma*dt/2.0*vz[i] + dt/(2.0*m)*0.6022*(harter_z[i]+angle_z[i]+vdW_z[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma*m)/dt)*ranr_z[i];
              new_z_coor[i]=z_coor[i]+new_vz[i]*dt;
              z_coor[i]=new_z_coor[i];
              vz[i]=new_vz[i];
          }

          harmonic_force(x_coor, y_coor, z_coor, bond_x, bond_y, bond_z, bond_r, bead_number, harmonic_constant, harter_x, harter_y, harter_z, d0);
          kratky_porod(x_coor, y_coor, z_coor, bead_number, angle_x, angle_y, angle_z, kbT, lk, sigma_lj);
          //angle(x_coor, y_coor, z_coor, bead_number, angle_constant, angle_x, angle_y, angle_z, theta0);
          lj_force(x_coor, y_coor, z_coor, vdW_x, vdW_y, vdW_z, bead_number, epsilon_lj, sigma_lj, rcut);
          //epigenetic_force(x_coor, y_coor, z_coor, vdW_x, vdW_y, vdW_z, bead_number, epsilon_lj, sigma_lj, epigenetic_identity, kbT, active_e, inactive_e);

          for (int i=0; i<bead_number; ++i){
              new_vx[i] = vx[i] - gamma*dt/2.0*vx[i] + dt/(2.0*m)*0.6022*(harter_x[i]+angle_x[i]+vdW_x[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma*m)/dt)*ranr_x[i];
              vx[i]=new_vx[i];          

              new_vy[i] = vy[i] - gamma*dt/2.0*vy[i] + dt/(2.0*m)*0.6022*(harter_y[i]+angle_y[i]+vdW_y[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma*m)/dt)*ranr_y[i];
              vy[i]=new_vy[i];          

              new_vz[i] = vz[i] - gamma*dt/2.0*vz[i] + dt/(2.0*m)*0.6022*(harter_z[i]+angle_z[i]+vdW_z[i]) + dt/(2.0*m)*sqrt((2.0*0.6022*kbT*gamma*m)/dt)*ranr_z[i];
              vz[i]=new_vz[i];          
          }

    }
    //for(int p=0; p<bead_number-1; ++p){
    //    myfile4 << p  << "   " << correl_angle[p]/float(bead_number-1-p)/float(N_steps-1-5000000) << endl;
    //}
//myfile1.close();
myfile3.close();
myfile6.close();
//myfile4.close();
}
return 0;
}
