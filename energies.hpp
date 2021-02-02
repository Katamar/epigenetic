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
    }
    return angleenergy;
}



double lj_energy(double *x_coor, double *y_coor, double *z_coor, double *vdW_x, double *vdW_y, double *vdW_z, int bead_number, double epsilon_lj, double sigma_lj){
    double total_lj=0;
    double rvector=0;
    for (int i=0; i<bead_number; ++i){
        for (int j=0; j<i; j++){
            if (j != i-1 && j != i+1 && j != i){
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
