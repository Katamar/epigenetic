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
            theta=acos(round(v1dotv2/(modv1*modv2)));      //theta is in radians, so should theta0 be
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
