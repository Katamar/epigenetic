double generateGaussianNoise(double mu, double sigma, int repetitions){
    double u1=0.5, u2=0.5, z0, z1, num0, num1;
    double two_pi = 2.0*3.14159265358979323846;
    int ctr=0;
    
    while(u1 != 0 && u1 !=1 && u2 != 0 && u2 != 1 && ctr<repetitions){
        double u1 = ((double)rand()/(double)RAND_MAX);
        double u2 = ((double)rand()/(double)RAND_MAX);
        if (u1==0) {u1=0.001;}
        if (u1==1) {u1=0.999;}
        z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
        z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
        num0 = z0 * sigma + mu;
        num1 = z1 * sigma + mu;
        ctr+=1;
    }
    return num0;
}
