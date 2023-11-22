//
//  implicit.c
//  
//
//  Created by ty on 2020/6/26.
//

#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

#include <stdarg.h>
#define INDEX(i,j) ((int) pow(2,m+1) - (int) pow(2,m-(j)+1)) + (i)/((int) pow(2,(j)))



void timestamp ();
int main (int argc, char *argv[]){
    int rank, nprocs;
    
    MPI_Status status;
    int mpi_error_code;
    
    mpi_error_code = MPI_Init (&argc, &argv);
    mpi_error_code = MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
    mpi_error_code = MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    
    double wall_time = MPI_Wtime ();
    if (rank == 0) timestamp ();
    
    int m = 13;
    int n = (int) round(log((double) nprocs)/log(2.));
    int nmesh = (int) pow(2,m)-1;
    int nghost = 2;
    int lenx = nmesh+2*nghost;
    double xmin = 2.0;
    double xmax = 10.0;
    double alpha_nu = 0.2;
    double alpha_d0 = 1.0;
    double alpha_dstart = 1000.0;
    double Calpha_dnull = 1.0e6;
    double tstart = 5.0*M_PI;
    double tCstart = 20.0*M_PI;
    double h0 = 0.1;    //@1AU
    int nout = 200;
    int dn = 5000;
    int nstep = dn*nout;
    double C_start = 0.2;
    double C0 = 0.2;
    double Cp0 = 0.2;
    double tpstart = 20.0*M_PI;
    double t = 0.0;
    double t0 = 10000.0;
    double omega0 = 0.2;
    double sigma_void = 1.0e-6;
    
    bool LOG = false;
    double xc = pow(omega0,-2./3.);
    double tau0 = 100.0;
    
    double Rstellar = 2.0;
    double mstellar = 1.0;
    double Bstellar = 1000.0;
    double omegastellar = 1.0;
    double Sigma0 = 10.0;    //at 1 AU kg/m^2 or 0.1 g/cm^2
    double beta = 100.0;
    double q = 0.25;
    double p = 1.0;
    
    const double G = 6.67259e-11;
    const double Rsolar = 6.957e8;
    const double msolar = 1.98847e30;
    const double r0 = 1.4959787e11;
    const double omegasolar = 2.865329e-6;
    const double Bgauss = 1.0e-4;
    const double year = 3600.0*24.0*365.0;
    const double mu0 = 4.0*M_PI*1.0E-7;
    double Rs = Rstellar*Rsolar;
    double ms = mstellar*msolar;
    double omegas = omegastellar*omegasolar;
    double Bs = Bstellar*Bgauss;
    double mdotd = 3.*M_PI*Sigma0*alpha_nu*pow(h0,2.0)*sqrt(G*ms)*sqrt(r0);
    double magmom = Bs*pow(Rs,3.0);
    double A = pow(magmom,2.0)/(mdotd*mu0*sqrt(G*ms*pow(Rs,7.0)));
    double omegastar_k = sqrt(G*ms/pow(Rs,3.));
    
    double Rc = pow(G*ms/pow(omegas,2.),1.0/3.0);
    
    double sigma0 = Sigma0*pow(Rs/r0,-p);
    sigma0 = sigma0*G*ms*pow(Rs,4.)*mu0/pow(magmom,2.);
    h0 = h0*pow(Rs/r0,q);
    double xx1 = pow(3.0*M_PI*alpha_nu*beta*h0*A/2.0,2./7.);
    double xx = pow(3.0*M_PI*alpha_nu*beta*h0*pow(xx1,q)*A/2.0,2./7.);
    double xd = 1.5*xx;
    double pmc = -2.9;
    
    double v0 = -1.5*alpha_nu*pow(h0,2.);
    
    double lx[lenx];
    double x[lenx];
    double dx[lenx-1];
    double dlx = (log10(xmax)-log10(xmin))/((double) lenx-1);
    double ddx = (xmax - xmin)/((double) lenx-1.);
    
    lx[0] = log10(xmin);
    for (int i=1; i < lenx;i++){
        lx[i] = lx[0] + i*dlx;
    }
    
    if (LOG){
        for (int i=0; i < lenx;i++){
            x[i] = pow(10.,lx[i]);
        }
        for (int i=0; i < lenx-1;i++){
            dx[i] = x[i]*(pow(10.,lx[i])-1.0);
        }
    }else{
        x[0] = xmin;
        for (int i=1; i < lenx;i++){
            x[i] = x[0] + i*ddx;
        }
        for (int i=0; i < lenx-1;i++){
            dx[i] = ddx;
        }
    }
    
    double precision = 0.001;
    double deltax = 0.1;
    double width = sqrt(-2.0*log(precision))*deltax;
    int nwidth = (int) (width/dx[0]) + 1;
    int nwind = (int) (2.*width+xc-xmin)/dx[0] + 1;
    
    double hdx = ddx/2.;
    double v[2][lenx];
    double omega[2][lenx];
    double sigma[2][lenx];
    double vi[lenx];
    double omegak[lenx];
    double B[lenx];
    double integ[lenx];
    double a [(int) pow(2,m+1)-2];
    double b [(int) pow(2,m+1)-2];
    double c [(int) pow(2,m+1)-2];
    double d [(int) pow(2,m+1)-2];
    double dtmin[nprocs];

    for (int i=0;i<lenx;i++){
        sigma[0][i] = sigma0/x[i];
        omega[0][i] = 0.0;
        v[0][i] = v0;
        omegak[i] = pow(x[i],-1.5);
    }
    
    integ[nghost] = 0.5*dx[0]*(exp(-pow(0.0*dx[0]/deltax,2.)/2.) + exp(-pow(0.5*dx[0]/deltax,2.)/2.))/2. + 0.5*deltax*sqrt(2.*M_PI);
    for (int i=1+nghost;i<nwidth+nghost;i++){
        integ[i] = integ[i-1] + dx[0]*(exp(-pow((i-0.5)*dx[0]/deltax,2.)/2.) + exp(-pow((i+0.5)*dx[0]/deltax,2.)/2.))/2.;
    }
    for (int i=nwidth+nghost;i<nmesh-nwidth+nghost;i++)
        integ[i] = deltax*sqrt(2.*M_PI);
    for (int i=0;i<nwidth;i++)
        integ[nmesh+nghost-1-i] = integ[nghost+i];
    
    double v_ff = -omegak[nghost]*xmin;
    FILE *output;
    char address[100];
    char dir[] = "/Users/ty/Desktop/magenetsphere/output/";
    if (rank==0){
        char fname[5] = "00000";
        strcpy(address, dir);
        strcat(address, fname);
        strcat(address, ".txt");
        output = fopen(address, "w");
        fprintf(output, "######TIME %f\n",t);
        for (int j=0;j<lenx;j++){
            fprintf(output, "%d %10.18f %10.18f %10.18f %10.18f\n",j,v[0][j],omega[0][j],sigma[0][j],B[j]);
        }
        fclose(output);
    }
    
    for (int i=0;i<nstep;i++){
        double alpha_dt;
        double C;
        double Cp;
        alpha_dt = 1./(1./alpha_dstart - (1./alpha_dstart-1./alpha_d0)/tstart*t);
        if (alpha_dt <= alpha_d0)
            alpha_dt = alpha_d0;
        if (t<=tstart)
            C = C_start;
        else if(t>=tCstart)
            C = C0;
        else
            C = (C_start*tCstart - C0*tstart + (C0 - C_start)*t)/(tCstart-tstart);
        
        Cp = Cp0 + (1.0-Cp0)/tpstart*t;
        if (Cp >= 1.0)
            Cp = 1.0;
        
        double dt1=1.0e10;
        double dt2;
        int lenproc = nmesh/nprocs + 1;
        int lenprocw = nwind/nprocs + 2;
        
        for (int j = nghost+lenproc*rank;j<nghost+lenproc*(rank+1) && j<nmesh+nghost;j++){
            if (dx[j]/fabs(v[0][j]) < dt1)
                dt1 = dx[j]/fabs(v[0][j]);
        }
        if (rank!=0){
            MPI_Send (&dt1, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        }else{
            for (int iproc = 1; iproc < nprocs; iproc++){
                MPI_Recv(&dt2, 1, MPI_DOUBLE, iproc, 1,
                         MPI_COMM_WORLD, &status);
                if (dt2 < dt1)
                    dt1 = dt2;
            }
        }
        
        double dt = C*dt1;
        MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        double d_force[lenx];
        double b_force[lenx];
        double d_flow[lenx];
        double wind[lenx];
        for (int j = nghost+lenproc*rank;j<nghost+lenproc*(rank+1) && j<nmesh+nghost;j++){
            double getalpha;
            if (x[j]<=xx)
                getalpha = 1.0;
            else if(x[j]<=xd)
                getalpha = pow(Calpha_dnull,(x[j]-xx)/(xd-xx));
            else
                getalpha = Calpha_dnull;
            
            double alpha_d = alpha_dt*getalpha;
            double gravity = 0.5*(pow(omega[0][j],2.) + 2.0*omega[0][j]/pow(x[j],1.5))*sigma[0][j]*x[j] + 0.5*(pow(omega[0][j+1],2.) + 2.0*omega[0][j+1]/pow(x[j+1],1.5))*sigma[0][j+1]*x[j+1];
            double pgradient = -pow(Cp*h0,2.)/sqrt(0.5*(x[j]+x[j+1]))*(sigma[0][j+1] - sigma[0][j])/(x[j+1] - x[j]);
            double Somega;
            if (0.5*(omega[0][j]+omegak[j]+omega[0][j+1]+omegak[j+1])>omega0)
                Somega = 1.;
            else
                Somega = -1.;
            
            double gc = fabs(0.5*(omega[0][j]+omegak[j]+omega[0][j+1]+omegak[j+1])-omega0)/(omega0*h0*pow(0.5*(x[j]+x[j+1]),0.25));
            double fc = gc/(1.+gc);
            double Pm = fc*v[0][j]/(alpha_d*h0*pow(0.5*(x[j]+x[j+1]),0.25)*0.5*(x[j]+x[j+1])*fabs(0.5*(omega[0][j]+omegak[j]+omega[0][j+1]+omegak[j+1])-omega0));
            if (Pm<pmc){
                Pm = pmc;
            }
            double grho = 3.*(3.+Pm)*(1.-pow(x[j]/xc,3.));
            double frho = -3.*(3.+Pm)/(1.+grho);
            double SU = (1.+Somega)*(1.-v[0][j]/fabs(v[0][j]))/4.;
            wind[j] = -2.*fc*SU*frho*v[0][j]/(alpha_d*pow(0.5*(x[j]+x[j+1]),6.5)*fabs(0.5*(omega[0][j]+omegak[j]+omega[0][j+1]+omegak[j+1])-omega0));
            //  wind[j] = 10.*wind[j];
            double Pr = -fc/sqrt(2.)/alpha_d/pow(0.5*(x[j]+x[j+1]),7)/fabs(0.5*(omega[0][j]+omegak[j]+omega[0][j+1]+omegak[j+1])-omega0);
            double Fd = - 2.*alpha_nu*pow(h0,2.0)*0.5*(sigma[0][j]+sigma[0][j+1])*0.5*(omega[0][j]+omegak[j]+omega[0][j+1]+omegak[j+1])*pow(0.5*(x[j]+x[j+1]),0.5);
            if (x[j]>xd){
                Pr = 0.0;
                wind[j] = 0.0;
            }
            B[j] = wind[j];
            double Fr = gravity + pgradient + 0.5*(Pr+Fd)*v[0][j];
            d_force[j] = Fr*dt;
            b_force[j] = -0.5*(Pr+Fd)*dt;
            /*
            if (rank != 0){
                MPI_Send (&wind[j], 1, MPI_DOUBLE, 0, j, MPI_COMM_WORLD);
            }else{
                for (int iproc = 1; iproc < nprocs; iproc++){
                    int k = j + lenproc*iproc;
                    if (k < nmesh+nghost){
                        MPI_Recv(&wind[k], 1, MPI_DOUBLE, iproc, k,
                                 MPI_COMM_WORLD, &status);
                    }
                }
            }*/
        }
        /*
        MPI_Barrier (MPI_COMM_WORLD);
        MPI_Bcast (wind, lenx, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Barrier (MPI_COMM_WORLD);

        for (int j = nghost+lenprocw*rank;j<nghost+lenprocw*(rank+1) && j<nwind+1;j++){
            B[j] = 0.0;
            if (j < nwidth + nghost){
                for (int k = nghost;k<nwidth+j;k++)
                    B[j] += exp(-pow((j-k)*dx[0]/deltax,2.)/2.)/integ[k]*wind[k]*dx[0];
            }else if (j>nmesh+nghost-nwidth-1){
                for (int k = j-nwidth+1;k<nmesh+nghost;k++)
                    B[j] += exp(-pow((j-k)*dx[0]/deltax,2.)/2.)/integ[k]*wind[k]*dx[0];
            }else{
                for (int k= j-nwidth+1;k<nwidth+j;k++)
                    B[j] += exp(-pow((j-k)*dx[0]/deltax,2.)/2.)/integ[k]*wind[k]*dx[0];
            }
            if (rank != 0){
                MPI_Send (&B[j], 1, MPI_DOUBLE, 0, j, MPI_COMM_WORLD);
            }else{
                for (int iproc = 1; iproc < nprocs; iproc++){
                    int k = j + lenprocw*iproc;
                    if (k < nwind+1){
                        MPI_Recv(&B[k], 1, MPI_DOUBLE, iproc, k,
                                 MPI_COMM_WORLD, &status);
                    }
                }
            }
        }
        MPI_Barrier (MPI_COMM_WORLD);
        MPI_Bcast (B, nwind+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Barrier (MPI_COMM_WORLD);
         */
        
        for (int j = nghost+lenproc*rank;j<nghost+lenproc*(rank+1) && j<nmesh+nghost;j++){
            double r = dt/(dx[j-1]+dx[j])/2.0;
            double s = alpha_nu*pow(h0,2.0)*dt/pow(dx[j],2.)/(0.5*(x[j]+x[j+1]));
            double vr = (dx[j+1]*v[0][j+1]+dx[j]*v[0][j])/(dx[j+1]+dx[j]);
            double vl = (dx[j]*v[0][j]+dx[j-1]*v[0][j-1])/(dx[j]+dx[j-1]);
            double a_advection;
            double b_advection;
            double c_advection;
            double d_advection;
            double a_diffusion;
            double b_diffusion;
            double c_diffusion;
            double d_diffusion;
            if (j == nghost){
                a_advection = 0.0;
                b_advection = 0.5*(sigma[0][j]+sigma[0][j+1])*(1. + r*(vr-vl));
                c_advection = r*0.5*(sigma[0][j]+sigma[0][j+1])*vr;
                d_advection = 0.5*(sigma[0][j]+sigma[0][j+1])*v[0][j] + r*0.5*(sigma[0][j]+sigma[0][j+1])*vl*v[0][j-1];
                a_diffusion = 0.0;
                b_diffusion = s*(sigma[0][j+1]*(omega[0][j+1]+omegak[j+1])*pow(x[j+1],3.5) + sigma[0][j]*(omega[0][j]+omegak[j])*pow(x[j],3.5));
                c_diffusion = -s*sigma[0][j+1]*(omega[0][j+1]+omegak[j+1])*pow(x[j+1],3.5);
                d_diffusion = s*(sigma[0][j+1]*(omega[0][j+1]+omegak[j+1])*pow(x[j+1],3.5)*v[0][j+1] - (sigma[0][j+1]*(omega[0][j+1]+omegak[j+1])*pow(x[j+1],3.5) + sigma[0][j]*(omega[0][j]+omegak[j])*pow(x[j],3.5))*v[0][j] + 2.0*sigma[0][j]*(omega[0][j]+omegak[j])*pow(x[j],3.5)*v[0][j-1]);
            }else if (j==nmesh+nghost-1){
                a_advection = -r*0.5*(sigma[0][j]+sigma[0][j+1])*vl;
                b_advection = 0.5*(sigma[0][j]+sigma[0][j+1])*(1. + r*(vr-vl));
                c_advection = 0.0;
                d_advection = 0.5*(sigma[0][j]+sigma[0][j+1])*v[0][j] - r*0.5*(sigma[0][j]+sigma[0][j+1])*vr*v[0][j+1];
                a_diffusion = -s*sigma[0][j]*(omega[0][j]+omegak[j])*pow(x[j],3.5);
                b_diffusion = s*(sigma[0][j+1]*(omega[0][j+1]+omegak[j+1])*pow(x[j+1],3.5) + sigma[0][j]*(omega[0][j]+omegak[j])*pow(x[j],3.5));
                c_diffusion = 0.0;
                d_diffusion = s*(2.0*sigma[0][j+1]*(omega[0][j+1]+omegak[j+1])*pow(x[j+1],3.5)*v[0][j+1] - (sigma[0][j+1]*(omega[0][j+1]+omegak[j+1])*pow(x[j+1],3.5) + sigma[0][j]*(omega[0][j]+omegak[j])*pow(x[j],3.5))*v[0][j] + sigma[0][j]*(omega[0][j]+omegak[j])*pow(x[j],3.5)*v[0][j-1]);
            }else{
                a_advection = -r*0.5*(sigma[0][j]+sigma[0][j+1])*vl;
                b_advection = 0.5*(sigma[0][j]+sigma[0][j+1])*(1. + r*(vr-vl));
                c_advection = r*0.5*(sigma[0][j]+sigma[0][j+1])*vr;
                d_advection = 0.5*(sigma[0][j]+sigma[0][j+1])*v[0][j];
                a_diffusion = -s*sigma[0][j]*(omega[0][j]+omegak[j])*pow(x[j],3.5);
                b_diffusion = s*(sigma[0][j+1]*(omega[0][j+1]+omegak[j+1])*pow(x[j+1],3.5) + sigma[0][j]*(omega[0][j]+omegak[j])*pow(x[j],3.5));
                c_diffusion = -s*sigma[0][j+1]*(omega[0][j+1]+omegak[j+1])*pow(x[j+1],3.5);
                d_diffusion = s*(sigma[0][j+1]*(omega[0][j+1]+omegak[j+1])*pow(x[j+1],3.5)*v[0][j+1] - (sigma[0][j+1]*(omega[0][j+1]+omegak[j+1])*pow(x[j+1],3.5) + sigma[0][j]*(omega[0][j]+omegak[j])*pow(x[j],3.5))*v[0][j] + sigma[0][j]*(omega[0][j]+omegak[j])*pow(x[j],3.5)*v[0][j-1]);
            }
           // b_force[j] += -0.5*B[j]*0.5*(x[j]+x[j+1])*dt;
            //d_force[j] += 0.5*B[j]*0.5*(x[j]+x[j+1])*v[0][j]*dt;
            a[INDEX(j-nghost+1,0)] = a_advection + a_diffusion;
            b[INDEX(j-nghost+1,0)] = b_advection + b_diffusion + b_force[j];
            c[INDEX(j-nghost+1,0)] = c_advection + c_diffusion;
            d[INDEX(j-nghost+1,0)] = d_advection + d_diffusion + d_force[j];
            
        }
        if (rank != 0){
            MPI_Send (&a[INDEX(1 + lenproc*rank,0)], 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
            MPI_Send (&b[INDEX(1 + lenproc*rank,0)], 1, MPI_DOUBLE, rank-1, 2, MPI_COMM_WORLD);
            MPI_Send (&c[INDEX(1 + lenproc*rank,0)], 1, MPI_DOUBLE, rank-1, 3, MPI_COMM_WORLD);
            MPI_Send (&d[INDEX(1 + lenproc*rank,0)], 1, MPI_DOUBLE, rank-1, 4, MPI_COMM_WORLD);
        }
        if (rank != nprocs-1){
            MPI_Recv(&a[INDEX(1 + lenproc*(rank+1),0)], 1, MPI_DOUBLE, rank+1, 1,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(&b[INDEX(1 + lenproc*(rank+1),0)], 1, MPI_DOUBLE, rank+1, 2,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(&c[INDEX(1 + lenproc*(rank+1),0)], 1, MPI_DOUBLE, rank+1, 3,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(&d[INDEX(1 + lenproc*(rank+1),0)], 1, MPI_DOUBLE, rank+1, 4,
                     MPI_COMM_WORLD, &status);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
        for (int j = 1; j < m; j++){
            
            int lenproc = (pow(2, m-j)-1)/ nprocs + 1;
            if (lenproc == 0) lenproc = 1;
            int s = (int) pow(2, j-1);
            
            for (int l = (1 + lenproc*rank)*2*s;
                 l < (1 + lenproc*(rank+1))*2*s
                 && l <= nmesh+1 - (int) pow(2, j);
                 l = l + 2*s){
                
                double e = - a[INDEX(l,j-1)] / b[INDEX(l-s,j-1)];
                double f = - c[INDEX(l,j-1)] / b[INDEX(l+s,j-1)];
                
                a[INDEX(l,j)] = e * a[INDEX(l-s,j-1)];
                c[INDEX(l,j)] = f * c[INDEX(l+s,j-1)];
                
                b[INDEX(l,j)] = b[INDEX(l,j-1)] +
                e * c[INDEX(l-s,j-1)] + f * a[INDEX(l+s,j-1)];
                
                d[INDEX(l,j)] = d[INDEX(l,j-1)] +
                e * d[INDEX(l-s,j-1)] + f * d[INDEX(l+s,j-1)];
            }
            if (lenproc != 1){
                if (rank != 0){
                    MPI_Send (&a[INDEX((1 + lenproc*rank) * 2*s,j)], 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
                    MPI_Send (&b[INDEX((1 + lenproc*rank) * 2*s,j)], 1, MPI_DOUBLE, rank-1, 2, MPI_COMM_WORLD);
                    MPI_Send (&c[INDEX((1 + lenproc*rank) * 2*s,j)], 1, MPI_DOUBLE, rank-1, 3, MPI_COMM_WORLD);
                    MPI_Send (&d[INDEX((1 + lenproc*rank) * 2*s,j)], 1, MPI_DOUBLE, rank-1, 4, MPI_COMM_WORLD);
                }
                if (rank != nprocs-1){
                    MPI_Recv(&a[INDEX((1 + lenproc*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, rank+1, 1,
                             MPI_COMM_WORLD, &status);
                    MPI_Recv(&b[INDEX((1 + lenproc*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, rank+1, 2,
                             MPI_COMM_WORLD, &status);
                    MPI_Recv(&c[INDEX((1 + lenproc*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, rank+1, 3,
                             MPI_COMM_WORLD, &status);
                    MPI_Recv(&d[INDEX((1 + lenproc*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, rank+1, 4,
                             MPI_COMM_WORLD, &status);
                }
            }else if(j<m-1){
                if (rank%2 == 1 && rank<pow(2, m-j)-1){
                    MPI_Send (&a[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2, 1, MPI_COMM_WORLD);
                    MPI_Send (&b[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2, 2, MPI_COMM_WORLD);
                    MPI_Send (&c[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2, 3, MPI_COMM_WORLD);
                    MPI_Send (&d[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2, 4, MPI_COMM_WORLD);
                }else if (rank != 0 && rank<pow(2, m-j)-2){
                    for (int k=0;k<2;k++){
                        MPI_Send (&a[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-k, 1, MPI_COMM_WORLD);
                        MPI_Send (&b[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-k, 2, MPI_COMM_WORLD);
                        MPI_Send (&c[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-k, 3, MPI_COMM_WORLD);
                        MPI_Send (&d[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-k, 4, MPI_COMM_WORLD);
                    }
                }else if (rank==pow(2, m-j)-2){
                    MPI_Send (&a[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-1, 1, MPI_COMM_WORLD);
                    MPI_Send (&b[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-1, 2, MPI_COMM_WORLD);
                    MPI_Send (&c[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-1, 3, MPI_COMM_WORLD);
                    MPI_Send (&d[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-1, 4, MPI_COMM_WORLD);
                }
                
                if (rank<pow(2, m-j-1)-1){
                    MPI_Recv (&a[INDEX((1 + rank) * 4*s,j)], 1, MPI_DOUBLE, 2*rank+1, 1, MPI_COMM_WORLD,&status);
                    MPI_Recv (&b[INDEX((1 + rank) * 4*s,j)], 1, MPI_DOUBLE, 2*rank+1, 2, MPI_COMM_WORLD,&status);
                    MPI_Recv (&c[INDEX((1 + rank) * 4*s,j)], 1, MPI_DOUBLE, 2*rank+1, 3, MPI_COMM_WORLD,&status);
                    MPI_Recv (&d[INDEX((1 + rank) * 4*s,j)], 1, MPI_DOUBLE, 2*rank+1, 4, MPI_COMM_WORLD,&status);
                    if (rank !=0 ){
                        for (int k=0;k<2;k++){
                            MPI_Recv(&a[INDEX((1+2*(rank+k)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+k), 1,
                                     MPI_COMM_WORLD, &status);
                            MPI_Recv(&b[INDEX((1+2*(rank+k)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+k), 2,
                                     MPI_COMM_WORLD, &status);
                            MPI_Recv(&c[INDEX((1+2*(rank+k)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+k), 3,
                                     MPI_COMM_WORLD, &status);
                            MPI_Recv(&d[INDEX((1+2*(rank+k)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+k), 4,
                                     MPI_COMM_WORLD, &status);
                        }
                    }
                    else{
                        MPI_Recv(&a[INDEX((1+2*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+1), 1,
                                 MPI_COMM_WORLD, &status);
                        MPI_Recv(&b[INDEX((1+2*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+1), 2,
                                 MPI_COMM_WORLD, &status);
                        MPI_Recv(&c[INDEX((1+2*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+1), 3,
                                 MPI_COMM_WORLD, &status);
                        MPI_Recv(&d[INDEX((1+2*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+1), 4,
                                 MPI_COMM_WORLD, &status);
                    }
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        
        
        int middle = (int) pow(2,m-1);
        if (rank==0)
            v[1][middle+nghost-1] = d[INDEX(middle,m-1)] / b[INDEX(middle,m-1)];
        
        
        for (int j = m-1; j > 0; j--){
            
            int lenproc = (pow(2, m-j+1)-1)/nprocs+1;
            if (lenproc == 0) lenproc = 1;
            
            int s = (int) pow(2, j-1);
            
            
            if (lenproc==1){
                int l = (1 + rank) * 2*s;
                int k = (1 + rank) * s;
                if (rank != 0&& rank <(int)  pow(2, m-j)-1){
                    MPI_Send (&v[1][l+nghost-1], 1, MPI_DOUBLE, 2*rank, 1, MPI_COMM_WORLD);
                    MPI_Send (&v[1][l+nghost-1], 1, MPI_DOUBLE, 2*rank+1, 2, MPI_COMM_WORLD);
                    MPI_Send (&v[1][l+nghost-1], 1, MPI_DOUBLE, 2*rank+2, 3, MPI_COMM_WORLD);
                }else if(rank==0 && rank <(int)  pow(2, m-j)-1){
                    MPI_Send (&v[1][l+nghost-1], 1, MPI_DOUBLE, 2*rank+1, 2, MPI_COMM_WORLD);
                    MPI_Send (&v[1][l+nghost-1], 1, MPI_DOUBLE, 2*rank+2, 3, MPI_COMM_WORLD);
                }
                if (rank != 0 && rank%2==0 && rank <(int)  pow(2, m-j+1)-2){
                    MPI_Recv(&v[1][k-s+nghost-1], 1, MPI_DOUBLE, rank/2-1, 3,
                             MPI_COMM_WORLD, &status);
                    MPI_Recv(&v[1][k+s+nghost-1], 1, MPI_DOUBLE, rank/2, 1,
                             MPI_COMM_WORLD, &status);
                }else if(rank%2==1 && rank <(int)  pow(2, m-j+1)-2){
                    MPI_Recv(&v[1][k+nghost-1], 1, MPI_DOUBLE, rank/2, 2,
                             MPI_COMM_WORLD, &status);
                }else if(rank == (int)  pow(2, m-j+1)-2){
                    MPI_Recv(&v[1][k-s+nghost-1], 1, MPI_DOUBLE, rank/2-1, 3,
                             MPI_COMM_WORLD, &status);
                }
                
                
                MPI_Barrier(MPI_COMM_WORLD);
                if (rank%2==0 && rank<(int) pow(2,m-j+1)-1){
                    if (k==s){
                        v[1][k+nghost-1] = (d[INDEX(k,j-1)] - c[INDEX(k,j-1)] * v[1][k+s+nghost-1])
                        / b[INDEX(k,j-1)];
                    }else if (k==nmesh+1-s){
                        v[1][k+nghost-1] = (d[INDEX(k,j-1)] -
                                          a[INDEX(k,j-1)] * v[1][k-s+nghost-1])
                        / b[INDEX(k,j-1)];
                    }else{
                        v[1][k+nghost-1] = (d[INDEX(k,j-1)] -
                                          a[INDEX(k,j-1)] * v[1][k-s+nghost-1] - c[INDEX(k,j-1)] * v[1][k+s+nghost-1])
                        / b[INDEX(k,j-1)];
                    }
                }
            }else{
                if(rank!=nprocs-1){
                    MPI_Send (&v[1][lenproc*(rank+1)*s+nghost-1], 1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD);
                }
                if(rank!=0){
                    MPI_Recv(&v[1][lenproc*rank*s+nghost-1], 1, MPI_DOUBLE, rank-1, 1,
                             MPI_COMM_WORLD, &status);
                }
                MPI_Barrier(MPI_COMM_WORLD);
                for (int l = (1 + lenproc*rank) * s;
                     (l < (1 + lenproc*(rank + 1)) * s) &&
                     (l <= nmesh+1 - (int) pow (2,j-1));
                     l = l + 2 * s){
                    
                    if (l==s){
                        v[1][l+nghost-1] = (d[INDEX(l,j-1)] - c[INDEX(l,j-1)] * v[1][l+s+nghost-1])
                        / b[INDEX(l,j-1)];
                    }else if (l==nmesh+1-s){
                        v[1][l+nghost-1] = (d[INDEX(l,j-1)] -
                                          a[INDEX(l,j-1)] * v[1][l-s+nghost-1])
                        / b[INDEX(l,j-1)];
                    }else{
                        v[1][l+nghost-1] = (d[INDEX(l,j-1)] -
                                          a[INDEX(l,j-1)] * v[1][l-s+nghost-1] - c[INDEX(l,j-1)] * v[1][l+s+nghost-1])
                        / b[INDEX(l,j-1)];
                    }
                }
            }
        }
        
        for (int j=0;j<nghost;j++){
            if (rank==0)
                v[1][j] = v[1][j+nghost];
            if (rank==nprocs-1)
                v[1][lenx-1-j] = v[1][lenx-1-nghost];
                //v[1][lenx-1-j] = v0;
        }
        
        if (rank != 0){
            MPI_Send (&v[1][lenproc*rank+nghost], 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
            for (int k=0;k<2;k++){
                MPI_Recv(&v[1][lenproc*rank-k+nghost-1], 1, MPI_DOUBLE, rank-1, 2+k,
                         MPI_COMM_WORLD, &status);
            }
        }
        if (rank != nprocs-1){
            MPI_Recv(&v[1][lenproc*(rank+1)+nghost], 1, MPI_DOUBLE, rank+1, 1,
                     MPI_COMM_WORLD, &status);
            for (int k=0;k<2;k++){
                MPI_Send (&v[1][lenproc*(rank+1)-k+nghost-1], 1, MPI_DOUBLE, rank+1, 2+k, MPI_COMM_WORLD);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
        for (int j = nghost+lenproc*rank;j<nghost+lenproc*(rank+1) && j<nmesh+nghost;j++){
            double getalpha;
            if (x[j]<=xx)
                getalpha = 1.0;
            else if(x[j]<=xd)
                getalpha = pow(Calpha_dnull,(x[j]-xx)/(xd-xx));
            else
                getalpha = Calpha_dnull;
            
            double alpha_d = alpha_dt*getalpha;
            double vh_new = (dx[j]*v[1][j]+dx[j-1]*v[1][j-1])/(dx[j]+dx[j-1]);
            double vh_old = (dx[j]*v[0][j]+dx[j-1]*v[0][j-1])/(dx[j]+dx[j-1]);
            double vh = 0.5*(vh_new+vh_old);
            double Somega;
            if (omega[0][j]+omegak[j]>omega0)
                Somega = 1.;
            else
                Somega = -1.;
            double gc = fabs(omega[0][j]+omegak[j]-omega0)/(omega0*h0*pow(x[j],0.25));
            double fc = gc/(1.+gc);
            double Pm = fc*vh/(alpha_d*h0*pow(x[j],0.25)*x[j]*fabs(omega[0][j]+omegak[j]-omega0));
            if (Pm<pmc){
                Pm = pmc;
            }
            double grho = 3.*(3.+Pm)*(1.-pow((x[j]/xc),3.));
            double frho = -3.*(3.+Pm)/(1.+grho);
            double SU = (1.+Somega)*(1.-vh/fabs(vh))/4.;
            wind[j] = -2.*fc*SU*frho*vh/(alpha_d*pow(x[j],6.5)*fabs(omega[0][j]+omegak[j]-omega0));
            //   wind[j] = 10.*wind[j];
            double torque = -(sqrt(2.)*fc*Somega)/(alpha_d*pow(x[j],6.)*M_PI);
            double Pphi = -2.*fc*SU*frho*vh/(alpha_d*pow(x[j],6.)*fabs(omega[0][j]+omegak[j]-omega0));
            double Fphi = torque + Pphi;
            if (x[j]>xd){
                Fphi = 0.0;
                wind[j] = 0.0;
            }
            B[j] = wind[j];
            double Fomega = omegak[j]*0.5*((sigma[0][j+1]*x[j+1]+sigma[0][j]*x[j])*v[1][j] - (sigma[0][j]*x[j]+sigma[0][j-1]*x[j-1])*v[1][j-1])/dx[j];
            d_force[j] = (Fphi+Fomega)*pow(x[j],2.)*dt;
            /*
            if (rank != 0){
                MPI_Send (&wind[j], 1, MPI_DOUBLE, 0, j, MPI_COMM_WORLD);
            }else{
                for (int iproc = 1; iproc < nprocs; iproc++){
                    int k = j + lenproc*iproc;
                    if (k < nmesh+nghost){
                        MPI_Recv(&wind[k], 1, MPI_DOUBLE, iproc, k,
                                 MPI_COMM_WORLD, &status);
                    }
                }
            }*/
        }
        /*
        MPI_Barrier (MPI_COMM_WORLD);
        MPI_Bcast (wind, lenx, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Barrier (MPI_COMM_WORLD);
        
        for (int j = nghost+lenprocw*rank;j<nghost+lenprocw*(rank+1) && j<nwind+1;j++){
            B[j] = 0.0;
            if (j < nwidth + nghost){
                for (int k = nghost;k<nwidth+j;k++)
                    B[j] += exp(-pow((j-k)*dx[0]/deltax,2.)/2.)/integ[k]*wind[k]*dx[0];
            }else if (j>nmesh+nghost-nwidth-1){
                for (int k = j-nwidth+1;k<nmesh+nghost;k++)
                    B[j] += exp(-pow((j-k)*dx[0]/deltax,2.)/2.)/integ[k]*wind[k]*dx[0];
            }else{
                for (int k= j-nwidth+1;k<nwidth+j;k++)
                    B[j] += exp(-pow((j-k)*dx[0]/deltax,2.)/2.)/integ[k]*wind[k]*dx[0];
            }
            if (rank != 0){
                MPI_Send (&B[j], 1, MPI_DOUBLE, 0, j, MPI_COMM_WORLD);
            }else{
                for (int iproc = 1; iproc < nprocs; iproc++){
                    int k = j + lenprocw*iproc;
                    if (k < nwind+1){
                        MPI_Recv(&B[k], 1, MPI_DOUBLE, iproc, k,
                                 MPI_COMM_WORLD, &status);
                    }
                }
            }
        }
        MPI_Barrier (MPI_COMM_WORLD);
        MPI_Bcast (B, nwind+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Barrier (MPI_COMM_WORLD);
         */
        
        for (int j = nghost+lenproc*rank;j<nghost+lenproc*(rank+1) && j<nmesh+nghost;j++){
            double r = dt/(dx[j-1]+dx[j])/2.0;
            double sl = alpha_nu*pow(h0,2.0)*dt/(dx[j-1]+dx[j])/2./dx[j-1];
            double sr = alpha_nu*pow(h0,2.0)*dt/(dx[j-1]+dx[j])/2./dx[j];
            double vr_new = (dx[j+1]*v[1][j+1]+dx[j]*v[1][j])/(dx[j+1]+dx[j]);
            double vl_new = (dx[j-1]*v[1][j-1]+dx[j-2]*v[1][j-2])/(dx[j-1]+dx[j-2]);
            double vr_old = (dx[j+1]*v[0][j+1]+dx[j]*v[0][j])/(dx[j+1]+dx[j]);
            double vl_old = (dx[j-1]*v[0][j-1]+dx[j-2]*v[0][j-2])/(dx[j-1]+dx[j-2]);
            double a_advection;
            double b_advection;
            double c_advection;
            double d_advection;
            double a_diffusion;
            double b_diffusion;
            double c_diffusion;
            double d_diffusion;
            if (j == nghost){
                a_advection = 0.0;
                b_advection = sigma[0][j]*pow(x[j],3.) - r*((sigma[0][j+1]*x[j+1]+sigma[0][j]*x[j])*v[1][j] - (sigma[0][j]*x[j]+sigma[0][j-1]*x[j-1])*v[1][j-1])*pow(x[j],2.);
                c_advection = r*sigma[0][j+1]*pow(x[j+1],3.)*vr_new - r*sigma[0][j-1]*pow(x[j-1],3.)*vl_new;
                d_advection = -r*sigma[0][j+1]*pow(x[j+1],3.)*vr_old*omega[0][j+1] + sigma[0][j]*pow(x[j],3.)*omega[0][j] + r*sigma[0][j-1]*pow(x[j-1],3.)*vl_old*omega[0][j-1] - r*sigma[0][j+1]*pow(x[j+1],1.5)*(vr_old+vr_new-2*v0) + r*sigma[0][j-1]*pow(x[j-1],1.5)*(vl_old+vl_new-2*v0) + r*((sigma[0][j+1]*x[j+1]+sigma[0][j]*x[j])*v[0][j] - (sigma[0][j]*x[j]+sigma[0][j-1]*x[j-1])*v[0][j-1])*omega[0][j]*pow(x[j],2.);
                a_diffusion = 0.0;
                b_diffusion = (sl*(pow(x[j-1],5.5)*sigma[0][j-1] + pow(x[j],5.5)*sigma[0][j]) + sr*(pow(x[j],5.5)*sigma[0][j] + pow(x[j+1],5.5)*sigma[0][j+1]))*(omega[0][j] + pow(x[j],-1.5));
                c_diffusion = -sr*(pow(x[j],5.5)*sigma[0][j] + pow(x[j+1],5.5)*sigma[0][j+1])*(omega[0][j+1] + pow(x[j+1],-1.5)) - sl*(pow(x[j-1],5.5)*sigma[0][j-1] + pow(x[j],5.5)*sigma[0][j])*(omega[0][j-1] + pow(x[j-1],-1.5));
                d_diffusion = sl*(pow(x[j-1],5.5)*sigma[0][j-1] + pow(x[j],5.5)*sigma[0][j])*omega[0][j-1]/pow(x[j-1],1.5) - (sl*(pow(x[j-1],5.5)*sigma[0][j-1] + pow(x[j],5.5)*sigma[0][j]) + sr*(pow(x[j],5.5)*sigma[0][j] + pow(x[j+1],5.5)*sigma[0][j+1]))*omega[0][j]/pow(x[j],1.5) + sr*(pow(x[j],5.5)*sigma[0][j] + pow(x[j+1],5.5)*sigma[0][j+1])*omega[0][j+1]/pow(x[j+1],1.5);
            } else if (j==nmesh+nghost-1){
                a_advection = -r*sigma[0][j-1]*pow(x[j-1],3.)*vl_new;
                b_advection = sigma[0][j]*pow(x[j],3.) - r*((sigma[0][j+1]*x[j+1]+sigma[0][j]*x[j])*v[1][j] - (sigma[0][j]*x[j]+sigma[0][j-1]*x[j-1])*v[1][j-1])*pow(x[j],2.);
                c_advection = 0.0;
                d_advection = -r*sigma[0][j+1]*pow(x[j+1],3.)*(vr_old+vr_new)*omega[0][j+1] + sigma[0][j]*pow(x[j],3.)*omega[0][j] + r*sigma[0][j-1]*pow(x[j-1],3.)*vl_old*omega[0][j-1] - r*sigma[0][j+1]*pow(x[j+1],1.5)*(vr_old+vr_new-2*v0) + r*sigma[0][j-1]*pow(x[j-1],1.5)*(vl_old+vl_new-2*v0) + r*((sigma[0][j+1]*x[j+1]+sigma[0][j]*x[j])*v[0][j] - (sigma[0][j]*x[j]+sigma[0][j-1]*x[j-1])*v[0][j-1])*omega[0][j]*pow(x[j],2.);
                a_diffusion = -sl*(pow(x[j-1],5.5)*sigma[0][j-1] + pow(x[j],5.5)*sigma[0][j])*(omega[0][j-1] + pow(x[j-1],-1.5));
                b_diffusion = (sl*(pow(x[j-1],5.5)*sigma[0][j-1] + pow(x[j],5.5)*sigma[0][j]) + sr*(pow(x[j],5.5)*sigma[0][j] + pow(x[j+1],5.5)*sigma[0][j+1]))*(omega[0][j] + pow(x[j],-1.5));
                c_diffusion = 0.0;
                d_diffusion = sl*(pow(x[j-1],5.5)*sigma[0][j-1] + pow(x[j],5.5)*sigma[0][j])*omega[0][j-1]/pow(x[j-1],1.5) - (sl*(pow(x[j-1],5.5)*sigma[0][j-1] + pow(x[j],5.5)*sigma[0][j]) + sr*(pow(x[j],5.5)*sigma[0][j] + pow(x[j+1],5.5)*sigma[0][j+1]))*omega[0][j]/pow(x[j],1.5) + sr*(pow(x[j],5.5)*sigma[0][j] + pow(x[j+1],5.5)*sigma[0][j+1])*omega[0][j+1]*(omega[0][j+1] + 2.*pow(x[j+1],-1.5));
            }else{
                a_advection = -r*sigma[0][j-1]*pow(x[j-1],3.)*vl_new;
                b_advection = sigma[0][j]*pow(x[j],3.) - r*((sigma[0][j+1]*x[j+1]+sigma[0][j]*x[j])*v[1][j] - (sigma[0][j]*x[j]+sigma[0][j-1]*x[j-1])*v[1][j-1])*pow(x[j],2.);
                c_advection = r*sigma[0][j+1]*pow(x[j+1],3.)*vr_new;
                d_advection = -r*sigma[0][j+1]*pow(x[j+1],3.)*vr_old*omega[0][j+1] + sigma[0][j]*pow(x[j],3.)*omega[0][j] + r*sigma[0][j-1]*pow(x[j-1],3.)*vl_old*omega[0][j-1] - r*sigma[0][j+1]*pow(x[j+1],1.5)*(vr_old+vr_new-2*v0) + r*sigma[0][j-1]*pow(x[j-1],1.5)*(vl_old+vl_new-2*v0) + r*((sigma[0][j+1]*x[j+1]+sigma[0][j]*x[j])*v[0][j] - (sigma[0][j]*x[j]+sigma[0][j-1]*x[j-1])*v[0][j-1])*omega[0][j]*pow(x[j],2.);
                a_diffusion = -sl*(pow(x[j-1],5.5)*sigma[0][j-1] + pow(x[j],5.5)*sigma[0][j])*(omega[0][j-1] + pow(x[j-1],-1.5));
                b_diffusion = (sl*(pow(x[j-1],5.5)*sigma[0][j-1] + pow(x[j],5.5)*sigma[0][j]) + sr*(pow(x[j],5.5)*sigma[0][j] + pow(x[j+1],5.5)*sigma[0][j+1]))*(omega[0][j] + pow(x[j],-1.5));
                c_diffusion = -sr*(pow(x[j],5.5)*sigma[0][j] + pow(x[j+1],5.5)*sigma[0][j+1])*(omega[0][j+1] + pow(x[j+1],-1.5));
                d_diffusion = sl*(pow(x[j-1],5.5)*sigma[0][j-1] + pow(x[j],5.5)*sigma[0][j])*omega[0][j-1]/pow(x[j-1],1.5) - (sl*(pow(x[j-1],5.5)*sigma[0][j-1] + pow(x[j],5.5)*sigma[0][j]) + sr*(pow(x[j],5.5)*sigma[0][j] + pow(x[j+1],5.5)*sigma[0][j+1]))*omega[0][j]/pow(x[j],1.5) + sr*(pow(x[j],5.5)*sigma[0][j] + pow(x[j+1],5.5)*sigma[0][j+1])*omega[0][j+1]/pow(x[j+1],1.5);
            }
            
            //b_force[j] = - 0.5*pow(x[j],3.)*B[j]*dt;
            //d_force[j] += 0.5*omega[0][j]*pow(x[j],3.)*B[j]*dt;
            a[INDEX(j-nghost+1,0)] = a_advection + a_diffusion;
            b[INDEX(j-nghost+1,0)] = b_advection + b_diffusion;
            c[INDEX(j-nghost+1,0)] = c_advection + c_diffusion;
            d[INDEX(j-nghost+1,0)] = d_advection + d_diffusion + d_force[j];
            
        }
        if (rank != 0){
            MPI_Send (&a[INDEX(1 + lenproc*rank,0)], 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
            MPI_Send (&b[INDEX(1 + lenproc*rank,0)], 1, MPI_DOUBLE, rank-1, 2, MPI_COMM_WORLD);
            MPI_Send (&c[INDEX(1 + lenproc*rank,0)], 1, MPI_DOUBLE, rank-1, 3, MPI_COMM_WORLD);
            MPI_Send (&d[INDEX(1 + lenproc*rank,0)], 1, MPI_DOUBLE, rank-1, 4, MPI_COMM_WORLD);
        }
        if (rank != nprocs-1){
            MPI_Recv(&a[INDEX(1 + lenproc*(rank+1),0)], 1, MPI_DOUBLE, rank+1, 1,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(&b[INDEX(1 + lenproc*(rank+1),0)], 1, MPI_DOUBLE, rank+1, 2,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(&c[INDEX(1 + lenproc*(rank+1),0)], 1, MPI_DOUBLE, rank+1, 3,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(&d[INDEX(1 + lenproc*(rank+1),0)], 1, MPI_DOUBLE, rank+1, 4,
                     MPI_COMM_WORLD, &status);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
        for (int j = 1; j < m; j++){
            
            int lenproc = (pow(2, m-j)-1)/ nprocs + 1;
            if (lenproc == 0) lenproc = 1;
            int s = (int) pow(2, j-1);
            
            for (int l = (1 + lenproc*rank)*2*s;
                 l < (1 + lenproc*(rank+1))*2*s
                 && l <= nmesh+1 - (int) pow(2, j);
                 l = l + 2*s){
                
                double e = - a[INDEX(l,j-1)] / b[INDEX(l-s,j-1)];
                double f = - c[INDEX(l,j-1)] / b[INDEX(l+s,j-1)];
                
                a[INDEX(l,j)] = e * a[INDEX(l-s,j-1)];
                c[INDEX(l,j)] = f * c[INDEX(l+s,j-1)];
                
                b[INDEX(l,j)] = b[INDEX(l,j-1)] +
                e * c[INDEX(l-s,j-1)] + f * a[INDEX(l+s,j-1)];
                
                d[INDEX(l,j)] = d[INDEX(l,j-1)] +
                e * d[INDEX(l-s,j-1)] + f * d[INDEX(l+s,j-1)];
            }
            if (lenproc != 1){
                if (rank != 0){
                    MPI_Send (&a[INDEX((1 + lenproc*rank) * 2*s,j)], 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
                    MPI_Send (&b[INDEX((1 + lenproc*rank) * 2*s,j)], 1, MPI_DOUBLE, rank-1, 2, MPI_COMM_WORLD);
                    MPI_Send (&c[INDEX((1 + lenproc*rank) * 2*s,j)], 1, MPI_DOUBLE, rank-1, 3, MPI_COMM_WORLD);
                    MPI_Send (&d[INDEX((1 + lenproc*rank) * 2*s,j)], 1, MPI_DOUBLE, rank-1, 4, MPI_COMM_WORLD);
                }
                if (rank != nprocs-1){
                    MPI_Recv(&a[INDEX((1 + lenproc*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, rank+1, 1,
                             MPI_COMM_WORLD, &status);
                    MPI_Recv(&b[INDEX((1 + lenproc*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, rank+1, 2,
                             MPI_COMM_WORLD, &status);
                    MPI_Recv(&c[INDEX((1 + lenproc*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, rank+1, 3,
                             MPI_COMM_WORLD, &status);
                    MPI_Recv(&d[INDEX((1 + lenproc*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, rank+1, 4,
                             MPI_COMM_WORLD, &status);
                }
            }else if(j<m-1){
                if (rank%2 == 1 && rank<pow(2, m-j)-1){
                    MPI_Send (&a[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2, 1, MPI_COMM_WORLD);
                    MPI_Send (&b[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2, 2, MPI_COMM_WORLD);
                    MPI_Send (&c[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2, 3, MPI_COMM_WORLD);
                    MPI_Send (&d[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2, 4, MPI_COMM_WORLD);
                }else if (rank != 0 && rank<pow(2, m-j)-2){
                    for (int k=0;k<2;k++){
                        MPI_Send (&a[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-k, 1, MPI_COMM_WORLD);
                        MPI_Send (&b[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-k, 2, MPI_COMM_WORLD);
                        MPI_Send (&c[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-k, 3, MPI_COMM_WORLD);
                        MPI_Send (&d[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-k, 4, MPI_COMM_WORLD);
                    }
                }else if (rank==pow(2, m-j)-2){
                    MPI_Send (&a[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-1, 1, MPI_COMM_WORLD);
                    MPI_Send (&b[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-1, 2, MPI_COMM_WORLD);
                    MPI_Send (&c[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-1, 3, MPI_COMM_WORLD);
                    MPI_Send (&d[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-1, 4, MPI_COMM_WORLD);
                }
                
                if (rank<pow(2, m-j-1)-1){
                    MPI_Recv (&a[INDEX((1 + rank) * 4*s,j)], 1, MPI_DOUBLE, 2*rank+1, 1, MPI_COMM_WORLD,&status);
                    MPI_Recv (&b[INDEX((1 + rank) * 4*s,j)], 1, MPI_DOUBLE, 2*rank+1, 2, MPI_COMM_WORLD,&status);
                    MPI_Recv (&c[INDEX((1 + rank) * 4*s,j)], 1, MPI_DOUBLE, 2*rank+1, 3, MPI_COMM_WORLD,&status);
                    MPI_Recv (&d[INDEX((1 + rank) * 4*s,j)], 1, MPI_DOUBLE, 2*rank+1, 4, MPI_COMM_WORLD,&status);
                    if (rank !=0 ){
                        for (int k=0;k<2;k++){
                            MPI_Recv(&a[INDEX((1+2*(rank+k)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+k), 1,
                                     MPI_COMM_WORLD, &status);
                            MPI_Recv(&b[INDEX((1+2*(rank+k)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+k), 2,
                                     MPI_COMM_WORLD, &status);
                            MPI_Recv(&c[INDEX((1+2*(rank+k)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+k), 3,
                                     MPI_COMM_WORLD, &status);
                            MPI_Recv(&d[INDEX((1+2*(rank+k)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+k), 4,
                                     MPI_COMM_WORLD, &status);
                        }
                    }
                    else{
                        MPI_Recv(&a[INDEX((1+2*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+1), 1,
                                 MPI_COMM_WORLD, &status);
                        MPI_Recv(&b[INDEX((1+2*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+1), 2,
                                 MPI_COMM_WORLD, &status);
                        MPI_Recv(&c[INDEX((1+2*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+1), 3,
                                 MPI_COMM_WORLD, &status);
                        MPI_Recv(&d[INDEX((1+2*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+1), 4,
                                 MPI_COMM_WORLD, &status);
                    }
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        
        
        if (rank==0)
            omega[1][middle+nghost-1] = d[INDEX(middle,m-1)] / b[INDEX(middle,m-1)];
        
        
        for (int j = m-1; j > 0; j--){
            
            int lenproc = (pow(2, m-j+1)-1)/nprocs+1;
            if (lenproc == 0) lenproc = 1;
            
            int s = (int) pow(2, j-1);
            
            if (lenproc==1){
                int l = (1 + rank) * 2*s;
                int k = (1 + rank) * s;
                if (rank != 0&& rank <(int)  pow(2, m-j)-1){
                    MPI_Send (&omega[1][l+nghost-1], 1, MPI_DOUBLE, 2*rank, 1, MPI_COMM_WORLD);
                    MPI_Send (&omega[1][l+nghost-1], 1, MPI_DOUBLE, 2*rank+1, 2, MPI_COMM_WORLD);
                    MPI_Send (&omega[1][l+nghost-1], 1, MPI_DOUBLE, 2*rank+2, 3, MPI_COMM_WORLD);
                }else if(rank==0 && rank <(int)  pow(2, m-j)-1){
                    MPI_Send (&omega[1][l+nghost-1], 1, MPI_DOUBLE, 2*rank+1, 2, MPI_COMM_WORLD);
                    MPI_Send (&omega[1][l+nghost-1], 1, MPI_DOUBLE, 2*rank+2, 3, MPI_COMM_WORLD);
                }
                if (rank != 0 && rank%2==0 && rank <(int)  pow(2, m-j+1)-2){
                    MPI_Recv(&omega[1][k-s+nghost-1], 1, MPI_DOUBLE, rank/2-1, 3,
                             MPI_COMM_WORLD, &status);
                    MPI_Recv(&omega[1][k+s+nghost-1], 1, MPI_DOUBLE, rank/2, 1,
                             MPI_COMM_WORLD, &status);
                }else if(rank%2==1 && rank <(int)  pow(2, m-j+1)-2){
                    MPI_Recv(&omega[1][k+nghost-1], 1, MPI_DOUBLE, rank/2, 2,
                             MPI_COMM_WORLD, &status);
                }else if(rank == (int)  pow(2, m-j+1)-2){
                    MPI_Recv(&omega[1][k-s+nghost-1], 1, MPI_DOUBLE, rank/2-1, 3,
                             MPI_COMM_WORLD, &status);
                }
                
                
                MPI_Barrier(MPI_COMM_WORLD);
                if (rank%2==0 && rank<(int) pow(2,m-j+1)-1){
                    if (k==s){
                        omega[1][k+nghost-1] = (d[INDEX(k,j-1)] - c[INDEX(k,j-1)] * omega[1][k+s+nghost-1])
                        / b[INDEX(k,j-1)];
                    }else if (k==nmesh+1-s){
                        omega[1][k+nghost-1] = (d[INDEX(k,j-1)] -
                                                a[INDEX(k,j-1)] * omega[1][k-s+nghost-1])
                        / b[INDEX(k,j-1)];
                    }else{
                        omega[1][k+nghost-1] = (d[INDEX(k,j-1)] -
                                                a[INDEX(k,j-1)] * omega[1][k-s+nghost-1] - c[INDEX(k,j-1)] * omega[1][k+s+nghost-1])
                        / b[INDEX(k,j-1)];
                    }
                }
            }else{
                if(rank!=nprocs-1){
                    MPI_Send (&omega[1][lenproc*(rank+1)*s+nghost-1], 1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD);
                }
                if(rank!=0){
                    MPI_Recv(&omega[1][lenproc*rank*s+nghost-1], 1, MPI_DOUBLE, rank-1, 1,
                             MPI_COMM_WORLD, &status);
                }
                MPI_Barrier(MPI_COMM_WORLD);
                for (int l = (1 + lenproc*rank) * s;
                     (l < (1 + lenproc*(rank + 1)) * s) &&
                     (l <= nmesh+1 - (int) pow (2,j-1));
                     l = l + 2 * s){
                    
                    if (l==s){
                        omega[1][l+nghost-1] = (d[INDEX(l,j-1)] - c[INDEX(l,j-1)] * omega[1][l+s+nghost-1])
                        / b[INDEX(l,j-1)];
                    }else if (l==nmesh+1-s){
                        omega[1][l+nghost-1] = (d[INDEX(l,j-1)] -
                                                a[INDEX(l,j-1)] * omega[1][l-s+nghost-1])
                        / b[INDEX(l,j-1)];
                    }else{
                        omega[1][l+nghost-1] = (d[INDEX(l,j-1)] -
                                                a[INDEX(l,j-1)] * omega[1][l-s+nghost-1] - c[INDEX(l,j-1)] * omega[1][l+s+nghost-1])
                        / b[INDEX(l,j-1)];
                    }
                }
            }
        }
        
        for (int j=0;j<nghost;j++){
            if (rank==0)
                omega[1][j] = omega[1][nghost+1];
            if (rank==nprocs-1)
                omega[1][lenx-1-j] = omega[1][lenx-1-nghost];
        }
        
        for (int j = nghost+lenproc*rank;j<nghost+lenproc*(rank+1) && j<nmesh+nghost;j++){
            double getalpha;
            if (x[j]<=xx)
                getalpha = 1.0;
            else if(x[j]<=xd)
                getalpha = pow(Calpha_dnull,(x[j]-xx)/(xd-xx));
            else
                getalpha = Calpha_dnull;
            
            double alpha_d = alpha_dt*getalpha;
            double vh_new = (dx[j]*v[1][j]+dx[j-1]*v[1][j-1])/(dx[j]+dx[j-1]);
            double vh_old = (dx[j]*v[0][j]+dx[j-1]*v[0][j-1])/(dx[j]+dx[j-1]);
            double vh = 0.5*(vh_new+vh_old);
            double Somega;
            double omega_mid = 0.5*(omega[1][j]+omega[0][j]);
            if (omega_mid+omegak[j]>omega0)
                Somega = 1.;
            else
                Somega = -1.;
            double gc = fabs(omega_mid+omegak[j]-omega0)/(omega0*h0*pow(x[j],0.25));
            double fc = gc/(1.+gc);
            double Pm = fc*vh/(alpha_d*h0*pow(x[j],0.25)*x[j]*fabs(omega_mid+omegak[j]-omega0));
            if (Pm<pmc){
                Pm = pmc;
            }
            double grho = 3.*(3.+Pm)*(1.-pow(x[j]/xc,3.));
            double frho = -3.*(3.+Pm)/(1.+grho);
            double SU = (1.+Somega)*(1.-vh/fabs(vh))/4.;
            wind[j] = -2.*fc*SU*frho*vh/(alpha_d*pow(x[j],6.5)*fabs(omega_mid+omegak[j]-omega0));
            //     wind[j] = 10.*wind[j];
            if (x[j]>xd)
                wind[j] = 0.0;
            B[j] = wind[j];
            if (rank != 0){
                MPI_Send (&wind[j], 1, MPI_DOUBLE, 0, j, MPI_COMM_WORLD);
            }else{
                for (int iproc = 1; iproc < nprocs; iproc++){
                    int k = j + lenproc*iproc;
                    if (k < nmesh+nghost){
                        MPI_Recv(&wind[k], 1, MPI_DOUBLE, iproc, k,
                                 MPI_COMM_WORLD, &status);
                    }
                }
            }
        }
        
        MPI_Barrier (MPI_COMM_WORLD);
        MPI_Bcast (wind, lenx, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Barrier (MPI_COMM_WORLD);
        
        for (int j = nghost+lenprocw*rank;j<nghost+lenprocw*(rank+1) && j<nwind+1;j++){
            B[j] = 0.0;
            if (j < nwidth + nghost){
                for (int k = nghost;k<nwidth+j;k++)
                    B[j] += exp(-pow((j-k)*dx[0]/deltax,2.)/2.)/integ[k]*wind[k]*dx[0];
            }else if (j>nmesh+nghost-nwidth-1){
                for (int k = j-nwidth+1;k<nmesh+nghost;k++)
                    B[j] += exp(-pow((j-k)*dx[0]/deltax,2.)/2.)/integ[k]*wind[k]*dx[0];
            }else{
                for (int k= j-nwidth+1;k<nwidth+j;k++)
                    B[j] += exp(-pow((j-k)*dx[0]/deltax,2.)/2.)/integ[k]*wind[k]*dx[0];
            }
            if (rank != 0){
                MPI_Send (&B[j], 1, MPI_DOUBLE, 0, j, MPI_COMM_WORLD);
            }else{
                for (int iproc = 1; iproc < nprocs; iproc++){
                    int k = j + lenprocw*iproc;
                    if (k < nwind+1){
                        MPI_Recv(&B[k], 1, MPI_DOUBLE, iproc, k,
                                 MPI_COMM_WORLD, &status);
                    }
                }
            }
        }
        MPI_Barrier (MPI_COMM_WORLD);
        MPI_Bcast (B, nwind+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Barrier (MPI_COMM_WORLD);
        
        for (int j = nghost+lenproc*rank;j<nghost+lenproc*(rank+1) && j<nmesh+nghost;j++){
            double a_advection;
            double b_advection;
            double c_advection;
            double d_advection;
            double r = dt/(dx[j-1]+dx[j])/2.0;
            double vr_new;
            double vl_new;
            double vr_old;
            double vl_old;
            if (v[0][j] <= 0.){
                vr_new = (dx[j+1]*v[1][j+1]+dx[j]*v[1][j])/(dx[j+1]+dx[j]);
                vl_new = (dx[j]*v[1][j]+dx[j-1]*v[1][j-1])/(dx[j]+dx[j-1]);
                vr_old = (dx[j+1]*v[0][j+1]+dx[j]*v[0][j])/(dx[j+1]+dx[j]);
                vl_old = (dx[j]*v[0][j]+dx[j-1]*v[0][j-1])/(dx[j]+dx[j-1]);
                if (j==nmesh+nghost-1){
                    a_advection = 0.0;
                    b_advection = x[j] - r*vl_new*x[j];
                    c_advection = 0.0;
                    d_advection = sigma[0][j]*x[j]+r*vl_old*x[j]*sigma[0][j]-r*vr_old*x[j+1]*sigma[0][j+1] - r*vr_new*x[j+1]*sigma[0][j+1];
                }else{
                    a_advection = 0.0;
                    b_advection = x[j] - r*vl_new*x[j];
                    c_advection = r*vr_new*x[j+1];
                    d_advection = sigma[0][j]*x[j]+r*vl_old*x[j]*sigma[0][j]-r*vr_old*x[j+1]*sigma[0][j+1];
                }
            }else{
                vr_new = (dx[j]*v[1][j]+dx[j-1]*v[1][j-1])/(dx[j]+dx[j-1]);
                vl_new = (dx[j-1]*v[1][j-1]+dx[j-2]*v[1][j-2])/(dx[j-1]+dx[j-2]);
                vr_old = (dx[j]*v[0][j]+dx[j-1]*v[0][j-1])/(dx[j]+dx[j-1]);
                vl_old = (dx[j-1]*v[0][j-1]+dx[j-2]*v[0][j-2])/(dx[j-1]+dx[j-2]);
                if (j==nghost){
                    a_advection = 0.0;
                    b_advection = x[j] + r*vr_new*x[j];
                    c_advection = 0.0;
                    d_advection = sigma[0][j]*x[j]+r*vl_old*x[j-1]*sigma[0][j-1]-r*vr_old*x[j]*sigma[0][j] + r*vl_new*x[j-1]*sigma[0][j-1];
                } else{
                    a_advection = -r*vl_new*x[j-1];
                    b_advection = x[j] + r*vr_new*x[j];
                    c_advection = 0.0;
                    d_advection = sigma[0][j]*x[j]+r*vl_old*x[j-1]*sigma[0][j-1]-r*vr_old*x[j]*sigma[0][j];
                }
            }
            d_flow[j] = B[j]*x[j]*dt;
            a[INDEX(j-nghost+1,0)] = a_advection;
            b[INDEX(j-nghost+1,0)] = b_advection;
            c[INDEX(j-nghost+1,0)] = c_advection;
            d[INDEX(j-nghost+1,0)] = d_advection + d_flow[j];
            
        }
        if (rank != 0){
            MPI_Send (&a[INDEX(1 + lenproc*rank,0)], 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
            MPI_Send (&b[INDEX(1 + lenproc*rank,0)], 1, MPI_DOUBLE, rank-1, 2, MPI_COMM_WORLD);
            MPI_Send (&c[INDEX(1 + lenproc*rank,0)], 1, MPI_DOUBLE, rank-1, 3, MPI_COMM_WORLD);
            MPI_Send (&d[INDEX(1 + lenproc*rank,0)], 1, MPI_DOUBLE, rank-1, 4, MPI_COMM_WORLD);
        }
        if (rank != nprocs-1){
            MPI_Recv(&a[INDEX(1 + lenproc*(rank+1),0)], 1, MPI_DOUBLE, rank+1, 1,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(&b[INDEX(1 + lenproc*(rank+1),0)], 1, MPI_DOUBLE, rank+1, 2,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(&c[INDEX(1 + lenproc*(rank+1),0)], 1, MPI_DOUBLE, rank+1, 3,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(&d[INDEX(1 + lenproc*(rank+1),0)], 1, MPI_DOUBLE, rank+1, 4,
                     MPI_COMM_WORLD, &status);
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        for (int j = 1; j < m; j++){
            
            int lenproc = (pow(2, m-j)-1)/ nprocs + 1;
            if (lenproc == 0) lenproc = 1;
            int s = (int) pow(2, j-1);
            
            for (int l = (1 + lenproc*rank)*2*s;
                 l < (1 + lenproc*(rank+1))*2*s
                 && l <= nmesh+1 - (int) pow(2, j);
                 l = l + 2*s){
                
                double e = - a[INDEX(l,j-1)] / b[INDEX(l-s,j-1)];
                double f = - c[INDEX(l,j-1)] / b[INDEX(l+s,j-1)];
                
                a[INDEX(l,j)] = e * a[INDEX(l-s,j-1)];
                c[INDEX(l,j)] = f * c[INDEX(l+s,j-1)];
                
                b[INDEX(l,j)] = b[INDEX(l,j-1)] +
                e * c[INDEX(l-s,j-1)] + f * a[INDEX(l+s,j-1)];
                
                d[INDEX(l,j)] = d[INDEX(l,j-1)] +
                e * d[INDEX(l-s,j-1)] + f * d[INDEX(l+s,j-1)];
            }
            if (lenproc != 1){
                if (rank != 0){
                    MPI_Send (&a[INDEX((1 + lenproc*rank) * 2*s,j)], 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
                    MPI_Send (&b[INDEX((1 + lenproc*rank) * 2*s,j)], 1, MPI_DOUBLE, rank-1, 2, MPI_COMM_WORLD);
                    MPI_Send (&c[INDEX((1 + lenproc*rank) * 2*s,j)], 1, MPI_DOUBLE, rank-1, 3, MPI_COMM_WORLD);
                    MPI_Send (&d[INDEX((1 + lenproc*rank) * 2*s,j)], 1, MPI_DOUBLE, rank-1, 4, MPI_COMM_WORLD);
                }
                if (rank != nprocs-1){
                    MPI_Recv(&a[INDEX((1 + lenproc*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, rank+1, 1,
                             MPI_COMM_WORLD, &status);
                    MPI_Recv(&b[INDEX((1 + lenproc*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, rank+1, 2,
                             MPI_COMM_WORLD, &status);
                    MPI_Recv(&c[INDEX((1 + lenproc*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, rank+1, 3,
                             MPI_COMM_WORLD, &status);
                    MPI_Recv(&d[INDEX((1 + lenproc*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, rank+1, 4,
                             MPI_COMM_WORLD, &status);
                }
            }else if(j<m-1){
                if (rank%2 == 1 && rank<pow(2, m-j)-1){
                    MPI_Send (&a[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2, 1, MPI_COMM_WORLD);
                    MPI_Send (&b[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2, 2, MPI_COMM_WORLD);
                    MPI_Send (&c[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2, 3, MPI_COMM_WORLD);
                    MPI_Send (&d[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2, 4, MPI_COMM_WORLD);
                }else if (rank != 0 && rank<pow(2, m-j)-2){
                    for (int k=0;k<2;k++){
                        MPI_Send (&a[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-k, 1, MPI_COMM_WORLD);
                        MPI_Send (&b[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-k, 2, MPI_COMM_WORLD);
                        MPI_Send (&c[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-k, 3, MPI_COMM_WORLD);
                        MPI_Send (&d[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-k, 4, MPI_COMM_WORLD);
                    }
                }else if (rank==pow(2, m-j)-2){
                    MPI_Send (&a[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-1, 1, MPI_COMM_WORLD);
                    MPI_Send (&b[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-1, 2, MPI_COMM_WORLD);
                    MPI_Send (&c[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-1, 3, MPI_COMM_WORLD);
                    MPI_Send (&d[INDEX((1 + rank) * 2*s,j)], 1, MPI_DOUBLE, rank/2-1, 4, MPI_COMM_WORLD);
                }
                
                if (rank<pow(2, m-j-1)-1){
                    MPI_Recv (&a[INDEX((1 + rank) * 4*s,j)], 1, MPI_DOUBLE, 2*rank+1, 1, MPI_COMM_WORLD,&status);
                    MPI_Recv (&b[INDEX((1 + rank) * 4*s,j)], 1, MPI_DOUBLE, 2*rank+1, 2, MPI_COMM_WORLD,&status);
                    MPI_Recv (&c[INDEX((1 + rank) * 4*s,j)], 1, MPI_DOUBLE, 2*rank+1, 3, MPI_COMM_WORLD,&status);
                    MPI_Recv (&d[INDEX((1 + rank) * 4*s,j)], 1, MPI_DOUBLE, 2*rank+1, 4, MPI_COMM_WORLD,&status);
                    if (rank !=0 ){
                        for (int k=0;k<2;k++){
                            MPI_Recv(&a[INDEX((1+2*(rank+k)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+k), 1,
                                     MPI_COMM_WORLD, &status);
                            MPI_Recv(&b[INDEX((1+2*(rank+k)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+k), 2,
                                     MPI_COMM_WORLD, &status);
                            MPI_Recv(&c[INDEX((1+2*(rank+k)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+k), 3,
                                     MPI_COMM_WORLD, &status);
                            MPI_Recv(&d[INDEX((1+2*(rank+k)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+k), 4,
                                     MPI_COMM_WORLD, &status);
                        }
                    }
                    else{
                        MPI_Recv(&a[INDEX((1+2*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+1), 1,
                                 MPI_COMM_WORLD, &status);
                        MPI_Recv(&b[INDEX((1+2*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+1), 2,
                                 MPI_COMM_WORLD, &status);
                        MPI_Recv(&c[INDEX((1+2*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+1), 3,
                                 MPI_COMM_WORLD, &status);
                        MPI_Recv(&d[INDEX((1+2*(rank+1)) * 2*s,j)], 1, MPI_DOUBLE, 2*(rank+1), 4,
                                 MPI_COMM_WORLD, &status);
                    }
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        
        
        if (rank==0)
            sigma[1][middle+nghost-1] = d[INDEX(middle,m-1)] / b[INDEX(middle,m-1)];
        
        for (int j = m-1; j > 0; j--){
            
            int lenproc = (pow(2, m-j+1)-1)/nprocs+1;
            if (lenproc == 0) lenproc = 1;
            
            int s = (int) pow(2, j-1);
            
            if (lenproc==1){
                int l = (1 + rank) * 2*s;
                int k = (1 + rank) * s;
                if (rank != 0&& rank <(int)  pow(2, m-j)-1){
                    MPI_Send (&sigma[1][l+nghost-1], 1, MPI_DOUBLE, 2*rank, 1, MPI_COMM_WORLD);
                    MPI_Send (&sigma[1][l+nghost-1], 1, MPI_DOUBLE, 2*rank+1, 2, MPI_COMM_WORLD);
                    MPI_Send (&sigma[1][l+nghost-1], 1, MPI_DOUBLE, 2*rank+2, 3, MPI_COMM_WORLD);
                }else if(rank==0 && rank <(int)  pow(2, m-j)-1){
                    MPI_Send (&sigma[1][l+nghost-1], 1, MPI_DOUBLE, 2*rank+1, 2, MPI_COMM_WORLD);
                    MPI_Send (&sigma[1][l+nghost-1], 1, MPI_DOUBLE, 2*rank+2, 3, MPI_COMM_WORLD);
                }
                if (rank != 0 && rank%2==0 && rank <(int)  pow(2, m-j+1)-2){
                    MPI_Recv(&sigma[1][k-s+nghost-1], 1, MPI_DOUBLE, rank/2-1, 3,
                             MPI_COMM_WORLD, &status);
                    MPI_Recv(&sigma[1][k+s+nghost-1], 1, MPI_DOUBLE, rank/2, 1,
                             MPI_COMM_WORLD, &status);
                }else if(rank%2==1 && rank <(int)  pow(2, m-j+1)-2){
                    MPI_Recv(&sigma[1][k+nghost-1], 1, MPI_DOUBLE, rank/2, 2,
                             MPI_COMM_WORLD, &status);
                }else if(rank == (int)  pow(2, m-j+1)-2){
                    MPI_Recv(&sigma[1][k-s+nghost-1], 1, MPI_DOUBLE, rank/2-1, 3,
                             MPI_COMM_WORLD, &status);
                }
                
                
                MPI_Barrier(MPI_COMM_WORLD);
                if (rank%2==0 && rank<(int) pow(2,m-j+1)-1){
                    if (k==s){
                        sigma[1][k+nghost-1] = (d[INDEX(k,j-1)] - c[INDEX(k,j-1)] * sigma[1][k+s+nghost-1])
                        / b[INDEX(k,j-1)];
                    }else if (k==nmesh+1-s){
                        sigma[1][k+nghost-1] = (d[INDEX(k,j-1)] -
                                                a[INDEX(k,j-1)] * sigma[1][k-s+nghost-1])
                        / b[INDEX(k,j-1)];
                    }else{
                        sigma[1][k+nghost-1] = (d[INDEX(k,j-1)] -
                                                a[INDEX(k,j-1)] * sigma[1][k-s+nghost-1] - c[INDEX(k,j-1)] * sigma[1][k+s+nghost-1])
                        / b[INDEX(k,j-1)];
                    }
                }
            }else{
                if(rank!=nprocs-1){
                    MPI_Send (&sigma[1][lenproc*(rank+1)*s+nghost-1], 1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD);
                }
                if(rank!=0){
                    MPI_Recv(&sigma[1][lenproc*rank*s+nghost-1], 1, MPI_DOUBLE, rank-1, 1,
                             MPI_COMM_WORLD, &status);
                }
                MPI_Barrier(MPI_COMM_WORLD);
                for (int l = (1 + lenproc*rank) * s;
                     (l < (1 + lenproc*(rank + 1)) * s) &&
                     (l <= nmesh+1 - (int) pow (2,j-1));
                     l = l + 2 * s){
                    
                    if (l==s){
                        sigma[1][l+nghost-1] = (d[INDEX(l,j-1)] - c[INDEX(l,j-1)] * sigma[1][l+s+nghost-1])
                        / b[INDEX(l,j-1)];
                    }else if (l==nmesh+1-s){
                        sigma[1][l+nghost-1] = (d[INDEX(l,j-1)] -
                                                a[INDEX(l,j-1)] * sigma[1][l-s+nghost-1])
                        / b[INDEX(l,j-1)];
                    }else{
                        sigma[1][l+nghost-1] = (d[INDEX(l,j-1)] -
                                                a[INDEX(l,j-1)] * sigma[1][l-s+nghost-1] - c[INDEX(l,j-1)] * sigma[1][l+s+nghost-1])
                        / b[INDEX(l,j-1)];
                    }
                }
            }
        }
        
        for (int j=0;j<nghost;j++){
            if (rank==0)
                sigma[1][j] = sigma[1][j+nghost];
            if (rank==nprocs-1)
                //sigma[1][lenx-1-j] = sigma0/x[lenx-1-j];
                sigma[1][lenx-1-j] = sigma[1][lenx-1-nghost]*x[lenx-1-nghost]/x[lenx-1-j];
        }
        
        lenproc = nmesh/nprocs + 1;
        if(rank!=0){
            MPI_Send (&omega[1][lenproc*rank+nghost], 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
            MPI_Send (&sigma[1][lenproc*rank+nghost], 1, MPI_DOUBLE, rank-1, 2, MPI_COMM_WORLD);
        }
        
        if(rank!=nprocs-1){
            MPI_Recv(&omega[1][lenproc*(rank+1)+nghost], 1, MPI_DOUBLE, rank+1, 1,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(&sigma[1][lenproc*(rank+1)+nghost], 1, MPI_DOUBLE, rank+1, 2,
                     MPI_COMM_WORLD, &status);
        }
        
        for (int j = lenproc*rank-nghost;j<lenproc*(rank+1)+1;j++){
            v[0][j+nghost] = v[1][j+nghost];
            omega[0][j+nghost] = omega[1][j+nghost];
            sigma[0][j+nghost] = sigma[1][j+nghost];
        }
        
        t = t + dt;
        MPI_Barrier (MPI_COMM_WORLD);
        if ((i+1)%dn == 0){
            for (int j = lenproc*rank;j<lenproc*(rank+1);j++){
                if (rank!=0){
                    MPI_Send (&v[0][j+nghost], 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
                    MPI_Send (&omega[0][j+nghost], 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
                    MPI_Send (&sigma[0][j+nghost], 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
                }else{
                    for (int iproc = 1; iproc < nprocs; iproc++){
                        int k = j + lenproc*iproc;
                        MPI_Recv(&v[0][k+nghost], 1, MPI_DOUBLE, iproc, 1,
                                 MPI_COMM_WORLD, &status);
                        MPI_Recv(&omega[0][k+nghost], 1, MPI_DOUBLE, iproc, 2,
                                 MPI_COMM_WORLD, &status);
                        MPI_Recv(&sigma[0][k+nghost], 1, MPI_DOUBLE, iproc, 3,
                                 MPI_COMM_WORLD, &status);
                    }
                }
                MPI_Barrier (MPI_COMM_WORLD);
            }
            if (rank==0){
                char fname[6];
                for (int j=0;j<5;j++){
                    fname[4-j] = (((i+1)/dn)%(int)pow(10,j+1))/(int)pow(10,j) + '0';
                }
                fname[5] = '\0';
                strcpy(address, dir);
                strcat(address, fname);
                strcat(address, ".txt");
                output = fopen(address, "w");
                fprintf(output, "######TIME %f\n",t);
                for (int j=0;j<lenx;j++){
                    fprintf(output, "%d %10.18f %10.18f %10.18f %10.18f %10.18f\n",j,v[0][j],omega[0][j],sigma[0][j],B[j],wind[j]);
                }
                fclose(output);
            }
        }
        MPI_Barrier (MPI_COMM_WORLD);
    }
    /*
     for (int j=0;j<lenx;j++){
     if(rank == 1)printf("%d %20.18f\n",j,sigma[1][j]);
     }*/
    
    wall_time = MPI_Wtime() - wall_time;
    
    if (rank == 0){
        printf ( "Running time = %f secs\n\n", wall_time );
        timestamp ();
        printf ("\n");
    }
    
    MPI_Finalize();
    
}


void timestamp (){
    static char time_buffer[40];
    const struct tm *tm;
    time_t now;
    
    now = time ( NULL );
    tm = localtime ( &now );
    
    strftime ( time_buffer, 40, "%d %B %Y %I:%M:%S %p", tm );
    
    printf ( "%s\n", time_buffer );
    
    return;
}



