#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <gsl/gsl_integration.h>

using namespace std;
#define pi 3.14159265358979323846
double y_muovert;
// double muovert1;
// double muovert2;

double rateIntegral(double qp, double qm, double v)
{
    double up = -log(qp);
    double um = -log(qm);
    double f, g, h; 
    double muovert1 = y_muovert; // y1=2 y2=-4 gives the right results from table
    double muovert2 = y_muovert;


    f = (pow(um, 3.0)*pow(1 - v, 2.0)*qp*qm*exp(-muovert1 - muovert2 + 2*v*um + 2*up))/
    (sqrt(up)*(exp(-muovert1 - muovert2) - qp*qp*qm*qm)*(expm1(-muovert1 - muovert2 + 2*v*um + 2*up)));
    //f = (pow(um, 3.0)*pow(1 - v, 2.0)*qp*qm)/(sqrt(up)*(exp(-muovert1 - muovert2) - qp*qp*qm*qm)*(-expm1(muovert1 + muovert2 - 2*v*um - 2*up)));


    //g = log(cosh(0.5*pow(pow(up, 0.5) + pow(um, 0.5), 2.0) - 0.5*muovert) / cosh(0.5*pow(pow(up, 0.5) - pow(um, 0.5), 2.0) - 0.5*muovert));
    g = log1p(exp(pow(sqrt(um) + sqrt(up), 2.0) - muovert1)) + log1p(exp(muovert2 - pow(sqrt(um) + sqrt(up), 2.0)))
    - log1p(exp(pow(sqrt(um) - sqrt(up), 2.0) - muovert1)) - log1p(exp(muovert2 - pow(sqrt(um) - sqrt(up), 2.0)));

    //h = log(cosh(0.5*pow(pow(v*um, 0.5) + pow(up, 0.5), 2.0) - 0.5*muovert) / cosh(0.5*pow(pow(v*um, 0.5) - pow(up, 0.5), 2.0) - 0.5*muovert));
    h = log1p(exp(pow(sqrt(v*um) + sqrt(up), 2.0) - muovert1)) + log1p(exp(muovert2 - pow(sqrt(v*um) + sqrt(up), 2.0)))
    - log1p(exp(pow(sqrt(v*um) - sqrt(up), 2.0) - muovert1)) - log1p(exp(muovert2 - pow(sqrt(v*um) - sqrt(up), 2.0)));//the term that takes the longest

    return f*g*h;
}

//3d integration routine
//************************************************************************************************************
static double xsav, ysav; //outer most integration variables 
static double (*nrfunc)(double, double, double);

double f3(double z, void *params)
{
    return (*nrfunc)(xsav,ysav,z); 
}

double f2(double y, void *params)
{
    double f3(double z, void *params);
    double z_min(double);
    double z_max(double);

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (5000);

    gsl_function integrand;
    integrand.function = &f3;

    double abs_error = 0.0;
    double rel_error = 1.0e-9;
    double result;
    double error;
    ysav = y;

    gsl_integration_qag (&integrand, 0.000, 1.0, abs_error, rel_error, 5000, 6, work_ptr, &result, &error);
    gsl_integration_workspace_free(work_ptr);

    return result;
}

double f1(double x, void *params) 
{
    double f2(double y, void *params);
    double y_min(double);
    double y_max(double);

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (5000);

    gsl_function integrand;
    integrand.function = &f2;

    double abs_error = 0.0;  
    double rel_error = 1.0e-9;  
    double result;  
    double error;
    xsav = x;

    gsl_integration_qag (&integrand, 0.000, 1.0, abs_error, rel_error, 5000, 6, work_ptr, &result, &error);
    gsl_integration_workspace_free(work_ptr);

    return result;
}

double quad3d(double (*func)(double, double, double), double x_min, double x_max) 
{
    double f1(double x, void *params);

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (5000);

    gsl_function integrand;
    integrand.function = &f1;

    double abs_error = 0.0; 
    double rel_error = 1.0e-9;  
    double result;  
    double error;

    nrfunc = func; 
    gsl_integration_qag (&integrand, x_min, x_max, abs_error, rel_error, 5000, 6, work_ptr, &result, &error);
    //read as: (integrand, lower bound, upper bound, absolute error, relative error, max # of subintervals, type of quadrature (see link below), workspace, result, error)
    //http://irtfweb.ifa.hawaii.edu/SoftwareDocs/gsl/gsl-ref-html/gsl-ref_16.html
    gsl_integration_workspace_free(work_ptr);
    return result;
}


//************************************************************************************************************
// Density Calculator 

double densityIntegral(double t, void *params)//integrand
{
    double y = ((double *)params)[0];
    return sqrt(-log(t))/(exp(-y)+t);
    //return y*t; //Solve 1=y*(1/2), so should give y=2
}

double integration_func(double muovert) 
{
    double x_min = 0.000;
    double x_max = 1.000;
    double params[] = { muovert };
    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (5000);

    gsl_function integrand;
    integrand.function = &densityIntegral;
    integrand.params = params;

    double abs_error = 0.0; 
    double rel_error = 1.0e-9;  
    double result;  
    double error;

    gsl_integration_qag (&integrand, x_min, x_max, abs_error, rel_error, 5000, 6, work_ptr, &result, &error);
    //read as: (integrand, lower bound, upper bound, absolute error, relative error, max # of subintervals, type of quadrature (see link below), workspace, result, error)
    //http://irtfweb.ifa.hawaii.edu/SoftwareDocs/gsl/gsl-ref-html/gsl-ref_16.html
    gsl_integration_workspace_free(work_ptr);
    return result;
}

//************************************************************************************************************


//Start of rootfinders

double Solve_Bisect(double nu, double(*func)(double), double x_min, double x_max, double tol, int *count)
{
    double x_mid, f_max, f_min, f_mid;
    double err; 
    int count_max;

    count_max = 200000; 
    *count += 1; //add 1 when BisectSolve is called

    x_mid = (x_min + x_max)/2.0; //calculate x_mid
    //printf("xmid=%e at iteration=%d \n", x_mid, *count);

        //Warn and exit
        if(*count > count_max)
        {
            printf("Solve_Bisect: Done %d iterations without convergence. \n", count_max);
            printf("Exiting. \n");
            exit(0); 
        }

    f_max = (*func)(x_max) - nu;
    f_min = (*func)(x_min) - nu; 

    if(f_max*f_min > 0.0) //we can't find a sol'n in this range
    {
            printf("No solution with x_min=%e and x_max=%e. Choose a different range. \n", x_min, x_max);
            printf("Exiting. \n");
        }

    f_mid = (*func)(x_mid) - nu; 

    //Calculates the error
    if(nu != 0.0) {err = fabsl(f_mid/nu);}
    else {err = fabsl(f_mid);}

    //if err<tol, we have a sol'n and the calculation ends.
    if(err < tol) {return x_mid;}

    if(f_max*f_mid < 0.0) //the soln is between x_max and x_mid
    {
        //call Solve_Bisect with the range (x_mid,x_max)
        return Solve_Bisect(nu, func, x_mid, x_max, tol, count);
    }
    else if(f_min*f_mid < 0.0) //the soln is between x_min and X_mid
    {
            //call Solve_Bisect with the range (x_min,x_mid)
        return Solve_Bisect(nu, func, x_min, x_mid, tol, count); //TD about recursion 
    }
    else //one of the factors is zero
    {
        if ((f_mid = 0.0)) 
            {return x_mid;} 
        else if ((f_max = 0.0)) 
            {return x_max;}
        else {return x_min;}
    }
}//BisectSolve

//************************************************************************************************************

int main (int argc, char **argv) 
{

    string Density = argv[1];
    string Mass = argv[2];
    double density = std::stod(Density);
    double mass = std::stod(Mass);
    double Temp = 6.0; //MeV 0.09 is for 10^9K
    double massconversion = 5.617977528e29; // convets kg to MeV
    double densityconversion = 197.3e-15*197.3e-15*197.3e-15;// converts 1/m^3 to MeV^3
    double tol = 1.0e-11;
    double axionrate;
    double f = 1.0; //Approx 1
    double mpi = 135.0; // MeV approx for neutral pion
    //double mu = pow(3*pi*pi*density*densityconversion, 1.0/3.0);
    //double g = 4.0e-6;

    FILE *output;
    string input_where_to_save = argv[3];
    string fname = input_where_to_save + "densityvsrate.dat";
    output = fopen(fname.c_str(), "w");

        int count = 0; //For Solve_Bisect

        double masstempfactor = (sqrt(2)/(pi*pi))*pow(mass*Temp*massconversion, 1.5);
        double lhs = density*densityconversion/masstempfactor;
        double SMsqd = (64*f*f*mass*mass*massconversion*massconversion)/(mpi*mpi*mpi*mpi);

        //Testing, uncomment here------------
        y_muovert = Solve_Bisect(lhs, integration_func, 92.0, 150.0, tol, &count);
        //y_muovert = (mu - mass*massconversion) / Temp; //For low temp
        axionrate = SMsqd*sqrt(mass*massconversion)*pow(Temp, 6.5)*quad3d(rateIntegral, 0.0, 1.0) / (pow(2, 5.5)*pow(pi, 7));
        //doesnt include axion coupling 
        

        printf("%le,%le,%le,%le\n", density*densityconversion, mass*massconversion, y_muovert, axionrate); //Main printf that I want

        fprintf(output, "Printing density (MeV), mass (MeV), y, 8.14511and R (MeV^5)\n");
        //Testing, uncomment here------------

    fclose(output);


        // fprintf(output1, "m^-3=%le kg=%le n=%le mass=%le\n lhs=%le massfac=%le y=%le R=%le\n\n",
        // density, mass, density*densityconversion, mass*massconversion, 
        // lhs, masstempfactor, y_muovert, axionrate);

    return 0;
}

