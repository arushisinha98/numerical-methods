//
//  main.cpp
//  Modeling Task 7
//
//  Created by Arushi Sinha on 21/5/19.
//  Copyright © 2019 Arushi Sinha. All rights reserved.
//
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <vector>
#include <string>

#include <math.h>

using namespace std;

// define and declare global variables
double L_x = 2000;
// domain size (x) = 2000 m
double L_z = 2000;
// domain size (y) = 2000 m
int nx = 401;
// 101 grid points in 0 <= x <= L_x
int nz = 401;
// 101 grid points in 0 <= y <= L_y
double dx = L_x/(nx-1);
// delta x
double dz = L_z/(nz-1);
// delta y
double r0 = 250;
// radius = 250 m
double dt = 0.5;
// delta t = 1s
double d_theta = 0.5;
// delta theta = 0.5 K
double theta_0 = 300;
// theta0 = 300 K
double g = 9.80665;
// gravitational acceleration = 9.80665 m/s^2
double T_total = 1500;
// T(total) = 1200 seconds
int nt = T_total/dt;
int save_frequency = 30;
// save full solution every 50 iterations

// (x01, y01)
double x0 = L_x/2;
double z0 = 260;

//  artificial diffusion/viscosity
double U = sqrt(2*r0*g*d_theta/theta_0);
double Re = 1500;
double K_x = (2*r0*U)/Re;
double K_z = K_x;

// formatting values
int width = 14; // width of fixed-length data fields
int precision = 6; // precision of data fields
// error bound values
double tolerance = 0.0000001;
double max_iterations = 2000;
// calculate alpha
double sigma = (1/(1+pow((dx/dz),2)))*(cos(M_PI/nx)+pow((dx/dz),2)*cos(M_PI/nz));
double alpha = 2/(1 + sqrt(1 - pow(sigma,2)));

// initialize potential temperature field
vector<vector<double>> initialize (vector<vector<double>>& u);
// prints all nodes
void print(vector<vector<double>> u, ofstream& output);
// integrates over all physical nodes (takes in field with ghost nodes)
double sum (vector<vector<double>> u);

double max_norm (vector<vector<double>> u);

double one_norm (vector<vector<double>> u);

double max (vector<vector<double>> u);

void Poisson (vector<vector<double>>& streamfunction, vector<vector<double>> f1, vector<vector<double>>& residual, double& R, double& V, double& S, double& epsilon);

vector<vector<double>> Jacobian (vector<vector<double>> f1,  vector<vector<double>> f2);

vector<vector<double>> EF_omega (vector<vector<double>> current, vector<vector<double>> J, vector<vector<double>> theta);

vector<vector<double>> EF_theta (vector<vector<double>> current, vector<vector<double>> J);

vector<vector<double>> AB3_omega (vector<vector<double>> current, vector<vector<double>> previous, vector<vector<double>> J, vector<vector<double>> J_1, vector<vector<double>> J_2, vector<vector<double>> T, vector<vector<double>> T_1, vector<vector<double>> T_2);

vector<vector<double>> AB3_theta (vector<vector<double>> current, vector<vector<double>> previous, vector<vector<double>> J, vector<vector<double>> J_1, vector<vector<double>> J_2);
// scheme used for diffusion for both omega and theta after time step
vector<vector<double>> FTCS (vector<vector<double>> u_this, vector<vector<double>> u);

double calculate_Cmax (vector<vector<double>> streamfunction);
// remove ghost nodes to calculate using physical node values
vector<vector<double>> remove_ghost_nodes (vector<vector<double>> u);
// add ghost nodes to calculate using five point stencil
vector<vector<double>> add_ghost_nodes (vector<vector<double>> u);

int main()
{
    cout << K_x << endl;
    vector<vector<double>> streamfunction(nz+2, vector<double>(nx+2));
    vector<vector<double>> residual(nz+2, vector<double>(nx+2));
    vector<vector<double>> omega(nz+2, vector<double>(nx+2));
    vector<vector<double>> theta(nz+2, vector<double>(nx+2));
    theta = initialize(theta);
    cout << "T: " << sum(theta) << endl;
    
    ofstream output1[2];
    ofstream output2[2];
    ofstream output3[2];
    ofstream outputA;
    
    outputA.open("Theta_0.csv", ifstream::out | ifstream::trunc);
    outputA.seekp(0, ios::beg);
    print(remove_ghost_nodes(theta),outputA);
    outputA.close();
    for (int j = 0; j < nz+2; j++)
    {
        for (int i = 0; i < nx+2; i++)
        {
            streamfunction[j][i] = 0.0;
            omega[j][i] = 0.0;
        }
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TIME STEP: n = 1
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Euler-Forward
    int loop = 1;
    // #1: calculate Jacobian (streamfunction, omega)
    vector<vector<double>> Jo_2 = Jacobian(streamfunction, omega);
    // #2: step omega in time
    // next_omega at time n+1
    // current_omega at time n
    // previous_omega at time n-1
    vector<vector<double>> previous_omega = omega;
    vector<vector<double>> current_omega = EF_omega(previous_omega, Jo_2, theta);
    
    /* output1[0].open("Omega_1.csv", ifstream::out | ifstream::trunc);
    output1[0].seekp(0, ios::beg);
    print(remove_ghost_nodes(current_omega),output1[0]);
    */
    // #3: calculate Jacobian (streamfunction, theta)
    vector<vector<double>> Jt_2 = Jacobian(streamfunction, theta);
    // #4: step theta in time
    // next_theta at time n+1
    // current_theta at time n
    // previous_theta at time n-1
    vector<vector<double>> previous_theta = theta;
    vector<vector<double>> current_theta = EF_theta(previous_theta, Jt_2);
    cout << "T: " << sum(current_theta) << endl;
    
    output2[0].open("Theta_1.csv", ifstream::out | ifstream::trunc);
    output2[0].seekp(0, ios::beg);
    print(remove_ghost_nodes(current_theta),output2[0]);
    
    // #5: update streamfunction
    // initialize streamfunction to 0 over domain
    // calculate residual for the initial epsilon
    for (int j = 2; j < nz; j++)
    {
        for (int i = 1; i < nx+1; i++)
        {
            residual[j][i] = (1.0/pow(dx,2))*(streamfunction[j][i-1] - 2.0*streamfunction[j][i] + streamfunction[j][i+1]) + (1.0/pow(dz,2))*(streamfunction[j-1][i] - 2.0*streamfunction[j][i] + streamfunction[j+1][i]) - current_omega[j][i];
        }
    }
    double R = max_norm(residual);
    double V = max_norm(current_omega);
    double S = one_norm(streamfunction);
    double epsilon = R/((2*((1/pow(dx,2))+(1/pow(dz,2)))*S) + V);
    int iterations = 0;
    while (epsilon > tolerance && iterations < max_iterations)
    {
        // run Poisson solver
        Poisson(streamfunction, current_omega, residual, R, V, S, epsilon);
        iterations++;
    }
    for (int j = 0; j < nz+2; j++)
    {
        for (int i = 0; i < nx+2; i++)
        {
            streamfunction[j][i] = (-1.0)*streamfunction[j][i];
        }
    }
    cout << loop << ": " << iterations << endl;
    
    /* output3[0].open("Streamfunction_1.csv", ifstream::out | ifstream::trunc);
    output3[0].seekp(0, ios::beg);
    print(remove_ghost_nodes(streamfunction),output3[0]);
    */
    outputA.open("Observations.csv", ifstream::out | ifstream::trunc);
    outputA.seekp(0, ios::beg);
    outputA << left << setw(width) << "n" << left << setw(width) << "C(max)" << left << setw(width) << "s(max)" <<  left << setw(width) << "T" << left << setw(width) << "Iterations" << endl;

    double Courant = calculate_Cmax(streamfunction);
    double Neumann = (K_x*2*dt/pow(dx,2))+(K_z*2*dt/pow(dz,2));
    double T = sum(current_theta);
    outputA << fixed << setprecision(precision) << left << setw(width) << loop << left << setw(width) << Courant << left << setw(width) << Neumann << left << setw(width) << T << left << setw(width) << iterations << endl;
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TIME STEP: n = 2
    ///////////////////////////////////////////////////////////////////////////////////////////////
    loop++;
    vector<vector<double>> Jo_1 = Jacobian(streamfunction, current_omega);
    vector<vector<double>> next_omega = EF_omega(current_omega, Jo_1, current_theta);
    
    /* output1[1].open("Omega_2.csv", ifstream::out | ifstream::trunc);
    output1[1].seekp(0, ios::beg);
    print(remove_ghost_nodes(next_omega),output1[1]);
     */

    vector<vector<double>> Jt_1 = Jacobian(streamfunction, current_theta);
    vector<vector<double>> next_theta = EF_theta(current_theta, Jt_1);
    cout << "T: " << sum(next_theta) << endl;
    
    output2[1].open("Theta_2.csv", ifstream::out | ifstream::trunc);
    output2[1].seekp(0, ios::beg);
    print(remove_ghost_nodes(next_theta),output2[1]);
    
    // initialize streamfunction to 0 over domain
    for (int j = 0; j < nz+2; j++)
    {
        for (int i = 0; i < nx+2; i++)
        {
            streamfunction[j][i] = 0.0;
        }
    }
    // calculate residual for epsilon
    for (int j = 2; j < nz; j++)
    {
        for (int i = 1; i < nx+1; i++)
        {
            residual[j][i] = (1/pow(dx,2))*(streamfunction[j][i-1] - 2*streamfunction[j][i] + streamfunction[j][i+1]) + (1/pow(dz,2))*(streamfunction[j-1][i] - 2*streamfunction[j][i] + streamfunction[j+1][i]) - next_omega[j][i];
        }
    }
    R = max_norm(residual);
    V = max_norm(next_omega);
    S = one_norm(streamfunction);
    epsilon = R/((2*((1/pow(dx,2))+(1/pow(dz,2)))*S) + V);
    iterations = 0;
    while (epsilon > tolerance && iterations < max_iterations)
    {
        // run Poisson solver
        Poisson(streamfunction, next_omega, residual, R, V, S, epsilon);
        iterations++;
    }
    for (int j = 0; j < nz+2; j++)
    {
        for (int i = 0; i < nx+2; i++)
        {
            streamfunction[j][i] = (-1.0)*streamfunction[j][i];
        }
    }
    cout << loop << ": " << iterations << endl;
    
    /* output3[1].open("Streamfunction_2.csv", ifstream::out | ifstream::trunc);
    output3[1].seekp(0, ios::beg);
    print(remove_ghost_nodes(streamfunction),output3[1]);
    */
    outputA << fixed << setprecision(precision) << left << setw(width) << loop << left << setw(width) << calculate_Cmax(streamfunction) << left << setw(width) << (K_x*dt/pow(dx,2))+(K_z*dt/pow(dz,2)) << left << setw(width) << setprecision(precision) << sum(next_theta) << left << setw(width) << iterations << endl;
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TIME STEP: n = 3  onwards
    ///////////////////////////////////////////////////////////////////////////////////////////////
    loop++;
    vector<vector<double>> Jo;
    vector<vector<double>> Jt;
    while (loop < nt)
    {
        Jo = Jacobian(streamfunction, next_omega);
        previous_omega = current_omega;
        current_omega = next_omega;
        next_omega = AB3_omega(current_omega, previous_omega, Jo, Jo_1, Jo_2, next_theta, current_theta, previous_theta);
        
        Jt = Jacobian(streamfunction, next_theta);
        previous_theta = current_theta;
        current_theta = next_theta;
        next_theta = AB3_theta(current_theta, previous_theta, Jt, Jt_1, Jt_2);
        cout << "T: " << sum(next_theta) << endl;
        
        Jo_2 = Jo_1;
        Jt_2 = Jt_1;
        Jo_1 = Jo;
        Jt_1 = Jt;
        
        for (int j = 0; j < nz+2; j++)
        {
            for (int i = 0; i < nx+2; i++)
            {
                streamfunction[j][i] = 0.0;
            }
        }
        for (int j = 2; j < nz; j++)
        {
            for (int i = 1; i < nx+1; i++)
            {
                residual[j][i] = (1/pow(dx,2))*(streamfunction[j][i-1] - 2*streamfunction[j][i] + streamfunction[j][i+1]) + (1/pow(dz,2))*(streamfunction[j-1][i] - 2*streamfunction[j][i] + streamfunction[j+1][i]) - next_omega[j][i];
            }
        }
        R = max_norm(residual);
        V = max_norm(next_omega);
        S = one_norm(streamfunction);
        epsilon = R/((2*((1/pow(dx,2))+(1/pow(dz,2)))*S) + V);
        iterations = 0;
        while (epsilon > tolerance && iterations < max_iterations)
        {
            // run Poisson solver
            Poisson(streamfunction, next_omega, residual, R, V, S, epsilon);
            iterations++;
        }
        for (int j = 0; j < nz+2; j++)
        {
            for (int i = 0; i < nx+2; i++)
            {
                streamfunction[j][i] = (-1.0)*streamfunction[j][i];
            }
        }
        cout << loop << ": " << iterations << endl;
        
        outputA << fixed << setprecision(precision) << left << setw(width) << loop << left << setw(width) << calculate_Cmax(streamfunction) << left << setw(width) << (K_x*dt/pow(dx,2))+(K_z*dt/pow(dz,2)) << left << setw(width) << setprecision(precision) << sum(next_theta) << left << setw(width) << iterations << endl;
        
        if (loop%save_frequency == 0)
        {
            ofstream output1;
            ofstream output2;
            ofstream output3;
            
            string header = "Omega_" + to_string(loop) + ".csv";
            /* output1.open(header, ifstream::out | ifstream::trunc);
            output1.seekp(0, ios::beg);
            print(remove_ghost_nodes(next_omega),output1);
            header.clear(); */
            
            header.clear();
            header = "Theta_" + to_string(loop) + ".csv";
            output2.open(header, ifstream::out | ifstream::trunc);
            output2.seekp(0, ios::beg);
            print(remove_ghost_nodes(next_theta),output2);
            
            /* header = "Streamfunction_" + to_string(loop) + ".csv";
            output3.open(header, ifstream::out | ifstream::trunc);
            output3.seekp(0, ios::beg);
            print(remove_ghost_nodes(streamfunction),output3);
            header.clear(); */
        }
        
        loop++;
    }
    return 0;
}

vector<vector<double>> initialize (vector<vector<double>>& u)
{
    // initial condition
    for (int col = 0; col < nz+2; col++)
    {
        for (int row = 0; row < nx+2; row++)
        {
            u[col][row] = theta_0;
            double x = static_cast<double>(dx*(row-1));
            double z = static_cast<double>(dz*(col-1));
            if (sqrt(pow((x - x0),2) + pow((z - z0),2)) <= r0)
            {
                u[col][row] += d_theta;
            }
        }
        // boundary conditions:
        // wall on x
        // periodic on x
        u[col][0] = u[col][nx-1];
        u[col][nx+1] = u[col][2];
    }
    return u;
}

void Poisson (vector<vector<double>>& streamfunction, vector<vector<double>> omega, vector<vector<double>>& residual, double& R, double& V, double& S, double& epsilon)
{
    for (int j = 2; j < nz; j++)
    {
        if (j == nz) // if periodic
        {
            streamfunction[nz+1] = streamfunction[2]; // update ghost at y = nz+1
        }
        for (int i = 1; i < nx+1; i++)
        {
            if (i == nx)
            {
                streamfunction[j][nx+1] = streamfunction[j][2]; // update ghost at x = nx+1
            }
            // calculate residual for interior nodes
            residual[j][i] = (1/pow(dx,2))*(streamfunction[j][i-1] - 2*streamfunction[j][i] + streamfunction[j][i+1]) + (1/pow(dz,2))*(streamfunction[j-1][i] - 2*streamfunction[j][i] + streamfunction[j+1][i]) - omega[j][i];
            // calculate streamfunction
            streamfunction[j][i] = streamfunction[j][i] + alpha*residual[j][i]/(2.0*((1/(pow(dx,2)))+(1/pow(dz,2))));
        }
        streamfunction[j][0] = streamfunction[j][nx-1]; // update ghost at x = 0
        streamfunction[j][1] = streamfunction[j][nx]; // enforce (x = 1) = (x = nx)
    }
    // streamfunction[0] = streamfunction[nz-1]; // update ghost at y = 0
    // streamfunction[1] = streamfunction[nz]; // enforce (y = 1) = (y = nz)
    
    for (int index = 0; index < nx+2; index++)
    {
        streamfunction[0][index] = 0.0;
        streamfunction[1][index] = 0.0;
        streamfunction[nz][index] = 0.0;
        streamfunction[nz+1][index] = 0.0;
    }
    
    R = max_norm(remove_ghost_nodes(residual));
    V = max_norm(remove_ghost_nodes(omega));
    S = one_norm(remove_ghost_nodes(streamfunction));
    epsilon = R/((2*((1/pow(dx,2))+(1/pow(dz,2)))*S) + V);
}

vector<vector<double>> Jacobian (vector<vector<double>> f1, vector<vector<double>> f2)
{
    double J1, J2, J3;
    vector<vector<double>> J(nz+2, vector<double>(nx+2));
    double A, B;
    for (int j = 2; j < nz; j++)
    {
        for (int i = 2; i < nx; i++)
        {
            J1 = (1.0/(2*dx))*(f1[j][i+1] - f1[j][i-1])*(1.0/(2*dz))*(f2[j+1][i] - f2[j-1][i]) - (1.0/(2*dz))*(f1[j+1][i] - f1[j-1][i])*(1.0/(2*dx))*(f2[j][i+1] - f2[j][i-1]);
            
            A = (1.0/(2*dx))*(f1[j][i+1]*(1.0/(2*dz))*(f2[j+1][i+1] - f2[j-1][i+1]) - f1[j][i-1]*(1.0/(2*dz))*(f2[j+1][i-1] - f2[j-1][i-1]));
            B = (1.0/(2*dz))*(f1[j+1][i]*(1.0/(2*dx))*(f2[j+1][i+1] - f2[j+1][i-1]) - f1[j-1][i]*(1.0/(2*dx))*(f2[j-1][i+1] - f2[j-1][i-1]));
            J2 = A - B;
            
            A = (1.0/(2*dz))*(f2[j+1][i]*(1.0/(2*dx))*(f1[j+1][i+1] - f1[j+1][i-1]) - f2[j-1][i]*(1.0/(2*dx))*(f1[j-1][i+1] - f1[j-1][i-1]));
            B = (1.0/(2*dx))*(f2[j][i+1]*(1.0/(2*dz))*(f1[j+1][i+1] - f1[j-1][i+1]) - f2[j][i-1]*(1.0/(2*dz))*(f1[j+1][i-1] - f1[j-1][i-1]));
            J3 = A - B;
            
            J[j][i] = (1.0/3.0)*(J1 + J2 + J3);
        }
    }
    // implement boundary counditions
    /* for (int index = 0; index < nz+2; index++)
    {
        J[index][0] = J[index][nx-1];
        J[index][nx+1] = J[index][2];
        J[index][1] = J[index][nx];
    } */
    return J;
}

vector<vector<double>> EF_omega (vector<vector<double>> current, vector<vector<double>> J, vector<vector<double>> theta)
{
    // calculate Jacobian (J) at current time-step (n)
    vector<vector<double>> next(nz+2, vector<double>(nx+2));
    for (int j = 1; j < nz+1; j++)
    {
        for (int i = 1; i < nx+1; i++)
        {
            next[j][i] = current[j][i] - dt*(J[j][i]) - (g/theta_0)*dt*(theta[j][i+1] - theta[j][i-1])/(2*dx) + K_x*(dt/pow(dx,2))*(current[j][i+1] - 2*current[j][i] + current[j][i-1]) + K_z*(dt/pow(dz,2))*(current[j+1][i] - 2*current[j][i] + current[j-1][i]);
        }
    }
    for (int index = 0; index < nz+2; index++)
    {
        // periodic boundary on x
        next[index][0] = next[index][nx-1];
        next[index][nx+1] = next[index][2];
        next[index][1] = next[index][nx]; // enforce (x = 1) = (x = nx)
    }
    for (int index = 0; index < nz+2; index++)
    {
        // wall on z
        next[0][index] = 0.0;
        next[1][index] = 0.0;
        next[nz][index] = 0.0;
        next[nz+1][index] = 0.0;
    }
    return next;
}

vector<vector<double>> EF_theta (vector<vector<double>> current, vector<vector<double>> J)
{
    // calculate Jacobian (J) at current time-step (n)
    // calculate Streamfunction (SF) at current time-step (n)
    vector<vector<double>> next(nz+2, vector<double>(nx+2));
    for (int j = 1; j < nz+1; j++)
    {
        for (int i = 1; i < nx+1; i++)
        {
            next[j][i] = current[j][i] - dt*(J[j][i]) + K_x*(dt/pow(dx,2))*(current[j][i+1] - 2*current[j][i] + current[j][i-1]) + K_z*(dt/pow(dz,2))*(current[j+1][i] - 2*current[j][i] + current[j-1][i]);
        }
    }
    for (int index = 0; index < nz+2; index++)
    {
        // periodic boundary on x
        next[index][0] = next[index][nx-1];
        next[index][nx+1] = next[index][2];
        next[index][1] = next[index][nx]; // enforce (x = 1) = (x = nx)
    }
    for (int index = 0; index < nz+2; index++)
    {
        // wall on z
        next[0][index] = theta_0;
        next[1][index] = theta_0;
        next[nz][index] = theta_0;
        next[nz+1][index] = theta_0;
    }
    return next;
}

vector<vector<double>> AB3_omega (vector<vector<double>> current, vector<vector<double>> previous, vector<vector<double>> J, vector<vector<double>> J_1, vector<vector<double>> J_2, vector<vector<double>> T, vector<vector<double>> T_1, vector<vector<double>> T_2)
{
    // calculate Jacobian (J) and theta (T) at current time-step (n)
    // calculate Jacobian (J_1) and theta (T_1) at n-1
    // calculate Jacobian (J_2) and theta (T_2) at n-2
    vector<vector<double>> next(nz+2, vector<double>(nx+2));
    for (int j = 2; j < nz; j++)
    {
        for (int i = 1; i < nx+1; i++)
        {
            next[j][i] = current[j][i] - (23.0/12.0)*((J[j][i])*dt + (g/theta_0)*dt*(T[j][i+1] - T[j][i-1])/(2*dx)) + (16.0/12.0)*((J_1[j][i])*dt + (g/theta_0)*dt*(T_1[j][i+1] - T_1[j][i-1])/(2*dx)) - (5.0/12.0)*((J_2[j][i])*dt + (g/theta_0)*dt*(T_2[j][i+1] - T_2[j][i-1])/(2*dx)) + K_x*(2*dt/pow(dx,2))*(previous[j][i+1] - 2*previous[j][i] + previous[j][i-1]) + K_z*(2*dt/pow(dz,2))*(previous[j+1][i] - 2*previous[j][i] + previous[j-1][i]);
        }
    }
    for (int index = 0; index < nz+2; index++)
    {
        // periodic boundary on x
        next[index][0] = next[index][nx-1];
        next[index][nx+1] = next[index][2];
        next[index][1] = next[index][nx]; // enforce (x = 1) = (x = nx)
    }
    for (int index = 0; index < nz+2; index++)
    {
        // wall on z
        next[0][index] = 0.0;
        next[1][index] = 0.0;
        next[nz][index] = 0.0;
        next[nz+1][index] = 0.0;
    }
    return next;
}

vector<vector<double>> AB3_theta (vector<vector<double>> current, vector<vector<double>> previous, vector<vector<double>> J, vector<vector<double>> J_1, vector<vector<double>> J_2)
{
    // calculate Jacobian (J) at current time-step (n)
    // calculate Jacobian (J_1) at n-1
    // calculate Jacobian (J_2) at n-2
    vector<vector<double>> next(nz+2, vector<double>(nx+2));
    for (int j = 2; j < nz; j++)
    {
        for (int i = 1; i < nx+1; i++)
        {
            next[j][i] = current[j][i] - (23.0/12.0)*(J[j][i])*dt + (16.0/12.0)*(J_1[j][i])*dt - (5.0/12.0)*(J_2[j][i])*dt + K_x*(2*dt/pow(dx,2))*(previous[j][i+1] - 2*previous[j][i] + previous[j][i-1]) + K_z*(2*dt/pow(dz,2))*(previous[j+1][i] - 2*previous[j][i] + previous[j-1][i]);
        }
    }
    for (int index = 0; index < nz+2; index++)
    {
        // periodic boundary on x
        next[index][0] = next[index][nx-1];
        next[index][nx+1] = next[index][2];
        next[index][1] = next[index][nx]; // enforce (x = 1) = (x = nx)
    }
    for (int index = 0; index < nz+2; index++)
    {
        // wall on z
        next[0][index] = theta_0;
        next[1][index] = theta_0;
        next[nz][index] = theta_0;
        next[nz+1][index] = theta_0;
    }
    return next;
}

vector<vector<double>> FTCS (vector<vector<double>> u_current, vector<vector<double>> u_previous)
{
    vector<vector<double>> u_next(nz+2, vector<double>(nx+2));
    for (int j = 2; j < nz; j++)
    {
        for (int i = 1; i < nx+1; i++)
        {
            u_next[j][i] = u_current[j][i] + K_x*(2*dt/pow(dx,2))*(u_previous[j][i+1] - 2*u_previous[j][i] + u_previous[j][i-1]) + K_z*(2*dt/pow(dz,2))*(u_previous[j+1][i] - 2*u_previous[j][i] + u_previous[j-1][i]);
        }
    }
    for (int index = 0; index < nz+2; index++)
    {
        // periodic boundary on x
        u_next[index][0] = u_next[index][nx-1];
        u_next[index][nx+1] = u_next[index][2];
        u_next[index][1] = u_next[index][nx]; // enforce (x = 1) = (x = nx)
    }
    return u_next;
}

double calculate_Cmax (vector<vector<double>> streamfunction)
{
    vector<vector<double>> u(nz+2, vector<double>(nx+2));
    vector<vector<double>> v(nz+2, vector<double>(nx+2));
    vector<double> C;
    for (int j = 2; j < nz; j++)
    {
        for (int i = 1; i < nx+1; i++)
        {
            u[j][i] = (streamfunction[j-1][i] - streamfunction[j+1][i])/(2*dz);
            v[j][i] = (streamfunction[j][i+1] - streamfunction[j][i-1])/(2*dx);
            C.push_back(abs(u[j][i])*(dt/dx) + abs(v[j][i])*(dt/dz));
        }
    }
    return *max_element(C.begin(),C.end());
}

double sum (vector<vector<double>> u)
{
    double sum = 0;
    for (int col = 2; col < u.size()-2; col++)
    {
        for (int row = 1; row < u[0].size()-1; row++)
        {
            sum += u[col][row];
        }
    }
    return sum;
}

double max_norm (vector<vector<double>> u)
{
    double max = abs(u[1][1]);
    for (int col = 2; col < u.size()-2; col++)
    {
        for (int row = 1; row < u[0].size()-1; row++)
        {
            if (abs(u[col][row]) > max)
            {
                max = abs(u[col][row]);
            }
        }
    }
    return max;
}

double one_norm (vector<vector<double>> u)
{
    double max = 0;
    for (int col = 2; col < u.size()-2; col++)
    {
        for (int row = 1; row < u[0].size()-1; row++)
        {
            max += abs(u[col][row]);
        }
    }
    return max;
}

double max (vector<vector<double>> u)
{
    double max = u[1][1];
    for (int col = 2; col < u.size()-2; col++)
    {
        for (int row = 1; row < u[0].size()-1; row++)
        {
            if (u[col][row] > max)
            {
                max = u[col][row];
            }
        }
    }
    return max;
}

vector<vector<double>> add_ghost_nodes (vector<vector<double>> u)
{
    vector<vector<double>> ghost(u.size()+2, vector<double> (u.size()+2));
    for (int j = 1; j < u.size()-1; j++)
    {
        for (int i = 1; i < u.size()-1; i++)
        {
            ghost[j+1][i+1] = u[j][i];
        }
    }
    // ghost nodes: L_y + 1 = 1
    ghost[u.size()+1] = ghost[2];
    // ghost nodes: -1 = L_y - 1
    ghost[0] = ghost[u.size()-1];
    for (int j = 0; j < u.size()+2; j++)
    {
        // ghost nodes: L_x + 1 = 1
        ghost[j][u[0].size()+1] = ghost[j][2];
        // ghost nodes: -1 = L_y - 1
        ghost[j][0] = ghost[j][u[0].size()-1];
    }
    return ghost;
}

vector<vector<double>> remove_ghost_nodes (vector<vector<double>> u)
{
    vector<vector<double>> original(u.size()-2, vector<double>(u[0].size()-2));
    for (int j = 1; j < u.size()-1; j++)
    {
        for (int i = 1; i < u[0].size()-1; i++)
        {
            original[j-1][i-1] = u[j][i];
        }
    }
    return original;
}

void print(vector<vector<double>> u, ofstream& output)
{
    // prints all nodes in domain (assuming there are no ghost nodes)
    for (int j = 0; j < nz; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            output << fixed << setprecision(precision) << left << setw(width) << u[j][i];
        }
        output << endl;
    }
}
