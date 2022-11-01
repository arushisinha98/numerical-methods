//
//  main.cpp
//  Modeling Task 2
//
//  Created by Arushi Sinha on 13/4/19.
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
double T_total = 200;
// total duration of observation (final time step)
double A = 1;
// amplitude
double lambda = 10;
// wavelength
// EXPERIMENT #1, #2, #3, #5: 50
// EXPERIMENT #4: 10
double c = 0.5;
// propagation velocity = 0.5 m/s
double dt = 2.0;
// delta t
// EXPERIMENT #1: 2
// EXPERIMENT #2: 0.25, 0.50, 0.75, 1.00
// EXPERIMENT #3: 4, 5, 8
// EXPERIMENT #4: 1, 2
// EXPERIMENT #5: 1
double L = 50;
// domain size = 50 m
int nx = 51;
// 51 grid points in 0 <= x <= L
double dx = L/(nx-1);
// delta x
double nt = T_total/dt;
// number of timesteps
double Courant = c*dt/dx;
// Courant's number
int save_frequency = 2;
// save full solution every 5 iterations
// save full solution every 2 iterations (EXPERIMENT #3 & #4 ONLY)

// formatting values
int width = 12; // width of fixed-length data fields
int precision = 6; // precision of data fields

// initialize u(x,t = 0)
vector<double> exact(vector<double>& u, double time);

// finding the error
// error(t) = max{|u_exact(x_i,t) - u_approx(x_i,t)|}
double find_error(vector<double> u, double time);

// finding the amplitude
// amplitude(t) = max{|u_approx(x_i,t)|}
double find_amplitude(vector<double> u);

// implement sign function for EXPERIMENT #5
int sign(double value);

int main()
{
    vector<double> u_exact(nx);
    u_exact = exact(u_exact, 0);
    double amplitude;
    double error;
    // vector of vector of full solutions to be saved and printed in file: "full_solution.txt"
    vector<vector<double>> save_solutions;
    vector<vector<double>> exact_solutions;
    ofstream output1;
    ofstream output2;
    ofstream output3;
    output1.open("summary_statistics.csv", ifstream::out | ifstream::trunc);
    output1.seekp(0, ios::beg);
    output1 << left << setw(width) << "Time" << left << setw(width) << "Amplitude" << left << setw(width) << "Error" << endl;
    
    // initialize vector u(x = i, t = n+1)
    vector<double> u_next(nx);
    // initialize vector u(x = i, t = n)
    vector<double> u_current = u_exact;

    for (int n = 0; n <= nt; n++)
    {
        // output time, amplitude, error
        u_exact = exact(u_exact, n*dt);
        amplitude = find_amplitude(u_current);
        error = find_error(u_current, n*dt);
        output1.seekp(0, ios::end);
        output1 << left << setw(width) << n*dt << left << setw(width) << amplitude << left << setw(width) << error << endl;
        // advance the solution for interior nodes
        for (int i = 1; i < nx-1; i++)
        {
            u_next[i] = u_current[i] - Courant*(u_current[i]-u_current[i-1]);
        }
        // handle boundary conditions
        // right boundary condition (backward in space)
        u_next[nx-1] = u_current[nx-1] - Courant*(u_current[nx-1]-u_current[nx-2]);
        // left boundary condition (cyclic boundary)
        u_next[0] = u_next[nx-1];
        
        if (n%save_frequency == 0)
        {
            save_solutions.push_back(u_current);
            exact_solutions.push_back(u_exact);
        }
        // update solution
        u_current = u_next;
    }
    output1.close();
    
    // output full solutions
    vector<string> header(nx);
    output2.open("full_solution.csv", ifstream::out | ifstream::trunc);
    output2.seekp(0, ios::beg);
    output3.open("exact_solution.csv", ifstream::out | ifstream::trunc);
    output3.seekp(0, ios::beg);
    for (int i = 0; i < nx; i++)
    {
        header[i] = "u(" + to_string(i) + ")";
        output2 << left << setw(width) << header[i];
        output3 << left << setw(width) << header[i];
        for (int j = 0; j < save_solutions.size(); j++)
        {
            output2.seekp(0, ios::end);
            output2 << left << setw(width) << save_solutions[j][i];
            output3.seekp(0, ios::end);
            output3 << left << setw(width) << exact_solutions[j][i];
        }
        output2.seekp(0, ios::end);
        output2 << endl;
        output3.seekp(0, ios::end);
        output3 << endl;
    }
    output2.close();
    output3.close();
    return 0;
}

vector<double> exact(vector<double>& u, double time)
{
    // initial condition
    // u(x,t = 0) = u0(x) = A sin(2*pi*x/lambda);
    vector<double> x(u.size());
    for (int i = 0; i < u.size(); i++)
    {
        x[i] = i*dx;
        u[i] = A*sin(2*M_PI*(x[i] - c*time)/lambda);
        // u[i] = sign(A*sin(2*M_PI*(x[i] - c*time)/lambda));
        // FOR EXPERIMENT #5 ONLY
    }
    return u;
}

double find_error(vector<double> u, double time)
{
    vector<double> difference(u.size());
    double u_exact;
    double x;
    for (int i = 0; i < u.size(); i++)
    {
        x = i*dx;
        u_exact = A*sin(2*M_PI*(x - c*time)/lambda);
        // u_exact = sign(A*sin(2*M_PI*(x - c*time)/lambda));
        // FOR EXPERIMENT #5 ONLY
        difference[i] = abs(u_exact - u[i]);
    }
    double error = *max_element(difference.begin(),difference.end());
    return error;
}

double find_amplitude(vector<double> u)
{
    vector<double> height(u.size());
    for (int i = 0; i < u.size(); i++)
    {
        height[i] = (abs(u[i]));
    }
    double amplitude = *max_element(height.begin(),height.end());
    return amplitude;
}

int sign(double value)
{
    if (value < 0)
        return -1;
    else if (value == 0)
        return 0;
    else
        return 1;
}
