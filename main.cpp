//
//  main.cpp
//  Modeling Task 3
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
double T_total = 2000;
// total duration of observation (final time step)
double A = 1;
// amplitude
double lambda = 50;
// wavelength
// EXPERIMENT #1, #2: 50
// EXPERIMENT #3: 10
// EXPERIMENT #4: 25
double c = 0.5;
// propagation velocity = 0.5 m/s
double dt = 2.0;
// delta t
// EXPERIMENT #1, #3: 2
// EXPERIMENT #2: 1.00, 2.50
// EXPERIMENT #4: 1
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
int save_frequency = 5;
// save full solution every 5 iterations

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

// implement sign function for square wave
int sign(double value);

// Forward-in-Time, Backward-in-Space (FTBS)
vector<double> FTBS (vector<double> u);

// Runge-Kutta 3
vector<double> RK3 (vector<double> u);

// Leap-frog
vector<double> leapfrog (vector<double> u_current, vector<double> u_previous);

// Adam's-Bashforth 3
vector<double> AB3 (vector<double> u_current, vector<double> u_previous1, vector<double> u_previous2);

int main()
{
    vector<double> u_exact(nx);
    u_exact = exact(u_exact, 0);
    double amplitude;
    double error;
    // vector of vector of full solutions to be saved and printed in file: "full_solution.txt"
    vector<vector<double>> save_solutions1; // for FTBS
    vector<vector<double>> save_solutions2; // for RK3
    vector<vector<double>> exact_solutions;
    ofstream output1[4]; // files for summary statistics
    ofstream output2[4]; // files for full solutions
    ofstream output3; // file for exact solution
    // output file for FTBS
    output1[0].open("FTBS_summary.csv", ifstream::out | ifstream::trunc);
    output1[0].seekp(0, ios::beg);
    // output file for RK3
    output1[1].open("RK3_summary.csv", ifstream::out | ifstream::trunc);
    output1[1].seekp(0, ios::beg);
    // output file for Leap-Frog
    output1[2].open("LF_summary.csv", ifstream::out | ifstream::trunc);
    output1[2].seekp(0, ios::beg);
    // output file for AB3
    output1[3].open("AB3_summary.csv", ifstream::out | ifstream::trunc);
    output1[3].seekp(0, ios::beg);
    for (int file = 0; file < 3; file++)
    {
        output1[file] << left << setw(width) << "Time" << left << setw(width) << "Amplitude" << left << setw(width) << "Error" << endl;
    }
    // initialize vector u(x = i, t = n+1)
    vector<double> u_next1(nx); // u_next for FTBS
    vector<double> u_next2(nx); // u_next for RK3
    // initialize vector u(x = i, t = n)
    vector<double> u_current1 = u_exact; // u_current for FTBS
    vector<double> u_current2 = u_exact; // u_current for RK3
    
    // loop for FTBS and RK3 (do not require saving previous timestep)
    for (int n = 0; n <= nt; n++)
    {
        // output time, amplitude, error for FTBS
        u_exact = exact(u_exact, n*dt);
        amplitude = find_amplitude(u_current1);
        error = find_error(u_current1, n*dt);
        output1[0].seekp(0, ios::end);
        output1[0] << left << setw(width) << n*dt << left << setw(width) << amplitude << left << setw(width) << error << endl;
        // output time, amplitude, error for RK3
        amplitude = find_amplitude(u_current2);
        error = find_error(u_current2, n*dt);
        output1[1].seekp(0, ios::end);
        output1[1] << left << setw(width) << n*dt << left << setw(width) << amplitude << left << setw(width) << error << endl;
        // advance u_next1 w/ FTBS
        u_next1 = FTBS(u_current1);
        // advance u_next2 w/ RK3
        u_next2 = RK3(u_current2);
        
        if (n%save_frequency == 0)
        {
            save_solutions1.push_back(u_current1); // save FTBS iteration
            save_solutions2.push_back(u_current2); // save RK3 iteration
            exact_solutions.push_back(u_exact); // save exact solution
        }
        // update solution
        u_current1 = u_next1; // update FTBS
        u_current2 = u_next2; // update RK3
    }
    output1[0].close();
    output1[1].close();
    
    // output full solutions
    vector<string> header(nx);
    output2[0].open("FTBS_full.csv", ifstream::out | ifstream::trunc);
    output2[0].seekp(0, ios::beg);
    output2[1].open("RK3_full.csv", ifstream::out | ifstream::trunc);
    output2[1].seekp(0, ios::beg);
    output3.open("exact_solution.csv", ifstream::out | ifstream::trunc);
    output3.seekp(0, ios::beg);
    for (int i = 0; i < nx; i++)
    {
        header[i] = "u(" + to_string(i) + ")";
        output2[0] << left << setw(width) << header[i];
        output2[1] << left << setw(width) << header[i];
        output3 << left << setw(width) << header[i];
        for (int j = 0; j < save_solutions1.size(); j++)
        {
            output2[0].seekp(0, ios::end);
            output2[0] << left << setw(width) << save_solutions1[j][i];
            output2[1].seekp(0, ios::end);
            output2[1] << left << setw(width) << save_solutions2[j][i];
            output3.seekp(0, ios::end);
            output3 << left << setw(width) << exact_solutions[j][i];
        }
        output2[0].seekp(0, ios::end);
        output2[0] << endl;
        output2[1].seekp(0, ios::end);
        output2[1] << endl;
        output3.seekp(0, ios::end);
        output3 << endl;
    }
    output2[0].close();
    output2[1].close();
    output3.close();
    
    vector<double> u_previous1(nx); // for AB3
    vector<double> u_previous2(nx); // for AB3
    vector<double> u_past(nx); // for Leap-Frog
    
    // reset exact solution
    u_exact = exact(u_exact, 0);
    // reset save_solutions
    save_solutions1.clear();
    save_solutions2.clear();
    // reset current solution
    u_current1 = u_exact; // u_current for Leap-Frog
    u_current2 = u_exact; // u_current for AB3
    
    save_solutions1.push_back(u_current1);
    save_solutions2.push_back(u_current2);
    
    for (int file = 2; file <= 3; file++)
    {
        amplitude = find_amplitude(u_current1);
        error = find_error(u_current1, 0*dt);
        output1[file].seekp(0, ios::end);
        output1[file] << left << setw(width) << 0*dt << left << setw(width) << amplitude << left << setw(width) << error << endl;
    }
    u_past = u_current1; // set to initial
    u_next1 = FTBS(u_current1); // advance first iteration for Leap-Frog using FTBS
    u_current1 = u_next1; // update current Leap-Frog w/ FTBS (first iteration)
    
    u_previous2 = u_current2; // set to initial
    u_previous1 = u_current2; // set to initial
    u_next2 = FTBS(u_current2); // advance first iteration for AB3 using FTBS
    u_current2 = u_next2; // update current AB3 w/ FTBS (first iteration)
    
    // output time, amplitude, error for Leap-Frog
    amplitude = find_amplitude(u_current1);
    error = find_error(u_current1, 1*dt);
    output1[2].seekp(0, ios::end);
    output1[2] << left << setw(width) << 1*dt << left << setw(width) << amplitude << left << setw(width) << error << endl;
    
    u_next1 = leapfrog(u_current1, u_past); // leap-frog advance
    u_past = u_current1;
    u_current1 = u_next1;
    
    // output time, amplitude, error for AB3
    amplitude = find_amplitude(u_current2);
    error = find_error(u_current2, 1*dt);
    output1[3].seekp(0, ios::end);
    output1[3] << left << setw(width) << 1*dt << left << setw(width) << amplitude << left << setw(width) << error << endl;
    
    u_next2 = FTBS(u_current2); // advance second iteration for AB3 using FTBS
    u_previous1 = u_current2; // update previous1
    u_current2 = u_next2;
    
    for (int n = 2; n <= nt; n++)
    {
        // output time, amplitude, error for Leap-Frog
        amplitude = find_amplitude(u_current1);
        error = find_error(u_current1, n*dt);
        output1[2].seekp(0, ios::end);
        output1[2] << left << setw(width) << n*dt << left << setw(width) << amplitude << left << setw(width) << error << endl;
        u_current1 = u_next1;
        // output time, amplitude, error for AB3
        amplitude = find_amplitude(u_current2);
        error = find_error(u_current2, n*dt);
        output1[3].seekp(0, ios::end);
        output1[3] << left << setw(width) << n*dt << left << setw(width) << amplitude << left << setw(width) << error << endl;
        // advance u_next1 w/ Leap-Frog
        u_next1 = leapfrog(u_current1, u_past);
        u_past = u_current1;
        // advance u_next2 w/ AB3
        u_next2 = AB3(u_current2, u_previous1, u_previous2);
        u_previous2 = u_previous1;
        u_previous1 = u_current2;
        u_current2 = u_next2;
        
        if (n%save_frequency == 0)
        {
            save_solutions1.push_back(u_current1); // save Leap-Frog iteration
            save_solutions2.push_back(u_current2); // save AB3 iteration
        }
    }
    output1[0].close();
    output1[1].close();
    
    // output full solutions
    output2[2].open("LF_full.csv", ifstream::out | ifstream::trunc);
    output2[2].seekp(0, ios::beg);
    output2[3].open("AB3_full.csv", ifstream::out | ifstream::trunc);
    output2[3].seekp(0, ios::beg);
    for (int i = 0; i < nx; i++)
    {
        header[i] = "u(" + to_string(i) + ")";
        output2[2] << left << setw(width) << header[i];
        output2[3] << left << setw(width) << header[i];
        output3 << left << setw(width) << header[i];
        for (int j = 0; j < save_solutions1.size(); j++)
        {
            output2[2].seekp(0, ios::end);
            output2[2] << left << setw(width) << save_solutions1[j][i];
            output2[3].seekp(0, ios::end);
            output2[3] << left << setw(width) << save_solutions2[j][i];
        }
        output2[2].seekp(0, ios::end);
        output2[2] << endl;
        output2[3].seekp(0, ios::end);
        output2[3] << endl;
    }
    output2[2].close();
    output2[3].close();
    
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
        // FOR SQUARE WAVE ONLY
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
        // FOR SQUARE WAVE
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

vector<double> FTBS (vector<double> u)
{
    // advance the solution for interior nodes
    vector<double> u_next(u.size());
    for (int i = 1; i < nx-1; i++)
    {
        u_next[i] = u[i] - Courant*(u[i] - u[i-1]);
    }
    // handle boundary conditions
    // right boundary condition (backward in space)
    u_next[nx-1] = u[nx-1] - Courant*(u[nx-1] - u[nx-2]);
    // left boundary condition (cyclic boundary)
    u_next[0] = u_next[nx-1];
    
    return u_next;
}

vector<double> RK3 (vector<double> u)
{
    // create first approximation
    // use Euler-Forward to advance dt/3
    vector<double> u1(u.size());
    for (int i = 1; i < u.size()-1; i++)
    {
        u1[i] = u[i] - (Courant/6)*(u[i+1] - u[i-1]);
    }
    u1[u.size()-1] = u[u.size()-1] - (Courant/6)*(u[1] - u[u.size()-2]);
    u1[0] = u1[u.size()-1];
    
    // create second approximation
    // advance dt/2 from original time
    vector<double> u2(u.size());
    for (int i = 1; i < u.size()-1; i++)
    {
        u2[i] = u[i] - (Courant/4)*(u1[i+1] - u1[i-1]);
    }
    u2[u.size()-1] = u[u.size()-1] - (Courant/4)*(u1[1] - u1[u.size()-2]);
    u2[0] = u2[u.size()-1];
    
    // create final approximation
    // advance dt from original time
    vector<double> u_next(u.size());
    for (int i = 1; i < u.size()-1; i++)
    {
        u_next[i] = u[i] - (Courant/2)*(u2[i+1] - u2[i-1]);
    }
    u_next[u.size()-1] = u[u.size()-1] - (Courant/2)*(u2[1] - u2[u.size()-2]);
    u_next[0] = u_next[u.size()-1];
    
    // return final approximation
    return u_next;
}

vector<double> leapfrog (vector<double> u_current, vector<double> u_previous)
{
    vector<double> u_next(u_current.size());
    for (int i = 1; i < u_current.size(); i++)
    {
        u_next[i] = u_previous[i] - Courant*(u_current[i+1]-u_current[i-1]);
    }
    u_next[nx-1] = u_previous[nx-1] - Courant*(u_current[1]-u_current[nx-2]);
    u_next[0] = u_next[nx-1];
    
    return u_next;
}

vector<double> AB3 (vector<double> u_current, vector<double> u_previous1, vector<double> u_previous2)
{
    vector<double> u_next(u_current.size());
    for (int i = 1; i < u_current.size()-1; i++)
    {
        u_next[i] = u_current[i] - (Courant/24)*(23*(u_current[i+1]-u_current[i-1]) - 16*(u_previous1[i+1]-u_previous1[i-1]) + 5*(u_previous2[i+1]-u_previous2[i-1]));
    }
    u_next[u_current.size()-1] = u_current[u_current.size()-1] - (Courant/24)*(23*(u_current[u_current.size()]-u_current[u_current.size()-2]) - 16*(u_previous1[u_current.size()]-u_previous1[u_current.size()-2]) + 5*(u_previous2[u_current.size()]-u_previous2[u_current.size()-2]));
    u_next[0] = u_next[u_current.size()-1];
    
    return u_next;
}
