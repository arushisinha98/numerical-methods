//
//  main.cpp
//  Modeling Task 1
//
//  Created by Arushi Sinha on 5/4/19.
//  Copyright © 2019 Arushi Sinha. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

// define and declare global variables
double t_total = 2000.0;
// total duration of observation (final time step)
double dt[2] = {10.0, 1.0};
// delta t (try: 10.0s, 1.0s)
double output_step = 1;
// number of time steps b/w outputs
double g = 9.81;
// gravitational constant
double theta = 290.0;
double dtheta_dz = 5.0/1000.0; // d(theta)/dz = 5 K/km = 0.005 K/m
double z0 = 100.0; // z0 = 100 m
double N_squared = g / theta * dtheta_dz;
double N = sqrt(N_squared);

// formatting values
int width = 12; // width of fixed-length data fields
int precision = 6; // precision of data fields


// Euler-Forward Function

// Leap-Frog Function

int main()
{
    // create files to store data
    ofstream file[3];
    file[0].open("EF_10s.txt", ifstream::out | ifstream::trunc);
    file[0].seekp(0, ios::beg);
    file[1].open("EF_1s.txt", ifstream:: out | ifstream::trunc);
    file[1].seekp(0, ios::beg);
    
    // loop to generate Euler-Forward dataset (10s iteration first followed by 1s iteration)
    for (int iteration = 0; iteration <= 1; iteration++)
    {
        // define and declare local variables
        double t_current = 0.0;
        double t_next = t_current;
        
        double z_current = z0; // initial condition: z(t = 0) = z0
        double z_next;
        double z_error = 0.0;
        double z_exact;
        
        double w_current = 0.0; // initial condition: w(t = 0) = 0
        double w_next;
        
        int nt = t_total / dt[iteration]; // total number of timesteps
        // headers: t   z(t) est    w(t) est    z(t) exact  error[z(t)]
        file[iteration] << left << setw(width) << "t" << left << setw(width) << "z(t) est" << left << setw(width) << "w(t) est" << left << setw(width) << "z(t) exact" << left << setw(width) << "error[z(t)]" << endl;
    
        // Euler-Forward equations:
        // z(n+1) = z(n) + dt * w(t(n))
        // w(n+1) = w(n) - dt * N^2 * z(t(n))
        for (int i = 0; i <= nt; i++)
        {
            t_next += dt[iteration];
            z_next = z_current + dt[iteration]*w_current;
            w_next = w_current - dt[iteration]*N_squared*z_current;
            z_exact = z0*cos(N*t_current);
            z_error = abs(z_exact - z_next);
        
            // output calculated variables to file
            // t_current, t_next
            // z_current, z_next
            // w_current, w_next
            file[iteration] << left << setw(width) << fixed << setprecision(precision) << t_current << left << setw(width) << z_current << left << setw(width) << w_current << left << setw(width) << z_exact << left << setw(width) << z_error << endl;
        
            // update current variables
            t_current = t_next;
            z_current = z_next;
            w_current = w_next;
        }
        file[iteration] << '\n' << endl;
    }
    file[0].close();
    file[1].close();
    
    file[0].open("LF_10s.txt", ifstream::in | ifstream::trunc);
    file[0].seekp(0, ios::beg);
    file[1].open("LF_1s.txt", ifstream::in | ifstream::trunc);
    file[1].seekp(0, ios::beg);
    
    // loop to generate Leap-Frog dataset (10s iteration first, followed by 1s iteration)
    for (int iteration = 0; iteration <= 1; iteration++)
    {
        // define and declare local variables
        double t_current = 0.0;
        double t_next = t_current;
        
        double z_current = z0; // initial condition: z(t = 0) = z0
        double z_error = 0.0;
        double z_exact;
        double z_next;
        double z_previous;
        
        double w_current = 0.0; // initial condition: w(t = 0) = 0
        double w_next;
        double w_previous;
        
        double nt = t_total / dt[iteration]; // total number of timesteps
        // headers: t   z(t) est    w(t) est    z(t) exact  error[z(t)]
        file[iteration] << left << setw(width) << "t" << left << setw(width) << "z(t) est" << left << setw(width) << "w(t) est" << left << setw(width) << "z(t) exact" << left << setw(width) << "error[z(t)]" << endl;
        
        t_next += dt[iteration];
        z_next = z_current + dt[iteration]*w_current;
        w_next = w_current - dt[iteration]*N_squared*z_current;
        z_exact = z0*cos(N*t_current);
        z_error = abs(z_exact - z_next);
        
        // output calculated variables to file
        // t_current, t_next
        // z_current, z_next
        // w_current, w_next
        file[iteration] << left << setw(width) << fixed << setprecision(precision) << t_current << left << setw(width) << z_current << left << setw(width) << w_current << left << setw(width) << z_exact << left << setw(width) << z_error << endl;
        
        t_current = t_next;
        z_previous = z_current;
        z_current = z_next;
        w_previous = w_current;
        w_current = w_next;
        
        // Leap-Frog equations:
        // z(n+1) = z(n-1) + 2 * dt * w(t(n))
        // w(n+1) = w(n-1) - 2 * dt * N^2 * z(t(n))
        for (int i = 0; i <= nt-1; i++)
        {
            t_next += dt[iteration];
            z_next = z_previous + 2*dt[iteration]*w_current;
            w_next = w_previous - 2*dt[iteration]*N_squared*z_current;
            z_exact = z0*cos(N*t_current);
            z_error = abs(z_exact - z_next);
            
            file[iteration] << left << setw(width) << fixed << setprecision(precision) << t_current << left << setw(width) << z_current << left << setw(width) << w_current << left << setw(width) << z_exact << left << setw(width) << z_error << endl;
            
            t_current = t_next;
            z_previous = z_current;
            z_current = z_next;
            w_previous = w_current;
            w_current = w_next;
        }
        file[iteration] << '\n' << endl;
    }
    file[0].close();
    file[1].close();
    
    return 0;
}
