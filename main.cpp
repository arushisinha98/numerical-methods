//
//  main.cpp
//  Modeling Task 4
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
double T_total = 400;
// total duration of observation (final time step)
double A = 1;
// amplitude
double c_x = 1.0;
// propagation velocity in x-direction = 1.0 m/s
double c_y = 1.0;
// propagation velocity in y-direction = 1.0 m/s
double L_x = 50;
// domain size (x) = 50 m
double L_y = 50;
// domain size (y) = 50 m
int nx = 51;
// 51 grid points in 0 <= x <= L_x
int ny = 51;
// 51 grid points in 0 <= y <= L_y
double dx = L_x/(nx-1);
// delta x
double dy = L_y/(ny-1);
// delta y
double K_x = 0.5;
// diffusivity in x
double K_y = 0.5;
// diffusivity in y
double r = 12.5;
// radius = 12.5 m
int save_frequency = 100;
// save full solution every 100 iterations

// (x0, y0)
double x_start = 25.0;
double y_start = 25.0;

double dt = 0.1;
// delta t
double nt = T_total/dt;
// number of timesteps

// formatting values
int width = 12; // width of fixed-length data fields
int precision = 6; // precision of data fields

// initialize u(x,y,t = 0)
vector<vector<double>> initialize(vector<vector<double>>& u);

// check conservation
double integrate(vector<vector<double>> u);

// find maximum value of u
double max_height(vector<vector<double>> u);

// add ghost nodes to iterate easily with boundary conditions
vector<vector<double>> add_ghost_nodes (vector<vector<double>> u);

// Euler Forward (used in first iteration for advection, before leap-frog)
vector<vector<double>> FTBS (vector<vector<double>> u_current);

// Leap-frog
vector<vector<double>> leapfrog (vector<vector<double>> u_current, vector<vector<double>> u_previous);

// Forward-in-Time, Centered-in-Space
vector<vector<double>> FTCS (vector<vector<double>> u_previous);

// generate cuts of solution
// first diagonal: y = x
// second diagonal: y = L_x - x
vector<vector<double>> diagonal_cut(vector<vector<double>> u);

// C(i,j) = A(i,j) + B(i,j)
vector<vector<double>> sum(vector<vector<double>> A, vector<vector<double>> B);

// C(i,j) = A(i,j) - B(i,j)
vector<vector<double>> subtract(vector<vector<double>> A, vector<vector<double>> B);

// to create a square wave
int sign(double value);

 int main()
{
    // save integral and maximum{u} every time-step
    ofstream output1;
    output1.open("control_output.csv", ifstream::out | ifstream::trunc);
    output1.seekp(0, ios::beg);
    output1 << left << setw(width) << "Time" << left << setw(width) << "Max{u(x,y)}" << left << setw(width) << "Integral" << endl;
    // save full solution every 10 time-steps
    ofstream output2;
    // save diagonal cuts of solution every 10 time-steps
    ofstream output3;
    
    vector<vector<double>> u_previous(ny, vector<double> (nx));
    u_previous = initialize(u_previous);
    double max_u = max_height(u_previous);
    double integral = integrate(u_previous);
    
    output1 << fixed << setprecision(precision) << left << setw(width) << 0 << left << setw(width) << max_u << left << setw(width) << integral << endl;
    
    output2.open("2D_field_0.csv", ifstream::out | ifstream::trunc);
    output2.seekp(0, ios::beg);
    // print with (0,0) in bottom left corner
    for (int j = ny-1; j >= 0; j--)
    {
        for (int i = 0; i < nx; i++)
        {
            output2 << left << setw(width) << u_previous[j][i];
        }
        output2 << endl;
    }
    output2.close();
    output3.open("Slice_0.csv", ifstream::out | ifstream::trunc);
    output3.seekp(0, ios::beg);
    vector<vector<double>> diagonal = diagonal_cut(u_previous);
    for (int j = 0; j <= 1; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            output3 << left << setw(width) << fixed << setprecision(precision) << diagonal[j][i];
        }
        output3 << endl;
    }
    output3.close();
    // PART 1: ADVECTION
    vector<vector<double>> u_current = FTBS(u_previous);
    vector<vector<double>> advection = leapfrog(u_current, u_previous);
    // PART 2: DIFFUSION
    vector<vector<double>> diffusion = FTCS(u_previous);
    vector<vector<double>> u_next;
    // PART 3: ADVECTION-DIFFUSION
    for (int n = 2; n <= nt; n++)
    {
        // u_next = advection; // EXPERIMENT #1
        // u_next = diffusion; // EXPERIMENT #2
        // EXPERIMENT #3
        u_next = sum(advection, diffusion);
        u_next = subtract(u_next, u_previous);
        max_u = max_height(u_next);
        integral = integrate(u_next);
        u_previous = u_current;
        u_current = u_next;
        advection = leapfrog(u_current, u_previous);
        diffusion = FTCS(u_previous);
        
        output1 << left << setw(width) << fixed << setprecision(precision) << n*dt << left << setw(width) << max_u << left << setw(width) << integral << endl;
        
        if (n%save_frequency == 0)
        {
            ofstream output2;
            ofstream output3;
            string header = "2D_field_" + to_string(n) + ".csv";
            output2.open(header, ifstream::out | ifstream::trunc);
            output2.seekp(0, ios::beg);
            // print with (0,0) in bottom left corner
            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    output2 << left << setw(width) << u_current[j][i]*1000000;
                }
                output2 << endl;
            }
            header.clear();
            header = "Slice_" + to_string(n) + ".csv";
            output3.open(header, ifstream::out | ifstream::trunc);
            output3.seekp(0, ios::beg);
            diagonal = diagonal_cut(u_current);
            for (int j = 0; j <= 1; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    output3 << left << setw(width) << fixed << setprecision(precision) << diagonal[j][i];
                }
                output3 << endl;
            }
            output2.close();
            output3.close();
        }
    }
    output1.close();
    output2.close();
    output3.close();
    return 0;
}

vector<vector<double>> initialize(vector<vector<double>>& u)
{
    // initial condition
    // u(x,y,t = 0) = u0(x,y) = (1/2) A [cos(d*pi) + 1];
    // d = min{1, (1/r)[(x-x0)^2 + (y-y0)^2]^(1/2)}
    double d;
    for (int col = 0; col < ny; col++)
    {
        for (int row = 0; row < nx; row++)
        {
            double x = static_cast<double>(row);
            double y = static_cast<double>(col);
            d = min(1.0,(1.0/r)*(pow(pow(x - x_start, 2) + pow(y - y_start, 2),(1.0/2.0))));
            u[col][row] = (1.0/2.0)*A*(cos(d*M_PI) + 1);
        }
    }
    return u;
}

double integrate(vector<vector<double>> u)
{
    double sum = 0;
    for (int y = 0; y < ny-1; y++)
    {
        for (int x = 0; x < nx-1; x++)
        {
            sum += u[y][x];
        }
    }
    return sum;
}

double max_height(vector<vector<double>> u)
{
    vector<double> height;
    for (int j = 0; j < u.size(); j++)
    {
        for (int i = 0; i < u[0].size(); i++)
        {
            height.push_back(u[j][i]);
        }
    }
    double max = *max_element(height.begin(),height.end());
    return max;
}

vector<vector<double>> FTBS (vector<vector<double>> u_current)
{
    vector<vector<double>> ghost = add_ghost_nodes(u_current);
    vector<vector<double>> u_next(ny, vector<double>(nx));
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            u_next[j][i] = ghost[j+1][i+1] - c_x*(dt/dx)*(ghost[j+1][i+1] - ghost[j+1][i]) - c_y*(dt/dy)*(ghost[j+1][i+1] - ghost[j][i+1]);
        }
    }
    return u_next;
}

vector<vector<double>> leapfrog (vector<vector<double>> u_current, vector<vector<double>> u_previous)
{
    vector<vector<double>> ghost = add_ghost_nodes(u_current);
    vector<vector<double>> u_next(ny, vector<double> (nx));
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            u_next[j][i] = u_previous[j][i] - c_x*(dt/dx)*(ghost[j+1][i+2] - ghost[j+1][i]) - c_y*(dt/dy)*(ghost[j+2][i+1] - ghost[j][i+1]);
        }
    }
    return u_next;
}

vector<vector<double>> FTCS (vector<vector<double>> u_previous)
{
    vector<vector<double>> ghost = add_ghost_nodes(u_previous);
    vector<vector<double>> u_next(ny, vector<double> (nx));
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            u_next[j][i] = u_previous[j][i] + K_x*(2*dt/pow(dx,2))*(ghost[j+1][i+2] - 2*ghost[j+1][i+1] + ghost[j+1][i]) + K_y*(2*dt/pow(dy,2))*(ghost[j+2][i+1] - 2*ghost[j+1][i+1] + ghost[j][i+1]);
        }
    }
    return u_next;
}

vector<vector<double>> add_ghost_nodes (vector<vector<double>> u)
{
    vector<vector<double>> ghost(ny+2, vector<double> (nx+2));
    for (int j = 0; j < u.size(); j++)
    {
        for (int i = 0; i < u.size(); i++)
        {
            ghost[j+1][i+1] = u[j][i];
        }
    }
    // ghost nodes: L_y + 1 = 1
    ghost[ny+1] = ghost[2];
    // ghost nodes: -1 = L_y - 1
    ghost[0] = ghost[ny-1];
    for (int j = 0; j < ny+2; j++)
    {
        // ghost nodes: L_x + 1 = 1
        ghost[j][nx+1] = ghost[j][2];
        // ghost nodes: -1 = L_y - 1
        ghost[j][0] = ghost[j][nx-1];
    }
    
    return ghost;
}

vector<vector<double>> diagonal_cut(vector<vector<double>> u)
{
    vector<double> first;
    vector<double> second;
    for (int j = 0; j < u.size(); j++)
    {
        first.push_back(u[j][j]);
        second.push_back(u[u.size()-1-j][j]);
    }
    vector<vector<double>> diagonals;
    diagonals.push_back(first);
    diagonals.push_back(second);
    return diagonals;
}

vector<vector<double>> sum(vector<vector<double>> A, vector<vector<double>> B)
{
    assert(A.size() == B.size());
    assert(A[0].size() == B[0].size());
    vector<vector<double>> sum(A.size(), vector<double>(A[0].size()));
    for (int j = 0; j < A.size(); j++)
    {
        for (int i = 0; i < A[0].size(); i++)
        {
            sum[j][i] = A[j][i] + B[j][i];
        }
    }
    return sum;
}

vector<vector<double>> subtract(vector<vector<double>> A, vector<vector<double>> B)
{
    assert(A.size() == B.size());
    assert(A[0].size() == B[0].size());
    vector<vector<double>> difference(A.size(), vector<double>(A[0].size()));
    for (int j = 0; j < A.size(); j++)
    {
        for (int i = 0; i < A[0].size(); i++)
        {
            difference[j][i] = A[j][i] - B[j][i];
        }
    }
    return difference;
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

