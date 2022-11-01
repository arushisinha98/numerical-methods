//
//  main.cpp
//  Modeling Task 5
//
//  Created by Arushi Sinha on 6/5/19.
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
double A = 1.0;
// amplitude
double r = 25.0;
// radius
double L_x = 100;
// domain size (x) = 100 m
double L_y = 100;
// domain size (y) = 100 m
int nx = 51;
// 51 grid points in 0 <= x <= L_x
int ny = 51;
// 51 grid points in 0 <= y <= L_y
double dx = L_x/(nx-1);
// delta x
double dy = L_y/(ny-1);
// delta y
int save_frequency = 10;
// save full solution every 10 iterations

// (x01, y01)
double x_start1 = 50.0;
double y_start1 = 25.0;

// (x02, y02)
double x_start2 = 50.0;
double y_start2 = 75.0;

// RUN A
// int initial_y = 1;
// int end_y = ny;
// RUN B
int initial_y = 2;
int end_y = ny-1;

// formatting values
int width = 14; // width of fixed-length data fields
int precision = 8; // precision of data fields
// error bound values
double tolerance = 0.0000001;
double max_iterations = 2000;

vector<vector<double>> initialize (vector<vector<double>>& u);

void print(vector<vector<double>> u, ofstream& output);

double sum (vector<vector<double>> u);

double max_norm (vector<vector<double>> u);

double one_norm (vector<vector<double>> u);

double max (vector<vector<double>> u);

// remove ghost nodes to calculate using physical node values
vector<vector<double>> remove_ghost_nodes (vector<vector<double>> u);

int main()
{
    double beta = dx/dy;
    double sigma = (1/(1+pow(beta,2)))*(cos(M_PI/nx)+pow(beta,2)*cos(M_PI/ny));
    double alpha = 2/(1 + sqrt(1 - pow(sigma,2)));
    
    vector<vector<double>> vorticity(ny+2, vector<double>(nx+2));
    vorticity = initialize(vorticity);
    // sum vorticity over physical nodes
    double s = 0;
    for (int j = 1; j < ny+1; j++)
    {
        for (int i = 1; i < nx+1; i++)
        {
            s += vorticity[j][i];
        }
    }
    s = s/(ny*nx);
    int counter = 0;
    while (s > tolerance)
    {
        // subtract average from each point
        for (int col = 0; col < ny+2; col++)
        {
            for (int row = 0; row < nx+2; row++)
            {
                vorticity[col][row] = vorticity[col][row] - s;
            }
        }
        s = 0;
        for (int j = 1; j < ny+1; j++)
        {
            for (int i = 1; i < nx+1; i++)
            {
                s += vorticity[j][i];
            }
        }
        s = s/(ny*nx);
        counter++;
    }
    
    ofstream vort;
    // output normalized vorticity
    vort.open("Normalized_vorticity.csv", ifstream::out | ifstream::trunc);
    vort.seekp(0, ios::beg);
    print(vorticity,vort);
    vort.close();
    
    vector<vector<double>> streamfunction(ny+2, vector<double>(nx+2));
    vector<vector<double>> residual(ny+2, vector<double>(nx+2));
    for (int j = 0; j < ny+2; j++)
    {
        for (int i = 0; i < nx+2; i++)
        {
            streamfunction[j][i] = 0.0;
            residual[j][i] = (-1.0)*vorticity[j][i];
        }
    }

    ofstream r;
    r.open("Initial_r.csv", ifstream::out | ifstream::trunc);
    r.seekp(0, ios::beg);
    print(residual,r);
    r.close();
    
    double R = max_norm(residual);
    double V = max_norm(vorticity);
    double S = one_norm(streamfunction);
    double epsilon = R/((2*((1/pow(dx,2))+(1/pow(dy,2)))*S) + V);
    
    int iterations = 0;
    
    ofstream output1;
    output1.open("control_output.csv", ifstream::out | ifstream::trunc);
    output1.seekp(0, ios::beg);
    output1 << left << setw(width) << "n" << left << setw(width) << "Epsilon" << left << setw(width) << "Max{R}" << endl;

    while (epsilon > tolerance && iterations < max_iterations)
    {
        for (int j = initial_y; j <= end_y; j++)
        {
            if (j == end_y)
            {
                for (int i = 1; i < nx+1; i++)
                {
                    residual[ny+1][i] = residual[2][i]; // update ghost at y = ny+1
                    streamfunction[ny+1][i] = streamfunction[2][i];
                }
            }
            for (int i = 1; i <= nx; i++)
            {
                if (i == nx)
                {
                    residual[j][nx+1] = residual[j][2]; // update ghost at x = nx+1
                    streamfunction[j][nx+1] = streamfunction[j][2];
                }
                // calculate residual for interior nodes
                residual[j][i] = (1/pow(dx,2))*(streamfunction[j][i-1] - 2*streamfunction[j][i] + streamfunction[j][i+1]) + (1/pow(dy,2))*(streamfunction[j-1][i] - 2*streamfunction[j][i] + streamfunction[j+1][i]) - vorticity[j][i];
                // calculate streamfunction
                streamfunction[j][i] = streamfunction[j][i] + alpha*residual[j][i]/(2*((1/(pow(dx,2)))+(1/pow(dy,2))));
            }
            residual[j][0] = residual[j][nx-1]; // update ghost at x = 0
            streamfunction[j][0] = streamfunction[j][nx-1];
            streamfunction[j][1] = streamfunction[j][nx]; // enforce (x = 1) = (x = 51)
        }
        for (int i = 1; i < nx+1; i++)
        {
            residual[0][i] = residual[ny-1][i]; // update ghost at y = 0
            streamfunction[0][i] = streamfunction[ny-1][i];
        }
        streamfunction[1] = streamfunction[ny]; // enforce (y = 1) = (y = 51)
        
        R = max_norm(remove_ghost_nodes(residual));
        V = max_norm(vorticity);
        S = one_norm(remove_ghost_nodes(streamfunction));
        epsilon = R/((2*((1/pow(dx,2))+(1/pow(dy,2)))*S) + V);
        
        output1 << left << setw(width) << fixed << setprecision(precision) << iterations << left << setw(width) << epsilon << left << setw(width) << max(residual) << endl;
        
        if (iterations%save_frequency == 0)
        {
            ofstream output3;
            ofstream output4;
            ofstream output5;
            string header = "Streamfunction_" + to_string(iterations) + ".csv";
            output3.open(header, ifstream::out | ifstream::trunc);
            output3.seekp(0, ios::beg);
            print(streamfunction,output3);
            header.clear();
            header = "u_" + to_string(iterations) + ".csv";
            output4.open(header, ifstream::out | ifstream::trunc);
            output4.seekp(0, ios::beg);
            for (int j = 1; j < ny+1; j++)
            {
                for (int i = 1; i < nx+1; i++)
                {
                    output4 << fixed << setprecision(precision) << left << setw(width) << (streamfunction[j-1][i] - streamfunction[j+1][i])/(2*dy);
                }
                output4 << endl;
            }
            header.clear();
            header = "v_" + to_string(iterations) + ".csv";
            output5.open(header, ifstream::out | ifstream::trunc);
            output5.seekp(0, ios::beg);
            for (int j = 1; j < ny+1; j++)
            {
                for (int i = 1; i < nx+1; i++)
                {
                    output5 << fixed << setprecision(precision) << left << setw(width) << (streamfunction[j][i+1] - streamfunction[j][i-1])/(2*dx);
                }
                output5 << endl;
            }
        }
        iterations++;
        cout << epsilon << endl;
    }
    vector<vector<double>> vort_num(ny+2, vector<double>(nx+2));
    for (int j = 1; j < ny+1; j++)
    {
        for (int i = 1; i < nx+1; i++)
        {
            vort_num[j][i] = (1/pow(dx,2))*(streamfunction[j][i-1] - 2*streamfunction[j][i] + streamfunction[j][i+1]) + (1/pow(dy,2))*(streamfunction[j-1][i] - 2*streamfunction[j][i] + streamfunction[j+1][i]);
        }
    }
    ofstream num;
    num.open("Numerical_vort.csv", ifstream::out | ifstream::trunc);
    num.seekp(0, ios::beg);
    for (int j = 1; j < ny+1; j++)
    {
        for (int i = 1; i < nx+1; i++)
        {
            num << left << setw(width) << vort_num[j][i];
        }
        num << endl;
    }
    return 0;
}

vector<vector<double>> initialize (vector<vector<double>>& u)
{
    // vorticity to be initialized contains ghost nodes
    
    // initial condition
    // vorticity(x,y) = (1/2) A [cos(d1*pi) + 1] + (1/2) A [cos(d2*pi) + 1];
    // d = min{1, (1/r)[(x-x0)^2 + (y-y0)^2]^(1/2)}
    double d1;
    double d2;
    for (int col = 1; col < u.size()-1; col++)
    {
        for (int row = 1; row < u[0].size()-1; row++)
        {
            double x = static_cast<double>(2*(row-1));
            double y = static_cast<double>(2*(col-1));
            d1 = min(1.0,(1.0/r)*(sqrt(pow(x - x_start1, 2) + pow(y - y_start1, 2))));
            d2 = min(1.0,(1.0/r)*(sqrt(pow(x - x_start2, 2) + pow(y - y_start2, 2))));
            u[col][row] = (1.0/2.0)*A*(cos(d1*M_PI) + 1) + (1.0/2.0)*A*(cos(d2*M_PI) + 1);
        }
    }
    // ghost nodes in y
    u[0] = u[50];
    u[52] = u[2];
    // ghost nodes in x
    for (int j = 0; j < u.size(); j++)
    {
        u[j][0] = u[j][50];
        u[j][52] = u[j][2];
    }
    ofstream vort;
    vort.open("Initial_vorticity.csv", ifstream::out | ifstream::trunc);
    vort.seekp(0, ios::beg);
    print(u,vort);
    vort.close();
    return u;
}

double sum (vector<vector<double>> u)
{
    double sum = 0;
    for (int col = 0; col < ny-1; col++)
    {
        for (int row = 0; row < nx-1; row++)
        {
            sum += u[col][row];
        }
    }
    return sum;
}

double max_norm (vector<vector<double>> u)
{
    double max = abs(u[0][0]);
    for (int col = 0; col < u.size(); col++)
    {
        for (int row = 0; row < u[0].size(); row++)
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
    for (int col = 0; col < u.size(); col++)
    {
        for (int row = 0; row < u[0].size(); row++)
        {
            max += abs(u[col][row]);
        }
    }
    return max;
}

double max (vector<vector<double>> u)
{
    double max = u[0][0];
    for (int col = 0; col < u.size(); col++)
    {
        for (int row = 0; row < u[0].size(); row++)
        {
            if (u[col][row] > max)
            {
                max = u[col][row];
            }
        }
    }
    return max;
}

vector<vector<double>> remove_ghost_nodes (vector<vector<double>> u)
{
    vector<vector<double>> original(ny, vector<double>(nx));
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
    for (int j = 1; j < ny+1; j++)
    {
        for (int i = 1; i < nx+1; i++)
        {
            output << left << setw(width) << u[j][i];
        }
        output << endl;
    }
}
