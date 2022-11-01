//
//  main.cpp
//  Modeling Task 6
//
//  Created by Arushi Sinha on 15/5/19.
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
double a = 180;
// radius of vortex = 180 km
double b = 600;
// distance b/w vortices = 600 km
double Gamma = 8*0.00001*M_PI*pow(a,2);

double beta = 2*pow(10,-8);
// coriolis
double L_x = 2000;
// domain size (x) = 2000 km
double L_y = 2000;
// domain size (y) = 2000 km
int nx = 101;
// 101 grid points in 0 <= x <= L_x
int ny = 101;
// 101 grid points in 0 <= y <= L_y
double dx = L_x/(nx-1);
// delta x
double dy = L_y/(ny-1);
// delta y
double dt = 500;
// delta t = 100s
double T_total = 20*24*60*60;
// T(total) = 20 days (convert to seconds)
int nt = T_total/dt;
int save_frequency = 50;
// save full solution every 50 iterations

// (x01, y01)
double x_start1 = (L_x/2)-(b/2);
double y_start1 = L_y/2;

// (x02, y02)
double x_start2 = (L_x/2)+(b/2);
double y_start2 = L_y/2;

// formatting values
int width = 14; // width of fixed-length data fields
int precision = 8; // precision of data fields
// error bound values
double tolerance = 0.0000001;
double max_iterations = 2000;
// calculate alpha
double sigma = (1/(1+pow((dx/dy),2)))*(cos(M_PI/nx)+pow((dx/dy),2)*cos(M_PI/ny));
double alpha = 2/(1 + sqrt(1 - pow(sigma,2)));

vector<vector<double>> initialize (vector<vector<double>>& u);
// prints only physical nodes
void print(vector<vector<double>> u, ofstream& output);

double sum (vector<vector<double>> u);

double max_norm (vector<vector<double>> u);

double one_norm (vector<vector<double>> u);

double max (vector<vector<double>> u);

void Poisson (vector<vector<double>>& streamfunction, vector<vector<double>> vorticity, vector<vector<double>>& residual, double& R, double& V, double& S, double& epsilon);

vector<vector<double>> Jacobian (vector<vector<double>> streamfunction,  vector<vector<double>> vorticity);

vector<vector<double>> EF (vector<vector<double>> current, vector<vector<double>> J);

vector<vector<double>> EF_coriolis (vector<vector<double>> current, vector<vector<double>> J, vector<vector<double>> SF);

vector<vector<double>> leapfrog (vector<vector<double>> previous, vector<vector<double>> J);

vector<vector<double>> AB3 (vector<vector<double>> current, vector<vector<double>> J, vector<vector<double>> J_1, vector<vector<double>> J_2);

vector<vector<double>> AB3_coriolis (vector<vector<double>> current, vector<vector<double>> J, vector<vector<double>> J_1, vector<vector<double>> J_2, vector<vector<double>> SF, vector<vector<double>> SF_1, vector<vector<double>> SF_2);

vector<vector<double>> u_velocity (vector<vector<double>> streamfunction);

vector<vector<double>> v_velocity (vector<vector<double>> streamfunction);

double calculate_KE (vector<vector<double>> u_velocity, vector<vector<double>> v_velocity);

double calculate_EN (vector<vector<double>> vorticity);

double calculate_Cmax (vector<vector<double>> u_velocity, vector<vector<double>> v_velocity);
// remove ghost nodes to calculate using physical node values
vector<vector<double>> remove_ghost_nodes (vector<vector<double>> u);
// add ghost nodes to calculate using five point stencil
vector<vector<double>> add_ghost_nodes (vector<vector<double>> u);

int main()
{
    cout << alpha << endl;
    vector<vector<double>> vorticity(ny+2, vector<double>(nx+2));
    vorticity = initialize(vorticity);
    vector<vector<double>> streamfunction(ny+2, vector<double>(nx+2));
    vector<vector<double>> residual(ny+2, vector<double>(nx+2));
    
    // initialize streamfunction to 0 over domain
    for (int j = 0; j < ny+2; j++)
    {
        for (int i = 0; i < ny+2; i++)
        {
            streamfunction[j][i] = 0.0;
        }
    }
    // calculate residual for the initial epsilon
    for (int j = 2; j < ny; j++)
    {
        for (int i = 2; i < nx; i++)
        {
            residual[j][i] = (1/pow(dx,2))*(streamfunction[j][i-1] - 2*streamfunction[j][i] + streamfunction[j][i+1]) + (1/pow(dy,2))*(streamfunction[j-1][i] - 2*streamfunction[j][i] + streamfunction[j+1][i]) - vorticity[j][i];
        }
    }
    
    double R = max_norm(remove_ghost_nodes(remove_ghost_nodes(residual)));
    double V = max_norm(remove_ghost_nodes(remove_ghost_nodes(vorticity)));
    double S = one_norm(remove_ghost_nodes(remove_ghost_nodes(streamfunction)));
    double epsilon = R/((2*((1/pow(dx,2))+(1/pow(dy,2)))*S) + V);
    
    int iterations = 0;
    while (epsilon > tolerance && iterations < max_iterations)
    {
        // run Poisson solver
        Poisson(streamfunction, vorticity, residual, R, V, S, epsilon);
        iterations++;
    }
    cout << "0: " << iterations << endl;
    vector<vector<double>> u = u_velocity(streamfunction);
    vector<vector<double>> v = v_velocity(streamfunction);
    // calculate conserved values
    // total kinetic energy
    double kinetic = calculate_KE(u,v);
    // total enstrophy
    double enstrophy = calculate_EN(vorticity);
    // Courant's number
    double Courant = calculate_Cmax(u,v);
    int loop = 0;
    ofstream output3;
    output3.open("Observations.csv", ifstream::out | ifstream::trunc);
    output3.seekp(0, ios::beg);
    output3 << left << setw(width) << "n" << left << setw(width) << "Kinetic" << left << setw(width) << "Enstrophy" << left << setw(width) << "Cmax" << left << setw(width) << "Iterations" << endl;
    output3 << fixed << setprecision(precision) << left << setw(width) << loop << left << setw(width) << kinetic << left << setw(width) << enstrophy << left << setw(width) << Courant << left << setw(width) << iterations << endl;

    ofstream output1[3];
    output1[0].open("Streamfunction_0.csv", ifstream::out | ifstream::trunc);
    output1[0].seekp(0, ios::beg);
    print(remove_ghost_nodes(streamfunction),output1[0]);
    
    ofstream output2[3];
    output2[0].open("Vorticity_0.csv", ifstream::out | ifstream::trunc);
    output2[0].seekp(0, ios::beg);
    print(remove_ghost_nodes(vorticity),output2[0]);
    
    ofstream outputA[3];
    outputA[0].open("u_0.csv", ifstream::out | ifstream::trunc);
    outputA[0].seekp(0, ios::beg);
    print(remove_ghost_nodes(u),outputA[0]);
    
    ofstream outputB[3];
    outputB[0].open("v_0.csv", ifstream::out | ifstream::trunc);
    outputB[0].seekp(0, ios::beg);
    print(remove_ghost_nodes(v),outputB[0]);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TIME STEP: n = 1
    ///////////////////////////////////////////////////////////////////////////////////////////////
    loop++;
    // Euler-Forward
    // #1: calculate Jacobian
    vector<vector<double>> J_2 = Jacobian(streamfunction, vorticity);
    vector<vector<double>> SF_2 = streamfunction;
    // #2: step vorticity in time
    // next at time n+1
    // current at time n
    // previous at time n-1
    vector<vector<double>> current = vorticity;
    vector<vector<double>> next = EF(current, J_2);
    // vector<vector<double>> next = EF_coriolis(current, J_2, SF_2);
    // #3: update streamfunction
    // initialize streamfunction to 0 over domain
    for (int j = 0; j < ny+2; j++)
    {
        for (int i = 0; i < ny+2; i++)
        {
            streamfunction[j][i] = 0.0;
        }
    }
    // calculate residual for epsilon
    for (int j = 2; j < ny; j++)
    {
        for (int i = 2; i < nx; i++)
        {
            residual[j][i] = (1/pow(dx,2))*(streamfunction[j][i-1] - 2*streamfunction[j][i] + streamfunction[j][i+1]) + (1/pow(dy,2))*(streamfunction[j-1][i] - 2*streamfunction[j][i] + streamfunction[j+1][i]) - next[j][i];
        }
    }
    R = max_norm(remove_ghost_nodes(remove_ghost_nodes(residual)));
    V = max_norm(remove_ghost_nodes(remove_ghost_nodes(next)));
    S = one_norm(remove_ghost_nodes(remove_ghost_nodes(streamfunction)));
    epsilon = R/((2*((1/pow(dx,2))+(1/pow(dy,2)))*S) + V);
    iterations = 0;
    while (epsilon > tolerance && iterations < max_iterations)
    {
        // run Poisson solver
        Poisson(streamfunction, next, residual, R, V, S, epsilon);
        iterations++;
    }
    cout << loop << ": " << iterations << endl;
    // #4: update variables
    u = u_velocity(streamfunction);
    v = v_velocity(streamfunction);
    kinetic = calculate_KE(u,v);
    enstrophy = calculate_EN(next);
    Courant = calculate_Cmax(u,v);
    output3 << fixed << setprecision(precision) << left << setw(width) << loop << left << setw(width) << kinetic << left << setw(width) << enstrophy << left << setw(width) << Courant << left << setw(width) << iterations << endl;
    
    output1[1].open("Streamfunction_1.csv", ifstream::out | ifstream::trunc);
    output1[1].seekp(0, ios::beg);
    print(remove_ghost_nodes(streamfunction),output1[1]);
    
    output2[1].open("Vorticity_1.csv", ifstream::out | ifstream::trunc);
    output2[1].seekp(0, ios::beg);
    print(remove_ghost_nodes(next),output2[1]);
    
    outputA[1].open("u_1.csv", ifstream::out | ifstream::trunc);
    outputA[1].seekp(0, ios::beg);
    print(remove_ghost_nodes(u),outputA[1]);
    
    outputB[1].open("v_1.csv", ifstream::out | ifstream::trunc);
    outputB[1].seekp(0, ios::beg);
    print(remove_ghost_nodes(v),outputB[1]);
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TIME STEP: n = 2
    ///////////////////////////////////////////////////////////////////////////////////////////////
    loop++;
    vector<vector<double>> J_1 = Jacobian(streamfunction, next);
    vector<vector<double>> SF_1 = streamfunction;
    // #2: step vorticity in time
    // next at time n+1
    // current at time n
    // previous at time n-1
    current = next;
    next = EF(current, J_1);
    // next = EF_coriolis(current, J_1, SF_1);
    // #3: update streamfunction
    // initialize streamfunction to 0 over domain
    for (int j = 0; j < ny+2; j++)
    {
        for (int i = 0; i < ny+2; i++)
        {
            streamfunction[j][i] = 0.0;
        }
    }
    // calculate residual for the initial epsilon
    for (int j = 2; j < ny; j++)
    {
        for (int i = 2; i < nx; i++)
        {
            residual[j][i] = (1/pow(dx,2))*(streamfunction[j][i-1] - 2*streamfunction[j][i] + streamfunction[j][i+1]) + (1/pow(dy,2))*(streamfunction[j-1][i] - 2*streamfunction[j][i] + streamfunction[j+1][i]) - next[j][i];
        }
    }
    R = max_norm(remove_ghost_nodes(remove_ghost_nodes(residual)));
    V = max_norm(remove_ghost_nodes(remove_ghost_nodes(next)));
    S = one_norm(remove_ghost_nodes(remove_ghost_nodes(streamfunction)));
    epsilon = R/((2*((1/pow(dx,2))+(1/pow(dy,2)))*S) + V);
    iterations = 0;
    while (epsilon > tolerance && iterations < max_iterations)
    {
        // run Poisson solver
        Poisson(streamfunction, next, residual, R, V, S, epsilon);
        iterations++;
    }
    cout << loop << ": " << iterations << endl;
    // #4: update variables
    u = u_velocity(streamfunction);
    v = v_velocity(streamfunction);
    kinetic = calculate_KE(u,v);
    enstrophy = calculate_EN(next);
    Courant = calculate_Cmax(u,v);
    output3 << fixed << setprecision(precision) << left << setw(width) << loop << left << setw(width) << kinetic << left << setw(width) << enstrophy << left << setw(width) << Courant << left << setw(width) << iterations << endl;
    
    output1[2].open("Streamfunction_2.csv", ifstream::out | ifstream::trunc);
    output1[2].seekp(0, ios::beg);
    print(remove_ghost_nodes(streamfunction),output1[2]);
    
    output2[2].open("Vorticity_2.csv", ifstream::out | ifstream::trunc);
    output2[2].seekp(0, ios::beg);
    print(remove_ghost_nodes(next),output2[2]);
    
    outputA[2].open("u_2.csv", ifstream::out | ifstream::trunc);
    outputA[2].seekp(0, ios::beg);
    print(remove_ghost_nodes(u),outputA[2]);
    
    outputB[2].open("v_2.csv", ifstream::out | ifstream::trunc);
    outputB[2].seekp(0, ios::beg);
    print(remove_ghost_nodes(v),outputB[2]);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TIME STEP: n = 3 onwards
    ///////////////////////////////////////////////////////////////////////////////////////////////
    loop++;
    vector<vector<double>> J;
    vector<vector<double>> SF;
    while (loop < nt)
    {
        J = Jacobian(streamfunction, next);
        SF = streamfunction;
        // #2: step vorticity in time
        // next at time n+1
        // current at time n
        // previous at time n-1
        current = next;
        next = AB3(current, J, J_1, J_2);
        // next = AB3_coriolis(current, J, J_1, J_2, SF, SF_1, SF_2);
        J_2 = J_1;
        SF_2 = SF_1;
        J_1 = J;
        SF_1 = SF;
        // #3: update streamfunction
        // initialize streamfunction to 0 over domain
        for (int j = 0; j < ny+2; j++)
        {
            for (int i = 0; i < ny+2; i++)
            {
                streamfunction[j][i] = 0.0;
            }
        }
        // calculate residual for the initial epsilon
        for (int j = 2; j < ny; j++)
        {
            for (int i = 2; i < nx; i++)
            {
                residual[j][i] = (1/pow(dx,2))*(streamfunction[j][i-1] - 2*streamfunction[j][i] + streamfunction[j][i+1]) + (1/pow(dy,2))*(streamfunction[j-1][i] - 2*streamfunction[j][i] + streamfunction[j+1][i]) - next[j][i];
            }
        }
        R = max_norm(remove_ghost_nodes(remove_ghost_nodes(residual)));
        V = max_norm(remove_ghost_nodes(remove_ghost_nodes(next)));
        S = one_norm(remove_ghost_nodes(remove_ghost_nodes(streamfunction)));
        epsilon = R/((2*((1/pow(dx,2))+(1/pow(dy,2)))*S) + V);
        iterations = 0;
        while (epsilon > tolerance && iterations < max_iterations)
        {
            // run Poisson solver
            Poisson(streamfunction, next, residual, R, V, S, epsilon);
            iterations++;
        }
        cout << loop << ": " << iterations << endl;
        // #4: update variables
        u = u_velocity(streamfunction);
        v = v_velocity(streamfunction);
        kinetic = calculate_KE(u,v);
        enstrophy = calculate_EN(next);
        Courant = calculate_Cmax(u,v);
        output3 << fixed << setprecision(precision) << left << setw(width) << loop << left << setw(width) << kinetic << left << setw(width) << enstrophy << left << setw(width) << Courant << left << setw(width) << iterations << endl;
        loop++;
        if (loop%save_frequency == 0)
        {
            ofstream output1;
            ofstream output2;
            ofstream outputA;
            ofstream outputB;
            
            string header = "Streamfunction_" + to_string(loop) + ".csv";
            output1.open(header, ifstream::out | ifstream::trunc);
            output1.seekp(0, ios::beg);
            print(remove_ghost_nodes(streamfunction),output1);
            
            header.clear();
            header = "Vorticity_" + to_string(loop) + ".csv";
            output2.open(header, ifstream::out | ifstream::trunc);
            output2.seekp(0, ios::beg);
            print(remove_ghost_nodes(next),output2);
            
            header.clear();
            header = "u_" + to_string(loop) + ".csv";
            outputA.open(header, ifstream::out | ifstream::trunc);
            outputA.seekp(0, ios::beg);
            print(remove_ghost_nodes(u),outputA);
            
            header.clear();
            header = "v_" + to_string(loop) + ".csv";
            outputB.open(header, ifstream::out | ifstream::trunc);
            outputB.seekp(0, ios::beg);
            print(remove_ghost_nodes(v),outputB);
            
            header.clear();
        }
    }
    return 0;
}

vector<vector<double>> initialize (vector<vector<double>>& u)
{
    // initial condition
    // vorticity(x,y)
    double A = Gamma/(M_PI*pow(a,2));
    for (int col = 2; col < ny; col++)
    {
        for (int row = 2; row < nx; row++)
        {
            double x = static_cast<double>(dx*(row-1));
            double y = static_cast<double>(dy*(col-1));
            u[col][row] = A*exp((-1)*(pow((x-x_start1),2) + pow((y-y_start1),2))/pow(a,2)) + A*exp((-1)*(pow((x-x_start2),2) + pow((y-y_start2),2))/pow(a,2));
        }
    }
    return u;
}

void Poisson (vector<vector<double>>& streamfunction, vector<vector<double>> vorticity, vector<vector<double>>& residual, double& R, double& V, double& S, double& epsilon)
{
    for (int j = 2; j < ny; j++)
    {
        if (j == ny-1)
        {
            for (int i = 2; i < nx; i++)
            {
                streamfunction[ny+1][i] = streamfunction[2][i]; // update ghost at y = ny+1
            }
        }
        for (int i = 2; i < nx; i++)
        {
            // calculate residual for interior nodes
            residual[j][i] = (1/pow(dx,2))*(streamfunction[j][i-1] - 2*streamfunction[j][i] + streamfunction[j][i+1]) + (1/pow(dy,2))*(streamfunction[j-1][i] - 2*streamfunction[j][i] + streamfunction[j+1][i]) - vorticity[j][i];
            // calculate streamfunction
            streamfunction[j][i] = streamfunction[j][i] + alpha*residual[j][i]/(2*((1/(pow(dx,2)))+(1/pow(dy,2))));
        }
        // streamfunction[j][0] = streamfunction[j][nx-1]; // update ghost at x = 0
        streamfunction[j][1] = streamfunction[j][nx]; // enforce (x = 1) = (x = nx) = 0
    }
    /* for (int i = 2; i < nx; i++)
    {
        streamfunction[0][i] = streamfunction[ny-1][i]; // update ghost at y = 0
    }*/
    streamfunction[1] = streamfunction[ny]; // enforce (y = 1) = (y = ny) = 0
    
    R = max_norm(remove_ghost_nodes(remove_ghost_nodes(residual)));
    V = max_norm(remove_ghost_nodes(remove_ghost_nodes(vorticity)));
    S = one_norm(remove_ghost_nodes(remove_ghost_nodes(streamfunction)));
    epsilon = R/((2*((1/pow(dx,2))+(1/pow(dy,2)))*S) + V);
}

vector<vector<double>> Jacobian (vector<vector<double>> streamfunction, vector<vector<double>> vorticity)
{
    double J1, J2, J3;
    vector<vector<double>> J(ny+2, vector<double>(nx+2));
    double A, B;
    for (int j = 2; j < ny; j++)
    {
        for (int i = 2; i < nx; i++)
        {
            J1 = (1.0/(2*dx))*(streamfunction[j][i+1] - streamfunction[j][i-1])*(1.0/(2*dy))*(vorticity[j+1][i] - vorticity[j-1][i]) - (1.0/(2*dy))*(streamfunction[j+1][i] - streamfunction[j-1][i])*(1.0/(2*dx))*(vorticity[j][i+1] - vorticity[j][i-1]);
            
            A = (1.0/(2*dx))*(streamfunction[j][i+1]*(1.0/(2*dy))*(vorticity[j+1][i+1] - vorticity[j-1][i+1]) - streamfunction[j][i-1]*(1.0/(2*dy))*(vorticity[j+1][i-1] - vorticity[j-1][i-1]));
            B = (1.0/(2*dy))*(streamfunction[j+1][i]*(1.0/(2*dx))*(vorticity[j+1][i+1] - vorticity[j+1][i-1]) - streamfunction[j-1][i]*(1.0/(2*dx))*(vorticity[j-1][i+1] - vorticity[j-1][i-1]));
            J2 = A - B;
            
            A = (1.0/(2*dy))*(vorticity[j+1][i]*(1.0/(2*dx))*(streamfunction[j+1][i+1] - streamfunction[j+1][i-1]) - vorticity[j-1][i]*(1.0/(2*dx))*(streamfunction[j-1][i+1] - streamfunction[j-1][i-1]));
            B = (1.0/(2*dx))*(vorticity[j][i+1]*(1.0/(2*dy))*(streamfunction[j+1][i+1] - streamfunction[j-1][i+1]) - vorticity[j][i-1]*(1.0/(2*dy))*(streamfunction[j+1][i-1] - streamfunction[j-1][i-1]));
            J3 = A - B;
            
            J[j][i] = (1.0/3.0)*(J1 + J2 + J3);
        }
    }
    // implement boundary counditions
    for (int index = 0; index < nx; index++)
    {
        J[index][0] = 0.0;
        J[0][index] = 0.0;
        J[index][nx] = 0.0;
        J[ny][index] = 0.0;
    }
    return J;
}

vector<vector<double>> EF (vector<vector<double>> current, vector<vector<double>> J)
{
    // calculate Jacobian (J) at current time-step (n)
    vector<vector<double>> next(ny+2, vector<double>(nx+2));
    for (int j = 2; j < ny; j++)
    {
        for (int i = 2; i < nx; i++)
        {
            next[j][i] = current[j][i] - dt*(J[j][i]);
        }
    }
    return next;
}

vector<vector<double>> EF_coriolis (vector<vector<double>> current, vector<vector<double>> J, vector<vector<double>> SF)
{
    // calculate Jacobian (J) at current time-step (n)
    // calculate Streamfunction (SF) at current time-step (n)
    vector<vector<double>> next(ny+2, vector<double>(nx+2));
    for (int j = 2; j < ny; j++)
    {
        for (int i = 2; i < nx; i++)
        {
            next[j][i] = current[j][i] - dt*(J[j][i]) - beta*(SF[j][i+1] - SF[j][i-1])/(2.0*dx);
        }
    }
    return next;
}

vector<vector<double>> leapfrog (vector<vector<double>> previous, vector<vector<double>> J)
{
    // calculate Jacobian (J) at previous time-step (n-1)
    vector<vector<double>> next(ny+2, vector<double>(nx+2));
    for (int j = 2; j < ny; j++)
    {
        for (int i = 2; i < nx; i++)
        {
            next[j][i] = previous[j][i] - 2.0*dt*(J[j][i]);
        }
    }
    return next;
}

vector<vector<double>> AB3 (vector<vector<double>> current, vector<vector<double>> J, vector<vector<double>> J_1, vector<vector<double>> J_2)
{
    // calculate Jacobian (J) at current time-step (n)
    // calculate Jacobian (J_1) at n-1
    // calculate Jacobian (J_2) at n-2
    vector<vector<double>> next(ny+2, vector<double>(nx+2));
    for (int j = 2; j < ny; j++)
    {
        for (int i = 2; i < nx; i++)
        {
            next[j][i] = current[j][i] - (23.0/12.0)*(J[j][i])*dt + (16.0/12.0)*(J_1[j][i])*dt - (5.0/12.0)*(J_2[j][i])*dt;
        }
    }
    return next;
}

vector<vector<double>> AB3_coriolis (vector<vector<double>> current, vector<vector<double>> J, vector<vector<double>> J_1, vector<vector<double>> J_2, vector<vector<double>> SF, vector<vector<double>> SF_1, vector<vector<double>> SF_2)
{
    // calculate Jacobian (J) at current time-step (n)
    // calculate Jacobian (J_1) at n-1
    // calculate Jacobian (J_2) at n-2
    // calculate Streamfunction (SF) at current time-step (n)
    // calculate Streamfunction (SF_1) at n-1
    // calculate Streamfunction (SF_2) at n-2
    vector<vector<double>> next(ny+2, vector<double>(nx+2));
    for (int j = 2; j < ny; j++)
    {
        for (int i = 2; i < nx; i++)
        {
            next[j][i] = current[j][i] + (23.0/12.0)*((-1.0)*(J[j][i])*dt - beta*(SF[j][i+1] - SF[j][i-1])/(2.0*dx)) - (16.0/12.0)*(((-1.0)*J_1[j][i])*dt - beta*(SF_1[j][i+1] - SF_1[j][i-1])/(2.0*dx)) + (5.0/12.0)*((-1.0)*(J_2[j][i])*dt - beta*(SF_2[j][i+1] - SF_2[j][i-1])/(2.0*dx));
        }
    }
    return next;
}

vector<vector<double>> u_velocity (vector<vector<double>> streamfunction)
{
    vector<vector<double>> u(ny+2, vector<double>(nx+2));
    for (int j = 2; j < ny; j++)
    {
        for (int i = 2; i < nx; i++)
        {
             u[j][i] = (streamfunction[j-1][i] - streamfunction[j+1][i])/(2*dy);
        }
    }
    return u;
}

vector<vector<double>> v_velocity (vector<vector<double>> streamfunction)
{
    vector<vector<double>> v(ny+2, vector<double>(nx+2));
    for (int j = 2; j < ny; j++)
    {
        for (int i = 2; i < nx; i++)
        {
            v[j][i] = (streamfunction[j][i+1] - streamfunction[j][i-1])/(2*dx);
        }
    }
    return v;
}

double calculate_KE (vector<vector<double>> u_velocity, vector<vector<double>> v_velocity)
{
    double KE = 0;
    for (int j = 2; j < ny; j++)
    {
        for (int i = 2; i < nx; i++)
        {
            KE += (1.0/2.0)*(pow(u_velocity[j][i],2)+pow(v_velocity[j][i],2));
        }
    }
    return KE;
}

double calculate_EN (vector<vector<double>> vorticity)
{
    double EN = 0;
    for (int j = 2; j < ny; j++)
    {
        for (int i = 2; i < nx; i++)
        {
            EN += pow(abs(vorticity[j][i]),2);
        }
    }
    return EN;
}

double calculate_Cmax (vector<vector<double>> u_velocity, vector<vector<double>> v_velocity)
{
    vector<double> C;
    for (int j = 2; j < ny; j++)
    {
        for (int i = 2; i < nx; i++)
        {
            C.push_back(abs(u_velocity[j][i])*(dt/dx) + abs(v_velocity[j][i])*(dt/dy));
        }
    }
    return *max_element(C.begin(),C.end());
}

double sum (vector<vector<double>> u)
{
    double sum = 0;
    for (int col = 1; col < u.size()-1; col++)
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
    double max = abs(u[0][0]);
    for (int col = 1; col < u.size()-1; col++)
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
    for (int col = 1; col < u.size()-1; col++)
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
    double max = u[0][0];
    for (int col = 1; col < u.size()-1; col++)
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
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            output << fixed << setprecision(precision) << left << setw(width) << u[j][i];
        }
        output << endl;
    }
}
