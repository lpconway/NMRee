#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

// Combination calculator

long combination(int N, int i)
{
    int a = 1;
    for(int loop = 1; loop < N+1; loop++)
    {
        a = a * loop;
    }

    int b = 1;
    for(int loop = 1; loop < i+1; loop++)
    {
        b = b * loop;
    }

    int c = 1;
    for(int loop = 1; loop < (N-i)+1; loop++)
    {
        c = c * loop;
    }
    return a/(b*c);
}

// Calculates the ratio of deuterium to total acidic solvent hydrogens

double solventfraction(double x[100], double d0, double p0, long N)
{
    double sum = 0;
    for(int i = 0; i < N + 1; i++)
    {
        sum = sum + (i*x[i]);
    }
    return (d0 - sum)/(d0 + p0);
}

// Calculates the expected integral based on the population of each species

double integral(long N, double xu, double x[100])
{
    double sum = 0;
    for(int i = 0; i < N + 1; i++)
    {
        sum = sum + (N - i) * x[i];
    }
    return N * xu + sum;
}

int main()
{
    double x[100]={};       // Fraction of substrate molecules with [index] deuterons
    double dx[100]={};      // Rate of change of fractions x
    double d0;              // Initial equivalents of acidic deuterons in solvent (relative to substrate)
    double p0;              // Initial equivalents of acidic protons in solvent (relative to substrate)
    int N;                  // Equivalents of exchangeable substrate hydrogens (relative to substrate)
    double xu = 1;          // Fraction of unracemised substrate molecules
    long steps = 1000;             // Number of time steps
    double k = 0.04;        // Rate constant (arbitrary units - choose so that the entire reaction is effectively covered)
    double ee0 = 100;       // Initial enantiomeric excess

    double k1[100]={};      // variables for Runge-Kutta method
    double k2[100]={};
    double k3[100]={};
    double k4[100]={};
    double sf;
    double rkx[100]={};

    // Take inputs for variables

    cout << "Exchangeable hydrogens (equiv):";
    cin >> N;
    cout << "\nSolvent deuterons (equiv):";
    cin >> d0;
    cout << "\nSolvent protons (equiv):";
    cin >> p0;
    cout << "\nInitial enantiomeric excess (%):";
    cin >> ee0;

    // Open output file

    ofstream output;
    output.open ("output.csv");

    // Create appropriate column headings for output file

    output << "Time step:,xu:";
    for (int i = 0; i < N + 1; i++)
    {
        output << ",x" << i << ":";
    }

    output << ",Integral:,%ee:\n";

    // Calculation of fractions

    for(long time = 0; time < steps + 1; time++)
    {

        // xu follows a simple exponential decay over time

        xu = exp(-k*time);

        // Output fractions, integral, and enantiomeric excess for this time step to file

        output << time << "," << xu;
        for(int index = 0; index < N + 1; index++)
        {
            output << "," << x[index];
        }

        output << "," << integral(N,xu,x) << "," << xu*ee0 << "\n";

        // Calculate sf

        sf = solventfraction(x, d0, p0, N);

        // Calculate k1

        for(long i = 0; i < N + 1; i++)
        {
            k1[i] = -k * x[i] + k * combination(N, i) * pow(sf, i) * pow((1 - sf), (N - i));
        }

        // Calculate new solvent fraction

        for(long i = 0; i < N + 1; i++)
        {
            rkx[i]=x[i] + 0.5 * k1[i];
        }

        sf = solventfraction(rkx, d0, p0, N);

        // Calculate k2

        for(long i = 0; i < N + 1; i++)
        {
            k2[i] = -k * rkx[i] + k * combination(N, i) * pow(sf, i) * pow((1 - sf), (N - i));
        }

        // Calculate k3 - same solvent fraction

        for(long i = 0; i < N + 1; i++)
        {
            rkx[i]=x[i] + 0.5 * k2[i];
        }

        for(long i = 0; i < N + 1; i++)
        {
            k3[i] = -k * rkx[i] + k * combination(N, i) * pow(sf, i) * pow((1 - sf), (N - i));
        }

        // Calculate new solvent fraction

        for(long i = 0; i < N + 1; i++)
        {
            rkx[i]=x[i] + k3[i];
        }

        sf = solventfraction(rkx, d0, p0, N);

        // Calculate k4

        for(long i = 0; i < N + 1; i++)
        {
            k4[i] = -k * rkx[i] + k * combination(N, i) * pow(sf, i) * pow((1 - sf), (N - i));
        }

        // Calculate new xi values

        for(long i = 0; i < N + 1; i++)
        {
            x[i] = x[i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
        }

    }

    // Close output file

    output.close();
    return 0;

}
