// Kernel density estimation by Tim Nugent (c) 2014

#ifndef KDE_HPP
#define KDE_HPP

#include "fftw3.h"
#include <cmath>
#include <iostream>
#include <map>
#include <stdint.h>
#include <stdlib.h>
#include <vector>

using namespace std;

class KDE
{

public:
    KDE() : extension(4){};
    ~KDE(){};
    void                   readData();
    void                   writeData();
    void                   initialization();
    void                   fft();
    pair< double, double > get_density_threshold(string type);
    pair< double, double > run(vector< double >& input, string type);

private:
    double        gauss_pdf(double x);
    void          filter();
    fftw_complex* bindist();
    vector< int > find_local_minima();
    double        pdf_linear_interpol(double v);
    double        get_min_mode();

private:
    double           min, max;
    vector< double > data_array;
    vector< double > kords, xords;
    vector< double > density, vec_x;
    double           bw;
    int              Count;
    int              N, n_user, n;
    unsigned int     extension, curr_var;
};

#endif