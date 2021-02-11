// Kernel density estimation by Tim Nugent (c) 2014
// Based on Philipp K. Janert's Perl module:
// http://search.cpan.org/~janert/Statistics-KernelEstimation-0.05
// Multivariate stuff from here:
// http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/ebooks/html/spm/spmhtmlnode18.html

#define _USE_MATH_DEFINES
#include "kde.hpp"
#include "fftw3.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <string>  // can delete

#include <spdlog/spdlog.h>

void KDE::initialization()
{
    Count  = 0;
    bw     = 0.1;
    n_user = 10000;
    n      = pow(2, ceil(log2(n_user)));
    kords.resize(n);
    xords.resize(n);
    density.resize(n_user);
    vec_x.resize(n_user);
}

void KDE::readData()
{
    string fname;
    int    data;
    fname = "testdata.csv";
    ifstream readfile(fname);
    while (readfile >> data)
    {  // assume the input has been preprocessed to eliminate NA or infinite items
        data_array.push_back(log10(
            data));  // we log10-transform the data (detransform would only take place after this module to save time)
    }
    readfile.close();
    N = data_array.size();
    filter();
}

void KDE::fft()
{
    int           i, Num = 2 * n;
    double        xlo  = min - 4 * bw;
    double        xhi  = max + 4 * bw;
    double        diff = 2 * (xhi - xlo) / (Num - 1);
    double        Temp;
    fftw_complex* in;
    fftw_complex *temp_out, *y_fft, *kords_fft;
    fftw_plan     p;
    // in        = ( fftw_complex* )fftw_malloc(sizeof(fftw_complex) * Num);
    y_fft     = ( fftw_complex* )fftw_malloc(sizeof(fftw_complex) * Num);
    kords_fft = ( fftw_complex* )fftw_malloc(sizeof(fftw_complex) * Num);
    temp_out  = ( fftw_complex* )fftw_malloc(sizeof(fftw_complex) * Num);

    in = bindist();  // "in" is actually in_y here.

    p = fftw_plan_dft_1d(Num, in, temp_out, FFTW_FORWARD,
                         FFTW_ESTIMATE);  // I cannot put it before "in = bindist()" but anyway.
    fftw_execute(p);                      // Calculate fft(y).

    memcpy(y_fft, temp_out, sizeof(fftw_complex) * Num);

    for (i = 0; i < n + 1; i++)
    {  // let in be initial in_kords
        Temp     = i * diff;
        in[i][0] = gauss_pdf(Temp);
        in[i][1] = 0;
    }
    for (i = n + 1; i < Num; i++)
    {
        in[i][0] = in[Num - i][0];
        in[i][1] = 0;
    }
    fftw_execute_dft(p, in, temp_out);  // calculate fft(kords)

    for (i = 0; i < Num; i++)
    {
        kords_fft[i][0] = temp_out[i][0] * y_fft[i][0] + temp_out[i][1] * y_fft[i][1];
        kords_fft[i][1] = temp_out[i][0] * y_fft[i][1]
                          - temp_out[i][1] * y_fft[i][0];  // let kords_fft equals to "fft(y)*Conj(fft(kords))"
    }

    p = fftw_plan_dft_1d(Num, kords_fft, temp_out, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(p);

    double diff_new = (xhi - xlo) / (n - 1);
    for (i = 0; i < n; i++)
    {
        Temp     = temp_out[i][0] / Num;
        kords[i] = 0 > Temp ? 0 : Temp;
        xords[i] = xlo + diff_new * i;
    }
    /*
    ofstream saveFile("testresult.txt");
    for (i = 0; i < Num; i++) {
        saveFile << kords[i] << endl;
    }
    saveFile.close();
    */
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(y_fft);
    fftw_free(temp_out);
    fftw_free(kords_fft);
}

void KDE::writeData()
{
    int      i = 0;
    ofstream saveFile("result.txt");
    double   x_increment = (max - min) / (n_user - 1);
    for (double x = min; x <= max; x += x_increment)
    {
        density[i] = pdf_linear_interpol(x);
        vec_x[i]   = x;
        saveFile << x << ' ' << density[i] << endl;
        i++;
    }
    saveFile.close();
}

pair< double, double > KDE::get_density_threshold(string type)
{
    auto minima = find_local_minima();

    int    x, local_min = 0;
    double abs_difference = 0.5;
    for (int i = minima.size() - 1; i >= 0; i--)
    {
        x = minima[i];
        if ((x >= 0.2 * n_user) & ((max - vec_x[x]) > abs_difference || (vec_x[x] < max / 2)))
        {
            local_min = x;
            break;
        }
    }

    double threshold = 0;
    if (local_min != 0)
    {
        threshold = pow(10, vec_x[local_min]);
    }
    else
    {
        threshold = 0;
    }

    double safety = 0;
    if (type == "bead" && (threshold > 100000 || threshold < 100))
    {
        safety = 500;
    }
    if (type == "jaccard" && (threshold > 0.5 || threshold < 0.000001))
    {
        safety = 0.005;
    }
    safety = safety > 0 ? safety : threshold;

    return { safety, threshold };
}

void KDE::filter()
{
    double threshold = get_min_mode() - 3;  // we have initially done the log10 transform, so '* 0.001' becomes '- 3'.
    while (!data_array.empty())
    {
        if (data_array.back() <= threshold)
        {
            data_array.pop_back();
        }
        else
        {
            break;
        }
    }
    min = data_array.back();
    max = data_array.front();
    //	cout << threshold << ' ' << min << ' ' << max << endl;
}

/*
int compare(const void * a, const void * b) {
    if (*(double*)a > *(double*)b) return -1;
    else if (*(double*)a < *(double*)b) return 1;
    else return 0;
}
*/

double KDE::get_min_mode()
{
    std::sort(data_array.begin(), data_array.end(), std::greater< double >());
    // cout << N << ' ' << data_array.size() << endl;
    int    count    = 1;
    int    maxcount = 1;
    double min_mode = data_array.back();
    double temp_var = min_mode;
    for (int i = data_array.size() - 2; i >= 0; i--)
    {
        double curr_var = data_array[i];
        if (curr_var == temp_var)
        {
            count += 1;
            if (maxcount <= count)
            {
                maxcount = count;
                min_mode = curr_var;
            }
        }
        else
            count = 1;
        temp_var = curr_var;
    }
    //	cout << min_mode << endl;
    return min_mode;
}

fftw_complex* KDE::bindist()
{
    double           w = 1.0 / N;
    vector< double > x = data_array;
    fftw_complex*    bindens;
    bindens      = ( fftw_complex* )fftw_malloc(sizeof(fftw_complex) * 2 * n);
    double xlo   = min - 4 * bw;
    double xhi   = max + 4 * bw;
    int    ixmin = 0, ixmax = n - 2;
    double xdelta = (xhi - xlo) / (n - 1);
    for (int i = 0; i < 2 * n; i++)
    {  // 2n=16384*2
        bindens[i][0] = 0;
        bindens[i][1] = 0;
    }
    for (int i = 0; i < N; i++)
    {
        double xpos = (x[i] - xlo) / xdelta;
        int    ix   = ( int )floor(xpos);
        double fx   = xpos - ix;
        if (ixmin <= ix && ix <= ixmax)
        {
            bindens[ix][0] += (1 - fx) * w;
            bindens[ix + 1][0] += fx * w;
        }
        //		else if (ix == -1) bindens[0] += fx * wi;  // since xlo < ixmin we don't need this now
        else if (ix == ixmax + 1)
            bindens[ix][0] += (1 - fx) * w;
    }

    return bindens;
}

double KDE::gauss_pdf(double x)
{
    double              m = 0, s = bw;
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double              z            = (x - m) / s;
    return exp(-0.5 * z * z) / s * inv_sqrt_2pi;
}

double KDE::pdf_linear_interpol(double v)
{
    double diff  = xords[1] - xords[0];
    double xlo   = min - 4 * bw;
    int    index = round((v - xlo) / diff);
    double d     = kords[index - 1]
               + (kords[index] - kords[index - 1]) * (v - xords[index - 1]) / (xords[index] - xords[index - 1]);
    return d;
}

vector< int > KDE::find_local_minima()
{
    vector< int > minima_temp;

    int flag = 0;
    if (density[1] - density[0] > 0)
        flag = 1;
    for (int i = 1; i < n_user - 1; i++)
    {
        if ((density[i] - density[i - 1]) * (density[i + 1] - density[i]) < 0)
        {
            minima_temp.push_back(i);
        }
    }

    vector< int > minima;
    if (minima_temp.size() > 2)
    {
        for (int i = 0; i < minima_temp.size() / 2; i++)
        {
            minima.push_back(minima_temp[2 * i + flag]);
        }
    }
    else
    {
        minima = minima_temp;
    }

    return minima;
}

pair< double, double > KDE::run(vector< double >& input, string type)
{
    // Prepare data and filter
    initialization();

    data_array.resize(input.size());
    std::transform(input.begin(), input.end(), data_array.begin(), [](double d) { return log10(d); });
    N = data_array.size();

    filter();

    // Do the FFT
    fft();

    // Get density threshold
    int    i           = 0;
    double x_increment = (max - min) / (n_user - 1);
    for (double x = min; x <= max; x += x_increment)
    {
        density[i] = pdf_linear_interpol(x);
        vec_x[i]   = x;
        i++;
    }

    auto res = get_density_threshold(type);

    return res;
}