/*
 * File: saturation.h
 * Created Data: 2020-6-16
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#pragma once

#include <mutex>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
using std::string;
using std::string_view;
using std::unordered_map;
using std::vector;

struct Metrics
{
    size_t reads;
    size_t reads_sat;
    size_t uniq;
    size_t median_genes;
};

class Saturation
{
public:
    Saturation();
    virtual ~Saturation();

    // Parse raw data to vectors
    virtual int addData([[maybe_unused]] std::unordered_map< std::string, std::unordered_map< std::string, int > >& raw)
    {
        return 0;
    }

    // Calculate sequencing saturation
    virtual int calculateSaturation(string& out_file);

    // Encode umi sequence to uint32
    unsigned int encodeUmi(string_view sv);

    // Encode gene name to uint32
    unsigned int encodeGene(string& gene);

    virtual string sample()
    {
        return "";
    }

    template < class T1, class T2 > Metrics saturation(unordered_map< T1, unordered_map< T2, int > >& data);

public:
    unordered_map< string, unsigned int > _ge2i;
    unsigned int                          _i;
    std::mutex                            _mutex;
    int                                   _base2i[128];

    size_t _nreads;

    char _sep;
    char _barcode_sep;
    int  _bin;

    vector< float > _samples;

    unsigned int nogene_idx;
};

class CoordinateBarcode : public Saturation
{
public:
    virtual int addData(std::unordered_map< std::string, std::unordered_map< std::string, int > >& raw);

    virtual string sample();

    struct ST
    {
        unsigned int b1;
        unsigned int b2;
        unsigned int ge;
        unsigned int umi;
    };
    vector< ST > _st;
};

class SequenceBarcode : public Saturation
{
public:
    virtual int addData(std::unordered_map< std::string, std::unordered_map< std::string, int > >& raw);

    virtual string sample();

    struct ST
    {
        string       bar;
        unsigned int ge;
        unsigned int umi;
    };
    vector< ST > _st;
};