/*
 * File: saturation.cpp
 * Created Data: 2020-6-16
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#include "saturation.h"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <fstream>
#include <limits>
#include <random>
#include <set>
#include <sstream>

#include "utils.h"

Saturation::Saturation()
{
    _ge2i.clear();
    _i = 0;

    _sep         = '|';
    _barcode_sep = '_';
    _bin         = 150;

    _nreads = 0;

    _samples = { 0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1 };

    _base2i['A'] = 0;
    _base2i['C'] = 1;
    _base2i['G'] = 2;
    _base2i['T'] = 3;
}

Saturation::~Saturation() {}

int Saturation::calculateSaturation(string& out_file)
{
    std::lock_guard< std::mutex > guard(_mutex);

    // Get NOGENE index
    nogene_idx = std::numeric_limits< unsigned int >::max();
    spdlog::debug("NOGENE index: {}", nogene_idx);
    if (_ge2i.count("NOGENE"))
    {
        nogene_idx = _ge2i["NOGENE"];
        spdlog::debug("NOGENE index: {}", nogene_idx);
    }

    spdlog::info("Calculate sequencing saturation");
    spdlog::debug("_nreads {}", _nreads);
    spdlog::info("physical memory used before calculating saturation:{}(GB)",
                 physical_memory_used_by_process() / (1024 * 1024));

    string sample_result = sample();

    spdlog::info("maximum physical memory used when calculating saturation:{}(GB)",
                 physical_memory_used_by_process() / (1024 * 1024));

    std::ofstream ofs(out_file, std::ofstream::out);
    if (ofs.is_open())
    {
        ofs << sample_result;
        ofs.close();
        spdlog::info("Success dump saturation file:{}", out_file);
    }
    else
    {
        spdlog::error("Error opening file:{}", out_file);
    }

    return 0;
}

unsigned int Saturation::encodeUmi(string_view sv)
{
    unsigned int res = 0;
    for (auto& s : sv)
    {
        res = (res * 4) + _base2i[( unsigned char )s];
    }
    return res;
}

unsigned int Saturation::encodeGene(string& gene)
{
    if (_ge2i.count(gene) == 0)
    {
        _ge2i[gene] = _i++;
    }
    return _ge2i[gene];
}

int CoordinateBarcode::addData(std::unordered_map< std::string, std::unordered_map< std::string, int > >& raw)
{
    std::lock_guard< std::mutex > guard(_mutex);

    for (const auto& [b, p] : raw)
    {
        for (const auto& [umi, count] : p)
        {
            if (count == 0)
                continue;

            string_view tmp(b);
            size_t      pos0    = tmp.find(_sep);
            string_view barcode = tmp.substr(0, pos0);
            string      gene    = string(tmp.substr(pos0 + 1));

            pos0            = barcode.find_last_of(_barcode_sep);
            string_view t   = barcode.substr(0, pos0);
            int         col = stoi(string(barcode.substr(pos0 + 1)));
            pos0            = t.find_last_of(_barcode_sep);
            int row         = stoi(string(t.substr(pos0 + 1)));

            ST st;
            st.b1  = col;
            st.b2  = row;
            st.umi = encodeUmi(string_view(umi));
            st.ge  = encodeGene(gene);
            _st.insert(_st.end(), count, st);

            _nreads += count;
        }
    }

    return 0;
}

template < class T1, class T2 > Metrics Saturation::saturation(unordered_map< T1, unordered_map< T2, int > >& data)
{
    Metrics                  metrics;
    size_t                   n_reads     = 0;
    size_t                   n_reads_sat = 0;
    size_t                   n_uniq      = 0;
    vector< int >            n_genes;
    std::set< unsigned int > genes;
    for (auto& b : data)
    {
        genes.clear();
        for (auto& p : b.second)
        {
            n_reads += p.second;
            unsigned int gene = ( unsigned int )(p.first >> 32);
            if (gene != nogene_idx)
            {
                genes.insert(gene);
                ++n_uniq;
                n_reads_sat += p.second;
            }
        }
        n_genes.push_back(genes.size());
    }
    if (n_reads == 0)
    {
        return metrics;
    }

    std::nth_element(n_genes.begin(), n_genes.begin() + n_genes.size() / 2, n_genes.end());

    metrics.reads        = n_reads;
    metrics.reads_sat    = n_reads_sat;
    metrics.uniq         = n_uniq;
    metrics.median_genes = *(n_genes.begin() + n_genes.size() / 2);

    return metrics;
}

string CoordinateBarcode::sample()
{
    std::stringstream ss;
    ss << "#sample bar_x bar_y1 bar_y2 bin_x bin_y1 bin_y2\n";

    std::random_device rd;
    std::mt19937       gen(rd());
    std::shuffle(_st.begin(), _st.end(), gen);

    unordered_map< unsigned long long, unordered_map< unsigned long long, int > > data;
    unordered_map< unsigned int, unordered_map< unsigned long long, int > >       data_bin;
    size_t                                                                        p = 0;
    for (size_t i = 1; i < _samples.size(); ++i)
    {
        spdlog::info("Saturation sample:{}", _samples[i]);
        ss << _samples[i] << " ";

        size_t size = size_t(_samples[i] * _nreads);
        for (; p < size; ++p)
        {
            ST                 st      = _st[p];
            unsigned long long barcode = (( unsigned long long )st.b1 << 32) + st.b2;
            if (data.count(barcode) == 0)
                data[barcode] = {};
            unsigned long long value = (( unsigned long long )st.ge << 32) + st.umi;
            ++data[barcode][value];

            unsigned int col         = st.b1 / _bin;
            unsigned int row         = st.b2 / _bin;
            unsigned int new_barcode = (row << 16) + col;
            if (data_bin.count(new_barcode) == 0)
                data_bin[new_barcode] = {};
            ++data_bin[new_barcode][value];
        }

        // Barcode
        Metrics metrics;
        metrics = saturation(data);
        if (metrics.reads == 0)
        {
            spdlog::warn("Invalid data: n_reads == 0");
            continue;
        }
        ss << metrics.reads / data.size() << " " << 1 - (metrics.uniq * 1.0 / metrics.reads_sat) << " "
           << metrics.median_genes << " ";

        // Bin
        metrics = saturation(data_bin);
        ss << metrics.reads / data_bin.size() << " " << 1 - (metrics.uniq * 1.0 / metrics.reads_sat) << " "
           << metrics.median_genes << std::endl;
    }

    return ss.str();
}

int SequenceBarcode::addData(std::unordered_map< std::string, std::unordered_map< std::string, int > >& raw)
{
    std::lock_guard< std::mutex > guard(_mutex);

    for (const auto& [b, p] : raw)
    {
        for (const auto& [umi, count] : p)
        {
            if (count == 0)
                continue;

            string_view tmp(b);
            size_t      pos0    = tmp.find(_sep);
            string      barcode = string(tmp.substr(0, pos0));
            string      gene    = string(tmp.substr(pos0 + 1));

            ST st;
            st.bar = barcode;
            st.umi = encodeUmi(string_view(umi));
            st.ge  = encodeGene(gene);
            _st.insert(_st.end(), count, st);

            _nreads += count;
        }
    }

    return 0;
}

string SequenceBarcode::sample()
{
    std::stringstream ss;
    ss << "#sample bar_x bar_y1 bar_y2\n";

    std::random_device rd;
    std::mt19937       gen(rd());
    std::shuffle(_st.begin(), _st.end(), gen);

    unordered_map< string, unordered_map< unsigned long long, int > > data;
    size_t                                                            p = 0;
    for (size_t i = 1; i < _samples.size(); ++i)
    {
        spdlog::info("Saturation sample:{}", _samples[i]);
        ss << _samples[i] << " ";

        size_t size = size_t(_samples[i] * _nreads);
        for (; p < size; ++p)
        {
            ST     st      = _st[p];
            string barcode = st.bar;
            if (data.count(barcode) == 0)
                data[barcode] = {};
            unsigned long long value = (( unsigned long long )st.ge << 32) + st.umi;
            ++data[barcode][value];
        }

        // Barcode
        Metrics metrics;
        metrics = saturation(data);
        if (metrics.reads == 0)
        {
            spdlog::warn("Invalid data: n_reads == 0");
            continue;
        }
        ss << metrics.reads / data.size() << " " << 1 - (metrics.uniq * 1.0 / metrics.reads_sat) << " "
           << metrics.median_genes << endl;
    }

    return ss.str();
}
