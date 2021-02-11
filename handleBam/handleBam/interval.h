/*
 * File: interval.h
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#pragma once

#include <string>

class Interval
{
public:
    Interval() {}
    int getStart()
    {
        return start;
    }
    int getEnd()
    {
        return end;
    }
    std::string& getContig()
    {
        return contig;
    }
    bool isNegativeStrand()
    {
        return negativeStrand;
    }

    void setContig(std::string& contig_)
    {
        contig = contig_;
    }
    void setStart(int start_)
    {
        start = start_;
    }
    void setEnd(int end_)
    {
        end = end_;
    }
    void setNegativeStrand(bool negativeStrand_)
    {
        negativeStrand = negativeStrand_;
    }

private:
    std::string contig;
    int         start;
    int         end;
    bool        negativeStrand;
    std::string name;
};
