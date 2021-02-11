/*
 * File: gtfRecord.h
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#pragma once

#include <string>

#include "interval.h"

class GTFRecord
{
public:
    GTFRecord() : geneVersion(-1), transcriptID(""), transcriptName(""), transcriptType("") {}

    std::string& getGeneName()
    {
        return geneName;
    }
    std::string& getFeatureType()
    {
        return featureType;
    }
    int getStart()
    {
        return interval.getStart();
    }
    int getEnd()
    {
        return interval.getEnd();
    }
    std::string& getContig()
    {
        return interval.getContig();
    }
    std::string& getGeneID()
    {
        return geneID;
    }
    int getGeneVersion()
    {
        return geneVersion;
    }
    bool isNegativeStrand()
    {
        return interval.isNegativeStrand();
    }
    std::string& getTranscriptType()
    {
        return transcriptType;
    }
    std::string& getTranscriptID()
    {
        return transcriptID;
    }
    std::string& getTranscriptName()
    {
        return transcriptName;
    }

    void setGeneName(std::string& geneName_)
    {
        geneName = geneName_;
    }
    void setGeneID(std::string& geneID_)
    {
        geneID = geneID_;
    }
    void setTranscriptName(std::string& transcriptName_)
    {
        transcriptName = transcriptName_;
    }
    void setTranscriptID(std::string& transcriptID_)
    {
        transcriptID = transcriptID_;
    }
    void setTranscriptType(std::string& transcriptType_)
    {
        transcriptType = transcriptType_;
    }
    void setFeatureType(std::string& featureType_)
    {
        featureType = featureType_;
    }
    void setInterval(int start, int end, std::string& contig, bool negativeStrand,
                     [[maybe_unused]] std::string name = "")
    {
        interval.setStart(start);
        interval.setEnd(end);
        interval.setContig(contig);
        interval.setNegativeStrand(negativeStrand);
    }

private:
    Interval    interval;
    std::string geneID;
    std::string geneName;
    std::string featureType;
    int         geneVersion;

    std::string transcriptID;
    std::string transcriptName;
    std::string transcriptType;
};
