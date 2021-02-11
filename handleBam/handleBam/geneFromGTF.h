/*
 * File: geneFromGTF.h
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "bamUtils.h"
#include "locusFunction.h"

using namespace std;

class Exon
{
public:
    Exon(int s, int e) : start(s), end(e) {}

    int start;
    int end;
};

bool compareExon(Exon& e1, Exon& e2);

class TranscriptFromGTF
{
public:
    TranscriptFromGTF(int transcriptionStart_, int transcriptionEnd_, int codingStart_, int codingEnd_, int numExons_,
                      std::string transcriptName_, std::string transcriptID_, std::string transcriptType_)
        : transcriptionStart(transcriptionStart_), transcriptionEnd(transcriptionEnd_), codingStart(codingStart_),
          codingEnd(codingEnd_), numExons(numExons_), transcriptName(transcriptName_), transcriptID(transcriptID_),
          transcriptType(transcriptType_)
    {
        // exons.resize(numExons);
    }
    TranscriptFromGTF() {}

    void addExons(std::vector< Exon >& exons_);

    int getTranscriptionStart()
    {
        return transcriptionStart;
    }
    int getTranscriptionEnd()
    {
        return transcriptionEnd;
    }
    int getCodingStart()
    {
        return codingStart;
    }
    int getCodingEnd()
    {
        return codingEnd;
    }
    const std::vector< Exon >& getExons() const
    {
        return exons;
    }

    void assignLocusFunction(int start, std::vector< LocusFunction >& locusFunctions) const;
    void assignLocusFunction(int start, LocusFunction& l, int len) const;
    void assignLocusFunction(int start, int len, pair< int, int >& max_cnts) const;
    bool inExon(int locus) const;
    bool inRange(int start, int end, int locus) const;

private:
    int transcriptionStart;
    int transcriptionEnd;
    int codingStart;
    int codingEnd;

    int numExons;
    int length;  // the number of bases in the transcript

    std::string transcriptName;
    std::string transcriptID;
    std::string transcriptType;

    std::vector< Exon > exons;
};

class GeneFromGTF
{
public:
    GeneFromGTF(std::string& contig_, int start_, int end_, bool isNegativeStrand_, std::string& geneName_,
                std::string featureType_, std::string geneID_, std::string transcriptType_, int geneVersion_)
        : contig(contig_), start(start_), end(end_), strand(isNegativeStrand_), geneName(geneName_),
          featureType(featureType_), geneID(geneID_), transcriptType(transcriptType_), geneVersion(geneVersion_)
    {
    }

    GeneFromGTF() {}

    bool isNegativeStrand() const
    {
        return strand;
    }
    const std::string& getContig() const
    {
        return contig;
    }
    const std::string& getName() const
    {
        return geneName;
    }
    int getStart()
    {
        return start;
    }
    int getEnd()
    {
        return end;
    }
    const std::unordered_map< std::string, TranscriptFromGTF >& getTranscripts() const
    {
        return transcripts;
    }

    TranscriptFromGTF* addTranscript(int transcriptionStart, int transcriptionEnd, int codingStart, int codingEnd,
                                     int numExons, std::string transcriptName, std::string transcriptID,
                                     std::string transcriptType);

private:
    std::string contig;

    // Construct interval.
    int         start;
    int         end;
    bool        strand;
    std::string geneName;

    std::string featureType;
    std::string geneID;
    std::string transcriptType;
    int         geneVersion;

    std::unordered_map< std::string, TranscriptFromGTF > transcripts;
};