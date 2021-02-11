/*
 * File: geneFromGTF.cpp
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#include <spdlog/spdlog.h>

#include "annotationException.h"
#include "geneFromGTF.h"

bool compareExon(Exon& e1, Exon& e2)
{
    return (e1.start < e2.start);
}

TranscriptFromGTF* GeneFromGTF::addTranscript(int transcriptionStart, int transcriptionEnd, int codingStart,
                                              int codingEnd, int numExons, std::string transcriptName,
                                              std::string transcriptID, std::string transcriptType)
{
    if (transcripts.count(transcriptName) != 0)
    {
        std::string error = "Transcript " + transcriptName + " for gene " + geneName + " appears more than once";
        throw AnnotationException(error);
    }
    else
    {
        {
            TranscriptFromGTF transcript(transcriptionStart, transcriptionEnd, codingStart, codingEnd, numExons,
                                         transcriptName, transcriptID, transcriptType);
            transcripts[transcriptName] = transcript;
            return &(transcripts[transcriptName]);
        }
    }
}

void TranscriptFromGTF::addExons(std::vector< Exon >& exons_)
{
    // spdlog::debug("exons num:{} exons_ num:{}", exons.size(), exons_.size());
    exons.assign(exons_.begin(), exons_.end());
    // spdlog::debug("exons num:{} exons_ num:{}", exons.size(), exons_.size());
    for (auto& e : exons)
        length += e.end - e.start + 1;
}

void TranscriptFromGTF::assignLocusFunction(int start, std::vector< LocusFunction >& locusFunctions) const
{
    int begin = std::max(start, transcriptionStart);
    int end   = std::min(int(start + locusFunctions.size() - 1), transcriptionEnd);
    // spdlog::debug("assignLocusFunction:{} {}", begin, end);
    // spdlog::debug("exons num:{}", exons.size());
    // for (auto& exon : exons)
    //     spdlog::debug("exons :{} {}", exon.start, exon.end);
    // spdlog::debug("start {} end {} coding {} {}", start, end, codingStart, codingEnd);
    for (int i = begin; i <= end; ++i)
    {
        if (locusFunctions[i - start] >= LocusFunction::CODING)
            break;

        LocusFunction locusFunction;
        if (inExon(i))
        {
            if (!inRange(codingStart, codingEnd, i))
                locusFunction = LocusFunction::UTR;
            else
                locusFunction = LocusFunction::CODING;
        }
        else
            locusFunction = LocusFunction::INTRONIC;
        if (locusFunction > locusFunctions[i - start])
            locusFunctions[i - start] = locusFunction;
        // Exit early to reduce unnecessary calculations
        if (locusFunction == LocusFunction::CODING)
            break;
    }
}

void TranscriptFromGTF::assignLocusFunction(int start, LocusFunction& l, int len) const
{
    int begin = std::max(start, transcriptionStart);
    int end   = std::min(int(start + len - 1), transcriptionEnd);
    for (int i = begin; i <= end; ++i)
    {
        LocusFunction locusFunction;
        if (inExon(i))
        {
            if (!inRange(codingStart, codingEnd, i))
                locusFunction = LocusFunction::UTR;
            else
                locusFunction = LocusFunction::CODING;
        }
        else
            locusFunction = LocusFunction::INTRONIC;
        if (locusFunction > l)
            l = locusFunction;
        // Exit early to reduce unnecessary calculations
        if (l == LocusFunction::CODING)
            break;
    }
}

void TranscriptFromGTF::assignLocusFunction(int start, int len, pair< int, int >& max_cnts) const
{
    int begin     = std::max(start, transcriptionStart);
    int end       = std::min(int(start + len - 1), transcriptionEnd);
    int exon_cnts = 0, intro_cnts = 0;
    for (int i = begin; i <= end; ++i)
    {
        if (inExon(i))
            ++exon_cnts;
        else
            ++intro_cnts;
    }
    if (exon_cnts > max_cnts.first)
        max_cnts = { exon_cnts, intro_cnts };
    else if (exon_cnts == max_cnts.first)
        max_cnts.second = max(max_cnts.second, intro_cnts);
}

bool TranscriptFromGTF::inExon(int locus) const
{
    for (auto& exon : exons)
    {
        if (exon.start > locus)
            return false;
        if (inRange(exon.start, exon.end, locus))
            return true;
    }
    return false;
}

inline bool TranscriptFromGTF::inRange(int start, int end, int locus) const
{
    return locus >= start && locus <= end;
}
