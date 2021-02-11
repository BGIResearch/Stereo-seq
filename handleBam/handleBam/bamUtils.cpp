/*
 * File: bamUtils.cpp
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#include "bamUtils.h"

/*
cigars's op:
BAM_CMATCH      0	M
BAM_CINS        1	I
BAM_CDEL        2	D
BAM_CREF_SKIP   3	N
BAM_CSOFT_CLIP  4	S
BAM_CHARD_CLIP  5	H
BAM_CPAD        6	P
BAM_CEQUAL      7	EQ/=
BAM_CDIFF       8	X
BAM_CBACK       9
*/

enum CIGAR_TYPE
{
    BAM_CMATCH,
    BAM_CINS,
    BAM_CDEL,
    BAM_CREF_SKIP,
    BAM_CSOFT_CLIP,
    BAM_CHARD_CLIP,
    BAM_CPAD,
    BAM_CEQUAL,
    BAM_CDIFF
};

std::vector< AlignmentBlock > getAlignmentBlocks(std::vector< std::pair< int, int > >& cigars, int alignmentStart)
{
    if (cigars.empty())
        return {};
    std::vector< AlignmentBlock > alignmentBlocks;
    int                           readBase = 1;
    int                           refBase  = alignmentStart;
    int                           len;
    for (auto& c : cigars)
    {
        switch (c.first)
        {
        case CIGAR_TYPE::BAM_CHARD_CLIP:
            break;  // ignore hard clips
        case CIGAR_TYPE::BAM_CPAD:
            break;  // ignore pads
        case CIGAR_TYPE::BAM_CSOFT_CLIP:
            readBase += c.second;
            break;  // soft clip read bases
        case CIGAR_TYPE::BAM_CREF_SKIP:
            refBase += c.second;
            break;  // reference skip
        case CIGAR_TYPE::BAM_CDEL:
            refBase += c.second;
            break;
        case CIGAR_TYPE::BAM_CINS:
            readBase += c.second;
            break;
        case CIGAR_TYPE::BAM_CMATCH:
        case CIGAR_TYPE::BAM_CEQUAL:
        case CIGAR_TYPE::BAM_CDIFF:
            len = c.second;
            alignmentBlocks.push_back(AlignmentBlock(readBase, refBase, len));
            readBase += len;
            refBase += len;
            break;
        default:
            break;
            // throw std::exception("Case statement didn't deal op: " + c.operation);
        }
    }
    return alignmentBlocks;
}

int getReferenceLength(std::vector< std::pair< int, int > >& cigars)
{
    int length = 0;
    for (auto& c : cigars)
    {
        switch (c.first)
        {
        case CIGAR_TYPE::BAM_CMATCH:
        case CIGAR_TYPE::BAM_CDEL:
        case CIGAR_TYPE::BAM_CREF_SKIP:
        case CIGAR_TYPE::BAM_CEQUAL:
        case CIGAR_TYPE::BAM_CDIFF:
            length += c.second;
            break;
        default:
            break;
        }
    }
    return length;
}