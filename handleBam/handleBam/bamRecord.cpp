/*
 * File: bamRecord.cpp
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#include <algorithm>
#include <exception>
#include <iostream>
#include <string>

#include "bamRecord.h"

int getQName(BamRecord b, std::string& qname)
{
    qname = bam_get_qname(b);
    return 0;
}

int getQual(BamRecord b)
{
    return b->core.qual;
}

std::vector< std::pair< int, int > > getCigar(BamRecord b)
{
    std::vector< std::pair< int, int > > cigars(b->core.n_cigar);

    uint32_t* cigar = bam_get_cigar(b);
    for (uint32_t i = 0; i < b->core.n_cigar; ++i)
    {
        int op    = bam_cigar_op(cigar[i]);
        int len   = bam_cigar_oplen(cigar[i]);
        cigars[i] = { op, len };
    }
    return cigars;
}

int getRefStart(BamRecord b)
{
    return b->core.pos;
}

bool getNegativeStrand(BamRecord b)
{
    return bam_is_rev(b);
}

int updateStrTags(BamRecord& b, std::string tag, std::string data)
{
    return bam_aux_update_str(b, tag.c_str(), data.size() + 1, data.c_str());
}

BamRecord createBamRecord()
{
    return bam_init1();
}

int destroyBamRecord(BamRecord b)
{
    bam_destroy1(b);
    return 0;
}

bool compareBamRecord(BamRecord b1, BamRecord b2)
{
    // if (getContig(b1) != getContig(b2))
    //     return false;
    if (getRefStart(b1) != getRefStart(b2))
        return false;
    if (getNegativeStrand(b1) != getNegativeStrand(b2))
        return false;
    if (b1->core.n_cigar != b2->core.n_cigar)
        return false;
    else
    {
        uint32_t* cigar1 = bam_get_cigar(b1);
        uint32_t* cigar2 = bam_get_cigar(b2);
        for (uint32_t i = 0; i < b1->core.n_cigar; ++i)
        {
            int op1  = bam_cigar_op(cigar1[i]);
            int len1 = bam_cigar_oplen(cigar1[i]);
            int op2  = bam_cigar_op(cigar2[i]);
            int len2 = bam_cigar_oplen(cigar2[i]);
            if (op1 != op2 || len1 != len2)
                return false;
        }
    }
    return true;
}

// Used for deduplication when no umi found
int getMarker(BamRecord b, std::string& marker)
{
    if (b->core.isize < 0)
    {
        marker = std::to_string(b->core.mpos - b->core.isize) + std::to_string(b->core.mpos);
    }
    else
    {
        marker = std::to_string(b->core.pos) + std::to_string(b->core.pos + b->core.isize);
    }
    return 0;
}

bool getTag(BamRecord b, const char tag[2], std::string& value)
{
    uint8_t* data = bam_aux_get(b, tag);
    if (data == NULL)
        return false;
    value = bam_aux2Z(data);
    return true;
}

bool getTagInt(BamRecord b, const char tag[2], int& value)
{
    uint8_t* data = bam_aux_get(b, tag);
    if (data == NULL)
        return false;
    value = bam_aux2i(data);
    return true;
}

void setQcFail(BamRecord b)
{
    b->core.flag |= FQCFAIL;
}

void setDuplication(BamRecord b)
{
    b->core.flag |= FDUP;
}

bool getQcFail(BamRecord b)
{
    return b->core.flag & FQCFAIL;
}

bool getDuplication(BamRecord b)
{
    return b->core.flag & FDUP;
}