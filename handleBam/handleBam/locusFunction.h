/*
 * File: locusFunction.h
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#pragma once

#include <algorithm>
#include <string>
#include <vector>

enum LocusFunction
{
    NONE,
    INTERGENIC,
    RIBOSOMAL,
    INTRONIC,
    UTR,
    CODING
};

static std::string LocusString[] = { "", "INTERGENIC", "RIBOSOMAL", "INTRONIC", "UTR", "EXONIC" };

/**
 * Summarize the locus functions that are at base-by-base level to a single annotation.
 *
 * @param locusFunctions
 * @param conservative If true, only return a LocusFunction if all locusFunctions are the same, otherwise return null.
 * If false, return the "best" annotation, where annotations like coding are preferred over intronic,intergenic.
 * @return locusFunction
 */
LocusFunction getLocusFunction(std::vector< LocusFunction >& locusFunctions, bool conservative);