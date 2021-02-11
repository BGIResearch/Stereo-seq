/*
 * File: locusFunction.cpp
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#include "locusFunction.h"

LocusFunction getLocusFunction(std::vector< LocusFunction >& locusFunctions, bool conservative)
{
    if (locusFunctions.empty())
        return LocusFunction::INTERGENIC;
    if (conservative)
    {
        LocusFunction f = locusFunctions[0];
        for (LocusFunction& lf : locusFunctions)
            if (lf != f)
                return LocusFunction::NONE;
        return f;
    }

    LocusFunction bestFunction = LocusFunction::INTERGENIC;
    for (LocusFunction& f : locusFunctions)
    {
        bestFunction = std::max(bestFunction, f);
    }
    return bestFunction;
}