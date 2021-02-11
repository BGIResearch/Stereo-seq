/*
 * File: geneBuilder.h
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "geneFromGTF.h"
#include "gtfRecord.h"

typedef std::unordered_map< std::string, std::vector< GTFRecord > >::const_iterator GTFIterator;
typedef std::unordered_map< int, std::vector< GTFRecord > >                         GTFMapForGeneVersion;
typedef std::unordered_map< std::string, std::vector< GTFRecord > >                 GTFMapForTranscript;
class GeneBuilder
{
public:
    GeneBuilder() {}
    std::vector< GeneFromGTF > makeGene(GTFIterator& iter);

private:
    std::vector< GeneFromGTF > makeGeneFromMultiVersionGTFRecords(std::vector< GTFRecord >& gtfRecords);
    bool                       geneInDiffChrs(std::vector< GTFRecord >& gtfRecords);
    GTFMapForGeneVersion       gatherByGeneVersion(std::vector< GTFRecord >& gtfRecords);
    GTFMapForTranscript        gatherByTransciptID(std::vector< GTFRecord >& gtfRecords);
    GeneFromGTF                makeGeneWithTranscriptsFromGTFRecords(std::vector< GTFRecord >& gtfRecords);
    std::vector< GeneFromGTF > makeGenesWithTranscriptsFromGTFRecords(std::vector< GTFRecord >& gtfRecords);
    GeneFromGTF                makeGeneFromGTFRecords(std::vector< GTFRecord >& gtfRecords);
    void                       validateGTFRecord(GTFRecord& record, GeneFromGTF& gene);
    void addTranscriptToGeneFromGTFRecords(GeneFromGTF& gene, std::vector< GTFRecord >& gtfRecords);
};
