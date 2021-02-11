/*
 * File: gtfReader.h
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#pragma once

#include <set>
#include <string>
#include <unordered_map>

#include "gtfRecord.h"

enum GTF_COLUMNS
{
    CHROMOSOME,
    SOURCE,
    FEATURE,
    START,
    END,
    SCORE,
    STRAND,
    FRAME,
    ATTRIBUTE
};

class GTFReader
{
public:
    GTFReader() {}
    int loadGTFFile(std::string gtf_file_, std::unordered_map< std::string, std::vector< GTFRecord > >& gtfMap);
    std::unordered_map< std::string, int >& getContigs()
    {
        return contigs;
    }

private:
    int  parseTabbedLine(std::string& line, GTFRecord& gtfRecord);
    int  parseAttributes(std::string& s);
    void splitString(const std::string& str, char delim, bool skip_empty);

    int parseTabbedLineGFF(std::string& line, GTFRecord& gtfRecord);
    int parseAttributesGFF(std::string& s, std::unordered_map< std::string, std::string >& attrs);

private:
    static const int                       GTF_TAB_NUM = 9;
    std::unordered_map< std::string, int > contigs;

    std::vector< std::string >             splitedLine;
    std::vector< std::string >             attrs;
    std::unordered_map< std::string, int > attrKeys;
    char                                   attrSep;
    std::string                            geneID, geneName, transID, transName;
};
