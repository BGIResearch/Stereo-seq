/*
 * File: gtfReader.cpp
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include <spdlog/spdlog.h>

#include "annotationException.h"
#include "gtfReader.h"

void GTFReader::splitString(const std::string& str, char delim, bool skip_empty)
{
    std::istringstream iss(str);
    int                i = 0;
    for (std::string item; getline(iss, item, delim);)
        if (skip_empty && item.empty())
            continue;
        else
            splitedLine[i++] = item;
}

// Read GTF file line by line, parse data what i want and gather by genename.
int GTFReader::loadGTFFile(std::string gtf_file, std::unordered_map< std::string, std::vector< GTFRecord > >& gtfMap)
{
    std::filesystem::path temp_file(gtf_file);
    if (temp_file.extension() == ".gtf")
        attrSep = ' ';
    else if (temp_file.extension() == ".gff")
        attrSep = '=';
    else
        throw AnnotationException("Invalid gtf/gff file format");

    splitedLine.resize(GTF_TAB_NUM);
    attrKeys.insert({ "gene_id", 0 });
    attrKeys.insert({ "gene_name", 1 });
    attrKeys.insert({ "transcript_id", 2 });
    attrKeys.insert({ "transcript_name", 3 });
    attrs.resize(attrKeys.size());
    // Open gtf file.
    std::ifstream input(gtf_file, std::ifstream::in);
    // Travel the record by line.
    std::string line;
    int         num = 0;
    while (std::getline(input, line))
    {
        GTFRecord gtfRecord;
        // Parse data and construct instance of GTFRecord.
        if (line[0] == '#')
            continue;
        if (attrSep == ' ')
            parseTabbedLine(line, gtfRecord);
        else
            parseTabbedLineGFF(line, gtfRecord);
        if (++num % 100000 == 0)
            spdlog::get("gtf")->info("read {:10d} GTF records. Last read position: {}:{}", num, gtfRecord.getContig(),
                                     gtfRecord.getStart());
        // Gather gtfRecord by genename.
        // Skip pseudo genes
        if (!gtfRecord.getGeneName().empty())
            gtfMap[gtfRecord.getGeneName()].push_back(std::move(gtfRecord));
    }
    return 0;
}

int GTFReader::parseTabbedLine(std::string& line, GTFRecord& gtfRecord)
{
    // Format:CHROMOSOME, SOURCE, FEATURE, START, END, SCORE, STRAND, FRAME, ATTRIBUTE
    splitString(line, '\t', true);
    // if (splitedLine.size() != GTF_TAB_NUM)
    // 	throw AnnotationException("GTF file format error!");
    gtfRecord.setInterval(stoi(splitedLine[GTF_COLUMNS::START]), stoi(splitedLine[GTF_COLUMNS::END]),
                          splitedLine[GTF_COLUMNS::CHROMOSOME], splitedLine[GTF_COLUMNS::STRAND] == "-");

    parseAttributes(splitedLine[GTF_COLUMNS::ATTRIBUTE]);
    gtfRecord.setFeatureType(splitedLine[GTF_COLUMNS::FEATURE]);
    gtfRecord.setGeneID(attrs.at(attrKeys["gene_id"]));
    gtfRecord.setGeneName(attrs.at(attrKeys["gene_name"]));
    if (gtfRecord.getFeatureType() != "gene")
    {
        std::string transcript_id = attrs.at(attrKeys["transcript_id"]);
        std::string transcript_name = attrs.at(attrKeys["transcript_name"]);
        gtfRecord.setTranscriptID(transcript_id);
        if (transcript_name.empty())
            transcript_name = transcript_id;
        gtfRecord.setTranscriptName(transcript_id);
        // TODO(fxzhao): when to use transcriptType and geneVersion?
        // gtfRecord.setTranscriptType();
    }

    return 0;
}

int GTFReader::parseTabbedLineGFF(std::string& line, GTFRecord& gtfRecord)
{
    // Format:CHROMOSOME, SOURCE, FEATURE, START, END, SCORE, STRAND, FRAME, ATTRIBUTE
    splitString(line, '\t', true);
    // if (splitedLine.size() != GTF_TAB_NUM)
    // 	throw AnnotationException("GTF file format error!");
    gtfRecord.setInterval(stoi(splitedLine[GTF_COLUMNS::START]), stoi(splitedLine[GTF_COLUMNS::END]),
                          splitedLine[GTF_COLUMNS::CHROMOSOME], splitedLine[GTF_COLUMNS::STRAND] == "-");

    std::unordered_map< std::string, std::string > attrs;
    parseAttributesGFF(splitedLine[GTF_COLUMNS::ATTRIBUTE], attrs);
    gtfRecord.setFeatureType(splitedLine[GTF_COLUMNS::FEATURE]);
    if (gtfRecord.getFeatureType() == "region")
    {
        geneID = geneName = transID = transName = "";
        return 0;
    }
    else if (gtfRecord.getFeatureType() == "gene")
    {
        if (attrs.count("ID"))
            geneID = attrs["ID"];
        if (attrs.count("Name"))
            geneName = attrs["Name"];
    }
    else if (gtfRecord.getFeatureType() == "mRNA")
    {
        if (attrs.count("ID"))
            transID = attrs["ID"];
        if (attrs.count("Name"))
            transName = attrs["Name"];
    }

    gtfRecord.setGeneID(geneID);
    gtfRecord.setGeneName(geneName);
    if (gtfRecord.getFeatureType() != "gene")
    {
        gtfRecord.setTranscriptID(transID);
        gtfRecord.setTranscriptName(transName);
        // TODO(fxzhao): when to use transcriptType and geneVersion?
        // gtfRecord.setTranscriptType();
    }

    return 0;
}

int GTFReader::parseAttributes(std::string& s)
{
    attrs.clear();
    attrs.resize(attrKeys.size());

    std::istringstream        iss(s);
    int                       i     = 0;
    decltype(attrKeys.size()) count = 0;
    std::string               key, value;
    for (std::string item; getline(iss, item, ' ');)
        if (item.empty())
            continue;
        else
        {
            if ((i++) % 2 == 0)
                key = item;
            else
            {
                if (attrKeys.find(key) == attrKeys.end())
                    continue;
                value = item;
                // Remove the quotes in begin and end of the string,
                // include comma in the end.
                // e.g. gene_name "DDX11L1"; level 2;
                if (value.size() > 3)
                    value = value.substr(1, value.size() - 3);
                else
                    value = value.substr(0, value.size() - 1);
                attrs[attrKeys[key]] = value;
                ++count;
                // When we get the information we want, return early.
                if (count == attrKeys.size())
                    return 0;
            }
        }
    return 0;
}

int GTFReader::parseAttributesGFF(std::string& s, std::unordered_map< std::string, std::string >& attrs)
{
    std::istringstream iss(s);
    std::string        key, value;
    for (std::string item; getline(iss, item, ';');)
        if (item.empty())
            continue;
        else
        {
            auto pos = item.find(attrSep);
            if (pos == std::string::npos)
                continue;
            key        = item.substr(0, pos);
            value      = item.substr(pos + 1);
            attrs[key] = value;
        }
    return 0;
}