/*
 * File: bamRecord.h
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#pragma once

#include <string>
using std::string;

#include <spdlog/spdlog.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

typedef bam1_t*    BamRecord;
typedef bam_hdr_t* BamHeader;

#define FQCFAIL BAM_FQCFAIL
#define FDUP BAM_FDUP

// Seperator of qname, e.g. 39_19_58583_28608@TTGGCGGGT@V300062757_T67L3C006R0070857078
// New version: V300062757_T67L3C006R0070857078|||CB:Z:39_19_58583_28608|||UR:Z:TTGGCGGGT|||UY:Z:E*5)<E0+(
static const string QNAME_SEP  = "|||";
static const int    PREFIX_LEN = 5;

int getQName(BamRecord b, std::string& qname);

int getQual(BamRecord b);

std::vector< std::pair< int, int > > getCigar(BamRecord b);

int getRefStart(BamRecord b);

bool getNegativeStrand(BamRecord b);

int updateStrTags(BamRecord& b, std::string tag, std::string data);

BamRecord createBamRecord();

int destroyBamRecord(BamRecord b);

bool compareBamRecord(BamRecord b1, BamRecord b2);

bool getTag(BamRecord b, const char tag[2], std::string& value);

bool getTagInt(BamRecord b, const char tag[2], int& value);

void setQcFail(BamRecord b);

void setDuplication(BamRecord b);

bool getQcFail(BamRecord b);

bool getDuplication(BamRecord b);

int getMarker(BamRecord b, std::string& marker);