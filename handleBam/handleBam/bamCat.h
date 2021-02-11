/*
 * File: bamCat.h
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#pragma once

#include <string>
#include <vector>

#include "htslib/bgzf.h"
#include "htslib/sam.h"

int bam_cat(std::vector< std::string > fn, bam_hdr_t* h, const char* outbam, char* arg_list, int no_pg);