/*
 * File: gzIO.h
 * Created Data: 2020-7-9
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#pragma once

#include <zlib.h>
//#include <gzguts.h>

#define cmpFile gzFile
#define cmpOpen(x) gzopen(x, "wb")
#define cmpClose(x) gzclose(x)
#define cmpFunc(x, y) gzputs(x, y)

#include <string>
using namespace std;

bool readline(gzFile f, string& l, int len = 128);