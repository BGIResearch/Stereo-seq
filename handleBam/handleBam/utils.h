/*
 * File: utils.h
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#pragma once

#include <string>
#include <vector>
using namespace std;

// For executing shell script
int exec_shell(const char* cmd, std::vector< std::string >& resvec);

// Get physical memory used by process in real time
// Only support linux
size_t physical_memory_used_by_process();

vector< string > split_str(const std::string& str, char delim = ' ', bool skip_empty = true);

// Return true: first file is older than second file
bool check_file_older(const std::string& first, const std::string& second);