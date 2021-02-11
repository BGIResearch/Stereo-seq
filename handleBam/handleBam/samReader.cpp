/*
 * File: samReader.cpp
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include <exception>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>
using namespace std;

#include <filesystem>
namespace fs = std::filesystem;

#include <spdlog/spdlog.h>

#include "samReader.h"
#include "types.h"
#include "utils.h"

// Class Implementations.
const int            _DEFAULT_HTS_BLOCK_SIZE = 128 * (1024 * 1024);
const int            QUAL_THRESHOLD          = 10;
static set< string > EXCLUDE_REFS{ "chrGL", "chrNC", "chrhs", "random", "chrU", "chrEK", "chrAQ" };
static char          SEP = '@';

// Common functions.
inline string parseQname(string& qname)
{
    string res;
    size_t pos0 = qname.find(SEP);
    res += qname.substr(0, pos0);
    size_t pos1 = qname.find(SEP, pos0 + 1);
    res += qname.substr(pos0 + 1, pos1 - pos0 - 1);
    return res;
}

// Class SamReader's functions.
SamReader::~SamReader()
{
    if (fp_)
    {
        Close();
    }
}

void check_index(const std::string& reads_path)
{
    // The htslib index data structure for our indexed BAM file. May be NULL if no
    // index was loaded.
    const int threads_num = 4;

    std::string bai_idx_path = reads_path + ".bai";
    std::string csi_idx_path = reads_path + ".csi";
    if (fs::exists(bai_idx_path) && !check_file_older(bai_idx_path, reads_path))
        return;
    if (fs::exists(csi_idx_path) && !check_file_older(csi_idx_path, reads_path))
        return;

    // Try Create bai index first
    if (!fs::exists(bai_idx_path) || check_file_older(bai_idx_path, reads_path))
    {
        int ret = sam_index_build3(reads_path.c_str(), bai_idx_path.c_str(), 0, threads_num);
        if (ret == 0)
            return;
    }

    // Bai-format limit < 512MB, try csi-format
    int ret = sam_index_build3(reads_path.c_str(), csi_idx_path.c_str(), 14, threads_num);
    if (ret != 0)
        throw std::runtime_error("Failed to create index for " + reads_path);
}

std::unique_ptr< SamReader > SamReader::FromFile(const std::string& reads_path)
{
    htsFile* fp = hts_open(reads_path.c_str(), "r");
    if (!fp)
    {
        char msg[128];
        sprintf(msg, "Could not open %s\n", reads_path.c_str());
        // printf(msg);
        std::runtime_error error(msg);
        throw std::exception(error);
    }

    // if (hts_set_opt(fp, HTS_OPT_BLOCK_SIZE, _DEFAULT_HTS_BLOCK_SIZE) != 0)
    // {
    //     char msg[128];
    //     sprintf(msg, "Failed to set HTS_OPT_BLOCK_SIZE\n");
    //     // printf(msg);
    //     std::runtime_error error(msg);
    //     throw std::exception(error);
    // }

    bam_hdr_t* header = sam_hdr_read(fp);
    if (header == nullptr)
    {
        char msg[128];
        sprintf(msg, "Couldn't parse header for %s\n", fp->fn);
        // printf(msg);
        std::runtime_error error(msg);
    }

    check_index(reads_path);

    hts_idx_t* idx = sam_index_load(fp, fp->fn);
    return std::unique_ptr< SamReader >(new SamReader(fp, header, idx));
}

int SamReader::QueryAll()
{
    bam1_t* b = bam_init1();
    // htsFile* out = hts_open("out.bam", "wb");
    // sam_hdr_write(out, header_);

    /*
    while (sam_read1(fp_, header_, b) >= 0)
    {
        if (b->core.qual < QUAL_THRESHOLD) continue;
        int r = sam_write1(out, header_, b);
    }
    */
    unordered_set< string > read_set;
    size_t                  total = 0, unique = 0;
    for (size_t i = 0; i < ref_.size(); ++i)
    {
        // cout<<"query ref:"<<ref_[i].first<<" "<<ref_[i].second<<endl;
        if (EXCLUDE_REFS.count(ref_[i].first) != 0)
            continue;
        hts_itr_t* iter = sam_itr_queryi(idx_, i, 0, ref_[i].second);
        if (iter == nullptr)
        {
            cout << "query unknown reference:" << ref_[i].first << endl;
            continue;
        }
        read_set.clear();
        while (sam_itr_next(fp_, iter, b) > 0)
        {
            // Filter reads which quality value less than quality threshold.
            if (b->core.qual < QUAL_THRESHOLD)
                continue;
            string qname  = bam_get_qname(b);
            string marker = parseQname(qname);
            if (read_set.count(marker) == 0)
            {
                ++unique;
                read_set.insert(marker);
                // sam_write1(out, header_, b);
            }
            ++total;
        }
    }
    cout << "Total Reads: " << total << "\nUnique Reads: " << unique << endl;
    cout << "Duplication rate: " << setiosflags(ios::fixed) << setprecision(2) << float(total - unique) * 100 / total
         << "%" << endl;

    bam_destroy1(b);
    // hts_close(out);
    // out = nullptr;

    return 0;
}

int SamReader::QueryOne(BamRecord b)
{
    for (size_t i = 0; i < ref_.size(); ++i)
    {
        if (EXCLUDE_REFS.count(ref_[i].first) != 0)
            continue;
        hts_itr_t* iter = sam_itr_queryi(idx_, i, 0, ref_[i].second);
        if (iter == nullptr)
        {
            spdlog::warn("query unknown reference:{}", ref_[i].first);
            continue;
        }
        if (sam_itr_next(fp_, iter, b) > 0)
        {
            break;
        }
    }
    return 0;
}

int SamReader::QueryAll(BamRecord b)
{
    if (sam_read1(fp_, header_, b) >= 0)
    {
        return true;
    }
    return false;
}

std::vector< std::pair< std::string, uint32 > > SamReader::getContigs()
{
    return ref_;
}

bool SamReader::QueryByContig(int tid)
{
    if (EXCLUDE_REFS.count(ref_[tid].first) != 0)
        return false;
    iter_ = sam_itr_queryi(idx_, tid, 0, ref_[tid].second);
    if (iter_ == nullptr || iter_->finished)
    {
        spdlog::debug("No reads for ref:{}", ref_[tid].first);
        return false;
    }
    return true;
}

bool SamReader::QueryByContigBE(int tid, const int beg, const int end, hts_itr_t*& iter)
{
    if (EXCLUDE_REFS.count(ref_[tid].first) != 0)
        return false;
    iter = sam_itr_queryi(idx_, tid, beg, end);
    if (iter == nullptr)
    {
        cout << "query unknown reference:" << ref_[tid].first << endl;
        return false;
    }
    // std::cout<<"query iter:"<<iter<<std::endl;
    // bam1_t* b = bam_init1();
    // int ret;
    // while ((ret = sam_itr_next(fp_, iter, b)) >= 0)
    // {
    //     spdlog::info("get");
    // }
    // bam_destroy1(b);

    return true;
}

bool SamReader::next(BamRecord b)
{
    // std::cout<<"iter_:"<<iter_<<" fp_:"<<fp_<< " b:"<< b<<std::endl;
    if (sam_itr_next(fp_, iter_, b) > 0)
    {
        return true;
    }
    return false;
}

std::string SamReader::refName(BamRecord b)
{
    const int tid = b->core.tid;
    if (tid >= 0)
        return header_->target_name[tid];
    else
        return "";
}

bool SamReader::next(BamRecord b, hts_itr_t*& iter)
{
    // std::cout<<"iter:"<<iter<<std::endl;
    int ret;
    if ((ret = sam_itr_next(fp_, iter, b)) > 0)
    {
        return true;
    }
    // std::cout<<"sam_itr_next ret:"<<ret<<std::endl;
    return false;
}

SamReader::SamReader(htsFile* fp, bam_hdr_t* header, hts_idx_t* idx) : fp_(fp), header_(header), idx_(idx)
{
    iter_ = nullptr;
    for (int i = 0; i < header->n_targets; ++i)
        ref_.push_back({ header_->target_name[i], header_->target_len[i] });
}

int SamReader::Close()
{
    hts_close(fp_);
    fp_ = nullptr;

    // free(idx_);
    if (idx_ != nullptr)
    {
        hts_idx_destroy(idx_);
        idx_ = nullptr;
    }
    if (iter_ != nullptr)
    {
        hts_itr_destroy(iter_);
        iter_ = nullptr;
    }

    bam_hdr_destroy(header_);
    header_ = nullptr;

    return 0;
}

BamHeader SamReader::getHeader()
{
    return header_;
}

void SamReader::setThreadPool(htsThreadPool* p)
{
    hts_set_opt(fp_, HTS_OPT_THREAD_POOL, p);
}