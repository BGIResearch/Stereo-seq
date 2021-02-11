/*
 * File: samWriter.h
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#pragma once

#include <fstream>

#include <spdlog/spdlog.h>

#include "samReader.h"

class SamWriter
{
public:
    SamWriter(string filename_) : filename(filename_), out(nullptr), header_(nullptr) {}
    ~SamWriter()
    {
        if (out != nullptr)
        {
            close();
        }
    }

    int close()
    {
        // if (header_ != nullptr)
        // {
        // 	bam_hdr_destroy(header_);
        // 	header_ = nullptr;
        // }
        if (out != nullptr)
        {
            hts_close(out);
            out = nullptr;
        }
        return 0;
    }

    int init(BamHeader header)
    {
        out                      = hts_open(filename.c_str(), "wb");
        header_                  = header;
        [[maybe_unused]] int ret = sam_hdr_write(out, header_);

        spdlog::debug("SamWriter init.");
        return 0;
    }

    int write(BamRecord& record)
    {
        [[maybe_unused]] int ret = sam_write1(out, header_, record);

        return 0;
    }

    void setThreadPool(htsThreadPool* p)
    {
        hts_set_opt(out, HTS_OPT_THREAD_POOL, p);
    }

private:
    string filename;

    // A pointer to the htslib file used to access the SAM/BAM data.
    htsFile* out;

    // A htslib header data structure obtained by parsing the header of this BAM.
    bam_hdr_t* header_;
};
