/*
 * File: samReader.h
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "bamRecord.h"
#include "timer.h"
#include "types.h"

class SamReader
{
public:
    // Creates a new SamReader reading from the BAM file reads_path
    static std::unique_ptr< SamReader > FromFile(const std::string& reads_path);

    ~SamReader();

    // DIsable assignment/copy operations.
    SamReader(const SamReader& other) = delete;
    SamReader& operator=(const SamReader&) = delete;

    // Gets all of the reads that overlap any bases in range.
    int QueryAll();

    // Query reads treat the bam as a whole contig
    int QueryAll(BamRecord b);

    // Get first bam record in order
    int QueryOne(BamRecord b);

    // Close the underlying resource descriptors.
    int Close();

    // Return <chromosome,index> pairs.
    std::vector< std::pair< std::string, uint32 > > getContigs();

    // Return header for construct other instance.
    BamHeader getHeader();

    // Query records by a given contig.
    bool QueryByContig(int tid);

    bool QueryByContigBE(int tid, const int beg, const int end, hts_itr_t*& iter);

    // Iterate all records from query range.
    // NOTICE: must call next() after call QueryByContig()
    bool next(BamRecord b);

    bool next(BamRecord b, hts_itr_t*& iter);

    void setThreadPool(htsThreadPool* p);

    // Parse the contig name of read
    std::string refName(BamRecord b);

private:
    // Private constructor; use FromFile to safely create a SamReader from a
    // file.
    SamReader(htsFile* fp, bam_hdr_t* header, hts_idx_t* idx);

    // A pointer to the htslib file used to access the SAM/BAM data.
    htsFile* fp_;

    // A htslib header data structure obtained by parsing the header of this BAM.
    bam_hdr_t* header_;

    // The htslib index data structure for our indexed BAM file.
    // May be NULL if no index was loaded.
    hts_idx_t* idx_;

    // Store reference name and length from header.
    std::vector< std::pair< std::string, uint32 > > ref_;

    // A pointer to the query result.
    hts_itr_t* iter_;

    Timer  timer;
    double nextTimes;
};
