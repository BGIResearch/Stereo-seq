/*
 * File: handleBam.h
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#pragma once

#include <string>
using std::string;
#include <filesystem>  // C++17 only
#include <mutex>
#include <queue>
#include <tuple>
namespace fs = std::filesystem;

#include <spdlog/spdlog.h>

#include "bamRecord.h"
#include "samReader.h"
#include "saturation.h"
#include "tagReadsWithGeneExon.h"

struct UmiConfig
{
    bool   on;        // true means do umi correction
    size_t min_num;   // the minimum umi numbers for correction
    int    mismatch;  // the maximum mismatch for finding duplication reads
};

struct BamConfig
{
    bool        save_lq;   // true means save low quality reads and set flag with BAM_FQCFAIL
    bool        save_dup;  // true means save duplicate reads and set flag with BAM_FDUP
    AnnoVersion anno_ver;  // choose annotation version
    UmiConfig   umi;       // structure for umi correction
};

struct UmiMetrics
{
    UmiMetrics() : uniq_barcode_gene_nums(0), umi_cnt_raw(0), umi_cnt_dedup(0), umi_mis_types{}, umi_mis_postions{} {}
    size_t uniq_barcode_gene_nums;
    size_t umi_cnt_raw;
    size_t umi_cnt_dedup;
    size_t umi_mis_types[64];
    size_t umi_mis_postions[64];
};

class HandleBam
{
public:
    HandleBam(std::vector< std::string >& input_bam_filenames_, string output_bam_filename_,
              string annotation_filename_, string metrics_filename_, int mapping_quality_threshold_, string exp_file_)
        : input_bam_filenames(input_bam_filenames_), output_bam_filename(output_bam_filename_),
          annotation_filename(annotation_filename_), metrics_filename(metrics_filename_),
          mapping_quality_threshold(mapping_quality_threshold_), exp_file(exp_file_), bFinish(false)
    {
        BASES_ENCODE['A'] = 0;
        BASES_ENCODE['C'] = 1;
        BASES_ENCODE['G'] = 2;
        BASES_ENCODE['T'] = 3;
    }

    ~HandleBam()
    {
        if (fs::exists(tmp_bam_path))
            fs::remove_all(tmp_bam_path);
        if (fs::exists(tmp_exp_path))
            fs::remove_all(tmp_exp_path);

        if (saturation)
        {
            delete saturation;
            saturation = nullptr;
        }
    }

    int doWork();

    std::tuple< int, int, int, int > processChromosome(std::string ctg, TagReadsWithGeneExon* tagReadsWithGeneExon);
    std::tuple< int, int, int, int > processChromosomeUmi(std::string ctg, TagReadsWithGeneExon* tagReadsWithGeneExon);
    std::tuple< int, int, int, int > processChromosomeWhole(std::string           ctg,
                                                            TagReadsWithGeneExon* tagReadsWithGeneExon);
    std::tuple< int, int, int, int > processChromosomeUmiWhole(std::string           ctg,
                                                               TagReadsWithGeneExon* tagReadsWithGeneExon);

    int createPath();

    void setBamConfig(bool save_lq, bool save_dup, int anno_mode);
    void setUmiConfig(bool on, int min_num, int mismatch);
    void setExtraConfig(std::string sat_file, bool filter_matrix, int cores, bool scrna);

private:
    int deDupUmi(std::unordered_map< std::string, std::unordered_map< std::string, int > >&         umi_mismatch,
                 std::unordered_map< std::string, std::unordered_map< std::string, std::string > >& umi_correct);
    int checkUmi();
    int umiDistance(const std::string& s1, const std::string& s2, vector< int >& types, vector< int >& positions);
    // Transform barcode gene expression file format to
    // matrix markert file format
    bool transform_txt2mtx();

private:
    std::vector< std::string > input_bam_filenames;
    string                     output_bam_filename;
    string                     annotation_filename;
    string                     metrics_filename;
    int                        mapping_quality_threshold;
    std::string                exp_file;

    std::vector< BamRecord > records;
    std::queue< int >        producer_queue, consumer_queue;
    std::mutex               producer_mutex, consumer_mutex;
    bool                     bFinish;

    std::unique_ptr< SamReader > samReader;

    fs::path tmp_bam_path, tmp_exp_path;

    BamConfig bam_config;
    UmiConfig umi_config;

    size_t barcode_len;
    size_t umi_len;

    int BASES_ENCODE[128];

    UmiMetrics umi_metrics;
    std::mutex metrics_mutex;

    std::string sat_file;
    Saturation* saturation;

    bool filter_matrix;

    int cpu_cores;

    bool scrna;
};
