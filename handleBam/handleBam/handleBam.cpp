/*
 * File: handleBam.cpp
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#include <exception>
#include <fstream>
#include <iomanip>
#include <set>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <unordered_set>

#include <htslib/thread_pool.h>

#include "annotationException.h"
#include "bamCat.h"
#include "bamRecord.h"
#include "density/kde.hpp"
#include "geneBuilder.h"
#include "gtfReader.h"
#include "gzIO.h"
#include "handleBam.h"
#include "intervalTree.h"
#include "samReader.h"
#include "samWriter.h"
#include "threadpool.h"
#include "timer.h"
#include "utils.h"

static std::set< std::string > EXCLUDE_REFS{ "chrGL", "chrNC", "chrhs", "random", "chrU", "chrEK", "chrAQ" };
static const char              GE_TAG[]    = "GE";
static const char              CB_TAG[]    = "CB";
static const char              UR_TAG[]    = "UR";
static const char              UY_TAG[]    = "UY";
static const char              HI_TAG[]    = "HI";
const unsigned short           MAX_THREADS = 24;

// CB/UR are in extra fields
//#define OLD_QNAME

static const int BASES_NUM = 4;

static const char BASES_DECODE[] = "ACGT";

constexpr int EXCESS_CONTIGS_NUM = 10000;

std::tuple< int, int, int, int > HandleBam::processChromosome(std::string           ctg,
                                                              TagReadsWithGeneExon* tagReadsWithGeneExon)
{
    Timer t;

    int       total = 0, filtered = 0, annotated = 0, unique = 0;
    BamRecord bamRecord = createBamRecord();

    std::unordered_set< std::string > read_set;
    // read_set.reserve(10*1000*1000);
    std::string                            marker;
    std::unordered_map< std::string, int > barcode_gene_exp;
    std::string                            ge_value, qname;
    int      hi_index;  // Query hit index, indicating the alignment record is the i-th one stored in SAM
    fs::path tmp_exp_file = tmp_exp_path / (ctg + ".txt");

    fs::path                     tmp_bam_file = tmp_bam_path / (ctg + ".bam");
    SamWriter                    samWriter(tmp_bam_file.string());
    std::unique_ptr< SamReader > sr = SamReader::FromFile(input_bam_filenames[0]);
    samWriter.init(sr->getHeader());

    for (auto& input_bam : input_bam_filenames)
    {
        std::unique_ptr< SamReader > sr      = SamReader::FromFile(input_bam);
        auto                         contigs = sr->getContigs();
        // Find index of contig
        int chr_id = -1;
        for (size_t i = 0; i < contigs.size(); ++i)
            if (contigs[i].first == ctg)
            {
                chr_id = i;
                break;
            }
        if (chr_id == -1 || !sr->QueryByContig(chr_id))
            continue;

        // Iterate over each record.
        while (true)
        {
            if (!sr->next(bamRecord))
                break;

            // Deduplication of STAR before process
            if (getTagInt(bamRecord, HI_TAG, hi_index) && hi_index != 1)
                continue;

            ++total;

            // Deduplication in each chromosome.
            char*       qname     = bam_get_qname(bamRecord);
            std::string barcode   = "";
            std::string umi       = "";
            std::string umi_score = "";
#ifdef OLD_QNAME
            getTag(bamRecord, CB_TAG, barcode);
            getTag(bamRecord, UR_TAG, umi);
            if (barcode.empty() || umi.empty())
                continue;
#else
            {
                std::string_view qname_view(qname);
                while (!qname_view.empty())
                {
                    size_t      pos = qname_view.find(QNAME_SEP);
                    std::string tmp;
                    if (pos != std::string::npos)
                    {
                        tmp        = qname_view.substr(0, pos);
                        qname_view = qname_view.substr(pos + QNAME_SEP.size());
                    }
                    else
                    {
                        tmp        = qname_view;
                        qname_view = "";
                    }
                    std::string prefix = tmp.substr(0, PREFIX_LEN);
                    if (prefix == "CB:Z:")
                    {
                        barcode = tmp.substr(PREFIX_LEN);
                    }
                    else if (prefix == "UR:Z:")
                    {
                        umi = tmp.substr(PREFIX_LEN);
                    }
                    else if (prefix == "UY:Z:")
                    {
                        umi_score = tmp.substr(PREFIX_LEN);
                    }
                }
            }
            int ret = 0;
            ret     = bam_aux_append(bamRecord, CB_TAG, 'Z', barcode.size() + 1, ( uint8_t* )barcode.c_str());
            if (ret != 0)
                spdlog::warn("bam_aux_append CB failed:{}", strerror(errno));
            if (!umi.empty())
            {
                ret = bam_aux_append(bamRecord, UR_TAG, 'Z', umi.size() + 1, ( uint8_t* )umi.c_str());
                if (ret != 0)
                    spdlog::warn("bam_aux_append UR failed:{}", strerror(errno));
            }
            if (!umi_score.empty())
            {
                ret = bam_aux_append(bamRecord, UY_TAG, 'Z', umi_score.size() + 1, ( uint8_t* )umi_score.c_str());
                if (ret != 0)
                    spdlog::warn("bam_aux_append UY failed:{}", strerror(errno));
            }

            qname = bam_get_qname(bamRecord);
            std::string_view qname_view(qname);
            size_t           pos = qname_view.find(QNAME_SEP);
            if (pos != std::string::npos)
            {
                memset(qname + pos, 0, qname_view.size() - pos);
                bamRecord->core.l_extranul += qname_view.size() - pos;
            }
#endif

            // Filter mapping quality.
            int score = getQual(bamRecord);
            if (score < mapping_quality_threshold)
            {
                // Save the reads that qc failed
                if (bam_config.save_lq)
                {
                    setQcFail(bamRecord);
                    samWriter.write(bamRecord);
                }
                continue;
            }
            ++filtered;

            tagReadsWithGeneExon->setAnnotation(bamRecord, ctg);
            if (!getTag(bamRecord, GE_TAG, ge_value))
            {
                samWriter.write(bamRecord);
                continue;
            }
            ++annotated;

            getMarker(bamRecord, marker);
            marker += barcode;
            if (read_set.count(marker) != 0)
            {
                // Save the reads that duplicate
                if (bam_config.save_dup)
                {
                    setDuplication(bamRecord);
                    samWriter.write(bamRecord);
                }
                continue;
            }
            read_set.insert(marker);
            ++unique;

            // Calculate barcode gene expression
            barcode_gene_exp[barcode + "\t" + ge_value]++;

            // Write disk of output bam data
            samWriter.write(bamRecord);
        }
    }
    samWriter.close();

    if (!barcode_gene_exp.empty())
    {
        std::ofstream exp_handle(tmp_exp_file, std::ofstream::out);
        if (!exp_handle.is_open())
        {
            std::string error = "Error opening file: " + tmp_exp_file.string();
            spdlog::error(error);
            throw std::runtime_error(error);
        }
        for (auto& p : barcode_gene_exp)
            exp_handle << p.first << "\t" << p.second << std::endl;
        exp_handle.close();
    }

    bam_destroy1(bamRecord);

    spdlog::info("chr:{} total:{} filtered:{} annotated:{} unique:{} time(s):{:.2f}", ctg, total, filtered, annotated,
                 unique, t.toc(1000));
    return make_tuple(total, filtered, annotated, unique);
}

std::tuple< int, int, int, int > HandleBam::processChromosomeWhole(std::string           ctg,
                                                                   TagReadsWithGeneExon* tagReadsWithGeneExon)
{
    Timer t;

    int       total = 0, filtered = 0, annotated = 0, unique = 0;
    BamRecord bamRecord = createBamRecord();

    std::unordered_set< std::string > read_set;
    // read_set.reserve(10*1000*1000);
    std::string                            marker;
    std::unordered_map< std::string, int > barcode_gene_exp;
    std::string                            ge_value, qname;
    int      hi_index;  // Query hit index, indicating the alignment record is the i-th one stored in SAM
    fs::path tmp_exp_file = tmp_exp_path / (ctg + ".txt");

    fs::path                     tmp_bam_file = tmp_bam_path / (ctg + ".bam");
    SamWriter                    samWriter(tmp_bam_file.string());
    std::unique_ptr< SamReader > sr = SamReader::FromFile(input_bam_filenames[0]);
    samWriter.init(sr->getHeader());

    for (auto& input_bam : input_bam_filenames)
    {
        std::unique_ptr< SamReader > sr = SamReader::FromFile(input_bam);

        // Iterate over each record.
        while (true)
        {
            if (!sr->QueryAll(bamRecord))
                break;

            // Deduplication of STAR before process
            if (getTagInt(bamRecord, HI_TAG, hi_index) && hi_index != 1)
                continue;

            ++total;

            // Deduplication in each chromosome.
            char*       qname     = bam_get_qname(bamRecord);
            std::string barcode   = "";
            std::string umi       = "";
            std::string umi_score = "";
#ifdef OLD_QNAME
            getTag(bamRecord, CB_TAG, barcode);
            getTag(bamRecord, UR_TAG, umi);
            if (barcode.empty() || umi.empty())
                continue;
#else
            {
                std::string_view qname_view(qname);
                while (!qname_view.empty())
                {
                    size_t      pos = qname_view.find(QNAME_SEP);
                    std::string tmp;
                    if (pos != std::string::npos)
                    {
                        tmp        = qname_view.substr(0, pos);
                        qname_view = qname_view.substr(pos + QNAME_SEP.size());
                    }
                    else
                    {
                        tmp        = qname_view;
                        qname_view = "";
                    }
                    std::string prefix = tmp.substr(0, PREFIX_LEN);
                    if (prefix == "CB:Z:")
                    {
                        barcode = tmp.substr(PREFIX_LEN);
                    }
                    else if (prefix == "UR:Z:")
                    {
                        umi = tmp.substr(PREFIX_LEN);
                    }
                    else if (prefix == "UY:Z:")
                    {
                        umi_score = tmp.substr(PREFIX_LEN);
                    }
                }
            }
            int ret = 0;
            ret     = bam_aux_append(bamRecord, CB_TAG, 'Z', barcode.size() + 1, ( uint8_t* )barcode.c_str());
            if (ret != 0)
                spdlog::warn("bam_aux_append CB failed:{}", strerror(errno));
            if (!umi.empty())
            {
                ret = bam_aux_append(bamRecord, UR_TAG, 'Z', umi.size() + 1, ( uint8_t* )umi.c_str());
                if (ret != 0)
                    spdlog::warn("bam_aux_append UR failed:{}", strerror(errno));
            }
            if (!umi_score.empty())
            {
                ret = bam_aux_append(bamRecord, UY_TAG, 'Z', umi_score.size() + 1, ( uint8_t* )umi_score.c_str());
                if (ret != 0)
                    spdlog::warn("bam_aux_append UY failed:{}", strerror(errno));
            }

            qname = bam_get_qname(bamRecord);
            std::string_view qname_view(qname);
            size_t           pos = qname_view.find(QNAME_SEP);
            if (pos != std::string::npos)
            {
                memset(qname + pos, 0, qname_view.size() - pos);
                bamRecord->core.l_extranul += qname_view.size() - pos;
            }
#endif

            // Filter mapping quality.
            int score = getQual(bamRecord);
            if (score < mapping_quality_threshold)
            {
                // Save the reads that qc failed
                if (bam_config.save_lq)
                {
                    setQcFail(bamRecord);
                    samWriter.write(bamRecord);
                }
                continue;
            }
            ++filtered;

            string ctg = sr->refName(bamRecord);
            tagReadsWithGeneExon->setAnnotation(bamRecord, ctg);
            if (!getTag(bamRecord, GE_TAG, ge_value))
            {
                samWriter.write(bamRecord);
                continue;
            }
            ++annotated;

            getMarker(bamRecord, marker);
            marker += barcode;

            if (read_set.count(marker) != 0)
            {
                // Save the reads that duplicate
                if (bam_config.save_dup)
                {
                    setDuplication(bamRecord);
                    samWriter.write(bamRecord);
                }
                continue;
            }
            read_set.insert(marker);
            ++unique;

            // Calculate barcode gene expression
            barcode_gene_exp[barcode + "\t" + ge_value]++;

            // Write disk of output bam data
            samWriter.write(bamRecord);
        }
    }
    samWriter.close();

    if (!barcode_gene_exp.empty())
    {
        std::ofstream exp_handle(tmp_exp_file, std::ofstream::out);
        if (!exp_handle.is_open())
        {
            std::string error = "Error opening file: " + tmp_exp_file.string();
            spdlog::error(error);
            throw std::runtime_error(error);
        }
        for (auto& p : barcode_gene_exp)
            exp_handle << p.first << "\t" << p.second << std::endl;
        exp_handle.close();
    }

    bam_destroy1(bamRecord);

    spdlog::info("chr:{} total:{} filtered:{} annotated:{} unique:{} time(s):{:.2f}", ctg, total, filtered, annotated,
                 unique, t.toc(1000));
    return make_tuple(total, filtered, annotated, unique);
}

std::tuple< int, int, int, int > HandleBam::processChromosomeUmi(std::string           ctg,
                                                                 TagReadsWithGeneExon* tagReadsWithGeneExon)
{
    Timer t;
    int   total = 0, filtered = 0, annotated = 0, unique = 0;

    BamRecord                                                                 bamRecord = createBamRecord();
    std::unordered_map< std::string, std::unordered_map< std::string, int > > umi_mismatch;

    std::string                                              marker;
    std::unordered_map< std::string, std::pair< int, int > > barcode_gene_exp;
    std::string                                              ge_value, qname;
    int hi_index;  // Query hit index, indicating the alignment record is the i-th one stored in SAM

    fs::path tmp_exp_file = tmp_exp_path / ctg;
    tmp_exp_file += ".txt";

    int index = 0;
    for (auto& input_bam : input_bam_filenames)
    {
        ++index;
        std::unique_ptr< SamReader > sr      = SamReader::FromFile(input_bam);
        auto                         contigs = sr->getContigs();
        // Find index of contig
        int chr_id = -1;
        for (size_t i = 0; i < contigs.size(); ++i)
            if (contigs[i].first == ctg)
            {
                chr_id = i;
                break;
            }
        if (chr_id == -1 || !sr->QueryByContig(chr_id))
            continue;

        fs::path  tmp_bam_file = tmp_bam_path / (ctg + "_" + to_string(index) + ".bam");
        SamWriter samWriter(tmp_bam_file.string());
        samWriter.init(sr->getHeader());

        // First read bam: fitler, set annotations, stat {barcode_gene: {umi: cnt}},
        // deduplication the first time, and write disk
        while (true)
        {
            if (!sr->next(bamRecord))
                break;
            // Deduplication of STAR before process
            if (getTagInt(bamRecord, HI_TAG, hi_index) && hi_index != 1)
                continue;

            ++total;

            char* qname = bam_get_qname(bamRecord);
            // Barcode length may be changed, e.g. 41_16_61982_24418 / 0_0_0_0
            // if (qname.size() <= barcode_len || qname[barcode_len] != QNAME_SEP)
            std::string barcode   = "";
            std::string umi       = "";
            std::string umi_score = "";
#ifdef OLD_QNAME
            getTag(bamRecord, CB_TAG, barcode);
            getTag(bamRecord, UR_TAG, umi);
            if (barcode.empty() || umi.empty())
                continue;
#else
            {
                std::string_view qname_view(qname);
                while (!qname_view.empty())
                {
                    size_t      pos = qname_view.find(QNAME_SEP);
                    std::string tmp;
                    if (pos != std::string::npos)
                    {
                        tmp        = qname_view.substr(0, pos);
                        qname_view = qname_view.substr(pos + QNAME_SEP.size());
                    }
                    else
                    {
                        tmp        = qname_view;
                        qname_view = "";
                    }
                    std::string prefix = tmp.substr(0, PREFIX_LEN);
                    if (prefix == "CB:Z:")
                    {
                        barcode = tmp.substr(PREFIX_LEN);
                    }
                    else if (prefix == "UR:Z:")
                    {
                        umi = tmp.substr(PREFIX_LEN);
                    }
                    else if (prefix == "UY:Z:")
                    {
                        umi_score = tmp.substr(PREFIX_LEN);
                    }
                }
            }

            // Cut flags from qname and paste them the extra fileds
            int ret = 0;
            ret     = bam_aux_append(bamRecord, CB_TAG, 'Z', barcode.size() + 1, ( uint8_t* )barcode.c_str());
            if (ret != 0)
                spdlog::warn("bam_aux_append CB failed:{}", strerror(errno));
            ret = bam_aux_append(bamRecord, UR_TAG, 'Z', umi.size() + 1, ( uint8_t* )umi.c_str());
            if (ret != 0)
                spdlog::warn("bam_aux_append UR failed:{}", strerror(errno));
            if (!umi_score.empty())
            {
                ret = bam_aux_append(bamRecord, UY_TAG, 'Z', umi_score.size() + 1, ( uint8_t* )umi_score.c_str());
                if (ret != 0)
                    spdlog::warn("bam_aux_append UY failed:{}", strerror(errno));
            }
            qname = bam_get_qname(bamRecord);
            std::string_view qname_view(qname);
            size_t           pos = qname_view.find(QNAME_SEP);
            if (pos != std::string::npos)
            {
                memset(qname + pos, 0, qname_view.size() - pos);
                bamRecord->core.l_extranul += qname_view.size() - pos;
            }
#endif
            // Filter mapping quality.
            int score = getQual(bamRecord);
            if (score < mapping_quality_threshold)
            {
                // For total reads in sequencing saturation
                ge_value        = "NOGENE";
                std::string key = barcode + "|" + ge_value;
                if (umi_mismatch.count(key) == 0)
                    umi_mismatch[key] = {};
                umi_mismatch[key][umi]++;

                // Save the reads that qc failed
                if (bam_config.save_lq)
                {
                    setQcFail(bamRecord);
                    samWriter.write(bamRecord);
                }
                continue;
            }
            ++filtered;

            // Discard reads that umi has 'N'
            if (umi.find('N') != std::string::npos)
                continue;

            // Set annotations, need the gene name for the next step
            tagReadsWithGeneExon->setAnnotation(bamRecord, ctg);

            // Calculate barcode gene expression
            if (getTag(bamRecord, GE_TAG, ge_value))
            {
                ++annotated;
                std::string key = barcode + "|" + ge_value;
                if (umi_mismatch.count(key) == 0)
                    umi_mismatch[key] = {};
                umi_mismatch[key][umi]++;
                if (umi_mismatch[key][umi] > 1)
                {
                    // Save the reads that duplicate
                    if (bam_config.save_dup)
                        setDuplication(bamRecord);
                    else
                        continue;
                }
                else
                    ++unique;
            }
            else
            {
                // For total reads in sequencing saturation
                ge_value        = "NOGENE";
                std::string key = barcode + "|" + ge_value;
                if (umi_mismatch.count(key) == 0)
                    umi_mismatch[key] = {};
                umi_mismatch[key][umi]++;

                // Do not save the reads with no gene name
                // continue;
            }

            // Write disk of output bam data
            samWriter.write(bamRecord);
        }
        samWriter.close();
    }

    if (total == 0)
    {
        // Free resources
        bam_destroy1(bamRecord);

        spdlog::info("chr:{} total:{} filtered:{} annotated:{} unique:{} time(s):{:.2f}", ctg, total, filtered,
                     annotated, unique, t.toc(1000));
        return make_tuple(total, filtered, annotated, unique);
    }

    // Calcluate which pattern of barcode_gene_umi should be duplicated
    std::unordered_map< std::string, std::unordered_map< std::string, std::string > > umi_correct;
    deDupUmi(umi_mismatch, umi_correct);

    fs::path                     inter_bam_file = tmp_bam_path / (ctg + ".bam");
    SamWriter                    samWriter2(inter_bam_file.string());
    std::unique_ptr< SamReader > sr2 = SamReader::FromFile(input_bam_filenames[0]);
    samWriter2.init(sr2->getHeader());
    for (size_t index = 1; index <= input_bam_filenames.size(); ++index)
    {
        fs::path tmp_bam_file = tmp_bam_path / (ctg + "_" + to_string(index) + ".bam");
        if (!fs::exists(tmp_bam_file))
            continue;  // maybe not exists

        std::unique_ptr< SamReader > sr2     = SamReader::FromFile(tmp_bam_file.string());
        auto                         contigs = sr2->getContigs();
        // Find index of contig
        int chr_id = -1;
        for (size_t i = 0; i < contigs.size(); ++i)
            if (contigs[i].first == ctg)
            {
                chr_id = i;
                break;
            }

        if (chr_id == -1 || !sr2->QueryByContig(chr_id))
            continue;

        // Second read bam: deduplication the second time, stat gene expression and write disk
        while (true)
        {
            if (!sr2->next(bamRecord))
                break;

            if ((bam_config.save_lq && getQcFail(bamRecord)) || (bam_config.save_dup && getDuplication(bamRecord)))
            {
                samWriter2.write(bamRecord);
                continue;
            }

            // Barcode length may be changed, e.g. 41_16_61982_24418 / 0_0_0_0
            // if (qname.size() <= barcode_len || qname[barcode_len] != QNAME_SEP)
            std::string barcode = "";
            std::string umi     = "";
            getTag(bamRecord, CB_TAG, barcode);
            getTag(bamRecord, UR_TAG, umi);

            // Calculate barcode gene expression
            if (getTag(bamRecord, GE_TAG, ge_value))
            {
                std::string key = barcode + "|" + ge_value;

                if (umi_mismatch[key][umi] == 0)
                {
                    --unique;
                    // Save the reads that duplicate
                    if (bam_config.save_dup)
                    {
                        int         ret     = 0;
                        std::string correct = umi_correct[key][umi];
                        ret = bam_aux_append(bamRecord, "UB", 'Z', correct.size() + 1, ( uint8_t* )correct.c_str());
                        if (ret != 0)
                            spdlog::warn("bam_aux_append UB failed:{}", strerror(errno));
                        setDuplication(bamRecord);
                    }
                    else
                        continue;
                }
                else
                {
                    // Calculate barcode gene expression
                    barcode_gene_exp[barcode + "\t" + ge_value].first++;
                    barcode_gene_exp[barcode + "\t" + ge_value].second += umi_mismatch[key][umi];
                }
            }

            // Write disk of output bam data
            samWriter2.write(bamRecord);
        }
    }
    samWriter2.close();

    if (!barcode_gene_exp.empty())
    {
        std::ofstream exp_handle(tmp_exp_file, std::ofstream::out);
        if (!exp_handle.is_open())
        {
            std::string error = "Error opening file: " + tmp_exp_file.string();
            spdlog::error(error);
            throw std::runtime_error(error);
        }
        if (scrna)
        {
            for (auto& p : barcode_gene_exp)
                exp_handle << p.first << "\t" << p.second.first << "\t" << p.second.second << std::endl;
        }
        else
        {
            for (auto& p : barcode_gene_exp)
                exp_handle << p.first << "\t" << p.second.first << std::endl;
        }

        exp_handle.close();
    }

    // std::unordered_map< std::string, int > tmp_map;
    // barcode_gene_exp.swap(tmp_map);
    // std::unordered_set< std::string > tmp_set;
    // read_set.swap(tmp_set);
    bam_destroy1(bamRecord);

    if (saturation)
    {
        saturation->addData(umi_mismatch);
    }

    spdlog::info("chr:{} total:{} filtered:{} annotated:{} unique:{} time(s):{:.2f}", ctg, total, filtered, annotated,
                 unique, t.toc(1000));
    return make_tuple(total, filtered, annotated, unique);
}

std::tuple< int, int, int, int > HandleBam::processChromosomeUmiWhole(std::string           ctg,
                                                                      TagReadsWithGeneExon* tagReadsWithGeneExon)
{
    Timer t;
    int   total = 0, filtered = 0, annotated = 0, unique = 0;

    BamRecord                                                                 bamRecord = createBamRecord();
    std::unordered_map< std::string, std::unordered_map< std::string, int > > umi_mismatch;

    std::string                                              marker;
    std::unordered_map< std::string, std::pair< int, int > > barcode_gene_exp;
    std::string                                              ge_value, qname;
    int hi_index;  // Query hit index, indicating the alignment record is the i-th one stored in SAM

    fs::path tmp_exp_file = tmp_exp_path / ctg;
    tmp_exp_file += ".txt";

    int index = 0;
    for (auto& input_bam : input_bam_filenames)
    {
        ++index;
        std::unique_ptr< SamReader > sr = SamReader::FromFile(input_bam);

        fs::path  tmp_bam_file = tmp_bam_path / (ctg + "_" + to_string(index) + ".bam");
        SamWriter samWriter(tmp_bam_file.string());
        samWriter.init(sr->getHeader());

        // First read bam: fitler, set annotations, stat {barcode_gene: {umi: cnt}},
        // deduplication the first time, and write disk
        while (true)
        {
            if (!sr->QueryAll(bamRecord))
                break;
            // Deduplication of STAR before process
            if (getTagInt(bamRecord, HI_TAG, hi_index) && hi_index != 1)
                continue;

            ++total;

            char* qname = bam_get_qname(bamRecord);
            // Barcode length may be changed, e.g. 41_16_61982_24418 / 0_0_0_0
            // if (qname.size() <= barcode_len || qname[barcode_len] != QNAME_SEP)
            std::string barcode   = "";
            std::string umi       = "";
            std::string umi_score = "";
#ifdef OLD_QNAME
            getTag(bamRecord, CB_TAG, barcode);
            getTag(bamRecord, UR_TAG, umi);
            if (barcode.empty() || umi.empty())
                continue;
#else
            {
                std::string_view qname_view(qname);
                while (!qname_view.empty())
                {
                    size_t      pos = qname_view.find(QNAME_SEP);
                    std::string tmp;
                    if (pos != std::string::npos)
                    {
                        tmp        = qname_view.substr(0, pos);
                        qname_view = qname_view.substr(pos + QNAME_SEP.size());
                    }
                    else
                    {
                        tmp        = qname_view;
                        qname_view = "";
                    }
                    std::string prefix = tmp.substr(0, PREFIX_LEN);
                    if (prefix == "CB:Z:")
                    {
                        barcode = tmp.substr(PREFIX_LEN);
                    }
                    else if (prefix == "UR:Z:")
                    {
                        umi = tmp.substr(PREFIX_LEN);
                    }
                    else if (prefix == "UY:Z:")
                    {
                        umi_score = tmp.substr(PREFIX_LEN);
                    }
                }
            }

            // Cut flags from qname and paste them the extra fileds
            int ret = 0;
            ret     = bam_aux_append(bamRecord, CB_TAG, 'Z', barcode.size() + 1, ( uint8_t* )barcode.c_str());
            if (ret != 0)
                spdlog::warn("bam_aux_append CB failed:{}", strerror(errno));
            ret = bam_aux_append(bamRecord, UR_TAG, 'Z', umi.size() + 1, ( uint8_t* )umi.c_str());
            if (ret != 0)
                spdlog::warn("bam_aux_append UR failed:{}", strerror(errno));
            if (!umi_score.empty())
            {
                ret = bam_aux_append(bamRecord, UY_TAG, 'Z', umi_score.size() + 1, ( uint8_t* )umi_score.c_str());
                if (ret != 0)
                    spdlog::warn("bam_aux_append UY failed:{}", strerror(errno));
            }
            qname = bam_get_qname(bamRecord);
            std::string_view qname_view(qname);
            size_t           pos = qname_view.find(QNAME_SEP);
            if (pos != std::string::npos)
            {
                memset(qname + pos, 0, qname_view.size() - pos);
                bamRecord->core.l_extranul += qname_view.size() - pos;
            }
#endif
            // Filter mapping quality.
            int score = getQual(bamRecord);
            if (score < mapping_quality_threshold)
            {
                // For total reads in sequencing saturation
                ge_value        = "NOGENE";
                std::string key = barcode + "|" + ge_value;
                if (umi_mismatch.count(key) == 0)
                    umi_mismatch[key] = {};
                umi_mismatch[key][umi]++;

                // Save the reads that qc failed
                if (bam_config.save_lq)
                {
                    setQcFail(bamRecord);
                    samWriter.write(bamRecord);
                }
                continue;
            }
            ++filtered;

            // Discard reads that umi has 'N'
            if (umi.find('N') != std::string::npos)
                continue;

            string ctg = sr->refName(bamRecord);
            // Set annotations, need the gene name for the next step
            tagReadsWithGeneExon->setAnnotation(bamRecord, ctg);

            // Calculate barcode gene expression
            if (getTag(bamRecord, GE_TAG, ge_value))
            {
                ++annotated;
                std::string key = barcode + "|" + ge_value;
                if (umi_mismatch.count(key) == 0)
                    umi_mismatch[key] = {};
                umi_mismatch[key][umi]++;
                if (umi_mismatch[key][umi] > 1)
                {
                    // Save the reads that duplicate
                    if (bam_config.save_dup)
                        setDuplication(bamRecord);
                    else
                        continue;
                }
                else
                    ++unique;
            }
            else
            {
                // For total reads in sequencing saturation
                ge_value        = "NOGENE";
                std::string key = barcode + "|" + ge_value;
                if (umi_mismatch.count(key) == 0)
                    umi_mismatch[key] = {};
                umi_mismatch[key][umi]++;

                // Do not save the reads with no gene name
                // continue;
            }

            // Write disk of output bam data
            samWriter.write(bamRecord);
        }
        samWriter.close();
    }

    if (total == 0)
    {
        // Free resources
        bam_destroy1(bamRecord);

        spdlog::info("total:{} filtered:{} annotated:{} unique:{} time(s):{:.2f}", ctg, total, filtered, annotated,
                     unique, t.toc(1000));
        return make_tuple(total, filtered, annotated, unique);
    }

    // Calcluate which pattern of barcode_gene_umi should be duplicated
    std::unordered_map< std::string, std::unordered_map< std::string, std::string > > umi_correct;
    deDupUmi(umi_mismatch, umi_correct);

    fs::path                     inter_bam_file = tmp_bam_path / (ctg + ".bam");
    SamWriter                    samWriter2(inter_bam_file.string());
    std::unique_ptr< SamReader > sr2 = SamReader::FromFile(input_bam_filenames[0]);
    samWriter2.init(sr2->getHeader());
    for (size_t index = 1; index <= input_bam_filenames.size(); ++index)
    {
        fs::path tmp_bam_file = tmp_bam_path / (ctg + "_" + to_string(index) + ".bam");
        if (!fs::exists(tmp_bam_file))
            continue;  // maybe not exists

        std::unique_ptr< SamReader > sr2 = SamReader::FromFile(tmp_bam_file.string());

        // Second read bam: deduplication the second time, stat gene expression and write disk
        while (true)
        {
            if (!sr2->QueryAll(bamRecord))
                break;

            if ((bam_config.save_lq && getQcFail(bamRecord)) || (bam_config.save_dup && getDuplication(bamRecord)))
            {
                samWriter2.write(bamRecord);
                continue;
            }

            // Barcode length may be changed, e.g. 41_16_61982_24418 / 0_0_0_0
            // if (qname.size() <= barcode_len || qname[barcode_len] != QNAME_SEP)
            std::string barcode = "";
            std::string umi     = "";
            getTag(bamRecord, CB_TAG, barcode);
            getTag(bamRecord, UR_TAG, umi);

            // Calculate barcode gene expression
            if (getTag(bamRecord, GE_TAG, ge_value))
            {
                std::string key = barcode + "|" + ge_value;

                if (umi_mismatch[key][umi] == 0)
                {
                    --unique;
                    // Save the reads that duplicate
                    if (bam_config.save_dup)
                    {
                        int         ret     = 0;
                        std::string correct = umi_correct[key][umi];
                        ret = bam_aux_append(bamRecord, "UB", 'Z', correct.size() + 1, ( uint8_t* )correct.c_str());
                        if (ret != 0)
                            spdlog::warn("bam_aux_append UB failed:{}", strerror(errno));
                        setDuplication(bamRecord);
                    }
                    else
                        continue;
                }
                else
                {
                    // Calculate barcode gene expression
                    barcode_gene_exp[barcode + "\t" + ge_value].first++;
                    barcode_gene_exp[barcode + "\t" + ge_value].second += umi_mismatch[key][umi];
                }
            }

            // Write disk of output bam data
            samWriter2.write(bamRecord);
        }
    }
    samWriter2.close();

    if (!barcode_gene_exp.empty())
    {
        std::ofstream exp_handle(tmp_exp_file, std::ofstream::out);
        if (!exp_handle.is_open())
        {
            std::string error = "Error opening file: " + tmp_exp_file.string();
            spdlog::error(error);
            throw std::runtime_error(error);
        }
        if (scrna)
        {
            for (auto& p : barcode_gene_exp)
                exp_handle << p.first << "\t" << p.second.first << "\t" << p.second.second << std::endl;
        }
        else
        {
            for (auto& p : barcode_gene_exp)
                exp_handle << p.first << "\t" << p.second.first << std::endl;
        }
        exp_handle.close();
    }

    // std::unordered_map< std::string, int > tmp_map;
    // barcode_gene_exp.swap(tmp_map);
    // std::unordered_set< std::string > tmp_set;
    // read_set.swap(tmp_set);
    bam_destroy1(bamRecord);

    if (saturation)
    {
        saturation->addData(umi_mismatch);
    }

    spdlog::info("chr:{} total:{} filtered:{} annotated:{} unique:{} time(s):{:.2f}", ctg, total, filtered, annotated,
                 unique, t.toc(1000));
    return make_tuple(total, filtered, annotated, unique);
}

// Calculate mismatch of two umis, and mismatch positions/types
int HandleBam::umiDistance(const std::string& s1, const std::string& s2, vector< int >& types, vector< int >& positions)
{
    types.clear();
    positions.clear();
    int distance = 0;
    for (size_t i = 0; i < s1.size(); ++i)
    {
        if (s1[i] != s2[i])
        {
            ++distance;
            types.push_back(BASES_ENCODE[( unsigned char )s1[i]] * BASES_NUM + BASES_ENCODE[( unsigned char )s2[i]]);
            positions.push_back(i);
        }
    }
    return distance;
}

inline bool compareBySecond(const pair< std::string, int >& p1, const pair< std::string, int >& p2)
{
    return p1.second > p2.second;
}

// Mark the duplicate umi through set cnt to 0 in {barcode_gene : {umi: cnt}}
int HandleBam::deDupUmi(std::unordered_map< std::string, std::unordered_map< std::string, int > >&         umi_mismatch,
                        std::unordered_map< std::string, std::unordered_map< std::string, std::string > >& umi_correct)
{
    if (umi_mismatch.empty())
        return 0;

    vector< pair< std::string, int > > array;
    int                                umi_total_nums = 0, umi_dedup_nums = 0;
    int                                umi_mis_types[64]    = { 0 };
    int                                umi_mis_postions[64] = { 0 };
    std::vector< int >                 types;
    std::vector< int >                 positions;

    for (auto& p : umi_mismatch)
    {
        umi_total_nums += p.second.size();
        umi_dedup_nums += p.second.size();
        if (p.second.size() < umi_config.min_num)
            continue;
        // Check gene name of 'NOGENE'
        std::string_view tmp = p.first;
        if (tmp.substr(tmp.find('|') + 1) == "NOGENE")
            continue;

        // spdlog::info("deDupUmi key:{}", p.first);
        array.clear();
        // Transform data from map to vector<pair> for sorting by value
        for (const auto& umi : p.second)
        {
            array.push_back({ umi.first, umi.second });
            // spdlog::info("{} {}", umi.first, umi.second);
        }
        // Sort the vector by cnt
        sort(array.begin(), array.end(), compareBySecond);
        // spdlog::info("sort by cnt:");
        // for (const auto& umi : array)
        // {
        //     spdlog::info("{} {}", umi.first, umi.second);
        // }

        // Pairwise comparison of all umis
        for (size_t i = 1; i < array.size(); ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                // Not this umi because it was already marked
                if (p.second[array[j].first] == 0)
                    continue;
                // Calculate distance of two umis
                if (umiDistance(array[i].first, array[j].first, types, positions) <= umi_config.mismatch)
                {
                    // spdlog::info("dedup {} {} {}", array[i].first, array[j].first, umiDistance(array[i].first,
                    // array[j].first));
                    p.second[array[j].first] += p.second[array[i].first];
                    p.second[array[i].first] = 0;
                    --umi_dedup_nums;

                    // Mark the correct umi
                    if (umi_correct.count(p.first) == 0)
                        umi_correct[p.first] = {};
                    umi_correct[p.first][array[i].first] = array[j].first;

                    // Calculate mismatch types and mismatch positions
                    for (auto& t : types)
                        umi_mis_types[t]++;
                    for (auto& p : positions)
                        umi_mis_postions[p]++;
                    break;
                }
            }
        }
        // spdlog::info("dedup by distance:");
        // for (const auto& umi : p.second)
        // {
        //    spdlog::info("{} {}", umi.first, umi.second);
        // }

        // break;
    }

    // Accumulate metrics of umi
    metrics_mutex.lock();

    for (size_t i = 0; i < umi_len; ++i)
        umi_metrics.umi_mis_postions[i] += umi_mis_postions[i];

    for (int i = 0; i < BASES_NUM * BASES_NUM; ++i)
        umi_metrics.umi_mis_types[i] += umi_mis_types[i];

    umi_metrics.uniq_barcode_gene_nums += umi_mismatch.size();
    umi_metrics.umi_cnt_raw += umi_total_nums;
    umi_metrics.umi_cnt_dedup += umi_dedup_nums;

    metrics_mutex.unlock();

    return 0;
}

int HandleBam::doWork()
{
    if (createPath() != 0)
        return -2;

    // Load annotations.
    TagReadsWithGeneExon tagReadsWithGeneExon(annotation_filename);
    if (tagReadsWithGeneExon.makeOverlapDetectorV2() != 0)
    {
        spdlog::error("Failed makeOverlapDetector!");
        return -1;
    }
    tagReadsWithGeneExon.setAnnoVersion(bam_config.anno_ver);

    // Open input bam file
    // There maybe more than one bam file, check they have same header
    samReader                                                     = SamReader::FromFile(input_bam_filenames[0]);
    std::vector< std::pair< std::string, unsigned int > > contigs = samReader->getContigs();
    {
        if (input_bam_filenames.size() > 1)
        {
            for (size_t i = 1; i < input_bam_filenames.size(); ++i)
            {
                std::unique_ptr< SamReader > reader      = SamReader::FromFile(input_bam_filenames[i]);
                auto                         tmp_contigs = reader->getContigs();
                if (tmp_contigs != contigs)
                {
                    spdlog::error("Different header of bam files: {} {}", input_bam_filenames[0],
                                  input_bam_filenames[i]);
                    return -4;
                }
            }
        }
    }
    spdlog::debug("Bam contigs num:{}", contigs.size());

    // Check if umi exists
#ifndef OLD_QNAME
    if (checkUmi() != 0)
        return -3;
#endif

    spdlog::info("Using threads num:{}", cpu_cores);

    Timer total_timer;

    size_t total = 0, filtered = 0, annotated = 0, unique = 0;

    // Using theadpool to accelerate process
    std::threadpool executor{ static_cast< unsigned short >(cpu_cores) };

    std::vector< std::future< std::tuple< int, int, int, int > > > results;
    // Iterate over each contig.
    if (cpu_cores == 1 || contigs.size() > EXCESS_CONTIGS_NUM)
    {
        string ctg = "whole";
        contigs.clear();
        contigs.push_back({ ctg, 0 });
        if (umi_config.on)
            results.emplace_back(
                executor.commit(std::bind(&HandleBam::processChromosomeUmiWhole, this, ctg, &tagReadsWithGeneExon)));
        else
            results.emplace_back(
                executor.commit(std::bind(&HandleBam::processChromosomeWhole, this, ctg, &tagReadsWithGeneExon)));
    }
    else
    {
        for (auto& [ctg, _] : contigs)
        {
            spdlog::debug("start query contig:{}", ctg);
            if (EXCLUDE_REFS.count(ctg) != 0)
                continue;
            // if (!samReader->QueryByContig(i))
            //     continue;
            // When to use special version for umi? According to user input parameter and check if umi exists in bam
            // file
            if (umi_config.on)
                results.emplace_back(
                    executor.commit(std::bind(&HandleBam::processChromosomeUmi, this, ctg, &tagReadsWithGeneExon)));
            else
                results.emplace_back(
                    executor.commit(std::bind(&HandleBam::processChromosome, this, ctg, &tagReadsWithGeneExon)));
        }
    }

    for (auto&& result : results)
    {
        int n1, n2, n3, n4;
        std::tie(n1, n2, n3, n4) = result.get();
        total += n1;
        filtered += n2;
        annotated += n3;
        unique += n4;
    }
    spdlog::debug("Process max memory(KB):{}", physical_memory_used_by_process());
    spdlog::info("Process time(s):{:.2f}", total_timer.toc(1000));

    // Merge bam files and gene expression files
    // std::string                cat_exp_cmd = "cat > " + exp_file;

    std::ofstream              ofs_exp(exp_file, std::ofstream::out);
    std::vector< std::string > bam_files;
    for (auto& [ctg, _] : contigs)
    {
        fs::path tmp_bam_file = tmp_bam_path / ctg;
        tmp_bam_file += ".bam";
        if (fs::exists(tmp_bam_file))
        {
            bam_files.push_back(tmp_bam_file.string());
        }
        fs::path tmp_exp_file = tmp_exp_path / ctg;
        tmp_exp_file += ".txt";
        if (fs::exists(tmp_exp_file))
        {
            std::ifstream ifs(tmp_exp_file, std::ifstream::in);
            ofs_exp << ifs.rdbuf();
            ifs.close();
            // cat_exp_cmd += " " + tmp_exp_file.string();
        }
    }
    ofs_exp.close();
    // cat_exp_cmd += " 2>&1";
    int                        cmd_rtn;
    std::vector< std::string > cmd_result;

    cmd_rtn = bam_cat(bam_files, nullptr, output_bam_filename.c_str(), nullptr, 0);
    if (cmd_rtn == 0)
        spdlog::info("Merge bam file success");
    else
        spdlog::info("Merge bam file fail, rtn:{}", cmd_rtn);
    spdlog::info("Merge bam and gene expression file time(s):{:.2f}", total_timer.toc(1000));

    // spdlog::debug("cat exp cmd: {}", cat_exp_cmd);
    // cmd_rtn = exec_shell(cat_exp_cmd.c_str(), cmd_result);
    // for (const auto& line : cmd_result)
    //     if (!line.empty())
    //         spdlog::error(line);
    // if (cmd_rtn == 0)
    //     spdlog::info("Merge expression file success");
    // else
    //     spdlog::info("Merge expression file fail, rtn:{}", cmd_rtn);

    // spdlog::info("Merge exp file time(s):{:.2f}", total_timer.toc(1000));

    float filter_rate        = total != 0 ? float(total - filtered) * 100 / total : 0;
    float fail_annotate_rate = filtered != 0 ? float(filtered - annotated) * 100 / filtered : 0;
    float dup_rate           = annotated != 0 ? float(annotated - unique) * 100 / annotated : 0;
    spdlog::info("Total reads:{} Pass filter reads:{} Annotated reads:{} Unique reads:{}", total, filtered, annotated,
                 unique);
    spdlog::info("Failed filter rate:{:.2f}%", filter_rate);
    spdlog::info("Failed annotate rate:{:.2f}%", fail_annotate_rate);
    spdlog::info("Duplication rate:{:.2f}%", dup_rate);

    // Dump metrics information to disk.
    std::string   anno_metrics = tagReadsWithGeneExon.dumpMetrics();
    std::ofstream ofs(metrics_filename, std::ofstream::out);
    if (ofs.is_open())
    {
        ofs.precision(2);
        ofs.setf(std::ios::fixed);
        ofs << "## "
               "FILTER & DEDUPLICATION "
               "METRICS\nTOTAL_READS\tPASS_FILTER\tANNOTATED_READS\tUNIQUE_READS\tFAIL_FILTER_RATE\tFAIL_ANNOTATE_"
               "RATE\tDUPLICATION_RATE\n";
        ofs << total << "\t" << filtered << "\t" << annotated << "\t" << unique << "\t" << filter_rate << "\t"
            << fail_annotate_rate << "\t" << dup_rate << std::endl;
        ofs << anno_metrics;
        if (umi_config.on)
        {
            ofs << "## UMI CORRECTIONS METRICS\n"
                   "BARCODE_GENE_NUM\tUMI_CNT_RAW\tUMI_CNT_DEDUP\tRAW_PCT DEDUP_PCT\n";
            ofs << umi_metrics.uniq_barcode_gene_nums << "\t" << umi_metrics.umi_cnt_raw << "\t"
                << umi_metrics.umi_cnt_dedup << "\t"
                << (umi_metrics.uniq_barcode_gene_nums != 0
                        ? umi_metrics.umi_cnt_raw * 100.0 / umi_metrics.uniq_barcode_gene_nums
                        : 0)
                << "\t"
                << (umi_metrics.uniq_barcode_gene_nums != 0
                        ? umi_metrics.umi_cnt_dedup * 100.0 / umi_metrics.uniq_barcode_gene_nums
                        : 0)
                << std::endl;

            size_t total_cnt = 0;
            for (size_t i = 0; i < umi_len; ++i)
                total_cnt += umi_metrics.umi_mis_postions[i];
            if (total_cnt == 0)
                total_cnt = 1;
            ofs << "## UMI MISMATCH POSITIONS METRICS\n"
                   "POSITION\tCNT\tPCT\n";
            for (size_t i = 0; i < umi_len; ++i)
                ofs << i + 1 << "\t" << umi_metrics.umi_mis_postions[i] << "\t"
                    << umi_metrics.umi_mis_postions[i] * 100.0 / total_cnt << std::endl;
            ofs << "## UMI MISMATCH TYPES METRICS\n"
                   "TYPE\tCNT\tPCT\n";
            for (int i = 0; i < BASES_NUM; ++i)
            {
                for (int j = 0; j < BASES_NUM; ++j)
                {
                    int k = i * BASES_NUM + j;
                    ofs << BASES_DECODE[i] << "_" << BASES_DECODE[j] << "\t" << umi_metrics.umi_mis_types[k] << "\t"
                        << umi_metrics.umi_mis_types[k] * 100.0 / total_cnt << std::endl;
                }
            }
        }
        ofs.close();

        spdlog::info("Success dump metrics file:{}", metrics_filename);
    }
    else
    {
        spdlog::error("Error opening file:{}", metrics_filename);
    }

    // Calculate sequencing saturation
    if (saturation)
    {
        saturation->calculateSaturation(sat_file);
    }

    if (scrna)
    {
        // Transform barcode gene expression file format to
        // matrix markert file format
        if (!transform_txt2mtx())
        {
            spdlog::warn("Failed transform txt to matrix market file!");
        }
        else
        {
            spdlog::info("Success transform txt to matrix market file!");
        }
    }

    spdlog::debug("Finish doWork");
    return 0;
}

int HandleBam::createPath()
{
    fs::path p = output_bam_filename;
    spdlog::debug("output_bam_filename:{}", p.string());
    p = p.parent_path();
    spdlog::debug("parent path:{}", p.string());
    tmp_bam_path = p;
    tmp_bam_path += "/_bam";
    spdlog::debug("tmp_bam_path:{}", tmp_bam_path.string());
    tmp_exp_path = p;
    tmp_exp_path += "/_exp";

    std::error_code ec;

    if (!fs::exists(tmp_bam_path) && !fs::create_directories(tmp_bam_path, ec))
    {
        spdlog::error("Failed create directories:{} error:{}", tmp_bam_path.string(), ec.message());
        return -1;
    }
    if (!fs::exists(tmp_exp_path) && !fs::create_directories(tmp_exp_path, ec))
    {
        spdlog::error("Failed create directories:{} error:{}", tmp_exp_path.string(), ec.message());
        return -1;
    }
    return 0;
}

void HandleBam::setBamConfig(bool save_lq, bool save_dup, int anno_mode)
{
    bam_config.save_lq  = save_lq;
    bam_config.save_dup = save_dup;
    bam_config.anno_ver = AnnoVersion(anno_mode);
}

void HandleBam::setUmiConfig(bool on, int min_num, int mismatch)
{
    umi_config.on       = on;
    umi_config.min_num  = min_num;
    umi_config.mismatch = mismatch;
}

// Check if umi exists? If true, then set length of barcode and umi
int HandleBam::checkUmi()
{
    barcode_len = 0;
    umi_len     = 0;

    BamRecord bamRecord = createBamRecord();
    samReader->QueryOne(bamRecord);
    std::string qname;
    getQName(bamRecord, qname);
    while (!qname.empty())
    {
        size_t      pos = qname.find(QNAME_SEP);
        std::string tmp;
        if (pos != std::string::npos)
        {
            tmp   = qname.substr(0, pos);
            qname = qname.substr(pos + QNAME_SEP.size());
        }
        else
        {
            tmp   = qname;
            qname = "";
        }

        std::string prefix = tmp.substr(0, PREFIX_LEN);
        // std::cout<<tmp<<std::endl;
        if (prefix == "CB:Z:")
        {
            barcode_len = tmp.size() - PREFIX_LEN;
        }
        else if (prefix == "UR:Z:")
        {
            umi_len = tmp.size() - PREFIX_LEN;
        }
    }

    if (barcode_len == 0)
    {
        spdlog::info("No barcode found.");
        return -1;
    }

    if (umi_len != 0)
    {
        spdlog::info("Barcode length:{} umi length:{}", barcode_len, umi_len);
    }
    else
    {
        if (umi_config.on)
        {
            spdlog::error("No umi found");
            return -1;
        }
    }

    bam_destroy1(bamRecord);
    return 0;
}

void HandleBam::setExtraConfig(std::string _sat_file, bool _filter_matrix, int _cores, bool _scrna)
{
    sat_file      = _sat_file;
    filter_matrix = _filter_matrix;
    cpu_cores     = _cores;
    scrna         = _scrna;

    saturation = nullptr;
    if (!sat_file.empty())
    {
        if (scrna)
        {
            saturation = new SequenceBarcode();
        }
        else
        {
            saturation = new CoordinateBarcode();
        }
    }
}

// Transform barcode gene expression file format to
// matrix markert file format
bool HandleBam::transform_txt2mtx()
{
    if (!fs::exists(exp_file))
        return false;
    const char sep = ' ';
    fs::path   output_path(exp_file);
    output_path = output_path.parent_path();

    // Load the raw data from disk
    ifstream                                                           ifs(exp_file);
    string                                                             line;
    unordered_map< string, unordered_map< string, pair< int, int > > > raw_data;
    while (std::getline(ifs, line))
    {
        vector< string > vec_s    = split_str(line, '\t');
        string           barcode  = vec_s[0];
        string           gene     = vec_s[1];
        int              umi_cnt  = stoi(vec_s[2]);
        int              read_cnt = vec_s.size() > 3 ? stoi(vec_s[3]) : 0;

        if (raw_data.count(barcode) == 0)
            raw_data[barcode] = {};
        auto& p = raw_data[barcode][gene];
        p.first += umi_cnt;
        p.second += read_cnt;
    }
    ifs.close();

    // Do the filter of low barcode count
    unordered_set< string > discard_barcodes;
    if (filter_matrix)
    {
        unordered_map< string, int > times;
        vector< double >             cnts;
        for (auto& p : raw_data)
        {
            int cnt = 0;
            for (auto& pp : p.second)
                cnt += pp.second.first;
            times[p.first] = cnt;
            cnts.push_back(cnt);
        }

        KDE    kde;
        auto   paras             = kde.run(cnts, "bead");
        double min_barcode_frags = paras.first;
        // double call_threshold = paras.second;
        spdlog::info("Barcode threshold of filtering matrix: {}", min_barcode_frags);

        for (auto& p : times)
        {
            if (p.second < min_barcode_frags)
                discard_barcodes.insert(p.first);
        }

        if (!discard_barcodes.empty())
        {
            fs::path ofs_exp(exp_file);
            if (ofs_exp.extension() == ".gz")
            {
                cmpFile out_exp;
                out_exp = cmpOpen(ofs_exp.c_str());
                for (auto& p : raw_data)
                {
                    if (discard_barcodes.count(p.first) != 0)
                        continue;

                    for (auto& pp : p.second)
                    {
                        string s = p.first + "\t" + pp.first + "\t" + to_string(pp.second.first);
                        if (pp.second.second != 0)
                            s += "\t" + to_string(pp.second.second);
                        s += "\n";
                        cmpFunc(out_exp, s.c_str());
                    }
                }
                cmpClose(out_exp);
            }
            else
            {
                ofstream out_exp(ofs_exp, std::ofstream::out);
                for (auto& p : raw_data)
                {
                    if (discard_barcodes.count(p.first) != 0)
                        continue;

                    for (auto& pp : p.second)
                    {
                        string s = p.first + "\t" + pp.first + "\t" + to_string(pp.second.first);
                        if (pp.second.second != 0)
                            s += "\t" + to_string(pp.second.second);
                        s += "\n";
                        out_exp << s;
                    }
                }
                out_exp.close();
            }
            spdlog::info("Dump filtered gene expression file");
        }
    }

    unordered_map< string, int >           map_barcodes, map_genes;
    unordered_map< string, int >::iterator it;
    vector< string >                       vec_barcodes, vec_genes, vec_matrix;
    // Travel all barcode-gene-cnt, calculate the order of barcodes, genes,
    // then the final matrix data
    for (auto& p : raw_data)
    {
        string barcode = p.first;
        if (discard_barcodes.count(barcode) != 0)
            continue;

        int i1, i2;
        it = map_barcodes.find(barcode);
        if (it == map_barcodes.end())
        {
            vec_barcodes.push_back(barcode);
            i1                    = vec_barcodes.size();
            map_barcodes[barcode] = i1;
        }
        else
        {
            i1 = it->second;
        }
        for (auto& pp : p.second)
        {
            string gene = pp.first;
            int    cnt  = pp.second.first;

            it = map_genes.find(gene);
            if (it == map_genes.end())
            {
                vec_genes.push_back(gene);
                i2              = vec_genes.size();
                map_genes[gene] = i2;
            }
            else
            {
                i2 = it->second;
            }
            vec_matrix.push_back(to_string(i1) + sep + to_string(i2) + sep + to_string(cnt));
        }
    }

    fs::path ofs_barcode(output_path / "barcodes.tsv.gz"), ofs_genes(output_path / "genes.tsv.gz"),
        ofs_matrix(output_path / "matrix.mtx.gz");
    cmpFile out_barcode;
    out_barcode = cmpOpen(ofs_barcode.c_str());
    for (auto& l : vec_barcodes)
    {
        string s = l + "\n";
        cmpFunc(out_barcode, s.c_str());
    }
    cmpClose(out_barcode);

    cmpFile out_genes;
    out_genes = cmpOpen(ofs_genes.c_str());
    for (auto& l : vec_genes)
    {
        string s = l + "\n";
        cmpFunc(out_genes, s.c_str());
    }
    cmpClose(out_genes);

    cmpFile out_matrix;
    out_matrix = cmpOpen(ofs_matrix.c_str());
    stringstream ss;
    ss << "%%MatrixMarket matrix coordinate real general\n%\n";
    ss << map_barcodes.size() << sep << map_genes.size() << sep << vec_matrix.size() << endl;
    cmpFunc(out_matrix, ss.str().c_str());
    for (auto& l : vec_matrix)
    {
        string s = l + "\n";
        cmpFunc(out_matrix, s.c_str());
    }
    cmpClose(out_matrix);

    return true;
}