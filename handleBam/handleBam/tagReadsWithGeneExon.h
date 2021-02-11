/*
 * File: tagReadsWithGeneExon.h
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#pragma once

#include <filesystem>
#include <mutex>
#include <set>
#include <string>
#include <tuple>
using std::string;
#include <atomic>
// namespace fs = std::filesystem;

#include <spdlog/spdlog.h>

// Alternative interval tree implemented base on red-black tree.
#include "ygg.hpp"

#include "bamRecord.h"
#include "bamUtils.h"
#include "geneBuilder.h"
#include "locusFunction.h"

enum AnnoVersion
{
    DROP_SEQ_V2,
    DROP_SEQ_V1,
    TENX,
};

using MyInterval = std::pair< int, int >;
template < class Node > class NodeTraits : public ygg::ITreeNodeTraits< Node >
{
public:
    using key_type = int;
    static int get_lower(const Node& node)
    {
        return node.lower;
    }
    static int get_upper(const Node& node)
    {
        return node.upper;
    }
    static int get_lower(const MyInterval& i)
    {
        return std::get< 0 >(i);
    }
    static int get_upper(const MyInterval& i)
    {
        return std::get< 1 >(i);
    }
};

class Node : public ygg::ITreeNodeBase< Node, NodeTraits< Node > >
{
public:
    // Node(int upper_, int lower_, GeneFremGTF gene_)
    // 	: upper(upper_),
    // 	lower(lower_),
    // 	value(gene_) {}

    int         upper;
    int         lower;
    GeneFromGTF value;
};

using MyTree = ygg::IntervalTree< Node, NodeTraits< Node > >;

class TagReadsWithGeneExon
{
public:
    TagReadsWithGeneExon(string annotation_filename_) : annotation_filename(annotation_filename_)
    {
        currContig = "";

        total_reads               = 0;
        reads_wrong_strand        = 0;
        ambiguous_reads_rejected  = 0;
        read_ambiguous_gene_fixed = 0;
        reads_right_strand        = 0;

        map_reads           = 0;
        exonic_reads        = 0;
        intronic_reads      = 0;
        intergenic_reads    = 0;
        transcriptome_reads = 0;
        multigene_reads     = 0;
        nogene_reads        = 0;

        ALLOW_MULTI_GENE_READS = false;
        USE_STRAND_INFO        = true;

        lastRecord = NULL;
    }
    ~TagReadsWithGeneExon()
    {
        if (lastRecord != NULL)
        {
            bam_destroy1(lastRecord);
            lastRecord = NULL;
        }
    }

    int makeOverlapDetector();
    int makeOverlapDetectorV2();

    int setAnnotation(BamRecord& record);
    int setAnnotation(BamRecord& record, const std::string& contig);

    void setContig(std::string& contig);

    std::string dumpMetrics();

    void setAnnoVersion(AnnoVersion version)
    {
        anno_ver = version;
    }

private:
    std::unordered_map< int, LocusFunction >
                       getLocusFunctionForReadByGene(std::vector< const GeneFromGTF* >& result,
                                                     std::vector< AlignmentBlock >&     alignmentBlock);
    bool               getAlignmentBlockonGeneExon(AlignmentBlock& b, const GeneFromGTF* gene, std::string contig);
    std::set< int >    getAlignmentBlockonGeneExon(std::vector< const GeneFromGTF* >&        result,
                                                   std::unordered_map< int, LocusFunction >& locusMap, AlignmentBlock& b,
                                                   std::string contig);
    std::vector< int > getGenesConsistentWithReadStrand(std::vector< const GeneFromGTF* >& result,
                                                        std::vector< int >& ids, bool recordNegative);
    std::pair< std::string, std::string > getCompoundNameAndStrand(std::vector< const GeneFromGTF* >& result,
                                                                   std::vector< int >&                ids);

    int setAnnotationTS(BamRecord& record, const std::string& contig);
    int setAnnotationTENX(BamRecord& record, const std::string& contig);

private:
    // Only used in query bam by contig, so the contigs are continuous.
    std::string currContig;

    std::unordered_map< std::string, MyTree > mytrees;
    std::vector< Node >                       nodes;

    std::unordered_map< std::string, int > contigs;

    string annotation_filename;

    // Annotation statics for Drop-seq
    std::atomic< int > total_reads;
    std::atomic< int > reads_right_strand;
    std::atomic< int > reads_wrong_strand;
    std::atomic< int > ambiguous_reads_rejected;
    std::atomic< int > read_ambiguous_gene_fixed;

    // Add statics metrics for TENX genomics
    std::atomic< int > map_reads;
    std::atomic< int > exonic_reads;
    std::atomic< int > intronic_reads;
    std::atomic< int > intergenic_reads;
    std::atomic< int > transcriptome_reads;
    std::atomic< int > multigene_reads;
    std::atomic< int > nogene_reads;

    // std::mutex stat_mutex;
    std::mutex query_mutex;

    bool ALLOW_MULTI_GENE_READS;
    bool USE_STRAND_INFO;

    std::string FUNCTION_TAG = "XF";
    std::string TAG          = "GE";
    std::string STRAND_TAG   = "GS";

    BamRecord                                                               lastRecord;
    std::tuple< std::string, std::string, std::string, int, int, int, int > lastCache;

    AnnoVersion anno_ver;
};
