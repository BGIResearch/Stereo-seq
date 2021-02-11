/*
 * File: tagReadsWithGeneExon.cpp
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#include <algorithm>
#include <fstream>
#include <sstream>

#include "tagReadsWithGeneExon.h"

#include "annotationException.h"
#include "bamRecord.h"
#include "gtfReader.h"
#include "timer.h"

// Default use drop-seq V2
//#define ANNO_V1

static const std::string RECORD_SEP = ",";

inline bool inRange(int start, int end, int locus)
{
    return (locus >= start && locus <= end);
}

struct CmpLocus
{
    bool operator()(const pair< int, int >& a, const pair< int, int >& b)
    {
        return a.first < b.first || (a.first == b.first && a.second < b.second);
    }
};

std::unordered_map< int, LocusFunction >
TagReadsWithGeneExon::getLocusFunctionForReadByGene(std::vector< const GeneFromGTF* >& result,
                                                    std::vector< AlignmentBlock >&     alignmentBlock)
{
    // Calculate LocusFunction for each overlapped Gene
    std::unordered_map< int, LocusFunction > locusMap;

    if (anno_ver == AnnoVersion::TENX)
    {
        vector< pair< int, int > > total_cnts(result.size());
        for (unsigned j = 0; j < result.size(); ++j)
        {
            const GeneFromGTF* gene     = result[j];
            pair< int, int >   max_cnts = { 0, 0 };  // <exon numbers, intron numbers>
            // For every block, find the confidently result
            for (auto& b : alignmentBlock)
            {
                pair< int, int > cnts = { 0, 0 };
                for (auto& p : gene->getTranscripts())
                {
                    p.second.assignLocusFunction(b.getReferenceStart(), b.getLength(), cnts);
                }
                max_cnts.first += cnts.first;
                max_cnts.second += cnts.second;
            }  // end Block

            int len = 0;
            for (auto& b : alignmentBlock)
                len += b.getLength();
            LocusFunction locusFunction(LocusFunction::INTERGENIC);
            if (max_cnts.first >= int(len * 0.5))
                locusFunction = LocusFunction::CODING;
            else if (max_cnts.second > 0)
                locusFunction = LocusFunction::INTRONIC;

            locusMap[j]   = locusFunction;
            total_cnts[j] = max_cnts;
        }  // end Gene
        // Pick the most confidently
        if (result.size() > 1)
        {
            std::vector< int > genes;
            for (auto& p : locusMap)
            {
                if (genes.empty())
                    genes.push_back(p.first);
                else
                {
                    if (p.second > locusMap[genes[0]])
                    {
                        // Find better gene
                        genes.clear();
                        genes.push_back(p.first);
                    }
                    else if (p.second == locusMap[genes[0]])
                        genes.push_back(p.first);
                }
            }

            std::unordered_map< int, LocusFunction > tmp_map;
            if (genes.size() == 1)
            {
                // Lucky, only one gene
                tmp_map[genes[0]] = locusMap[genes[0]];
            }
            else
            {
                // Pick one gene using the overlap numbers of exon and intro
                auto iter    = std::max_element(total_cnts.begin(), total_cnts.end(), CmpLocus());
                int  pos     = iter - total_cnts.begin();
                tmp_map[pos] = locusMap[pos];
            }
            tmp_map.swap(locusMap);
        }
    }
    else
    {
        for (unsigned j = 0; j < result.size(); ++j)
        {
            const GeneFromGTF* gene = result[j];
            LocusFunction      locusFunction(LocusFunction::INTERGENIC);
            for (auto& b : alignmentBlock)
            {
                for (auto& p : gene->getTranscripts())
                {
                    p.second.assignLocusFunction(b.getReferenceStart(), locusFunction, b.getLength());
                    // Exit early for reducing unnecessary calculations
                    if (locusFunction == LocusFunction::CODING)
                        break;
                }
                // Exit early for reducing unnecessary calculations
                if (locusFunction == LocusFunction::CODING)
                    break;
            }  // end Block
            locusMap[j] = locusFunction;
        }  // end Gene
    }
    return locusMap;
}

inline bool intersect(int s1, int e1, int s2, int e2)
{
    return (s1 <= s2 && s2 <= e1) || (s1 <= e2 && e2 <= e1) || (s2 <= s1 && e1 <= e2);
}

bool TagReadsWithGeneExon::getAlignmentBlockonGeneExon(AlignmentBlock& b, const GeneFromGTF* gene, std::string contig)
{
    for (auto& p : gene->getTranscripts())
    {
        for (auto& e : p.second.getExons())
        {
            if (contig == gene->getContig()
                && intersect(b.getReferenceStart(), b.getReferenceStart() + b.getLength() - 1, e.start, e.end))
                return true;
        }
    }
    return false;
}

std::set< int > TagReadsWithGeneExon::getAlignmentBlockonGeneExon(std::vector< const GeneFromGTF* >&        result,
                                                                  std::unordered_map< int, LocusFunction >& locusMap,
                                                                  AlignmentBlock& b, std::string contig)
{
    std::set< int > geneSet;

    for (auto& pair : locusMap)
    {
        int geneId = pair.first;
        if (getAlignmentBlockonGeneExon(b, result[geneId], contig))
            geneSet.insert(geneId);
    }
    return geneSet;
}

inline bool readAnnotationMatchStrand(bool annoNegative, bool recordNegative)
{
    return annoNegative == recordNegative;
}

std::vector< int > TagReadsWithGeneExon::getGenesConsistentWithReadStrand(std::vector< const GeneFromGTF* >& result,
                                                                          std::vector< int >& ids, bool recordNegative)
{
    std::vector< int > sameStrand;
    std::vector< int > oppositeStrand;
    // std::get<3>(lastCache) = 0;
    // std::get<4>(lastCache) = 0;
    // std::get<5>(lastCache) = 0;
    // std::get<6>(lastCache) = 0;

    for (auto& id : ids)
    {
        bool annoNegative = result[id]->isNegativeStrand();
        bool strandCheck  = readAnnotationMatchStrand(annoNegative, recordNegative);
        if (strandCheck)
            sameStrand.push_back(id);
        else
            oppositeStrand.push_back(id);
    }

    if (sameStrand.size() == 0 && oppositeStrand.size() > 0)
    {
        reads_wrong_strand++;
        // std::get<3>(lastCache) = 1;
        return {};
    }
    if (sameStrand.size() > 1)
    {
        ambiguous_reads_rejected++;
        // std::get<4>(lastCache) = 1;
        return {};
    }
    // otherwise, the read is unambiguously assigned to a gene on the correct strand - the sameStrandSize must be 1
    // as it's not 0 and not > 1.
    if (oppositeStrand.size() > 0)
    {
        read_ambiguous_gene_fixed++;
        // std::get<5>(lastCache) = 1;
    }

    reads_right_strand++;
    // std::get<6>(lastCache) = 1;
    return sameStrand;
}

std::pair< std::string, std::string >
TagReadsWithGeneExon::getCompoundNameAndStrand(std::vector< const GeneFromGTF* >& result, std::vector< int >& ids)
{
    std::string geneName, geneStrand;
    for (auto& id : ids)
    {
        std::string name = result[id]->getName();
        if (geneName.empty())
            geneName = name;
        else
            geneName += RECORD_SEP + name;

        std::string strand = result[id]->isNegativeStrand() ? "-" : "+";
        if (geneStrand.empty())
            geneStrand = strand;
        else
            geneStrand += RECORD_SEP + strand;
        // spdlog::debug("id:{} name:{} strand:{}", id, name, strand);
    }
    return { geneName, geneStrand };
}

void TagReadsWithGeneExon::setContig(std::string& contig)
{
    currContig = contig;
}

// Speed through caching the same result
// Deprecated
int TagReadsWithGeneExon::setAnnotation(BamRecord& record)
{
    // Note, make sure the input bam is sorted, otherwise it will
    // reduce the speed of program.
    if (compareBamRecord(record, lastRecord))
    {
        // spdlog::debug("compareBamRecord true.");
        std::string locusFunction, name, strand;
        int         i1, i2, i3, i4;
        std::tie(locusFunction, name, strand, i1, i2, i3, i4) = lastCache;
        if (!locusFunction.empty())
            updateStrTags(record, FUNCTION_TAG, locusFunction);
        if (!name.empty() && !strand.empty())
        {
            updateStrTags(record, TAG, name);
            updateStrTags(record, STRAND_TAG, strand);
        }
        ++total_reads;
        if (i1 == 1)
            ++reads_wrong_strand;
        if (i2 == 1)
            ++ambiguous_reads_rejected;
        if (i3 == 1)
            ++read_ambiguous_gene_fixed;
        if (i4 == 1)
            ++reads_right_strand;
        return 0;
    }
    bam_copy1(lastRecord, record);

    // spdlog::debug("setAnnotation");
    std::string contig = currContig;

    // spdlog::debug("setAnnotation contig:{}", contig);
    std::vector< std::pair< int, int > > cigars = getCigar(record);

    // Change begin position from 0-based to 1-based
    int beginPos = getRefStart(record) + 1;
    // spdlog::debug("beginPos:{}", beginPos);
    int queryBegin = beginPos;
    int queryEnd   = queryBegin + getReferenceLength(cigars) - 1;
    // spdlog::debug("setAnnotation query:{} {} - {}", contig, queryBegin, queryEnd);
    // GTree::interval_vector overlapped = trees[contig].findOverlapping(queryBegin, queryEnd);
    MyInterval                        query_range{ queryBegin, queryEnd };
    std::vector< const GeneFromGTF* > result;
    for (const auto& o : mytrees[contig].query(query_range))
    {
        // spdlog::debug("{} {}", o.lower, o.upper);
        result.push_back(&(o.value));
    }
    // spdlog::debug("setAnnotation query results num:{}", overlapped.size());

    // std::vector<GeneFromGTF> result(overlapped.size());
    // for (int i = 0; i < overlapped.size(); ++i)
    // 	result[i] = std::move(overlapped[i].value);

    // Set annotations and record the metrics.
    std::vector< AlignmentBlock >            alignmentBlock = getAlignmentBlocks(cigars, beginPos);
    std::unordered_map< int, LocusFunction > locusMap       = getLocusFunctionForReadByGene(result, alignmentBlock);

    // Get Genes which overlapped the exons
    std::set< int > exonsForRead;
    for (auto& b : alignmentBlock)
    {
        std::set< int > temp = getAlignmentBlockonGeneExon(result, locusMap, b, contig);
        if (anno_ver == AnnoVersion::DROP_SEQ_V2)
            exonsForRead.insert(temp.begin(), temp.end());
        else if (anno_ver == AnnoVersion::DROP_SEQ_V1)
        {
            // if result is not empty and blockGenes isn't empty, intersect the set and set the result as this new
            // set.
            if (exonsForRead.size() > 0 && temp.size() > 0)
            {
                if (!ALLOW_MULTI_GENE_READS)
                {
                    // Intersection of two sets
                    for (auto& g : exonsForRead)
                    {
                        if (temp.count(g) == 0)
                            exonsForRead.erase(g);
                    }
                }
                else
                    exonsForRead.insert(temp.begin(), temp.end());
            }
            else
                // if blockGenes is populated and you're here, then result is empty, so set result to these results
                exonsForRead = temp;
        }
        else
        {
        }
    }
    if (anno_ver == AnnoVersion::DROP_SEQ_V2)
    {
        if (!ALLOW_MULTI_GENE_READS && exonsForRead.size() > 1)
            exonsForRead.clear();
    }

    // Determine the final LocusFunction
    std::vector< int > genes;
    for (auto& geneId : exonsForRead)
    {
        LocusFunction f = locusMap[geneId];
        if (f == LocusFunction::CODING || f == LocusFunction::UTR)
            genes.push_back(geneId);
    }
    std::vector< LocusFunction > allPassingFunction;
    if (USE_STRAND_INFO)
    {
        // constrain gene exons to read strand.
        genes = getGenesConsistentWithReadStrand(result, genes, getNegativeStrand(record));
        if (anno_ver == AnnoVersion::DROP_SEQ_V2)
        {
            // only retain functional map entries that are on the correct strand.
            for (auto& pair : locusMap)
            {
                int  geneId       = pair.first;
                bool annoNegative = result[geneId]->isNegativeStrand();
                bool strandCheck  = readAnnotationMatchStrand(annoNegative, getNegativeStrand(record));
                if (strandCheck)
                    allPassingFunction.push_back(pair.second);
            }
        }
    }
    if (anno_ver == AnnoVersion::DROP_SEQ_V2)
    {
        if (!USE_STRAND_INFO)
        {
            for (auto& pair : locusMap)
            {
                allPassingFunction.push_back(pair.second);
            }
        }
        // if strand tag is used, only add locus function values for passing genes.
        for (auto& g : genes)
            allPassingFunction.push_back(locusMap[g]);
    }
    else if (anno_ver == AnnoVersion::DROP_SEQ_V1)
    {
        allPassingFunction.clear();
        for (auto& pair : locusMap)
        {
            allPassingFunction.push_back(pair.second);
        }
    }
    else
    {
    }

    LocusFunction f = getLocusFunction(allPassingFunction, false);

    if (genes.size() > 1 && !ALLOW_MULTI_GENE_READS)
        spdlog::error("There should only be 1 gene assigned to a read for DGE purposes.");

    if (f != LocusFunction::NONE)
    {
        [[maybe_unused]] int ret = updateStrTags(record, FUNCTION_TAG, LocusString[int(f)]);
        std::get< 0 >(lastCache) = LocusString[int(f)];
        // spdlog::debug("updateStrTags {} {} {}", FUNCTION_TAG, LocusString[int(f)], ret);
    }
    else
        std::get< 0 >(lastCache) = "";

    auto p = getCompoundNameAndStrand(result, genes);
    if (!p.first.empty() && !p.second.empty())
    {
        updateStrTags(record, TAG, p.first);
        updateStrTags(record, STRAND_TAG, p.second);
        std::get< 1 >(lastCache) = p.first;
        std::get< 2 >(lastCache) = p.second;
    }
    else
    {
        // TODO(fxzhao): should i erase the tag?
        std::get< 1 >(lastCache) = "";
        std::get< 2 >(lastCache) = "";
    }

    return 0;
}

int TagReadsWithGeneExon::setAnnotation(BamRecord& record, const std::string& contig)
{
    total_reads++;
    if (anno_ver != AnnoVersion::TENX)
        return setAnnotationTS(record, contig);
    else
        return setAnnotationTENX(record, contig);
}

// Thread safe version of Drop-seq
int TagReadsWithGeneExon::setAnnotationTS(BamRecord& record, const std::string& contig)
{
    // spdlog::debug("setAnnotation");
    // std::string contig = currContig;

    // spdlog::debug("setAnnotation contig:{}", contig);
    std::vector< std::pair< int, int > > cigars = getCigar(record);

    // Change begin position from 0-based to 1-based
    int beginPos = getRefStart(record) + 1;
    // spdlog::debug("beginPos:{}", beginPos);
    int queryBegin = beginPos;
    int queryEnd   = queryBegin + getReferenceLength(cigars) - 1;
    // spdlog::debug("setAnnotation query:{} {} - {}", contig, queryBegin, queryEnd);
    // GTree::interval_vector overlapped = trees[contig].findOverlapping(queryBegin, queryEnd);
    MyInterval                        query_range{ queryBegin, queryEnd };
    std::vector< const GeneFromGTF* > result;

    // query_mutex.lock();
    // const auto& tree = mytrees[contig];
    // if (mytrees.find(contig) != mytrees.end())
    if (mytrees.count(contig) != 0)
    {
        const auto& tree = mytrees.at(contig);
        // query_mutex.unlock();
        const auto& res = tree.query(query_range);
        for (const auto& o : res)
        {
            // spdlog::debug("{} {}", o.lower, o.upper);
            result.push_back(&(o.value));
        }
    }

    if (result.empty())
    {
        nogene_reads++;
        return 0;
    }
    // spdlog::debug("setAnnotation query results num:{}", overlapped.size());

    // std::vector<GeneFromGTF> result(overlapped.size());
    // for (int i = 0; i < overlapped.size(); ++i)
    // 	result[i] = std::move(overlapped[i].value);

    // Set annotations and record the metrics.
    std::vector< AlignmentBlock >            alignmentBlock = getAlignmentBlocks(cigars, beginPos);
    std::unordered_map< int, LocusFunction > locusMap       = getLocusFunctionForReadByGene(result, alignmentBlock);

    // Get Genes which overlapped the exons
    std::set< int > exonsForRead;
    for (auto& b : alignmentBlock)
    {
        std::set< int > temp = getAlignmentBlockonGeneExon(result, locusMap, b, contig);
        if (anno_ver == AnnoVersion::DROP_SEQ_V2)
        {
            exonsForRead.insert(temp.begin(), temp.end());
        }
        else if (anno_ver == AnnoVersion::DROP_SEQ_V1)
        {
            // if result is not empty and blockGenes isn't empty, intersect the set and set the result as this new
            // set.
            if (exonsForRead.size() > 0 && temp.size() > 0)
            {
                if (!ALLOW_MULTI_GENE_READS)
                {
                    // Intersection of two sets
                    for (auto& g : exonsForRead)
                    {
                        if (temp.count(g) == 0)
                            exonsForRead.erase(g);
                    }
                }
                else
                    exonsForRead.insert(temp.begin(), temp.end());
            }
            else
                // if blockGenes is populated and you're here, then result is empty, so set result to these results
                exonsForRead = temp;
        }
    }
    if (anno_ver == AnnoVersion::DROP_SEQ_V2)
    {
        if (!ALLOW_MULTI_GENE_READS && exonsForRead.size() > 1)
            exonsForRead.clear();
    }

    // Determine the final LocusFunction
    std::vector< int > genes;
    for (auto& geneId : exonsForRead)
    {
        LocusFunction f = locusMap[geneId];
        if (f == LocusFunction::CODING || f == LocusFunction::UTR)
            genes.push_back(geneId);
    }
    std::vector< LocusFunction > allPassingFunction;
    if (USE_STRAND_INFO)
    {
        // constrain gene exons to read strand.
        genes = getGenesConsistentWithReadStrand(result, genes, getNegativeStrand(record));
        if (anno_ver == AnnoVersion::DROP_SEQ_V2)
        {
            // only retain functional map entries that are on the correct strand.
            for (auto& pair : locusMap)
            {
                int  geneId       = pair.first;
                bool annoNegative = result[geneId]->isNegativeStrand();
                bool strandCheck  = readAnnotationMatchStrand(annoNegative, getNegativeStrand(record));
                if (strandCheck)
                    allPassingFunction.push_back(pair.second);
            }
        }
    }
    if (anno_ver == AnnoVersion::DROP_SEQ_V2)
    {
        if (!USE_STRAND_INFO)
        {
            for (auto& pair : locusMap)
            {
                allPassingFunction.push_back(pair.second);
            }
        }
        // if strand tag is used, only add locus function values for passing genes.
        for (auto& g : genes)
            allPassingFunction.push_back(locusMap[g]);
    }
    else if (anno_ver == AnnoVersion::DROP_SEQ_V1)
    {
        allPassingFunction.clear();
        for (auto& pair : locusMap)
        {
            allPassingFunction.push_back(pair.second);
        }
    }
    else
    {
    }

    LocusFunction f = getLocusFunction(allPassingFunction, false);

    if (genes.size() > 1 && !ALLOW_MULTI_GENE_READS)
        spdlog::error("There should only be 1 gene assigned to a read for DGE purposes.");

    if (f != LocusFunction::NONE)
    {
        [[maybe_unused]] int ret = updateStrTags(record, FUNCTION_TAG, LocusString[int(f)]);
        // std::get<0>(lastCache) = LocusString[int(f)];
        // spdlog::debug("updateStrTags {} {} {}", FUNCTION_TAG, LocusString[int(f)], ret);
    }
    // else
    // std::get<0>(lastCache) = "";

    std::pair< std::string, std::string > p = getCompoundNameAndStrand(result, genes);
    if (!p.first.empty() && !p.second.empty())
    {
        updateStrTags(record, TAG, p.first);
        updateStrTags(record, STRAND_TAG, p.second);
        // std::get<1>(lastCache) = p.first;
        // std::get<2>(lastCache) = p.second;
    }
    else
    {
        // TODO(fxzhao): should i erase the tag?
        // std::get<1>(lastCache) = "";
        // std::get<2>(lastCache) = "";
    }

    return 0;
}

// TENX version
int TagReadsWithGeneExon::setAnnotationTENX(BamRecord& record, const std::string& contig)
{
    // if (record->core.flag == 0)
    //     return 0;
    bool confidently = getQual(record) >= 255;
    if (confidently)
        map_reads++;

    std::vector< std::pair< int, int > > cigars = getCigar(record);

    // Change begin position from 0-based to 1-based
    int                               beginPos   = getRefStart(record) + 1;
    int                               queryBegin = beginPos;
    int                               queryEnd   = queryBegin + getReferenceLength(cigars) - 1;
    MyInterval                        query_range{ queryBegin, queryEnd };
    std::vector< const GeneFromGTF* > result;

    if (mytrees.count(contig) != 0)
    {
        const auto& tree = mytrees.at(contig);
        const auto& res  = tree.query(query_range);
        for (const auto& o : res)
        {
            result.push_back(&(o.value));
        }
    }

    if (result.empty())
    {
        nogene_reads++;
        return 0;
    }
    // Ignore multi-genes
    // if (result.size() > 1)
    // {
    //     multigene_reads++;
    //     return 0;
    // }

    // Set annotations and record the metrics.
    std::vector< AlignmentBlock >            alignmentBlock = getAlignmentBlocks(cigars, beginPos);
    std::unordered_map< int, LocusFunction > locusMap       = getLocusFunctionForReadByGene(result, alignmentBlock);

    // Pick one gene from multi-genes
    std::vector< int > genes;
    for (auto& p : locusMap)
    {
        genes.push_back(p.first);
    }
    // if (genes.size() > 1)
    // {
    //     // multigene_reads++;
    //     // return 0;
    //     std::vector< int > temp;
    //     for (auto& g : genes)
    //     {
    //         // Filter by same strand
    //         if (readAnnotationMatchStrand(result[g]->isNegativeStrand(), getNegativeStrand(record)))
    //             temp.push_back(g);
    //     }
    //     temp.swap(genes);
    // }
    if (genes.empty())
    {
        intergenic_reads++;
        return 0;
    }
    // if (genes.size() > 1 && !ALLOW_MULTI_GENE_READS)
    //     genes.resize(1);

    bool annoNegative = result[genes[0]]->isNegativeStrand();
    bool strandCheck  = readAnnotationMatchStrand(annoNegative, getNegativeStrand(record));
    if (strandCheck)
        reads_right_strand++;
    else
        reads_wrong_strand++;

    LocusFunction        f   = locusMap[genes[0]];
    [[maybe_unused]] int ret = updateStrTags(record, FUNCTION_TAG, LocusString[int(f)]);

    if (confidently)
    {
        if (f == LocusFunction::CODING)
        {
            exonic_reads++;
            if (strandCheck)
                transcriptome_reads++;
        }
        else if (f == LocusFunction::INTERGENIC)
            intergenic_reads++;
        else if (f == LocusFunction::INTRONIC)
            intronic_reads++;
    }

    // Only dump gene name when locus is Exon or Intro
    if (f == LocusFunction::INTERGENIC)
        genes.clear();
    std::pair< std::string, std::string > p = getCompoundNameAndStrand(result, genes);
    if (!p.first.empty() && !p.second.empty())
    {
        updateStrTags(record, TAG, p.first);
        updateStrTags(record, STRAND_TAG, p.second);
    }

    return 0;
}

// int TagReadsWithGeneExon::makeOverlapDetector()
// {
// 	spdlog::debug("Begin makeOverlapDetector");

// 	GTFReader gtfReader;
// 	std::unordered_map<std::string, std::vector<GTFRecord>> gatherByGeneName;
// 	try
// 	{
// 		gtfReader.loadGTFFile(annotation_filename, gatherByGeneName);
// 	}
// 	catch (AnnotationException& e)
// 	{
// 		spdlog::error(e.what());
// 		return -1;
// 	}
// 	catch (std::exception& e)
// 	{
// 		spdlog::error(e.what());
// 		return -1;
// 	}

// 	// Store key-value from chromosome to index, like "chr1":0, "chr2":1 ...
// 	contigs = gtfReader.getContigs();

// 	std::unordered_map<std::string, GTree::interval_vector> intervals;
// 	GeneBuilder geneBuilder;
// 	GTFIterator gtfIterator = gatherByGeneName.begin();
// 	for (; gtfIterator != gatherByGeneName.end(); ++gtfIterator)
// 	{
// 		try
// 		{
// 			GeneFromGTF gene = geneBuilder.makeGene(gtfIterator);
// 			intervals[gene.getContig()].push_back(GTree::interval(gene.getStart(), gene.getEnd(), gene));
// 		}
// 		catch (AnnotationException& e)
// 		{
// 			spdlog::debug("{} {}", e.what(), "-- skipping");
// 			continue;
// 		}

// 	}
// 	int numGene = 0;
// 	for (auto& p : intervals)
// 	{
// 		numGene += p.second.size();
// 		trees[p.first] = GTree(std::move(p.second), 16, 1);
// 	}
// 	spdlog::info("Gene count:{}", numGene);
// 	spdlog::debug("End makeOverlapDetector");

// 	lastRecord = createBamRecord();
// 	lastCache = std::make_tuple("", "", "", 0, 0, 0, 0);

// 	return 0;
// }

int TagReadsWithGeneExon::makeOverlapDetectorV2()
{
    spdlog::debug("Begin makeOverlapDetector");

    Timer                                                       load_timer;
    GTFReader                                                   gtfReader;
    std::unordered_map< std::string, std::vector< GTFRecord > > gatherByGeneName;
    try
    {
        gtfReader.loadGTFFile(annotation_filename, gatherByGeneName);
    }
    catch (AnnotationException& e)
    {
        spdlog::error(e.what());
        return -1;
    }
    catch (std::exception& e)
    {
        spdlog::error(e.what());
        return -1;
    }
    spdlog::info("loadGTFFile time(s):{:.2f}", load_timer.toc(1000));

    int         numGene = 0;
    GeneBuilder geneBuilder;
    GTFIterator gtfIterator = gatherByGeneName.begin();

    for (; gtfIterator != gatherByGeneName.end(); ++gtfIterator)
    {
        try
        {
            std::vector< GeneFromGTF > genes = geneBuilder.makeGene(gtfIterator);
            for (auto& gene : genes)
            {
                Node node;
                node.lower = gene.getStart();
                node.upper = gene.getEnd();
                node.value = gene;
                nodes.push_back(std::move(node));
                ++numGene;
            }
        }
        catch (AnnotationException& e)
        {
            spdlog::debug("{} {}", e.what(), "-- skipping");
            continue;
        }
    }
    // spdlog::debug("mytrees size:{}", mytrees.size());
    // for (auto& p : mytrees)
    // 	spdlog::debug("mytrees key:{}", p.first);
    for (auto& node : nodes)
        mytrees[node.value.getContig()].insert(node);
    spdlog::info("makeOverlapDetector time(s):{:.2f}", load_timer.toc(1000));

    spdlog::info("Gene count:{}", numGene);
    spdlog::debug("End makeOverlapDetector");

    lastRecord = createBamRecord();
    lastCache  = std::make_tuple("", "", "", 0, 0, 0, 0);

    return 0;
}

std::string TagReadsWithGeneExon::dumpMetrics()
{
    spdlog::info("TOTAL READS [{}] CORRECT_STRAND [{}]  WRONG_STRAND [{}] AMBIGUOUS_STRAND_FIXED [{}] AMBIGUOUS "
                 "REJECTED READS [{}]",
                 total_reads, reads_right_strand, reads_wrong_strand, read_ambiguous_gene_fixed,
                 ambiguous_reads_rejected);
    std::stringstream ss;
    ss << "## ANNOTATION METRICS\n";
    if (anno_ver != AnnoVersion::TENX)
    {
        ss << "TOTAL_READS\tREADS_WRONG_STRAND\tREADS_RIGHT_STRAND\tREAD_AMBIGUOUS_GENE_FIXED\tAMBIGUOUS_"
              "READS_REJECTED\n";
        ss << total_reads << "\t" << reads_wrong_strand << "\t" << reads_right_strand << "\t"
           << read_ambiguous_gene_fixed << "\t" << ambiguous_reads_rejected << std::endl;
    }
    else
    {
        ss << "TOTAL_READS\tMAP\tEXONIC\tINTRONIC\tINTERGENIC\tTRANSCRIPTOME\tANTISENSE\n";
        ss << total_reads << "\t" << map_reads << "\t" << exonic_reads << "\t" << intronic_reads << "\t"
           << intergenic_reads + nogene_reads << "\t" << transcriptome_reads << "\t" << reads_wrong_strand << std::endl;
        ss.flags(std::ios::fixed);
        ss.precision(1);
        ss << 100.0 << "\t" << map_reads * 100.0 / total_reads << "\t" << exonic_reads * 100.0 / total_reads << "\t"
           << intronic_reads * 100.0 / total_reads << "\t" << (intergenic_reads + nogene_reads) * 100.0 / total_reads
           << "\t" << transcriptome_reads * 100.0 / total_reads << "\t" << reads_wrong_strand * 100.0 / total_reads
           << "\t" << std::endl;
    }

    return ss.str();
}
