/*
 * File: geneBuilder.cpp
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#include <climits>

#include <algorithm>
#include <set>

#include "annotationException.h"
#include "geneBuilder.h"

std::vector< GeneFromGTF > GeneBuilder::makeGene(GTFIterator& iter)
{
    std::vector< GTFRecord > recordsForGene = iter->second;
    return makeGeneFromMultiVersionGTFRecords(recordsForGene);
}

std::vector< GeneFromGTF > GeneBuilder::makeGeneFromMultiVersionGTFRecords(std::vector< GTFRecord >& gtfRecords)
{
    GTFMapForGeneVersion gtfRecordsByGeneVersion = gatherByGeneVersion(gtfRecords);
    int                  maxVersion              = -1;
    for (auto& p : gtfRecordsByGeneVersion)
        maxVersion = std::max(maxVersion, p.first);
    return makeGenesWithTranscriptsFromGTFRecords(gtfRecordsByGeneVersion[maxVersion]);
}

GTFMapForGeneVersion GeneBuilder::gatherByGeneVersion(std::vector< GTFRecord >& gtfRecords)
{
    GTFMapForGeneVersion res;
    for (auto& record : gtfRecords)
    {
        res[record.getGeneVersion()].push_back(record);
    }
    return res;
}

GTFMapForTranscript GeneBuilder::gatherByTransciptID(std::vector< GTFRecord >& gtfRecords)
{
    GTFMapForTranscript res;
    for (auto& record : gtfRecords)
    {
        // Romove featureType == gene when makeing transcripts.
        if (record.getFeatureType() == "gene")
            continue;
        if (record.getTranscriptID() == "")
        {
            std::string error = "GTFRecord does not have transcriptID: " + record.getGeneName();
            throw AnnotationException(error);
        }
        res[record.getTranscriptID()].push_back(record);
    }
    return res;
}

bool GeneBuilder::geneInDiffChrs(std::vector< GTFRecord >& gtfRecords)
{
    std::string chr(gtfRecords[0].getContig());
    for (auto& r : gtfRecords)
    {
        if (r.getContig() != chr)
            return true;
    }
    return false;
}

std::vector< GeneFromGTF > GeneBuilder::makeGenesWithTranscriptsFromGTFRecords(std::vector< GTFRecord >& gtfRecords)
{
    // Check is there same gene in different chromosomes
    if (!geneInDiffChrs(gtfRecords))
    {
        return { makeGeneWithTranscriptsFromGTFRecords(gtfRecords) };
    }
    else
    {
        std::vector< GeneFromGTF > res;
        // Split input gtf records into multi fragments
        std::vector< GTFRecord > temp;
        for (auto& r : gtfRecords)
        {
            if (!temp.empty() && r.getContig() != temp.front().getContig())
            {
                res.push_back(makeGeneWithTranscriptsFromGTFRecords(temp));
                temp.clear();
            }
            temp.push_back(r);
        }
        if (!temp.empty())
        {
            res.push_back(makeGeneWithTranscriptsFromGTFRecords(temp));
            temp.clear();
        }
        return res;
    }
}

GeneFromGTF GeneBuilder::makeGeneWithTranscriptsFromGTFRecords(std::vector< GTFRecord >& gtfRecords)
{
    GeneFromGTF gene = makeGeneFromGTFRecords(gtfRecords);

    GTFMapForTranscript gtfLinesByTranscript = gatherByTransciptID(gtfRecords);
    for (auto& p : gtfLinesByTranscript)
        addTranscriptToGeneFromGTFRecords(gene, p.second);

    if (gene.getTranscripts().size() == 0)
    {
        std::string error = "No transcript in GTF for gene " + gene.getName();
        throw AnnotationException(error);
    }
    return gene;
}

std::string setToString(std::set< std::string >& ss, std::string sep = ", ")
{
    std::string res;
    for (auto& s : ss)
    {
        if (!res.empty())
            res += sep;
        res += s;
    }
    return res;
}
GeneFromGTF GeneBuilder::makeGeneFromGTFRecords(std::vector< GTFRecord >& gtfRecords)
{
    auto& geneRecord = gtfRecords[0];

    // Figure out the extend of the gene.
    int                     start = INT_MAX;
    int                     end   = INT_MIN;
    std::set< std::string > geneIds, chromosomes;
    for (auto& r : gtfRecords)
    {
        start = std::min(start, r.getStart());
        end   = std::max(end, r.getEnd());
        geneIds.insert(r.getGeneID());
        chromosomes.insert(r.getContig());
    }
    // Keep the same genes in different chromosomes
    // if (chromosomes.size() > 1)
    // {
    //     std::string error = "Chromosome disagreement(" + setToString(chromosomes) + ") in GTF file for gene "
    //                         + geneRecord.getGeneName();
    //     throw AnnotationException(error);
    // }
    if (geneIds.size() > 1)
    {
        std::string error = "Multiple gene IDs for gene " + geneRecord.getGeneName() + ": " + setToString(geneIds);
        throw AnnotationException(error);
    }
    GeneFromGTF gene(geneRecord.getContig(), start, end, geneRecord.isNegativeStrand(), geneRecord.getGeneName(),
                     geneRecord.getFeatureType(), geneRecord.getGeneID(), geneRecord.getTranscriptType(),
                     geneRecord.getGeneVersion());

    for (auto& r : gtfRecords)
        validateGTFRecord(r, gene);
    return gene;
}

void GeneBuilder::validateGTFRecord(GTFRecord& record, GeneFromGTF& gene)
{
    // TODO(fxzhao): maybe i can ignore those compare.
    if (gene.getContig() != record.getContig())
    {
        std::string error = "Strand disagreement in GTF file for gene " + gene.getName();
        throw AnnotationException(error);
    }
    if (gene.getContig() != record.getContig())
    {
        std::string error = "Chromosome disagreement(" + gene.getContig() + " != " + record.getContig()
                            + ") in GTF file for gene " + gene.getName();
        throw AnnotationException(error);
    }
    if (record.getFeatureType() == "gene")
    {
        if (record.getStart() != gene.getStart() || record.getEnd() != gene.getEnd())
            throw AnnotationException("gene GTFRecord != GeneFromGTF");
    }
}

// Conversion from 0-based half-open to 1-based inclusive
// intervals is done here.
void GeneBuilder::addTranscriptToGeneFromGTFRecords(GeneFromGTF& gene, std::vector< GTFRecord >& gtfRecords)
{
    auto&       geneName              = gene.getName();
    auto&       lineOne               = gtfRecords[0];
    auto&       transcriptName        = lineOne.getTranscriptName();
    auto&       transcriptID          = lineOne.getTranscriptID();
    std::string transcriptDescription = geneName + ":" + transcriptName;
    auto&       transcriptType        = lineOne.getTranscriptType();

    std::vector< Exon > exons;
    int                 transcriptionStart = INT_MAX;
    int                 transcriptionEnd   = INT_MIN;
    int                 codingStart        = INT_MAX;
    int                 codingEnd          = INT_MIN;
    for (auto& r : gtfRecords)
    {
        std::string featureType = r.getFeatureType();
        int         start       = r.getStart();
        int         end         = r.getEnd();
        if (featureType == "exon")
        {
            Exon e(start, end);
            exons.push_back(e);
            transcriptionStart = std::min(transcriptionStart, start);
            transcriptionEnd   = std::max(transcriptionEnd, end);
        }
        else if (featureType == "CDS")
        {
            codingStart = std::min(codingStart, start);
            codingEnd   = std::max(codingEnd, end);
        }
    }

    std::sort(exons.begin(), exons.end(), compareExon);
    if (codingStart == INT_MAX)
        codingStart = transcriptionStart;
    if (codingEnd == INT_MIN)
        codingEnd = transcriptionEnd;

    TranscriptFromGTF* transcript = gene.addTranscript(transcriptionStart, transcriptionEnd, codingStart, codingEnd,
                                                       exons.size(), transcriptName, transcriptID, transcriptType);
    for (size_t i = 0; i < exons.size(); ++i)
    {
        auto& e = exons[i];
        if (e.start > e.end)
        {
            std::string error = "Exon has 0 or negative extent for " + transcriptDescription;
            throw AnnotationException(error);
        }
        if (i > 0 && exons[i - 1].end > e.start)
        {
            std::string error = "Exons overlap for " + transcriptDescription;
            throw AnnotationException(error);
        }
    }

    transcript->addExons(exons);
}
