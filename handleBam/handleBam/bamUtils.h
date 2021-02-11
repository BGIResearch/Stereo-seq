/*
 * File: bamUtils.h
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#pragma once

#include <vector>

/**
 * Represents the contiguous alignment of a subset of read bases to a reference
 * sequence. Simply put an alignment block tells you that read bases from
 * readStart are aligned to the reference (matching or mismatching) from
 * referenceStart for length bases.
 */
class AlignmentBlock
{
public:
    /** Constructs a new alignment block with the supplied read and ref starts and length. */
    AlignmentBlock() {}

    AlignmentBlock(int readStart_, int referenceStart_, int length_)
        : readStart(readStart_), referenceStart(referenceStart_), length(length_)
    {
    }

    /** The first, 1-based, base in the read that is aligned to the reference reference. */
    int getReadStart()
    {
        return readStart;
    }

    /** The first, 1-based, position in the reference to which the read is aligned. */
    int getReferenceStart()
    {
        return referenceStart;
    }

    /** The number of contiguous bases aligned to the reference. */
    int getLength()
    {
        return length;
    }

private:
    int readStart;
    int referenceStart;
    int length;
};

/**
 * Given a Cigar, Returns blocks of the sequence that have been aligned directly to the
 * reference sequence. Note that clipped portions, and inserted and deleted bases (vs. the reference)
 * are not represented in the alignment blocks.
 *
 * @param cigar          The cigar containing the alignment information
 * @param alignmentStart The start (1-based) of the alignment
 * @param cigarTypeName  The type of cigar passed - for error logging.
 * @return List of alignment blocks
 */
std::vector< AlignmentBlock > getAlignmentBlocks(std::vector< std::pair< int, int > >& cigars, int alignmentStart);

int getReferenceLength(std::vector< std::pair< int, int > >& cigars);