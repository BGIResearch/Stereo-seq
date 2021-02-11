/*
 * File: bamCat.cpp
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#include <string.h>

#include "bamCat.h"
#include "spdlog/spdlog.h"

/*
 * Copy from bam_cat.c of samtools, and modify it.
 */

#define BUF_SIZE 0x10000

#define GZIPID1 31
#define GZIPID2 139

#define BGZF_EMPTY_BLOCK_SIZE 28

int bam_cat(std::vector< std::string > fn, bam_hdr_t* h, const char* outbam, [[maybe_unused]] char* arg_list,
            [[maybe_unused]] int no_pg)
{
    BGZF *              fp, *in = NULL;
    uint8_t*            buf = NULL;
    uint8_t             ebuf[BGZF_EMPTY_BLOCK_SIZE];
    const int           es = BGZF_EMPTY_BLOCK_SIZE;
    decltype(fn.size()) i;

    fp = strcmp(outbam, "-") ? bgzf_open(outbam, "w") : bgzf_fdopen(fileno(stdout), "w");
    if (fp == 0)
    {
        spdlog::error("bam_cat fail to open output file:{}", outbam);
        return -1;
    }
    if (h)
    {
        // if (!no_pg && sam_hdr_add_pg(h, "samtools",
        //                              "VN", samtools_version(),
        //                              arg_list ? "CL": NULL,
        //                              arg_list ? arg_list : NULL,
        //                              NULL))
        //     goto fail;

        if (bam_hdr_write(fp, h) < 0)
        {
            spdlog::error("bam_cat couldn't write header");
            goto fail;
        }
    }

    buf = ( uint8_t* )malloc(BUF_SIZE);
    if (!buf)
    {
        fprintf(stderr, "[%s] Couldn't allocate buffer\n", __func__);
        goto fail;
    }
    for (i = 0; i < fn.size(); ++i)
    {
        bam_hdr_t* old;
        int        len, j;

        in = strcmp(fn[i].c_str(), "-") ? bgzf_open(fn[i].c_str(), "r") : bgzf_fdopen(fileno(stdin), "r");
        if (in == 0)
        {
            spdlog::error("bam_cat fail to open file: {}", fn[i]);
            goto fail;
        }
        if (in->is_write)
            return -1;

        old = bam_hdr_read(in);
        if (old == NULL)
        {
            fprintf(stderr, "[%s] ERROR: couldn't read header for '%s'.\n", __func__, fn[i].c_str());
            goto fail;
        }
        if (h == 0 && i == 0)
        {
            // if (!no_pg && sam_hdr_add_pg(old, "samtools",
            //                              "VN", samtools_version(),
            //                              arg_list ? "CL": NULL,
            //                              arg_list ? arg_list : NULL,
            //                              NULL))
            //     goto fail;

            if (bam_hdr_write(fp, old) < 0)
            {
                spdlog::error("bam_cat couldn't write header");
                goto fail;
            }
        }

        if (in->block_offset < in->block_length)
        {
            if (bgzf_write(fp, ( char* )in->uncompressed_block + in->block_offset, in->block_length - in->block_offset)
                < 0)
                goto write_fail;
            if (bgzf_flush(fp) != 0)
                goto write_fail;
        }

        j = 0;
        while ((len = bgzf_raw_read(in, buf, BUF_SIZE)) > 0)
        {
            if (len < es)
            {
                int diff = es - len;
                if (j == 0)
                {
                    fprintf(stderr, "[%s] ERROR: truncated file?: '%s'.\n", __func__, fn[i].c_str());
                    goto fail;
                }
                if (bgzf_raw_write(fp, ebuf, len) < 0)
                    goto write_fail;

                memcpy(ebuf, ebuf + len, diff);
                memcpy(ebuf + diff, buf, len);
            }
            else
            {
                if (j != 0)
                {
                    if (bgzf_raw_write(fp, ebuf, es) < 0)
                        goto write_fail;
                }
                len -= es;
                memcpy(ebuf, buf + len, es);
                if (bgzf_raw_write(fp, buf, len) < 0)
                    goto write_fail;
            }
            j = 1;
        }

        /* check final gzip block */
        {
            const uint8_t  gzip1 = ebuf[0];
            const uint8_t  gzip2 = ebuf[1];
            const uint32_t isize = *(( uint32_t* )(ebuf + es - 4));
            if (((gzip1 != GZIPID1) || (gzip2 != GZIPID2)) || (isize != 0))
            {
                fprintf(stderr, "[%s] WARNING: Unexpected block structure in file '%s'.", __func__, fn[i].c_str());
                fprintf(stderr, " Possible output corruption.\n");
                if (bgzf_raw_write(fp, ebuf, es) < 0)
                    goto write_fail;
            }
        }
        bam_hdr_destroy(old);
        bgzf_close(in);
        in = NULL;
    }
    free(buf);
    if (bgzf_close(fp) < 0)
    {
        fprintf(stderr, "[%s] Error on closing '%s'.\n", __func__, outbam);
        return -1;
    }
    return 0;

write_fail:
    fprintf(stderr, "[%s] Error writing to '%s'.\n", __func__, outbam);
fail:
    if (in)
        bgzf_close(in);
    if (fp)
        bgzf_close(fp);
    free(buf);
    return -1;
}