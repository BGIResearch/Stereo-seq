/*
 * File: main.cpp
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI-Research
 */

#include "handleBam.h"
#include "utils.h"

#include <ctime>

#include <chrono>
#include <filesystem>
#include <iostream>
#include <string>
namespace fs = std::filesystem;

#include <CLI11.hpp>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

static const std::string version = "1.0.1";

int main(int argc, char** argv)
{
    Timer timer;
    // Parse the command line parameters.
    CLI::App app{ "HandleBam: mapping quality filter, deduplication, "
                  "set annotation, stat gene expression data." };
    app.footer("HandleBam version: " + version);
    app.get_formatter()->column_width(40);

    string input_bam, output_bam, annotation_file, metrics_file, exp_file;
    // Required parameters
    app.add_option("-I,-i", input_bam, "Input bam filename or file list separated by comma")->required();
    app.add_option("-O,-o", output_bam, "Output bam filename")->required();
    app.add_option("-A,-a", annotation_file, "Input annotation filename")->check(CLI::ExistingFile)->required();
    app.add_option("-S,-s", metrics_file, "Output summary filename")->required();
    app.add_option("-E,-e", exp_file, "Output barcode gene expression filename")->required();
    // Optional parameters
    int mapping_quality_threshold = 10;
    app.add_option("-Q,-q", mapping_quality_threshold, "Set mapping quality threshold, default 10")
        ->check(CLI::PositiveNumber);
    int cpu_cores = std::thread::hardware_concurrency();
    ;
    app.add_option("-C,-c", cpu_cores, "Set cpu cores, default detect")->check(CLI::PositiveNumber);

    bool save_low_quality = false;
    app.add_flag("--save_lq", save_low_quality, "Save low quality reads, default false");
    bool save_duplicate = false;
    app.add_flag("--save_dup", save_duplicate, "Save duplicate reads, default false");
    int annotation_mode = 2;
    app.add_option("--anno_mode", annotation_mode,
                   "Select annotation mode, default 2\n0->'v1'\n1->'v2'\n2->'v3'")
        ->check(CLI::Range(0, 2));
    bool umi_on      = false;
    auto umi_option  = app.add_flag("--umi_on", umi_on, "Open umi correction, default off");
    int  umi_min_num = 5;
    app.add_option("--umi_min_num", umi_min_num, "Minimum umi number for correction, default 5")
        ->check(CLI::PositiveNumber);
    int umi_mismatch = 1;
    app.add_option("--umi_mismatch", umi_mismatch, "Maximum mismatch for umi correction, default 1")
        ->check(CLI::PositiveNumber);
    std::string saturation_file = "";
    app.add_option("--sat_file", saturation_file, "Output sequencing saturation file, default None")->needs(umi_option);
    bool scrna            = false;
    auto scrna_option     = app.add_flag("--scrna,--scRNA,--SCRNA", scrna, "Set scRNA mode, default false");
    bool no_filter_matrix = false;
    app.add_flag("--no_filter_matrix", no_filter_matrix, "Not filter the gene expression matrix, default false")
        ->needs(scrna_option);

    CLI11_PARSE(app, argc, argv);

    // Check the input bam files
    vector< string > bam_lists = split_str(input_bam, ',');
    if (bam_lists.empty())
    {
        std::cerr << "Invalid parameter of -I: " << input_bam << std::endl;
        exit(-1);
    }
    else
    {
        for (auto& f : bam_lists)
        {
            if (!fs::exists(f))
            {
                std::cerr << "Not exists bam file: " << f << std::endl;
                exit(-1);
            }
        }
    }

    // Set the default logger to file logger.
    std::time_t        t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::ostringstream ostr;
    ostr << "HandleBam_" << std::put_time(std::localtime(&t), "%Y%m%d_%H%M%S") << ".log";
    try
    {
        // auto file_logger = spdlog::basic_logger_mt("main", "logs/" + ostr.str());
        auto file_sink      = std::make_shared< spdlog::sinks::basic_file_sink_mt >("logs/" + ostr.str());
        auto main_logger    = std::make_shared< spdlog::logger >("main", file_sink);
        auto gtf_logger     = std::make_shared< spdlog::logger >("gtf", file_sink);
        auto process_logger = std::make_shared< spdlog::logger >("process", file_sink);
        spdlog::register_logger(main_logger);
        spdlog::register_logger(gtf_logger);
        spdlog::register_logger(process_logger);
        spdlog::set_default_logger(process_logger);
    }
    catch (const spdlog::spdlog_ex& ex)
    {
        std::cerr << "Log init failed: " << ex.what() << std::endl;
    }
    spdlog::set_level(spdlog::level::info);  // Set global log level.
    // spdlog::cfg::load_argv_levels(argc, argv); // Set log level from command line
    spdlog::flush_on(spdlog::level::info);
    spdlog::set_pattern("%Y-%m-%d %H:%M:%S.%e %L %n: %v");

    spdlog::get("main")->info("{} INPUT={} OUTPUT={} SUMMARY={} TAG={} ANNOTATION_FILE={} SAVE_LOW_QUALITY={} "
                              "SAVE_DUPLICATE={} ANNOTATION_MODE={} UMI_ON={} UMI_MIN_NUM={} UMI_MISMATCH={} "
                              "SAT_FILE={} NO_FILTER_MATRIX={} CPU_CORES={} SCRNA={}",
                              argv[0], input_bam, output_bam, metrics_file, "GE", annotation_file, save_low_quality,
                              save_duplicate, annotation_mode, umi_on, umi_min_num, umi_mismatch, saturation_file,
                              no_filter_matrix, cpu_cores, scrna);
    HandleBam handleBam(bam_lists, output_bam, annotation_file, metrics_file, mapping_quality_threshold, exp_file);
    handleBam.setBamConfig(save_low_quality, save_duplicate, annotation_mode);
    handleBam.setUmiConfig(umi_on, umi_min_num, umi_mismatch);
    handleBam.setExtraConfig(saturation_file, !no_filter_matrix, cpu_cores, scrna);
    try
    {
        handleBam.doWork();
    }
    catch (std::exception& e)
    {
        spdlog::error("Error: {}", e.what());
    }
    catch (...)
    {
        spdlog::error("Unknown error");
    }

    spdlog::get("main")->info("HandleBam done. Elapsed time(s):{:.2f}", timer.toc(1000));

    return 0;
}
