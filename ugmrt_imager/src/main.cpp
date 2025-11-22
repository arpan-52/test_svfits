/**
 * @file main.cpp
 * @brief uGMRT Fast Imager - Command line interface
 *
 * Usage:
 *   ugmrt_imager -u param_file.txt -a antsamp.hdr [options]
 *
 * Options:
 *   -u FILE    Parameter file (svfits_par.txt format)
 *   -a FILE    Antenna/sampler header file
 *   -B FILE    Bulletin A file for DUT1 corrections (optional)
 *   -o FILE    Output FITS image filename
 *   -n NX,NY   Grid size (default: 512,512)
 *   -c CELL    Cell size in arcseconds (default: 1.0)
 *   -d DEVICE  Device: cuda, openmp, serial (default: cuda)
 *   -t THREADS Number of CPU threads for reading (default: 4)
 *   -h         Show help
 */

#include "ugmrt_imager.hpp"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <getopt.h>

void print_usage(const char* prog) {
    std::cout << "uGMRT Fast Imager v1.0\n\n";
    std::cout << "Usage: " << prog << " -u param_file -a antsamp_file [options]\n\n";
    std::cout << "Required:\n";
    std::cout << "  -u FILE    Parameter file (svfits_par.txt format)\n";
    std::cout << "  -a FILE    Antenna/sampler header file\n\n";
    std::cout << "Optional:\n";
    std::cout << "  -B FILE    Bulletin A file for DUT1 corrections\n";
    std::cout << "  -o FILE    Output FITS image (default: burst_image.fits)\n";
    std::cout << "  -n NX,NY   Grid size (default: 512,512)\n";
    std::cout << "  -c CELL    Cell size in arcseconds (default: 1.0)\n";
    std::cout << "  -s SUPPORT CF support in pixels (default: 7)\n";
    std::cout << "  -d DEVICE  Device: cuda, openmp, serial (default: cuda)\n";
    std::cout << "  -t THREADS CPU threads for reading (default: 4)\n";
    std::cout << "  -T THRESH  Flagging threshold in MAD (default: 5.0)\n";
    std::cout << "  --no-bandpass    Disable bandpass correction\n";
    std::cout << "  --no-baseline    Disable baseline subtraction\n";
    std::cout << "  --no-flag        Disable RFI flagging\n";
    std::cout << "  -h, --help       Show this help\n\n";
    std::cout << "Example:\n";
    std::cout << "  " << prog << " -u svfits_par.txt -a antsamp.hdr -o burst.fits -d cuda\n";
}

int main(int argc, char* argv[]) {
    ugmrt::ImagerConfig config;

    // Default values
    config.imaging.output_fits = "burst_image.fits";

    // Long options
    static struct option long_options[] = {
        {"no-bandpass", no_argument, nullptr, 1},
        {"no-baseline", no_argument, nullptr, 2},
        {"no-flag", no_argument, nullptr, 3},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int opt;
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "u:a:B:o:n:c:s:d:t:T:h", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'u':
                config.param_file = optarg;
                break;
            case 'a':
                config.antsamp_file = optarg;
                break;
            case 'B':
                config.bulletin_a = optarg;
                break;
            case 'o':
                config.imaging.output_fits = optarg;
                break;
            case 'n': {
                int nx, ny;
                if (sscanf(optarg, "%d,%d", &nx, &ny) == 2) {
                    config.grid.nx = nx;
                    config.grid.ny = ny;
                } else if (sscanf(optarg, "%d", &nx) == 1) {
                    config.grid.nx = config.grid.ny = nx;
                }
                break;
            }
            case 'c':
                config.grid.cell_size_asec = atof(optarg);
                break;
            case 's':
                config.grid.cf_support = atoi(optarg);
                break;
            case 'd':
                config.device = optarg;
                break;
            case 't':
                config.processing.num_threads = atoi(optarg);
                break;
            case 'T':
                config.processing.flag_threshold = atof(optarg);
                break;
            case 1:  // --no-bandpass
                config.processing.do_bandpass = false;
                break;
            case 2:  // --no-baseline
                config.processing.do_baseline = false;
                break;
            case 3:  // --no-flag
                config.processing.do_flag = false;
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    // Check required arguments
    if (config.param_file.empty() || config.antsamp_file.empty()) {
        std::cerr << "Error: -u and -a options are required\n\n";
        print_usage(argv[0]);
        return 1;
    }

    // Print configuration
    std::cout << "========================================\n";
    std::cout << "  uGMRT Fast Imager v1.0\n";
    std::cout << "========================================\n";
    std::cout << "Configuration:\n";
    std::cout << "  Parameter file: " << config.param_file << "\n";
    std::cout << "  Antenna file: " << config.antsamp_file << "\n";
    std::cout << "  Output: " << config.imaging.output_fits << "\n";
    std::cout << "  Grid: " << config.grid.nx << "x" << config.grid.ny << "\n";
    std::cout << "  Cell size: " << config.grid.cell_size_asec << " arcsec\n";
    std::cout << "  Device: " << config.device << "\n";
    std::cout << "  Threads: " << config.processing.num_threads << "\n";
    std::cout << "  Bandpass: " << (config.processing.do_bandpass ? "yes" : "no") << "\n";
    std::cout << "  Baseline: " << (config.processing.do_baseline ? "yes" : "no") << "\n";
    std::cout << "  Flagging: " << (config.processing.do_flag ? "yes" : "no");
    if (config.processing.do_flag) {
        std::cout << " (threshold: " << config.processing.flag_threshold << " MAD)";
    }
    std::cout << "\n";
    std::cout << "========================================\n\n";

    // Create and run imager
    ugmrt::UGMRTImager imager(config);

    if (!imager.initialize()) {
        std::cerr << "Failed to initialize imager\n";
        return 1;
    }

    if (!imager.run()) {
        std::cerr << "Failed to run imager\n";
        return 1;
    }

    // Print statistics
    auto stats = imager.get_stats();
    std::cout << "\n========================================\n";
    std::cout << "Statistics:\n";
    std::cout << "  Total visibilities: " << stats.total_visibilities << "\n";
    std::cout << "  Flagged: " << stats.flagged_visibilities
              << " (" << 100.0 * stats.flagged_visibilities / stats.total_visibilities << "%)\n";
    std::cout << "  Gridded: " << stats.gridded_visibilities << "\n";
    std::cout << "========================================\n";

    std::cout << "\nDone!\n";
    return 0;
}
