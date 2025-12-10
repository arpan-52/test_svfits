#include "config.hpp"
#include "utils.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

Config Config::fromFile(const std::string& filename) {
    Config config;
    config.setDefaults();

    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open config file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        config.parseLine(line);
    }

    config.validate();
    config.computeDerived();

    return config;
}

void Config::setDefaults() {
    // Defaults already set in class definition
}

void Config::parseLine(const std::string& line) {
    if (line.empty() || line[0] == '*' || line[0] == '#') return;

    std::string key, value;
    if (!parseKeyValue(line, key, value)) return;

    try {
        if (key == "NFILE") nfiles = std::stoi(value);
        else if (key == "PATH") path = value;
        else if (key == "INPUT") {
            input_files = split(value, ',');
        }
        else if (key == "FITS") fits_output = value;
        else if (key == "FREQ_SET") {
            auto parts = split(value, ':');
            if (parts.size() >= 3) {
                freq_start = std::stod(parts[0]);
                freq_end = std::stod(parts[1]);
                n_channels = std::stoi(parts[2]);
            }
        }
        else if (key == "BURST_NAME") burst_name = value;
        else if (key == "BURST_MJD") burst_mjd = std::stod(value);
        else if (key == "BURST_DM") burst_dm = std::stod(value);
        else if (key == "BURST_INTWD") burst_intwd = std::stod(value);
        else if (key == "BURST_FREQ") burst_freq = std::stod(value);
        else if (key == "RA_APP") ra_app = std::stod(value);
        else if (key == "DEC_APP") dec_app = std::stod(value);
        else if (key == "RA_MEAN") ra_mean = std::stod(value);
        else if (key == "DEC_MEAN") dec_mean = std::stod(value);
        else if (key == "EPOCH") epoch = std::stod(value);
        else if (key == "ANTMASK") antmask = std::stoul(value);
        else if (key == "DO_BAND") do_bandpass = (std::stoi(value) != 0);
        else if (key == "DO_BASE") do_baseline = (std::stoi(value) != 0);
        else if (key == "DO_FLAG") do_flag = (std::stoi(value) != 0);
        else if (key == "THRESH") flag_threshold = std::stof(value);
        else if (key == "NUM_THREADS") num_threads = std::stoi(value);
        else if (key == "IATUTC") iatutc = std::stof(value);
    } catch (const std::exception& e) {
        std::cerr << "Warning: Failed to parse " << key << " = " << value << std::endl;
    }
}

void Config::validate() {
    if (input_files.size() != static_cast<size_t>(nfiles)) {
        std::cerr << "Warning: Number of input files (" << input_files.size()
                  << ") doesn't match NFILE (" << nfiles << ")" << std::endl;
    }
    if (burst_mjd == 0.0) {
        throw std::runtime_error("BURST_MJD not set");
    }
}

void Config::computeDerived() {
    switch (ref_freq_mode) {
        case RefFreqMode::BURST_FREQ:
            ref_freq = burst_freq;
            break;
        case RefFreqMode::HIGHEST_FREQ:
            ref_freq = freq_end;
            break;
        case RefFreqMode::MANUAL:
            break;
    }

    std::cout << "Reference frequency: " << ref_freq / 1e6 << " MHz" << std::endl;
}
