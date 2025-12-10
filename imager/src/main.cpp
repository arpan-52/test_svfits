#include "config.hpp"
#include "reader.hpp"
#include "kernel.hpp"
#include "gridder.hpp"
#include "cleaner.hpp"
#include "fits_output.hpp"
#include <iostream>
#include <exception>

int main(int argc, char** argv) {
    try {
        std::cout << "========================================" << std::endl;
        std::cout << "  uGMRT Burst Imager with W-Projection" << std::endl;
        std::cout << "========================================" << std::endl;

        // 1. Load configuration
        std::string config_file = (argc > 1) ? argv[1] : "svfits_par.txt";
        std::cout << "\n[1] Loading configuration from " << config_file << std::endl;
        Config config = Config::fromFile(config_file);

        // 2. Initialize reader
        std::cout << "\n[2] Initializing reader..." << std::endl;
        Reader reader(config);
        reader.initialize();

        // 3. Read visibilities
        std::cout << "\n[3] Reading visibilities..." << std::endl;
        auto visibilities = reader.readAllVisibilities();
        std::cout << "Loaded " << visibilities.size() << " visibilities" << std::endl;

        if (visibilities.empty()) {
            std::cerr << "No visibilities loaded!" << std::endl;
            return 1;
        }

        // 4. Generate W-kernels
        std::cout << "\n[4] Generating W-projection kernels..." << std::endl;
        KernelGenerator kernel_gen(config);
        auto kernel = kernel_gen.generateWProjection(reader.getMaxW());

        // 5. Initialize gridder
        std::cout << "\n[5] Initializing gridder..." << std::endl;
        Gridder gridder(config, kernel);

        // 6. Grid visibilities
        std::cout << "\n[6] Gridding visibilities..." << std::endl;
        auto grid = gridder.createGrid();
        gridder.gridVisibilities(grid, visibilities);
        gridder.normalize(grid);

        // 7. Compute PSF
        std::cout << "\n[7] Computing PSF..." << std::endl;
        auto psf = gridder.computePSF(visibilities);

        // 8. FFT to image domain
        std::cout << "\n[8] Transforming to image domain..." << std::endl;
        gridder.fft(grid);
        gridder.fftshift(grid);

        // 9. CLEAN deconvolution
        std::cout << "\n[9] Running CLEAN..." << std::endl;
        Cleaner cleaner(config);
        auto components = cleaner.clean(grid, psf);
        auto restored = cleaner.restore(components, grid);

        // 10. Write FITS output
        std::cout << "\n[10] Writing FITS output..." << std::endl;
        writeFITS(config.fits_output, restored, config);

        std::cout << "\n========================================" << std::endl;
        std::cout << "  Imaging Complete!" << std::endl;
        std::cout << "  Output: " << config.fits_output << std::endl;
        std::cout << "========================================" << std::endl;

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}
