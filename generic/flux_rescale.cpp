#include <phypp.hpp>

void print_help();

int phypp_main(int argc, char* argv[]) {
    vec1s  files;                 // 1D spectra to analyze
    bool   verbose = false;       // print progress in the standard output
    double min_snr = 3.0;         // minimum S/N on the rescaling factor
    vec1u exclude_hdu;            // HDUs for which rescaling should not be applied

    bool help = false;

    read_args(argc, argv, arg_list(files, verbose, exclude_hdu, help));

    if (help) {
        print_help();
        return 0;
    }

    if (files.empty()) {
        error("no spectrum to stack (files=...)");
        return 1;
    }

    // Read and stack fluxes
    vec1d mflx, merr;
    vec1d cflx, cerr;
    vec1s bands;
    for (uint_t i : range(files)) {
        vec1d tflx, terr, tcflx, tcerr;

        fits::read_table(file::remove_extension(files[i])+"_broadband.fits",
            "flux", tflx, "flux_err", terr, "catalog_flux", tcflx, "catalog_flux_err", tcerr,
            "bands", bands
        );

        if (mflx.empty()) {
            mflx = tflx;
            merr = terr;
            cflx = tcflx;
            cerr = tcerr;
        } else {
            mflx += tflx;
            merr = sqrt(sqr(merr) + sqr(terr));
            cflx += tcflx;
            cerr = sqrt(sqr(cerr) + sqr(tcerr));
        }
    }

    // Compute rescaling
    vec1d fratio(mflx.size());
    vec1d eratio(mflx.size());
    for (uint_t b : range(mflx)) {
        fratio[b] = cflx[b]/mflx[b];
        eratio[b] = sqrt(sqr(cerr[b]/mflx[b]) + sqr(cflx[b]*merr[b]/sqr(mflx[b])));

        if (verbose) {
            note(bands[b], " ", fratio[b], " +/- ", eratio[b]);
        }
    }

    vec1u idb = where(!is_finite(eratio));
    fratio[idb] = 1.0;
    eratio[idb] = 1e20;

    auto p = optimal_mean(fratio, eratio);
    if (verbose) {
        note(" --> ", p.first, " +/- ", p.second, " (S/N=", p.first/p.second, ")");
    }

    // Apply rescaling
    for (uint_t i : range(files)) {
        file::copy(files[i], file::remove_extension(files[i])+"_rescaled.fits");
        fits::image img(file::remove_extension(files[i])+"_rescaled.fits");

        if (p.first/p.second >= min_snr) {
            for (uint_t u : range(img.hdu_count())) {
                if (is_any_of(u, exclude_hdu)) continue;

                img.reach_hdu(u);
                if (!img.image_dims().empty()) {
                    vec1d data;
                    img.read(data);
                    data *= p.first;
                    img.update(data);
                }
            }
        }
    }

    return 0;
}

void print_help() {
    using namespace terminal_format;

    print("flux_rescale v1.0");
    print("usage: flux_rescale files=[file1.fits,file2.fits,...] [options]");
    print("");
    print("Main parameters:");
    paragraph("The files listed in 'files=[...]' must be valid 1D spectra FITS files, for example "
        "created by 'extract2d', and one for each source in a given slit. It is expected that the "
        "fluxes are located in the first HDU, and uncertainties in the second HDU. For each file in "
        "this list, there must be a corresponding <file>_broadband.fits file in the same directory, "
        "which contains the measured synthetic fluxes (e.g., obtained with 'get_fluxes'). The "
        "program will then rescale all the data in the file so that the synthetic broadband fluxes "
        "agree as best as possible with their measured fluxes in the catalog. The rescaling factor "
        "is determined as an average value, and is assumed to be the same for all the sources whose "
        "spectra re listed in 'files=[...]'.");

    print("");
    print("Output format:");
    paragraph("The program will create a file called <file>_rescaled.fits for each spectrum in "
        "the list 'files=[...]'. This file has the same format as the original file, but with the "
        "flux rescaled.");

    print("");
    print("Options:");
    bullet("files=[...]", "Must be a list of files. This is the list of 1D spectra that will be "
        "rescaled, one spectrum per source. The way this tool is meant to be used, this list should "
        "correspond to all the sources detected in a given slit for one OB (or one exposure).");
    bullet("min_snr=...", "Must be a number. This defines the minimum allowed S/N on the rescaling "
        "factor. If the S/N falls below this threshold, no rescaling is attempted. Default is 3.");
    bullet("exclude_hdu=[...]", "Must be a list of integers. The HDUs listed here will not be "
        "affected by the rescaling (e.g., useful for flags)");
    bullet("verbose", "Set this flag to print the derived rescaling factors in the terminal.");
}
