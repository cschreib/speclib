#include <phypp.hpp>

vec1d cure_error(vec1d err, vec1d err_formal) {
    vec1u idl = where(err > 0);
    if (!idl.empty()) {
        err_formal[idl] *= median(err[idl]/err_formal[idl]);
    }

    // Keep the bootstrap unless it is too small for some reason
    idl = where(err < 0.7*err_formal);
    err[idl] = err_formal[idl];

    return err;
}

void print_help();

int phypp_main(int argc, char* argv[]) {
    std::string out;

    vec1s  files;           // 1D spectra to stack
    double rescale = 1.0;   // rescaling factor to apply to flux and errors
    vec1u  rebin = {1};     // spectral binnings (1: no binning, x: bin 'x' pixels into one)
    bool   verbose = false; // print progress in the standard output
    vec1s  filters;         // filters to use to build broadband fluxes
    std::string filter_db = data_dir+"fits/filter-db/db.dat";

    // Weighting scheme for stacking spectra and spectral binning
    // "optimal": inverse variance weighting, "uniform": no weighting
    std::string stack_weight = "optimal";
    std::string stack_op = "mean";

    bool do_sigma_clip = true;          // enable/disable sigma clipping of outliers
    double sigma_clip_threshold = 5.0;  // significance threshold for rejecting outliers
    uint_t sigma_clip_width = 1;        // width (in pixels) of the wavelength bin in which to define outliers

    bool save_cubes = false;            // save the "cube" of 1D spectra prior to stacking them

    bool help = false;

    read_args(argc, argv, arg_list(
        out, files, stack_weight, stack_op, filters,
        do_sigma_clip, sigma_clip_threshold, sigma_clip_width, save_cubes, rebin,
        rescale, verbose, filter_db, help
    ));

    if (help) {
        print_help();
        return 0;
    }

    file::mkdir(file::get_directory(out));

    if (stack_weight != "optimal" && stack_weight != "uniform") {
        error("unknown weighting scheme '", stack_weight, "' for stacking");
        return 1;
    }

    if (files.empty()) {
        error("no spectrum to stack (files=...)");
        return 1;
    }

    vec<1,filter_t> fbb;
    if (!filters.empty()) {
        auto db = read_filter_db(filter_db);

        if (!get_filters(db, filters, fbb)) {
            return 1;
        }

        rebin.push_back(npos);
    }

    vec2d cflx, cerr;

    // Read spectra
    vec1d lam;
    vec1s exp_band(files.size());
    vec1d exp_seeing(files.size());
    vec1d exp_offset(files.size());
    for (uint_t i : range(files)) {
        vec1d tflx, terr;

        fits::input_image iimg(files[i]);
        iimg.reach_hdu(1);
        iimg.read(tflx);

        vec1d tlam = astro::build_axis(astro::wcs(iimg.read_header()), 0, astro::axis_unit::wave_um);
        if (cflx.empty()) {
            lam = tlam;
        } else {
            if (count(abs(tlam - lam)/mean(lam) > 1e-6) > 0) {
                error("mismatch in wavelength range, cannot stack");
                return 1;
            }
        }

        iimg.read_keyword("BAND", exp_band[i]);
        iimg.read_keyword("GWIDTH", exp_seeing[i]);
        exp_seeing[i] *= 2.335; // sigma to FWHM (NB: pixel units, not arcsec!)
        if (verbose) {
            note(files[i]);
            note(" seeing: ", exp_seeing[i], " pixels, band: ", exp_band[i]);
        }

        iimg.read_keyword("GPOS", exp_offset[i]);

        iimg.reach_hdu(2);
        iimg.read(terr);

        if (cflx.empty()) {
            cflx.resize(tflx.dims, files.size());
            cerr.resize(cflx.dims);
        }

        tflx *= rescale;
        terr *= rescale;

        cflx(_,i) = tflx;
        cerr(_,i) = terr;
    }

    // Define weights
    vec2d cwei;
    if (stack_weight == "optimal") {
        // Inverse variance weighting
        cwei = 1/sqr(cerr);
        // Adjust global normalization to avoid numerical errors
        cwei /= max(cwei[where(is_finite(cwei))]);
    } else if (stack_weight == "uniform") {
        cwei = replicate(1.0, cflx.dims);
    }

    // Apply sigma clipping (if asked)
    vec2b crej(cflx.dims);
    if (do_sigma_clip) {
        uint_t d1 = (sigma_clip_width-1)/2;
        uint_t d2 = sigma_clip_width-1 - d1;

        for (uint_t l : range(cflx.dims[0])) {
            // Define wavelength region
            uint_t l0 = (l > d1 ?              l - d1 : 0);
            uint_t l1 = (l < cflx.dims[0]-d2 ? l + d2 : cflx.dims[0]-1);

            // First compute the weighted median (which we assume is unbiased)
            vec2d med(l1-l0+1, cflx.dims[1]);
            for (uint_t ll : range(med.dims[0])) {
                med(ll,_) = weighted_median(cflx(l0+ll,_), cwei(l0+ll,_));
            }

            // Compute the absolute residuals
            vec2d res = abs(cflx(l0-_-l1,_) - med);

            // Compute the RMS of these using the MAD (which we assume is unbiased)
            double rms = 1.48*median(res);

            // Select significant outliers in this wavelength element
            crej(l,_) = res(d1,_) > sigma_clip_threshold*rms;
        }

        if (verbose) {
            uint_t nfinite = count(is_finite(cflx));
            uint_t cnt = count(crej);
            note(cnt, "/", nfinite, " elements sigma clipped (",
                round(10.0*100.0*cnt/float(nfinite))/10.0, "%)");
        }
    }

    // Determine information about exposures
    vec1d exp_fclip(cflx.dims[1]);
    vec1d exp_wei(cflx.dims[1]);
    for (uint_t i : range(cflx.dims[1])) {
        vec1u idl = where(is_finite(cflx(_,i)));
        exp_fclip[i] = fraction_of(crej(idl,i));
        exp_wei[i] = median(cwei(idl,i));
    }

    // Rescale weights for display purposes (max of 1)
    for (auto b : unique_values(exp_band)) {
        vec1u idb = where(exp_band == b);
        exp_wei[idb] /= max(exp_wei[idb]);
    }

    // Define function to save this data into the FITS headers
    auto write_exposures = [&](fits::output_image& oimg) {
        oimg.write_keyword("NEXP", cflx.dims[1]);
        for (uint_t i : range(cflx.dims[1])) {
            std::string base = "E"+align_right(to_string(i), 4, '0');
            oimg.write_keyword(base+"SRC", file::get_basename(files[i]));
            oimg.write_keyword(base+"BND", exp_band[i], "observing band");
            oimg.write_keyword(base+"SNG", round(1000*exp_seeing[i])/1000.0, "seeing (arcsec)");
            oimg.write_keyword(base+"OFF", round(1000*exp_offset[i])/1000.0, "slit offset (pixels)");
            oimg.write_keyword(base+"WEI", round(1000*exp_wei[i])/1000.0, "average weight");
            oimg.write_keyword(base+"REJ", round(1000*exp_fclip[i])/1000.0, "fraction of rejected pixels");
        }
    };

    // Save cubes (if asked)
    if (save_cubes) {
        fits::output_image oimg(file::remove_extension(out)+"_cubes.fits");
        oimg.write_empty();

        auto write_wcs = [&]() {
            oimg.write_keyword("CTYPE1", "WAVE");
            oimg.write_keyword("CUNIT1", "um");
            oimg.write_keyword("CRPIX1", 1.0);
            oimg.write_keyword("CRVAL1", lam[0]);
            oimg.write_keyword("CDELT1", lam[1]-lam[0]);
            oimg.write_keyword("CTYPE2", "EPOCH");
            oimg.write_keyword("CRPIX2", 1.0);
            oimg.write_keyword("CRVAL2", 1.0);
            oimg.write_keyword("CDELT2", 1.0);

            write_exposures(oimg);
        };

        oimg.reach_hdu(1);
        oimg.write(transpose(cflx));
        write_wcs();
        oimg.reach_hdu(2);
        oimg.write(transpose(cerr));
        write_wcs();
        oimg.reach_hdu(3);
        oimg.write(transpose(cwei));
        write_wcs();
        oimg.reach_hdu(4);
        oimg.write(transpose(crej));
        write_wcs();
    }

    // Down-weight bad pixels
    {
        vec1u idb = where(!is_finite(cflx) || !is_finite(cerr) || !is_finite(cwei) || crej);
        cwei[idb] = 0; cflx[idb] = 0; cerr[idb] = 0;
    }

    // Resample the filters to the grid of the spectra
    for (auto& f : fbb) {
        uint_t i0 = where_first(lam > f.lam.front());
        uint_t i1 = where_last(lam < f.lam.back());

        if (i0 != npos && i1 != npos) {
            f.res = interpolate(f.res, f.lam, lam);
            f.res[_-i0] = 0;
            f.res[i1-_] = 0;
            f.lam = lam;
        } else {
            // Filter is not covered
            f.res.clear();
            f.lam.clear();
        }
    }

    auto stack_elements = [&](const vec2d& fdata, const vec2d& edata, const vec2d& wdata,
        uint_t l0, uint_t l1, double& flx, double& err, double& errb) {

        double wei = 0;
        err = 0;
        errb = 0;

        vec1d iflx(fdata.dims[1]);
        vec1d iwei(fdata.dims[1]);

        bool noflx = false;
        if (stack_op == "median") {
            flx = weighted_median(fdata.safe(l0-_-l1,_), wdata.safe(l0-_-l1,_));
            for (uint_t i : range(fdata.dims[1])) {
                iflx.safe[i] = weighted_median(fdata.safe(l0-_-l1,i), wdata.safe(l0-_-l1,i));
            }
            noflx = true;
        } else {
            flx = 0;
        }

        for (uint_t l = l0; l <= l1; ++l)
        for (uint_t i : range(fdata.dims[1])) {
            wei += wdata.safe(l,i);
            err += sqr(wdata.safe(l,i)*edata.safe(l,i));
            iwei.safe[i] += wdata.safe(l,i);

            if (!noflx) {
                flx += wdata.safe(l,i)*fdata.safe(l,i);
                iflx.safe[i] += wdata.safe(l,i)*fdata.safe(l,i);
            }
        }

        err = sqrt(err)/wei;

        if (!noflx) {
            flx /= wei;
            iflx /= iwei;
        }

        // For sample-estimated uncertainty, be sure to compute it on the binned data
        // rather than on all the data at once
        wei = 0;
        uint_t npt = 0;
        for (uint_t i : range(fdata.dims[1])) {
            if (iwei.safe[i] > 0) {
                ++npt;
                wei += iwei.safe[i];
                errb += sqr(iwei.safe[i]*(iflx.safe[i] - flx));
            }
        }

        errb = sqrt(errb*(npt/(npt-1.0)))/wei;
    };

    auto smooth_spectrum = [&](const vec2d& fcube, const vec2d& ecube, const vec2d& wcube,
        uint_t bwidth, vec1d& flx, vec1d& err, vec1d& errb) {

        if (bwidth != npos) {
            uint_t d1 = (bwidth-1)/2;
            uint_t d2 = bwidth-1 - d1;

            for (uint_t l : range(flx)) {
                uint_t l0 = (l > d1 ?            l - d1 : 0);
                uint_t l1 = (l < flx.size()-d2 ? l + d2 : flx.size()-1);

                stack_elements(fcube, ecube, wcube, l0, l1, flx[l], err[l], errb[l]);
            }
        } else {
            // Broadband fluxes, use a trick: only set one element of the spectrum to
            // the flux value, and ignore the rest since we do not want to "smooth" really
            for (uint_t l : range(fbb)) {
                if (fbb[l].res.empty()) {
                    flx[l] = dnan;
                    err[l] = dnan;
                    errb[l] = dnan;
                    continue;
                }

                vec2d tw = wcube*transpose(replicate(fbb[l].res, wcube.dims[1]));
                stack_elements(fcube, ecube, tw, 0, wcube.dims[0]-1, flx[l], err[l], errb[l]);
            }
        }
    };

    // Define function to write a 1D spectrum to the disk
    auto write_spectrum = [&](uint_t bwidth, std::string filename,
        const vec1d& l, const vec1d& flx, const vec1d& err, const vec1d& errb) {

        fits::output_image oimg(filename);
        oimg.write_empty();

        auto write_wcs = [&]() {
            oimg.write_keyword("BINNING", bwidth);
            oimg.write_keyword("CTYPE1", "WAVE");
            oimg.write_keyword("CUNIT1", "um");
            oimg.write_keyword("CRPIX1", 1.0);
            oimg.write_keyword("CRVAL1", l[0]);
            oimg.write_keyword("CDELT1", l[1]-l[0]);

            write_exposures(oimg);
        };

        oimg.reach_hdu(1);
        oimg.write(flx);
        write_wcs();

        if (!errb.empty()) {
            vec1f best_err = cure_error(errb, err);

            oimg.reach_hdu(2);
            oimg.write(best_err);
            write_wcs();

            oimg.reach_hdu(3);
            oimg.write(errb);
            write_wcs();
        } else {
            oimg.reach_hdu(2);
            oimg.write(err);
            write_wcs();

            oimg.reach_hdu(3);
            oimg.write_empty();
        }

        oimg.reach_hdu(4);
        oimg.write(err);
        write_wcs();

        if (verbose) note("wrote ", filename);
    };

    // Binning, stack, bootstrap and save files
    for (uint_t bwidth : rebin) {
        // Smooth with a boxcar of that width and stack
        vec1d flx(cflx.dims[0]);
        vec1d err(cflx.dims[0]);
        vec1d errb(cflx.dims[0]);

        smooth_spectrum(cflx, cerr, cwei, bwidth, flx, err, errb);

        std::string ofile = out;
        if (bwidth != 1) {
            ofile = file::remove_extension(out)+"_s"+to_string(bwidth)+".fits";
        }

        if (bwidth != npos) {
            write_spectrum(bwidth, ofile, lam, flx, err, errb);

            if (bwidth != 1) {
                // Then bin by picking regularly spaced elements in the smoothed spectrum
                vec1u idb = bwidth*uindgen(uint_t(ceil(flx.size()/float(bwidth)))) + (bwidth-1)/2;
                idb = idb[where(idb < flx.size())];
                flx = flx[idb];
                err = err[idb];
                errb = errb[idb];
                vec1d tlam = lam[idb];

                ofile = file::remove_extension(out)+"_b"+to_string(bwidth)+".fits";
                write_spectrum(bwidth, ofile, tlam, flx, err, errb);
            }
        } else {
            // Only keep the two elements we care about
            vec1u idf = uindgen(fbb.size());
            flx = flx[idf];
            err = err[idf];
            errb = errb[idf];
            vec1d tlam;
            for (auto& f : fbb) {
                tlam.push_back(f.rlam);
            }

            ofile = file::remove_extension(out)+"_broadband.fits";
            fits::write_table(ofile,
                "lambda", tlam, "flux", flx, "flux_err", err, "flux_err_bstrap", errb,
                "bands", filters
            );
        }
    }

    return 0;
}

void print_help() {
    using namespace terminal_format;

    print("stack1d v1.0");
    print("usage: stack1d files=[file1.fits,file2.fits,...] out=... [options]");
    print("");
    print("Main parameters:");
    paragraph("The files listed in 'files=[...]' must be valid 1D spectra FITS files, for example "
        "created by 'extract2d'. It is expected that the fluxes are located in the first HDU, and "
        "uncertainties in the second HDU. These spectra will be combined into a stacked spectrum, "
        "with accurate uncertainties determined from the variance between exposures. The 'out=...' "
        "parameter specifies the file name of the output stacked spectrum.");

    print("");
    print("Output format:");
    paragraph("Each output spectrum will contain the stacked flux in the first FITS extension, "
        "and the best estimate of the uncertainty in the second extension. The third extension "
        "will contain the bootstrap uncertainty estimate, while the fourth extension will "
        "contain the formal uncertainty estimate (from standard error propagation).");
    paragraph("If the 'filters=[...]' option is set (see below), the program will also create a "
        "file called '<out>_broadband.fits', which contains the stacked broadband fluxes. This "
        "file is a FITS binary table (column oriented). It contains the following columns: 'lambda'"
        "is the central wavelength of the filter, 'flux' is the flux, 'flux_err' is the best "
        "estimate of the uncertainty, 'flux_err_bstrap' is the bootstrap uncertainty, and 'bands' "
        "is the code name of the filter.");

    print("");
    print("Options:");
    bullet("files=[...]", "Must be a list of files. This is the list of 1D spectra that will be "
        "stacked.");
    bullet("out=...", "Must be a string. This sets the name of the file in which to store the "
        "stacked spectrum. It can include directories, which will be created if they do not exist. "
        "All other output files will use this file name as a base. For example, setting "
        "'out=somedir/spectrum.fits' will save all files in the 'somedir' directory, and their "
        "names will all start with 'spectrum...'.");
    bullet("rescale=...", "Must be a number. This is a global flux rescaling factor that is "
        "applied to both fluxes and uncertainties, to account for global transmission issues. "
        "Default is 1, meaning no rescaling is performed.");
    bullet("rebin=[...]", "Must be a list of integers. This defines different values of binning "
        "to apply to the stacked spectrum. Because the uncertainty in binned flux may not scale "
        "like pure Gaussian noise, the spectra are binned before being stacked, and the "
        "uncertainties for the binned spectra are computed from the variance in binned fluxes. "
        "For each value 'X' in this list, the program will create two files: <out>_bX.fits, "
        "which contains the stacked binned spectrum, and <out>_sX.fits, which contains the "
        "spectrum smoothed with a boxcar window of width 'X'. A value of one will produce an "
        "unbinned spectrum (simply saved as <out>.fits). Default is no rebinning.");
    bullet("filters=[...]", "Must be a list of strings. Each value must be the name of a filter "
        "form the filter database (see 'filter_db' below). For each filter, the program will "
        "compute the synthetic broadband flux, with the corresponding uncertainty. If omitted, "
        "no broadband flux will be computed.");
    bullet("filter_db=...", "Must be the path to a 'filter.db' file, which contains a list of "
        "filters. Each filter must be listed on a separate line, with the format '<name>=<path>', "
        "where 'name' is the code name of the filter (e.g., 'subaru-B' for the Subaru B band), "
        "and where 'path' is the path to the filter response curve. These response curves must "
        "be either FITS tables (columns: 'LAM' for the wavelength in microns, and 'RES' for the "
        "response, normalized to unit integral), or ASCII tables with two columns (wavelength "
        "and response, with the same units as for the FITS file).");
    bullet("stack_weight=...", "Must be either 'optimal' or 'uniform'. This defines the weighting "
        "used when stacking the 1D spectra. Using 'optimal' weighting will weight spectral data "
        "by inverse variance, which provides the best S/N. Using 'uniform' weighting will weight "
        "all data equally. Default is 'optimal'.");
    bullet("stack_op=...", "Must be either 'mean' or 'median'. This defines the method used for "
        "stacking the 1D spectra. Default is 'mean'. Both methods make use of the weights "
        "specified in 'stack_weight'.");
    bullet("do_sigma_clip", "Set this flag to enable sigma-clipping when stacking the data. "
        "If enabled, for each spectral element of the final spectrum, the data points that would "
        "enter the stack are compared to one another, and data points with strong deviations are "
        "excluded from the stack. Enabled by default.");
    bullet("sigma_clip_threshold", "Must be a number. When 'do_sigma_clip' is enabled, this sets "
        "the significance of the deviation required for a spectral element to be excluded. Default "
        "is 5 sigma.");
    bullet("sigma_clip_width", "Must be an integer. When 'do_sigma_clip' is enabled, this sets "
        "the width of the spectral window over which the background level and noise amplitude "
        "are computed to determine the pixel deviations. Default is 1 element.");
}
