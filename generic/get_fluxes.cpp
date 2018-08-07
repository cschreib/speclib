#include <phypp.hpp>

void print_help();

int phypp_main(int argc, char* argv[]) {
    vec1s  files;                 // 1D spectra to analyze
    bool   verbose = false;       // print progress in the standard output
    vec1s  filters;               // filters to use to build broadband fluxes
    std::string filter_db = data_dir+"fits/filter-db/db.dat";

    bool do_sigma_clip = true;          // enable/disable sigma clipping of outliers
    double sigma_clip_threshold = 5.0;  // significance threshold for rejecting outliers
    uint_t sigma_clip_width = 1;        // width (in pixels) of the wavelength bin in which to define outliers

    // Flux catalog used to rescale exposures to the right total flux
    std::string catalog_file;
    std::string catalog_id;

    bool help = false;

    read_args(argc, argv, arg_list(
        files, filters, name(catalog_file, "catalog"), catalog_id, filter_db,
        do_sigma_clip, sigma_clip_threshold, sigma_clip_width, verbose, help
    ));

    if (help) {
        print_help();
        return 0;
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
    }

    vec1d blam(fbb.size());
    for (uint_t f : range(fbb)) {
        blam[f] = fbb[f].rlam;
    }

    vec2d cflx, cerr;

    // Read spectra
    vec1d lam;
    for (uint_t i : range(files)) {
        vec1d tflx, terr;

        fits::input_image iimg(files[i]);
        iimg.reach_hdu(1);
        iimg.read(tflx);

        if (cflx.empty()) {
            double crpix, cdelt, crval;
            if (iimg.read_keyword("CRPIX1", crpix) && iimg.read_keyword("CRVAL1", crval) &&
                iimg.read_keyword("CDELT1", cdelt)) {
                lam = cdelt*(findgen(tflx.size())+1 - crpix) + crval;
            } else {
                std::string ctype;
                if (!iimg.read_keyword("CTYPE1", ctype)) {
                    error("could not read WCS of wavelength axis");
                    return 1;
                }

                if (begins_with(ctype, "TAB")) {
                    // Tabulated axis, read from other extensions
                    uint_t lowext = npos, upext = npos;
                    std::string axis = erase_begin(ctype, "TAB");
                    if (iimg.read_keyword(axis+"LOWEXT", lowext) &&
                        iimg.read_keyword(axis+"UPEXT", upext)) {
                        vec1d xl, xu;
                        iimg.reach_hdu(lowext);
                        iimg.read(xl);
                        iimg.reach_hdu(upext);
                        iimg.read(xu);
                        lam = 0.5*(xl + xu);
                        iimg.reach_hdu(1);
                    }
                } else {
                    error("could not read WCS of wavelength axis");
                    return 1;
                }
            }

            std::string cunit;
            if (iimg.read_keyword("CUNIT1", cunit)) {
                cunit = to_lower(cunit);
                double conv = 1.0;
                if (cunit == "angstrom") {
                    conv = 1e-4;
                } else if (cunit == "nm") {
                    conv = 1e-3;
                } else if (cunit == "um" || cunit == "micron") {
                    conv = 1.0;
                } else if (cunit == "mm") {
                    conv = 1e3;
                } else if (cunit == "cm") {
                    conv = 1e4;
                } else if (cunit == "m") {
                    conv = 1e6;
                } else {
                    error("unrecognized wavelength/frequency unit '", cunit, "'");
                    return false;
                }

                lam *= conv;
            }
        }

        iimg.reach_hdu(2);
        iimg.read(terr);

        if (cflx.empty()) {
            cflx.resize(tflx.dims, files.size());
            cerr.resize(cflx.dims);
        }

        cflx(_,i) = tflx;
        cerr(_,i) = terr;
    }

    // Define weights
    vec2d cwei = 1/sqr(cerr);
    // Adjust global normalization to avoid numerical errors
    cwei /= max(cwei[where(is_finite(cwei))]);

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

    // Down-weight bad pixels
    {
        vec1u idb = where(!is_finite(cflx) || !is_finite(cerr) || !is_finite(cwei) || crej);
        cwei[idb] = 0; cflx[idb] = 0; cerr[idb] = 0;
    }

    // Resample the filters to the grid of the spectra
    for (auto& f : fbb) {
        f.res = interpolate(f.res, f.lam, lam);
        uint_t i0 = where_first(lam > f.lam.front());
        uint_t i1 = where_last(lam < f.lam.back());
        if (i0 == npos || i1 == npos) {
            error("filter is not covered by spectrum (", mean(f.lam), " vs. ", min(lam), " - ", max(lam), ")");
            return 1;
        }

        f.res[_-i0] = 0;
        f.res[i1-_] = 0;
        f.lam = lam;
    }

    // Read catalog fluxes
    vec2f cflux, cflux_err;
    vec1s csid, cbands;
    uint_t cid = npos;

    if (!catalog_file.empty()) {
        fits::input_table itbl(catalog_file);
        itbl.read_columns("flux", cflux, "flux_err", cflux_err, "bands", cbands);
        if (!itbl.read_column("id", csid)) {
            vec1u cuid;
            if (itbl.read_column("id", cuid)) {
                csid = to_string_vector(cuid);
            } else {
                error("could not find column 'ID' in ", catalog_file);
                return 1;
            }
        }

        cid = where_first(csid == catalog_id);
        if (cid == npos) {
            error("could not find source with ID=", catalog_id, " in ", catalog_file);
            return 1;
        }
    }

    // Compute broad band fluxes for each exposure
    for (uint_t i : range(cflx.dims[1])) {
        vec1f bflx(filters.size());
        vec1f berr(filters.size());
        vec1f cat_flx = replicate(fnan, filters.size());
        vec1f cat_err = replicate(fnan, filters.size());

        for (uint_t b : range(filters)) {
            vec1f tw = cwei(_,i)*fbb[b].res;
            double w = total(tw);

            bflx[b] = total(tw*cflx(_,i))/w;
            berr[b] = sqrt(total(sqr(tw*cerr(_,i))))/w;

            if (!catalog_file.empty()) {
                uint_t ib = where_first(cbands == filters[b]);
                if (ib != npos) {
                    cat_flx[b] = uJy2cgs(fbb[b].rlam, cflux(cid,ib));
                    cat_err[b] = uJy2cgs(fbb[b].rlam, cflux_err(cid,ib));
                }
            }
        }

        if (!catalog_file.empty()) {
            fits::write_table(file::remove_extension(files[i])+"_broadband.fits",
                "flux", bflx, "flux_err", berr, "bands", filters, "lambda", blam,
                "catalog_flux", cat_flx, "catalog_flux_err", cat_err
            );
        } else {
            fits::write_table(file::remove_extension(files[i])+"_broadband.fits",
                "flux", bflx, "flux_err", berr, "bands", filters, "lambda", blam
            );
        }
    }

    return 0;
}

void print_help() {
    using namespace terminal_format;

    print("get_fluxes v1.0");
    print("usage: get_fluxes files=[file1.fits,file2.fits,...] filters=[...] catalog=... catalog_id=... [options]");
    print("");
    print("Main parameters:");
    paragraph("The files listed in 'files=[...]' must be valid 1D spectra FITS files, for example "
        "created by 'extract2d'. It is expected that the fluxes are located in the first HDU, and "
        "uncertainties in the second HDU. Synthetic broadband fluxes will be produced for each "
        "spectrum in this list. The spectra in the list are expected to eventually be stacked into "
        "a single 1D spectrum (for example with 'stack1d'). The reason why the entire list of 1D "
        "spectra must be provided here is that this program will use the other spectra to perform a "
        "sigma clipping, to ensure the synthetic broadband fluxes are not affected by bad pixels or "
        "other artifacts.");
    paragraph("The program will compute broad band fluxes for all filters listed in 'filters=[...]'. "
        "These must be valid entries in the filter data base (see 'filter_db' below), and there "
        "must exist an observed flux for this object in each of these filters in the provided "
        "flux catalog. The path to this flux catalog must be given in 'catalog=...'. The catalog "
        "must be a FITS table with columns 'ID' (to identify the galaxy, set 'catalog_id=...' to "
        "the corresponding ID in the catalog), 'FLUX' (2D column with dimensions [BAND,GALAXY]), "
        ", 'FLUX_ERR' (corresponding flux uncertainties), and 'BANDS' (list of filters with names "
        "drawn from the filter data base). Fluxes are expected to be given in microJansky.");

    print("");
    print("Output format:");
    paragraph("The program will create a file called <file>_broadband.fits for each spectrum in "
        "the list 'files=[...]'. This file is a FITS binary table (column oriented). It contains "
        "the following columns: 'lambda' is the central wavelength of the filter, 'flux' is the "
        "synthetic flux, 'flux_err' is the best estimate of the uncertainty, 'catalog_flux' is "
        "the corresponding flux in the catalog, 'catalog_flux_err' is the associated uncertainty, "
        "and 'bands' is the code name of the filter.");

    print("");
    print("Options:");
    bullet("files=[...]", "Must be a list of files. This is the list of 1D spectra that will be "
        "analyzed.");
    bullet("filters=[...]", "Must be a list of strings. Each value must be the name of a filter "
        "form the filter database (see 'filter_db' below). For each filter, the program will "
        "compute the synthetic broadband flux, with the corresponding uncertainty, and extract "
        "the associated observed broadband flux from the catalog.");
    bullet("filter_db=...", "Must be the path to a 'filter.db' file, which contains a list of "
        "filters. Each filter must be listed on a separate line, with the format '<name>=<path>', "
        "where 'name' is the code name of the filter (e.g., 'subaru-B' for the Subaru B band), "
        "and where 'path' is the path to the filter response curve. These response curves must "
        "be either FITS tables (columns: 'LAM' for the wavelength in microns, and 'RES' for the "
        "response, normalized to unit integral), or ASCII tables with two columns (wavelength "
        "and response, with the same units as for the FITS file).");
    bullet("catalog_file=...", "Must be the path to a FITS table containing the fluxes of this "
        "galaxy (among others) as observed in broadband filters from other surveys. See above "
        "for the expected format. If omitted, the program will not use the catalog at all and "
        "only compute the synthetic photometry.");
    bullet("catalog_id=...", "Must be an integer. It must correspond to the value of the column "
        "'ID' in the 'catalog_file' corresponding to the galaxy being reduced.");
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
