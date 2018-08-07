#include <phypp.hpp>

void print_help();

int phypp_main(int argc, char* argv[]) {
    std::string flux_file;
    std::string error_file;
    std::string telluric_gaps;
    std::string out;
    uint_t flux_hdu = 1;
    uint_t error_hdu = 2;
    bool help = false;

    read_args(argc, argv, arg_list(
        flux_file, error_file, flux_hdu, error_hdu, telluric_gaps, out, help
    ));

    if (help) {
        print_help();
        return 0;
    }

    if (error_file.empty()) {
        error_file = flux_file;
    }

    if (flux_file.empty()) {
        error("please provide the name of the telluric star spectrum file in flux_file=...");
        return 1;
    }

    if (out.empty()) {
        error("please provide the name of the output file in out=...");
        return 1;
    }

    file::mkdir(file::get_directory(out));

    // Build synthetic "photometry"
    fits::header hdr;
    vec1f lambda;
    vec1f flux, flux_err;
    vec1d tflx, terr, tlam; {
        fits::input_image fimg(flux_file);
        fimg.reach_hdu(flux_hdu);
        fimg.read(tflx);

        hdr = fimg.read_header();
        fits::setkey(hdr, "CTYPE1", "'WAVE'");
        tlam = astro::build_axis(astro::wcs(hdr), 0, astro::axis_unit::wave_um);

        note("lambda: ", min(tlam), " to ", max(tlam));

        fits::input_image eimg(error_file);
        eimg.reach_hdu(error_hdu);
        eimg.read(terr);
    } {
        vec1f lam_low, lam_up;
        ascii::read_table(telluric_gaps, lam_low, lam_up);

        vec1u idg = where(lam_low > min(tlam) && lam_up < max(tlam));
        lam_low = lam_low[idg];
        lam_up = lam_up[idg];

        lambda.resize(idg.size());
        flux.resize(idg.size());
        flux_err.resize(idg.size());
        for (uint_t l : range(lam_low)) {
            vec1u idl = where(tlam >= lam_low[l] && tlam <= lam_up[l]);

            vec1d f = tflx[idl];
            vec1d e = terr[idl];
            vec1d w = 1.0/sqr(terr[idl]);

            vec1u idb = where(!is_finite(f) || !is_finite(e) || !is_finite(w));
            f[idb] = 0; e[idb] = 0; w[idb] = 0;

            w /= max(w);

            flux[l] = total(f*w)/total(w);
            flux_err[l] = sqrt(total(sqr(e*w)))/total(w);
            lambda[l] = total(tlam[idl]*w)/total(w);
        }
    }

    if (flux.empty()) {
        error("no usable data in input spectrum");
        return 1;
    }

    vec1d baseline = interpolate(flux, lambda, tlam);

    fits::write_table(file::remove_extension(out)+"_intermediate.fits",
        "lam", tlam, "flx", tflx, "err", terr, "baseline", baseline, "blam", lambda, "bflx", flux, "berr", flux_err
    );

    vec1d cor = baseline/tflx;
    cor = max(cor, 1.0);

    fits::output_image oimg(out);
    oimg.reach_hdu(1);
    oimg.write(cor);
    oimg.write_header(hdr);

    return 0;
}

void print_help() {
    using namespace terminal_format;

    print("get_telluric v1.0 (xshooter)");
    print("usage: get_telluric flux_file=... telluric_gaps=... out=... [options]");
    print("");
    print("Main parameters:");
    paragraph("The program will analyze the spectrum of a telluric standard star given in "
        "'flux_file'. It will use the regions devoid of telluric absorption defined in "
        "'telluric_gaps' to determine the telluric correction, which will be saved in 'out'. No "
        "knowledge of the type of the star is needed.");

    print("");
    print("Output format:");
    paragraph("The program will create a file at the path specified in 'out', which will be a "
        "FITS 1D spectrum containing a single extension. This extension will hold the telluric "
        "correction to apply as a function of wavelength. This file can be used directly by "
        "'extract2d' to perform the telluric correction ('telluric=...' options).");

    print("");
    print("Options:");
    bullet("flux_file=...", "Must be a file. This is the 1D spectrum of the telluric standard star, "
        "as produced by the X-SHOOTER pipeline (TELL_SLIT_FLUX_MERGE1D_xxx.fits). This file must "
        "contain the observed flux of the star.");
    bullet("error_file=...", "Must be a file. It must contain the uncertainties on the flux of the "
        "standard star. The default is to be the same file as 'flux_file'.");
    bullet("flux_hdu=...", "Must be an integer. It sets the FITS HDU from which to read the flux "
        "of the standard star inside the file 'flux_file'. Default is 0 (primary HDU).");
    bullet("error_hdu=...", "Must be an integer. It sets the FITS HDU from which to read the flux "
        "uncertainties of the standard star inside the file 'error_file'. Default is 1 (first "
        "extension).");
    bullet("telluric_gaps=...", "Must be the path to a file listing the gaps in telluric "
        "absorption. Such a file is provided with this program ('telluric_regions.dat'). The file "
        "must be an ASCII table containing two floating point numbers on each line. Each pair of "
        "numbers defines a narrow wavelength window which is free of atmospheric absorption. The "
        "average flux in each window will be calculated, and used to interpolate the intrinsic "
        "spectrum of the standard star.");
    bullet("out=...", "Specifies the path of the output file containing the telluric correction.");
}
