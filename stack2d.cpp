#include <phypp.hpp>

vec2d cure_error(vec2d err, vec2d err_formal) {
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

    vec1s files;
    vec1s specs;
    uint_t errext = 1;
    double sigma_clip_threshold = 5.0;
    bool do_sigma_clip = true;
    bool verbose = false;
    uint_t rebin = 0;

    // Weighting scheme for stacking spectra and spectral binning
    // "optimal": inverse variance weighting, "uniform": no weighting
    std::string stack_weight = "optimal";

    bool help = false;

    read_args(argc, argv, arg_list(out, files, specs, errext, stack_weight, do_sigma_clip,
        sigma_clip_threshold, verbose, rebin, help));

    if (help) {
        print_help();
        return 0;
    }

    if (stack_weight != "optimal" && stack_weight != "uniform") {
        error("unknown weighting scheme '", stack_weight, "' for stacking");
        return 1;
    }

    if (files.empty()) {
        error("no spectrum to stack (files=...)");
        return 1;
    }

    if (!specs.empty() && specs.size() != files.size()) {
        error("mismatch in specs=... (", specs.size(), " files) and files=... (", files.size(), " files)");
        return 1;
    }

    vec1d rms_renorm = replicate(1.0, files.size());
    if (!specs.empty()) {
        for (uint_t i : range(specs)) {
            fits::input_image iimg(specs[i]);
            iimg.reach_hdu(errext);

            vec1d err;
            iimg.read(err);

            vec1u idg = where(is_finite(err) && err > 0);
            rms_renorm[i] = median(err[idg]);
        }

        rms_renorm /= min(rms_renorm);

        if (verbose) {
            for (uint_t i : range(files)) {
                print(files[i], ": ", rms_renorm[i]);
            }
        }
    }

    vec3d cflx, cerr;

    // Read spectra
    vec1d lam;
    fits::header hdr;
    for (uint_t i : range(files)) {
        vec2d tflx, terr;

        fits::input_image iimg(files[i]);
        iimg.reach_hdu(1);
        iimg.read(tflx);

        if (cflx.empty()) {
            hdr = iimg.read_header();
            lam = astro::build_axis(astro::wcs(hdr), 1, astro::axis_unit::wave_um);
        }

        iimg.reach_hdu(2);
        iimg.read(terr);

        if (cflx.empty()) {
            cflx.resize(tflx.dims, files.size());
            cerr.resize(cflx.dims);
        }

        cflx(_,_,i) = tflx;
        cerr(_,_,i) = terr*rms_renorm[i];
    }

    // Define weights
    vec3d cwei;
    if (stack_weight == "optimal") {
        cwei = 1/sqr(cerr);
    } else if (stack_weight == "uniform") {
        cwei = replicate(1.0, cflx.dims);
    }

    // Apply sigma clipping (if asked)
    if (do_sigma_clip) {
        for (uint_t p : range(cflx.dims[0]))
        for (uint_t l : range(cflx.dims[1])) {
            cflx(p,l,where(!sigma_clip(cflx(p,l,_), sigma_clip_threshold))) = dnan;
        }
    }

    // Down-weight bad pixels
    {
        vec1u idb = where(!is_finite(cflx) || !is_finite(cerr) || !is_finite(cwei));
        cwei[idb] = 0; cflx[idb] = 0; cerr[idb] = 0;
    }

    // Stack spectra
    vec2d flx(cflx.dims[0], cflx.dims[1]);
    vec2d err(cflx.dims[0], cflx.dims[1]);
    vec2d errb(cflx.dims[0], cflx.dims[1]);
    for (uint_t p : range(cflx.dims[0]))
    for (uint_t l : range(cflx.dims[1])) {
        double wei = total(cwei.safe(p,l,_));
        flx.safe(p,l) = total(cflx.safe(p,l,_)*cwei.safe(p,l,_))/wei;
        err.safe(p,l) = sqrt(total(sqr(cerr.safe(p,l,_)*cwei.safe(p,l,_))))/wei;
        errb.safe(p,l) = sqrt(cflx.dims[2]/(cflx.dims[2]-1.0)*total(sqr((cflx.safe(p,l,_) - flx.safe(p,l))*cwei.safe(p,l,_))))/wei;
    }

    vec2d best_err = cure_error(errb, err);

    // Save stacked spectrum
    {
        fits::output_image oimg(out);
        oimg.write_empty();
        oimg.reach_hdu(1);
        oimg.write(flx);
        oimg.write_header(hdr);
        oimg.reach_hdu(2);
        oimg.write(best_err);
        oimg.write_header(hdr);
        oimg.reach_hdu(3);
        oimg.write(errb);
        oimg.write_header(hdr);
        oimg.reach_hdu(4);
        oimg.write(err);
        oimg.write_header(hdr);

        if (verbose) note("wrote ", out);
    }

    if (rebin > 1) {
        // Stack spectra
        uint_t nlam = floor(cflx.dims[1]/float(rebin));

        uint_t d1 = (rebin-1)/2;
        uint_t d2 = rebin-1 - d1;

        flx.resize(cflx.dims[0], nlam);
        err.resize(cflx.dims[0], nlam);
        errb.resize(cflx.dims[0], nlam);
        for (uint_t p : range(cflx.dims[0]))
        for (uint_t ll : range(nlam)) {
            uint_t l = ll*rebin + (rebin-1)/2;
            uint_t l0 = (l > d1 ?              l - d1 : 0);
            uint_t l1 = (l < cflx.dims[1]-d2 ? l + d2 : cflx.dims[1]-1);

            double wei = total(cwei.safe(p,l0-_-l1,_));
            flx.safe(p,ll) = total(cflx.safe(p,l0-_-l1,_)*cwei.safe(p,l0-_-l1,_))/wei;
            err.safe(p,ll) = sqrt(total(sqr(cerr.safe(p,l0-_-l1,_)*cwei.safe(p,l0-_-l1,_))))/wei;
            errb.safe(p,ll) = sqrt(cflx.dims[2]/(cflx.dims[2]-1.0)*total(sqr((cflx.safe(p,l0-_-l1,_) - flx.safe(p,ll))*cwei.safe(p,l0-_-l1,_))))/wei;
        }

        best_err = cure_error(errb, err);

        // Update header
        double dl1 = lam[1]-lam[0];
        double dl2 = dl1*rebin;
        fits::setkey(hdr, "CTYPE1", "WAVE");
        fits::setkey(hdr, "CUNIT1", "um");
        fits::setkey(hdr, "CRPIX1", 1.0);
        // fits::setkey(hdr, "CRVAL1", lam[0]+0.5*(dl2 - dl1));
        fits::setkey(hdr, "CRVAL1", lam[(rebin-1)/2]);
        fits::setkey(hdr, "CDELT1", dl2);

        // Save stacked spectrum
        {
            out = file::remove_extension(out)+"_b"+to_string(rebin)+".fits";
            fits::output_image oimg(out);
            oimg.write_empty();
            oimg.reach_hdu(1);
            oimg.write(flx);
            oimg.write_header(hdr);
            oimg.reach_hdu(2);
            oimg.write(best_err);
            oimg.write_header(hdr);
            oimg.reach_hdu(3);
            oimg.write(errb);
            oimg.write_header(hdr);
            oimg.reach_hdu(4);
            oimg.write(err);
            oimg.write_header(hdr);

            if (verbose) note("wrote ", out);
        }
    }

    return 0;
}

void print_help() {
    using namespace terminal_format;

    print("stack2d v1.0");
    print("usage: stack2d files=[file1.fits,file2.fits,...] specs=[spec1.fits,spec2.fits,...] "
        "out=... [options]");
    print("");
    print("Main parameters:");
    paragraph("The files listed in 'files=[...]' must be valid 2D spectra FITS files, for example "
        "created by 'extract2d'. It is expected that the fluxes are located in the first HDU, and "
        "uncertainties in the second HDU. The 2D spectra will be combined into a stacked 2D "
        "spectrum, using the same weighting as would be used to stack the 1D spectra listed in "
        "specs=[...] (e.g., as if doing the stacking with 'stack1d'). Accurate uncertainties are "
        "determined from the variance between exposures. The 'out=...' parameter specifies the "
        "file name of the output stacked spectrum.");

    print("");
    print("Output format:");
    paragraph("Each output spectrum will contain the stacked flux in the first FITS extension, "
        "and the best estimate of the uncertainty in the second extension. The third extension "
        "will contain the bootstrap uncertainty estimate, while the fourth extension will "
        "contain the formal uncertainty estimate (from standard error propagation).");

    print("");
    print("Options:");
    bullet("files=[...]", "Must be a list of files. This is the list of 2D spectra that will be "
        "stacked.");
    bullet("specs=[...]", "Must be a list of files. This is the list of 1D spectra, extracted from "
        "the 2D spectra listed in 'files=[...]' (in the same order). These 1D spectra are used "
        "to define the weighting when stacking the 2D spectra, in order to build a stacked 2D "
        "spectrum that would be most representative for the stacked 1D spectrum.");
    bullet("out=...", "Must be a string. This sets the name of the file in which to store the "
        "stacked spectrum. It can include directories, which will be created if they do not exist. "
        "All other output files will use this file name as a base. For example, setting "
        "'out=somedir/spectrum.fits' will save all files in the 'somedir' directory, and their "
        "names will all start with 'spectrum...'.");
    bullet("rebin=...", "Must be an integer. This defines the binning to apply to the stacked "
        "spectrum. Because the uncertainty in binned flux may not scale like pure Gaussian noise, "
        "the spectra are binned before being stacked, and the uncertainties for the binned spectra "
        "are computed from the variance in binned fluxes. For a value 'X', the program will saved "
        "the binned spectrum in the file <out>_bX.fits (in addition to the unbinned spectrum). "
        "Default is 0, which does not produce a re-binned spectrum.");
    bullet("stack_weight=...", "Must be either 'optimal' or 'uniform'. This defines the weighting "
        "used when stacking the 1D spectra. Using 'optimal' weighting will weight spectral data "
        "by inverse variance, which provides the best S/N. Using 'uniform' weighting will weight "
        "all data equally. Default is 'optimal'.");
    bullet("do_sigma_clip", "Set this flag to enable sigma-clipping when stacking the data. "
        "If enabled, for each spectral element of the final spectrum, the data points that would "
        "enter the stack are compared to one another, and data points with strong deviations are "
        "excluded from the stack. Enabled by default.");
    bullet("sigma_clip_threshold", "Must be a number. When 'do_sigma_clip' is enabled, this sets "
        "the significance of the deviation required for a spectral element to be excluded. Default "
        "is 5 sigma.");
}
