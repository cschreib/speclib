#include <phypp.hpp>
#include <phypp/math/mpfit.hpp>

void print_help();

int phypp_main(int argc, char* argv[]) {
    std::string out_dir;           // output directory
    std::string name;              // name of the target in the spectrum
    std::string instrument;        // name of instrument
    vec1s sources;                 // names of the multiple sources in the spectrum
    uint_t dpix = 10;              // spatial width of stacked 2D
    uint_t dpixfit = 6;            // spatial width of the Gaussian fit area within the 2D spectrum
    uint_t crop_border = 0;        // number of pixels to remove from spectrum borders before stacking
    bool   no_ohflag = false;      // set to true to keep OH line residuals in the spectrum
    bool   no_gapflag = false;     // set to true to keep transmission gaps in the spectrum
    double flag_begin = fnan;      // set to flag the first wavelenth elements until this wavelenth
    double flag_end = fnan;        // set to flag the last wavelenth elements after this wavelenth
    float  oh_line_width = 100.0;  // velocity width of OH lines (km/s)
    std::string suffix;            // suffix to append to the source ID
    float  offset = 0.0;           // offset along the slit wrt the expected source position
    bool   discard_wcs_pos = false; // do not use WCS slit position, use half of slit width as center
    bool   unit_flux = false;      // set input 2d spectrum = 1 (to derive effective telluric correction)
    float  rms_rescale = 1.0;      // RMS rescaling factor for each epoch (applied error before stacking)
    bool   median_background = false; // subtract the median of off-source pixels at each wavelength
    bool   verbose = false;        // print intermediate step status in the terminal
    std::string telluric;          // path to a FITS table with telluric + flux correction
    bool   auto_mask = false;      // automatically mask bad pixels in the 2D spectrum
    uint_t auto_mask_dl = 12;      // width of the spectral region to use for auto masking (larger = may mask lines!)
    float  auto_mask_sigma = 4;    // significance threshold for masking pixel outliers
    std::string flux_file;         // matching pattern for finding flux frames
    std::string error_file;        // matching pattern for finding uncertainty frames
    std::string flag_file;         // matching pattern for finding flag frames
    uint_t flux_hdu = 1;           // HDU of data containing flux extension
    uint_t error_hdu = 2;          // HDU of data containing uncertainty extension
    uint_t flag_hdu = npos;        // HDU of data containing quality flag extension
    uint_t max_flag = 0;           // Maximum value allowed for the 'flag' of a pixel to consider it a good pixel

    float seeing = 0.0;                // estimated Gaussian sigma for point sources (in pixels)
    vec1f default_gauss_width = {1.4}; // sigma of the Gaussian model (if fit to the data fails)
    vec1f default_gauss_pos = {0.0};   // position of the Gaussian model (if fit to the data fails)
    bool gauss_fit_background = false; // fit the Gaussian model over a flat background
    bool gauss_fit_background_linear = false; // fit the Gaussian model over a linear background
    bool gauss_nofit = false;          // don't fit for the Gaussian model, just use the default values
    bool gauss_nofit_width = false;    // don't fit for the Gaussian model widths, just use the default values
    float gauss_min_snr = 6.0;         // minimum peak SNR of the collapsed profile to do the Gaussian fit

    double lmin = dnan;           // minimum wavelength of the stacked spectrum [um]
    double lmax = dnan;           // maximum wavelength of the stacked spectrum [um]

    bool help = false;

    read_args(argc, argv, arg_list(instrument, name, flux_file, error_file, flux_hdu, error_hdu,
        dpix, dpixfit, no_ohflag, no_gapflag, oh_line_width, offset, suffix, out_dir, verbose, unit_flux,
        sources, default_gauss_width, gauss_fit_background, gauss_fit_background_linear,
        gauss_nofit, gauss_nofit_width, gauss_min_snr, lmin, lmax, rms_rescale, telluric,
        auto_mask, auto_mask_dl, auto_mask_sigma, default_gauss_pos, median_background,
        crop_border, flux_hdu, error_hdu, discard_wcs_pos, flag_begin, flag_end,
        flag_file, flag_hdu, max_flag, seeing, help
    ));

    if (help) {
        print_help();
        return 0;
    }

    if (error_file.empty()) {
        error_file = flux_file;
    }
    if (flag_file.empty()) {
        flag_file = flux_file;
    }

    if (name.empty()) {
        error("please provide the name of the source in name=...");
    }

    name += suffix;

    if (instrument.empty()) {
        error("please provide the name of the instrument in instrument=...");
        return 1;
    }

    if (flux_file.empty()) {
        error("please provide the file name of the flux frame in flux_file=...");
        return 1;
    }

    bool full_coverage = false;
    if (!is_finite(lmin) || !is_finite(lmax)) {
        if (verbose) note("using full wavelength coverage of input spectrum");
        full_coverage = true;
    }

    if (dpixfit > dpix) {
        error("fit aperture for 1D spectrum (dpixfit=", dpixfit, ") must be equal or smaller than the width "
            "of the 2D spectrum (dpix=", dpix, ")");
        return 1;
    }

    if (sources.size() > 1 && default_gauss_width.size() == 1.0) {
        default_gauss_width = replicate(default_gauss_width[0], sources.size());
    }
    if (default_gauss_width.size() > 1 && default_gauss_pos.size() == 1.0) {
        default_gauss_pos = replicate(default_gauss_pos[0], default_gauss_width.size());
    }
    if (default_gauss_pos.size() > 1 && default_gauss_width.size() == 1.0) {
        default_gauss_width = replicate(default_gauss_width[0], default_gauss_pos.size());
    }
    if (default_gauss_pos.size() > 1 && sources.size() != default_gauss_pos.size()) {
        error("must provide names for each gaussian component in sources=...");
        return 1;
    }

    if (gauss_fit_background_linear) {
        gauss_fit_background = true;
    }

    out_dir = file::directorize(out_dir);
    file::mkdir(out_dir);

    // Strongest OH line database
    vec1d lines;

    if (instrument == "MOSFIRE") {
        lines = {
            // H band
            1.4564, 1.4606, 1.4666, 1.4699, 1.4784, 1.4832, 1.4864, 1.4888,
            1.4909, 1.4932, 1.4803, 1.5054, 1.5057, 1.5069, 1.5088, 1.5188,
            1.5241, 1.5287, 1.5332, 1.5395, 1.5433, 1.5505, 1.5515, 1.5541,
            1.5544, 1.5570, 1.5598, 1.5632, 1.5656, 1.5702, 1.5832, 1.5848,
            1.5869, 1.5973, 1.6031, 1.6080, 1.6129, 1.6195, 1.6236, 1.6317,
            1.6352, 1.6388, 1.6415, 1.6444, 1.6477, 1.6503, 1.6555, 1.6611,
            1.6690, 1.6693, 1.6706, 1.6709, 1.6733, 1.6841, 1.6904, 1.6955,
            1.7009, 1.7079, 1.7124, 1.7211, 1.7249, 1.7283, 1.7330, 1.7350,
            1.7358, 1.7386, 1.7428, 1.7450, 1.7505, 1.7651, 1.7669, 1.7696,
            1.7810, 1.7879, 1.7993, 1.8065, 1.8117, 1.8125, 1.8205,

            // K band
            1.9250, 1.9350, 1.9399, 1.9642, 1.9701, 1.9753, 1.9774, 1.9841,
            2.0010, 2.0034, 2.0278, 2.0342, 2.0414, 2.0500, 2.0566, 2.0731,
            2.0862, 2.0909, 2.1177, 2.1233, 2.1252, 2.1509, 2.1541, 2.1710,
            2.1804, 2.1874, 2.1957, 2.2125
        };
    } else if (instrument == "X-SHOOTER") {
        lines = {
            // VIS
            0.53408, 0.53430, 0.53456, 0.56242, 0.61481, 0.63003, 0.66230,
            0.76243, 0.76282, 0.79137, 0.80944, 0.83446, 0.83992, 0.84302,
            0.84654, 0.88271, 0.88859, 0.89197, 0.89581, 0.89881, 0.90011,
            0.93166, 0.93238, 0.93435, 0.93539, 0.93582, 0.93757, 0.94402,
            0.94768, 0.94818, 0.95013, 0.95024, 0.95194, 0.95531, 0.95671,
            0.96993, 0.97198, 0.97910, 0.98721, 0.99147, 1.00130,

            // J band
            1.03731, 1.04184, 1.04691, 1.08315, 1.09236, 1.09725, 1.10271,
            1.11630, 1.12173, 1.12218, 1.12248, 1.12344, 1.12527, 1.12758,
            1.12947, 1.13310, 1.13445, 1.13496, 1.13577, 1.14060, 1.14369,
            1.14405, 1.14636, 1.15239, 1.15356, 1.15548, 1.15887, 1.16478,
            1.20279, 1.21194, 1.22262, 1.22838, 1.23483, 1.24200, 1.28322,
            1.29024, 1.30182, 1.30491, 1.30818, 1.31247, 1.31535, 1.32330,
            1.33212, 1.34748, 1.34811, 1.35003,

            // Gap JH
            1.35144, 1.35279, 1.35405, 1.35468, 1.35540, 1.35561, 1.35672,
            1.35684, 1.35738, 1.35801, 1.35840, 1.35894, 1.35912, 1.35945,
            1.36041, 1.36092, 1.36110, 1.36152, 1.36170, 1.36209, 1.36260,
            1.36320, 1.36362, 1.36380, 1.36404, 1.36440, 1.36488, 1.36533,
            1.36644, 1.36692, 1.36794, 1.36812, 1.36854, 1.36953, 1.36983,
            1.37037, 1.37085, 1.37103, 1.37142, 1.37166, 1.37187, 1.37202,
            1.37262, 1.37307, 1.37457, 1.37556, 1.37574, 1.37667, 1.37838,
            1.37982, 1.38006, 1.38036, 1.38051, 1.38150, 1.38165, 1.38198,
            1.38213, 1.38234, 1.38270, 1.38312, 1.38354, 1.38450, 1.38477,
            1.38549, 1.38603, 1.38738, 1.38768, 1.38795, 1.38951, 1.39044,
            1.39170, 1.39212, 1.39230, 1.39251, 1.39353, 1.39425, 1.39440,
            1.39566, 1.39626, 1.39641, 1.39683, 1.39905, 1.39938, 1.39974,
            1.39992, 1.40016, 1.40040, 1.40076, 1.40115, 1.40136, 1.40271,
            1.40385, 1.40427, 1.40505, 1.40631, 1.40718, 1.40868, 1.40934,
            1.41192, 1.41321,

            // H band
            1.41402, 1.41450, 1.41498, 1.41836, 1.41859, 1.41879, 1.41909,
            1.42026, 1.42236, 1.42293, 1.42356, 1.43082, 1.43136, 1.43208,
            1.43397, 1.43498, 1.43546, 1.43546, 1.43691, 1.43736, 1.43916,
            1.44108, 1.44255, 1.44486, 1.44528, 1.44664, 1.45008, 1.45032,
            1.45062, 1.45155, 1.45192, 1.45494, 1.45599, 1.46004, 1.46280,
            1.46484, 1.46613, 1.46937, 1.46983, 1.47033, 1.47198, 1.47348,
            1.47531, 1.47675, 1.47792, 1.47966, 1.48005, 1.48023, 1.48293,
            1.48611, 1.48839, 1.49054, 1.49280, 1.49609, 1.50499, 1.50516,
            1.50613, 1.50651, 1.50846, 1.51080, 1.51833, 1.52370, 1.52838,
            1.53105, 1.53285, 1.53448, 1.53727, 1.53912, 1.54281, 1.54347,
            1.54632, 1.54597, 1.54965, 1.54846, 1.55055, 1.55136, 1.55361,
            1.55421, 1.55658, 1.55934, 1.56273, 1.56516, 1.56801, 1.56954,
            1.56984, 1.57072, 1.57152, 1.57206, 1.57291, 1.57341, 1.57426,
            1.57474, 1.57566, 1.57624, 1.57772, 1.57998, 1.58036, 1.58280,
            1.58292, 1.58398, 1.58439, 1.58654, 1.59684, 1.60266, 1.60734,
            1.60761, 1.61149, 1.61220, 1.61250, 1.61554, 1.62181, 1.61904,
            1.62312, 1.62417, 1.62481, 1.62636, 1.62983, 1.63033, 1.63128,
            1.63375, 1.63470, 1.63564, 1.63842, 1.64103, 1.64376, 1.64433,
            1.64473, 1.64725, 1.64745, 1.64979, 1.65255, 1.65495, 1.65820,
            1.66068, 1.66850, 1.66885, 1.66992, 1.67043, 1.67282, 1.67506,
            1.67597, 1.68363, 1.68996, 1.69037, 1.69503, 1.70043, 1.70518,
            1.70736, 1.71192, 1.71847, 1.72062, 1.72394, 1.72437, 1.72533,
            1.72785, 1.72992, 1.73265, 1.73468, 1.73559, 1.73790, 1.73825,
            1.74225, 1.74453, 1.74979, 1.75011, 1.75248, 1.76456, 1.76487,
            1.76548, 1.76606, 1.76672, 1.76804, 1.76844, 1.76940, 1.77271,
            1.77300, 1.78074, 1.78256, 1.78437, 1.78758, 1.79652, 1.79883,
            1.79907, 1.79994, 1.80063, 1.80552, 1.80603, 1.80657, 1.80696,
            1.81002, 1.81059, 1.81155, 1.81236, 1.81475, 1.81515, 1.81581,
            1.81660, 1.81812, 1.81956, 1.82070, 1.82133, 1.82208,

            // K band
            1.91883, 1.92102, 1.92324, 1.92342, 1.92423, 1.92678, 1.92813,
            1.92984, 1.93056, 1.93356, 1.93596, 1.93617, 1.93758, 1.93893,
            1.94040, 1.94100, 1.94316, 1.94424, 1.95273, 1.95882, 1.96743,
            1.96962, 1.97667, 1.97934, 1.98648, 2.00016, 2.00046, 2.00568,
            2.02707, 2.03346, 2.04072, 2.05572, 2.09694, 2.09718, 2.09754,
            2.09775, 2.09874, 2.09955, 2.09991, 2.10006, 2.10036, 2.10096,
            2.10114
        };
    } else {
        error("unknown instrument '", instrument, "'");
        return 1;
    }

    vec1f gaps_low;
    vec1f gaps_up;

    if (instrument == "X-SHOOTER") {
        gaps_low = {0.75979, 1.3500487, 1.8221314};
        gaps_up  = {0.76109, 1.4147945, 1.9186203};
    } else if (instrument == "MOSFIRE") {
        // Nothing to flag, single band spectra
    } else {
        error("unknown instrument '", instrument, "'");
        return 1;
    }

    vec1d tell_cor, tell_lam;
    if (!telluric.empty()) {
        fits::input_image iimg(telluric);
        iimg.reach_hdu(1);
        iimg.read(tell_cor);
        fits::header hdr = iimg.read_header();
        fits::setkey(hdr, "CTYPE1", "'WAVE'");
        tell_lam = astro::build_axis(astro::wcs(hdr), 0, astro::axis_unit::wave_um);
    } else {
        note("assuming telluric correction is already applied or not necessary (else use telluric=...)");
    }

    double exptime = 0;
    double ps0, aspix;
    vec1d lam;
    vec2d tflx, terr; {
        // Read flux
        fits::input_image iimg(flux_file);
        iimg.reach_hdu(flux_hdu);
        iimg.read(tflx);

        fits::header hdr = iimg.read_header();
        fits::setkey(hdr, "CTYPE1", "'WAVE'");
        astro::wcs w(hdr);

        // Read exposure time
        if (instrument == "X-SHOOTER") {
            uint_t nframe = 0;
            float dit = 0.0;
            if (fits::getkey(hdr, "ESO PRO DATANCOM", nframe) &&
                fits::getkey(hdr, "HIERARCH ESO DET DIT", dit)) {
                exptime = nframe*dit;
            }
        } else if (instrument == "MOSFIRE") {
            uint_t nframe = 0;
            float dit = 0.0;
            if (fits::getkey(hdr, "NUMFRM", nframe) && fits::getkey(hdr, "TRUITIME", dit)) {
                exptime = 2*nframe*dit;
            }
        }

        // Build wavelength axis
        lam = astro::build_axis(w, 1, astro::axis_unit::wave_um);
        if (verbose) note("lambda min: ", min(lam), ", max: ", max(lam));

        // Get vertical position of source
        double lam0, sp, sp1;
        xy2ad(w, tflx.dims[1]/2, tflx.dims[0]/2, lam0, sp);
        xy2ad(w, tflx.dims[1]/2, tflx.dims[0]/2+1, lam0, sp1);
        double pl0;
        ad2xy(w, lam0, 0.0, pl0, ps0);
        ps0 -= 1.0;
        aspix = sp1 - sp;

        if (instrument == "MOSFIRE") {
            aspix = 0.18;
        }

        if (discard_wcs_pos) {
            ps0 = tflx.dims[0]/2;
        }

        // Apply manual offset if needed
        ps0 += offset;

        if (verbose) {
            note("extracting source row ", ps0);
        }

        // Read uncertainty
        fits::input_image ierr(error_file);
        ierr.reach_hdu(error_hdu);
        ierr.read(terr);
    }

    double seeing_pix = seeing/(2.335*aspix); // convert seeing from FWHM in arcsec to sigma in pixels
    default_gauss_pos /= aspix;   // convert arcsec to pixels
    default_gauss_width /= aspix; // convert arcsec to pixels

    if (!full_coverage) {
        double cdelt = lam[1]-lam[0];

        if (lmin > lam.back() || lmax < lam[0]) {
            error("requested wavelength range is not covered by these data (asked ",
                lmin, " to ", lmax, ", got ", lam[0], " to ", lam.back(), ")");
            return 1;
        }

        // Apply min/max lam
        if (lmin < lam[0] - 0.5*cdelt) {
            uint_t nadd = ceil((lam[0] - 0.5*cdelt - lmin)/cdelt);
            prepend(lam, lam[0] + (dindgen(nadd) - nadd)*cdelt);
            prepend<0>(tflx, replicate(dnan, nadd, tflx.dims[1]));
            prepend<0>(terr, replicate(dnan, nadd, tflx.dims[1]));
        } else if (lmin > lam[0] + 0.5*cdelt) {
            uint_t p = lower_bound(lam, lmin-0.5*cdelt);
            lam = lam[p-_];
            tflx = tflx(_,p-_);
            terr = terr(_,p-_);
        }

        if (lmax > lam.back() + 0.5*cdelt) {
            uint_t nadd = ceil((lmax - (lam.back() + 0.5*cdelt))/cdelt);
            append(lam, lam.back() + (dindgen(nadd) + 1)*cdelt);
            append<0>(tflx, replicate(dnan, nadd, tflx.dims[1]));
            append<0>(terr, replicate(dnan, nadd, tflx.dims[1]));
        } else if (lmax < lam.back() - 0.5*cdelt) {
            uint_t p = upper_bound(lam, lmax+0.5*cdelt);
            lam = lam[_-p];
            tflx = tflx(_,_-p);
            terr = terr(_,_-p);
        }
    }

    uint_t tnpix = lam.dims[1];

    std::string band;
    if (mean(lam) > 1.8 && mean(lam) < 2.3) {
        band = "K";
    } else if (mean(lam) > 1.4 && mean(lam) < 1.8) {
        band = "H";
    } else if (mean(lam) > 1.1 && mean(lam) < 1.4) {
        band = "J";
    } else if (mean(lam) > 0.95 && mean(lam) < 1.05) {
        band = "Y";
    } else if (mean(lam) > 0.83 && mean(lam) < 0.95) {
        band = "Z";
    } else if (mean(lam) > 0.70 && mean(lam) < 0.83) {
        band = "I";
    } else if (mean(lam) > 0.58 && mean(lam) < 0.70) {
        band = "R";
    } else if (mean(lam) > 0.51 && mean(lam) < 0.58) {
        band = "V";
    } else if (mean(lam) > 0.40 && mean(lam) < 0.51) {
        band = "B";
    } else if (mean(lam) > 0.30 && mean(lam) < 0.4) {
        band = "U";
    }

    if (band == "H" && min(lam) < 1.2 && max(lam) > 2.1) {
        band = "NIR";
    } else if (band == "I" && min(lam) < 0.6 && max(lam) > 0.9) {
        band = "VIS";
    } else if (band == "B" && min(lam) < 0.35 && max(lam) > 0.52) {
        band = "UVB";
    }

    if (band.empty()) {
        warning("could not identify band: ", min(lam), " to ", max(lam), " um");
        band = "UNK";
    }

    // Subtract background
    if (median_background) {
        if (verbose) note("removing background from median in slice");

        vec1b mask = replicate(true, tflx.dims[0]);
        for (uint_t i : range(default_gauss_width)) {
            double d0 = ps0 + default_gauss_pos[i] - 1.5*default_gauss_width[i];
            double d1 = ps0 + default_gauss_pos[i] + 1.5*default_gauss_width[i];
            uint_t i0 = max(0, round(d0));
            uint_t i1 = min(mask.size()-1, round(d1));
            mask[i0-_-i1] = false;
        }

        vec1u idm = where(mask);
        for (uint_t l : range(tflx.dims[1])) {
            tflx.safe(_,l) -= median(tflx.safe(idm,l));
        }
    }

    // Mask bad pixels automatically (if asked)
    if ((auto_mask || flag_hdu != npos) && verbose) note("masking bad pixels");

    if (auto_mask) {
        uint_t auto_mask_ndl = auto_mask_dl/2;
        uint_t auto_mask_pdl = auto_mask_dl - auto_mask_ndl;

        vec2b mask(tflx.dims);
        for (uint_t l : range(tflx.dims[1])) {
            uint_t l0 = (l < auto_mask_ndl ? 0 : l - auto_mask_ndl);
            uint_t l1 = (l >= tflx.dims[1]-auto_mask_pdl ? tflx.dims[1]-1 : l + auto_mask_pdl);
            for (uint_t x : range(tflx.dims[0])) {
                double m = median(tflx.safe(x,l0-_-l1));
                double s = max(1.48*mad(tflx.safe(x,l0-_-l1)), 1.48*mad(tflx.safe(_,l)));
                mask.safe(x,l) = abs(tflx.safe(x,l) - m)/s > auto_mask_sigma;
            }
        }

        mask = mask_inflate(mask, 1);

        // Mask anomalously low errors
        for (uint_t l : range(terr.dims[1])) {
            double m = median(terr.safe(_,l));
            mask.safe(where(terr.safe(_,l)/m < 0.5),l) = true;
        }

        vec1u idm = where(mask);
        tflx[idm] = dnan;
        terr[idm] = dnan;
    }

    // Read and apply flag, if asked
    if (flag_hdu != npos) {
        vec2u tflag;

        fits::input_image iflag(flag_file);
        iflag.reach_hdu(flag_hdu);
        iflag.read(tflag);

        vec1u idb = where(tflag > max_flag);
        tflx[idb] = dnan;
        terr[idb] = dnan;
    }

    // Flag out strong OH residuals (if asked)
    vec1b flagged(lam.dims);
    if (!no_ohflag) {
        if (verbose) note("flagging OH residuals");
        for (uint_t l : range(lines)) {
            flagged = flagged || (lam >= lines[l]*(1.0 - oh_line_width/3e5) &&
                                  lam <= lines[l]*(1.0 + oh_line_width/3e5));
        }
    }

    if (!no_gapflag) {
        if (verbose) note("flagging transmission gaps");
        for (uint_t l : range(gaps_low)) {
            flagged = flagged || (lam >= gaps_low[l] && lam <= gaps_up[l]);
        }
    }

    if (is_finite(flag_begin)) {
        flagged = flagged || lam <= flag_begin;
    }

    if (is_finite(flag_end)) {
        flagged = flagged || lam >= flag_end;
    }

    // Flag out spectrum borders (if asked)
    if (crop_border != 0) {
        vec1d t1d = tflx(tflx.dims[0]/2,_);
        uint_t tlmin = where_first(is_finite(t1d));
        if (tlmin != npos) {
            tlmin = min(tnpix-1, tlmin+crop_border);
            flagged[_-tlmin] = true;
        }

        uint_t tlmax = where_last(is_finite(t1d));
        if (tlmax != npos) {
            tlmax -= min(crop_border, tlmax);
            flagged[tlmax-_] = true;
        }
    }

    for (uint_t l : where(flagged)) {
        tflx.safe(_,l) = dnan;
        terr.safe(_,l) = dnan;
    }

    // Truncate the 2D spectrum around the region of interest
    // NB: do it after "auto_mask" to make sure the masking uses as much statistics as possible
    int_t ips0 = round(ps0);
    if (ips0 < 0 || uint_t(ips0) > tflx.dims[0]-1) {
        error("source is outside of the 2D spectrum (Y coordinate: ", ips0, ")");
        return 1;
    }

    if (ips0 < int_t(dpix) || ips0 > int_t(tflx.dims[0])-1-int_t(dpix)) {
        vec2d oflx = tflx;
        vec2d oerr = terr;
        tflx = replicate(dnan, 2*dpix+1, oflx.dims[1]);
        terr = replicate(dnan, 2*dpix+1, oflx.dims[1]);

        uint_t i0 = max(ips0-int_t(dpix), 0);
        uint_t i1 = min(ips0+int_t(dpix), int_t(oflx.dims[0])-1);
        uint_t o0 = dpix + (int_t(i0) - ips0);
        uint_t o1 = dpix + (int_t(i1) - ips0);

        tflx(o0-_-o1,_) = oflx(i0-_-i1,_);
        terr(o0-_-o1,_) = oerr(i0-_-i1,_);
    } else {
        tflx = tflx((ips0-dpix)-_-(ips0+dpix),_);
        terr = terr((ips0-dpix)-_-(ips0+dpix),_);
    }

    // Apply rescaling factor (if asked)
    terr *= rms_rescale;

    // Make flux = 1 if asking for effective telluric correction
    if (unit_flux) {
        tflx[_] = 1.0;
    }

    // Apply telluric correction (if asked)
    if (!telluric.empty()) {
        if (verbose) note("apply telluric correction");

        vec1d tc = interpolate(tell_cor, tell_lam, lam);
        for (uint_t p : range(tflx.dims[0])) {
            tflx.safe(p,_) *= tc;
            terr.safe(p,_) *= tc;
        }
    }

    // Mask out the points with weird values
    {
        vec1u idb = where(!is_finite(tflx) || terr <= 0 || !is_finite(terr));
        tflx[idb] = fnan; terr[idb] = fnan;
    }

    if (verbose) note("write 2D spectrum");

    {
        fits::output_image oimg2d(out_dir+"stacked_"+name+"_spec2d.fits");

        auto write_wcs_2d = [&]() {
            oimg2d.write_keyword("CTYPE1", "WAVE");
            oimg2d.write_keyword("CUNIT1", "um");
            oimg2d.write_keyword("CRPIX1", 1.0);
            oimg2d.write_keyword("CRVAL1", lam[0]);
            oimg2d.write_keyword("CDELT1", lam[1]-lam[0]);

            oimg2d.write_keyword("CTYPE2", "SLIT");
            oimg2d.write_keyword("CUNIT2", "arcsec");
            oimg2d.write_keyword("CRPIX2", dpix+1);
            oimg2d.write_keyword("CRVAL2", 0.0);
            oimg2d.write_keyword("CDELT2", aspix);

            oimg2d.write_keyword("BAND", band);
            oimg2d.write_keyword("INSTR", instrument);
            oimg2d.write_keyword("EXPTIME", exptime);
            oimg2d.write_keyword("SEEING", seeing);
        };

        oimg2d.write_empty();
        oimg2d.reach_hdu(1);
        oimg2d.write(tflx);
        write_wcs_2d();
        oimg2d.reach_hdu(2);
        oimg2d.write(terr);
        write_wcs_2d();
    }

    // Extract 1D spectrum by fitting one or several gaussian profiles
    if (verbose) note("extract 1D spectrum");

    // Only fit a portion of the 2D spectrum
    vec2d tflxfit, terrfit;
    if (dpixfit < dpix) {
        tflxfit = tflx((dpix-dpixfit)-_-(dpix+dpixfit),_);
        terrfit = terr((dpix-dpixfit)-_-(dpix+dpixfit),_);
    } else {
        tflxfit = tflx;
        terrfit = terr;
    }

    // First build the collapsed profile
    vec1d sprof(tflxfit.dims[0]);
    vec1d serror(tflxfit.dims[0]);
    for (uint_t p : range(sprof)) {
        vec1d f = tflxfit(p,_);
        vec1d e = terrfit(p,_);
        vec1d w = 1/sqr(e);

        vec1u idb = where(!is_finite(w) || !is_finite(f) || !is_finite(e) || !sigma_clip(f/e, 10.0));
        w[idb] = 0; e[idb] = 0; f[idb] = 0;

        w /= median(w[where(w > 0)]);

        double tw = total(w);
        sprof[p] = total(f*w)/tw;
        serror[p] = sqrt(total(sqr(e*w)))/tw;
    }

    uint_t np = sprof.size();
    double psnr = sprof[np/2]/serror[np/2];
    double msnr = max(sprof/serror);
    bool badsnr = (psnr < gauss_min_snr || !is_finite(psnr));
    if (verbose) {
        note("peak profile S/N in ", band, " is ", psnr);
        note("max profile S/N in ", band, " is ", msnr);
    }

    vec1d sx = findgen(np) - np/2;
    serror /= max(sprof);
    sprof /= max(sprof);

    vec1u idbp = where(!is_finite(sprof) || serror <= 0.0);
    serror[idbp] = 1e9;
    sprof[idbp] = 0;

    // Increase weights towards the center
    serror *= sqr(sx/(sprof.size()/2.0)) + 1.0;

    uint_t nsrc = sources.size();
    if (nsrc == 0) nsrc = 1;

    // Fit it with a Gaussian
    auto model_func = [&](vec1d x, vec1d p) {
        vec1d f = replicate(p[p.size()-2], x.size());
        f += p[p.size()-1]*x;

        for (uint_t s : range(nsrc)) {
            f += integrate_gauss(x-0.5, x+0.5,
                p[1*nsrc+s], sqrt(sqr(p[2*nsrc+s]) + sqr(seeing_pix)), p[s]
            );
        }

        return f;
    };

    vec1d sp = replicate(1.0, nsrc);
    append(sp, default_gauss_pos);
    append(sp, default_gauss_width);
    sp.push_back(0.0);
    sp.push_back(0.0);

    mpfit_result fitres;
    if (!gauss_nofit && !badsnr) {
        if (verbose) {
            note("starting parameters: ", sp);
        }

        mpfit_options opts = mpfit_options(sp.size());
        if (gauss_nofit_width) {
            opts.frozen[(2*nsrc)-_-(3*nsrc-1)] = true;
        }
        if (!gauss_fit_background) {
            opts.frozen[sp.size()-2] = true;
        }
        if (!gauss_fit_background_linear) {
            opts.frozen[sp.size()-1] = true;
        }

        bool converged = false;
        bool first = true;
        do {
            fitres = mpfitfun(sprof, serror, sx, model_func, sp, opts);

            if (!fitres.success) {
                fitres.params = sp; // default values...
                warning("fit of collapsed spatial profile in ", band,
                    " failed, using default values: ", fitres.params);
            } else {
                converged = true;
                for (uint_t s : range(nsrc)) {
                    // Make sure the width is positive
                    // (mpfit will allow it to go negative because it is always squared in the model)
                    fitres.params[2*nsrc+s] = abs(fitres.params[2*nsrc+s]);

                    if (fitres.params[2*nsrc+s]/fitres.errors[2*nsrc+s] < 3.0) {
                        converged = false;
                        vec1u ids = {s, 1*nsrc+s, 2*nsrc+s};
                        vec1d oparams = fitres.params[ids];
                        fitres.params[ids] = sp[ids]; // default values...
                        warning("fit of collapsed spatial profile for source",
                            (sources.empty() ? "" : " "+sources[s]), " in ", band,
                            " has low S/N, using default values: ", fitres.params[ids]);
                        note("fitted values were: ", oparams, " (uncertainties: ", fitres.errors[ids], ")");
                    }
                }
            }

            // On the second iteration, force convergence
            if (!first) converged = true;

            if (!converged && first) {
                if (verbose) note("trying a two-pass fit fixing shapes first...");

                first = false;

                // Fix the widths and positions to at least get amplitude right
                vec1b old_freeze = opts.frozen;
                opts.frozen[(1*nsrc)-_-(3*nsrc-1)] = true;

                fitres = mpfitfun(sprof, serror, sx, model_func, fitres.params, opts);

                // Now refit
                opts.frozen = old_freeze;
                sp[_-(nsrc-1)] = fitres.params[_-(nsrc-1)];
            }
        } while (!converged);

        if (verbose) note("model for collapsed spatial profile in ", band, ": ", fitres.params,
            " (uncertainties: ", fitres.errors, ")");
    } else {
        fitres.params = sp; // default values...
        if (verbose && is_finite(psnr)) {
            note("using default model for the profile in ", band, " : ", sp);
        }
    }

    // Eliminate background from model
    fitres.params[fitres.params.size()-2] = 0.0;
    fitres.params[fitres.params.size()-1] = 0.0;

    // Save source profiles
    {
        fits::output_image oimgp(out_dir+"stacked_"+name+"_profile.fits");
        oimgp.reach_hdu(1);
        oimgp.write(sprof);
        oimgp.write_keyword("NSRC", nsrc);
        oimgp.reach_hdu(2);
        oimgp.write(serror);

        vec1d oamp = fitres.params[_-(nsrc-1)];

        for (uint_t s : range(nsrc)) {
            fitres.params[_-(nsrc-1)] = 0.0;
            fitres.params[s] = oamp[s];
            vec1d model = model_func(sx, fitres.params);
            oimgp.reach_hdu(3+s);
            oimgp.write(model);
        }
    }

    // Define models for fit
    uint_t nmodel = nsrc;
    if (gauss_fit_background) {
        ++nmodel;
        if (gauss_fit_background_linear) {
            ++nmodel;
        }
    }

    vec2f models(nmodel, tflxfit.dims[0]);

    if (gauss_fit_background) {
        models(nsrc+0,_) = 1.0;
        if (gauss_fit_background_linear) {
            models(nsrc+1,_) = sx;
        }
    }

    for (uint_t s : range(nsrc)) {
        // Build models of unit integral
        fitres.params[_-(nsrc-1)] = 0.0;
        fitres.params[s] = 1.0;
        models(s,_) = model_func(sx, fitres.params);
    }

    vec2d flx1d = replicate(dnan, nsrc, tflxfit.dims[1]);
    vec2d err1d = replicate(dnan, nsrc, tflxfit.dims[1]);
    vec1d bg1d, bg1dx;
    if (gauss_fit_background) {
        bg1d = replicate(dnan, tflxfit.dims[1]);
        if (gauss_fit_background_linear) {
            bg1dx = replicate(dnan, tflxfit.dims[1]);
        }
    }

    // Fit each wavelength slice
    for (uint_t l : range(lam)) {
        vec1d ttflx = tflxfit.safe(_,l);
        vec1d tterr = terrfit.safe(_,l);
        vec1u idg = where(is_finite(ttflx) && is_finite(tterr) && tterr > 0);

        if (idg.empty()) continue;

        // Re-normalize by median error to avoid numerical accuracy issues with big exponents
        double m = median(tterr[idg]);
        ttflx /= m;
        tterr /= m;

        linfit_result res = linfit_pack(ttflx.safe[idg], tterr.safe[idg], models.safe(_,idg));

        if (res.success) {
            flx1d.safe(_,l) = m*res.params.safe[_-(nsrc-1)];
            err1d.safe(_,l) = m*res.errors.safe[_-(nsrc-1)];
            if (gauss_fit_background) {
                if (gauss_fit_background_linear) {
                    bg1d.safe[l] = m*res.params[res.params.size()-2];
                    bg1dx.safe[l] = m*res.params[res.params.size()-1];
                } else {
                    bg1d.safe[l] = m*res.params[res.params.size()-1];
                }
            }
        }
    }

    // Save spectra
    for (uint_t s : range(nsrc)) {
        std::string tsuffix;
        if (!sources.empty()) tsuffix = "_"+sources[s];
        fits::output_image oimg(out_dir+"stacked_"+name+tsuffix+".fits");

        auto write_wcs = [&]() {
            oimg.write_keyword("CTYPE1", "WAVE");
            oimg.write_keyword("CUNIT1", "um");
            oimg.write_keyword("CRPIX1", 1.0);
            oimg.write_keyword("CRVAL1", lam[0]);
            oimg.write_keyword("CDELT1", lam[1]-lam[0]);

            oimg.write_keyword("GWIDTH", fitres.params[2*nsrc+s]*aspix);
            oimg.write_keyword("GPOS",   fitres.params[1*nsrc+s]*aspix);
            oimg.write_keyword("EXPTIME", exptime);
            oimg.write_keyword("SEEING", seeing);

            oimg.write_keyword("BAND", band);
            oimg.write_keyword("INSTR", instrument);
        };

        oimg.reach_hdu(1);
        oimg.write(flx1d(s,_));
        write_wcs();
        oimg.reach_hdu(2);
        oimg.write(err1d(s,_));
        write_wcs();
    }

    auto write2d = [&](std::string ofile, const vec2d& oflx, const vec2d& oerr) {
        fits::output_image oimg2d(ofile);

        auto write_wcs_2d = [&]() {
            oimg2d.write_keyword("CTYPE1", "WAVE");
            oimg2d.write_keyword("CUNIT1", "um");
            oimg2d.write_keyword("CRPIX1", 1.0);
            oimg2d.write_keyword("CRVAL1", lam[0]);
            oimg2d.write_keyword("CDELT1", lam[1]-lam[0]);

            oimg2d.write_keyword("CTYPE2", "SLIT");
            oimg2d.write_keyword("CUNIT2", "arcsec");
            oimg2d.write_keyword("CRPIX2", dpix+1);
            oimg2d.write_keyword("CRVAL2", 0.0);
            oimg2d.write_keyword("CDELT2", aspix);

            oimg2d.write_keyword("BAND", band);
            oimg2d.write_keyword("INSTR", instrument);
            oimg2d.write_keyword("EXPTIME", exptime);
            oimg2d.write_keyword("SEEING", seeing);
        };

        oimg2d.write_empty();
        oimg2d.reach_hdu(1);
        oimg2d.write(oflx);
        write_wcs_2d();
        oimg2d.reach_hdu(2);
        oimg2d.write(oerr);
        write_wcs_2d();
    };

    // Build background model
    sx = findgen(2*dpix+1) - dpix;
    vec2d bg(tflx.dims);
    if (gauss_fit_background) {
        if (gauss_fit_background_linear) {
            for (uint_t l : range(lam)) {
                bg.safe(_,l) = bg1d.safe[l] + sx*bg1dx.safe[l];
            }
        } else {
            for (uint_t l : range(lam)) {
                bg.safe(_,l) = bg1d.safe[l];
            }
        }
    }

    // Subtract (does nothing if no background fitted)
    tflx -= bg;

    // Save bg sub residual
    if (gauss_fit_background) {
        write2d(out_dir+"stacked_"+name+"_spec2d_bgsub.fits", tflx, terr);
    }

    // Save full residual
    vec2d rflx = tflx;
    for (uint_t s : range(nsrc)) {
        // Build models of unit integral
        fitres.params[_-(nsrc-1)] = 0.0;
        fitres.params[s] = 1.0;
        vec1d model = model_func(sx, fitres.params);

        for (uint_t l : range(lam)) {
            rflx.safe(_,l) -= model*flx1d.safe(s,l);
        }
    }

    write2d(out_dir+"stacked_"+name+"_spec2d_residual.fits", rflx, terr);

    if (!sources.empty()) {
        for (uint_t s : range(nsrc)) {
            rflx = tflx;
            for (uint_t os : range(nsrc)) {
                if (os == s) continue;

                // Build models of unit integral
                fitres.params[_-(nsrc-1)] = 0.0;
                fitres.params[os] = 1.0;
                vec1d model = model_func(sx, fitres.params);

                for (uint_t l : range(lam)) {
                    rflx.safe(_,l) -= model*flx1d.safe(os,l);
                }
            }

            std::string tsuffix = "_"+sources[s];
            write2d(out_dir+"stacked_"+name+tsuffix+"_spec2d_residual.fits", rflx, terr);
        }
    }

    return 0;
}

void print_help() {
    using namespace terminal_format;

    print("extract2d v1.0");
    print("usage: extract2d flux_file=<spectrum.fits> instrument=... [options]");
    print("");
    print("Main parameters:");
    paragraph("'spectrum.fits' must be a valid 2D spectrum file, i.e., a FITS file containing "
        "a 2D image. The spectral axis (wavelength or frequency) must be the X axis, and the "
        "units must be adequately determined in the WCS information in the FITS header. You must "
        "also specify from which instrument the spectrum was taken, as this will set a number "
        "of internal parameters or tweaks that are necessary for this instrument. Curently "
        "supported instruments are X-SHOOTER and MOSFIRE.");

    print("");
    print("Output format:");
    paragraph("Each output 1D spectrum will contain the extracted flux in the first FITS "
        "extension, and the formal uncertainty in the second extension. The program will create one "
        "1D spectrum per source specified in source=[...]. WCS keywords are set to define the "
        "spectral axis units and format. The fitting parameters of the source's profile are also "
        "saved as keywords ('GPOS': position, 'GWIDTH': width).");
    paragraph("The program will also 2D spectra, extracted around the reference position and with "
        "all the rescaling and flagging applied, for inspection. The format is the same as the 1D "
        "spectra. In addition, the program will save residual 2D spectra, with all sources removed, "
        "or with only one source removed, for inspection.");

    print("");
    print("Specifying input data:");
    bullet("flux_file=...", "Path to the FITS file containing the flux data.");
    bullet("error_file=...", "Path to the FITS file containing the flux uncertainty data. If not "
        "specified, it will be the same file as given in 'flux_file' (the uncertainties are then "
        "read from a different FITS extension, see 'error_hdu' below).");
    bullet("flag_file=...", "Path to the FITS file containing the pixel flags. If not specified, "
        "it will be the same file as given in 'flux_file' (the uncertainties are then read from a "
        "different FITS extension, see 'flag_hdu' below).");
    bullet("flux_hdu=...", "Must be a number. Defines the FITS \"header data unit\" (HDU, also "
        "known as \"FITS extension\") from which to read the flux data in 'flux_file'. Default "
        "is 1 (first extension), which is the ESO convention.");
    bullet("error_hdu=...", "Must be a number. Defines the FITS \"header data unit\" (HDU, also "
        "known as \"FITS extension\") from which to read the flux uncertainty data in "
        "'error_file'. Default is 2 (second extension), which is the ESO convention.");
    bullet("flag_hdu=...", "Must be a number. Defines the FITS \"header data unit\" (HDU, also "
        "known as \"FITS extension\") from which to read the pixel flags data in "
        "'flags_file'. Default is -1, which means there are no flags to read.");
    bullet("instrument=...", "Must be the name of the instrument which was used to take the 2D "
        "spectrum. Accepted values: X-SHOOTER, MOSFIRE.");

    print("");
    print("Reduction and flagging options:");
    bullet("max_flag=...", "Must be a number. If pixel flags are provided (see 'flag_hdu'), this "
        "sets the maximum allowed flag value for valid pixels. If the flag value is larger, the "
        "corresponding pixel will be excluded in the final spectrum (possibly creating NaN pixel).");
    bullet("rms_rescale=...", "Must be a number. This number will be multiplied to the "
        "input uncertainties, which can be used to manually correct under-estimated or "
        "over-estimated uncertainty spectra. This can have an impact especially if 'auto_mask' is "
        "set (see below).");
    bullet("auto_mask", "Set this flag if you want to automatically identify bad pixels in the "
        "data. The algorithm is described in the Appendix of Schreiber et al. (2018), and "
        "essentially locates pixels whose value are larger than X sigma from the expected noise "
        "around them. This algorithm expects that the flux of the target and its emission lines "
        "are either significantly fainter than the sky noise, or significantly extended in the "
        "spectral direction (i.e., that the emission lines are wider than one pixel). This "
        "flagging is not enabled by default.");
    bullet("auto_mask_dl=...", "If 'auto_mask' is set, this defines the width (in pixels) in the "
        "spectral dimension over which the local noise and background will be estimated. The larger "
        "this value, the more likely emission lines will get flagged out (which is bad). The "
        "smaller the value, the more noisy the flagging will be. Default is 12 pixels, which is "
        "adequate for R=3000 spectra (e.g., MOSFIRE), while 30 pixels is more adequate for high "
        "resolution spectra (e.g., X-SHOOTER).");
    bullet("auto_mask_sigma=...", "If 'auto_mask' is set, this sets the maximum allowed deviation "
        "of pixels before they are flagged. Default is 4 sigma.");
    bullet("crop_border=...", "Must be an integer. This defines how many pixels are ignored at the "
        "beginning and end of the 2D spectrum, which can be useful to remove bad edges. Default is "
        "0, which means no pixel is excluded.");
    bullet("flag_begin=...", "Must be a number. This defines the value of the spectral axis below "
        "which all pixels will be flagged. This is an alternative to 'crop_border'. Default is NaN, "
        "which means no pixel is flagged.");
    bullet("flag_end=...", "Must be a number. This defines the value of the spectral axis above "
        "which all pixels will be flagged. This is an alternative to 'crop_border'. Default is NaN, "
        "which means no pixel is flagged.");
    bullet("lmin=...", "Must be a number. It sets the minimum point on the spectral axis for which "
        "to extract data. The difference between this and 'flag_begin' is that 'flag_begin' will "
        "flag the pixels (set to NaN) and keep them in the output data, while 'lmin' will truncate "
        "the output data. Default is NaN, which means use the full range of the input data.");
    bullet("lmax=...", "Must be a number. It sets the minimum point on the spectral axis for which "
        "to extract data. The difference between this and 'flag_end' is that 'flag_end' will "
        "flag the pixels (set to NaN) and keep them in the output data, while 'lmax' will truncate "
        "the output data. Default is NaN, which means use the full range of the input data.");
    bullet("no_ohflag", "Set this flag to disable flagging of OH line residuals. By default, a "
        "pre-determined list of OH lines will be masked from the spectrum, since the residuals are "
        "often imperfect and larger than expected from the noise (e.g., because of atmospheric or "
        "telescope response variations between the source and sky positions). Set this to zero "
        "to disable OH line flagging.");
    bullet("oh_line_width=...", "This sets the expected width of OH lines on the spectral axis, "
        "in km/s (c*Delta_lambda/lambda). Default is 100 km/s, which is adequate for R=3000 "
        "spectra (e.g., MOSFIRE), while 60 km/s is a better value for high resolution spectra "
        "(e.g., X-SHOOTER).");
    bullet("no_gapflag", "Set this flag to disable flagging transmission gaps in between "
        "atmospheric transmission windows. This is only relevant for instruments with wide band "
        "capabilities, such as X-SHOOTER, KMOS, or SINFONI, and not for single-band instruments "
        "like MOSFIRE.");
    bullet("telluric=...", "Provide the path to a 1D FITS spectrum containing the telluric "
        "correction. The values in this 1D spectrum will be divided from the input 2D spectrum. "
        "If not provided, no telluric correction is performed. This can also be used to convert "
        "the 2D spectrum to an absolute flux reference.");
    bullet("median_background", "Set this flag to subtract the median of each spectral slice in "
        "the 2D spectrum before the extraction. This is necessary if sky subtraction is imperfect "
        "(or absent), and will only work correctly if the source is not filling all the slit in "
        "the spatial direction (there must be enough sky pixels). This is disabled by default. "
        "This uses the default profiles of the sources set in 'default_gauss_pos' and "
        "'default_gauss_width' (see below) to mask the source before computing the background.");

    print("");
    print("Spatial axis options:");
    bullet("offset=...", "Must be a number. This sets a global shift of the astrometry in the "
        "spatial axis, to correct, e.g., for telescope pointing accuracy or inaccurate input "
        "coordinates. Must be given in pixels. For a fixed expected sky position, this number is "
        "added to the corresponding pixel coordinate in the 2D spectrum (this effectively \"moves "
        "the spectrum down\" when the value is positive). Default is zero, and no shift is applied.");
    bullet("discard_wcs_pos", "Set this flag to ignore the WCS data for the spatial direction, and "
        "force the center of the slit as reference position (0 arcsec offset). Can be useful in "
        "cases the WCS data is incorrect or missing.");
    bullet("dpix=...", "Must be an integer. It specifies the half-size of the extraction window "
        "along the spatial dimension (total of 2*dpix+1 pixels will be extracted around the "
        "reference position, see above). Default is 10 pixels, which means a total of 21 pixels "
        "in the spatial direction.");
    bullet("seeing=...", "Must be a number. The instrument PSF is assumed to be seeing-limited, "
        "and Gaussian. This number sets the width of the PSF (FWHM) in arcsec. Default is zero, "
        "which assumes there is no PSF.");

    print("");
    print("Source extraction options:");
    bullet("name=...", "Must be a string. This gives a name to the extracted source (or group of "
        "sources), which will be added to the output file name. It must be provided.");
    bullet("sources=[...]", "Must be a list of strings. This gives a name to each source that will "
        "be extracted. It can be left empty if there is only one source.");
    bullet("dpixfit=...", "Must be an integer. It specifies the half-size of the fitting region on "
        "the spatial direction, in pixels (total of 2*dpixfit+1 pixels). Default is 6. This can be "
        "used to effectively ignore pixels far away from the source(s).");
    bullet("default_gauss_pos=[...]", "Must be a list of numbers or a single number. This sets the "
        "starting spatial position (along the slit) of each source, in arcsec, and relative to "
        "the reference spatial position (see 'offset' and 'discard_wcs_pos' above). This list "
        "must contain as many values as the 'sources' list. The value will be used as starting "
        "point for the fitting algorithm, and as fallback value in case the fit failed or when the "
        "S/N is too low.");
    bullet("default_gauss_width=[...]", "Must be a list of numbers or a single number. This sets "
        "the starting spatial width (along the slit) of each source, which is defined as the sigma "
        "of a Gaussian profile, in arcsec. If 'seeing' is set to a non-zero value, the seeing will "
        "be added in quadrature to this width to model the source profile. This list must contain "
        "as many values as the 'sources' list. The value will be used as starting point for the "
        "fitting algorithm, and as fallback value in case the fit failed or when the S/N is too "
        "low.");
    bullet("gauss_min_snr=...", "Must be a number. This defines the minimum allowed S/N to perform "
        "a fit of the source profile. Below this S/N, the fitting will be skipped and the default "
        "values listed in 'default_gauss_pos' and 'default_gauss_width' will be used. The S/N is "
        "computed as the S/N of the central pixel in the spectrally-averaged 2D spectrum. Default "
        "value is 6.");
    bullet("gauss_nofit", "Set this flag to skip the fitting stage, and directly use the default "
        "source profiles listed in 'default_gauss_pos' and 'default_gauss_width'.");
    bullet("gauss_nofit_width", "Set this flag to only fit for the source positions, while keeping "
        "the source widths frozen to their default values.");
    bullet("gauss_fit_background", "Set this flag to allow a constant level of background flux when "
        "fitting the source profile, and later on when fitting each spectral slice.");
    bullet("gauss_fit_background_linear", "If 'gauss_fit_background' is set, set this flag to "
        "model the background as a linear function, rather than constant pedestal.");

    print("");
    print("Output options:");
    bullet("out_dir=...", "Must be a string. This specifies the output directory, in which all "
        "the data products will be saved. Default is current directory. The program will create "
        "the directory if it does not exist.");
    bullet("verbose", "Set this flag to print the progress of the program in the terminal. "
        "Useful for figuring out what is going on, or for debugging.");
}
