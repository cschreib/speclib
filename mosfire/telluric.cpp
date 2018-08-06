#include <phypp.hpp>

void wavelength_smooth(vec1d& flx, vec1d& err, uint_t smooth) {
    vec1d owei = 1.0/sqr(err);
    vec1d oflx = flx;
    vec1d oerr = err;

    vec1u idb = where(!is_finite(oflx) || !is_finite(owei) || oerr == 0.0);
    oflx[idb] = 0;
    oerr[idb] = 0;
    owei[idb] = 0;

    uint_t d1 = (smooth-1)/2;
    uint_t d2 = smooth-1 - d1;

    for (uint_t l : range(flx)) {
        uint_t i0 = (l > d1 ?            l - d1 : 0);
        uint_t i1 = (l < flx.size()-d2 ? l + d2 : flx.size()-1);

        double twei = total(owei[i0-_-i1]);
        flx[l] = total(owei[i0-_-i1]*oflx[i0-_-i1])/twei;
        err[l] = sqrt(total(sqr(owei[i0-_-i1]*oerr[i0-_-i1])))/twei;
    }
}

int phypp_main(int argc, char* argv[]) {
    if (argc < 2) {
        print("usage:");
        print("  telluric star=[ascii_star_spectrum.dat] files=[files.fits...] filters=[filter_dir]");
        return 0;
    }

    std::string star_spec;
    vec1s files;
    std::string filter_dir;
    bool verbose = false;
    bool no_baseline = false;
    bool no_throughput = false;
    uint_t smooth = 0;
    uint_t lambda_padding = 0;

    read_args(argc, argv, arg_list(
        name(star_spec, "star"), files, name(filter_dir, "filters"), no_baseline,
        no_throughput, verbose, smooth, lambda_padding
    ));

    if (star_spec.empty()) {
        error("need star=... parameter to point to an ASCII file with the star's theoretical spectrum");
        note("this spectrum must be in units of erg/s/cm2/A, rescaled to the broad band flux of the star");
        return 1;
    }

    if (files.empty()) {
        note("no spectrum to correct");
        return 0;
    }

    if (filter_dir.empty()) {
        error("missing filters=... parameter");
        error("please using this parameter to specify the path to the MOSFIRE filters");
        return 1;
    }

    if (verbose) {
        note("computing correction for ", files.size(), " spectra");
    }

    std::string suffix;
    if (smooth > 1) {
        note("smoothing spectra by ", smooth, " resolution elements");
        suffix += "_s"+to_string(smooth);
    }

    // Read the star's intrinsic spectrum
    vec1d tlam, tsed;
    ascii::read_table(star_spec, ascii::find_skip(star_spec), tlam, tsed);

    // Read the MOSFIRE filter's throughput
    vec1s filter_band;
    vec<1,filter_t> filters; {
        filter_dir = file::directorize(filter_dir);
        vec1s ffiles = file::list_files(filter_dir+"*.dat");
        for (uint_t i : range(ffiles)) {
            filter_t f;
            std::string band = file::remove_extension(ffiles[i]);
            ascii::read_table(filter_dir+ffiles[i], ascii::find_skip(filter_dir+ffiles[i]),
                f.lam, f.res
            );

            vec1u ids = sort(f.lam);
            f.lam = f.lam[ids];
            f.res = f.res[ids];
            f.res /= max(f.res);

            filter_band.push_back(band);
            filters.push_back(f);
        }
    }

    // Read the average transmission at Mauna Kea
    vec1d omk_lam, omk_cor;
    ascii::read_table(filter_dir+"/mauna_kea_avg.dat", ascii::find_skip(filter_dir+"mauna_kea_avg.dat"),
        omk_lam, omk_cor
    );

    // Read the observed spectrum (just the first one, to get wavelength grid)
    vec1d oflx, oerr;
    fits::input_image iimg(files[0]);
    iimg.reach_hdu(1);
    iimg.read(oflx);

    double crpix = 0, crval = 0, cdelt = 0;
    if (!iimg.read_keyword("CRPIX1", crpix) ||
        !iimg.read_keyword("CRVAL1", crval) ||
        !iimg.read_keyword("CDELT1", cdelt)) {
        error("could not read WCS");
        return 1;
    }

    vec1d olam = (dindgen(oflx.size()) + 1 - crpix)*cdelt + crval;

    // First rebin the template and throughputs to the resolution of the star spectrum
    if (verbose) note("rebinning intrinsic spectrum");
    vec1d sflx(oflx.size());
    for (uint_t i : range(olam)) {
        sflx[i] = integrate(tlam, tsed, olam[i]-0.5*cdelt, olam[i]+0.5*cdelt)/cdelt;
    }

    // Smooth it
    if (smooth > 1) {
        vec1d terr = replicate(1.0, sflx.dims);
        wavelength_smooth(sflx, terr, smooth);
    }

    if (verbose) note("rebinning filters");
    for (uint_t b : range(filters)) {
        vec1f nres(oflx.size());
        nres[_] = dnan;

        vec1u idg = where(is_finite(1/filters[b].res));
        double milam = min(filters[b].lam[idg]);
        double malam = max(filters[b].lam[idg]);

        for (uint_t i : range(olam)) {
            if (olam[i]-0.5*cdelt > milam && olam[i]+0.5*cdelt < malam) {
                nres[i] = integrate(filters[b].lam, filters[b].res, olam[i]-0.5*cdelt, olam[i]+0.5*cdelt)/cdelt;
            }
        }

        filters[b].lam.clear();
        filters[b].res = nres;
    }

    if (verbose) note("rebinning & filtering Mauna Kea average transmission");
    vec1d mk_cor(oflx.size());
    for (uint_t i : range(olam)) {
        double k = integrate(omk_lam, omk_cor, olam[i]-0.5*cdelt, olam[i]+0.5*cdelt)/cdelt;
        if (i > 0) mk_cor[i-1] += 0.25*k;
        mk_cor[i] += 0.5*k;
        if (i < olam.size()-1) mk_cor[i+1] += 0.25*k;
    }

    mk_cor = 1/mk_cor;

    // Smooth it
    if (smooth > 1) {
        vec1d terr = replicate(1.0, sflx.dims);
        wavelength_smooth(mk_cor, terr, smooth);
    }

    if (!no_baseline) {
        mk_cor -= 1.0;
    }

    // Compute correction for each spectrum
    if (verbose) note("computing raw corrections");
    vec2d ccor(oflx.size(), files.size());
    vec2d cbas(oflx.size(), files.size());
    vec2d cerr(oflx.size(), files.size());
    vec2d ccov(oflx.size(), files.size());
    vec1s band(files.size());

    for (uint_t i : range(files)) {
        // Read the observed spectrum
        fits::input_image tiimg(files[i]);
        tiimg.reach_hdu(1);
        tiimg.read(oflx);
        tiimg.reach_hdu(2);
        tiimg.read(oerr);

        // Apply padding if requested
        if (lambda_padding > 0) {
            uint_t i0 = where_first(is_finite(oflx) && is_finite(oerr) && oerr > 0);
            oflx[_-(i0+lambda_padding-1)] = dnan;
            oerr[_-(i0+lambda_padding-1)] = dnan;
            uint_t i1 = where_last(is_finite(oflx) && is_finite(oerr) && oerr > 0);
            oflx[(i1-(lambda_padding-1))-_] = dnan;
            oerr[(i1-(lambda_padding-1))-_] = dnan;
        }

        // Smooth it
        if (smooth > 1) {
            wavelength_smooth(oflx, oerr, smooth);
        }

        // Compute the correction
        vec1f tcor = sflx/oflx;
        vec1f terr = oerr*abs(sflx/sqr(oflx));
        vec1f tcov = 1/sqr(oerr);
        vec1b good = is_finite(tcor) && is_finite(terr);
        vec1b inband = olam >= min(olam[where(good)]) && olam <= max(olam[where(good)]);
        vec1u idb = where(!good && inband);
        vec1u idg = where(good);

        // Identify the band
        double mlam = mean(olam[where(inband)]);
        if (mlam < 1.5 || mlam > 2.4) {
            error("this code is only written for H or K band");
            error("got an exposure with mean wavelength ", mlam, " um");
            return 1;
        }

        band[i] = (mlam > 1.9 ? "K" : "H");

        // Fill the gaps
        auto p = interpolate(tcor[idg], terr[idg], olam[idg], olam[idb]);
        tcor[idb] = p.first;
        terr[idb] = p.second;
        tcov[idb] = interpolate(tcov[idg], olam[idg], olam[idb]);

        // Remove throughput contribution, if possible
        if (!no_throughput) {
            uint_t idf = where_first(filter_band == band[i]);
            if (idf != npos) {
                tcor *= filters[idf].res;
                terr *= filters[idf].res;
            } else {
                warning("no throughput available for filter '", band[i], "'");
            }
        }

        vec1d tcornb = tcor;
        if (!no_baseline) {
            // // Fit a polynomial to the curve to determine the baseline, masking largest features
            // // This baseline is mostly residuals in the response of the grism that were not
            // // properly accounted for by the advertised throughput
            // vec1d fc = tcor;
            // vec1d fe = terr;

            // // Mask regions of strong telluric absorption, because we only want the baseline for now
            // vec1b mask;
            // if (band[i] == "K") {
            //     mask = abs(olam - 1.96) < 0.01 || abs(olam - 2.01) < 0.015 || abs(olam - 2.06) < 0.01  || olam < 1.94 || olam > 2.375;
            // } else {
            //     mask = abs(olam - 1.575) < 0.01 || olam < 1.46 || olam > 1.78;
            // }

            // // Down-weight the bad values and masked regions
            // vec1u idb = where(fc <= 0 || !is_finite(fc) || fe <= 0 || !is_finite(fe) || mask);
            // fe[idb] = median(fe[where(is_finite(fe) && fe > 0)])*1e9;
            // fc[idb] = 1.0;

            // vec1u idf = (band[i] == "K" ? where(olam > 1.9) : where(olam < 1.9));

            // // Define baseline model
            // const double lref = (band[i] == "K" ? 2.15 : 1.6);
            // vec2d xfit(porder+1, fc.size());
            // xfit(0,_) = 1.0;
            // for (uint_t k : range(porder)) {
            //     xfit(k+1,_) = pow(olam - lref, k+1);
            // }

            // // Do the fit
            // auto res = linfit_pack(fc[idf], fe[idf], xfit(_,idf));
            // if (!res.success) {
            //     warning("fit of the baseline failed");
            //     note("  fitting '", files[i], "'");
            //     note("  assuming a constant baseline of 1");
            //     res.params[_] = 0;
            //     res.params[0] = 1;
            // }

            // // Compute the baseline and subtract it from the spectrum
            // // What remains should be mostly zero everywhere except for regions of telluric absorption
            // // which will have positive values
            // vec1f baseline(olam.size());
            // for (uint_t k : range(res.params)) {
            //     baseline += res.params[k]*xfit(k,_);
            // }

            // tcor -= baseline;


            // Compute average of chunks without strong features
            vec1d bchunkl, bchunku;
            double milam, malam;
            if (band[i] == "K") {
                milam = 1.90;
                malam = 2.43;
                bchunkl = {1.980, 2.030, 2.090, 2.170, 2.220, 2.290, 2.330, 2.380};
                bchunku = {1.995, 2.045, 2.130, 2.190, 2.240, 2.310, 2.350, 2.400};
            } else {
                milam = 1.43;
                malam = 1.85;
                bchunkl = {1.520, 1.585, 1.625, 1.670, 1.770};
                bchunku = {1.550, 1.595, 1.632, 1.710, 1.783};
            }

            vec1d clam = replicate(dnan, bchunku.size());
            vec1d cavg = replicate(dnan, bchunku.size());
            for (uint_t c : range(clam)) {
                vec1u idbl = where(is_finite(tcor) && is_finite(terr) && terr > 0 &&
                    olam >= bchunkl[c] && olam <= bchunku[c]);

                // Ignore this chunk if there is too little good data
                if (!idbl.empty() && max(olam[idbl])-min(olam[idbl]) > 0.009) {
                    vec1d wei = 1.0/sqr(terr[idbl]);
                    cavg[c] = weighted_median(tcor[idbl], wei);
                    clam[c] = weighted_median(olam[idbl], wei);
                }
            }

            vec1u idbg = where(is_finite(cavg) && cavg > 0.0);
            clam = clam[idbg];
            cavg = cavg[idbg];

            // Build a baseline from it using a spline
            vec1d baseline;
            bool iterate = true;
            bool first = true;
            while (iterate) {
                if (cavg.size() >= 3) {
                    double cv = median(cavg);
                    baseline = cv*interpolate_3spline(cavg/cv, clam, olam);
                } else if (cavg.size() >= 1) {
                    baseline = replicate(median(cavg), tcor.size());
                } else {
                    error("no finite data available in this exposure");
                    error(files[i]);
                    return 1;
                }

                if (count(baseline < 0 && olam > milam && olam < malam) > 0) {
                    if (cavg.size() >= 3 && first) {
                        warning("baseline went negative!");
                        warning("  analyzing '", files[i], "'");
                        warning("will remove end points and try again, but check the slit star's spectrum");
                        first = false;

                        uint_t i0 = where_first(olam > milam);
                        if (i0 == npos) i0 = 0;
                        if (baseline[i0] < 0) {
                            cavg[0] = cavg[1];
                        }

                        uint_t i1 = where_first(olam > malam);
                        if (i1 == npos) i1 = baseline.size()-1;
                        if (baseline[i1] < 0) {
                            cavg[cavg.size()-1] = cavg[cavg.size()-2];
                        }
                    } else {
                        error("baseline went negative!");
                        error("  analyzing '", files[i], "'");
                        std::string fname = "/tmp/failed_fits_0.fits";
                        uint_t ntry = 0;
                        while (file::exists(fname) && ntry < 1000) {
                            ++ntry;
                            fname = "/tmp/failed_fits_"+to_string(ntry)+".fits";
                        }

                        fits::write_table(fname, ftable(olam, baseline, clam, cavg));
                        note("  data is saved in '"+fname+"' for inspection");
                        return 1;
                    }
                } else {
                    iterate = false;
                }
            }

            // Save it for later
            cbas(_,i) = baseline;

            // And subtract it from the curve
            tcor /= baseline;
            terr /= baseline;
            tcor -= 1.0;
        }

        // Store it
        ccor(_,i) = tcor;
        cerr(_,i) = terr;
        ccov(_,i) = tcov;

        // Save it
        fits::output_image oimg(file::remove_extension(files[i])+"_correction"+suffix+".fits");
        oimg.reach_hdu(0);
        oimg.write_empty();

        auto write_wcs = [&]() {
            oimg.write_keyword("CTYPE1", "WAVE");
            oimg.write_keyword("CUNIT1", "um");
            oimg.write_keyword("CDELT1", olam[1]-olam[0]);
            oimg.write_keyword("CRPIX1", 1.0);
            oimg.write_keyword("CRVAL1", olam[0]);
        };

        oimg.reach_hdu(1);
        oimg.write(tcornb);
        write_wcs();

        oimg.reach_hdu(2);
        oimg.write(tcor);
        write_wcs();

        oimg.reach_hdu(3);
        oimg.write(terr);
        write_wcs();

        oimg.reach_hdu(4);
        oimg.write(vec1d{cbas(_,i)});
        write_wcs();
    }

    if (verbose) note("build template curve");

    // const uint_t porder = 2;

    vec2d cbcor(oflx.size(), files.size());
    vec2d cbwei(oflx.size(), files.size());
    vec2d cbcov(oflx.size(), files.size());
    for (uint_t i : range(files)) {
        vec1f tcor = ccor(_,i);
        vec1f terr = cerr(_,i);
        vec1f tcov = ccov(_,i);

        // Normalize to a common reference amplitude to build the template
        // double amp;
        // if (band[i] == "K") {
        //     amp = median(abs(tcor[where(olam > 2.1 && olam < 2.3)]));
        // } else {
        //     amp = median(abs(tcor[where(olam > 1.6 && olam < 1.7)]));
        // }

        // tcor /= amp;
        // terr /= amp;

        // Save the curve for later stacking
        vec1d twei = 1/sqr(terr);
        vec1u idb = where(!is_finite(tcor) || !is_finite(terr) || terr <= 0);
        tcor[idb] = 0;
        twei[idb] = 0;
        tcov[idb] = 0;

        cbcor(_,i) = tcor;
        cbwei(_,i) = twei;
        cbcov(_,i) = tcov;
    }

    // Stack the template telluric absorption spectrum
    vec1d tpl_cor = partial_total(1, cbcor*cbwei)/partial_total(1, cbwei);
    vec1d tpl_core = sqrt(partial_total(1, sqr(cbcor - transpose(replicate(tpl_cor, cbcor.dims[1])))*cbwei)/partial_total(1, cbwei));
    // Compute the coverage fraction
    // (% of exposures in that band that cover a wavelength element)
    vec1d tpl_cov = partial_total(1, cbcov)/partial_max(1, cbcov)/partial_total(1, cbcov > 0);

    // Use the Mauna Kea average to fill the edges and gaps
    vec1b tpl_edge(tpl_cor.dims); {
        vec1d bl = {1.4, 1.85};
        vec1d bu = {1.85, 2.45};

        for (uint_t b : range(bl)) {
            // Flag out the uncovered and poorly covered edges
            vec1b inband = olam >= bl[b] && olam <= bu[b];
            vec1u idc = where(inband && is_finite(tpl_cov));
            if (idc.empty()) continue;
            double medcov = median(tpl_cov[idc]);
            double rmscov = 1.48*mad(tpl_cov[idc]);

            vec1b gcor = is_finite(tpl_cor) && tpl_cor > -1.0 && tpl_core > 0 && tpl_cov > medcov - 3*rmscov;
            uint_t il = where_first(gcor && inband);
            if (il != npos) {
                // Take a sample close to the edge and compute the average
                // and subtract the expected average from MK
                vec1u idl = where(gcor && olam >= olam[il] && olam < olam[il]*1.02);
                double avg = median(tpl_cor[idl]) - median(mk_cor[idl]);
                double avge = median(tpl_core[idl]/tpl_cor[idl]);

                // Fill the gap with MK
                vec1u idb = where(inband && olam < olam[il]);
                tpl_cor[idb] = mk_cor[idb] + avg;
                tpl_core[idb] = tpl_cor[idb]*avge;
                tpl_edge[idb] = true;
            }

            uint_t iu = where_last(gcor && inband);
            if (iu != npos) {
                // Take a sample close to the edge and compute the average
                // and subtract the expected average from MK
                vec1u idl = where(gcor && olam <= olam[iu] && olam > olam[iu]*0.98);
                double avg = median(tpl_cor[idl]) - median(mk_cor[idl]);
                double avge = median(tpl_core[idl]/tpl_cor[idl]);

                // Fill the gap with MK
                vec1u idb = where(inband && olam > olam[iu]);
                tpl_cor[idb] = mk_cor[idb] + avg;
                tpl_core[idb] = tpl_cor[idb]*avge;
                tpl_edge[idb] = true;
            }
        }
    }

    tpl_cor[where(!is_finite(tpl_cor) || tpl_cor < -0.1)] = 0;

    // Save master curve
    {
        fits::output_image oimg(file::get_directory(files[0])+"template_correction"+suffix+".fits");
        oimg.reach_hdu(0);
        oimg.write_empty();

        auto write_wcs = [&]() {
            oimg.write_keyword("CTYPE1", "WAVE");
            oimg.write_keyword("CUNIT1", "um");
            oimg.write_keyword("CDELT1", olam[1]-olam[0]);
            oimg.write_keyword("CRPIX1", 1.0);
            oimg.write_keyword("CRVAL1", olam[0]);
        };

        oimg.reach_hdu(1);
        oimg.write(tpl_cor);
        write_wcs();
        oimg.reach_hdu(2);
        oimg.write(tpl_core);
        write_wcs();
        oimg.reach_hdu(3);
        oimg.write(tpl_cov);
        write_wcs();
        oimg.reach_hdu(4);
        oimg.write(tpl_edge);
        write_wcs();
        oimg.reach_hdu(5);
        oimg.write(mk_cor);
        write_wcs();
        oimg.reach_hdu(6);
        oimg.write(sflx);
        write_wcs();
    }

    if (verbose) note("fit each curve with the master template");

    // Fit each curve with the template
    for (uint_t i : range(files)) {
        vec1f tcor = ccor(_,i);
        vec1f terr = cerr(_,i);

        // Renormalize the curve before fitting to avoid numerical imprecision
        // vec1u idg = where(is_finite(tcor) && is_finite(terr) && terr > 0);
        // double conv = median(abs(tcor[idg]));
        // tcor /= conv;
        // terr /= conv;

        // Define wavelength coverage
        double minlam = min(olam);
        double maxlam = max(olam);

        // Cut the template into independent regions
        vec1d chunks;
        if (band[i] == "K") {
            chunks = {1.951, 1.961, 1.98, 1.991, 2.009, 2.035, 2.06, 2.09, 2.2, 2.3};
        } else {
            chunks = {1.565, 1.585, 1.615, 1.66, 1.67};
        }

        // Make sure the regions end at the spectrum's boundaries
        {
            if (minlam < chunks[0]) {
                prepend(chunks, {minlam});
            }
            if (maxlam > chunks.back()) {
                append(chunks, {maxlam});
            }

            // This removes chunks completely, which is not good at all
            // uint_t i0 = lower_bound(minlam, chunks);
            // if (i0 != npos) {
            //     if (i0 == chunks.size()-1) {
            //         chunks.clear();
            //     } else {
            //         chunks = chunks[(i0+1)-_];
            //     }
            // }

            // i0 = upper_bound(maxlam, chunks);
            // if (i0 != npos) {
            //     if (i0 == 0) {
            //         chunks.clear();
            //     } else {
            //         chunks = chunks[_-(i0-1)];
            //     }
            // }

            // prepend(chunks, {minlam});
            // append(chunks, {maxlam});
        }

        const uint_t nchunk = chunks.size()-1;

        // Select good points
        vec1u idg = where(is_finite(tcor) && is_finite(terr) && terr > 0);

        // Filter bad points
        // vec1u idb = where(!is_finite(tcor) || !is_finite(terr) || terr <= 0);

        // tcor[idb] = 1.0;

        // tcor[idb] = 0.0;
        // terr[idb] = median(terr[where(is_finite(terr) && terr > 0)])*1e9;

        // Define the fit model
        const double lref = (band[i] == "K" ? 2.15 : 1.6);
        vec2d x, xe;
        vec1b chunk_covered;

        // Baseline
        uint_t first_chunk = 0;
        // if (!no_baseline) {
        //     x.resize(porder+1+nchunk, tcor.size());
        //     xe.resize(porder+1+nchunk, tcor.size());

        //     x(0,_) = 1.0;
        //     for (uint_t k : range(porder)) {
        //         x(k+1,_) = pow(olam - lref, k+1);
        //     }

        //     first_chunk = porder+1;

        //     // x.resize(nchunk*3, tcor.size());

        //     // for (uint_t k : range(nchunk)) {
        //     //     vec1d tmp = replicate(1.0, tcor.size());
        //     //     uint_t l0 = lower_bound(chunks[k], olam);
        //     //     if (l0 != npos) tmp[_-l0] = 0;
        //     //     uint_t l1 = upper_bound(chunks[k+1], olam);
        //     //     if (l1 != npos) tmp[l1-_] = 0;
        //     //     x(2*k+0,_) = tmp;
        //     //     x(2*k+1,_) = tmp*(olam - 0.5*(chunks[k]+chunks[k+1]));
        //     // }

        //     // first_chunk = 2*nchunk;
        // } else {
            x.resize(nchunk, tcor.size());
            xe.resize(nchunk, tcor.size());
        // }

        // Template absorption
        for (uint_t k : range(nchunk)) {
            vec1d tmp = tpl_cor;
            vec1d tmpe = tpl_core;

            auto b = bounds(olam, chunks[k], chunks[k+1]);
            if (b[0] != npos) {
                tmp[_-b[0]] = 0;
                tmpe[_-b[0]] = 0;
            }
            if (b[1] != npos) {
                tmp[b[1]-_] = 0;
                tmpe[b[1]-_] = 0;
            }

            x(k+first_chunk,_) = tmp;
            xe(k+first_chunk,_) = tmpe;
        }

        // Select properly covered chunks
        vec1u idcov = where(partial_total(1, sqr(x(_,idg))) > 0);

        if (idcov.empty()) {
            error("no covered chunk in '"+files[i]+"'");
            return 1;
        }

        auto res = linfit_pack(tcor[idg], terr[idg], x(idcov,idg));
        if (!res.success) {
            error("fit of the absorption template failed");
            note("  fitting '", files[i], "'");
            std::string fname = "/tmp/failed_fits_0.fits";
            uint_t ntry = 0;
            while (file::exists(fname) && ntry < 1000) {
                ++ntry;
                fname = "/tmp/failed_fits_"+to_string(ntry)+".fits";
            }

            fits::write_table(fname, ftable(tcor, terr, x, res.cov));
            note("  invalid fitted data is saved in '"+fname+"' for inspection");
            return 1;
        }

        // Enlarge res.params to original chunks
        if (idcov.size() != x.dims[0]) {
            vec1d tmp = res.params;
            res.params.resize(x.dims[0]);
            res.params[_] = dnan;
            res.params[idcov] = tmp;

            // Set amplitude of un-covered chunks to the average of the other chunks
            double avgamp = mean(res.params[idcov[where(idcov >= first_chunk)]]);
            vec1u idn = where(!is_finite(res.params));
            // res.params[idn] = 0.0; // for baseline
            // res.params[idn[where(idn >= first_chunk)]] = avgamp;
            res.params[idn] = avgamp;

            tmp = res.errors;
            res.errors.resize(x.dims[0]);
            res.errors[_] = 0;
            res.errors[idcov] = tmp;
        }

        // Adjust value of fit parameters
        // Amplitude cannot go negative or very low. When this happens, this is probably because the
        // curve has low S/N or poor coverage, so we use the average absorption shape
        for (uint_t k : range(nchunk)) {
            if (res.params[k+first_chunk] < 0.2) {
                res.params[k+first_chunk] = 1.0;
                res.errors[k+first_chunk] = 0.0;
            }
        }

        // Build the best fit model
        vec1d model(tcor.size());
        vec1d modele(tcor.size());
        // vec1d baseline(tcor.size());
        for (uint_t k : range(x.dims[0])) {
            model += res.params[k]*x(k,_);
            modele += sqr(res.params[k]*xe(k,_));
            // if (k < first_chunk) {
            //     baseline += res.params[k]*x(k,_);
            // }
        }

        modele = sqrt(modele);

        // Compute local chi2
        vec1d resid = (tcor - model)/terr;
        vec1d chi2(nchunk);
        for (uint_t k : range(nchunk)) {
            uint_t l0 = lower_bound(olam, chunks[k]);
            if (l0 == npos) l0 = 0;
            uint_t l1 = upper_bound(olam, chunks[k+1]);
            if (l1 == npos) l1 = olam.size()-1;

            vec1d lres = resid[l0-_-l1];
            lres = lres[where(is_finite(lres))];
            if (!lres.empty()) {
                chi2[k] = mean(sqr(lres))/lres.size();
            }
        }

        // Re-normalize to original amplitude
        // model *= conv;
        // modele *= conv;
        // baseline *= conv;

        // Re-add baseline
        vec1d baseline = cbas(_,i);
        if (!no_baseline) {
            model += 1.0;
            model *= baseline;
            modele *= baseline;
        }

        // Re-apply throughput, if used
        if (!no_throughput) {
            uint_t idf = where_first(filter_band == band[i]);
            if (idf != npos) {
                model /= filters[idf].res;
                modele /= filters[idf].res;
                baseline /= filters[idf].res;
            }
        }

        // Write model to disk
        fits::output_image oimg(file::remove_extension(files[i])+"_correction"+suffix+"_fit.fits");
        oimg.reach_hdu(0);
        oimg.write_empty();

        auto write_wcs = [&]() {
            oimg.write_keyword("CTYPE1", "WAVE");
            oimg.write_keyword("CUNIT1", "um");
            oimg.write_keyword("CDELT1", olam[1]-olam[0]);
            oimg.write_keyword("CRPIX1", 1.0);
            oimg.write_keyword("CRVAL1", olam[0]);
        };

        oimg.reach_hdu(1);
        oimg.write(model);
        write_wcs();

        // oimg.write_keyword("GNORM", conv);

        for (uint_t k : range(res.params)) {
            std::string type = (k < first_chunk && !no_baseline ? "baseline" : "feature");
            oimg.write_keyword("FITTYP"+align_right(to_string(k), 2, '0'), type);
        }

        for (uint_t k : range(res.params)) {
            oimg.write_keyword("FITPAR"+align_right(to_string(k), 2, '0'), res.params[k]);
            oimg.write_keyword("FITPER"+align_right(to_string(k), 2, '0'), res.errors[k]);
        }

        for (uint_t k : range(first_chunk, res.params.size())) {
            oimg.write_keyword("FITLLO"+align_right(to_string(k), 2, '0'), chunks[k-first_chunk]);
            oimg.write_keyword("FITLUP"+align_right(to_string(k), 2, '0'), chunks[(k+1)-first_chunk]);
        }

        for (uint_t k : range(first_chunk, res.params.size())) {
            oimg.write_keyword("FITCHI"+align_right(to_string(k), 2, '0'), chi2[k-first_chunk]);
        }

        oimg.reach_hdu(2);
        oimg.write(modele);
        write_wcs();

        oimg.reach_hdu(3);
        oimg.write(baseline);
        write_wcs();
    }

    if (verbose) note("done");

    return 0;
}
