#include <phypp.hpp>
#include <phypp/astro/template_fit.hpp>

int phypp_main(int argc, char* argv[]) {
    uint_t sid = npos;
    std::string catalog;
    std::string library;
    double lambda_min = 0.8;
    double lambda_max = 3.0;

    read_args(argc, argv, arg_list(name(sid, "id"), catalog, library, lambda_min, lambda_max));

    vec1u id;
    vec2f flux, flux_err;
    vec1s bands;
    fits::read_table(catalog, ftable(id, flux, flux_err, bands));

    sid = where_first(id == sid);
    if (sid == npos) {
        error("could not find source ID");
        return 1;
    }

    // Load filter response curves
    std::string filter_db_file = data_dir+"fits/filter-db/db.dat";
    if (!file::exists(filter_db_file)) {
        filter_db_file = system_var("PHYPP_FILTERS_PATH", ".")+"/db.dat";
        if (!file::exists(filter_db_file)) {
            error("cannot find filter database");
            note("looked at:");
            note(" - ", data_dir+"fits/filter-db/db.dat");
            note(" - ", filter_db_file);
            return 1;
        }
    }

    auto fdb = read_filter_db(filter_db_file);
    filter_bank_t fbank(bands.size());
    for (uint_t b : range(bands)) {
        get_filter(fdb, bands[b], fbank[b]);
    }

    // Load stellar SED library
    struct {
        vec2f lam, sed;
    } tlib;

    fits::read_table(library, "lam", tlib.lam, "sed", tlib.sed);

    float milam = min(tlib.lam);
    float malam = max(tlib.lam);

    vec1f lambda(fbank.size());
    vec1b empty_filter(fbank.size());
    for (uint_t i : range(fbank)) {
        lambda[i] = fbank[i].rlam;

        // Truncate filters to match the stellar library range
        vec1u idin = where(fbank[i].lam > milam && fbank[i].lam < malam);
        fbank[i].res = fbank[i].res[idin];
        fbank[i].lam = fbank[i].lam[idin];
        empty_filter[i] = idin.empty();
    }

    // Find chi2 for each template
    auto seed = make_seed(42);
    vec1u ifit = where(flux_err(sid,_) > 0 && is_finite(flux(sid,_)) &&
        lambda > lambda_min && lambda < lambda_max && !empty_filter);

    print("using ", ifit.size(), "/", flux.dims[1], " photometric points");

    template_fit_params tparams;
    tparams.ulim = false;
    tparams.nsim = 0;
    tparams.renorm = true;
    tparams.lib_obs = true;

    flux_err(sid,ifit) = max(flux_err(sid,ifit), 0.1*flux(sid,ifit));

    auto res = template_fit(tlib, seed,
        0.0, 1.0, vec1d{flux(sid,ifit)}, vec1d{flux_err(sid,ifit)}, fbank[ifit], tparams);

    print("best chi2: ", res.chi2[res.bfit]);

    fits::write_table("best_fit_slit_star.fits",
        "bfit", res.bfit, "chi2", res.chi2, "amp", res.amp,
        "flux_model", res.flux, "ifit", ifit,
        "flux", uJy2cgs(lambda, flux(sid,_)), "flux_err", uJy2cgs(lambda, flux_err(sid,_)),
        "lambda", lambda
    );

    vec1d lam = tlib.lam(res.bfit,_);
    vec1d flx = res.amp[res.bfit]*tlib.sed(res.bfit,_);
    flx = uJy2cgs(lam, flx);

    vec1s hdr = {"lambda[um]", "flambda[cgs]"};
    ascii::write_table_hdr("best_fit_slit_star.dat", 22, hdr, lam, flx);

    return 0;
}
