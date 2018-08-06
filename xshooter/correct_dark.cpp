#include <phypp.hpp>

vec2b make_mask1(const vec2d& wave) {
    // Make a first subtraction excluding slits
    vec2b mask = wave > 0;
    for (uint_t i : range(mask.dims[1])) {
        double ty0 = max(-1.2*i + 600.0, 0.8*i - 1150);
        if (ty0 > 0) {
            mask(_-uint_t(ty0),i) = false;
        }
    }

    return mask;
}

vec2f make_sky(const vec2f& data, const vec2d& wave, const vec2b& mask) {
    segment_output sdo;
    vec2u seg = segment(mask, sdo);
    note("found ", sdo.id.size(), " slits");

    vec2d bg(seg.dims);
    foreach_segment(seg, sdo.origin, [&](uint_t sid, vec1u ids) {
        vec1d w = wave.safe[ids];
        vec1d d = data.safe[ids];
        vec1d b = bg.safe[ids];

        double mi = min(w);
        double ma = max(w);
        vec2d lbin = make_bins(rgen(mi, ma, 2048));

        histogram(w, lbin, [&](uint_t l, vec1u idl) {
            if (idl.size() > 1) {
                double t = median(d.safe[idl]);
                b.safe[idl] = t;
            }
        });

        bg.safe[ids] = b;
    });

    return bg;
}

vec2b make_mask2(const vec2f& data) {
    vec2b mask = abs(data) > 10;
    mask = mask_inflate(!mask_inflate(!mask_inflate(mask, 1), 2), 1);
    return mask;
}

void subtract_pattern(vec2f& data, vec2f& cor, const vec2b& mask) {
    // Second pass with larger data
    // Subtract X then Y patterns
    for (uint_t ix : range(data.dims[1])) {
        vec1u idg = where(!mask(_,ix));
        if (idg.size() > 10) {
            double d = median(data.safe(idg,ix));
            data(_,ix) -= d;
            cor(_,ix) += d;
        }
    }

    for (uint_t iy : range(data.dims[0])) {
        vec1u idg = where(!mask(iy,_));
        if (idg.size() > 10) {
            double d = median(data.safe(iy,idg));
            data(iy,_) -= d;
            cor(iy,_) += d;
        }
    }
}

vec2d from_pre(vec2d data) {
    vec2d ret = replicate(dnan, 1100, 2048);
    data = transpose(data);
    data = flip_x(data);

    const uint_t y0 = 20;
    ret(y0-_-(y0+data.dims[0]-1),_-(data.dims[1]-1)) = data;

    return ret;
}

void print_help();

int phypp_main(int argc, char* argv[]) {
    if (argc <= 1) {
        print_help();
        return 0;
    }

    std::string indir = argv[1];

    if (indir == "help") {
        print_help();
        return 0;
    }

    indir = file::directorize(indir);
    std::string odir = file::directorize(argv[2]);
    file::mkdir(odir);

    vec1s files = file::list_files(indir+"*-nir.fits");
    inplace_sort(files); // make sure they are sorted by date

    // Read wavelength map
    vec2d wave = from_pre(fits::read(argv[3]));

    // Read frames
    vec3f cube;
    vec1f offsets;
    for (uint_t i : range(files)) {
        file::copy(indir+files[i], odir+files[i]);

        fits::input_image iimg(odir+files[i]);
        vec2f tmp;
        iimg.read(tmp);

        double yoff;
        iimg.read_keyword("HIERARCH ESO SEQ CUMOFF Y", yoff);
        offsets.push_back(yoff);

        if (cube.empty()) {
            cube.resize(files.size(), tmp.dims);
        }

        cube(i,_,_) = tmp;
    }

    // Subtract mean offset to make sure we have some + and some -
    offsets -= mean(offsets);

    if (files.size() % 2 != 0 || count(offsets > 0) % 2 != 0) {
        error("needs an even number of frames with +/- Y to do NOD subtraction");
        return 1;
    }

    for (uint_t i : range(files.size()/2)) {
        uint_t ia = 2*i+0, ib = 2*i+1;

        if (offsets[ia]*offsets[0] < 0) {
            // Handle ABBA sequence
            std::swap(ia, ib);
        }

        print(odir+files[ia]);

        if (offsets[ia]*offsets[ib] > 0) {
            error("pair is not a NOD pair");
            return 1;
        }

        std::string qc_dir = odir+"/dark_qc/"+file::get_basename(files[ia])+"/";
        file::mkdir(qc_dir);

        vec2f a = cube(ia,_,_);
        vec2f b = cube(ib,_,_);

        // Subtract
        vec2f sum = a - b;

        fits::write(qc_dir+"sub_orig.fits", sum);

        // First mask the slits
        vec2b mask = make_mask1(wave);
        fits::write(qc_dir+"mask1.fits", mask);

        // Do a first subtraction
        vec2f cor(sum.dims);
        subtract_pattern(sum, cor, mask);

        // Now subtract sky
        vec2d sky = make_sky(sum, wave, mask);
        fits::write(qc_dir+"sky.fits", sky);
        sum -= sky;

        // Mask strong residuals
        mask = make_mask2(sum);
        fits::write(qc_dir+"mask2.fits", mask);

        // Second subtraction
        subtract_pattern(sum, cor, mask);
        sum += sky;

        fits::write(qc_dir+"sub_cor.fits", sum);
        fits::write(qc_dir+"cor.fits", cor);

        // Subtract pattern from first frame
        a -= cor;

        // Save (only A, leave B untouched)
        fits::image oimg(odir+files[ia]);
        oimg.update(a);
    }

    return 0;
}

void print_help() {
    using namespace terminal_format;

    print("correct_dark v1.0 (xshooter)");
    print("usage: correct_dark <input_dir> <output_dir> <wavelength_solution_file>");
    print("");
    print("Main parameters:");
    paragraph("The first argument must be a backup directory containing X-SHOOTER raw data for one "
        "OB. The program will analyze the data in this directory, and save the dark-corrected "
        "frames in the directory specified as second argument. The third argument must be a "
        "wavelength calibration file produced by the X-SHOOTER pipeline for this OB "
        "(WAVE_MAP_NIR.fits).");
    paragraph("WARNING: the program assumes exposures will be subtracted in pairs in the final "
        "reduction (ABBA or ABAB pattern), and will not support STARE reductions.");
}
