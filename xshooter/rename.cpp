#include <phypp.hpp>

void print_help();

int phypp_main(int argc, char* argv[]) {
    if (argc < 2) {
        print_help();
        return 0;
    }

    std::string dir = file::directorize(argv[1]);

    vec1s files = file::list_files(dir+"XSHOO*.fits");
    print(files.size(), " FITS files to analyze");

    for (auto& file : files) {
        fits::header hdr = fits::read_header_hdu(dir+file, 0);

        std::string obj;
        if (!fits::getkey(hdr, "OBJECT", obj)) {
            error("missing OBJECT keyword in '", file, "'");
            continue;
        }

        obj = trim(to_lower(replace(obj, ",", "-")));

        if (obj == "object") {
            obj = "acq";
        } else if (obj == "dark" || obj == "bias" || obj == "lamp-fmtchk" || obj == "lamp-orderdef" ||
            obj == "lamp-flat" || obj == "lamp-dflat" || obj == "lamp-qflat" ||
            obj == "lamp-dorderdef" || obj == "lamp-qorderdef" || obj == "lamp-wave" ||
            obj == "lamp-afc" || obj == "std-flux" || obj == "std-telluric") {
            // nothing to do
        } else {
            // This must be a science frame
            // Make sure
            std::string tpl;
            if (!fits::getkey(hdr, "HIERARCH ESO TPL ID", tpl)) {
                error("missing HIERARCH ESO TPL ID keyword in '", file, "'");
                continue;
            }

            tpl = trim(to_lower(tpl));

            if (!begins_with(tpl, "xshooter_slt_obs")) {
                warning("unknown frame type '", obj, "' with template '", tpl, "'");
                continue;
            }

            obj = "sci";
        }

        std::string arm;
        if (!fits::getkey(hdr, "HIERARCH ESO SEQ ARM", arm)) {
            error("missing HIERARCH ESO ESO SEQ keyword in '", file, "'");
            continue;
        }

        arm = trim(to_lower(arm));

        if (arm != "agc") {
            obj = obj+"-"+arm;
        }

        if (!ends_with(file, obj+".fits")) {
            spawn("mv "+dir+file+" "+dir+file::remove_extension(file)+"-"+obj+".fits");
        }
    }

    return 0;
}

void print_help() {
    using namespace terminal_format;

    print("rename v1.0 (xshooter)");
    print("usage: rename <directory>");
    print("");
    print("Main parameters:");
    paragraph("'directory' must contain an abitrary number of X-SHOOTER raw data frames in "
        "FITS format. This program will rename these files one by one, appending a suffix "
        "to the file's name to allow you to identify the file's purpose.");
}
