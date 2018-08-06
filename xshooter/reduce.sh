# ---------------------------
# join() function
#
# example:
#    LIST=(a b c d)
#    `join , ${LIST[@]}`
# outputs: "a,b,c,d"
# ---------------------------

join() {
    # $1 is separator (example: ,)
    # $2... are the elements to join
    local sep=${1} IFS=
    local join_ret=${2}
    shift 2 || shift $(($#))
    join_ret+="${*/#/$sep}"
    echo ${join_ret}
}

# MODIFY the path to the directory containing the speclib tools
SPECLIB_DIR=../../xshooter-scripts

# MODIFY this to reduce only one arm instead of the three
ARMS=(UVB VIS NIR)

# MODIFY here the default seeing (FWHM in arcsec) for each ARM
DEFAULT_SEEING=(0.5 0.5 0.5)

# MODIFY this to match the OBs you have, or to reduce only one OB
OBS=(A B C D E F G H I)

# MODIFY this if your data is at a different location
DATA_DIR=raw

# MODIFY this to save the 1D spectra in a different directory
OUT_DIR=reduced

# MODIFY the name of your target
# (unimportant unless you're dealing with several targets in the same program)
TARGET=7329

# MODIFY the list of sources to extract for this target.
# A "target" is one slit pointing, and a "source" is an actual galaxy seen in this slit.
# Therefore there can be more than one "source" for a given "target".
#  - SOURCES: name of the sources (arbitrary, you choose how you name them).
#  - SOURCES_ID: IDs of these sources in the photometry catalog. This will be used for flux
#    rescaling later.
#  - SOURCES_POS: the position offsets of the sources (in arcsec) along the slit.
#  - SOURCES_WIDTH: the radii of the sources (sigma of Gaussian in arcsec), before seeing.
SOURCES=(blobA blobB)
SOURCES_ID=(7329 7365)
SOURCES_POS=(-1.16 1.16)
SOURCES_WIDTH=(0.133 0.133)

# MODIFY here the options passed to 'extract2d' (see 'extract2d help')
EXTRACT_OPTIONS="gauss_nofit"

# MODIFY the rebinning to use to produce stacked spectra.
# (you can list as many as you like, but be sure to include "1" to produce an unbinned spectrum)
REBIN_STACKS=[1,30,100]

# MODIFY the rebinning to use to produce stacked 2D spectra (can only set one value)
# The programs will also produce an unbinned version automatically.
REBIN_STACKS2D=30

# MODIFY here the path to the filter database
FILTER_DB=${SPECLIB_DIR}/filter-db/db.dat

# MODIFY the list of broadband filters to use for rescaling to photometry
# Must be listed in the FILTER_DB.
RESCALING_FILTERS=[ukidss_H]

# MODIFY the path to the photometric catalog containing the fluxes
RESCALING_CATALOG=/home/cschreib/data/fits/uds/uds_zfourge_dr1.fits

# MODIFY the list of filters for which to produce stacked synthetic photometry (for checks).
# Must be listed in the FILTER_DB.
STACK_FILTERS=[ukidss_H,fourstar_H1,fourstar_H2,f160w]

# MODIFY enable/disable components of the reduction
DO_TELLURIC=0 # estimate telluric correction
DO_EXTRACT=1  # extract 1D spectra
DO_RESCALE=0  # rescale 1D spectra to photometry
DO_STACK=1    # stack 1D spectra

# Disable telluric correction in UVB arm?
# (should be disabled, as standard stars are typically not observed with UVB arm so you
# cannot compute the telluric correction anyway)
NO_TELLURIC_UVB=1

for ((a=0; a<${#ARMS[@]}; a++)); do
    ARM=${ARMS[$a]}

    # 1) Extract 1D spectra for each OB
    for ((o=0; o<${#OBS[@]}; o++)); do
        OB=${OBS[$o]}

        # 1.1) Compute telluric correction
        if [ ${DO_TELLURIC} -eq 1 ]; then
            if [ ${ARM} != "UVB" ] || [ ${NO_TELLURIC_UVB} -eq 0 ]; then
                echo "note: processing telluric star for OB ${OB} and arm ${ARM}"
                ophy++ ${SPECLIB_DIR}/xshooter/get_telluric \
                    flux_file=${DATA_DIR}/${OB}/telluric/TELL_SLIT_FLUX_MERGE1D_${ARM}.fits \
                    flux_hdu=0 error_hdu=1 out=${OUT_DIR}/${OB}/telluric_${ARM}.fits \
                    telluric_gaps=${SPECLIB_DIR}/xshooter/telluric_regions.dat
            fi
        fi

        # 1.2) Extract spectra and apply telluric correction
        if [ ${DO_EXTRACT} -eq 1 ]; then
            EPOCHS=(`ls ${DATA_DIR}/${OB}/ | grep -E "${ARM}$"`)
            echo "note: found ${#EPOCHS[@]} epochs for OB ${OB}"

            OFFSETS=()
            SEEINGS=()

            # MODIFY. Here you can specify the offsets and seeing for each OB & ARM (separately).
            # This is optional. If you don't do anything, no offset will be used, and the seeing
            # will be set to the value of DEFAULT_SEEING.
            # For example:
            if [ ${OB} == "A" ] && [ ${ARM} == "NIR" ]; then
                # Set seeing in OB A and NIR arm to 0.7 arcsec, and an offset of -1.2 pixels
                SEEINGS=(0.7)
                OFFSETS=(-1.2)
            fi

            for ((e=0; e<${#EPOCHS[@]}; e++)); do
                EPOCH=${EPOCHS[$e]}

                if [ ${#OFFSETS[@]} -eq 0 ]; then
                    OFFSET=0
                elif [ ${#OFFSETS[@]} -eq 1 ]; then
                    OFFSET=${OFFSETS[0]}
                else
                    OFFSET=${OFFSETS[$e]}
                fi

                if [ ${#SEEINGS[@]} -eq 0 ]; then
                    SEEING=${DEFAULT_SEEING[$a]}
                elif [ ${#SEEINGS[@]} -eq 1 ]; then
                    SEEING=${SEEINGS[0]}
                else
                    SEEING=${SEEINGS[$e]}
                fi

                echo "note: analyzing epoch ${EPOCH}"

                EXTRA_OPTS=""
                if [ ${ARM} != "UVB" ] || [ ${NO_TELLURIC_UVB} -eq 0 ]; then
                    EXTRA_OPTS="${EXTRA_OPTS} telluric=${OUT_DIR}/${OB}/telluric_${ARM}.fits"
                fi

                if [ ${ARM} == "UVB" ]; then
                    EXTRA_OPTS="${EXTRA_OPTS} flag_begin=0.314"
                fi
                if [ ${ARM} == "VIS" ]; then
                    EXTRA_OPTS="${EXTRA_OPTS} flag_end=1.01 flag_begin=0.56"
                fi
                if [ ${ARM} == "NIR" ]; then
                    EXTRA_OPTS="${EXTRA_OPTS} flag_end=1.91 flag_begin=1.036"
                fi

                ophy++ ${SPECLIB_DIR}/extract2d \
                    flux_file=${DATA_DIR}/${OB}/${EPOCH}/SCI_SLIT_FLUX_MERGE2D_${ARM}.fits \
                    flux_hdu=0 error_hdu=1 name=${TARGET} out_dir=${OUT_DIR}/${OB} suffix=_${EPOCH} \
                    dpix=37 dpixfit=12 oh_line_width=60 verbose auto_mask auto_mask_dl=30 \
                    offset=${OFFSET} default_gauss_width=[`join , ${SOURCES_WIDTH[@]}`] \
                    default_gauss_pos=[`join , ${SOURCES_POS[@]}`] \
                    sources=[`join , ${SOURCES[@]}`] instrument="X-SHOOTER" \
                    ${EXTRA_OPTS} discard_wcs_pos ${EXTRACT_OPTIONS}
            done
        fi
    done

    # 2) Compute rescaling to total flux for each OB and stack
    if [ ${DO_STACK} -eq 1 ]; then
        if [ ${DO_RESCALE} -eq 1 ]; then
            # 2.1) Compute rescaling to total flux
            for ((i=0; i<${#SOURCES[@]}; i++)); do
                SOURCE=${SOURCES[$i]}
                CATID=${SOURCES_ID[$i]}

                echo "note: computing broadband fluxes for ${SOURCE} ${ARM}"

                FILES=()
                for ((o=0; o<${#OBS[@]}; o++)); do
                    OB=${OBS[$o]}
                    EPOCHS=(`ls ${DATA_DIR}/${OB}/ | grep -E "${ARM}$"`)
                    for ((e=0; e<${#EPOCHS[@]}; e++)); do
                        EPOCH=${EPOCHS[$e]}
                        FILES=("${FILES[@]}" ${OUT_DIR}/${OB}/stacked_${TARGET}_${EPOCH}_${SOURCE}.fits)
                    done
                done

                echo "note: found ${#FILES[@]} exposures"

                ophy++ ${SPECLIB_DIR}/get_fluxes files=[`join , ${FILES[@]}`] \
                    catalog=${RESCALING_CATALOG} catalog_id=${CATID} \
                    verbose filters=${RESCALING_FILTERS}
            done

            # 2.2) Rescale
            for ((o=0; o<${#OBS[@]}; o++)); do
                OB=${OBS[$o]}
                EPOCHS=(`ls ${DATA_DIR}/${OB}/ | grep -E "${ARM}$"`)
                for ((e=0; e<${#EPOCHS[@]}; e++)); do
                    EPOCH=${EPOCHS[$e]}
                    echo "note: rescaling flux for ${EPOCH}"

                    FILES=()
                    for ((i=0; i<${#SOURCES[@]}; i++)); do
                        SOURCE=${SOURCES[$i]}
                        FILES=("${FILES[@]}" ${OUT_DIR}/${OB}/stacked_${TARGET}_${EPOCH}_${SOURCE}.fits)
                    done

                    ophy++ ${SPECLIB_DIR}/flux_rescale files=[`join , ${FILES[@]}`] verbose
                done
            done
        fi

        # 2.3) Stack
        for ((i=0; i<${#SOURCES[@]}; i++)); do
            SOURCE=${SOURCES[$i]}
            CATID=${SOURCES_ID[$i]}

            echo "note: stacking ${SOURCE} ${ARM}"

            FILES=()
            FILESR=()
            FILES2DR=()
            FILES2D=()
            for ((o=0; o<${#OBS[@]}; o++)); do
                OB=${OBS[$o]}
                EPOCHS=(`ls ${DATA_DIR}/${OB}/ | grep -E "${ARM}$"`)
                for ((e=0; e<${#EPOCHS[@]}; e++)); do
                    EPOCH=${EPOCHS[$e]}
                    FILES=("${FILES[@]}" ${OUT_DIR}/${OB}/stacked_${TARGET}_${EPOCH}_${SOURCE}.fits)
                    if [ ${DO_RESCALE} -eq 1 ]; then
                        FILESR=("${FILESR[@]}" ${OUT_DIR}/${OB}/stacked_${TARGET}_${EPOCH}_${SOURCE}_rescaled.fits)
                    else
                        FILESR=("${FILESR[@]}" ${OUT_DIR}/${OB}/stacked_${TARGET}_${EPOCH}_${SOURCE}.fits)
                    fi
                    FILES2DR=("${FILES2DR[@]}" ${OUT_DIR}/${OB}/stacked_${TARGET}_${EPOCH}_${SOURCE}_spec2d_residual.fits)
                    FILES2D=("${FILES2D[@]}" ${OUT_DIR}/${OB}/stacked_${TARGET}_${EPOCH}_spec2d_bgsub.fits)
                done
            done

            echo "note: found ${#FILES[@]} exposures"

            ophy++ ${SPECLIB_DIR}/stack1d files=[`join , ${FILES[@]}`] \
                out=${OUT_DIR}/${SOURCE}/stacked_${ARM}.fits rebin=${REBIN_STACKS} \
                verbose filters=${STACK_FILTERS}

            if [ ${DO_RESCALE} -eq 1 ]; then
                ophy++ ${SPECLIB_DIR}/stack1d files=[`join , ${FILESR[@]}`] \
                    out=${OUT_DIR}/${SOURCE}/stacked_${ARM}_rescaled.fits rebin=${REBIN_STACKS} \
                    verbose filters=${STACK_FILTERS}
            fi

            ophy++ ${SPECLIB_DIR}/stack2d files=[`join , ${FILES2DR[@]}`] specs=[`join , ${FILESR[@]}`] \
                out=${OUT_DIR}/${SOURCE}/stacked_${ARM}_spec2d_sub.fits rebin=${REBIN_STACKS2D}

            ophy++ ${SPECLIB_DIR}/stack2d files=[`join , ${FILES2D[@]}`] specs=[`join , ${FILESR[@]}`] \
                out=${OUT_DIR}/${SOURCE}/stacked_${ARM}_spec2d.fits rebin=${REBIN_STACKS2D}
        done
    fi
done
