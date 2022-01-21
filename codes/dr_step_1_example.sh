for f in /apollo/jung/bandersen/noaa_class_dr_2020_10_10
do
 newName=${f/GCRSO-*\_/GCRSO_}
 mv -i "$f" "$newName"
done

for filename in /apollo/jung/bandersen/noaa_class_dr_2020_10_10/*.h5; do
 [ -f "$filename" ] || continue
 cp "$filename" "${filename//GCRSO_/}"
done

for file in /apollo/jung/bandersen/noaa_class_dr_2020_10_10/_j01_*.h5
 do
  mv "$file" "SCRIF$file"
 done

mkdir /apollo/jung/bandersen/raw_dr_output_2020_10_10
cd /apollo/jung/bandersen/raw_dr_output_2020_10_10

export HS_RET_DIR=/apollo/jung/bandersen/CSPP_UW_HSRTV_2_0
source $HS_RET_DIR/env/uw_hs_l2.bash_env
run_HSRTV.scr 3 /apollo/jung/bandersen/noaa_class_dr_2020_10_10


for f in /apollo/jung/bandersen/full_dr_output_2020_10_10/CrIS_FSR_j01_*19.atm_prof_rtv.h5
do
 mv -r "${f}" "${f/19./.}"
done
