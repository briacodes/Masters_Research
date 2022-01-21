cd /apollo/jung/bandersen/CloudDetection/src
rm -f *.o *.mod ../Cloud_and_Aerosol_Detect
make
cd /home/bandersen/iris-home
mkdir ecmwf_workspace
cd ecmwf_workspace
rm -f cloud_and_aerosol_detection_input.dat
##### Change to the date of interest
ln -s /home/bandersen/iris-home/data/ecmwf_cd_input_20200821_n20_0700-0800.txt cloud_and_aerosol_detection_input.dat
rm -f reference_output.dat
ln -s /home/bandersen/iris-home/data/ecmwf_cd.output reference_output.dat
ln -s /apollo/jung/bandersen/CloudDetection/namelist/CRIS_AERDET.NL .
ln -s /apollo/jung/bandersen/CloudDetection/namelist/CRIS_CLDDET.NL .
/apollo/jung/bandersen/CloudDetection/Cloud_and_Aerosol_Detect
rm -f cloud_and_aerosol_detection_input.dat
rm -f CRIS_AERDET.NL
rm -f CRIS_CLDDET.NL
rm -f reference_output.dat
rm -f result.txt
##### Change to the date of interest
mv cloud_and_aerosol_detection_output.dat /home/bandersen/iris-home/data/ecmwf_cd_output_20200821_0700_0800_n20.txt.newest
cd /home/bandersen/iris-home
rm -d ecmwf_workspace
exit

