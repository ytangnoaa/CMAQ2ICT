"plot-flight-cmaq.py" is the python script for quickview plotting, and "comp-dc8-1.index" is its
index file for plotting. The default data is stored at ../data. To use it, one needs to install
the python packages, such cartopy etc, and type

plot-flight-cmaq.py  $start_flight_index  $end_flight_index $output_file

e.g. plot-flight-cmaq.py  1   21  dc8-1-21.pdf
