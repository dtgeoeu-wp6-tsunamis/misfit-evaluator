 # Misfit evaluator for the DT-GEO WP6 tsunami workflow. 

As of now, the misfit evaluator only includes a Python script for the Samos test case as well as a src directory that includes class definitions for the misfit evaluator measurement types (e.g. tide gauges, GNSS)

This software is still under development and validation (10/2024)

# Instructions for this misfit evaluator

This script can be run from the command line and needs some additional arguments.

How to run the script: 
  python misfit_samos.py --data_path data_path --gauge_path gauge_path (--syntethic_gauge bool --statistical_misfit bool --gauge_POI gauge_POI_list --arrival_time_percentage percentage --plot_arrivaltimes boolean)

Arguments that need/can to be provided:
  * --data_path data_path                               Path to the directory where the scenario data is stored. 
  * --gauge_path gauge_path                         Path to the directory where the gauge data is stored.
  * --syntethic_gauge                                     Handle whether the gauge data is synthetic or real gauge data. True/False; default: False 
  * --statistical_misfit                                     Handle whether the synthetic gauge data will be statistical (misfit analysis at all POIs) or not. True/False; default: False
  * --gauge_POI                                              List of the POIs closest to the gauges. Especially relevant for real gauge data.
  * --arrival_time_percentage percentage      Percentage of the maximum wave height at which the arrival time is picked. Default is 10% (=0.1).
  * --plot_arrivaltimes boolean                       Optional handle to trigger plotting of the gauge data with picked arrival times. True/False; default: False


Note(s): 
The scenario data should be located in a directory where each scenario has an additional individual directory within 'data_path'. The scenario data itself then should be given by a netCDF file by the name of 'out_ts.nc'.

The gauge data should be located in a directory given by 'gauge_path' and can either be synthetic or real gauge data (synthetic_gauge handle). In the case of real gauge data, we assume that ASCII files with times and waveheight are given. The directory should ONLY include the gauge ASCII files!
The gauge POI list should be ordered exactly like the ASCII files (alphabetically, upper case first, then lowercase names).

For the synthetic data that is situated on Leonardo, the gauges (data for the first 23 indices) are used if gauge data misfit is calculated vs. the statistical misfit where all POIs are compared to the syntethic results.

Examples for the script usage are:
Gauge data misfit: 
python misfit_synthetic.py --data_path data_path --gauge_path gauge_path --synthetic_gauge True  --statistical_misfit False --gauge_POI 1980 3801 3800 2115 1755 3337 2202 2199 2200 2210 2193 2172 2172 167 2167 3622 1712 3614 214 2175 2175 2029 0

where --gauge_POI gives the POIs closest to the gauges for comparison.

Statistical misfit:
python misfit_synthetic.py --data_path data_path --gauge_path gauge_path --synthetic_gauge True  --statistical_misft True

# Required Python packages
To run the misfit evaluator, the following Python packages have to be installed on the system. 

Standard Packages:
  * argparse
  * matplotlib
  * numpy
  * os
  * time 
  
Non-standard packages:
  * netCDF4
