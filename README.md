# Misfit evaluator for the DT-GEO WP6 tsunami workflow. 

As of now, the misfit evaluator only includes a Python script for the Samos test case as well as a src directory that includes class definitions for the misfit evaluator measurement types (e.g. tide gauges, GNSS)

# Instructions for this misfit evaluator

This script can be run from the command line and needs some additional arguments.

How to run the script: 
  python misfit_samos.py --data_path data_path (--arrival_time_percentage percentage --plot_arrivaltimes boolean --plot_arrivaltimes_scenario integer --remove_gauges boolean)

Arguments that need/can to be provided:
  * --data_path data_path                               Path to the directory where both the scenario data and the gauge data is stored. 
  * --arrival_time_percentage percentage      Percentage of the maximum wave height at which the arrival time is picked. Default is 10% (=0.1).
  * --plot_arrivaltimes boolean                       Optional handle to trigger plotting of the gauge data with picked arrival times. True/False; default: False
  * --plot_arrivaltimes_scenario integer         Optional handle to trigger plotting the comparison of gauge data with scenario data for this scenario index with arrival times. Integer; default: None
  * --remove_gauges boolean                        Optional handle to remove badly located gauges relative to closest POIs. Removes kos1, kos2, syro and NOA03 tide gauges. True/False; default: False


Note: The script assumes that the data is stored in a folder (= the data path given by the input 'data_path').
The scenario data should be located in a directory called 'data' (where each scenario has an additional individual directory) within 'data_path'.
The gauge data should be located in a directory called 'gauges' within 'data_path'.

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
