"""
Misfit evaluator for the Samos test case 

Based on the GITEWS approach. See the misfit draft for more detail.

*** Instructions for this module ***

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

"""

# Generic modules that are needed 
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
from src.TypeMeasurement import TypeMeasurement 
from netCDF4 import Dataset
import time

# Define arguments the script has to be called with
parser = argparse.ArgumentParser(
    description="Misfit evaluator for the Samos test case"
)
parser.add_argument(
    "--data_path",
    help="Path to the directory where both the scenario data and the gauge data is stored.",
    required=True,
)
parser.add_argument(
    "--arrival_time_percentage",
    help="Percentage of the maximum wave height at which the arrival time is picked. Default is 10% (=0.1).",
    default=0.1,
)
parser.add_argument("--plot_arrivaltimes", 
    help="Optional handle to trigger plotting of the gauge data with picked arrival times. True/False; default: False", 
    default=False)
parser.add_argument("--plot_arrivaltimes_scenario", 
    help="Optional handle to trigger plotting the comparison of gauge data with scenario data for this scenario index with arrival times. Integer; default: None", 
    default=None)
parser.add_argument("--remove_gauges", 
    help="Optional handle to remove badly located gauges relative to closest POIs. Removes kos1, kos2, syro and NOA03 tide gauges. True/False; default: False",
    default=False)

args = parser.parse_args()
scenario_data_path = args.data_path
plot_arrivaltimes = args.plot_arrivaltimes
remove_gauges = args.remove_gauges
arrtime_percentage = float(args.arrival_time_percentage)
plot_scenario = args.plot_arrivaltimes_scenario
if (not(plot_scenario) == None):
  plot_arrivaltimes = True
  plot_scenario = int(plot_scenario)

# Some functions that allow for better readability
def get_scenario_waves(N, data_path, data_folders, indices):
  """
  Function to read the scenario data and output them in an array
  """
  wave_data = []
  for scenario in range(N):
    ncfile = os.path.join(data_path, data_folders[scenario], 'out_ts.nc')
    ds = Dataset(ncfile, 'r', format='NETCDF4')
    scenario_wave_amplitude = ds.variables['eta'][:]
    scenario_min_height = ds.variables['min_height'][indices]
    scenario_max_height = ds.variables['max_height'][indices]
    scenario_POIs_wave = scenario_wave_amplitude[:, indices]
    wave_data.append(scenario_POIs_wave)
    
  scenario_time = ds.variables['time'][:]
  return wave_data, scenario_time, scenario_min_height, scenario_max_height

def get_PTF_waveheights(ngauges, gauge_data):
  """
  Function to save the minimum and maximum PTF waveheight.
  """
  minimum_waveheight = np.zeros(ngauges)
  maximum_waveheight = np.zeros(ngauges)
  for gauge in range(ngauges):
    minimum_waveheight[gauge] = np.min(gauge_data[gauge])
    maximum_waveheight[gauge] = np.max(gauge_data[gauge])
  return minimum_waveheight, maximum_waveheight
  
def get_PTF_time_indices(ngauges, gauge_times, scenario_time):
  """
  Function to compare gauge times and choose the closest PTF scenario time.
  """
  
  PTF_indices = []
  for gauge in range(ngauges):
    time_array = np.array(gauge_times[gauge])
    sub_indices = np.zeros(len(time_array), dtype=int)
    for index in range(len(time_array)):
      time_value = time_array[index]
      time_diff = np.abs(scenario_time - time_value)
      sub_indices[index] = np.argmin(time_diff)
    PTF_indices.append(sub_indices)
  return PTF_indices

#------------------------------------------------------------------------------------------------------------------------------------------------
"""
This part loads the points of interest (POIs) and gives context about the gauge stations.

The following stations are given (see SamosMap as well, for reference):

1: gokce            30 s sampling
2: plom             60 s sampling
3: bodru            30 s sampling
4: kos1 and kos2 (they have the same coordinates)   60 s sampling
5: marma          30 s sampling
6: NOA03          60 s sampling
7: NOA04          60 s sampling
8: NOA10          60 s sampling
9: syro              60 s sampling

Some of these POIs cane be removed since the gauge is located badly relative to the POI (4, 6 and 9)

POIs closest to the gauges (selected from with another script)
gauge_POI = [507, 119, 147, 146, 157, 630, 650, 641, 547]

with the following longitude coordinates:
lon_tg = [25.8936, 26.317, 27.424, 27.333, 28.3848, 26.92, 25.743, 25.15, 24.941] 
and the latitude coordinates:
lat_tg = [40.2314, 38.977, 37.032, 36.906, 36.8380, 35.42, 35.009, 35.35, 37.438] 

# POIs can be loaded from the pois.npy file
# The POIs contain the following 'pois_coords', 'pois_labels' and 'pois_index' 
pois = np.load(os.path.join(scenario_data_path,"pois.npy"),allow_pickle=True).item()
[:]
# Get coordinates and indices
POIs_coordinates = pois["pois_coords"]
POIs_indices = pois["pois_index"]

"""

# These are the index sorted gauges 
gauge_POI = [119, 146, 146, 147, 157, 507, 547, 630, 641, 650]

ngauges = len(gauge_POI)

gauge_list = ['plom_res.dat', 'kos1_res.dat', 'kos2_res.dat', 'bodru.rad.rmn.txt', 
                      'marma.rad.rmn.txt', 'gokce.rad.rmn.txt', 'syro.pr1.rmn.txt', 
                      'NOA03.rad.rmn.txt', 'NOA10.rad.rmn.txt', 'NOA04.rad.rmn.txt']
# Sampling rates for the gauges. Note: HySEA sampling rate is about 30 seconds 
sampling = [60., 60., 60., 30., 30., 30., 60., 60., 60., 60.]

# Load reference data (= gauge data)
gaugedir = os.path.join(scenario_data_path, 'gauges')

time_list = []
data_list = []
for gauge in range(ngauges):
  time_list.append(np.loadtxt(os.path.join(gaugedir, gauge_list[gauge]))[:,0])
  data_list.append(np.loadtxt(os.path.join(gaugedir, gauge_list[gauge]))[:,1])
  

#------------------------------------------------------------------------------------------------------------------------------------------------
"""
Fetch PTF scenario data
"""

# Get data folder names
data_path = os.path.join(scenario_data_path, 'data')
data_folders = []
data_folders += [each for each in sorted(os.listdir(data_path))]
N = len(data_folders)

# Read and save the scenario results at the gauge POIs and time data
print("Reading scenario data.")
start = time.time()
scenario_results, scenario_time, scenario_min_height, scenario_max_height = get_scenario_waves(N, data_path, data_folders, gauge_POI)
stop = time.time()    
print(f"Reading the data took {stop - start} s.\n")

#------------------------------------------------------------------------------------------------------------------------------------------------
"""
Pick arrival times (they have to be used in the time range picking)
"""

# Get time range and set up time arrays for the gauge data
time_range = [np.min(scenario_time), np.max(scenario_time)]
N_scenario_time = len(scenario_time)

# Overwrite gauge times and data so that they match the scenario time range
for gauge in range(ngauges):
  data_list[gauge] = data_list[gauge][(time_list[gauge] >= time_range[0]) &
                                                          (time_list[gauge] <= time_range[1])]
  time_list[gauge] = time_list[gauge][(time_list[gauge] >= time_range[0]) &
                                                          (time_list[gauge] <= time_range[1])]
                                                          
# Get PTF min/max waveheight
min_waveheight, max_waveheight = get_PTF_waveheights(ngauges, data_list)

# Arrival times for the gauges
arrival_times = np.zeros(ngauges)
for gauge in range(ngauges):
  absmax_waveheight = max(np.abs(min_waveheight[gauge]), np.abs(max_waveheight[gauge]))
  arrtime_trigger = absmax_waveheight * arrtime_percentage
  tmp_data = data_list[gauge][(time_list[gauge] >= time_range[0])]
  arrtime_index = np.argmax(np.array(np.abs(tmp_data)) >= arrtime_trigger)
  arrival_times[gauge] = time_list[gauge][arrtime_index]

# Arrival times for the scenarios
scenario_arrival_times = np.zeros((N, ngauges))
for scenario in range(N):
  for gauge in range(ngauges):
    current_wave = scenario_results[scenario][:,gauge]
    scenario_absmaxwaveheight_calc = np.max(np.abs(current_wave))
    scenario_absmaxwaveheight = max(scenario_absmaxwaveheight_calc,
                                                      max(np.abs(scenario_min_height[gauge]),
                                                              np.abs(scenario_max_height[gauge])))
    arrtime_trigger = scenario_absmaxwaveheight * arrtime_percentage
    arrtime_trigger = scenario_absmaxwaveheight_calc * arrtime_percentage
    arrtime_index = np.argmax(np.array(np.abs(current_wave)) >= arrtime_trigger)
    scenario_arrival_times[scenario, gauge] = scenario_time[arrtime_index]

stop = time.time()    
print(f"Calculating the arrival times took {stop - start} s.\n")
#------------------------------------------------------------------------------------------------------------------------------------------------

# Cut off gauge times and data that lie outside the arrival time range 
gauge_cut_times = []
gauge_cut_data = []
for gauge in range(ngauges):
  gauge_tmp_times = time_list[gauge][time_list[gauge] >= arrival_times[gauge]]
  gauge_tmp_data = data_list[gauge][time_list[gauge] >= arrival_times[gauge]]  
  gauge_cut_times.append(gauge_tmp_times)
  gauge_cut_data.append(gauge_tmp_data)

# Get indices which compare the gauge time with the closest PTF scenario time
gauge_indices = get_PTF_time_indices(ngauges, gauge_cut_times, scenario_time)

# Cut off scenario data that lies outside the arrival time range 
scenario_results_cut = []
for scenario in range(N):
  for gauge in range(ngauges):
    scenario_results_tmp = scenario_results[scenario][scenario_time >= scenario_arrival_times[scenario, gauge], gauge]
    scenario_results_cut.append(scenario_results_tmp)

# Get maximum index for wave comparison (gauge and scenario data have window that needs to be matched for comparison)
PTF_maxindex = np.zeros((N, ngauges),dtype=int)
for gauge in range(ngauges):
  max_index_length = len(gauge_cut_data[gauge])
  gauge_index = gauge_indices[gauge]
  # determine index increment for gauge (is either 1 or 2; gauge_data always has increment 1, so the gauge index will be stored with an increment of 1)
  index_dt = gauge_index[1]-gauge_index[0]  
  len_gauge = len(gauge_index)    
  for scenario in range(N):
    # get lenght of cut scenario results and check if length has to be changed due to index increment
    len_scenario_old = len(scenario_results_cut[gauge + (ngauges-1)*scenario])
    new_scenario_index = np.arange(0, len_scenario_old, index_dt, dtype=int)
    len_scenario = len(new_scenario_index)
    
    # calculate highest possible index = lowest length of timeseries
    min_len = min(len_scenario, max_index_length, len_gauge)
    PTF_maxindex[scenario, gauge] = min_len
    #print(min_len)
    
#------------------------------------------------------------------------------------------------------------------------------------------------
"""
Gauge setup
"""

# Set up GAUGES type 
coords_gauges = [[26.317, 38.977], # plom
                             [27.333, 36.906], # kos1/2
                             [27.424, 37.032], # bodru
                             [28.3848, 36.8380], # marma
                             [25.8936, 40.2314], # gokce
                             [24.941, 37.438], # syro
                             [26.92, 35.42], # NOA03
                             [25.15, 35.35], # NOA10
                             [25.743, 35.009]] # NOA04
weigths_gauges = np.ones(ngauges) / ngauges
range_gauges = [np.min(min_waveheight), np.max(max_waveheight)]

GAUGES = TypeMeasurement(ngauges, coords_gauges, weigths_gauges, range_gauges, measuretypeweight = 0.5)
GAUGES.scale_measured_data()

if (remove_gauges): # This removes gauges 1, 2, 6 and 7 (kos1, kos2, syro and NOA03) which may not be useful due to their position relative to the POIs
  weigths_gauges[[1,2,6,7]] = 0.
  GAUGES.renormalize_stationweights()

#------------------------------------------------------------------------------------------------------------------------------------------------
"""
Arrival times setup
"""
print("Calculating arrival times.\n")
start = time.time()

# Set up ARRTIME (= arrival times) type at each gauge (will use the same coordinates)
# The arrival time will be based on the first occurrence of a certain percentage of the maximum wave height
weights_arrtime = np.ones(ngauges) / ngauges
range_arrtime = time_range

ARRTIME = TypeMeasurement(ngauges, coords_gauges, weights_arrtime, range_arrtime, measuretypeweight = 0.5)
ARRTIME.scale_measured_data()

if (remove_gauges): # This removes gauges 1, 2, 6 and 7 (kos1, kos2, syro and NOA03) which may not be useful due to their position relative to the POIs
  weights_arrtime[[1,2,6,7]] = 0.
  ARRTIME.renormalize_stationweights()




#------------------------------------------------------------------------------------------------------------------------------------------------
"""
Plotting arrival times (and gauge data)
"""

if (plot_arrivaltimes):
  for gauge in range(ngauges):
    if (plot_scenario is None):
      plt.figure()
      plt.plot(gauge_times[gauge], gauge_data[gauge])
      plt.axvline(x = arrival_times[gauge], color = 'b')
      plt.title(f"Gauge data and arrival time (blue line) for {gauge_list[gauge]}")
      plt.xlabel('Time after earthquake = model time [s]')
      plt.ylabel('Wave height [m]')
    else:
      fig, ax = plt.subplots(2)
      fig.suptitle(f"Gauge data vs scenario and arrival times (blue/red line) for {gauge_list[gauge]}")
      ax[0].plot(gauge_times[gauge], gauge_data[gauge])
      ax[0].axvline(x = arrival_times[gauge], color = 'b', label='gauge')
      ax[0].axvline(x = scenario_arrival_times[plot_scenario,gauge], color = 'r', label='scenario')
      ax[1].plot(scenario_time, scenario_results[plot_scenario][:,gauge])
      plt.ylabel('Wave height [m]')
      ax[1].axvline(x = arrival_times[gauge], color = 'b', label='gauge')
      ax[1].axvline(x = scenario_arrival_times[plot_scenario,gauge], color = 'r', label='scenario')
      ax[0].legend()
      ax[1].legend()
      plt.ylabel('Wave height [m]')
      plt.xlabel('Time after earthquake = model time [s]')
  plt.show()

#------------------------------------------------------------------------------------------------------------------------------------------------
"""
Calculate distances for gauge data
"""

print("Calculating gauge distances and norms.")
start = time.time()

# This loop structure is faster than calulating the scaled distances/norms in one big loop
gauges_dist = []
for gauge in range(ngauges):
  for scenario in range(N):
    max_index = PTF_maxindex[scenario, gauge]
    current_scenario_data = scenario_results_cut[scenario][0:max_index][gauge]
    current_gauge_data = gauge_cut_data[gauge][0:max_index]
    # Each entry is a len(indices) array of distances; Results will be normalized
    gauges_dist.append( np.abs(current_gauge_data - current_scenario_data))
  
gauges_norms = np.zeros((N, ngauges))
idx = 0
for gauge in range(ngauges):
  for scenario in range(N):
    current_dist_gauge = gauges_dist[idx]
    gauges_norms[scenario, gauge] = GAUGES.scaling_func(np.linalg.norm(current_dist_gauge))
    idx += 1
    
stop = time.time()    
print(f"Calculating gauge distances and norms took {stop - start} s.\n")

#------------------------------------------------------------------------------------------------------------------------------------------------
"""
Calculate distances for arrival times
"""

print("Calculating arrival time distances and norms.")
start = time.time()

arrtime_dist_scaled = np.zeros((N, ngauges))
for scenario in range(N):
  for gauge in range(ngauges):
    arrtime_dist = np.abs(scenario_arrival_times[scenario, gauge] - arrival_times[gauge])
    arrtime_dist_scaled[scenario, gauge] = ARRTIME.scaling_func(arrtime_dist)
  
stop = time.time()    
print(f"Calculating arrival time distances and norms took {stop - start} s.\n")

#------------------------------------------------------------------------------------------------------------------------------------------------
"""
Calculate misfit
"""

print("Calculating misfit.")
start = time.time()

# Calculate misfit for each type (= sum of the norms and each stationweight)
GAUGES_misfit = np.dot(gauges_norms, GAUGES.stationweights)

ARRTIME_misfit = np.dot(arrtime_dist_scaled, ARRTIME.stationweights)

# Calculate total misfit (sum of misfits for each type and the typeweights)
MISFIT = GAUGES_misfit * GAUGES.measuretypeweight + \
               ARRTIME_misfit * ARRTIME.measuretypeweight

stop = time.time()    
print(f"Calculating the misfit took {stop - start} s.\n")

#------------------------------------------------------------------------------------------------------------------------------------------------
"""
Plot misfit
"""

fig = plt.figure()
fig.subplots_adjust(top=0.88)
if (remove_gauges):
  titlestring = "some gauges removed"
else:
  titlestring = "all gauges"
  
fig.suptitle(f"Misfit for the Samos test case with {titlestring}" )
plt.subplot(2,2,1)
plt.plot(ARRTIME_misfit, '.')
plt.title('Arrival time misfit')
plt.xlabel('Scenario number')
plt.ylabel('Normalized arrival time misfit')

plt.subplot(2,2,2)
plt.plot(GAUGES_misfit, '.')
plt.title('Gauge misfit')
plt.xlabel('Scenario number')
plt.ylabel('Normalized gauge misfit')

plt.subplot(2,2,3)
plt.plot(MISFIT, '.')
plt.title('Total misfit')
plt.xlabel('Scenario number')
plt.ylabel('Normalized total misfit')
#fig.set_size_inches(18.5, 10.5)
plt.tight_layout()
plt.show()
