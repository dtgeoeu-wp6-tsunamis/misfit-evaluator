"""
Misfit evaluator  

Based on the GITEWS approach. See the misfit draft for more detail.

*** Instructions for this module ***

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
    help="Path to the directory where the scenario data is stored.",
    required=True,
)
parser.add_argument(
    "--gauge_path",
    help="Path to the directory where the gauge data is stored.",
    required=True,
)
parser.add_argument(
    "--synthetic_gauge",
    help="Handle whether the gauge data is synthetic or real gauge data. True/False; default: False",
    default=False,
    type=bool,
)
parser.add_argument(
    "--statistical_misfit",
    help="Handle whether the synthetic gauge data will be statistical (misfit analysis at all POIs) or not. True/False; default: False" ,
    default=False,
    type=bool,
)
parser.add_argument(
    "--gauge_POI",
    nargs='+',
    type=int,
    help="List of the POIs closest to the gauges. Especially relevant for real gauge data.",
    default=None,
)
parser.add_argument(
    "--arrival_time_percentage",
    help="Percentage of the maximum wave height at which the arrival time is picked. Default is 10% (=0.1).",
    default=0.1,
)
parser.add_argument("--plot_arrivaltimes", 
    help="Optional handle to trigger plotting of the gauge data with picked arrival times. True/False; default: False", 
    default=False,
    type=bool,
)

args = parser.parse_args()
scenario_data_path = args.data_path
gauge_data_path = args.gauge_path
synthetic_gauge = args.synthetic_gauge
statistical_misfit = args.statistical_misfit
gauge_POI = args.gauge_POI
plot_arrivaltimes = args.plot_arrivaltimes
arrtime_percentage = float(args.arrival_time_percentage)

#------------------------------------------------------------------------------------------------------------------------------------------------
# Some functions that allow for better readability

def remove_masked(data):
  """
  Function to remove the masked entries from the synthetic data and save the indices that have been removed
  """
  maskedarray_indices = []
  tmp_data = []
  data_len = np.shape(data)[0]
  time_line = np.shape(data)[1]
  for idx in range(data_len):
    if np.ma.is_masked(data[idx]):
      maskedarray_indices.append(idx)
    else: 
      tmp_data.append(data[idx])
  nMask = len(maskedarray_indices)
  cleaned_data = np.zeros((data_len-nMask, time_line))
  for idx in range(data_len-nMask):
    cleaned_data[idx,:] = tmp_data[idx]
  return cleaned_data, maskedarray_indices



def get_gauge_synthetic(gauge_data_path, gauge_POI, statistical_misfit):
  """
  Function to read the syntethic gauge data and output them in an array
  """
  ncfile = os.path.join(gauge_data_path, 'grid_A_ts.nc')
  ds = Dataset(ncfile, 'r', format='NETCDF4')
  gauge_time = ds.variables['time'][:]
 
  # Get time and wave height data
  time_data = ds.variables['time'][:]
  tmp_wave_data = ds.variables['eta'][:]  
 
  # get number of gauges (which is equal to number of grid points for the synthetic data if not otherwise specified)
  idx_remove = 23
  if (statistical_misfit is True):
    if (gauge_POI is None):
    # The first data points (23) need to be removed because they are observation and not PTF data
    
      ngauges = ds.dimensions['grid_npoints'].size - idx_remove
      gauge_POI = list(range(ngauges))
    
      tmp2_wave_data = np.transpose(tmp_wave_data[:, idx_remove:ngauges+idx_remove])
    
      # The whole syntethic data has many masked entries that need to be removed
      wave_data, maskedarray_indices = remove_masked(tmp2_wave_data)
    else: 
      ngauges = len(gauge_POI)
      wave_data, maskedarray_indices = remove_masked(np.transpose(tmp_wave_data[:, gauge_POI]))
  else:
    wave_data = np.transpose(tmp_wave_data[:,0:idx_remove])
    if (gauge_POI is None):
      gauge_POI = list(range(idx_remove))
    maskedarray_indices = []

  # Calculate new length for gauges
  ngauges = np.shape(wave_data)[0]
  return ngauges, gauge_POI, time_data, wave_data, maskedarray_indices



def get_scenario_waves(N, data_path, data_folders, indices):
  """
  Function to read the scenario data and output them in an array
  """
  # Check index size (mainly for synthetic gauge data)
  ncfile = os.path.join(data_path, data_folders[0], 'out_ts.nc')
  ds = Dataset(ncfile, 'r', format='NETCDF4')
  max_index = ds.dimensions['grid_npoints'].size
  if (indices[-1] > max_index): 
    indices = indices[0:max_index]
  
  wave_data = []
  scenario_min_height = []
  scenario_max_height = []
  for scenario in range(N):
    if (scenario % 250 == 0):
      print(f"Fetching scenario {scenario} out of {N} scenarios.")
    ncfile = os.path.join(data_path, data_folders[scenario], 'out_ts.nc')
    ds = Dataset(ncfile, 'r', format='NETCDF4')
    scenario_wave_amplitude = ds.variables['eta'][:]
    scenario_min_height.append(ds.variables['min_height'][indices])
    scenario_max_height.append(ds.variables['max_height'][indices])
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
  nan_indices = []
  for gauge in range(ngauges):
    minimum_waveheight[gauge] = np.min(gauge_data[gauge])
    maximum_waveheight[gauge] = np.max(gauge_data[gauge])
    if (np.isnan(maximum_waveheight[gauge]) or np.isnan(minimum_waveheight[gauge])):
      print("Gauge entry is a masked array. Setting min/max waveheight to zero.")
      minimum_waveheight[gauge] = 0.
      maximum_waveheight[gauge] = 0.
  return minimum_waveheight, maximum_waveheight, nan_indices
  
  
  
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
# Load gauge data
print("Reading gauge data ")
start = time.time()

if (synthetic_gauge):
  # Synthetic gauge data
  ngauges, gauge_POI, time_list, data_list, maskedarray_indices = get_gauge_synthetic(gauge_data_path, gauge_POI, statistical_misfit)
else:
  # Real gauge data
  # Determine number of gauges from POI list
  if (gauge_POI is None):
    raise ValueError("Real gauge data chosen. In this case, gauge POIs have to be specified!")
  else:
    ngauges = len(gauge_POI)

  # Get gauge ASCII file names
  gauge_list = []
  gauge_list += [each for each in sorted(os.listdir(gauge_data_path))]

  # Read time and waveheight data
  time_list = []
  data_list = []
  gaugedir = os.path.join(gauge_data_path)
  for gauge in range(ngauges):
    time_list.append(np.loadtxt(os.path.join(gaugedir, gauge_list[gauge]))[:,0])
    data_list.append(np.loadtxt(os.path.join(gaugedir, gauge_list[gauge]))[:,1])
stop = time.time()    
print(f"Reading the data took {stop - start} s.\n")

#------------------------------------------------------------------------------------------------------------------------------------------------
"""
Fetch PTF scenario data
"""

# Get data folder names
data_folders = []
data_folders += [each for each in sorted(os.listdir(scenario_data_path))]
N = len(data_folders)

# Read and save the scenario results at the gauge POIs and time data
print("Reading scenario data.")
start = time.time()
scenario_results, scenario_time, scenario_min_height, scenario_max_height = get_scenario_waves(N, scenario_data_path, data_folders, gauge_POI)
stop = time.time()    
print(f"Reading the data took {stop - start} s.\n")

# Get time range and set up time arrays for the gauge data
time_range = [np.min(scenario_time), np.max(scenario_time)]
N_scenario_time = len(scenario_time)

# Remove necessary indices if syntethic data is used
if (synthetic_gauge is True):
  # Determine length of masked array indices
  masked_length = len(maskedarray_indices)
  
  # Initialize new arrays
  tmp_scenario_results = []
  tmp_scenario_min_height = []
  tmp_scenario_max_height = []
  # Delete masked array indices from PTF data
  for scenario in range(N):
    tmp_scenario_min_height.append(np.delete(scenario_min_height[scenario], maskedarray_indices).tolist())
    tmp_scenario_max_height.append(np.delete(scenario_max_height[scenario], maskedarray_indices).tolist())
    for time_idx in range(N_scenario_time):
      tmp_scenario_results.append(np.delete(scenario_results[scenario][time_idx], maskedarray_indices).tolist())
  
  # Overwrite old data
  scenario_results = np.reshape(tmp_scenario_results, (N, N_scenario_time, ngauges))
  scenario_min_height = tmp_scenario_min_height
  scenario_max_height = tmp_scenario_max_height
  
#------------------------------------------------------------------------------------------------------------------------------------------------
"""
Pick arrival times (they have to be used in the time range picking)
"""

# Overwrite gauge times and data so that they match the scenario time range
if (not synthetic_gauge):
  for gauge in range(ngauges):
    data_list[gauge] = data_list[gauge][(time_list[gauge] >= time_range[0]) &
                                                            (time_list[gauge] <= time_range[1])]
    time_list[gauge] = time_list[gauge][(time_list[gauge] >= time_range[0]) &
                                                            (time_list[gauge] <= time_range[1])]
# Or vice-versa if synthetic data is used
else:
  N_new_time = len(scenario_time[scenario_time <= np.max(time_list)])
  new_scenario_results = np.zeros((N, N_new_time, ngauges))
  for scenario in range(N):
    for gauge in range(ngauges):
      new_scenario_results[scenario][:, gauge] = scenario_results[scenario][0:N_new_time, gauge]
      scenario_time = scenario_time[scenario_time <= np.max(time_list)]      
  scenario_results = new_scenario_results
  gauge_Ntime = len(time_list)
  time_list = np.ones((ngauges, gauge_Ntime)) * time_list
  
# Get PTF min/max waveheight (and indices with NaNs for the synthetic data)
min_waveheight, max_waveheight, nan_indices = get_PTF_waveheights(ngauges, data_list)

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
                                                      max(np.abs(scenario_min_height[scenario][gauge]),
                                                              np.abs(scenario_max_height[scenario][gauge])))
    arrtime_trigger = scenario_absmaxwaveheight * arrtime_percentage
    #arrtime_trigger = scenario_absmaxwaveheight_calc * arrtime_percentage
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
PTF_maxindex = np.zeros((N, ngauges), dtype=int)
for gauge in range(ngauges):
  max_index_length = len(gauge_cut_data[gauge])
  gauge_index = gauge_indices[gauge]
  # determine index increment for gauge (is either 1 or 2; gauge_data always has increment 1, so the gauge index will be stored with an increment of 1)
  index_dt = gauge_index[1]-gauge_index[0]  
  len_gauge = len(gauge_index)    
  for scenario in range(N):
    # get length of cut scenario results and check if length has to be changed due to index increment
    len_scenario_old = len(scenario_results_cut[gauge + (ngauges-1)*scenario])
    new_scenario_index = np.arange(0, len_scenario_old, index_dt, dtype=int)
    len_scenario = len(new_scenario_index)
    
    # calculate highest possible index = lowest length of timeseries
    min_len = min(len_scenario, max_index_length, len_gauge)
    PTF_maxindex[scenario, gauge] = min_len
    
    
#------------------------------------------------------------------------------------------------------------------------------------------------
"""
Gauge setup
"""

# Remove index 0 and 21 from misfit if syntethic data is used (because for both, the closest station is over 1° away). Index 3 and 22 are removed because the data is bad.
if (synthetic_gauge):
  if (not statistical_misfit):
    min_waveheight[[0,3,4,21,22]] = 0.
    max_waveheight[[0,3,4,21,22]] = 0.

# Set up GAUGES type 
# dummy coordinates (because they are not needed for this type of measurement)
coords_gauges = np.ones(ngauges)
weigths_gauges = np.ones(ngauges) / ngauges

range_gauges = [np.min(min_waveheight), np.max(max_waveheight)]
GAUGES = TypeMeasurement(ngauges, coords_gauges, weigths_gauges, range_gauges, measuretypeweight = 0.5)
GAUGES.scale_measured_data()


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

#------------------------------------------------------------------------------------------------------------------------------------------------
# Remove index 0 and 21 from misfit if syntethic data is used (because for both, the closest station is over 1° away). Index 3 and 22 are removed because the data is bad.

if (synthetic_gauge):
  if (not statistical_misfit):
    weigths_gauges[[0,3,4,21,22]] = 0.
    GAUGES.renormalize_stationweights()

    weights_arrtime[[0,3,4,21,22]] = 0.
    ARRTIME.renormalize_stationweights()
    print('Removing bad gauge data.')
#------------------------------------------------------------------------------------------------------------------------------------------------
"""
Plotting arrival times (and gauge data)
"""

if (plot_arrivaltimes):
  for gauge in range(ngauges):
    plt.figure()
    plt.plot(gauge_times[gauge], gauge_data[gauge])
    plt.axvline(x = arrival_times[gauge], color = 'b')
    plt.title(f"Gauge data and arrival time (blue line) for {gauge_list[gauge]}")
    plt.xlabel('Time after earthquake = model time [s]')
    plt.ylabel('Wave height [m]')
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
    current_scenario_data = scenario_results_cut[gauge + (ngauges-1)*scenario][0:max_index]
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
Print misfit
"""

print("The misfits are:\n")
print('Arrival time:\n', ARRTIME_misfit)
np.savetxt('arrtime_misfit.txt', ARRTIME_misfit)
print('Gauges:\n', GAUGES_misfit)
np.savetxt('gauge_misfit.txt', GAUGES_misfit)
print('Total:\n', MISFIT)
np.savetxt('total_misfit.txt', MISFIT)

#------------------------------------------------------------------------------------------------------------------------------------------------
"""
Plot misfit
"""

fig = plt.figure()
fig.subplots_adjust(top=0.88)
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

