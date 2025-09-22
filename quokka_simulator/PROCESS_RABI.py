# After running RUN_RABI.bash to execute the C code that simulates the Rabi experiments, which prints out check-qubit measurements during the rounds of stabiliser measurements and the data qubit measurements at the end, run this PROCESS_RABI.py to do the classical post-processing of the measurements, namely what the trajectory is, whether it landed in |0>_L or |1>_L etc.

import time
from datetime import datetime
import glob
import os
import numpy as np
import re
from collections import defaultdict
import sys

print(datetime.now().strftime('%H:%M:%S'))
start = time.time()
debut = time.process_time()
def read_output_file(filepath):
    with open(filepath, 'r') as file:
        return file.read()

def parse_output(output, num_sweeps):
    lines = output.splitlines()
    
    stabiliser_measurements = []
    data_qubit_measurements = []
    
    for i, line in enumerate(lines):
        if "Stabiliser measurements" in line:
            stabiliser_measurements.append(lines[i+1:i+1+num_sweeps])
        elif "Data qubit measurements" in line:
            data_qubit_measurements.append(lines[i+1])
    return stabiliser_measurements, data_qubit_measurements

def binary_add_modulo_2(str1, str2):
    if len(str1) != len(str2):
        raise ValueError("Input strings must have the same length")
    
    binary_sum = ''
    for i in range(len(str1)):
        binary_sum += str((int(str1[i]) + int(str2[i])) % 2)
    return binary_sum

def find_detection_events(stabiliser_measurements):
    all_detection_events = []
    at_least_one_detection_event = []
    for rep in stabiliser_measurements:
        rep_detection_events = []
        error_detected = False
        for i in range(len(rep[:-1])):
            sweep1 = rep[i]
            sweep2 = rep[i+1]
            detection_events = binary_add_modulo_2(sweep1, sweep2)
            if '1' in detection_events: error_detected = True
            rep_detection_events.append(detection_events)
        all_detection_events.append(rep_detection_events)
        at_least_one_detection_event.append(error_detected)
    return all_detection_events, at_least_one_detection_event

def measure_ZL(distance, data_qubit_measurements):
    ZL_measurements = []
    for i in range(len(data_qubit_measurements)):
        thisrep = data_qubit_measurements[i]
        firstrow = thisrep[0:distance]
        binarysum = (int(firstrow[0]) + int(firstrow[1])) % 2
        if binarysum == 0:
            ZLeigenstate = '|0>_L'
        else:
            ZLeigenstate = '|1>_L'
        ZL_measurements.append(ZLeigenstate)
    return ZL_measurements

def extract_parameters_from_filename(filename):
    regex_distance = r"d(\d+)_"
    regex_reuse = r"reuse"
    regex_unrotation = r"unrot"
    regex_phystheta = r"theta(-?[\d.]+)π"
    regex_p_error = r"p=([\d.]+)"
    regex_physphi = r"phi(-?[\d.]+)π"
    regex_num_sweeps = r"sweeps(\d+)_"
    regex_unique_id = r"id(\d+)"

    distance = int(re.search(regex_distance, filename).group(1))
    reuse_aux = bool(re.search(regex_reuse, filename))
    rotation = 'unrot' if re.search(regex_unrotation, filename) else 'rot'
    phystheta_match = re.search(regex_phystheta, filename)
    physthetacoef = float(phystheta_match.group(1)) if phystheta_match else None
    phystheta = float(phystheta_match.group(1)) * np.pi if phystheta_match else None
    p_error = float(re.search(regex_p_error, filename).group(1))
    physphi_match = re.search(regex_physphi, filename)
    physphi = float(physphi_match.group(1)) * np.pi if physphi_match else None
    physphicoef = float(physphi_match.group(1)) if physphi_match else None
    num_sweeps = int(re.search(regex_num_sweeps, filename).group(1))

    return distance, reuse_aux, rotation, physthetacoef, phystheta, p_error, physphicoef, physphi, num_sweeps

grouped_files = defaultdict(list)

dist = sys.argv[1]

file_list = glob.glob(f"results_of_run_rabi/unprocessed_data_from_module/d{dist}/*.dat")

print(f"Appending data to files in results_of_run_rabi/processed/d{dist}")

for filepath in file_list:
    filename = os.path.basename(filepath)
    params = extract_parameters_from_filename(filename)
    grouped_files[params].append(filepath)

for params, filepaths in grouped_files.items():
    distance, reuse_aux, rotation, physthetacoef, phystheta, p_error, physphicoef, physphi, num_sweeps = params

    # Construct the processed filename following the original pattern but without the unique ID
    processed_filename = f"results_of_run_rabi/processed/d{distance}/p={p_error:.4f}_d{distance}_{'reuse' if reuse_aux else ''}_{rotation}_theta{physthetacoef}π_phi{physphicoef}π_sweeps{num_sweeps}.dat"

    
    # Ensure the directory exists
    os.makedirs(os.path.dirname(processed_filename), exist_ok=True)
    
    # Initialize the output file
    with open(processed_filename, "w") as f:
        pass

    for filepath in filepaths:
        filename = os.path.basename(filepath)
        # print(f"Reading file {filename}")
        output = read_output_file(filepath)
        
        stabiliser_measurements, data_qubit_measurements = parse_output(output, num_sweeps)

        detection_events_by_sweep, rep_error_flags = find_detection_events(stabiliser_measurements)
        ZL_measurements = measure_ZL(distance, data_qubit_measurements)

        num_digits = len(str(len(rep_error_flags)))

        prefix_pattern = r"results_of_run_rabi/unprocessed/d\d+/"  # Regex pattern
        # Remove the prefix using regex
        filename_stripped = re.sub(prefix_pattern, "", filename, 1)

        print(f"Processing {filename_stripped}")


        with open(processed_filename, "a") as f:
            for i in range(len(rep_error_flags)):
                print(f"Rep {i:0{num_digits}d}", end="\r")
                print("Stabiliser measurements:", file=f)
                for bla in stabiliser_measurements[i]: print(bla, file=f)
                print("Detection events:", file=f)
                for bla in detection_events_by_sweep[i]: print(bla, file=f)
                if rep_error_flags[i]:
                    print("At least one detection", file=f)
                else:
                    print("No detections", file=f)
                print(ZL_measurements[i], file=f)
                print("", file=f)

    print(f" - Results appended to {processed_filename}")

cputime = (time.process_time() - debut)
executiontime = (time.time() - start)
print("")
print("")
print("CPU time = ", str(cputime), "s")
print("Real time = ", str(executiontime), "s")
