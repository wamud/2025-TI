import time

start = time.time()
debut = time.process_time()

import subprocess
from concurrent.futures import ThreadPoolExecutor
import sys
import datetime
import numpy as np

def run_executable(distance,rotation, phystheta, physphi, p_error, num_sweeps, num_reps,reuse_aux):
    executable = f"./d{distance}module_{rotation}" if reuse_aux == False else f"./d{distance}module_{rotation}_reuse"
    try:
        cmd = [executable, str(distance), str(phystheta), str(physphi), str(p_error), str(num_sweeps), str(num_reps)]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True)
        
        # print any errors that arise from executable:
        if result.stderr:
            print("Error:", result.stderr)
        return result.stdout

    except subprocess.CalledProcessError as e:
        print("CalledProcessError:", e.stderr)
        print("CalledProcessError code:", e.returncode)
        return e.stderr

def run_parallel_executions(distance,rotation, physthetacoef, physphicoef, p_error, num_sweeps, num_reps,num_workers,reuse_aux):
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(run_executable, distance,rotation, physthetacoef, physphicoef, p_error, num_sweeps, num_reps,reuse_aux) for _ in range(num_workers)]
        results = [future.result() for future in futures]
    return results

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
            detection_events = binary_add_modulo_2(sweep1,sweep2)
            if '1' in detection_events: error_detected = True
            rep_detection_events.append(detection_events)
        all_detection_events.append(rep_detection_events)
        at_least_one_detection_event.append(error_detected)
    return all_detection_events, at_least_one_detection_event

def measure_ZL(distance,data_qubit_measurements):
    ZL_measurements = []
    for i in range(len(data_qubit_measurements)):
        # each list entry is one rep of rabi experiment
        thisrep = data_qubit_measurements[i]
        firstrow = thisrep[0:distance]
        binarysum = (int(firstrow[0]) + int(firstrow[1])) % 2
        if binarysum == 0:
            ZLeigenstate = '|0>_L'
        else:
            ZLeigenstate = '|1>_L'
        ZL_measurements.append(ZLeigenstate)
    return ZL_measurements


# Run the Rabi Experiment:

# Parameters

reuse_aux = True
rotation = 'unrot'

distance = 3

physphicoef = 0
physphi = physphicoef*np.pi

num_sweeps = distance
num_reps = 10  # number of repetitions

num_workers = 1  # I thought putting multiple workers here would run in parallel but it's doing it in serial, someting to do with global python lock?

for p_error in [0]:
    for physthetacoef in [-1.3]:
        phystheta = physthetacoef*np.pi 

        print(f"Running Rabi experiment with p_error={p_error}, d{distance}{rotation}, θ={physthetacoef}π, φ={physphicoef}π, sweeps_per_shot={num_sweeps}, num_reps={num_reps}")
        


        # To actually run it:
        results = run_parallel_executions(distance,rotation, physthetacoef, physphicoef, p_error, num_sweeps, num_reps, num_workers,reuse_aux)

        



        stabiliser_measurements = []
        data_qubit_measurements = []

        for output in results: # parse each worker's output
            stab_measurements, dataq_measurements = parse_output(output, num_sweeps)
            stabiliser_measurements.extend(stab_measurements)
            data_qubit_measurements.extend(dataq_measurements)

        # Stabiliser measurements is now a list of lists, with each list containing the stabiliser measurements (with num_sweeps sweeps) from a rep of the rabi experiment. 
        # The corresponding entry in data_qubit_measurements contains the measurement results of the data qubits at the end of each rep.

        detections, error_detected = find_detection_events(stabiliser_measurements)

        ZL_measurements = measure_ZL(distance, data_qubit_measurements)

        # Get current time to the nearest millisecond for file_name
        current_time = datetime.datetime.now().strftime("%S%f")[:-3]

        if reuse_aux == True:
            filename = f"results_of_run_rabi/processed/d{distance}/p={p_error:.4f}_d{distance}_reuse_{rotation}_theta{physthetacoef:.3f}π_phi{physphicoef:.3f}π_sweeps{num_sweeps}_id{current_time}.dat"
        else:
            filename = f"results_of_run_rabi/processed/d{distance}/p={p_error:.4f}_d{distance}_{rotation}_theta{physthetacoef:.3f}π_phi{physphicoef:.3f}π_sweeps{num_sweeps}_reps{num_reps}_id{current_time}.dat"

        ## Going to print out my results to a file in the following format:
        #  n (num_sweeps) is num of stab. measurement rounds per rep.

        #                                                      line(s):
        # Stabiliser measurements:                              0
        # (n lines o/msrmnts - is traj. if no errors)  1 to n   
        # Detection events:                                     n + 1
        # (n-1) lines of detection events              2n 
        # At least one detection / no detections                2n + 1
        # |1>_L / |0>_L                                          2n + 2
        # (empty line)                                          2n + 3

        # => There is 2n + 4 lines per rep of a rabi experiment, first line is 0-th line, last is (2n + 3)th line

        with open(filename, "w") as f:
            for i in range(len(error_detected)):
                print("Stabiliser measurements:",file=f)
                for bla in stabiliser_measurements[i]: print(bla,file=f)
                print("Detection events:",file=f)
                for bla in detections[i]: print(bla,file=f)
                if error_detected[i]:
                    print("At least one detection",file=f)
                else:
                    print("No detections",file=f)
                print(ZL_measurements[i],file=f)
                print("",file=f)

        print(f" - Results printed to {filename}")



cputime = (time.process_time() - debut)
executiontime = (time.time() - start)
print("")
print("")
print("CPU time = ",str(cputime),"s")
print("Real time = ",str(executiontime),"s")