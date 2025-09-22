import time
from datetime import datetime

print(datetime.now().strftime('%H:%M:%S'))
print("")
start = time.time()
debut = time.process_time()

import glob
import os
import numpy as np
import matplotlib.pyplot as plt
#import qutip
import sys
import math
from colorsys import hls_to_rgb
import re
from collections import defaultdict
from scipy.optimize import curve_fit




dist = sys.argv[1]


# New: try group by theta bins:
try_group_by_bins = False


if try_group_by_bins == True:
    num_bins = 20
    cutoff = 5





def find_coefficients(θ,φ):
    # We will be subbing values of a=cos(θ/2) and b=e^(iφ)*sin(θ/2) (from the initial state of n rotated qubits (a|0> + b|1>)^n)
    # into the analytical equations for aL and bL in filea and fileb. Numpy introduces a lot of numerical errors though.
    # For example, sin(2*np.pi / 2) ≠ 0. Hence need to do a bit of tinkering here to get more exact values for a and b:

    #Here are some exact values for theta:

    if θ == 2*np.pi:
        cosθon2 = -1
        sinθon2 = 0
    elif θ == 0:
        cosθon2 = 1
        sinθon2 = 0
    elif θ == np.pi:
        cosθon2 = 0
        sinθon2 = 1
    elif θ == 0.5*np.pi:
        c = 0.7071067811865475
        cosθon2 = c
        sinθon2 = c

    elif θ == 1.5*np.pi:
        c = 0.7071067811865475
        cosθon2 = -c
        sinθon2 = c
    elif θ == 7*np.pi/8:
        cosθon2 = 0.195090322016128
        sinθon2 = 0.980785280403230

    elif θ == np.pi/8:
        cosθon2 = 0.980785280403230
        sinθon2 = 0.195090322016128
    else:
        cosθon2 = np.cos(θ/2)
        sinθon2 = np.sin(θ/2)

    #And some exact values for phi:

    if φ == np.pi:
        expterm = -1
    elif φ == 0:
        expterm = 1
    elif φ == 2*np.pi:
        expterm = 1
    elif φ == np.pi/2:
        expterm = 1j
    elif φ == 3*np.pi/2:
        expterm = -1j
    else:
        expterm = np.exp(1j*φ)

    a = cosθon2
    b = expterm * sinθon2

    return a,b

#region calcLogicalAlphaBeta
# Takes a particular trajectory (binary string of stabiliser measurement results; see end of this comment for explanation of their ordering),
# the contents of two text files which contain symbolic expressions for αL and βL for that trajectory calculated by mathematica, the initial values
# of a and b and returns the values for αL and βL as well as a boolean 'flag'. If the trajectory is impossible to reach, implying an error in the measurements,
# the flag is True, else False.
# 
# ( Order of trajectory: binary string that is stabiliser measurement results: 
#   We list X-type then Z-type stabilisers which are labelled left-to-right, top-to-bottom. 
#   For the rotated surface code we imagine a Z-type stabilier in the top-left corner of the bulk and a weight-2 X-type stabiliser directly above it.
#   For the unrotated code X-type boundaries are the top and bottom, Z-type the left and right. )
#endregion
def calcLogicalAlphaBeta(traj,contenta,contentb,a,b):

    aLeqn = contenta[traj] # This is the equation for alphaL (in terms of alpha (a) and beta (b) physical corresponding to the trajectory 
    bLeqn = contentb[traj] 

    αLprime = eval(aLeqn) # This subs in a and b physical as calculated by find_coefficients for this theta and phi
    βLprime = eval(bLeqn)

    threshold = 1e-6
    if abs(αLprime) < threshold:
        αLprime = 0

    if abs(βLprime) < threshold:
        βLprime = 0

    norm = np.sqrt(abs(αLprime)**2 + abs(βLprime)**2)

    flag = False
    
    if norm == 0:
        # print(f"      Unreachable trajectory {traj}")
        flag = True
        return 0, 0, flag
    else:
        αL = αLprime/norm
        βL = βLprime/norm
        return αL, βL, flag


#region calcLogicalCoords
# This function uses αL and βL from the logical state αL|0_L〉+ βL|1_L〉to find 
# θL and φL with the state in the form cos(θL/2)|0_L〉+ e^(iφL)sin(θL/2)|1_L〉
# Additionally it finds the spherical coords x,y,z to plot on qutip bloch sphere.
#endregion
def calcLogicalCoords(αL,βL):

    if αL == βL == 0:
        print("Unreachable trajectory")
        return 0,0,0,0

    # Now convert this aL and bL to logical theta, logical phi and x,y,z:
    
    aLr = αL.real
    aLi = αL.imag

    r_alpha = math.sqrt((aLr**2) + (aLi**2))
    if r_alpha >  1: # Occasionally comes out as 1.0000000000000002 so round to 1
        r_alpha = 1
    theta = 2 * np.arccos(r_alpha) #THES ES LOGICAL THETA!!

    r_beta = math.sqrt((βL.real**2) + (βL.imag**2))
    theta2 = 2*np.arcsin(r_beta) # should be the same as theta, but just in case

    if abs(theta - theta2) > 1e-6:
        print("Theta and theta2 are not the same, theta = ",theta/np.pi,"π, theta2 = ",theta2/np.pi,"π")

    #Manually set some values to avoid numerical error.
    #If logical theta is 0, for example, then the logical coordinate is just (x,y,z)=(0,0,1).

    if theta == 0 or theta == 2*np.pi:
        x = 0
        y = 0
        z = 1
        phi = 0

    elif theta == np.pi:
        x = 0
        y = 0
        z = -1
        phi = 0
    
    else:

        #Now to figure out logical phi, via each θ from aL and bL.
        # First convert the imaginary numbers aL and bL to their polar coords
        # aL = (r_a)e^iθ_1 , bL = (r_b)e^iθ_2 
        # Also find phi:
        # φ_L = θ_2 - θ_1

        # We must take into account the quadrants of aL and bL to find the correct polar angle θ, which is measured from the x-axis.
        '''
        S   |   A
        ____|____
            |   
        T   |   C
        
        '''

        if aLr > 0 and aLi >= 0:
            aL_theta = math.atan(abs(aLi)/abs(aLr))
        elif aLr < 0 and aLi >= 0:
            aL_theta = np.pi - math.atan(abs(aLi)/abs(aLr))
        elif aLr < 0 and aLi < 0:
            aL_theta = np.pi + math.atan(abs(aLi)/abs(aLr))
        elif aLr > 0 and aLi < 0:
            aL_theta = 2*np.pi - math.atan(abs(aLi)/abs(aLr))
        elif aLr == 0 and aLi > 0:
            aL_theta = np.pi/2
        elif aLr == 0 and aLi < 0:
            aL_theta = 3*np.pi/2
        elif aLr == 0 and aLi ==0:
            aL_theta = 0
        else:
            print(f"Bug in code: it does not account for this αL = {αL} in finding θ_L so have skipped this trajectory")
            print(f"aLi = {aLi}\naLr = {aLr}")
            aL_theta = 0
            with open(f"aL bugs, φ = {φ/np.pi}π , θ = {θ/np.pi}π.txt", "a") as filea:
                filea.write(f"Trajectory = {traj}\naLi = {aLi}\naLr = {aLr}\n")
            error = True



        bLr = βL.real
        bLi = βL.imag

        if bLr > 0 and bLi >= 0:
            bL_theta = math.atan(abs(bLi)/abs(bLr))
        elif bLr < 0 and bLi >= 0:
            bL_theta = np.pi - math.atan(abs(bLi)/abs(bLr))
        elif bLr < 0 and bLi < 0:
            bL_theta = np.pi + math.atan(abs(bLi)/abs(bLr))
        elif bLr > 0 and bLi < 0:
            bL_theta = 2*np.pi - math.atan(abs(bLi)/abs(bLr))
        elif bLr == 0 and bLi > 0:
            bL_theta = np.pi/2
        elif bLr == 0 and bLi < 0:
            bL_theta = 3*np.pi/2
        elif bLr == 0 and bLi ==0:
            bL_theta = 0
        else:
            print(f"Bug in code: it does not account for this bL = {βL}, so have skipped this trajectory")
            print(f"bLi = {bLi}\nbLr = {bLr}")
            bL_theta = 0
            with open(f"bL bugs, φ = {phi/np.pi}π , θ = {th/np.pi}π.txt", "a") as fileb:
                fileb.write(f"Trajectory = {traj}\nbLi = {bLi}\nbLr = {bLr}\n")
        error = True


    # print("φ = ",φ/np.pi,"π , θ = ",θ/np.pi,"π")

        # r_beta = math.sqrt((bLr**2) + (bLi**2))   #not used

        # use global phase to get a real only value for amplitude associated with state |0> on bloch sphere
        phi = bL_theta - aL_theta #Want 0 ≤ ph < 2π

        #The current range that it's producing is 0 < th <= π and -2π <= phi <= 2π
        #Convert to domains for Rabi experiment:
            # -π/2 ≤ ph < π/2  (visually for great circles near XZ plane this makes more sense, rather than using 0 <= phi < π and having points with minutely negative y's being given a theta value of theta = 2π - theta on the Rabi curve)
            # 0 ≤ th < 2π (Rabi curve is |a_L|^2 vs. theta so want 0 ≤ th < 2π)
        pi = np.pi
        if phi == 2*np.pi or phi == -2*np.pi:
            phi = 0 
        if -2*np.pi < phi < -np.pi/2:
            phi = 2*np.pi + phi 
            #no need to change θ 
            #should now have -π/2 ≤ φ ≤ 2π 
        if 3*np.pi/2 <= phi < 2*np.pi:
            phi = phi - 2*np.pi
            #should now have -π/2 ≤ φ < 3π/2  
        if -np.pi/2 <= phi < np.pi/2:
            pass
        if np.pi/2 <= phi <= np.pi: #2nd quadrant inclusive
            phi = phi - np.pi 
            #this moves it from 2nd quadrant to -π/2 ≤ φ ≤ 0 
            theta = 2*np.pi - theta #reach the same point with an increased theta
        if np.pi < phi < 3*np.pi/2 :
            phi = phi - np.pi
            theta = 2*np.pi - theta
            #should now have -π/2 ≤ φ < π/2  
        if phi < -np.pi/2 or phi >= np.pi/2:
            print("φ out of expected domain, φ = ",phi/np.pi,"π")

        x = math.sin(theta)*math.cos(phi)
        y = math.sin(theta)*math.sin(phi)
        z = math.cos(theta)

    phiL = phi
    thetaL = theta  #(J'aurais dû l'appeler comme ça au début)

    return thetaL, x,y,z


def rainbow_color_stops(n, end):
    return [ hls_to_rgb((end * i/(n-1)), 0.5, 1) for i in range(n) ][::-1]

def group_by_theta(trajectorycounts, logicalzeros, thetaLarray, analyticalprobarray, xLp, yLp, zLp):
    # This function sums the trajectorycounts and logicalzerocount for different trajectories under the same theta (if they have the same theta)
    # and calculates the mean for xLp, yLp, and zLp entries corresponding to each unique thetaL.

    duplicates = defaultdict(list)

    for i, thetaL in enumerate(thetaLarray):
        
        # Equate equal thetaL's:
            # accumulate positions of thetaL to the same number:

        # Group equal thetaL values: 

        duplicates[thetaL].append(i)  # duplicates is a dictionary with keys as thetaL and values as the positions in the array where that thetaL occurs. Where any two thetaL's are exactly equal (up to python's floating-point precision), the positions will be grouped under the same key, so one thetaL could have multiple positions in the thetaLarray which is ordered by trajectory number.

    dict = {key: value for key, value in duplicates.items()}

    newthetas = list(dict.keys()) # newthetas is a list of the keys of the dictionary duplicates, which are the unique thetaL values
    positionslist = list(dict.values())

    newprobs = np.empty(len(newthetas)) 
    newprobs[:] = np.nan

    newanalyticalprobs = np.empty(len(newthetas))
    newanalyticalprobs[:] = np.nan

    newxLp = np.empty(len(newthetas))
    newyLp = np.empty(len(newthetas))
    newzLp = np.empty(len(newthetas))

    # For each theta, group (sum) its trajectory count and logical zero count
    for i in range(len(newthetas)):
        theta_positions = positionslist[i]
        trajectorycountsum = np.sum(trajectorycounts[theta_positions])
        logical0sum = np.sum(logicalzeros[theta_positions])

        # Calculate newprobs
        with np.errstate(divide='ignore', invalid='ignore'): # (ignore divide by 0 warning (traj never occurred))
            newprobs[i] = logical0sum/trajectorycountsum

        # Aggregate analyticalprobarray entries corresponding to each unique theta
        # They should all be the same so let's double check that:
        reference_prob = analyticalprobarray[theta_positions][0]
        
        if not np.isnan(reference_prob): # (will be nan if the trajectory was never landed in, and hence the thetaL never calculated)
            for prob in analyticalprobarray[theta_positions]:
                if not (reference_prob - 0.001 <= prob <= reference_prob + 0.001):
                    print("Not all probs for theta = ", newthetas[i], " are within ±0.001 of each other")
        
        analyticalprobsum = np.mean(analyticalprobarray[theta_positions])
        newanalyticalprobs[i] = analyticalprobsum

        # Calculate mean positions for xLp, yLp, zLp
        newxLp[i] = np.mean(xLp[theta_positions])
        newyLp[i] = np.mean(yLp[theta_positions])
        newzLp[i] = np.mean(zLp[theta_positions])

    newthetas = np.array(newthetas) # convert to numpy array
    
    return newthetas, newprobs, newanalyticalprobs, newxLp, newyLp, newzLp

# New: 

def group_by_theta_bins(trajectorycounts, logicalzeros, thetaLarray, analyticalprobarray, xLp, yLp, zLp, num_bins):
    # Define bins between 0 and 2π
    bins = np.linspace(0, 2*np.pi, num_bins+1)
    bin_indices = np.digitize(thetaLarray*np.pi, bins) - 1  # -1 to make indices 0-based

    # Prepare output arrays
    newthetas = np.zeros(num_bins)
    newprobs = np.full(num_bins, np.nan)
    newanalyticalprobs = np.full(num_bins, np.nan)
    newxLp = np.full(num_bins, np.nan)
    newyLp = np.full(num_bins, np.nan)
    newzLp = np.full(num_bins, np.nan)
    counts = np.zeros(num_bins, dtype=int)

    # Aggregate values for each bin
    for i in range(num_bins):
        positions = np.where(bin_indices == i)[0]
        if len(positions) == 0:
            continue
        counts[i] = len(positions)
        trajectory_sum = np.sum(trajectorycounts[positions])
        logical0_sum = np.sum(logicalzeros[positions])
        with np.errstate(divide='ignore', invalid='ignore'):
            newprobs[i] = logical0_sum / trajectory_sum
        newanalyticalprobs[i] = np.mean(analyticalprobarray[positions])
        newxLp[i] = np.mean(xLp[positions])
        newyLp[i] = np.mean(yLp[positions])
        newzLp[i] = np.mean(zLp[positions])
        newthetas[i] = np.mean(thetaLarray[positions])

    return newthetas, newprobs, newanalyticalprobs, newxLp, newyLp, newzLp, counts



def extract_parameters_from_filename_v2(filename):
    # Define regular expressions to extract parameters
    regex_distance = r"d(\d+)_"
    regex_reuse = r"reuse"  # Just checking if 'reuse' is present
    regex_unrotation = r"unrot"  # Search explicitly for 'unrot'
    regex_phystheta = r"theta(-?[\d.]+)π"
    regex_p_error = r"p=([\d.]+)"
    regex_physphi = r"phi(-?[\d.]+)π"
    regex_num_sweeps = r"sweeps(\d+)"

    # Extract parameters using regular expressions
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

    # Return extracted parameters as a tuple
    return distance, reuse_aux, rotation, physthetacoef, phystheta, p_error, physphicoef, physphi, num_sweeps


# # # PLOTTER CODE:
file_list = glob.glob(f"results_of_run_rabi/processed/d{dist}/*.dat") # will plot graphs for all files in this directory (I have added a way to combine files)

# file_list = glob.glob(f"results_of_run_rabi/processed/d{dist}/p=0.0000_d3_reuse_unrot_theta0.700π_phi0.000π_sweeps3_reps_500000_id550539030118980394.dat")

for i,filename in enumerate(file_list):

    # Remove prefix:
    full_string = filename
    prefix_pattern = r"results_of_run_rabi/processed/d\d+/"  # Regex pattern
    # Remove the prefix using regex
    stripped_string = re.sub(prefix_pattern, "", full_string, 1)  # The '1' limits the replacement to the first occurrence


    print("\n",stripped_string)
    distance, reuse_aux, rotation, physthetacoef, phystheta, p_error, physphicoef, physphi, num_sweeps = extract_parameters_from_filename_v2(stripped_string)

    file = open(filename,"r")   
    lines = file.readlines()

    θ = phystheta
    φ = physphi
    n = num_sweeps

    reuse_str = ", reuse aux" if reuse_aux else ""
    print(f"d{distance} {rotation}{reuse_str}, θ = {physthetacoef}π, φ = {physphicoef}π, p = {p_error}, num_rounds = {num_sweeps}")

    numlines = len(lines)
    if numlines == 0:
        print(f"     Empty file: {filename}")
        continue

    # Structure of data file:
    # (where n (num_sweeps) is num of stab. measurement rounds per rep.):

    #                                                      line(s):
    # Stabiliser measurements:                              0
    # (n lines o/msrmnts - is traj. if no errors)           1 to n
    # Detection events:                                     n + 1
    # (n-1) lines of detection events                       2n 
    # At least one detection / no detections                2n + 1
    # |1>_L / |0>_L                                         2n + 2
    # (empty line)                                          2n + 3

    # => There are 2n + 4 lines per rep of a rabi experiment, first line is 0-th line, last is (2n + 3)th line


    lines_per_block = 2*n + 4  # (lines per rep)
    num_blocks = round(numlines/lines_per_block)
    samples = num_blocks
    n_digits = math.floor(math.log10(samples)+1) # finds number of digits in 'samples', to be used for formatting the self-overwriting print statement


    # Open these files w/ eqns for α_L and β_L created with files in logical_coefs folder. E.g."ant-transversal analytical D3 v5.nb"
    filea = open(f"logical_coefs/aL_D{distance}_{rotation}ated.txt","r") 
    fileb = open(f"logical_coefs/bL_D{distance}_{rotation}ated.txt","r")
    contenta = filea.readlines()
    contentb = fileb.readlines()

    # Replace ^ (put in by Mathematica) with ** (python exponentiation operator)
    # Also remove { and } characters if they're there:
    contenta = [line.replace('^', '**').replace('{', '').replace('}', '') for line in contenta]
    contentb = [line.replace('^', '**').replace('{', '').replace('}', '') for line in contentb]



    a,b = find_coefficients(θ,φ)  # actually this was only making a difference because I wasn't rounding near-zero results, then their explosion when they're normalised to a state would be affected by slight differences in the exact values.


    dist = distance
    nq = dist**2 + (dist-1)**2  # number of data qubits (unrotated surface code)
    ns = nq - 1                 # number of stabilisers
    nt = 2**ns                  # number of trajectories

    colours = rainbow_color_stops(nt, 0.9)


    ## The way I'm defining these arrays is only sustainable up to distance 5, after that it's too much memory and I should instead only save the trajectories that appear rather than leaving empty spaces for all of the possible trajectories. This is what I was previously doing for the x,y,z coordinates (see the commented out 'append' code below) and I should be able to do for the trajectory counts, just need to make sure the actual index of the trajectory isn't used for anything else.
    trajectorycounts = np.zeros(2**ns, dtype = int)  
    logicalzeros = np.zeros(2**ns, dtype = int)
    thetaLarray = np.empty(2**ns)
    thetaLarray[:] = np.nan
    analyticalprobarray = np.empty(2**ns)
    analyticalprobarray[:] = np.nan

    xLp = np.empty(2**ns)
    xLp[:] = np.nan
    yLp = np.empty(2**ns)
    yLp[:] = np.nan
    zLp = np.empty(2**ns)
    zLp[:] = np.nan


    ## POST-SELECTION:
    ## Let's start reading the memory experiment results and post-select for ones which saw no errors:

    n_selected = 0  # count amount that are post-selected

    for k in range (num_blocks):   
        # extract the line which says if there were any detection events:
        any_detections = lines[2*n + 1 + k*lines_per_block].strip("\n") # strip removes the \n symbol from each line

        # Begin post-selection:

        if any_detections != "No detections": # != is does not equal
            # implies file reads "At least one detection event" so don't post-select this run
            continue

        traj_binary = lines[1 + k*lines_per_block].strip("\n") # trajectory in binary string
        traj = int(traj_binary,2)                              # trajectory in integer form
        
        trajectorycounts[traj] += 1 # adds 1 to relevant trajectory count, is keeping count of how many times each trajectory appears.

        # For the trajectory from each sample, I will need its corresponding thetaL to later plot the rabi curve. As my method of finding thetaL requires finding x,y,z in the Bloch sphere anyway, I am also going to find each trajectory's logical bloch sphere coordinates as it arrives. Furthermore I will save its |αL|^2 vs. thetaL to build up the analytical rabi curve. These calculations will only be done once for each trajectory.

        ZL_measurement_result = lines[2*n+2 + k*lines_per_block].strip("\n")

        if trajectorycounts[traj] >= 2: # (trajectory has been landed in before)
            n_selected += 1

            if '0' in ZL_measurement_result:
                logicalzeros[traj] += 1 # adds 1 to count of logicalzeros at index corresponding to traj.

            print("Trajectory = ",str(traj_binary).zfill(4),", logical 0 count = ",str(logicalzeros[traj]).zfill(n_digits-1),"/",str(trajectorycounts[traj]).zfill(n_digits-1),", sample = ", str(k).zfill(n_digits-1),end = '\r')
            continue 

        # Else, it's the first time landed in so need to calculate logical coords (and check that it is possible to land in this trajectory) before adding to trajectory count

        αL, βL, flag = calcLogicalAlphaBeta(traj,contenta,contentb,a,b)
        
        if flag == True:
            # an unreachable trajectory (given initial theta and phi) was reached by the simulation, indicating an error in the very first syndrome measurement

            trajectorycounts[traj] = 0 # reset that trajectory's count

            # alternatively could keep a record of the unreachable trajectories which are often appearing and not re-calculate that they are unreachable every time, however I will add this later only if it will save a substantial amount of comp. time because I then need to account for the fact that these trajectories have non-zero counts (which might be used later for processing and plotting)

            continue # go to next iteration of rabi experiment coz this one's errored
        
        # Is a reachable trajectory so add to count:
        analyticalprobarray[traj] = round(abs(αL)**2,10)
        n_selected += 1

        if '0' in ZL_measurement_result:
            logicalzeros[traj] += 1 # adds 1 to count of logicalzeros at index corresponding to traj.

        print("Trajectory = ",str(traj).zfill(4),", logical 0 count = ",str(logicalzeros[traj]).zfill(n_digits-1),"/",str(trajectorycounts[traj]).zfill(n_digits-1),", sample = ", str(k).zfill(n_digits-1),end = '\r')


        # Now calc. position on logical bloch sphere: 
        thetaL,xL,yL,zL = calcLogicalCoords(αL,βL) # (comment out and calculat below instead just for grouped thetas rather than every theta (where some thetas appear more than once as different trajectories have the same thetaL))


        #Save analytical thetaL and points for trajectory which appeared in simulator. These are ordered by trajectory number, any blanks means the trajectory didn't appear in the simulation.

        thetaLarray[traj] = round(thetaL/np.pi,10)   # (will be plotted in units of π)
        xLp[traj] = xL
        yLp[traj] = yL
        zLp[traj] = zL


        # # save the points (these will be in order they appeared, not in trajectory order)
        # xLp = np.append(xLp,[xL])  # define new array as old array with xL appended. (unlike appending to lists)
        # yLp = np.append(yLp,[yL])
        # zLp = np.append(zLp,[zL])

    if n_selected == 0:
        print("No post-selected runs")
        continue # continues to next file
    

    # Above 'for' loop has created:
    # - trajectorycounts: how many times simulation landed in each trajectory w no detections
    # - logicalzeros: how many times --> ︱0 〉_L for each traj. given no detections 
    # - thetaLarray and analyticalprobarray, ordered by trajectory with analytical values for thetaL (in units of π) and abs(α)^2 for the trajectories which appeared in the simulation
    # - x,y,z points for visualising the analytical trajectories on the Bloch sphere

    # Calculate and save simulated probabilities:
    with np.errstate(divide='ignore', invalid='ignore'):  #(ignore divide by zero warning, which occurs when a trajctory was never landed in - scatter plot wil ignore those points anyway)
        problogical0 = logicalzeros/trajectorycounts
        trajprobs = trajectorycounts/n_selected   # creates list w the probs of landing in each traj., given no detection events


    ## Print out data (note that you can get nan values for the thetaL of a particular trajectory if it was never landed in because it was never calculated)

    # for j in range(len(trajectorycounts)):
    #     print("")
    #     print(f"j = {j}")
    #     print(f"thetaL = {thetaLarray[j]}")
    #     print(f"trajectorycount = {trajectorycounts[j]}")
    #     print(f"analyticalprob = {analyticalprobarray[j]}")
    #     print(f"problogical0 = {problogical0[j]}")
    #     print(f"trajprob = {trajprobs[j]}")


    # Each trajectory has a thetaL value. For trajectories with exactly the same thetaL, I will group their counts. groupedthetas is now an array of unique thetaL values and groupedprobs is the average of the probabilities of the grouped trajectories. 

    # original:
    if try_group_by_bins == False:
        groupedthetas,groupedprobs,groupedanalyticalprobs,groupedxLp,groupedyLp,groupedzLp = group_by_theta(trajectorycounts,logicalzeros,thetaLarray,analyticalprobarray,xLp,yLp,zLp)

    # New: try group by bins:
    else:
        groupedthetas, groupedprobs, groupedanalyticalprobs, groupedxLp, groupedyLp, groupedzLp, counts = group_by_theta_bins(trajectorycounts, logicalzeros, thetaLarray, analyticalprobarray, xLp, yLp, zLp, num_bins)



    ## Save the calculated probabilities, gropued thetas, trajectory counts etc.
    file_string = "results_of_plot_rabi/trajcounts_probs_thetas/d{}/d{},θ={}π, φ={}π, p={} (post-selected {} of {}).npz".format(distance,distance,format(θ/np.pi, '.3f'),format(φ/np.pi, '.3f'),p_error,n_selected,samples)
    np.savez_compressed(file_string, 
    a = trajectorycounts, b = trajprobs, c = problogical0, d = thetaLarray,e=analyticalprobarray,f=groupedthetas,g=groupedanalyticalprobs,h=groupedprobs)

    print("")
    # # Print out trajctory counts:
    # for i,count in enumerate(trajectorycounts):
    #     traj = str(i).zfill(ns) 
    #     print(f"{traj} occurred {count} times")

    # # Post-selecting for trajectories which occurred above a threshold amount of times:
    # threshold = 6000
    # for i,count in enumerate(trajectorycounts):
    #     if count < threshold:
    #         trajectorycounts[i] = 0
    #         logicalzeros[i] = 0


    # PLOTTING:

    titlefontsize = 25
    tickfontsize = 20
    labelfontsize = 25
    dotsize_analytical = 60  # (for scatter plot)
    dotsize_simulation = 60  #(for scatter plot)

    # Plot ANALYTICAL Rabi curve:

    # add cos^2 curve:
    x = np.arange(0,2*np.pi,0.01*np.pi)   # start,stop,step
    y = (np.cos(x/2))**2 

    plt.figure(figsize=(9.6,7.2))
    plt.tight_layout()
    plt.plot(x/np.pi,y, c = 'C0',linewidth = 0.5)    #the cos^2 curve
    plt.xlabel('$θ_L$',fontsize = labelfontsize)
    plt.ylabel('P(0 logical)',fontsize = labelfontsize)
    plt.title('d{}{}, θ={}π, φ={}π, p={})'.format(distance,rotation,format(θ/np.pi, '.3f'),format(φ/np.pi, '.3f'),p_error,n_selected,num_blocks),fontsize = titlefontsize)
    steps = [0,1/4,1/2,3*1/4,1,5*1/4,3*1/2,7*1/4,2*1]
    plt.xticks( steps, ['0', 'π/4', 'π/2', '3π/4', 'π', '5π/4', '3π/2', '7π/4', '2π'], fontsize = tickfontsize)
    plt.yticks(fontsize = tickfontsize)

    # add ANALYTICAL values of trajectories which appeared in simulated experiment (groupedanalyticalprobs is the analytical values)
    plt.scatter(groupedthetas,groupedanalyticalprobs, s = dotsize_analytical, color = colours[:len(groupedthetas)],label = 'analytical')
    plt.legend(fontsize=tickfontsize)
    plt.tight_layout()


    plt.savefig("results_of_plot_rabi/rabi curves/d{}/θ={}π, φ={}π, p={} analytical.png".format(distance,format(round(θ/np.pi,4), '.4f'),format(round(φ/np.pi,4), '.4f'),p_error,n_selected,samples))
    plt.clf()
    plt.close()


    # Bloch sphere analytical graph (visualising the analytical values on bloch sphere)

   # B = qutip.Bloch()
   # B.figsize = [5,6.2]
   # B.make_sphere()
   # B.zlabel = ["θ = {}π  , φ = {}π \n ∣0⟩ᴸ".format(format(round(θ/np.pi,2), '.2f'),format(round(φ/np.pi,2), '.2f')),'  ∣1⟩ᴸ'] 
   # B.zlpos = [1.3,-1.25]
   # pnts = [groupedxLp, groupedyLp, groupedzLp]
   # B.point_size = [dotsize_analytical/1.5]
   # B.point_marker = ["o"]
   # B.point_color = ['#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255)) for r,g,b in colours[:len(groupedthetas)]]
   # B.add_points(pnts,'m')
   # B.render()
   # B.fig.savefig("results_of_plot_rabi/bloch spheres/d{}/θ={}π, φ={}π, p={} (post-selected {} of {}).png".format(distance,format(round(θ/np.pi,4), '.4f'),format(round(φ/np.pi,4), '.4f'),p_error,n_selected,samples))
   # B.fig.savefig("results_of_plot_rabi/bloch sphere.png")

    # Could add functionality to the plots to actually match up colours. Would be possible by adding each point in for loops but probs a faster way too ...


    #plt.clf()    #get rid of the analytical points (saved under curve - analytical)
    #plt.close()


    # Plot SIMULATED Rabi curve (with results from simulation)

    plt.figure(figsize=(9.6,7.2))
    x = np.arange(0,2*np.pi,0.01*np.pi) 
    plt.tight_layout()

    # Plot simulated points:    
    # ('groupedprobs' contains the probabilities calculated from simulation):
    ## original:
    
    if try_group_by_bins == False:
        plt.scatter(groupedthetas,groupedprobs,s = dotsize_simulation,label = 'simulated',color = colours[:len(groupedthetas)])  #automatically doesn't plot the points which were never landed in (will have nan values for theta)


    # New: try group by bins:
    else:
        valid_bins = counts >= cutoff
        plt.scatter(groupedthetas[valid_bins],groupedprobs[valid_bins],s = dotsize_simulation,label = 'simulated',color = colours[:len(groupedthetas[valid_bins])])  #automatically doesn't plot the points which were never landed in (will have nan values for theta)
    

    def cos_curve(x, A):
        return A * np.cos(x) + 0.5  # A Rabi curve, centred at y = 0.5
    
    ## Fit cos curve to data:

    ## Filter out NaN values from groupedthetas and groupedprobs to be able to do fit
    # valid_indices = ~np.isnan(groupedthetas) & ~np.isnan(groupedprobs)
    # groupedthetas_valid = groupedthetas[valid_indices]
    # groupedprobs_valid = groupedprobs[valid_indices]\    
    # print(groupedthetas_valid, groupedprobs_valid)
    # popt, pcov = curve_fit(cos_curve, groupedthetas_valid, groupedprobs_valid,p0 = [1])
    # A = popt[0]  # amplitude
    # visibility = 2*A
    # print("Visibility = ",visibility)  

    ## Above fit is terrible for d2. Will just plot a cos curve for now

    plt.plot(x/np.pi,cos_curve(x,0.5), c = 'C0',linewidth = 0.5,label = 'theoretical curve')

    # Format graph:
    plt.xlabel('$θ_L$',fontsize = labelfontsize)
    plt.ylabel('P(0 logical)',fontsize = labelfontsize)
    steps = [0,1/4,1/2,3*1/4,1,5*1/4,3*1/2,7*1/4,2*1]
    plt.xticks( steps, ['0', 'π/4', 'π/2', '3π/4', 'π', '5π/4', '3π/2', '7π/4', '2π'], fontsize = tickfontsize)
    plt.yticks(fontsize = tickfontsize)
    plt.title('d{}{}, θ={}π, φ={}π, p={} \n (post-selected {}/{})'.format(distance,rotation,format(θ/np.pi, '.3f'),format(φ/np.pi, '.3f'),p_error,n_selected,samples),fontsize = titlefontsize)
    plt.legend(fontsize=tickfontsize)    #gives a legend to whichever ones you've labelled
    plt.tight_layout()
    plt.savefig("results_of_plot_rabi/rabi curves/d{}/θ={}π, φ={}π, p={} (post-selected {} of {}).png".format(distance,format(round(θ/np.pi,4), '.4f'),format(round(φ/np.pi,4), '.4f'),p_error,n_selected,samples))
    plt.savefig("results_of_plot_rabi/rabi curve.png")
    plt.clf()
    plt.close()

    filea.close()
    fileb.close()
    file.close()

print("")
print(datetime.now().strftime('%H:%M:%S'))
cputime = (time.process_time() - debut)
executiontime = (time.time() - start)
print("CPU time = ",str(round(cputime,5)),"s")
print("Real time = ",str(round(executiontime,5)),"s")
