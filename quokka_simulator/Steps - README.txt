Background (skip if just want steps):

If all of the data qubits in the surface code are rotated to a particular θ and φ before the first round of stabiliser measurements, the logical state after the first round will be a random but heralded logical state somewhere on the logical Bloch sphere. In the absence of errors, the stabiliser measurement results (in this case we call this the "trajectory") can be used to determine exactly what the logical state is. Performing this multiple times enables enough logical states with arbitrary values of θ_L to be arrived at to form a logical Rabi curve, the visibility of which can be used to diagnose underlying error rates of a system and indicate the fidelity with which these logical states are arrived at or 'injected' into the surface code.



Steps:

0. Navigate into quokka_simulator directory

1. Run the ./RUN_RABI.bash after inserting in the bash file the particular distance module, parameters and number of parallel instances you want:
    e.g. ./rabi_simulator_modules/d3rabi 3 0.25 0 0 2 100 > /dev/null &
    would run the distance 3 module, physthetat = 0.25pi, physphi = 0, p = 0, sweeps = 2, reps = 100.

    The only thing the modules print out is qubit measurements - ancillas during the rounds and the top row of data qubits at the end. Post-processing is then done in the python script next.

    (Note I'm using the bash file to run parallel instances because while the file RUN_RABI.py was written to be a nice wrapper that any simulator could be plugged into, it can't run parallel instances due to global python lock)

2. Run python3 PROCESS_RABI.py <distance>  
    
    (this was written with the idea that the only output that will be coming from a quantum computer will be the measurement results, processing of detection events and ZL measurements needs to be done classicaly.)

2. Run python3 PLOT_RABI.py <distance>
    This divides into theta bins and has a cutoff for minimum number of data points to be plotted



Contents:

'Simulator' -- contains the quokka state-simulator files

Results_of_PLOT_ANALYTICAL_RABI has all the analytically calculated Bloch spheres and Rabi curves etc.

results_of_plot_rabi is plots of the results of the simulations