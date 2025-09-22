5/9/2025

To visualise potential injected states go to TI/2024/MAIN_RABI_SIMULATORS/Results_of_PLOT_ANALYTICAL_RABI/display_blochs-without both.html

While I worked with Alan on writing a new transversal injection simulator, which is in the git repository Alan-Robertson/transversal_injection (cloned to my OneDrive TI/2024n5/transversal_injection) and managed to do the simulation / calculation (of every trajectory and its probability) of transversal injection to rotated surface code distance 5, the next step of adding errors to the final state (finding the probability mass function of what errors should be applied to the final state to simulate it having passed through a circuit by monte carlo sampling of many runs of the stabiliser extraction circuit in stim to build the PMF) will be done after just having a complete workflow using this quokka simulator.

So this workflow
- simulate transversal injection 
        - rotate all data qubits by same rotation (to non-eigenstate of stabs)
        - apply sydrome extraction circuit -- projects state to an eigenstate of the stabilisers, somewhere on the logical bloch sphere
- build Rabi curve from all the projected states
- do this for varying p (noise) values
- write function that fits a Rabi curve to the data and finds it visibility
- plot visibility versus p value.
- do this for d2 and 3 before mid-October. Also start running 4.


See README.txt within quokka_simulator for how to run it using that.