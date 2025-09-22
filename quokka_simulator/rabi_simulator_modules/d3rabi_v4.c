//v4: made changes to work with the whole RUN_RABI workflow, sending output to results_of_run_rabi > unprocessed_data_from_module etc. like my d2_module_unrotated_reuse.c file.
// (Basically I was trying to write code which reused a single auxiliary but applied idling errors as if it was separate auxiliaries. The results that came out of that looked jumbled, and should now be superceded by my C simulator I wrote with Alan once I get the error PMF's. So now I shall just get a fast and loose complete workflow going, which performs all the stabiliser measurements in serial and uses a single auxiliary qubit (note this assumes long-range connections and creates more idling errors.))

// Also changed Z checks to just use CNOTs going from data qubit to auxiliary in zero state.



//include a bunch of stuff from the simulator (TheQ / the quokka) folder:

#include <stdlib.h> 
#include <time.h>
#include "stdio.h"
#include <math.h>
#include "Simulator/sim.h"
#include "Simulator/norm.h"
#include <string.h>
#include <unistd.h>
#include<stdbool.h>
//#include <stdint.h>

// This file runs a transversal injection Rabi experiment for 'reps' number of samples.
// It rotates n data qubits, measures all the stabilisers, measures ZL and reports if the state was |0>_L or |1>_L
// It prints each round of stabiliser results, then the detection events detected therein, then the final measurement of ZL to a file


int main(int argc, char *argv[]){  
    double p_error, physphi, phystheta, physthetacoef,physphicoef;
    ds_Register reg; 
    int distance, numsweeps, reps; // reps is number of repetitions/samples to take 

    sscanf(argv[1],"%i",&distance); 
    sscanf(argv[2],"%le",&physthetacoef); 
    sscanf(argv[3],"%le",&physphicoef);
    sscanf(argv[4],"%le",&p_error);
    sscanf(argv[5],"%i",&numsweeps);
    sscanf(argv[6],"%i",&reps); 

    if (distance != 3) {
        //// printf("This code is currently only for distance 3 \n");
        return 1;
    }

    phystheta = physthetacoef * M_PI;
    physphi = physphicoef * M_PI;

    bool TI = true; // are we doing transversal injection? 
    bool faulty_reset_and_measurement = true; // set to true if you want to simulate faulty reset and measurement. 

    //printf("phi %f , theta %f , p %f, d %i, sweeps %i,  \n", physphi,phystheta,p_error,distance,numsweeps);

// Simulation of unrot. d=3: auxindex data qubits (numbered 0 to 12) and  1 auxiliary (numbered auxindex)
// (Normally 12 auxiliarys but I'm just using one and resetting it)

    int dataqubits, auxqubits, auxindex, qubits, num_stabs;

    dataqubits = pow(distance,2) + pow(distance - 1, 2); // unrotated
    auxqubits = 1; 
    qubits = dataqubits + auxqubits;
    num_stabs = dataqubits - 1;
    auxindex = dataqubits;


    // Set up print to file:
    
    // Get nanoseconds for multiplying by pid to be a unique identifier in filename
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    unsigned long long nanoseconds = ts.tv_sec * 1000000000LL + ts.tv_nsec;

    char filename[300];
    sprintf(filename, "results_of_run_rabi/unprocessed_data_from_module/d%d/p=%.4f_d%d_reuse_unrot_theta%.3fπ_phi%.3fπ_sweeps%i_reps_%i_id%llu.dat", distance,p_error, distance, physthetacoef, physphicoef, numsweeps, reps, nanoseconds * getpid());
    
    FILE *output = fopen(filename, "a");

    if (output == NULL) {
        perror("Hi! Failed to open file -- make sure the required distance directory exists and that you're running this from RUN_RABI.bash");
        return 1; // or handle error appropriately
    }

    // Loop over 'reps' number of shots of the experiment (prepare transversally injected state, measured num_sweeps rounds of stab measurements, measure in logical basis):
    for(int ii = 0; ii < reps; ii++){


        fprintf(output,"Rep. %i: \n",ii);
        //// printf("Rep. %i: \n",ii);

        // Getting 'random' seed for initialising simulator:
        struct timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        unsigned long long nanoseconds = ts.tv_sec * 1000000000LL + ts.tv_nsec;

        ds_initialize_simulator(nanoseconds * getpid());
        
        reg = ds_create_register(qubits, p_error, 0); 
        ds_set_state(reg, 0, 1, 0); // set to |0000...000>


        // simulate faulty reset:
        if (p_error != 0){
            if (faulty_reset_and_measurement == true){
                for(int j = 0; j < (qubits); j++){
                    ds_xerr(reg, j); // probabilistically apply an x error to each qubit
                    }
            }
        }

        // Do transversal injection - i.e. rotate all the dataqubits

        if (TI == true){
            for(int j = 0; j < (dataqubits); j++){
                ds_yrot(reg, j, -phystheta, 1); // 3rd argument adds a time step if 1 and models errors
                
                // ds_lerr(reg, j, 1); // gate error // Gozde had it turned off in her simulations up to now (April 2024)
                // Though this adds a time step, so is doing all the gates in serial in this for loop anyway -- might as well just put a 1 in ds_yrot then -- same effect. Todo: make parallel rotations.
                
            }
            
            if (physphicoef != 0){
                for(int j = 0; j < (dataqubits); j++){
                    ds_zrot(reg, j, -physphi, 1);   //-ve angle as yrot and zrot go opposite direction to usual defn on Bloch sphere
                }
            }
        }

            // ds_print(reg);
            //printf("Above shows: \n state in binary : its coefficient : the number zero \n");




        // PARITY CHECKS (i.e. stabiliser measurements):

        //Parity checks. I.e. hadamard on auxiliary then controlled (by auxiliary) stabiliser then hadamard on auxiliary (just turns into CNOTS for Z stabs)

        int X1, X2, X3, X4, X5, X6, Z1, Z2, Z3, Z4, Z5, Z6; //will store trajectory 
        int bit[1], syndrome[numsweeps][num_stabs], change[numsweeps-1][num_stabs];

        bit[0] = auxindex; //The 13th is the auxiliary qubit which will be measured

        // Loop over 'numsweeps' rounds of stabilizer measurements:
        for(int sweep = 0; sweep < numsweeps; sweep++) { 
            
            // Parity check K1 =  X X I X I I I I I I I I I  

                ds_Hadamard(reg,auxindex,1);

                //Apply controlled XXIXIIIIIIIII

                    static const int indicesK1[] = { 0,1,3 };

                    for (int i = 0; i != sizeof(indicesK1) / sizeof(indicesK1[0]); ++i) {
                        const int k = indicesK1[i];
                        ds_cnot(reg,auxindex,k,1); //ds_cnot(ds_Register reg, int qcont, int qtarg, int time);
                    }

                ds_Hadamard(reg,auxindex,1);
                
                syndrome[sweep][0] = ds_measure(reg,1,bit);  //probabilistic measurement (rather than forced w ds_set_measure)

                //Reset auxiliary if syndrome[sweep][0] stab was measured to one:
                if (syndrome[sweep][0] == 1){ ds_X(reg, auxindex, 1);}


            // Parity check K2 = I X X I X I I I I I I I I  

                ds_Hadamard(reg,auxindex,1);

                //Apply controlled IXXIXIIIIIIII

                    static const int  indicesK2[] = { 1,2,4 };

                    for (int i = 0; i != sizeof( indicesK2) / sizeof( indicesK2[0]); ++i) {
                        const int k =  indicesK2[i];
                        ds_cnot(reg,auxindex,k,1); //ds_cnot(ds_Register reg, int qcont, int qtarg, int time);
                    }

                ds_Hadamard(reg,auxindex,1);

                syndrome[sweep][1] = ds_measure(reg,1,bit);

                //Reset auxiliary if syndrome[sweep][1] stab was measured to one:
                if (syndrome[sweep][1] == 1){ ds_X(reg, auxindex, 1);}

                
            // Parity check K3 =  I I I X I X X I X I I I I 

                ds_Hadamard(reg,auxindex,1);

                //Apply controlled K3

                    static const int  indicesK3[] = { 3,5,6,8 };

                    for (int i = 0; i != sizeof( indicesK3) / sizeof( indicesK3[0]); ++i) {
                        const int k =  indicesK3[i];
                        ds_cnot(reg,auxindex,k,1); //ds_cnot(ds_Register reg, int qcont, int qtarg, int time);
                    }

                ds_Hadamard(reg,auxindex,1);

                syndrome[sweep][2] = ds_measure(reg,1,bit);

                //Reset auxiliary if syndrome[sweep][2] stab was measured to one:
                if (syndrome[sweep][2] == 1){ ds_X(reg, auxindex, 1);}

            // Parity check K4 =   I I I I X I X X I X I I I 

                ds_Hadamard(reg,auxindex,1);

                //Apply controlled K4

                    static const int  indicesK4[] = { 4,6,7,9 };

                    for (int i = 0; i != sizeof( indicesK4) / sizeof( indicesK4[0]); ++i) {
                        const int k =  indicesK4[i];
                        ds_cnot(reg,auxindex,k,1); //ds_cnot(ds_Register reg, int qcont, int qtarg, int time);
                    }

                ds_Hadamard(reg,auxindex,1);
                
                syndrome[sweep][3] = ds_measure(reg,1,bit);

                //Reset auxiliary if syndrome[sweep][3] stab was measured to one:
                if (syndrome[sweep][3] == 1){ ds_X(reg, auxindex, 1);}


            // Parity check K5 =   IIIIIIIIXIXXI 

                ds_Hadamard(reg,auxindex,1);

                //Apply controlled K5

                    static const int  indicesK5[] = { 8,10,11};

                    for (int i = 0; i != sizeof( indicesK5) / sizeof( indicesK5[0]); ++i) {
                        const int k =  indicesK5[i];
                        ds_cnot(reg,auxindex,k,1); //ds_cnot(ds_Register reg, int qcont, int qtarg, int time);
                    }

                ds_Hadamard(reg,auxindex,1);

                syndrome[sweep][4] = ds_measure(reg,1,bit);

                //Reset auxiliary if syndrome[sweep][4] stab was measured to one:
                if (syndrome[sweep][4] == 1){ ds_X(reg, auxindex, 1);}


            // Parity check K6 =   IIIIIIIIIXIXX 

                ds_Hadamard(reg,auxindex,1);

                //Apply controlled K6

                    static const int  indicesK6[] = { 9, 11,12 };

                    for (int i = 0; i != sizeof( indicesK6) / sizeof( indicesK6[0]); ++i) {
                        const int k =  indicesK6[i];
                        ds_cnot(reg,auxindex,k,1); //ds_cnot(ds_Register reg, int qcont, int qtarg, int time);
                    }

                ds_Hadamard(reg,auxindex,1);
                
                syndrome[sweep][5] = ds_measure(reg,1,bit);

                //Reset auxiliary if syndrome[sweep][5] stab was measured to one:
                if (syndrome[sweep][5] == 1){ ds_X(reg, auxindex, 1);}




            // Z-CHECKS !!!



            // Parity check K7 = ZIIZIZIIIIIII

                //Apply controlled K7 

                    static const int indicesK7[] = { 0,3,5 };

                    for (int i = 0; i != sizeof(indicesK7) / sizeof(indicesK7[0]); ++i) {
                        const int k = indicesK7[i];
                        
                        int control = k;
                        int target = auxindex;

                        ds_cnot(reg, control, target, 1);
                    }

                syndrome[sweep][6] = ds_measure(reg,1,bit);

                // If measurement was 1 reset the auxiliary
                if (syndrome[sweep][6] == 1){ ds_X(reg, auxindex, 1);}
                

            // Parity check K8 = IZIZZIZIIIIII


                //Apply controlled K8

                    static const int indicesK8[] = { 1,3,4,6 };

                    for (int i = 0; i != sizeof(indicesK8) / sizeof(indicesK8[0]); ++i) {
                        const int k = indicesK8[i];
                        
                        int control = k;
                        int target = auxindex;

                        ds_cnot(reg, control, target, 1);
                    }


                syndrome[sweep][7] = ds_measure(reg,1,bit);


                if (syndrome[sweep][7] == 1){ ds_X(reg, auxindex, 1);}


            // Parity check K9 = IIZIZIIZIIIII



                //Apply controlled K9

                    static const int indicesK9[] = { 2,4,7 };

                    for (int i = 0; i != sizeof(indicesK9) / sizeof(indicesK9[0]); ++i) {
                        const int k = indicesK9[i];
                        
                        int control = k;
                        int target = auxindex;

                        ds_cnot(reg, control, target, 1);
                    }


                syndrome[sweep][8] = ds_measure(reg,1,bit);


                if (syndrome[sweep][8] == 1){ ds_X(reg, auxindex, 1);}



            // Parity check K10 = IIIIIZIIZIZII


                //Apply controlled K10

                    static const int indicesK10[] = { 5,8,10 };

                    for (int i = 0; i != sizeof(indicesK10) / sizeof(indicesK10[0]); ++i) {
                        const int k = indicesK10[i];
                        
                        int control = k;
                        int target = auxindex;

                        ds_cnot(reg, control, target, 1);
                    }


                syndrome[sweep][9] = ds_measure(reg,1,bit);


                if (syndrome[sweep][9] == 1){ ds_X(reg, auxindex, 1);}


            // Parity check K11 = IIIIIIZIZZIZI


                //Apply controlled K11

                    static const int indicesK11[] = { 6,8,9,11 };

                    for (int i = 0; i != sizeof(indicesK11) / sizeof(indicesK11[0]); ++i) {
                        const int k = indicesK11[i];
                        
                        int control = k;
                        int target = auxindex;

                        ds_cnot(reg, control, target, 1);
                    }


                syndrome[sweep][10] = ds_measure(reg,1,bit);


                if (syndrome[sweep][10] == 1){ ds_X(reg, auxindex, 1);}


            // Parity check K12 = IIIIIIIZIZIIZ


                //Apply controlled K12

                    static const int indicesK12[] = { 7,9,12 };

                    for (int i = 0; i != sizeof(indicesK12) / sizeof(indicesK12[0]); ++i) {
                        const int k = indicesK12[i];
                        
                        int control = k;
                        int target = auxindex;

                        ds_cnot(reg, control, target, 1);
                    }


                syndrome[sweep][11] = ds_measure(reg,1,bit);


                if (syndrome[sweep][11] == 1){ ds_X(reg, auxindex, 1);}
        }



        // All the stabiliser results are now saved in the matrix 'syndrome'. Let's print them out:

        // Stabiliser ordering: I'm imagining qubits labelled from 0 in rows going left to right, top to bottom. Stabilisers are labelled in reading order too. I'm currently measuring stabilisers in the same way (in the unrotated code they are in nice rows), however we've been reporting trajectories with X stabiliser results (even rows) first, then Z, so will need to reorder to report this.

        // First d-1 stabs are x, next d stabs are Z, next d-1 are X etc. There are d rows of d-1 X stabs and d-1 rows of d Z stabs, implying d(d-1) X stabs and d(d-1) Z stabs:


        fprintf(output,"Stabiliser measurements:\n");

        for (int sweep = 0; sweep < numsweeps; sweep++) {
            
            // Print out X-type stabiliser measurements first (see explanation above)
            for (int x = 0; x < 2 * distance * (distance - 1); x += distance + (distance - 1)) {
                for (int xi = 0; xi < distance - 1; xi++) {
                    int index = x + xi;
                    if (index < num_stabs) { // Ensure within bounds
                        fprintf(output, "%i", syndrome[sweep][index]);
                        //// printf("%i", syndrome[sweep][index]);
                        
                        //printf("(i=%i)", index); // check it's actually printing correctly
                    }
                }
            }

            // Print out Z-type stabilisers next (see explanation above)
            for (int z = distance - 1; z < 2 * distance * (distance - 1); z += distance + (distance - 1)) {
                for (int zi = 0; zi < distance; zi++) {
                    int index = z + zi;
                    if (index < num_stabs) { // Ensure within bounds
                        fprintf(output, "%i", syndrome[sweep][index]);
                        //// printf("%i", syndrome[sweep][index]);
                        //printf("(i=%i)", index); // check it's actually printing correctly 
                    }
                }
            }

            // New line at the end of each sweep
            fprintf(output, "\n");
            //// printf("\n");
        }


        // Now let's measure all the data qubits (which can be used to get our ZL measurement). 
        // No matter the trajectory (i.e. eigenvalues of the stabilisers) we will be saying that 0_L is defined as the +1 eigenstate of Z's along the top row. (Transforming ZL to a different chain requires multiplying by stabilisers with their eigenvalue. If one or more of them has a -1 given the particular trajectory this means you might have a mix of -ve and +ve horizontal chains of ZL, so we just stay consistent and say 0_L is a +1 eigenvector of the top chain.

        int j;

        //// printf( "Data qubit measurements:\n");
        fprintf(output,"Data qubit measurements:\n");
        
        // Let's measure just the top row of data qubits to measure ZL. (Skip auxiliary qubits)
        
        // If using all auxiliaries: 
        // for (j=0; j<2*distance-1; j++){  // when running stab. measurements in parallel, each row of qubits has 2d - 1 qubits (incl. data and auxiliary qubits)

        
        // When running stab. measurements in serial with one auxiliary, each row only has d data qubits (There are still 2d -1 rows however):
        for (j = 0; j < distance; j++){  // just measuring the first d data qubits for ZL. todo: change this to measure all data qubits if using that information, in decoding for example.
            
            // If using all auxiliaries, skip each auxiliary:
            //   if (j%2==0)  {  //when j is even, measure the qubit 

            // // When re-using one auxiliary and measuring all data qubits, don't need to skip any:
            bit[0] = j; // list of qubits to measure
            int measurement = ds_measure(reg,1,bit); // recall: ds_measure(reg,nq2m,*lq2m)
            //// printf( "%d",measurement);
            fprintf(output, "%d",measurement);

        }

        //// printf("\n\n");
        fprintf(output,"\n\n");

        ds_destroy_register(reg);

        }  // end of for loop of all the samples

        
        

    //// printf("\n pid: %i\n\n",getpid());

    fclose(output);
    
    return 0;



}

;
/*-----------------------end-------------------------------*/



    


   

