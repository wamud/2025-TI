// This file runs many shots (repetitions) of a TI Rabi experiment, used in python wrapper RUN_RABI.py. To 'transversally inject' a random but heralded state it rotates all the data qubits then measures all the stabilisers. It repeats the stab. measurements numsweeps times then measures all the data qubits.
// The processing finding detection events and calculating whether it ended up in a ±1 eigenstate of ZL is done in python, so that this simulator can be replaced like a module.

// It has manual idling errors rather than relying on time step count. Also manually doing gate errors.

// Output is each round of stabiliser measurement results then all the data qubits measurement results


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

void perform_stabiliser_measurement(ds_Register reg, int auxindex, int stab_type, bool idling_errors, int pre_aux_idling, int stab_weight, int *datas, int *preceding_data_idlings, int *aux_idlings, int *final_data_idlings) {
    // Time step 1: Perform Hadamard or idling on the auxiliary qubit
    ds_Hadamard_or_idling_on_aux(reg, auxindex, stab_type, idling_errors);

    // Any aux idling errors preceding first CNOT of this stabiliser:
    if (pre_aux_idling != 0) {
        ds_aux_idling(reg, auxindex, pre_aux_idling, idling_errors);
    }

    // Perform the sequence of CNOTs with appropriate idling errors
    for (int cnot_number = 0; cnot_number < stab_weight; cnot_number++) {
        int k = datas[cnot_number]; // get data qubit index
        int num_preceding = preceding_data_idlings[cnot_number]; // number of preceding idling errors on this data qubit
        int aux_idling = aux_idlings[cnot_number]; // number of succeeding auxiliary idling errors
        int num_final_data_idling = final_data_idlings[cnot_number]; // If final CNOT on this data qubit, number of idling errors after

        // Do any preceding data idling errors
        ds_preceding_data_idling(reg, k, num_preceding, idling_errors);
        
        // Perform CNOT
        ds_my_CNOT(reg, stab_type, k, auxindex);
        
        // Do any succeeding auxiliary idling errors
        ds_aux_idling(reg, auxindex, aux_idling, idling_errors);

        // Do any final idling errors on data qubit if last CNOT on this data qubit before hadamard step
        ds_data_idling(reg, k, num_final_data_idling, idling_errors);
    }

    // Time step 6: Perform Hadamard or idling on the auxiliary qubit again
    ds_Hadamard_or_idling_on_aux(reg, auxindex, stab_type, idling_errors);
}



int main(int argc, char *argv[]){  
    double p_error, physphi, phystheta, physthetacoef,physphicoef;
    ds_Register reg; 
    int distance, numsweeps, reps; // reps is number of repetitions/samples to take 
    bool faulty_reset_and_measurement, TI, idling_errors;

    sscanf(argv[1],"%i",&distance); 
    sscanf(argv[2],"%le",&physthetacoef); 
    sscanf(argv[3],"%le",&physphicoef);
    sscanf(argv[4],"%le",&p_error);
    sscanf(argv[5],"%i",&numsweeps);
    sscanf(argv[6],"%i",&reps); 

    phystheta = physthetacoef * M_PI;
    physphi = physphicoef * M_PI;

    TI = true; // are we doing transversal injection? 
    idling_errors = true; 
    faulty_reset_and_measurement = false; // set to true if you want to simulate faulty reset and measurement. Our previous results did not have this set to true. If want to compare to stim files from 'compare the pair' need to change errors to be p/15 for two-qubit gates.


    // printf("theta %f,phi %f , p %f, d %i, sweeps %i,  \n", phystheta,physphi,p_error,distance,numsweeps);

    // Simulation of unrot. d=2: 5 data qubits (numbered 0 to 4) and  1 auxiliary (numbered 5)
    // (Normally 4 auxiliarys but I'm just using one and resetting it)

    int dataqubits, auxqubits, qubits, num_stabs,auxindex,j;


    if (distance != 2) {
        printf("This code is currently only for distance 2. \n");
        return 0;
    }

    dataqubits = pow(distance,2) + pow(distance-1,2); // unrotated
    auxqubits = 1; // re-using one auxiliary (need to change ds_update still)
    qubits = dataqubits + auxqubits;
    num_stabs = dataqubits-1;
    auxindex = dataqubits; // this will be the index for the auxiliary qubit. With 13 data qubits, the auxiliary is the 14th. Or just the 13th (=num_dataqubits) when you count the 0-th one too, as we are doing.

    // I'm imagining qubits labelled from 0 in rows reading left to right, top to bottom. Stabilisers are labelled in the same way: left to right, top to bottom.

    // Set up print to file:

    // Get nanoseconds for multiplying by pid to be a unique identifier in filename
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    unsigned long long nanoseconds = ts.tv_sec * 1000000000LL + ts.tv_nsec;

    char filename[300];
    sprintf(filename, "results_of_run_rabi/unprocessed_data_from_module/d%d/p=%.4f_d%d_reuse_unrot_theta%.3fπ_phi%.3fπ_sweeps%i_reps_%i_id%llu.dat", distance,p_error, distance, physthetacoef, physphicoef, numsweeps, reps, nanoseconds * getpid());
    
    FILE *output = fopen(filename, "a");

    if (output == NULL) {
        perror("Failed to open file");
        return 1; // or handle error appropriately
    }

      // Counts for logical 1 and logical 0: (for testing, this post-processing will be done in python wrapper)
    int count_0 = 0;
    int count_1=0;

    // Loop over 'reps' number of shots of the experiment (prepare transversally injected state, measured num_sweeps rounds of stab measurements, measure in logical basis):
    for(int ii = 0; ii < reps; ii++){  

        fprintf(output,"Rep. %i: \n",ii);
        printf("Rep. %i: \n",ii);

        // Getting 'random' seed for initialising simulator:
        struct timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        unsigned long long nanoseconds = ts.tv_sec * 1000000000LL + ts.tv_nsec;

        ds_initialize_simulator(nanoseconds * getpid());
    
        reg = ds_create_register(qubits, p_error, 0); 
        ds_set_state(reg, 0, 1, 0); 

        // simulate faulty reset:
        if (faulty_reset_and_measurement == true){
            for(j = 0; j < (qubits); j++){
                ds_xerr(reg, j); // probabilistically apply an x error to each qubit
                }
        }

        if (TI == true){
            // Do transversal injection (TI) - i.e. rotate all the dataqubits before the parity checks:
            for(j = 0; j < (dataqubits); j++){
                ds_yrot(reg, j, -phystheta, 0); // '0' in 3rd argument to turn off step count
                ds_lerr(reg, j, 1); // gate error // Gozde had it turned off in her simulations up to now (April 2024)
                ds_zrot(reg, j, -physphi, 0);   //-ve angle as yrot and zrot go opposite direction to usual defn on Bloch sphere
                ds_lerr(reg, j, 1); // gate error // Gozde had it turned off in her simulations up to now (April 2024)
                }
        }



        // PARITY CHECKS (i.e. stabiliser measurements):

        // These involve hadamard on auxiliary then controlled (by auxiliary) stabiliser then hadamard on auxiliary (just turns into CNOTS from data to aux for Z stabs)

        int bit[1];  // intializes an array 'bit' that can hold one integer element. This will be the list of qubits to measure.
        int k; // use to label the data qubit index
        int stab_index,stab_type; //stab_type = 0 for X-type, 1 for Z-type 
        int syndrome[numsweeps][num_stabs];
        int change[numsweeps-1][num_stabs];
        int stab_weight;
        int cnot_number, num_preceding, aux_idling, pre_aux_idling,num_final_data_idling;
        
        
        // Initialise the arrays for different weight (2,3 or 4) stabilisers:
        int datas2[2], preceding_data_idlings2[2], aux_idlings2[2],final_data_idlings2[2];
        int datas3[3], preceding_data_idlings3[3], aux_idlings3[3],final_data_idlings3[3];
        int datas4[4], preceding_data_idlings4[4], aux_idlings4[4],final_data_idlings4[4];

        bit[0] = auxindex;

        // Loop over 'numsweeps' rounds of stabiliser measurements:

        // Stabilisers S0 = XXXII, S1=ZIZZI, S2=IZZIZ, S3=IIXXX but reported S0,S3,S1,S2 (X-types then Z-types)

        for(int sweep = 0; sweep < numsweeps; sweep++) { 
        
            // Before we measure the stabilisers, let's put idling errors on all the data qubits (which occur while the hadamard is occuring on half the auxiliary qubits)
            if (idling_errors) {
                for(int j = 0; j < (dataqubits); j++){
                    ds_lerr(reg, j,1); 
                    }
            }

            // Now stabiliser measurements (parity checks):
            
            // Parity check S0 = X X X I I  so data qubits in order 2,1,0 
                                        
                stab_type = 0; //X-type
                stab_index = 0;
                stab_weight = 3;

                // Aux index idling preceding first sweep (during transversal rotations):
                if (sweep == 0 && idling_errors && TI){ds_idling_during_TI(reg, auxindex);}

                // Any aux idling errors preceding first CNOT of this stabiliser:
                pre_aux_idling = 0;

                // Data qubits to perform CNOTs on
                datas3[0] = 2;
                datas3[1] = 1; 
                datas3[2] = 0;
                
                // How many idling errors immediately precede that data qubit's CNOT
                preceding_data_idlings3[0] = 0;
                preceding_data_idlings3[1] = 1;
                preceding_data_idlings3[2] = 2;
            
                // If last CNOT on this data qubit (if not set all to zero), how many idling errors between it and the hadamard step?
                final_data_idlings3[0] = 0;
                final_data_idlings3[1] = 0;
                final_data_idlings3[2] = 0;

                // Number of aux idling errors immediately following CNOT (max 1 (2) in unrot (rot))
                aux_idlings3[0] = 0;
                aux_idlings3[1] = 0;
                aux_idlings3[2] = 1;


                // Perform stabiliser measurement: (Hadamards if there, CNOTs and appropriate idling errors)
                 perform_stabiliser_measurement(reg, auxindex, stab_type, idling_errors, pre_aux_idling, stab_weight, datas3, preceding_data_idlings3, aux_idlings3, final_data_idlings3);

                // time step 7 (measure & reset in 1 timestep (like stim) but both faulty)
                    if (faulty_reset_and_measurement) {ds_xerr(reg,auxindex);}
                    syndrome[sweep][stab_index] = ds_measure(reg,1,bit);  // measurement. recall: int ds_measure(ds_Register reg, int nq2m, int *lq2m). So 'bit' is list of qubits to measure.

                    // Reset auxiliary to 0 if measured to 1:
                        if (syndrome[sweep][stab_index] == 1){ ds_X(reg, auxindex, 1);}
                    if (faulty_reset_and_measurement == true){ds_xerr(reg,auxindex);}


            // Parity check S1 = Z I Z Z I so data qubits 3,2,0
                
                stab_index = 1;
                stab_type = 1; //Z-type
                stab_weight = 3;
                
                // The setup information below is read off of the stim circuit timeline-svg with CNOT order 2130 and with idling errors inserted:

                // Any aux idling errors between hadamard time step and first CNOT?
                pre_aux_idling = 0;
                
                // Data qubits to perform CNOTs on
                datas3[0] = 3;
                datas3[1] = 2; 
                datas3[2] = 0;

                // How many idling errors immediately precede that data qubit's CNOT
                preceding_data_idlings3[0] = 0;
                preceding_data_idlings3[1] = 0;
                preceding_data_idlings3[2] = 0;

            
                // If last CNOT on this data qubit (if not then set all to zero), how many idling errors between it and the hadamard step?
                final_data_idlings3[0] = 0;
                final_data_idlings3[1] = 0;
                final_data_idlings3[2] = 1;
            
                // Number of aux idling errors following CNOT (max 1 (2) in unrot (rot))
                aux_idlings3[0] = 0;
                aux_idlings3[1] = 1;
                aux_idlings3[2] = 0;

                // Aux idling preceding first sweep (during transversal rotations):
                if (sweep == 0 && idling_errors && TI){ds_idling_during_TI(reg, auxindex);}
                
                // Perform stabiliser measurement: (Hadamards if there, CNOTs and appropriate idling errors)
                    perform_stabiliser_measurement(reg, auxindex, stab_type, idling_errors, pre_aux_idling, stab_weight, datas3, preceding_data_idlings3, aux_idlings3, final_data_idlings3);
                
                // time step 7 (measure & reset in 1 timestep (like stim) but both faulty)
                    if (faulty_reset_and_measurement) {ds_xerr(reg,auxindex);}
                    syndrome[sweep][stab_index] = ds_measure(reg,1,bit);  // measurement. recall: int ds_measure(ds_Register reg, int nq2m, int *lq2m). So 'bit' is list of qubits to measure.

                // Reset auxiliary to 0 if measured to 1:
                        if (syndrome[sweep][stab_index] == 1){ ds_X(reg, auxindex, 1);}
                    if (faulty_reset_and_measurement == true){ds_xerr(reg,auxindex);}


            
            // Parity check S2 = I Z Z I Z so data qubits 4,2,1 
                
                stab_index = 2;
                stab_type = 1; //Z-type
                stab_weight = 3;
            
                // Any aux idling errors between hadamard time step and first CNOT?
                pre_aux_idling = 0;
                
                // Data qubits to perform CNOTs on, any idling errors immediately preceding that CNOT (apart from aux idling during hadamard), any aux idling following it (apart from idling during hadamard):
                datas3[0] = 4;  preceding_data_idlings3[0] = 0; final_data_idlings3[0] = 0; aux_idlings3[0] = 1; 
                datas3[1] = 2;  preceding_data_idlings3[1] = 0; final_data_idlings3[1] = 0; aux_idlings3[1] = 0; 
                datas3[2] = 1;  preceding_data_idlings3[2] = 1; final_data_idlings3[2] = 1; aux_idlings3[2] = 0; 

                // Aux idling preceding first sweep (during transversal rotations):
                if (sweep == 0 && idling_errors && TI){ds_idling_during_TI(reg, auxindex);}

                // Perform stabiliser measurement:
                    perform_stabiliser_measurement(reg, auxindex, stab_type, idling_errors, pre_aux_idling, stab_weight, datas3, preceding_data_idlings3, aux_idlings3, final_data_idlings3);


                // time step 7 (measure & reset in 1 timestep (like stim) but both faulty)
                    if (faulty_reset_and_measurement) {ds_xerr(reg,auxindex);}
                    syndrome[sweep][stab_index] = ds_measure(reg,1,bit);  // measurement. recall: int ds_measure(ds_Register reg, int nq2m, int *lq2m). So 'bit' is list of qubits to measure.

                    // Reset auxiliary to 0 if measured to 1:
                        if (syndrome[sweep][stab_index] == 1){ ds_X(reg, auxindex, 1);}
                    if (faulty_reset_and_measurement == true){ds_xerr(reg,auxindex);}

            // Parity check S3 = I I X X X so data qubits 4,3,2
                
                stab_index = 3;
                stab_type = 0; //0 is X-type, 1 is Z-type
                stab_weight = 3;
            
                // Any aux idling errors between hadamard time step and first CNOT?
                pre_aux_idling = 1;
                
                // Data qubits to perform CNOTs on, any idling errors immediately preceding that CNOT (apart from aux idling during hadamard), any aux idling following it (apart from idling during hadamard):
                datas3[0] = 4;  preceding_data_idlings3[0] = 0;  final_data_idlings3[0] = 2; aux_idlings3[0] = 0;
                datas3[1] = 3;  preceding_data_idlings3[1] = 1;  final_data_idlings3[1] = 1; aux_idlings3[1] = 0;
                datas3[2] = 2;  preceding_data_idlings3[2] = 0;  final_data_idlings3[2] = 0; aux_idlings3[2] = 0;
                
                // Aux idling preceding first sweep (during transversal rotations):
                if (sweep == 0 && idling_errors && TI){ds_idling_during_TI(reg, auxindex);}

                // Perform stabiliser measurement:
                    perform_stabiliser_measurement(reg, auxindex, stab_type, idling_errors, pre_aux_idling, stab_weight, datas3, preceding_data_idlings3, aux_idlings3, final_data_idlings3);

                // time step 7 (measure & reset in 1 timestep (like stim) but both faulty)
                    if (faulty_reset_and_measurement) {ds_xerr(reg,auxindex);}
                    syndrome[sweep][stab_index] = ds_measure(reg,1,bit);  // measurement. recall: int ds_measure(ds_Register reg, int nq2m, int *lq2m). So 'bit' is list of qubits to measure.

                    // Reset auxiliary to 0 if measured to 1:
                        if (syndrome[sweep][stab_index] == 1){ ds_X(reg, auxindex, 1);}
                    if (faulty_reset_and_measurement == true){ds_xerr(reg,auxindex);}

            // STABILISERS MEASURED!

            // Finally, idling on data qubits during time step 7 (while auxiliaries are measured):
            if (idling_errors == true){
                for(int j = 0; j < (dataqubits); j++){
                    ds_lerr(reg, j,1); 
                    }
            }


        } // end of numsweeps rounds of stabiliser measurements loop

        
        
        // All the stabiliser results are now saved in the matrix 'syndrome'. Let's print them out.
            
    // I'm imagining qubits labelled from 0 in rows going left to right, top to bottom. Stabilisers are labelled in the same way: left to right, top to bottom. I'm currently measuring stabilisers in the same way (in the unrotated code they are in nice rows), however we've been reporting trajectories with X stabiliser results (even rows) first, then Z, so will need to reorder to report this.

    // First d-1 stabs are x, next d stabs are Z, next d-1 are X etc. There are d rows of d-1 X stabs and d-1 rows of d Z stabs, implying d(d-1) X stabs and d(d-1) Z stabs:

        fprintf(output,"Stabiliser measurements:\n");
        printf("Stabiliser measurements:\n");

        for (int sweep = 0; sweep < numsweeps; sweep++) {
            
            // Print out X-type stabiliser measurements first (see explanation above)
            for (int x = 0; x < 2 * distance * (distance - 1); x += distance + (distance - 1)) {
                for (int xi = 0; xi < distance - 1; xi++) {
                    int index = x + xi;
                    if (index < num_stabs) { // Ensure within bounds
                        fprintf(output, "%i", syndrome[sweep][index]);
                        printf("%i", syndrome[sweep][index]);
                        
                        // printf("(i=%i)", index); // check it's actually printing correctly
                    }
                }
            }

            // Print out Z-type stabilisers next (see explanation above)
            for (int z = distance - 1; z < 2 * distance * (distance - 1); z += distance + (distance - 1)) {
                for (int zi = 0; zi < distance; zi++) {
                    int index = z + zi;
                    if (index < num_stabs) { // Ensure within bounds
                        fprintf(output, "%i", syndrome[sweep][index]);
                        printf("%i", syndrome[sweep][index]);
                        // printf("(i=%i)", index); // check it's actually printing correctly 
                    }
                }
            }

            // New line at the end of each sweep
            fprintf(output, "\n");
            printf("\n");
        }


            
    // // Previous code which just measured the stabilisers in ascending order:
    //     for(int sweep = 0; sweep < numsweeps; sweep++){
    //         for(int j = 0; j < num_stabs; j++){
    //             fprintf(output,"%i",syndrome[sweep][j]);
    //             printf("%i",syndrome[sweep][j]);
    //         }
    //         fprintf(output,"\n");
    //         printf("\n");
    //     }
        
            
            
        // Now data qubit measurements

        int sum_of_measurement=0; // used logical measurement
        int j;

        printf( "Data qubit measurements:\n");
        fprintf(output,"Data qubit measurements:\n");
        
        // for (j=0; j<2*distance-1; j++){  // when running stab. measurements in parallel, each row would have 2d - 1 qubits (incl. data and auxiliary qubits):

    // Unrotated: when running stab. measurements in serial and re-using the auxiliary, each row only has d data qubits. There are still 2d-1 rows however.
        for (j=0; j<dataqubits; j++){ 
            // // When not re-using auxiliary and just measuring top row, need to skip the auxiliary qubits:
            //   if (j%2==0)  {  //when j is even, measure the qubit 

            // // When re-using auxiliary and measuring all data qubits, don't need to skip any:
            bit[0]=j;
            int measurement = ds_measure(reg,1,bit); // recall: ds_measure(reg,nq2m,*lq2m)
            printf( "%d",measurement);
            fprintf(output, "%d",measurement);
            
            if (j<distance){ // only count the first d qubits for logical Z measurement
                sum_of_measurement = measurement + sum_of_measurement;
            }
            
            
        }
        
        // // See if logical 0 or not (done in python wrapper but can do here for testing):
            
        // if(sum_of_measurement%2 == 1){ // if odd number of ones, then logical 1
        //     count_1=count_1+1;
        //     }
        // else{
        //         count_0=count_0+1;
        //     }
        

        printf("\n\n");
        fprintf(output,"\n\n");

        ds_destroy_register(reg);


        }   // end of for loop going over shots / reps
    
    
    
    // printf("\nFinal count: \n︱0_L〉: %d   ,    ︱1_L〉: %d\n\n",  count_0, count_1); 
    printf("\n pid: %i\n\n",getpid());

    
    return 0;

}

;
/*-----------------------end-------------------------------*/



    


   

