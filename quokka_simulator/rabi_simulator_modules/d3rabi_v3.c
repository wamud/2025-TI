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
    double p_error, physphi, phystheta;
    ds_Register reg; 
    int distance, numsweeps, reps; // reps is number of repetitions/samples to take 

    sscanf(argv[1],"%le",&phystheta); 
    sscanf(argv[2],"%le",&physphi);
    sscanf(argv[3],"%le",&p_error);
    sscanf(argv[4],"%i",&distance); 
    sscanf(argv[5],"%i",&numsweeps);
    sscanf(argv[6],"%i",&reps); 



    printf("phi %f , theta %f , p %f, d %i, sweeps %i,  \n", physphi,phystheta,p_error,distance,numsweeps);

// Simulation of unrot. d=3: 13 data qubits (numbered 0 to 12) and  1 auxiliary (numbered 13)
// (Normally 12 auxiliarys but I'm just using one and resetting it)

    int dataqubits, auxqubits, qubits, num_stabs;

    dataqubits = pow(distance,2) + pow(distance -1,2); // unrotated
    auxqubits = 1; // re-using one auxiliary (need to change ds_update still)
    qubits = dataqubits + auxqubits;
    num_stabs = dataqubits-1;


    // Set up print to file:

    int pid;
    pid = getpid(); // to differentiate the files and combine them later

    char filename[100];

    // LEARN HOW TO USE PI AND ROUNDING FORMAT IN C SO CAN EXPRESS ANGLE VALUES IN FILENAME IN TERMS OF PI ROTATION. E.G. IN PYTHON ITS ("Results from FULLRABI_v4/data/φ={}π, θ={}π ({} samples)".format(format(φ/np.pi, '.3f'),format(θ/np.pi, '.3f'),samples),

    sprintf(filename, "Results/data/d%d_p%.4f_theta%.3f_phi%.3f_sweeps%i_pid%i.dat",distance, p_error, phystheta, physphi, numsweeps,pid);
    FILE *output = fopen(filename, "a");

    for(int ii = 0; ii < reps; ii++){

    printf("%i \n",ii);
    
    ds_initialize_simulator((unsigned) time(NULL)*getpid());  //////////// put a random number generator here instead!
    
    reg = ds_create_register(qubits, p_error, 0); 
    ds_set_state(reg, 0, 1, 0); 

    // Do transversal injection - i.e. rotate all the dataqubits
        for(int j = 0; j < (dataqubits); j++){
            ds_yrot(reg, j, -phystheta, 0);
            ds_zrot(reg, j, -physphi, 0);   //-ve angle as yrot and zrot go opposite direction to how theta and phi specified on Bloch sphere
        }
        // ds_print(reg);
        // printf("Above shows: \n state in binary : its coefficient : the number zero \n");


    //Parity checks. I.e. hadamard on auxiliary then controlled (by auxiliary) stabiliser then hadamard on auxiliary (just turns into CNOTS for Z stabs)

    int X1, X2, X3, X4, X5, X6, Z1, Z2, Z3, Z4, Z5, Z6; //will store trajectory 
    int bit[1], syndrome[numsweeps][num_stabs], change[numsweeps-1][num_stabs];

    bit[0] = 13; //The 13th is the auxiliary qubit which will be measured

    // Loop over 'numsweeps' rounds of stabilizer measurements:
    for(int sweep = 0; sweep < numsweeps; sweep++) { 
        
        // Parity check K1 =  X X I X I I I I I I I I I  

            ds_Hadamard(reg,13,1);

            //Apply controlled XXIXIIIIIIIII

                static const int indicesK1[] = { 0,1,3 };

                for (int i = 0; i != sizeof(indicesK1) / sizeof(indicesK1[0]); ++i) {
                    const int k = indicesK1[i];
                    ds_cnot(reg,13,k,1); //ds_cnot(ds_Register reg, int qcont, int qtarg, int time);
                }

            ds_Hadamard(reg,13,1);
            
            syndrome[sweep][0] = ds_measure(reg,1,bit);  //probabilistic measurement (rather than forced w ds_set_measure)

            //Reset auxiliary if syndrome[sweep][0] stab was measured to one:
            if (syndrome[sweep][0] == 1){ ds_X(reg, 13, 1);}


        // Parity check K2 = I X X I X I I I I I I I I  

            ds_Hadamard(reg,13,1);

            //Apply controlled IXXIXIIIIIIII

                static const int  indicesK2[] = { 1,2,4 };

                for (int i = 0; i != sizeof( indicesK2) / sizeof( indicesK2[0]); ++i) {
                    const int k =  indicesK2[i];
                    ds_cnot(reg,13,k,1); //ds_cnot(ds_Register reg, int qcont, int qtarg, int time);
                }

            ds_Hadamard(reg,13,1);

            syndrome[sweep][1] = ds_measure(reg,1,bit);

            //Reset auxiliary if syndrome[sweep][1] stab was measured to one:
            if (syndrome[sweep][1] == 1){ ds_X(reg, 13, 1);}

            
        // Parity check K3 =  I I I X I X X I X I I I I 

            ds_Hadamard(reg,13,1);

            //Apply controlled K3

                static const int  indicesK3[] = { 3,5,6,8 };

                for (int i = 0; i != sizeof( indicesK3) / sizeof( indicesK3[0]); ++i) {
                    const int k =  indicesK3[i];
                    ds_cnot(reg,13,k,1); //ds_cnot(ds_Register reg, int qcont, int qtarg, int time);
                }

            ds_Hadamard(reg,13,1);

            syndrome[sweep][2] = ds_measure(reg,1,bit);

            //Reset auxiliary if syndrome[sweep][2] stab was measured to one:
            if (syndrome[sweep][2] == 1){ ds_X(reg, 13, 1);}

        // Parity check K4 =   I I I I X I X X I X I I I 

            ds_Hadamard(reg,13,1);

            //Apply controlled K4

                static const int  indicesK4[] = { 4,6,7,9 };

                for (int i = 0; i != sizeof( indicesK4) / sizeof( indicesK4[0]); ++i) {
                    const int k =  indicesK4[i];
                    ds_cnot(reg,13,k,1); //ds_cnot(ds_Register reg, int qcont, int qtarg, int time);
                }

            ds_Hadamard(reg,13,1);
            
            syndrome[sweep][3] = ds_measure(reg,1,bit);

            //Reset auxiliary if syndrome[sweep][3] stab was measured to one:
            if (syndrome[sweep][3] == 1){ ds_X(reg, 13, 1);}


        // Parity check K5 =   IIIIIIIIXIXXI 

            ds_Hadamard(reg,13,1);

            //Apply controlled K5

                static const int  indicesK5[] = { 8,10,11};

                for (int i = 0; i != sizeof( indicesK5) / sizeof( indicesK5[0]); ++i) {
                    const int k =  indicesK5[i];
                    ds_cnot(reg,13,k,1); //ds_cnot(ds_Register reg, int qcont, int qtarg, int time);
                }

            ds_Hadamard(reg,13,1);

            syndrome[sweep][4] = ds_measure(reg,1,bit);

            //Reset auxiliary if syndrome[sweep][4] stab was measured to one:
            if (syndrome[sweep][4] == 1){ ds_X(reg, 13, 1);}


        // Parity check K6 =   IIIIIIIIIXIXX 

            ds_Hadamard(reg,13,1);

            //Apply controlled K6

                static const int  indicesK6[] = { 9, 11,12 };

                for (int i = 0; i != sizeof( indicesK6) / sizeof( indicesK6[0]); ++i) {
                    const int k =  indicesK6[i];
                    ds_cnot(reg,13,k,1); //ds_cnot(ds_Register reg, int qcont, int qtarg, int time);
                }

            ds_Hadamard(reg,13,1);
            
            syndrome[sweep][5] = ds_measure(reg,1,bit);

            //Reset auxiliary if syndrome[sweep][5] stab was measured to one:
            if (syndrome[sweep][5] == 1){ ds_X(reg, 13, 1);}


        // Parity check K7 = ZIIZIZIIIIIII

            ds_Hadamard(reg,13,1);

            //Apply controlled K7

                static const int indicesK7[] = { 0,3,5 };

                for (int i = 0; i != sizeof(indicesK7) / sizeof(indicesK7[0]); ++i) {
                    const int k = indicesK7[i];
                    
                    int control = 13;
                    int target = k;

                    //The CZ:
                    ds_Hadamard(reg, target, 1);
                    ds_cnot(reg, control, target, 1);
                    ds_Hadamard(reg, target, 1);
                }

            ds_Hadamard(reg,13,1);

            // ds_set_measure(reg, 1, bit, Z1);

            syndrome[sweep][6] = ds_measure(reg,1,bit);

            // If measurement was 1 reset the auxiliary
            if (syndrome[sweep][6] == 1){ ds_X(reg, 13, 1);}
            


        
        // Parity check K8 = IZIZZIZIIIIII

            ds_Hadamard(reg,13,1);

            //Apply controlled K8

                static const int indicesK8[] = { 1,3,4,6 };

                for (int i = 0; i != sizeof(indicesK8) / sizeof(indicesK8[0]); ++i) {
                    const int k = indicesK8[i];
                    
                    int control = 13;
                    int target = k;

                    //The CZ:
                    ds_Hadamard(reg, target, 1);
                    ds_cnot(reg, control, target, 1);
                    ds_Hadamard(reg, target, 1);
                }

            ds_Hadamard(reg,13,1);

            syndrome[sweep][7] = ds_measure(reg,1,bit);


            if (syndrome[sweep][7] == 1){ ds_X(reg, 13, 1);}


        // Parity check K9 = IIZIZIIZIIIII

            ds_Hadamard(reg,13,1);


            //Apply controlled K9

                static const int indicesK9[] = { 2,4,7 };

                for (int i = 0; i != sizeof(indicesK9) / sizeof(indicesK9[0]); ++i) {
                    const int k = indicesK9[i];
                    
                    int control = 13;
                    int target = k;

                    //The CZ:
                    ds_Hadamard(reg, target, 1);
                    ds_cnot(reg, control, target, 1);
                    ds_Hadamard(reg, target, 1);
                }

            ds_Hadamard(reg,13,1);

            syndrome[sweep][8] = ds_measure(reg,1,bit);


            if (syndrome[sweep][8] == 1){ ds_X(reg, 13, 1);}



        // Parity check K10 = IIIIIZIIZIZII

            ds_Hadamard(reg,13,1);

            //Apply controlled K10

                static const int indicesK10[] = { 5,8,10 };

                for (int i = 0; i != sizeof(indicesK10) / sizeof(indicesK10[0]); ++i) {
                    const int k = indicesK10[i];
                    
                    int control = 13;
                    int target = k;

                    //The CZ:
                    ds_Hadamard(reg, target, 1);
                    ds_cnot(reg, control, target, 1);
                    ds_Hadamard(reg, target, 1);
                }

            ds_Hadamard(reg,13,1);

            syndrome[sweep][9] = ds_measure(reg,1,bit);


            if (syndrome[sweep][9] == 1){ ds_X(reg, 13, 1);}


        // Parity check K11 = IIIIIIZIZZIZI

            ds_Hadamard(reg,13,1);

            //Apply controlled K11

                static const int indicesK11[] = { 6,8,9,11 };

                for (int i = 0; i != sizeof(indicesK11) / sizeof(indicesK11[0]); ++i) {
                    const int k = indicesK11[i];
                    
                    int control = 13;
                    int target = k;

                    //The CZ:
                    ds_Hadamard(reg, target, 1);
                    ds_cnot(reg, control, target, 1);
                    ds_Hadamard(reg, target, 1);
                }

            ds_Hadamard(reg,13,1);

            syndrome[sweep][10] = ds_measure(reg,1,bit);


            if (syndrome[sweep][10] == 1){ ds_X(reg, 13, 1);}


        // Parity check K12 = IIIIIIIZIZIIZ

            ds_Hadamard(reg,13,1);

            //Apply controlled K12

                static const int indicesK12[] = { 7,9,12 };

                for (int i = 0; i != sizeof(indicesK12) / sizeof(indicesK12[0]); ++i) {
                    const int k = indicesK12[i];
                    
                    int control = 13;
                    int target = k;

                    //The CZ:
                    ds_Hadamard(reg, target, 1);
                    ds_cnot(reg, control, target, 1);
                    ds_Hadamard(reg, target, 1);
                }

            ds_Hadamard(reg,13,1);

            syndrome[sweep][11] = ds_measure(reg,1,bit);


            if (syndrome[sweep][11] == 1){ ds_X(reg, 13, 1);}
    }

    // All the stabiliser results are now saved in the matrix 'syndrome'. Let's print them out:



    fprintf(output,"Stabiliser measurements:\n");

    for(int sweep = 0; sweep < numsweeps; sweep++){
        for(int j = 0; j < num_stabs; j++){
            fprintf(output,"%i",syndrome[sweep][j]);
        }
        fprintf(output,"\n");
    }
    
    // The detection event is 1 if corresponding syndrome change from round to round:
    fprintf(output,"Detection events: \n");
    
    bool flag = false; // will flag if any detection event at all (for post-selection purposes)
    for(int i = 0; i < numsweeps - 1; i++){
        for(int j = 0; j < num_stabs; j++){
            change[i][j] = (syndrome[i][j] + syndrome[i+1][j]) %2;
            if (change[i][j] == 1){flag = true;}
            fprintf(output,"%i",change[i][j]);
        }
        fprintf(output,"\n");
    }

    //Parity check ZLogical = ZZZIIIIIIIIII
    //Need to change this to be a classical parity check

    
       ds_Hadamard(reg,13,1);

        //Apply controlled ZL

            for(int j = 0; j < 3; j++){

                int control = 13;
                int target = j;

                //The CZ:
                ds_Hadamard(reg, target, 1);
                ds_cnot(reg, control, target, 1);
                ds_Hadamard(reg, target, 1);
            }

        ds_Hadamard(reg,13,1);

        int ZLresult;

        ZLresult = ds_measure(reg, 1, bit);

        // printf("ZL = ");
        // printf("%d\n",ZLresult);

    // post-select. Flag turned to true if any detection events
    if (flag == false){
        fprintf(output,"No detections\n");
    }
    else 
        fprintf(output,"At least one detection\n");

    // Count logical 0 (yup if logical 0, nope if logical 1)

    if (ZLresult == 0){
        fprintf(output,"|0>_L \n\n"); 
        }
        else fprintf(output,"|1>_L \n\n");

    ds_destroy_register(reg);

    }  // end of for loop of all the samples
    
    fclose(output);
    
    return 0;

}

;
/*-----------------------end-------------------------------*/



    


   

