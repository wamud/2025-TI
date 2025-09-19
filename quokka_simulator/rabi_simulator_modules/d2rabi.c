
// This seems to be an analytical rabi simulator in that it sets the stabiliser measurement reusults (projects the state) rather than randoly projects

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


int main(int argc, char *argv[]){  
    double p_error, angle, X1, X2, Z1, Z2; 
    ds_Register reg; 
    int n; 
    n = 6;
    sscanf(argv[1],"%lf",&angle);
    sscanf(argv[2],"%lf",&p_error);
    sscanf(argv[3],"%lf",&X1);
    sscanf(argv[4],"%lf",&X2);
    sscanf(argv[5],"%lf",&Z1);
    sscanf(argv[6],"%lf",&Z2);
   

    ds_initialize_simulator((unsigned) time(NULL)*getpid());
    
    reg = ds_create_register(n, p_error, 0);  
    ds_set_state(reg, 0, 1, 0); 


    // Do transversal injection - i.e. rotate all the qubits
        for(int j = 0; j < (n-1); j++){
            ds_yrot(reg, j, -angle, 0); //-ve angle as yrot opposite?
        }

    //Parity checks. I.e. controlled (by ancilla) stabiliser nested in hadamards on ancilla

    int bit[1];
    bit[0] = 5; //This is the ancilla qubit which will be measured

    // Parity check K1 = XXXII

        ds_Hadamard(reg,5,0);

        //Apply controlled XXXII

            for(int j = 0; j < 3; j++){
                ds_cnot(reg,5,j,1); //ds_cnot(ds_Register reg, int qcont, int qtarg, int time);
            }

        ds_Hadamard(reg,5,0);

        // This is setting the measurement
        ds_set_measure(reg, 1, bit, X1); //double ds_set_measure(ds_Register reg, int nq2m, int *lq2m, int val)
        
        //Reset ancilla if was measured to one:
        if (X1 == 1){ ds_X(reg, 5, 1);}
    
    
    // Parity check K2 = IIXXX 

        ds_Hadamard(reg,5,0);

        //Apply controlled IIXXX
            
            for(int j = 2; j < 5; j++){
                ds_cnot(reg,5,j,1); 
            }

        ds_Hadamard(reg,5,0);

        ds_set_measure(reg, 1, bit, X2);

        if (X2 == 1){ ds_X(reg, 5, 1);}


    // Parity check K3 = ZIZZI

        ds_Hadamard(reg,5,0);

        //Apply controlled ZIZZI

            static const int indicesK3[] = { 0,2,3 };

            for (int i = 0; i != sizeof(indicesK3) / sizeof(indicesK3[0]); ++i) {
                const int k = indicesK3[i];
                
                int control = 5;
                int target = k;

                //The CZ:
                ds_Hadamard(reg, target, 0);
                ds_cnot(reg, control, target, 1);
                ds_Hadamard(reg, target, 0);
            }

        ds_Hadamard(reg,5,0);

        ds_set_measure(reg, 1, bit, Z1);

        if (Z1 == 1){ ds_X(reg, 5, 1);}
    
    //Parity check K4 = IZZIZ

       ds_Hadamard(reg,5,0);

        //Apply controlled IZZIZ

            static const int indicesK4[] = { 1,2,4 };

            for (int i = 0; i != sizeof(indicesK4) / sizeof(indicesK4[0]); ++i) {
                const int k = indicesK4[i];
                
                int control = 5;
                int target = k;

                //The CZ:
                ds_Hadamard(reg, target, 0);
                ds_cnot(reg, control, target, 1);
                ds_Hadamard(reg, target, 0);
            }

        ds_Hadamard(reg,5,0);

        ds_set_measure(reg, 1, bit, Z2);

        if (Z2 == 1){ ds_X(reg, 5, 1);}

    //Parity check ZL = ZZIII

    
       ds_Hadamard(reg,5,0);

        //Apply controlled ZZIII

            for(int j = 0; j < 2; j++){

                int control = 5;
                int target = j;

                //The CZ:
                ds_Hadamard(reg, target, 0);
                ds_cnot(reg, control, target, 1);
                ds_Hadamard(reg, target, 0);
            }

        ds_Hadamard(reg,5,0);

        int result;

        result = ds_measure(reg, 1, bit);

    // Count logical 0:

    if (result == 0){
        printf("yup, logical zero \n");
        }
        
        else { printf("nope \n");
    }
    

    ds_destroy_register(reg);
    
    return 0;

}

;
/*-----------------------end-------------------------------*/



    


   

