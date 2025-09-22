//v5 -- added printing to a buffer rather than just printing to file every step to speed it up hopefully (avoid memory allocation bottlenecks)


//v4: made changes to work with the whole RUN_RABI workflow, sending output to results_of_run_rabi > unprocessed_data_from_module etc. like my d2_module_unrotated_reuse.c file.
// (Basically I was trying to write code which reused a single auxiliary but applied idling errors as if it was separate auxiliaries. The results that came out of that looked jumbled, and should now be superceded by my C simulator I wrote with Alan once I get the error PMF's. So now I shall just get a fast and loose complete workflow going, which performs all the stabiliser measurements in serial and uses a single auxiliary qubit (note this assumes long-range connections and creates more idling errors.))

// Also changed Z checks to just use CNOTs going from data qubit to auxiliary in zero state.

// 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include <time.h>
#include "Simulator/sim.h"
#include "Simulator/norm.h"

#define REPS_PER_BUFFER 1000
#define ESTIMATED_REP_SIZE 100   
#define BUFFER_SIZE (REPS_PER_BUFFER * ESTIMATED_REP_SIZE)

int main(int argc, char *argv[]) {
    double p_error, physphi, phystheta, physthetacoef, physphicoef;
    ds_Register reg;
    int distance, numsweeps, reps;

    sscanf(argv[1], "%i", &distance);
    sscanf(argv[2], "%le", &physthetacoef);
    sscanf(argv[3], "%le", &physphicoef);
    sscanf(argv[4], "%le", &p_error);
    sscanf(argv[5], "%i", &numsweeps);
    sscanf(argv[6], "%i", &reps);

    if (distance != 3) return 1;

    phystheta = physthetacoef * M_PI;
    physphi = physphicoef * M_PI;

    int dataqubits = pow(distance, 2) + pow(distance - 1, 2);
    int auxqubits = 1;
    int qubits = dataqubits + auxqubits;
    int num_stabs = dataqubits - 1;
    int auxindex = dataqubits;

    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    unsigned long long nanoseconds = ts.tv_sec * 1000000000LL + ts.tv_nsec;

    char filename[300];
    sprintf(filename, "results_of_run_rabi/unprocessed_data_from_module/d%d/p=%.4f_d%d_reuse_unrot_theta%.3fπ_phi%.3fπ_sweeps%i_reps_%i_id%llu.dat",
            distance, p_error, distance, physthetacoef, physphicoef, numsweeps, reps, nanoseconds * getpid());

    FILE *output = fopen(filename, "a");
    if (!output) { perror("Failed to open file"); return 1; }

    char *buffer = malloc(BUFFER_SIZE);
    if (!buffer) { perror("Failed to allocate buffer"); return 1; }
    size_t pos = 0;

    for (int ii = 0; ii < reps; ii++) {
        char temp[ESTIMATED_REP_SIZE];
        size_t temp_pos = 0;

        temp_pos += snprintf(temp + temp_pos, ESTIMATED_REP_SIZE - temp_pos, "Rep. %i:\n", ii);

        clock_gettime(CLOCK_REALTIME, &ts);
        nanoseconds = ts.tv_sec * 1000000000LL + ts.tv_nsec;
        ds_initialize_simulator(nanoseconds * getpid());
        reg = ds_create_register(qubits, p_error, 0);
        ds_set_state(reg, 0, 1, 0);

        bool faulty_reset_and_measurement = true;
        if (faulty_reset_and_measurement) {
            for (int j = 0; j < qubits; j++) ds_xerr(reg, j);
        }

        bool TI = true;
        if (TI) {
            for (int j = 0; j < dataqubits; j++) ds_yrot(reg, j, -phystheta, 1);
            if (physphicoef != 0) {
                for (int j = 0; j < dataqubits; j++) ds_zrot(reg, j, -physphi, 1);
            }
        }

        int bit[1], syndrome[numsweeps][num_stabs];
        bit[0] = auxindex;

        // ----------------- Stabiliser measurements -----------------
        for (int sweep = 0; sweep < numsweeps; sweep++) {
            // K1 = X X I X I I I I I I I I I
            ds_Hadamard(reg, auxindex, 1);
            int indicesK1[] = {0,1,3};
            for (int i = 0; i < sizeof(indicesK1)/sizeof(indicesK1[0]); i++) ds_cnot(reg, auxindex, indicesK1[i], 1);
            ds_Hadamard(reg, auxindex, 1);
            syndrome[sweep][0] = ds_measure(reg, 1, bit);
            if (syndrome[sweep][0] == 1) ds_X(reg, auxindex, 1);

            // K2 = I X X I X I I I I I I I I
            ds_Hadamard(reg, auxindex, 1);
            int indicesK2[] = {1,2,4};
            for (int i = 0; i < sizeof(indicesK2)/sizeof(indicesK2[0]); i++) ds_cnot(reg, auxindex, indicesK2[i], 1);
            ds_Hadamard(reg, auxindex, 1);
            syndrome[sweep][1] = ds_measure(reg, 1, bit);
            if (syndrome[sweep][1] == 1) ds_X(reg, auxindex, 1);

            // K3 = I I I X I X X I X I I I I
            ds_Hadamard(reg, auxindex, 1);
            int indicesK3[] = {3,5,6,8};
            for (int i = 0; i < sizeof(indicesK3)/sizeof(indicesK3[0]); i++) ds_cnot(reg, auxindex, indicesK3[i], 1);
            ds_Hadamard(reg, auxindex, 1);
            syndrome[sweep][2] = ds_measure(reg, 1, bit);
            if (syndrome[sweep][2] == 1) ds_X(reg, auxindex, 1);

            // K4 = I I I I X I X X I X I I I
            ds_Hadamard(reg, auxindex, 1);
            int indicesK4[] = {4,6,7,9};
            for (int i = 0; i < sizeof(indicesK4)/sizeof(indicesK4[0]); i++) ds_cnot(reg, auxindex, indicesK4[i], 1);
            ds_Hadamard(reg, auxindex, 1);
            syndrome[sweep][3] = ds_measure(reg, 1, bit);
            if (syndrome[sweep][3] == 1) ds_X(reg, auxindex, 1);

            // K5 = I I I I I I I I X I X X I
            ds_Hadamard(reg, auxindex, 1);
            int indicesK5[] = {8,10,11};
            for (int i = 0; i < sizeof(indicesK5)/sizeof(indicesK5[0]); i++) ds_cnot(reg, auxindex, indicesK5[i], 1);
            ds_Hadamard(reg, auxindex, 1);
            syndrome[sweep][4] = ds_measure(reg, 1, bit);
            if (syndrome[sweep][4] == 1) ds_X(reg, auxindex, 1);

            // K6 = I I I I I I I I I X I X X
            ds_Hadamard(reg, auxindex, 1);
            int indicesK6[] = {9,11,12};
            for (int i = 0; i < sizeof(indicesK6)/sizeof(indicesK6[0]); i++) ds_cnot(reg, auxindex, indicesK6[i], 1);
            ds_Hadamard(reg, auxindex, 1);
            syndrome[sweep][5] = ds_measure(reg, 1, bit);
            if (syndrome[sweep][5] == 1) ds_X(reg, auxindex, 1);

            // Z-checks K7-K12
            int indicesK7[] = {0,3,5};
            for (int i = 0; i < sizeof(indicesK7)/sizeof(indicesK7[0]); i++) ds_cnot(reg, indicesK7[i], auxindex, 1);
            syndrome[sweep][6] = ds_measure(reg, 1, bit);
            if (syndrome[sweep][6] == 1) ds_X(reg, auxindex, 1);

            int indicesK8[] = {1,3,4,6};
            for (int i = 0; i < sizeof(indicesK8)/sizeof(indicesK8[0]); i++) ds_cnot(reg, indicesK8[i], auxindex, 1);
            syndrome[sweep][7] = ds_measure(reg, 1, bit);
            if (syndrome[sweep][7] == 1) ds_X(reg, auxindex, 1);

            int indicesK9[] = {2,4,7};
            for (int i = 0; i < sizeof(indicesK9)/sizeof(indicesK9[0]); i++) ds_cnot(reg, indicesK9[i], auxindex, 1);
            syndrome[sweep][8] = ds_measure(reg, 1, bit);
            if (syndrome[sweep][8] == 1) ds_X(reg, auxindex, 1);

            int indicesK10[] = {5,8,10};
            for (int i = 0; i < sizeof(indicesK10)/sizeof(indicesK10[0]); i++) ds_cnot(reg, indicesK10[i], auxindex, 1);
            syndrome[sweep][9] = ds_measure(reg, 1, bit);
            if (syndrome[sweep][9] == 1) ds_X(reg, auxindex, 1);

            int indicesK11[] = {6,8,9,11};
            for (int i = 0; i < sizeof(indicesK11)/sizeof(indicesK11[0]); i++) ds_cnot(reg, indicesK11[i], auxindex, 1);
            syndrome[sweep][10] = ds_measure(reg, 1, bit);
            if (syndrome[sweep][10] == 1) ds_X(reg, auxindex, 1);

            int indicesK12[] = {7,9,12};
            for (int i = 0; i < sizeof(indicesK12)/sizeof(indicesK12[0]); i++) ds_cnot(reg, indicesK12[i], auxindex, 1);
            syndrome[sweep][11] = ds_measure(reg, 1, bit);
            if (syndrome[sweep][11] == 1) ds_X(reg, auxindex, 1);
        }

        // ----------------- Writing to buffer -----------------
        temp_pos += snprintf(temp + temp_pos, ESTIMATED_REP_SIZE - temp_pos, "Stabiliser measurements:\n");
        for (int sweep = 0; sweep < numsweeps; sweep++) {
            for (int x = 0; x < 2*distance*(distance-1); x += distance+(distance-1)) {
                for (int xi = 0; xi < distance-1; xi++) {
                    int index = x+xi;
                    if (index < num_stabs) temp_pos += snprintf(temp + temp_pos, ESTIMATED_REP_SIZE - temp_pos, "%i", syndrome[sweep][index]);
                }
            }
            for (int z = distance-1; z < 2*distance*(distance-1); z += distance+(distance-1)) {
                for (int zi = 0; zi < distance; zi++) {
                    int index = z+zi;
                    if (index < num_stabs) temp_pos += snprintf(temp + temp_pos, ESTIMATED_REP_SIZE - temp_pos, "%i", syndrome[sweep][index]);
                }
            }
            temp_pos += snprintf(temp + temp_pos, ESTIMATED_REP_SIZE - temp_pos, "\n");
        }

        temp_pos += snprintf(temp + temp_pos, ESTIMATED_REP_SIZE - temp_pos, "Data qubit measurements:\n");
        for (int j = 0; j < distance; j++) {
            bit[0] = j;
            int measurement = ds_measure(reg, 1, bit);
            temp_pos += snprintf(temp + temp_pos, ESTIMATED_REP_SIZE - temp_pos, "%d", measurement);
        }
        temp_pos += snprintf(temp + temp_pos, ESTIMATED_REP_SIZE - temp_pos, "\n\n");

        ds_destroy_register(reg);

        if (pos + temp_pos >= BUFFER_SIZE) {
            fwrite(buffer, 1, pos, output);
            pos = 0;
        }
        memcpy(buffer + pos, temp, temp_pos);
        pos += temp_pos;
    }

    if (pos > 0) fwrite(buffer, 1, pos, output);

    free(buffer);
    fclose(output);

    return 0;
}
