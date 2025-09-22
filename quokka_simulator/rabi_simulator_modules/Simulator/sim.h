#ifndef SIM_H
#define SIM_H

#include <stdio.h>

typedef struct {
   double x,y ;
} ds_Complex;

typedef struct {
   ds_Complex *state;
   int *steps;
   int nq, nc;
   int *n_errors;
   double err, sigma;
} ds_Register;

extern double ds_Pi, ds_Pio2, ds_Pio4, ds_Pio8, ds_root2_2;
extern ds_Register ds_reg;

ds_Register ds_create_register(int nq_L, double err_L, double sigma_L);
void ds_destroy_register(ds_Register reg);
void ds_equate_registers(ds_Register reg1, ds_Register reg2);
void ds_initialize_simulator(long ds_seed);
void ds_clearreg(ds_Register reg);
void ds_set_state(ds_Register reg, int n, double x, double y);
int ds_query_state(ds_Register reg, int n, double tol);
void ds_print(ds_Register reg);
void ds_update(ds_Register reg, int q1, int q2);
void ds_global_update(ds_Register reg);

// Ant's additions [


void ds_idling_during_TI(ds_Register reg, int aux_index);
void ds_Hadamard_or_idling_on_aux(ds_Register reg,int auxindex,int stab_type,int idling_errors);
void ds_preceding_data_idling(ds_Register reg, int data_qubit_index, int num_preceding_data_idling_errors,int idling_errors);
void ds_data_idling(ds_Register reg, int data_qubit_index, int num_idling_errors,int idling_errors);

void ds_my_CNOT(ds_Register reg, int stab_type, int k, int auxindex);
void ds_aux_idling(ds_Register reg,int auxindex, int aux_idling, int idling_errors);


void ds_xerr(ds_Register reg, int q);
void ds_depolarising(ds_Register reg, int q);
void ds_print_reg(ds_Register reg);

int ds_get_qubit_step(ds_Register reg, int qubit_index);
void ds_print_all_qubit_steps(ds_Register reg);
void ds_print_qubit_step(ds_Register reg, int qubit_index); //Ant: prints qubit's reg.steps
void ds_set_qubit_step(ds_Register reg, int qubit_index, int step_value);

// ] 

ds_Complex ds_eitheta(double theta);
ds_Complex ds_add(ds_Complex z1, ds_Complex z2);
ds_Complex ds_multiply(ds_Complex z1, ds_Complex z2);
ds_Complex ds_zstarz(ds_Complex z1, ds_Complex z2);
double ds_modsq(ds_Complex z);
double ds_inner_product(ds_Register reg1, ds_Register reg2);

void ds_esigx(ds_Complex *zPtr1, ds_Complex *zPtr2, double theta);
void ds_esigy(ds_Complex *zPtr1, ds_Complex *zPtr2, double theta);
void ds_esigz(ds_Complex *zPtr1, ds_Complex *zPtr2, double theta);
void ds_unitary(ds_Complex *zPtr1, ds_Complex *zPtr2,
             double alpha, double beta, double theta);

void ds_one_qubit_indices(int n, int q, int *iPtr, int *jPtr);
void ds_controlled_indices(int n, int q1, int q2, int *iPtr, int *jPtr);
void ds_swap_indices(int n, int q1, int q2, int *iPtr, int *jPtr);
void ds_two_qubit_indices(int n, int q1, int q2, int *i1Ptr, int *i2Ptr, int *i3Ptr, int*i4Ptr);

void ds_xrot(ds_Register reg, int q, double theta, int time);
void ds_yrot(ds_Register reg, int q, double theta, int time);
void ds_zrot(ds_Register reg, int q, double theta, int time);

void ds_X(ds_Register reg, int q, int time);
void ds_Z(ds_Register reg, int q, int time);
void ds_XZ(ds_Register reg, int q, int time);

void ds_cnot(ds_Register reg, int qcont, int qtarg, int time);
void ds_swap(ds_Register reg, int q1, int q2, int time);
void ds_lerr(ds_Register reg, int q, int time);
void ds_cphase(ds_Register reg, int qcont, int qtarg, double theta, int time);
void ds_Hadamard(ds_Register reg, int q, int time);

int ds_measure(ds_Register reg, int nq2m, int *lq2m);

double ds_set_measure(ds_Register reg, int nq2m, int *lq2m, int val);

#endif
