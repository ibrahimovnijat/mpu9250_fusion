#ifndef KALMAN_H_
#define KALMAN_H_

#include "arm_math.h"
#include "arm_const_structs.h"

#include <malloc.h>
#include <stdlib.h>

#define dt_kalman   0.0025f		// kalman filter sampling time
#define Pi			 3.14159f
#define deg2rad		 Pi / 180.0f
#define rad2deg 	 180.0f / Pi


arm_matrix_instance_f32 A;
arm_matrix_instance_f32 A_trans;
arm_matrix_instance_f32 B;
arm_matrix_instance_f32 H;
arm_matrix_instance_f32 H_trans;
arm_matrix_instance_f32 U;
arm_matrix_instance_f32 Z;

arm_matrix_instance_f32 P;
arm_matrix_instance_f32 P_prev;
arm_matrix_instance_f32 X;
arm_matrix_instance_f32 X_hat;
arm_matrix_instance_f32 X_hat_prev;

arm_matrix_instance_f32 ident_mat;

arm_matrix_instance_f32 K;
arm_matrix_instance_f32 Q;
arm_matrix_instance_f32 R;

//initialize temperary matrices to be used in the filter recursively to hold temperary results
arm_matrix_instance_f32 temp_mat1;
arm_matrix_instance_f32 temp_mat2;
arm_matrix_instance_f32 temp_mat3;
arm_matrix_instance_f32 temp_mat4;
arm_matrix_instance_f32 temp_mat5;

float32_t* null_arr1;
float32_t* null_arr2;
float32_t* null_arr3;
float32_t* null_arr4;
float32_t* null_arr5;


void initKalman(void);

void printMat(arm_matrix_instance_f32 mat);

void kalman_filter(float32_t gx, float32_t gy, float32_t roll_accel, float32_t pitch_accel, float32_t* angles);

void measure_accelAngles(float* roll_out, float* pitch_out, float* accel_data);

void test_func(void);

#endif /* KALMAN_H_ */
