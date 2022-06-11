/*
 * kalman.c
 *
 *  Created on: Oct 23, 2018
 *      Author: Nijat
 */


#include <stdio.h>
#include "kalman.h"

static float32_t X_arr[2] = {0, 0};
static float32_t X_hat_arr[2] = {0, 0};
static float32_t X_hat_prev_arr[2] = {0, 0};

static float32_t A_arr[4] = {1, 0, 0, 1};
static float32_t B_arr[4] = {dt_kalman, 0, 0, dt_kalman};
static float32_t H_arr[4] = {1, 0, 0, 1};
static float32_t A_trans_arr[4] = {1, 0, 0, 1};
static float32_t H_trans_arr[4] = {1, 0, 0, 1};
static float32_t ident_mat_arr[4] = {1.0, 0.0, 0.0, 1.0};

static float32_t Q_arr[4] = {0.00001, 0.0, 0.0, 0.00001};
static float32_t R_arr[4] = {0.1, 0.0, 0.0, 0.1};
static float32_t P_arr[4] = {0.5, 0.0, 0.0, 0.5};
static float32_t P_prev_arr[4] = {0.5, 0.0, 0.0, 0.5};
static float32_t K_arr[4] = {0.1, 0.1, 0.1, 0.1};



void printMat(arm_matrix_instance_f32 mat)
{
	for (uint16_t i = 0; i < mat.numRows; i++){
		for (uint16_t j = 0; j < mat.numCols; j++){
			printf("%.1f ", *mat.pData);
			mat.pData++;
		}
	}
	printf("\n");
}

void initKalman(void)
{
	arm_mat_init_f32(&A, 2, 2, (float32_t*)A_arr);
	arm_mat_init_f32(&B, 2, 2, (float32_t*)B_arr);
	arm_mat_init_f32(&H, 2, 2, (float32_t*)H_arr);
	arm_mat_init_f32(&A_trans, 2, 2, (float32_t*)A_trans_arr);
	arm_mat_init_f32(&H_trans, 2, 2, (float32_t*)H_trans_arr);

	arm_mat_init_f32(&X, 2, 1, (float32_t*)X_arr);
	arm_mat_init_f32(&X_hat, 2, 1, (float32_t*)X_hat_arr);
	arm_mat_init_f32(&X_hat_prev, 2, 1, (float32_t*)X_hat_prev_arr);

	arm_mat_init_f32(&P, 2, 2, (float32_t*)P_arr);
	arm_mat_init_f32(&P_prev, 2, 2, (float32_t*)P_prev_arr);
	arm_mat_init_f32(&Q, 2, 2, (float32_t*)Q_arr);
	arm_mat_init_f32(&R, 2, 2, (float32_t*)R_arr);
	arm_mat_init_f32(&K, 2, 2, (float32_t*)K_arr);

	arm_mat_init_f32(&ident_mat, 2, 2, (float32_t*)ident_mat_arr);
}
//

void kalman_filter(float32_t gx, float32_t gy, float32_t roll_accel, float32_t pitch_accel, float32_t* angles)
{
	/*----------------------------- TIME UPDATE ---------------------------------*/

	// Step 1. Calculate the priori state matrix, initial state estimate from system dynamics
	null_arr1 = (float32_t*)calloc(2, sizeof(float32_t));
	null_arr2 = (float32_t*)calloc(2, sizeof(float32_t));

	arm_mat_init_f32(&temp_mat1, 2, 1, null_arr1);
	arm_mat_init_f32(&temp_mat2, 2, 1, null_arr2);
	arm_mat_mult_f32(&A, &X_hat_prev, &temp_mat1);
	float32_t u_arr[2] = {gx, gy};

	arm_mat_init_f32(&U, 2, 1, (float32_t*)u_arr);
	arm_mat_mult_f32(&B, &U, &temp_mat2);

	arm_mat_add_f32(&temp_mat1, &temp_mat2, &X);	//estimate priori states
	free(null_arr1); free(null_arr2);

	//Step 2. Calculate error covariance matrix
	null_arr1 = (float32_t*)calloc(4, sizeof(float32_t));
	null_arr2 = (float32_t*)calloc(4, sizeof(float32_t));

	arm_mat_init_f32(&temp_mat1, 2, 2, null_arr1);
	arm_mat_init_f32(&temp_mat2, 2, 2, null_arr2);

	arm_mat_mult_f32(&A, &P_prev, &temp_mat1);
	//arm_mat_trans_f32(&A, &temp_mat2);
	arm_mat_mult_f32(&temp_mat1, &A_trans, &temp_mat2);
	arm_mat_add_f32(&temp_mat2, &Q, &P);
	free(null_arr1); free(null_arr2);

	/*---------------------------- MEASUREMENT UPDATE ------------------------------*/
	//Step 3. Calculate Kalman gain
	null_arr1 = (float32_t*)calloc(4, sizeof(float32_t));
	null_arr2 = (float32_t*)calloc(4, sizeof(float32_t));
	null_arr3 =	(float32_t*)calloc(4, sizeof(float32_t));
	null_arr4 = (float32_t*)calloc(4, sizeof(float32_t));
	null_arr5 = (float32_t*)calloc(4, sizeof(float32_t));

	arm_mat_init_f32(&temp_mat1, 2, 2, null_arr1);
	arm_mat_init_f32(&temp_mat2, 2, 2, null_arr2);
	arm_mat_init_f32(&temp_mat3, 2, 2, null_arr3);
	arm_mat_init_f32(&temp_mat4, 2, 2, null_arr4);
	arm_mat_init_f32(&temp_mat5, 2, 2, null_arr5);

	arm_mat_mult_f32(&P, &H_trans, &temp_mat1);

	arm_mat_mult_f32(&H, &P, &temp_mat2);
	arm_mat_mult_f32(&temp_mat2, &H_trans, &temp_mat3);
	arm_mat_add_f32(&temp_mat3, &R, &temp_mat4);
	arm_mat_inverse_f32(&temp_mat4, &temp_mat5);
	arm_mat_mult_f32(&temp_mat1, &temp_mat5, &K);
	free(null_arr1); free(null_arr2); free(null_arr3); free(null_arr4); free(null_arr5);

	//Step 4. Calculate posteriori state matrix
	float32_t out_arr[2] = {roll_accel, pitch_accel};
	arm_mat_init_f32(&Z, 2, 1, (float32_t*)out_arr);

	null_arr1 = (float32_t*)calloc(2, sizeof(float32_t));
	null_arr2 = (float32_t*)calloc(2, sizeof(float32_t));
	null_arr3 = (float32_t*)calloc(2, sizeof(float32_t));

	arm_mat_init_f32(&temp_mat1, 2, 1, null_arr1);
	arm_mat_init_f32(&temp_mat2, 2, 1, null_arr2);
	arm_mat_init_f32(&temp_mat3, 2, 1, null_arr3);

	arm_mat_mult_f32(&H, &X, &temp_mat1);
	arm_mat_sub_f32(&Z, &temp_mat1, &temp_mat2);
	arm_mat_mult_f32(&K, &temp_mat2, &temp_mat3);
	arm_mat_add_f32(&X, &temp_mat3, &X_hat);
	free(null_arr1); free(null_arr2); free(null_arr3);

	//Step 5. Update error covariance matrix
	null_arr1 = (float32_t*)calloc(4, sizeof(float32_t));
	null_arr2 = (float32_t*)calloc(4, sizeof(float32_t));

	arm_mat_init_f32(&temp_mat1, 2, 2, null_arr1);
	arm_mat_init_f32(&temp_mat2, 2, 2, null_arr2);

	arm_mat_mult_f32(&K, &H, &temp_mat1);
	arm_mat_sub_f32(&ident_mat, &temp_mat1, &temp_mat2);
	arm_mat_mult_f32(&temp_mat2, &P, &P_prev);
	free(null_arr1); free(null_arr2);

	//update state vector
	arm_copy_f32(X_hat.pData, X_hat_prev.pData, 2);
	angles[0] = *X_hat.pData;
	X_hat.pData++;
	angles[1] = *X_hat.pData;
	X_hat.pData--;
}

void measure_accelAngles(float32_t* roll_out, float32_t* pitch_out, float32_t* accel_data)
{
	*roll_out  = (float32_t) atan2(accel_data[1], (float32_t) sqrt(accel_data[0]*accel_data[0] + accel_data[2]*accel_data[2])) * rad2deg;
	*pitch_out = (float32_t) atan2(accel_data[0], (float32_t) sqrt(accel_data[1]*accel_data[1] + accel_data[2]*accel_data[2])) * rad2deg;
}

