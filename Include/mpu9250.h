#ifndef MPU9250_H_
#define MPU9250_H_

#include "stm32f4xx_hal.h"
#include <stdbool.h>
#include <stdint.h>
#include "serialCommProtocols.h"

/*------------- MPU defines ----------------*/
#define I2C_MST_CTRL                0x24
#define I2C_SLV0_ADDR               0x25
#define I2C_SLV0_REG                0x26
#define I2C_SLV0_CTRL               0x27
#define I2C_SLV0_DO                 0x63

#define EXT_SENS_DATA_00            0x49
#define EXT_SENS_DATA_01            0x4A
#define EXT_SENS_DATA_02            0x4B
#define EXT_SENS_DATA_03            0x4C
#define EXT_SENS_DATA_04            0x4D
#define EXT_SENS_DATA_05            0x4E
#define EXT_SENS_DATA_06            0x4F
#define EXT_SENS_DATA_07            0x50
#define EXT_SENS_DATA_08            0x51

#define MPU_SLAVE_ADDR              0xD0      			// MPU I2C address (AD0 is low) but 0xD2 if AD0 is high
#define I2C_DATA_READ               0x01      			// From MPU to micro-controller
#define I2C_MAG_DATA_READ           0x80        		// From slave (Magnetometer) to MPU
/*----------------------------------------------*/
#define PWR_MGMT_1                  0x6B
#define PWR_MGMT_2            		 0x6C
#define CONFIG                   	 0x1A
#define SAMPLE_DIVIDE_RATE          0x19
#define GYRO_CONFIG              	 0x1B
#define ACCEL_CONFIG                0x1C
#define ACCEL_CONFIG_2              0x1D
#define INT_PIN_CONFIG              0x37
#define INT_ENABLE                  0x38
#define INT_STATUS                  0x3A
#define MOT_DETECT_CTRL             0x69
#define WOM_THRESH                  0x1F
#define LP_ACCEL_ODR                0x1E

#define FIFO_EN                		 0x23
#define FIFO_COUNTH             	 0x72
#define FIFO_COUNTL             	 0x73
#define FIFO_R_W                 	 0x74

#define ACCEL_XOUT_H                0x3B
#define ACCEL_YOUT_H                0x3D
#define ACCEL_ZOUT_H                0x3F
#define GYRO_XOUT_H                 0x43
#define GYRO_YOUT_H             	 0x45
#define GYRO_ZOUT_H             	 0x47
#define TEMP_OUT_H                  0x41

#define ACCEL_XOUT_L                0x3C
#define ACCEL_YOUT_L                0x3E
#define ACCEL_ZOUT_L                0x40
#define GYRO_XOUT_L             	 0x44
#define GYRO_YOUT_L             	 0x46
#define GYRO_ZOUT_L             	 0x48
#define TEMP_OUT_L                  0x42

#define USER_CTRL                   0x6A
#define SIGNAL_PATH_RESET           0x68
#define I2C_MST_CTRL                0x24
#define MPU_WHO_AM_I                0x75
/*---------------------------------------------*/
#define ACCEL_SENSITIVITY_2G        16384.0f
#define ACCEL_SENSITIVITY_4G        8192.0f
#define ACCEL_SENSITIVITY_8G        4096.0f
#define ACCEL_SENSITIVITY_16G       2048.0f

#define GYRO_SENSITIVITY_250        131.0f
#define GYRO_SENSITIVITY_500        65.5f
#define GYRO_SENSITIVITY_1000       32.8f
#define GYRO_SENSITIVITY_2000       16.4f

/*--------------- Variable declarations -------------------*/
extern volatile float accelerometerData[3], gyroscopeData[3];
extern float accelerometerBias[3], gyroscopeBias[3]; // bias values of IMU connected to SPI
extern float prev_accelData[3], prev_gyroData[3]; 	// for gyro and accel low pass filter
extern float totalAcc[3], totalGyr[3];

extern SPI_HandleTypeDef hspi2;
extern I2C_HandleTypeDef hi2c3;

/*-------------- Function declarations -------------------*/
void whoAmI_MPU9250_I2C(void);
void whoAmI_MPU9250_SPI(void);

void initMPU9250_I2C(void);
void initMPU9250_SPI(void);

bool isIMUDataReady(void);

void readAccelerometerData_I2C(float* accelerometerData);
void readGyroscopeData_I2C(float* gyroscopeData);

void readAccelerometerData_SPI(float* accelerometerData);
void readGyroscopeData_SPI(float* gyroscopeData);

void readAccelerometerData_SPI_calibr(float* accelerometerData);
void readGyroscopeData_SPI_calibr(float* gyroscopeData);

void mpuCalibration_I2C(float* accelerometerBias, float* gyroscopeData);
void mpuCalibration_SPI(float* accelerometerBias, float* gyroscopeData);

#endif /* MPU9250_H_ */
