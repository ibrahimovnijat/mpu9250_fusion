#include "mpu9250.h"



void whoAmI_MPU9250_I2C(void)
{
	uint8_t whoAmI_regValue;
    I2C_Receive(&hi2c3, MPU_SLAVE_ADDR, MPU_WHO_AM_I, &whoAmI_regValue, 1);
    if (whoAmI_regValue == 0x71 || whoAmI_regValue == 0x73)
    {
    	printf("Connected to IMU..");
	    printf("MPU address is 0x%X ..\n", whoAmI_regValue);
    }
    else
    {
    	printf("Not connected to IMU, Please, reset the device.. \n");
    	while(1);
    }
}

void whoAmI_MPU9250_SPI(void)
{
	uint8_t whoAmI_regValue;
    SPI_Receive(&hspi2, MPU_WHO_AM_I, &whoAmI_regValue, 1);
    if (whoAmI_regValue == 0x71 || whoAmI_regValue == 0x73)
    {
    	printf("Connected to IMU..");
    	printf("MPU address is 0x%X ..\n", whoAmI_regValue);
    }
    else
    {
    	printf("Not connected to IMU, Please, reset the device.. \n");
    	while(1);
    }
}

void initMPU9250_I2C(void)
{
	uint8_t data[2] = {0x00};
	uint8_t temp_read_data = {0x00};

	//select the power source
	data[0] = PWR_MGMT_1; data[1] = 0x80; //reset MPU
	I2C_Transmit(&hi2c3, MPU_SLAVE_ADDR, (uint8_t*)&data, 2);
	HAL_Delay(500);
	// select best available clock source - PLL if ready
	data[1] = 0x01;
	I2C_Transmit(&hi2c3, MPU_SLAVE_ADDR, (uint8_t*)&data, 2);
	HAL_Delay(10);
	//enable accelerometer and gyroscope
	data[0] = PWR_MGMT_2; data[1] = 0x00;
	I2C_Transmit(&hi2c3, MPU_SLAVE_ADDR, (uint8_t*)&data, 2);
	HAL_Delay(10);

	// internal sampling rate, bandwitdh(cutoff) freq,
	data[0] = CONFIG; data[1] = 0x01; // 184MHz cutoff freq, 1KHz internal sampling, 2.9 ms delay
	I2C_Transmit(&hi2c3, MPU_SLAVE_ADDR, (uint8_t*)&data, 2);
	HAL_Delay(10);

    data[0] = SAMPLE_DIVIDE_RATE; data[1] = 0x00;	// internal sampling is the same as ODR freq, 1Khz sampling rate
    I2C_Transmit(&hi2c3, MPU_SLAVE_ADDR, (uint8_t*)&data, 2);
    HAL_Delay(10);

	//configre gyroscope at +=500dps rate
	data[0] = GYRO_CONFIG;
	I2C_Receive(&hi2c3, MPU_SLAVE_ADDR, GYRO_CONFIG, &temp_read_data, 1);
	data[1] = temp_read_data & ~0x02; //clear fchoice (DLPF) bits (to enable sampling rate config)
	data[1] = data[1] & ~0x18; // clear Gyro FS Scale bits
	data[1] = data[1] | 0x08; //clear gyro scale bits -> set to +- 500dps
	I2C_Transmit(&hi2c3, MPU_SLAVE_ADDR, (uint8_t*)&data, 2);
	HAL_Delay(100);

	//configure accelerometer
	data[0] = ACCEL_CONFIG;
	I2C_Receive(&hi2c3, MPU_SLAVE_ADDR, ACCEL_CONFIG, &temp_read_data, 1);
	data[1] = temp_read_data & ~0x18;  //clear Accel full scale bits
	data[1]  = data[1] | 0x08;		// set to +- 4G
    I2C_Transmit(&hi2c3, MPU_SLAVE_ADDR, (uint8_t*)&data, 2);
    HAL_Delay(100);

    //configure accelerometer sampling rate
    data[0] = ACCEL_CONFIG_2;
    I2C_Receive(&hi2c3, MPU_SLAVE_ADDR, ACCEL_CONFIG_2, &temp_read_data, 1);
    data[1] = temp_read_data & ~0x08;  //clear accel fchoice (DLPS) bit
    data[1] = data[1] | 0x01;  		// 184Mhz cutoff freq, 1KHz internal sampling, 5.8 ms delay
    I2C_Transmit(&hi2c3, MPU_SLAVE_ADDR, (uint8_t*)&data, 2);
    HAL_Delay(100);

    // configure data ready interrupt
    data[0] = INT_PIN_CONFIG; //data[1] = 0x22;
    I2C_Receive(&hi2c3, MPU_SLAVE_ADDR, ACCEL_CONFIG_2, &temp_read_data, 1);
    data[1] = temp_read_data & ~ 0x01;	// clear interrupt after INT_STATUS is read
    I2C_Transmit(&hi2c3, MPU_SLAVE_ADDR, (uint8_t*)&data, 2);

    data[0] = INT_ENABLE; data[1] = 0x01;	//enable data ready interrupt
    I2C_Transmit(&hi2c3, MPU_SLAVE_ADDR, (uint8_t*)&data, 2);
    HAL_Delay(10);
}

void initMPU9250_SPI(void)
{
	uint8_t data_spi = 0x00;
	uint8_t temp_read_data = {0x00};

	//select the power source
	data_spi = 0x80; //reset MPU
	SPI_Transmit(&hspi2, PWR_MGMT_1, data_spi, 1);
	HAL_Delay(500);

	// select best available clock source - PLL if ready
	data_spi = 0x01;
	SPI_Transmit(&hspi2, PWR_MGMT_1, data_spi, 1);
	HAL_Delay(10);

	//enable accelerometer and gyroscope
	data_spi = 0x00;
	SPI_Transmit(&hspi2, PWR_MGMT_2, data_spi, 1 );
	HAL_Delay(10);

	// internal sampling rate, bandwitdh(cutoff) freq,
	data_spi = 0x01; // // 184MHz cutoff freq, 1KHz internal sampling, 2.9 ms delay
	SPI_Transmit(&hspi2, CONFIG, data_spi, 1 );
	HAL_Delay(10);

	data_spi = 0x00;	// internal sampling is the same as ODR freq, 1Khz sampling rate
	SPI_Transmit(&hspi2, SAMPLE_DIVIDE_RATE, data_spi, 1 );
    HAL_Delay(10);

	//configure gyroscope at +=500dps rate
    SPI_Receive( &hspi2, GYRO_CONFIG, &temp_read_data, 1 );
    data_spi = temp_read_data & ~0x02; //clear fchoice (DLPF) bits (to enable sampling rate config)
    data_spi = data_spi & ~0x18; // clear Gyro FS Scale bits
    //data_spi = data_spi | 0x08;//clear gyro scale bits -> set to +- 500dps
    data_spi = data_spi | 0x10; //set gyro range to +- 1000dps
	SPI_Transmit(&hspi2, GYRO_CONFIG, data_spi, 1 );
	HAL_Delay(100);

	//configure accelerometer
    SPI_Receive( &hspi2, ACCEL_CONFIG, &temp_read_data, 1 );
    data_spi = temp_read_data & ~0x18; //clear accel full scale bits
    //data_spi = data_spi | 0x08; // set to +- 4G
    data_spi = data_spi | 0x10; //set accel range to +- 8G
	SPI_Transmit(&hspi2, ACCEL_CONFIG, data_spi, 1 );
    HAL_Delay(100);

    //configure accelerometer sampling rate
    SPI_Receive( &hspi2, ACCEL_CONFIG_2, &temp_read_data, 1 );
    data_spi = temp_read_data & ~0x08;  //clear accel fchoice (DLPS) bit
    data_spi = data_spi | 0x01;  // 184Mhz cutoff freq, 1KHz internal sampling, 5.8 ms delay
	SPI_Transmit(&hspi2, ACCEL_CONFIG_2, data_spi, 1 );
    HAL_Delay(100);
}

bool isIMUDataReady(void)
{
	uint8_t mpuDataReadyStatus = 0x00;
	I2C_Receive(&hi2c3, MPU_SLAVE_ADDR, INT_STATUS, (uint8_t*)&mpuDataReadyStatus, 1); //if new mpu
	return (mpuDataReadyStatus & 0x01);
}

void readAccelerometerData_I2C(float* accelerometerData)
{
	uint8_t raw_accel_data[6] = {0};
	int16_t accel[3] = {0};
	uint8_t baseAddress = ACCEL_XOUT_H;

	I2C_Receive(&hi2c3, MPU_SLAVE_ADDR, baseAddress, (uint8_t*)&raw_accel_data, 6);

	accel[0] = (int16_t) ((raw_accel_data[0] << 8) | raw_accel_data[1]);
	accel[1] = (int16_t) ((raw_accel_data[2] << 8) | raw_accel_data[3]);
	accel[2] = (int16_t) ((raw_accel_data[4] << 8) | raw_accel_data[5]);

	accelerometerData[0] = (float) accel[0] / ACCEL_SENSITIVITY_8G;
	accelerometerData[1] = (float) accel[1] / ACCEL_SENSITIVITY_8G;
	accelerometerData[2] = (float) accel[2] / ACCEL_SENSITIVITY_8G;

	accelerometerData[0] -= accelerometerBias[0];
	accelerometerData[1] -= accelerometerBias[1];
	accelerometerData[2] -= accelerometerBias[2];

	//convert x and z directions to match North-East-Down (NED) coordinate system
	accelerometerData[0] = -accelerometerData[0];
	accelerometerData[2] = -accelerometerData[2];

}

void readGyroscopeData_I2C(float* gyroscopeData)
{
	uint8_t raw_gyro_data[6] = {0};
	uint8_t baseAddress = GYRO_XOUT_H;
	int16_t gyro[3] = {0};

	I2C_Receive(&hi2c3, MPU_SLAVE_ADDR, baseAddress, (uint8_t*)&raw_gyro_data, 6);

	gyro[0] = (int16_t) ((raw_gyro_data[0] << 8) | raw_gyro_data[1]);
	gyro[1] = (int16_t) ((raw_gyro_data[2] << 8) | raw_gyro_data[3]);
	gyro[2] = (int16_t) ((raw_gyro_data[4] << 8) | raw_gyro_data[5]);

	gyroscopeData[0] = (float) gyro[0] / GYRO_SENSITIVITY_1000;
	gyroscopeData[1] = (float) gyro[1] / GYRO_SENSITIVITY_1000;
	gyroscopeData[2] = (float) gyro[2] / GYRO_SENSITIVITY_1000;

	gyroscopeData[0] -= gyroscopeBias[0];
	gyroscopeData[1] -= gyroscopeBias[1];
	gyroscopeData[2] -= gyroscopeBias[2];

	//convert x and z directions to match North-East-Down (NED) coordinate system
	gyroscopeData[0] = -gyroscopeData[0];
	gyroscopeData[2] = -gyroscopeData[2];

}

void readAccelerometerData_SPI_calibr(float* accelerometerData)
{
	uint8_t raw_accel_data[6] = {0};
	int16_t accel[3] = {0};
	uint8_t baseAddress = ACCEL_XOUT_H;

	SPI_Receive(&hspi2, baseAddress, (uint8_t*)&raw_accel_data, 6);
	accel[0] = (int16_t) ((raw_accel_data[0] << 8) | raw_accel_data[1]);
	accel[1] = (int16_t) ((raw_accel_data[2] << 8) | raw_accel_data[3]);
	accel[2] = (int16_t) ((raw_accel_data[4] << 8) | raw_accel_data[5]);

	accelerometerData[0] = (float) accel[0] / ACCEL_SENSITIVITY_8G;
	accelerometerData[1] = (float) accel[1] / ACCEL_SENSITIVITY_8G;
	accelerometerData[2] = (float) accel[2] / ACCEL_SENSITIVITY_8G;
}

void readAccelerometerData_SPI(float* accelerometerData)
{
	uint8_t raw_accel_data[6] = {0};
	int16_t accel[3] = {0};
	uint8_t baseAddress = ACCEL_XOUT_H;

	SPI_Receive(&hspi2, baseAddress, (uint8_t*)&raw_accel_data, 6);

	accel[0] = (int16_t) ((raw_accel_data[0] << 8) | raw_accel_data[1]);
	accel[1] = (int16_t) ((raw_accel_data[2] << 8) | raw_accel_data[3]);
	accel[2] = (int16_t) ((raw_accel_data[4] << 8) | raw_accel_data[5]);

	accelerometerData[0] = (float) accel[0] / ACCEL_SENSITIVITY_8G;
	accelerometerData[1] = (float) accel[1] / ACCEL_SENSITIVITY_8G;
	accelerometerData[2] = (float) accel[2] / ACCEL_SENSITIVITY_8G;

	accelerometerData[0] -= accelerometerBias[0];
	accelerometerData[1] -= accelerometerBias[1];
	accelerometerData[2] -= accelerometerBias[2];
	//convert x and z directions to match North-East-Down (NED) coordinate system
	accelerometerData[0] = -accelerometerData[0];
	accelerometerData[2] = -accelerometerData[2];
}

void readGyroscopeData_SPI_calibr(float* gyroscopeData)
{
	uint8_t raw_gyro_data[6] = {0};
	uint8_t baseAddress = GYRO_XOUT_H;
	int16_t gyro[3] = {0};

	SPI_Receive(&hspi2, baseAddress, (uint8_t*)&raw_gyro_data, 6);

	gyro[0] = (int16_t) ((raw_gyro_data[0] << 8) | raw_gyro_data[1]);
	gyro[1] = (int16_t) ((raw_gyro_data[2] << 8) | raw_gyro_data[3]);
	gyro[2] = (int16_t) ((raw_gyro_data[4] << 8) | raw_gyro_data[5]);

	gyroscopeData[0] = (float) gyro[0] / GYRO_SENSITIVITY_1000;
	gyroscopeData[1] = (float) gyro[1] / GYRO_SENSITIVITY_1000;
	gyroscopeData[2] = (float) gyro[2] / GYRO_SENSITIVITY_1000;
}

void readGyroscopeData_SPI(float* gyroscopeData)
{
	uint8_t raw_gyro_data[6] = {0};
	uint8_t baseAddress = GYRO_XOUT_H;
	int16_t gyro[3] = {0};

	SPI_Receive(&hspi2, baseAddress, (uint8_t*)&raw_gyro_data, 6);

	gyro[0] = (int16_t) ((raw_gyro_data[0] << 8) | raw_gyro_data[1]);
	gyro[1] = (int16_t) ((raw_gyro_data[2] << 8) | raw_gyro_data[3]);
	gyro[2] = (int16_t) ((raw_gyro_data[4] << 8) | raw_gyro_data[5]);

	gyroscopeData[0] = (float) gyro[0] / GYRO_SENSITIVITY_1000;
	gyroscopeData[1] = (float) gyro[1] / GYRO_SENSITIVITY_1000;
	gyroscopeData[2] = (float) gyro[2] / GYRO_SENSITIVITY_1000;

	gyroscopeData[0] -= gyroscopeBias[0];
	gyroscopeData[1] -= gyroscopeBias[1];
	gyroscopeData[2] -= gyroscopeBias[2];

	//convert x and z directions to match North-East-Down (NED) coordinate system
	gyroscopeData[0] = -gyroscopeData[0];
	gyroscopeData[2] = -gyroscopeData[2];
}


void mpuCalibration_I2C(float* accelerometerBias, float* gyroscopeBias)
{
	float tempAccelData[3] = {0.0f}, tempGyroData[3] = {0.0f};
	float tempAccelDataSum[3] = {0.0f}, tempGyroDataSum[3] = {0.0f};
	int bufferSize = 0;

	printf("IMU calibration...\nWaiting for temperature stabilization..\n");
	HAL_Delay(5000);
	printf("Calibrating accelerometer. Please, wait...\n");

	do
	{
		if (isIMUDataReady())
		{
			readAccelerometerData_I2C((float*)&tempAccelData);
			tempAccelDataSum[0] += tempAccelData[0];
			tempAccelDataSum[1] += tempAccelData[1];
			tempAccelDataSum[2] += tempAccelData[2];
			bufferSize += 1;
			HAL_Delay(7);
		}
	}
	while(bufferSize < 2000);

	for (int i = 0; i < 3; i++)
	{
		accelerometerBias[i] = (float) tempAccelDataSum[i] / bufferSize;	//calculate the mean of data
	}
	accelerometerBias[2] -= 1.0f;
	printf("Accelerometer calibration done...\n");
	HAL_Delay(2000);
	bufferSize = 0;
	printf("Calibrating gyroscope. Please, wait...\n");
	do
	{
		if(isIMUDataReady())
		{
			readGyroscopeData_I2C((float*)&tempGyroData);
			tempGyroDataSum[0] += tempGyroData[0];
			tempGyroDataSum[1] += tempGyroData[1];
			tempGyroDataSum[2] += tempGyroData[2];
			bufferSize += 1;
			HAL_Delay(7);
		}
	}
	while (bufferSize < 2000);

	for (int i = 0; i < 3; i++)
	{
		gyroscopeBias[i] = (float) tempGyroDataSum[i] / bufferSize;
	}

	printf("Gyroscope calibration done..\nAcceleromter bias:\n");
	printf("%.5f  %.5f  %.5f\n", accelerometerBias[0], accelerometerBias[1], accelerometerBias[2]);
	printf("Gyroscope bias: \n");
	printf("%.5f  %.5f  %.5f\n", gyroscopeBias[0], gyroscopeBias[1], gyroscopeBias[2]);
}

void mpuCalibration_SPI(float* accelerometerBias, float* gyroscopeBias)
{
	float tempAccelData[3] = {0.0f}, tempGyroData[3] = {0.0f};
	float tempAccelDataSum[3] = {0.0f}, tempGyroDataSum[3] = {0.0f};
	uint16_t bufferSize = 0;

	printf("IMU calibration...\nWaiting for temperature stabilization..\n");
	HAL_Delay(500);
	printf("Calibrating accelerometer. Please, wait...\n");

	do
	{
		readAccelerometerData_SPI_calibr((float*)&tempAccelData);
		for (uint8_t i = 0; i < 3; i++)
		{
			tempAccelDataSum[i] += tempAccelData[i];
		}
		bufferSize += 1;
		HAL_Delay(8);
	}
	while (bufferSize < 500);

	for(uint8_t i = 0; i < 3; i++)
	{
		accelerometerBias[i] = (float) (tempAccelDataSum[i] / (float)bufferSize);
	}
	accelerometerBias[2] -= 1.0f;
	printf("Accelerometer calibration done...\nCalibrating gyroscope. Please, wait...\n");
	HAL_Delay(500);
	bufferSize = 0;

	do
	{
		readGyroscopeData_SPI_calibr((float*)&tempGyroData);
		for (uint8_t i = 0; i < 3; i++)
		{
			tempGyroDataSum[i] += tempGyroData[i];
		}
		bufferSize += 1;
		HAL_Delay(8);
	}
	while(bufferSize < 500);

	for (uint8_t i = 0; i < 3; i++)
	{
		gyroscopeBias[i] = (float)(tempGyroDataSum[i] / (float) bufferSize);
	}

	printf("Gyroscope calibration done..\nAcceleromter bias:\n");
	printf("%.5f  %.5f  %.5f\n Gyroscope bias: \n", accelerometerBias[0], accelerometerBias[1], accelerometerBias[2]);
	printf("%.5f  %.5f  %.5f\n", gyroscopeBias[0], gyroscopeBias[1], gyroscopeBias[2]);
}
