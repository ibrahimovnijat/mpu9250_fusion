#include "serialCommProtocols.h"


void I2C_Transmit(I2C_HandleTypeDef* hi2c, uint8_t DeviceAddr, uint8_t* WriteData, uint16_t size)
{
   if (HAL_I2C_Master_Transmit(hi2c, (uint16_t)DeviceAddr, WriteData, size, (uint32_t)HAL_MAX_DELAY) != HAL_OK)
   {
      printf("I2C master transmit error \n");
   }
}

void I2C_Receive(I2C_HandleTypeDef* hi2c, uint8_t DeviceAddr, uint8_t ReadAddr, uint8_t* ReadData, uint16_t size)
{
   if (HAL_I2C_Master_Transmit(hi2c, (uint16_t)DeviceAddr, (uint8_t*)&ReadAddr, 1, (uint32_t)HAL_MAX_DELAY) != HAL_OK)
   {
      printf("I2C master transmit error (receive) \n");
   }

   DeviceAddr = DeviceAddr | I2C_DATA_READ;

   if (HAL_I2C_Master_Receive(hi2c, (uint16_t)DeviceAddr, ReadData, size, (uint32_t)HAL_MAX_DELAY) != HAL_OK)
   {
      printf("I2C master receive error \n");
   }
}

void SPI_Receive( SPI_HandleTypeDef* hspi, uint8_t ReadAddr, uint8_t* ReadData, uint16_t size)
{
	CS_ON;
	ReadAddr = ReadAddr | SPI_DATA_READ;
    HAL_SPI_Transmit(hspi, &ReadAddr, 1, (uint32_t)HAL_MAX_DELAY);
    HAL_SPI_Receive(hspi, ReadData, size, (uint32_t)HAL_MAX_DELAY);
	CS_OFF;
}

void SPI_Transmit(SPI_HandleTypeDef* hspi,  uint8_t WriteReg, uint8_t WriteData, uint16_t size )
{
	CS_ON;
    HAL_SPI_Transmit(hspi,  &WriteReg, 1, (uint32_t)HAL_MAX_DELAY);
    HAL_SPI_Transmit(hspi,  &WriteData, size, (uint32_t)HAL_MAX_DELAY);
    CS_OFF;
}
