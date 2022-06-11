
#ifndef SERIALCOMMPROTOCOLS_H_
#define SERIALCOMMPROTOCOLS_H_

#include "stm32f4xx_hal.h"

#define I2C_DATA_READ 	0x01
#define SPI_DATA_READ	0x80
#define CS_ON		    { HAL_GPIO_WritePin(GPIOB, GPIO_PIN_12, GPIO_PIN_RESET); }
#define CS_OFF			{ HAL_GPIO_WritePin(GPIOB, GPIO_PIN_12, GPIO_PIN_SET);   }

/*----------------- Variable declarations ---------------------*/

/*----------------- Function declarations ---------------------*/
void I2C_Transmit(I2C_HandleTypeDef* hi2c, uint8_t DeviceAddr, uint8_t* WriteData, uint16_t size);
void I2C_Receive(I2C_HandleTypeDef* hi2c, uint8_t DeviceAddr, uint8_t ReadAddr, uint8_t* ReadData, uint16_t size);

void SPI_Receive( SPI_HandleTypeDef* hspi, uint8_t ReadAddr, uint8_t* ReadData, uint16_t size);
void SPI_Transmit(SPI_HandleTypeDef* hspi,  uint8_t WriteReg, uint8_t WriteData, uint16_t size);


#endif /* SERIALCOMMPROTOCOLS_H_ */
