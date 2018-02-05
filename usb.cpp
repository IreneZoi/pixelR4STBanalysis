// usb.cpp

#ifndef _WIN32
#include <libusb-1.0/libusb.h>
#include <cstring>
#endif

#include "profiler.h"
#include "rpc_error.h"
#include "usb.h"
#include <iostream> // cout
#include <sys/time.h> // gettimeofday, timeval

#include <ftdi.h>
#include <unistd.h> // sleep

const char* CUSB::GetErrorMsg(int error)
{
  switch (error)
    {
    case FT_OK:                          return "ok";
    case FT_INVALID_HANDLE:              return "invalid handle";
    case FT_DEVICE_NOT_FOUND:            return "device not found";
    case FT_DEVICE_NOT_OPENED:           return "device not opened";
    case FT_IO_ERROR:                    return "io error";
    case FT_INSUFFICIENT_RESOURCES:      return "insufficient resource";
    case FT_INVALID_PARAMETER:           return "invalid parameter";
    case FT_INVALID_BAUD_RATE:           return "invalid baud rate";
    case FT_DEVICE_NOT_OPENED_FOR_ERASE: return "device not opened for erase";
    case FT_DEVICE_NOT_OPENED_FOR_WRITE: return "device not opened for write";
    case FT_FAILED_TO_WRITE_DEVICE:      return "failed to write device";
    case FT_EEPROM_READ_FAILED:          return "eeprom read failed";
    case FT_EEPROM_WRITE_FAILED:         return "eeprom write failed";
    case FT_EEPROM_ERASE_FAILED:         return "eeprom erase failed";
    case FT_EEPROM_NOT_PRESENT:          return "eeprom not present";
    case FT_EEPROM_NOT_PROGRAMMED:       return "eeprom not programmed";
    case FT_INVALID_ARGS:                return "invalid args";
    case FT_NOT_SUPPORTED:               return "not supported";
    case FT_OTHER_ERROR:                 return "other error";
    }
  return "unknown error";
}


bool CUSB::EnumFirst(unsigned int &nDevices)
{

  //ftStatus = FT_ListDevices(&enumCount,
  //			    NULL,FT_LIST_NUMBER_ONLY|FT_OPEN_BY_SERIAL_NUMBER);

   enumCount = ftdi_usb_find_all(ftdi, &ftdiDevices, 0x0403, 0x6014);
   if ( 0 >= enumCount )
   {
	  nDevices = enumCount = enumPos = 0;
	  return false;
   }

   nDevices = enumCount;
   enumPos = 0;
   return true;
} // EnumFirst


bool CUSB::EnumNext(char name[])
{
  if (enumPos >= enumCount) return false;

  ftStatus = FT_ListDevices((PVOID)enumPos, name, FT_LIST_BY_INDEX);
  if (ftStatus != FT_OK)
    {
      enumCount = enumPos = 0;
      return false;
    }

  enumPos++;
  return true;
}

bool CUSB::Enum(char name[], unsigned int pos)
{
  enumPos=pos;
  if (enumPos >= enumCount) return false;
  ftStatus = FT_ListDevices((PVOID)enumPos, name, FT_LIST_BY_INDEX);
  if (ftStatus != FT_OK)
    {
      enumCount = enumPos = 0;
      return false;
    }

  return true;
}


bool CUSB::Open(char serialNumber[])
{
  if (isUSB_open) { ftStatus = FT_DEVICE_NOT_OPENED; return false; }

  m_posR = m_sizeR = m_posW = 0;

  //Open the first device with given vendor and product IDs
  ftStatus = ftdi_usb_open(ftdi, 0x0403, 0x6014);
  if (ftStatus != 0) {
	 printf("could not open the specified FTDI device\n");
	 return EXIT_FAILURE;
  }

  //Reset the FTDI chip
  ftStatus = ftdi_usb_reset(ftdi);
  if (ftStatus != 0) {
	 ftdi_usb_close(ftdi);
	 printf("unable to RESET the the chip\n");
	 return EXIT_FAILURE;
  } // if

  //Bitmode Reset
  ftStatus = ftdi_set_bitmode(ftdi, 0xFF, BITMODE_RESET);
  if (ftStatus != 0) {
	 ftdi_usb_close(ftdi);
	 printf("unable to RESET the FTDI bitmode\n");
	 return EXIT_FAILURE;
  }
  usleep(5000);

  //Set the FTDI chip to 245 Synchronous FIFO mode
  ftStatus = ftdi_set_bitmode(ftdi, 0xFF, BITMODE_SYNCFF);
  if (ftStatus != 0) {
	 ftdi_usb_close(ftdi);
	 printf("unable to set the FTDI chip to synchronous FIFO mode\n");
	 return EXIT_FAILURE;
  }

  //Set Latency Timer
  //Time that the FTDI chip keeps data in the internal buffer if buffer is not full yet
  ftStatus = ftdi_set_latency_timer(ftdi, 2);
  if (ftStatus != 0) {
	 ftdi_usb_close(ftdi);
	 printf("unable to set FTDI Latency Timer\n");
	 return EXIT_FAILURE;
  }

  //Set USB read buffer chunk size. Default is 4096
  //libftdi limits the read buffer chunksize to 16384 if a higher value is set
  ftStatus = ftdi_read_data_set_chunksize(ftdi, USBREADBUFFERSIZE);
  if (ftStatus != 0) {
	 ftdi_usb_close(ftdi);
	 printf("unable to set FTDI read chunk size \n");
	 return EXIT_FAILURE;
  }
  ftStatus |= ftdi_write_data_set_chunksize(ftdi, USBWRITEBUFFERSIZE);
  if (ftStatus != 0) {
	 ftdi_usb_close(ftdi);
	 printf("unable to set FTDI write chunk size \n");
	 return EXIT_FAILURE;
  }

  //Set FTDI Flow Control
  ftStatus = ftdi_setflowctrl(ftdi, SIO_RTS_CTS_HS);
  if (ftStatus != 0) {
	 ftdi_usb_close(ftdi);
	 printf("unable to set FTDI Flow Control\n");
	 return EXIT_FAILURE;
  }

  //Clear read and write buffers on the chip and interal read buffer
  ftStatus = ftdi_usb_purge_buffers(ftdi);
  if (ftStatus != 0) {
	 printf("could not purge FTDI buffers\n");
	 return EXIT_FAILURE;
  }

  //Set FTDI timeouts
  //These settings are passsed to the underlying libusb.
  //However, at least READ timeout is not respected. I believe the reason is
  //the the FTDI chip sends 2 status bytes every 40 ms, even in the absense of data.
  ftdi->usb_read_timeout = 4000;
  ftdi->usb_write_timeout = 1000;

  isUSB_open = true;
  return true;
}


void CUSB::Close()
{
  if (!isUSB_open) return;
  ftdi_usb_close(ftdi);
  isUSB_open = 0;
}


void CUSB::Write(const void *buffer, unsigned int bytesToWrite)
{ PROFILING
    if (!isUSB_open) throw CRpcError(CRpcError::WRITE_ERROR);

  DWORD k=0;
  for (k=0; k < bytesToWrite; k++)
    {
      if (m_posW >= USBWRITEBUFFERSIZE) Flush();
      m_bufferW[m_posW++] = ((unsigned char*)buffer)[k];
    }
}


void CUSB::Flush()
{ PROFILING
  DWORD bytesWritten;
  DWORD bytesToWrite = m_posW;
  m_posW = 0;

  if (!isUSB_open) throw CRpcError(CRpcError::WRITE_ERROR);

  if (!bytesToWrite) return;

  bytesWritten = ftdi_write_data(ftdi, m_bufferW, bytesToWrite);

  if (bytesWritten <= 0) throw CRpcError(CRpcError::WRITE_ERROR);
  if (bytesWritten != bytesToWrite) { ftStatus = FT_IO_ERROR; throw CRpcError(CRpcError::WRITE_ERROR); }
}


void CUSB::Read( void *buffer,  unsigned int bytesToRead ) {
  m_sizeR = 0;
  DWORD n = bytesToRead;

  if( n > USBREADBUFFERSIZE )
	n = USBREADBUFFERSIZE; // 4096

timeval tv;
gettimeofday( &tv, NULL );
long s0 = tv.tv_sec;
long u0 = tv.tv_usec;
  while(m_sizeR < bytesToRead) {
    m_posR = 0;

    DWORD bytesRead = ftdi_read_data(ftdi, m_bufferR, n); 

     for( DWORD i = m_sizeR; i < bytesRead; ++i ) {
	   ((unsigned char*)buffer)[i] = m_bufferR[m_posR++];
     } // for i
    m_sizeR += bytesRead;
  }// end while
gettimeofday( &tv, NULL );
long s1 = tv.tv_sec;
long u1 = tv.tv_usec;
//gettimeofday( &tv, NULL );
//long s2 = tv.tv_sec;
//long u2 = tv.tv_usec;
std::cout << "read " << m_sizeR << " in " << s1 - s0 + ( u1 - u0 ) * 1e-6 << " s" << std::endl;
//std::cout << "copied " << bytesRead << " in " << s2 - s1 + ( u2 - u1 ) * 1e-6 << " s" << std::endl;

}

void CUSB::Clear()
{ PROFILING
    if (!isUSB_open) return;

  ftdi_usb_purge_buffers(ftdi);
  m_posR = m_sizeR = 0;
  m_posW = 0;
}


