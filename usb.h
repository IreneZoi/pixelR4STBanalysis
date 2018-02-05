// usb.h
//
// Author: Aliakbar Ebrahimi, UHH
//
// Class provides basic functionalities for the USB interface using libFTDI1

#ifndef USB_H
#define USB_H

#ifdef _WIN32
#include "windows.h"
#include "win32\FTD2XX.h"
#else
#include "linux/ftd2xx.h"
#endif

#include "rpc_io.h"

#include <ftdi.h>

#define USBWRITEBUFFERSIZE  65536
#define USBREADBUFFERSIZE  200*65536 

class CUSB : public CRpcIo
{
  // libftdi1 stuff
  struct ftdi_context *ftdi;
  struct ftdi_device_list *ftdiDevices;

  bool isUSB_open;
  FT_HANDLE ftHandle;
  FT_STATUS ftStatus;

  DWORD enumPos, enumCount;

  DWORD m_posW;
  unsigned char m_bufferW[USBWRITEBUFFERSIZE];

  DWORD m_posR, m_sizeR;
  unsigned char m_bufferR[USBREADBUFFERSIZE];

 public:
  CUSB()
    {
      m_posR = m_sizeR = m_posW = 0;
      isUSB_open = false;
      ftHandle = 0; ftStatus = 0;
      enumPos = enumCount = 0;
      // libFTDI stuff
      ftdi = ftdi_new();
    }
  ~CUSB() { /* Close(); */ }
  int GetLastError() { return ftStatus; }
  static const char* GetErrorMsg(int error);
  bool EnumFirst(unsigned int &nDevices);
  bool EnumNext(char name[]);
  bool Enum(char name[], unsigned int pos);
  bool Open(char serialNumber[]);
  void Close();
  bool Connected() { return isUSB_open; };


  void Write(const void *buffer, unsigned int size);
  void Flush();
  void Clear();
  void Read(void *buffer, unsigned int size);
};

#endif
