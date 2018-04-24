/*---------------------------------------------------------------------
 *
 *  filename:    pixel_dtb.h
 *
 *  description: PSI46 testboard API for DTB
 *	author:      Beat Meier
 *	date:        1.12.2015
 *	rev:
 *
 *---------------------------------------------------------------------
 */

#pragma once

#include "profiler.h"
#include "rpc.h"
#include "config.h"

#ifdef _WIN32
#include "pipe.h"
#endif

#include "usb.h"


class CTestboard
{
  RPC_DEFS
    RPC_THREAD

#ifdef _WIN32
    CPipeClient pipe;
#endif
  CUSB usb;

 public:
  CRpcIo& GetIo() { return *rpc_io; }

  CTestboard() { RPC_INIT rpc_io = &usb; }
  ~CTestboard() { RPC_EXIT }

  // === RPC ==============================================================

  // Don't change the following two entries
  RPC_EXPORT uint16_t GetRpcVersion();
  RPC_EXPORT int32_t  GetRpcCallId(string &callName);

  RPC_EXPORT void GetRpcTimestamp(stringR &ts);

  RPC_EXPORT int32_t GetRpcCallCount();
  RPC_EXPORT bool    GetRpcCallName(int32_t id, stringR &callName);

  bool RpcLink(bool verbose = true);

  // === DTB connection ====================================================

  bool EnumFirst(unsigned int &nDevices) { return usb.EnumFirst(nDevices); };
  bool EnumNext(string &name);
  bool Enum(unsigned int pos, string &name);

  bool FindDTB(string &usbId);
  bool Open(string &name, bool init=true); // opens a connection
  void Close();				// closes the connection to the testboard

#ifdef _WIN32
  bool OpenPipe(const char *name) { return pipe.Open(name); }
  void ClosePipe() { pipe.Close(); }
#endif

  bool IsConnected() { return usb.Connected(); }
  const char * ConnectionError()
  { return usb.GetErrorMsg(usb.GetLastError()); }

  void Flush() { rpc_io->Flush(); }
  void Clear() { rpc_io->Clear(); }

  // === DTB identification ===============================================

  RPC_EXPORT void GetInfo(stringR &info);
  RPC_EXPORT uint16_t GetBoardId();	// reads the board number

  RPC_EXPORT void GetHWVersion(stringR &version);
  RPC_EXPORT uint16_t GetFWVersion();
  RPC_EXPORT uint16_t GetSWVersion();

  // === DTB service ======================================================

  // --- upgrade
  RPC_EXPORT uint16_t UpgradeGetVersion();
  RPC_EXPORT uint8_t  UpgradeStart(uint16_t version);
  RPC_EXPORT uint8_t  UpgradeData(string &record);
  RPC_EXPORT uint8_t  UpgradeError();
  RPC_EXPORT void     UpgradeErrorMsg(stringR &msg);
  RPC_EXPORT void     UpgradeExec(uint16_t recordCount);

  // === DTB functions ====================================================

  RPC_EXPORT void Init();

  RPC_EXPORT void Welcome();
  RPC_EXPORT void SetLed(uint8_t x);

  // === Clock, Timing ====================================================

  RPC_EXPORT void cDelay(uint16_t clocks);
  RPC_EXPORT void uDelay(uint16_t us);
  void mDelay(uint16_t ms);

  // --- select ROC/Module clock source

#define CLK_SRC_INT  0
#define CLK_SRC_EXT  1
  RPC_EXPORT void SetClockSource(uint8_t source);

  // --- check if external clock is present

  RPC_EXPORT bool IsClockPresent();

  // --- Signal Delay -----------------------------------------------------

#define SIG_CLK 0
#define SIG_CTR 1
#define SIG_SDA 2
#define SIG_TIN 3

  RPC_EXPORT void Sig_SetLevel(uint8_t signal, uint8_t level);
  RPC_EXPORT void Sig_SetOffset(uint8_t offset);
  RPC_EXPORT void Sig_SetLVDS();
  RPC_EXPORT void Sig_SetLCDS();
  RPC_EXPORT void Sig_SetRdaToutDelay(uint8_t delay);

  // --- digital signal probe ---------------------------------------------

#define PROBE_OFF             0
#define PROBE_ADC_DAVAIL     12
#define PROBE_HOLD           11
#define PROBE_CAL_PULSE      10
#define PROBE_CAL_ENA         9
#define PROBE_LA              8
#define PROBE_READOUT_LINE    7
#define PROBE_FIRST_LINE      6
#define PROBE_MEASURE         5
#define PROBE_SEQ_END         4
#define PROBE_SEQ_START       3
#define PROBE_SEQ_ABORT       2
#define PROBE_SEQ_RUNNING     1

  RPC_EXPORT void SignalProbeD1(uint8_t signal);
  RPC_EXPORT void SignalProbeD2(uint8_t signal);

  // --- analog signal probe ----------------------------------------------

#define PROBEA_TIN     0
#define PROBEA_SDATA1  1
#define PROBEA_SDATA2  2
#define PROBEA_CTR     3
#define PROBEA_CLK     4
#define PROBEA_SDA     5
#define PROBEA_TOUT    6
#define PROBEA_OFF     7

#define GAIN_1   0
#define GAIN_2   1
#define GAIN_3   2
#define GAIN_4   3

  RPC_EXPORT void SignalProbeA1(uint8_t signal);
  RPC_EXPORT void SignalProbeA2(uint8_t signal);
  RPC_EXPORT void SignalProbeADC(uint8_t signal, uint8_t gain = 0);

  // --- ROC/Module power VD/VA -------------------------------------------

  RPC_EXPORT void Pon();
  RPC_EXPORT void Poff();

  RPC_EXPORT void _SetVD(uint16_t mV);
  RPC_EXPORT void _SetVA(uint16_t mV);
  RPC_EXPORT void _SetID(uint16_t uA100);
  RPC_EXPORT void _SetIA(uint16_t uA100);

  RPC_EXPORT uint16_t _GetVD();
  RPC_EXPORT uint16_t _GetVA();
  RPC_EXPORT uint16_t _GetID();
  RPC_EXPORT uint16_t _GetIA();

  //void SetVA(double V) { _SetVA(uint16_t(V*1000)); }  // set VA voltage
  //void SetVD(double V) { _SetVD(uint16_t(V*1000)); }  // set VD voltage
  //void SetIA(double A) { _SetIA(uint16_t(A*10000)); }  // set VA current limit
  //void SetID(double A) { _SetID(uint16_t(A*10000)); }  // set VD current limit

  double GetVA() { return _GetVA()/1000.0; }   // get VA voltage in V
  double GetVD() { return _GetVD()/1000.0; }	 // get VD voltage in V
  double GetIA() { return _GetIA()/10000.0; }  // get VA current in A
  double GetID() { return _GetID()/10000.0; }  // get VD current in A

  RPC_EXPORT uint16_t _GetVD_Cap();
  RPC_EXPORT uint16_t _GetVD_Reg();
  RPC_EXPORT uint16_t _GetVDAC_Reg();

  double GetVD_Cap(){ return _GetVD_Cap()/1000.0;}
  double GetVD_Reg(){ return _GetVD_Reg()/1000.0;}
  double GetVDAC_Reg(){ return _GetVDAC_Reg()/1000.0;}


  RPC_EXPORT void HVon();
  RPC_EXPORT void HVoff();

  RPC_EXPORT uint8_t GetStatus();

  // --- data aquisition --------------------------------------------------

  RPC_EXPORT uint32_t Daq_Open(uint32_t buffersize = 10000000);
  RPC_EXPORT void Daq_Close();
  RPC_EXPORT void Daq_Start();
  RPC_EXPORT void Daq_Stop();
  RPC_EXPORT void Daq_MemReset();
  RPC_EXPORT uint32_t Daq_GetSize();

  // --- Read data from DAQ FIFO
  // Meaning if return bits
#define DAQ_FIFO_OVFL 4
#define DAQ_MEM_OVFL  2
#define DAQ_STOPPED   1

  RPC_EXPORT uint8_t Daq_Read(HWvectorR<uint16_t> &data,
			      uint32_t blocksize = 65536);

  RPC_EXPORT uint8_t Daq_Read(HWvectorR<uint16_t> &data,
			      uint32_t blocksize, uint32_t &availsize);

  // Ethernet test functions

  RPC_EXPORT void Ethernet_Send(string &message);
  RPC_EXPORT uint32_t Ethernet_RecvPackets();

  // --- ROC4Sens adapter test functions ----------------------------------

  RPC_EXPORT void r4s_SetHoldPos(uint8_t t);
  RPC_EXPORT void r4s_AdcDelay(uint8_t t);
  RPC_EXPORT void r4s_SetSequence(vector<uint32_t> &prog);

  void r4s_SetSeqReadout(int ext = 0); // FW 0.8
  void r4s_SetSeqReadCol(int ext = 0); // FW 0.7, 0.8
  void r4s_SetSeqCalScan(int ext = 0); // row-wise
  void r4s_SetSeqCalScanCol(int ext = 0); // Column-wise

  void r4s_SetPixCal(uint8_t x, uint8_t y);
  void r4s_Set2PixCal(uint8_t x, uint8_t y);
  void r4s_Set4PixCal(uint8_t x, uint8_t y);

  RPC_EXPORT void r4s_SetRegX(vector<uint32_t> &shr_x);
  RPC_EXPORT void r4s_SetRegY(vector<uint32_t> &shr_y);

  RPC_EXPORT void r4s_Enable(uint8_t t);
  RPC_EXPORT void r4s_Start();
  RPC_EXPORT void r4s_Stop();
  RPC_EXPORT bool r4s_Running();

  RPC_EXPORT void r4s_SetVdig(uint16_t mV);
  RPC_EXPORT void r4s_SetVana(uint16_t mV);
  RPC_EXPORT void r4s_SetVcal(uint16_t mV);
  RPC_EXPORT void r4s_SetRgsh(uint16_t mV);
  RPC_EXPORT void r4s_SetRgpr(uint16_t mV);
  RPC_EXPORT void r4s_SetVref(uint16_t mV);
  RPC_EXPORT void r4s_SetVaux1(uint16_t mV);
  RPC_EXPORT void r4s_SetVaux2(uint16_t mV);

};
