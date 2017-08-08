/* -------------------------------------------------------------
 *
 *  file:        command.cpp
 *
 *  description: command line interpreter for Chip/Wafer tester
 *
 *  author:      Beat Meier
 *  modified:    31.8.2007
 *
 *  rev:
 *
 * -------------------------------------------------------------
 */

#include "cmd.h"

// =======================================================================
//  connection, communication, startup commands
// =======================================================================

CMD_PROC(scan)
{
  CTestboard *tb = new CTestboard;
  string name;
  vector<string> devList;
  unsigned int nDev, nr;

  try
    {
      if( !tb->EnumFirst(nDev)) throw int(1);
      for( nr=0; nr<nDev; nr++)
	{
	  if( !tb->EnumNext(name)) throw int(2);
	  if( name.size() < 4) continue;
	  if( name.compare(0, 4, "DTB_") == 0) devList.push_back(name);
	}
    }
  catch (int e)
    {
      switch (e)
	{
	case 1: printf("Cannot access the USB driver\n"); break;
	case 2: printf("Cannot read name of connected device\n"); break;
	}
      delete tb;
      return;
    }

  if( devList.size() == 0)
    {
      printf("no DTB connected\n");
      return;
    }

  for( nr = 0; nr < devList.size(); nr++)
    try
      {
	printf("%10s: ", devList[nr].c_str());
	if( !tb->Open(devList[nr],false))
	  {
	    printf("DTB in use\n");
	    continue;
	  }

	unsigned int bid = tb->GetBoardId();
	printf("DTB Id %u\n", bid);
	tb->Close();
      }
    catch (...)
      {
	printf("DTB not identifiable\n");
	tb->Close();
      }

  delete tb;
}

CMD_PROC(open)
{
  if( tb.IsConnected())
    {
      printf("Already connected to DTB.\n");
      return;
    }

  string usbId;
  char name[80];
  if( PAR_IS_STRING(name,79)) usbId = name;
  else if( !tb.FindDTB(usbId)) return;

  bool status = tb.Open(usbId, false);

  if( !status)
    {
      printf("USB error: %s\nCould not connect to DTB %s\n", tb.ConnectionError(), usbId.c_str());
      return;
    }
  printf("DTB %s opened\n", usbId.c_str());

  string info;
  tb.GetInfo(info);
  printf("--- DTB info-------------------------------------\n"
	 "%s"
	 "-------------------------------------------------\n", info.c_str());
}

CMD_PROC(close)
{
  tb.Close();
}


CMD_PROC(rpclink)
{
  if( tb.RpcLink()) printf("ok\n");
}


CMD_PROC(welcome)
{
  tb.Welcome();
  DO_FLUSH
    }


CMD_PROC(setled)
{
  int value;
  PAR_INT(value, 0, 0x3f);
  tb.SetLed(value);
  DO_FLUSH
    }

CMD_PROC(log)
{
  char s[256];
  PAR_STRINGEOL(s,255);
  Log.printf("%s\n", s);
}

bool UpdateDTB(const char *filename)
{
  fstream src;

  if( tb.UpgradeGetVersion() == 0x0100)
    {
      // open file
      src.open(filename);
      if( !src.is_open())
	{
	  printf("ERROR UPGRADE: Could not open \"%s\"!\n", filename);
	  return false;
	}

      // check if upgrade is possible
      printf("Start upgrading DTB.\n");
      if( tb.UpgradeStart(0x0100) != 0)
	{
	  string msg;
	  tb.UpgradeErrorMsg(msg);
	  printf("ERROR UPGRADE: %s!\n", msg.data());
	  return false;
	}

      // download data
      printf("Download running ...\n");
      string rec;
      uint16_t recordCount = 0;
      while (true)
	{
	  getline(src, rec);
	  if( src.good())
	    {
	      if( rec.size() == 0) continue;
	      recordCount++;
	      if( tb.UpgradeData(rec) != 0)
		{
		  string msg;
		  tb.UpgradeErrorMsg(msg);
		  printf("ERROR UPGRADE: %s!\n", msg.data());
		  return false;
		}
	    }
	  else if( src.eof()) break;
	  else
	    {
	      printf("ERROR UPGRADE: Error reading \"%s\"!\n", filename);
	      return false;
	    }
	}

      if( tb.UpgradeError() != 0)
	{
	  string msg;
	  tb.UpgradeErrorMsg(msg);
	  printf("ERROR UPGRADE: %s!\n", msg.data());
	  return false;
	}

      // write EPCS FLASH
      printf("DTB download complete.\n");
      tb.mDelay(200);
      printf("FLASH write start (LED 1..4 on)\n"
	     "DO NOT INTERUPT DTB POWER !\n"
	     "Wait till LEDs goes off.\n"
	     "Restart the DTB.\n");
      tb.UpgradeExec(recordCount);
      tb.Flush();

      return true;
    }

  printf("ERROR UPGRADE: Could not upgrade this DTB version!\n");
  return false;
}


CMD_PROC(upgrade)
{
  char filename[256];
  PAR_STRING(filename, 255);
  UpdateDTB(filename);
}

CMD_PROC(rpcinfo)
{
  string name, call, ts;

  tb.GetRpcTimestamp(ts);
  int version = tb.GetRpcVersion();
  int n = tb.GetRpcCallCount();

  printf("--- DTB RPC info ----------------------------------------\n");
  printf("RPC version:     %i.%i\n", version/256, version & 0xff);
  printf("RPC timestamp:   %s\n", ts.c_str());
  printf("Number of calls: %i\n", n);
  printf("Function calls:\n");
  for( int i = 0; i < n; i++)
    {
      tb.GetRpcCallName(i, name);
      rpc_TranslateCallName(name, call);
      printf("%5i: %s\n", i, call.c_str());
    }
}

CMD_PROC(info)
{
  string s;
  tb.GetInfo(s);
  printf("--- DTB info ------------------------------------\n%s"
	 "-------------------------------------------------\n", s.c_str());
}

CMD_PROC(ver)
{
  string hw;
  tb.GetHWVersion(hw);
  int fw = tb.GetFWVersion();
  int sw = tb.GetSWVersion();
  printf("%s: FW=%i.%02i SW=%i.%02i\n", hw.c_str(), fw/256, fw%256, sw/256, sw%256);
}

CMD_PROC(version)
{
  string hw;
  tb.GetHWVersion(hw);
  int fw = tb.GetFWVersion();
  int sw = tb.GetSWVersion();
  printf("%s: FW=%i.%02i SW=%i.%02i\n", hw.c_str(), fw/256, fw%256, sw/256, sw%256);
}

CMD_PROC(boardid)
{
  int id = tb.GetBoardId();
  printf("\nBoard Id = %i\n", id);
}

CMD_PROC(init)
{
  tb.Init();
  DO_FLUSH
    }

CMD_PROC(flush)
{
  tb.Flush();
}

CMD_PROC(clear)
{
  tb.Clear();
}




// =======================================================================
//  delay commands
// =======================================================================

CMD_PROC(udelay)
{
  int del;
  PAR_INT(del, 0, 1000);
  if( del) tb.uDelay(del);
  DO_FLUSH
    }

CMD_PROC(mdelay)
{
  int ms;
  PAR_INT(ms,1,10000)
    tb.mDelay(ms);
}


// =======================================================================
//  test board commands
// =======================================================================

CMD_PROC(clksrc)
{
  int source;
  PAR_INT(source, 0, 1);
  tb.SetClockSource(source);
  DO_FLUSH
    }

CMD_PROC(clkok)
{
  if( tb.IsClockPresent())
    printf("clock ok\n");
  else
    printf("clock missing\n");
}


CMD_PROC(clklvl)
{
  int lvl;
  PAR_INT(lvl,0,15);
  tb.Sig_SetLevel(SIG_CLK, lvl);
  DO_FLUSH
    }

CMD_PROC(sdalvl)
{
  int lvl;
  PAR_INT(lvl,0,15);
  tb.Sig_SetLevel(SIG_SDA, lvl);
  DO_FLUSH
    }

CMD_PROC(ctrlvl)
{
  int lvl;
  PAR_INT(lvl,0,15);
  tb.Sig_SetLevel(SIG_CTR, lvl);
  DO_FLUSH
    }

CMD_PROC(tinlvl)
{
  int lvl;
  PAR_INT(lvl,0,15);
  tb.Sig_SetLevel(SIG_TIN, lvl);
  DO_FLUSH
    }

CMD_PROC(lvds)
{
  tb.Sig_SetLVDS();
  DO_FLUSH
    }

CMD_PROC(lcds)
{
  tb.Sig_SetLCDS();
  DO_FLUSH
    }


CMD_PROC(pon)
{
  tb.Pon();
  DO_FLUSH
    }

CMD_PROC(poff)
{
  tb.Poff();
  DO_FLUSH
    }

CMD_PROC(va)
{
  int value;
  PAR_INT(value, 0, 4000);
  tb._SetVA(value);
  DO_FLUSH
    }

CMD_PROC(vd)
{
  int value;
  PAR_INT(value, 0, 4000);
  tb._SetVD(value);
  DO_FLUSH
    }


CMD_PROC(ia)
{
  int value;
  PAR_INT(value, 0, 1200);
  tb._SetIA(value*10);
  DO_FLUSH
    }

CMD_PROC(id)
{
  int value;
  PAR_INT(value, 0, 1200);
  tb._SetID(value*10);
  DO_FLUSH
    }

CMD_PROC(getva)
{
  double v = tb.GetVA();
  printf("\n VA = %1.3fV\n", v);
}

CMD_PROC(getvd)
{
  double v = tb.GetVD();
  printf("\n VD = %1.3fV\n", v);
}

CMD_PROC(getia)
{
  double i = tb.GetIA();
  printf("\n IA = %1.1fmA\n", i*1000.0);
}

CMD_PROC(getid)
{
  double i = tb.GetID();
  printf("\n ID = %1.1fmA\n", i*1000.0);
}


CMD_PROC(hvon)
{
  tb.HVon();
  DO_FLUSH
    }

CMD_PROC(hvoff)
{
  tb.HVoff();
  DO_FLUSH
    }


CMD_PROC(status)
{
  uint8_t status = tb.GetStatus();
  printf("SD card detect: %c\n", (status&8) ? '1' : '0');
  printf("CRC error:      %c\n", (status&4) ? '1' : '0');
  printf("Clock good:     %c\n", (status&2) ? '1' : '0');
  printf("CLock present:  %c\n", (status&1) ? '1' : '0');
}

CMD_PROC(d1)
{
  int sig;
  PAR_INT(sig, 0, 255)
    tb.SignalProbeD1(sig);
  DO_FLUSH
    }

CMD_PROC(d2)
{
  int sig;
  PAR_INT(sig, 0, 255)
    tb.SignalProbeD2(sig);
  DO_FLUSH
    }

CMD_PROC(a1)
{
  int sig;
  PAR_INT(sig, 0, 7)
    tb.SignalProbeA1(sig);
  DO_FLUSH
    }

CMD_PROC(a2)
{
  int sig;
  PAR_INT(sig, 0, 7)
    tb.SignalProbeA2(sig);
  DO_FLUSH
    }

CMD_PROC(probeadc)
{
  int sig;
  PAR_INT(sig, 0, 7)
    tb.SignalProbeADC(sig);
  DO_FLUSH
    }


// === DAQ ==================================================================

CMD_PROC(dopen)
{
  int buffersize;
  PAR_INT(buffersize, 0, 60000000);
  buffersize = tb.Daq_Open(buffersize);
  printf("%i words allocated for data buffer %i\n", buffersize);
  if( buffersize == 0) printf("error\n");
}

CMD_PROC(dclose)
{
  tb.Daq_Close();
  DO_FLUSH
    }


CMD_PROC(dstart)
{
  tb.Daq_Start();
  DO_FLUSH
    }

CMD_PROC(dstop)
{
  tb.Daq_Stop();
  DO_FLUSH
    }

CMD_PROC(dsize)
{
  unsigned int size = tb.Daq_GetSize();
  printf("size = %u\n", size);
}



CMD_PROC(dread)
{
  uint32_t words_remaining = 0;
  vector<uint16_t> data;
  tb.Daq_Read(data, 256, words_remaining);
  int size = data.size();
  printf("#samples: %i remaining: %i\n", size, int(words_remaining));

  for( int i=0; i<size; i++)
    {
      printf(" %04X", int(data[i]));
      if( i%10 == 9) printf("\n");
    }
  printf("\n");
}



// === ROC4sens ===================================================

CMD_PROC(cal)
{
  int x, y;
  PAR_INT(x, 0, 154)
    PAR_INT(y, 0, 159)
    tb.r4s_SetPixCal(x, y);
  DO_FLUSH
    }

CMD_PROC(hold)
{
  int t;
  PAR_INT(t, 0, 255)
    tb.r4s_SetHoldPos(t);
  DO_FLUSH
    }

CMD_PROC(adcdel)
{
  int t;
  PAR_INT(t, 0, 31)
    tb.r4s_AdcDelay(t);
  DO_FLUSH
    }

CMD_PROC(daqena)
{
  int flags;
  PAR_INT(flags, 0, 3)
    tb.r4s_Enable(flags);
  DO_FLUSH
    }

CMD_PROC(go)
{
  tb.r4s_Start();
  bool running = tb.r4s_Running();
  if( running) printf("running\n"); else printf("not running\n");
  DO_FLUSH
    }

CMD_PROC(goloop)
{
  while (!keypressed())
    {
      tb.r4s_Start();
      tb.mDelay(100);
    }
}

CMD_PROC(running)
{
  bool running = tb.r4s_Running();
  if( running) printf("true\n"); else printf("false\n");
}

CMD_PROC(vcal)
{
  int mV;
  PAR_INT(mV, 0, 2500)
    tb.r4s_SetVcal(mV);
  DO_FLUSH
    }

CMD_PROC(rgsh)
{
  int mV;
  PAR_INT(mV, 0, 2500)
    tb.r4s_SetRgsh(mV);
  DO_FLUSH
    }

CMD_PROC(rgpr)
{
  int mV;
  PAR_INT(mV, 0, 2500)
    tb.r4s_SetRgpr(mV);
  DO_FLUSH
    }

CMD_PROC(vref)
{
  int mV;
  PAR_INT(mV, 0, 2500)
    tb.r4s_SetVref(mV);
  DO_FLUSH
    }

CMD_PROC(vaux1)
{
  int mV;
  PAR_INT(mV, 0, 2500)
    tb.r4s_SetVaux1(mV);
  DO_FLUSH
    }

CMD_PROC(vaux2)
{
  int mV;
  PAR_INT(mV, 0, 2500)
    tb.r4s_SetVaux2(mV);
  DO_FLUSH
    }

CMD_PROC(vaux3)
{
  int mV;
  PAR_INT(mV, 0, 2500)
    tb.r4s_SetVaux3(mV);
  DO_FLUSH
    }
