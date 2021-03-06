// psi46_tb.cpp

#include "pixel_dtb.h"
#include <stdio.h>

#ifndef _WIN32
#include <unistd.h>
#include <iostream>
#endif


bool CTestboard::RpcLink(bool verbose)
{
  bool error = false;
  for( unsigned short i = 2; i < rpc_cmdListSize; ++i ) {
    try {
      rpc_GetCallId(i);
    }
    catch (CRpcError &e) {
      e.SetFunction(0);
      if( verbose) {
	if( !error) printf("\nMissing DTB functions:\n");
	std::string fname(rpc_cmdName[i]);
	std::string fname_pretty;
	rpc_TranslateCallName(fname, fname_pretty);
	printf("%s\n", fname_pretty.c_str());
      }
      error = true;
    }
  }
  return !error;
}


bool CTestboard::EnumNext(string &name)
{
  char s[64];
  if( !usb.EnumNext(s)) return false;
  name = s;
  return true;
}


bool CTestboard::Enum(unsigned int pos, string &name)
{
  char s[64];
  if( !usb.Enum(s, pos)) return false;
  name = s;
  return true;
}


bool CTestboard::FindDTB(string &usbId)
{
  string name;
  vector<string> devList;
  unsigned int nDev;
  unsigned int nr;

  try {
    if( !EnumFirst(nDev)) throw int(1);
    for( nr=0; nr<nDev; nr++) {
      if( !EnumNext(name)) continue;
      if( name.size() < 4) continue;
      if( name.compare(0, 4, "DTB_") == 0) devList.push_back(name);
    }
  }
  catch (int e) {
    switch (e)
      {
      case 1: printf("Cannot access the USB driver\n"); return false;
      default: return false;
      }
  }

  if( devList.size() == 0) {
    printf("No DTB connected.\n");
    return false;
  }

  if( devList.size() == 1) {
    usbId = devList[0];
    return true;
  }

  // If more than 1 connected device list them
  printf("\nConnected DTBs:\n");
  for( nr=0; nr<devList.size(); nr++) {
    printf("%2u: %s", nr, devList[nr].c_str());
    if( Open(devList[nr], false)) {
      try {
	unsigned int bid = GetBoardId();
	printf("  BID=%2u\n", bid);
      }
      catch (...) {
	printf("  Not identifiable\n");
      }
      Close();
    }
    else printf(" - in use\n");
  }

  printf("Please choose DTB (0-%u): ", (nDev-1));
  char choice[8];
  fgets(choice, 8, stdin);
  sscanf (choice, "%d", &nr);
  if( nr >= devList.size() ) {
    nr = 0;
    printf("No DTB opened\n");
    return false;
  }

  usbId = devList[nr];
  return true;
}


bool CTestboard::Open(string &usbId, bool init)
{
  rpc_Clear();
  if( !usb.Open(&(usbId[0]))) return false;

  if( init) Init();
  return true;
}


void CTestboard::Close()
{
  //	if( usb.Connected()) Daq_Close();
  usb.Close();
  rpc_Clear();
}


void CTestboard::mDelay(uint16_t ms)
{
  Flush();
#ifdef _WIN32
  Sleep(ms);			// Windows
#else
  usleep(ms*1000);	// Linux
#endif
}


void CTestboard::r4s_SetPixCal( uint8_t x, uint8_t y )
{
  vector <uint32_t> regx(5); // initialized with zero
  vector <uint32_t> regy(5);

  unsigned int wpos, bpos;

  wpos = x / 32;
  bpos = 1 << (x % 32);
  if( wpos < 5 )
    regx[wpos] = bpos;

  wpos = y / 32;
  bpos = 1 << (y % 32);
  if( wpos < 5 )
    regy[wpos] = bpos;

  r4s_SetRegX(regx);
  r4s_SetRegY(regy);
}


void CTestboard::r4s_Set2PixCal( uint8_t x, uint8_t y )
{
  vector <uint32_t> regx(5); // initialized with zero
  vector <uint32_t> regy(5);

  unsigned int wpos, bpos;

  wpos = x / 32;
  bpos = 3 << (x % 32); // 3 = b_11 = two adjacent col
  if( wpos < 5 )
    regx[wpos] = bpos;

  wpos = y / 32;
  bpos = 1 << (y % 32);
  if( wpos < 5 )
    regy[wpos] = bpos;

  r4s_SetRegX(regx);
  r4s_SetRegY(regy);
}


void CTestboard::r4s_Set4PixCal( uint8_t x, uint8_t y )
{
  vector <uint32_t> regx(5); // initialized with zero
  vector <uint32_t> regy(5);

  unsigned int wpos, bpos;

  wpos = x / 32;
  bpos = 15 << (x % 32); // 15 = b_1111 = four adjacent col
  if( wpos < 5 )
    regx[wpos] = bpos;

  wpos = y / 32;
  bpos = 1 << (y % 32);
  if( wpos < 5 )
    regy[wpos] = bpos;

  r4s_SetRegX(regx);
  r4s_SetRegY(regy);
}

//------------------------------------------------------------------------------
// sequence building blocks:
// see ~/psi/r4s/src/FPGA/roc4sens_firmware/dtb/seq_main.v
/*
   0: stop
   1: seq_resetx
   2: seq_resety
   3: seq_loadx (from regx)
   4: seq_loady (from regy)
   5: seq_measure
   6: seq_readline
   7: seq_firstline
   8: seq_nextline
FW 0.7:
   9: seq_firstcol
   A: seq_nextcol
   B: seq_readcol
   C: seq_exttrg
   D: seq_calib
*/
//------------------------------------------------------------------------------
void CTestboard::r4s_SetSeqReadout( int ext )
{
  cout << "r4s_SetSeqReadout called with ext " << ext << endl;

  vector<uint32_t> prog(42); // read from right to left:
  if( ext ) // external trigger flag
    prog[ 0] = 0x21CD4321; // 21=clear SR x, y, 43=load cal D=cal C=trg 21=clear SR
  else
    prog[ 0] = 0xf2154321; // 21=clear SR x, y, 43=load cal, 5=meas, 21=clear SR
  prog[ 1] = 0x68686867; // 7=start row 0
  prog[ 2] = 0x68686868; // 8=next row, 6=read entire row (all cols)
  prog[ 3] = 0x68686868;
  prog[ 4] = 0x68686868;
  prog[ 5] = 0x68686868;
  prog[ 6] = 0x68686868;
  prog[ 7] = 0x68686868;
  prog[ 8] = 0x68686868;
  prog[ 9] = 0x68686868;
  prog[10] = 0x68686868; // 40
  prog[11] = 0x68686868;
  prog[12] = 0x68686868;
  prog[13] = 0x68686868;
  prog[14] = 0x68686868;
  prog[15] = 0x68686868;
  prog[16] = 0x68686868;
  prog[17] = 0x68686868;
  prog[18] = 0x68686868;
  prog[19] = 0x68686868;
  prog[20] = 0x68686868; // 80
  prog[21] = 0x68686868;
  prog[22] = 0x68686868;
  prog[23] = 0x68686868;
  prog[24] = 0x68686868;
  prog[25] = 0x68686868;
  prog[26] = 0x68686868;
  prog[27] = 0x68686868;
  prog[28] = 0x68686868;
  prog[29] = 0x68686868;
  prog[30] = 0x68686868; // 120
  prog[31] = 0x68686868;
  prog[32] = 0x68686868;
  prog[33] = 0x68686868;
  prog[34] = 0x68686868;
  prog[35] = 0x68686868;
  prog[36] = 0x68686868;
  prog[37] = 0x68686868;
  prog[38] = 0x68686868;
  prog[39] = 0x68686868;
  prog[40] = 0x68686868; // 160
  prog[41] = 0x08686868; // 163
  r4s_SetSequence(prog);
}

//------------------------------------------------------------------------------
void CTestboard::r4s_SetSeqReadCol( int ext ) // FW 0.7
{
  vector<uint32_t> prog(41);
  if( ext ) // external trigger flag
    prog[ 0] = 0x21CD4321; // clear_SR cal trg clear_SR
  else
    prog[ 0] = 0xf2154321; // clear SR, cal hold, clear SR,
  prog[ 1] = 0xbababab9; // 9 = start col 0
  prog[ 2] = 0xbabababa; // a = next column, b = read entire column (all rows)
  prog[ 3] = 0xbabababa;
  prog[ 4] = 0xbabababa;
  prog[ 5] = 0xbabababa;
  prog[ 6] = 0xbabababa;
  prog[ 7] = 0xbabababa;
  prog[ 8] = 0xbabababa;
  prog[ 9] = 0xbabababa;
  prog[10] = 0xbabababa; // 40
  prog[11] = 0xbabababa;
  prog[12] = 0xbabababa;
  prog[13] = 0xbabababa;
  prog[14] = 0xbabababa;
  prog[15] = 0xbabababa;
  prog[16] = 0xbabababa;
  prog[17] = 0xbabababa;
  prog[18] = 0xbabababa;
  prog[19] = 0xbabababa;
  prog[20] = 0xbabababa; // 80
  prog[21] = 0xbabababa;
  prog[22] = 0xbabababa;
  prog[23] = 0xbabababa;
  prog[24] = 0xbabababa;
  prog[25] = 0xbabababa;
  prog[26] = 0xbabababa;
  prog[27] = 0xbabababa;
  prog[28] = 0xbabababa;
  prog[29] = 0xbabababa;
  prog[30] = 0xbabababa; // 120
  prog[31] = 0xbabababa;
  prog[32] = 0xbabababa;
  prog[33] = 0xbabababa;
  prog[34] = 0xbabababa;
  prog[35] = 0xbabababa;
  prog[36] = 0xbabababa;
  prog[37] = 0xbabababa;
  prog[38] = 0xbabababa;
  prog[39] = 0xbabababa;
  prog[40] = 0x00000aba; // 157
  r4s_SetSequence(prog);
}

//------------------------------------------------------------------------------
void CTestboard::r4s_SetSeqCalScan( int ext )
{
  vector<uint32_t> prog(103);
  prog[ 0] = 0x86153721; // 1=resetx 2=resety 7=first 3=loadx 5=meas 1=resetx 6=read 8=next
  prog[ 1] = 0x15386153; prog[ 2] = 0x38615386; prog[ 3] = 0x61538615; prog[ 4] = 0x53861538;	prog[ 5] = 0x86153861;
  prog[ 6] = 0x15386153; prog[ 7] = 0x38615386; prog[ 8] = 0x61538615; prog[ 9] = 0x53861538;	prog[10] = 0x86153861;
  prog[11] = 0x15386153; prog[12] = 0x38615386; prog[13] = 0x61538615; prog[14] = 0x53861538;	prog[15] = 0x86153861;
  prog[16] = 0x15386153; prog[17] = 0x38615386; prog[18] = 0x61538615; prog[19] = 0x53861538;	prog[20] = 0x86153861;
  prog[21] = 0x15386153; prog[22] = 0x38615386; prog[23] = 0x61538615; prog[24] = 0x53861538;	prog[25] = 0x86153861;
  prog[26] = 0x15386153; prog[27] = 0x38615386; prog[28] = 0x61538615; prog[29] = 0x53861538;	prog[30] = 0x86153861;
  prog[31] = 0x15386153; prog[32] = 0x38615386; prog[33] = 0x61538615; prog[34] = 0x53861538;	prog[35] = 0x86153861;
  prog[36] = 0x15386153; prog[37] = 0x38615386; prog[38] = 0x61538615; prog[39] = 0x53861538;	prog[40] = 0x86153861;
  prog[41] = 0x15386153; prog[42] = 0x38615386; prog[43] = 0x61538615; prog[44] = 0x53861538;	prog[45] = 0x86153861;
  prog[46] = 0x15386153; prog[47] = 0x38615386; prog[48] = 0x61538615; prog[49] = 0x53861538;	prog[50] = 0x86153861;
  prog[51] = 0x15386153; prog[52] = 0x38615386; prog[53] = 0x61538615; prog[54] = 0x53861538;	prog[55] = 0x86153861;
  prog[56] = 0x15386153; prog[57] = 0x38615386; prog[58] = 0x61538615; prog[59] = 0x53861538;	prog[60] = 0x86153861;
  prog[61] = 0x15386153; prog[62] = 0x38615386; prog[63] = 0x61538615; prog[64] = 0x53861538;	prog[65] = 0x86153861;
  prog[66] = 0x15386153; prog[67] = 0x38615386; prog[68] = 0x61538615; prog[69] = 0x53861538;	prog[70] = 0x86153861;
  prog[71] = 0x15386153; prog[72] = 0x38615386; prog[73] = 0x61538615; prog[74] = 0x53861538;	prog[75] = 0x86153861;
  prog[76] = 0x15386153; prog[77] = 0x38615386; prog[78] = 0x61538615; prog[79] = 0x53861538;	prog[80] = 0x86153861;
  prog[81] = 0x15386153; prog[82] = 0x38615386; prog[83] = 0x61538615; prog[84] = 0x53861538;	prog[85] = 0x86153861;
  prog[86] = 0x15386153; prog[87] = 0x38615386; prog[88] = 0x61538615; prog[89] = 0x53861538;	prog[90] = 0x86153861;
  prog[91] = 0x15386153; prog[92] = 0x38615386; prog[93] = 0x61538615; prog[94] = 0x53861538;	prog[95] = 0x86153861;
  prog[96] = 0x15386153; prog[97] = 0x38615386; prog[98] = 0x61538615; prog[99] = 0x53861538;	prog[100]= 0x86153861;
  prog[101]= 0x15386153; prog[102]= 0x00000086;  // limit 1023
  r4s_SetSequence(prog);

  vector <uint32_t> calx( 5, 0xffffffff ); // activate all col cal: 5x32 = 160
  cout << "calx.size " << calx.size() << endl;

  r4s_SetRegX(calx);

}

//------------------------------------------------------------------------------
void CTestboard::r4s_SetSeqCalScanCol( int ext )
{
  vector<uint32_t> prog(99);
  prog[ 0] = 0xAB254921; // 1=resetx 2=resety 9=first col 4=loady 5=meas B=read col A=next col
  prog[ 1] = 0x254AB254; prog[ 2] = 0x4AB254AB; prog[ 3] = 0xB254AB25; prog[ 4] = 0x54AB254A;	prog[ 5] = 0xAB254AB2;
  prog[ 6] = 0x254AB254; prog[ 7] = 0x4AB254AB; prog[ 8] = 0xB254AB25; prog[ 9] = 0x54AB254A;	prog[10] = 0xAB254AB2;
  prog[11] = 0x254AB254; prog[12] = 0x4AB254AB; prog[13] = 0xB254AB25; prog[14] = 0x54AB254A;	prog[15] = 0xAB254AB2;
  prog[16] = 0x254AB254; prog[17] = 0x4AB254AB; prog[18] = 0xB254AB25; prog[19] = 0x54AB254A;	prog[20] = 0xAB254AB2;
  prog[21] = 0x254AB254; prog[22] = 0x4AB254AB; prog[23] = 0xB254AB25; prog[24] = 0x54AB254A;	prog[25] = 0xAB254AB2;
  prog[26] = 0x254AB254; prog[27] = 0x4AB254AB; prog[28] = 0xB254AB25; prog[29] = 0x54AB254A;	prog[30] = 0xAB254AB2;
  prog[31] = 0x254AB254; prog[32] = 0x4AB254AB; prog[33] = 0xB254AB25; prog[34] = 0x54AB254A;	prog[35] = 0xAB254AB2;
  prog[36] = 0x254AB254; prog[37] = 0x4AB254AB; prog[38] = 0xB254AB25; prog[39] = 0x54AB254A;	prog[40] = 0xAB254AB2;
  prog[41] = 0x254AB254; prog[42] = 0x4AB254AB; prog[43] = 0xB254AB25; prog[44] = 0x54AB254A;	prog[45] = 0xAB254AB2;
  prog[46] = 0x254AB254; prog[47] = 0x4AB254AB; prog[48] = 0xB254AB25; prog[49] = 0x54AB254A;	prog[50] = 0xAB254AB2;
  prog[51] = 0x254AB254; prog[52] = 0x4AB254AB; prog[53] = 0xB254AB25; prog[54] = 0x54AB254A;	prog[55] = 0xAB254AB2;
  prog[56] = 0x254AB254; prog[57] = 0x4AB254AB; prog[58] = 0xB254AB25; prog[59] = 0x54AB254A;	prog[60] = 0xAB254AB2;
  prog[61] = 0x254AB254; prog[62] = 0x4AB254AB; prog[63] = 0xB254AB25; prog[64] = 0x54AB254A;	prog[65] = 0xAB254AB2;
  prog[66] = 0x254AB254; prog[67] = 0x4AB254AB; prog[68] = 0xB254AB25; prog[69] = 0x54AB254A;	prog[70] = 0xAB254AB2;
  prog[71] = 0x254AB254; prog[72] = 0x4AB254AB; prog[73] = 0xB254AB25; prog[74] = 0x54AB254A;	prog[75] = 0xAB254AB2;
  prog[76] = 0x254AB254; prog[77] = 0x4AB254AB; prog[78] = 0xB254AB25; prog[79] = 0x54AB254A;	prog[80] = 0xAB254AB2;
  prog[81] = 0x254AB254; prog[82] = 0x4AB254AB; prog[83] = 0xB254AB25; prog[84] = 0x54AB254A;	prog[85] = 0xAB254AB2;
  prog[86] = 0x254AB254; prog[87] = 0x4AB254AB; prog[88] = 0xB254AB25; prog[89] = 0x54AB254A;	prog[90] = 0xAB254AB2;
  prog[91] = 0x254AB254; prog[92] = 0x4AB254AB; prog[93] = 0xB254AB25; prog[94] = 0x54AB254A;	prog[95] = 0xAB254AB2;
  prog[96] = 0x254AB254; prog[97] = 0x4AB254AB; prog[98] = 0x0000AB25;

  r4s_SetSequence(prog);

  vector <uint32_t> caly( 5, 0xffffffff ); // activate all col cal: 5x32 = 160
  cout << "caly.size " << caly.size() << endl;

  r4s_SetRegY(caly);

}

//------------------------------------------------------------------------------
// A sequence to read one single column. Column should be selected by r4s_SetPixCal(x, y)
void CTestboard::r4s_SetSeqSingleCol( int ext )
{
  vector<uint32_t> prog(1);
  prog[ 0] = 0x0b254321; // 1=resetx 2=resety 3=loadx 4=loady 5=meas b=read col

  r4s_SetSequence(prog);

}
