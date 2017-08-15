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
//#include "datalink.h"
#include <iostream> // cout

#include <TFile.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TProfile.h>
#include <TProfile2D.h>

using namespace std;

#define IMG_WIDTH 157
#define IMG_HEIGHT 163
#define R4S_VALUE_OR 4096

//------------------------------------------------------------------------------
class R4sImg
{
  int * data;
public:
  R4sImg() : data(0) {}

  void Clear() { if( data ) delete[] data; data = 0; }

  bool CreateRaw( const vector<uint16_t> &rawdata ); // fills data

  int Get(int x, int y) { return data[y*IMG_WIDTH + x]; } // one pix

  void Print( unsigned int count ); // one row

  void Save( const string &filename ); // map.txt
};

//------------------------------------------------------------------------------
bool R4sImg::CreateRaw( const vector<uint16_t> &rawdata )
{
  Clear();
  if( rawdata.size() < IMG_WIDTH * IMG_HEIGHT ) return false;

  data = new int[IMG_WIDTH * IMG_HEIGHT];

  for( int pos = 0; pos < IMG_WIDTH * IMG_HEIGHT; ++pos ) {

    int value = rawdata[pos];

    if( value & 0x1000 ) // ADC overrange
      value = R4S_VALUE_OR;
    else if( value & 0x0800 ) // negative
      value -= 0x1000;

    data[pos] = value;
  }

  return true;
};

//------------------------------------------------------------------------------
void R4sImg::Print( unsigned int count )
{
  if( !data) {
    printf( "<<empty>>\n" );
    return;
  }

  if( count > IMG_WIDTH * IMG_HEIGHT)
    count = IMG_WIDTH * IMG_HEIGHT;

  unsigned i0 = 8*IMG_WIDTH;
  cout << "row " << i0/IMG_WIDTH << " from top" << endl;

  for( unsigned int i = i0; i < i0+count; ++i ) {
    int v = data[i];
    if( v == R4S_VALUE_OR )
      printf( "   or" );
    else
      printf( " %4i", v );
  }
  printf( "\n" );

  tb.Daq_Close();
}

//------------------------------------------------------------------------------
void R4sImg::Save( const string &filename )
{
  FILE *f = fopen( filename.c_str(), "wt" );

  TFile * histoFile = new TFile( "one.root", "RECREATE" );
  TProfile2D mapxy( "mapxy",
		    "ADC map;col;row;<ADC>",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );

  if( !f ) return;

  for( int y = IMG_HEIGHT-1; y >= 0; --y ) { // start top row

    for( int x = 0; x < IMG_WIDTH; ++x ) {

      int v = Get( x, y );

      if( v == R4S_VALUE_OR )
	fprintf( f, "   or" );
      else
	fprintf( f, " %4i", v );

      mapxy.Fill( x+0.5, y+0.5, v );

    }

    fputs( "\n", f );

  } // rows

  fclose(f);

  histoFile->Write();
  histoFile->Close();
  cout << histoFile->GetName() << endl;
}

//------------------------------------------------------------------------------
void ReadImage(R4sImg &map)
{
  tb.Daq_Open(50000);

  // prepare ADC:

  tb.SignalProbeADC( PROBEA_SDATA1, GAIN_1 );
  //	tb.r4s_AdcDelay(0);
  tb.r4s_Enable(3); // slow readout
  tb.uDelay(100);

  // take data:

  tb.Daq_Start();

  tb.r4s_Start();
  tb.uDelay(3000);

  tb.Daq_Stop();

  // stop ADC:

  tb.r4s_Enable(0);

  // read buffer:

  vector<uint16_t> data;
  unsigned int ret = tb.Daq_Read(data);

  printf("--- status = %u; n = %u\n", ret, (unsigned int)(data.size()));

  tb.Daq_Close();

  map.CreateRaw( data );
}

//------------------------------------------------------------------------------
CMD_PROC(getimg)
{
  R4sImg map;
  ReadImage(map);
  map.Print(200);
  map.Save("map.txt");
}

//------------------------------------------------------------------------------
CMD_PROC(getped)
{
  TFile * histoFile = new TFile( "ped.root", "RECREATE" );
  TProfile2D pedxy( "pedxy",
		    "pedestal map;col;row;<ADC>",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );
  TH2I hsuma( "suma",
	      "suma map;col;row;sum ADC",
	      IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT );
  TH2I hsumaa( "sumaa",
	       "sumaa map;col;row;sum ADC^2",
	       IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT );
  TProfile2D rmsxy( "rmsxy",
		    "rms map;col;row;RMS",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, 0, 2222 );
  TH1I hrms( "rms",
	     "rms;rms;pixels",
	     100, 0, 100 );

  R4sImg map;

  tb.r4s_SetSeqReadout(); // pedestal

  unsigned N = 99;
  for( unsigned i = 0; i < N; ++i ) { // events

    cout << "event " << i << endl;

    ReadImage( map );

    for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

      for( int x = 0; x < IMG_WIDTH; ++x ) {

	int a = map.Get( x, y );
	pedxy.Fill( x+0.5, y+0.5, a );
	hsuma.Fill( x+0.5, y+0.5, a );
	hsumaa.Fill( x+0.5, y+0.5, a*a );

      }

  } // events

  for( int y = IMG_HEIGHT-1; y >= 0; --y ) { // start top row

    for( int x = 0; x < IMG_WIDTH; ++x ) {
      double suma = hsuma.GetBinContent( x+1, y+1 ); // bins start at 1
      double sumaa = hsumaa.GetBinContent( x+1, y+1 ); // bins start at 1
      double rms = sqrt( sumaa/N - suma/N*suma/N );
      cout << " " << rms;
      rmsxy.Fill( x+0.5, y+0.5, rms );
      hrms.Fill( rms );
    }
    cout << endl;
  }
  cout << "mean rms " << hrms.GetMean() << endl;
  histoFile->Write();
  histoFile->Close();
  cout << histoFile->GetName() << endl;

}

//------------------------------------------------------------------------------
CMD_PROC(getcal) // ph-ped map
{
  TFile * histoFile = new TFile( "cal.root", "RECREATE" );
  TProfile2D pedxy( "pedxy",
		    "pedestal map;col;row;<ADC>",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );
  TH2I hsuma( "suma",
	      "suma map;col;row;sum ADC",
	      IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT );
  TH2I hsumaa( "sumaa",
	       "sumaa map;col;row;sum ADC^2",
	       IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT );
  TProfile2D rmsxy( "rmsxy",
		    "rms map;col;row;RMS",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, 0, 2222 );
  TH1I hrms( "rms",
	     "rms;rms;pixels",
	     100, 0, 100 );

  R4sImg map;

  tb.r4s_SetSeqReadout(); // pedestal

  unsigned Np = 99;

  for( unsigned i = 0; i < Np; ++i ) { // events

    cout << "event " << i << endl;

    ReadImage( map );

    for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

      for( int x = 0; x < IMG_WIDTH; ++x ) {

	int a = map.Get( x, y );
	pedxy.Fill( x+0.5, y+0.5, a );
	hsuma.Fill( x+0.5, y+0.5, a );
	hsumaa.Fill( x+0.5, y+0.5, a*a );

      }

  } // events

  for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

    for( int x = 0; x < IMG_WIDTH; ++x ) {
      double suma = hsuma.GetBinContent( x+1, y+1 ); // bins start at 1
      double sumaa = hsumaa.GetBinContent( x+1, y+1 ); // bins start at 1
      double rms = sqrt( sumaa/Np - suma/Np*suma/Np );
      rmsxy.Fill( x+0.5, y+0.5, rms );
      hrms.Fill( rms );
    }

  cout << "mean rms " << hrms.GetMean() << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  tb.r4s_SetSeqCalScan(); // Cal

  TProfile2D phxy( "phxy",
		   "pulse height;col;row;<ADC>",
		   IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );

  unsigned Nc = 10;
  for( unsigned i = 0; i < Nc; ++i ) { // events

    cout << "event " << i << endl;

    ReadImage( map );

    for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

      for( int x = 0; x < IMG_WIDTH; ++x )
	phxy.Fill( x+0.5, y+0.5, map.Get( x, y ) );


  } // events

  TProfile2D calxy( "calxy",
		    "cal;col;row;<ADC> - ped",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );
  TH1I hcal( "cal",
	     "cal;cal [ADC-ped];pixels",
	     100, -2000, 2000 );

  for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

    for( int x = 0; x < IMG_WIDTH; ++x ) {
      double a = phxy.GetBinContent( x+1, y+1 ); // bins start at 1
      double p = pedxy.GetBinContent( x+1, y+1 ); // bins start at 1
      calxy.Fill( x+0.5, y+0.5, a-p );
      hcal.Fill( a-p );
    }
  
  histoFile->Write();
  histoFile->Close();
  cout << histoFile->GetName() << endl;

}

//------------------------------------------------------------------------------
CMD_PROC(scancal) // scan Vcal
{
  TFile * histoFile = new TFile( "scancal.root", "RECREATE" );
  TProfile2D pedxy( "pedxy",
		    "pedestal map;col;row;<ADC>",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );
  TH2I hsuma( "suma",
	      "suma map;col;row;sum ADC",
	      IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT );
  TH2I hsumaa( "sumaa",
	       "sumaa map;col;row;sum ADC^2",
	       IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT );
  TProfile2D rmsxy( "rmsxy",
		    "rms map;col;row;RMS",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, 0, 2222 );
  TH1I hrms( "rms",
	     "rms;rms;pixels",
	     100, 0, 100 );

  R4sImg map;

  tb.r4s_SetSeqReadout(); // pedestal

  unsigned Np = 99;

  for( unsigned i = 0; i < Np; ++i ) { // events

    cout << "event " << i << endl;

    ReadImage( map );

    for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

      for( int x = 0; x < IMG_WIDTH; ++x ) {

	int a = map.Get( x, y );
	pedxy.Fill( x+0.5, y+0.5, a );
	hsuma.Fill( x+0.5, y+0.5, a );
	hsumaa.Fill( x+0.5, y+0.5, a*a );

      }

  } // events

  for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

    for( int x = 0; x < IMG_WIDTH; ++x ) {
      double suma = hsuma.GetBinContent( x+1, y+1 ); // bins start at 1
      double sumaa = hsumaa.GetBinContent( x+1, y+1 ); // bins start at 1
      double rms = sqrt( sumaa/Np - suma/Np*suma/Np );
      rmsxy.Fill( x+0.5, y+0.5, rms );
      hrms.Fill( rms );
    }

  cout << "mean rms " << hrms.GetMean() << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  tb.r4s_SetSeqCalScan(); // Cal

  TProfile phvscal( "phvscal", "PH vs Vcal pix 77 88;Vcal [mV];PH-ped [ADC]", 40, 25, 2025 );

  for( unsigned vcal = 50; vcal < 2002; vcal += 50 ) { // [mV]

    cout << "Vcal " << vcal << endl;

    tb.r4s_SetVcal(vcal);

    TProfile2D phxy( Form( "phxy_v%i", vcal ),
		     Form( "pulse height Vcal %i;col;row;<ADC>", vcal ),
		     IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );

    unsigned Nc = 10;
    for( unsigned i = 0; i < Nc; ++i ) { // events

      cout << "event " << i << endl;

      ReadImage( map );

      for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

      for( int x = 0; x < IMG_WIDTH; ++x )
	phxy.Fill( x+0.5, y+0.5, map.Get( x, y ) );


    } // events

    TProfile2D calxy( Form( "calxy_v%i", vcal ),
		      Form( "cal Vcal %i;col;row;<ADC> - ped", vcal ),
		      IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );
    TH1I hcal( Form( "cal_v%i", vcal ),
	       Form( "cal Vcal %i;cal [ADC-ped];pixels", vcal ),
	       100, -2000, 2000 );

    for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

      for( int x = 0; x < IMG_WIDTH; ++x ) {
	double a = phxy.GetBinContent( x+1, y+1 ); // bins start at 1
	double p = pedxy.GetBinContent( x+1, y+1 ); // bins start at 1
	calxy.Fill( x+0.5, y+0.5, p-a );
	hcal.Fill( p-a );
	if( x == 77 && y == 88 )
	  phvscal.Fill( vcal, p-a );
      }

    phxy.Write();
    calxy.Write();
    hcal.Write();

  } // Vcal

  histoFile->Write();
  histoFile->Close();
  cout << histoFile->GetName() << endl;

}

//------------------------------------------------------------------------------
CMD_PROC(seqreadout)
{
  tb.r4s_SetSeqReadout();
  DO_FLUSH;
}

//------------------------------------------------------------------------------
CMD_PROC(seqreadcol)
{
  tb.r4s_SetSeqReadCol();
  DO_FLUSH;
}

//------------------------------------------------------------------------------
CMD_PROC(seqcalscan)
{
  tb.r4s_SetSeqCalScan();
  DO_FLUSH;
}

/*
//------------------------------------------------------------------------------
CMD_PROC(gui)
{
  printf("Connect to GUI ...");
  CDataLink link(1024, 65536);
  //	ShellExecute(NULL, "open", "file.exe", NULL, NULL, SW_SHOWDEFAULT);

  system ("start roc4sens_view\\bin\\Release\\roc4sens_view.exe");
  //	system ("start roc4sens_view.exe");

  link.Run();
  printf("Connection to GUI closed.\n");
}
*/
