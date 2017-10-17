
// Daniel Pitzl (DESY) Aug 2017
// .raw to .root

// r2r -n   100 raw000005.txt
// r2r -n  1340 raw000005.txt
// r2r -n 12500 raw000005.txt

#include <stdlib.h> // atoi
#include <fstream> // files
#include <iostream> // cout
#include <iomanip> // setw
#include <string> // strings
#include <sstream> // stringstream

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>

using namespace std;

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  cout << "main " << argv[0] << " called with " << argc << " arguments" << endl;

  if( argc == 1 ) {
    cout << "give file name" << endl;
    return 1;
  }

  // file name = last argument:

  string evFileName( argv[argc-1] );
  cout << "try to open  " << evFileName;

  ifstream evFile( argv[argc-1] );

  if( !evFile ) {
    cout << " : failed " << endl;
    return 2;
  }

  cout << " : succeed " << endl;

  // further arguments:

  int Nev = 2;

  for( int i = 1; i < argc; i++ ) {

    if( !strcmp( argv[i], "-n" ) )
      Nev = atoi( argv[++i] );

  } // argc

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  TFile* histoFile = new TFile( "r2r.root", "RECREATE" );

  // book histos:

  TH1I hadc( "adc", "ADC;ADC [ADC];pixels", 4097, -2048.5, 2048.5 );

  TH2D msump( "sump", "map running pedestal;col;row;ped [ADC]", 155, -0.5, 154.5, 160, -0.5, 159.5 );

  TH2I msum( "sum", "map entries;col;row;entries", 155, -0.5, 154.5, 160, -0.5, 159.5 );
  TH2I msuma( "suma", "map sum A;col;row;ADC", 155, -0.5, 154.5, 160, -0.5, 159.5 );
  TH2I msumaa( "sumaa", "map sum AA;col;row;ADC^{2}", 155, -0.5, 154.5, 160, -0.5, 159.5 );

  TH1I hped( "ped", "pedestals;PH [ADC];pixels", 4097, -2048.5, 2048.5 );

  TH1I hph( "ph", "PH;ADC-PED [ADC];pixels", 100, -100, 900 );

  TH1I hrms( "rms", "RMS;RMS [ADC];pixels", 100, 0, 20 );
  TH2D mrms( "mrms", "map RMS;col;row;RMS", 155, -0.5, 154.5, 160, -0.5, 159.5 );

  TH1I hsg( "sg", "PH/RMS;(ADC-PED)/RMS;pixels", 100, -10, 90 );

  TH1I hpht( "pht", "PH above cut;ADC-PED [ADC];pixels above cut", 100, -100, 900 );
  TH1I hsgt( "sgt", "PH/RMS above cut;(ADC-PED)/RMS;pixels above cut", 100, -10, 90 );
  TH2I mhit( "mhit", "map hits;col;row;hits", 155, -0.5, 154.5, 160, -0.5, 159.5 );
  TProfile2D mpht( "mpht", "map hit PH;col;row;hit <PH> [ADC]", 155, -0.5, 154.5, 160, -0.5, 159.5,-9999, 9999 );
  TH1I hphc1( "phc1", "PH with 1st cal;ADC-PED [ADC];1st pixel with cal", 100, -100, 900 );
  TH1I hphc2( "phc2", "PH with 2nd cal;ADC-PED [ADC];2nd pixel with cal", 100, -100, 900 );
  TH1I hfphc2( "fphc2", "fraction PH2;PH2/PH;2nd pixel with cal", 100, 0, 1 );

  TH1I hdph( "dph", "DPH;DPH [ADC];pixels", 500, -500, 500 );
  TH1I hdpht( "dpht", "DPH hits;DPH [ADC];pixel hits", 500, -0, 500 );
  TProfile2D mdph( "mdph", "map hit DPH;col;row;<DPH> [ADC]", 154, -0.5, 153.5, 160, -0.5, 159.5, 0, 9999 );
  TH2I mdht( "mdht", "map diff hits;col;row;diff hits", 154, -0.5, 153.5, 160, -0.5, 159.5 );

  TH1I hnht( "nht", "signal pixels per event;signal pixels;events", 101, -0.5, 100.5 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // event loop:

  int nev = 0;

  const int M = 100; // pedestal averaging range
  const double f = ( M - 1.0 ) / M; // update fraction

  // Read file by lines:

  string evnum;
  getline( evFile, evnum ); // read one line into string

  while( evFile.good() && ! evFile.eof() && nev < Nev ) {

    bool ldb = 0;
    if( nev < 10 ) ldb = 1;

    if( ldb ) cout << endl << "event " << nev << endl;
    if( nev < 100 && nev% 10 == 0 ) cout << "event " << nev << endl;
    if(              nev%100 == 0 ) cout << "event " << nev << endl;

    int npx = 0;
    int nht = 0;
    TH2I mev1( "ev", "map event;col;row;ADC", 155, -0.5, 154.5, 160, -0.5, 159.5 );
    double ph1 = 0;

    for( int row = 0; row < 163; ++row ) {

      string row_of_Pixels;
      getline( evFile, row_of_Pixels ); // read one line into string

      istringstream is( row_of_Pixels ); // tokenize string

      int col = 0;
      do {

	int adc;
	is >> adc;

	if( row < 160 && col < 155 ) { // suppress fake pixels

	  ++npx;

	  hadc.Fill( adc );

	  mev1.Fill( col, row, adc );

	  int ibin = msump.FindBin( col, row );

	  int nped = msum.GetBinContent( ibin );
	  double sump = msump.GetBinContent( ibin );

	  double ped = adc;
	  if( nped > 0 )
	    ped = sump/nped;

	  hped.Fill( ped );

	  double ph = ped-adc; // positive
	  if( nped > 0 )
	    hph.Fill( ph );

	  double avg = adc;
	  double rms = 10;

	  if( nped > 4 ) {
	    double suma = msuma.GetBinContent( ibin );
	    double sumaa = msumaa.GetBinContent( ibin );
	    avg = suma/nped;
	    rms = sqrt(sumaa/nped - avg*avg);
	    if( nped == M-1 ) { // plot once: final noise
	      hrms.Fill( rms );
	      mrms.Fill( col, row, rms );
	    }
	    hsg.Fill( ph/rms );
	  }

	  double cut = 0;
	  if(      nped > 0 ) cut = ped - 5*rms;
	  else if( nev > 1 ) cut = -2048;

	  if( adc > cut ) { // no (negative) signal

	    if( nped < M ) {
	      msum.Fill( col, row );
	      msuma.Fill( col, row, adc );
	      msumaa.Fill( col, row, adc*adc ); // for RMS
	      msump.Fill( col, row, adc ); // initial accumulated ped
	    }
	    else {
	      sump = adc + f*sump;
	      msump.SetBinContent( ibin, sump ); // update ped
	    }

	  }
	  else // hit
	    if( nped >= M ) {
	      //cout << nev << " " << col << " " << row << " " << ph << endl;
	      ++nht;
	      hpht.Fill( ph );
	      hsgt.Fill( ph/rms );
	      mhit.Fill( col, row );
	      mpht.Fill( col, row, ph );
	    }

	  // pix with cal:
	  if( (nev-100)/155 == row && (nev-100)%155 == col ) {
	    hphc1.Fill( ph );
	  }
	  else if( (nev-100)/155 == row && (nev-100)%155+1 == col && col%32 > 0 ) {
	    hphc2.Fill( ph );
	    hfphc2.Fill( (double) ph / (ph+ph1) );
	    if( ph < 40 )
	      cout << "small PH1: " << nev << " " << col << " " << row << " " << ph << endl;
	  }

	  ph1 = ph; // remember (for 2-col clus)

	} // valid pix

	++col;

      }
      while( ! is.eof() ); // one line = row

    } // rows

    if( ldb )
      cout
	<< "read " << npx << " valid pixels"
	<< ", found " << nht << " signals"
	<< endl;

    hnht.Fill( nht );

    // forward difference:

    if( nev > 99 ) {

      for( int col = 0; col < 154; ++col ) // one less

	for( int row = 0; row < 160; ++row ) {

	  int ibin = msump.FindBin( col, row );

	  int nped = msum.GetBinContent( ibin );
	  double sump = msump.GetBinContent( ibin );
	  double ped = sump/nped;

	  int adc = mev1.GetBinContent( ibin );
	  double ph = ped - adc; // positive

	  int jbin = msum.FindBin( col+1, row );
	  nped = msum.GetBinContent( jbin );
	  sump = msump.GetBinContent( jbin );
	  ped = sump/nped;

	  adc = mev1.GetBinContent( jbin );
	  double qh = ped - adc; // positive

	  int dph = qh-ph;
	  hdph.Fill( dph ); // RMS 3
	  if( fabs(dph) > 5*3 ) {
	    hdpht.Fill( dph );
	    mdht.Fill( col, row );
	    mdph.Fill( col, row, dph );
	  }

	} // row

    } // nev

    getline( evFile, evnum ); // read one line into string

    ++nev;

  } // while events

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cout << "done"
       << endl << "events " << nev
       << endl << histoFile->GetName()
       << endl;

  histoFile->Write();
  histoFile->Close();

  return 0;
}
