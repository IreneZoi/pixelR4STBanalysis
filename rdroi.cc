
// Daniel Pitzl (DESY) Sep 2017
// read region-of-interest data

// rdroi -n 2500 roi000102.txt

#include <stdlib.h> // atoi
#include <iostream> // cout
#include <iomanip> // setw
#include <string> // strings
#include <sstream> // stringstream
#include <fstream> // files
#include <vector>
#include <cmath>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>

using namespace std;

struct pixel {
  int col;
  int row;
  double ph;
};

struct cluster {
  vector <pixel> vpix;
  int size;
  int sum;
  double col,row;
};

// globals:

pixel pb[155*160]; // pixel block for clustering
int fNHit; // for clustering

// ----------------------------------------------------------------------
vector<cluster> getClus()
{
  // returns clusters with local coordinates
  // decodePixels should have been called before to fill pixel buffer pb 
  // simple clusterization
  // cluster search radius fCluCut ( allows fCluCut-1 empty pixels)

  const int fCluCut = 1; // clustering: 1 = no gap (15.7.2012)
  //const int fCluCut = 2;

  vector <cluster> v;
  if( fNHit == 0 ) return v;

  int * gone = new int[fNHit];

  for( int i = 0; i < fNHit; ++i )
    gone[i] = 0;

  int seed = 0;

  while( seed < fNHit ) {

    // start a new cluster

    cluster c;
    c.vpix.push_back( pb[seed] );
    gone[seed] = 1;

    // let it grow as much as possible:

    int growing;
    do{
      growing = 0;
      for( int i = 0; i < fNHit; ++i ) {
        if( !gone[i] ){ // unused pixel
          for( unsigned int p = 0; p < c.vpix.size(); ++p ) { // vpix in cluster so far
            int dr = c.vpix.at(p).row - pb[i].row;
            int dc = c.vpix.at(p).col - pb[i].col;
            if( (   dr>=-fCluCut) && (dr<=fCluCut) 
		&& (dc>=-fCluCut) && (dc<=fCluCut) ) {
              c.vpix.push_back(pb[i]);
	      gone[i] = 1;
              growing = 1;
              break; // important!
            }
          } // loop over vpix
        } // not gone
      } // loop over all pix
    }
    while( growing );

    // added all I could. determine position and append it to the list of clusters:

    c.sum = 0;
    c.size = 0;
    c.col = 0;
    c.row = 0;

    for( vector<pixel>::iterator p = c.vpix.begin();  p != c.vpix.end();  ++p ) {
      double PH = p->ph;
      c.sum += PH;
      c.col += (*p).col*PH;
      c.row += (*p).row*PH;
    }

    c.size = c.vpix.size();

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    c.col /= c.sum;
    c.row /= c.sum;

    v.push_back(c); // add cluster to vector

    // look for a new seed = used pixel:

    while( (++seed < fNHit) && gone[seed] );

  } // while over seeds

  // nothing left,  return clusters

  delete gone;
  return v;
}

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

  int Nev = 2000*1000;

  for( int i = 1; i < argc; i++ ) {

    if( !strcmp( argv[i], "-n" ) )
      Nev = atoi( argv[++i] );

  } // argc

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  TFile* histoFile = new TFile( "rdroi.root", "RECREATE" );

  // book histos:

  TH1I hph( "ph", "PH;ADC-PED [ADC];pixels", 1000, -100, 900 );
  TH1I hnpx( "npx", "PH pixels per event;PH pixels;events", 50, 0.5, 50.5 );
  TH2I * hpxmap = new TH2I( "pxmap", "pixel map, PH > 30;col;row;PH pixels",
			    155, -0.5, 154.5, 160, -0.5, 159.5 );
  TH1I hdph( "dph", "dPH;#DeltaPH [ADC];pixel", 1000, -100, 900 );

  TH1D hncl( "ncl", "cluster per event;cluster;events", 20, 0.5,  20.5 );
  TH2I * hclmap = new TH2I( "clmap", "cluster map;col;row;clusters",
			    155, -0.5, 154.5, 160, -0.5, 159.5 );
  TH1I hclsz( "clsz", "cluster size;cluster size [pixels];clusters", 20, 0.5, 20.5 );
  TH1I hclph( "clph", "cluster PH;cluster ph [ADC];clusters", 200, 0, 1000 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Read file by lines:

  string START {"START"};
  string hd;

  while( hd != START ) {
    getline( evFile, hd ); // read one line into string
    cout << "  " << hd << endl;
  }

  string F {"F"}; // filled flag

  string evseed;
  getline( evFile, evseed ); // read one line into string
  int nev = 0;

  while( evFile.good() && ! evFile.eof() && nev < Nev ) {

    istringstream iss( evseed ); // tokenize string
    int iev;
    iss >> iev;

    cout << "ev " << iev;

    int mpx = 0; // event pixels

    string filled;
    iss >> filled;
    if( filled == F ) {

      string roi;
      getline( evFile, roi );
      istringstream css( roi ); // tokenize string

      int npx = 0;
      vector <pixel> vpx;
      vpx.reserve(35);

      while( ! css.eof() ) {

	int col;
	int row;
	double ph;
	css >> col;
	css >> row;
	css >> ph;

	pixel px { col, row, ph };
	vpx.push_back(px);

	++npx;

	if( ph > 30 )
	  hpxmap->Fill( col, row );

      } // roi px

      // columns-wise common mode correction:

      for( unsigned ipx = 0; ipx < vpx.size(); ++ipx ) {

	int col4 = vpx[ipx].col;
	int row4 = vpx[ipx].row;
	double ph4 = vpx[ipx].ph;

	int row1 = row4;
	int row7 = row4;
	double ph1 = ph4;
	double ph7 = ph4;

	for( unsigned jpx = 0; jpx < vpx.size(); ++jpx ) {

	  if( jpx == ipx ) continue;
	  if( vpx[jpx].col != col4 ) continue; // want same column

	  int jrow = vpx[jpx].row;

	  if( jrow < row1 ) {
	    row1 = jrow;
	    ph1 = vpx[jpx].ph;
	  }

	  if( jrow > row7 ) {
	    row7 = jrow;
	    ph7 = vpx[jpx].ph;
	  }

	} // jpx

	if( row4 == row1 ) continue; // Randpixel
	if( row4 == row7 ) continue;

	double dph;
	if( row4 - row1 < row7 - row4 )
	  dph = ph4 - ph1;
	else
	  dph = ph4 - ph7;

	hph.Fill( ph4 );
	hdph.Fill( dph );

	if( dph > 24 ) {
	  pb[mpx].col = col4;
	  pb[mpx].row = row4;
	  pb[mpx].ph = dph;
	  ++mpx;
	}

      } // ipx

      cout << " roi " << npx << ", hits " << mpx << endl;

    } // filled

    else
      cout << "  empty" << endl;

    // clustering:

    hnpx.Fill( mpx );

    fNHit = mpx; // for cluster search

    vector <cluster> vcl = getClus();	    

    hncl.Fill( vcl.size() );
    if( vcl.size() )
      cout << "  clusters " << vcl.size() << endl;

    for( unsigned icl = 0; icl < vcl.size(); ++ icl ) {

      hclmap->Fill( vcl[icl].col, vcl[icl].row );

      if( vcl[icl].row > 3.5 && vcl[icl].row < 155.5 ) { // Wed 13.9. runs 163-169

	if( vcl[icl].size < 6 )
	  hclph.Fill( vcl[icl].sum );

	if( vcl[icl].sum > 55 )
	  hclsz.Fill( vcl[icl].size );

      }

      // debug:

      if( vcl[icl].row > 155.5 ) {
	cout << "cl " << icl << " at col " <<  vcl[icl].col << ", row " << vcl[icl].row << endl;
	for( int ipx = 0; ipx < vcl[icl].size; ++ipx )
	  cout << "  " << vcl[icl].vpix[ipx].col
	       << ", " << vcl[icl].vpix[ipx].row
	       << ": " << vcl[icl].vpix[ipx].ph
	    ;
	cout << endl;
      }

    } // icl

    ++nev;

    getline( evFile, evseed ); // read ahead

  } // while events

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cout << "done " << evFileName
       << endl << "events " << nev
       << endl << histoFile->GetName()
       << endl;

  histoFile->Write();
  histoFile->Close();

  return 0;
}
