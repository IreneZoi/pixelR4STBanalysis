
// Daniel Pitzl (DESY) Sep 2017
// 2 x R4S: A and B

// zweiAB 338

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
    cout << "give run number" << endl;
    return 1;
  }

  // file name = last argument:

  string runnum( argv[argc-1] );
  int run = atoi( argv[argc-1] );

  string Afile = "107/roi000" + runnum + ".txt";
  cout << "try to open  " << Afile;
  ifstream Astream( Afile.c_str() );
  if( !Astream ) {
    cout << " : failed " << endl;
    return 2;
  }
  cout << " : succeed " << endl;

  string Bfile = "108/roi000" + runnum + ".txt";
  cout << "try to open  " << Bfile;
  ifstream Bstream( Bfile.c_str() );
  if( !Bstream ) {
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

  TFile* histoFile = new TFile( "zweiAB.root", "RECREATE" );

  // book histos:

  TH1I hphA( "phA", "A PH;ADC-PED [ADC];pixels", 200, -100, 900 );
  TH1I hnpxA( "npxA", "A PH pixels per event;PH pixels;events", 50, 0.5, 50.5 );
  TH2I * hpxmapA = new TH2I( "pxmapA", "A pixel map, PH > 30;col;row;PH pixels",
			     155, -0.5, 154.5, 160, -0.5, 159.5 );

  TH1I hphB( "phB", "B PH;ADC-PED [ADC];pixels", 200, -100, 900 );
  TH1I hnpxB( "npxB", "B PH pixels per event;PH pixels;events", 50, 0.5, 50.5 );
  TH2I * hpxmapB = new TH2I( "pxmapB", "B pixel map, PH > 30;col;row;PH pixels",
			     155, -0.5, 154.5, 160, -0.5, 159.5 );

  TH1I hnclA( "nclA", "A cluster per event;cluster;events", 20, 0.5,  20.5 );
  TH2I * hclmapA = new TH2I( "clmapA", "A cluster map;col;row;clusters",
			    155, -0.5, 154.5, 160, -0.5, 159.5 );
  TH1I hclszA( "clszA", "A cluster size;cluster size [pixels];clusters", 20, 0.5, 20.5 );
  TH1I hclphA( "clphA", "A cluster PH;cluster ph [ADC];clusters", 200, 0, 1000 );

  TH1I hnclB( "nclB", "B cluster per event;cluster;events", 20, 0.5,  20.5 );
  TH2I * hclmapB = new TH2I( "clmapB", "B cluster map;col;row;clusters",
			    155, -0.5, 154.5, 160, -0.5, 159.5 );
  TH1I hclszB( "clszB", "B cluster size;cluster size [pixels];clusters", 20, 0.5, 20.5 );
  TH1I hclphB( "clphB", "B cluster PH;cluster ph [ADC];clusters", 200, 0, 1000 );

  TH2I hxxBA( "xxBA", "B vs A;row A;row B;clusters", 160, -0.5, 159.5, 160, -0.5, 159.5 );
  TH2I hyyBA( "yyBA", "B vs A;col A;col B;clusters", 155, -0.5, 154.5, 155, -0.5, 154.5 );

  TH1I hdxBA( "dxBA", "Bx-Ax;x-x [px];cluster pairs", 200, -50, 50 );
  TH1I hdyBA( "dyBA", "By-Ay;y-y [px];cluster pairs", 200, -50, 50 );

  TH2I hdxvsev( "dxvsev", "Bx-Ax vs events;events;#Deltax [px];clusters",
		100, 0, 10000, 100, -50, 50 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Read file by lines:

  string START {"START"};
  string hd;

  while( hd != START ) {
    getline( Astream, hd ); // read one line into string
    cout << "  " << hd << endl;
  }

  hd.clear();
  while( hd != START ) {
    getline( Bstream, hd ); // read one line into string
    cout << "  " << hd << endl;
  }

  string F {"F"}; // filled flag

  int nev = 0;

  while( Astream.good() && ! Astream.eof() &&
	 Bstream.good() && ! Bstream.eof() &&
	 nev < Nev ) {

    string evseed;
    getline( Astream, evseed ); // read one line into string

    istringstream issA( evseed ); // tokenize string

    int iev;
    issA >> iev;
    int mpx = 0; // event pixels

    cout << "A ev " << iev;

    string filled;
    issA >> filled;
    if( filled == F ) {

      string roi;
      getline( Astream, roi );
      istringstream css( roi ); // tokenize string

      int npx = 0;

      while( ! css.eof() ) {

	int col;
	int row;
	double ph;
	css >> col;
	css >> row;
	css >> ph;

	hphA.Fill( ph );

	if( ph > 30 )
	  hpxmapA->Fill( col, row );

	if( ph > 30 && row > 4.5 && row < 154.5 ) { // for r108
	  pb[mpx].col = col;
	  pb[mpx].row = row;
	  pb[mpx].ph = ph;
	  ++mpx;
	}

	++npx;

      } // roi px

      cout << " roi " << npx << ", hits " << mpx << endl;

    } // filled

    else
      cout << "  empty" << endl;

    // clustering:

    hnpxA.Fill( mpx );

    fNHit = mpx; // for cluster search

    vector <cluster> vclA = getClus();	    

    hnclA.Fill( vclA.size() );
    if( vclA.size() )
      cout << "  clusters " << vclA.size() << endl;

    for( unsigned icl = 0; icl < vclA.size(); ++ icl ) {

      hclmapA->Fill( vclA[icl].col, vclA[icl].row );

      if( vclA[icl].row > 3.5 && vclA[icl].row < 155.5 ) { // Wed 13.9. runs 163-169

	if( vclA[icl].size < 6 )
	  hclphA.Fill( vclA[icl].sum );

	if( vclA[icl].sum > 55 )
	  hclszA.Fill( vclA[icl].size );

      }

    } // icl

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    getline( Bstream, evseed ); // read one line into string

    istringstream issB( evseed ); // tokenize string

    issB >> iev;
    mpx = 0; // event pixels

    cout << "B ev " << iev;

    issB >> filled;
    if( filled == F ) {

      string roi;
      getline( Bstream, roi );
      istringstream css( roi ); // tokenize string

      int npx = 0;

      while( ! css.eof() ) {

	int col;
	int row;
	double ph;
	css >> col;
	css >> row;
	css >> ph;

	hphB.Fill( ph );

	if( ph > 30 )
	  hpxmapB->Fill( col, row );

	if( ph > 30 && row > 4.5 && row < 154.5 ) { // for r108
	  pb[mpx].col = col;
	  pb[mpx].row = row;
	  pb[mpx].ph = ph;
	  ++mpx;
	}

	++npx;

      } // roi px

      cout << " roi " << npx << ", hits " << mpx << endl;

    } // filled

    else
      cout << "  empty" << endl;

    // clustering:

    hnpxB.Fill( mpx );

    fNHit = mpx; // for cluster search

    vector <cluster> vclB = getClus();	    

    hnclB.Fill( vclB.size() );
    if( vclB.size() )
      cout << "  clusters " << vclB.size() << endl;

    for( unsigned icl = 0; icl < vclB.size(); ++ icl ) {

      hclmapB->Fill( vclB[icl].col, vclB[icl].row );

      if( vclB[icl].row > 3.5 && vclB[icl].row < 155.5 ) { // Wed 13.9. runs 163-169

	if( vclB[icl].size < 6 )
	  hclphB.Fill( vclB[icl].sum );

	if( vclB[icl].sum > 55 )
	  hclszB.Fill( vclB[icl].size );

      }

    } // icl

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // B-A cluster correlations:

    for( vector<cluster>::iterator cB = vclB.begin(); cB != vclB.end(); ++cB ) {

      double xB = cB->row; // rot90 Dreimaster
      double yB = cB->col; // down

      for( vector<cluster>::iterator cA = vclA.begin(); cA != vclA.end(); ++cA ) {

	double xA = cA->row;
	double yA = cA->col;

	hxxBA.Fill( xA, xB );
	hyyBA.Fill( yA, yB );

	double dx = xB - xA;
	double dy = yB - yA;
	hdxBA.Fill( dx );
	hdyBA.Fill( dy );

	hdxvsev.Fill( nev, dx );

      } // clusters

    } // cl

    ++nev;

  } // while events

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cout << "done " << Afile
       << endl << "events " << nev
       << endl << histoFile->GetName()
       << endl;

  histoFile->Write();
  histoFile->Close();

  return 0;
}
