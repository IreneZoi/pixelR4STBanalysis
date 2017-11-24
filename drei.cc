
// Daniel Pitzl (DESY) Sep 2017
// 3 x R4S, 25x100 on rot90 PCB

// drei 338

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
  double q;
};

struct cluster {
  vector <pixel> vpix;
  int size; // [px]
  int sum; // [ADC]
  double q; // [ke]
  double col, row; // [px]
};

// globals:

pixel pb[155*160]; // pixel block for clustering
int fNHit; // for clustering

//------------------------------------------------------------------------------
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
      double ph = p->ph;
      c.sum += ph;
      double q = p->q;
      c.q += q;
      //c.col += (*p).col*ph;
      //c.row += (*p).row*ph;
      c.col += (*p).col*q;
      c.row += (*p).row*q;
    }

    c.size = c.vpix.size();

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    //c.col /= c.sum;
    //c.row /= c.sum;
    c.col /= c.q;
    c.row /= c.q;

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

  string Afile = "A/roi000" + runnum + ".txt";
  cout << "try to open  " << Afile;
  ifstream Astream( Afile.c_str() );
  if( !Astream ) {
    cout << " : failed " << endl;
    return 2;
  }
  cout << " : succeed " << endl;

  string Bfile = "B/roi000" + runnum + ".txt";
  cout << "try to open  " << Bfile;
  ifstream Bstream( Bfile.c_str() );
  if( !Bstream ) {
    cout << " : failed " << endl;
    return 2;
  }
  cout << " : succeed " << endl;

  string Cfile = "C/roi000" + runnum + ".txt";
  cout << "try to open  " << Cfile;
  ifstream Cstream( Cfile.c_str() );
  if( !Cstream ) {
    cout << " : failed " << endl;
    return 2;
  }
  cout << " : succeed " << endl;

  // further arguments:

  int Nev = 20*1000*1000;

  for( int i = 1; i < argc; i++ ) {

    if( !strcmp( argv[i], "-n" ) )
      Nev = atoi( argv[++i] );

  } // argc

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // gains:

  double Ap0[155][160]; // Fermi
  double Ap1[155][160];
  double Ap2[155][160];
  double Ap3[155][160];

  //string gainA{ "A/r113-scancal-tb21-0921.dat"}; // lots of negative q
  string gainA{ "A/r113-scancal-tb21-0923.dat"};
  if( run >= 423 ) gainA = "A/r113-scancal-tb21-0923.dat";
  if( run >= 430 ) gainA = "A/r112-scancal-tb21-0925.dat";
  if( run >= 439 ) gainA = "A/r113-scancal-tb21-0923.dat";

  double Ake = 0.039; // Landau peak at 11 ke

  if( run >= 423 ) Ake = 0.0396; // Landau peak at 11 ke
  if( run >= 432 ) Ake = 0.039; // r112

  ifstream gainFileA( gainA );

  if( ! gainFileA ) {
    cout << "gain file for A not found" << endl;
    return 1;
  }

  while( ! gainFileA.eof() ) {

    int icol;
    int irow;
    gainFileA >> icol;
    gainFileA >> irow;
    gainFileA >> Ap0[icol][irow];
    gainFileA >> Ap1[icol][irow];
    gainFileA >> Ap2[icol][irow];
    gainFileA >> Ap3[icol][irow];

  } // while

  // B:

  double Bp0[155][160]; // Fermi
  double Bp1[155][160];
  double Bp2[155][160];
  double Bp3[155][160];

  string gainB{ "B/r108-scancal-tb21-0921.dat" };

  //if( run >= 423 ) gainB = "B/r108-scancal-tb21-0923-hold24.dat";
  if( run >= 423 ) gainB = "B/r108-scancal-tb21-0923-hold25.dat";

  //if( run >= 432 ) gainB = "C/r110-scancal-tb21-0921.dat"; // old, test
  if( run >= 432 ) gainB = "B/r110-scancal-tb21-0925-hold25.dat"; // 12.6 ke
  //if( run >= 432 ) gainB = "B/r110-scancal-tb21-0925-hold26.dat"; // 12.0 ke

  if( run >= 444 ) gainB = "B/r110-scancal-tb21-0928-hold24.dat";

  double Bke = 0.0276; // Landau peak at 11 ke

  if( run >= 423 ) Bke = 0.026;

  //if( run >= 432 ) Bke = 0.039; // old, test
  if( run >= 432 ) Bke = 0.0326; // r110

  if( run >= 443 ) Bke = 0.036; // cmspixel-daq

  ifstream gainFileB( gainB );

  if( ! gainFileB ) {
    cout << "gain file for B not found" << endl;
    return 1;
  }

  while( ! gainFileB.eof() ) {

    int icol;
    int irow;
    gainFileB >> icol;
    gainFileB >> irow;
    gainFileB >> Bp0[icol][irow];
    gainFileB >> Bp1[icol][irow];
    gainFileB >> Bp2[icol][irow];
    gainFileB >> Bp3[icol][irow];

  } // while

  // C:

  double Cp0[155][160]; // Fermi
  double Cp1[155][160];
  double Cp2[155][160];
  double Cp3[155][160];

  string gainC{ "C/r110-scancal-tb21-0921.dat" };

  if( run >= 423 ) gainC = "C/r110-scancal-tb21-0923-hold25.dat";

  //if( run >= 432 ) gainC = "C/r114-scancal-tb21-0925-hold24.dat"; // 14.9
  if( run >= 432 ) gainC = "C/r114-scancal-tb21-0925-hold25.dat"; // 15.2
  //if( run >= 432 ) gainC = "C/r114-scancal-tb21-0925-hold26.dat"; // 14.9

  double Cke = 0.0366; // Landau peak at 11 ke
  if( run >= 432 ) Cke = 0.028; // Landau peak at 11 ke

    ifstream gainFileC( gainC );

  if( ! gainFileC ) {
    cout << "gain file for C not found" << endl;
    return 1;
  }

  while( ! gainFileC.eof() ) {

    int icol;
    int irow;
    gainFileC >> icol;
    gainFileC >> irow;
    gainFileC >> Cp0[icol][irow];
    gainFileC >> Cp1[icol][irow];
    gainFileC >> Cp2[icol][irow];
    gainFileC >> Cp3[icol][irow];

  } // while

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // run 392:

  double alignxA =  0.092; // [mm] same sign as dx
  double alignyA =  0.100; // [mm] same sign as dy
  double alignfA =  0.0085; // [rad] same sign dxvsy
  //double txA = 0;
  //double tyA = 0;

  double alignxC = -0.466; // [mm] same sign as dx
  double alignyC = -0.050; // [mm] same sign as dy
  double alignfC =  0.0105; // [rad] same sign dxvsy
  //double txC = 0;
  //double tyC = 0;

  if( run >= 393 ) { // tilt 10 deg

    alignxA =  0.39; // [mm] same sign as dx
    alignyA = -0.20; // [mm] same sign as dy
    alignfA =  0.008; // [rad] same sign dxvsy

    alignxC = -0.73; // [mm] same sign as dx
    alignyC =  0.25; // [mm] same sign as dy
    alignfC =  0.01; // [rad] same sign dxvsy
  }

  if( run >= 402 ) { // tilt 11 deg

    alignxA =  0.405; // [mm] same sign as dx
    alignyA = -0.30; // [mm] same sign as dy
    alignfA =  0.008; // [rad] same sign dxvsy

    alignxC = -0.76; // [mm] same sign as dx
    alignyC =  0.30; // [mm] same sign as dy
    alignfC =  0.01; // [rad] same sign dxvsy
  }

  if( run >= 403 ) { // tilt 12 deg

    alignxA =  0.405; // [mm] same sign as dx
    alignyA = -0.30; // [mm] same sign as dy
    alignfA =  0.008; // [rad] same sign dxvsy

    alignxC = -0.76; // [mm] same sign as dx
    alignyC =  0.30; // [mm] same sign as dy
    alignfC =  0.01; // [rad] same sign dxvsy

  }

  if( run >= 404 ) { // tilt 13 deg

    alignxA =  0.423; // [mm] same sign as dx
    alignyA = -0.30; // [mm] same sign as dy
    alignfA =  0.008; // [rad] same sign dxvsy

    alignxC = -0.766; // [mm] same sign as dx
    alignyC =  0.30; // [mm] same sign as dy
    alignfC =  0.01; // [rad] same sign dxvsy

  }

  if( run >= 405 ) { // tilt 10 deg Thu night

    alignxA =  0.408; // [mm] same sign as dx
    alignyA = -0.300; // [mm] same sign as dy
    alignfA =  0.0083; // [rad] same sign dxvsy

    alignxC = -0.758; // [mm] same sign as dx
    alignyC =  0.350; // [mm] same sign as dy
    alignfC =  0.0104; // [rad] same sign dxvsy

  }

  if( run >= 406 ) { // turn 8

    alignxA =  0.414; // [mm] same sign as dx
    alignyA = -0.25; // [mm] same sign as dy
    alignfA =  0.0075; // [rad] same sign dxvsy

    alignxC = -0.766; // [mm] same sign as dx
    alignyC =  0.30; // [mm] same sign as dy
    alignfC =  0.010; // [rad] same sign dxvsy

  }

  if( run >= 407 ) { // turn 7

    alignxA =  0.411; // [mm] same sign as dx
    alignyA = -0.25; // [mm] same sign as dy
    alignfA =  0.008; // [rad] same sign dxvsy

    alignxC = -0.768; // [mm] same sign as dx
    alignyC =  0.30; // [mm] same sign as dy
    alignfC =  0.010; // [rad] same sign dxvsy

  }

  if( run >= 408 ) {

    alignxA =  0.418; // [mm] same sign as dx
    alignyA = -0.20; // [mm] same sign as dy
    alignfA =  0.0075; // [rad] same sign dxvsy

    alignxC = -0.776; // [mm] same sign as dx
    alignyC =  0.25; // [mm] same sign as dy
    alignfC =  0.0079; // [rad] same sign dxvsy

  }

  if( run >= 409 ) {

    alignxA =  0.413; // [mm] same sign as dx
    alignyA = -0.20; // [mm] same sign as dy
    alignfA =  0.0075; // [rad] same sign dxvsy

    alignxC = -0.783; // [mm] same sign as dx
    alignyC =  0.25; // [mm] same sign as dy
    alignfC =  0.0079; // [rad] same sign dxvsy

  }

  if( run >= 423 ) { // tilt 11 deg Fri night

    alignxA =  0.325; // [mm] same sign as dx
    alignyA = -0.20; // [mm] same sign as dy
    alignfA =  0.0086; // [rad] same sign dxvsy

    alignxC = -0.678; // [mm] same sign as dx
    alignyC =  0.25; // [mm] same sign as dy
    alignfC =  0.0106; // [rad] same sign dxvsy

  }

  if( run >= 432 ) { // Sun 24.9.2017

    alignxA = -0.073; // [mm] same sign as dx
    alignyA =  0.05;  // [mm] same sign as dy
    alignfA = -0.0012; // [rad] same sign dxvsy

    alignxC =  0.005; // [mm] same sign as dx
    alignyC = -0.15;  // [mm] same sign as dy
    alignfC = -0.0013; // [rad] same sign dxvsy

  }

  if( run >= 443 ) { // Wed 27.9.2017

    alignxA = -0.458; // [mm] same sign as dx
    alignyA =  0.00;  // [mm] same sign as dy
    alignfA = -0.008; // [rad] same sign dxvsy

    alignxC =  0.112; // [mm] same sign as dx
    alignyC = -0.20;  // [mm] same sign as dy
    alignfC = -0.007; // [rad] same sign dxvsy

  }

  if( run >= 444 ) { // Thu 28.9.2017

    alignxA = -0.440; // [mm] same sign as dx
    alignyA =  0.00;  // [mm] same sign as dy
    alignfA = -0.006; // [rad] same sign dxvsy

    alignxC =  0.094; // [mm] same sign as dx
    alignyC = -0.25;  // [mm] same sign as dy
    alignfC = -0.004; // [rad] same sign dxvsy

  }

  if( run >= 446 ) { // Thu 28.9.2017 turn 10

    alignxA = -0.465; // [mm] same sign as dx
    alignyA =  0.00;  // [mm] same sign as dy
    alignfA = -0.006; // [rad] same sign dxvsy

    alignxC =  0.077; // [mm] same sign as dx
    alignyC = -0.25;  // [mm] same sign as dy
    alignfC = -0.004; // [rad] same sign dxvsy

  }

  if( run >= 447 ) { // Thu 28.9.2017 turn 9

    alignxA = -0.461; // [mm] same sign as dx
    alignyA =  0.00;  // [mm] same sign as dy
    alignfA = -0.006; // [rad] same sign dxvsy

    alignxC =  0.081; // [mm] same sign as dx
    alignyC = -0.25;  // [mm] same sign as dy
    alignfC = -0.004; // [rad] same sign dxvsy

  }

  double cfA = cos(alignfA);
  double sfA = sin(alignfA);
  double cfC = cos(alignfC);
  double sfC = sin(alignfC);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  TFile* histoFile = new TFile( Form( "drei-r%i.root", run ), "RECREATE" );

  // book histos:

  TH1I hphA( "phA", "A PH;ADC-PED [ADC];pixels", 500, -100, 900 );
  TH1I hdphA( "dphA", "A #DeltaPH;#DeltaPH [ADC];pixels", 500, -100, 900 );
  TH1I hnpxA( "npxA", "A PH pixels per event;PH pixels;events", 50, 0.5, 50.5 );
  TH2I * hpxmapA = new TH2I( "pxmapA", "A pixel map, PH > cut;col;row;PH pixels",
			     155, -0.5, 154.5, 160, -0.5, 159.5 );

  TH1I hphB( "phB", "B PH;ADC-PED [ADC];pixels", 500, -100, 900 );
  TH1I hdphB( "dphB", "B #DeltaPH;#DeltaPH [ADC];pixels", 500, -100, 900 );
  TH1I hnpxB( "npxB", "B PH pixels per event;PH pixels;events", 50, 0.5, 50.5 );
  TH2I * hpxmapB = new TH2I( "pxmapB", "B pixel map, PH > cut;col;row;PH pixels",
			     155, -0.5, 154.5, 160, -0.5, 159.5 );

  TH1I hphC( "phC", "C PH;ADC-PED [ADC];pixels", 500, -100, 900 );
  TH1I hdphC( "dphC", "C #DeltaPH;#DeltaPH [ADC];pixels", 500, -100, 900 );
  TH1I hnpxC( "npxC", "C PH pixels per event;PH pixels;events", 50, 0.5, 50.5 );
  TH2I * hpxmapC = new TH2I( "pxmapC", "C pixel map, PH > cut;col;row;PH pixels",
			     155, -0.5, 154.5, 160, -0.5, 159.5 );

  TH1D hnclA( "nclA", "A cluster per event;A clusters;events", 20, 0.5,  20.5 );
  TH2I * hclmapA = new TH2I( "clmapA", "A cluster map;col;row;A clusters",
			    155, -0.5, 154.5, 160, -0.5, 159.5 );
  TH1I hclszA( "clszA", "A cluster size;cluster size [pixels];A clusters", 20, 0.5, 20.5 );
  TH1I hclphA( "clphA", "A cluster PH;cluster ph [ADC];A clusters", 200, 0, 1000 );
  TH1I hclqA( "clqA", "A cluster charge;cluster charge [ke];A clusters", 100, 0, 50 );

  TH1D hnclB( "nclB", "B cluster per event;A clusters;events", 20, 0.5,  20.5 );
  TH2I * hclmapB = new TH2I( "clmapB", "B cluster map;col;row;B clusters",
			    155, -0.5, 154.5, 160, -0.5, 159.5 );
  TH1I hclszB( "clszB", "B cluster size;cluster size [pixels];B clusters", 20, 0.5, 20.5 );
  TH1I hclphB( "clphB", "B cluster PH;cluster ph [ADC];B clusters", 200, 0, 1000 );
  TH1I hclqB( "clqB", "B cluster charge;cluster charge [ke];B clusters", 100, 0, 50 );

  TH1D hnclC( "nclC", "C cluster per event;C cluster;events", 20, 0.5,  20.5 );
  TH2I * hclmapC = new TH2I( "clmapC", "C cluster map;col;row;C clusters",
			    155, -0.5, 154.5, 160, -0.5, 159.5 );
  TH1I hclszC( "clszC", "C cluster size;cluster size [pixels];C clusters", 20, 0.5, 20.5 );
  TH1I hclphC( "clphC", "C cluster PH;cluster ph [ADC];C clusters", 200, 0, 1000 );
  TH1I hclqC( "clqC", "C cluster charge;cluster charge [ke];C clusters", 100, 0, 50 );

  // correlations:

  TH2I hxxAB( "xxAB", "B vs A;row A;row B;clusters", 320, -4, 4, 320, -4, 4 );
  TH2I hyyAB( "yyAB", "B vs A;col A;col B;clusters",  80, -4, 4,  80, -4, 4 );

  //double f = 0.5;
  double f = 0.1; // aligned

  TH1I hdxAB( "dxAB", "Bx-Ax;x-x [mm];cluster pairs", 200, -1, 1 );
  TH1I hdyAB( "dyAB", "By-Ay;y-y [mm];cluster pairs", 200, -1, 1 );
  TProfile dxvsxAB( "dxvsxAB", "dx vs x A-B;y [mm];<dx> [mm]", 320, -4, 4, -f, f );
  TProfile dxvsyAB( "dxvsyAB", "dx vs y A-B;y [mm];<dx> [mm]",  80, -4, 4, -f, f );

  TH2I hdxvsev( "dxvsev", "Bx-Ax vs events;events;#Deltax [px];clusters",
		100, 0, 10000, 100, -f, f );

  TH2I hxxCB( "xxCB", "C vs B;row B;row C;clusters", 320, -4, 4, 320, -4, 4 );
  TH2I hyyCB( "yyCB", "C vs B;col B;col C;clusters",  80, -4, 4,  80, -4, 4 );

  TH1I hdxCB( "dxCB", "Cx-Bx;x-x [mm];cluster pairs", 200, -1, 1 );
  TH1I hdyCB( "dyCB", "Cy-By;y-y [mm];cluster pairs", 200, -1, 1 );
  TProfile dxvsxCB( "dxvsxCB", "dx vs x C-B;y [mm];<dx> [mm]", 320, -4, 4, -f, f );
  TProfile dxvsyCB( "dxvsyCB", "dx vs y C-B;y [mm];<dx> [mm]",  80, -4, 4, -f, f );

  // triplets:

  TH2I hxxCA( "xxCA", "C vs A;row A;row C;clusters", 320, -4, 4, 320, -4, 4 );
  TH2I hyyCA( "yyCA", "C vs A;col A;col C;clusters",  80, -4, 4,  80, -4, 4 );

  TH1I hdxCA( "dxCA", "Cx-Ax;x-x [mm];cluster pairs", 200, -1, 1 );
  TH1I hdyCA( "dyCA", "Cy-Ay;y-y [mm];cluster pairs", 200, -1, 1 );

  TProfile dxvsxCA( "dxvsxCA", "dx vs x C-A;y [mm];<dx> [mm]", 320, -4, 4, -f, f );
  TProfile dxvsyCA( "dxvsyCA", "dx vs y C-A;y [mm];<dx> [mm]",  80, -4, 4, -f, f );
  TH1I hdyCAc( "dyCAc", "Cy-Ay, cut dx;y-y [mm];cluster pairs", 200, -1, 1 );
  TProfile dyvsyCA( "dyvsyCA", "dy vs y C-A;y [mm];<dy> [mm]",  80, -4, 4, -f, f );

  TH1I hdx3( "dx3", "triplet dx;dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I hdy3( "dy3", "triplet dy;dy [mm];triplets", 200, -1, 1 );
  TProfile dx3vsx( "dx3vsx", "dx vs x;x [mm];<dx3> [mm]", 320, -4, 4, -0.5, 0.5 );
  TProfile dx3vsy( "dx3vsy", "dx vs y;y [mm];<dx3> [mm]",  80, -4, 4, -0.5, 0.5 );
  TProfile dx3vsxm( "dx3vsxm", "dx vs x mod 100 um;x mod 100 [#mum];<dx3> [mm]",
		    100, 0, 100, -0.5, 0.5 );
  TH1I hdx3m( "dx3m", "triplet dx, x < 0;dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I hdx3p( "dx3p", "triplet dx, x > 0;dx [mm];triplets", 500, -0.5, 0.5 );

  TProfile madx3vsx( "madx3vsx", "MAD(dx3) vs x;x [mm];MAD dx3 [mm]", 320, -4, 4, 0, 0.1 );
  TProfile madx3vsdx( "madx3vsdx", "MAD(dx3) vs dx C-A;C-A dx [#mum];MAD dx3 [mm]",
		      100, -100, 100, 0, 0.1 );
  TProfile madx3vsy( "madx3vsy", "MAD(dx3) vs y;y [mm];MAD dx3 [mm]",  80, -4, 4, 0, 0.1 );
  TProfile madx3vsxm( "madx3vsxm", "MAD(dx3) vs xmod;x mod 100 [#mum];MAD dx3 [mm]",
		      100, 0, 100, 0, 0.1 );
  TProfile etavsxmB3( "etavsxmB3", "eta vs xmod;x mod 100 [#mum];B <eta>",
		      100, 0, 100, -1.1, 1.1 );
  TProfile madx3vseta( "madx3vseta", "MAD(dx3) vs eta;eta;MAD dx3 [mm]",
		      100, -1, 1, 0, 0.1 );

  TH1I hclszA3( "clszA3", "A cluster size on tracks;cluster size [pixels];Aclusters on tracks",
		40, 0.5, 40.5 );
  TH1I hclphA3( "clphA3", "A cluster PH on tracks;cluster ph [ADC];A clusters on tracks",
		200, 0, 1000 );
  TH1I hclqA3( "clqA3", "A cluster charge on tracks;cluster charge [ke];A clusters on tracks",
	       160, 0, 80 );

  TH1I hclszB3( "clszB3", "B cluster size on tracks;cluster size [pixels];B clusters on tracks",
		40, 0.5, 40.5 );
  TH1I hclphB3( "clphB3", "B cluster PH on tracks;cluster ph [ADC];B clusters on tracks",
		200, 0, 1000 );
  TH1I hclqB3( "clqB3", "B cluster charge on tracks;cluster charge [ke];B clusters on tracks",
	       160, 0, 80 );
  TH1I hclqB3n( "clqB3n",
		"B cluster charge on tracks, npx < 4;cluster charge [ke];B clusters on tracks, npx < 4",
		160, 0, 80 );

  TH1I hclszC3( "clszC3", "C cluster size on tracks;cluster size [pixels];C clusters on tracks",
		40, 0.5, 40.5 );
  TH1I hclphC3( "clphC3", "C cluster PH on tracks;cluster ph [ADC];C clusters on tracks",
		200, 0, 1000 );
  TH1I hclqC3( "clqC3", "C cluster charge on tracks;cluster charge [ke];C clusters on tracks",
	       160, 0, 80 );

  TProfile nrowvsxmB3( "nrowvsxmB3",
		       "B rows vs xmod;x mod 100 [#mum];<B cluster size [rows]>",
		       100, 0, 100, 0.5, 10.5 );
  TProfile clqvsxmB3( "clqvsxmB3",
		      "B cluster charge vs xmod;x mod 100 [#mum];<B cluster charge [ke]>",
		      100, 0, 100, 0, 50 );

  TProfile madx3vsq( "madx3vsq", "MAD(dx) vs Q;B cluster charge [ke];MAD dx [mm]",
		      50, 0, 50, 0, 0.1 );
  TProfile madx3vsn( "madx3vsn", "MAD(dx) vs cluster size;B cluster size [pixels];MAD dx [mm]",
		      20, 0.5, 20.5, 0, 0.1 );
  TH1I hdx3c( "dx3c", "triplet dx, cut dy;dx [mm];triplets",
	      500, -0.25, 0.25 );

  TH1I hdx3c1( "dx3c1", "triplet dx, cut dy, npx 1;dx [mm];triplets, B npx 1",
	       500, -0.25, 0.25 );
  TH1I hdx3c2( "dx3c2", "triplet dx, cut dy, npx 2;dx [mm];triplets, B npx 2",
	       500, -0.25, 0.25 );
  TH1I hdx3c3( "dx3c3", "triplet dx, cut dy, npx 3;dx [mm];triplets, B npx 3",
	       500, -0.25, 0.25 );
  TH1I hdx3c4( "dx3c4", "triplet dx, cut dy, npx 4;dx [mm];triplets, B npx 4",
	       500, -0.25, 0.25 );
  TH1I hdx3c5( "dx3c5", "triplet dx, cut dy, npx 5;dx [mm];triplets, B npx 5",
	       500, -0.25, 0.25 );
  TH1I hdx3c6( "dx3c6", "triplet dx, cut dy, npx 6;dx [mm];triplets, B npx 6",
	       500, -0.25, 0.25 );
  TH1I hdx3c7( "dx3c7", "triplet dx, cut dy, npx > 6;dx [mm];triplets, B npx > 6",
	       500, -0.25, 0.25 );

  TH1I hdx3ct( "dx3ct", "triplet dx, cut dy, tx;dx [mm];triplets",
	       500, -0.25, 0.25 );
  TH1I hdx3cq( "dx3cq", "triplet dx, Landau peak;dx [mm];Landau peak triplets",
	       500, -0.25, 0.25 );
  TH1I hdx3cq3( "dx3cq3", "triplet dx, 3 Landau peak;dx [mm];3 Landau peak triplets",
		500, -0.25, 0.25 );
  TH1I hdx3cq3t( "dx3cq3t",
		 "triplet dx, 3 Landau peak, forward;dx [mm];3 Landau peak forward triplets",
		 500, -0.25, 0.25 );
  TH1I hdx3cq4t( "dx3cq4t",
		 "triplet dx, 3 Landau peak, forward;dx [mm];3 Landau peak forward triplets",
		 500, -0.25, 0.25 );
  TH1I hdx3cq3t2( "dx3cq3t2",
		  "triplet dx, 3 Landau peak, forward, 2-px;dx [mm];3 Landau peak forward triplets",
		  500, -0.25, 0.25 );

  TH1I hetaB3( "etaB3", "B cluster eta;eta;B 2-pix clusters on tracks",
	       100, -1, 1 );
  TH1I hpxqA3( "pxqA3", "A pixel charge;pixel charge [ke];A pixels on tracks",
	       100, 0, 20 );
  TH1I hpxqB3( "pxqB3", "B pixel charge;pixel charge [ke];B pixels on tracks",
	       100, 0, 20 );
  TH1I hpxqC3( "pxqC3", "C pixel charge;pixel charge [ke];C pixels on tracks",
	       100, 0, 20 );
  TH1I hpxq1stB3( "pxq1stB3", "B 1st pixel charge;pixel charge [ke];B `st pixels on tracks",
	       100, 0, 20 );
  TH1I hpxq2ndB3( "pxq2ndB3", "B 2nd pixel charge;pixel charge [ke];B `st pixels on tracks",
	       100, 0, 20 );

  TProfile effvsdxy( "effvsdxy",
		     "DUT efficiency vs triplet dxy;xy match radius [mm];DUT efficiency",
		     1000, 0, 5, -0.1, 1.1 );

  TProfile2D * effvsxy =
    new TProfile2D( "effvsxy",
		    "DUT efficiency map;x [mm];y[mm];DUT efficiency",
		    80, -4, 4, 80, -4, 4, -0.1, 1.1 );
  TProfile effvsx( "effvsx", "eff vs x;x [mm];DUT efficiency",
		   320, -4, 4, -0.1, 1.1 );
  TProfile effvsy( "effvsy", "eff vs y;y [mm];DUT efficiency",
		   80, -4, 4, -0.1, 1.1 );
  TProfile effvsxm( "effvsxm", "eff vs x mod 100;x mod 100 [#mum];DUT efficiency",
		    100, 0, 100, -0.1, 1.1 );

  TProfile effvsev( "effvsev", "eff vs time;trigger;DUT efficiency",
		    200, 0, 2000*1000, -0.1, 1.1 );
  TProfile effvsiev( "effvsiev", "eff vs event;event mod 100;DUT efficiency",
		    100, -0.5, 99.5, -0.1, 1.1 );
  TProfile effvsmpxA( "effvsmpxA", "eff vs occupancy A;occupancy A [pixels];DUT efficiency",
		    50, 0.5, 50.5, -0.1, 1.1 );
  TProfile effvsqA( "effvsqA", "eff vs charge A;cluster charge A [ke];DUT efficiency",
		    100, 0, 100, -0.1, 1.1 );
  TProfile effvstxy( "effvstxy", "eff vs angle;dxy CA [mm];DUT efficiency",
		    100, 0, 0.2, -0.1, 1.1 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Read file by lines:

  string START {"START"};
  string hd; // header

  while( hd != START ) {
    getline( Astream, hd ); // read one line into string
    cout << "  " << hd << endl;
  }

  hd.clear();
  while( hd != START ) {
    getline( Bstream, hd ); // read one line into string
    cout << "  " << hd << endl;
  }

  hd.clear();
  while( hd != START ) {
    getline( Cstream, hd ); // read one line into string
    cout << "  " << hd << endl;
  }

  string F {"F"}; // filled flag
  string A {"A"}; // added flag

  int nev = 0;
  int maxpx = 0;
  bool ldb = 0;

  while( Astream.good() && ! Astream.eof() &&
	 Bstream.good() && ! Bstream.eof() &&
	 Cstream.good() && ! Cstream.eof() &&
	 nev < Nev ) {

    string evseed;
    getline( Astream, evseed ); // read one line into string

    istringstream iss( evseed ); // tokenize string
    int iev;
    iss >> iev;

    if( ldb )
      cout << "ev " << iev;
    else
      if( iev%1000 == 0 )
	cout << "  ev " << iev << flush;

    int mpxA = 0; // event pixels

    string filledA;
    iss >> filledA;
    if( filledA == F ) {

      string roi;
      getline( Astream, roi );
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

	hphA.Fill( ph );

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

	hdphA.Fill( dph ); // sig 2.7

	//if( dph > 20 ) { // old adapter
	//if( dph > 16 ) { // 
	if( dph > 12 ) { // 

	  hpxmapA->Fill( col4, row4 );

 	  pb[mpxA].col = (col4+1)/2; // 100 um

	  if( col4%2 ) 
	    pb[mpxA].row = 2*row4 + 0;
	  else
	    pb[mpxA].row = 2*row4 + 1;

	  pb[mpxA].ph = dph;

	  // r4scal.C

	  double U = ( dph - Ap3[col4][row4] ) / Ap2[col4][row4];

	  if( U >= 1 )
	    U = 0.9999999; // avoid overflow

	  double vcal = Ap0[col4][row4] - Ap1[col4][row4] * log( (1-U)/U ); // inverse Fermi

	  pb[mpxA].q = Ake*vcal;

	  ++mpxA;
	}

      } // ipx

      if( ldb ) cout << " roi " << npx << ", hits " << mpxA << endl;

    } // filled

    else
      if( ldb ) cout << "  empty" << endl;

    // clustering:

    hnpxA.Fill( mpxA );

    if( mpxA > maxpx ) {
      cout << " A maxpx " << mpxA << " in " << nev;
      maxpx = mpxA;
    }

    if( mpxA > 50 ) mpxA = 0; // speed

    fNHit = mpxA; // for cluster search

    vector <cluster> vclA = getClus();	    

    hnclA.Fill( vclA.size() );
    if( vclA.size() )
      if( ldb ) cout << "  clusters " << vclA.size() << endl;

    for( unsigned icl = 0; icl < vclA.size(); ++ icl ) {

      hclmapA->Fill( vclA[icl].col, vclA[icl].row );

      if( vclA[icl].size < 6 ) {
	hclphA.Fill( vclA[icl].sum );
	hclqA.Fill( vclA[icl].q );
      }

      if( vclA[icl].sum > 55 )
	hclszA.Fill( vclA[icl].size );

    } // icl

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // B:

    getline( Bstream, evseed ); // read one line into string

    istringstream issB( evseed ); // tokenize string

    issB >> iev;

    if( ldb ) cout << "B ev " << iev;

    int mpxB = 0; // event pixels

    string filledB;
    issB >> filledB;
    if( filledB == F ) {

      string roi;
      getline( Bstream, roi );
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

	hphB.Fill( ph );

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

	// simulate n-bit ADC: run 449

	//int nb = 256; // 446 dx3cq3t  4.10
	//int nb = 128; // 446 dx3cq3t  4.12
	//int nb = 64; // 446 dx3cq3t  4.16
	//int nb = 32; // 446 dx3cq3t  4.27
	//int nb = 16; // 446 dx3cq3t  4.63
	//int nb =  8; // 446 dx3cq3t  5.5
	//int nb =  4; // 446 dx3cq3t  5.94
	//int nb =  2; // 446 dx3cq3t  6.5 ( 25/sqrt(12) = 7.2 )

	//dph = int( dph/400*nb) * 400/nb;

	//if( dph > 400 ) dph = 400; // overflow bit

	hdphB.Fill( dph ); // sig 2.7

	// r4scal.C

	double U = ( dph - Bp3[col4][row4] ) / Bp2[col4][row4];

	if( U >= 1 )
	  U = 0.9999999; // avoid overflow

	double vcal = Bp0[col4][row4] - Bp1[col4][row4] * log( (1-U)/U ); // inverse Fermi

	double q = Bke*vcal;

	//if( dph > 20 ) { // sig 4.4
	//if( dph > 16 ) { // 449 dx3cq3t sig 4.27
	//if( dph > 12 ) { // 449 dx3cq3t 449 sig 4.24
	//if( dph > 10 ) { // sig 4.3, 434 dx3cq3t 4.9

	//if( q > 0.8 ) { // 446 dx3cq3t  5.2
	//if( q > 1.0 ) { // 446 dx3cq3t  4.3
	if( q > 1.2 ) { // 446 dx3cq3t  4.11
	//if( q > 1.4 ) { // 446 dx3cq3t  4.14
	//if( q > 1.7 ) { // 446 dx3cq3t  4.25
	//if( q > 2.0 ) { // 446 dx3cq3t  4.41
	//if( q > 2.5 ) { // 446 dx3cq3t  4.71
	//if( q > 3.0 ) { // 446 dx3cq3t  5.03
	//if( q > 3.5 ) { // 446 dx3cq3t  5.39
	//if( q > 4.0 ) { // 446 dx3cq3t  5.78
	//if( q > 4.5 ) { // 446 dx3cq3t  6.25
	//if( q > 5.0 ) { // 446 dx3cq3t  6.46
	//if( q > 6.0 ) { // 446 dx3cq3t  7.02

	  //q = 5; // binary above threshold

	  hpxmapB->Fill( col4, row4 );

 	  pb[mpxB].col = (col4+1)/2; // 100 um

	  if( col4%2 ) 
	    pb[mpxB].row = 2*row4 + 0;
	  else
	    pb[mpxB].row = 2*row4 + 1;

	  pb[mpxB].ph = dph;

	  pb[mpxB].q = q;

	  ++mpxB;
	}

      } // ipx

      if( ldb ) cout << " roi " << npx << ", hits " << mpxB << endl;

    } // filled

    else
      if( ldb ) cout << "  empty" << endl;

    // clustering:

    hnpxB.Fill( mpxB );

    if( mpxB > maxpx ) {
      cout << " B maxpx " << mpxB << " in " << nev;
      maxpx = mpxB;
    }

    fNHit = mpxB; // for cluster search

    vector <cluster> vclB = getClus();	    

    hnclB.Fill( vclB.size() );
    if( vclB.size() )
      if( ldb ) cout << "  clusters " << vclB.size() << endl;

    for( unsigned icl = 0; icl < vclB.size(); ++ icl ) {

      hclmapB->Fill( vclB[icl].col, vclB[icl].row );

      if( vclB[icl].size < 6 ) {
	hclphB.Fill( vclB[icl].sum );
	hclqB.Fill( vclB[icl].q );
      }

      if( vclB[icl].sum > 55 )
	hclszB.Fill( vclB[icl].size );

    } // icl

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // C:

    getline( Cstream, evseed ); // read one line into string

    istringstream issC( evseed ); // tokenize string

    issC >> iev;

    if( ldb ) cout << "C ev " << iev;

    int mpxC = 0; // event pixels

    string filledC;
    issC >> filledC;
    if( filledC == F ) {

      string roi;
      getline( Cstream, roi );
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

	hphC.Fill( ph );

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

	hdphC.Fill( dph ); // sig 2.7

	//if( dph > 20 ) { // sig 4.6
	//if( dph > 16 ) { // sig 4.5
	if( dph > 12 ) { // sig 4.5

	  hpxmapC->Fill( col4, row4 );

 	  pb[mpxC].col = (col4+1)/2; // 100 um

	  if( col4%2 ) 
	    pb[mpxC].row = 2*row4 + 0;
	  else
	    pb[mpxC].row = 2*row4 + 1;

	  pb[mpxC].ph = dph;

	  // r4scal.C

	  double U = ( dph - Cp3[col4][row4] ) / Cp2[col4][row4];

	  if( U >= 1 )
	    U = 0.9999999; // avoid overflow

	  double vcal = Cp0[col4][row4] - Cp1[col4][row4] * log( (1-U)/U ); // inverse Fermi

	  pb[mpxC].q = Cke*vcal;

	  ++mpxC;
	}

      } // ipx

      if( ldb ) cout << " roi " << npx << ", hits " << mpxC << endl;

    } // filled

    else
      if( ldb ) cout << "  empty" << endl;

    // clustering:

    hnpxC.Fill( mpxC );

    if( mpxC > maxpx ) {
      cout << " C maxpx " << mpxC << " in " << nev;
      maxpx = mpxC;
    }

    if( mpxC > 50 ) mpxC = 0; // speed

    fNHit = mpxC; // for cluster search

    vector <cluster> vclC = getClus();	    

    hnclC.Fill( vclC.size() );
    if( vclC.size() )
      if( ldb ) cout << "  clusters " << vclC.size() << endl;

    for( unsigned icl = 0; icl < vclC.size(); ++ icl ) {

      hclmapC->Fill( vclC[icl].col, vclC[icl].row );

      if( vclC[icl].size < 6 ) {
	hclphC.Fill( vclC[icl].sum );
	hclqC.Fill( vclC[icl].q );
      }

      if( vclC[icl].sum > 55 )
	hclszC.Fill( vclC[icl].size );

    } // icl

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // A-B cluster correlations:

    for( vector<cluster>::iterator cA = vclA.begin(); cA != vclA.end(); ++cA ) {

      double xA = cA->row*0.025 - 4.0 - alignxA; // rot90 Dreimaster
      double yA = cA->col*0.100 - 3.9 - alignyA; // 100 um px

      double xAr = xA*cfA - yA*sfA;
      double yAr = xA*sfA + yA*cfA;

      for( vector<cluster>::iterator cB = vclB.begin(); cB != vclB.end(); ++cB ) {

	double xB = cB->row*0.025 - 4.0;
	double yB = cB->col*0.100 - 3.9;

	hxxAB.Fill( xAr, xB );
	hyyAB.Fill( yAr, yB );

	double dx = xAr - xB;
	double dy = yAr - yB;
	hdxAB.Fill( dx );
	hdyAB.Fill( dy );
	dxvsxAB.Fill( xB, dx );
	dxvsyAB.Fill( yB, dx );

	hdxvsev.Fill( nev, dx );

      } // clusters

    } // cl

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // B-C cluster correlations:

    for( vector<cluster>::iterator cC = vclC.begin(); cC != vclC.end(); ++cC ) {

      double xC = cC->row*0.025 - 4.0 - alignxC; // rot90 Dreimaster
      double yC = cC->col*0.100 - 3.9 - alignyC; // down

      double xCr = xC*cfC - yC*sfC;
      double yCr = xC*sfC + yC*cfC;

      for( vector<cluster>::iterator cB = vclB.begin(); cB != vclB.end(); ++cB ) {

	double xB = cB->row*0.025 - 4.0;
	double yB = cB->col*0.100 - 3.9;

	hxxCB.Fill( xB, xCr );
	hyyCB.Fill( yB, yCr );

	double dx = xCr - xB;
	double dy = yCr - yB;
	hdxCB.Fill( dx );
	hdyCB.Fill( dy );
	dxvsxCB.Fill( xB, dx );
	dxvsyCB.Fill( yB, dx );

      } // clusters

    } // cl

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // A-C cluster correlations:

    for( vector<cluster>::iterator cA = vclA.begin(); cA != vclA.end(); ++cA ) {

      double xA = cA->row*0.025 - 4.0 - alignxA; // rot90 Dreimaster
      double yA = cA->col*0.100 - 3.9 - alignyA; // down

      double xAr = xA*cfA - yA*sfA;
      double yAr = xA*sfA + yA*cfA;

      for( vector<cluster>::iterator cC = vclC.begin(); cC != vclC.end(); ++cC ) {

	double xC = cC->row*0.025 - 4.0 - alignxC;
	double yC = cC->col*0.100 - 3.9 - alignyC;

	double xCr = xC*cfC - yC*sfC;
	double yCr = xC*sfC + yC*cfC;

	hxxCA.Fill( xAr, xCr );
	hyyCA.Fill( yAr, yCr );

	double dxCA = xCr - xAr;
	double dyCA = yCr - yAr;
	double dxyCA = sqrt( dxCA*dxCA + dyCA*dyCA );
	hdxCA.Fill( dxCA );
	hdyCA.Fill( dyCA );

	if( fabs( dxCA ) > 0.10 ) continue; // includes beam divergence: +-3 sigma

	hdyCAc.Fill( dyCA );
	dyvsyCA.Fill( yAr, dyCA );

	if( fabs( dyCA ) > 0.25 ) continue;

	dxvsyCA.Fill( yAr, dxCA );
	dxvsxCA.Fill( xAr, dxCA ); // linear trend in run 392, 403: acceptance and beam divergence ?

	double xavg = 0.5 * ( xAr + xCr );
	double yavg = 0.5 * ( yAr + yCr );
	double xmod = fmod( xavg + 4, 0.1 ); // [mm] 0..0.1

	int eff[999] = {0};

	for( vector<cluster>::iterator cB = vclB.begin(); cB != vclB.end(); ++cB ) {

	  double xB = cB->row*0.025 - 4.0;
	  double yB = cB->col*0.100 - 3.9;

	  double eta = -2;
	  if( cB->size == 2 ) {
	    double q0 = cB->vpix[0].q;
	    double q1 = cB->vpix[1].q;
	    eta = (q1-q0)/(q1+q0);
	  }

	  int rowmin = 999;
	  int rowmax = 0;
	  for( int ipx = 0; ipx < cB->size; ++ipx ) {
	    int row = cB->vpix[ipx].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	  }
	  int nrowB = rowmax-rowmin+1;

	  // triplet residual:

	  double dx3 = xB - xavg;
	  if( run == 447 || run == 449 )
	    dx3 -= 0.00076*xavg;
	  double dy3 = yB - yavg;
	  double dxy = sqrt( dx3*dx3 + dy3*dy3 );

	  hdx3.Fill( dx3 );
	  hdy3.Fill( dy3 );

	  if( fabs( dy3 ) < 0.15 && yavg > -3.5 ) { // cut on y, look at x, see madxvsy

	    hdx3c.Fill( dx3 );
	    if( cB->size == 1 )
	      hdx3c1.Fill( dx3 ); // r447 4.4
	    if( cB->size == 2 )
	      hdx3c2.Fill( dx3 ); // r447 4.4
	    if( cB->size == 3 )
	      hdx3c3.Fill( dx3 ); // r447 4.8
	    if( cB->size == 4 )
	      hdx3c4.Fill( dx3 ); // r447 6.5
	    if( cB->size == 5 )
	      hdx3c5.Fill( dx3 ); // r447 16.6
	    if( cB->size == 6 )
	      hdx3c6.Fill( dx3 ); // r447 24.5
	    if( cB->size > 6 )
	      hdx3c7.Fill( dx3 ); // r447 39.9

	    dx3vsx.Fill( xB, dx3 ); // turn
	    if( xB < 0 )
	      hdx3m.Fill( dx3 );
	    else
	      hdx3p.Fill( dx3 );
	    dx3vsy.Fill( yB, dx3 ); // rot
	    dx3vsxm.Fill( xmod*1E3, dx3 );

	    if( fabs( dxCA ) < 0.050 ) { // track angle
	      hdx3ct.Fill( dx3 );
	      madx3vsq.Fill( cB->q, fabs(dx3) );
	      madx3vsn.Fill( cB->size, fabs(dx3) );
	    }

	    if( cB->q > 8 && cB->q < 17 ) {

	      hdx3cq.Fill( dx3 );

	      if( cA->q > 6 && cA->q < 18 &&
		  cC->q > 6 && cC->q < 18 ) {

		hdx3cq3.Fill( dx3 );

		madx3vsdx.Fill( dxCA*1E3, fabs(dx3) ); // dxCA

		if( fabs( dxCA ) < 0.050 ) { // track angle

		  hdx3cq3t.Fill( dx3 ); // 447 4.27 um

		  madx3vsx.Fill( xB, fabs(dx3) );
		  madx3vsy.Fill( yB, fabs(dx3) );
		  madx3vsxm.Fill( xmod*1E3, fabs(dx3) );
		  if( cB->size == 2 ) {
		    etavsxmB3.Fill( xmod*1E3, eta ); // sine
		    madx3vseta.Fill( eta, fabs(dx3) ); // flat
		    hdx3cq3t2.Fill( dx3 ); // 447 4.25 um
		  }

		  if( cB->q > 9 && cB->q < 14 )
		    hdx3cq4t.Fill( dx3 ); // 447 4.22 um

		} // angle

	      } // Qa, qC

	    } // qB

	  } // cut dy

	  if( fabs( dx3 ) < 0.07 && fabs( dy3 ) < 0.15 ) { // hit on track

	    hclszA3.Fill( cA->size );
	    hclszB3.Fill( cB->size );
	    hclszC3.Fill( cC->size );

	    hclphA3.Fill( cA->sum );
	    hclphB3.Fill( cB->sum );
	    hclphC3.Fill( cC->sum );

	    hclqA3.Fill( cA->q );
	    hclqB3.Fill( cB->q );
	    hclqC3.Fill( cC->q );
	    if( cB->size < 4 )
	      hclqB3n.Fill( cB->q );

	    nrowvsxmB3.Fill( xmod*1E3, nrowB );
	    clqvsxmB3.Fill( xmod*1E3, cB->q );

	    if( cB->size == 2 )
	      hetaB3.Fill( eta );

	    for( int ipx = 0; ipx < cA->size; ++ipx )
	      hpxqA3.Fill( cA->vpix[ipx].q );
	    for( int ipx = 0; ipx < cB->size; ++ipx )
	      hpxqB3.Fill( cB->vpix[ipx].q );
	    for( int ipx = 0; ipx < cC->size; ++ipx )
	      hpxqC3.Fill( cC->vpix[ipx].q );

	    if( cB->size == 2 ) {
	      hpxq1stB3.Fill( cB->vpix[0].q );
	      hpxq2ndB3.Fill( cB->vpix[1].q ); // identical
	    }

	    // task: store track

	  } // linked

	  for( int iw = 1; iw < 999; ++iw )
	    if( dxy < iw*0.020 )
	      eff[iw] = 1; // eff

	} // clusters B

	if( fabs( dyCA ) > 0.12 ) continue; // clean reference "tracks"

	if( filledB == A ) continue; // padded event

	effvsxy->Fill( xavg, yavg, eff[50] );

	if( yavg > -3.7 && yavg < 3.5 )
	  effvsx.Fill( xavg, eff[50] );

	if( xavg > -3.3 && xavg < 3.6 )
	  effvsy.Fill( yavg, eff[50] );

	if( xavg > -3.3 && xavg < 3.6 &&
	    yavg > -3.7 && yavg < 3.5 ) {

	  for( int iw = 1; iw < 999; ++iw )
	    effvsdxy.Fill( iw*0.020+0.005, eff[iw] );

	  effvsxm.Fill( xmod*1E3, eff[50] ); // bias dot
	  effvsev.Fill( nev, eff[50] );
	  effvsiev.Fill( nev%100, eff[50] );
	  effvsmpxA.Fill( mpxA, eff[50] );
	  effvsqA.Fill( cA->q, eff[50] );
	  effvstxy.Fill( dxyCA, eff[50] ); // flat

	} // fiducial x, y

      } // clusters A

    } // cl C

    // task: track correlations, intersects

    ++nev;

  } // while events

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( !ldb ) cout << endl;

  cout << "done " << Afile
       << endl << "events " << nev
       << endl << histoFile->GetName()
       << endl;

  histoFile->Write();
  histoFile->Close();

  return 0;
}
