
// Daniel Pitzl (DESY) Oct 2017
// R4S edge-on

// edges B/roi000612.txt
// edges -c 110  B/roi000631.txt
// edges -c 114  B/roi000634.txt
// edges -c 109  B/roi000637.txt
// edges -c 109  B/roi000648.txt
// edges -c 109  B/roi000657.txt  rot 13.5  2732k

// edges -c 117 -f B/roi000662.txt

// edges -c 1092  B/roi000689.txt  rot 4

// edges -f -c 118  B/roi000814.txt # edge-on 50x50

// edges -f -c 144  B/roi000946.txt # shallow

// edges -c 124  B/roi001388.txt # irrad

// edges -c 108  B/roi001491.txt

// edges -c 136  B/roi001655.txt

#include <cstdlib> // atoi
#include <iostream> // cout
#include <iomanip> // setw
#include <string> // strings
#include <sstream> // stringstream
#include <fstream> // files
#include <vector>
#include <set>
#include <cmath>
#include <time.h> // clock_gettime

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>

using namespace std;

class pixel {
 public:
  int col;
  int row;
  double ph;
  double q;
};

class pixsrt {
public:
  bool operator() ( const pixel p1,  const pixel p2 )
  {
    return p1.row < p2.row;
  }
};

class cluster {
public:
  vector <pixel> vpix;
  int size;
  double sum;
  double q;
  double col, row;
};

class clusrt {
public:
  bool operator() ( const cluster c1,  const cluster c2 )
  {
    return c1.row < c2.row;
  }
};

//------------------------------------------------------------------------------
unsigned digit_value( char c )
{
  return unsigned( c - '0' ); // negatives get flipped to large positives
}

//------------------------------------------------------------------------------
int fast_atoi( const char * p )
{
  bool neg = false;
  if( *p == '-' ) {
    neg = true;
    ++p;
  }
  //cout << " " << *p << flush;
  int x = digit_value( *p ); // first digit, must be present
  unsigned d;
  while( ( d = digit_value( *++p ) ) <= 9 ) {
    x = x*10 + d;
    //cout << " " << *p << flush;
  }
  if( neg )
    x = -x;

  return x;
}

//------------------------------------------------------------------------------
float fast_atof( string f )
{
  string DOT{"."};
  int idot = f.find( DOT, 0 );
  string pre = f.substr( 0, idot );
  string dec = f.substr( idot+1 );
  float x = fast_atoi( pre.c_str() );
  string SIGN{"-"};
  if( f.substr( 0, 1 ) == SIGN )
    x -= fast_atoi( dec.c_str() ) / pow( 10, dec.size() );
  else
    x += fast_atoi( dec.c_str() ) / pow( 10, dec.size() );
  //cout << "  " << f << " (" << pre << "." << dec << " " << x << ")";
  return x;
}

//------------------------------------------------------------------------------
vector <cluster> getClus( vector <pixel> pb, int fCluCut = 1 ) // 1 = no gap
{
  // returns clusters with local coordinates
  // next-neighbour topological clustering (allows fCluCut-1 empty pixels)

  vector <cluster> v;
  if( pb.size() == 0 ) return v;

  int * gone = new int[pb.size()] {0};

  unsigned seed = 0;

  while( seed < pb.size() ) {

    // start a new cluster

    cluster c;
    c.vpix.push_back( pb[seed] );
    gone[seed] = 1;

    // let it grow as much as possible:

    int growing;
    do{
      growing = 0;
      for( unsigned i = 0; i < pb.size(); ++i ) {
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
    c.q = 0;
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

    while( ( ++seed < pb.size() ) && gone[seed] );

  } // while over seeds

  // nothing left,  return clusters

  delete gone;
  return v;
}

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  timespec ts;
  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s0 = ts.tv_sec; // seconds since 1.1.1970
  long f0 = ts.tv_nsec; // nanoseconds

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

  // B/roi000634.txt

  int run = stoi( evFileName.substr( 5, 6 ) );
  cout << evFileName << "  " << evFileName.substr( 5, 6 ) << "  " << run << endl;

  // further arguments:

  int Nev = 20*1000*1000;
  bool fifty = 0;
  int chip = 108;

  for( int i = 1; i < argc; i++ ) {

    if( !strcmp( argv[i], "-c" ) )
      chip = atoi( argv[++i] );

    if( !strcmp( argv[i], "-n" ) )
      Nev = atoi( argv[++i] );

    if( !strcmp( argv[i], "-f" ) )
      fifty = 1;

  } // argc

  // gain:

  double ke = 0.036; // edge-on Landau peak at 7.5 ke in 100 um

  double p0[155][160]; // Fermi
  double p1[155][160];
  double p2[155][160];
  double p3[155][160];

  string gain{ "B/r108-scancal-tb24-1024.dat" };

  if( chip == 108 )
    gain = "B/c108-scancal2-tb21-2018-02-24-ia125-hold24.dat";

  if( chip == 109 )
    gain = "B/r109-scancal-tb24-1025.dat";

  if( chip == 1092 ) {
    gain = "B/r109-scancal-tb24-1027.dat";
    ke = 0.033; // run 689 edges4 at 30 ke
  }

  if( chip == 110 )
    //gain = "B/r110-scancal-tb21-0928-hold24.dat";
    gain = "B/r110-scancal-tb24-1025.dat";

  if( chip == 114 )
    gain = "B/r114-scancal-tb24-1025.dat";

  if( chip == 117 )
    gain = "B/r117-scancal-tb24-1027.dat";

  if( chip == 118 ) {
    gain = "B/r118-scancal-tb21-1105.dat";
    ke = 0.033; // edge-on Landau peak at 30 ke in 400 um for c118
  }

  if( chip == 124 ) {
    gain = "B/c124i-scancal2-tb21-icy-ia150-hold28-2018-02-21.dat";
    ke = 0.035; // default
  }

  if( chip == 136 ) {
    gain = "B/scm136i-scancal2-tb21-icy-pr800-sh600-ia125-2018-03-08-hold20.dat";
    ke = 0.035; // default
  }

  if( chip == 144 )
    //gain = "B/r144-scancal-tb21-1109.dat";
    //gain = "B/r144-scancal-tb21-1110.dat";
    gain = "B/r144-scancal-tb21-1205.dat";

  if( chip == 332 )
    gain = "B/r332-scancal-tb21-1030.dat";

  ifstream gainFile( gain );

  if( ! gainFile ) {
    cout << "gain file not found" << endl;
    return 1;
  }

  cout << "gain " << gain << endl;

  while( ! gainFile.eof() ) {

    int icol;
    int irow;
    gainFile >> icol;
    gainFile >> irow;
    gainFile >> p0[icol][irow];
    gainFile >> p1[icol][irow];
    gainFile >> p2[icol][irow];
    gainFile >> p3[icol][irow];

  } // while

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  TFile * histoFile = new TFile( Form( "edges-%i.root", run ), "RECREATE" );

  // book histos:

  TH1I hph( "ph", "PH;ADC-PED [ADC];pixels", 1000, -100, 900 );
  TH1I hdph( "dph", "dPH;#DeltaPH [ADC];pixel", 1000, -100, 900 );
  TH1I hpxph( "pxph", "pixel PH;pixel ph [ADC];pixels in clusters", 200, 0, 400 );
  TH1I hpxq( "pxq", "pixel charge;pixel charge [ke];pixels in clusters", 100, 0, 25 );

  int nbx =  78;
  int nby = 320;
  if( fifty ) {
    nbx = 155;
    nby = 160;
  }

  double pitch = 25; // [um]
  if( fifty )
    pitch = 50;

  TH2I * hpxmap = new TH2I( "pxmap", "pixel map, dph > cut;col;row;pixels above cut",
			    nbx, -0.5, nbx-0.5, nby, -0.5, nby-0.5 );
  TH1I hnpx( "npx", "PH pixels per event;PH pixels;events", 200, 0.5, 200.5 );

  TH1D hncl( "ncl", "cluster per event;cluster;events", 20, 0.5,  20.5 );
  TH2I * hclmap = new TH2I( "clmap", "cluster map;col;row;clusters",
			    155, -0.5, 154.5, 160, -0.5, 159.5 );
  TH1I hclph( "clph", "cluster PH;cluster ph [ADC];clusters", 200, 0, 1000 );
  TH1I hclsz( "clsz", "cluster size;cluster size [pixels];clusters", 200, 0.5, 200.5 );

  TH1I hncol( "ncol", "cluster cols;cluster size [cols];clusters", nbx+5, 0.5, nbx+5.5 );
  TH1I hnrow( "nrow", "cluster rows;cluster size [rows];clusters", nby+5, 0.5, nby+5.5 );
  TH1I haspect( "aspect", "rows/cols;cluster aspectonal [rows/cols];clusters", 100, 0, 2 );
  TH1I hcolsz( "colsz", "column size;column size [rows];columns", 11, -0.5, 10.5 );
  TH1I hpixph( "pixph", "pixel ADC PH;pixel ph [ADC];pixels in clusters", 200, 0, 400 );
  TH1I hcolph( "colph", "column ADC PH;column ph [ADC];columns", 200, 0, 800 );
  TH1I hcolq1( "colq1", "column charge;column charge [ke];columns", 200, 0, 100 );
  TH1I hcolq2( "colq2", "column charge;column charge [ke];columns", 100, 0, 20 );
  TProfile qvscol( "qvscol", "column charge;column;<column charge> [ke]",
		   nbx, -0.5, nbx-0.5, 0, 99 );

  TH1D hslp( "slp", "track slope;track angle [pixels];tracks", 400, -2, 2 );
  TProfile slpvsncol( "slpvsncol", "slope vs track length;track length [columns];<slope> [pixels]",
		    nbx, 0.5, nbx+0.5, -2, 2 );

  TH2I * hq0q1 = new TH2I( "q0q1", "charge correlation;q0 [ke];q1 [ke];column pairs",
			   100, 0, 25, 100, 0, 25 );
  TH2I * hq0q2 = new TH2I( "q0q2", "charge correlation;q0 [ke];q2 [ke];column pairs",
			   100, 0, 25, 100, 0, 25 );
  TH1I hcolq0( "colq0", "normal column charge;normal column charge [ke];columns", 200, 0, 50 );

  TH1D hdy( "dy", "triplet residual;triplet residual [um];triplets",
	    200, -100, 100 );
  TH1D hdyc( "dyc", "triplet residual;triplet residual [um];triplets",
	    200, -100, 100 );
  TH1D hdyc2( "dyc2", "triplet residual, 1-2-rows;triplet residual [um];1-2-row triplets",
	      200, -100, 100 );

  TProfile madyvsq( "madyvsq", "mad y vs q;column charge [ke];mad y [#mum]",
		    100, 0, 50, 0, 100 );
  TProfile nrowvsq( "nrowvsq", "rows vs q;column charge [ke];<rows>",
		    100, 0, 50, 0, 50 );

  TProfile madyvsx( "madyvsx", "mad y vs x;x [columns];mad y [#mum]",
		    nbx, -0.5, nbx-0.5, 0, 100 );
  TProfile madyvsy( "madyvsy", "mad y vs y;y [rows];mad y [#mum]",
		    nby, -0.5, nby-0.5, 0, 100 );
  TProfile dyvsym( "dyvsym", "dy vs ymod;y mod 50 [#mum];<#Deltay> [#mum]",
		   100, 0, 50, -50, 50 );
  TProfile madyvsym( "madyvsym", "mad y vs ymod;y mod 50 [#mum];mad y [#mum]",
		     100, 0, 50, 0, 100 );
  TProfile nrowvsym( "nrowvsym", "rows vs ymod;y mod 50 [#mum];<rows>",
		     100, 0, 50, 0, 10 );

  TH1D hdycq( "dycq", "triplet residual;triplet residual [um];Q peak triplets",
	     200, -100, 100 );

  TH1D hdycqx1( "dycqx1", "triplet residual, col < 2;triplet residual [um];col < 2 Q peak triplets",
		200, -100, 100 );
  TH1D hdycqx2( "dycqx2",
		"triplet residual, 2 < col < 30;triplet residual [um];2 < col < 30 Q peak triplets",
		200, -100, 100 );
  TH1D hdycqx3( "dycqx3",
		"triplet residual, 30 < col < 50;triplet residual [um];30 < col < 50 Q peak triplets",
		200, -100, 100 );
  TH1D hdycqx4( "dycqx4", "triplet residual, 50 < col;triplet residual [um];50 < col Q peak triplets",
		200, -100, 100 );

  TH1D hdycqm1( "dycqm1", "triplet residual, ymod < 7;triplet residual [um];ymod < 7 Q peak triplets",
		200, -100, 100 );
  TH1D hdycqm2( "dycqm2",
		"triplet residual, 7 < ymod < 8;triplet residual [um];7 < ymod < 8 Q peak triplets",
		200, -100, 100 );
  TH1D hdycqm3( "dycqm3",
		"triplet residual, 8 < ymod < 10;triplet residual [um];8 < ymod < 10 Q peak triplets",
		200, -100, 100 );
  TH1D hdycqm4( "dycqm4",
		"triplet residual, 10 < ymod < 11;triplet residual [um];10 < ymod < 11 Q peak triplets",
		200, -100, 100 );
  TH1D hdycqm5( "dycqm5",
		"triplet residual, 11 < ymod < 12;triplet residual [um];11 < ymod < 12 Q peak triplets",
		200, -100, 100 );
  TH1D hdycqm6( "dycqm6",
		"triplet residual, 12 < ymod < 13;triplet residual [um];12 < ymod < 13 Q peak triplets",
		200, -100, 100 );

  TH1D hdycq3( "dycq3", "triplet residual;triplet residual [um];Q peak triplets",
	       200, -50, 50 );
  TH1D hdycq4( "dycq4", "triplet residual;triplet residual [um];Q peak triplets",
	       200, -50, 50 );
  TH1D heta( "eta", "charge sharing;eta;2-row clusters",
	     200, -1, 1 );
  TProfile madyvseta( "madyvseta", "mad y vs eta;eta;mad dy [#mum]",
		      80, -0.8, 0.8, 0, 100 );
  TH1D hdyc2eta( "dyc2eta", "triplet residual;triplet residual [um];low eta triplets",
	      200, -100, 100 );

  TH1I hpixphn( "pixphn", "pixel new ADC PH;pixel ph [ADC];pixels in clusters", 200, 0, 200 );
  TH1I hcolphn( "colphn", "column new ADC PH;column ph [ADC];columns", 200, 0, 200 );

  TH1D hdu( "du", "nonet residual;nonet residual [um];nonets",
	    200, -50, 50 );
  TH1D hduc( "duc", "nonet residual;nonet residual [um];nonets",
	    200, -50, 50 );
  TH1D hduc2( "duc2", "nonet residual, 1-2-rows;nonet residual [um];1-2-row nonets",
	      200, -50, 50 );

  TProfile maduvsq( "maduvsq", "mad y vs q;column charge [ke];mad y [#mum]",
		    100, 0, 50, 0, 100 );

  TProfile maduvsx( "maduvsx", "mad y vs x;x [columns];mad y [#mum]",
		    nbx, -0.5, nbx-0.5, 0, 100 );
  TProfile maduvsy( "maduvsy", "mad y vs y;y [rows];mad y [#mum]",
		    nby, -0.5, nby-0.5, 0, 100 );
  TProfile duvsum( "duvsum", "du vs umod;y mod 50 [#mum];<#Deltay> [#mum]",
		   100, 0, 50, -50, 50 );
  TProfile maduvsum( "maduvsum", "mad y vs umod;y mod 50 [#mum];mad y [#mum]",
		     100, 0, 50, 0, 100 );
  TProfile nrowvsum( "nrowvsum", "rows vs umod;y mod 50 [#mum];<rows>",
		     100, 0, 50, 0, 10 );

  TH1D hducq( "ducq", "nonet residual;nonet residual [um];Q peak nonets",
	     200, -50, 50 );
  TH1D hducq2( "ducq2", "nonet residual;nonet residual [um];Q peak nonets",
	       200, -50, 50 );

  TH1D hducq3( "ducq3", "nonet residual;nonet residual [um];Q peak nonets",
	       200, -50, 50 );
  TH1D hducq4( "ducq4", "nonet residual;nonet residual [um];Q peak nonets",
	       200, -50, 50 );
  TProfile maduvseta( "maduvseta", "mad y vs eta;eta;mad du [#mum]",
		      100, -1, 1, 0, 100 );
  TH1D hduc2eta( "duc2eta", "nonet residual;nonet residual [um];low eta nonets",
	      200, -50, 50 );

  //const double wt = atan(1.0) / 45.0; // pi/180 deg

  double qmn =  5; // cut
  double qmx = 11; // cut
  if( chip == 124 ) {
    qmn = 1;
    qmx = 6;
  }

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

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s1 = ts.tv_sec; // seconds since 1.1.1970
  long f1 = ts.tv_nsec; // nanoseconds
  cout << "setup time " << s1 - s0 + ( f1 - f0 ) * 1e-9 << " s" << endl;

  while( evFile.good() && ! evFile.eof() && nev < Nev ) {

    istringstream iss( evseed ); // tokenize string
    int iev;
    iss >> iev;

    if( iev%1000 == 0 ) cout << " " << iev << flush;

    vector <pixel> pb; // for clustering

    string filled;
    iss >> filled;

    if( filled == F ) {

      string roi;
      getline( evFile, roi );

      vector <pixel> vpx;
      vpx.reserve(35);

      string BLANK{" "};
      size_t start = 0;
      size_t gap = 0;

      while( gap < roi.size()-1 ) { // data have trailing blank

	gap = roi.find( BLANK, start );
	string s1( roi.substr( start, gap - start ) );
	//int col = stoi(s1);
	//int col = atoi( s1.c_str() ); // 4% faster
	int col = fast_atoi( s1.c_str() ); // another 6% faster
	start = gap + BLANK.size();

	gap = roi.find( BLANK, start );
	string s2( roi.substr( start, gap - start ) );
	//int row = stoi(s2);
	//int row = atoi( s2.c_str() );
	int row = fast_atoi( s2.c_str() );
	start = gap + BLANK.size();

	gap = roi.find( BLANK, start );
	string s3( roi.substr( start, gap - start ) );
	//double ph = stod(s3);
	//double ph = atof(s3.c_str());
	double ph = fast_atof(s3);
	start = gap + BLANK.size();

	pixel px { col, row, ph, ph };
	vpx.push_back(px); // comment out = no clustering

      } // roi

      /* slower
	 istringstream css( roi ); // tokenize string
	 while( ! css.eof() ) { // one line = one event

	 int col;
	 int row;
	 double ph;
	 css >> col;
	 css >> row;
	 css >> ph;

	 pixel px { col, row, ph, ph };
	 vpx.push_back(px);

	 } // roi px
      */
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // column-wise common mode correction: (takes time)

      if( vpx.size() > 999 ) vpx.clear();

      set <pixel,pixsrt> colpx[155]; // per column, sorted along row

      for( unsigned ipx = 0; ipx < vpx.size(); ++ipx ) {

	int col = vpx[ipx].col;
	int row = vpx[ipx].row;
	double ph = vpx[ipx].ph;
	double q = vpx[ipx].q;
	pixel px { col, row, ph, q };
	colpx[col].insert(px); // sorted along row
      }

      for( unsigned col = 0; col < 155; ++col ) {

	if( colpx[col].size() < 2 ) continue;

	auto px1 = colpx[col].begin();
	auto px7 = colpx[col].end(); --px7; // last
				       
	int row1 = px1->row;
	int row7 = px7->row;
	double ph1 = px1->ph;
	double ph7 = px7->ph;

	auto px4 = px1; ++px4;
	for( ; px4 != px7; ++px4 ) { // between 1 and 7, exclusively

	  int col4 = px4->col;
	  int row4 = px4->row;
	  double ph4 = px4->ph;

	  double dph;
	  if( row4 - row1 < row7 - row4 )
	    dph = ph4 - ph1;
	  else
	    dph = ph4 - ph7;

	  hph.Fill( ph4 );
	  hdph.Fill( dph );

	  // r4scal.C

	  double U = ( dph - p3[col4][row4] ) / p2[col4][row4];

	  if( U >= 1 )
	    U = 0.9999999; // avoid overflow

	  double vcal = p0[col4][row4] - p1[col4][row4] * log( (1-U)/U ); // inverse Fermi

	  double q = ke*vcal;

	  //double dphcut = 12; // 648 dyc 1.90
	  double dphcut = 10; // 648 duc 1.84
	  if( chip == 124 )
	    dphcut = 30; // gain_2 noisy irrad
	    //dphcut = 40; // gain_2 noisy irrad
	  if( dph > dphcut ) {
	    //if( q > 1.0 ) {
	    //if( q > 0.8 ) { // noisy

	    pixel px;

	    if( fifty ) {
	      px.col = col4;
	      px.row = row4;
	    }
	    else{
	      px.col = (col4+1)/2; // 100 um
	      if( col4%2 ) 
		px.row = 2*row4 + 0;
	      else
		px.row = 2*row4 + 1;
	    }

	    px.ph = dph;
	    px.q = q;
	    pb.push_back(px);

	    hpxph.Fill( dph );
	    hpxq.Fill( q );
	    hpxmap->Fill( px.col, px.row );

	  } // dph

	} // p4

      } // cols

      //cout << "ev " << iev << " roi " << vpx.size() << ", hits " << pb.size() << endl;

    } // filled

    //else cout << "  empty" << endl;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // clusters:

    hnpx.Fill( pb.size() );

    vector <cluster> vcl = getClus(pb); // clustering

    hncl.Fill( vcl.size() );
    //if( vcl.size() ) cout << "  clusters " << vcl.size();

    for( unsigned icl = 0; icl < vcl.size(); ++ icl ) {

      hclmap->Fill( vcl[icl].col, vcl[icl].row );

      //cout << " size " << vcl[icl].size;

      if( vcl[icl].size < 6 )
	hclph.Fill( vcl[icl].sum );

      if( vcl[icl].sum > 55 )
	hclsz.Fill( vcl[icl].size );

      // pixels in the cluster:

      int colmin = 999;
      int colmax = 0;
      int rowmin = 999;
      int rowmax = 0;
      vector < pixel > colvec[155]; // 155 col max, array of vectors

      for ( int ipx = 0; ipx < vcl[icl].size; ++ipx ) {

	int col = vcl[icl].vpix[ipx].col;
	if( col < colmin ) colmin = col;
	if( col > colmax ) colmax = col;

	int row = vcl[icl].vpix[ipx].row;
	if( row < rowmin ) rowmin = row;
	if( row > rowmax ) rowmax = row;

	colvec[col].push_back( vcl[icl].vpix[ipx] ); // add px to col

      } // px

      int ncol = colmax - colmin + 1;
      hncol.Fill( ncol );

      int nrow = rowmax - rowmin + 1;
      hnrow.Fill( nrow );

      double aspect = (double)nrow / (double)ncol;

      if( ncol > 5 ) { // skip noise artifact from roi columns = 5

	haspect.Fill( aspect ); // peak at 1

	double yvec[155];
	double qvec[155];

	for( int col = colmin; col <= colmax; ++col ) {

	  hcolsz.Fill( colvec[col].size() );

	  if( colvec[col].size() == 0 ) continue; // gap in cluster ?

	  double sumph = 0;
	  double sumq = 0;
	  double sumqy = 0;

	  for( unsigned ir = 0; ir < colvec[col].size(); ++ir ) { // rows (pixels) in this col
	    double ph = colvec[col][ir].ph;
	    sumph += ph;
	    hpixph.Fill( ph );
	    double q = colvec[col][ir].q;
	    sumq += q;
	    sumqy += q * colvec[col][ir].row;
	  } // rows

	  if( col > colmin && col < colmax ) {
	    hcolph.Fill( sumph );
	    hcolq1.Fill( sumq );
	    hcolq2.Fill( sumq );
	  }
	  qvscol.Fill( col, sumq );
	  qvec[col] = sumq;
	  yvec[col] = sumqy / sumq; // weighted average y coordinate [px]

	} // cols

	double x0 = colmin;
	double x9 = colmax;
	double y0 = yvec[colmin];
	double y9 = yvec[colmax];
	double slp = (y9-y0) / (x9-x0);
	hslp.Fill( slp ); // peak at -1
	slpvsncol.Fill( ncol, slp ); // flat
	double turn = atan( 0.25*slp ); // signed
	if( fifty )
	  turn = atan( slp ); // signed
	double norm = cos(turn); // positive

	// triplets:

	for( int col = colmin+1; col < colmax-2; ++col ) { // skip 1st and lst

	  if( colvec[col+0].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+1].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+2].size() == 0 ) continue; // gap in cluster ?

	  double y0 = yvec[col+0];
	  double y1 = yvec[col+1];
	  double y2 = yvec[col+2];
	  double ay = 0.5*(y0+y2);
	  double dy = y1 - ay;
	  double ymod = fmod( (ay+0.5)*25, 50 );
	  double ymod25 = fmod( (ay+0.5)*25, 25 );

	  hdy.Fill( dy*pitch ); // [um]	    dphcut = 40; // gain_2 noisy irrad


	  hq0q1->Fill( qvec[col+0], qvec[col+1] ); // Landau correlations
	  hq0q2->Fill( qvec[col+0], qvec[col+2] );
	  hcolq0.Fill( qvec[col+1]*norm ); // peak at 7.5 ke

	  if( qvec[col+0]*norm > qmn &&
	      qvec[col+2]*norm > qmn &&
	      qvec[col+0]*norm < qmx &&
	      qvec[col+2]*norm < qmx ) {

	    hdyc.Fill( dy*pitch ); // [um] 2.64

	    if( colvec[col+1].size() <= 2 )
	      hdyc2.Fill( dy*pitch ); // [um] 2.20

	    double q1 = qvec[col+1]; // [ke]

	    madyvsq.Fill( q1*norm, abs(dy)*pitch );
	    nrowvsq.Fill( q1*norm, colvec[col+1].size() );

	    if( q1*norm > qmn && q1*norm < qmx ) {

	      hdycq.Fill( dy*pitch ); // [um] 2.15

	      dyvsym.Fill( ymod, dy*pitch );
	      madyvsx.Fill( col+1, abs(dy)*pitch );
	      madyvsy.Fill( ay, abs(dy)*pitch );
	      madyvsym.Fill( ymod, abs(dy)*pitch );
	      nrowvsym.Fill( ymod, colvec[col+1].size() );

	      if(      col+1 <  2 )
		hdycqx1.Fill( dy*pitch );
	      else if( col+1 < 30 )
		hdycqx2.Fill( dy*pitch ); // 2.2
	      else if( col+1 < 50 )
		hdycqx3.Fill( dy*pitch ); // 2.1
	      else
		hdycqx4.Fill( dy*pitch ); // 2.3

	      if(      ymod25 <  7 )
		hdycqm1.Fill( dy*pitch ); // 1.91
	      else if( ymod25 <  8 )
		hdycqm2.Fill( dy*pitch );
	      else if( ymod25 < 10 )
		hdycqm3.Fill( dy*pitch ); // 3.2
	      else if( ymod25 < 11 )
		hdycqm4.Fill( dy*pitch );
	      else if( ymod25 < 12 )
		hdycqm5.Fill( dy*pitch );
	      else if( ymod25 < 13 )
		hdycqm6.Fill( dy*pitch );
	      else if( ymod25 < 14 )
		hdycqm5.Fill( dy*pitch );
	      else if( ymod25 < 15 )
		hdycqm4.Fill( dy*pitch );
	      else if( ymod25 < 17 )
		hdycqm3.Fill( dy*pitch );
	      else if( ymod25 < 18 )
		hdycqm2.Fill( dy*pitch );
	      else
		hdycqm1.Fill( dy*pitch );

	    } // q

	    if( q1*norm > 6 && q1*norm < 9 &&
		qvec[col+0]*norm > 6 &&
		qvec[col+2]*norm > 6 &&
		qvec[col+0]*norm < 9 &&
		qvec[col+2]*norm < 9 )
	      hdycq3.Fill( dy*pitch ); // [um] 1.80

	    if( q1*norm > 6.5 && q1*norm < 8 &&
		qvec[col+0]*norm > 6.5 &&
		qvec[col+2]*norm > 6.5 &&
		qvec[col+0]*norm < 8 &&
		qvec[col+2]*norm < 8 )
	      hdycq4.Fill( dy*pitch ); // [um] 1.80

	    if( colvec[col+1].size() == 2 ) {

	      double a0 = colvec[col+1][0].ph;
	      double a1 = colvec[col+1][1].ph;
	      double eta = (a1-a0)/(a1+a0);
	      heta.Fill( eta );
	      madyvseta.Fill( eta, abs(dy)*pitch ); // flat

	      if( abs(eta) < 0.4 )
		hdyc2eta.Fill( dy*pitch ); // [um] 2.85

	    } // 2-rows

	  } // cuts

	} // loop cols

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// nonets: 0+1+2+3 + 5+6+7+8 vs 4: hdu

	// simulate ADC bits:

	//int nb = pow( 2, 10 ); // 648 duc  
	int nb = pow( 2,  8 ); // 648 duc  1.76
	//int nb = pow( 2,  7 ); // 648 duc  1.77
	//int nb = pow( 2,  6 ); // 648 duc  1.79
	//int nb = pow( 2,  5 ); // 648 duc  1.86
	//int nb = pow( 2,  4 ); // 648 duc  2.045
	//int nb = pow( 2,  3 ); // 648 duc  2.77
	//int nb = pow( 2,  2 ); // 648 duc  4.94 RMS
	//int nb = pow( 2,  1 ); // 648 duc  5.46 RMS

	for( int col = colmin+1; col < colmax-8; ++col ) { // skip 1st and lst

	  if( colvec[col+0].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+1].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+2].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+3].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+4].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+5].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+6].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+7].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+8].size() == 0 ) continue; // gap in cluster ?

	  double u0 = yvec[col+0];
	  double u1 = yvec[col+1];
	  double u2 = yvec[col+2];
	  double u3 = yvec[col+3];

	  int col4 = col+4;
	  double u4 = yvec[col4];

	  // ADC:

	  double sumph = 0;
	  double sumq = 0;
	  double sumqy = 0;

	  for( unsigned ir = 0; ir < colvec[col4].size(); ++ir ) { // rows (pixels) in this col

	    double ph = colvec[col4][ir].ph;
	    int row4 = colvec[col4][ir].row;

	    int ia = ph/200 * nb;
	    if( ia > nb-1 ) ia = nb-1; // overflow
	    ph = ia * 200.0 / nb;
	    sumph += ph;

	    hpixphn.Fill( ph );

	    double U = ( ph - p3[col4][row4] ) / p2[col4][row4];
	    if( U >= 1 )
	      U = 0.9999999; // avoid overflow
	    double vcal = p0[col4][row4] - p1[col4][row4] * log( (1-U)/U ); // inverse Fermi
	    double q = ke*vcal;
	    sumq += q;
	    sumqy += q * row4;

	  } // rows

	  hcolphn.Fill( sumph );

	  u4 = sumqy / sumq; // overwrite !!

	  double u5 = yvec[col+5];
	  double u6 = yvec[col+6];
	  double u7 = yvec[col+7];
	  double u8 = yvec[col+8];

	  double au = 0.125 * ( u0+u1+u2+u3 + u5+u6+u7+u8 );
	  double du = u4 - au;
	  double umod = fmod( (au+0.5)*25, 50 );

	  hdu.Fill( du*pitch ); // [um]

	  if( qvec[col+0]*norm > qmn &&
	      qvec[col+1]*norm > qmn &&
	      qvec[col+2]*norm > qmn &&
	      qvec[col+3]*norm > qmn &&
	      qvec[col+5]*norm > qmn &&
	      qvec[col+6]*norm > qmn &&
	      qvec[col+7]*norm > qmn &&
	      qvec[col+8]*norm > qmn &&
	      qvec[col+0]*norm < qmx &&
	      qvec[col+1]*norm < qmx &&
	      qvec[col+2]*norm < qmx &&
	      qvec[col+3]*norm < qmx &&
	      qvec[col+5]*norm < qmx &&
	      qvec[col+6]*norm < qmx &&
	      qvec[col+7]*norm < qmx &&
	      qvec[col+8]*norm < qmx ) {

	    hduc.Fill( du*pitch ); // [um] 648 1.84

	    if( colvec[col+4].size() <= 2 )
	      hduc2.Fill( du*pitch ); // [um] 648 1.83

	    double q4 = qvec[col+4]; // [ke]

	    maduvsq.Fill( q4*norm, abs(du)*pitch );

	    if( q4*norm >qmn && q4*norm < qmx ) {

	      hducq.Fill( du*pitch ); // [um] 648 1.72

	      duvsum.Fill( umod, du*pitch );
	      maduvsx.Fill( col+4, abs(du)*pitch );
	      maduvsy.Fill( au, abs(du)*pitch );
	      maduvsum.Fill( umod, abs(du)*pitch );
	      nrowvsum.Fill( umod, colvec[col+4].size() );

	    } // q

	    if( q4*norm > 6 && q4 < 9*norm &&
		qvec[col+3]*norm > 6 &&
		qvec[col+5]*norm > 6 &&
		qvec[col+3]*norm < 9 &&
		qvec[col+5]*norm < 9 )
	      hducq3.Fill( du*pitch ); // [um] 648 1.61

	    if( q4*norm > 6.5 && q4 *norm< 8 &&
		qvec[col+3]*norm > 6.5 &&
		qvec[col+5]*norm > 6.5 &&
		qvec[col+3]*norm < 8 &&
		qvec[col+5]*norm < 8 )
	      hducq4.Fill( du*pitch ); // [um] 648 1.50

	    if( colvec[col+4].size() == 2 ) {

	      double a0 = colvec[col+4][0].ph;
	      double a1 = colvec[col+4][1].ph;
	      double eta = (a1-a0)/(a1+a0);
	      maduvseta.Fill( eta, abs(du)*pitch ); // flat

	      if( abs(eta) < 0.4 )
		hduc2eta.Fill( du*pitch ); // [um] 2.85

	    } // 2-rows

	  } // cuts

	} // loop cols nonet

      } // long

      //if( vcl.size() ) cout << endl;

    } // cl

    ++nev;

    getline( evFile, evseed ); // read ahead

  } // while events

  //cout << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s9 = ts.tv_sec; // seconds since 1.1.1970
  long f9 = ts.tv_nsec; // nanoseconds

  cout << "done " << evFileName
       << endl << "events " << nev
       << " (time " << s9 - s1 + ( f9 - f1 ) * 1e-9 << " s)"
       << endl << histoFile->GetName()
       << endl;

  histoFile->Write();
  histoFile->Close();

  return 0;
}
