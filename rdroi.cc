
// Daniel Pitzl (DESY) Sep 2017
// read region-of-interest data

// rdroi -n 2500 roi000102.txt
// rdroi A/roi000476.txt  # shallow copy from cmshannonb
// rdroi A/roi000477.txt  # shallow copy
// rdroi A/roi000480.txt  # shallow

// rdroi x/roi010014.txt  # X-ray
// rdroi -f x/roi010026.txt  # Sr90

// rdroi -f B/roi000798.txt  # edge-on

// rdroi -f B/roi000972.txt

// rdroi -n 546000 C/roi001019.txt

// rdroi -f B/roi001047.txt
// rdroi B/roi001525.txt
// rdroi -f B/roi001553.txt
// rdroi -f B/roi001637.txt  133i

#include <stdlib.h> // atoi
#include <iostream> // cout
#include <iomanip> // setw
#include <string> // strings
#include <sstream> // stringstream
#include <fstream> // files
#include <vector>
#include <cmath>
#include <sys/time.h> // gettimeofday, timeval

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
  int size;
  int sum;
  double q;
  double col, row;
};

//------------------------------------------------------------------------------
vector<cluster> getClus( vector <pixel> pb, int fCluCut = 1 ) // 1 = no gap
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

  //std::ios_base::sync_with_stdio( false ); // faster read? no

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

  int Nev = 2000*1000;
  bool fifty = 0;

  for( int i = 1; i < argc; i++ ) {

    if( !strcmp( argv[i], "-n" ) )
      Nev = atoi( argv[++i] );

    if( !strcmp( argv[i], "-f" ) )
      fifty = 1;

  } // argc

  // gain:

  double p0[155][160]; // Fermi
  double p1[155][160];
  double p2[155][160];
  double p3[155][160];

  //string gain{ "x/r111-scancal-x.dat" };
  //string gain{ "B/c108-scancal2-tb21-2018-02-24-ia125-hold24.dat" };
  //string gain{ "B/c128i-scancal2-tb21-icy-2018-02-23-ia115-hold20.dat"};
  //string gain{ "B/c109-scancal2-pr900-sh670-2018-02-25-hold20b.dat" };
  string gain{ "B/scm133i-scancal2-tb21-icy-pr800-sh600-ia115-2018-03-07-hold20.dat" };

  double ke = 0.036;

  ifstream gainFile( gain );

  bool haveGain = 0;

  if( gainFile ) {

    haveGain = 1;

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

  } // gain
  else
    cout << "missing " << gain << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  TFile * histoFile = new TFile( Form( "rdroi-%i.root", run ), "RECREATE" );

  // book histos:

  TH1I hph( "ph", "PH;ADC-PED [ADC];pixels", 1000, -100, 900 );
  TH1I hdph( "dph", "dPH;#DeltaPH [ADC];pixel", 1000, -100, 900 );
  TH1I hdph20( "dph20", "dPH;#DeltaPH [ADC];pixel", 1000, -100, 900 );

  int nbx =  78;
  int nby = 320;
  if( fifty ) {
    nbx = 155;
    nby = 160;
  }
  TH2I * hpxmap = new TH2I( "pxmap", "pixel map, dph > cut;col;row;pixels above cut",
			    nbx, -0.5, nbx-0.5, nby, -0.5, nby-0.5 );
  TH1I hnpx( "npx", "PH pixels per event;PH pixels;events", 80, 0.5, 80.5 );

  TH1D hncl( "ncl", "cluster per event;cluster;events", 20, 0.5,  20.5 );
  TH2I * hclmap = new TH2I( "clmap", "cluster map;col;row;clusters",
			    155, -0.5, 154.5, 160, -0.5, 159.5 );
  TH1I hclsz( "clsz", "cluster size;cluster size [pixels];clusters", 80, 0.5, 80.5 );
  TH1I hclph( "clph", "cluster PH;cluster ph [ADC];clusters", 200, 0, 1000 );
  TH1I hpxph( "pxph", "pixel PH;pixel ph [ADC];pixels in clusters", 200, 0, 400 );
  TH1I hpxph35( "pxph35", "pixel PH, long clusters;pixel ph [ADC];pixels in long clusters",
		200, 0, 400 );
  TH1I hg200( "g200", "200mV/PH;pixel gain [mv/ADC];pixels in clusters", 100, 0, 2 );
  TH1I hpxq( "pxq", "pixel charge;pixel charge [ke];pixels in clusters", 100, 0, 20 );
  TH1I hncol( "ncol", "cluster cols;cluster size [cols];clusters", 90, 0.5, 90.5 );
  TH1I hnrow( "nrow", "cluster rows;cluster size [rows];clusters", 90, 0.5, 90.5 );
  TH1I hdiag( "diag", "rows/cols;cluster diagonal [rows/cols];clusters", 100, 0, 2 );
  TH1I hcolph( "colph", "column PH;column ph [ADC];columns", 100, 0, 500 );
  TH1I hclq( "clq", "cluster charge;cluster charge [ke];cluster", 100, 0, 20 );
  TH1I hclq1( "clq1", "cluster charge;cluster charge [ke];cluster", 200, 0, 100 );
  TH1I hclcol( "clcol", "cluster column;column;clusters", nbx, -0.5, nbx-0.5 );
  TH1I hclrow( "clrow", "cluster row;row;clusters", nby, -0.5, nby-0.5 );
  TProfile clqvscol( "clqvscol", "Q vs columnnnnn;column;<cluster charge> [ke]", nbx, -0.5, nbx-0.5 );

	   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Read file by lines:

  string START {"START"};
  string hd;

  while( hd != START ) {
    getline( evFile, hd ); // read one line into string
    cout << "  " << hd << endl;
  }

  timeval tv;
  gettimeofday( &tv, NULL );
  long s0 = tv.tv_sec; // seconds since 1.1.1970
  long u0 = tv.tv_usec; // microseconds

  string F {"F"}; // filled flag
  string BLANK{" "};

  string evseed;
  getline( evFile, evseed ); // read one line into string
  int nev = 0;

  while( evFile.good() && ! evFile.eof() && nev < Nev ) {

    istringstream iss( evseed ); // tokenize string
    int iev;
    iss >> iev;

    if( iev%1000 == 0 ) cout << " " << iev << flush;

    vector <pixel> pb; // for clustering

    string filled;
    iss >> filled;

    if( filled == F ) {

      int npx = 0;
      vector <pixel> vpx;
      vpx.reserve(35);

      string roi;
      getline( evFile, roi );

      //cout << roi << " (" << roi.size() << ")" << endl;

      size_t start = 0;
      size_t gap = 0;
      while( gap < roi.size()-1 ) { // data have trailing blank

        gap = roi.find( BLANK, start );
	string s1( roi.substr( start, gap - start ) );
	//cout << " " << s1 << "(" << gap << ")";
	//int col = stoi(s1);
	int col = atoi(s1.c_str()); // 5% faster
	start = gap + BLANK.size();

        gap = roi.find( BLANK, start );
	string s2( roi.substr( start, gap - start ) );
	//cout << " " << s2 << "(" << gap << ")";
	//int row = stoi(s2);
	int row = atoi(s2.c_str());
	start = gap + BLANK.size();

        gap = roi.find( BLANK, start );
	string s3( roi.substr( start, gap - start ) );
	//cout << " " << s1 << "(" << gap << ")";
	//double ph = stod(s3);
	double ph = atof(s3.c_str());
	start = gap + BLANK.size();

	pixel px { col, row, ph, ph };
	vpx.push_back(px); // comment out = no clustering
	++npx;

      }
      //cout << endl;
      /*
      istringstream css( roi ); // tokenize string
      while( ! css.eof() ) { // slower

	int col;
	int row;
	double ph;
	css >> col;
	css >> row;
	css >> ph;

	pixel px { col, row, ph, ph };
	vpx.push_back(px);
	++npx;

      } // roi stream
      */
      // columns-wise common mode correction:

      if( vpx.size() > 999 ) vpx.clear();

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
	if( col4 > 19 ) hdph20.Fill( dph ); // skip noisy region

	//if( dph > 12 ) { // gain_1
	if( dph > 24 ) { // gain_2
	  //if( dph > 40 ) { // gain_2 irrad

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

	  // r4scal.C

	  if( haveGain ) {

	    double U = ( dph - p3[col4][row4] ) / p2[col4][row4];

	    if( U >= 1 )
	      U = 0.9999999; // avoid overflow

	    double vcal = p0[col4][row4] - p1[col4][row4] * log( (1-U)/U ); // inverse Fermi
	    /*
	    double t200 = ( 200 - p0[col4][row4] ) / p1[col4][row4];
	    double ph200 = p3[col4][row4] + p2[col4][row4] / ( 1 + exp(-t200) );
	    double g200 = 200 / ph200;
	    vcal = dph*g200;
	    hg200.Fill( g200 );
	    */
	    double q = ke*vcal;

	    hpxq.Fill( q );

	    px.q = q;

	  }
	  else
	    px.q = dph;

	  pb.push_back(px);

	  hpxmap->Fill( px.col, px.row );

	} // dph

      } // ipx

      //cout << "ev " << iev << " roi " << npx << ", hits " << pb.size() << endl;

    } // filled

    //else cout << "  empty" << endl;

    // clustering:

    hnpx.Fill( pb.size() );

    vector <cluster> vcl = getClus(pb);

    hncl.Fill( vcl.size() );
    //if( vcl.size() ) cout << "  clusters " << vcl.size();

    for( unsigned icl = 0; icl < vcl.size(); ++ icl ) {

      hclmap->Fill( vcl[icl].col, vcl[icl].row );

      //cout << " size " << vcl[icl].size;

      if( vcl[icl].col > 15 ) { // 109 noisy
	hclsz.Fill( vcl[icl].size );
	hclph.Fill( vcl[icl].sum );
	hclq.Fill( vcl[icl].q );
	hclq1.Fill( vcl[icl].q );
      }
      if( vcl[icl].q > 4 ) {
	hclcol.Fill( vcl[icl].col );
	hclrow.Fill( vcl[icl].row );
	clqvscol.Fill( vcl[icl].col, vcl[icl].q );
      }

      // pixels in clusters:

      int colmin = 999;
      int colmax = 0;
      int rowmin = 999;
      int rowmax = 0;
      double colph[155] { 0 };

      for ( int ipx = 0; ipx < vcl[icl].size; ++ipx ) {

	int col = vcl[icl].vpix[ipx].col;
	if( col < colmin ) colmin = col;
	if( col > colmax ) colmax = col;

	int row = vcl[icl].vpix[ipx].row;
	if( row < rowmin ) rowmin = row;
	if( row > rowmax ) rowmax = row;

	colph[col] += vcl[icl].vpix[ipx].ph;
	hpxph.Fill( vcl[icl].vpix[ipx].ph );

	if( vcl[icl].size > 31 && vcl[icl].size < 39 ) // shallow road
	  hpxph35.Fill( vcl[icl].vpix[ipx].ph );

      } // px

      int ncol = colmax - colmin + 1;
      hncol.Fill( ncol );

      int nrow = rowmax - rowmin + 1;
      hnrow.Fill( nrow );

      double diag = (double)nrow / (double)ncol;

      if( ncol > 4 ) {
	hdiag.Fill( diag ); // peak at 1

	for( int col = colmin; col <= colmax; ++col )
	  hcolph.Fill( colph[col] );
      }

      //if( vcl.size() ) cout << endl;

    } // cl

    ++nev;

    getline( evFile, evseed ); // read ahead

  } // while events

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  gettimeofday( &tv, NULL );
  long s9 = tv.tv_sec; // seconds since 1.1.1970
  long u9 = tv.tv_usec; // microseconds

  cout << endl << "done " << evFileName
       << endl << nev << " events"
       << " in " << s9 - s0 + ( u9 - u0 ) * 1e-6 << " s"
       << endl << histoFile->GetName()
       << endl;

  histoFile->Write();
  histoFile->Close();

  cout << "have gain " << haveGain << endl;

  return 0;
}
