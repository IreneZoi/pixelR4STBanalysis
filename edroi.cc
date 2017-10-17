
// Daniel Pitzl (DESY) Sep 2017
// read region-of-interest data: event display

// edroi A/roi000530.txt # 50x50 tilt 0
// edroi A/roi000531.txt # 50x50 tilt 0
// edroi A/roi000540.txt # 50x50 tilt 0

// edroi A/roi000534.txt # til 20

#include <stdlib.h> // atoi
#include <iostream> // cout
#include <iomanip> // setw
#include <string> // strings
#include <sstream> // stringstream
#include <fstream> // files
#include <vector>
#include <cmath>
#include <sys/ioctl.h>

#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2.h>

using namespace std;

struct pixel {
  int col;
  int row;
  double ph;
};

struct cluster {
  vector <pixel> vpix;
  int size;
  double sum;
  double col,row;
};

// globals:

pixel pb[155*160]; // pixel block for clustering
int fNHit; // for clustering

// ----------------------------------------------------------------------
bool kbhit()
{
  int byteswaiting;
  ioctl(0, FIONREAD, &byteswaiting);
  return byteswaiting > 0;
}

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
  // set styles:

  gStyle->SetTextFont(62); // 62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align
  gStyle->SetTitleX( 0.15 ); // global title
  gStyle->SetTitleY( 0.98 ); // global title

  gStyle->SetTitleFont( 62, "XYZ" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 1.4, "y" );
  gStyle->SetTitleOffset( 1.7, "z" );

  gStyle->SetLabelFont( 62, "XYZ" );

  gStyle->SetLabelOffset( 0.022, "xyz" );

  gStyle->SetTickLength( -0.02, "xyz" ); // tick marks outside

  gStyle->SetLineWidth(1);// frames
  gStyle->SetHistLineColor(4); // 4=blau
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistFillColor(5); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle(1001); // 1001 = solid

  gStyle->SetFrameLineWidth(2);

  // statistics box:

  gStyle->SetOptStat(10);
  gStyle->SetStatFont(42); // 42 = Helvetica normal
  gStyle->SetStatBorderSize(1); // no 'shadow'
  gStyle->SetStatX(0.82);
  gStyle->SetStatY(0.92);

  //gStyle->SetPalette(55); // sunset
  //gStyle->SetPalette(56); // white to red
  //gStyle->SetPalette(73); // blue
  gStyle->SetPalette(90); // green to magenta
  //gStyle->SetPalette(109); // sea blue to magenta
  gStyle->SetNumberContours(32); // -20..300

  TApplication app( "app", 0, 0 );
  TCanvas c1( "c1", "pixel event display", 900, 800 ); // square
  c1.SetTopMargin( 0.12 );
  c1.SetRightMargin( 0.18 );

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
  bool more = 1;
  string q{"q"};

  while( evFile.good() && ! evFile.eof() && more ) {

    istringstream iss( evseed ); // tokenize string
    int iev;
    iss >> iev;
    cout << "ev " << iev;

    string filled;
    iss >> filled;

    if( filled == F ) {

      string roi;
      getline( evFile, roi );
      istringstream css( roi ); // tokenize string

      int npx = 0;
      vector <pixel> vpx;
      vpx.reserve(35);

      cout << " px";

      while( ! css.eof() ) {

	int col;
	int row;
	double ph;
	css >> col;
	css >> row;
	css >> ph;

	++npx;

	pixel px { col, row, ph };
	vpx.push_back(px);

      } // roi px

      int mpx = 0;
      double sumph = 0;

      TH2I hpxmap( "pxmap",
		   Form( "pixel map %i;col;row;PH [ADC]", iev ),
		   155, -0.5, 154.5, 160, -0.5, 159.5 );
      hpxmap.SetMinimum(-20);
      hpxmap.SetMaximum(300);

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

	hpxmap.Fill( col4, row4, dph );

	//cout << " " << col << " " << row << " " << ph;

	if( dph > 12 ) {
	  pb[mpx].col = col4;
	  pb[mpx].row = row4;
	  pb[mpx].ph = dph;
	  ++mpx;
	  sumph += dph;
	}

      } // roi px

      cout << ", roi size " << npx
	   << ", hits " << mpx
	   << ", sumdph " << sumph
	   << endl;

      // clustering:

      fNHit = mpx; // for cluster search

      vector <cluster> vcl = getClus();	    

      if( vcl.size() )
	cout << "  clusters " << vcl.size() << endl;

      for( unsigned icl = 0; icl < vcl.size(); ++ icl ) {

	cout << "  size " << vcl[icl].size
	     << " at " << vcl[icl].col
	     << ", " << vcl[icl].row
	     << endl;

	//if( mpx > 11 ) { // delta-rays
	//if( sumph > 555 ) { // delta-rays
	//if( vcl[icl].size > 11 ) { // delta-rays
	//if( vcl[icl].size > 4 ) { // delta-rays
	if( vcl[icl].size > 20 ) { // delta-rays
	  //if( vcl[icl].sum > 888 ) { // delta-rays

	  hpxmap.Draw( "colz" );
	  c1.Update();

	  cout << "enter any key, q to stop" << endl;

	  while( !kbhit() )
	    gSystem->ProcessEvents(); // ROOT

	  string any;
	  cin >> any;
	  if( any == q )
	    more = 0;

	} // show

      } // icl

    } // filled

    else
      cout << "  empty" << endl;

    ++nev;

    getline( evFile, evseed ); // read ahead

  } // while events

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cout << "done " << evFileName
       << endl;

  return 0;

}
