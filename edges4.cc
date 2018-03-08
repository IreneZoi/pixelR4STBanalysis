
// Daniel Pitzl (DESY) Oct 2017
// R4S edge-on, combined 4 cols

// edges4 -c 1092 B/roi000689.txt  rot 3.4 = atan(25/400) = 4-col
// edges4 -f -c 118 B/roi000808.txt  50x50
// edges4 -c 124 B/roi001427.txt

#include <cstdlib> // atoi
#include <iostream> // cout
#include <iomanip> // setw
#include <string> // strings
#include <sstream> // stringstream
#include <fstream> // files
#include <vector>
#include <set>
#include <cmath>

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

struct cluster {
  vector <pixel> vpix;
  int size;
  double sum;
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

  double p0[155][160]; // Fermi
  double p1[155][160];
  double p2[155][160];
  double p3[155][160];

  string gain{ "B/r108-scancal-tb24-1024.dat" };

  double ke = 0.036; // edge-on Landau peak at 7.5 ke in 100 um

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

  if( chip == 118 )
    gain = "B/r118-scancal-tb21-1105.dat";

  if( chip == 124 ) {
    gain = "B/c124i-scancal2-tb21-icy-ia150-hold28-2018-02-21.dat";
    ke = 0.035; // default
  }

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

  TFile * histoFile = new TFile( Form( "edges4-%i.root", run ), "RECREATE" );

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

  double leng = 1; // 100 um
  if( fifty )
    leng = 0.5; // 50 um

  TH2I * hpxmap = new TH2I( "pxmap", "pixel map, dph > cut;col;row;pixels above cut",
			    nbx, -0.5, nbx-0.5, nby, -0.5, nby-0.5 );
  TH1I hnpx( "npx", "PH pixels per event;PH pixels;events", 80, 0.5, 80.5 );

  TH1I hncl( "ncl", "cluster per event;cluster;events", 20, 0.5,  20.5 );
  TH2I * hclmap = new TH2I( "clmap", "cluster map;col;row;clusters",
			    155, -0.5, 154.5, 160, -0.5, 159.5 );
  TH1I hclsz( "clsz", "cluster size;cluster size [pixels];clusters", 200, 0.5, 200.5 );
  TH1I hclph( "clph", "cluster PH;cluster ph [ADC];clusters", 200, 0, 1000 );
  TH1I hncol( "ncol", "cluster cols;cluster size [cols];clusters", nbx+5, 0.5, nbx+5.5 );
  TH1I hnrow( "nrow", "cluster rows;cluster size [rows];clusters", nby+5, 0.5, nby+5.5 );
  TH1I hrowq( "rowq", "row charge;row charge [ke];rows", 100, 0, 50 );
  TH1I haspect( "aspect", "rows/cols;cluster aspect ratio [rows/cols];clusters", 150, 0, 1.5 );

  TH2I * hq0q1 = new TH2I( "q0q1", "charge correlation;q0 [ke];q1 [ke];column pairs",
			   100, 0, 25, 100, 0, 25 );
  TH2I * hq0q2 = new TH2I( "q0q2", "charge correlation;q0 [ke];q2 [ke];column pairs",
			   100, 0, 25, 100, 0, 25 );
  TH2I * hq0q9 = new TH2I( "q0q9", "charge correlation;q0 [ke];q9 [ke];column pairs",
			   100, 0, 25, 100, 0, 25 );

  TH1I hq3any( "q3any", "3-column charge;3-column charge [ke];3-columns", 150, 0, 150 );
  TH1I hq3min( "q3min", "min 3-of-4 column charge;min 3-of-4 column charge [ke];3-columns",
	       150, 0, 150 );

  TH1I hcolsz( "colsz", "column size;column size [rows];4-columns", 11, -0.5, 10.5 );
  TH1I hcolph( "colph", "column PH;4-column pulse height [ADC];4-columns", 150, 0, 3000 );
  TH1I hcolq( "colq", "column charge;4-column charge [ke];4-columns", 150, 0, 150 );
  TProfile qvscol( "qvscol", "column charge;column;<4-column charge> [ke]",
		   nbx, -0.5, nbx-0.5, 0, 99 );
  TH1I hncut( "ncut", "size cut cols;columns above size cut;4-clusters", 5, -0.5, 4.5 );

  TH1I hdy( "dy", "4-column triplet residual;triplet residual [#mum];4-column triplets",
	    200, -50, 50 );
  TH1I hdyc( "dyc", "4-column triplet residual;triplet residual [#mum];4-column triplets",
	    200, -50, 50 );
  TH1I hduc( "duc", "4-column triplet Moyal residual;triplet Moyal residual [#mum];4-column triplets",
	    200, -50, 50 );
  TH1I hdvc( "dvc",
	     "4-column triplet truncated residual;triplet truncated residual [#mum];4-column triplets",
	    200, -50, 50 );
  TH1I hdyc2( "dyc2",
	      "4-column triplet residual, 1-2-rows;triplet residual [#mum];1-2-row 4-column triplets",
	      200, -50, 50 );

  TProfile nrowvsym( "nrowvsym", "rows vs ymod;y mod 25 [#mum];<rows>",
		     50, 0, pitch, 0, 10 );
  TProfile nrowvsq( "nrowvsq", "rows vs q;4-column charge [ke];<rows>",
		    100, 0, 100, 0, 50 );
  TProfile qvsym( "qvsym", "charge vs ymod;y mod 25 [#mum];<charge> [ke]",
		     50, 0, pitch, 0, 99 );
  TProfile madyvsn( "madyvsn", "MAD y vs nrow;4-column size [rows];MAD y [#mum]",
		    20, 0.5, 20.5, 0, 100 );

  TProfile madyvsx( "madyvsx", "MAD y vs x;x [columns];MAD y [#mum]",
		    nbx, -0.5, nbx-0.5, 0, 100 );
  TProfile madyvsy( "madyvsy", "MAD y vs y;y [rows];MAD y [#mum]",
		    nby, -0.5, nby-0.5, 0, 100 );

  TProfile dyvsym( "dyvsym", "dy vs ymod;y mod 25 [#mum];<#Deltay> [#mum]",
		   50, 0, pitch, -50, 50 );
  TProfile madyvsym( "madyvsym", "MAD y vs ymod;y mod 25 [#mum];MAD #Deltay [#mum]",
		     50, 0, pitch, 0, 100 );

  TProfile madyvsq( "madyvsq", "MAD y vs q;4-column charge [ke];MAD #Deltay [#mum]",
		    100, 0, 100, 0, 100 );
  TProfile maduvsq( "maduvsq", "MAD Moyal #Deltay vs q;4-column charge [ke];MAD Moyal #Deltay [#mum]",
		    100, 0, 100, 0, 100 );
  TProfile madvvsq( "madvvsq", "MAD truncated #Deltay vs q;4-column charge [ke];MAD truncated #Deltay [#mum]",
		    100, 0, 100, 0, 100 );

  TH1I hdycq2( "dycq2", "triplet residual;triplet residual [#mum];Q peak triplets",
	       200, -50, 50 );

  TH1I hdycq3( "dycq3", "triplet residual;triplet residual [#mum];Q peak triplets",
	       200, -50, 50 );
  TH1I hdycq4( "dycq4", "triplet residual;triplet residual [#mum];Q peak triplets",
	       200, -50, 50 );
  TH1I hdycq22( "dycq22", "triplet residual;triplet residual [#mum];Q < 22 ke",
	       200, -50, 50 );
  TH1I hdycq32( "dycq32", "triplet residual;triplet residual [#mum];22 < Q < 32 ke",
	       200, -50, 50 );
  TH1I hdycq36( "dycq36", "triplet residual;triplet residual [#mum];32 < Q < 36 ke",
	       200, -50, 50 );
  TH1I hdycq40( "dycq40", "triplet residual;triplet residual [#mum];36 < Q < 40 ke",
	       200, -50, 50 );
  TH1I hdycq50( "dycq50", "triplet residual;triplet residual [#mum];40 < Q < 50 ke",
	       200, -50, 50 );
  TH1I hdycq60( "dycq60", "triplet residual;triplet residual [#mum];50 < Q < 60 ke",
	       200, -50, 50 );
  TH1I hdycq80( "dycq80", "triplet residual;triplet residual [#mum];60 < Q < 80 ke",
	       200, -100, 100 );
  TH1I hdycq99( "dycq99", "triplet residual;triplet residual [#mum];80 ke < Q",
	       200, -100, 100 );

  TH1I hducq32( "ducq32", "triplet truncated residual;triplet truncated residual [#mum];22 < Q < 32 ke",
	       200, -50, 50 );
  TH1I hducq36( "ducq36", "triplet truncated residual;triplet truncated residual [#mum];32 < Q < 36 ke",
	       200, -50, 50 );
  TH1I hducq40( "ducq40", "triplet truncated residual;triplet truncated residual [#mum];36 < Q < 40 ke",
	       200, -50, 50 );
  TH1I hducq50( "ducq50", "triplet truncated residual;triplet truncated residual [#mum];40 < Q < 50 ke",
	       200, -50, 50 );
  TH1I hducq60( "ducq60", "triplet truncated residual;triplet truncated residual [#mum];50 < Q < 60 ke",
	       200, -50, 50 );
  TH1I hducq80( "ducq80", "triplet truncated residual;triplet truncated residual [#mum];60 < Q < 80 ke",
	       200, -100, 100 );
  TH1I hducq99( "ducq99", "triplet truncated residual;triplet truncated residual [#mum];80 ke < Q",
	       200, -100, 100 );

  TH1D heta( "eta", "charge sharing;eta;2-row clusters",
	     200, -1, 1 );
  TProfile etavsym( "etavsym", "eta vs ymod;y mod 25 [#mum];<eta>",
		     50, 0, pitch, -1, 1 );
  TProfile madyvseta( "madyvseta", "mad y vs eta;eta;MAD #Deltay [#mum]",
		      100, -1, 1, 0, 100 );
  TH1D hdyc2eta( "dyc2eta", "triplet residual;triplet residual [#mum];high eta triplets",
	      200, -50, 50 );
  TH1D hdyc2eta3( "dyc2eta3", "triplet residual;triplet residual [#mum];high eta triplets",
	      200, -50, 50 );

  TH1I hslp( "slp", "track slope;track angle [pixels];tracks", 400, -2, 2 );
  TProfile slpvsncol( "slpvsncol", "slope vs track length;track length [columns];<slope> [pixels]",
		    nbx, 0.5, nbx+0.5, -2, 2 );

  TH1I hq[156];
  for( int lng = 1; lng <= 155; ++lng )
    hq[lng] = TH1I( Form( "q%i", lng ),
		    Form( "q/l for l = %i;charge/length [ke/columns];clusters", lng ),
		    100, 0, 50 );

  //const double qwid = 1.0; // Moyal width in 100 um

  double qmn = 22; // cut
  double qmx = 40; // cut
  if( chip == 124 ) {
    qmn =  6;
    qmx = 20;
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
      istringstream css( roi ); // tokenize string

      int npx = 0;
      vector <pixel> vpx;
      vpx.reserve(35);

      while( ! css.eof() ) { // one line = one event

	int col;
	int row;
	double ph;
	css >> col;
	css >> row;
	css >> ph;

	pixel px { col, row, ph, ph };
	vpx.push_back(px);

	++npx;

      } // roi px

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // column-wise common mode correction:

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
    // clustering:

    hnpx.Fill( pb.size() );

    vector <cluster> vcl = getClus(pb);

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
      vector < pixel > colvec[155]; // 155 col max, array of vectors of pixels
      double qvec[155]{0};

      for ( int ipx = 0; ipx < vcl[icl].size; ++ipx ) {

	int col = vcl[icl].vpix[ipx].col;
	if( col < colmin ) colmin = col;
	if( col > colmax ) colmax = col;

	int row = vcl[icl].vpix[ipx].row;
	if( row < rowmin ) rowmin = row;
	if( row > rowmax ) rowmax = row;

	colvec[col].push_back( vcl[icl].vpix[ipx] ); // add px to col

	qvec[col] += vcl[icl].vpix[ipx].q;

      } // px

      int ncol = colmax - colmin + 1;
      hncol.Fill( ncol );

      int nrow = rowmax - rowmin + 1;
      hnrow.Fill( nrow );

      if( ncol > 13 ) { // 1+4+4+4+1

	double aspect = (double)nrow / (double)ncol;
	haspect.Fill( aspect ); // should be 0.25

	double x0 = colmin+2;
	double x9 = colmax-2;
	double y0 = rowmax; // negative slp
	double y9 = rowmin;

	// 4-col triplets:

	for( int col = colmin+1; col < colmax-11; col += 4 ) { // 0 1+2+3+4 5+6+7+8 9+10+11+12 13

	  if( colvec[col+0].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+1].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+2].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+3].size() == 0 ) continue; // gap in cluster ?

	  if( colvec[col+4].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+5].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+6].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+7].size() == 0 ) continue; // gap in cluster ?

	  if( colvec[col+8].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+9].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+10].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+11].size() == 0 ) continue; // gap in cluster ?

	  hq0q1->Fill( qvec[col+0], qvec[col+1] );
	  hq0q2->Fill( qvec[col+0], qvec[col+2] );
	  hq0q9->Fill( qvec[col+0], qvec[col+9] );

	  unsigned imaxA = 0;
	  double qmaxA = qvec[col+0];
	  double sumq3 = 0;
	  for( unsigned i = 1; i < 4; ++ i ) {
	    sumq3 += qvec[col+i];
	    if( qvec[col+i] > qmaxA ) {
	      qmaxA = qvec[col+i];
	      imaxA = i;
	    }
	  }
	  hq3any.Fill( sumq3 );

	  double sumq3min = 0;
	  for( unsigned i = 0; i < 4; ++ i ) {
	    if( i == imaxA ) continue;
	    sumq3min += qvec[col+i];
	  }
	  hq3min.Fill( sumq3min );

	  double qA = 0;
	  double qyA = 0;
	  rowmin = 999;
	  rowmax = 0;
	  double rowqA[320] = { 0 };

	  for( unsigned ic = 0; ic <= 3; ++ic ) {

	    for( unsigned ir = 0; ir < colvec[col+ic].size(); ++ir ) {
	      int row = colvec[col+ic][ir].row;
	      if( row < rowmin ) rowmin = row;
	      if( row > rowmax ) rowmax = row;
	      double q = colvec[col+ic][ir].q;
	      qA += q;
	      qyA += q * row;
	      rowqA[row] += q;
	    }

	  } // ic

	  int nrowA = rowmax-rowmin+1;
	  for( int ir = rowmin; ir <= rowmax; ++ir )
	    hrowq.Fill( rowqA[ir] );
	  qvscol.Fill( col+0, qA );
	  qvscol.Fill( col+1, qA );
	  qvscol.Fill( col+2, qA );
	  qvscol.Fill( col+3, qA );
	  double yA = qyA / qA;
	  if( col == colmin+1 )
	    y0 = yA;

	  double etaA = 0;
	  if( nrowA == 2 ) {
	    double a0 = rowqA[rowmin];
	    double a1 = rowqA[rowmax];
	    etaA = (a1-a0)/(a1+a0);
	  }

	  // B = 4+5+6+7

	  double phB = 0;
	  double qB = 0;
	  double qyB = 0;
	  double qxB = 0;
	  double qxyB = 0;
	  double qtB = 0;
	  double qtyB = 0;
	  int ncut = 0;
	  rowmin = 999;
	  rowmax = 0;
	  double rowqB[320] = { 0 };

	  for( unsigned ic = 4; ic <= 7; ++ic ) {

	    for( unsigned ir = 0; ir < colvec[col+ic].size(); ++ir ) {

	      int row = colvec[col+ic][ir].row;
	      if( row < rowmin ) rowmin = row;
	      if( row > rowmax ) rowmax = row;
	      phB += colvec[col+ic][ir].ph;
	      double q = colvec[col+ic][ir].q;
	      rowqB[row] += q;

	      qB += q;
	      qyB += q * row;

	      if( colvec[col+ic].size() < 3 ) {
		qxB += q;
		qxyB += q * row;
	      }

	      if( qvec[col+ic] < 16 ) {
		qtB += q;
		qtyB += q * row;
	      }

	    } // ir

	    if( colvec[col+ic].size() > 2 ) ++ncut;

	  } // ic

	  int nrowB = rowmax-rowmin+1;
	  hcolsz.Fill( nrowB );
	  for( int ir = rowmin; ir <= rowmax; ++ir )
	    hrowq.Fill( rowqB[ir] );
	  hcolph.Fill( phB );
	  hcolq.Fill( qB );
	  qvscol.Fill( col+4, qB );
	  qvscol.Fill( col+5, qB );
	  qvscol.Fill( col+6, qB );
	  qvscol.Fill( col+7, qB );
	  hncut.Fill( ncut );
	  double yB = qyB / qB; // 688 dyc 1.32
	  double uB = qxyB / qxB; // truncated
	  double vB = qtyB / qtB; // truncated

	  qB = 0;
	  qyB = 0;
	  for( int row = rowmin; row <= rowmax; ++row ) {
	    // simulate threshold: colq peak at 31 ke
	    /*
	    if( rowqB[row] < 15 ) continue; // 688 dyc2 RMS 6.46
	    if( rowqB[row] < 14 ) continue; // 688 dyc2 RMS 6.08
	    if( rowqB[row] < 13 ) continue; // 688 dyc2 RMS 5.63
	    if( rowqB[row] < 12 ) continue; // 688 dyc2 RMS 5.10
	    if( rowqB[row] < 11 ) continue; // 688 dyc2 RMS 4.57
	    if( rowqB[row] < 10 ) continue; // 688 dyc2 RMS 4.06
	    if( rowqB[row] <  9 ) continue; // 688 dyc2 RMS 3.59
	    if( rowqB[row] <  8 ) continue; // 688 dyc2 RMS 3.16
	    if( rowqB[row] <  7 ) continue; // 688 dyc2 RMS 2.77
	    if( rowqB[row] <  6 ) continue; // 688 dyc2 RMS 2.46
	    if( rowqB[row] <  5 ) continue; // 688 dyc2 RMS 2.22
	    if( rowqB[row] <  4 ) continue; // 688 dyc2 RMS 2.03
	    if( rowqB[row] <  3 ) continue; // 688 dyc2 RMS 1.88
	    if( rowqB[row] < 2.5 ) continue; // 688 dyc2 RMS 1.83
	    if( rowqB[row] <  2 ) continue; // 688 dyc2 RMS 1.79
	    if( rowqB[row] < 1.5 ) continue; // 688 dyc2 RMS 1.77
	    if( rowqB[row] < 1.2 ) continue; // 688 dyc2 RMS 1.76
	    if( rowqB[row] <  1 ) continue; // 688 dyc2 RMS 1.76
	    */
	    qB += rowqB[row];
	    qyB += row * rowqB[row];
	  }
	  yB = qyB / qB; // dyc2 1.33

	  double etaB = 2;
	  if( nrowB == 2 ) {
	    double a0 = rowqB[rowmin];
	    double a1 = rowqB[rowmax];
	    etaB = (a1-a0)/(a1+a0);
	  }

	  // C = 8+9+10+11

	  double qC = 0;
	  double qyC = 0;
	  rowmin = 999;
	  rowmax = 0;
	  double rowqC[320] = { 0 };

	  for( unsigned ic = 8; ic <= 11; ++ic ) {

	    for( unsigned ir = 0; ir < colvec[col+ic].size(); ++ir ) {
	      int row = colvec[col+ic][ir].row;
	      if( row < rowmin ) rowmin = row;
	      if( row > rowmax ) rowmax = row;
	      double q = colvec[col+ic][ir].q;
	      qC += q;
	      qyC += q * row;
	      rowqC[row] += q;
	    }

	  } // ic

	  int nrowC = rowmax-rowmin+1;
	  for( int ir = rowmin; ir <= rowmax; ++ir )
	    hrowq.Fill( rowqC[ir] );
	  qvscol.Fill( col+8, qC );
	  qvscol.Fill( col+9, qC );
	  qvscol.Fill( col+10, qC );
	  qvscol.Fill( col+11, qC );
	  double yC = qyC / qC;
	  if( col+11 == colmax-1 || col+11 == colmax-2 || col+11 == colmax-3 || col+11 == colmax-4 ) {
	    x9 = col+10;
	    y9 = yC;
	  }

	  double etaC = 0;
	  if( nrowC == 2 ) {
	    double a0 = rowqC[rowmin];
	    double a1 = rowqC[rowmax];
	    etaC = (a1-a0)/(a1+a0);
	  }

	  // triplet residual: B vs A+C

	  double ay = 0.5*(yA+yC);
	  double dy = yB - ay;
	  double du = uB - ay;
	  double dv = vB - ay;
	  double ymod = fmod( (ay+0.5)*pitch, pitch );

	  hdy.Fill( dy*pitch ); // [um]

	  if( nrowA <= 2 && // suppress delta rays
	      nrowC <= 2 ) {

	    hdyc.Fill( dy*pitch ); // [um] 687: 1.37  688: 1.31  689: 1.30
	    hduc.Fill( du*pitch ); // [um] 687: 1.37  688: 1.31  689: 1.30
	    hdvc.Fill( dv*pitch ); // [um] 687: 1.37  688: 1.31  689: 1.30

	    nrowvsym.Fill( ymod, nrowB );
	    nrowvsq.Fill( qB, nrowB );
	    qvsym.Fill( ymod, qB ); // flat
	    madyvsn.Fill( nrowB, abs(dy)*pitch );

	    if( nrowB <= 2 ) {

	      hdyc2.Fill( dy*pitch ); // [um] 687: 1.33  688: 1.33  689: 1.31

	      madyvsx.Fill( col+6, abs(dy)*pitch );
	      madyvsy.Fill( ay, abs(dy)*pitch );
	      dyvsym.Fill( ymod, dy*pitch );
	      madyvsym.Fill( ymod, abs(dy)*pitch );

	    } // 2-row B

	    madyvsq.Fill( qB, abs(dy)*pitch );
	    maduvsq.Fill( qB, abs(du)*pitch );
	    madvvsq.Fill( qB, abs(dv)*pitch );

	    if( qA > qmn*leng && qA < qmx*leng &&
		qB > qmn*leng && qB < qmx*leng &&
		qC > qmn*leng && qC < qmx*leng )
	      hdycq2.Fill( dy*pitch ); // [um] 687: 1.23  688: 1.25  689: 1.15

	    if( qA > qmn*leng && qA < 36*leng &&
		qB > qmn*leng && qB < 36*leng &&
		qC > qmn*leng && qC < 36*leng )
	      hdycq3.Fill( dy*pitch ); // [um] 687: 1.12  688: 1.16  689: 1.08

	    if( qA > qmn*leng && qA < 32*leng && // peak at 32 ke
		qB > qmn*leng && qB < 32*leng &&
		qC > qmn*leng && qC < 32*leng )
	      hdycq4.Fill( dy*pitch ); // [um] 687: 1.02  688: 1.08  689: 0.99

	    if(      qB*leng < 22 )
	      hdycq22.Fill( dy*pitch );
	    else if( qB*leng < 32 ) {
	      hdycq32.Fill( dy*pitch );
	      hducq32.Fill( du*pitch );
	    }
	    else if( qB*leng < 36 ) {
	      hdycq36.Fill( dy*pitch );
	      hducq36.Fill( du*pitch );
	    }
	    else if( qB*leng < 40 ) {
	      hdycq40.Fill( dy*pitch );
	      hducq40.Fill( du*pitch );
	    }
	    else if( qB*leng < 50 ) {
	      hdycq50.Fill( dy*pitch );
	      hducq50.Fill( du*pitch );
	    }
	    else if( qB*leng < 60 ) {
	      hdycq60.Fill( dy*pitch );
	      hducq60.Fill( du*pitch );
	    }
	    else if( qB*leng < 80 ) {
	      hdycq80.Fill( dy*pitch );
	      hducq80.Fill( du*pitch );
	    }
	    else {
	      hdycq99.Fill( dy*pitch );
	      hducq99.Fill( du*pitch );
	    }

	    if( nrowB == 2 ) {

	      heta.Fill( etaB );
	      etavsym.Fill( ymod, etaB );
	      madyvseta.Fill( etaB, abs(dy)*pitch ); // 

	      if( abs(etaB) > 0.6 )
		hdyc2eta.Fill( dy*pitch ); // [um] 1.05 um

	      if( abs(etaA) > 0.6 && abs(etaB) > 0.6 && abs(etaC) > 0.6 )
		hdyc2eta3.Fill( dy*pitch ); // [um] 0.90 um

	    } // 2-rows

	  } // cuts A,C

	} // col

	double slp = (y9-y0) / (x9-x0);
	hslp.Fill( slp );
	slpvsncol.Fill( ncol, slp ); // flat

	// Landau vs length:

	int mxl = colmax-colmin-1; // skip 1st and lst, > 11

	for( int lng = 1; lng <= mxl; ++lng ) {

	  for( int col1 = colmin+1; col1 <= colmax-lng; col1 += lng ) { // step along track

	    double sumq = 0;

	    for( int col = col1; col < col1+lng; ++col ) { // columns in lng

	      for( unsigned ir = 0; ir < colvec[col].size(); ++ir )
		sumq += colvec[col][ir].q;

	    } // col in lng

	    hq[lng].Fill( sumq/lng );

	  } // col1

	} // lng

      } // long

      //if( vcl.size() ) cout << endl;

    } // cl

    ++nev;

    getline( evFile, evseed ); // read ahead

  } // while events

  cout << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cout << "done " << evFileName
       << endl << "events " << nev
       << endl << histoFile->GetName()
       << endl;

  histoFile->Write();
  histoFile->Close();

  return 0;
}
