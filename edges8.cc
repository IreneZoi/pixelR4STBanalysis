
// Daniel Pitzl (DESY) Nov 2017
// R4S edge-on, combined 8 cols 50x50

// edges8 -f -c 118 B/roi000798.txt  103k
// edges8 -f -c 118 B/roi000802.txt  160k
// edges8 -f -c 118 B/roi000814.txt 2570k

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

  string gain;

  if( chip == 118 )
    gain = "B/r118-scancal-tb21-1105.dat";

  double ke = 0.033; // edge-on Landau peak at 30 ke in 400 um for c118

  ifstream gainFile( gain );

  if( ! gainFile ) {
    cout << "gain file not found" << endl;
    return 1;
  }

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

  TFile * histoFile = new TFile( Form( "edges8-%i.root", run ), "RECREATE" );

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
  TH2I * hpxmap = new TH2I( "pxmap", "pixel map, dph > cut;col;row;pixels above cut",
			    nbx, -0.5, nbx-0.5, nby, -0.5, nby-0.5 );
  TH1I hnpx( "npx", "PH pixels per event;PH pixels;events", 80, 0.5, 80.5 );

  TH1D hncl( "ncl", "cluster per event;cluster;events", 20, 0.5,  20.5 );
  TH2I * hclmap = new TH2I( "clmap", "cluster map;col;row;clusters",
			    155, -0.5, 154.5, 160, -0.5, 159.5 );
  TH1I hclsz( "clsz", "cluster size;cluster size [pixels];clusters", 80, 0.5, 80.5 );
  TH1I hclph( "clph", "cluster PH;cluster ph [ADC];clusters", 200, 0, 1000 );
  TH1I hncol( "ncol", "cluster cols;cluster size [cols];clusters", nbx+5, 0.5, nbx+5.5 );
  TH1I hnrow( "nrow", "cluster rows;cluster size [rows];clusters", nby+5, 0.5, nby+5.5 );
  TH1I hrowq( "rowq", "row charge;row charge [ke];rows", 100, 0, 50 );
  TH1I haspect( "aspect", "rows/cols;cluster aspect ratio [rows/cols];clusters", 150, 0, 1.5 );
  TH1I hcolsz( "colsz", "column size;column size [rows];columns", 11, -0.5, 10.5 );
  TH1I hcolph( "colph", "column PH;column pulse height [ADC];columns", 150, 0, 3000 );
  TH1I hcolq( "colq", "column charge;column charge [ke];columns", 150, 0, 150 );
  TH1I hcolq2( "colq2", "2-row column charge;column charge [ke];2-row columns", 150, 0, 150 );
  TProfile qvscol( "qvscol", "column charge;column;<column charge> [ke]",
		   nbx, -0.5, nbx-0.5, 0, 99 );

  TH1D hdy( "dy", "triplet residual;triplet residual [um];triplets",
	    200, -50, 50 );
  TH1D hdyc( "dyc", "triplet residual;triplet residual [um];triplets",
	    200, -50, 50 );
  TH1D hdyc2( "dyc2", "triplet residual, 1-2-rows;triplet residual [um];1-2-row triplets",
	      200, -50, 50 );

  TProfile nrowvsym( "nrowvsym", "rows vs ymod;y mod 50 [#mum];<rows>",
		     100, 0, 50, 0, 10 );
  TProfile nrowvsq( "nrowvsq", "rows vs q;column charge [ke];<rows>",
		    100, 0, 100, 0, 50 );
  TProfile madyvsn( "madyvsn", "MAD y vs nrow;column size [rows];MAD y [#mum]",
		    20, 0.5, 20.5, 0, 100 );

  TProfile madyvsx( "madyvsx", "MAD y vs x;x [columns];MAD y [#mum]",
		    nbx, -0.5, nbx-0.5, 0, 100 );
  TProfile madyvsy( "madyvsy", "MAD y vs y;y [rows];MAD y [#mum]",
		    nby, -0.5, nby-0.5, 0, 100 );

  TProfile dyvsym( "dyvsym", "dy vs ymod;y mod 50 [#mum];<#Deltay> [#mum]",
		   100, 0, 50, -50, 50 );
  TProfile madyvsym( "madyvsym", "MAD y vs ymod;y mod 50 [#mum];MAD y [#mum]",
		     100, 0, 50, 0, 100 );

  TProfile madyvsq( "madyvsq", "MAD y vs q;column charge [ke];MAD y [#mum]",
		    100, 0, 100, 0, 100 );

  TH1D hdycq2( "dycq2", "triplet residual;triplet residual [um];Q peak triplets",
	       200, -50, 50 );

  TH1D hdycq3( "dycq3", "triplet residual;triplet residual [um];Q peak triplets",
	       200, -50, 50 );
  TH1D hdycq4( "dycq4", "triplet residual;triplet residual [um];Q peak triplets",
	       200, -50, 50 );

  TH1D hslp( "slp", "track slope;track angle [pixels];tracks", 250, -2, 0.5 );
  TProfile slpvsncol( "slpvsncol", "slope vs track length;track length [columns];<slope> [pixels]",
		    nbx, 0.5, nbx+0.5, -2, 1 );

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

	  //if( dph > 12 ) { // 648 dyc 1.90
	  if( dph > 10 ) { // 648 duc 1.84
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

      if( ncol > 13 ) { // 1+4+4+4+1

	double aspect = (double)nrow / (double)ncol;
	haspect.Fill( aspect ); // should be 0.25

	double x0 = colmin+1;
	double x9 = colmax-1;
	double y0 = rowmax; // negative slp
	double y9 = rowmin;

	// 8-col triplets:

	for( int col = colmin+1; col < colmax-23; col += 8 ) { // 0  1+2+3+4+5+6+7+8  9+10+11+12+13+14+15+16  17+18+19+20+21+22+23+24  25

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
	  if( colvec[col+12].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+13].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+14].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+15].size() == 0 ) continue; // gap in cluster ?

	  if( colvec[col+16].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+17].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+18].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+19].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+20].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+21].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+22].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+23].size() == 0 ) continue; // gap in cluster ?

	  double phA = 0;
	  double qA = 0;
	  double qyA = 0;
	  rowmin = 999;
	  rowmax = 0;
	  double rowqA[160] = { 0 };

	  for( unsigned ir = 0; ir < colvec[col+0].size(); ++ir ) {
	    phA += colvec[col+0][ir].ph;
	    double q = colvec[col+0][ir].q;
	    qA += q;
	    int row = colvec[col+0][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyA += q * row;
	    rowqA[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+1].size(); ++ir ) {
	    phA += colvec[col+1][ir].ph;
	    double q = colvec[col+1][ir].q;
	    qA += q;
	    int row = colvec[col+1][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyA += q * row;
	    rowqA[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+2].size(); ++ir ) {
	    phA += colvec[col+2][ir].ph;
	    double q = colvec[col+2][ir].q;
	    qA += q;
	    int row = colvec[col+2][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyA += q * row;
	    rowqA[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+3].size(); ++ir ) {
	    phA += colvec[col+3][ir].ph;
	    double q = colvec[col+3][ir].q;
	    qA += q;
	    int row = colvec[col+3][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyA += q * row;
	    rowqA[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+4].size(); ++ir ) {
	    phA += colvec[col+4][ir].ph;
	    double q = colvec[col+4][ir].q;
	    qA += q;
	    int row = colvec[col+4][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyA += q * row;
	    rowqA[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+5].size(); ++ir ) {
	    phA += colvec[col+5][ir].ph;
	    double q = colvec[col+5][ir].q;
	    qA += q;
	    int row = colvec[col+5][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyA += q * row;
	    rowqA[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+6].size(); ++ir ) {
	    phA += colvec[col+6][ir].ph;
	    double q = colvec[col+6][ir].q;
	    qA += q;
	    int row = colvec[col+6][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyA += q * row;
	    rowqA[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+7].size(); ++ir ) {
	    phA += colvec[col+7][ir].ph;
	    double q = colvec[col+7][ir].q;
	    qA += q;
	    int row = colvec[col+7][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyA += q * row;
	    rowqA[row] += q;
	  }
	  int nrowA = rowmax-rowmin+1;
	  hcolsz.Fill( nrowA );
	  for( int ir = rowmin; ir <= rowmax; ++ir )
	    hrowq.Fill( rowqA[ir] );
	  hcolph.Fill( phA );
	  hcolq.Fill( qA ); // peak at 30 ke
	  if( nrowA <= 2 )
	    hcolq2.Fill( qA ); // peak at 30 ke
	  qvscol.Fill( col+0, qA );
	  qvscol.Fill( col+1, qA );
	  qvscol.Fill( col+2, qA );
	  qvscol.Fill( col+3, qA );
	  qvscol.Fill( col+4, qA );
	  qvscol.Fill( col+5, qA );
	  qvscol.Fill( col+6, qA );
	  qvscol.Fill( col+7, qA );
	  double yA = qyA / qA;
	  if( col == colmin+1 )
	    y0 = yA;

	  // B = 8..15

	  double qB = 0;
	  double qyB = 0;
	  rowmin = 999;
	  rowmax = 0;
	  double rowqB[160] = { 0 };

	  for( unsigned ir = 0; ir < colvec[col+8].size(); ++ir ) {
	    double q = colvec[col+8][ir].q;
	    qB += q;
	    int row = colvec[col+8][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyB += q * row;
	    rowqB[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+9].size(); ++ir ) {
	    double q = colvec[col+9][ir].q;
	    qB += q;
	    int row = colvec[col+9][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyB += q * row;
	    rowqB[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+10].size(); ++ir ) {
	    double q = colvec[col+10][ir].q;
	    qB += q;
	    int row = colvec[col+10][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyB += q * row;
	    rowqB[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+11].size(); ++ir ) {
	    double q = colvec[col+11][ir].q;
	    qB += q;
	    int row = colvec[col+11][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyB += q * row;
	    rowqB[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+12].size(); ++ir ) {
	    double q = colvec[col+12][ir].q;
	    qB += q;
	    int row = colvec[col+12][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyB += q * row;
	    rowqB[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+13].size(); ++ir ) {
	    double q = colvec[col+13][ir].q;
	    qB += q;
	    int row = colvec[col+13][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyB += q * row;
	    rowqB[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+14].size(); ++ir ) {
	    double q = colvec[col+14][ir].q;
	    qB += q;
	    int row = colvec[col+14][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyB += q * row;
	    rowqB[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+15].size(); ++ir ) {
	    double q = colvec[col+15][ir].q;
	    qB += q;
	    int row = colvec[col+15][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyB += q * row;
	    rowqB[row] += q;
	  }
	  int nrowB = rowmax-rowmin+1;
	  hcolsz.Fill( nrowB );
	  for( int ir = rowmin; ir <= rowmax; ++ir )
	    hrowq.Fill( rowqB[ir] );
	  hcolq.Fill( qB );
	  if( nrowB <= 2 )
	    hcolq2.Fill( qB ); // peak at 30 ke
	  qvscol.Fill( col+8, qB );
	  qvscol.Fill( col+9, qB );
	  qvscol.Fill( col+10, qB );
	  qvscol.Fill( col+11, qB );
	  qvscol.Fill( col+12, qB );
	  qvscol.Fill( col+13, qB );
	  qvscol.Fill( col+14, qB );
	  qvscol.Fill( col+15, qB );
	  double yB = qyB / qB;
	  qB = 0;
	  qyB = 0;
	  for( int ir = rowmin; ir <= rowmax; ++ir ) {
	    qB += rowqB[ir];
	    qyB += ir * rowqB[ir];
	  }
	  yB = qyB / qB;
 
	  // C = 16..23

	  double phC = 0;
	  double qC = 0;
	  double qyC = 0;
	  rowmin = 999;
	  rowmax = 0;
	  double rowqC[160] = { 0 };

	  for( unsigned ir = 0; ir < colvec[col+16].size(); ++ir ) {
	    phC += colvec[col+16][ir].ph;
	    double q = colvec[col+16][ir].q;
	    qC += q;
	    int row = colvec[col+16][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyC += q * row;
	    rowqC[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+17].size(); ++ir ) {
	    phC += colvec[col+17][ir].ph;
	    double q = colvec[col+17][ir].q;
	    qC += q;
	    int row = colvec[col+17][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyC += q * row;
	    rowqC[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+18].size(); ++ir ) {
	    phC += colvec[col+18][ir].ph;
	    double q = colvec[col+18][ir].q;
	    qC += q;
	    int row = colvec[col+18][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyC += q * row;
	    rowqC[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+19].size(); ++ir ) {
	    phC += colvec[col+19][ir].ph;
	    double q = colvec[col+19][ir].q;
	    qC += q;
	    int row = colvec[col+19][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyC += q * row;
	    rowqC[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+20].size(); ++ir ) {
	    phC += colvec[col+20][ir].ph;
	    double q = colvec[col+20][ir].q;
	    qC += q;
	    int row = colvec[col+20][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyC += q * row;
	    rowqC[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+21].size(); ++ir ) {
	    phC += colvec[col+21][ir].ph;
	    double q = colvec[col+21][ir].q;
	    qC += q;
	    int row = colvec[col+21][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyC += q * row;
	    rowqC[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+22].size(); ++ir ) {
	    phC += colvec[col+22][ir].ph;
	    double q = colvec[col+22][ir].q;
	    qC += q;
	    int row = colvec[col+22][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyC += q * row;
	    rowqC[row] += q;
	  }
	  for( unsigned ir = 0; ir < colvec[col+23].size(); ++ir ) {
	    phC += colvec[col+23][ir].ph;
	    double q = colvec[col+23][ir].q;
	    qC += q;
	    int row = colvec[col+23][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyC += q * row;
	    rowqC[row] += q;
	  }
	  int nrowC = rowmax-rowmin+1;
	  hcolsz.Fill( nrowC );
	  for( int ir = rowmin; ir <= rowmax; ++ir )
	    hrowq.Fill( rowqC[ir] );
	  hcolph.Fill( phC );
	  hcolq.Fill( qC ); // peak at 31 ke
	  if( nrowC <= 2 )
	    hcolq2.Fill( qC ); // peak at 30 ke
	  qvscol.Fill( col+16, qC );
	  qvscol.Fill( col+17, qC );
	  qvscol.Fill( col+18, qC );
	  qvscol.Fill( col+19, qC );
	  qvscol.Fill( col+20, qC );
	  qvscol.Fill( col+21, qC );
	  qvscol.Fill( col+22, qC );
	  qvscol.Fill( col+23, qC );
	  double yC = qyC / qC;
	  if( col+23 == colmax-1 || col+23 == colmax-2 || col+23 == colmax-3 || col+23 == colmax-4 ||
	      col+23 == colmax-5 || col+23 == colmax-6 || col+23 == colmax-7 || col+23 == colmax-8 ) {
	    x9 = col+23;
	    y9 = yC;
	  }

	  // triplet residual: B vs A+C

	  double ay = 0.5*(yA+yC);
	  double dy = yB - ay;
	  double ymod = fmod( (ay+0.5)*50, 50 );

	  hdy.Fill( dy*50 ); // [um]

	  if( nrowA <= 2 && // suppress delta rays
	      nrowC <= 2 ) {

	    hdyc.Fill( dy*50 ); // [um] 687: 1.37  688: 1.31  689: 1.30

	    nrowvsym.Fill( ymod, nrowB );
	    nrowvsq.Fill( qB, nrowB );
	    madyvsn.Fill( nrowB, abs(dy)*50 );

	    if( nrowB <= 2 ) {

	      hdyc2.Fill( dy*50 ); // [um] 687: 1.33  688: 1.33  689: 1.31

	      madyvsx.Fill( col+10, abs(dy)*50 );
	      madyvsy.Fill( ay, abs(dy)*50 );
	      dyvsym.Fill( ymod, dy*50 );
	      madyvsym.Fill( ymod, abs(dy)*50 );

	    } // 2-row

	    madyvsq.Fill( qB, abs(dy)*50 );

	    if( qA > 22 && qA < 40 &&
		qB > 22 && qB < 40 &&
		qC > 22 && qC < 40 )
	      hdycq2.Fill( dy*50 ); // [um] 687: 1.23  688: 1.25  689: 1.22

	    if( qA > 22 && qA < 36 &&
		qB > 22 && qB < 36 &&
		qC > 22 && qC < 36 )
	      hdycq3.Fill( dy*50 ); // [um] 687: 1.12  688: 1.16  689: 1.13

	    if( qA > 22 && qA < 32 &&
		qB > 22 && qB < 32 &&
		qC > 22 && qC < 32 )
	      hdycq4.Fill( dy*50 ); // [um] 687: 1.02  688: 1.08  689: 1.05

	  } // cuts

	} // col

	double slp = (y9-y0) / (x9-x0);
	hslp.Fill( slp );
	slpvsncol.Fill( ncol, slp ); // flat

      } // long

      //if( vcl.size() ) cout << endl;

    } // cl

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
