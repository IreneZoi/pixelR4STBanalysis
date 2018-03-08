
// Daniel Pitzl (DESY) Oct 2017
// R4S edge-on, combined 2 cols

// edges2 -c 109 B/roi000681.txt  rot 7 = atan(25/200) = 2-col
// edges2 -c 124 B/roi001386.txt  rot 7 = atan(25/200) = 2-col
// edges2 -c 124 B/roi001430.txt

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

  double ke = 0.036; // edge-on Landau peak at 7.5 ke in 100 um

  double p0[155][160]; // Fermi
  double p1[155][160];
  double p2[155][160];
  double p3[155][160];

  string gain{ "B/r108-scancal-tb24-1024.dat" };

  if( chip == 109 )
    gain = "B/r109-scancal-tb24-1025.dat";

  if( chip == 1092 )
    gain = "B/r109-scancal-tb24-1027.dat";

  if( chip == 110 )
    //gain = "B/r110-scancal-tb21-0928-hold24.dat";
    gain = "B/r110-scancal-tb24-1025.dat";

  if( chip == 114 )
    gain = "B/r114-scancal-tb24-1025.dat";

  if( chip == 117 )
    gain = "B/r117-scancal-tb24-1027.dat";

  if( chip == 124 ) {
    gain = "B/c124i-scancal2-tb21-icy-ia150-hold28-2018-02-21.dat";
    ke = 0.035; // default
  }

  ifstream gainFile( gain );

  if( ! gainFile ) {
    cout << "gain file not found" << endl;
    return 1;
  }

  cout << "using gain " << gain << endl;

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

  TFile * histoFile = new TFile( Form( "edges2-%i.root", run ), "RECREATE" );

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

  TH1D hncl( "ncl", "cluster per event;cluster;events", 20, 0.5,  20.5 );
  TH2I * hclmap = new TH2I( "clmap", "cluster map;col;row;clusters",
			    155, -0.5, 154.5, 160, -0.5, 159.5 );
  TH1I hclsz( "clsz", "cluster size;cluster size [pixels];clusters", 80, 0.5, 80.5 );
  TH1I hclph( "clph", "cluster PH;cluster ph [ADC];clusters", 200, 0, 1000 );
  TH1I hncol( "ncol", "cluster cols;cluster size [cols];clusters", nbx+5, 0.5, nbx+5.5 );
  TH1I hnrow( "nrow", "cluster rows;cluster size [rows];clusters", nby+5, 0.5, nby+5.5 );
  TH1I haspect( "aspect", "rows/cols;cluster aspect ratio [rows/cols];clusters", 150, 0, 1.5 );
  TH1I hcolsz( "colsz", "column size;column size [rows];columns", 11, -0.5, 10.5 );
  TH1I hcolq( "colq", "2-column charge;2-column charge [ke];2-columns", 160, 0, 80 );
  TProfile qvscol( "qvscol", "2-column charge;column;<2-column charge> [ke]",
		   nbx, -0.5, nbx-0.5, 0, 99 );
  TH1D hy0( "hy0", "first y;y0 [pixels];tracks", nby, -0.5, nby-0.5 );
  TH1D hy9( "hy9", "last y;y9 [pixels];tracks", nby, -0.5, nby-0.5 );

  TH1D hdy( "dy", "triplet residual;triplet residual [um];triplets",
	    200, -50, 50 );
  TH1D hdyc( "dyc", "triplet residual;triplet residual [um];triplets",
	    200, -50, 50 );
  TH1D hdyc2( "dyc2", "triplet residual, 1-2-rows;triplet residual [um];1-2-row triplets",
	      200, -50, 50 );

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

  TH1D hdycq2( "dycq2", "triplet residual;triplet residual [um];Q peak triplets",
	       200, -50, 50 );

  TH1D hdycq3( "dycq3", "triplet residual;triplet residual [um];Q peak triplets",
	       200, -50, 50 );
  TH1D hdycq4( "dycq4", "triplet residual;triplet residual [um];Q peak triplets",
	       200, -50, 50 );

  TH1D hslp( "slp", "track slope;track angle [pixels];tracks", 400, -2, 2 );
  TProfile slpvsncol( "slpvsncol", "slope vs track length;track length [columns];<slope> [pixels]",
		    nbx, 0.5, nbx+0.5, -2, 2 );

  double qmn = 10; // cut
  double qmx = 22; // cut
  if( chip == 124 ) {
    qmn =  3;
    qmx = 10;
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

      while( ! css.eof() ) {

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

      if( ncol > 7 ) {

	double aspect = (double)nrow / (double)ncol;
	haspect.Fill( aspect ); // should be 0.5

	double x0 = colmin+1;
	double x9 = colmax-1;
	double y0 = rowmax; // negative slp
	double y9 = rowmin;

	//cout << "  cols " << ncol << ":";

	// 2-col triplets:

	for( int col = colmin+1; col < colmax-5; col += 2 ) { // 0 1+2 3+4 5+6 7

	  if( colvec[col+0].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+1].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+2].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+3].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+4].size() == 0 ) continue; // gap in cluster ?
	  if( colvec[col+5].size() == 0 ) continue; // gap in cluster ?

	  //cout << " " << col;

	  double qA = 0;
	  double qyA = 0;
	  rowmin = 999;
	  rowmax = 0;

	  for( unsigned ir = 0; ir < colvec[col+0].size(); ++ir ) {
	    double q = colvec[col+0][ir].q;
	    qA += q;
	    int row = colvec[col+0][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyA += q * row;
	  }
	  for( unsigned ir = 0; ir < colvec[col+1].size(); ++ir ) {
	    double q = colvec[col+1][ir].q;
	    qA += q;
	    int row = colvec[col+1][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyA += q * row;
	  }
	  int nrowA = rowmax-rowmin+1;
	  qvscol.Fill( col+0, qA );
	  qvscol.Fill( col+1, qA );
	  double yA = qyA / qA;
	  if( col == colmin+1 ) {
	    y0 = yA;
	    hy0.Fill( y0 );
	  }

	  // B = 2+3

	  double qB = 0;
	  double qyB = 0;
	  rowmin = 999;
	  rowmax = 0;

	  for( unsigned ir = 0; ir < colvec[col+2].size(); ++ir ) {
	    double q = colvec[col+2][ir].q;
	    qB += q;
	    int row = colvec[col+2][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyB += q * row;
	  }
	  for( unsigned ir = 0; ir < colvec[col+3].size(); ++ir ) {
	    double q = colvec[col+3][ir].q;
	    qB += q;
	    int row = colvec[col+3][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyB += q * row;
	  }
	  int nrowB = rowmax-rowmin+1;
	  hcolsz.Fill( nrowB );
	  hcolq.Fill( qB );
	  qvscol.Fill( col+2, qB );
	  qvscol.Fill( col+3, qB );
	  double yB = qyB / qB;

	  // C = 4+5

	  double qC = 0;
	  double qyC = 0;
	  rowmin = 999;
	  rowmax = 0;

	  for( unsigned ir = 0; ir < colvec[col+4].size(); ++ir ) {
	    double q = colvec[col+4][ir].q;
	    qC += q;
	    int row = colvec[col+4][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyC += q * row;
	  }
	  for( unsigned ir = 0; ir < colvec[col+5].size(); ++ir ) {
	    double q = colvec[col+5][ir].q;
	    qC += q;
	    int row = colvec[col+5][ir].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    qyC += q * row;
	  }
	  int nrowC = rowmax-rowmin+1;
	  qvscol.Fill( col+4, qC );
	  qvscol.Fill( col+5, qC );
	  double yC = qyC / qC;
	  if( col+5 == colmax-1 || col+5 == colmax-2 ) {
	    x9 = col+5;
	    y9 = yC;
	    hy9.Fill( y9 );
	  }

	  // triplet residual:

	  double ay = 0.5*(yA+yC);
	  double dy = yB - ay;
	  double ymod = fmod( (ay+0.5)*pitch, 50 );

	  hdy.Fill( dy*pitch ); // [um]

	  if( nrowA <= 2 &&
	      nrowC <= 2 ) {

	    hdyc.Fill( dy*pitch ); // [um] 1.64

	    nrowvsym.Fill( ymod, nrowB );
	    nrowvsq.Fill( qB, nrowB );

	    if( nrowB <= 2 ) {

	      hdyc2.Fill( dy*pitch ); // [um] 1.60

	      dyvsym.Fill( ymod, dy*pitch );
	      madyvsx.Fill( col+3, abs(dy)*pitch );
	      madyvsy.Fill( ay, abs(dy)*pitch );
	      madyvsym.Fill( ymod, abs(dy)*pitch );

	      madyvsq.Fill( qB, abs(dy)*pitch );

	      if( qA > qmn*leng && qA < qmx*leng &&
		  qB > qmn*leng && qB < qmx*leng &&
		  qC > qmn*leng && qC < qmx*leng )
		hdycq2.Fill( dy*pitch ); // [um] 1.46

	      if( qA > 11*leng && qA < 18*leng &&
		  qB > 11*leng && qB < 18*leng &&
		  qC > 11*leng && qC < 18*leng )
		hdycq3.Fill( dy*pitch ); // [um] 1.35

	      if( qA > 12*leng && qA < 16*leng &&
		  qB > 12*leng && qB < 16*leng &&
		  qC > 12*leng && qC < 16*leng )
		hdycq4.Fill( dy*pitch ); // [um] 1.27

	    } // cuts

	  } // cuts

	} // loop cols

	//cout << " max " << colmax << endl;

	double slp = (y9-y0) / (x9-x0);
	hslp.Fill( slp ); // peak at -1
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
