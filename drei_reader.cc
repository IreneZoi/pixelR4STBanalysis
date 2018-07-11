
// Daniel Pitzl (DESY) Sep 2017
// Jan 2018: version for openMP (faster)
// 3 x R4S, 25x100 or 50x50 on rot90 PCB

// drei 338
// drei -l 99999 866
// drei -p 1.2 876
// drei -p 6 885
// drei -p 4.8 998
// drei -p 6.0 -l 81000 1010
// drei -p 6.0 -l 546000 1019
// drei -p 6.0 1020  (2.1M)
// drei -p 5.6 -f 1024  (50x50)

// drei 1757
// drei -f 1842  ************* -f is the option for 50x50 **********************
// drei -f 1894

#include <cstdlib> // atoi
#include <iostream> // cout
#include <iomanip> // setw
#include <string> // strings
#include <sstream> // stringstream
#include <fstream> // files
#include <vector>
#include <list>
#include <cmath>
#include <time.h> // clock_gettime
#include <sched.h> // getcpu
#include <sys/resource.h>
#include <unistd.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>


using namespace std;

struct evInfo {
  uint64_t evtime;
  bool skip;
  string filled;
};

struct pixel {
  int col;
  int row;
  double ph;
  double q;
};

struct cluster {
  vector <pixel> vpix;
  int size; // [px]
  double sum; // [ADC]
  double q; // [ke]
  double col, row; // [px]
  bool iso;
};

const int A{0};
const int B{1};
const int C{2};
string PN[]{"A","B","C"};

double p0[3][155][160]; // Fermi
double p1[3][155][160];
double p2[3][155][160];
double p3[3][155][160];
double ke[3];

TProfile phvsprev[3];
TProfile dphvsprev[3];
TH1I hph[3];
TH1I hdph[3];
TH1I hnpx[3];
TH1I hnht[3];
TH2I * hpxmap[3];
TH1I hncl[3];
TH2I * hclmap[3];
TH1I hclsz[3];
TH1I hclph[3];
TH1I hclq[3];

list < evInfo > infoA;
list < evInfo > infoB;
list < evInfo > infoC;

#define halfSensorX 4.0
#define halfSensorY 3.9
#define ACspacing 40 //[mm]


double nSigmaTolerance = 3;
double beamDivergence = 0.001; //1 mrad
double straightTracks = ACspacing*beamDivergence*nSigmaTolerance;


/*
//------------------------------------------------------------------------------
unsigned digit_value( char c )
{
return unsigned( c - '0' ); // negatives get flipped to large positives
}

//------------------------------------------------------------------------------
int fast_atoi( const char * p )
{
bool neg = false;
if( *p == '-') {
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
*/
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

    // start a new cluster:

    cluster c;
    c.vpix.push_back( pb[seed] );
    gone[seed] = 1;

    // let it grow as much as possible:

    int growing;
    do{
      growing = 0;
      for( unsigned i = 0; i < pb.size(); ++i ) {
        if( !gone[i] ) { // unused pixel
          for( unsigned int p = 0; p < c.vpix.size(); ++p ) { // vpix in cluster so far
            int dr = c.vpix.at(p).row - pb[i].row;
            int dc = c.vpix.at(p).col - pb[i].col;
            if( ( dr >= -fCluCut ) && ( dr <= fCluCut ) &&
		( dc >= -fCluCut ) && ( dc <= fCluCut ) ) {
              c.vpix.push_back(pb[i]);
	      gone[i] = 1;
              growing = 1;
              break; // p, important!
            }
          } // p loop over vpix
        } // not gone
      } // i loop over all pix
    }
    while( growing );

    // added all I could. determine position and append it to the list of clusters:

    c.sum = 0;
    c.q = 0;
    c.size = 0;
    c.col = 0;
    c.row = 0;
    c.iso = 1;

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
list < vector < cluster > > oneplane( int plane, string runnum, unsigned Nev, bool fifty )
{
  int run = stoi( runnum );

  list < vector < cluster > > evlist;
  string datapath = "/mnt/pixeldata/";
  string Xfile;
  if( plane == A ) {
    Xfile = datapath+"a/roi000" + runnum + ".txt";
    if( run > 999 )
      Xfile = datapath+"a/roi00" + runnum + ".txt";
  }
  if( plane == B ) {
    Xfile = datapath+"b/roi000" + runnum + ".txt";
    if( run > 999 )
      Xfile = datapath+"b/roi00" + runnum + ".txt";
  }
  if( plane == C ) {
    Xfile = datapath+"c/roi000" + runnum + ".txt";
    if( run > 999 )
      Xfile = datapath+"c/roi00" + runnum + ".txt";
  }
  // cout << "FOK " << access(Xfile.c_str(),F_OK) << endl;
  // cout << "ROK " << access(Xfile.c_str(),R_OK) << endl;
  cout << "try to open  " << Xfile;
  ifstream Xstream( Xfile.c_str(), ifstream::in );
  cout << " is opening good? " << Xstream.good() << endl;
  if( !Xstream ) {
    cout << " : failed " << endl;
    return evlist;
  }
  cout << " : succeed " << endl;

  string START {"START"};
  string hd; // header

  while( hd != START ) {
    getline( Xstream, hd ); // read one line into string
    cout << "  " << hd << endl;
  }

  string F {"F"}; // filled flag
  string E {"E"}; // empty flag
  string BLANK{" "};

  bool ldb = 0;
  uint64_t prevtime = 0;

  while( Xstream.good() && ! Xstream.eof() &&
	 evlist.size() < Nev ) {

    string evseed;
    getline( Xstream, evseed ); // read one line into string

    istringstream iss( evseed ); // tokenize string
    int iev;
    iss >> iev;

    if( iev%1000 == 0 )
      cout << "  " << PN[plane] << " " << iev << flush;

    string filled;
    iss >> filled;

    int iblk; // event block number: 100, 200, 300...
    iss >> iblk;

    uint64_t evtime = prevtime;

    if( filled == F )
      iss >> evtime; // from run 456

    else if( filled == E )
      iss >> evtime; // from run 456

    vector <pixel> pb; // for clustering
    vector <cluster> vcl;
    evInfo evinf;
    evinf.evtime = evtime;
    evinf.filled = filled;
    evinf.skip = 0;
    prevtime = evtime;

    if( filled == F ) {

      string roi;
      getline( Xstream, roi );

      size_t start = 0;
      size_t gap = 0;
      unsigned ng = 0; // 3 per pix
      string BLANK{" "};
      while( gap < roi.size()-1 ) { // data have trailing blank
	gap = roi.find( BLANK, start );
	start = gap + BLANK.size();
	++ng;
      }
      hnpx[plane].Fill( ng/3 );

      vector <pixel> vpx;

      if( ng/3 < 400 ) {

	vpx.reserve(ng/3);

	size_t start = 0;
	size_t gap = 0;
	while( gap < roi.size()-1 ) { // data have trailing blank

	  gap = roi.find( BLANK, start );
	  string s1( roi.substr( start, gap - start ) );
	  //cout << " " << s1 << "(" << gap << ")";
	  //int col = stoi(s1);
	  int col = atoi( s1.c_str() ); // 4% faster
	  //int col = fast_atoi( s1.c_str() ); // another 6% faster
	  start = gap + BLANK.size();

	  gap = roi.find( BLANK, start );
	  string s2( roi.substr( start, gap - start ) );
	  //cout << " " << s2 << "(" << gap << ")";
	  //int row = stoi(s2);
	  int row = atoi( s2.c_str() );
	  //int row = fast_atoi( s2.c_str() );
	  start = gap + BLANK.size();

	  gap = roi.find( BLANK, start );
	  string s3( roi.substr( start, gap - start ) );
	  //cout << " " << s1 << "(" << gap << ")";
	  //double ph = stod(s3);
	  double ph = atof(s3.c_str());
	  hph[plane].Fill( ph );
	  start = gap + BLANK.size();

	  pixel px { col, row, ph, ph };
	  vpx.push_back(px); // comment out = no clustering

	}
	/* slower
	  istringstream css( roi ); // tokenize string
	  while( ! css.eof() ) {
	  int col;
	  int row;
	  double ph;
	  css >> col;
	  css >> row;
	  css >> ph;

	  pixel px { col, row, ph, ph };
	  vpx.push_back(px);

	  hph[plane].Fill( ph );

	  } // roi px
	*/
      } // size
      else {
	evinf.skip = 1;
	cout << " (" << iev << ": B ROI " << ng/3 << " skip)";
      }

      // column-wise common mode correction:

      double phprev = 0;
      double dphprev = 0;

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

	if( row4 == row1 ) {
	  phprev = ph1;
	  continue; // Randpixel
	}
	if( row4 == row7 ) continue;

	phvsprev[plane].Fill( phprev, ph4 );

	if( run >= 1873 && run <= 1905 ) { // fresh
	  if( plane == A )
	    ph4 -= 0.17*phprev; // Tsunami
	  if( plane == B )
	    ph4 -= 0.14*phprev; // Tsunami
	  if( plane == C )
	    ph4 -= 0.15*phprev; // Tsunami
	}

	phprev = vpx[ipx].ph; // original ph4

	double dph;
	if( row4 - row1 < row7 - row4 )
	  dph = ph4 - ph1;
	else
	  dph = ph4 - ph7;
 
	hdph[plane].Fill( dph ); // sig 2.7

	dphvsprev[plane].Fill( dphprev, dph );
	dphprev = dph;

	double dphcut = 12; // sig 4.5

	if( run >= 1747 ) {

	  dphcut = 40; // gain_2

	  if( plane == B )
	    dphcut = 30; // 136i
	  //dphcut = 40; // 136i

	}

	if( run >= 1757 ) { // fresh

	  //dphcut = 30; // gain_2 1767 dx3cq3 3.39
	  dphcut = 40; // gain_2 1767 dx3cq3 3.18
	  //dphcut = 50; // gain_2 1767 dx3cq3 3.31
	  //dphcut = 60; // gain_2 1767 dx3cq3 

	  if( plane == B )
	    //dphcut = 30; // 1783 dx3cq3 3.10
	    dphcut = 40; // 1783 dx3cq3 3.03
	  //dphcut = 50; // 1767 dx3cq3 3.31

	}

	if( run >= 1789 ) {

	  dphcut = 30; // gain_2, 1820: dx3cq3 5.07

	  if( plane == B )
	    //dphcut = 20; // 130i, 1820: dx3cq3 5.41
	    dphcut = 30; // 130i, 1820: dx3cq3 5.07
	  //dphcut = 40; // 130i, 1820: dx3cq3 5.2307

	}

	if( run >= 1823 ) { // fresh

	  //dphcut = 20; // gain_2 1840: dx3cq3 4.6
	  //dphcut = 30; // gain_2 1840: dx3cq3 4.04
	  //dphcut = 40; // gain_2 1840: dx3cq3 3.91
	  //dphcut = 50; // gain_2 1840: dx3cq3 3.99

	  //dphcut = 20; // gain_2 1840: dx3cq3 3.92
	  dphcut = 30; // gain_2 1840: dx3cq3 3.88
	  if( plane == B )
	    dphcut = 40; // gain_2 1840: dx3cq3 3.88
	}

	if( run >= 1842 ) {

	  dphcut = 20; // gain_2 7.5
	  //dphcut = 30; // gain_2 7.5
	  //dphcut = 40; // gain_2 7.6

	  if( plane == B )
	    //dphcut = 15; // 133i 8.8
	    //dphcut = 20; // 133i 7.7
	    dphcut = 25; // 133i 7.5
	  //dphcut = 30; // 133i 7.6
	}

	if( run >= 1865 ) { // fresh

	  //dphcut = 15; // gain_2, 1894 dx3cq3 4.85 Tsunami corrected
	  dphcut = 20; // gain_2, 1894 dx3cq3 4.65 Tsunami corrected
	  //dphcut = 25; // gain_2, 1894 dx3cq3 4.74 Tsunami
	  //dphcut = 30; // gain_2, 

	}

	if( dph > dphcut ) {

	  hpxmap[plane]->Fill( col4, row4 );

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
	  /*
	    int ia = dph/300 * nb; // Landau peak at 200 ADC
	    if( ia > nb-1 ) ia = nb-1; // overflow
	    dph = ia * 300.0/nb; // bits
	  */
	  px.ph = dph;

	  // r4scal.C

	  double U = ( dph - p3[plane][col4][row4] ) / p2[plane][col4][row4];

	  if( U >= 1 )
	    U = 0.9999999; // avoid overflow

	  double vcal = p0[plane][col4][row4] - p1[plane][col4][row4] * log( (1-U)/U ); // inverse Fermi

	  px.q = ke[plane]*vcal;

	  pb.push_back(px);

	}

      } // ipx

      //if( ldb ) cout << " roi " << vpx.size() << ", hits " << pb.size() << endl;

    } // filled

    // clustering:

    hnht[plane].Fill( pb.size() );

    if( pb.size() > 50 ) pb.clear(); // speed

    vcl = getClus(pb);

    hncl[plane].Fill( vcl.size() );
    if( vcl.size() )
      if( ldb ) cout << "  clusters " << vcl.size() << endl;

    for( unsigned icl = 0; icl < vcl.size(); ++icl ) {

      hclmap[plane]->Fill( vcl[icl].col, vcl[icl].row );

      hclph[plane].Fill( vcl[icl].sum );
      hclq[plane].Fill( vcl[icl].q );

      if( vcl[icl].sum > 55 ) hclsz[plane].Fill( vcl[icl].size );

      // cluster isolation:

      for( unsigned jcl = icl+1; jcl < vcl.size(); ++jcl ) {

	bool done = 0;

	for( unsigned ipx = 0; ipx < vcl[icl].vpix.size(); ++ipx ) {
	  for( unsigned jpx = 0; jpx < vcl[jcl].vpix.size(); ++jpx )

	    if( fabs( vcl[icl].vpix[ipx].col - vcl[jcl].vpix[jpx].col ) < 3 &&
		fabs( vcl[icl].vpix[ipx].row - vcl[jcl].vpix[jpx].row ) < 3 ) {
	      if( vcl[icl].q < vcl[jcl].q ) // Thu 22.3.2018
		vcl[icl].iso = 0; // flag smaller cluster
	      else
		vcl[jcl].iso = 0;
	      done = 1;
	      break; // jpx
	    }
	  if( done ) break; // ipx

	} // ipx

      } // jcl

    } // icl

    evlist.push_back(vcl);
    if( plane == A )
      infoA.push_back(evinf);
    else if( plane == B )
      infoB.push_back(evinf);
    else if( plane == C )
      infoC.push_back(evinf);

  } // events

  cout << endl;
  cout << "done " << Xfile << ", read " << evlist.size() << " events" << endl;

  return evlist;

} // oneplane

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{

  //std::ios_base::sync_with_stdio( false ); // faster read ?

  cout << "main " << argv[0] << " called with " << argc << " arguments" << endl;

  if( argc == 1 ) {
    cout << "give run number" << endl;
    return 1;
  }

  // file name = last argument:

  string runnum( argv[argc-1] );
  int run = atoi( argv[argc-1] );

  // further arguments:

  int Nev = 20*1000*1000;
  double beamEnergy = 5.6; // [GeV]
  double beamDivergenceScaled = 5/beamEnergy; // 5sigma/energy
  bool fifty = 0;

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      Nev = atoi( argv[++i] );

    if( !strcmp( argv[i], "-p" ) )
      beamEnergy = atof( argv[++i] );

    if( !strcmp( argv[i], "-f" ) )
      fifty = 1;

  } // argc

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // gains:

  //string gainA{ "A/r113-scancal-tb21-0921.dat"}; // lots of negative q
  string gainA{ "/home/cmspix/A/r113-scancal-tb21-0923.dat"};
  if( run >=  423 ) gainA = "A/r113-scancal-tb21-0923.dat";
  if( run >=  430 ) gainA = "A/r112-scancal-tb21-0925.dat";
  if( run >=  439 ) gainA = "A/r113-scancal-tb21-0923.dat";
  if( run >=  866 ) gainA = "A/r109-scancal-tb21-1112.dat";
  if( run >=  998 ) gainA = "A/r146-scancal-tb21-1208.dat";
  if( run >= 1010 ) gainA = "A/r163-scancal-tb21-1209.dat";
  if( run >= 1024 ) gainA = "A/r158-scancal-tb21-1210.dat";
  if( run >= 1037 ) gainA = "A/r152-scancal-tb21-1211.dat";
  if( run >= 1037 ) gainA = "A/r152-scancal-tb21-1211.dat";
  if( run >= 1747 ) gainA = "A/scm108-scancal-tb21-drei-2018-03-11-hold20.dat";
  if( run >= 1757 ) gainA = "A/scm146-scancal2-2018-03-12-hold20.dat";
  if( run >= 1784 ) gainA = "A/scm146-scancal2-drei-pr650-sh630-hold17.dat";
  if( run >= 1787 ) gainA = "B/scm148-scancal2-drei-pr650-sh630-2018-03-13-hold17.dat";
  if( run >= 1823 ) gainA = "/home/cmspix/r4sclient/B/scm148-scancal2-drei-pr650-sh630-2018-03-13-hold17.dat";
  if( run >= 1842 ) gainA = "A/scm152-scancal2-drei-2018-03-16-hold20.dat";
  if( run >= 1865 ) gainA = "A/scm152-scancal2-drei-warm-2018-03-16-hold20.dat";
  if( run >= 1872 ) gainA = "/home/cmspix/r4sclient/A/scm160-scancal2-drei-warm-2018-03-17-hold20.dat";

  ke[A] = 0.039; // Landau peak at 11 ke
  if( run >= 423 ) ke[A] = 0.0396; // Landau peak at 11 ke
  if( run >= 432 ) ke[A] = 0.039; // r112
  if( run >= 866 ) ke[A] = 0.038; // r109
  if( run >= 998 ) ke[A] = 0.0354; // r146
  if( run >= 1010 ) ke[A] = 0.0356; // r163
  if( run >= 1024 ) ke[A] = 0.0285; // r158 thicker deep diff at 11 ke
  if( run >= 1037 ) ke[A] = 0.0283; // r152 thicker deep diff at 11 ke
  if( run >= 1747 ) ke[A] = 0.035;
  //if( run >= 1757 ) ke[A] = 0.0374; // 11 ke dphcut 30
  if( run >= 1757 ) ke[A] = 0.0390; // 11 ke dphcut 40
  if( run >= 1764 ) ke[A] = 0.0386; // 11 ke turn  7
  if( run >= 1765 ) ke[A] = 0.0384; // 11 ke turn  8
  if( run >= 1766 ) ke[A] = 0.0382; // 11 ke turn  9
  if( run >= 1767 ) ke[A] = 0.0379; // 11 ke turn 10
  if( run >= 1768 ) ke[A] = 0.0375; // 11 ke turn 11
  if( run >= 1769 ) ke[A] = 0.0371; // 11 ke turn 12
  if( run >= 1770 ) ke[A] = 0.0369; // 11 ke turn 13
  if( run >= 1771 ) ke[A] = 0.0365; // 11 ke turn 14
  if( run >= 1772 ) ke[A] = 0.0361; // 11 ke turn 15
  if( run >= 1773 ) ke[A] = 0.0355; // 11 ke turn 16
  if( run >= 1774 ) ke[A] = 0.0350; // 11 ke turn 17
  if( run >= 1775 ) ke[A] = 0.0344; // 11 ke turn 18
  if( run >= 1776 ) ke[A] = 0.0337; // 11 ke turn 19
  if( run >= 1777 ) ke[A] = 0.0331; // 11 ke turn 20
  if( run >= 1778 ) ke[A] = 0.0345; // 11 ke turn
  if( run >= 1779 ) ke[A] = 0.0379; // 11 ke turn
  if( run >= 1780 ) ke[A] = 0.0375; // 11 ke turn
  if( run >= 1781 ) ke[A] = 0.0381; // 11 ke turn
  if( run >= 1782 ) ke[A] = 0.0381; // 11 ke turn
  if( run >= 1784 ) ke[A] = 0.0400; // 11 ke 146
  if( run >= 1787 ) ke[A] = 0.0513; // 11 ke 148
  if( run >= 1789 ) ke[A] = 0.0505; // 11 ke 148 with cold 130i
  if( run >= 1823 ) ke[A] = 0.0511; // 11 ke 148 with warm 108
  if( run >= 1842 ) ke[A] = 0.0353; // 11 ke 152
  if( run >= 1865 ) ke[A] = 0.0380; // 11 ke 152
  if( run >= 1872 ) ke[A] = 0.0605; // 11 ke 160 dph>25
  if( run >= 1872 ) ke[A] = 0.0624; // 11 ke 160 Tsunami

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
    gainFileA >> p0[A][icol][irow];
    gainFileA >> p1[A][icol][irow];
    gainFileA >> p2[A][icol][irow];
    gainFileA >> p3[A][icol][irow];

  } // while

  // B:

  string gainB{ "B/r108-scancal-tb21-0921.dat" };
  //if( run >= 423 ) gainB = "B/r108-scancal-tb21-0923-hold24.dat";
  if( run >=  423 ) gainB = "B/r108-scancal-tb21-0923-hold25.dat";
  if( run >=  430 ) gainB = "A/r117-scancal-tb21-1005.dat";
  if( run >=  432 ) gainB = "B/r110-scancal-tb21-0925-hold25.dat"; // 12.6 ke
  //if( run >= 432 ) gainB = "B/r110-scancal-tb21-0925-hold26.dat"; // 12.0 ke
  if( run >=  444 ) gainB = "B/r110-scancal-tb21-0928-hold24.dat";
  if( run >=  866 ) gainB = "B/r148-scancal-tb21-1112.dat";
  if( run >=  998 ) gainB = "B/r150-scancal-tb21-1208.dat";
  if( run >= 1010 ) gainB = "B/r146-scancal-tb21-1208.dat";
  if( run >= 1024 ) gainB = "B/r152-scancal-tb21-1210.dat";
  if( run >= 1037 ) gainB = "B/r160-scancal-tb21-2017-12-11.dat";
  if( run >= 1747 ) gainB = "B/scm136i-scancal2-tb21-icy-drei-2018-03-11-hold20.dat";
  if( run >= 1757 ) gainB = "B/scm148-scancal2-2018-03-12-hold16.dat";
  if( run >= 1782 ) gainB = "B/scm148-scancal2-2018-03-12-hold24.dat";
  if( run >= 1784 ) gainB = "B/scm148-scancal2-drei-pr650-sh630-2018-03-13-hold17.dat";
  if( run >= 1787 ) gainB = "A/scm146-scancal2-drei-pr650-sh630-hold17.dat";
  if( run >= 1789 ) gainB = "B/scm130i-scancal2-drei-icy-pr800-sh600-ia119-2018-03-13-hold20.dat";
  if( run >= 1798 ) gainB = "B/scm130i-scancal2-drei-icy-pr800-sh600-ia125-2018-03-13-hold16.dat";
  //  if( run >= 1823 ) gainB = "/home/cmspix/r4sclient/B/c108-scancal2-tb21-2018-02-24-ia125-hold24.dat";
  if( run >= 1823 ) gainB = "/home/cmspix/r4sclient/B/scm108-scancal2-drei-2018-3-15-hold24.dat";
  if( run >= 1842 ) gainB = "B/scm133-scancal2-drei-icy-2018-3-16-hold20.dat";
  if( run >= 1865 ) gainB = "B/scm102-scancal2-drei-warm-2018-03-16-hold20.dat";
  if( run >= 1872 ) gainB = "/home/cmspix/r4sclient/B/scm159-scancal2-drei-warm-2018-03-17-pr650-sh700-hold20.dat";
  //  if( run >= 1872 ) gainB = "/home/cmspix/r4sclient/A/scm159-scancal1-tb21-pr900-sh670-ia125-vb120-hold24.dat";

  ke[B] = 0.0276; // Landau peak at 11 ke
  if( run >= 423 ) ke[B] = 0.026;
  if( run >= 432 ) ke[B] = 0.0326; // r110
  if( run >= 443 ) ke[B] = 0.036; // cmspixel-daq
  if( run >= 866 ) ke[B] = 0.029; // c148
  if( run >= 998 ) ke[B] = 0.0264; // c150 clcut 2
  if( run >= 1010 ) ke[B] = 0.033; // r146
  if( run >= 1024 ) ke[B] = 0.0228; // r152 thicker deep diff at 11 ke
  if( run >= 1037 ) ke[B] = 0.029; // r160
  if( run >= 1747 ) ke[B] = 0.035; // default
  //if( run >= 1757 ) ke[B] = 0.0289; // 11 ke 148 dphcut 30
  if( run >= 1757 ) ke[B] = 0.0307; // 11 ke 148 dphcut 40
  if( run >= 1763 ) ke[B] = 0.0300; // 11 ke turn
  if( run >= 1764 ) ke[B] = 0.0293; // 11 ke turn
  if( run >= 1765 ) ke[B] = 0.0288; // 11 ke turn
  if( run >= 1766 ) ke[B] = 0.0281; // 11 ke turn
  if( run >= 1767 ) ke[B] = 0.0274; // 11 ke turn
  if( run >= 1768 ) ke[B] = 0.0267; // 11 ke turn
  if( run >= 1769 ) ke[B] = 0.0263; // 11 ke turn
  if( run >= 1770 ) ke[B] = 0.0260; // 11 ke turn
  if( run >= 1771 ) ke[B] = 0.0255; // 11 ke turn
  if( run >= 1772 ) ke[B] = 0.0251; // 11 ke turn
  if( run >= 1773 ) ke[B] = 0.0245; // 11 ke turn
  if( run >= 1774 ) ke[B] = 0.0240; // 11 ke turn
  if( run >= 1775 ) ke[B] = 0.0234; // 11 ke turn
  if( run >= 1776 ) ke[B] = 0.0227; // 11 ke turn
  if( run >= 1777 ) ke[B] = 0.0223; // 11 ke turn
  if( run >= 1778 ) ke[B] = 0.0235; // 11 ke turn
  if( run >= 1779 ) ke[B] = 0.0269; // 11 ke turn
  if( run >= 1780 ) ke[B] = 0.0266; // 11 ke turn
  if( run >= 1781 ) ke[B] = 0.0276; // 11 ke turn
  if( run >= 1782 ) ke[B] = 0.0307; // 11 ke
  if( run >= 1783 ) ke[B] = 0.0307; // 11 ke
  if( run >= 1784 ) ke[B] = 0.0505; // 11 ke 148
  if( run >= 1787 ) ke[B] = 0.0402; // 11 ke 146
  if( run >= 1789 ) ke[B] = 0.035; // default
  if( run >= 1823 ) ke[B] = 0.0385; // 11 ke 108
  if( run >= 1842 ) ke[B] = 0.035; // default
  if( run >= 1865 ) ke[B] = 0.038; // 11 ke 102
  if( run >= 1872 ) ke[B] = 0.0377; // 11 ke 159 dph>25
  if( run >= 1872 ) ke[B] = 0.0401; // 11 ke 159 Tsunami

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
    gainFileB >> p0[B][icol][irow];
    gainFileB >> p1[B][icol][irow];
    gainFileB >> p2[B][icol][irow];
    gainFileB >> p3[B][icol][irow];

  } // while

  // C:

  string gainC{ "C/r110-scancal-tb21-0921.dat" };
  if( run >= 423 ) gainC = "C/r110-scancal-tb21-0923-hold25.dat";
  //if( run >= 432 ) gainC = "C/r114-scancal-tb21-0925-hold24.dat"; // 14.9
  if( run >= 432 ) gainC = "C/r114-scancal-tb21-0925-hold25.dat"; // 15.2
  //if( run >= 432 ) gainC = "C/r114-scancal-tb21-0925-hold26.dat"; // 14.9
  if( run >= 866 ) gainC = "C/r110-scancal-tb21-1112.dat";
  if( run >= 998 ) gainC = "C/r148-scancal-tb21-1208.dat";
  if( run >= 1024 ) gainC = "C/r159-scancal-tb21-1210.dat";
  if( run >= 1747 ) gainC = "C/scm109-scancal2-tb21-drei-2018-03-11-hold20.dat";
  if( run >= 1757 ) gainC = "C/scm109-scancal2-tb21-drei-2018-03-12-hold20.dat";
  if( run >= 1784 ) gainC = "C/scm109-scancal2-drei-pr650-sh630-ia125-2018-03-13-hold20.dat";
  if( run >= 1823 ) gainC = "/home/cmspix/r4sclient/C/scm109-scancal2-drei-pr650-sh630-ia125-2018-03-13-hold20.dat";
  if( run >= 1842 ) gainC = "C/scm159-scancal2-drei-2018-3-16-hold20.dat";
  if( run >= 1865 ) gainC = "C/scm159-scancal2-drei-warm-2018-03-16-hold20.dat";
  if( run >= 1872 ) gainC = "/home/cmspix/r4sclient/C/scm102-scanhold-drei-warm-2018-03-17-hold20.dat";

  ke[C] = 0.0366; // Landau peak at 11 ke
  if( run >=  432 ) ke[C] = 0.028; // Landau peak at 11 ke
  if( run >=  866 ) ke[C] = 0.034; // c110
  if( run >=  998 ) ke[C] = 0.041; // c148
  if( run >= 1024 ) ke[C] = 0.039; // c159
  if( run >= 1747 ) ke[C] = 0.040; // c109
  if( run >= 1757 ) ke[C] = 0.0393; // c109
  if( run >= 1758 ) ke[C] = 0.0398; // c109
  if( run >= 1768 ) ke[C] = 0.0393; // c109 turn
  if( run >= 1769 ) ke[C] = 0.0391; // c109 turn
  if( run >= 1770 ) ke[C] = 0.0388; // c109 turn
  if( run >= 1771 ) ke[C] = 0.0386; // c109 turn
  if( run >= 1772 ) ke[C] = 0.0383; // c109 turn
  if( run >= 1773 ) ke[C] = 0.0379; // c109 turn
  if( run >= 1774 ) ke[C] = 0.0374; // c109 turn
  if( run >= 1775 ) ke[C] = 0.0370; // c109 turn
  if( run >= 1776 ) ke[C] = 0.0364; // c109 turn
  if( run >= 1777 ) ke[C] = 0.0359; // c109 turn
  if( run >= 1778 ) ke[C] = 0.0367; // c109 turn
  if( run >= 1779 ) ke[C] = 0.0392; // c109 turn
  if( run >= 1780 ) ke[C] = 0.0392; // c109 turn
  if( run >= 1781 ) ke[C] = 0.0395; // c109 turn
  if( run >= 1782 ) ke[C] = 0.0395; // c109 turn
  if( run >= 1784 ) ke[C] = 0.0422; // c109 turn
  if( run >= 1823 ) ke[C] = 0.0463; // c109 with warm 108
  if( run >= 1842 ) ke[C] = 0.0404; // c159
  if( run >= 1865 ) ke[C] = 0.042; // c159
  if( run >= 1872 ) ke[C] = 0.0416; // c102 dph>25
  if( run >= 1872 ) ke[C] = 0.0429; // c102 Tsunami

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
    gainFileC >> p0[C][icol][irow];
    gainFileC >> p1[C][icol][irow];
    gainFileC >> p2[C][icol][irow];
    gainFileC >> p3[C][icol][irow];

  } // while

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string alignpath = "align/";
  int aligniteration = 0;

  double alignxA = 0.0; // [mm] same sign as dx
  double alignyA = 0.0; // [mm] same sign as dy
  double alignfA = 0.0; // [rad] same sign dxvsy

  double alignxC = 0.0; // [mm] same sign as dx
  double alignyC = 0.0; // [mm] same sign as dy
  double alignfC = 0.0; // [rad] same sign dxvsy

  string alignFileName = alignpath+"align_" + runnum + ".dat";

  ifstream alignFile( alignFileName );

  cout << endl;

  if( alignFile.bad() || ! alignFile.is_open() ) {
    cout << "no " << alignFileName << ", will bootstrap" << endl;
  }
  else {

    cout << "read alignment from " << alignFileName << endl;

    string HASH( "#" );
    string ITER( "iteration" );
    string ALXA( "alignxA" );
    string ALYA( "alignyA" );
    string ALFA( "alignfA" );
    string ALXC( "alignxC" );
    string ALYC( "alignyC" );
    string ALFC( "alignfC" );

    while( ! alignFile.eof() ) {

      string line;
      getline( alignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == HASH ) // comments start with #
	continue;

      if( tag == ITER )
	tokenizer >> aligniteration;

      double val;
      tokenizer >> val;
      if(      tag == ALXA )
	alignxA = val;
      else if( tag == ALYA )
	alignyA = val;
      else if( tag == ALFA )
	alignfA = val;
      else if( tag == ALXC )
	alignxC = val;
      else if( tag == ALYC )
	alignyC = val;
      else if( tag == ALFC )
	alignfC = val;

      // anything else on the line and in the file gets ignored

    } // while getline

    alignFile.close();

  } // alignFile

  double cfA = cos(alignfA);
  double sfA = sin(alignfA);
  double cfC = cos(alignfC);
  double sfC = sin(alignfC);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  TFile * histoFile = new TFile( Form( "/home/zoiirene/Output/drei-r%i.root", run ), "RECREATE" );

  // book histos:

  int nbx =  80;
  int nby = 320;
  if( fifty ) {
    nbx = 160;
    nby = 160;
  }

  for( unsigned ipl = 0; ipl < 3; ++ipl ) {

    phvsprev[ipl] = TProfile( Form( "phvsprev%s", PN[ipl].c_str() ),
			      Form( "%s Tsunami;previous PH [ADC];%s <PH> [ADC]",
				    PN[ipl].c_str(), PN[ipl].c_str() ),
			      80, 0, 800, -999, 1999 );
    dphvsprev[ipl] = TProfile( Form( "dphvsprev%s", PN[ipl].c_str() ),
			      Form( "%s Tsunami;previous #DeltaPH [ADC];%s <#DeltaPH> [ADC]",
				    PN[ipl].c_str(), PN[ipl].c_str() ),
			       80, 0, 800, -999, 1999 );

    hph[ipl] = TH1I( Form( "ph%s", PN[ipl].c_str() ),
		     Form("%s PH;ADC-PED [ADC];%s pixels",
			  PN[ipl].c_str(), PN[ipl].c_str() ),
		     500, -100, 900 );
    hdph[ipl] = TH1I( Form( "dph%s", PN[ipl].c_str() ),
		      Form( "%s #DeltaPH;#DeltaPH [ADC];%s pixels",
			    PN[ipl].c_str(), PN[ipl].c_str() ),
		      500, -100, 900 );
    hnpx[ipl] = TH1I( Form( "npx%s", PN[ipl].c_str() ),
		      Form( "%s ROI pixels per event;ROI pixels;%s events",
			    PN[ipl].c_str(), PN[ipl].c_str() ),
		      200, 0, 1000 );
    hnht[ipl] = TH1I( Form( "nht%s", PN[ipl].c_str() ),
		      Form( "%s pixel hits per event;pixel hits;%s events",
			    PN[ipl].c_str(), PN[ipl].c_str() ),
		      50, 0.5, 50.5 );
    hpxmap[ipl] = new TH2I( Form( "pxmap%s", PN[ipl].c_str() ),
			    Form( "%s pixel map, PH > cut;col;row;%s PH pixels",
				  PN[ipl].c_str(), PN[ipl].c_str() ),
			    155, -0.5, 154.5, 160, -0.5, 159.5 );

    hncl[ipl] = TH1I( Form( "ncl%s", PN[ipl].c_str() ),
		      Form( "%s cluster per event;clusters;%s events",
			    PN[ipl].c_str(), PN[ipl].c_str() ),
		      21, -0.5, 20.5 );
    hclmap[ipl] = new TH2I( Form( "clmap%s", PN[ipl].c_str() ),
			    Form( "%s cluster map;col;row;%s clusters",
				  PN[ipl].c_str(), PN[ipl].c_str() ),
			    nbx, 0, nbx, nby, 0, nby );
    hclsz[ipl] = TH1I( Form( "clsz%s", PN[ipl].c_str() ),
		       Form( "%s cluster size;cluster size [pixels];%s clusters",
			     PN[ipl].c_str(), PN[ipl].c_str() ),
		       20, 0.5, 20.5 );
    hclph[ipl] = TH1I( Form( "clph%s", PN[ipl].c_str() ),
		       Form( "%s cluster PH;cluster ph [ADC];%s clusters",
			     PN[ipl].c_str(), PN[ipl].c_str() ),
		       200, 0, 1000 );
    hclq[ipl] = TH1I( Form( "clq%s", PN[ipl].c_str() ),
		      Form( "%s cluster charge;cluster charge [ke];%s clusters",
			    PN[ipl].c_str(), PN[ipl].c_str() ),
		      100, 0, 50 );

  } // ipl

  TH1I hdt( "dt", "time between events;log_{10}(#Deltat [s]);events", 100, -4, 1 );
  TH1I hddtAB( "ddtAB", "dtA - dtB;A-B #delta#Deltat [clocks];events", 100, -1000, 1000 );
  TH1I hddtCB( "ddtCB", "dtC - dtB;C-B #delta#Deltat [clocks];events", 100, -1000, 1000 );
  TH1I hddtCA( "ddtCA", "dtC - dtA;C-A #delta#Deltat [clocks];events", 100, -1000, 1000 );
  TProfile ddtvsdtAB( "ddtvsdtAB",
		      "A-B time lag vs intervall;log_{10}(#Deltat [s]);<#delta#Deltat> [clocks]",
		      100, -4, 1,-1e99, 1e99 );

  // correlations:

  TH1I hxA( "xA", "x A;x [mm];clusters A", 100, -5, 5 );
  TH1I hyA( "yA", "y A;y [mm];clusters A", 100, -5, 5 ); 
  TH1I hxAi( "xAi", "x A isolated;x [mm];isolated clusters A", 100, -5, 5 );
  TH1I hyAi( "yAi", "y A isolated;y [mm];isolated clusters A", 100, -5, 5 ); 
  TH1I hclqAi( "clqAi", "A isolated cluster charge;cluster charge [ke];A isolatewd clusters",
	       100, 0, 50 );

  TH2I * hxxAB = new TH2I( "xxAB", "B vs A;row A;row B;clusters", 320, -4, 4, 320, -4, 4 );
  TH2I * hyyAB = new TH2I( "yyAB", "B vs A;col A;col B;clusters",  80, -4, 4,  80, -4, 4 );

  double f = 0.5;
  //double f = 0.1; // aligned

  TH1I hdxAB( "dxAB", "Bx-Ax;x-x [mm];cluster pairs", 800, -2, 2 );
  TH1I hdyAB( "dyAB", "By-Ay;y-y [mm];cluster pairs", 400, -2, 2 );
  TProfile dxvsxAB( "dxvsxAB", "dx vs x A-B;x [mm];<dx> [mm]", 320, -4, 4, -f, f );
  TProfile dxvsyAB( "dxvsyAB", "dx vs y A-B;y [mm];<dx> [mm]",  80, -4, 4, -f, f );

  TH2I * hdxvsev = new
    TH2I( "dxvsev", "Bx-Ax vs events;events;#Deltax [px];clusters",
	  100, 0, 10000, 100, -f, f );

  TProfile nmvsevAB( "nmvsevAB", "AB matches vs time;time [events];AB matches",
		     3100, 0, 3100*1000, -1, 99 );

  TH1I hxB( "xB", "x B;x [mm];clusters B", 100, -5, 5 );
  TH1I hyB( "yB", "y B;y [mm];clusters B", 100, -5, 5 ); 
  TH1I hxBi( "xBi", "x B isolated;x [mm];isolated clusters B", 100, -5, 5 );
  TH1I hyBi( "yBi", "y B isolated;y [mm];isolated clusters B", 100, -5, 5 ); 
  TH1I hclqBi( "clqBi", "B isolated cluster charge;cluster charge [ke];B isolatewd clusters",
	       100, 0, 50 );

  TH2I * hxxCB = new TH2I( "xxCB", "C vs B;row B;row C;clusters", 320, -4, 4, 320, -4, 4 );
  TH2I * hyyCB = new TH2I( "yyCB", "C vs B;col B;col C;clusters",  80, -4, 4,  80, -4, 4 );

  TH1I hdxCB( "dxCB", "Cx-Bx;x-x [mm];cluster pairs", 800, -2, 2 );
  TH1I hdyCB( "dyCB", "Cy-By;y-y [mm];cluster pairs", 400, -2, 2 );
  TProfile dxvsxCB( "dxvsxCB", "dx vs x C-B;x [mm];<dx> [mm]", 320, -4, 4, -f, f );
  TProfile dxvsyCB( "dxvsyCB", "dx vs y C-B;y [mm];<dx> [mm]",  80, -4, 4, -f, f );
  TProfile nmvsevCB( "nmvsevCB", "CB matches vs time;time [events];CB matches",
		     3100, 0, 3100*1000, -1, 99 );

  // triplets:

  TH1I hxC( "xC", "x C;x [mm];clusters C", 100, -5, 5 );
  TH1I hyC( "yC", "y C;y [mm];clusters C", 100, -5, 5 ); 
  TH1I hxCi( "xCi", "x C isolated;x [mm];isolated clusters C", 100, -5, 5 );
  TH1I hyCi( "yCi", "y C isolated;y [mm];isolated clusters C", 100, -5, 5 ); 
  TH1I hclqCi( "clqCi", "C isolated cluster charge;cluster charge [ke];C isolatewd clusters",
	       100, 0, 50 );

  TH2I * hxxCA = new TH2I( "xxCA", "C vs A;row A;row C;clusters", 320, -4, 4, 320, -4, 4 );
  TH2I * hyyCA = new TH2I( "yyCA", "C vs A;col A;col C;clusters",  80, -4, 4,  80, -4, 4 );

  TH1I hdxCA( "dxCA", "Cx-Ax;x-x [mm];cluster pairs", 400, -1, 1 );
  TH1I hdyCA( "dyCA", "Cy-Ay;y-y [mm];cluster pairs", 200, -1, 1 );

  TProfile dxvsxCA( "dxvsxCA", "dx vs x C-A;x [mm];<dx> [mm]", 320, -4, 4, -f, f );
  TProfile dxvsyCA( "dxvsyCA", "dx vs y C-A;y [mm];<dx> [mm]",  80, -4, 4, -f, f );
  TH1I hdyCAc( "dyCAc", "Cy-Ay, cut dx;y-y [mm];cluster pairs", 200, -1, 1 );
  TProfile dyvsyCA( "dyvsyCA", "dy vs y C-A;y [mm];<dy> [mm]",  80, -4, 4, -f, f );
  TProfile nmvsevCA( "nmvsevCA", "CA matches vs time;time [events];CA matches",
		     3100, 0, 3100*1000, -1, 99 );

  TH1I hdx3( "dx3", "triplet dx;dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I hdy3( "dy3", "triplet dy;dy [mm];triplets", 200, -1, 1 );

  TH1I hdx3c( "dx3c", "triplet dx, cut dy;dx [mm];triplets", 500, -0.25, 0.25 );
  TH1I hdx3ci( "dx3ci", "triplet dx, cut dy, isolated;dx [mm];isolated triplets", 500, -0.25, 0.25 );

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

  TH1I hdx3m( "dx3m", "triplet dx, x < 0;dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I hdx3p( "dx3p", "triplet dx, x > 0;dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I hdx3ct( "dx3ct", "triplet dx, cut dy, tx;dx [mm];triplets",
	       500, -0.25, 0.25 );
  TProfile madx3vsq( "madx3vsq", "MAD(dx) vs Q;B cluster charge [ke];MAD dx [mm]",
		     100, 0, 100, 0, 0.1 );
  TProfile madx3vsn( "madx3vsn", "MAD(dx) vs cluster size;B cluster size [pixels];MAD dx [mm]",
		     20, 0.5, 20.5, 0, 0.1 );

  TH1I hdx3cq( "dx3cq", "triplet dx, Landau peak;dx [mm];Landau peak triplets",
	       500, -0.25, 0.25 );
  TH1I hdx3cqi( "dx3cqi", "triplet dx, Landau peak, isolated;dx [mm];isolated Landau peak triplets",
		500, -0.25, 0.25 );
  TH1I hdx3cq3( "dx3cq3", "triplet dx, 3 Landau peak;dx [mm];Landau peak triplets",
		500, -0.25, 0.25 );
  TH1I hdx3cq3i( "dx3cq3i", "triplet dx, 3 Landau peak, isolated;dx [mm];isolated Landau peak triplets",
		 500, -0.25, 0.25 );

  TProfile dx3vsev( "dx3vsev", "dx3 vs time;trigger;<dx3> [mm]",
		    310, 0, 3100*1000, -0.5, 0.5 );

  TProfile dx3vsx( "dx3vsx", "dx vs x;x [mm];<dx3> [mm]", 320, -4, 4, -0.5, 0.5 );
  TProfile dx3vsy( "dx3vsy", "dx vs y;y [mm];<dx3> [mm]",  80, -4, 4, -0.5, 0.5 );
  TProfile dx3vsxm( "dx3vsxm", "dx vs x mod 50 um;x mod 50 [#mum];<dx3> [mm]",
		    50, 0, 50, -0.5, 0.5 );

  TProfile madx3vsdx( "madx3vsdx", "MAD(dx3) vs dx C-A;C-A dx [#mum];MAD dx3 [mm]",
		      100, -100, 100, 0, 0.1 );

  TH1I hdx3cq3t( "dx3cq3t",
		 "triplet dx, 3 Landau peak, forward;dx [mm];Landau peak forward triplets",
		 500, -0.25, 0.25 );
  TProfile madx3vsx( "madx3vsx", "MAD(dx3) vs x;x [mm];MAD dx3 [mm]", 320, -4, 4, 0, 0.1 );
  TProfile madx3vsy( "madx3vsy", "MAD(dx3) vs y;y [mm];MAD dx3 [mm]",  80, -4, 4, 0, 0.1 );
  TProfile madx3vsxm( "madx3vsxm", "MAD(dx3) vs xmod;x mod 50 [#mum];MAD dx3 [mm]",
		      50, 0, 50, 0, 0.1 );

  TProfile etavsxmB3( "etavsxmB3", "eta vs xmod;x mod 50 [#mum];B <eta>",
		      50, 0, 50, -1.1, 1.1 );
  TProfile madx3vseta( "madx3vseta", "MAD(dx3) vs eta;eta;MAD dx3 [mm]",
		       100, -1, 1, 0, 0.1 );
  TH1I hdx3cq3t2( "dx3cq3t2",
		  "triplet dx, 3 Landau peak, forward, 2-px;dx [mm];Landau peak forward triplets",
		  500, -0.25, 0.25 );

  TH2I * hclmapB3 = new
    TH2I( "clmapB3", "linked cluster map B;col;row;B clusters on tracks",
	  80, 0, 80, 320, 0, 320 );

  TH1I hxA3( "xA3", "x A linked;x [mm];A clusters on tracks", 100, -5, 5 );
  TH1I hyA3( "yA3", "y A linked;y [mm];A clusters on tracks", 100, -5, 5 ); 
  TH1I hxB3( "xB3", "x B linked;x [mm];B clusters on tracks", 100, -5, 5 );
  TH1I hyB3( "yB3", "y B linked;y [mm];B clusters on tracks", 100, -5, 5 ); 
  TH1I hxC3( "xC3", "x C linked;x [mm];C clusters on tracks", 100, -5, 5 );
  TH1I hyC3( "yC3", "y C linked;y [mm];C clusters on tracks", 100, -5, 5 ); 

  TH1I hclszA3( "clszA3", "A cluster size on tracks;cluster size [pixels];Aclusters on tracks",
		40, 0.5, 40.5 );
  TH1I hclphA3( "clphA3", "A cluster PH on tracks;cluster ph [ADC];A clusters on tracks",
		200, 0, 1000 );
  TH1I hclqA3( "clqA3", "A cluster charge on tracks;cluster charge [ke];A clusters on tracks",
	       160, 0, 80 );

  TH1I hclszB3( "clszB3", "B cluster size on tracks;cluster size [pixels];B clusters on tracks",
		40, 0.5, 40.5 );
  TH1I hncolB3( "ncolB3", "B cluster size on tracks;cluster size [columns];B clusters on tracks",
		20, 0.5, 20.5 );
  TH1I hnrowB3( "nrowB3", "B cluster size on tracks;cluster size [rows];B clusters on tracks",
		20, 0.5, 20.5 );
  TH1I hclphB3( "clphB3", "B cluster PH on tracks;cluster ph [ADC];B clusters on tracks",
		200, 0, 1000 );
  TH1I hclqB3( "clqB3", "B cluster charge on tracks;cluster charge [ke];B clusters on tracks",
	       160, 0, 80 );
  TH1I hclqB3i( "clqB3i", "B cluster charge on tracks;cluster charge [ke];B clusters on tracks",
	       80, 0, 20 );
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
		       "B rows vs xmod;x mod 50 [#mum];<B cluster size [rows]>",
		       50, 0, 50, 0.5, 10.5 );
  TProfile clqvsxmB3( "clqvsxmB3",
		      "B cluster charge vs xmod;x mod 50 [#mum];<B cluster charge [ke]>",
		      50, 0, 50, 0, 50 );

  TH1I hetaA3( "etaA3", "A cluster eta;eta;A 2-pix clusters on tracks",
	       100, -1, 1 );
  TH1I hetaB3( "etaB3", "B cluster eta;eta;B 2-pix clusters on tracks",
	       100, -1, 1 );
  TH1I hetaC3( "etaC3", "C cluster eta;eta;C 2-pix clusters on tracks",
	       100, -1, 1 );

  TH1I hpxqA3( "pxqA3", "A pixel charge;pixel charge [ke];A pixels on tracks",
	       100, 0, 20 );
  TH1I hpxqB3( "pxqB3", "B pixel charge;pixel charge [ke];B pixels on tracks",
	       100, 0, 20 );
  TH1I hpxqC3( "pxqC3", "C pixel charge;pixel charge [ke];C pixels on tracks",
	       100, 0, 20 );

  TH1I hpxpA3( "pxpA3", "A pixel PH;pixel PH [ADC];A pixels on tracks", 250, 0, 500 );
  TH1I hpxpB3( "pxpB3", "B pixel PH;pixel PH [ADC];B pixels on tracks", 250, 0, 500 );
  TH1I hpxpC3( "pxpC3", "C pixel PH;pixel PH [ADC];C pixels on tracks", 250, 0, 500 );

  TH1I hpxq1stB3( "pxq1stB3", "B 1st pixel charge;pixel charge [ke];B `st pixels on tracks",
		  100, 0, 20 );
  TH1I hpxq2ndB3( "pxq2ndB3", "B 2nd pixel charge;pixel charge [ke];B `st pixels on tracks",
		  100, 0, 20 );

  TProfile effvsdxy( "effvsdxy",
		     "DUT efficiency vs triplet dxy;xy match radius [mm];DUT efficiency",
		     1000, 0, 10, -0.1, 1.1 );

  TProfile2D * effvsxy =
    new TProfile2D( "effvsxy",
		    "DUT efficiency map;x [mm];y[mm];DUT efficiency",
		    80, -4, 4, 80, -4, 4, -0.1, 1.1 );
  TProfile effvsx( "effvsx", "eff vs x;x [mm];DUT efficiency",
		   320, -4, 4, -0.1, 1.1 );
  TProfile effvsy( "effvsy", "eff vs y;y [mm];DUT efficiency",
		   80, -4, 4, -0.1, 1.1 );
  TProfile effvsxm( "effvsxm", "eff vs x mod 50;x mod 50 [#mum];DUT efficiency",
		    50, 0, 50, -0.1, 1.1 );

  TProfile effvsev( "effvsev", "eff vs time;trigger;DUT efficiency",
		    3100, 0, 3100*1000, -0.1, 1.1 );
  TProfile effvsiev( "effvsiev", "eff vs event;event mod 200;DUT efficiency",
		     100, -0.5, 199.5, -0.1, 1.1 );
  TProfile effvsmpxA( "effvsmpxA", "eff vs occupancy A;occupancy A [pixels];DUT efficiency",
		      50, 0.5, 50.5, -0.1, 1.1 );
  TProfile effvsqA( "effvsqA", "eff vs charge A;cluster charge A [ke];DUT efficiency",
		    100, 0, 100, -0.1, 1.1 );
  TProfile effvstxy( "effvstxy", "eff vs angle;dxy CA [mm];DUT efficiency",
		     100, 0, 0.2, -0.1, 1.1 );

  const double log10 = log(10);
  string ADD {"A"}; // added flag

  double qL =  9;
  double qR = 15;

  double qLB = qL;
  double qRB = qR;

  if( run >= 1789 && run <= 1822 ) { // 130i peak at 4
    qLB = 2;
    qRB = 7;
  }

  if( run >= 1842 && run <= 1864 ) { // 133i peak at 5
    qLB = 3;
    qRB = 8;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  timespec ts;
  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s0 = ts.tv_sec; // seconds since 1.1.1970
  long f0 = ts.tv_nsec; // nanoseconds

  clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ); // sum of threads time
  time_t t0 = ts.tv_sec; // seconds
  long n0 = ts.tv_nsec; // nanoseconds

  double zeit1 = 0;
  double zeit2 = 0;

  // simulate n-bit ADC: run 884

  //int nb = pow( 2, 8 ); // 885 dx3cq3  3.04
  //int nb = pow( 2, 7 ); // 885 dx3cq3  3.05
  //int nb = pow( 2, 6 ); // 885 dx3cq3  3.07
  //int nb = pow( 2, 5 ); // 885 dx3cq3  3.11
  //int nb = pow( 2, 4 ); // 885 dx3cq3  3.29
  //int nb = pow( 2, 3 ); // 885 dx3cq3  3.895
  //int nb = pow( 2, 2 ); // 885 dx3cq3  4.95
  //int nb = pow( 2, 1 ); // 885 dx3cq3  6.78

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Read and process the run for each plane in parallel:

  list < vector < cluster > > evlistA;
  list < vector < cluster > > evlistB;
  list < vector < cluster > > evlistC;

  //#pragma omp sections // test, not parallel
#pragma omp parallel sections
  {
#pragma omp section
    {
      //cout << "A " << sched_getcpu() << endl << flush; // changes

      evlistA = oneplane( A, runnum, Nev, fifty );

    } //omp section

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // B:

#pragma omp section
    {
      //cout << "B " << sched_getcpu() << endl << flush; // different from A

      evlistB = oneplane( B, runnum, Nev, fifty );

    } // omp section

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // C:

#pragma omp section
    {
      //cout << "C " << sched_getcpu() << endl << flush; // different from A and B

      evlistC = oneplane( C, runnum, Nev, fifty );

    } // omp section

  } // omp parallel sections

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s2 = ts.tv_sec; // seconds since 1.1.1970
  long f2 = ts.tv_nsec; // nanoseconds
  zeit1 += s2 - s0 + ( f2 - f0 ) * 1e-9; // read and cluster

  clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ); // sum threads
  time_t t9 = ts.tv_sec; // seconds
  long n9 = ts.tv_nsec; // nanoseconds

  cout << "time " << zeit1 << " s"
       << " (sum threads " << t9 - t0 + ( n9 - n0 ) * 1e-9 << " s)"
       << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // loop over events, correlate planes:

  auto evA = evlistA.begin();
  auto evB = evlistB.begin();
  auto evC = evlistC.begin();

  auto evinfoA = infoA.begin();
  auto evinfoB = infoB.begin();
  auto evinfoC = infoC.begin();

  unsigned iev = 0;
  uint64_t prevtimeA = 0;
  uint64_t prevtimeB = 0;
  uint64_t prevtimeC = 0;

  cout << endl << "tracking ev" << flush;

  for( ; evinfoB != infoB.end() && evA != evlistA.end() && evB != evlistB.end() && evC != evlistC.end();
       ++evinfoA, ++evinfoB, ++evinfoC, ++evA, ++evB, ++evC ) {

    vector <cluster> vclA = *evA;
    vector <cluster> vclB = *evB;
    vector <cluster> vclC = *evC;

    ++iev;
    if( iev%10000 == 0 )
      cout << " " << iev << flush;

    if( evinfoA->filled == ADD ) cout << endl << "ev " << iev << " added A" << endl;
    if( evinfoB->filled == ADD ) cout << endl << "ev " << iev << " added B" << endl;
    if( evinfoC->filled == ADD ) cout << endl << "ev " << iev << " added C" << endl;

    uint64_t dtA = evinfoA->evtime - prevtimeA;
    prevtimeA = evinfoA->evtime;

    uint64_t dtB = evinfoB->evtime - prevtimeB;
    prevtimeB = evinfoB->evtime;
    if( iev > 1 && dtB > 0 ) // added events have same time
      hdt.Fill( log(dtB/40e6) / log10 ); // MHz -> s

    uint64_t dtC = evinfoC->evtime - prevtimeC;
    prevtimeC = evinfoC->evtime;

    long ddtAB = dtA - dtB;
    long ddtCB = dtC - dtB;
    long ddtCA = dtC - dtA;
    if( iev > 1 && dtB > 0 ) {
      hddtAB.Fill( ddtAB );
      hddtCB.Fill( ddtCB );
      hddtCA.Fill( ddtCA );
      ddtvsdtAB.Fill( log(dtB/40e6) / log10, ddtAB );
    }
    if( abs( ddtAB ) > 1999 )
      cout << endl << "ev " << iev
	   << " dtA " << dtA
	   << ", dtB " << dtB << " (" << dtB/40E6 << "s)"
	   << ", ddtAB " << ddtAB
	   << endl;
    if( abs( ddtCB ) > 1999 )
      cout << endl << "ev " << iev
	   << " dtC " << dtC
	   << ", dtB " << dtB << " (" << dtB/40E6 << "s)"
	   << ", ddtCB " << ddtCB
	   << endl;
    if( abs( ddtCA ) > 1999 )
      cout << endl << "ev " << iev
	   << " dtC " << dtC
	   << ", dtA " << dtA << " (" << dtA/40E6 << "s)"
	   << ", ddtCA " << ddtCA
	   << endl;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // A-B cluster correlations:

    double ptchc = 0.100; // [mm] col size
    double ptchr = 0.025; // [mm] row size
    if( fifty ) {
      ptchc = 0.050;
      ptchr = 0.050;
    }
    int nm = 0;

    for( vector<cluster>::iterator cA = vclA.begin(); cA != vclA.end(); ++cA ) {

      double xA = cA->row*ptchr - halfSensorX - alignxA; // rot90 Dreimaster
      double yA = cA->col*ptchc - halfSensorY - alignyA; // 100 um px
      if( fifty ) {
	xA = cA->col*ptchc - halfSensorY - alignxA; // straight
	yA = cA->row*ptchr - halfSensorX - alignyA; // PCB
      }

      double xAr = xA*cfA - yA*sfA;
      double yAr = xA*sfA + yA*cfA;

      hxA.Fill( xAr );
      hyA.Fill( yAr );
      if( cA->iso ) {
	hxAi.Fill( xAr );
	hyAi.Fill( yAr );
	hclqAi.Fill( cA->q );
      }

      for( vector<cluster>::iterator cB = vclB.begin(); cB != vclB.end(); ++cB ) {

	double xB = cB->row*ptchr - halfSensorX;
	double yB = cB->col*ptchc - halfSensorY;
	if( run == 431 ) {
	  xB = cB->row*ptchr - halfSensorX; // rot90
	  yB = cB->col*ptchc - halfSensorY; // PCB
	}
	if( fifty ) {
	  xB = cB->col*ptchc - halfSensorY; // straight
	  yB = cB->row*ptchr - halfSensorX; // PCB
	}

	hxxAB->Fill( xAr, xB );
	hyyAB->Fill( yAr, yB );

	double dx = xAr - xB;
	double dy = yAr - yB;

	if( cA->q > qL  && cA->q < qR &&
	    cB->q > qLB && cB->q < qRB &&
	    cA->iso && cB->iso ) {

	  hdxAB.Fill( dx );
	  hdyAB.Fill( dy );
	  dxvsxAB.Fill( xB, dx );
	  dxvsyAB.Fill( yB, dx );

	}

	hdxvsev->Fill( iev, dx );

	if( fabs( dx ) < straightTracks * beamDivergenceScaled + 0.020 &&
	    fabs( dy ) < straightTracks * beamDivergenceScaled + 0.100 )
	  ++nm;

      } // clusters

    } // cl

    nmvsevAB.Fill( iev, nm );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // B-C cluster correlations:

    nm = 0;

    for( vector<cluster>::iterator cB = vclB.begin(); cB != vclB.end(); ++cB ) {

      double xB = cB->row*ptchr - halfSensorX;
      double yB = cB->col*ptchc - halfSensorY;
      if( run == 431 ) {
	xB = cB->row*ptchr - halfSensorX; // rot90
	yB = cB->col*ptchc - halfSensorY; // PCB
      }
      if( fifty ) {
	xB = cB->col*ptchc - halfSensorY; // straight
	yB = cB->row*ptchr - halfSensorX; // PCB
      }

      hxB.Fill( xB );
      hyB.Fill( yB );
      if( cB->iso ) {
	hxBi.Fill( xB );
	hyBi.Fill( yB );
	hclqBi.Fill( cB->q );
      }

      for( vector<cluster>::iterator cC = vclC.begin(); cC != vclC.end(); ++cC ) {

	double xC = cC->row*ptchr - halfSensorX - alignxC; // rot90 Dreimaster
	double yC = cC->col*ptchc - halfSensorY - alignyC; // down
	if( fifty ) {
	  xC = cC->col*ptchc - halfSensorY - alignxC; // straight
	  yC = cC->row*ptchr - halfSensorX - alignyC; // PCB
	}

	double xCr = xC*cfC - yC*sfC;
	double yCr = xC*sfC + yC*cfC;

	hxxCB->Fill( xB, xCr );
	hyyCB->Fill( yB, yCr );

	double dx = xCr - xB;
	double dy = yCr - yB;

	if( cC->q > qL  && cC->q < qR &&
	    cB->q > qLB && cB->q < qRB &&
	    cC->iso && cB->iso ) {

	  hdxCB.Fill( dx );
	  hdyCB.Fill( dy );
	  dxvsxCB.Fill( xB, dx );
	  dxvsyCB.Fill( yB, dx );

	}

	if( fabs( dx ) < straightTracks * beamDivergenceScaled + 0.020 &&
	    fabs( dy ) < straightTracks * beamDivergenceScaled + 0.100 )
	  ++nm;

      } // clusters B

    } // cl C

    nmvsevCB.Fill( iev, nm );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // A-C cluster correlations:

    nm = 0;

    for( vector<cluster>::iterator cC = vclC.begin(); cC != vclC.end(); ++cC ) {

      double xC = cC->row*ptchr - halfSensorX - alignxC;
      double yC = cC->col*ptchc - halfSensorY - alignyC;
      if( fifty ) {
	xC = cC->col*ptchc - halfSensorY - alignxC; // straight
	yC = cC->row*ptchr - halfSensorX - alignyC; // PCB
      }

      double xCr = xC*cfC - yC*sfC;
      double yCr = xC*sfC + yC*cfC;

      hxC.Fill( xCr );
      hyC.Fill( yCr );
      if( cC->iso ) {
	hxCi.Fill( xCr );
	hyCi.Fill( yCr );
	hclqCi.Fill( cC->q );
      }

      double etaC = -2;
      if( cC->size == 2 ) {
	double q0 = cC->vpix[0].q;
	double q1 = cC->vpix[1].q;
	etaC = (q1-q0)/(q1+q0);
      }

      for( vector<cluster>::iterator cA = vclA.begin(); cA != vclA.end(); ++cA ) {

	double xA = cA->row*ptchr - halfSensorX - alignxA; // rot90 Dreimaster
	double yA = cA->col*ptchc - halfSensorY - alignyA; // down
	if( fifty ) {
	  xA = cA->col*ptchc - halfSensorY - alignxA; // straight
	  yA = cA->row*ptchr - halfSensorX - alignyA; // PCB
	}

	double xAr = xA*cfA - yA*sfA;
	double yAr = xA*sfA + yA*cfA;

	hxxCA->Fill( xAr, xCr );
	hyyCA->Fill( yAr, yCr );

	double dxCA = xCr - xAr;
	double dyCA = yCr - yAr;
	double dxyCA = sqrt( dxCA*dxCA + dyCA*dyCA );
	hdxCA.Fill( dxCA );
	hdyCA.Fill( dyCA );

	if( fabs( dxCA ) > straightTracks * beamDivergenceScaled + 0.02 ) continue; // includes beam divergence: +-5 sigma

	hdyCAc.Fill( dyCA );
	dyvsyCA.Fill( yAr, dyCA );

	if( fabs( dyCA ) > straightTracks * beamDivergenceScaled + 0.1 ) continue; // [mm]

	++nm;

	dxvsyCA.Fill( yAr, dxCA );
	dxvsxCA.Fill( xAr, dxCA ); // linear trend in run 392, 403: acceptance and beam divergence ?

	double xavg = 0.5 * ( xAr + xCr );
	double yavg = 0.5 * ( yAr + yCr );

	double xmod = fmod( xavg + 8, 0.05 ); // [mm] 0..0.05
	if( fifty )
	  xmod = fmod( xavg + 8.025, 0.05 ); // [mm] 0..0.05

	double etaA = -2;
	if( cA->size == 2 ) {
	  double q0 = cA->vpix[0].q;
	  double q1 = cA->vpix[1].q;
	  etaA = (q1-q0)/(q1+q0);
	}

	int eff[999] = {0};

	for( vector<cluster>::iterator cB = vclB.begin(); cB != vclB.end(); ++cB ) {

	  double xB = cB->row*ptchr - halfSensorX;
	  double yB = cB->col*ptchc - halfSensorY;
	  if( run == 431 ) {
	    xB = cB->row*ptchr - halfSensorX; // rot90
	    yB = cB->col*ptchc - halfSensorY; // PCB
	  }
	  if( fifty ) {
	    xB = cB->col*ptchc - halfSensorY; // straight
	    yB = cB->row*ptchr - halfSensorX; // PCB
	  }

	  double etaB = -2;
	  if( cB->size == 2 ) {
	    double q0 = cB->vpix[0].q;
	    double q1 = cB->vpix[1].q;
	    etaB = (q1-q0)/(q1+q0);
	  }

	  int rowmin = 999;
	  int rowmax = 0;
	  int colmin = 999;
	  int colmax = 0;
	  for( int ipx = 0; ipx < cB->size; ++ipx ) {
	    int row = cB->vpix[ipx].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;
	    int col = cB->vpix[ipx].col;
	    if( col < colmin ) colmin = col;
	    if( col > colmax ) colmax = col;
	  }
	  int nrowB = rowmax-rowmin+1;
	  int ncolB = colmax-colmin+1;

	  // triplet residual:

	  double dx3 = xB - xavg;

	  if( run == 447 || run == 449 )
	    dx3 -= 0.00076*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 866 )
	    dx3 += 0.0011*xavg; // from -dx3vsx
	  if( run == 871 )
	    dx3 += 0.0018*xavg; // from -dx3vsx
	  if( run == 872 )
	    dx3 += 0.0018*xavg; // from -dx3vsx
	  if( run == 873 )
	    dx3 += 0.0017*xavg; // from -dx3vsx
	  if( run == 875 )
	    dx3 += 0.0013*xavg; // from -dx3vsx
	  if( run == 876 )
	    dx3 += 0.0011*xavg; // from -dx3vsx
	  if( run == 877 )
	    dx3 += 0.0011*xavg; // from -dx3vsx
	  if( run == 878 )
	    dx3 += 0.0013*xavg; // from -dx3vsx
	  if( run == 879 )
	    dx3 += 0.0010*xavg; // from -dx3vsx
	  if( run == 880 )
	    dx3 += 0.0012*xavg; // from -dx3vsx
	  if( run == 881 )
	    dx3 += 0.0011*xavg; // from -dx3vsx
	  if( run == 882 )
	    dx3 += 0.0011*xavg; // from -dx3vsx
	  if( run == 883 )
	    dx3 += 0.0013*xavg; // from -dx3vsx
	  if( run == 884 )
	    dx3 += 0.0010*xavg; // from -dx3vsx
	  if( run == 885 )
	    dx3 += 0.0010*xavg; // from -dx3vsx
	  if( run == 998 )
	    dx3 -= 0.0005*xavg; // from -dx3vsx
	  if( run == 999 )
	    dx3 -= 0.0005*xavg; // from -dx3vsx
	  if( run == 1000 )
	    dx3 -= 0.0010*xavg; // from -dx3vsx
	  if( run == 1006 )
	    dx3 -= 0.0005*xavg; // from -dx3vsx
	  if( run >= 1010 && run <= 1020 )
	    dx3 -= 0.00086*xavg; // from -dx3vsx
	  if( run == 1025 )
	    dx3 -= 0.00084*xavg; // from -dx3vsx
	  if( run == 1026 )
	    dx3 -= 0.00082*xavg; // from -dx3vsx
	  if( run == 1027 )
	    dx3 -= 0.00072*xavg; // from -dx3vsx
	  if( run >= 1028 && run <= 1035 )
	    dx3 -= 0.00082*xavg; // from -dx3vsx
	  if( run >= 1036 && run <= 1045 )
	    dx3 -= 0.00085*xavg; // from -dx3vsx
	  if( run == 1046 )
	    dx3 += 0.00036*xavg; // from -dx3vsx
	  if( run == 1047 )
	    dx3 -= 0.00082*xavg; // from -dx3vsx
	  if( run == 1764 )
	    dx3 -= 0.00067*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1765 )
	    dx3 -= 0.00085*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1766 )
	    dx3 -= 0.00106*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1767 )
	    dx3 -= 0.00119*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1768 )
	    dx3 -= 0.00139*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1769 )
	    dx3 -= 0.00152*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1770 )
	    dx3 -= 0.00165*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1771 )
	    dx3 -= 0.00191*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1772 )
	    dx3 -= 0.00206*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1773 )
	    dx3 -= 0.00206*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1774 )
	    dx3 -= 0.00261*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1775 )
	    dx3 -= 0.00293*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1776 )
	    dx3 -= 0.00316*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1777 )
	    dx3 -= 0.00337*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1778 )
	    dx3 -= 0.00031*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1779 )
	    dx3 += 0.00049*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1780 )
	    dx3 -= 0.00055*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1781 )
	    dx3 += 0.00027*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1782 )
	    dx3 += 0.00028*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1783 )
	    dx3 += 0.00028*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1784 )
	    dx3 += 0.00028*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1785 )
	    dx3 += 0.00028*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1786 )
	    dx3 -= 0.00001*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1787 )
	    dx3 -= 0.00117*xavg; // from -dx3vsx.Fit("pol1")

	  if( run == 1789 )
	    dx3 -= 0.00118*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1789 )
	    dx3 -= 0.00213*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1799 )
	    dx3 -= 0.00214*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1800 )
	    dx3 -= 0.00218*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1801 )
	    dx3 -= 0.00225*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1802 )
	    dx3 -= 0.00240*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1803 )
	    dx3 -= 0.00240*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1804 )
	    dx3 -= 0.00240*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1805 )
	    dx3 -= 0.00240*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1806 )
	    dx3 -= 0.00240*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1807 )
	    dx3 -= 0.00240*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1808 )
	    dx3 -= 0.00240*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1809 )
	    dx3 -= 0.00240*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1810 )
	    dx3 -= 0.00220*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1811 )
	    dx3 -= 0.00220*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1812 )
	    dx3 -= 0.00220*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1813 )
	    dx3 -= 0.00220*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1814 )
	    dx3 -= 0.00220*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1815 )
	    dx3 -= 0.00220*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1816 )
	    dx3 -= 0.00230*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1817 )
	    dx3 -= 0.00230*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1818 )
	    dx3 -= 0.00240*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1819 )
	    dx3 -= 0.00250*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1820 )
	    dx3 -= 0.00220*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1822 )
	    dx3 -= 0.00218*xavg; // from -dx3vsx.Fit("pol1")

	  if( run == 1825 )
	    dx3 += 0.00043*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1826 )
	    dx3 += 0.00036*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1827 )
	    dx3 -= 0.00064*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1828 )
	    dx3 -= 0.00064*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1829 )
	    dx3 -= 0.00064*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1830 )
	    dx3 -= 0.00064*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1831 )
	    dx3 -= 0.00064*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1832 )
	    dx3 -= 0.00064*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1833 )
	    dx3 -= 0.00064*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1834 )
	    dx3 -= 0.00064*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1835 )
	    dx3 -= 0.00064*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1836 )
	    dx3 -= 0.00064*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1837 )
	    dx3 -= 0.00064*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1838 )
	    dx3 -= 0.00064*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1839 )
	    dx3 -= 0.00064*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1840 )
	    dx3 -= 0.00054*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1843 )
	    dx3 -= 0.00200*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1845 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1846 )
	    dx3 -= 0.00251*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1847 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1848 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1849 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1850 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1851 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1852 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1853 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1854 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1855 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1856 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1857 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1858 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1859 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1860 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1861 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1862 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1863 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1864 )
	    dx3 -= 0.00238*xavg; // from -dx3vsx.Fit("pol1")

	  if( run == 1865 )
	    dx3 -= 0.00049*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1866 )
	    dx3 -= 0.00110*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1868 )
	    dx3 -= 0.00155*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1869 )
	    dx3 -= 0.00063*xavg; // from -dx3vsx.Fit("pol1")

	  if( run == 1873 )
	    dx3  = 0.000000018*xavg; // from -dx3vsx.Fit("pol1") - opposite of the slope from a linear fit has to be added as correction
	  if( run == 1874 )
	    dx3 -= 0.00067*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1875 )
	    dx3 -= 0.00098*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1876 )
	    dx3 -= 0.00137*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1877 )
	    dx3 -= 0.00163*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1878 )
	    dx3 -= 0.00170*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1879 )
	    dx3 -= 0.00178*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1880 )
	    dx3 -= 0.00200*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1881 )
	    dx3 -= 0.00224*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1882 )
	    dx3 -= 0.00242*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1883 )
	    dx3 -= 0.00257*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1884 )
	    dx3 -= 0.00275*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1885 )
	    dx3 -= 0.00300*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1886 )
	    dx3 -= 0.00323*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1887 )
	    dx3 -= 0.00341*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1888 )
	    dx3 -= 0.00367*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1889 )
	    dx3 -= 0.00374*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1890 )
	    dx3 -= 0.00400*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1891 )
	    dx3 -= 0.00088*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1892 )
	    dx3 -= 0.00068*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1893 )
	    dx3 -= 0.00068*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1894 )
	    dx3 -= 0.00068*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1895 )
	    dx3 -= 0.00068*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1896 )
	    dx3 -= 0.00068*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1897 )
	    dx3 -= 0.00068*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1898 )
	    dx3 -= 0.00068*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1899 )
	    dx3 -= 0.00068*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1900 )
	    dx3 -= 0.00068*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1901 )
	    dx3 -= 0.00068*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1902 )
	    dx3 -= 0.00068*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1903 )
	    dx3 -= 0.00068*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1904 )
	    dx3 -= 0.00068*xavg; // from -dx3vsx.Fit("pol1")
	  if( run == 1905 )
	    dx3 -= 0.00068*xavg; // from -dx3vsx.Fit("pol1")

	  double dy3 = yB - yavg;
	  double dxy = sqrt( dx3*dx3 + dy3*dy3 );

	  hdx3.Fill( dx3 );
	  hdy3.Fill( dy3 );

	  if( fabs( dy3 ) < straightTracks * beamDivergenceScaled + 0.05 ) { // cut on y, look at x, see madxvsy

	    hdx3c.Fill( dx3 );

	    if( dx3 > 0.04 && dx3 < -0.06 ) { // side lobe
	      cout << endl;
	      cout << "x: " << xAr << ", " << xB << ", " << xCr << ", dx3 " << dx3 << endl;
	      cout << "A:";
	      for( unsigned icl = 0; icl < vclA.size(); ++icl )
		cout << " (" << vclA[icl].col << ", " << vclA[icl].row << ", " << vclA[icl].q << ")";
	      cout << endl;
	      cout << "B:";
	      for( unsigned icl = 0; icl < vclB.size(); ++icl )
		cout << " (" << vclB[icl].col << ", " << vclB[icl].row << ", " << vclB[icl].q << ")";
	      cout << endl;
	      cout << "C:";
	      for( unsigned icl = 0; icl < vclC.size(); ++icl )
		cout << " (" << vclC[icl].col << ", " << vclC[icl].row << ", " << vclC[icl].q << ")";
	      cout << endl;
	    }

	    if( cA->iso && cB->iso && cC->iso ) {

	      hdx3ci.Fill( dx3 );

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

	      if( xB < 0 )
		hdx3m.Fill( dx3 );
	      else
		hdx3p.Fill( dx3 );

	      if( fabs( dxCA ) < straightTracks * beamDivergenceScaled ) { // track angle
		hdx3ct.Fill( dx3 );
		madx3vsq.Fill( cB->q, fabs(dx3) );
		madx3vsn.Fill( cB->size, fabs(dx3) );
	      }

	    } // iso

	    if( cB->q > qLB && cB->q < qRB ) {

	      hdx3cq.Fill( dx3 );

	      if( cA->iso && cB->iso && cC->iso )
		hdx3cqi.Fill( dx3 );
	      
	      if( cA->q > qL && cA->q < qR &&
		  cC->q > qL && cC->q < qR ) {

		hdx3cq3.Fill( dx3 );

		if( cA->iso && cB->iso && cC->iso )
		  hdx3cq3i.Fill( dx3 );

		dx3vsev.Fill( iev, dx3 );

		dx3vsx.Fill( xB, dx3 ); // turn
		dx3vsy.Fill( yB, dx3 ); // rot
		dx3vsxm.Fill( xmod*1E3, dx3 );

		madx3vsdx.Fill( dxCA*1E3, fabs(dx3) ); // dxCA

		if( fabs( dxCA ) < straightTracks * beamDivergenceScaled ) { // track angle

		  hdx3cq3t.Fill( dx3 ); // 447 4.27 um

		  madx3vsx.Fill( xB, fabs(dx3) );
		  madx3vsy.Fill( yB, fabs(dx3) );
		  madx3vsxm.Fill( xmod*1E3, fabs(dx3) );
		  if( cB->size == 2 ) {
		    etavsxmB3.Fill( xmod*1E3, etaB ); // sine
		    madx3vseta.Fill( etaB, fabs(dx3) ); // flat
		    hdx3cq3t2.Fill( dx3 ); // 447 4.25 um
		  }

		} // angle

	      } // Qa, qC

	    } // qB

	  } // cut dy

	  if( fabs( dx3 ) < 0.07 && // hit on track
	      fabs( dy3 ) < 0.15 &&
	      cA->iso && cB->iso && cC->iso) {

	    hclmapB3->Fill( cB->col, cB->row );

	    hxA3.Fill( xAr );
	    hyA3.Fill( yAr );
	    hxB3.Fill( xB  );
	    hyB3.Fill( yB  );
	    hxC3.Fill( xCr );
	    hyC3.Fill( yCr );

	    hclszA3.Fill( cA->size );
	    hclszB3.Fill( cB->size );
	    hncolB3.Fill( ncolB );
	    hnrowB3.Fill( nrowB );
	    hclszC3.Fill( cC->size );

	    hclphA3.Fill( cA->sum );
	    hclphB3.Fill( cB->sum );
	    hclphC3.Fill( cC->sum );

	    hclqA3.Fill( cA->q );
	    hclqB3.Fill( cB->q );
	    hclqB3i.Fill( cB->q );
	    hclqC3.Fill( cC->q );
	    if( cB->size < 4 )
	      hclqB3n.Fill( cB->q );

	    if( fifty )
	      nrowvsxmB3.Fill( xmod*1E3, ncolB );
	    else
	      nrowvsxmB3.Fill( xmod*1E3, nrowB );

	    clqvsxmB3.Fill( xmod*1E3, cB->q );

	    if( cA->size == 2 )
	      hetaA3.Fill( etaA );

	    if( cB->size == 2 )
	      hetaB3.Fill( etaB );

	    if( cC->size == 2 )
	      hetaC3.Fill( etaC );

	    for( int ipx = 0; ipx < cA->size; ++ipx ) {
	      hpxpA3.Fill( cA->vpix[ipx].ph );
	      hpxqA3.Fill( cA->vpix[ipx].q );
	    }
	    for( int ipx = 0; ipx < cB->size; ++ipx ) {
	      hpxpB3.Fill( cB->vpix[ipx].ph );
	      hpxqB3.Fill( cB->vpix[ipx].q );
	    }
	    for( int ipx = 0; ipx < cC->size; ++ipx ) {
	      hpxpC3.Fill( cC->vpix[ipx].ph );
	      hpxqC3.Fill( cC->vpix[ipx].q );
	    }

	    if( cB->size == 2 ) {
	      hpxq1stB3.Fill( cB->vpix[0].q );
	      hpxq2ndB3.Fill( cB->vpix[1].q ); // identical
	    }

	    // task: store track

	  } // linked, iso

	  for( int iw = 1; iw < 999; ++iw )
	    if( dxy < iw*0.010 )
	      eff[iw] = 1; // eff

	} // clusters B

	if( fabs( dyCA ) > straightTracks*beamDivergenceScaled ) continue; // clean reference "tracks"

	if( evinfoB->filled == ADD ) continue; // padded event in B
	if( evinfoB->skip ) continue; // fat event in B

	if( cA->iso == 0 ) continue;
	if( cC->iso == 0 ) continue;

	effvsxy->Fill( xavg, yavg, eff[50] );

	if( yavg > -3.7 && yavg < 3.5 )
	  effvsx.Fill( xavg, eff[50] );

	if( xavg > -3.3 && xavg < 3.6 )
	  effvsy.Fill( yavg, eff[50] );

	if( xavg > -3.3 && xavg < 3.6 &&
	    yavg > -3.7 && yavg < 3.5 ) {

	  for( int iw = 1; iw < 999; ++iw )
	    effvsdxy.Fill( iw*0.010+0.005, eff[iw] );

	  effvsxm.Fill( xmod*1E3, eff[50] ); // bias dot
	  effvsev.Fill( iev, eff[50] );
	  effvsiev.Fill( iev%200, eff[50] );
	  //effvsmpxA.Fill( pbA.size(), eff[50] );
	  effvsqA.Fill( cA->q, eff[50] );
	  effvstxy.Fill( dxyCA, eff[50] ); // flat

	} // fiducial x, y

      } // clusters A

    } // cl C

    nmvsevCA.Fill( iev, nm );

    // task: track correlations, intersects
  } // events

  cout << endl;

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s3 = ts.tv_sec; // seconds since 1.1.1970
  long f3 = ts.tv_nsec; // nanoseconds
  zeit2 += s3 - s2 + ( f3 - f2 ) * 1e-9; // read and cluster
  cout << "time " << zeit2 << " s" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s9 = ts.tv_sec; // seconds since 1.1.1970
  long f9 = ts.tv_nsec; // nanoseconds

  cout << "full time " << s9 - s0 + ( f9 - f0 ) * 1e-9 << " s"
       << " (read and cluster " << zeit1 << " s, tracking " << zeit2 << " s)"
       << endl;

  rusage usage;
  getrusage( RUSAGE_SELF, &usage );
  cout << endl << "resource usage:"
       << endl << "time " << usage.ru_utime.tv_sec + usage.ru_utime.tv_usec*1E-6 << " s"
       << endl << "max resident set size " << usage.ru_maxrss << " kB"
       << endl << "file inputs " << usage.ru_inblock
       << endl << "file ouputs " << usage.ru_oublock
       << endl << "page faults without I/O " << usage.ru_minflt
       << endl << "page faults with I/O " << usage.ru_majflt
       << endl << "voluntary context switches " << usage.ru_nvcsw
       << endl << "involuntary context switches " << usage.ru_nivcsw
       << endl;

  histoFile->Write();
  histoFile->Close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // alignment fits:

  double newalignxA = alignxA;

  cout << endl << hdxAB.GetTitle() << " entries " << hdxAB.GetEntries() << endl;

  if( hdxAB.GetEntries() > 999 ) {

    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -10, 10 );
    double xpk = hdxAB.GetBinCenter( hdxAB.GetMaximumBin() );
    fgp0->SetParameter( 0, hdxAB.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, xpk );
    fgp0->SetParameter( 2, hdxAB.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, hdxAB.GetBinContent( hdxAB.FindBin(xpk-1) ) ); // BG
    hdxAB.Fit( "fgp0", "q", "", xpk-1, xpk+1 ); // fit range around peak
    cout << "  Fit Gauss + BG:"
	 << endl << "  area " << fgp0->GetParameter(0)
	 << endl << "  mean " << fgp0->GetParameter(1)
	 << endl << "  sigm " << fgp0->GetParameter(2)
	 << endl << "  offs " << fgp0->GetParameter(3)
	 << endl;
    newalignxA += fgp0->GetParameter(1);
  }
  else
    cout << "  not enough" << endl;

  // y:

  double newalignyA = alignyA;

  cout << endl << hdyAB.GetTitle() << " entries " << hdyAB.GetEntries() << endl;
  if( hdyAB.GetEntries() > 999 ) {
    cout << "  y correction " << hdyAB.GetMean() << endl;
    newalignyA += hdyAB.GetMean();
  }
  else
    cout << "  not enough" << endl;

  // dxvsy -> -rot

  double newalignfA = alignfA;

  cout << endl << dxvsyAB.GetTitle() << " entries " << dxvsyAB.GetEntries() << endl;
  if( aligniteration > 0 && dxvsyAB.GetEntries() > 999 ) {
    dxvsyAB.Fit( "pol1", "q", "", -3, 3 );
    TF1 * fdxvsy = dxvsyAB.GetFunction( "pol1" );
    cout << "  extra rot " << fdxvsy->GetParameter(1) << endl;
    newalignfA += fdxvsy->GetParameter(1);
  }
  else
    cout << "  not enough" << endl;

  // C:

  double newalignxC = alignxC;

  cout << endl << hdxCB.GetTitle() << " entries " << hdxCB.GetEntries() << endl;
  if( hdxCB.GetEntries() > 999 ) {

    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -10, 10 );
    double xpk = hdxCB.GetBinCenter( hdxCB.GetMaximumBin() );
    fgp0->SetParameter( 0, hdxCB.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, xpk );
    fgp0->SetParameter( 2, hdxCB.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, hdxCB.GetBinContent( hdxCB.FindBin(xpk-1) ) ); // BG
    hdxCB.Fit( "fgp0", "q", "", xpk-1, xpk+1 ); // fit range around peak
    cout << "  Fit Gauss + BG:"
	 << endl << "  area " << fgp0->GetParameter(0)
	 << endl << "  mean " << fgp0->GetParameter(1)
	 << endl << "  sigm " << fgp0->GetParameter(2)
	 << endl << "  offs " << fgp0->GetParameter(3)
	 << endl;
    newalignxC += fgp0->GetParameter(1);
  }
  else
    cout << "  not enough" << endl;

  // y:

  double newalignyC = alignyC;

  cout << endl << hdyCB.GetTitle() << " entries " << hdyCB.GetEntries() << endl;
  if( hdyCB.GetEntries() > 999 ) {
    cout << "  y correction " << hdyCB.GetMean() << endl;
    newalignyC += hdyCB.GetMean();
  }
  else
    cout << "  not enough" << endl;

  // dxvsy -> -rot

  double newalignfC = alignfC;

  cout << endl << dxvsyCB.GetTitle() << " entries " << dxvsyCB.GetEntries() << endl;
  if( aligniteration > 0 && dxvsyCB.GetEntries() > 999 ) {
    dxvsyCB.Fit( "pol1", "q", "", -3, 3 );
    TF1 * fdxvsy = dxvsyCB.GetFunction( "pol1" );
    cout << "  extra rot " << fdxvsy->GetParameter(1) << endl;
    newalignfC += fdxvsy->GetParameter(1);
  }
  else
    cout << "  not enough" << endl;

  ++aligniteration;
  cout << endl
       << "for " << alignFileName << endl
       << "iteration " << aligniteration << endl
       << "alignxA " << setw(11) << newalignxA << endl
       << "alignyA " << setw(11) << newalignyA << endl
       << "alignfA " << setw(11) << newalignfA << endl
       << "alignxC " << setw(11) << newalignxC << endl
       << "alignyC " << setw(11) << newalignyC << endl
       << "alignfC " << setw(11) << newalignfC << endl
    ;

  // cout << "update alignment file? (y/n)" << endl;
  // string ans;
  // cin >> ans;
  // string YES{"y"};
  // if( ans == YES ) {

  //   ofstream alignFile( alignFileName );

  //   alignFile << "# alignment for run " << run << endl;
  //   alignFile << "iteration " << aligniteration << endl;
  //   alignFile << "alignxA " << setw(11) << newalignxA << endl;
  //   alignFile << "alignyA " << setw(11) << newalignyA << endl;
  //   alignFile << "alignfA " << setw(11) << newalignfA << endl;
  //   alignFile << "alignxC " << setw(11) << newalignxC << endl;
  //   alignFile << "alignyC " << setw(11) << newalignyC << endl;
  //   alignFile << "alignfC " << setw(11) << newalignfC << endl;

  //   alignFile.close();
  // }

  cout << endl << histoFile->GetName() << endl;

  cout << endl;

  return 0;
}
