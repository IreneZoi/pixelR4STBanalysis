
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
// ****************************use flag: '-p beamEnergy' to change the beam energy! ***********************************
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
#include "./drei_reader.h"


#define r4sRows 160
#define r4sColumns 155

#define halfSensorX 4.0
#define halfSensorY 3.9
#define ACspacing 40 //[mm] distance between planes A and C in the dreimaster


using namespace std;
bool PRINT = false;

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

double p0[DreiMasterPlanes][r4sColumns][r4sRows]; // Fermi
double p1[DreiMasterPlanes][r4sColumns][r4sRows];
double p2[DreiMasterPlanes][r4sColumns][r4sRows];
double p3[DreiMasterPlanes][r4sColumns][r4sRows];
double ke[DreiMasterPlanes];


list < evInfo > infoA;
list < evInfo > infoB;
list < evInfo > infoC;


double nSigmaTolerance = 3;
double beamDivergence = 0.001; //1 mrad
double straightTracks = ACspacing*beamDivergence*nSigmaTolerance;


//functions definition
vector<cluster> getClus( vector <pixel> pb, int fCluCut = 1 ); // 1 = no gap
list < vector < cluster > > oneplane( int plane, string runnum, unsigned Nev, bool fifty );
void getGain( string gainfile, double (*p0)[r4sColumns][r4sRows], double (*p1)[r4sColumns][r4sRows], double (*p2)[r4sColumns][r4sRows], double (*p3)[r4sColumns][r4sRows], int plane);

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
  double beamEnergy = -1.0; // [GeV]
  double beamDivergenceScaled = 5/beamEnergy; // 5sigma/energy
  bool fifty = 0;
  int alignversion = 1; 

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      Nev = atoi( argv[++i] );

    // if( !strcmp( argv[i], "-p" ) )
    //   beamEnergy = atof( argv[++i] );

    // if( !strcmp( argv[i], "-f" ) )
    //   fifty = 1;

    if( !strcmp( argv[i], "-a" ) )
      alignversion = atoi( argv[++i] );

  } // argc

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //ALIGN
  
  string alignpath = "align/";
  int aligniteration = 0;

  double alignxA = 0.0; // [mm] same sign as dx
  double alignyA = 0.0; // [mm] same sign as dy
  double alignfA = 0.0; // [rad] same sign dxvsy

  double alignxC = 0.0; // [mm] same sign as dx
  double alignyC = 0.0; // [mm] same sign as dy
  double alignfC = 0.0; // [rad] same sign dxvsy
  string gainA = "a ";
  string gainB = "b ";
  string gainC = "c ";
  double keA = 0.0;
  double keB = 0.0;
  double keC = 0.0;
  int pitch = 0; 
  double ptchc = 0.0;// 0.100; // [mm] col size
  double ptchr = 0.0;//0.025; // [mm] row size
  
  string alignFileName = "0";
  if(alignversion == 1)  
    alignFileName = alignpath+"align_" + runnum + ".dat";
  if(alignversion == 2)  
    alignFileName = alignpath+"align_v2_" + runnum + ".dat";

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
    string GAINA( "gainA" );
    string GAINB( "gainB" );
    string GAINC( "gainC" );
    string KEA( "keA" );
    string KEB( "keB" );
    string KEC( "keC" );
    string BEAMENERGY( "beamEnergy" );
    string PITCH( "pitch" );
    
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

      if( tag == ALXA )
	tokenizer >>	alignxA; 
      else if( tag == ALYA )
	tokenizer >> 	alignyA;
      else if( tag == ALFA )
	tokenizer >> 	alignfA;
      else if( tag == ALXC )
	tokenizer >> 	alignxC;
      else if( tag == ALYC )
	tokenizer >> 	alignyC;
      else if( tag == ALFC )
	tokenizer >> 	alignfC;
      else if( tag == GAINA )
	tokenizer >> 	gainA;
      else if( tag == GAINB )
	tokenizer >> 	gainB;
      else if( tag == GAINC )
	tokenizer >> 	gainC;
      else if( tag == KEA )
	tokenizer >> 	keA;
      else if( tag == KEB )
	tokenizer >> 	keB;
      else if( tag == KEC )
	tokenizer >> 	keC;
      else if( tag == BEAMENERGY )
	tokenizer >> 	beamEnergy;
      else if( tag == PITCH )
	tokenizer >> pitch;


      // anything else on the line and in the file gets ignored

    } // while getline

    alignFile.close();

  } // alignFile

  double cfA = cos(alignfA);
  double sfA = sin(alignfA);
  double cfC = cos(alignfC);
  double sfC = sin(alignfC);

  if(PRINT)  cout << "Gains " <<  gainA << " " << gainB << " " << gainC << endl;
  if(PRINT)  cout << " alignfA " << alignfA << endl;
  

  ke[A] = keA; // Landau peak at 11 ke
  ke[B] = keB; // Landau peak at 11 ke
  ke[C] = keC; // Landau peak at 11 ke


  getGain( gainA, p0, p1, p2, p3, 0);
  getGain( gainB, p0, p1, p2, p3, 1);
  getGain( gainC, p0, p1, p2, p3, 2);

  if(pitch == 25)
    {
      ptchc = 0.100; // [mm] col size
      ptchr = 0.025; // [mm] row size
    }
  if(pitch == 50)
    {
      fifty = true;
      ptchc = 0.050;
      ptchr = 0.050;
    }


  
  // (re-)create root file:

  TFile * histoFile = new TFile( Form( "/home/zoiirene/Output/drei-r%i_irene.root", run ), "RECREATE" );

  // book histos:

  int nbx =  80; //number of bins x
  int nby = 320; //number of bins y
  if( fifty ) {
    nbx = 160;
    nby = 160;
  }

  for( unsigned ipl = 0; ipl < DreiMasterPlanes; ++ipl ) {

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
			    r4sColumns, -0.5, 154.5, r4sRows, -0.5, 159.5 );

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


  const double log10 = log(10);
  string ADD {"A"}; // added flag

  //selection for Landau peak
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

	  if( fabs( dy3 ) < straightTracks * beamDivergenceScaled + 0.05 ) { // cut on y, look at x, see madx3vsy

	    hdx3c.Fill( dx3 );

	    if( (cB->q <= qLB || cB->q >= qRB) && ( cA->q <= qL || cA->q >= qR) && ( cC->q <= qL || cC->q >= qR))
	      hdx3nocq3.Fill( dx3);

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

	    if( cB->iso )
	      hdx3ci.Fill( dx3 );
 
	    if( cA->iso && cC->iso )
	      hdx3cii.Fill( dx3 );

	    if( cA->iso && cB->iso && cC->iso ) {

	      hdx3ciii.Fill( dx3 );

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



//------------------------------------------------------------------------------
vector<cluster> getClus( vector <pixel> pb, int fCluCut ) // 1 = no gap
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
}//getClus


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
	 
  while(Xstream.good() && ! Xstream.eof() && evlist.size() < Nev ) {

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

      if( vcl[icl].sum > 55 ) hclsz[plane].Fill( vcl[icl].size );//55 is historical adc, for noise suppression

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



void getGain( string gainfile, double (*p0)[r4sColumns][r4sRows], double (*p1)[r4sColumns][r4sRows], double (*p2)[r4sColumns][r4sRows], double (*p3)[r4sColumns][r4sRows], int plane)
{
  ifstream gainFile( gainfile );

  if( ! gainFile ) {
    cout << "gain file " << gainfile << " not found" << endl;
  }
  
  while( ! gainFile.eof() ) {

    int icol;
    int irow;
    gainFile >> icol;
    gainFile >> irow;
    gainFile >> p0[plane][icol][irow];
    gainFile >> p1[plane][icol][irow];
    gainFile >> p2[plane][icol][irow];
    gainFile >> p3[plane][icol][irow];
    
  } // while

}
