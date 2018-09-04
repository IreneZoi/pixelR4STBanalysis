
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
  //  bool fifty = 0;
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
  //read from ALIGN file, now containing also calibration, pitch, beam energy, etc.
  
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
  double beamEnergy = -1.0; // [GeV]
  double qL =  0;
  double qR = 0;
  double qLB = 0;
  double qRB = 0;
  double Tsunami[DreiMasterPlanes];
  double dphcut[DreiMasterPlanes];
  
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
    string QL( "qL" );
    string QR( "qR" );
    string QLB( "qLB" );
    string QRB( "qRB" );
    string TSUNAMIA( "TsunamiA" );    
    string TSUNAMIB( "TsunamiB" );    
    string TSUNAMIC( "TsunamiC" );
    string DPHCUTA( "dphcutA" );    
    string DPHCUTB( "dphcutB" );    
    string DPHCUTC( "dphcutC" );
    
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
	tokenizer >>  beamEnergy;
      else if( tag == PITCH )
	tokenizer >> pitch;
      else if( tag == QL )
	tokenizer >> qL;
      else if( tag == QR )
	tokenizer >> qR;
      else if( tag == QLB )
	tokenizer >> qLB;
      else if( tag == QRB )
	tokenizer >> qRB;
      else if( tag == TSUNAMIA )
	tokenizer >> Tsunami[0];
      else if( tag == TSUNAMIB )
	tokenizer >> Tsunami[1];
      else if( tag == TSUNAMIC )
	tokenizer >> Tsunami[2];
      else if( tag == DPHCUTA )
	tokenizer >> dphcut[0];
      else if( tag == DPHCUTB )
	tokenizer >> dphcut[1];
      else if( tag == DPHCUTC )
	tokenizer >> dphcut[2];

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

  double beamDivergenceScaled = 5/beamEnergy; // 5sigma/energy
  //---------------------------------------------------------------------------------- 
  
  // (re-)create root file:

  TFile * histoFile = new TFile( Form( "/home/zoiirene/Output/drei-r%i_irene.root", run ), "RECREATE" );

  // book histos:
  if(PRINT) cout << "***** going to book hists ***********" << endl;
  bookHists();
  if(PRINT)   cout << hdxAB->GetTitle() << " entries " << hdxAB->GetEntries() << endl;

  const double log10 = log(10);
  string ADD {"A"}; // added flag

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
      if(PRINT) cout << "A " << sched_getcpu() << endl << flush; // changes
      if(PRINT) cout << "***** Plane A ***********" << endl;

      evlistA = oneplane( A, runnum, Nev, fifty, Tsunami[A],dphcut[A] );

    } //omp section

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // B:

#pragma omp section
    {
      if(PRINT) cout << "B " << sched_getcpu() << endl << flush; // different from A
      if(PRINT) cout << "***** Plane B ***********" << endl;

      evlistB = oneplane( B, runnum, Nev, fifty, Tsunami[B],dphcut[B]  );

    } // omp section

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // C:

#pragma omp section
    {
      if(PRINT)cout << "C " << sched_getcpu() << endl << flush; // different from A and B
      if(PRINT) cout << "***** Plane C ***********" << endl;

      evlistC = oneplane( C, runnum, Nev, fifty, Tsunami[C],dphcut[C]  );

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

  if(PRINT) cout << "***** CLUSTERS CORRELATION ***********" << endl;

  cout << endl << "tracking ev" << flush;

  ///////     Big Loop on the event ///////
  for( ; evinfoB != infoB.end() && evA != evlistA.end() && evB != evlistB.end() && evC != evlistC.end();       ++evinfoA, ++evinfoB, ++evinfoC, ++evA, ++evB, ++evC )
    {

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
      
      uint64_t dtC = evinfoC->evtime - prevtimeC;
      prevtimeC = evinfoC->evtime;

      if( iev > 1 && dtB > 0 ) // added events have same time
	hdt->Fill( log(dtB/40e6) / log10 ); // MHz -> s
      
      
      long ddtAB = dtA - dtB;
      long ddtCB = dtC - dtB;
      long ddtCA = dtC - dtA;
      if( iev > 1 && dtB > 0 )
	{
	  hddtAB->Fill( ddtAB );
	  hddtCB->Fill( ddtCB );
	  hddtCA->Fill( ddtCA );
	  ddtvsdtAB->Fill( log(dtB/40e6) / log10, ddtAB );
	}
      
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ////////////////                A-B cluster correlations:    //////////////////////////
      
      int nm = 0;


      
      for( vector<cluster>::iterator cA = vclA.begin(); cA != vclA.end(); ++cA )
	{
	  if(PRINT) cout << " cluster cA" << endl;


	  // double xA = cA->row*ptchr - halfSensorX - alignxA; // rot90 Dreimaster
	  // double yA = cA->col*ptchc - halfSensorY - alignyA; // 100 um px
	  // if( fifty )
	  //   {
	  //     xA = cA->col*ptchc - halfSensorY - alignxA; // straight
	  //     yA = cA->row*ptchr - halfSensorX - alignyA; // PCB
	  //   }
	  // if(PRINT) cout <<  " cA->row " << cA->row << " ptchr " << ptchr << " halfSensorX " << halfSensorX <<  " alignxA " << alignxA << endl;; // rot90 Dreimaster
	
	  double xA = xcoordinate(0, cA, alignxA, ptchc, ptchr);
	  double yA = ycoordinate(0, cA, alignyA, ptchc, ptchr);
	  double xAr = xA*cfA - yA*sfA;
	  double yAr = xA*sfA + yA*cfA;
	  //	if(PRINT) cout << " xAr "<< xAr << endl;

	  hxA->Fill( xAr );
	  hyA->Fill( yAr );
	  if( cA->iso )
	    {
	      hxAi->Fill( xAr );
	      hyAi->Fill( yAr );
	      hclqAi->Fill( cA->q );
	    }
	
	  for( vector<cluster>::iterator cB = vclB.begin(); cB != vclB.end(); ++cB )
	    {
	  
	      // double xB = cB->row*ptchr - halfSensorX;
	      // double yB = cB->col*ptchc - halfSensorY;
	      // if( run == 431 )
	      // 	{
	      // 	  xB = cB->row*ptchr - halfSensorX; // rot90
	      // 	  yB = cB->col*ptchc - halfSensorY; // PCB
	      // 	}
	      // if( fifty )
	      // 	{
	      // 	  xB = cB->col*ptchc - halfSensorY; // straight
	      // 	  yB = cB->row*ptchr - halfSensorX; // PCB
	      // 	}
	    
	      double xB = xcoordinate(1, cB, 0, ptchc, ptchr);
	      double yB = ycoordinate(1, cB, 0, ptchc, ptchr);

	      hxxAB->Fill( xAr, xB );
	      hyyAB->Fill( yAr, yB );
	      
	      double dx = xAr - xB;
	      double dy = yAr - yB;
	      
	      if( cA->q > qL  && cA->q < qR && cB->q > qLB && cB->q < qRB && cA->iso && cB->iso )
		{
		  
		  hdxAB->Fill( dx );
		  if(PRINT)   cout << hdxAB->GetTitle() << " entries " << hdxAB->GetEntries() << endl;
		  
		  hdyAB->Fill( dy );
		  dxvsxAB->Fill( xB, dx );
		  dxvsyAB->Fill( yB, dx );
		
		}
	  
	      hdxvsev->Fill( iev, dx );
	    
	      if( fabs( dx ) < straightTracks * beamDivergenceScaled + 0.020 && fabs( dy ) < straightTracks * beamDivergenceScaled + 0.100 )
		++nm;
	      
	    } // clusters
	  
	} // cl
      
      nmvsevAB->Fill( iev, nm );
      
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // B-C cluster correlations:
      
      nm = 0;
      
      for( vector<cluster>::iterator cB = vclB.begin(); cB != vclB.end(); ++cB )
	{
	
	  double xB = cB->row*ptchr - halfSensorX;
	  double yB = cB->col*ptchc - halfSensorY;
	  if( run == 431 )
	    {
	      xB = cB->row*ptchr - halfSensorX; // rot90
	      yB = cB->col*ptchc - halfSensorY; // PCB
	    }
	  if( fifty )
	    {
	      xB = cB->col*ptchc - halfSensorY; // straight
	      yB = cB->row*ptchr - halfSensorX; // PCB
	    }
	  
	  hxB->Fill( xB );
	  hyB->Fill( yB );
	  if( cB->iso )
	    {
	      hxBi->Fill( xB );
	      hyBi->Fill( yB );
	      hclqBi->Fill( cB->q );
	    }

	  for( vector<cluster>::iterator cC = vclC.begin(); cC != vclC.end(); ++cC )
	    {
	  
	      double xC = cC->row*ptchr - halfSensorX - alignxC; // rot90 Dreimaster
	      double yC = cC->col*ptchc - halfSensorY - alignyC; // down
	      if( fifty )
		{
		  xC = cC->col*ptchc - halfSensorY - alignxC; // straight
		  yC = cC->row*ptchr - halfSensorX - alignyC; // PCB
		}

	      double xCr = xC*cfC - yC*sfC;
	      double yCr = xC*sfC + yC*cfC;
	      
	      hxxCB->Fill( xB, xCr );
	      hyyCB->Fill( yB, yCr );
	      
	      double dx = xCr - xB;
	      double dy = yCr - yB;
	  
	  if( cC->q > qL  && cC->q < qR && cB->q > qLB && cB->q < qRB && cC->iso && cB->iso )
	    {
	    
	      hdxCB->Fill( dx );
	      hdyCB->Fill( dy );
	      dxvsxCB->Fill( xB, dx );
	      dxvsyCB->Fill( yB, dx );
	      
	    }
	  
	  if( fabs( dx ) < straightTracks * beamDivergenceScaled + 0.020 && fabs( dy ) < straightTracks * beamDivergenceScaled + 0.100 )
	    ++nm;
	  
	    } // clusters B
	  
	} // cl C
      
      nmvsevCB->Fill( iev, nm );
      
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // A-C cluster correlations:
      
      nm = 0;
      
      for( vector<cluster>::iterator cC = vclC.begin(); cC != vclC.end(); ++cC )
	{
	
	  double xC = cC->row*ptchr - halfSensorX - alignxC;
	  double yC = cC->col*ptchc - halfSensorY - alignyC;
	  if( fifty )
	    {
	      xC = cC->col*ptchc - halfSensorY - alignxC; // straight
	      yC = cC->row*ptchr - halfSensorX - alignyC; // PCB
	    }

	  double xCr = xC*cfC - yC*sfC;
	  double yCr = xC*sfC + yC*cfC;
	
	  hxC->Fill( xCr );
	  hyC->Fill( yCr );
	  if( cC->iso )
	    {
	      hxCi->Fill( xCr );
	      hyCi->Fill( yCr );
	      hclqCi->Fill( cC->q );
	    }
      
	  double etaC = -2;
	  if( cC->size == 2 )
	    {
	      double q0 = cC->vpix[0].q;
	      double q1 = cC->vpix[1].q;
	      etaC = (q1-q0)/(q1+q0);
	    }
	
	  for( vector<cluster>::iterator cA = vclA.begin(); cA != vclA.end(); ++cA )
	    {
	
	      double xA = cA->row*ptchr - halfSensorX - alignxA; // rot90 Dreimaster
	      double yA = cA->col*ptchc - halfSensorY - alignyA; // down
	      if( fifty )
		{
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
	      hdxCA->Fill( dxCA );
	      hdyCA->Fill( dyCA );


	      if(PRINT) cout << " tracks condition " << fabs( dxCA ) << " vs "  << straightTracks * beamDivergenceScaled + 0.02 << endl;
	      if( fabs( dxCA ) > straightTracks * beamDivergenceScaled + 0.02 ) continue; // includes beam divergence: +-5 sigma
	      if(PRINT) cout << " dyCA " << dyCA << endl;
	      hdyCAc->Fill( dyCA );
	      dyvsyCA->Fill( yAr, dyCA );
	      
	      if( fabs( dyCA ) > straightTracks * beamDivergenceScaled + 0.1 ) continue; // [mm]
	      
	      ++nm;
	  
	      dxvsyCA->Fill( yAr, dxCA );
	      
	      dxvsxCA->Fill( xAr, dxCA ); // linear trend in run 392, 403: acceptance and beam divergence ?
	      
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
	  
	      for( vector<cluster>::iterator cB = vclB.begin(); cB != vclB.end(); ++cB )
		{
	    
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
		  for( int ipx = 0; ipx < cB->size; ++ipx )
		    {
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

	  hdx3->Fill( dx3 );
	  hdy3->Fill( dy3 );

	  if( fabs( dy3 ) < straightTracks * beamDivergenceScaled + 0.05 ) { // cut on y, look at x, see madx3vsy

	    hdx3c->Fill( dx3 );

	    if( (cB->q <= qLB || cB->q >= qRB) && ( cA->q <= qL || cA->q >= qR) && ( cC->q <= qL || cC->q >= qR))
	      hdx3nocq3->Fill( dx3);

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
	      hdx3ci->Fill( dx3 );
 
	    if( cA->iso && cC->iso )
	      hdx3cii->Fill( dx3 );

	    if( cA->iso && cB->iso && cC->iso ) {

	      hdx3ciii->Fill( dx3 );

	      if( cB->size == 1 )
		hdx3c1->Fill( dx3 ); // r447 4.4
	      if( cB->size == 2 )
		hdx3c2->Fill( dx3 ); // r447 4.4
	      if( cB->size == 3 )
		hdx3c3->Fill( dx3 ); // r447 4.8
	      if( cB->size == 4 )
		hdx3c4->Fill( dx3 ); // r447 6.5
	      if( cB->size == 5 )
		hdx3c5->Fill( dx3 ); // r447 16.6
	      if( cB->size == 6 )
		hdx3c6->Fill( dx3 ); // r447 24.5
	      if( cB->size > 6 )
		hdx3c7->Fill( dx3 ); // r447 39.9

	      if( xB < 0 )
		hdx3m->Fill( dx3 );
	      else
		hdx3p->Fill( dx3 );

	      if( fabs( dxCA ) < straightTracks * beamDivergenceScaled ) { // track angle
		hdx3ct->Fill( dx3 );
		madx3vsq->Fill( cB->q, fabs(dx3) );
		madx3vsn->Fill( cB->size, fabs(dx3) );
	      }

	    } // iso

	    
	    if( cB->q > qLB && cB->q < qRB ) {

	      hdx3cq->Fill( dx3 );

	      if( cA->iso && cB->iso && cC->iso )
		hdx3cqi->Fill( dx3 );
	      
	      if( cA->q > qL && cA->q < qR &&
		  cC->q > qL && cC->q < qR ) {

		hdx3cq3->Fill( dx3 );

		if( cA->iso && cB->iso && cC->iso )
		  hdx3cq3i->Fill( dx3 );

		dx3vsev->Fill( iev, dx3 );

		dx3vsx->Fill( xB, dx3 ); // turn
		dx3vsy->Fill( yB, dx3 ); // rot
		dx3vsxm->Fill( xmod*1E3, dx3 );

		madx3vsdx->Fill( dxCA*1E3, fabs(dx3) ); // dxCA

		if( fabs( dxCA ) < straightTracks * beamDivergenceScaled ) { // track angle

		  hdx3cq3t->Fill( dx3 ); // 447 4.27 um

		  madx3vsx->Fill( xB, fabs(dx3) );
		  madx3vsy->Fill( yB, fabs(dx3) );
		  madx3vsxm->Fill( xmod*1E3, fabs(dx3) );
		  if( cB->size == 2 ) {
		    etavsxmB3->Fill( xmod*1E3, etaB ); // sine
		    madx3vseta->Fill( etaB, fabs(dx3) ); // flat
		    hdx3cq3t2->Fill( dx3 ); // 447 4.25 um
		  }

		} // angle

	      } // Qa, qC

	    } // qB

	  } // cut dy

	  if( fabs( dx3 ) < 0.07 && // hit on track
	      fabs( dy3 ) < 0.15 &&
	      cA->iso && cB->iso && cC->iso) {

	    hclmapB3->Fill( cB->col, cB->row );

	    hxA3->Fill( xAr );
	    hyA3->Fill( yAr );
	    hxB3->Fill( xB  );
	    hyB3->Fill( yB  );
	    hxC3->Fill( xCr );
	    hyC3->Fill( yCr );

	    hclszA3->Fill( cA->size );
	    hclszB3->Fill( cB->size );
	    hncolB3->Fill( ncolB );
	    hnrowB3->Fill( nrowB );
	    hclszC3->Fill( cC->size );

	    hclphA3->Fill( cA->sum );
	    hclphB3->Fill( cB->sum );
	    hclphC3->Fill( cC->sum );

	    hclqA3->Fill( cA->q );
	    hclqB3->Fill( cB->q );
	    hclqB3i->Fill( cB->q );
	    hclqC3->Fill( cC->q );
	    if( cB->size < 4 )
	      hclqB3n->Fill( cB->q );

	    if( fifty )
	      nrowvsxmB3->Fill( xmod*1E3, ncolB );
	    else
	      nrowvsxmB3->Fill( xmod*1E3, nrowB );

	    clqvsxmB3->Fill( xmod*1E3, cB->q );

	    if( cA->size == 2 )
	      hetaA3->Fill( etaA );

	    if( cB->size == 2 )
	      hetaB3->Fill( etaB );

	    if( cC->size == 2 )
	      hetaC3->Fill( etaC );

	    for( int ipx = 0; ipx < cA->size; ++ipx ) {
	      hpxpA3->Fill( cA->vpix[ipx].ph );
	      hpxqA3->Fill( cA->vpix[ipx].q );
	    }
	    for( int ipx = 0; ipx < cB->size; ++ipx ) {
	      hpxpB3->Fill( cB->vpix[ipx].ph );
	      hpxqB3->Fill( cB->vpix[ipx].q );
	    }
	    for( int ipx = 0; ipx < cC->size; ++ipx ) {
	      hpxpC3->Fill( cC->vpix[ipx].ph );
	      hpxqC3->Fill( cC->vpix[ipx].q );
	    }

	    if( cB->size == 2 ) {
	      hpxq1stB3->Fill( cB->vpix[0].q );
	      hpxq2ndB3->Fill( cB->vpix[1].q ); // identical
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
	  effvsx->Fill( xavg, eff[50] );

	if( xavg > -3.3 && xavg < 3.6 )
	  effvsy->Fill( yavg, eff[50] );

	if( xavg > -3.3 && xavg < 3.6 &&
	    yavg > -3.7 && yavg < 3.5 ) {

	  for( int iw = 1; iw < 999; ++iw )
	    effvsdxy->Fill( iw*0.010+0.005, eff[iw] );

	  effvsxm->Fill( xmod*1E3, eff[50] ); // bias dot
	  effvsev->Fill( iev, eff[50] );
	  effvsiev->Fill( iev%200, eff[50] );
	  //effvsmpxA->Fill( pbA.size(), eff[50] );
	  effvsqA->Fill( cA->q, eff[50] );
	  effvstxy->Fill( dxyCA, eff[50] ); // flat

	} // fiducial x, y

      } // clusters A

    } // cl C

    nmvsevCA->Fill( iev, nm );

    // task: track correlations, intersects

    
    } // events

  cout << endl;
  if(PRINT) cout << " after big loop on events!!!************************" << endl;
  if(PRINT)cout << hdxAB->GetTitle() << " entries " << hdxAB->GetEntries() << endl;
  if(PRINT)cout << hdt->GetTitle() << " entries " << hdt->GetEntries() << endl;


  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s3 = ts.tv_sec; // seconds since 1.1.1970
  long f3 = ts.tv_nsec; // nanoseconds
  zeit2 += s3 - s2 + ( f3 - f2 ) * 1e-9; // read and cluster
  cout << "time " << zeit2 << " s" << endl;
  if(PRINT)cout << hdxAB->GetTitle() << " entries " << hdxAB->GetEntries() << endl;
  if(PRINT)cout << hdt->GetTitle() << " entries " << hdt->GetEntries() << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s9 = ts.tv_sec; // seconds since 1.1.1970
  long f9 = ts.tv_nsec; // nanoseconds

  cout << "full time " << s9 - s0 + ( f9 - f0 ) * 1e-9 << " s"
       << " (read and cluster " << zeit1 << " s, tracking " << zeit2 << " s)"
       << endl;
  if(PRINT)cout << hdxAB->GetTitle() << " entries " << hdxAB->GetEntries() << endl;
  if(PRINT)cout << hdt->GetTitle() << " entries " << hdt->GetEntries() << endl;

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
  if(PRINT) cout << " going to write the hist file" << endl;
  histoFile->Write();
  if(PRINT) cout << " wrote the hist file" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // alignment fits:

  double newalignxA = alignxA;
  if(PRINT) cout << " newalignxA " << newalignxA << " alignxA " << alignxA << endl;

  cout << hdxAB->GetTitle() << " entries " << hdxAB->GetEntries() << endl;
  cout << hdt->GetTitle() << " entries " << hdt->GetEntries() << endl;

  if( hdxAB->GetEntries() > 999 ) {

    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -10, 10 );
    double xpk = hdxAB->GetBinCenter( hdxAB->GetMaximumBin() );
    fgp0->SetParameter( 0, hdxAB->GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, xpk );
    fgp0->SetParameter( 2, hdxAB->GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, hdxAB->GetBinContent( hdxAB->FindBin(xpk-1) ) ); // BG
    hdxAB->Fit( "fgp0", "q", "", xpk-1, xpk+1 ); // fit range around peak
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

  cout << endl << hdyAB->GetTitle() << " entries " << hdyAB->GetEntries() << endl;
  if( hdyAB->GetEntries() > 999 ) {
    cout << "  y correction " << hdyAB->GetMean() << endl;
    newalignyA += hdyAB->GetMean();
  }
  else
    cout << "  not enough" << endl;

  // dxvsy -> -rot

  double newalignfA = alignfA;

  cout << endl << dxvsyAB->GetTitle() << " entries " << dxvsyAB->GetEntries() << endl;
  if( aligniteration > 0 && dxvsyAB->GetEntries() > 999 ) {
    dxvsyAB->Fit( "pol1", "q", "", -3, 3 );
    TF1 * fdxvsy = dxvsyAB->GetFunction( "pol1" );
    cout << "  extra rot " << fdxvsy->GetParameter(1) << endl;
    newalignfA += fdxvsy->GetParameter(1);
  }
  else
    cout << "  not enough" << endl;

  // C:

  double newalignxC = alignxC;

  cout << endl << hdxCB->GetTitle() << " entries " << hdxCB->GetEntries() << endl;
  if( hdxCB->GetEntries() > 999 ) {

    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -10, 10 );
    double xpk = hdxCB->GetBinCenter( hdxCB->GetMaximumBin() );
    fgp0->SetParameter( 0, hdxCB->GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, xpk );
    fgp0->SetParameter( 2, hdxCB->GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, hdxCB->GetBinContent( hdxCB->FindBin(xpk-1) ) ); // BG
    hdxCB->Fit( "fgp0", "q", "", xpk-1, xpk+1 ); // fit range around peak
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

  cout << endl << hdyCB->GetTitle() << " entries " << hdyCB->GetEntries() << endl;
  if( hdyCB->GetEntries() > 999 ) {
    cout << "  y correction " << hdyCB->GetMean() << endl;
    newalignyC += hdyCB->GetMean();
  }
  else
    cout << "  not enough" << endl;

  // dxvsy -> -rot

  double newalignfC = alignfC;

  cout << endl << dxvsyCB->GetTitle() << " entries " << dxvsyCB->GetEntries() << endl;
  if( aligniteration > 0 && dxvsyCB->GetEntries() > 999 ) {
    dxvsyCB->Fit( "pol1", "q", "", -3, 3 );
    TF1 * fdxvsy = dxvsyCB->GetFunction( "pol1" );
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
  histoFile->Close();
  if(PRINT) cout << " closed the hist file" << endl;

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

  while( seed < pb.size() )
    {

      // start a new cluster:      
      cluster c;
      c.vpix.push_back( pb[seed] );
      gone[seed] = 1;
      
      // let it grow as much as possible:

      int growing;
      do{
	growing = 0;
	for( unsigned i = 0; i < pb.size(); ++i )
	  {
	    if( !gone[i] )
	      { // unused pixel
		for( unsigned int p = 0; p < c.vpix.size(); ++p )
		  { // vpix in cluster so far
		    int dr = c.vpix.at(p).row - pb[i].row;
		    int dc = c.vpix.at(p).col - pb[i].col;
		    if( ( dr >= -fCluCut ) && ( dr <= fCluCut ) && ( dc >= -fCluCut ) && ( dc <= fCluCut ) )//if this condition is satisfied the new pixel c.vpix.at(p) and the old one pb[i] are adjacent
		      {
			c.vpix.push_back(pb[i]); //add adjacent pixel to cluster
			gone[i] = 1; //pixel has been used
			growing = 1; //cluster is growing
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
      
      for( vector<pixel>::iterator p = c.vpix.begin();  p != c.vpix.end();  ++p )
	{
	  double ph = p->ph;
	  c.sum += ph; //cluster charge in ph units
	  double q = p->q;
	  c.q += q; //cluster charge
	  //c.col += (*p).col*ph;
	  //c.row += (*p).row*ph;
	  c.col += (*p).col*q; //pixel charge is a weight on the pixel position
	  c.row += (*p).row*q;
	}

      c.size = c.vpix.size();

      if(PRINT) cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

      //cluster coordinates is the average of the pixel coordianates 
      c.col /= c.q;
      c.row /= c.q;
      
      v.push_back(c); // add cluster to vector
      
      // look for a new seed = used pixel:

      while( ( ++seed < pb.size() ) && gone[seed] );
      
    }// while over seeds

  // nothing left,  return clusters

  delete gone;
  return v;
}//getClus


list < vector < cluster > > oneplane( int plane, string runnum, unsigned Nev, bool fifty, double Tsunami, double dphcut )
{
  int run = stoi( runnum );

  list < vector < cluster > > evlist;
  string datadir;
  string Xfile;
  string runpath = "roi000";
  if(run > 999)
    runpath = "roi00";
  if(plane == A) datadir="a/";
  if(plane == B) datadir="b/";
  if(plane == C) datadir="c/";
  
  Xfile = datapath+datadir+runpath + runnum + ".txt";
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

  //empty and filled event flag
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
    iss >> evtime; // from run 456

    vector <pixel> pb; // for clustering
    vector <cluster> vcl;
    evInfo evinf;
    evinf.evtime = evtime;
    evinf.filled = filled;
    evinf.skip = 0;
    prevtime = evtime;

    //do the analysis only for filled events
    if( filled == F )
      {
	if(PRINT)  cout << "Analayzing filled data" << endl; 
	string roi;
	getline( Xstream, roi );
	if(PRINT)  cout << "roi: " << roi << endl;
	if(PRINT) cout << "roi size: " << roi.size() << endl;
	size_t start = 0;
	size_t gap = 0;
	unsigned ng = 0; // 3 entries per pixel (col, row and ph)
	string BLANK{" "};
	while( gap < roi.size()-1 )
	  { // data have trailing blank
	    gap = roi.find( BLANK, start );
	    start = gap + BLANK.size();
	    ++ng;
	  }
	hnpx[plane].Fill( ng/3 );
	
	vector <pixel> vpx;
	
	if( ng/3 < 400 )
	  { //if there are more than 400 hits, it is a noisy event and takes a lot!
	  
	    vpx.reserve(ng/3); //pixels in the event
	  
	    size_t start = 0;
	    size_t gap = 0;
	    while( gap < roi.size()-1 )
	      { // data have trailing blank
	    
		//getting pixel column
		gap = roi.find( BLANK, start );
		string s1( roi.substr( start, gap - start ) );
		if(PRINT) cout << " " << s1 << "(" << gap << ")"<<endl;
		int col = atoi( s1.c_str() ); // 4% faster
		start = gap + BLANK.size();
	    
		//getting pixel row
		gap = roi.find( BLANK, start );
		string s2( roi.substr( start, gap - start ) );
		if(PRINT) cout << " " << s2 << "(" << gap << ")";
		int row = atoi( s2.c_str() );
		start = gap + BLANK.size();

		//getting pixel puls height
		gap = roi.find( BLANK, start );
		string s3( roi.substr( start, gap - start ) );
		if(PRINT) cout << " " << s1 << "(" << gap << ")";
		double ph = atof(s3.c_str());
		start = gap + BLANK.size();

		hph[plane].Fill( ph );

		pixel px { col, row, ph, ph };
		vpx.push_back(px); // comment out = no clustering
	    
	      }//while 
	  } // size (less than 400 hits!)
      else
	{
	  evinf.skip = 1;
	  cout << " (" << iev << ": B ROI " << ng/3 << " skip)";
	}
	
	// column-wise common mode correction:
	// 4 = central pixel, 1 lower rows, 7 upper rows
	double phprev = 0; //pulse height of previuos pixel
	double dphprev = 0; //ph difference

	for( unsigned ipx = 0; ipx < vpx.size(); ++ipx )
	  {

	    int col4 = vpx[ipx].col;
	    int row4 = vpx[ipx].row;
	    double ph4 = vpx[ipx].ph;
	    
	    int row1 = row4;
	    int row7 = row4;
	    double ph1 = ph4;
	    double ph7 = ph4;
	    
	    for( unsigned jpx = 0; jpx < vpx.size(); ++jpx )
	      {
		
		if( jpx == ipx ) continue;
		if( vpx[jpx].col != col4 ) continue; // want same column
		
		int jrow = vpx[jpx].row;
	      
		if( jrow < row1 )
		{
		  row1 = jrow;
		  ph1 = vpx[jpx].ph;
		}

		if( jrow > row7 )
		  {
		    row7 = jrow;
		    ph7 = vpx[jpx].ph;
		  }
	      
	      } // jpx

	    if( row4 == row1 ) //this is possible only when jpx=ipx => nothing happened in the previous loop => ph1 is still storing the ph4 of the previous ipx
	      {
		phprev = ph1;   
		continue; // Randpixel
	      }
	    if( row4 == row7 ) continue;

	    phvsprev[plane].Fill( phprev, ph4 ); //looking at this plot I can see if there are correlections = "tsunami effect"

	    if(PRINT) cout << " ph4 before " << ph4;
	    ph4 -= Tsunami*phprev; // Tsunami
	    if(PRINT) cout << " ph4 after " << ph4 << endl;
	    
	    phprev = vpx[ipx].ph; // original ph4

	    double dph;
	    
	    //find closest row to define the difference
	    if( row4 - row1 < row7 - row4 )
	      dph = ph4 - ph1;
	    else
	      dph = ph4 - ph7;
 
	    hdph[plane].Fill( dph ); // sig 2.7

	    dphvsprev[plane].Fill( dphprev, dph );
	    dphprev = dph;


	    if( dph > dphcut ) 
	      {

		hpxmap[plane]->Fill( col4, row4 );

		pixel px;
		
		if( fifty )
		  {
		    px.col = col4;
		    px.row = row4;
		  }
		else
		  {
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

	      }//cut on dph
	
	  } // ipx
	
	//if( ldb ) cout << " roi " << vpx.size() << ", hits " << pb.size() << endl;
	
      } // filled

    // clustering:

    hnht[plane].Fill( pb.size() );

    if( pb.size() > 50 ) pb.clear(); // speed because it deletes events with more than 50 hits

    vcl = getClus(pb); // this function is doing the clustering

    hncl[plane].Fill( vcl.size() );
    if( vcl.size() )
      if( ldb ) cout << "  clusters " << vcl.size() << endl;

    for( unsigned icl = 0; icl < vcl.size(); ++icl )
      {
	//filling some hists about clusters
	hclmap[plane]->Fill( vcl[icl].col, vcl[icl].row );
	hclph[plane].Fill( vcl[icl].sum );
	hclq[plane].Fill( vcl[icl].q );
	if( vcl[icl].sum > 55 ) hclsz[plane].Fill( vcl[icl].size );//55 is historical adc, for noise suppression

      // cluster isolation:

      for( unsigned jcl = icl+1; jcl < vcl.size(); ++jcl )
	{

	  bool done = 0;
	  
	  for( unsigned ipx = 0; ipx < vcl[icl].vpix.size(); ++ipx ) {

	  for( unsigned jpx = 0; jpx < vcl[jcl].vpix.size(); ++jpx )
	    if( fabs( vcl[icl].vpix[ipx].col - vcl[jcl].vpix[jpx].col ) < 3 &&
		fabs( vcl[icl].vpix[ipx].row - vcl[jcl].vpix[jpx].row ) < 3 )
	      {
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

double xcoordinate(int plane, vector<cluster>::iterator c, double align, double pitchc, double pitchr)
{
  double variable;
  
  
  variable = c->row*pitchr - halfSensorX - align;
  // if( run == 431 && plane == 1 )
  //   variable = c->row*ptchr - halfSensorX; // rot90
  if( fifty )
    variable = c->col*pitchc - halfSensorY - align; // straight
  
  return variable;
}
      
double ycoordinate(int plane, vector<cluster>::iterator c, double align, double pitchc, double pitchr)
{
  double variable;
  
   variable =  c->col*pitchc - halfSensorY - align;
  // if( run == 431 && plane == 1  )
  //   variable =  c->col*ptchc - halfSensorY; // PCB
  if( fifty )
    variable = c->row*pitchr - halfSensorX - align; // PCB

  return variable;
}

void bookHists()
{

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

  hdt = new TH1I( "dt", "time between events;log_{10}(#Deltat [s]);events", 100, -4, 1 );
  hddtAB = new TH1I( "ddtAB", "dtA - dtB;A-B #delta#Deltat [clocks];events", 100, -1000, 1000 );
    hddtCB = new TH1I( "ddtCB", "dtC - dtB;C-B #delta#Deltat [clocks];events", 100, -1000, 1000 );
    hddtCA = new TH1I( "ddtCA", "dtC - dtA;C-A #delta#Deltat [clocks];events", 100, -1000, 1000 );
    ddtvsdtAB = new  TProfile( "ddtvsdtAB",		      "A-B time lag vs intervall;log_{10} = new  (#Deltat [s]);<#delta#Deltat> [clocks]",		      100, -4, 1,-1e99, 1e99 );

    // correlations:

     hxA = new TH1I( "xA", "x A;x [mm];clusters A", 100, -5, 5 );
     hyA = new  TH1I( "yA", "y A;y [mm];clusters A", 100, -5, 5 ); 
   hxAi = new  TH1I( "xAi", "x A isolated;x [mm];isolated clusters A", 100, -5, 5 );
   hyAi = new  TH1I( "yAi", "y A isolated;y [mm];isolated clusters A", 100, -5, 5 ); 
   hclqAi = new  TH1I( "clqAi", "A isolated cluster charge;cluster charge [ke];A isolatewd clusters",
	       100, 0, 50 );

   hxxAB = new TH2I( "xxAB", "B vs A;row A;row B;clusters", 320, -4, 4, 320, -4, 4 );
   hyyAB = new TH2I( "yyAB", "B vs A;col A;col B;clusters",  80, -4, 4,  80, -4, 4 );

   hdxAB = new  TH1I( "dxAB", "Bx-Ax;x-x [mm];cluster pairs", 800, -2, 2 );
   hdyAB = new  TH1I( "dyAB", "By-Ay;y-y [mm];cluster pairs", 400, -2, 2 );
   dxvsxAB = new  TProfile( "dxvsxAB", "dx vs x A-B;x [mm];<dx> [mm]", 320, -4, 4, -f, f );
 dxvsyAB = new  TProfile( "dxvsyAB", "dx vs y A-B;y [mm];<dx> [mm]",  80, -4, 4, -f, f );

  hdxvsev = new    TH2I( "dxvsev", "Bx-Ax vs events;events;#Deltax [px];clusters",	  100, 0, 10000, 100, -f, f );

   nmvsevAB = new  TProfile( "nmvsevAB", "AB matches vs time;time [events];AB matches",		     3100, 0, 3100*1000, -1, 99 );

   hxB = new  TH1I( "xB", "x B;x [mm];clusters B", 100, -5, 5 );
 hyB = new  TH1I( "yB", "y B;y [mm];clusters B", 100, -5, 5 ); 
   hxBi = new  TH1I( "xBi", "x B isolated;x [mm];isolated clusters B", 100, -5, 5 );
 hyBi = new  TH1I( "yBi", "y B isolated;y [mm];isolated clusters B", 100, -5, 5 ); 
   hclqBi = new  TH1I( "clqBi", "B isolated cluster charge;cluster charge [ke];B isolatewd clusters",100, 0, 50 );

  hxxCB = new TH2I( "xxCB", "C vs B;row B;row C;clusters", 320, -4, 4, 320, -4, 4 );
  hyyCB = new TH2I( "yyCB", "C vs B;col B;col C;clusters",  80, -4, 4,  80, -4, 4 );

hdxCB = new  TH1I( "dxCB", "Cx-Bx;x-x [mm];cluster pairs", 800, -2, 2 );
 hdyCB = new  TH1I( "dyCB", "Cy-By;y-y [mm];cluster pairs", 400, -2, 2 );
   dxvsxCB = new  TProfile( "dxvsxCB", "dx vs x C-B;x [mm];<dx> [mm]", 320, -4, 4, -f, f );
 dxvsyCB = new  TProfile( "dxvsyCB", "dx vs y C-B;y [mm];<dx> [mm]",  80, -4, 4, -f, f );
 nmvsevCB = new  TProfile( "nmvsevCB", "CB matches vs time;time [events];CB matches",
		     3100, 0, 3100*1000, -1, 99 );

  // triplets:

   hxC = new  TH1I( "xC", "x C;x [mm];clusters C", 100, -5, 5 );
 hyC = new  TH1I( "yC", "y C;y [mm];clusters C", 100, -5, 5 ); 
 hxCi = new  TH1I( "xCi", "x C isolated;x [mm];isolated clusters C", 100, -5, 5 );
 hyCi = new  TH1I( "yCi", "y C isolated;y [mm];isolated clusters C", 100, -5, 5 ); 
 hclqCi = new  TH1I( "clqCi", "C isolated cluster charge;cluster charge [ke];C isolatewd clusters",
	       100, 0, 50 );

  hxxCA = new TH2I( "xxCA", "C vs A;row A;row C;clusters", 320, -4, 4, 320, -4, 4 );
  hyyCA = new TH2I( "yyCA", "C vs A;col A;col C;clusters",  80, -4, 4,  80, -4, 4 );

 hdxCA = new  TH1I( "dxCA", "Cx-Ax;x-x [mm];cluster pairs", 400, -1, 1 );
 hdyCA = new  TH1I( "dyCA", "Cy-Ay;y-y [mm];cluster pairs", 200, -1, 1 );

 dxvsxCA = new  TProfile( "dxvsxCA", "dx vs x C-A;x [mm];<dx> [mm]", 320, -4, 4, -f, f );
 dxvsyCA = new  TProfile( "dxvsyCA", "dx vs y C-A;y [mm];<dx> [mm]",  80, -4, 4, -f, f );
  hdyCAc = new  TH1I( "dyCAc", "Cy-Ay, cut dx;y-y [mm];cluster pairs", 200, -1, 1 );
 dyvsyCA = new  TProfile( "dyvsyCA", "dy vs y C-A;y [mm];<dy> [mm]",  80, -4, 4, -f, f );
 nmvsevCA = new  TProfile( "nmvsevCA", "CA matches vs time;time [events];CA matches",
		     3100, 0, 3100*1000, -1, 99 );

  hdx3 = new  TH1I( "dx3", "triplet dx;dx [mm];triplets", 500, -0.5, 0.5 );
  hdy3 = new  TH1I( "dy3", "triplet dy;dy [mm];triplets", 200, -1, 1 );

  hdx3c = new  TH1I( "dx3c", "triplet dx, cut dy;dx [mm];triplets", 500, -0.25, 0.25 );
  hdx3ci = new  TH1I( "dx3ci", "triplet dx, cut dy, isolated;dx [mm];isolated triplets", 500, -0.25, 0.25 );
  hdx3cii = new  TH1I( "dx3cii", "triplet dx, cut dy, isolated;dx [mm];isolated triplets", 500, -0.25, 0.25 );
  hdx3ciii = new  TH1I( "dx3ciii", "triplet dx, cut dy, isolated;dx [mm];isolated triplets", 500, -0.25, 0.25 );

  hdx3c1 = new  TH1I( "dx3c1", "triplet dx, cut dy, npx 1;dx [mm];triplets, B npx 1",
	       500, -0.25, 0.25 );
  hdx3c2 = new  TH1I( "dx3c2", "triplet dx, cut dy, npx 2;dx [mm];triplets, B npx 2",
	       500, -0.25, 0.25 );
  hdx3c3 = new TH1I( "dx3c3", "triplet dx, cut dy, npx 3;dx [mm];triplets, B npx 3",
	       500, -0.25, 0.25 );
  hdx3c4 = new TH1I( "dx3c4", "triplet dx, cut dy, npx 4;dx [mm];triplets, B npx 4",
	       500, -0.25, 0.25 );
  hdx3c5 = new TH1I( "dx3c5", "triplet dx, cut dy, npx 5;dx [mm];triplets, B npx 5",
	       500, -0.25, 0.25 );
  hdx3c6 = new TH1I( "dx3c6", "triplet dx, cut dy, npx 6;dx [mm];triplets, B npx 6",
	       500, -0.25, 0.25 );
  hdx3c7 = new TH1I( "dx3c7", "triplet dx, cut dy, npx > 6;dx [mm];triplets, B npx > 6",
	       500, -0.25, 0.25 );

  hdx3m = new TH1I( "dx3m", "triplet dx, x < 0;dx [mm];triplets", 500, -0.5, 0.5 );
  hdx3p = new TH1I( "dx3p", "triplet dx, x > 0;dx [mm];triplets", 500, -0.5, 0.5 );
  hdx3ct = new TH1I( "dx3ct", "triplet dx, cut dy, tx;dx [mm];triplets",
	       500, -0.25, 0.25 );
 madx3vsq = new TProfile( "madx3vsq", "MAD = new  (dx) vs Q;B cluster charge [ke];MAD dx [mm]",
		     100, 0, 100, 0, 0.1 );
 madx3vsn = new TProfile( "madx3vsn", "MAD = new  (dx) vs cluster size;B cluster size [pixels];MAD dx [mm]",
		     20, 0.5, 20.5, 0, 0.1 );

  hdx3cq = new TH1I( "dx3cq", "triplet dx, Landau peak;dx [mm];Landau peak triplets",
	       500, -0.25, 0.25 );
  hdx3cqi = new TH1I( "dx3cqi", "triplet dx, Landau peak, isolated;dx [mm];isolated Landau peak triplets",
		500, -0.25, 0.25 );
  hdx3cq3 = new TH1I( "dx3cq3", "triplet dx, 3 Landau peak;dx [mm];Landau peak triplets",
		500, -0.25, 0.25 );
  hdx3nocq3 = new TH1I( "dx3nocq3", "triplet dx, no 3 Landau peak;dx [mm];no Landau peak triplets",
		500, -0.25, 0.25 );
  hdx3cq3i = new TH1I( "dx3cq3i", "triplet dx, 3 Landau peak, isolated;dx [mm];isolated Landau peak triplets",
		 500, -0.25, 0.25 );

 dx3vsev = new TProfile( "dx3vsev", "dx3 vs time;trigger;<dx3> [mm]",
		    310, 0, 3100*1000, -0.5, 0.5 );

 dx3vsx = new TProfile( "dx3vsx", "dx vs x;x [mm];<dx3> [mm]", 320, -4, 4, -0.5, 0.5 );
 dx3vsy = new TProfile( "dx3vsy", "dx vs y;y [mm];<dx3> [mm]",  80, -4, 4, -0.5, 0.5 );
 dx3vsxm = new TProfile( "dx3vsxm", "dx vs x mod 50 um;x mod 50 [#mum];<dx3> [mm]",
		    50, 0, 50, -0.5, 0.5 );

madx3vsdx = new TProfile( "madx3vsdx", "MAD = new  (dx3) vs dx C-A;C-A dx [#mum];MAD dx3 [mm]",
		      100, -100, 100, 0, 0.1 );

  hdx3cq3t = new TH1I( "dx3cq3t",
		 "triplet dx, 3 Landau peak, forward;dx [mm];Landau peak forward triplets",
		 500, -0.25, 0.25 );
 madx3vsx = new TProfile( "madx3vsx", "MAD = new  (dx3) vs x;x [mm];MAD dx3 [mm]", 320, -4, 4, 0, 0.1 );
 madx3vsy = new TProfile( "madx3vsy", "MAD = new  (dx3) vs y;y [mm];MAD dx3 [mm]",  80, -4, 4, 0, 0.1 );
 madx3vsxm = new TProfile( "madx3vsxm", "MAD = new  (dx3) vs xmod;x mod 50 [#mum];MAD dx3 [mm]",		      50, 0, 50, 0, 0.1 );

 etavsxmB3 = new TProfile( "etavsxmB3", "eta vs xmod;x mod 50 [#mum];B <eta>",
		      50, 0, 50, -1.1, 1.1 );
 madx3vseta = new TProfile( "madx3vseta", "MAD = new  (dx3) vs eta;eta;MAD dx3 [mm]",
		       100, -1, 1, 0, 0.1 );
  hdx3cq3t2 = new TH1I( "dx3cq3t2",
		  "triplet dx, 3 Landau peak, forward, 2-px;dx [mm];Landau peak forward triplets",
		  500, -0.25, 0.25 );

  hclmapB3 = new  TH2I( "clmapB3", "linked cluster map B;col;row;B clusters on tracks",
	  80, 0, 80, 320, 0, 320 );

  hxA3 = new TH1I( "xA3", "x A linked;x [mm];A clusters on tracks", 100, -5, 5 );
  hyA3 = new TH1I( "yA3", "y A linked;y [mm];A clusters on tracks", 100, -5, 5 ); 
  hxB3 = new TH1I( "xB3", "x B linked;x [mm];B clusters on tracks", 100, -5, 5 );
  hyB3 = new TH1I( "yB3", "y B linked;y [mm];B clusters on tracks", 100, -5, 5 ); 
  hxC3 = new TH1I( "xC3", "x C linked;x [mm];C clusters on tracks", 100, -5, 5 );
  hyC3 = new TH1I( "yC3", "y C linked;y [mm];C clusters on tracks", 100, -5, 5 ); 

  hclszA3 = new TH1I( "clszA3", "A cluster size on tracks;cluster size [pixels];Aclusters on tracks",
		40, 0.5, 40.5 );
  hclphA3 = new TH1I( "clphA3", "A cluster PH on tracks;cluster ph [ADC];A clusters on tracks",
		200, 0, 1000 );
  hclqA3 = new TH1I( "clqA3", "A cluster charge on tracks;cluster charge [ke];A clusters on tracks",
	       160, 0, 80 );

  hclszB3 = new TH1I( "clszB3", "B cluster size on tracks;cluster size [pixels];B clusters on tracks",
		40, 0.5, 40.5 );
  hncolB3 = new TH1I( "ncolB3", "B cluster size on tracks;cluster size [columns];B clusters on tracks",
		20, 0.5, 20.5 );
  hnrowB3 = new TH1I( "nrowB3", "B cluster size on tracks;cluster size [rows];B clusters on tracks",
		20, 0.5, 20.5 );
  hclphB3 = new TH1I( "clphB3", "B cluster PH on tracks;cluster ph [ADC];B clusters on tracks",
		200, 0, 1000 );
  hclqB3 = new TH1I( "clqB3", "B cluster charge on tracks;cluster charge [ke];B clusters on tracks",
	       160, 0, 80 );
  hclqB3i = new TH1I( "clqB3i", "B cluster charge on tracks;cluster charge [ke];B clusters on tracks",
	       80, 0, 20 );
  hclqB3n = new TH1I( "clqB3n",
		"B cluster charge on tracks, npx < 4;cluster charge [ke];B clusters on tracks, npx < 4",
		160, 0, 80 );

  hclszC3 = new TH1I( "clszC3", "C cluster size on tracks;cluster size [pixels];C clusters on tracks",
		40, 0.5, 40.5 );
  hclphC3 = new TH1I( "clphC3", "C cluster PH on tracks;cluster ph [ADC];C clusters on tracks",
		200, 0, 1000 );
  hclqC3 = new TH1I( "clqC3", "C cluster charge on tracks;cluster charge [ke];C clusters on tracks",
	       160, 0, 80 );

 nrowvsxmB3 = new TProfile( "nrowvsxmB3",
		       "B rows vs xmod;x mod 50 [#mum];<B cluster size [rows]>",
		       50, 0, 50, 0.5, 10.5 );
   clqvsxmB3 = new TProfile( "clqvsxmB3",
		      "B cluster charge vs xmod;x mod 50 [#mum];<B cluster charge [ke]>",
		      50, 0, 50, 0, 50 );

  hetaA3 = new TH1I( "etaA3", "A cluster eta;eta;A 2-pix clusters on tracks",
	       100, -1, 1 );
  hetaB3 = new TH1I( "etaB3", "B cluster eta;eta;B 2-pix clusters on tracks",
	       100, -1, 1 );
  hetaC3 = new TH1I( "etaC3", "C cluster eta;eta;C 2-pix clusters on tracks",
	       100, -1, 1 );

  hpxqA3 = new TH1I( "pxqA3", "A pixel charge;pixel charge [ke];A pixels on tracks",
	       100, 0, 20 );
  hpxqB3 = new TH1I( "pxqB3", "B pixel charge;pixel charge [ke];B pixels on tracks",
	       100, 0, 20 );
   hpxqC3 = new TH1I( "pxqC3", "C pixel charge;pixel charge [ke];C pixels on tracks",
	       100, 0, 20 );

   hpxpA3 = new TH1I( "pxpA3", "A pixel PH;pixel PH [ADC];A pixels on tracks", 250, 0, 500 );
   hpxpB3 = new TH1I( "pxpB3", "B pixel PH;pixel PH [ADC];B pixels on tracks", 250, 0, 500 );
   hpxpC3 = new TH1I( "pxpC3", "C pixel PH;pixel PH [ADC];C pixels on tracks", 250, 0, 500 );

  hpxq1stB3 = new TH1I( "pxq1stB3", "B 1st pixel charge;pixel charge [ke];B `st pixels on tracks",		  100, 0, 20 );
  hpxq2ndB3 = new TH1I( "pxq2ndB3", "B 2nd pixel charge;pixel charge [ke];B `st pixels on tracks",		  100, 0, 20 );

  effvsdxy = new TProfile( "effvsdxy",		     "DUT efficiency vs triplet dxy;xy match radius [mm];DUT efficiency",		     1000, 0, 10, -0.1, 1.1 );

effvsxy =
    new TProfile2D( "effvsxy",
		    "DUT efficiency map;x [mm];y[mm];DUT efficiency",
		    80, -4, 4, 80, -4, 4, -0.1, 1.1 );
 effvsx = new TProfile( "effvsx", "eff vs x;x [mm];DUT efficiency",
		   320, -4, 4, -0.1, 1.1 );
 effvsy = new TProfile( "effvsy", "eff vs y;y [mm];DUT efficiency",
		   80, -4, 4, -0.1, 1.1 );
 effvsxm = new TProfile( "effvsxm", "eff vs x mod 50;x mod 50 [#mum];DUT efficiency",
		    50, 0, 50, -0.1, 1.1 );

 effvsev = new TProfile( "effvsev", "eff vs time;trigger;DUT efficiency",
		    3100, 0, 3100*1000, -0.1, 1.1 );
 effvsiev = new TProfile( "effvsiev", "eff vs event;event mod 200;DUT efficiency",
		     100, -0.5, 199.5, -0.1, 1.1 );
effvsmpxA = new TProfile( "effvsmpxA", "eff vs occupancy A;occupancy A [pixels];DUT efficiency",
		      50, 0.5, 50.5, -0.1, 1.1 );
 effvsqA = new TProfile( "effvsqA", "eff vs charge A;cluster charge A [ke];DUT efficiency",
		    100, 0, 100, -0.1, 1.1 );
   effvstxy = new TProfile("effvstxy", "eff vs angle;dxy CA [mm];DUT efficiency",
		     100, 0, 0.2, -0.1, 1.1 );

}
