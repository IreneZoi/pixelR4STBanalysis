
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

bool PRINT = false;
bool DOALIGNMENT = false;


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
  double dx3corr;
  
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
    string DX3C( "dx3c" );
    
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
      else if( tag == DX3C )
	tokenizer >> dx3corr;

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
  TString fileName;
  fileName.Form( "/home/zoiirene/Output/drei-r%i_irene.root", run );
  TFile * histoFile = new TFile( fileName, "RECREATE" ); //Form( "/home/zoiirene/Output/drei-r%i_irene.root", run ), "RECREATE" );

  // book histos:
  if(PRINT) cout << "***** going to book hists ***********" << endl;
  bookHists();
  histoMap raw =  bookControlHists("raw",histoFile);
  histoMap dxCAcut =  bookControlHists("dxCAcut",histoFile);

  histoMap beforeCorrections =   bookControlHists("beforeCorrections",histoFile);

  if(PRINT) cout << "hists booking whit new method"<<endl;
  histoMap nocuts =   bookControlHists("nocuts",histoFile);
  if(PRINT) cout << " booked! "  <<endl;


  if(PRINT)
    {
      bool found = false;

      auto it = nocuts.find("xA_nocuts");
      if(it != nocuts.end())
	{
	  found = true;
	  std::cout<<it->first<<" :: "<<it->second->GetEntries()<<std::endl;
	  it++;
	}
      
      cout << "hists cuts " << found << endl;
    }

  histoMap straightTracksY =   bookControlHists("straightTracksY",histoFile);
  histoMap straightTracksY_isoAandC =   bookControlHists("straightTracksY_isoAandC",histoFile);
  histoMap straightTracksY_isoAandCandB =   bookControlHists("straightTracksY_isoAandCandB",histoFile);
  histoMap straightTracksY_isoAandCandB_straightTracksX =   bookControlHists("straightTracksY_isoAandCandB_straightTracksX",histoFile);
  histoMap straightTracksY_isoAandCandB_chargerAandC  =   bookControlHists("straightTracksY_isoAandCandB_chargerAandC",histoFile);
  histoMap straightTracksY_isoAandCandB_chargerAandCandB  =   bookControlHists("straightTracksY_isoAandCandB_chargerAandCandB",histoFile);
  histoMap straightTracksY_isoAandCandB_chargeAandC  =   bookControlHists("straightTracksY_isoAandCandB_chargeAandC",histoFile);
  histoMap straightTracksY_isoAandCandB_chargeAandCandB  =   bookControlHists("straightTracksY_isoAandCandB_chargeAandCandB",histoFile);
  histoMap hitsOnTrack  =   bookControlHists("hitsOnTrack",histoFile);




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
	  double xA = xcoordinate(0, cA, alignxA, ptchc, ptchr);
	  double yA = ycoordinate(0, cA, alignyA, ptchc, ptchr);

	  double xAr = xA*cfA - yA*sfA;
	  double yAr = xA*sfA + yA*cfA;

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
	
	  double xB = xcoordinate(1, cB, 0, ptchc, ptchr);
	  double yB = ycoordinate(1, cB, 0, ptchc, ptchr);
	  
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
	      double xC = xcoordinate(2, cC, alignxC, ptchc, ptchr);
	      double yC = ycoordinate(2, cC, alignyC, ptchc, ptchr);

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

      ///////        A-C cluster correlations: ///////////
      
      nm = 0;
      
      for( vector<cluster>::iterator cC = vclC.begin(); cC != vclC.end(); ++cC )
	{
	
	  double xC = xcoordinate(2, cC, alignxC, ptchc, ptchr);
	  double yC = ycoordinate(2, cC, alignyC, ptchc, ptchr);
	  
	  double xCr = xC*cfC - yC*sfC;
	  double yCr = xC*sfC + yC*cfC;
	
      
	  double etaC = eta(cC);

	  for( vector<cluster>::iterator cA = vclA.begin(); cA != vclA.end(); ++cA )
	    {
	
              double xA = xcoordinate(0, cA, alignxA, ptchc, ptchr);
	      double yA = ycoordinate(0, cA, alignyA, ptchc, ptchr);
	      
	      double xAr = xA*cfA - yA*sfA;
	      double yAr = xA*sfA + yA*cfA;
	      
	      hxxCA->Fill( xAr, xCr );
	      hyyCA->Fill( yAr, yCr );
	  
	      double dxCA = xCr - xAr;
	      double dyCA = yCr - yAr;
	      
	      double dxyCA = sqrt( dxCA*dxCA + dyCA*dyCA );
	      hdxCA->Fill( dxCA );
	      hdyCA->Fill( dyCA );

	      if(PRINT) cout << "#################    EVENT " << iev << endl;
   
	      // NB until I initialize cB, I am using cA instead!!!!
	      fillControlHists(raw,"raw",0,0,cA,cA,cC,0,0,0,iev,0,0,xAr,yAr,xCr,yCr,dxCA,0,0,etaC,histoFile,fileName);


	      if(PRINT) cout << " tracks condition " << fabs( dxCA ) << " vs "  << straightTracks * beamDivergenceScaled + 0.02 << endl;

	      if( fabs( dxCA ) > straightTracks * beamDivergenceScaled + 0.02 ) continue; // includes beam divergence: +-5 sigma

	      fillControlHists(dxCAcut,"dxCAcut",0,0,cA,cA,cC,0,0,0,iev,0,0,xAr,yAr,xCr,yCr,dxCA,0,0,etaC,histoFile,fileName);


	      if(PRINT) cout << " dyCA " << dyCA << endl;
	      hdyCAc->Fill( dyCA );
	      dyvsyCA->Fill( yAr, dyCA );
	      
	      if( fabs( dyCA ) > straightTracks * beamDivergenceScaled + 0.1 ) continue; // [mm]
	      
	      ++nm;
	  
	      dxvsyCA->Fill( yAr, dxCA );
	      
	      dxvsxCA->Fill( xAr, dxCA ); // linear trend in run 392, 403: acceptance and beam divergence ?
	      
	      double xavg = 0.5 * ( xAr + xCr );
	      double yavg = 0.5 * ( yAr + yCr );
	      
	      //	      double xmod = fmod( xavg + 8, 0.05 ); // [mm] 0..0.05
	      double xmod = fmod( xavg + 8.0125, 0.025 ); // [mm] 0..0.05
	      if( fifty )
		xmod = fmod( xavg + 8.025, 0.05 ); // [mm] 0..0.05
	      
	      double etaA = eta(cA);

	      int eff[999] = {0};
	  
	      for( vector<cluster>::iterator cB = vclB.begin(); cB != vclB.end(); ++cB )
		{
	    
		  double xB = xcoordinate(1, cB, 0, ptchc, ptchr);
		  double yB = ycoordinate(1, cB, 0, ptchc, ptchr);

		  double etaB = eta(cB);
	    
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
		  double dy3 = yB - yavg;

		  fillControlHists(beforeCorrections,"beforeCorrections",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName);

		  dx3 = dx3 - dx3corr*xavg; // from -dx3vsx.Fit("pol1")

		  double dxy = sqrt( dx3*dx3 + dy3*dy3 );
		  fillControlHists(nocuts,"nocuts",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName);
	
		  if(PRINT) cout << "nocuts" << endl;
		  
		  
		  if( fabs( dy3 ) < straightTracks * beamDivergenceScaled + 0.05 )
		    { // cut on y, look at x, see madx3vsy
		    
		      fillControlHists(straightTracksY,"straightTracksY",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName);

		      if( cA->iso && cC->iso )
			{
			  fillControlHists(straightTracksY_isoAandC,"straightTracksY_isoAandC",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName);
			  if( cB->iso )
			    {
			      fillControlHists(straightTracksY_isoAandCandB,"straightTracksY_isoAandCandB",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName);
			      if(PRINT) cout << "isoAandCandB done" << endl;


			      if( fabs( dxCA ) < straightTracks * beamDivergenceScaled )
				{ // track angle
				  fillControlHists(straightTracksY_isoAandCandB_straightTracksX,"straightTracksY_isoAandCandB_straightTracksX",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName);
				}//fabs( dxCA ) < straightTracks * beamDivergenceScaled 

			      
			      if( (cA->q >= qR) && (cC->q >= qR))
				{
				  if(PRINT) cout << "cA->q " << cA->q << " >= qR " << qR << " && cC->q " << cC->q << "  >= qR " << endl;

				  if(PRINT) cout << "straightTracksY_isoAandCandB_chargerAandC going to be done" << endl;
			      
				  fillControlHists(straightTracksY_isoAandCandB_chargerAandC,"straightTracksY_isoAandCandB_chargerAandC",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName);
				  if(PRINT) cout << "straightTracksY_isoAandCandB_chargerAandC going done" << endl;
			      
				  if( ( cB->q >= qRB) && ( cA->q >= qR) && ( cC->q >= qR))
				    {
				      
				      fillControlHists(straightTracksY_isoAandCandB_chargerAandCandB,"straightTracksY_isoAandCandB_chargerAandCandB",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName);
				    }// ( cB->q >= qRB) && ( cA->q >= qR) && ( cC->q >= qR)

				  if( ( cA->q <= qL || cA->q >= qR) && ( cC->q <= qL || cC->q >= qR))
				    {
				      fillControlHists(straightTracksY_isoAandCandB_chargeAandC,"straightTracksY_isoAandCandB_chargeAandC",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName);

				      if( (cB->q <= qLB || cB->q >= qRB) && ( cA->q <= qL || cA->q >= qR) && ( cC->q <= qL || cC->q >= qR))
					{					  
					  fillControlHists(straightTracksY_isoAandCandB_chargeAandCandB,"straightTracksY_isoAandCandB_chargeAandCandB",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName);
					}//(cB->q <= qLB || cB->q >= qRB) && ( cA->q <= qL || cA->q >= qR) && ( cC->q <= qL || cC->q >= qR)

				    }//( cA->q <= qL || cA->q >= qR) && ( cC->q <= qL || cC->q >= qR)
			
				}//(cA->q >= qR) && (cC->q >= qR)

			    }//iso B
			}//iso A and C

    
		    } // cut dy (fabs( dy3 ) < straightTracks * beamDivergenceScaled + 0.05)
		  if(PRINT) cout << " done with hist filling my way " << endl;

		  
		  if(PRINT) cout << " filling hists on track " << endl;

		  if( fabs( dx3 ) < 0.07 && fabs( dy3 ) < 0.15 && cA->iso && cB->iso && cC->iso) // hit on track
		    {

		      if(PRINT) cout << "  it is on track " << endl;
		      
		      fillControlHists(hitsOnTrack,"hitsOnTrack",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName);
				    


		    } // linked, iso (hit on track)
		  
		  if(PRINT) cout << "  end of selections!  " << endl;
		  
		  for( int iw = 1; iw < 999; ++iw )
		    if( dxy < iw*0.010 )
		      eff[iw] = 1; // eff
		  if(PRINT) cout << "  hist filling 14 " << endl;

		} // clusters B
	      if(PRINT) cout << "  hist filling 15 " << endl;

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
  /////////    end of the big loop on the events!!! //////////////////
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

  if(DOALIGNMENT)
    {
      double newalignxA = alignxA;
      if(PRINT) cout << " newalignxA " << newalignxA << " alignxA " << alignxA << endl;
      
      if(PRINT) cout << hdxAB->GetTitle() << " entries " << hdxAB->GetEntries() << endl;
      
      if( hdxAB->GetEntries() > 999 )
	{
	  
	  newalignxA += alignx(hdxAB);
	}
      else
	cout << "  not enough statistics" << endl;
      
      // y:

      double newalignyA = alignyA;
      
      if(PRINT) cout << endl << hdyAB->GetTitle() << " entries " << hdyAB->GetEntries() << endl;
      if( hdyAB->GetEntries() > 999 ) {
	cout << "  y correction " << hdyAB->GetMean() << endl;
	newalignyA += aligny(hdyAB);
      }
      else
	cout << "  not enough statistics" << endl;
      
      // dxvsy -> -rot
      
      double newalignfA = alignfA;
      
      if(PRINT) cout << endl << dxvsyAB->GetTitle() << " entries " << dxvsyAB->GetEntries() << endl;
      if( aligniteration > 0 && dxvsyAB->GetEntries() > 999 ) {
	newalignfA += alignangle(dxvsyAB);
      }
      else
	cout << "  not enough" << endl;
      
      // C:

      double newalignxC = alignxC;
  
      if(PRINT) cout << endl << hdxCB->GetTitle() << " entries " << hdxCB->GetEntries() << endl;
      if( hdxCB->GetEntries() > 999 ) {   
	newalignxC += alignx(hdxCB);
      }
      else
	cout << "  not enough" << endl;

      // y:
      
      double newalignyC = alignyC;

      if(PRINT)      cout << endl << hdyCB->GetTitle() << " entries " << hdyCB->GetEntries() << endl;
      if( hdyCB->GetEntries() > 999 ) {
	cout << "  y correction " << hdyCB->GetMean() << endl;
	newalignyC += aligny(hdyCB);
      }
      else
	cout << "  not enough" << endl;
      
      // dxvsy -> -rot

      double newalignfC = alignfC;
      
      cout << endl << dxvsyCB->GetTitle() << " entries " << dxvsyCB->GetEntries() << endl;
      if( aligniteration > 0 && dxvsyCB->GetEntries() > 999 ) {
	newalignfC += alignangle(dxvsyCB);
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

      cout << "update alignment file? (y/n)" << endl;
      string ans;
      cin >> ans;
      string YES{"y"};
      if( ans == YES ) {
      
        ofstream alignFile( alignFileName );
      
        alignFile << "# alignment for run " << run << endl;
        alignFile << "iteration " << aligniteration << endl;
        alignFile << "alignxA " << setw(11) << newalignxA << endl;
        alignFile << "alignyA " << setw(11) << newalignyA << endl;
        alignFile << "alignfA " << setw(11) << newalignfA << endl;
        alignFile << "alignxC " << setw(11) << newalignxC << endl;
        alignFile << "alignyC " << setw(11) << newalignyC << endl;
        alignFile << "alignfC " << setw(11) << newalignfC << endl;
	alignFile << "gainA " << setw(11) << gainA << endl;
	alignFile << "gainB " << setw(11) << gainB << endl;
	alignFile << "gainC " << setw(11) << gainC << endl;
	alignFile << "keA " << setw(11) << keA << endl;
	alignFile << "keB " << setw(11) << keB << endl;
	alignFile << "keC " << setw(11) << keC << endl;
	alignFile << "beamEnergy " << setw(11) << beamEnergy << endl;
	alignFile << "pitch " << setw(11) << pitch << endl;
	alignFile << "qL " << setw(11) << qL << endl;
	alignFile << "qR " << setw(11) << qR << endl;
	alignFile << "qLB " << setw(11) << qLB << endl;
	alignFile << "qRB " << setw(11) << qRB << endl;
	alignFile << "TsunamiA " << setw(11) << Tsunami[A] << endl;
	alignFile << "TsunamiB " << setw(11) << Tsunami[B] << endl;
	alignFile << "TsunamiC " << setw(11) << Tsunami[C] << endl;
	alignFile << "dphcutA " << setw(11) << dphcut[A] << endl;
	alignFile << "dphcutB " << setw(11) << dphcut[B] << endl;
	alignFile << "dphcutC " << setw(11) << dphcut[C] << endl;
	alignFile << "dx3c " << setw(11) << dx3corr << endl;
        alignFile.close();
      }
    }//DOALIGNMENT
  cout << endl << histoFile->GetName() << endl;
  histoFile->Close();
  if(PRINT) cout << " closed the hist file" << endl;

  cout << endl;

  return 0;
}


/////// end of main code



//------------------------------------------------------------------------------

double alignx(TH1I * h)
{
  TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -10, 10 );
  double xpk = h->GetBinCenter( h->GetMaximumBin() );
  fgp0->SetParameter( 0, h->GetMaximum() ); // amplitude
  fgp0->SetParameter( 1, xpk ); //mean
  fgp0->SetParameter( 2, h->GetBinWidth(1) ); // sigma
  fgp0->SetParameter( 3, h->GetBinContent( hdxAB->FindBin(xpk-1) ) ); // BG
  h->Fit( "fgp0", "q", "", xpk-1, xpk+1 ); // fit range around peak
  cout << "  Fit Gauss + BG:"
       << endl << "  area " << fgp0->GetParameter(0)
       << endl << "  mean " << fgp0->GetParameter(1)
       << endl << "  sigm " << fgp0->GetParameter(2)
       << endl << "  offs " << fgp0->GetParameter(3)
       << endl;
  return fgp0->GetParameter(1);
}


double aligny(TH1I * h)
{
  return h->GetMean();
}

double alignangle(TProfile * h)
{
  h->Fit( "pol1", "q", "", -3, 3 );
  TF1 * fdxvsy = h->GetFunction( "pol1" );
  return fdxvsy->GetParameter(1);
}


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
	  if(c.q < 0)	  cout << "negative cluster charge " << c.q << endl;
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

	    //if(PRINT) cout << " ph4 before " << ph4;
	    ph4 -= Tsunami*phprev; // Tsunami
	    //if(PRINT) cout << " ph4 after " << ph4 << endl;
	    
	    phprev = vpx[ipx].ph; // original ph4

	    double dph;
	    
	    //find closest row to define the difference
	    // if( row4 - row1 < row7 - row4 )
	    //   dph = ph4 - ph1;
	    // else
	    //   dph = ph4 - ph7;

	    dph = ph4 - (ph1 + ph7)/2;  // Finn's common mode, try and give feedback
	    
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
		// subtract Vcal offset:

		double U0 = -p3[plane][col4][row4] / p2[plane][col4][row4]; // dph = 0
		double v0 = p0[plane][col4][row4] - p1[plane][col4][row4] * log( (1-U0)/U0 ); // inverse Fermi

		double q = ke[plane] * ( vcal - v0 );
		
		px.q = q;//ke[plane]*vcal;
		
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
	  double isolationCut = 0.3; // [mm]
	  double pitchc = 0.100; // [mm]
	  double pitchr = 0.025; // [mm]
	  if(fifty)
	    {
	      pitchc = 0.050;
	      pitchr = 0.050;
	    }
	  for( unsigned ipx = 0; ipx < vcl[icl].vpix.size(); ++ipx ){
	  for( unsigned jpx = 0; jpx < vcl[jcl].vpix.size(); ++jpx )
	    {
	      double dcol = fabs( vcl[icl].vpix[ipx].col - vcl[jcl].vpix[jpx].col)*pitchc;
	      double drow = fabs( vcl[icl].vpix[ipx].row - vcl[jcl].vpix[jpx].row )*pitchr; 
	      if(  dcol  < isolationCut && drow < isolationCut )
		{
		  if( vcl[icl].q < vcl[jcl].q ) // Thu 22.3.2018
		    vcl[icl].iso = 0; // flag smaller cluster
		  else
		    vcl[jcl].iso = 0;
		  done = 1;
		  break; // jpx
		}
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

double eta(vector<cluster>::iterator c)
{
  double eta = -2;
  if( c->size == 2 )
    {
      double q0 = c->vpix[0].q;
      double q1 = c->vpix[1].q;
      eta = (q1-q0)/(q1+q0);
    }
  return eta;
}

	  
// ******** xcoordinate and ycoordinate functions:
//to go from the row/column of the cluster to a (x,y)  coordinate on the sensor. being (0,0) the center of the sensor and keeping the alignment into account
double xcoordinate(int plane, vector<cluster>::iterator c, double align, double pitchc, double pitchr)
{
  double variable;
  variable = c->row*pitchr - halfSensorX - align; //[mm]
  // if( run == 431 && plane == 1 )
  //   variable = c->row*ptchr - halfSensorX; // rot90
  if( fifty )
    variable = c->col*pitchc - halfSensorY - align; // straight
  
  return variable;
}
      
double ycoordinate(int plane, vector<cluster>::iterator c, double align, double pitchc, double pitchr)
{
  double variable;
  
  variable =  c->col*pitchc - halfSensorY - align; //[mm]
  // if( run == 431 && plane == 1  )
  //   variable =  c->col*ptchc - halfSensorY; // PCB
  if( fifty )
    variable = c->row*pitchr - halfSensorX - align; // PCB

  return variable;
}
//****************

histoMap  bookControlHists(TString selection, TFile * histofile)
{
  if(PRINT) cout << " boking my hists " << selection << endl;
  TDirectory *cdtof = histofile->mkdir(selection);
  cdtof->cd(); 
  int nbx =  80; //number of bins x
  int nby = 320; //number of bins y
  if( fifty ) {
    nbx = 160;
    nby = 160;
  }
  

  histoMap mapOfHists;

  TH1I * hxA = new TH1I("xA","x Ar;x [mm];clusters A", 100, -5, 5 );
  TH1I * hyA = new TH1I("yA","y Ar;y [mm];clusters A", 100, -5, 5 ); 
  TH1I * hxB = new TH1I("xB","x Br;x [mm];clusters B", 100, -5, 5 );
  TH1I * hyB = new TH1I("yB","y Br;y [mm];clusters B", 100, -5, 5 ); 
  TH1I * hxC = new TH1I("xC","x Cr;x [mm];clusters C", 100, -5, 5 );
  TH1I * hyC = new TH1I("yC","y Cr;y [mm];clusters C", 100, -5, 5 ); 
  
  TH1I * hdx3 = new TH1I("dx3", "triplet dx; dx [mm];triplets", 500, -0.5, 0.5 );// "dx", "x A " +selection+ " ;x [mm];clusters A", 100, -5, 5 );
  TH1I * hdy3 = new TH1I("dy3", "triplet dy; dy [mm];triplets", 200, -1., 1. );// "dx", "x A " +selection+ " ;x [mm];clusters A", 100, -5, 5 );

  TH1I * hclsizeA = new TH1I("clsizeA", "A cluster size "+selection+";cluster size [pixels];A clusters on tracks", 40, 0.5, 40.5 ); 
  TH1I * hclsizeB = new TH1I("clsizeB", "B cluster size "+selection+";cluster size [pixels];B clusters on tracks", 40, 0.5, 40.5 ); 
  TH1I * hclsizeC = new TH1I("clsizeC", "C cluster size "+selection+";cluster size [pixels];C clusters on tracks", 40, 0.5, 40.5 ); 

  TH1I * hnrowB = new TH1I("nrowB", "B number of rows "+selection+";number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 ); 
  TH1I * hncolB = new TH1I("ncolB", "B number of cols "+selection+";number of cols [pixels];B clusters on tracks", 40, 0.5, 40.5 ); 

  TH2I * hclmapA = new TH2I("clmapA","A cluster map;col;row;A clusters", nbx, 0, nbx, nby, 0, nby );
  TH2I * hclmapB = new TH2I("clmapB","B cluster map;col;row;B clusters", nbx, 0, nbx, nby, 0, nby );
  TH2I * hclmapC = new TH2I("clmapC","C cluster map;col;row;C clusters", nbx, 0, nbx, nby, 0, nby );

  TH1I * hclchargeA = new TH1I("clchargeA", "A cluster charge "+selection+";cluster charge [ke];A clusters on tracks", 160, 0., 80. ); 
  TH1I * hclchargeB = new TH1I("clchargeB", "B cluster charge "+selection+";cluster charge [ke];B clusters on tracks", 160, 0., 80. ); 
  TH1I * hclchargeC = new TH1I("clchargeC", "C cluster charge "+selection+";cluster charge [ke];C clusters on tracks", 160, 0., 80. ); 

  TH1I * hclphA = new TH1I("clphA", "A cluster ph "+selection+";cluster ph [ADC];A clusters on tracks", 200, 0., 1000. ); 
  TH1I * hclphB = new TH1I("clphB", "B cluster ph "+selection+";cluster ph [ADC];B clusters on tracks", 200, 0., 1000. ); 
  TH1I * hclphC = new TH1I("clphC", "C cluster ph "+selection+";cluster ph [ADC];C clusters on tracks", 200, 0., 1000. ); 

  TH1I * hetaA = new TH1I("etaA", "A eta "+selection+";eta;A 2-pix clusters on tracks",            100, -1, 1 ); 
  TH1I * hetaB = new TH1I("etaB", "B eta "+selection+";eta;B 2-pix clusters on tracks",            100, -1, 1 ); 
  TH1I * hetaC = new TH1I("etaC", "C eta "+selection+";eta;C 2-pix clusters on tracks",            100, -1, 1 ); 


  TH1I * hpxphA = new TH1I( "pxphA", "A pixel PH;pixel PH [ADC];A pixels on tracks", 250, 0, 500 ); 
  TH1I * hpxchargeA = new TH1I( "pxchargeA", "A pixel charge;pixel charge [ke];A pixels on tracks",             100, 0, 20 );     
  TH1I * hpxphB = new TH1I( "pxphB", "B pixel PH;pixel PH [ADC];B pixels on tracks", 250, 0, 500 ); 
  TH1I * hpxchargeB = new TH1I( "pxchargeB", "B pixel charge;pixel charge [ke];B pixels on tracks",             100, 0, 20 );     
  TH1I * hpxphC = new TH1I( "pxphC", "C pixel PH;pixel PH [ADC];C pixels on tracks", 250, 0, 500 ); 
  TH1I * hpxchargeC = new TH1I( "pxchargeC", "C pixel charge;pixel charge [ke];C pixels on tracks",             100, 0, 20 );     

  TH1I * hpxchargeB1st = new TH1I( "pxchargeB1st", "B 1st pixel charge;pixel charge [ke];B pixels on tracks",             100, 0, 20 );     
  TH1I * hpxchargeB2nd = new TH1I( "pxchargeB2nd", "B 2nd pixel charge;pixel charge [ke];B pixels on tracks",             100, 0, 20 );     
		    
  TProfile * hnrowvsxmB3;
  if(fifty)
    {
      hnrowvsxmB3 = new TProfile( "nrowvsxmB3","B rows vs xmod @ res;x mod 50 [#mum];<B cluster size [rows]>",50, 0, 50, 0.5, 10.5 );
    }
  else
    {
      hnrowvsxmB3 = new TProfile( "nrowvsxmB3","B rows vs xmod @ res;x mod 25 [#mum];<B cluster size [rows]>",25, 0, 25, 0.5, 10.5 );
    }
  

  TProfile * clqvsxmB3= new TProfile( "clqvsxmB3",                "B cluster charge vs xmod;x mod 50 [#mum];<B cluster charge [ke]>",                     50, 0, 50, 0, 50 );

  TProfile * dx3vsev   = new TProfile( "dx3vsev", "dx3 vs time;trigger;<dx3> [mm]",                310, 0, 3100*1000, -0.5, 0.5 ); 
  TProfile * dx3vsx    = new TProfile( "dx3vsx", "dx vs x;x [mm];<dx3> [mm]", 320, -4, 4, -0.5, 0.5 );//turn 
  TProfile * dx3vsy    = new TProfile( "dx3vsy", "dx vs y;y [mm];<dx3> [mm]",  80, -4, 4, -0.5, 0.5 );// rot
  TProfile * dx3vsxm   = new TProfile( "dx3vsxm", "dx vs x mod 50 um;x mod 50 [#mum];<dx3> [mm]",                  50, 0, 50, -0.5, 0.5 ); 
  TProfile * madx3vsdx = new TProfile( "madx3vsdx", "MAD;// (dx3) vs dx C-A;C-A dx [#mum];MAD dx3 [mm]",                   100, -100, 100, 0, 0.1 );// dxCA
  TProfile * madx3vsx  = new TProfile( "madx3vsx", "MAD;// (dx3) vs x;x [mm];MAD dx3 [mm]", 320, -4, 4, 0, 0.1 );
  TProfile * madx3vsy  = new TProfile( "madx3vsy", "MAD;// (dx3) vs y;y [mm];MAD dx3 [mm]",  80, -4, 4, 0, 0.1 );
  TProfile * madx3vsxm;

  if( fifty)
    madx3vsxm = new TProfile( "madx3vsxm", "MAD = new  (dx3) vs xmod;x mod 50 [#mum];MAD dx3 [mm]",                    50, 0, 50, 0, 0.1 );
  else
    madx3vsxm = new TProfile( "madx3vsxm", "MAD = new  (dx3) vs xmod;x mod 25 [#mum];MAD dx3 [mm]",                    25, 0, 25, 0, 0.1 );

  TProfile * madx3vsq = new TProfile( "madx3vsq", "MAD;// (dx) vs Q;B cluster charge [ke];MAD dx [mm]", 100, 0, 100, 0, 0.1 );
  TProfile * madx3vsn = new TProfile( "madx3vsn", "MAD;// (dx) vs cluster size;B cluster size [pixels];MAD dx [mm]", 20, 0.5, 20.5, 0, 0.1 );
  
  TProfile * etavsxmB3 = new TProfile( "etavsxmB3", "eta vs xmod;x mod 50 [#mum];B <eta>",                 50, 0, 50, -1.1, 1.1 );//sine
  TProfile * madx3vseta= new TProfile( "madx3vseta", "MAD;// (dx3) vs eta;eta;MAD dx3 [mm]",                     100, -1, 1, 0, 0.1 );   //flat

  TH1I * hdx3_clsizeB1 = new TH1I("dx3_clsizeB1", "triplet dx_clsizeB1; dx [mm];triplets", 500, -0.5, 0.5 );// "dx", "x A " +selection+ " ;x [mm];clusters A", 100, -5, 5 );
  TH1I * hdx3_clsizeB2 = new TH1I("dx3_clsizeB2", "triplet dx_clsizeB2; dx [mm];triplets", 500, -0.5, 0.5 );// "dx", "x A " +selection+ " ;x [mm];clusters A", 100, -5, 5 );
  TH1I * hdx3_clsizeB3 = new TH1I("dx3_clsizeB3", "triplet dx_clsizeB3; dx [mm];triplets", 500, -0.5, 0.5 );// "dx", "x A " +selection+ " ;x [mm];clusters A", 100, -5, 5 );
  TH1I * hdx3_clsizeB4 = new TH1I("dx3_clsizeB4", "triplet dx_clsizeB4; dx [mm];triplets", 500, -0.5, 0.5 );// "dx", "x A " +selection+ " ;x [mm];clusters A", 100, -5, 5 );
  TH1I * hdx3_clsizeB5 = new TH1I("dx3_clsizeB5", "triplet dx_clsizeB5; dx [mm];triplets", 500, -0.5, 0.5 );// "dx", "x A " +selection+ " ;x [mm];clusters A", 100, -5, 5 );
  TH1I * hdx3_clsizeB6 = new TH1I("dx3_clsizeB6", "triplet dx_clsizeB6; dx [mm];triplets", 500, -0.5, 0.5 );// "dx", "x A " +selection+ " ;x [mm];clusters A", 100, -5, 5 );
  TH1I * hdx3_clsizeB7m = new TH1I("dx3_clsizeB7m", "triplet dx_clsizeB7m; dx [mm];triplets", 500, -0.5, 0.5 );// "dx", "x A " +selection+ " ;x [mm];clusters A", 100, -5, 5 );
  
  if(PRINT) cout << "going to insert first hist" << endl;

  mapOfHists.insert(std::make_pair("xA",hxA));
  mapOfHists.insert(std::make_pair("yA",hyA));
  mapOfHists.insert(std::make_pair("xB",hxB));
  mapOfHists.insert(std::make_pair("yB",hyB));
  mapOfHists.insert(std::make_pair("xC",hxC));
  mapOfHists.insert(std::make_pair("yC",hyC));
  
  mapOfHists.insert(std::make_pair("dx3",hdx3));
  mapOfHists.insert(std::make_pair("dy3",hdy3));

  mapOfHists.insert(std::make_pair("clsizeA",hclsizeA));
  mapOfHists.insert(std::make_pair("clsizeB",hclsizeB));
  mapOfHists.insert(std::make_pair("clsizeC",hclsizeC));

  mapOfHists.insert(std::make_pair("nrowB",hnrowB));
  mapOfHists.insert(std::make_pair("ncolB",hncolB));

  mapOfHists.insert(std::make_pair("clmapA",hclmapA));
  mapOfHists.insert(std::make_pair("clmapB",hclmapB));
  mapOfHists.insert(std::make_pair("clmapC",hclmapC));

  mapOfHists.insert(std::make_pair("clchargeA",hclchargeA));
  mapOfHists.insert(std::make_pair("clchargeB",hclchargeB));
  mapOfHists.insert(std::make_pair("clchargeC",hclchargeC));

  mapOfHists.insert(std::make_pair("clphA",hclphA));
  mapOfHists.insert(std::make_pair("clphB",hclphB));
  mapOfHists.insert(std::make_pair("clphC",hclphC));

  mapOfHists.insert(std::make_pair("etaA",hetaA));
  mapOfHists.insert(std::make_pair("etaB",hetaB));
  mapOfHists.insert(std::make_pair("etaC",hetaC));

  mapOfHists.insert(std::make_pair("pxphA",hpxphA));
  mapOfHists.insert(std::make_pair("pxchargeA",hpxchargeA));
  mapOfHists.insert(std::make_pair("pxphB",hpxphB));
  mapOfHists.insert(std::make_pair("pxchargeB",hpxchargeB));
  mapOfHists.insert(std::make_pair("pxphC",hpxphC));
  mapOfHists.insert(std::make_pair("pxchargeC",hpxchargeC));

  mapOfHists.insert(std::make_pair("pxchargeB1st",hpxchargeB1st));
  mapOfHists.insert(std::make_pair("pxchargeB2nd",hpxchargeB2nd));
  
  mapOfHists.insert(std::make_pair("nrowvsxmB3",hnrowvsxmB3));

  mapOfHists.insert(std::make_pair("clqvsxmB3",clqvsxmB3));

  mapOfHists.insert(std::make_pair("dx3vsev",dx3vsev));
  mapOfHists.insert(std::make_pair("dx3vsx",dx3vsx));
  mapOfHists.insert(std::make_pair("dx3vsy",dx3vsy));
  mapOfHists.insert(std::make_pair("dx3vsxm",dx3vsxm));
  mapOfHists.insert(std::make_pair("madx3vsdx",madx3vsdx));
  mapOfHists.insert(std::make_pair("madx3vsx",madx3vsx));
  mapOfHists.insert(std::make_pair("madx3vsy",madx3vsy));
  mapOfHists.insert(std::make_pair("madx3vsxm",madx3vsxm));
  mapOfHists.insert(std::make_pair("etavsxmB3",etavsxmB3));
  mapOfHists.insert(std::make_pair("madx3vseta",madx3vseta));
  mapOfHists.insert(std::make_pair("madx3vsq",madx3vsq));
  mapOfHists.insert(std::make_pair("madx3vsn",madx3vsn));

  mapOfHists.insert(std::make_pair("dx3_clsizeB1",hdx3_clsizeB1));
  mapOfHists.insert(std::make_pair("dx3_clsizeB2",hdx3_clsizeB2));
  mapOfHists.insert(std::make_pair("dx3_clsizeB3",hdx3_clsizeB3));
  mapOfHists.insert(std::make_pair("dx3_clsizeB4",hdx3_clsizeB4));
  mapOfHists.insert(std::make_pair("dx3_clsizeB5",hdx3_clsizeB5));
  mapOfHists.insert(std::make_pair("dx3_clsizeB6",hdx3_clsizeB6));
  mapOfHists.insert(std::make_pair("dx3_clsizeB7m",hdx3_clsizeB7m));

  
  return mapOfHists;
}


void  fillControlHists(histoMap mapOfHists, TString selection, double dx3, double dy3, vector<cluster>::iterator clusterA, vector<cluster>::iterator clusterB, vector<cluster>::iterator clusterC, int nrowB, int ncolB,double xmod, unsigned iev,double xB, double yB,double xAr, double yAr, double xCr, double yCr,double dxCA,double etaA, double etaB, double etaC,TFile * histofile, TString fileName)
{

  if(PRINT) std::cout << "********** filling hists of "<< selection << endl;
  
  TDirectory *cdtof = (TDirectory *)histofile->Get(selection);
  cdtof->cd(fileName+":"+selection);

  if(PRINT)  cout << " directory " << fileName << ":"<< selection << " found "<< endl;
  if(PRINT)cout << " cl B " << clusterB->size << endl;
  if(PRINT) cout << "filling in map" << endl;

	
  
  auto search =  mapOfHists.find("xA");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(xAr);
  } else {
    if(PRINT) std::cout << "Not found "<< "xA" << endl;
  }

  search =  mapOfHists.find("yA");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(yAr);
  } else {
    if(PRINT) std::cout << "Not found "<< "yAr" << endl;
  }

  search =  mapOfHists.find("xB");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(xB);
  } else {
    if(PRINT) std::cout << "Not found "<< "xB" << endl;
  }

  search =  mapOfHists.find("yB");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(yB);
  } else {
    if(PRINT) std::cout << "Not found "<< "yB" << endl;
  }

  search =  mapOfHists.find("xC");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(xCr);
  } else {
    if(PRINT) std::cout << "Not found "<< "xC" << endl;
  }

  search =  mapOfHists.find("yC");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(yCr);
  } else {
    if(PRINT) std::cout << "Not found "<< "yCr" << endl;
  }


  
  search =  mapOfHists.find("dx3");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3" << endl;
  }


  search =  mapOfHists.find("dy3");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(dy3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dy3" << endl;
  }

  search =  mapOfHists.find("clsizeA");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterA->size);
  } else {
    if(PRINT) std::cout << "Not found "<< "clsizeA" << endl;
  }

  search =  mapOfHists.find("clsizeB");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterB->size);
  } else {
    if(PRINT) std::cout << "Not found "<< "clsizeB" << endl;
  }

    search =  mapOfHists.find("clsizeC");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterC->size);
  } else {
    if(PRINT) std::cout << "Not found "<< "clsizeC" << endl;
  }


  search =  mapOfHists.find("nrowB");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(nrowB);
  } else {
    if(PRINT) std::cout << "Not found "<< "nrowB" << endl;
  }

  search =  mapOfHists.find("ncolB");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(ncolB);
  } else {
    if(PRINT) std::cout << "Not found "<< "ncolB" << endl;
  }

  search =  mapOfHists.find("clmapA");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterA->col, clusterA->row);
  } else {
    if(PRINT) std::cout << "Not found "<< "clmapA" << endl;
  }

  search =  mapOfHists.find("clmapB");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterB->col, clusterB->row);
  } else {
    if(PRINT) std::cout << "Not found "<< "clmapB" << endl;
  }

  search =  mapOfHists.find("clmapC");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterC->col, clusterC->row);
  } else {
    if(PRINT) std::cout << "Not found "<< "clmapC" << endl;
  }


  search =  mapOfHists.find("clchargeA");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterA->q);
  } else {
    if(PRINT) std::cout << "Not found "<< "clchargeA" << endl;
  }

    search =  mapOfHists.find("clchargeB");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterB->q);
  } else {
    if(PRINT) std::cout << "Not found "<< "clchargeB" << endl;
  }

    search =  mapOfHists.find("clchargeC");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterC->q);
  } else {
    if(PRINT) std::cout << "Not found "<< "clchargeC" << endl;
  }

  search =  mapOfHists.find("clphA");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterA->sum);
  } else {
    if(PRINT) std::cout << "Not found "<< "clphA" << endl;
  }

    search =  mapOfHists.find("clphB");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterB->sum);
  } else {
    if(PRINT) std::cout << "Not found "<< "clphB" << endl;
  }

  search =  mapOfHists.find("clphC");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterC->sum);
  } else {
    if(PRINT) std::cout << "Not found "<< "clphC" << endl;
  }

  search =  mapOfHists.find("pxphA");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    for( int ipx = 0; ipx < clusterA->size; ++ipx )     search->second->Fill( clusterA->vpix[ipx].ph );
  } else {
    if(PRINT) std::cout << "Not found "<< "pxphA" << endl;
  }

  search =  mapOfHists.find("pxchargeA");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    for( int ipx = 0; ipx < clusterA->size; ++ipx )     search->second->Fill( clusterA->vpix[ipx].q );
  } else {
    if(PRINT) std::cout << "Not found "<< "pxchargeA" << endl;
  }

  search =  mapOfHists.find("pxphB");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    for( int ipx = 0; ipx < clusterB->size; ++ipx )     search->second->Fill( clusterB->vpix[ipx].ph );
  } else {
    if(PRINT) std::cout << "Not found "<< "pxphB" << endl;
  }

  search =  mapOfHists.find("pxchargeB");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    for( int ipx = 0; ipx < clusterB->size; ++ipx )     search->second->Fill( clusterB->vpix[ipx].q );
  } else {
    if(PRINT) std::cout << "Not found "<< "pxchargeB" << endl;
  }

  search =  mapOfHists.find("pxphC");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    for( int ipx = 0; ipx < clusterC->size; ++ipx )     search->second->Fill( clusterC->vpix[ipx].ph );
  } else {
    if(PRINT) std::cout << "Not found "<< "pxphC" << endl;
  }

  search =  mapOfHists.find("pxchargeC");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    for( int ipx = 0; ipx < clusterC->size; ++ipx )     search->second->Fill( clusterC->vpix[ipx].q );
  } else {
    if(PRINT) std::cout << "Not found "<< "pxchargeC" << endl;
  }


  search =  mapOfHists.find("etaA");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    if(clusterA->size ==2)
      search->second->Fill(etaA);
  } else {
    if(PRINT) std::cout << "Not found "<< "etaA" << endl;
  }

  search =  mapOfHists.find("etaB");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size ==2)
      search->second->Fill(etaB);
  } else {
    if(PRINT) std::cout << "Not found "<< "etaB" << endl;
  }

  search =  mapOfHists.find("etaC");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    if(clusterC->size ==2)
      search->second->Fill(etaC);
  } else {
    if(PRINT) std::cout << "Not found "<< "etaA" << endl;
  }
		    
  search =  mapOfHists.find("pxchargeB1st");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size ==2)
      search->second->Fill( clusterB->vpix[0].q );
  } else {
    if(PRINT) std::cout << "Not found "<< "pxchargeB1st" << endl;
  }

  search =  mapOfHists.find("pxchargeB2nd");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size ==2)
      search->second->Fill( clusterB->vpix[1].q );
  } else {
    if(PRINT) std::cout << "Not found "<< "pxchargeB2nd" << endl;
  }




  search =  mapOfHists.find("nrowvsxmB3");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';

  if( fifty )
    search->second->Fill( xmod*1E3, ncolB );
  else
    search->second->Fill( xmod*1E3, nrowB );

  } else {
    if(PRINT) std::cout << "Not found "<< "nrowvsxmB3" << endl;
  }

  search =  mapOfHists.find("clqvsxmB3");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill( xmod*1E3, clusterB->q );
  } else {
    if(PRINT) std::cout << "Not found "<< "clqvsxmB3" << endl;
  }

  search =  mapOfHists.find("dx3vsev");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill( iev, dx3 );
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3vsev" << endl;
  }

  search =  mapOfHists.find("dx3vsx");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill( xB, dx3 ); // turn 
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3vsx" << endl;
  }

  search =  mapOfHists.find("dx3vsy");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill( yB, dx3 ); // rot
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3vsy" << endl;
  }

  search =  mapOfHists.find("dx3vsxm");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill( xmod*1E3, dx3 ); // rot
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3vsxm" << endl;
  }
  search =  mapOfHists.find("madx3vsdx");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill( dxCA*1E3, fabs(dx3) ); // dxCA
  } else {
    if(PRINT) std::cout << "Not found "<< "madx3vsdx" << endl;
  }

  search =  mapOfHists.find("madx3vsx");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill( xB, fabs(dx3) );
  } else {
    if(PRINT) std::cout << "Not found "<< "madx3vsx" << endl;
  }

  search =  mapOfHists.find("madx3vsy");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill( yB, fabs(dx3) );
      } else {
    if(PRINT) std::cout << "Not found "<< "madx3vsy" << endl;
  }

  search =  mapOfHists.find("madx3vsxm");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(xmod*1E3, fabs(dx3) );
      } else {
    if(PRINT) std::cout << "Not found "<< "madx3vsxm" << endl;
  }

  search =  mapOfHists.find("madx3vsq");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterB->q, fabs(dx3) );
      } else {
    if(PRINT) std::cout << "Not found "<< "madx3vsq" << endl;
  }

  search =  mapOfHists.find("madx3vsn");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterB->size, fabs(dx3) );
      } else {
    if(PRINT) std::cout << "Not found "<< "madx3vsn" << endl;
  }

  search =  mapOfHists.find("etavsxmB3");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end() && clusterB->size == 2) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill( xmod*1E3, etaB ); // sine 
      } else {
    if(PRINT) std::cout << "Not found "<< "etavsxmB3" << endl;
  }
  
  search =  mapOfHists.find("madx3vseta");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end() && clusterB->size == 2) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    search->second->Fill( etaB, fabs(dx3) ); // flat 
      } else {
    if(PRINT) std::cout << "Not found "<< "madx3vseta" << endl;
  }
  

  search =  mapOfHists.find("dx3_clsizeB1");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size ==1)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clsizeB1" << endl;
  }

  search =  mapOfHists.find("dx3_clsizeB2");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size ==2)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clsizeB2" << endl;
  }

  search =  mapOfHists.find("dx3_clsizeB3");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size ==3)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clsizeB3" << endl;
  }

  search =  mapOfHists.find("dx3_clsizeB4");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size ==4)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clsizeB4" << endl;
  }
  search =  mapOfHists.find("dx3_clsizeB5");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size ==5)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clsizeB5" << endl;
  }
  search =  mapOfHists.find("dx3_clsizeB6");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size ==6)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clsizeB6" << endl;
  }
  search =  mapOfHists.find("dx3_clsizeB7m");
  if(PRINT) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(PRINT) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size >6)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clsizeB7m" << endl;
  }

  if(PRINT) std::cout << "filled all hists of "<< selection << endl;
  
}//fillControlHists


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


