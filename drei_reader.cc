

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
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "./drei_reader.h"

bool PRINT = false;
bool DOALIGNMENT = false;
bool DPHCUT = false;
bool DEBUG = false;
bool DO1CL = false;
TString method = "mine";
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
  if(alignversion == 3)  
    alignFileName = alignpath+"align_v3_" + runnum + ".dat";

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
  if(DPHCUT)
    {
      int i_dph = (int) dphcut[1];
      //      stringstream stream;
      //stream << fixed << setprecision(3) << dphcut[1];
      //string s = stream.str();
      fileName.Form( "/home/zoiirene/Output/drei-r%i_irene_dphcut%i.root", run,i_dph);
    }
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
  histoMap closest3M = bookControlHists("closest3M",histoFile);
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
  /* oneplane function: for each plane                                                                                                                                                                                                                                        
- Read data file & assemble pixel hits                                                                                                                                                                                                                                        
- Skip noisy events with more than 400 hits                                                                                                                                                                                                                                    
- Tsunami correction                                                                                                                                                                                                                                                           
- Common mode correction (column-wise)                                                                                                                                                                                                                                         
  - in roc coordinates                                                                                                                                                                                                                                                         
  - we read along one column                                                                                                                                                                                                                                                   
  - a roi is 7 pixels long                                                                                                                                                                                                                                                    
  - we find the pixel with the highest and lowest row index in one column                                                                                                                                                                                                     
  - we assume there is no significant charge in these pixels                                                                                                                                                                                                                   
  - we take the average pulse height of these two pixels and subtract it from those in between  dph_i = ph_i - (ph_up + ph_low)/2                                                                                                                                              
- Continue analysis if dph> dphcut                                                                                                                                                                                                                                             
  - for each sample we use the cut yielding the best resolution                                                                                                                                                                                                                
  - 25 vs 50                                                                                                                                                                                                                                                                   
  - Conversion from ADC to ke (updated for Vcal offset correction - as Finn)                                                                                                                                                                                                   
- Clustering for events with less than 50 hits (speeding)                                                                                                                                                                                                                      
  - From a seed pixel add to it all adjacent pixels with a hit                                                                                                                                                                                                                 
  - Evaluating the cluster isolation                                                                                                                                                                                                                                           
  */
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
      if(PRINT) cout << "#################    EVENT " << iev << endl;

      vector <cluster> vclA = *evA;
      vector <cluster> vclB = *evB;
      vector <cluster> vclC = *evC;

      vector<closest> AMatchC ( evlistA.size(), closest() );
      vector<closest> CMatchA ( evlistC.size(), closest() );
      vector<closest> AMatchB ( evlistA.size(), closest() );
      vector<closest> BMatchA ( evlistB.size(), closest() );
      vector<closest> CMatchB ( evlistC.size(), closest() );
      vector<closest> BMatchC ( evlistB.size(), closest() );
      
      
      ++iev;
      if(DEBUG)      cout << " event " << iev << " cl A " << vclA.size() << " B " << vclB.size() << " C " << vclC.size() << endl;
     

      if(DO1CL){
	cout << "########################## 1 cluster per plane!!! ################# "<< endl;
	if(vclA.size()!=1 || vclB.size()!=1 || vclC.size() !=1){
	  cout << " more or less  than 1 clusters on at least one plane! " << endl;
	  continue;
	}
      }

      nclustA = vclA.size();
      nclustB = vclB.size();
      nclustC = vclC.size();


      
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
      

      ///////////////              Matching for closest
      cout << " ########## Matching for closest  ########## "<< endl;
      ///////        A-C cluster correlations: ///////////
      if(PRINT) cout << "entering AC correlation loop " << endl;
      for( vector<cluster>::iterator cC = vclC.begin(); cC != vclC.end(); ++cC ){
      	
	  double xC = xcoordinate(2, cC, alignxC, ptchc, ptchr);
	  double yC = ycoordinate(2, cC, alignyC, ptchc, ptchr);
	  
	  double xCr = xC*cfC - yC*sfC;
	  double yCr = xC*sfC + yC*cfC;
	
      
	  for( vector<cluster>::iterator cA = vclA.begin(); cA != vclA.end(); ++cA )   {
	
              double xA = xcoordinate(0, cA, alignxA, ptchc, ptchr);
	      double yA = ycoordinate(0, cA, alignyA, ptchc, ptchr);
	      
	      double xAr = xA*cfA - yA*sfA;
	      double yAr = xA*sfA + yA*cfA;
	      
	      double dxCA = xCr - xAr;
	      double dyCA = yCr - yAr;
	      
	      double dxyCA = sqrt( dxCA*dxCA + dyCA*dyCA );


	      //finn closest: update if closer
	      if( dxyCA < AMatchC.at(icA).distance ){
		AMatchC.at(icA).index = icC;
		AMatchC.at(icA).distance = dxyCA;
		//                AMatchC.at(iA).used = true;
		cout << "index cA " << icA << " index cC " << icC << " distance AC " << AMatchC.at(icA).distance << " dxyCA " << dxyCA << endl;
	      }

	      if( dxyCA < CMatchA.at(icC).distance ){
		CMatchA.at(icC).index = icA;
		CMatchA.at(icC).distance = dxyCA;
		cout << "index cA " << icA << " index cC " << icC << " distance CA " << CMatchA.at(icC).distance << " dxyCA " << dxyCA << endl;
	      }



	      
   
	      double xavg = 0.5 * ( xAr + xCr );
	      double yavg = 0.5 * ( yAr + yCr );
	      
	      for( vector<cluster>::iterator cB = vclB.begin(); cB != vclB.end(); ++cB )	{
	    
		  double xB = xcoordinate(1, cB, 0, ptchc, ptchr);
		  double yB = ycoordinate(1, cB, 0, ptchc, ptchr);

		  double dx3 = xB - xavg;
		  double dy3 = yB - yavg;

		  double dxy = (dx3*dx3+dy3*dy3);
		  //finn  closest:  update if closer

		  if(AMatchB.at(icA).used==false){
		    AMatchB.at(icA).index = icB;
		    AMatchB.at(icA).distance = dxy;
		    AMatchB.at(icA).used ==true;
		  }
		  else if(AMatchB.at(icA).used==true){
		    if(AMatchB.at(icA).distance> dxy){
		      AMatchB.at(icA).index = icB;
		      AMatchB.at(icA).distance = dxy;
		    }
		  }
			  
		  cout << "index cA " << icA << " index B " << icB << " distance A-B " << AMatchB.at(icA).distance << " dxy " << dxy << endl;
		  

		  if(BMatchA.at(icB).used==false){
		    BMatchA.at(icB).index = icA;
		    BMatchA.at(icB).distance = dxy;
		    BMatchA.at(icB).used=true;
		  }
		  else if(BMatchA.at(icB).used==true){
		    if(BMatchA.at(icB).distance> dxy){
		      BMatchA.at(icB).index = icB;
		      BMatchA.at(icB).distance = dxy;
		    }
		  }    
		  cout << "index cA " << icA << " index B " << icB << " distance B-A " << BMatchA.at(icB).distance << " dxy " << dxy << endl;
		  

		  if(CMatchB.at(icC).used==false){		    
		    CMatchB.at(icC).index = icB;
		    CMatchB.at(icC).distance = dxy;
		    CMatchB.at(icC).used =true;
		  }
		  else if(CMatchB.at(icC).used==true){
		    if(CMatchB.at(icC).distance> dxy){
		      CMatchB.at(icC).index = icB;
		      CMatchB.at(icC).distance = dxy;
		      }
		  }
		  cout << "index cC " << icC << " index B " << icB << " distance C-B " << CMatchB.at(icC).distance << " dxy " << dxy << endl;
		  
		  if(BMatchC.at(icB).used==false){
		      
		    BMatchC.at(icB).index = icC;
		    BMatchC.at(icB).distance = dxy;
		    BMatchC.at(icB).used=true;
		  }
		  else if(BMatchC.at(icB).used==true){
		    if(BMatchC.at(icB).distance> dxy){
		      BMatchC.at(icB).index = icB;
		      BMatchC.at(icB).distance = dxy;
		    }
		  }
		  
		  cout << "index cC " << icC << " index B " << icB << " distance B-C " << BMatchC.at(icB).distance << " dxy " << dxy << endl;
		  

		  icB++;
	      }//clB
	      icA++;   
	  }//clA
	  icC++;
      }//clC

      cout << " #########     summary  " << iev << endl;

      for( vector<cluster>::iterator cC = vclC.begin(); cC != vclC.end(); ++cC ){
	
	for( vector<cluster>::iterator cA = vclA.begin(); cA != vclA.end(); ++cA ){

	  for( vector<cluster>::iterator cB = vclB.begin(); cB != vclB.end(); ++cB ){

	    

	    cout << " index final IA " << IA << " AMatchC index " <<  AMatchC.at(IA).index << " distance dxyCA " << AMatchC.at(IA).distance << endl;
	    cout << " index final IC " << IC << " CMatchA index " <<  CMatchA.at(IC).index << " distance dxyCA " << CMatchA.at(IC).distance << endl;

	    cout << " index final IA " << IA << " AMatchB index " <<  AMatchB.at(IA).index << " distance dxy " << AMatchB.at(IA).distance << endl;
	    cout << " index final IB " << IB << " BMatchA index " <<  BMatchA.at(IB).index << " distance dxy " << BMatchA.at(IB).distance << endl;

	    cout << " index final IC " << IC << " CMatchB index " <<  CMatchB.at(IC).index << " distance dxy " << CMatchB.at(IC).distance << endl;
	    cout << " index final IB " << IB << " BMatchC index " <<  BMatchC.at(IB).index << " distance dxy " << BMatchC.at(IB).distance << endl;
	    IB++;
	  }
	  IA++;
	}
	IC++;
      }

      
      cout << " ########## DONE Matching for closest  ########## "<< endl;


      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ////////////////                A-B cluster correlations:    //////////////////////////
      if(PRINT) cout << "entering AB correlation loop " << endl;
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


	      double dxy = sqrt( dx*dx + dy*dy );
	      /*
	      //finn  closest:  update if closer
	      if( dxy < AMatchB.at(icA).distance ){
		AMatchB.at(icA).index = icB;
		AMatchB.at(icA).distance = dxy;
		cout << "index cA " << icA << " index cB " << icB << " distance A-B " << AMatchB.at(icA).distance << " dxy " << dxy << endl;
	      }

	      if( dxy < BMatchA.at(icB).distance ){
		BMatchA.at(icB).index = icA;
		BMatchA.at(icB).distance = dxy;
		cout << "index cA " << icA << " index cB " << icB << " distance B-A " << BMatchA.at(icB).distance << " dxy " << dxy << endl;
	      }
	      cout << " index A " << icA << " index cB " << icB <<  " distance dxy " << dxy << endl;

	      */
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
	      //icB++;
	    } // clusters
	  //icA++;
	} // cl
      
      nmvsevAB->Fill( iev, nm );
      
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // B-C cluster correlations:
      if(PRINT) cout << "entering BC correlation loop " << endl;
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

	      double dxy = sqrt( dx*dx + dy*dy );
	      /*
	      //finn  closest:  update if closer
	      if( dxy < CMatchB.at(icC).distance ){
		CMatchB.at(icC).index = icB2;
		CMatchB.at(icC).distance = dxy;
		cout << "index cC " << icC << " index cB2 " << icB2 << " distance C-B " << CMatchB.at(icC).distance << " dxy " << dxy << endl;
	      }

	      if( dxy < BMatchC.at(icB2).distance ){
		BMatchC.at(icB2).index = icC;
		BMatchC.at(icB2).distance = dxy;
		cout << "index cC " << icC << " index cB2 " << icB2 << " distance B-C " << BMatchC.at(icB2).distance << " dxy " << dxy << endl;
	      }
	      cout << " index cC " << icC << " index cB2 " << icB2 <<  " distance dxy " << dxy << endl;
	      */
	      
	      
	  if( cC->q > qL  && cC->q < qR && cB->q > qLB && cB->q < qRB && cC->iso && cB->iso )
	    {
	    
	      hdxCB->Fill( dx );
	      hdyCB->Fill( dy );
	      dxvsxCB->Fill( xB, dx );
	      dxvsyCB->Fill( yB, dx );
	      
	    }
	  
	  if( fabs( dx ) < straightTracks * beamDivergenceScaled + 0.020 && fabs( dy ) < straightTracks * beamDivergenceScaled + 0.100 )
	    ++nm;
	  //	  icC++;
	    } // clusters C
	  //icB2++;
	} // cl B
      
      nmvsevCB->Fill( iev, nm );
      
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      // ---------------------------------------
      ///////        A-C cluster correlations: ///////////
      if(PRINT) cout << "entering AC correlation loop " << endl;
      nm = 0;
      //      int iC =0;
      for( vector<cluster>::iterator cC = vclC.begin(); cC != vclC.end(); ++cC )
	{
	
	  double xC = xcoordinate(2, cC, alignxC, ptchc, ptchr);
	  double yC = ycoordinate(2, cC, alignyC, ptchc, ptchr);
	  
	  double xCr = xC*cfC - yC*sfC;
	  double yCr = xC*sfC + yC*cfC;
	
      
	  double etaC = eta(cC);
	  //          int iA=0;
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
	      /*
	      //finn closest: update if closer
	      if( dxyCA < AMatchC.at(iA).distance ){
		AMatchC.at(iA).index = iC;
		AMatchC.at(iA).distance = dxyCA;
		//                AMatchC.at(iA).used = true;
		cout << "index A " << iA << " index C " << iC << " distance AC " << AMatchC.at(iA).distance << " dxyCA " << dxyCA << endl;
	      }

	      if( dxyCA < CMatchA.at(iC).distance ){
		CMatchA.at(iC).index = iA;
		CMatchA.at(iC).distance = dxyCA;
		cout << "index A " << iA << " index C " << iC << " distance CA " << CMatchA.at(iC).distance << " dxyCA " << dxyCA << endl;
	      }

              cout << " index A " << iA << " index C " << iC << " distance dxyCA " << dxyCA << endl;

	      */
   
	      // NB until I initialize cB, I am using cA instead!!!!
	      fillControlHists(raw,"raw",0,0,cA,cA,cC,0,0,0,iev,0,0,xAr,yAr,xCr,yCr,dxCA,0,0,etaC,histoFile,fileName,hclphAiii,hclphBiii,hclphCiii,hclqAiii,hclqBiii,hclqCiii);


	      if(PRINT) cout << " tracks condition " << fabs( dxCA ) << " vs "  << straightTracks * beamDivergenceScaled + 0.02 << endl;

	      if( fabs( dxCA ) > straightTracks * beamDivergenceScaled + 0.02 ) continue; // includes beam divergence: +-5 sigma

	      fillControlHists(dxCAcut,"dxCAcut",0,0,cA,cA,cC,0,0,0,iev,0,0,xAr,yAr,xCr,yCr,dxCA,0,0,etaC,histoFile,fileName,hclphAiii,hclphBiii,hclphCiii,hclqAiii,hclqBiii,hclqCiii);


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

	      
	      //              int iB = 0; 
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

                  dx3vsx->Fill( xB, dx3 ); // turn
		  dx3vsy->Fill( yB, dx3 );
		  
		  
		  fillControlHists(beforeCorrections,"beforeCorrections",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName,hclphAiii,hclphBiii,hclphCiii,hclqAiii,hclqBiii,hclqCiii);

		  dx3 = dx3 - dx3corr*xavg; // from -dx3vsx.Fit("pol1")


		  double dxy = sqrt( dx3*dx3 + dy3*dy3 );

		  fillControlHists(nocuts,"nocuts",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName,hclphAiii,hclphBiii,hclphCiii,hclqAiii,hclqBiii,hclqCiii);
		  if(PRINT) cout << "nocuts" << endl;

		  /*
		  //finn  closest:  update if closer
		  if( dxy < AMatchB.at(iA).distance ){
		    AMatchB.at(iA).index = iB;
		    AMatchB.at(iA).distance = dxy;
		    cout << "index A " << iA << " index B " << iB << " distance A-B " << AMatchB.at(iA).distance << " dxy " << dxy << endl;
		  }

		  if( dxy < BMatchA.at(iB).distance ){
		    BMatchA.at(iB).index = iA;
		    BMatchA.at(iB).distance = dxy;
		    cout << "index A " << iA << " index B " << iB << " distance B-A " << BMatchA.at(iB).distance << " dxy " << dxy << endl;
		  }
		  cout << " index A " << iA << " index B " << iB <<  " distance dxy " << dxy << endl;

		  if( dxy < CMatchB.at(iC).distance ){
		    CMatchB.at(iC).index = iB;
		    CMatchB.at(iC).distance = dxy;
		    cout << "index C " << iC << " index B " << iB << " distance C-B " << CMatchB.at(iC).distance << " dxy " << dxy << endl;
		  }

		  if( dxy < BMatchC.at(iB).distance ){
		    BMatchC.at(iB).index = iC;
		    BMatchC.at(iB).distance = dxy;
		    cout << "index C " << iC << " index B " << iB << " distance B-C " << BMatchC.at(iB).distance << " dxy " << dxy << endl;
		  }
		  cout << " index C " << iC << " index B " << iB <<  " distance dxy " << dxy << endl;
*/

	
		    if( fabs( dy3 ) < straightTracks * beamDivergenceScaled ) //+ 0.05 )
		      { // cut on y, look at x, see madx3vsy
		    
			fillControlHists(straightTracksY,"straightTracksY",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName,hclphAiii,hclphBiii,hclphCiii,hclqAiii,hclqBiii,hclqCiii);
			
			if(AMatchC.at(iA).index == iC && CMatchA.at(iC).index == iA && AMatchB.at(iA).index == iB && BMatchA.at(iB).index == iA  && CMatchB.at(iC).index == iB && BMatchC.at(iB).index == iC){
			    //cout << "event "<< iev << " used distance AMatchC " << AMatchC.at(iA).distance << " CMatchA " << CMatchA.at(iC).distance << " ACMatchB " << ACMatchB.at(iA).distance << " BMatchAC " <<  BMatchAC.at(iB).distance << endl;
			    fillControlHists(closest3M,"closest3M",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName,hclphAiii,hclphBiii,hclphCiii,hclqAiii,hclqBiii,hclqCiii);

			    
			    if( cA->iso && cC->iso )
			      {
				fillControlHists(straightTracksY_isoAandC,"straightTracksY_isoAandC",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName,hclphAiii,hclphBiii,hclphCiii,hclqAiii,hclqBiii,hclqCiii);
				if( cB->iso )
				  {
				    fillControlHists(straightTracksY_isoAandCandB,"straightTracksY_isoAandCandB",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName,hclphAiii,hclphBiii,hclphCiii,hclqAiii,hclqBiii,hclqCiii);
				    if(PRINT) cout << "isoAandCandB done" << endl;
				

				    if( fabs( dxCA ) < straightTracks * beamDivergenceScaled )
				      { // track angle
					fillControlHists(straightTracksY_isoAandCandB_straightTracksX,"straightTracksY_isoAandCandB_straightTracksX",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName,hclphAiii,hclphBiii,hclphCiii,hclqAiii,hclqBiii,hclqCiii);
					cout << "event "<< iev << " used distance dxyCA " << dxyCA << " dxy " << dxy   << endl; //" iA " << iA << " iB " << iB << " iC " << iC << endl;
			    

				    
					if(PRINT) cout << " filling hists for resolution studies! " << endl;		  
					hclqAiii->Fill(cA->q);
					hclqBiii->Fill(cB->q); 
					hclqCiii->Fill(cC->q); 
					hclphAiii->Fill(cA->sum);
					hclphBiii->Fill(cB->sum);
					hclphCiii->Fill(cC->sum);
					if(PRINT) cout << "done  filling hists for resolution studies! " << endl;		  
					
					dx3tree = dx3;
				  
					clqAiii = cA->q;
					clqBiii = cB->q;
					clqCiii = cC->q;
					clphAiii = cA->sum;
					clphBiii = cB->sum;
					clphCiii = cC->sum;
					nrowBtree = nrowB;
					evt = iev;
					dxyCAtree=dxyCA;
					dxytree=dxy;
					charge_res->Fill();
					if(PRINT) cout << "done  filling the tree " << endl;
				
				  }//fabs( dxCA ) < straightTracks * beamDivergenceScaled 
				
			      }//iso B
			  }//iso CA
			}//closest
			
		      }//straight tracks 

		  iB++;
		}//clB
	      iA++;
	    }//clA
	  iC++;
	}//clC
    }//first big loop on events

  charge_res->Print();
  /*
  for( ; evinfoB != infoB.end() && evA != evlistA.end() && evB != evlistB.end() && evC != evlistC.end();       ++evinfoA, ++evinfoB, ++evinfoC, ++evA, ++evB, ++evC ){
    vector <cluster> vclA = *evA;
    vector <cluster> vclB = *evB;
    vector <cluster> vclC = *evC;

    ++iev;
    if( iev%10000 == 0 )
      cout << " " << iev << flush;

    if( evinfoA->filled == ADD ) cout << endl << "ev " << iev << " added A" << endl;
    if( evinfoB->filled == ADD ) cout << endl << "ev " << iev << " added B" << endl;
    if( evinfoC->filled == ADD ) cout << endl << "ev " << iev << " added C" << endl;
    
      for( vector<cluster>::iterator cC = vclC.begin(); cC != vclC.end(); ++cC ){
      for( vector<cluster>::iterator cA = vclA.begin(); cA != vclA.end(); ++cA ){
	      for( vector<cluster>::iterator cB = vclB.begin(); cB != vclB.end(); ++cB ){
	      if( fabs( dy3 ) < straightTracks * beamDivergenceScaled + 0.05 )
		    { 
		      
		      if( cA->iso && cC->iso )
			{
			  if( cB->iso )
			    {
			      if( (cA->q <= qR) && (cC->q <= qR))
				{
				  if(PRINT) cout << "cA->q " << cA->q << " >= qR " << qR << " && cC->q " << cC->q << "  >= qR " << endl;

				  if(PRINT) cout << "straightTracksY_isoAandCandB_chargerAandC going to be done" << endl;
			      
				  fillControlHists(straightTracksY_isoAandCandB_chargerAandC,"straightTracksY_isoAandCandB_chargerAandC",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName,hclphAiii,hclphBiii,hclphCiii,hclqAiii,hclqBiii,hclqCiii);
				  if(PRINT) cout << "straightTracksY_isoAandCandB_chargerAandC going done" << endl;
			      
				  if( ( cB->q <= qRB) && ( cA->q <= qR) && ( cC->q <= qR))
				    {
				      
				      fillControlHists(straightTracksY_isoAandCandB_chargerAandCandB,"straightTracksY_isoAandCandB_chargerAandCandB",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName,hclphAiii,hclphBiii,hclphCiii,hclqAiii,hclqBiii,hclqCiii);
				    }// ( cB->q >= qRB) && ( cA->q >= qR) && ( cC->q >= qR)

				  if( ( cA->q >= qL || cA->q <= qR) && ( cC->q >= qL || cC->q <= qR))
				    {
				      fillControlHists(straightTracksY_isoAandCandB_chargeAandC,"straightTracksY_isoAandCandB_chargeAandC",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName,hclphAiii,hclphBiii,hclphCiii,hclqAiii,hclqBiii,hclqCiii);

				      if( (cB->q >= qLB || cB->q <= qRB) && ( cA->q >= qL || cA->q <= qR) && ( cC->q >= qL || cC->q <= qR))
					{					  
					  fillControlHists(straightTracksY_isoAandCandB_chargeAandCandB,"straightTracksY_isoAandCandB_chargeAandCandB",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName,hclphAiii,hclphBiii,hclphCiii,hclqAiii,hclqBiii,hclqCiii);
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
		      
		      fillControlHists(hitsOnTrack,"hitsOnTrack",dx3,dy3,cA,cB,cC,nrowB,ncolB,xmod,iev,xB,yB,xAr,yAr,xCr,yCr,dxCA,etaA,etaB,etaC,histoFile,fileName,hclphAiii,hclphBiii,hclphCiii,hclqAiii,hclqBiii,hclqCiii);
				    

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
  */
  /////////    end of the big loop on the events!!! //////////////////
  cout << endl;
  if(PRINT) cout << " after big loop on events!!!************************" << endl;
  if(PRINT)cout << hdxAB->GetTitle() << " entries " << hdxAB->GetEntries() << endl;
  if(PRINT)cout << hdt->GetTitle() << " entries " << hdt->GetEntries() << endl;
  cout << "tree entries " << charge_res->GetEntries() << endl;  
  TH1I * hdx3tree = new TH1I("hdx3tree", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
  charge_res->Draw("dx3tree>>hdx3tree","","goff");
  hdx3tree = (TH1I*)gDirectory->Get("hdx3tree");
  cout << hdx3tree->GetTitle() << " entries " << hdx3tree->GetEntries() << endl;
  hdx3tree->Write();

  for(int i =0; i<hdx3tree->GetEntries(); i++)    {    
      hdx3tree2->SetBinContent(i+1,hdx3tree->GetBinContent(i+1));
    }
  cout << hdx3tree2->GetTitle() << " entries " << hdx3tree2->GetEntries() << endl;
  hdx3tree2->Write();


  // preparation to find the 90 % of the Landau with the lowest charge in all three 3M planes
  
  double integral[DreiMasterPlanes];
  double integralPH[DreiMasterPlanes];

  double integral90[DreiMasterPlanes];
  double integralPH90[DreiMasterPlanes];
  
  int high90[DreiMasterPlanes];
  int highPH90[DreiMasterPlanes];

  TH1I * hclq[DreiMasterPlanes];
  hclq[0]=hclqAiii;
  hclq[1]=hclqBiii;
  hclq[2]=hclqCiii;

  TH1I * hclph[DreiMasterPlanes];
  hclph[0]=hclphAiii;
  hclph[1]=hclphBiii;
  hclph[2]=hclphCiii;

  
  for(int j =0; j < DreiMasterPlanes; j++)   {
    if(PRINT)   cout << " plane " << j << endl;
    integral[j] = hclq[j]->Integral(0,hclq[j]->GetNbinsX()+1);
    if(PRINT)      cout << " integral " << integral[j] << " entries " << hclq[j]->GetEntries()<<  endl;
    
    integralPH[j] = hclph[j]->Integral(0,hclph[j]->GetNbinsX()+1);
    if(PRINT)   cout << " integralPH " << integralPH[j] << endl;
    
    integral90[j] = 0.9*integral[j]; //careful!! you should take into account the peak at low value! Hist is now filled only for isolated clusters and the effect is reduced
    if(DEBUG)   cout << " integral90 " << integral90[j] << endl;

    integralPH90[j] = 0.9*integralPH[j];

    high90[j] = 0;
    highPH90[j] = 0;
    
    int i = 0;
    while(integral[j]>integral90[j])    {
      cout << " while "<< i << endl;
      integral[j] = hclq[j]->Integral(1,hclq[j]->GetNbinsX()-i);
      high90[j] = hclq[j]->GetBinCenter(hclq[j]->GetNbinsX()-i);
      cout << " integral " << integral[j] << " high " << high90[j] << endl;
      i++;
    }
    cout << " integral90 " << integral90[j] << endl;
    i=0;
    while(integralPH[j]>integralPH90[j])  {
      integralPH[j] = hclph[j]->Integral(1,hclph[j]->GetNbinsX()-i);
      highPH90[j] = hclph[j]->GetBinCenter(hclph[j]->GetNbinsX()-i);
      i++;
    } 
  }

  // plot dx3 with above conditions
  //90 % Landau
  TString qA;
  qA.Form("%d",high90[0]);
  TString qB;
  qB.Form("%d",high90[1]);
  TString qC;
  qC.Form("%d",high90[2]);
  
  TH1D * hdx3treeq = new TH1D("hdx3treeq", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
  charge_res->Draw("dx3tree>>hdx3treeq","clqAiiitree<"+qA+"&&clqBiiitree<"+qB+"&&clqCiiitree<"+qC,"goff");
  hdx3treeq = (TH1D*)gDirectory->Get("hdx3treeq");
  cout << hdx3treeq->GetTitle() << " entries " << hdx3treeq->GetEntries() << " mean " << hdx3treeq->GetMean() << " 20bin " << hdx3treeq->GetBinContent(20)<<endl;
  hdx3treeq->Write();
  
  for(int i =0; i<hdx3treeq->GetEntries(); i++)    {
      hdx3_clchargeABC90evR->SetBinContent(i+1,hdx3treeq->GetBinContent(i+1));
  }
  cout << hdx3_clchargeABC90evR->GetTitle() << " entries " << hdx3_clchargeABC90evR->GetEntries() << "mean " << hdx3_clchargeABC90evR->GetMean() << endl;
  hdx3_clchargeABC90evR->Write();


  TString phA;
  phA.Form("%d",highPH90[0]);
  TString phB;
  phB.Form("%d",highPH90[1]);
  TString phC;
  phC.Form("%d",highPH90[2]);
  
  
  TH1I * hdx3treeph = new TH1I("hdx3treeph", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
  charge_res->Draw("dx3tree>>hdx3treeph","clphAiiitree<"+phA+"&&clphBiiitree<"+phB+"&&clphCiiitree<"+phC,"goff");
  hdx3treeph = (TH1I*)gDirectory->Get("hdx3treeph");
  cout << hdx3treeph->GetTitle() << " entries " << hdx3treeph->GetEntries() << endl;
  hdx3treeph->Write();
  
  for(int i =0; i<hdx3treeph->GetEntries(); i++)  {
    hdx3_clphABC90evR->SetBinContent(i+1,hdx3treeph->GetBinContent(i+1));
  }
  cout << hdx3_clphABC90evR->GetTitle() << " entries " << hdx3_clphABC90evR->GetEntries() << endl;
  hdx3_clphABC90evR->Write();

  TH1I * hdxyCAtreeph = new TH1I("hdxyCAtreeph", "triplet dxyCA ; dx [mm];triplets", 100,0., 0.3 );
  charge_res->Draw("dxyCAtree>>hdxyCAtreeph","clphAiiitree<"+phA+"&&clphBiiitree<"+phB+"&&clphCiiitree<"+phC,"goff");
  hdxyCAtreeph = (TH1I*)gDirectory->Get("hdxyCAtreeph");
  cout << hdxyCAtreeph->GetTitle() << " entries " << hdxyCAtreeph->GetEntries() << endl;
  hdxyCAtreeph->Write();
  
  for(int i =0; i<hdxyCAtreeph->GetEntries(); i++)  {
    hdxyCA_clphABC90evR->SetBinContent(i+1,hdxyCAtreeph->GetBinContent(i+1));
  }
  cout << hdxyCA_clphABC90evR->GetTitle() << " entries " << hdxyCA_clphABC90evR->GetEntries() << endl;
  hdxyCA_clphABC90evR->Write();

  TH1I * hdxytreeph = new TH1I("hdxytreeph", "triplet dxy ; dx [mm];triplets", 100, 0., 0.3 );
  charge_res->Draw("dxytree>>hdxytreeph","clphAiiitree<"+phA+"&&clphBiiitree<"+phB+"&&clphCiiitree<"+phC,"goff");
  hdxytreeph = (TH1I*)gDirectory->Get("hdxytreeph");
  cout << hdxytreeph->GetTitle() << " entries " << hdxytreeph->GetEntries() << endl;
  hdxytreeph->Write();
  
  for(int i =0; i<hdxytreeph->GetEntries(); i++)  {
    hdxy_clphABC90evR->SetBinContent(i+1,hdxytreeph->GetBinContent(i+1));
  }
  cout << hdxy_clphABC90evR->GetTitle() << " entries " << hdxy_clphABC90evR->GetEntries() << endl;
  hdxy_clphABC90evR->Write();


  
  double dlow, dhigh;
  TString ss_low,ss_high;

  cout << "charge 95 " << endl;

  getPercentRange(hdx3treeq, &(dlow),&(dhigh), 0.95,method);

  ss_low.Form("%f",dlow);
  ss_high.Form("%f",dhigh);

  TString ss_low_q95= ss_low;
  TString ss_high_q95 = ss_high;
  
  cout << " low and hig bin centers strings " << ss_low << " " << ss_high << endl;
  TH1D * hdx3treeq95 = new TH1D("hdx3treeq95", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
  charge_res->Draw("dx3tree>>hdx3treeq95","clqAiiitree<"+qA+"&&clqBiiitree<"+qB+"&&clqCiiitree<"+qC+"&&dx3tree<"+ss_high+"&&dx3tree>"+ss_low,"goff");
  hdx3treeq95 = (TH1D*)gDirectory->Get("hdx3treeq95");
  cout << hdx3treeq95->GetTitle() << " entries " << hdx3treeq95->GetEntries() << " mean " << hdx3treeq95->GetMean() << endl; 
  hdx3treeq95->Write();

  for(int i =0; i<hdx3treeq95->GetEntries(); i++)    {
    hdx3_clchargeABC90evR95->SetBinContent(i+1,hdx3treeq95->GetBinContent(i+1));
  }
  cout << hdx3_clchargeABC90evR95->GetTitle() << " entries " << hdx3_clchargeABC90evR95->GetEntries() << "mean " << hdx3_clchargeABC90evR95->GetMean() << endl;
  hdx3_clchargeABC90evR95->Write();

  TH1I * hnrowBtreeq95 = new TH1I("hnrowBtreeq95", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
  charge_res->Draw("nrowBtree>>hnrowBtreeq95","clqAiiitree<"+qA+"&&clqBiiitree<"+qB+"&&clqCiiitree<"+qC+"&&dx3tree<"+ss_high+"&&dx3tree>"+ss_low,"goff");
  hnrowBtreeq95 = (TH1I*)gDirectory->Get("hnrowBtreeq95");
  cout << hnrowBtreeq95->GetTitle() << " entries " << hnrowBtreeq95->GetEntries() << " mean " << hnrowBtreeq95->GetMean() << endl;
  hnrowBtreeq95->Write();

  for(int i =0; i<hnrowBtreeq95->GetEntries(); i++)    {
    hnrowB_q95->SetBinContent(i+1,hnrowBtreeq95->GetBinContent(i+1));
  }
  cout << hnrowB_q95->GetTitle() << " entries " << hnrowB_q95->GetEntries() << "mean " << hnrowB_q95->GetMean() << endl;
  hnrowB_q95->Write();

  TH1I * hnrowBtreeq5 = new TH1I("hnrowBtreeq5", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
  charge_res->Draw("nrowBtree>>hnrowBtreeq5","clqAiiitree<"+qA+"&&clqBiiitree<"+qB+"&&clqCiiitree<"+qC+"&&(dx3tree>"+ss_high+"||dx3tree<"+ss_low+")","goff");
  hnrowBtreeq5 = (TH1I*)gDirectory->Get("hnrowBtreeq5");
  cout << hnrowBtreeq5->GetTitle() << " entries " << hnrowBtreeq5->GetEntries() << " mean " << hnrowBtreeq5->GetMean() << endl;
  hnrowBtreeq5->Write();

  for(int i =0; i<hnrowBtreeq5->GetEntries(); i++)    {
    hnrowB_q5->SetBinContent(i+1,hnrowBtreeq5->GetBinContent(i+1));
  }
  cout << hnrowB_q5->GetTitle() << " entries " << hnrowB_q5->GetEntries() << "mean " << hnrowB_q5->GetMean() << endl;
  hnrowB_q5->Write();


  
  TH1I * hclBtreeq95 = new TH1I( "hclBtreeq95", "B cluster charge of 95% ;cluster charge [ke]; B clusters",            100, 0, 50 );  
  charge_res->Draw("clqBiiitree>>hclBtreeq95","clqAiiitree<"+qA+"&&clqBiiitree<"+qB+"&&clqCiiitree<"+qC+"&&dx3tree<"+ss_high+"&&dx3tree>"+ss_low,"goff");
  hclBtreeq95 = (TH1I*)gDirectory->Get("hclBtreeq95");
  cout << hclBtreeq95->GetTitle() << " entries " << hclBtreeq95->GetEntries() << " mean " << hclBtreeq95->GetMean() << endl;
  hclBtreeq95->Write();

  for(int i =0; i<hclBtreeq95->GetEntries(); i++)    {
    hclB_q95->SetBinContent(i+1,hclBtreeq95->GetBinContent(i+1));
  }
  cout << hclB_q95->GetTitle() << " entries " << hclB_q95->GetEntries() << "mean " << hclB_q95->GetMean() << endl;
  hclB_q95->Write();

  TH1I * hclBtreeq5 = new TH1I( "hclBtreeq5", "B cluster charge of 95% ;cluster charge [ke]; B clusters",            100, 0, 50 );  
  charge_res->Draw("clqBiiitree>>hclBtreeq5","clqAiiitree<"+qA+"&&clqBiiitree<"+qB+"&&clqCiiitree<"+qC+"&&(dx3tree>"+ss_high+"||dx3tree<"+ss_low+")","goff");
  hclBtreeq5 = (TH1I*)gDirectory->Get("hclBtreeq5");
  cout << hclBtreeq5->GetTitle() << " entries " << hclBtreeq5->GetEntries() << " mean " << hclBtreeq5->GetMean() << endl;
  hclBtreeq5->Write();

  for(int i =0; i<hclBtreeq5->GetEntries(); i++)    {
    hclB_q5->SetBinContent(i+1,hclBtreeq5->GetBinContent(i+1));
  }
  cout << hclB_q5->GetTitle() << " entries " << hclB_q5->GetEntries() << "mean " << hclB_q5->GetMean() << endl;
  hclB_q5->Write();


  cout << "ph 95 " << endl;

  getPercentRange(hdx3treeph, &(dlow),&(dhigh), 0.95,method);

  ss_low.Form("%f",dlow);
  ss_high.Form("%f",dhigh);
  TString ss_low_ph95= ss_low;
  TString ss_high_ph95 = ss_high;
  
  TH1D * hdx3treeph95 = new TH1D("hdx3treeph95", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
  charge_res->Draw("dx3tree>>hdx3treeph95","clphAiiitree<"+phA+"&&clphBiiitree<"+phB+"&&clphCiiitree<"+phC+"&&dx3tree<"+ss_high+"&&dx3tree>"+ss_low,"goff");
  hdx3treeph95 = (TH1D*)gDirectory->Get("hdx3treeph95");
  cout << hdx3treeph95->GetTitle() << " entries " << hdx3treeph95->GetEntries() << " mean " << hdx3treeph95->GetMean() << endl; 
  hdx3treeph95->Write();

  for(int i =0; i<hdx3treeph95->GetEntries(); i++)    {
    hdx3_clphABC90evR95->SetBinContent(i+1,hdx3treeph95->GetBinContent(i+1));
  }
  cout << hdx3_clphABC90evR95->GetTitle() << " entries " << hdx3_clphABC90evR95->GetEntries() << "mean " << hdx3_clphABC90evR95->GetMean() << endl;
  
  cout << "TO COMPARE resolution " <<  hdx3_clphABC90evR95->GetRMS() * 1000 << " ± " <<  hdx3_clphABC90evR95->GetRMSError() * 1000 << endl;
  
  hdx3_clphABC90evR95->Write();

  TH1I * hnrowBtreeph95 = new TH1I("hnrowBtreeph95", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
  charge_res->Draw("nrowBtree>>hnrowBtreeph95","clphAiiitree<"+phA+"&&clphBiiitree<"+phB+"&&clphCiiitree<"+phC+"&&dx3tree<"+ss_high+"&&dx3tree>"+ss_low,"goff");
  hnrowBtreeph95 = (TH1I*)gDirectory->Get("hnrowBtreeph95");
  cout << hnrowBtreeph95->GetTitle() << " entries " << hnrowBtreeph95->GetEntries() << " mean " << hnrowBtreeph95->GetMean() << endl;
  hnrowBtreeph95->Write();

  for(int i =0; i<hnrowBtreeph95->GetEntries(); i++)    {
    hnrowB_ph95->SetBinContent(i+1,hnrowBtreeph95->GetBinContent(i+1));
  }
  cout << hnrowB_ph95->GetTitle() << " entries " << hnrowB_ph95->GetEntries() << "mean " << hnrowB_ph95->GetMean() << endl;
  hnrowB_ph95->Write();

  TH1I * hnrowBtreeph5 = new TH1I("hnrowBtreeph5", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
  charge_res->Draw("nrowBtree>>hnrowBtreeph5","clphAiiitree<"+phA+"&&clphBiiitree<"+phB+"&&clphCiiitree<"+phC+"&&(dx3tree>"+ss_high+"||dx3tree<"+ss_low+")","goff");
  hnrowBtreeph5 = (TH1I*)gDirectory->Get("hnrowBtreeph5");
  cout << hnrowBtreeph5->GetTitle() << " entries " << hnrowBtreeph5->GetEntries() << " mean " << hnrowBtreeph5->GetMean() << endl;
  hnrowBtreeph5->Write();

  for(int i =0; i<hnrowBtreeph5->GetEntries(); i++)    {
    hnrowB_ph5->SetBinContent(i+1,hnrowBtreeph5->GetBinContent(i+1));
  }
  cout << hnrowB_ph5->GetTitle() << " entries " << hnrowB_ph5->GetEntries() << "mean " << hnrowB_ph5->GetMean() << endl;
  hnrowB_ph5->Write();


  
  TH1I * hclBtreeph95 = new TH1I( "hclBtreeph95", "B cluster charge of 95% ;cluster charge [ke]; B clusters",            200, 0, 1000 );  
  charge_res->Draw("clphBiiitree>>hclBtreeph95","clphAiiitree<"+phA+"&&clphBiiitree<"+phB+"&&clphCiiitree<"+phC+"&&dx3tree<"+ss_high+"&&dx3tree>"+ss_low,"goff");
  hclBtreeph95 = (TH1I*)gDirectory->Get("hclBtreeph95");
  cout << hclBtreeph95->GetTitle() << " entries " << hclBtreeph95->GetEntries() << " mean " << hclBtreeph95->GetMean() << endl;
  hclBtreeph95->Write();

  for(int i =0; i<hclBtreeph95->GetEntries(); i++)    {
    hclB_ph95->SetBinContent(i+1,hclBtreeph95->GetBinContent(i+1));
  }
  cout << hclB_ph95->GetTitle() << " entries " << hclB_ph95->GetEntries() << "mean " << hclB_ph95->GetMean() << endl;
  hclB_ph95->Write();

  TH1I * hdxyCAtreeph95 = new TH1I( "hdxyCAtreeph95", "B cluster charge of 95% ;cluster charge [ke]; B clusters",            100, 0.,0.3 );
  charge_res->Draw("dxyCAtree>>hdxyCAtreeph95","clphAiiitree<"+phA+"&&clphBiiitree<"+phB+"&&clphCiiitree<"+phC+"&&dx3tree<"+ss_high+"&&dx3tree>"+ss_low,"goff");
  hdxyCAtreeph95 = (TH1I*)gDirectory->Get("hdxyCAtreeph95");
  cout << hdxyCAtreeph95->GetTitle() << " entries " << hdxyCAtreeph95->GetEntries() << " mean " << hdxyCAtreeph95->GetMean() << endl;
  hdxyCAtreeph95->Write();

  for(int i =0; i<hdxyCAtreeph95->GetEntries(); i++)    {
    hdxyCA_ph95->SetBinContent(i+1,hdxyCAtreeph95->GetBinContent(i+1));
  }
  cout << hdxyCA_ph95->GetTitle() << " entries " << hdxyCA_ph95->GetEntries() << "mean " << hdxyCA_ph95->GetMean() << endl;
  hdxyCA_ph95->Write();
  
  TH1I * hdxytreeph95 = new TH1I( "hdxytreeph95", "B cluster charge of 95% ;cluster charge [ke]; B clusters",            100, 0, 0.3 );
  charge_res->Draw("dxytree>>hdxytreeph95","clphAiiitree<"+phA+"&&clphBiiitree<"+phB+"&&clphCiiitree<"+phC+"&&dx3tree<"+ss_high+"&&dx3tree>"+ss_low,"goff");
  hdxytreeph95 = (TH1I*)gDirectory->Get("hdxytreeph95");
  cout << hdxytreeph95->GetTitle() << " entries " << hdxytreeph95->GetEntries() << " mean " << hdxytreeph95->GetMean() << endl;
  hdxytreeph95->Write();

  for(int i =0; i<hdxytreeph95->GetEntries(); i++)    {
    hdxy_ph95->SetBinContent(i+1,hdxytreeph95->GetBinContent(i+1));
  }
  cout << hdxy_ph95->GetTitle() << " entries " << hdxy_ph95->GetEntries() << "mean " << hdxy_ph95->GetMean() << endl;
  hdxy_ph95->Write();


  TH1I * hclBtreeph5 = new TH1I( "hclBtreeph5", "B cluster charge of 95% ;cluster charge [ke]; B clusters",            200, 0, 1000 );  
  charge_res->Draw("clphBiiitree>>hclBtreeph5","clphAiiitree<"+phA+"&&clphBiiitree<"+phB+"&&clphCiiitree<"+phC+"&&(dx3tree>"+ss_high+"||dx3tree<"+ss_low+")","goff");
  hclBtreeph5 = (TH1I*)gDirectory->Get("hclBtreeph5");
  cout << hclBtreeph5->GetTitle() << " entries " << hclBtreeph5->GetEntries() << " mean " << hclBtreeph5->GetMean() << endl;
  hclBtreeph5->Write();

  for(int i =0; i<hclBtreeph5->GetEntries(); i++)    {
    hclB_ph5->SetBinContent(i+1,hclBtreeph5->GetBinContent(i+1));
  }
  cout << hclB_ph5->GetTitle() << " entries " << hclB_ph5->GetEntries() << "mean " << hclB_ph5->GetMean() << endl;
  hclB_ph5->Write();


  cout << "charge 99 " << endl;

  getPercentRange(hdx3treeq, &(dlow),&(dhigh), 0.9999,method);

  ss_low.Form("%f",dlow);
  ss_high.Form("%f",dhigh);
  TString ss_low_q99= ss_low;
  TString ss_high_q99 = ss_high;
  
  cout << " low and hig bin centers strings " << ss_low << " " << ss_high << endl;
  TH1D * hdx3treeq99 = new TH1D("hdx3treeq99", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
  charge_res->Draw("dx3tree>>hdx3treeq99","clqAiiitree<"+qA+"&&clqBiiitree<"+qB+"&&clqCiiitree<"+qC+"&&dx3tree<"+ss_high+"&&dx3tree>"+ss_low,"goff");
  hdx3treeq99 = (TH1D*)gDirectory->Get("hdx3treeq99");
  cout << hdx3treeq99->GetTitle() << " entries " << hdx3treeq99->GetEntries() << " mean " << hdx3treeq99->GetMean() << endl; 
  hdx3treeq99->Write();

  for(int i =0; i<hdx3treeq99->GetEntries(); i++)    {
    hdx3_clchargeABC90evR99->SetBinContent(i+1,hdx3treeq99->GetBinContent(i+1));
  }
  cout << hdx3_clchargeABC90evR99->GetTitle() << " entries " << hdx3_clchargeABC90evR99->GetEntries() << "mean " << hdx3_clchargeABC90evR99->GetMean() << endl;
  hdx3_clchargeABC90evR99->Write();

  TH1I * hnrowBtreeq99 = new TH1I("hnrowBtreeq99", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
  charge_res->Draw("nrowBtree>>hnrowBtreeq99","clqAiiitree<"+qA+"&&clqBiiitree<"+qB+"&&clqCiiitree<"+qC+"&&dx3tree<"+ss_high+"&&dx3tree>"+ss_low,"goff");
  hnrowBtreeq99 = (TH1I*)gDirectory->Get("hnrowBtreeq99");
  cout << hnrowBtreeq99->GetTitle() << " entries " << hnrowBtreeq99->GetEntries() << " mean " << hnrowBtreeq99->GetMean() << endl;
  hnrowBtreeq99->Write();

  for(int i =0; i<hnrowBtreeq99->GetEntries(); i++)    {
    hnrowB_q99->SetBinContent(i+1,hnrowBtreeq99->GetBinContent(i+1));
  }
  cout << hnrowB_q99->GetTitle() << " entries " << hnrowB_q99->GetEntries() << "mean " << hnrowB_q99->GetMean() << endl;
  hnrowB_q99->Write();

  TH1I * hnrowBtreeq1 = new TH1I("hnrowBtreeq1", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
  charge_res->Draw("nrowBtree>>hnrowBtreeq1","clqAiiitree<"+qA+"&&clqBiiitree<"+qB+"&&clqCiiitree<"+qC+"&&(dx3tree>"+ss_high+"||dx3tree<"+ss_low+")","goff");
  hnrowBtreeq1 = (TH1I*)gDirectory->Get("hnrowBtreeq5");
  cout << hnrowBtreeq5->GetTitle() << " entries " << hnrowBtreeq1->GetEntries() << " mean " << hnrowBtreeq1->GetMean() << endl;
  hnrowBtreeq1->Write();

  for(int i =0; i<hnrowBtreeq1->GetEntries(); i++)    {
    hnrowB_q1->SetBinContent(i+1,hnrowBtreeq1->GetBinContent(i+1));
  }
  cout << hnrowB_q1->GetTitle() << " entries " << hnrowB_q1->GetEntries() << "mean " << hnrowB_q1->GetMean() << endl;
  hnrowB_q1->Write();


  
  TH1I * hclBtreeq99 = new TH1I( "hclBtreeq99", "B cluster charge of 99% ;cluster charge [ke]; B clusters",            100, 0, 50 );  
  charge_res->Draw("clqBiiitree>>hclBtreeq99","clqAiiitree<"+qA+"&&clqBiiitree<"+qB+"&&clqCiiitree<"+qC+"&&dx3tree<"+ss_high+"&&dx3tree>"+ss_low,"goff");
  hclBtreeq99 = (TH1I*)gDirectory->Get("hclBtreeq99");
  cout << hclBtreeq99->GetTitle() << " entries " << hclBtreeq99->GetEntries() << " mean " << hclBtreeq99->GetMean() << endl;
  hclBtreeq99->Write();

  for(int i =0; i<hclBtreeq99->GetEntries(); i++)    {
    hclB_q99->SetBinContent(i+1,hclBtreeq99->GetBinContent(i+1));
  }
  cout << hclB_q99->GetTitle() << " entries " << hclB_q99->GetEntries() << "mean " << hclB_q99->GetMean() << endl;
  hclB_q99->Write();

  TH1I * hclBtreeq1 = new TH1I( "hclBtreeq1", "B cluster charge of 99% ;cluster charge [ke]; B clusters",            100, 0, 50 );  
  charge_res->Draw("clqBiiitree>>hclBtreeq1","clqAiiitree<"+qA+"&&clqBiiitree<"+qB+"&&clqCiiitree<"+qC+"&&(dx3tree>"+ss_high+"||dx3tree<"+ss_low+")","goff");
  hclBtreeq1 = (TH1I*)gDirectory->Get("hclBtreeq1");
  cout << hclBtreeq1->GetTitle() << " entries " << hclBtreeq1->GetEntries() << " mean " << hclBtreeq1->GetMean() << endl;
  hclBtreeq1->Write();

  for(int i =0; i<hclBtreeq1->GetEntries(); i++)    {
    hclB_q1->SetBinContent(i+1,hclBtreeq1->GetBinContent(i+1));
  }
  cout << hclB_q1->GetTitle() << " entries " << hclB_q1->GetEntries() << "mean " << hclB_q1->GetMean() << endl;
  hclB_q1->Write();


  cout << "ph 99 " << endl;

  getPercentRange(hdx3treeph, &(dlow),&(dhigh), 0.9999,method);

  ss_low.Form("%f",dlow);
  ss_high.Form("%f",dhigh);
  TString ss_low_ph99= ss_low;
  TString ss_high_ph99 = ss_high;
  
  TH1D * hdx3treeph99 = new TH1D("hdx3treeph99", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
  charge_res->Draw("dx3tree>>hdx3treeph99","clphAiiitree<"+phA+"&&clphBiiitree<"+phB+"&&clphCiiitree<"+phC+"&&dx3tree<"+ss_high+"&&dx3tree>"+ss_low,"goff");
  hdx3treeph99 = (TH1D*)gDirectory->Get("hdx3treeph99");
  cout << hdx3treeph99->GetTitle() << " entries " << hdx3treeph99->GetEntries() << " mean " << hdx3treeph99->GetMean() << endl; 
  hdx3treeph99->Write();

  for(int i =0; i<hdx3treeph99->GetEntries(); i++)    {
    hdx3_clphABC90evR99->SetBinContent(i+1,hdx3treeph99->GetBinContent(i+1));
  }
  cout << hdx3_clphABC90evR99->GetTitle() << " entries " << hdx3_clphABC90evR99->GetEntries() << "mean " << hdx3_clphABC90evR99->GetMean() << endl;
  cout << "TO COMPARE resolution " <<  hdx3_clphABC90evR99->GetRMS() * 1000 << " ± " <<  hdx3_clphABC90evR99->GetRMSError() * 1000 << endl;
  hdx3_clphABC90evR99->Write();

  TH1I * hnrowBtreeph99 = new TH1I("hnrowBtreeph99", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
  charge_res->Draw("nrowBtree>>hnrowBtreeph99","clphAiiitree<"+phA+"&&clphBiiitree<"+phB+"&&clphCiiitree<"+phC+"&&dx3tree<"+ss_high+"&&dx3tree>"+ss_low,"goff");
  hnrowBtreeph99 = (TH1I*)gDirectory->Get("hnrowBtreeph99");
  cout << hnrowBtreeph99->GetTitle() << " entries " << hnrowBtreeph99->GetEntries() << " mean " << hnrowBtreeph99->GetMean() << endl;
  hnrowBtreeph99->Write();

  for(int i =0; i<hnrowBtreeph99->GetEntries(); i++)    {
    hnrowB_ph99->SetBinContent(i+1,hnrowBtreeph99->GetBinContent(i+1));
  }
  cout << hnrowB_ph99->GetTitle() << " entries " << hnrowB_ph99->GetEntries() << "mean " << hnrowB_ph99->GetMean() << endl;
  hnrowB_ph99->Write();


  TH1I * hdxyCAtreeph99 = new TH1I( "hdxyCAtreeph99", "B cluster charge of 95% ;cluster charge [ke]; B clusters",            100, 0, 0.3 );
  charge_res->Draw("dxyCAtree>>hdxyCAtreeph99","clphAiiitree<"+phA+"&&clphBiiitree<"+phB+"&&clphCiiitree<"+phC+"&&dx3tree<"+ss_high+"&&dx3tree>"+ss_low,"goff");
  hdxyCAtreeph99 = (TH1I*)gDirectory->Get("hdxyCAtreeph99");
  cout << hdxyCAtreeph99->GetTitle() << " entries " << hdxyCAtreeph99->GetEntries() << " mean " << hdxyCAtreeph99->GetMean() << endl;
  hdxyCAtreeph99->Write();

  for(int i =0; i<hdxyCAtreeph99->GetEntries(); i++)    {
    hdxyCA_ph99->SetBinContent(i+1,hdxyCAtreeph99->GetBinContent(i+1));
  }
  cout << hdxyCA_ph99->GetTitle() << " entries " << hdxyCA_ph99->GetEntries() << "mean " << hdxyCA_ph99->GetMean() << endl;
  hdxyCA_ph99->Write();
  
  TH1I * hdxytreeph99 = new TH1I( "hdxytreeph99", "B cluster charge of 95% ;cluster charge [ke]; B clusters",            100, 0, 0.3 );
  charge_res->Draw("dxytree>>hdxytreeph99","clphAiiitree<"+phA+"&&clphBiiitree<"+phB+"&&clphCiiitree<"+phC+"&&dx3tree<"+ss_high+"&&dx3tree>"+ss_low,"goff");
  hdxytreeph99 = (TH1I*)gDirectory->Get("hdxytreeph99");
  cout << hdxytreeph99->GetTitle() << " entries " << hdxytreeph99->GetEntries() << " mean " << hdxytreeph99->GetMean() << endl;
  hdxytreeph99->Write();

  for(int i =0; i<hdxytreeph99->GetEntries(); i++)    {
    hdxy_ph99->SetBinContent(i+1,hdxytreeph99->GetBinContent(i+1));
  }
  cout << hdxy_ph99->GetTitle() << " entries " << hdxy_ph99->GetEntries() << "mean " << hdxy_ph99->GetMean() << endl;
  hdxy_ph99->Write();

  
  TH1I * hnrowBtreeph1 = new TH1I("hnrowBtreeph1", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
  charge_res->Draw("nrowBtree>>hnrowBtreeph1","clphAiiitree<"+phA+"&&clphBiiitree<"+phB+"&&clphCiiitree<"+phC+"&&(dx3tree>"+ss_high+"||dx3tree<"+ss_low+")","goff");
  hnrowBtreeph1 = (TH1I*)gDirectory->Get("hnrowBtreeph1");
  cout << hnrowBtreeph1->GetTitle() << " entries " << hnrowBtreeph1->GetEntries() << " mean " << hnrowBtreeph1->GetMean() << endl;
  hnrowBtreeph1->Write();

  for(int i =0; i<hnrowBtreeph1->GetEntries(); i++)    {
    hnrowB_ph1->SetBinContent(i+1,hnrowBtreeph1->GetBinContent(i+1));
  }
  cout << hnrowB_ph1->GetTitle() << " entries " << hnrowB_ph1->GetEntries() << "mean " << hnrowB_ph1->GetMean() << endl;
  hnrowB_ph1->Write();


  
  TH1I * hclBtreeph99 = new TH1I( "hclBtreeph99", "B cluster charge of 99% ;cluster charge [ke]; B clusters",            200, 0, 1000 );  
  charge_res->Draw("clphBiiitree>>hclBtreeph99","clphAiiitree<"+phA+"&&clphBiiitree<"+phB+"&&clphCiiitree<"+phC+"&&dx3tree<"+ss_high+"&&dx3tree>"+ss_low,"goff");
  hclBtreeph99 = (TH1I*)gDirectory->Get("hclBtreeph99");
  cout << hclBtreeph99->GetTitle() << " entries " << hclBtreeph99->GetEntries() << " mean " << hclBtreeph99->GetMean() << endl;
  hclBtreeph99->Write();

  for(int i =0; i<hclBtreeph99->GetEntries(); i++)    {
    hclB_ph99->SetBinContent(i+1,hclBtreeph99->GetBinContent(i+1));
  }
  cout << hclB_ph99->GetTitle() << " entries " << hclB_ph99->GetEntries() << "mean " << hclB_ph99->GetMean() << endl;
  hclB_ph99->Write();

  TH1I * hclBtreeph1 = new TH1I( "hclBtreeph1", "B cluster charge of 99% ;cluster charge [ke]; B clusters",            200, 0, 1000 );  
  charge_res->Draw("clphBiiitree>>hclBtreeph1","clphAiiitree<"+phA+"&&clphBiiitree<"+phB+"&&clphCiiitree<"+phC+"&&(dx3tree>"+ss_high+"||dx3tree<"+ss_low+")","goff");
  hclBtreeph1 = (TH1I*)gDirectory->Get("hclBtreeph1");
  cout << hclBtreeph1->GetTitle() << " entries " << hclBtreeph1->GetEntries() << " mean " << hclBtreeph1->GetMean() << endl;
  hclBtreeph1->Write();

  for(int i =0; i<hclBtreeph1->GetEntries(); i++)    {
    hclB_ph1->SetBinContent(i+1,hclBtreeph1->GetBinContent(i+1));
  }
  cout << hclB_ph1->GetTitle() << " entries " << hclB_ph1->GetEntries() << "mean " << hclB_ph1->GetMean() << endl;
  hclB_ph1->Write();








  

  
  //charge_res->Print();

  cout << " Cluster charge cut for 90 % events A: " << qA << " B: " << qB << " C: " << qC << endl;
  cout << " 95 % residual is between " << ss_low_q95 << " and " << ss_high_q95 << endl; 
  cout << " 99 % residual is between " << ss_low_q99 << " and " << ss_high_q99 << endl; 

  cout << " Cluster ph cut for 90 % events A: " << phA << " B: " << phB << " C: " << phC << endl;
  cout << " 95 % residual is between " << ss_low_ph95 << " and " << ss_high_ph95 << endl; 
  cout << " 99 % residual is between " << ss_low_ph99 << " and " << ss_high_ph99 << endl; 





  
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
  histoFile->Write("",TObject::kOverwrite);
  if(PRINT) cout << " wrote the hist file" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // alignment fits:

  if(DOALIGNMENT)
    {
      cout << "**************** alignment iteration " << aligniteration << endl;
      cout << "******* ALIGNA A ********" << endl;
      double newalignxA = alignxA;
      if(PRINT) cout << " newalignxA " << newalignxA << " alignxA " << alignxA << endl;
      
      if(PRINT) cout << hdxAB->GetTitle() << " entries " << hdxAB->GetEntries() << endl;
      
      if( hdxAB->GetEntries() > 999 )
	{
      cout << "******* A x ********" << endl;
	  
      newalignxA += alignx(hdxAB,"A",runnum,aligniteration);
	}
      else
	cout << "  not enough statistics" << endl;
      
      // y:

      double newalignyA = alignyA;
      
      if(PRINT) cout << endl << hdyAB->GetTitle() << " entries " << hdyAB->GetEntries() << endl;
      if( hdyAB->GetEntries() > 999 ) {
	//cout << "  y correction " << hdyAB->GetMean() << endl;
	cout << "******* A y ********" << endl;

	newalignyA += aligny(hdyAB,"A",runnum,aligniteration);
      }
      else
	cout << "  not enough statistics" << endl;
      
      // dxvsy -> -rot
      
      double newalignfA = alignfA;
      
      if(PRINT) cout << endl << dxvsyAB->GetTitle() << " entries " << dxvsyAB->GetEntries() << endl;
      if( aligniteration > 0 && dxvsyAB->GetEntries() > 999 ) {
	cout << "******* A angle ********" << endl;
	newalignfA += alignangle(dxvsyAB,"A",runnum,aligniteration);
      }
      else
	cout << "  not enough" << endl;
      
      // C:

      cout << "******* ALIGNA C ********" << endl;
      double newalignxC = alignxC;
  
      if(PRINT) cout << endl << hdxCB->GetTitle() << " entries " << hdxCB->GetEntries() << endl;
      if( hdxCB->GetEntries() > 999 ) {
	cout << "******* C x ********" << endl;
	newalignxC += alignx(hdxCB,"C",runnum,aligniteration);
      }
      else
	cout << "  not enough" << endl;

      // y:
      
      double newalignyC = alignyC;

      if(PRINT)      cout << endl << hdyCB->GetTitle() << " entries " << hdyCB->GetEntries() << endl;
      if( hdyCB->GetEntries() > 999 ) {
	//	cout << "  y correction " << hdyCB->GetMean() << endl;
	cout << "******* C y ********" << endl;
	newalignyC += aligny(hdyCB,"C",runnum,aligniteration);
      }
      else
	cout << "  not enough" << endl;
      
      // dxvsy -> -rot

      double newalignfC = alignfC;
      
      cout << endl << dxvsyCB->GetTitle() << " entries " << dxvsyCB->GetEntries() << endl;
      if( aligniteration > 0 && dxvsyCB->GetEntries() > 999 ) {
	cout << "******* C angle ********" << endl;
	newalignfC += alignangle(dxvsyCB,"C",runnum,aligniteration);
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

void getPercentRange(TH1 * h,double * dlow, double * dhigh, double percent, TString method="mine"){
  cout << " initial entries " << h->GetEntries() << endl;
  cout << " percent "<< percent << " -> " <<  percent*h->GetEntries() << endl;
  int low = 0;
  int high = 0;


  if(method =="mine"){
    double tolerance = (h->GetBinContent(0)+h->GetBinContent(h->GetNbinsX()+1))/h->GetEntries(); // since we are not working with continuous quantities but with binned hists, the difference between the integral and the 95% it may not be zero, so we ask it to be lower than 15% (before I was using 1.1% but since I changed to use the full range of hists, the overflow bin plays a too bigger role when stats is low and the 1.1% is not reached. 
    cout << " bin 0 " << h->GetBinContent(0) << " overflow " << h->GetBinContent(h->GetNbinsX()+1) << " entries " << h->GetEntries() << endl;
    cout << " RMS method charge based with tolerance " << tolerance << " with over/under flow would be " << (h->GetBinContent(0)+h->GetBinContent(h->GetNbinsX()+1))/h->GetEntries() << endl;
    double maximum = h->GetMaximum();
    int maxbin = h->GetMaximumBin();
    cout <<  " maxbin " << maxbin << " at " << h->GetBinCenter(maxbin)<<endl;
    cout << "intital sigma = " << h->GetRMS() * 1000 << " sigmaerr = " << h->GetRMSError() * 1000 << endl;
    double Integral = h->Integral(0,h->GetNbinsX()+1);
    cout << " integral " << Integral << " entries " << h->GetEntries() << endl;
    double integral95 = percent*Integral;
    cout << " integra "<< percent << " " << integral95 << endl;
    int i = 0;
    
    for(int i =0; i<h->GetNbinsX()/2; i++)
      {
	low = maxbin-i;
	high = maxbin+i;
      
	Integral = h->Integral(low,high);
	cout << " integral " << Integral << "low " <<low << " high " << high << endl;
	cout << " while "<< i << " fabs(integral-integral95)/integral95 " << fabs(Integral-integral95)/integral95 << endl;
	
	//      if(fabs(Integral-integral95)/integral95 < 0.011 || Integral>integral95)//integral>integral95)
	if(fabs(Integral-integral95)/integral95 < tolerance || Integral>integral95)//integral>integral95)
	  break;
	cout << "final integral " << Integral << "low " <<low <<" high " << high << endl;

  

      }
  }
  if(method == "simon")
    {
      double integral = h->GetEntries();
      int maxbin = h->GetMaximumBin();
      cout << "entries=" << integral << " maxbin=" << maxbin << endl;

      double subrange_integral = h->GetBinContent(maxbin);
      int bin = 0;
      while(subrange_integral < percent*integral) {
	bin++;
	// Add one bin to the left:
	subrange_integral += h->GetBinContent(maxbin-bin);
	// Add one bin to the right:
	subrange_integral += h->GetBinContent(maxbin+bin);
	cout << "subrange " << (maxbin-bin) << "-" << (maxbin+bin) << ": entries=" << subrange_integral << endl;
      }
      cout << "subrange " << (maxbin-bin) << "-" << (maxbin+bin) << " now has " << subrange_integral << " entries, this is " << (100.0*subrange_integral)/integral << "%" << endl;

      // Correct by overshoot bin:
      subrange_integral -= h->GetBinContent(maxbin+bin);
      subrange_integral -= h->GetBinContent(maxbin-bin);
      bin--;

      low = maxbin-bin;
      high = maxbin+bin;
      cout << "subrange " << (maxbin-bin) << "-" << (maxbin+bin) << " now has " << subrange_integral << " entries, this is " << (100.0*subrange_integral)/integral << "%" << endl;


    }

  h->GetXaxis()->SetRange(low,high); //to restrict range to bins binlow to binhigh

  
  double sigma = h->GetRMS() * 1000;
  double sigmaerr = h->GetRMSError() * 1000;
  cout << " resolution " << sigma << " ± " << sigmaerr << endl;

  cout << "low and hig bin centers strings " << h->GetBinCenter(low) << " " << h->GetBinCenter(high) << endl;
  *dlow  = h->GetBinCenter(low);
  *dhigh = h->GetBinCenter(high);
  

}



double alignx(TH1I * h, TString plane,TString run,int iteration)
{
  TString iter;
  iter.Form("%d",iteration);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas("c","c",700,700);
  c->cd();
  c->SetFrameFillStyle(1000);
  c->SetFrameFillColor(0);
  gPad->SetTicks(1,1);

  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.025);
  h->GetXaxis()->SetTitleSize(0.035);
  h->GetXaxis()->SetTitleOffset(0.8);
  h->GetXaxis()->SetTitleFont(42);
  
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.025);
  h->GetYaxis()->SetTitleSize(0.035);
  h->GetYaxis()->SetTitleOffset(1.4);
  h->GetYaxis()->SetTitleFont(42);

  h->SetMarkerStyle(20);
  h->Draw("PZ");
      
  //important part
  TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -10, 10 );
  double xpk = h->GetBinCenter( h->GetMaximumBin() );
  fgp0->SetParameter( 0, h->GetMaximum() ); // amplitude
  fgp0->SetParameter( 1, xpk ); //mean
  fgp0->SetParameter( 2, h->GetBinWidth(1) ); // sigma
  fgp0->SetParameter( 3, h->GetBinContent( hdxAB->FindBin(xpk-1) ) ); // BG
  fgp0->SetParName(0, "amplitude");
  fgp0->SetParName(1, "mean");
  fgp0->SetParName(2, "sigma");
  fgp0->SetParName(3, "BG");
  h->Fit( "fgp0", "qI", "", xpk-1, xpk+1 ); // fit range around peak
  fgp0->SetLineColor(kRed);
  fgp0->SetLineWidth(2);
  fgp0->Draw("same");
  
  cout << "  Fit Gauss + BG:"
       << endl << "  area " << fgp0->GetParameter(0)
       << endl << "  mean " << fgp0->GetParameter(1)
       << endl << "  mean error " << fgp0->GetParError(1)
       << endl << "  sigm " << fgp0->GetParameter(2)
       << endl << "  offs " << fgp0->GetParameter(3)
       << endl;

  h->GetXaxis()->SetRangeUser(-0.5,0.5);
  c->Update();
  gStyle->SetOptFit(1111);
  TString outputDir="/home/zoiirene/Output/alignment/";

  //  TString outputDir = outputDir+"";
  TString outputFile = outputDir+"alignx_"+plane+"_run"+run+"_iteration_"+iter;
  c->SaveAs(outputFile+".eps");
  c->SaveAs(outputFile+".pdf");
  c->SaveAs(outputFile+".png");
  c->SaveAs(outputFile+".root");
  
  return fgp0->GetParameter(1);
}


double aligny(TH1I * h, TString plane,TString run,int iteration)
{
  TString iter;
  iter.Form("%d",iteration);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1);

  TCanvas * c = new TCanvas("c","c",700,700);
  c->cd();
  c->SetFrameFillStyle(1000);
  c->SetFrameFillColor(0);
  gPad->SetTicks(1,1);

  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.025);
  h->GetXaxis()->SetTitleSize(0.035);
  h->GetXaxis()->SetTitleOffset(0.8);
  h->GetXaxis()->SetTitleFont(42);

  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.025);
  h->GetYaxis()->SetTitleSize(0.035);
  h->GetYaxis()->SetTitleOffset(1.4);
  h->GetYaxis()->SetTitleFont(42);

  h->SetLineWidth(2);

  //important part
  cout << "  Getting hist mean:"
       << endl << "  mean " << h->GetMean()
       << endl << "  mean error " << h->GetMeanError()
       << endl;
  double mean = h->GetMean();
  h->GetXaxis()->SetRangeUser(-0.5,0.5);
  h->Draw("hist");


  TString outputDir="/home/zoiirene/Output/alignment/";

  //  TString outputDir = outputDir+"";
  TString outputFile = outputDir+"aligny_"+plane+"_run"+run+"_iteration_"+iter;
  c->SaveAs(outputFile+".eps");
  c->SaveAs(outputFile+".pdf");
  c->SaveAs(outputFile+".png");
  c->SaveAs(outputFile+".root");

  
  return mean;



  
}

double alignangle(TProfile * h, TString plane,TString run,int iteration)
{
  TString iter;
  iter.Form("%d",iteration);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas("c","c",700,700);
  c->cd();
  c->SetFrameFillStyle(1000);
  c->SetFrameFillColor(0);
  gPad->SetTicks(1,1);


  h->SetTitle("");
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.025);
  h->GetXaxis()->SetTitleSize(0.035);
  h->GetXaxis()->SetTitleOffset(0.8);
  h->GetXaxis()->SetTitleFont(42);
  //  h->GetXaxis()->SetRangeUser(-0.5,0.5);

  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.025);
  h->GetYaxis()->SetTitleSize(0.035);
  h->GetYaxis()->SetTitleOffset(1.4);
  h->GetYaxis()->SetTitleFont(42);

  h->SetMarkerStyle(20);
  h->Draw("PZ");

  //important part
  h->Fit( "pol1", "q", "", -3, 3 );
  TF1 * fdxvsy = h->GetFunction( "pol1" );
  cout << "  Linear Fit:"
       << endl << "  slope " << fdxvsy->GetParameter(1)
       << endl << "  slope error " << fdxvsy->GetParError(1)
       << endl << "  offs " << fdxvsy->GetParameter(0)
       << endl;

  fdxvsy->SetLineColor(kRed);
  fdxvsy->SetLineWidth(2);
  fdxvsy->Draw("same");
  
  
  c->Update();
  gStyle->SetOptFit(1111);
  TString outputDir="/home/zoiirene/Output/alignment/";

  //  TString outputDir = outputDir+"";
  TString outputFile = outputDir+"alignf_"+plane+"_run"+run+"_iteration_"+iter;
  c->SaveAs(outputFile+".eps");
  c->SaveAs(outputFile+".pdf");
  c->SaveAs(outputFile+".png");
  c->SaveAs(outputFile+".root");

  
  
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

      //if(PRINT) cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

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

/* oneplane function: for each plane 
- Read data file & assemble pixel hits
- Skip noisy events with more than 400 hits
- Tsunami correction
- Common mode correction (column-wise)
  - in roc coordinates
  - we read along one column
  - a roi is 7 pixels long
  - we find the pixel with the highest and lowest row index in one column
  - we assume there is no significant charge in these pixels
  - we take the average pulse height of these two pixels and subtract it from those in between  dph_i = ph_i - (ph_up + ph_low)/2
- Continue analysis if dph> dphcut
  - for each sample we use the cut yielding the best resolution
  - 25 vs 50
  - Conversion from ADC to ke (updated for Vcal offset correction - as Finn)
- Clustering for events with less than 50 hits (speeding)
  - From a seed pixel add to it all adjacent pixels with a hit
  - Evaluating the cluster isolation
*/
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
	//if(PRINT)  cout << "Analayzing filled data" << endl; 
	string roi;
	getline( Xstream, roi );
	//if(PRINT)  cout << "roi: " << roi << endl;
	//if(PRINT) cout << "roi size: " << roi.size() << endl;
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
		//if(PRINT) cout << " " << s1 << "(" << gap << ")"<<endl;
		int col = atoi( s1.c_str() ); // 4% faster
		start = gap + BLANK.size();
	    
		//getting pixel row
		gap = roi.find( BLANK, start );
		string s2( roi.substr( start, gap - start ) );
		//if(PRINT) cout << " " << s2 << "(" << gap << ")";
		int row = atoi( s2.c_str() );
		start = gap + BLANK.size();

		//getting pixel puls height
		gap = roi.find( BLANK, start );
		string s3( roi.substr( start, gap - start ) );
		//if(PRINT) cout << " " << s1 << "(" << gap << ")";
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
  
  /*
  TH1I * hdx3_clchargeB2e = new TH1I("dx3_clchargeB2e", "triplet dx_clchargeB2e; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clchargeB2e5e = new TH1I("dx3_clchargeB2e5e", "triplet dx_clchargeB2e5e; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clchargeB5e8e = new TH1I("dx3_clchargeB5e8e", "triplet dx_clchargeB5e8e; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clchargeB8e10e = new TH1I("dx3_clchargeB8e10e", "triplet dx_clchargeB8e10e; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clchargeB10e13e = new TH1I("dx3_clchargeB10e13e", "triplet dx_clchargeB10e13e; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clchargeB13e16e = new TH1I("dx3_clchargeB13e16e", "triplet dx_clchargeB13e16e; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clchargeB16e22e = new TH1I("dx3_clchargeB16e22e", "triplet dx_clchargeB16e22e; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clchargeB22e40e = new TH1I("dx3_clchargeB22e40e", "triplet dx_clchargeB22e40e; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clchargeB40e = new TH1I("dx3_clchargeB40e", "triplet dx_clchargeB40e; dx [mm];triplets", 500, -0.5, 0.5 );

  TH1I * hdx3_clphB50adc = new TH1I("dx3_clphB50adc", "triplet dx_clphB50adc; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clphB50adc100adc = new TH1I("dx3_clphB50adc100adc", "triplet dx_clphB50adc100adc; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clphB100adc150adc = new TH1I("dx3_clphB100adc150adc", "triplet dx_clphB100adc150adc; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clphB150adc200adc = new TH1I("dx3_clphB150adc200adc", "triplet dx_clphB150adc200adc; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clphB200adc250adc = new TH1I("dx3_clphB200adc250adc", "triplet dx_clphB200adc250adc; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clphB250adc300adc = new TH1I("dx3_clphB250adc300adc", "triplet dx_clphB250adc300adc; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clphB300adc400adc = new TH1I("dx3_clphB300adc400adc", "triplet dx_clphB300adc400adc; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clphB400adc500adc = new TH1I("dx3_clphB400adc500adc", "triplet dx_clphB400adc500adc; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clphB500adc600adc = new TH1I("dx3_clphB500adc600adc", "triplet dx_clphB500adc600adc; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clphB600adc700adc = new TH1I("dx3_clphB600adc700adc", "triplet dx_clphB600adc700adc; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clphB700adc = new TH1I("dx3_clphB700adc", "triplet dx_clphB700adc; dx [mm];triplets", 500, -0.5, 0.5 );
  */

  /*
  TH1I * hdx3_clchargeAC90evR = new TH1I("dx3_clchargeAC90evR ", "triplet dx_clchargeAC90evR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 90% events in Landau (only high tail)
  TH1I * hdx3_clphAC90evR = new TH1I("dx3_clphAC90evR ", "triplet dx_clphAC90evR ; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clchargeAC80evR = new TH1I("dx3_clchargeAC80evR ", "triplet dx_clchargeAC80evR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 80% events in Landau (only high tail)
  TH1I * hdx3_clphAC80evR = new TH1I("dx3_clphAC80evR ", "triplet dx_clphAC80evR ; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clchargeAC70evR = new TH1I("dx3_clchargeAC70evR ", "triplet dx_clchargeAC70evR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 70% events in Landau (only high tail)
  TH1I * hdx3_clphAC70evR = new TH1I("dx3_clphAC70evR ", "triplet dx_clphAC70evR ; dx [mm];triplets", 500, -0.5, 0.5 );


  TH1I * hdx3_clchargeAC10hR = new TH1I("dx3_clchargeAC10hR ", "triplet dx_clchargeAC10hR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 10% height in Landau (only high tail)
  TH1I * hdx3_clphAC10hR = new TH1I("dx3_clphAC10hR ", "triplet dx_clphAC10hR ; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clchargeAC20hR = new TH1I("dx3_clchargeAC20hR ", "triplet dx_clchargeAC20hR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 20% height in Landau (only high tail)
  TH1I * hdx3_clphAC20hR = new TH1I("dx3_clphAC20hR ", "triplet dx_clphAC20hR ; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clchargeAC30hR = new TH1I("dx3_clchargeAC30hR ", "triplet dx_clchargeAC30hR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 30% height in Landau (only high tail)
  TH1I * hdx3_clphAC30hR = new TH1I("dx3_clphAC30hR ", "triplet dx_clphAC30hR ; dx [mm];triplets", 500, -0.5, 0.5 );
  */
  TH1I * hdx3_clchargeABC90evR = new TH1I("dx3_clchargeABC90evR ", "triplet dx_clchargeABC90evR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 90% events in Landau (only high tail)
  TH1I * hdx3_clphABC90evR = new TH1I("dx3_clphABC90evR ", "triplet dx_clphABC90evR ; dx [mm];triplets", 500, -0.5, 0.5 );
  /*
  TH1I * hdx3_clchargeABC80evR = new TH1I("dx3_clchargeABC80evR ", "triplet dx_clchargeABC80evR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 80% events in Landau (only high tail)
  TH1I * hdx3_clphABC80evR = new TH1I("dx3_clphABC80evR ", "triplet dx_clphABC80evR ; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clchargeABC70evR = new TH1I("dx3_clchargeABC70evR ", "triplet dx_clchargeABC70evR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 70% events in Landau (only high tail)
  TH1I * hdx3_clphABC70evR = new TH1I("dx3_clphABC70evR ", "triplet dx_clphABC70evR ; dx [mm];triplets", 500, -0.5, 0.5 );

  TH1I * hdx3_clchargeABC10hR = new TH1I("dx3_clchargeABC10hR ", "triplet dx_clchargeABC10hR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 10% height in Landau (only high tail)
  TH1I * hdx3_clphABC10hR = new TH1I("dx3_clphABC10hR ", "triplet dx_clphABC10hR ; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clchargeABC20hR = new TH1I("dx3_clchargeABC20hR ", "triplet dx_clchargeABC20hR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 20% height in Landau (only high tail)
  TH1I * hdx3_clphABC20hR = new TH1I("dx3_clphABC20hR ", "triplet dx_clphABC20hR ; dx [mm];triplets", 500, -0.5, 0.5 );
  TH1I * hdx3_clchargeABC30hR = new TH1I("dx3_clchargeABC30hR ", "triplet dx_clchargeABC30hR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 30% height in Landau (only high tail)
  TH1I * hdx3_clphABC30hR = new TH1I("dx3_clphABC30hR ", "triplet dx_clphABC30hR ; dx [mm];triplets", 500, -0.5, 0.5 );
  */





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

  /*
  mapOfHists.insert(std::make_pair("dx3_clchargeB2e",hdx3_clchargeB2e));
  mapOfHists.insert(std::make_pair("dx3_clchargeB2e5e",hdx3_clchargeB2e5e));
  mapOfHists.insert(std::make_pair("dx3_clchargeB5e8e",hdx3_clchargeB5e8e));
  mapOfHists.insert(std::make_pair("dx3_clchargeB8e10e",hdx3_clchargeB8e10e));
  mapOfHists.insert(std::make_pair("dx3_clchargeB10e13e",hdx3_clchargeB10e13e));
  mapOfHists.insert(std::make_pair("dx3_clchargeB13e16e",hdx3_clchargeB13e16e));
  mapOfHists.insert(std::make_pair("dx3_clchargeB16e22e",hdx3_clchargeB16e22e));
  mapOfHists.insert(std::make_pair("dx3_clchargeB22e40e",hdx3_clchargeB22e40e));
  mapOfHists.insert(std::make_pair("dx3_clchargeB40e",hdx3_clchargeB40e));
  mapOfHists.insert(std::make_pair("dx3_clphB50adc",hdx3_clphB50adc));
  mapOfHists.insert(std::make_pair("dx3_clphB50adc100adc",hdx3_clphB50adc100adc));
  mapOfHists.insert(std::make_pair("dx3_clphB100adc150adc",hdx3_clphB100adc150adc));
  mapOfHists.insert(std::make_pair("dx3_clphB150adc200adc",hdx3_clphB150adc200adc));
  mapOfHists.insert(std::make_pair("dx3_clphB200adc250adc",hdx3_clphB200adc250adc));
  mapOfHists.insert(std::make_pair("dx3_clphB250adc300adc",hdx3_clphB250adc300adc));
  mapOfHists.insert(std::make_pair("dx3_clphB300adc400adc",hdx3_clphB300adc400adc));
  mapOfHists.insert(std::make_pair("dx3_clphB400adc500adc",hdx3_clphB400adc500adc));
  mapOfHists.insert(std::make_pair("dx3_clphB500adc600adc",hdx3_clphB500adc600adc));
  mapOfHists.insert(std::make_pair("dx3_clphB600adc700adc",hdx3_clphB600adc700adc));
  mapOfHists.insert(std::make_pair("dx3_clphB700adc",hdx3_clphB700adc));
  */

  /*
  mapOfHists.insert(std::make_pair("dx3_clchargeAC90evR",hdx3_clchargeAC90evR));
  mapOfHists.insert(std::make_pair("dx3_clphAC90evR",hdx3_clphAC90evR));
  mapOfHists.insert(std::make_pair("dx3_clchargeAC80evR",hdx3_clchargeAC80evR));
  mapOfHists.insert(std::make_pair("dx3_clphAC80evR",hdx3_clphAC80evR));
  mapOfHists.insert(std::make_pair("dx3_clchargeAC70evR",hdx3_clchargeAC70evR));
  mapOfHists.insert(std::make_pair("dx3_clphAC70evR",hdx3_clphAC70evR));


  mapOfHists.insert(std::make_pair("dx3_clchargeAC10hR",hdx3_clchargeAC10hR));
  mapOfHists.insert(std::make_pair("dx3_clphAC10hR",hdx3_clphAC10hR));
  mapOfHists.insert(std::make_pair("dx3_clchargeAC20hR",hdx3_clchargeAC20hR));
  mapOfHists.insert(std::make_pair("dx3_clphAC20hR",hdx3_clphAC20hR));
  mapOfHists.insert(std::make_pair("dx3_clchargeAC30hR",hdx3_clchargeAC30hR));
  mapOfHists.insert(std::make_pair("dx3_clphAC30hR",hdx3_clphAC30hR));
  
  
  mapOfHists.insert(std::make_pair("dx3_clchargeABC90evR",hdx3_clchargeABC90evR));
  mapOfHists.insert(std::make_pair("dx3_clphABC90evR",hdx3_clphABC90evR));
  
  mapOfHists.insert(std::make_pair("dx3_clchargeABC80evR",hdx3_clchargeABC80evR));
  mapOfHists.insert(std::make_pair("dx3_clphABC80evR",hdx3_clphABC80evR));
  mapOfHists.insert(std::make_pair("dx3_clchargeABC70evR",hdx3_clchargeABC70evR));
  mapOfHists.insert(std::make_pair("dx3_clphABC70evR",hdx3_clphABC70evR));

  mapOfHists.insert(std::make_pair("dx3_clchargeABC10hR",hdx3_clchargeABC10hR));
  mapOfHists.insert(std::make_pair("dx3_clphABC10hR",hdx3_clphABC10hR));
  mapOfHists.insert(std::make_pair("dx3_clchargeABC20hR",hdx3_clchargeABC20hR));
  mapOfHists.insert(std::make_pair("dx3_clphABC20hR",hdx3_clphABC20hR));
  mapOfHists.insert(std::make_pair("dx3_clchargeABC30hR",hdx3_clchargeABC30hR));
  mapOfHists.insert(std::make_pair("dx3_clphABC30hR",hdx3_clphABC30hR));

  */


  
  return mapOfHists;
}


void  fillControlHists(histoMap mapOfHists, TString selection, double dx3, double dy3, vector<cluster>::iterator clusterA, vector<cluster>::iterator clusterB, vector<cluster>::iterator clusterC, int nrowB, int ncolB,double xmod, unsigned iev,double xB, double yB,double xAr, double yAr, double xCr, double yCr,double dxCA,double etaA, double etaB, double etaC,TFile * histofile, TString fileName, TH1I *hclphA,TH1I *hclphB,TH1I *hclphC, TH1I *hclqA, TH1I *hclqB, TH1I *hclqC)
{

  if(PRINT) std::cout << "********** filling hists of "<< selection << endl;
  
  TDirectory *cdtof = (TDirectory *)histofile->Get(selection);
  if(PRINT)  cout << "looking for directory " << fileName << ":"<< selection << endl;

  cdtof->cd(fileName+":"+selection);

  if(PRINT)  cout << " directory " << fileName << ":"<< selection << " found "<< endl;
  if(PRINT)cout << " cl B " << clusterB->size << endl;
  if(PRINT) cout << "filling in map" << endl;

  if(PRINT)cout << " hclph B " << hclphB->GetEntries() << endl;
  if(PRINT)cout << " hclq B " << hclqB->GetEntries() << endl;

  TH1I * hclq[DreiMasterPlanes];
  hclq[0]=hclqA;
  hclq[1]=hclqB;
  hclq[2]=hclqC;

  TH1I * hclph[DreiMasterPlanes];
  hclph[0]=hclphA;
  hclph[1]=hclphB;
  hclph[2]=hclphC;


  
  double dyCA = yCr - yAr;
  double dxyCA = sqrt( dxCA*dxCA + dyCA*dyCA );

  auto  search =  mapOfHists.find("dxyCA");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(dxyCA);
  } else {
    if(PRINT) std::cout << "Not found "<< "dxyCA" << endl;
  }

  
  double dxy = sqrt( dx3*dx3 + dy3*dy3 );

  search =  mapOfHists.find("dxy");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(dxy);
  } else {
    if(PRINT) std::cout << "Not found "<< "dxy" << endl;
  }


  
  search =  mapOfHists.find("xA");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(xAr);
  } else {
    if(PRINT) std::cout << "Not found "<< "xA" << endl;
  }

  search =  mapOfHists.find("yA");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(yAr);
  } else {
    if(PRINT) std::cout << "Not found "<< "yAr" << endl;
  }

  search =  mapOfHists.find("xB");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(xB);
  } else {
    if(PRINT) std::cout << "Not found "<< "xB" << endl;
  }

  search =  mapOfHists.find("yB");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(yB);
  } else {
    if(PRINT) std::cout << "Not found "<< "yB" << endl;
  }

  search =  mapOfHists.find("xC");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(xCr);
  } else {
    if(PRINT) std::cout << "Not found "<< "xC" << endl;
  }

  search =  mapOfHists.find("yC");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(yCr);
  } else {
    if(PRINT) std::cout << "Not found "<< "yCr" << endl;
  }


  
  search =  mapOfHists.find("dx3");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3" << endl;
  }


  search =  mapOfHists.find("dy3");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(dy3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dy3" << endl;
  }

  search =  mapOfHists.find("clsizeA");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterA->size);
  } else {
    if(PRINT) std::cout << "Not found "<< "clsizeA" << endl;
  }

  search =  mapOfHists.find("clsizeB");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterB->size);
  } else {
    if(PRINT) std::cout << "Not found "<< "clsizeB" << endl;
  }

    search =  mapOfHists.find("clsizeC");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterC->size);
  } else {
    if(PRINT) std::cout << "Not found "<< "clsizeC" << endl;
  }


  search =  mapOfHists.find("nrowB");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(nrowB);
  } else {
    if(PRINT) std::cout << "Not found "<< "nrowB" << endl;
  }

  search =  mapOfHists.find("ncolB");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(ncolB);
  } else {
    if(PRINT) std::cout << "Not found "<< "ncolB" << endl;
  }

  search =  mapOfHists.find("clmapA");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterA->col, clusterA->row);
  } else {
    if(PRINT) std::cout << "Not found "<< "clmapA" << endl;
  }

  search =  mapOfHists.find("clmapB");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterB->col, clusterB->row);
  } else {
    if(PRINT) std::cout << "Not found "<< "clmapB" << endl;
  }

  search =  mapOfHists.find("clmapC");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterC->col, clusterC->row);
  } else {
    if(PRINT) std::cout << "Not found "<< "clmapC" << endl;
  }


  search =  mapOfHists.find("clchargeA");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterA->q);
  } else {
    if(PRINT) std::cout << "Not found "<< "clchargeA" << endl;
  }

    search =  mapOfHists.find("clchargeB");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterB->q);
  } else {
    if(PRINT) std::cout << "Not found "<< "clchargeB" << endl;
  }

    search =  mapOfHists.find("clchargeC");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterC->q);
  } else {
    if(PRINT) std::cout << "Not found "<< "clchargeC" << endl;
  }

  search =  mapOfHists.find("clphA");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterA->sum);
  } else {
    if(PRINT) std::cout << "Not found "<< "clphA" << endl;
  }

    search =  mapOfHists.find("clphB");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterB->sum);
  } else {
    if(PRINT) std::cout << "Not found "<< "clphB" << endl;
  }

  search =  mapOfHists.find("clphC");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterC->sum);
  } else {
    if(PRINT) std::cout << "Not found "<< "clphC" << endl;
  }

  search =  mapOfHists.find("pxphA");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    for( int ipx = 0; ipx < clusterA->size; ++ipx )     search->second->Fill( clusterA->vpix[ipx].ph );
  } else {
    if(PRINT) std::cout << "Not found "<< "pxphA" << endl;
  }

  search =  mapOfHists.find("pxchargeA");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    for( int ipx = 0; ipx < clusterA->size; ++ipx )     search->second->Fill( clusterA->vpix[ipx].q );
  } else {
    if(PRINT) std::cout << "Not found "<< "pxchargeA" << endl;
  }

  search =  mapOfHists.find("pxphB");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    for( int ipx = 0; ipx < clusterB->size; ++ipx )     search->second->Fill( clusterB->vpix[ipx].ph );
  } else {
    if(PRINT) std::cout << "Not found "<< "pxphB" << endl;
  }

  search =  mapOfHists.find("pxchargeB");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    for( int ipx = 0; ipx < clusterB->size; ++ipx )     search->second->Fill( clusterB->vpix[ipx].q );
  } else {
    if(PRINT) std::cout << "Not found "<< "pxchargeB" << endl;
  }

  search =  mapOfHists.find("pxphC");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    for( int ipx = 0; ipx < clusterC->size; ++ipx )     search->second->Fill( clusterC->vpix[ipx].ph );
  } else {
    if(PRINT) std::cout << "Not found "<< "pxphC" << endl;
  }

  search =  mapOfHists.find("pxchargeC");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    for( int ipx = 0; ipx < clusterC->size; ++ipx )     search->second->Fill( clusterC->vpix[ipx].q );
  } else {
    if(PRINT) std::cout << "Not found "<< "pxchargeC" << endl;
  }


  search =  mapOfHists.find("etaA");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterA->size ==2)
      search->second->Fill(etaA);
  } else {
    if(PRINT) std::cout << "Not found "<< "etaA" << endl;
  }

  search =  mapOfHists.find("etaB");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size ==2)
      search->second->Fill(etaB);
  } else {
    if(PRINT) std::cout << "Not found "<< "etaB" << endl;
  }

  search =  mapOfHists.find("etaC");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterC->size ==2)
      search->second->Fill(etaC);
  } else {
    if(PRINT) std::cout << "Not found "<< "etaA" << endl;
  }
		    
  search =  mapOfHists.find("pxchargeB1st");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size ==2)
      search->second->Fill( clusterB->vpix[0].q );
  } else {
    if(PRINT) std::cout << "Not found "<< "pxchargeB1st" << endl;
  }

  search =  mapOfHists.find("pxchargeB2nd");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size ==2)
      search->second->Fill( clusterB->vpix[1].q );
  } else {
    if(PRINT) std::cout << "Not found "<< "pxchargeB2nd" << endl;
  }




  search =  mapOfHists.find("nrowvsxmB3");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

  if( fifty )
    search->second->Fill( xmod*1E3, ncolB );
  else
    search->second->Fill( xmod*1E3, nrowB );

  } else {
    if(PRINT) std::cout << "Not found "<< "nrowvsxmB3" << endl;
  }

  search =  mapOfHists.find("clqvsxmB3");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill( xmod*1E3, clusterB->q );
  } else {
    if(PRINT) std::cout << "Not found "<< "clqvsxmB3" << endl;
  }

  search =  mapOfHists.find("dx3vsev");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill( iev, dx3 );
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3vsev" << endl;
  }

  search =  mapOfHists.find("dx3vsx");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill( xB, dx3 ); // turn 
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3vsx" << endl;
  }

  search =  mapOfHists.find("dx3vsy");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill( yB, dx3 ); // rot
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3vsy" << endl;
  }

  search =  mapOfHists.find("dx3vsxm");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill( xmod*1E3, dx3 ); // rot
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3vsxm" << endl;
  }
  search =  mapOfHists.find("madx3vsdx");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill( dxCA*1E3, fabs(dx3) ); // dxCA
  } else {
    if(PRINT) std::cout << "Not found "<< "madx3vsdx" << endl;
  }

  search =  mapOfHists.find("madx3vsx");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill( xB, fabs(dx3) );
  } else {
    if(PRINT) std::cout << "Not found "<< "madx3vsx" << endl;
  }

  search =  mapOfHists.find("madx3vsy");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill( yB, fabs(dx3) );
      } else {
    if(PRINT) std::cout << "Not found "<< "madx3vsy" << endl;
  }

  search =  mapOfHists.find("madx3vsxm");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(xmod*1E3, fabs(dx3) );
      } else {
    if(PRINT) std::cout << "Not found "<< "madx3vsxm" << endl;
  }

  search =  mapOfHists.find("madx3vsq");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterB->q, fabs(dx3) );
      } else {
    if(PRINT) std::cout << "Not found "<< "madx3vsq" << endl;
  }

  search =  mapOfHists.find("madx3vsn");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    search->second->Fill(clusterB->size, fabs(dx3) );
      } else {
    if(PRINT) std::cout << "Not found "<< "madx3vsn" << endl;
  }

  search =  mapOfHists.find("etavsxmB3");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end() ) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size == 2)
      search->second->Fill( xmod*1E3, etaB ); // sine 
  } else {
    if(PRINT) std::cout << "Not found "<< "etavsxmB3" << endl;
  }
  
  search =  mapOfHists.find("madx3vseta");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end() ) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size == 2)
      search->second->Fill( etaB, fabs(dx3) ); // flat 
  } else {
    if(PRINT) std::cout << "Not found "<< "madx3vseta" << endl;
  }
  

  search =  mapOfHists.find("dx3_clsizeB1");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size ==1)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clsizeB1" << endl;
  }

  search =  mapOfHists.find("dx3_clsizeB2");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size ==2)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clsizeB2" << endl;
  }

  search =  mapOfHists.find("dx3_clsizeB3");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size ==3)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clsizeB3" << endl;
  }

  search =  mapOfHists.find("dx3_clsizeB4");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size ==4)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clsizeB4" << endl;
  }
  search =  mapOfHists.find("dx3_clsizeB5");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size ==5)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clsizeB5" << endl;
  }
  search =  mapOfHists.find("dx3_clsizeB6");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size ==6)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clsizeB6" << endl;
  }
  search =  mapOfHists.find("dx3_clsizeB7m");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->size >6)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clsizeB7m" << endl;
  }




  

  /*
  search =  mapOfHists.find("dx3_clchargeB2e");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->q <= 2)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeB2e" << endl;
  }

  search =  mapOfHists.find("dx3_clchargeB2e5e");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->q > 2 && clusterB->q <=5)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeB2e5e" << endl;
  }

  search =  mapOfHists.find("dx3_clchargeB5e8e");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->q > 5 && clusterB->q <=8)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeB5e8e" << endl;
  }

  search =  mapOfHists.find("dx3_clchargeB8e10e");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->q > 8 && clusterB->q <=10)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeB8e10e" << endl;
  }

  search =  mapOfHists.find("dx3_clchargeB10e13e");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->q > 10 && clusterB->q <=13)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeB10e13e" << endl;
  }

  search =  mapOfHists.find("dx3_clchargeB13e16e");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->q > 13 && clusterB->q <=16)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeB13e16e" << endl;
  }

  search =  mapOfHists.find("dx3_clchargeB16e22e");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->q > 16 && clusterB->q <=22)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeB16e22e" << endl;
  }

  search =  mapOfHists.find("dx3_clchargeB22e40e");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->q > 22 && clusterB->q <=40)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeB22e40e" << endl;
  }

  search =  mapOfHists.find("dx3_clchargeB40e");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->q > 40 )
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeB40e" << endl;
  }

  search =  mapOfHists.find("dx3_clphB50adc");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->sum <= 50)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphB50adc" << endl;
  }

  search =  mapOfHists.find("dx3_clphB50adc100adc");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->sum > 50 && clusterB->sum <=100)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphB50adc100adc" << endl;
  }

  search =  mapOfHists.find("dx3_clphB100adc150adc");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->sum > 100 && clusterB->sum <=150)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphB100adc150adc" << endl;
  }

  search =  mapOfHists.find("dx3_clphB150adc200adc");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->sum > 150 && clusterB->sum <=200)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphB150adc200adc" << endl;
  }

  search =  mapOfHists.find("dx3_clphB200adc250adc");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->sum > 200 && clusterB->sum <=250)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphB200adc250adc" << endl;
  }

  search =  mapOfHists.find("dx3_clphB250adc300adc");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->sum > 250 && clusterB->sum <=300)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphB250adc300adc" << endl;
  }

  search =  mapOfHists.find("dx3_clphB300adc400adc");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->sum > 300 && clusterB->sum <=400)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphB300adc400adc" << endl;
  }

  search =  mapOfHists.find("dx3_clphB400adc500adc");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->sum > 400 && clusterB->sum <=500)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphB400adc500adc" << endl;
  }

  search =  mapOfHists.find("dx3_clphB500adc600adc");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->sum > 500 && clusterB->sum <=600)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphB500adc600adc" << endl;
  }

  search =  mapOfHists.find("dx3_clphB600adc700adc");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->sum > 600 && clusterB->sum <=700)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphB600adc700adc" << endl;
  }

  search =  mapOfHists.find("dx3_clphB700adc");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';
    if(clusterB->sum > 700)
      search->second->Fill(dx3);
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphB700adc" << endl;
  }
  
  
  double integral[DreiMasterPlanes];
  double integralPH[DreiMasterPlanes];

  double integral90[DreiMasterPlanes];
  double integral80[DreiMasterPlanes];
  double integral70[DreiMasterPlanes];
  double integralPH90[DreiMasterPlanes];
  double integralPH80[DreiMasterPlanes];
  double integralPH70[DreiMasterPlanes];

  int high90[DreiMasterPlanes];
  int high80[DreiMasterPlanes];
  int high70[DreiMasterPlanes];
  double highPH90[DreiMasterPlanes];
  double highPH80[DreiMasterPlanes];
  double highPH70[DreiMasterPlanes];

  for(int j =0; j < DreiMasterPlanes; j++)
    {
      if(DEBUG)   cout << " plane " << j << endl;
      integral[j] = hclq[j]->Integral();
      if(DEBUG)   cout << " integral " << integral[j] << endl;

      integralPH[j] = hclph[j]->Integral();
      if(DEBUG)   cout << " integralPH " << integralPH[j] << endl;

      integral90[j] = 0.9*hclq[j]->Integral(); //careful!! you should take into account the peak at low value! Hist is now filled only for isolated clusters and the effect is reduced
      if(DEBUG)   cout << " integral90 " << integral90[j] << endl;
      integral80[j] = 0.8*hclq[j]->Integral(); //careful!! you should take into account the peak at low value! Hist is now filled only for isolated clusters and the effect is reduced
      integral70[j] = 0.7*hclq[j]->Integral(); //careful!! you should take into account the peak at low value! Hist is now filled only for isolated clusters and the effect is reduced

      integralPH90[j] = 0.9*hclph[j]->Integral();
      integralPH80[j] = 0.8*hclph[j]->Integral();
      integralPH70[j] = 0.7*hclph[j]->Integral();

      high90[j] = 0;
      high80[j] = 0;
      high70[j] = 0;
      highPH90[j] = 0;
      highPH80[j] = 0;
      highPH70[j] = 0;


      int i = 0;
      while(integral[j]>integral90[j])
	{
	  //      cout << " while "<< i << endl;
	  integral[j] = hclq[j]->Integral(1,hclq[j]->GetNbinsX()-i);
	  high90[j] = hclq[j]->GetNbinsX()-i;
	  //cout << " integral " << integral << " low " << low << " high " << high << endl;
	  i++;
	}
      if(DEBUG)    cout << "final integral 90 " << integral[j] << " high " << high90[j] << endl;
      i = 0;
      while(integral[j]>integral80[j])
	{
	  //      cout << " while "<< i << endl;
	  integral[j] = hclq[j]->Integral(1,hclq[j]->GetNbinsX()-i);
	  high80[j] = hclq[j]->GetNbinsX()-i;
	  //cout << " integral " << integral << " low " << low << " high " << high << endl;
	  i++;
	}
      if(DEBUG)    cout << "final integral 80 " << integral[j] << " high " << high80[j] << endl;
      i = 0;
      while(integral[j]>integral70[j])
	{
	  //      cout << " while "<< i << endl;
	  integral[j] = hclq[j]->Integral(1,hclq[j]->GetNbinsX()-i);
	  high70[j] = hclq[j]->GetNbinsX()-i;
	  //cout << " integral " << integral << " low " << low << " high " << high << endl;
	  i++;
	}
      if(DEBUG)    cout << "final integral 70 " << integral[j] << " high " << high70[j] << endl;

      i=0;
      while(integralPH[j]>integralPH90[j])
	{
	  //      cout << " while "<< i << endl;
	  integralPH[j] = hclph[j]->Integral(1,hclph[j]->GetNbinsX()-i);
	  highPH90[j] = hclph[j]->GetNbinsX()-i;
	  //cout << " integral " << integral << " low " << low << " high " << high << endl;
	  i++;
	}
      if(PRINT)  cout << "final integralPH 90 " << integralPH[j]  << " high " << highPH90[j] << endl;

      i=0;
      while(integralPH[j]>integralPH80[j])
	{
	  //      cout << " while "<< i << endl;
	  integralPH[j] = hclph[j]->Integral(1,hclph[j]->GetNbinsX()-i);
	  highPH80[j] = hclph[j]->GetNbinsX()-i;
	  //cout << " integral " << integral << " low " << low << " high " << high << endl;
	  i++;
	}
      if(DEBUG)  cout << "final integralPH 80 " << integralPH[j]  << " high " << highPH80[j] << endl;

      i=0;
      while(integralPH[j]>integralPH70[j])
	{
	  //      cout << " while "<< i << endl;
	  integralPH[j] = hclph[j]->Integral(1,hclph[j]->GetNbinsX()-i);
	  highPH70[j] = hclph[j]->GetNbinsX()-i;
	  //cout << " integral " << integral << " low " << low << " high " << high << endl;
	  i++;
	}
      if(DEBUG)  cout << "final integralPH 70 " << integralPH[j]  << " high " << highPH70[j] << endl;

    }

  /*
  search =  mapOfHists.find("dx3_clchargeAC90evR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    if(clusterA->q <= hclq[0]->GetBinCenter(high90[0]) && clusterC->q <= hclq[2]->GetBinCenter(high90[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeAC90evR" << endl;
  }  

  
  search =  mapOfHists.find("dx3_clchargeABC90evR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    if(clusterB->q <= hclq[1]->GetBinCenter(high90[1]) && clusterA->q <= hclq[0]->GetBinCenter(high90[0]) && clusterC->q <= hclq[2]->GetBinCenter(high90[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeABC90evR" << endl;
  }


  
  search =  mapOfHists.find("dx3_clchargeAC80evR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    if(clusterA->q <= hclq[0]->GetBinCenter(high80[0]) && clusterC->q <= hclq[2]->GetBinCenter(high80[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeAC80evR" << endl;
  }
  search =  mapOfHists.find("dx3_clchargeABC80evR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    if(clusterB->q <= hclq[1]->GetBinCenter(high80[1]) && clusterA->q <= hclq[0]->GetBinCenter(high80[0]) && clusterC->q <= hclq[2]->GetBinCenter(high80[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeABC80evR" << endl;
  }

  search =  mapOfHists.find("dx3_clchargeAC80evR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    if(clusterA->q <= hclq[0]->GetBinCenter(high80[0]) && clusterC->q <= hclq[2]->GetBinCenter(high80[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeAC80evR" << endl;
  }
  search =  mapOfHists.find("dx3_clchargeABC70evR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    if(clusterB->q <= hclq[1]->GetBinCenter(high70[1]) && clusterA->q <= hclq[0]->GetBinCenter(high70[0]) && clusterC->q <= hclq[2]->GetBinCenter(high70[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeABC70evR" << endl;
  }



  search =  mapOfHists.find("dx3_clphAC90evR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    if( clusterA->sum <= hclph[0]->GetBinCenter(highPH90[0]) && clusterC->sum <= hclph[2]->GetBinCenter(highPH90[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphAC90evR" << endl;
  }

  
  search =  mapOfHists.find("dx3_clphABC90evR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    if(clusterB->sum <= hclph[1]->GetBinCenter(highPH90[1]) && clusterA->sum <= hclph[0]->GetBinCenter(highPH90[0]) && clusterC->sum <= hclph[2]->GetBinCenter(highPH90[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphABC90evR" << endl;
  }


  
  search =  mapOfHists.find("dx3_clphAC80evR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    if( clusterA->sum <= hclph[0]->GetBinCenter(highPH80[0]) && clusterC->sum <= hclph[2]->GetBinCenter(highPH80[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphAC80evR" << endl;
  }

  search =  mapOfHists.find("dx3_clphABC80evR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    if(clusterB->sum <= hclph[1]->GetBinCenter(highPH80[1]) && clusterA->sum <= hclph[0]->GetBinCenter(highPH80[0]) && clusterC->sum <= hclph[2]->GetBinCenter(highPH80[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphABC80evR" << endl;
  }

  search =  mapOfHists.find("dx3_clphAC70evR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    if( clusterA->sum <= hclph[0]->GetBinCenter(highPH70[0]) && clusterC->sum <= hclph[2]->GetBinCenter(highPH70[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphAC70evR" << endl;
  }

  search =  mapOfHists.find("dx3_clphABC70evR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    if(clusterB->sum <= hclph[1]->GetBinCenter(highPH70[1]) && clusterA->sum <= hclph[0]->GetBinCenter(highPH70[0]) && clusterC->sum <= hclph[2]->GetBinCenter(highPH70[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphABC70evR" << endl;
  }


  //use a fit!!! - probably not needed after iso requirements
  double maximum[DreiMasterPlanes];
  int maximumBin[DreiMasterPlanes];
  double maximum10[DreiMasterPlanes];
  int bin10[DreiMasterPlanes];
  double maximum20[DreiMasterPlanes];
  int bin20[DreiMasterPlanes];
  double maximum30[DreiMasterPlanes];
  int bin30[DreiMasterPlanes];

  double maximumPH[DreiMasterPlanes];
  int maximumBinPH[DreiMasterPlanes];
  double maximum10PH[DreiMasterPlanes];
  int bin10PH[DreiMasterPlanes];
  double maximum20PH[DreiMasterPlanes];
  int bin20PH[DreiMasterPlanes];
  double maximum30PH[DreiMasterPlanes];
  int bin30PH[DreiMasterPlanes];



  
  for(int j =0; j < DreiMasterPlanes; j++)
    {
      if(PRINT)       cout << " plane " << j << endl;
      maximum[j] = hclq[j]->GetMaximum();
      if(PRINT)       cout << " maximum " << maximum[j] << endl;
      maximumBin[j] = hclq[j]->GetMaximumBin();
      if(PRINT) 	cout << " maximumBin " << maximumBin[j] << endl;

      maximumPH[j] = hclph[j]->GetMaximum();
      if(PRINT)       cout << " maximum " << maximumPH[j] << endl;
      maximumBinPH[j] = hclph[j]->GetMaximumBin();
      if(PRINT) 	cout << " maximumBin " << maximumBinPH[j] << endl;

      
      maximum10[j] = 0.1*hclq[j]->GetMaximum();
      maximum20[j] = 0.2*hclq[j]->GetMaximum();
      maximum30[j] = 0.3*hclq[j]->GetMaximum();
      maximum10PH[j] = 0.1*hclph[j]->GetMaximum();
      maximum20PH[j] = 0.2*hclph[j]->GetMaximum();
      maximum30PH[j] = 0.3*hclph[j]->GetMaximum();

   
      hclq[j]->GetBinWithContent(maximum10[j],bin10[j],maximumBin[j],hclq[j]->GetNbinsX(),11);
      if(PRINT) cout << " bin on 10 % " << bin10[j] << endl;
   
      hclq[j]->GetBinWithContent(maximum20[j],bin20[j],maximumBin[j],hclq[j]->GetNbinsX(),11);
      if(PRINT) cout << " bin on 20 % " << bin20[j] << endl;
   
      hclq[j]->GetBinWithContent(maximum30[j],bin30[j],maximumBin[j],hclq[j]->GetNbinsX(),11);
      if(PRINT) cout << " bin on 30 % " << bin30[j] << endl;

   
      hclph[j]->GetBinWithContent(maximum10PH[j],bin10PH[j],maximumBinPH[j],hclph[j]->GetNbinsX(),11);
      if(PRINT) cout << " bin on 10 % " << bin10PH[j] << endl;
   
      hclph[j]->GetBinWithContent(maximum20PH[j],bin20PH[j],maximumBinPH[j],hclph[j]->GetNbinsX(),11);
      if(PRINT) cout << " bin on 20 % " << bin20PH[j] << endl;

     
      hclph[j]->GetBinWithContent(maximum30PH[j],bin30PH[j],maximumBinPH[j],hclph[j]->GetNbinsX(),11);
      if(PRINT) cout << " bin on 30 % " << bin30PH[j] << endl;
      
    }

  search =  mapOfHists.find("dx3_clchargeAC10hR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    
    if( clusterA->q <= hclq[0]->GetBinCenter(bin10[0]) && clusterC->q <= hclq[2]->GetBinCenter(bin10[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeAC10hR" << endl;
  }

  search =  mapOfHists.find("dx3_clchargeABC10hR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    
    if(clusterB->q <= hclq[1]->GetBinCenter(bin10[1]) && clusterA->q <= hclq[0]->GetBinCenter(bin10[0]) && clusterC->q <= hclq[2]->GetBinCenter(bin10[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeABC10hR" << endl;
  }

  search =  mapOfHists.find("dx3_clchargeAC20hR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    
    if( clusterA->q <= hclq[0]->GetBinCenter(bin20[0]) && clusterC->q <= hclq[2]->GetBinCenter(bin20[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeAC20hR" << endl;
  }
  search =  mapOfHists.find("dx3_clchargeABC20hR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    
    if(clusterB->q <= hclq[1]->GetBinCenter(bin20[1]) && clusterA->q <= hclq[0]->GetBinCenter(bin20[0]) && clusterC->q <= hclq[2]->GetBinCenter(bin20[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeABC20hR" << endl;
  }

  search =  mapOfHists.find("dx3_clchargeAC30hR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    
    if( clusterA->q <= hclq[0]->GetBinCenter(bin30[0]) && clusterC->q <= hclq[2]->GetBinCenter(bin30[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeAC30hR" << endl;
  }
  search =  mapOfHists.find("dx3_clchargeABC30hR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    
    if(clusterB->q <= hclq[1]->GetBinCenter(bin30[1]) && clusterA->q <= hclq[0]->GetBinCenter(bin30[0]) && clusterC->q <= hclq[2]->GetBinCenter(bin30[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clchargeABC30hR" << endl;
  }



  search =  mapOfHists.find("dx3_clphABC10hR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    if(clusterB->sum <= hclph[1]->GetBinCenter(bin10[1]) && clusterA->sum <= hclph[0]->GetBinCenter(bin10[0]) && clusterC->sum <= hclph[2]->GetBinCenter(bin10[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphABC10hR" << endl;
  }

  search =  mapOfHists.find("dx3_clphAC10hR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    if(clusterA->sum <= hclph[0]->GetBinCenter(bin10[0]) && clusterC->sum <= hclph[2]->GetBinCenter(bin10[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphAC10hR" << endl;
  }


  search =  mapOfHists.find("dx3_clphABC20hR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    if(clusterB->sum <= hclph[1]->GetBinCenter(bin20[1]) && clusterA->sum <= hclph[0]->GetBinCenter(bin20[0]) && clusterC->sum <= hclph[2]->GetBinCenter(bin20[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphABC20hR" << endl;
  }

  search =  mapOfHists.find("dx3_clphAC20hR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    if(clusterA->sum <= hclph[0]->GetBinCenter(bin20[0]) && clusterC->sum <= hclph[2]->GetBinCenter(bin20[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphAC20hR" << endl;
  }

  search =  mapOfHists.find("dx3_clphABC30hR");
  if(DEBUG) cout << "search" << endl;
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    if(clusterB->sum <= hclph[1]->GetBinCenter(bin30[1]) && clusterA->sum <= hclph[0]->GetBinCenter(bin30[0]) && clusterC->sum <= hclph[2]->GetBinCenter(bin30[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphABC30hR" << endl;
  }

  search =  mapOfHists.find("dx3_clphAC30hR");
  if(DEBUG) cout << "search" << endl; 
  if (search !=  mapOfHists.end()) {
    if(DEBUG) std::cout << "Found " << search->first  << '\n';

    if(clusterA->sum <= hclph[0]->GetBinCenter(bin30[0]) && clusterC->sum <= hclph[2]->GetBinCenter(bin30[2]))
       search->second->Fill(dx3);
       
  } else {
    if(PRINT) std::cout << "Not found "<< "dx3_clphAC30hR" << endl;
  }
  */

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


  //hdx3b = new TH1I("hdx3", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
  // hdx3tree = new TH1I("hdx3tree", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
 hdx3_clchargeABC90evR = new TH1I("dx3_clchargeABC90evR ", "triplet dx_clchargeABC90evR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 90% events in Landau (only high tail)
 hdx3_clphABC90evR = new TH1I("dx3_clphABC90evR ", "triplet dx_clphABC90evR ; dx [mm];triplets", 500, -0.5, 0.5 );
 hdxyCA_clphABC90evR  = new TH1I("dxyCA_clphABC90evR","dxyCA_clphABC90evR ; dxyCA [mm]; triplets",100, 0.,0.3);
 hdxy_clphABC90evR = new TH1I("dxy_clphABC90evR","dxy_clphABC90evR ; dxy [mm]; triplets",100, 0.,0.3);
 
 hdx3_clchargeABC90evR95 = new TH1I("dx3_clchargeABC90evR95 ", "triplet dx_clchargeABC90evR95 ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 90% events in Landau (only high tail)
 hdx3_clphABC90evR95 = new TH1I("dx3_clphABC90evR95 ", "triplet dx_clphABC90evR95 ; dx [mm];triplets", 500, -0.5, 0.5 );
 hdx3_clchargeABC90evR99 = new TH1I("dx3_clchargeABC90evR99 ", "triplet dx_clchargeABC90evR99 ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 90% events in Landau (only high tail)
 hdx3_clphABC90evR99 = new TH1I("dx3_clphABC90evR99 ", "triplet dx_clphABC90evR99 ; dx [mm];triplets", 500, -0.5, 0.5 );
 
 
 

 hnrowB_ph95 = new TH1I("hnrowB_ph95", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
 hnrowB_q95 = new TH1I("hnrowB_q95", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
 hclB_ph95 = new TH1I("hclB_ph95 ", "triplet dx_clphABC90evR RMS 95; dx [mm];triplets", 200, 0, 1000 );
 hclB_q95 = new TH1I("hclB_q95 ", "triplet dx_clchargeABC90evR RMS 95; dx [mm];triplets", 100, 0, 50 );;
 hdxyCA_ph95  = new TH1I("dxyCA_ph95","dxyCA_ph95 ; dxyCA [mm]; triplets",100, 0.,0.3);
 hdxy_ph95 = new TH1I("dxy_ph95","dxy_ph95 ; dxy [mm]; triplets",100, 0.,0.3);
 hnrowB_ph5 = new TH1I("hnrowB_ph5", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
 hnrowB_q5 = new TH1I("hnrowB_q5", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
 hclB_ph5 = new TH1I("hclB_ph5 ", "triplet dx_clphABC90evR RMS out 5; dx [mm];triplets", 200, 0, 1000 );
 hclB_q5 = new TH1I("hclB_q5 ", "triplet dx_clchargeABC90evR RMS out 5; dx [mm];triplets", 100, 0, 50 );;
 
 hnrowB_ph99 = new TH1I("hnrowB_ph99", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
 hnrowB_q99 = new TH1I("hnrowB_q99", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
 hclB_ph99 = new TH1I("hclB_ph99 ", "triplet dx_clphABC90evR RMS 99; dx [mm];triplets", 200, 0, 1000 );
 hclB_q99 = new TH1I("hclB_q99 ", "triplet dx_clchargeABC90evR RMS 99; dx [mm];triplets", 100, 0, 50 );;
 hdxyCA_ph99  = new TH1I("dxyCA_ph99","dxyCA_ph99 ; dxyCA [mm]; triplets",100, 0.,0.3);
 hdxy_ph99 = new TH1I("dxy_ph99","dxy_ph99 ; dxy [mm]; triplets",100, 0.,0.3);
 hnrowB_ph1 = new TH1I("hnrowB_ph1", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
 hnrowB_q1 = new TH1I("hnrowB_q1", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
 hclB_ph1 = new TH1I("hclB_ph1 ", "triplet dx_clphABC90evR RMS out 5; dx [mm];triplets", 200, 0, 1000 );
 hclB_q1 = new TH1I("hclB_q1 ", "triplet dx_clchargeABC90evR RMS out 5; dx [mm];triplets", 100, 0, 50 );;
 
 /*
 hdx3_clchargeABC85evR = new TH1I("dx3_clchargeABC85evR ", "triplet dx_clchargeABC85evR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 85% events in Landau (only high tail)
 hdx3_clphABC85evR = new TH1I("dx3_clphABC85evR ", "triplet dx_clphABC85evR ; dx [mm];triplets", 500, -0.5, 0.5 );
 hdx3_clchargeABC85evR95 = new TH1I("dx3_clchargeABC85evR95 ", "triplet dx_clchargeABC85evR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 85% events in Landau (only high tail)
 hdx3_clphABC85evR95 = new TH1I("dx3_clphABC85evR95 ", "triplet dx_clphABC85evR ; dx [mm];triplets", 500, -0.5, 0.5 );

 hnrowB_ph95_85 = new TH1I("hnrowB_ph95_85", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
 hnrowB_q95_85 = new TH1I("hnrowB_q95_85", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
 hclB_ph95_85 = new TH1I("hclB_ph95_85 ", "triplet dx_clphABC90evR RMS 95; dx [mm];triplets", 200, 0, 1000 );
 hclB_q95_85 = new TH1I("hclB_q95_85 ", "triplet dx_clchargeABC90evR RMS 95; dx [mm];triplets", 100, 0, 50 );;
 hnrowB_ph5_85 = new TH1I("hnrowB_ph5_85", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
 hnrowB_q5_85 = new TH1I("hnrowB_q5_85", "B number of rows ;number of rows [pixels];B clusters on tracks", 40, 0.5, 40.5 );
 hclB_ph5_85 = new TH1I("hclB_ph5_85 ", "triplet dx_clphABC90evR RMS out 5; dx [mm];triplets", 200, 0, 1000 );
 hclB_q5_85 = new TH1I("hclB_q5_85 ", "triplet dx_clchargeABC90evR RMS out 5; dx [mm];triplets", 100, 0, 50 );;
 */

 
 hdx3tree2 = new TH1I("hdx3tree2", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );

 charge_res = new TTree("charge_res","landau and residuals");
 charge_res->Branch("dx3tree",&dx3tree);
 charge_res->Branch("clqAiiitree",&clqAiii);
 charge_res->Branch("clqBiiitree",&clqBiii);
 charge_res->Branch("clqCiiitree",&clqCiii);
 charge_res->Branch("clphAiiitree",&clphAiii);
 charge_res->Branch("clphBiiitree",&clphBiii);
 charge_res->Branch("clphCiiitree",&clphCiii);
 charge_res->Branch("nrowBtree",&nrowBtree);
 charge_res->Branch("dxyCAtree",&dxyCAtree);
 charge_res->Branch("dxytree",&dxytree);
 charge_res->Branch("evt",&evt);
 charge_res->Branch("nclustA",&nclustA);
 charge_res->Branch("nclustB",&nclustB);
 charge_res->Branch("nclustC",&nclustC);
  
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
  hclqAi = new  TH1I( "clqAi", "A isolated cluster charge;cluster charge [ke];A isolatewd clusters",100, 0, 50 );

  hxxAB = new TH2I( "xxAB", "B vs A;row A;row B;clusters", 320, -4, 4, 320, -4, 4 );
  hyyAB = new TH2I( "yyAB", "B vs A;col A;col B;clusters",  80, -4, 4,  80, -4, 4 );

  hdxAB = new  TH1I( "dxAB", "Bx-Ax;x-x [mm];cluster pairs", 800, -2, 2 );
  hdyAB = new  TH1I( "dyAB", "By-Ay;y-y [mm];cluster pairs", 400, -2, 2 );
  dxvsxAB = new  TProfile( "dxvsxAB", "dx vs x A-B;x [mm];<dx> [mm]", 320, -4, 4, -f, f );
  dxvsyAB = new  TProfile( "dxvsyAB", "dx vs y A-B;y [mm];<dx> [mm]",  80, -4, 4, -f, f );

  hdxyCA =  new  TH1I( "dxyCA","dxyCA; dxyCA [mm], events",100,-0.3,0.3);
  hdxy =  new  TH1I( "dxy","dxy; dxy [mm], events",100,-0.3,0.3);

  
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
  nmvsevCB = new  TProfile( "nmvsevCB", "CB matches vs time;time [events];CB matches",3100, 0, 3100*1000, -1, 99 );

  hxC = new  TH1I( "xC", "x C;x [mm];clusters C", 100, -5, 5 );
  hyC = new  TH1I( "yC", "y C;y [mm];clusters C", 100, -5, 5 ); 
  hxCi = new  TH1I( "xCi", "x C isolated;x [mm];isolated clusters C", 100, -5, 5 );
  hyCi = new  TH1I( "yCi", "y C isolated;y [mm];isolated clusters C", 100, -5, 5 ); 
  hclqCi = new  TH1I( "clqCi", "C isolated cluster charge;cluster charge [ke];C isolatewd clusters",	       100, 0, 50 );

  hxxCA = new TH2I( "xxCA", "C vs A;row A;row C;clusters", 320, -4, 4, 320, -4, 4 );
  hyyCA = new TH2I( "yyCA", "C vs A;col A;col C;clusters",  80, -4, 4,  80, -4, 4 );

  hdxCA = new  TH1I( "dxCA", "Cx-Ax;x-x [mm];cluster pairs", 400, -1, 1 );
  hdyCA = new  TH1I( "dyCA", "Cy-Ay;y-y [mm];cluster pairs", 200, -1, 1 );
  
  dxvsxCA = new  TProfile( "dxvsxCA", "dx vs x C-A;x [mm];<dx> [mm]", 320, -4, 4, -f, f );
  dxvsyCA = new  TProfile( "dxvsyCA", "dx vs y C-A;y [mm];<dx> [mm]",  80, -4, 4, -f, f );
  hdyCAc = new  TH1I( "dyCAc", "Cy-Ay, cut dx;y-y [mm];cluster pairs", 200, -1, 1 );
  dyvsyCA = new  TProfile( "dyvsyCA", "dy vs y C-A;y [mm];<dy> [mm]",  80, -4, 4, -f, f );
  nmvsevCA = new  TProfile( "nmvsevCA", "CA matches vs time;time [events];CA matches",3100, 0, 3100*1000, -1, 99 );


  dx3vsx    = new TProfile( "dx3vsx", "dx vs x;x [mm];<dx3> [mm]", 320, -4, 4, -0.5, 0.5 );//turn
  dx3vsy    = new TProfile( "dx3vsy", "dx vs y;y [mm];<dx3> [mm]",  80, -4, 4, -0.5, 0.5 );// rot
  

  hclqAiii = new  TH1I( "clqAiii", "A isolated cluster charge;cluster charge [ke];A isolated clusters",100, 0, 50 );
  hclqBiii = new  TH1I( "clqBiii", "B isolated cluster charge;cluster charge [ke];B isolated clusters",100, 0, 50 );
  hclqCiii = new  TH1I( "clqCiii", "C isolated cluster charge;cluster charge [ke];C isolated clusters",100, 0, 50 );
  hclphAiii = new  TH1I( "clphAiii", "A isolated cluster charge;cluster PH [ADC];A isolated clusters",  200, 0, 1000 );
  hclphBiii = new  TH1I( "clphBiii", "B isolated cluster charge;cluster PH [ADC];B isolated clusters",  200, 0, 1000 );
  hclphCiii = new  TH1I( "clphCiii", "C isolated cluster charge;cluster PH [ADC];C isolated clusters",  200, 0, 1000 );


  
  effvsdxy = new TProfile( "effvsdxy",		     "DUT efficiency vs triplet dxy;xy match radius [mm];DUT efficiency",		     1000, 0, 10, -0.1, 1.1 );

  effvsxy =    new TProfile2D( "effvsxy",		    "DUT efficiency map;x [mm];y[mm];DUT efficiency",		    80, -4, 4, 80, -4, 4, -0.1, 1.1 );
  effvsx = new TProfile( "effvsx", "eff vs x;x [mm];DUT efficiency",		   320, -4, 4, -0.1, 1.1 );
  effvsy = new TProfile( "effvsy", "eff vs y;y [mm];DUT efficiency",		   80, -4, 4, -0.1, 1.1 );
  effvsxm = new TProfile( "effvsxm", "eff vs x mod 50;x mod 50 [#mum];DUT efficiency",		    50, 0, 50, -0.1, 1.1 );

  effvsev = new TProfile( "effvsev", "eff vs time;trigger;DUT efficiency",		    3100, 0, 3100*1000, -0.1, 1.1 );
  effvsiev = new TProfile( "effvsiev", "eff vs event;event mod 200;DUT efficiency",		     100, -0.5, 199.5, -0.1, 1.1 );
  effvsmpxA = new TProfile( "effvsmpxA", "eff vs occupancy A;occupancy A [pixels];DUT efficiency",		      50, 0.5, 50.5, -0.1, 1.1 );
  effvsqA = new TProfile( "effvsqA", "eff vs charge A;cluster charge A [ke];DUT efficiency",		    100, 0, 100, -0.1, 1.1 );
  effvstxy = new TProfile("effvstxy", "eff vs angle;dxy CA [mm];DUT efficiency",		     100, 0, 0.2, -0.1, 1.1 );

  //  hdx3 = new TH1I("dx3", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );



  
}


