#define DreiMasterPlanes 3
bool fifty = false;
#define r4sRows 160
#define r4sColumns 155

#define halfSensorX 4.0
#define halfSensorY 3.9
#define ACspacing 40 //[mm] distance between planes A and C in the dreimaster

using namespace std;
/* bool PRINT = false; */
/* bool DOALIGNMENT = false; */

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
string datapath = "/mnt/pixeldata/";
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
list < vector < cluster > > oneplane( int plane, string runnum, unsigned Nev, bool fifty, double Tsunami, double dphcut);
void getGain( string gainfile, double (*p0)[r4sColumns][r4sRows], double (*p1)[r4sColumns][r4sRows], double (*p2)[r4sColumns][r4sRows], double (*p3)[r4sColumns][r4sRows], int plane);
void bookHists();
double xcoordinate(int plane, vector<cluster>::iterator c, double align, double pitchc, double pitchr); 
double ycoordinate(int plane, vector<cluster>::iterator c, double align, double pitchc, double pichr); 
double eta(vector<cluster>::iterator c); 
double alignx(TH1I * h);
double aligny(TH1I * h);
double alignangle(TProfile * h);


TProfile phvsprev[DreiMasterPlanes];
TProfile dphvsprev[DreiMasterPlanes];
TH1I hph[DreiMasterPlanes];
TH1I hdph[DreiMasterPlanes];
TH1I hnpx[DreiMasterPlanes];
TH1I hnht[DreiMasterPlanes];
TH2I * hpxmap[DreiMasterPlanes];
TH1I hncl[DreiMasterPlanes];
TH2I * hclmap[DreiMasterPlanes];
TH1I hclsz[DreiMasterPlanes];
TH1I hclph[DreiMasterPlanes];
TH1I hclq[DreiMasterPlanes];



TH1I * hdt; // = new TH1I( "dt", "time between events;log_{10}(#Deltat [s]);events", 100, -4, 1 );
TH1I * hddtAB;//( "ddtAB", "dtA - dtB;A-B #delta#Deltat [clocks];events", 100, -1000, 1000 );
TH1I * hddtCB;// ( "ddtCB", "dtC - dtB;C-B #delta#Deltat [clocks];events", 100, -1000, 1000 );
TH1I * hddtCA;// ( "ddtCA", "dtC - dtA;C-A #delta#Deltat [clocks];events", 100, -1000, 1000 );
TProfile * ddtvsdtAB;// ( "ddtvsdtAB", "A-B time lag vs intervall;log_{10};// (#Deltat [s]);<#delta#Deltat> [clocks]", 100, -4, 1,-1e99, 1e99 );

// correlations:

TH1I * hxA;// ( "xA", "x A;x [mm];clusters A", 100, -5, 5 );
TH1I * hyA;// ( "yA", "y A;y [mm];clusters A", 100, -5, 5 ); 
TH1I * hxAi;// ( "xAi", "x A isolated;x [mm];isolated clusters A", 100, -5, 5 );
TH1I * hyAi;// ( "yAi", "y A isolated;y [mm];isolated clusters A", 100, -5, 5 ); 
TH1I * hclqAi;// ( "clqAi", "A isolated cluster charge;cluster charge [ke];A isolatewd clusters", 	       100, 0, 50 );

TH2I * hxxAB;// = new TH2I;// ( "xxAB", "B vs A;row A;row B;clusters", 320, -4, 4, 320, -4, 4 );
TH2I * hyyAB;// = new TH2I;// ( "yyAB", "B vs A;col A;col B;clusters",  80, -4, 4,  80, -4, 4 );

double f = 0.5;
//double f = 0.1; // aligned

TH1I * hdxAB;// ( "dxAB", "Bx-Ax;x-x [mm];cluster pairs", 800, -2, 2 );
TH1I * hdyAB;// ( "dyAB", "By-Ay;y-y [mm];cluster pairs", 400, -2, 2 );
TProfile * dxvsxAB;// ( "dxvsxAB", "dx vs x A-B;x [mm];<dx> [mm]", 320, -4, 4, -f, f );
TProfile * dxvsyAB;// ( "dxvsyAB", "dx vs y A-B;y [mm];<dx> [mm]",  80, -4, 4, -f, f );

TH2I * hdxvsev;// ( "dxvsev", "Bx-Ax vs events;events;#Deltax [px];clusters",100, 0, 10000, 100, -f, f );

TProfile * nmvsevAB;// ( "nmvsevAB", "AB matches vs time;time [events];AB matches", 3100, 0, 3100*1000, -1, 99 );

TH1I * hxB;// ( "xB", "x B;x [mm];clusters B", 100, -5, 5 );
TH1I * hyB;// ( "yB", "y B;y [mm];clusters B", 100, -5, 5 ); 
TH1I * hxBi;// ( "xBi", "x B isolated;x [mm];isolated clusters B", 100, -5, 5 );
TH1I * hyBi;// ( "yBi", "y B isolated;y [mm];isolated clusters B", 100, -5, 5 ); 
TH1I * hclqBi;// ( "clqBi", "B isolated cluster charge;cluster charge [ke];B isolatewd clusters", 100, 0, 50 );

TH2I * hxxCB = new TH2I;// ( "xxCB", "C vs B;row B;row C;clusters", 320, -4, 4, 320, -4, 4 );
TH2I * hyyCB = new TH2I;// ( "yyCB", "C vs B;col B;col C;clusters",  80, -4, 4,  80, -4, 4 );

TH1I * hdxCB;// ( "dxCB", "Cx-Bx;x-x [mm];cluster pairs", 800, -2, 2 );
TH1I * hdyCB;// ( "dyCB", "Cy-By;y-y [mm];cluster pairs", 400, -2, 2 );
TProfile * dxvsxCB;// ( "dxvsxCB", "dx vs x C-B;x [mm];<dx> [mm]", 320, -4, 4, -f, f );
TProfile * dxvsyCB;// ( "dxvsyCB", "dx vs y C-B;y [mm];<dx> [mm]",  80, -4, 4, -f, f );
TProfile * nmvsevCB;// ( "nmvsevCB", "CB matches vs time;time [events];CB matches",3100, 0, 3100*1000, -1, 99 );

// triplets:

TH1I * hxC;// ( "xC", "x C;x [mm];clusters C", 100, -5, 5 );
TH1I * hyC;// ( "yC", "y C;y [mm];clusters C", 100, -5, 5 ); 
TH1I * hxCi;// ( "xCi", "x C isolated;x [mm];isolated clusters C", 100, -5, 5 );
TH1I * hyCi;// ( "yCi", "y C isolated;y [mm];isolated clusters C", 100, -5, 5 ); 
TH1I * hclqCi;// ( "clqCi", "C isolated cluster charge;cluster charge [ke];C isolatewd clusters",100, 0, 50 );

TH2I * hxxCA = new TH2I;// ( "xxCA", "C vs A;row A;row C;clusters", 320, -4, 4, 320, -4, 4 );
TH2I * hyyCA = new TH2I;// ( "yyCA", "C vs A;col A;col C;clusters",  80, -4, 4,  80, -4, 4 );

TH1I * hdxCA;// ( "dxCA", "Cx-Ax;x-x [mm];cluster pairs", 400, -1, 1 );
TH1I * hdyCA;// ( "dyCA", "Cy-Ay;y-y [mm];cluster pairs", 200, -1, 1 );

TProfile * dxvsxCA;// ( "dxvsxCA", "dx vs x C-A;x [mm];<dx> [mm]", 320, -4, 4, -f, f );
TProfile * dxvsyCA;// ( "dxvsyCA", "dx vs y C-A;y [mm];<dx> [mm]",  80, -4, 4, -f, f );
TH1I * hdyCAc;// ( "dyCAc", "Cy-Ay, cut dx;y-y [mm];cluster pairs", 200, -1, 1 );
TProfile * dyvsyCA;// ( "dyvsyCA", "dy vs y C-A;y [mm];<dy> [mm]",  80, -4, 4, -f, f );
TProfile * nmvsevCA;// ( "nmvsevCA", "CA matches vs time;time [events];CA matches", 3100, 0, 3100*1000, -1, 99 );

TH1I * hdx3;// ( "dx3", "triplet dx;dx [mm];triplets", 500, -0.5, 0.5 );
TH1I * hdy3;// ( "dy3", "triplet dy;dy [mm];triplets", 200, -1, 1 );

TH1I * hdx3c;// ( "dx3c", "triplet dx, cut dy;dx [mm];triplets", 500, -0.25, 0.25 );
TH1I * hdx3ci;// ( "dx3ci", "triplet dx, cut dy, isolated;dx [mm];isolated triplets", 500, -0.25, 0.25 );
TH1I * hdx3cii;// ( "dx3cii", "triplet dx, cut dy, isolated;dx [mm];isolated triplets", 500, -0.25, 0.25 );
TH1I * hdx3ciii;// ( "dx3ciii", "triplet dx, cut dy, isolated;dx [mm];isolated triplets", 500, -0.25, 0.25 );
TH1I * hdx3ciiiqr;// ( "dx3ciiiqr", "triplet dx, cut dy, isolated, qB < qR;dx [mm];isolated triplets", 500, -0.25, 0.25 );
TH1I * hdx3ciiiq;// ( "dx3ciiiq", "triplet dx, cut dy, isolated, qL < qB < qR;dx [mm];isolated triplets", 500, -0.25, 0.25 );
TH1I * hdx3ciiiqr3;// ( "dx3ciiiqr3", "triplet dx, cut dy, isolated, q3 < qR;dx [mm];isolated triplets", 500, -0.25, 0.25 );
TH1I * hdx3ciiiq3;// ( "dx3ciiiq3", "triplet dx, cut dy, isolated, qL < q3 < qR;dx [mm];isolated triplets", 500, -0.25, 0.25 );



TH1I * hdx3c1;// ( "dx3c1", "triplet dx, cut dy, npx 1;dx [mm];triplets, B npx 1",500, -0.25, 0.25 );
TH1I * hdx3c2;// ( "dx3c2", "triplet dx, cut dy, npx 2;dx [mm];triplets, B npx 2",500, -0.25, 0.25 );
TH1I * hdx3c3;// ( "dx3c3", "triplet dx, cut dy, npx 3;dx [mm];triplets, B npx 3",500, -0.25, 0.25 );
TH1I * hdx3c4;// ( "dx3c4", "triplet dx, cut dy, npx 4;dx [mm];triplets, B npx 4",500, -0.25, 0.25 );
TH1I * hdx3c5;// ( "dx3c5", "triplet dx, cut dy, npx 5;dx [mm];triplets, B npx 5",500, -0.25, 0.25 );
TH1I * hdx3c6;// ( "dx3c6", "triplet dx, cut dy, npx 6;dx [mm];triplets, B npx 6",500, -0.25, 0.25 );
TH1I * hdx3c7;// ( "dx3c7", "triplet dx, cut dy, npx > 6;dx [mm];triplets, B npx > 6",500, -0.25, 0.25 );

TH1I * hdx3m;// ( "dx3m", "triplet dx, x < 0;dx [mm];triplets", 500, -0.5, 0.5 );
TH1I * hdx3p;// ( "dx3p", "triplet dx, x > 0;dx [mm];triplets", 500, -0.5, 0.5 );
TH1I * hdx3ct;// ( "dx3ct", "triplet dx, cut dy, tx;dx [mm];triplets",500, -0.25, 0.25 );
TProfile * madx3vsq;// ( "madx3vsq", "MAD;// (dx) vs Q;B cluster charge [ke];MAD dx [mm]", 100, 0, 100, 0, 0.1 );
TProfile * madx3vsn;// ( "madx3vsn", "MAD;// (dx) vs cluster size;B cluster size [pixels];MAD dx [mm]", 20, 0.5, 20.5, 0, 0.1 );

TH1I * hdx3cq;// ( "dx3cq", "triplet dx, Landau peak;dx [mm];Landau peak triplets",500, -0.25, 0.25 );
TH1I * hdx3cqi;// ( "dx3cqi", "triplet dx, Landau peak, isolated;dx [mm];isolated Landau peak triplets",500, -0.25, 0.25 );
TH1I * hdx3cq3;// ( "dx3cq3", "triplet dx, 3 Landau peak;dx [mm];Landau peak triplets",		500, -0.25, 0.25 );
TH1I * hdx3nocq3;// ( "dx3nocq3", "triplet dx, no 3 Landau peak;dx [mm];no Landau peak triplets",	500, -0.25, 0.25 );
TH1I * hdx3cq3i;// ( "dx3cq3i", "triplet dx, 3 Landau peak, isolated;dx [mm];isolated Landau peak triplets",		 500, -0.25, 0.25 );

TProfile * dx3vsev;// ( "dx3vsev", "dx3 vs time;trigger;<dx3> [mm]",		    310, 0, 3100*1000, -0.5, 0.5 );

TProfile * dx3vsx;// ( "dx3vsx", "dx vs x;x [mm];<dx3> [mm]", 320, -4, 4, -0.5, 0.5 );
TProfile * dx3vsy;// ( "dx3vsy", "dx vs y;y [mm];<dx3> [mm]",  80, -4, 4, -0.5, 0.5 );
TProfile * dx3vsxm;// ( "dx3vsxm", "dx vs x mod 50 um;x mod 50 [#mum];<dx3> [mm]",		    50, 0, 50, -0.5, 0.5 );

TProfile * madx3vsdx;// ( "madx3vsdx", "MAD;// (dx3) vs dx C-A;C-A dx [#mum];MAD dx3 [mm]",		      100, -100, 100, 0, 0.1 );

TH1I * hdx3cq3t;// ( "dx3cq3t",		 "triplet dx, 3 Landau peak, forward;dx [mm];Landau peak forward triplets",		 500, -0.25, 0.25 );
TProfile * madx3vsx;// ( "madx3vsx", "MAD;// (dx3) vs x;x [mm];MAD dx3 [mm]", 320, -4, 4, 0, 0.1 );
TProfile * madx3vsy;// ( "madx3vsy", "MAD;// (dx3) vs y;y [mm];MAD dx3 [mm]",  80, -4, 4, 0, 0.1 );
TProfile * madx3vsxm;// ( "madx3vsxm", "MAD;// (dx3) vs xmod;x mod 50 [#mum];MAD dx3 [mm]",	      50, 0, 50, 0, 0.1 );

TProfile * etavsxmB3;// ( "etavsxmB3", "eta vs xmod;x mod 50 [#mum];B <eta>",		      50, 0, 50, -1.1, 1.1 );
  TProfile * madx3vseta;// ( "madx3vseta", "MAD;// (dx3) vs eta;eta;MAD dx3 [mm]",		       100, -1, 1, 0, 0.1 );
TH1I * hdx3cq3t2;// ( "dx3cq3t2",	  "triplet dx, 3 Landau peak, forward, 2-px;dx [mm];Landau peak forward triplets",		  500, -0.25, 0.25 );

TH2I * hclmapB3;// ( "clmapB3", "linked cluster map B;col;row;B clusters on tracks",	  80, 0, 80, 320, 0, 320 );

TH1I * hxA3;// ( "xA3", "x A linked;x [mm];A clusters on tracks", 100, -5, 5 );
TH1I * hyA3;// ( "yA3", "y A linked;y [mm];A clusters on tracks", 100, -5, 5 ); 
TH1I * hxB3;// ( "xB3", "x B linked;x [mm];B clusters on tracks", 100, -5, 5 );
TH1I * hyB3;// ( "yB3", "y B linked;y [mm];B clusters on tracks", 100, -5, 5 ); 
TH1I * hxC3;// ( "xC3", "x C linked;x [mm];C clusters on tracks", 100, -5, 5 );
TH1I * hyC3;// ( "yC3", "y C linked;y [mm];C clusters on tracks", 100, -5, 5 ); 

TH1I * hclszAiii;// ( "clszA3", "A cluster size on tracks;cluster size [pixels];Aclusters on tracks",		40, 0.5, 40.5 );
TH1I * hclphAiii;// ( "clphA3", "A cluster PH on tracks;cluster ph [ADC];A clusters on tracks",		200, 0, 1000 );
TH1I * hclqAiii;// ( "clqA3", "A cluster charge on tracks;cluster charge [ke];A clusters on tracks",	       160, 0, 80 );
TH1I * hclszBiii;// ( "clszB3", "B cluster size on tracks;cluster size [pixels];B clusters on tracks",		40, 0.5, 40.5 );
TH1I * hncolBiii;// ( "ncolB3", "B cluster size on tracks;cluster size [columns];B clusters on tracks",		20, 0.5, 20.5 );
TH1I * hnrowBiii;// ( "nrowB3", "B cluster size on tracks;cluster size [rows];B clusters on tracks",		20, 0.5, 20.5 );
TH1I * hclphBiii;// ( "clphB3", "B cluster PH on tracks;cluster ph [ADC];B clusters on tracks",		200, 0, 1000 );
TH1I * hclqBiii;// ( "clqB3", "B cluster charge on tracks;cluster charge [ke];B clusters on tracks",	       160, 0, 80 );
TH1I * hclszCiii;// ( "clszC3", "C cluster size on tracks;cluster size [pixels];C clusters on tracks",		40, 0.5, 40.5 );
TH1I * hclphCiii;// ( "clphC3", "C cluster PH on tracks;cluster ph [ADC];C clusters on tracks",		200, 0, 1000 );
TH1I * hclqCiii;// ( "clqC3", "C cluster charge on tracks;cluster charge [ke];C clusters on tracks",	       160, 0, 80 );




TH1I * hclszA3;// ( "clszA3", "A cluster size on tracks;cluster size [pixels];Aclusters on tracks",		40, 0.5, 40.5 );
TH1I * hclphA3;// ( "clphA3", "A cluster PH on tracks;cluster ph [ADC];A clusters on tracks",		200, 0, 1000 );
TH1I * hclqA3;// ( "clqA3", "A cluster charge on tracks;cluster charge [ke];A clusters on tracks",	       160, 0, 80 );

TH1I * hclszB3;// ( "clszB3", "B cluster size on tracks;cluster size [pixels];B clusters on tracks",		40, 0.5, 40.5 );
TH1I * hncolB3;// ( "ncolB3", "B cluster size on tracks;cluster size [columns];B clusters on tracks",		20, 0.5, 20.5 );
TH1I * hnrowB3;// ( "nrowB3", "B cluster size on tracks;cluster size [rows];B clusters on tracks",		20, 0.5, 20.5 );
TH1I * hclphB3;// ( "clphB3", "B cluster PH on tracks;cluster ph [ADC];B clusters on tracks",		200, 0, 1000 );
TH1I * hclqB3;// ( "clqB3", "B cluster charge on tracks;cluster charge [ke];B clusters on tracks",	       160, 0, 80 );
TH1I * hclqB3i;// ( "clqB3i", "B cluster charge on tracks;cluster charge [ke];B clusters on tracks",	       80, 0, 20 );
TH1I * hclqB3n;// ( "clqB3n",		"B cluster charge on tracks, npx < 4;cluster charge [ke];B clusters on tracks, npx < 4",		160, 0, 80 );

TH1I * hclszC3;// ( "clszC3", "C cluster size on tracks;cluster size [pixels];C clusters on tracks",		40, 0.5, 40.5 );
TH1I * hclphC3;// ( "clphC3", "C cluster PH on tracks;cluster ph [ADC];C clusters on tracks",		200, 0, 1000 );
TH1I * hclqC3;// ( "clqC3", "C cluster charge on tracks;cluster charge [ke];C clusters on tracks",	       160, 0, 80 );

TProfile * nrowvsxmB3;// ( "nrowvsxmB3",		       "B rows vs xmod;x mod 50 [#mum];<B cluster size [rows]>",		       50, 0, 50, 0.5, 10.5 );
  TProfile * clqvsxmB3;// ( "clqvsxmB3",		      "B cluster charge vs xmod;x mod 50 [#mum];<B cluster charge [ke]>",		      50, 0, 50, 0, 50 );

TH1I * hetaA3;// ( "etaA3", "A cluster eta;eta;A 2-pix clusters on tracks",	       100, -1, 1 );
TH1I * hetaB3;// ( "etaB3", "B cluster eta;eta;B 2-pix clusters on tracks",	       100, -1, 1 );
TH1I * hetaC3;// ( "etaC3", "C cluster eta;eta;C 2-pix clusters on tracks",	       100, -1, 1 );

TH1I * hpxqA3;// ( "pxqA3", "A pixel charge;pixel charge [ke];A pixels on tracks",	       100, 0, 20 );
TH1I * hpxqB3;// ( "pxqB3", "B pixel charge;pixel charge [ke];B pixels on tracks",	       100, 0, 20 );
TH1I * hpxqC3;// ( "pxqC3", "C pixel charge;pixel charge [ke];C pixels on tracks",	       100, 0, 20 );

TH1I * hpxpA3;// ( "pxpA3", "A pixel PH;pixel PH [ADC];A pixels on tracks", 250, 0, 500 );
TH1I * hpxpB3;// ( "pxpB3", "B pixel PH;pixel PH [ADC];B pixels on tracks", 250, 0, 500 );
TH1I * hpxpC3;// ( "pxpC3", "C pixel PH;pixel PH [ADC];C pixels on tracks", 250, 0, 500 );
TH1I * hpxq1stB3;// ( "pxq1stB3", "B 1st pixel charge;pixel charge [ke];B `st pixels on tracks",100, 0, 20 );
TH1I * hpxq2ndB3;// ( "pxq2ndB3", "B 2nd pixel charge;pixel charge [ke];B `st pixels on tracks",		  100, 0, 20 );

TProfile * effvsdxy;// ( "effvsdxy",		     "DUT efficiency vs triplet dxy;xy match radius [mm];DUT efficiency",		     1000, 0, 10, -0.1, 1.1 );

TProfile2D * effvsxy;// =    new TProfile2D;// ( "effvsxy",		    "DUT efficiency map;x [mm];y[mm];DUT efficiency",		    80, -4, 4, 80, -4, 4, -0.1, 1.1 );
TProfile * effvsx;// ( "effvsx", "eff vs x;x [mm];DUT efficiency",		   320, -4, 4, -0.1, 1.1 );
TProfile * effvsy;// ( "effvsy", "eff vs y;y [mm];DUT efficiency",		   80, -4, 4, -0.1, 1.1 );
TProfile * effvsxm;// ( "effvsxm", "eff vs x mod 50;x mod 50 [#mum];DUT efficiency",		    50, 0, 50, -0.1, 1.1 );

TProfile * effvsev;// ( "effvsev", "eff vs time;trigger;DUT efficiency",		    3100, 0, 3100*1000, -0.1, 1.1 );
TProfile * effvsiev;// ( "effvsiev", "eff vs event;event mod 200;DUT efficiency",		     100, -0.5, 199.5, -0.1, 1.1 );
TProfile * effvsmpxA;// ( "effvsmpxA", "eff vs occupancy A;occupancy A [pixels];DUT efficiency",		      50, 0.5, 50.5, -0.1, 1.1 );
TProfile * effvsqA;// ( "effvsqA", "eff vs charge A;cluster charge A [ke];DUT efficiency",		    100, 0, 100, -0.1, 1.1 );
TProfile * effvstxy;// ( "effvstxy", "eff vs angle;dxy CA [mm];DUT efficiency",		     100, 0, 0.2, -0.1, 1.1 );
