#define DreiMasterPlanes 3
bool fifty = false;
#define r4sRows 160
#define r4sColumns 155

#define halfSensorX 4.0 //[mm]
#define halfSensorY 3.9 //[mm]
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

typedef std::map<TString, TH1*> histoMap;
histoMap bookControlHists(TString selection, TFile * histofile);
void fillControlHists(histoMap mapOfHists, TString selection, double dx3, double dy3, vector<cluster>::iterator clusterA, vector<cluster>::iterator clusterB, vector<cluster>::iterator clusterC, int ncolB,int nrowB,double xmod,unsigned iev,double xB, double yB,double xAr, double yAr, double xCr, double yCr, double dxCA,double etaA, double etaB, double etaC, TFile * histofile, TString fileName, TH1I *hclph, TH1I *hclq);
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

//for cluster charge cuts
TH1I * hclqAiii;// ( "clqAi", "A isolated cluster charge;cluster charge [ke];A isolatewd clusters", 	       100, 0, 50 );
TH1I * hclqBiii;// ( "clqAi", "A isolated cluster charge;cluster charge [ke];A isolatewd clusters", 	       100, 0, 50 );
TH1I * hclqCiii;// ( "clqAi", "A isolated cluster charge;cluster charge [ke];A isolatewd clusters", 	       100, 0, 50 );
TH1I * hclphAiii;// ( "clqAi", "A isolated cluster charge;cluster charge [ke];A isolatewd clusters", 	       100, 0, 50 );
TH1I * hclphBiii;// ( "clqAi", "A isolated cluster charge;cluster charge [ke];A isolatewd clusters", 	       100, 0, 50 );
TH1I * hclphCiii;// ( "clqAi", "A isolated cluster charge;cluster charge [ke];A isolatewd clusters", 	       100, 0, 50 );


// triplets:
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
