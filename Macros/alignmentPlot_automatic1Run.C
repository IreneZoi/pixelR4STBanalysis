#include <stdio.h>
#include "TCanvas.h"
#include "TF1.h"
#include <TTreeFormula.h>
#include <TMath.h>
#include <Rtypes.h>
#include <TROOT.h>
#include <TVectorT.h>
#include <TMap.h>
#include <TAxis.h>
#include <sstream>
#include "fileHandler.h"

#define Iterations 8
#define Thresholds 1

bool print=true;
using namespace std;


void alignmentPlot_automatic1Run()
{

  TString inputDir="../align/";
  TString outputDir="/home/zoiirene/Output/alignment/";

  int Iteration[Iterations];
   
  TString Run = "2743";
  TString ss_dphcuts[Thresholds] = {"12"};

  Float_t alignxA[Thresholds][Iterations];
  Float_t alignxAerror[Thresholds][Iterations]; 
  Float_t alignyA[Thresholds][Iterations]; 
  Float_t alignyAerror[Thresholds][Iterations]; 
  Float_t alignfA[Thresholds][Iterations]; 
  Float_t alignfAerror[Thresholds][Iterations]; 
  Float_t alignxC[Thresholds][Iterations]; 
  Float_t alignxCerror[Thresholds][Iterations]; 
  Float_t alignyC[Thresholds][Iterations]; 
  Float_t alignyCerror[Thresholds][Iterations]; 
  Float_t alignfC[Thresholds][Iterations]; 
  Float_t alignfCerror[Thresholds][Iterations]; 
  Float_t dx3vsx[Thresholds][Iterations];
  Float_t dx3vsxError[Thresholds][Iterations];




  
  for(int l = 0; l< Thresholds; l++){
      cout << "dphcut  "<< ss_dphcuts[l] << endl;
	alignxA[l][0] = 0;
	alignxAerror[l][0] = 0; 
	alignyA[l][0] = 0; 
	alignyAerror[l][0] = 0; 
	alignfA[l][0] = 0; 
	alignfAerror[l][0] = 0; 
	alignxC[l][0] = 0; 
	alignxCerror[l][0] = 0; 
	alignyC[l][0] = 0; 
	alignyCerror[l][0] = 0; 
	alignfC[l][0] = 0; 
	alignfCerror[l][0] = 0; 
	dx3vsx[l][0] = 0;
	dx3vsxError[l][0] = 0;
	Iteration[0]=0;
      for(int i = 1; i< Iterations; i++){
	alignxA[l][i] = 0;
	alignxAerror[l][i] = 0; 
	alignyA[l][i] = 0; 
	alignyAerror[l][i] = 0; 
	alignfA[l][i] = 0; 
	alignfAerror[l][i] = 0; 
	alignxC[l][i] = 0; 
	alignxCerror[l][i] = 0; 
	alignyC[l][i] = 0; 
	alignyCerror[l][i] = 0; 
	alignfC[l][i] = 0; 
	alignfCerror[l][i] = 0; 
	dx3vsx[l][i] = 0;
	dx3vsxError[l][i] = 0;

	Iteration[i]=i;

	TString inter;
	inter.Form("%d",i-1);
	TString filename      = inputDir+"align_run"+Run+"_dphcutB"+ss_dphcuts[0]+".000000_iteration_"+inter+".txt";
	cout << filename << endl;
	ifstream stream(filename);
	
	if(!stream.is_open())     {
	    cout << " File " << filename << " not opened" << endl;
	}
	else{
	  TString dummy;
	  Float_t al,er;
	  string dummyLine;
	  getline(stream, dummyLine);
	  
	  while (stream >> dummy >> al >> er){
	    cout << dummy << " " << al << " " << er << endl;
	    if( dummy =="alignxA"){
	      alignxA[l][i] = al;
	      alignxAerror[l][i] = er;
	    }
	    if( dummy =="alignyA"){
	      alignyA[l][i] = al;
	      alignyAerror[l][i] = er;
	    }
	    if( dummy =="alignfA"){
	      alignfA[l][i] = al*1000;
	      alignfAerror[l][i] = er*1000;
	      }
	    if( dummy =="alignxC"){
	      alignxC[l][i] = al;
	      alignxCerror[l][i] = er;
	    }
	    if( dummy =="alignyC"){
	      alignyC[l][i] = al;
	      alignyCerror[l][i] = er;
	    }
	    if( dummy =="alignfC"){
	      alignfC[l][i] = al*1000;
	      alignfCerror[l][i] = er*1000;
	    }
	    if( dummy =="gamma"){
	      dx3vsx[l][i] =al*1000;
	      dx3vsxError[l][i]=er*1000;
	    }
	  } //line scanning  
	}//file scanning
      }//iteratins
  }//thresholds


  cout << " test content "<<endl;
  for(int l = 0; l< Thresholds; l++){
    cout << "dphcut  "<< ss_dphcuts[l] << endl;
    for(int i = 0; i< Iterations; i++){
      cout << i << " "<< alignxA[l][i] << endl;
    }
  }
  
  if(print) cout << "Plotting nrow big scan "  << endl;

  DrawTGraphWithError(Iterations, Iteration, alignxA[0], alignxAerror[0], "alignxA",Run, "thesis", "error1","f_{A}^{x} [mm]", 0,Iterations,-0.5,0.5, "alignment","Iterations");
  DrawTGraphWithError(Iterations, Iteration, alignxA[0], alignxAerror[0], "alignxA",Run, "thesis", "zoom","f_{A}^{x} [mm]", 0,Iterations,-0.0007,0.0003, "alignment","Iterations");
  DrawTGraphWithError(Iterations, Iteration, alignyA[0], alignyAerror[0], "alignyA",Run, "thesis", "error1","f_{A}^{y} [mm]", 0,Iterations,-0.2,0.1, "alignment","Iterations");
  DrawTGraphWithError(Iterations, Iteration, alignyA[0], alignyAerror[0], "alignyA",Run, "thesis", "zoom","f_{A}^{y} [mm]", 0,Iterations,-0.004,0.004, "alignment","Iterations");
  DrawTGraphWithError(Iterations, Iteration, alignfA[0], alignfAerror[0], "alignfA",Run, "thesis", "error1","#alpha_{A} [mrad]", 0,Iterations,-3,3, "alignment","Iterations");
  DrawTGraphWithError(Iterations, Iteration, alignfA[0], alignfAerror[0], "alignfA",Run, "thesis", "zoom","#alpha_{A} [mrad]", 0,Iterations,-0.7,0.3, "alignment","Iterations");
    
  DrawTGraphWithError(Iterations, Iteration, alignxC[0], alignxCerror[0], "alignxC",Run, "thesis", "error1","f_{C}^{x} [mm]", 0,Iterations,-0.05,0.25, "alignment","Iterations");
  DrawTGraphWithError(Iterations, Iteration, alignxC[0], alignxCerror[0], "alignxC",Run, "thesis", "zoom","f_{C}^{x} [mm]", 1,Iterations,-0.0002,0.001, "alignment","Iterations");
  DrawTGraphWithError(Iterations, Iteration, alignyC[0], alignyCerror[0], "alignyC",Run, "thesis", "error1","f_{C}^{y} [mm]", 0,Iterations,-0.4,0.1, "alignment","Iterations");
  DrawTGraphWithError(Iterations, Iteration, alignyC[0], alignyCerror[0], "alignyC",Run, "thesis", "zoom","f_{C}^{y} [mm]", 0,Iterations,-0.01,0.003, "alignment","Iterations");
  DrawTGraphWithError(Iterations, Iteration, alignfC[0], alignfCerror[0], "alignfC",Run, "thesis", "error1","#alpha_{C} [mrad]", 0,Iterations,-1,5, "alignment","Iterations");
  DrawTGraphWithError(Iterations, Iteration, alignfC[0], alignfCerror[0], "alignfC",Run, "thesis", "zoom","#alpha_{C} [mrad]", 0,Iterations,-0.2,0.2, "alignment","Iterations");

  DrawTGraphWithError(Iterations, Iteration, dx3vsx[0], dx3vsxError[0], "dx3vsx",Run, "thesis", "error1","#gamma [mrad]", 0,Iterations,-0.1,1.5, "alignment","Iterations");

}//resolution 
