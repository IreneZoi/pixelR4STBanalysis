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

#define Iterations 10 //12
#define Thresholds 3

bool print=true;
using namespace std;


void alignmentPlot_automatic()
{

  TString inputDir="../align/";
  TString outputDir="/home/zoiirene/Output/alignment/";

  int Iteration[Iterations];
  TString Label[Thresholds]={"6 %","12 %","20 %"}; 
  TString Run = "2801";
  TString ss_dphcuts[Thresholds] = {"7","15","24"};

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
	TString filename      = inputDir+"align_run"+Run+"_dphcutB"+ss_dphcuts[l]+".000000_iteration_"+inter+".txt";
	cout << filename << endl;
	ifstream stream(filename);

	std::string line;
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
	      alignfA[l][i] = al;
	      alignfAerror[l][i] = er;
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
	      alignfC[l][i] = al;
	      alignfCerror[l][i] = er;
	    }
	    if( dummy =="gamma"){
	      dx3vsx[l][i] =al;
	      dx3vsxError[l][i]=er;
	    }
          } //line scanning                                                                                                                                                                                                                  
	}//file scanning
      }//iteratins
  }//thresholds
  int lastit = 9;
  ofstream myfile;
  myfile.open ("compare_align_2801.txt");
  myfile << "\\begin{tabular}{c|ccc}\n";
  myfile << "\\multicolumn{1}{c|}{Parameter} &\\multicolumn{3}{c}{Offline Threshold [\\%]} \\\\ \n";
  myfile << "Parameter & 6  & 12 & 20 \\\\ \n";
  myfile << "\\hline \n";
  myfile << " f$_{A}^{x}$ [$\\upmu$m] & " << std::setprecision(1) <<alignxA[0][lastit]*1000 << "$\\pm $ "<< std::setprecision(1) <<alignxAerror[0][lastit]*1000 << " & " << std::setprecision(1) << alignxA[1][lastit]*1000 << "$\\pm $ "<< std::setprecision(1) <<alignxAerror[1][lastit]*1000 << " & " << std::setprecision(1) << alignxA[2][lastit]*1000 <<"$\\pm $ " << std::setprecision(1)<<alignxAerror[2][lastit]*1000 << " \\\\ \n";
  myfile << " f$_{A}^{y}$ [$\\upmu$m] & " << std::setprecision(1) <<alignyA[0][lastit]*1000 << "$\\pm $ "<< std::setprecision(1) <<alignyAerror[0][lastit]*1000 << " & " << std::setprecision(1) << alignyA[1][lastit]*1000 << "$\\pm $ "<< std::setprecision(1) <<alignyAerror[1][lastit]*1000 << " & " << std::setprecision(1) << alignyA[2][lastit]*1000 <<"$\\pm $ " << std::setprecision(1)<<alignyAerror[2][lastit]*1000 << " \\\\ \n";
  myfile << " $\\alpha_{A}$ [mrad] & " << std::setprecision(1) <<alignfA[0][lastit]*1000 << "$\\pm $ "<< std::setprecision(1) <<alignfAerror[0][lastit]*1000 << " & " << std::setprecision(1) << alignfA[1][lastit]*1000 << "$\\pm $ "<< std::setprecision(1) <<alignfAerror[1][lastit]*1000 << " & " << std::setprecision(1) << alignfA[2][lastit]*1000 <<"$\\pm $ " << std::setprecision(1)<<alignfAerror[2][lastit]*1000 << " \\\\ \n";
  myfile << " f$_{C}^{x}$ [$\\upmu$m] & " << std::setprecision(1) <<alignxC[0][lastit]*1000 << "$\\pm $ "<< std::setprecision(1) <<alignxCerror[0][lastit]*1000 << " & " << std::setprecision(1) << alignxC[1][lastit]*1000 << "$\\pm $ "<< std::setprecision(1) <<alignxCerror[1][lastit]*1000 << " & " << std::setprecision(1) << alignxC[2][lastit]*1000 <<"$\\pm $ " << std::setprecision(1)<<alignxCerror[2][lastit]*1000 << " \\\\ \n";
  myfile << " f$_{C}^{y}$ [$\\upmu$m] & " << std::setprecision(1) <<alignyC[0][lastit]*1000 << "$\\pm $ "<< std::setprecision(1) <<alignyCerror[0][lastit]*1000 << " & " << std::setprecision(1) << alignyC[1][lastit]*1000 << "$\\pm $ "<< std::setprecision(1) <<alignyCerror[1][lastit]*1000 << " & " << std::setprecision(1) << alignyC[2][lastit]*1000 <<"$\\pm $ " << std::setprecision(1)<<alignyCerror[2][lastit]*1000 << " \\\\ \n";
  myfile << " $\\alpha_{C}$ [mrad] & " << std::setprecision(1) <<alignfC[0][lastit]*1000 << "$\\pm $ "<< std::setprecision(1) <<alignfCerror[0][lastit]*1000 << " & " << std::setprecision(1) << alignfC[1][lastit]*1000 << "$\\pm $ "<< std::setprecision(1) <<alignfCerror[1][lastit]*1000 << " & " << std::setprecision(1) << alignfC[2][lastit]*1000 <<"$\\pm $ " << std::setprecision(1)<<alignfCerror[2][lastit]*1000 << " \\\\ \n";
  myfile << " $\\gamma$ [mrad] & " << std::setprecision(2) <<dx3vsx[0][lastit]*1000 << "$\\pm $ "<< std::setprecision(1) <<dx3vsxError[0][lastit]*1000 << " & " << std::setprecision(2) << dx3vsx[1][lastit]*1000 << "$\\pm $ "<< std::setprecision(1) <<dx3vsxError[1][lastit]*1000 << " & " << std::setprecision(1) << dx3vsx[2][lastit]*1000 <<"$\\pm $ " << std::setprecision(1)<<dx3vsxError[2][lastit]*1000 << " \\\\ \n";
 myfile << "\\end{tabular}\n";
  myfile.close();
  
  // FIXME : add other distributions and print uncertainties in table!
  if(print) cout << "Plotting nrow big scan "  << endl;

  cout << " alignxA " << (alignxA[0][10]-alignxA[1][10])/alignxA[1][10]*100  << " % vs low and " << (alignxA[2][10]-alignxA[1][10])/alignxA[1][10]*100 << " % vs high" << endl;
  DrawMoreTGraphWithError(Iterations, Iteration, alignxA[0], alignxAerror[0],alignxA[1], alignxAerror[1],alignxA[2], alignxAerror[2],"alignxA",Run, Label, "error1","f_{A}^{x} [mm]", 0,Iterations,-1.2,0.1, "alignment","Iterations");
  DrawMoreTGraphWithError(Iterations, Iteration, alignxA[0], alignxAerror[0],alignxA[1], alignxAerror[1],alignxA[2], alignxAerror[2],"alignxA",Run, Label, "zoom","f_{A}^{x} [mm]", 0,Iterations,-0.0003,0.0003, "alignment","Iterations");
  DrawMoreTGraphWithError(Iterations, Iteration, alignyA[0], alignyAerror[0],alignyA[1], alignyAerror[1],alignyA[2], alignyAerror[2],"alignyA",Run, Label, "error1","f_{A}^{y} [mm]", 0,Iterations,-0.25,0.01, "alignment","Iterations");
  DrawMoreTGraphWithError(Iterations, Iteration, alignyA[0], alignyAerror[0],alignyA[1], alignyAerror[1],alignyA[2], alignyAerror[2],"alignyA",Run, Label, "zoom","f_{A}^{y} [mm]", 1,Iterations,-0.002,0.002, "alignment","Iterations");
  DrawMoreTGraphWithError(Iterations, Iteration, alignfA[0], alignfAerror[0],alignfA[1], alignfAerror[1],alignfA[2], alignfAerror[2],"alignfA",Run, Label, "error1","#alpha_{A} [rad]", 0,Iterations,-0.003,0.003, "alignment","Iterations");
  DrawMoreTGraphWithError(Iterations, Iteration, alignfA[0], alignfAerror[0],alignfA[1], alignfAerror[1],alignfA[2], alignfAerror[2],"alignfA",Run, Label, "zoom","#alpha_{A} [rad]", 0,Iterations,-0.0003,0.003, "alignment","Iterations");
    
  DrawMoreTGraphWithError(Iterations, Iteration, alignxC[0], alignxCerror[0],alignxC[1], alignxCerror[1],alignxC[2], alignxCerror[2], "alignxC",Run, Label, "error1","f_{C}^{x} [mm]", 0,Iterations,-1.,0.05, "alignment","Iterations");
  DrawMoreTGraphWithError(Iterations, Iteration, alignxC[0], alignxCerror[0],alignxC[1], alignxCerror[1],alignxC[2], alignxCerror[2], "alignxC",Run, Label, "zoom","f_{C}^{x} [mm]", 1,Iterations,-0.0003,0.0015, "alignment","Iterations");
  DrawMoreTGraphWithError(Iterations, Iteration, alignyC[0], alignyCerror[0],alignyC[1], alignyCerror[1],alignyC[2], alignyCerror[2], "alignyC",Run, Label, "error1","f_{C}^{y} [mm]", 0,Iterations,-0.85,0.02, "alignment","Iterations");
  DrawMoreTGraphWithError(Iterations, Iteration, alignyC[0], alignyCerror[0],alignyC[1], alignyCerror[1],alignyC[2], alignyCerror[2], "alignyC",Run, Label, "zoom","f_{C}^{y} [mm]", 0,Iterations,-0.006,0.002, "alignment","Iterations");
  DrawMoreTGraphWithError(Iterations, Iteration, alignfC[0], alignfCerror[0],alignfC[1], alignfCerror[1],alignfC[2], alignfCerror[2],"alignfC",Run, Label, "error1","#alpha_{C} [rad]", 0,Iterations,-0.005,0.001, "alignment","Iterations");
  DrawMoreTGraphWithError(Iterations, Iteration, alignfC[0], alignfCerror[0],alignfC[1], alignfCerror[1],alignfC[2], alignfCerror[2],"alignfC",Run, Label, "zoom","#alpha_{C} [rad]", 0,Iterations,-0.0005,0.0003, "alignment","Iterations");

  DrawMoreTGraphWithError(Iterations, Iteration, dx3vsx[0], dx3vsxError[0],dx3vsx[1], dx3vsxError[1],dx3vsx[2], dx3vsxError[2], "dx3vsx",Run, Label, "error1","#gamma [rad]", 0,Iterations,-0.0005,0.005, "alignment","Iterations");

}//resolution 
