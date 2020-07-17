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

#define Iterations 6

bool print=true;
using namespace std;


void alignmentPlot_2743()
{


  TString outputDir="/home/zoiirene/Output/alignment/";

  TString detectorA="146";
  TString detectorB="148";
  TString detectorC="163";
  TString labelA="FDB150Y_2_R4S100x25-Y2_2, thr 12 ADC";
  TString labelB="FDB150P_12_R4S100x25-P4_1, 120V, thr 12 ADC";
  TString labelC="FDB150Y_2_R4S100x25-Y6_1, thr 12 ADC";
  TString info = "beam energy 5.6 GeV";




  int Iteration[Iterations]={0,1,2,3,4,5};//,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
   
  TString Run = "2743";


  Float_t alignxA[Iterations];
  Float_t alignxAerror[Iterations]; 
  alignxA[0] =  0.303518;
  alignxAerror[0] = 0.000106696; 
  alignxA[1] = 1.33904e-06;
  alignxAerror[1] = 0.000106658; 
  alignxA[2] = -0.000540149;
  alignxAerror[2] = 0.000101909; 
  alignxA[3] = -1.42524e-05;
  alignxAerror[3] = 0.000101902; 
  alignxA[4] = 6.26049e-07;
  alignxAerror[4] = 0.000101883;
  alignxA[5] = -1.30954e-06;
  alignxAerror[5] = 0.000101894; 

  Float_t alignyA[Iterations]; 
  Float_t alignyAerror[Iterations]; 
  alignyA[0] = -0.114475;
  alignyAerror[0] = 0.0010123; 
  alignyA[1] = -0.00145694;
  alignyAerror[1] = 0.00101139; 
  alignyA[2] = 0.000870129;
  alignyAerror[2] = 0.00101055; 
  alignyA[3] = 9.20373e-06;
  alignyAerror[3] = 0.00101054;
  
  alignyA[4] = 1.94264e-07;
  alignyAerror[4] = 0.00101054;
  alignyA[5] = 1.96467e-07;
  alignyAerror[5] = 0.00101054; 

  Float_t alignfA[Iterations]; 
  Float_t alignfAerror[Iterations]; 
  alignfA[0] = -0.00237107;
  alignfAerror[0] = 0.000155721; 
  alignfA[1] = -0.00217713;
  alignfAerror[1] = 0.000152344; 
  alignfA[2] = -2.56616e-05;
  alignfAerror[2] = 0.000152344; 
  alignfA[3] = -7.06145e-08;
  alignfAerror[3] = 0.000152345; 
  alignfA[4] = -8.02762e-10;
  alignfAerror[4] = 0.000152345;
  alignfA[5] = -8.02762e-10;
  alignfAerror[5] = 0.000152345; 

  Float_t alignxC[Iterations]; 
  Float_t alignxCerror[Iterations]; 
  alignxC[0] = 0.230829;
  alignxCerror[0] = 0.000114747; 
  alignxC[1] = -1.91685e-05;
  alignxCerror[1] = 0.000115324; 
  alignxC[2] = 0.000761364;
  alignxCerror[2] = 0.000100212; 
  alignxC[3] = 4.74637e-06;
  alignxCerror[3] = 9.97433e-05;
  alignxC[4] = 2.15038e-07;
  alignxCerror[4] = 9.95396e-05;
  alignxC[5] = 2.15038e-07;
  alignxCerror[5] = 9.95396e-05; 

  Float_t alignyC[Iterations]; 
  Float_t alignyCerror[Iterations]; 
  alignyC[0] = -0.297691;
  alignyCerror[0] = 0.0012313; 
  alignyC[1] = -0.00442152;
  alignyCerror[1] = 0.00123523; 
  alignyC[2] = -0.00179324;
  alignyCerror[2] = 0.0012294; 
  alignyC[3] = -2.50259e-05;
  alignyCerror[3] = 0.00122941; 
  alignyC[4] = -8.27905e-08;
  alignyCerror[4] = 0.00122941;
  alignyC[5] = -8.27905e-08;
  alignyCerror[5] = 0.00122941; 

  Float_t alignfC[Iterations]; 
  Float_t alignfCerror[Iterations]; 
  alignfC[0] = 0.00376352;
  alignfCerror[0] = 0.000174656; 
  alignfC[1] = 0.00414759;
  alignfCerror[1] = 0.000169144; 
  alignfC[2] =  5.49711e-05;
  alignfCerror[2] = 0.000169901; 
  alignfC[3] = 8.3873e-08;
  alignfCerror[3] = 0.000169904; 
  alignfC[4] = 3.99351e-09;
  alignfCerror[4] = 0.000169904;
  alignfC[5] = 3.99351e-09;
  alignfCerror[5] = 0.000169904; 

  Float_t dx3vsx[Iterations];
  Float_t dx3vsxError[Iterations];
  dx3vsx[0] = 0.00110718;
  dx3vsxError[0] = 9.95385e-05;
  dx3vsx[1] = 0.00123891;
  dx3vsxError[1] = 9.76266e-05;
  dx3vsx[2] = 0.00122139;
  dx3vsxError[2] = 9.84464e-05;
  dx3vsx[3] = 0.00122484;
  dx3vsxError[3] = 9.84865e-05;
  dx3vsx[4] = 1.59129e-08;
  dx3vsxError[4] = 9.8478e-05;
  dx3vsx[5] = 3.31324e-11;
  dx3vsxError[5] = 9.8478e-05;


  if(print) cout << "Plotting nrow big scan "  << endl;

  DrawTGraphWithError(Iterations, Iteration, alignxA, alignxAerror, "alignxA",Run, "4it", "error1","xAr-xB", 0,Iterations,-0.5,0.5, "alignment","iteration");
  DrawTGraphWithError(Iterations, Iteration, alignxA, alignxAerror, "alignxA",Run, "4it", "zoom","xAr-xB", 0,Iterations,-0.0007,0.0003, "alignment","iteration");
  DrawTGraphWithError(Iterations, Iteration, alignyA, alignyAerror, "alignyA",Run, "4it", "error1","yAr-yB", 0,Iterations,-0.2,0.1, "alignment","iteration");
  DrawTGraphWithError(Iterations, Iteration, alignyA, alignyAerror, "alignyA",Run, "4it", "zoom","yAr-yB", 0,Iterations,-0.004,0.004, "alignment","iteration");
  DrawTGraphWithError(Iterations, Iteration, alignfA, alignfAerror, "alignfA",Run, "4it", "error1","slope dxAB vs yB", 0,Iterations,-0.003,0.003, "alignment","iteration");
  DrawTGraphWithError(Iterations, Iteration, alignfA, alignfAerror, "alignfA",Run, "4it", "zoom","slope dxAB vs yB", 0,Iterations,-0.0007,0.0003, "alignment","iteration");
    
  DrawTGraphWithError(Iterations, Iteration, alignxC, alignxCerror, "alignxC",Run, "4it", "error1","xCr-xB", 0,Iterations,-0.1,0.3, "alignment","iteration");
  DrawTGraphWithError(Iterations, Iteration, alignxC, alignxCerror, "alignxC",Run, "4it", "zoom","xCr-xB", 0,Iterations,-0.0002,0.001, "alignment","iteration");
  DrawTGraphWithError(Iterations, Iteration, alignyC, alignyCerror, "alignyC",Run, "4it", "error1","yCr-xB", 0,Iterations,-0.4,0.1, "alignment","iteration");
  DrawTGraphWithError(Iterations, Iteration, alignyC, alignyCerror, "alignyC",Run, "4it", "zoom","yCr-xB", 0,Iterations,-0.01,0.003, "alignment","iteration");
  DrawTGraphWithError(Iterations, Iteration, alignfC, alignfCerror, "alignfC",Run, "4it", "error1","slope dxCB vs yB", 0,Iterations,-0.005,0.005, "alignment","iteration");
  DrawTGraphWithError(Iterations, Iteration, alignfC, alignfCerror, "alignfC",Run, "4it", "zoom","slope dxCB vs yB", 0,Iterations,-0.0007,0.0003, "alignment","iteration");

  DrawTGraphWithError(Iterations, Iteration, dx3vsx, dx3vsxError, "dx3vsx",Run, "4it", "error1","slope dx3 vs xB", 0,Iterations,-0.002,0.002, "alignment","iteration");

}//resolution 
