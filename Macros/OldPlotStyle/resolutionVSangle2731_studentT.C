#include <stdio.h>
#include "TCanvas.h"
#include "TF1.h"
#include <TTreeFormula.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <Rtypes.h>
#include <TROOT.h>
#include <TVectorT.h>
#include <TMap.h>
#include <TAxis.h>
#include <sstream>
#include "fileHandler.h"

#define Angles 29 //32
bool print=true;
#define dphcuts  1
#define comparisons 1
using namespace std;
TString outputDir = "/home/zoiirene/Output/Plots/";
  
void resolutionVSangle2731_studentT(TString function = "RMS")
{

  TString detectorA="146";
  TString detectorB="148";
  TString detectorC="163";
  TString labelA="Pspray_default_FDB, thr 12 ADC";
  TString labelB="Pstop_RD53Apads_FDB, 120V, thr 12 ADC";
  TString labelC="Pspray_RD53Apads_FDB, thr 12 ADC";
  TString info = "beam energy 5.6 GeV";
  
  
  //  Double_t Angle[Angles]={-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,9.5,10,10.5,11,11.5,12,13,14,15,16,17,18,19,20,21,22,23,24};
  Double_t Angle[Angles]={-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
  Double_t AngleError[Angles];
  Double_t TanAngle[Angles];
  Double_t TanAngleError[Angles];
  Double_t Resolution[Angles];
  Double_t Ncol[Angles];
  Double_t NcolError[Angles];
  Double_t ResolutionError[Angles];
  TString ss_Angle[Angles];
  
  TH1F * h_res;//[Angles];

  Int_t run[Angles]={2731,2732,2733,2734,2735,2736,2737,2738,2739,2740,2741,2743,2744,2745,2746,2747,2748,2749,2750,2751,2752,2753,2754,2755,2756,2757,2758,2759,2760};//,2764,2766,2767};
  TString Run[Angles];
  ostringstream strs[Angles];

  float angleunit = 1.25;
  if(print) cout << "Getting files "  << endl;
  for(int i=0; i<Angles; i++)
    {
      if(print) cout << "Run: "<< i <<" " << run[i] << endl;
      Run[i].Form("%d",run[i]);
      if(print) cout << "Run: " << Run << endl;
      ss_Angle[i].Form("%d",i);
      Angle[i]=Angle[i]*angleunit; //grad
      TanAngle[i]=TMath::Tan(Angle[i]*TMath::Pi()/((Double_t)180.));
      cout << " Angle " <<  Angle[i] << " TanAngle " <<  TanAngle[i] << endl;
      AngleError[i]=1; //grad
      TanAngleError[i]=(1+TMath::Power(TMath::Tan(AngleError[i]*TMath::Pi()/((Double_t) 180.)),2))*AngleError[i]*TMath::Pi()/((Double_t) 180.); //grad

      if(print) cout << "Angle " << i<< ": " << Angle[i] << " grads" << endl;

    }

  MapTH1 res_map;
  int i_dphcut[dphcuts] = {12};
  TString ss_dphcut[dphcuts] = {"12"};

  TString Label[comparisons] = {"straightTracksY_isoAandCandB_straightTracksX"};
  TString Hist = "dx3_clchargeABC90evR"; 
  bool dphcut = true;

  GetHists(&res_map,Angles, dphcuts, comparisons,dphcut, Run, i_dphcut, Label, Hist, h_res);

 
  if(print) cout << "Fit:"  << endl;
  
  ofstream myfile;
  myfile.open ("/home/zoiirene/Output/TextFiles/Ascan_rms_"+detectorB+"_dphcutB"+ss_dphcut[0]+"_"+function+".txt");
  myfile << "A " << detectorA << "\n";
  myfile << labelA << "\n";
  myfile << "B " << detectorB<< "\n";
  myfile << labelB << "\n";
  myfile << "C " << detectorC<< "\n";
  myfile << labelC << "\n";
  myfile << info << "\n";
  myfile << "Angle RMS95(um) Error\n";

  for(int i=0; i<Angles; i++)
    {
      auto it2 = res_map.find(std::make_pair(Run[i]+"_"+ss_dphcut[0]+"_"+Label[0],Hist));
      if(it2  != res_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)                   cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;

	  
	  FitTH1(it2->second, &(Resolution[i]), &(ResolutionError[i]), ss_Angle[i], detectorA, detectorB, detectorC, function );
	  if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;
	  myfile << " " << Angle[i] << " " << Resolution[i] << " " << ResolutionError[i] << "\n";
	}
    }

 

  myfile.close();

  if(print) cout << "Plotting Res vs angle"  << endl;

  DrawTGraphWithErrorDouble(Angles, Angle, Resolution, ResolutionError, Hist,"2743", Label[0], "angleScan_dphcutB"+ss_dphcut[0]+"_rms", "RMS 95% [#mum]",-7.,29.,0.,8., "angleScan", "Angle [degrees]"  );

  ofstream myfile2;
  myfile2.open ("/home/zoiirene/Output/TextFiles/Ascan_res_"+detectorB+"_dphcutB"+ss_dphcut[0]+"_"+function+".txt");
  myfile2 << "A " << detectorA << "\n";
  myfile2 << labelA << "\n";
  myfile2 << "B " << detectorB<< "\n";
  myfile2 << labelB << "\n";
  myfile2 << "C " << detectorC<< "\n";
  myfile2 << labelC << "\n";
  myfile2 << info << "\n";
  myfile2 << "Angle Resolution(um) Error\n";

  for(int i=0; i<Angles; i++)
    {

      if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;
      ExtractRes(&(Resolution[i]),&(ResolutionError[i]));
      if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;
      myfile2 << " " << Angle[i] << " " << Resolution[i] << " " << ResolutionError[i] << "\n";
    }
  myfile2.close();
    

  DrawTGraphWithErrorDouble(Angles, Angle, Resolution, ResolutionError, Hist,"2743", Label[0], "angleScan_dphcutB"+ss_dphcut[0]+"_res", "Resolution [#mum]",-7.,29.,0.,8., "angleScan", "Angle [degrees]"  );
  




  Hist = "nrowB";
  GetHists(&res_map,Angles, dphcuts, comparisons,dphcut, Run, i_dphcut, Label, Hist, h_res);

  ofstream myfile3;
  myfile3.open ("/home/zoiirene/Output/TextFiles/Ascan_clsizeB_"+detectorB+"_dphcutB"+ss_dphcut[0]+"_"+function+".txt");
  myfile3 << "A " << detectorA << "\n";
  myfile3 << labelA << "\n";
  myfile3 << "B " << detectorB<< "\n";
  myfile3 << labelB << "\n";
  myfile3 << "C " << detectorC<< "\n";
  myfile3 << labelC << "\n";
  myfile3 << info << "\n";
  myfile3 << "Angle MeanNrowB Error\n";

  Double_t NrowB[Angles];
  Double_t NrowBerror[Angles];
  for(int i=0; i<Angles; i++)
    {
      auto it2 = res_map.find(std::make_pair(Run[i]+"_"+ss_dphcut[0]+"_"+Label[0],Hist));
      if(it2  != res_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;

	  

	  NrowB[i] = it2->second->GetMean();
	  NrowBerror[i] = it2->second->GetMeanError();

	  if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> nrowB: " << NrowB[i] << " and err: " << NrowBerror[i] << endl;
	  myfile3 << " " << Angle[i] << " " << NrowB[i]<< " " << NrowBerror[i] << "\n";
	}
      
    }
  myfile3.close();
  
  DrawTGraphWithErrorDouble(Angles, TanAngle, NrowB, NrowBerror, Hist,"2743", Label[0], "angleScan_dphcutB"+ss_dphcut[0]+"_nrowB_vs_tanangle", "nrowB",-7.,29.,0.,5., "angleScan", "Tan angle"  );
  DrawTGraphWithErrorDouble(Angles, Angle, NrowB, NrowBerror, Hist,"2743", Label[0], "angleScan_dphcutB"+ss_dphcut[0]+"_nrowB_vs_angle", "nrowB",-7.,29.,0.,5., "angleScan", "Angle [degrees]"  );
  
  Hist = "clsizeB";

  GetHists(&res_map,Angles, dphcuts, comparisons,dphcut, Run, i_dphcut, Label, Hist, h_res);



  for(int i=0; i<Angles; i++)
    {
      auto it2 = res_map.find(std::make_pair(Run[i]+"_"+ss_dphcut[0]+"_"+Label[0],Hist));
      if(it2  != res_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;

	  DrawHist(it2->second, Run[i]+"_"+ss_Angle[i], ss_dphcut[0], Label[0], Hist,0.,10.,0., 50000.)  ;

	}

    }
  



}//resolution 

