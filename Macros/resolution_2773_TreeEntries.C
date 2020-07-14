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

#define Angles 1 //32
bool print=true;
#define dphcuts  1
#define comparisons 1
using namespace std;
TString outputDir = "/home/zoiirene/Output/Plots/";
  
void resolution_2773_TreeEntries(TString function = "RMSself")
{

  TString detectorA[Angles]={"148"};
  TString detectorB[Angles]={"146"};
  TString detectorC[Angles]={"163"};
  TString labelA[Angles]={"Pstop_RD53Apads_FDB, thr 11 ADC"};
  TString labelB[Angles]={"Pspray_default_FDB, thr 13 ADC"};
  TString labelC[Angles]={"Pspray_RD53Apads_FDB, thr 14 ADC"};
  TString info = "beam energy 5.6 GeV";
  
  
  //  Double_t Angle[Angles]={-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,9.5,10,10.5,11,11.5,12,13,14,15,16,17,18,19,20,21,22,23,24};
  Double_t Angle[Angles]={15};
  Double_t Percentage[Angles];
  Double_t AngleError[Angles];
  Double_t TanAngle[Angles];
  Double_t TanAngleError[Angles];
  Double_t Resolution[Angles];
  Double_t Ncol[Angles];
  Double_t NcolError[Angles];
  Double_t ResolutionError[Angles];
  TString ss_Angle[Angles];
  TString name = "closest_A11C14";
  TH1F * h_res;//[Angles];

  Int_t run[Angles]={2773}; 
  TString Run[Angles];
  ostringstream strs[Angles];

  float angleunit = 1.25;
  if(print) cout << "Getting files "  << endl;
  TString pitch = "25";


  for(int i=0; i<Angles; i++)
    {

      if(print) cout << "Run: "<< i <<" " << run[i] << endl;
      Run[i].Form("%d",run[i]);
      if(print) cout << "Run: " << Run[i] << endl;
      ss_Angle[i].Form("%d",i);
      //      Angle[i]=Angle[i]*angleunit; //grad
      TanAngle[i]=TMath::Tan(Angle[i]*TMath::Pi()/((Double_t)180.));
      cout << " Angle " <<  Angle[i] << " TanAngle " <<  TanAngle[i] << endl;
      AngleError[i]=1; //grad
      TanAngleError[i]=(1+TMath::Power(TMath::Tan(AngleError[i]*TMath::Pi()/((Double_t) 180.)),2))*AngleError[i]*TMath::Pi()/((Double_t) 180.); //grad

      if(print) cout << "Angle " << i<< ": " << Angle[i] << " grads" << endl;

    }

 
  MapTH1 res_map;
  std::map<std::pair<TString, TString>, TH1F *>::iterator it;
  int i_dphcut[dphcuts] = {13};
  TString ss_dphcut[dphcuts] = {"13"};

  TString Label[comparisons] = {"straightTracksY_isoAandCandB_straightTracksX"};
  TString Hist = "dx3_clchargeABC90evR"; 
  bool dphcut = true;

  //GetHists(&res_map,Angles, dphcuts, comparisons,dphcut, Run, i_dphcut, Label, Hist, h_res);
  TH1F * hdx3_clchargeABC90evR; // = new TH1I("dx3_clchargeABC90evR ", "triplet dx_clchargeABC90evR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 90% events in Landau (only high tail)
  res_map =    GetCheckHists(&res_map, Angles, dphcuts,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,name);
  
  for ( it = res_map.begin(); it != res_map.end(); it++ )	{
    cout << it->first.first << " " << it->first.second << " " <<  it->second->GetEntries() << endl;
  }
  
  
 
  if(print) cout << "Fit:"  << endl;

  ofstream myfile[Angles];
  for(int i=0; i<Angles; i++){

      myfile[i].open ("/home/zoiirene/Output/TextFiles/Ascan_rmsTree_"+detectorB[i]+"_dphcutB"+ss_dphcut[i]+"_"+function+"_"+name+".txt");
      myfile[i] << "A " << detectorA[i] << "\n";
      myfile[i] << labelA[i] << "\n";
      myfile[i] << "B " << detectorB[i]<< "\n";
      myfile[i] << labelB[i] << "\n";
      myfile[i] << "C " << detectorC[i]<< "\n";
      myfile[i] << labelC[i] << "\n";
      myfile[i] << info << " " << Label[0] << "\n";
      myfile[i] << "Angle RMS95(um) Error\n";

      auto it2 = res_map.find(std::make_pair(Run[i]+"_"+ss_dphcut[i]+"_"+pitch,Hist));
      if(it2  != res_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)                   cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;

	  
	  FitTH1(it2->second, &(Resolution[i]), &(ResolutionError[i]), ss_Angle[i], detectorA[i], detectorB[i], detectorC[i], function, &(Percentage[i]) );
	  if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;
	  myfile[i] << " " << Angle[i] << " " << Resolution[i] << " " << ResolutionError[i] << "\n";
	}
  myfile[i].close();

  }

 


  if(print) cout << "Plotting Res vs angle"  << endl;

//  DrawTGraphWithErrorDouble(Angles, Angle, Resolution, ResolutionError, Hist,"2743", Label[0], "angleScan_dphcutB"+ss_dphcut[0]+"_rmsTree", "RMS 95% [#mum]",-7.,29.,0.,8., "angleScan", "Angle [degrees]"  );


  ofstream myfile2[Angles];
  for(int i=0; i<Angles; i++){
  
  myfile2[i].open ("/home/zoiirene/Output/TextFiles/Ascan_resTree_"+detectorB[i]+"_dphcutB"+ss_dphcut[i]+"_"+function+"_"+name+".txt");
  myfile2[i] << "A " << detectorA[i] << "\n";
  myfile2[i] << labelA[i] << "\n";
  myfile2[i] << "B " << detectorB[i]<< "\n";
  myfile2[i] << labelB[i] << "\n";
  myfile2[i] << "C " << detectorC[i]<< "\n";
  myfile2[i] << labelC[i] << "\n";
  myfile2[i] << info << " " << Label[0] << "\n";
  myfile2[i] << "Angle Resolution(um) Error\n";

 
      if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;
      ExtractRes(&(Resolution[i]),&(ResolutionError[i]));
      if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;
      myfile2[i] << " " << Angle[i] << " " << Resolution[i] << " " << ResolutionError[i] << "\n";
  myfile2[i].close();
    }
    

//DrawTGraphWithErrorDouble(Angles, Angle, Resolution, ResolutionError, Hist,"2743", Label[0], "angleScan_dphcutB"+ss_dphcut[0]+"_resTree", "Resolution [#mum]",-7.,29.,0.,8., "angleScan", "Angle [degrees]"  );
  


/*

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
*/



}//resolution 

