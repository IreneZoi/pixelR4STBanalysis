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
#define comparisons 4
#define dphcuts 1
using namespace std;
TString outputDir = "/home/zoiirene/Output/Plots/";
  
void resolutionVSangle2731_LandauRMStest()
{

  TString detectorA="146";
  TString detectorB="148";
  TString detectorC="163";
  TString labelA="Pspray_default_FDB, thr 12 ADC";
  TString labelB="Pstop_RD53Apads_FDB, 120V, thr 12 ADC";
  TString labelC="Pspray_RD53Apads_FDB, thr 12 ADC";
  TString info = "beam energy 5.6 GeV";
  

  TString inputDir="/home/zoiirene/Output/";
  TString label = "testABC";  
  //  Double_t Angle[Angles]={-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,9.5,10,10.5,11,11.5,12,13,14,15,16,17,18,19,20,21,22,23,24};
  Double_t Angle[Angles]={-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
  Double_t AngleError[Angles];
  Double_t TanAngle[Angles];
  Double_t TanAngleError[Angles];
  Double_t Resolution[Angles][comparisons];
  Double_t Ncol[Angles];
  Double_t NcolError[Angles];
  Double_t ResolutionError[Angles][comparisons];
  TString ss_Angle[Angles];
  
  TH1F * h_res;//[Angles];

  Int_t run[Angles]={2731,2732,2733,2734,2735,2736,2737,2738,2739,2740,2741,2743,2744,2745,2746,2747,2748,2749,2750,2751,2752,2753,2754,2755,2756,2757,2758,2759,2760};//,2764,2766,2767};
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
      Angle[i]=Angle[i]*angleunit; //grad
      TanAngle[i]=TMath::Tan(Angle[i]*TMath::Pi()/((Double_t)180.));
      cout << " Angle " <<  Angle[i] << " TanAngle " <<  TanAngle[i] << endl;
      AngleError[i]=1; //grad
      TanAngleError[i]=(1+TMath::Power(TMath::Tan(AngleError[i]*TMath::Pi()/((Double_t) 180.)),2))*AngleError[i]*TMath::Pi()/((Double_t) 180.); //grad

      if(print) cout << "Angle " << i<< ": " << Angle[i] << " grads" << endl;

    }

 
  MapTH1 res_map;
  std::map<std::pair<TString, TString>, TH1F *>::iterator it;
  //int i_dphcut[dphcuts] = {10, 12,  15, 17, 20, 22, 25, 27, 30, 32, 35, 37, 40};
  //TString ss_dphcut[dphcuts] = {"10", "12", "15",  "17",  "20", "22", "25", "27", "30", "32", "35", "37", "40"};
  int i_dphcut[dphcuts] = {22}; 
  TString ss_dphcut[dphcuts] = {"22"};

  TString Label[1] = {"straightTracksY_isoAandCandB_straightTracksX"};

  //TString Hist = "dx3_clchargeABC90evR";



  //LANDAU 90 - RMS
  TString Hist = "dx3_clphABC90evR";
  TString function = "RMS95";
  //  if(function=="RMS95") Hist = "dx3_clphABC90evR95";


  bool dphcut = true;
  cout << " hdx3_clchargeABC90evR"<< endl;
  //GetHists(&res_map,Angles, dphcuts, comparisons,dphcut, Run, i_dphcut, Label, Hist, h_res);
  TH1F * hdx3_clchargeABC90evR; // = new TH1I("dx3_clchargeABC90evR ", "triplet dx_clchargeABC90evR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 90% events in Landau (only high tail)
  //hdx3_clchargeABC90evR =
  res_map = GetCheckHists(&res_map, Angles, dphcuts,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,"testABC");

  cout << " hdx3_clchargeABC90evR done"<< endl;

  for ( it = res_map.begin(); it != res_map.end(); it++ )	{
    cout << it->first.first << " " << it->first.second << " " <<  it->second->GetEntries() << endl;
    cout << it->first.first << " " << it->first.second << " " <<  it->second->GetMean() << endl;
    cout << it->first.first << " " << it->first.second << " " <<  it->second->GetRMS() << endl;

  }
  
  cout << " hdx3_clchargeABC90evR printed"<< endl;
  /*
 
  if(print) cout << "Fit:"  << endl;
  ofstream myfile;
  myfile.open ("/home/zoiirene/Output/TextFiles/Ascan_rmsTree_"+detectorB+"_dphcutB"+ss_dphcut[0]+"_"+function+"_testThrABC_landau90RMS.txt");
  myfile << "A " << detectorA << "\n";
  myfile << labelA << "\n";
  myfile << "B " << detectorB<< "\n";
  myfile << labelB << "\n";
  myfile << "C " << detectorC<< "\n";
  myfile << labelC << "\n";
  myfile << info << " " << Label[0] << "\n";
  myfile << "Angle RMS95(um) Error\n";
  
  for(int i=0; i<Angles; i++){
    auto it2 = res_map.find(std::make_pair(Run[i]+"_"+ss_dphcut[0]+"_"+pitch,Hist));
    if(it2  != res_map.end()){
      if(print)               cout << " found map " << endl;
      if(print)                   cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
      if(print)                   cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetMean() << endl;
      if(print)                   cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetRMS() << endl;

      
      
      FitTH1(it2->second, &(Resolution[i][0]), &(ResolutionError[i][0]), ss_Angle[i], detectorA, detectorB, detectorC, function );
      if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i][0] << " and res err: " << ResolutionError[i][0] << endl;
      myfile << " " << Angle[i] << " " << Resolution[i][0] << " " << ResolutionError[i][0] << "\n";
    }
  }
  
  myfile.close();
    
  ofstream myfile2;
  myfile2.open ("/home/zoiirene/Output/TextFiles/Ascan_resTree_"+detectorB+"_dphcutB"+ss_dphcut[0]+"_"+function+"_testThrABC_landau90RMS.txt");
  myfile2 << "A " << detectorA << "\n";
  myfile2 << labelA << "\n";
  myfile2 << "B " << detectorB<< "\n";
  myfile2 << labelB << "\n";
  myfile2 << "C " << detectorC<< "\n";
  myfile2 << labelC << "\n";
  myfile2 << info << " " << Label[0] << "\n";
  myfile2 << "Angle Resolution(um) Error\n";
    
  for(int i=0; i<Angles; i++)      {
	
    if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i][0] << " and res err: " << ResolutionError[i][0] << endl;
    ExtractRes(&(Resolution[i][0]),&(ResolutionError[i][0]));
    if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i][0] << " and res err: " << ResolutionError[i][0] << endl;
    myfile2 << " " << Angle[i] << " " << Resolution[i][0] << " " << ResolutionError[i][0] << "\n";
  }
  myfile2.close();
    


  */


  //LANDAU 90 - RMS95
  Hist = "dx3_clphABC90evR95";
  function = "RMS95";


  dphcut = true;
  cout << " hdx3_clchargeABC90evR95"<< endl;
  res_map = GetCheckHists(&res_map, Angles, dphcuts,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,"testABC");

  cout << " hdx3_clchargeABC90evR95 done"<< endl;

  for ( it = res_map.begin(); it != res_map.end(); it++ )	{
    cout << it->first.first << " " << it->first.second << " " <<  it->second->GetEntries() << endl;
    cout << it->first.first << " " << it->first.second << " " <<  it->second->GetMean() << endl;
    cout << it->first.first << " " << it->first.second << " " <<  it->second->GetRMS() << endl;

  }
  
  cout << " hdx3_clchargeABC90evR5 printed"<< endl;
  
  /*
  if(print) cout << "Fit:"  << endl;
  ofstream secondfile;
  secondfile.open ("/home/zoiirene/Output/TextFiles/Ascan_rmsTree_"+detectorB+"_dphcutB"+ss_dphcut[0]+"_"+function+"_testThrABC_landau90RMS95.txt");
  secondfile << "A " << detectorA << "\n";
  secondfile << labelA << "\n";
  secondfile << "B " << detectorB<< "\n";
  secondfile << labelB << "\n";
  secondfile << "C " << detectorC<< "\n";
  secondfile << labelC << "\n";
  secondfile << info << " " << Label[0] << "\n";
  secondfile << "Angle RMS95(um) Error\n";
  
  for(int i=0; i<Angles; i++){
    auto it2 = res_map.find(std::make_pair(Run[i]+"_"+ss_dphcut[0]+"_"+pitch,Hist));
    if(it2  != res_map.end()){
      if(print)               cout << " found map " << endl;
      if(print)                   cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
      if(print)                   cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetMean() << endl;
      if(print)                   cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetRMS() << endl;

      
      
      FitTH1(it2->second, &(Resolution[i][1]), &(ResolutionError[i][1]), ss_Angle[i], detectorA, detectorB, detectorC, function );
      if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i][1] << " and res err: " << ResolutionError[i][1] << endl;
      secondfile << " " << Angle[i] << " " << Resolution[i][1] << " " << ResolutionError[i][1] << "\n";
    }
  }
  
  secondfile.close();
    
  ofstream secondfile2;
  secondfile2.open ("/home/zoiirene/Output/TextFiles/Ascan_resTree_"+detectorB+"_dphcutB"+ss_dphcut[0]+"_"+function+"_testThrABC_landau90RMS95.txt");
  secondfile2 << "A " << detectorA << "\n";
  secondfile2 << labelA << "\n";
  secondfile2 << "B " << detectorB<< "\n";
  secondfile2 << labelB << "\n";
  secondfile2 << "C " << detectorC<< "\n";
  secondfile2 << labelC << "\n";
  secondfile2 << info << " " << Label[0] << "\n";
  secondfile2 << "Angle Resolution(um) Error\n";
    
  for(int i=0; i<Angles; i++)      {
	
    if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i][1] << " and res err: " << ResolutionError[i][1] << endl;
    ExtractRes(&(Resolution[i][1]),&(ResolutionError[i][1]));
    if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i][1] << " and res err: " << ResolutionError[i][1] << endl;
    secondfile2 << " " << Angle[i] << " " << Resolution[i][1] << " " << ResolutionError[i][1] << "\n";
  }
  secondfile2.close();
  */




  //Landau all - RMS  
  cout << "****************   Landau all - RMS   *********************" << endl;
  TH1I * hdx3tree[Angles];

  for(int i=0; i<Angles; i++)      {
    TString inputfile = inputDir+"drei-r"+Run[i]+"_irene_dphcutB"+ss_dphcut[0]+"_"+label+".root";
    cout << inputfile << endl;
    TFile * file = new TFile(inputfile);
    TTree * tree = (TTree*)file->Get("charge_res");
    tree->Print();

    TH1I *    hdx3treeAll = new TH1I("hdx3treeAll", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
    tree->Draw("dx3tree>>hdx3treeAll","","goff");
    hdx3tree[i] = (TH1I*)gDirectory->Get("hdx3treeAll");
  }

    
 
  if(print) cout << "Fit:"  << endl;
  ofstream thirdfile;
  thirdfile.open ("/home/zoiirene/Output/TextFiles/Ascan_rmsTree_"+detectorB+"_dphcutB"+ss_dphcut[0]+"_"+function+"_testThrABC_landauAllRMS.txt");
  thirdfile << "A " << detectorA << "\n";
  thirdfile << labelA << "\n";
  thirdfile << "B " << detectorB<< "\n";
  thirdfile << labelB << "\n";
  thirdfile << "C " << detectorC<< "\n";
  thirdfile << labelC << "\n";
  thirdfile << info << " " << Label[0] << "\n";
  thirdfile << "Angle RMS95(um) Error\n";
  
  for(int i=0; i<Angles; i++){

    if(print)                   cout << hdx3tree[i]->GetEntries() << endl;
    if(print)                   cout << hdx3tree[i]->GetMean() << endl;
    if(print)                   cout << hdx3tree[i]->GetRMS() << endl;

      
      
      FitTH1(hdx3tree[i], &(Resolution[i][2]), &(ResolutionError[i][2]), ss_Angle[i], detectorA, detectorB, detectorC, function );
      if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i][2] << " and res err: " << ResolutionError[i][2] << endl;
      thirdfile << " " << Angle[i] << " " << Resolution[i][2] << " " << ResolutionError[i][2] << "\n";
    }
 
  
  thirdfile.close();
    
  ofstream thirdfile2;
  thirdfile2.open ("/home/zoiirene/Output/TextFiles/Ascan_resTree_"+detectorB+"_dphcutB"+ss_dphcut[0]+"_"+function+"_testThrABC_landauAllRMS.txt");
  thirdfile2 << "A " << detectorA << "\n";
  thirdfile2 << labelA << "\n";
  thirdfile2 << "B " << detectorB<< "\n";
  thirdfile2 << labelB << "\n";
  thirdfile2 << "C " << detectorC<< "\n";
  thirdfile2 << labelC << "\n";
  thirdfile2 << info << " " << Label[0] << "\n";
  thirdfile2 << "Angle Resolution(um) Error\n";
    
  for(int i=0; i<Angles; i++)      {
	
    if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i][2] << " and res err: " << ResolutionError[i][2] << endl;
    ExtractRes(&(Resolution[i][2]),&(ResolutionError[i][2]));
    if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i][2] << " and res err: " << ResolutionError[i][2] << endl;
    thirdfile2 << " " << Angle[i] << " " << Resolution[i][2] << " " << ResolutionError[i][2] << "\n";
  }
  thirdfile2.close();




  //Landau all - RMS95  
  cout << "Landau all - RM95S   " << endl;
  for(int i=0; i<Angles; i++)      {
    double tolerance = (hdx3tree[i]->GetBinContent(0)+hdx3tree[i]->GetBinContent(hdx3tree[i]->GetNbinsX()+1))/hdx3tree[i]->GetEntries(); // since we are not working with continuous quantities but with binned hists, the difference between the integ ral   and the 95% it may not be zero, so we ask it to be lower than 15% (before I was using 1.1% but since I changed to use the full range of hists, the overflow bin plays a too bigger role when stats is low and the 1.1% is not reached.

    cout << " RMS method charge based with tolerance " << tolerance << endl;
    double maximum = hdx3tree[i]->GetMaximum();
    int maxbin = hdx3tree[i]->GetMaximumBin();
    cout <<  " maxbin " << maxbin << " at " << hdx3tree[i]->GetBinCenter(maxbin)<<endl;
    cout << "intital sigma = " << hdx3tree[i]->GetRMS() * 1000 << " sigmaerr = " << hdx3tree[i]->GetRMSError() * 1000 << endl;
    double Integral = hdx3tree[i]->Integral(0,hdx3tree[i]->GetNbinsX()+1);
    cout << " integral " << Integral << " entries " << hdx3tree[i]->GetEntries() << endl;
    double integral95 = 0.95*Integral;
    cout << " integral95 " << integral95 << endl;
    double low = 0;
    double high = 0;
    for(int l =0; i<hdx3tree[i]->GetNbinsX()/2; l++)  {
      low = maxbin-l;
      high = maxbin+l;
      
      Integral = hdx3tree[i]->Integral(low,high);
      cout << " integral " << Integral << " low " << low << " high " << high << endl;
      cout << " while "<< l << " fabs(integral-integral95)/integral95 " << fabs(Integral-integral95)/integral95 << endl;
      
      //      if(fabs(Integral-integral95)/integral95 < 0.011 || Integral>integral95)//integral>integral95)
      if(fabs(Integral-integral95)/integral95 < tolerance || Integral>integral95)//integral>integral95)
	break;
      
      
    }
    cout << "final integral " << Integral << " low " << low <<" high " << high << endl;
    hdx3tree[i]->GetXaxis()->SetRange(low,high);
								       
  }

    
 
  if(print) cout << "Fit:"  << endl;
  ofstream fourthfile;
  fourthfile.open ("/home/zoiirene/Output/TextFiles/Ascan_rmsTree_"+detectorB+"_dphcutB"+ss_dphcut[0]+"_"+function+"_testThrABC_landauAllRMS95.txt");
  fourthfile << "A " << detectorA << "\n";
  fourthfile << labelA << "\n";
  fourthfile << "B " << detectorB<< "\n";
  fourthfile << labelB << "\n";
  fourthfile << "C " << detectorC<< "\n";
  fourthfile << labelC << "\n";
  fourthfile << info << " " << Label[0] << "\n";
  fourthfile << "Angle RMS95(um) Error\n";
  
  for(int i=0; i<Angles; i++){

    if(print)                   cout << hdx3tree[i]->GetEntries() << endl;
    if(print)                   cout << hdx3tree[i]->GetMean() << endl;
    if(print)                   cout << hdx3tree[i]->GetRMS() << endl;

      
      
      FitTH1(hdx3tree[i], &(Resolution[i][3]), &(ResolutionError[i][3]), ss_Angle[i], detectorA, detectorB, detectorC, function );
      if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i][3] << " and res err: " << ResolutionError[i][3] << endl;
      fourthfile << " " << Angle[i] << " " << Resolution[i][3] << " " << ResolutionError[i][3] << "\n";
    }
  
  
  fourthfile.close();
    
  ofstream fourthfile2;
  fourthfile2.open ("/home/zoiirene/Output/TextFiles/Ascan_resTree_"+detectorB+"_dphcutB"+ss_dphcut[0]+"_"+function+"_testThrABC_landauAllRMS95.txt");
  fourthfile2 << "A " << detectorA << "\n";
  fourthfile2 << labelA << "\n";
  fourthfile2 << "B " << detectorB<< "\n";
  fourthfile2 << labelB << "\n";
  fourthfile2 << "C " << detectorC<< "\n";
  fourthfile2 << labelC << "\n";
  fourthfile2 << info << " " << Label[0] << "\n";
  fourthfile2 << "Angle Resolution(um) Error\n";
    
  for(int i=0; i<Angles; i++)      {
	
    if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i][3] << " and res err: " << ResolutionError[i][3] << endl;
    ExtractRes(&(Resolution[i][3]),&(ResolutionError[i][3]));
    if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i][3] << " and res err: " << ResolutionError[i][3] << endl;
    fourthfile2 << " " << Angle[i] << " " << Resolution[i][3] << " " << ResolutionError[i][3] << "\n";
  }
  fourthfile2.close();




  
}//resolution 

