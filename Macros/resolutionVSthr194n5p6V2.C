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

#define BiasVoltages 1
#define Angles 1
bool print=true;
using namespace std;
#define dphcuts  9
#define comparisons 1


void resolutionVSthr194n5p6V2(TString func = "RMSself")
{

  TString detectorA="148";
  TString detectorB="194i";
  TString detectorC="150";
  TString labelA="FDB150P_12_R4S100x25-P4_1, thr 12 ADC";
  TString labelB="FDB150P_12_R4S100x25-P1_1, thr 18 ADC";
  TString labelC="FDB150P_12_R4S100x25-P1_3, thr 12 ADC";
  TString info  ="beam momentum 5.6 GeV";  
  
  
  double BiasVoltage[BiasVoltages];
  BiasVoltage[0]=800;//V  
  double BiasVoltageError[BiasVoltages];
  BiasVoltageError[0]=0;//V   
  double Angle[Angles]={12}; 
  TString ss_Angle[Angles] = {"12"}; ;//={10};//,12};//,14,16,18,20};
  
  Double_t Resolution[Angles][dphcuts];
  Double_t Percentage[Angles][dphcuts];
  Double_t ResolutionError[Angles][dphcuts];
  TString ss_BiasVoltage[BiasVoltages] = {"800"};
  
  TString inputDir="/home/zoiirene/Output/TextFiles/";
  TString outputDir="/home/zoiirene/Output/Plots/";
  TString inputfile;

  TH1F * h_res;

  // Int_t run[BiasVoltages]={3763,3764,3765,3766,3767,3768,3769,3770};
  Int_t run = 3839 ;
  TString Run[Angles][BiasVoltages];
  Run[0][0] = "3839";
  
  MapTH1 res_map[Angles];
  std::map<std::pair<TString, TString>, TH1F *>::iterator it;
  TString pitch = "25";
  int i_dphcut[dphcuts] = { 6 ,  10 ,  15 ,  19 ,  24 ,  29 ,  32 ,  40 ,  48 }; //{15};
  Double_t d_dphcut[dphcuts] = {15};
  Double_t d_dphcuterr[dphcuts] = {0};
  TString ss_dphcut[dphcuts] = { "6", "10", "15", "19", "24", "29", "32", "40", "48"};;//{"15"};

  TString Label[comparisons] = {"straightTracksY_isoAandCandB_straightTracksX"};
  TString Hist = "dx3_clphABC90evR";
  bool dphcut = true;
  TString label = "thrScan_A12C13"; 
  TH1F * hdx3_clchargeABC90evR[Angles];
  for(int j = 0; j < Angles; j++)
    {
      res_map[j] =  GetCheckHists(&(res_map[j]), Angles, dphcuts,dphcut, Run[j],ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR[j],true,label);
      //      GetHists(&(res_map[j]),BiasVoltages, dphcuts, comparisons,dphcut, Run[j], i_dphcut, Label, Hist, h_res);
    }
  

  


  if(print) cout << "*************    RMS: *************"  << endl;
  
  for (int l = 0;l<dphcuts;l++){
    ofstream myfile;
    myfile.open ("/home/zoiirene/Output/TextFiles/Ascan_"+detectorB+"_5p6_"+func+"_"+ss_Angle[0]+"deg_rmsTree_"+ss_dphcut[l]+"_"+label+".txt");
    myfile << "A " << detectorA << "\n";
    myfile << labelA << "\n";
    myfile << "B " << detectorB<< "\n";
    myfile << labelB << "\n";
    myfile << "C " << detectorC<< "\n";
    myfile << labelC << "\n";
  
    myfile << "Angle(deg RMS(um) Error\n";
      
    for(int i=0; i<Angles; i++)	{
      auto it2 = res_map[0].find(std::make_pair(Run[0][0]+"_"+ss_dphcut[l]+"_"+pitch,Hist));
      if(it2  != res_map[0].end()){
	if(print)               cout << " found map " << endl;
	if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
	
	FitTH1(it2->second, &(Resolution[0][i]), &(ResolutionError[0][i]), ss_Angle[0]+"_"+ss_BiasVoltage[0]+"_"+ss_dphcut[l], detectorA, detectorB, detectorC,func,&(Percentage[0][i]) );
	if(print) cout << "thr " << i<< ": " << ss_dphcut[i] << " ADC -> Resolution: " << Resolution[0][i] << " and res err: " << ResolutionError[0][i] << " with percentage "<< Percentage[0][i] <<  endl;
	myfile << ss_Angle[i] << " " << Resolution[0][i] << " " << ResolutionError[0][i] << "\n";
      }
      else{
	cout << " not found " << Run[0][0]+"_"+ss_dphcut[l]+"_"+Label[0] << " hist " << Hist << endl;
      }
    }
    
    myfile.close();
    DrawTGraphWithErrorDouble(dphcuts, d_dphcut, Resolution[0], ResolutionError[0], Hist,Run[0][0], Label[0], "AScan_5p6_"+ss_Angle[0]+"_rmsTree", "RMS 95% [#mum]",0.,7.,0.,5., "AScan_"+ss_Angle[0], "dphcut [ADC]"  );

  
  
    //////////////          RESOLUTION ///////////////////////
    if(print) cout << "*************    RMS: *************"  << endl;
    double freshres = 3.89; // Ascan_resTree_148_dphcutB12_RMSself_closest_A13C14_bestnonirr.txt  29/06/2020 deg 12.5
    double freshres_err = 0.03;
    ofstream myfileUnf;
	myfileUnf.open ("/home/zoiirene/Output/TextFiles/AScanTreeCorr_"+detectorB+"_5p6_"+func+"_"+ss_BiasVoltage[0]+"_bestAngle_res_dphcut"+ss_dphcut[l]+"_"+label+".txt");
	myfileUnf << "A " << detectorA << "\n";
	myfileUnf << labelA << "\n";
	myfileUnf << "B " << detectorB<< "\n";
	myfileUnf << labelB << "\n";
	myfileUnf << "C " << detectorC<< "\n";
	myfileUnf << labelC << "\n";
	
	myfileUnf << "Ascan(deg) RES(um) Error\n";

	for(int i = 0; i < Angles; i++)
	  {
      
   	  
	    if(print) cout << "thr " << i<< ": " << ss_dphcut[i] << " ADC -> Resolution: " << Resolution[0][i] << " and res err: " << ResolutionError[0][i] << endl;
	    ExtractRes(&(Resolution[0][i]),&(ResolutionError[0][i]),true, freshres, freshres_err);
	    if(print) cout << "thr " << i<< ": " << ss_dphcut[i] << " ADC -> Resolution: " << Resolution[0][i] << " and res err: " << ResolutionError[0][i] << endl;
	    myfileUnf << ss_Angle[0] << " " << Resolution[0][i] << " " << ResolutionError[0][i] << "\n";
	  }
	DrawTGraphWithErrorDouble(dphcuts, d_dphcut, Resolution[0], ResolutionError[0], Hist,Run[0][0], Label[0], "AScan_5p6_"+ss_Angle[0]+"_resTreeCorr", "RMS 95% [#mum]",0.,7.,0.,5., "thrScan_"+ss_Angle[0], "dphcut [ADC]"  );
 
    
    myfileUnf.close();
  
  }//dphcuts
  

  ////// cluster size
  /*
  Hist = "hnrowB_ph95";
  TH1F * hnrowB_clphABC90evR[Angles]; // = new TH1I("dx3_clchargeABC90evR ", "triplet dx_clchargeABC90evR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 90% events in Landau (only high tail)
  for(int j = 0; j < Angles; j++)
    hnrowB_clphABC90evR[j] =   GetCheckHists(&(res_map[j]), BiasVoltages, dphcuts,dphcut, Run[j],ss_dphcut,pitch, Hist,hnrowB_clphABC90evR[j],true,""+label+"");
  */
  /*
  Hist = "nrowB";
  cout << " Hist "<< Hist << endl;  
   for(int j = 0; j < Angles; j++)    {
     GetHists(&(res_map[j]),BiasVoltages, dphcuts, comparisons,dphcut, Run[j], i_dphcut, Label, Hist, h_res,true,label);
  }
  
   ofstream myfile3; 
   myfile3.open ("/home/zoiirene/Output/TextFiles/ThrScan_clsizeB_"+detectorB+"_"+ss_BiasVoltage[0]+"_bestAngle_5p6_"+label+".txt");
   myfile3 << "A " << detectorA << "\n";
   myfile3 << labelA << "\n";
   myfile3 << "B " << detectorB<< "\n";
   myfile3 << labelB << "\n";
   myfile3 << "C " << detectorC<< "\n";
   myfile3 << labelC << "\n";
   myfile3 << info << "\n";
   myfile3 << "Angle MeanNrowB Error\n";

   double NrowB[Angles][dphcuts];
   double  NrowBerror[Angles][dphcuts];
   

   for(int i=0; i<dphcuts; i++)     {
     auto it2 = res_map[0].find(std::make_pair(Run[0][0]+"_"+ss_dphcut[i]+"_"+Label[0],Hist));
     //auto it2 = res_map[0].find(std::make_pair(Run[0][0]+"_"+ss_dphcut[i]+"_"+pitch[0],Hist));
     if(it2  != res_map[0].end())
       {
	 if(print)  cout << " found map " << endl;
	 if(print)  cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetMean() << endl;
	 NrowB[0][i] = it2->second->GetMean();
	 NrowBerror[0][i] = it2->second->GetMeanError();
	 myfile3 << " " << i_dphcut[0] << " " << NrowB[0][i]<< " " << NrowBerror[0][i] << "\n";
       }else{
       cout << "not found " << Run[0][0]+"_"+ss_dphcut[i]+"_"+Label[0] << " hist  " << Hist << endl;
     }
   }
   
  myfile3.close();
  */

  /*
  TCanvas *csza[Angles];/// = new TCanvas("csza", "csza", 1500, 900);

  double NrowB[Angles][BiasVoltages];
  double  NrowBerror[Angles][BiasVoltages];



  TLegend* leg2a = new TLegend(0.7,0.35,0.85,0.85);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.04);

  
  for(int j = 0; j < Angles; j++)
    {
      csza[j] = new TCanvas("csza", "csza", 1500, 900);  
      gPad->SetTicks(1,1);
      gROOT->SetStyle("Plain");
      gStyle->SetPadGridX(0);
      gStyle->SetPadGridY(0);
      gStyle->SetPalette(1);

      //gStyle->SetOptStat(1110);
      gStyle->SetOptTitle(0);
      
      for(int i=0; i<BiasVoltages; i++)
	{
	  gStyle->SetOptStat(0);
	  auto it2 = res_map[j].find(std::make_pair(Run[j][i]+"_"+Label[0],Hist));
	  if(it2  != res_map[j].end())
	    {
	      if(print)  cout << " found map " << endl;
	      if(print)  cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetMean() << endl;
	      NrowB[j][i] = it2->second->GetMean();
	      NrowBerror[j][i] = it2->second->GetMeanError();
	      myfile3 << " " << Angle[j] << " " << NrowB[j][i]<< " " << NrowBerror[j][i] << "\n";
	    }
	  color = i+1;
	  if(i+1 ==5) color = 95;
	  if(i+1 ==10) color = 29;

	  it2->second->SetLineColor(color);
	  it2->second->SetLineWidth(1);
	  it2->second->SetLineStyle(1);
	  //	  it2->second->SetMarkerColor(color);
	  //it2->second->SetMarkerSize(1.5);
	  //it2->second->SetMarkerStyle(21);
	  it2->second->GetXaxis()->SetRangeUser(0.,5.);
	  it2->second->SetMinimum(0.);
	  it2->second->SetMaximum(40000.);

	  //  it2->second->SetTitle("Option ACP example");
	  it2->second->GetXaxis()->SetTitle("cluster size B");
	  it2->second->GetYaxis()->SetTitle("Entries");
	  //it2->second->SetFillStyle(3003);
	  //it2->second->SetFillColor(kRed-8);
	  //      leg2->AddEntry(it2->second,ss_BiasVoltage[j], "lp");
	  if(j==0)          leg2a->AddEntry(it2->second,ss_BiasVoltage[i], "l");

	  it2->second->Draw("histsame");
	  //	  it2->second->Draw("CPEsame");

	  
	}
      leg2a->SetNColumns(2);
      leg2a->SetTextSize(0.03);
      leg2a->Draw("same");

      name = outputDir+"compare_angle"+ss_Angle[j]+"_bias_Scan_clustersize_"+detectorB+"_"+Label[0]+"_"+Hist;
      csza[j]->SaveAs(name+".eps");
      csza[j]->SaveAs(name+".pdf");
      csza[j]->SaveAs(name+".png");
      csza[j]->SaveAs(name+".root");
    }


  
  TCanvas *csz = new TCanvas("csz", "csz", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  //  gStyle->SetOptStat();
  //gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);


  TGraph * grB[BiasVoltages];
  for(int j=0; j<BiasVoltages; j++)
    {
       
      grB[j]  = new TGraph(Angles, Angle, NrowB[j]);
      color = j+1;
      if(j+1 ==5) color = 95;
      if(j+1 ==10) color = 29;

      grB[j]->SetLineColor(color);
      grB[j]->SetLineWidth(1);
      grB[j]->SetLineStyle(2);
      grB[j]->SetMarkerColor(color);
      grB[j]->SetMarkerSize(1.5);
      grB[j]->SetMarkerStyle(21);
      grB[j]->GetXaxis()->SetRangeUser(8.,22.);
      grB[j]->SetMinimum(0.);
      grB[j]->SetMaximum(4.);

      //  grB[j]->SetTitle("Option ACP example");
      grB[j]->GetXaxis()->SetTitle("Angle [deg]");
      grB[j]->GetYaxis()->SetTitle("Mean cluster size B");
      //grB[j]->SetFillStyle(3003);
      //grB[j]->SetFillColor(kRed-8);
      //      leg2->AddEntry(grB[j],ss_BiasVoltage[j], "lp");


      if(j==0)   grB[j]->Draw("ACPE");
      else grB[j]->Draw("CPEsame");
	
    }

  
  leg2->SetTextSize(0.03);
  leg2->Draw("same");
  
  name = outputDir+"compare_angle_bias_Scan_clustersize_"+detectorB+"_"+Label[0]+"_"+Hist;
  csz->SaveAs(name+".eps");
  csz->SaveAs(name+".pdf");
  csz->SaveAs(name+".png");
  csz->SaveAs(name+".root");






  */



  
}//resolution 

