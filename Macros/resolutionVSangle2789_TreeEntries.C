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

#define Angles 28
bool print=true;
#define dphcuts 9
#define comparisons 1

using namespace std;
TString outputDir = "/home/zoiirene/Output/Plots/testEntries/";
TString inputDir = "/home/zoiirene/Output/TextFiles/";
  
void resolutionVSangle2789_TreeEntries(TString function = "RMSself"){

  TString detectorA="148";
  TString detectorB="120i";
  TString detectorC="163";
  TString labelA="Pstop_RD53Apads_FDB, thr 12 ADC";
  TString labelB="Pstop_default_FTH, 800V, thr 15 ADC";
  TString labelC="Pspray_RD53Apads_FDB, thr 12 ADC";
  TString info = "beam energy 5.6 GeV";
  TString testname="beamdiv_A12C15";//"closest_A12C15_bestproton";
  TString testnameUnf="RMSselfiso600_beamdiv_A13C14";//"RMSself6sig_closest_A13C14_bestnonirr";
  
  Double_t Angle[Angles];
  Double_t AngleError[Angles];
  Double_t TanAngle[Angles];
  Double_t TanAngleError[Angles];
  Double_t Resolution[Angles];
  Double_t RMS[Angles];
  Double_t Ncol[Angles];
  Double_t NcolError[Angles];
  Double_t ResolutionError[Angles];
  Double_t RMSError[Angles];
  Double_t Percentage[Angles];
  Double_t Max[Angles];
  Double_t Min[Angles];
  TString ss_Angle[Angles];
  
  TH1F * h_res;//[Angles];

  Int_t run=2789;
  TString Run[Angles];
  ostringstream strs[Angles];
  float angleunit = 1.25;  
  TString pitch = "25";
  

  if(print) cout << "Getting files "  << endl;
  for(int i=0; i<11; i++)
    {
	  Run[i].Form("%d",run+i);
	  if(print) cout << "Run: " << Run << endl;
    }
  for(int i=11; i<Angles; i++)
    {
	  Run[i].Form("%d",run+i+1);
	  if(print) cout << "Run: " << Run << endl;
    }



  for(int i=0; i<Angles; i++)
    {

      Angle[i] = -4+i;
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
  int i_dphcut[dphcuts] = {5, 8, 11, 15, 19,  23, 25, 31, 38}; //{15};
  TString ss_dphcut[dphcuts] ={"5", "8", "11", "15", "19", "23", "25", "31", "38"}; // {"15"};

  TString Label[comparisons] = {"straightTracksY_isoAandCandB_straightTracksX"};
  TString Hist = "dx3_clphABC90evR";
  bool dphcut = true;

  TH1F * hdx3_clchargeABC90evR; // = new TH1I("dx3_clchargeABC90evR ", "triplet dx_clchargeABC90evR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 90% events in Landau (only high tail)
  //hdx3_clchargeABC90evR
  res_map =    GetCheckHists(&res_map, Angles, dphcuts,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,testname);

  for ( it = res_map.begin(); it != res_map.end(); it++ )
    {

      cout << it->first.first << " " << it->first.second << " " <<  it->second->GetEntries() << endl;
    }
  //GetHists(&res_map,Angles, dphcuts, comparisons,dphcut, Run, i_dphcut, Label, Hist, h_res);



  for(int k = 0; k<dphcuts; k ++){  


  if(print) cout << "Fit:"  << endl;
  
  ofstream myfile;
  myfile.open ("/home/zoiirene/Output/TextFiles/AscanTree_"+detectorB+"_dphcutB"+ss_dphcut[k]+"_"+function+"_"+testname+".txt");
  myfile << "A " << detectorA << "\n";
  myfile << labelA << "\n";
  myfile << "B " << detectorB<< "\n";
  myfile << labelB << "\n";
  myfile << "C " << detectorC<< "\n";
  myfile << labelC << "\n";
  myfile << info << "\n";
  myfile << "Angle Resolution(um) Error\n";


  ofstream myentries;
  myentries.open ("/home/zoiirene/Output/TextFiles/Ascan_Entries_"+detectorB+"_dphcutB"+ss_dphcut[k]+"_"+function+"_"+testname+".txt");
  myentries << "A " << detectorA << "\n";
  myentries << labelA << "\n";
  myentries << "B " << detectorB<< "\n";
  myentries << labelB << "\n";
  myentries << "C " << detectorC<< "\n";
  myentries << labelC << "\n";
  myentries << info << "\n";
  myentries << "Angle Entries\n";

  ofstream myperc;
  myperc.open ("/home/zoiirene/Output/TextFiles/Ascan_Percentage_"+detectorB+"_dphcutB"+ss_dphcut[k]+"_"+function+"_"+testname+".txt");
  myperc << "A " << detectorA << "\n";
  myperc << labelA << "\n";
  myperc << "B " << detectorB<< "\n";
  myperc << labelB << "\n";
  myperc << "C " << detectorC<< "\n";
  myperc << labelC << "\n";
  myperc << info << "\n";
  myperc << "Angle Percentage\n";


  //  function="RMS95";
  for(int i=0; i<Angles; i++)
    {

      auto it2 = res_map.find(std::make_pair(Run[i]+"_"+ss_dphcut[k]+"_"+pitch,Hist));
      if(it2  != res_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)                   cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;

	  FitTH1(it2->second, &(RMS[i]), &(RMSError[i]), ss_Angle[i]+"_dphcutB"+ss_dphcut[k]+"_", detectorA, detectorB, detectorC, function, &(Percentage[i]),&(Min[i]),&(Max[i]));
	  if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << RMS[i] << " and res err: " << ResolutionError[i] << endl;
	  myfile << Angle[i] << " "  << RMS[i] << " " << RMSError[i] << "\n";
	  myentries << Angle[i] << " "  << it2->second->GetEntries() << "\n";
	  myperc << Angle[i] << " "  << Percentage[i] << "\n";
	}
    }

 

  myfile.close();
  myentries.close();
  myperc.close();
  if(print) cout << "Plotting Res vs angle"  << endl;



  DrawTGraphWithErrorDouble(Angles, Angle, RMS, RMSError, Hist,"2789", Label[0], "angleScan_dphcutB"+ss_dphcut[k]+"_rmsTree", "RMS 95% [#mum]",-7.,29.,0.,8., "angleScan", "Angle [degrees]"  );

  cout << "########     resolution unfolding       ###########" << endl;
  
  double freshres = 3.8;
  double freshres_err = 0.01;
  //if(print) cout << "scaling residual to get resolution " <<  conversion <<endl;
  //function="RMS";
  
  ofstream myfile2;
  myfile2.open ("/home/zoiirene/Output/TextFiles/Ascan_resTree_"+detectorB+"_dphcutB"+ss_dphcut[k]+"_"+function+"_"+testname+".txt");
  myfile2 << "A " << detectorA << "\n";
  myfile2 << labelA << "\n";
  myfile2 << "B " << detectorB<< "\n";
  myfile2 << labelB << "\n";
  myfile2 << "C " << detectorC<< "\n";
  myfile2 << labelC << "\n";
  myfile2 << info << "\n";
  myfile2 << "Angle RMS(um) Error\n";


  for(int i=0; i<Angles; i++)
    {

      Resolution[i] = RMS[i];
      ResolutionError[i] = RMSError[i];
      if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> RMS: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;

      ExtractRes(&(Resolution[i]), &(ResolutionError[i]), true,freshres,freshres_err);
      if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;
      myfile2 << " " << Angle[i] << " " << Resolution[i] << " " << ResolutionError[i] << "\n";
    }
  myfile2.close();


  DrawTGraphWithErrorDouble(Angles, Angle, Resolution, ResolutionError, Hist,"2789", Label[0], "angleScan_dphcutB"+ss_dphcut[k]+"_resTree", "Resolution [#mum]",-7.,29.,0.,8., "angleScan", "Angle [degrees]"  );


  cout << "######## Angle dependent    resolution unfolding       ###########" << endl;

 TString  filenames = "Ascan_resTree_148_dphcutB12_"+testnameUnf+".txt"; 
 filenames      = inputDir+filenames;
 cout << filenames << endl;
 ifstream stream(filenames);
 double angle;
 double res,reserror;
 double freshresA[Angles]; 
 double freshresA_err[Angles];
 
 std::string line;
 if(!stream.is_open())
   {
     cout << " File " << filenames << " not opened" << endl;
   }
 else
   {
     for(int k =0; k<8;k++)
       {
	 std::getline(stream,line);
	 if(print) cout << k << " line " << line << endl;
       }
     //      while(!stream.eof())
     for(int j= 0; j < Angles; j++)
       {
	 stream  >> angle >> res >> reserror ;
	 if(print)  cout << "line " << j+7 << endl;

	 freshresA[j] = res;
	 freshresA_err[j] = reserror;
	 if(print)      cout  << angle << " " << freshresA[j] << " " << freshresA_err[j] << endl;
       }//angles
   }// file scanning
	  

 
  
   ofstream myfile2c;
   TString filenameR = "/home/zoiirene/Output/TextFiles/Ascan_resTreeCorr_"+detectorB+"_dphcutB"+ss_dphcut[k]+"_"+function+"_"+testname+"_"+testnameUnf+".txt";
   cout << filenameR << endl;
   myfile2c.open (filenameR);
   myfile2c << "A " << detectorA << "\n";
   myfile2c << labelA << "\n";
   myfile2c << "B " << detectorB<< "\n";
   myfile2c << labelB << "\n";
   myfile2c << "C " << detectorC<< "\n";
   myfile2c << labelC << "\n";
   myfile2c << info << "\n";
   myfile2c << "Angle RMS(um) Error\n";
   

  for(int i=0; i<Angles; i++)
    {

      Resolution[i] = RMS[i];
      ResolutionError[i] = RMSError[i];
      if(print) cout << "before" <<endl;
      if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> RMS: " << RMS[i] << " and res err: " << RMSError[i] << endl;
      if(print) cout << "fresh" <<endl;
      if(print)      cout  << i << " " << freshresA[i] << " " << freshresA_err[i] << endl;
      ExtractRes(&(Resolution[i]), &(ResolutionError[i]), true,freshresA[i],freshresA_err[i]);
      if(print) cout << "after" <<endl;
      if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;
      myfile2c << " " << Angle[i] << " " << Resolution[i] << " " << ResolutionError[i] << "\n";
    }
  myfile2c.close();


  DrawTGraphWithErrorDouble(Angles, Angle, Resolution, ResolutionError, Hist,"2789", Label[0], "angleScan_dphcutB"+ss_dphcut[k]+"_resTreeCorr", "Resolution [#mum]",-7.,29.,0.,8., "angleScan", "Angle [degrees]"  );

  }//dphcuts

  cout << "making other plots " <<  endl;

  /*
  Hist = "hnrowB_ph95";
  TH1F * hnrowB_clphABC90evR; // = new TH1I("dx3_clchargeABC90evR ", "triplet dx_clchargeABC90evR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 90% events in Landau (only high tail)
  hnrowB_clphABC90evR =    GetCheckHists(&res_map, Angles, dphcuts,dphcut, Run,ss_dphcut,pitch, Hist,hnrowB_clphABC90evR,true,"testEntries");
  */


  Hist = "nrowB_clphABC90evR";
  res_map =    GetCheckHists(&res_map, Angles, dphcuts,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,testname);

  for (int k = 0; k<dphcuts; k++){
    ofstream myfile3;
    myfile3.open ("/home/zoiirene/Output/TextFiles/Ascan_clsizeB_"+detectorB+"_dphcutB"+ss_dphcut[k]+"_"+function+"_"+testname+"_ABC90.txt");
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
	auto it2 = res_map.find(std::make_pair(Run[i]+"_"+ss_dphcut[k]+"_"+pitch,Hist));
	//auto it2 = res_map.find(std::make_pair(Run[i]+"_"+ss_dphcut[k]+"_"+pitch,Hist));
	if(it2  != res_map.end())
	  {
	    if(print)               cout << " found map " << endl;
	    if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
	    cout << " Getting clusters " << endl;


	    NrowB[i] = it2->second->GetMean();
	    NrowBerror[i] = it2->second->GetMeanError();

	    if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> nrowB: " << NrowB[i] << " and err: " << NrowBerror[i] << endl;
	    myfile3 << " " << Angle[i] << " " << NrowB[i]<< " " << NrowBerror[i] << "\n";
	  }

      }
    myfile3.close();

    DrawTGraphWithErrorDouble(Angles, TanAngle, NrowB, NrowBerror, Hist,"2743", Label[0], "angleScan_dphcutB"+ss_dphcut[k]+"_nrowB_vs_tanangle_ABC90", "nrowB",-7.,29.,0.,5., "angleScan", "Tan angle"  );
    DrawTGraphWithErrorDouble(Angles, Angle, NrowB, NrowBerror, Hist,"2743", Label[0], "angleScan_dphcutB"+ss_dphcut[k]+"_nrowB_vs_angle_ABC90", "nrowB",-7.,29.,0.,5., "angleScan", "Angle [degrees]"  );
  }
  


  /////





  Hist = "nrowB";

  GetHists(&res_map,Angles, dphcuts, comparisons,dphcut, Run, i_dphcut, Label, Hist, h_res,true,testname);
  for(int k = 0; k<dphcuts; k ++){
  ofstream myfile4;
  myfile4.open ("/home/zoiirene/Output/TextFiles/Ascan_clsizeB_"+detectorB+"_dphcutB"+ss_dphcut[k]+"_"+"_"+function+"_"+testname+".txt");
  myfile4 << "A " << detectorA << "\n";
  myfile4 << labelA << "\n";
  myfile4 << "B " << detectorB<< "\n";
  myfile4 << labelB << "\n";
  myfile4 << "C " << detectorC<< "\n";
  myfile4 << labelC << "\n";
  myfile4 << info << "\n";
  myfile4 << "Angle MeanNrowB Error\n";
  
  Double_t NrowB[Angles];
  Double_t NrowBerror[Angles];
  for(int i=0; i<Angles; i++)
    {
      auto it2 = res_map.find(std::make_pair(Run[i]+"_"+ss_dphcut[k]+"_"+Label[0],Hist));
      //auto it2 = res_map.find(std::make_pair(Run[i]+"_"+ss_dphcut[k]+"_"+pitch[0],Hist));
      if(it2  != res_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;



	  NrowB[i] = it2->second->GetMean();
	  NrowBerror[i] = it2->second->GetMeanError();
	  myfile4 << " " << Angle[i] << " " << NrowB[i]<< " " << NrowBerror[i] << "\n";
	  if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> nrowB: " << NrowB[i] << " and err: " << NrowBerror[i] << endl;

	}

    }
  myfile4.close();
  
  DrawTGraphWithErrorDouble(Angles, TanAngle, NrowB, NrowBerror, Hist,"2789", Label[0], "angleScan_dphcutB"+ss_dphcut[k]+"_nrowB_vs_tanangle", "nrowB",-7.,29.,0.,5., "angleScan", "Tan angle"  );
  DrawTGraphWithErrorDouble(Angles, Angle, NrowB, NrowBerror, Hist,"2789", Label[0], "angleScan_dphcutB"+ss_dphcut[k]+"_nrowB_vs_angle", "nrowB",-7.,29.,0.,5., "angleScan", "Angle [degrees]"  );
  }//dphcuts
  /*
  Hist = "clsizeB";
  GetHists(&res_map,Angles, dphcuts, comparisons,dphcut, Run, i_dphcut, Label, Hist, h_res);

  for(int i=0; i<Angles; i++)
    {
      auto it2 = res_map.find(std::make_pair(Run[i]+"_"+ss_dphcut[k]+"_"+Label[0],Hist));
      if(it2  != res_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;

	  DrawHist(it2->second, Run[i]+"_"+ss_Angle[i], ss_dphcut[k], Label[0], Hist,0.,10.,0., 50000.)  ;

	}

    }
  */
  
  
  /*
  //*******************size plot
  TCanvas *c3 = new TCanvas("c3", "Mean Number of columns vs angle", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TGraph* NcolPlot = new TGraphErrors(Angles,Angle,Ncol,AngleError,NcolError);

  NcolPlot->SetTitle(" ");
  NcolPlot->GetYaxis()->SetTitle("Mean Number of columns");
  NcolPlot->GetXaxis()->SetTitle("Angle [degrees]");
  NcolPlot->SetMarkerSize(2.5);
  NcolPlot->SetLineColor(kBlue);
  NcolPlot->SetMarkerColor(kBlue);
  NcolPlot->SetMarkerStyle(20);
  NcolPlot->SetLineWidth(2);
  NcolPlot->SetLineStyle(2);
  NcolPlot->GetXaxis()->SetRangeUser(0.,29.);
  NcolPlot->GetYaxis()->SetRangeUser(0.,3.);
  NcolPlot->Draw("AEPL");

  TLegend* leg3 = new TLegend(0.65,0.4,0.75,0.5);
  leg3->SetLineColor(0);
  leg3->SetTextSize(0.03);
  leg3->AddEntry(NcolPlot,"DUT "+detectorB , "ep");
  leg3->Draw();

  
  c3->SaveAs(outputDir+"Ncol_vs_angle_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".eps");
  c3->SaveAs(outputDir+"Ncol_vs_angle_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".root");

  //*******************size plot vs tang

  for(int i=0; i<Angles; i++)
    cout << " Angle " <<  Angle[i] << " TanAngle " <<  TanAngle[i] << endl;
       

  TCanvas *c5 = new TCanvas("c5", "Number of columns vs tan angle", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TGraph* NcolTanPlot = new TGraphErrors(Angles,TanAngle,Ncol,TanAngleError,NcolError);

  NcolTanPlot->SetTitle(" ");
  NcolTanPlot->GetYaxis()->SetTitle("Mean Number of columns");
  NcolTanPlot->GetXaxis()->SetTitle("Tan(Angle)");
  NcolTanPlot->SetMarkerSize(2.5);
  NcolTanPlot->SetLineColor(kBlue);
  NcolTanPlot->SetMarkerColor(kBlue);
  NcolTanPlot->SetMarkerStyle(20);
  NcolTanPlot->SetLineWidth(2);
  NcolTanPlot->SetLineStyle(2);
  //  NcolTanPlot->GetXaxis()->SetRangeUser(0.017,0.018);
  NcolTanPlot->GetXaxis()->SetLimits(0.,0.6);
  NcolTanPlot->GetYaxis()->SetRangeUser(0.,3.);
  NcolTanPlot->Draw("AEPL");
  NcolTanPlot->Fit("pol1");
  NcolTanPlot->Draw("ALPEsame");

  TLegend* leg4 = new TLegend(0.65,0.4,0.75,0.5);
  leg4->SetLineColor(0);
  leg4->SetTextSize(0.03);
  leg4->AddEntry(NcolTanPlot,"DUT "+detectorB , "ep");
  leg4->Draw();

  
  c5->SaveAs(outputDir+"Ncol_vs_tanangle_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".eps");
  c5->SaveAs(outputDir+"Ncol_vs_tanangle_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".root");

  /*
  TCanvas *c = new TCanvas("c", "resolution vs inverse beam energy", 1500, 900);
  gPad->SetTicks(1,1);
  TGraph* resolutionPlotInv = new TGraphErrors(Angles,BeamEnergyInverse,Resolution,BeamEnergyError,ResolutionError);

  resolutionPlotInv->SetTitle(" ");
  resolutionPlotInv->GetYaxis()->SetTitle("Resolution [#mum]");
  resolutionPlotInv->GetXaxis()->SetTitle("1/Beam Energy [GeV^{-1}]");
  resolutionPlotInv->SetMarkerSize(2.5);
  resolutionPlotInv->SetLineColor(2);
  resolutionPlotInv->SetMarkerColor(2);
  resolutionPlotInv->SetMarkerStyle(20);
  resolutionPlotInv->SetLineWidth(2);
  resolutionPlotInv->SetLineStyle(2);
  resolutionPlotInv->GetYaxis()->SetRangeUser(0.,20.);
  resolutionPlotInv->GetXaxis()->SetLimits(0.,1.);

  //  TF1 *fit1 = new TF1("fit1","pol1", 0.15, 1.255); 
  TF1 *fit1 = new TF1("fit1","sqrt([0]*[0]+([1]*x)*([1]*x))", 0., 0.8); 
  fit1->SetLineColor(kBlue);
  fit1->SetParameter(0,3);
  fit1->SetParameter(1,15);
  fit1->SetParName(0,"#sigma_{hit}");
  fit1->SetParName(1,"#sigma_{MS}");
  
  resolutionPlotInv->Fit("fit1","R");
  resolutionPlotInv->Draw("AEP");

  // TPaveStats * ps = (TPaveStats *)fit1->FindObject("stats");
  // ps->SetX1NDC(0.2);
  // ps->SetX2NDC(0.2);
  gStyle->SetStatX(0.5);
  gStyle->SetStatY(0.85);
  
  TLegend* leg = new TLegend(0.65,0.4,0.75,0.6);
  leg->SetLineColor(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(resolutionPlotInv,"DUT "+detectorB , "ep");
  //  leg->AddEntry(fit1,"p_{0}+x*p_{1}" , "l");
  leg->AddEntry(fit1,"#sqrt{#sigma_{hit}^{2}+(#sigma_{MS}/p)^{2}}" , "l");
  leg->Draw();

  c->SaveAs("/home/irene/Documents/CMS_material/DESY_TB/Plot/Res_vs_InvEnergy_newFit_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".png");  
  c->SaveAs("/home/irene/Documents/CMS_material/DESY_TB/Plot/Res_vs_InvEnergy_newFit_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".eps");
  c->SaveAs("/home/irene/Documents/CMS_material/DESY_TB/Plot/Res_vs_InvEnergy_newFit_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".pdf");


  TCanvas *c4 = new TCanvas("c4", "resolution vs inverse beam energy squared", 1500, 900);
  gPad->SetTicks(1,1);

  TGraph* resolutionPlotInvGeorg = new TGraphErrors(Angles,BeamEnergyInverseSquare,ResolutionSquare,BeamEnergyError,ResolutionErrorSquare);

  resolutionPlotInvGeorg->SetTitle(" ");
  resolutionPlotInvGeorg->GetYaxis()->SetTitle("Resolution^{2} [#mum^{2}]");
  resolutionPlotInvGeorg->GetXaxis()->SetTitle("1/(Beam Energy)^{2} [GeV^{-2}]");
  resolutionPlotInvGeorg->SetMarkerSize(2.5);
  resolutionPlotInvGeorg->SetLineColor(2);
  resolutionPlotInvGeorg->SetMarkerColor(2);
  resolutionPlotInvGeorg->SetMarkerStyle(20);
  resolutionPlotInvGeorg->SetLineWidth(2);
  resolutionPlotInvGeorg->SetLineStyle(2);
  //  resolutionPlotInvGeorg->GetYaxis()->SetRangeUser(0.,40.);
  resolutionPlotInvGeorg->GetXaxis()->SetLimits(0.,1.);

  //  TF1 *fit1 = new TF1("fit1","pol1", 0.15, 1.255); 
  TF1 *fit3 = new TF1("fit3","pol1", 0., 0.8); 
  fit3->SetLineColor(kBlue);
  //  fit1->SetParameter(0,3);
  //fit1->SetParameter(1,15);
  fit3->SetParName(0,"#sigma_{hit}^{2}");
  fit3->SetParName(1,"#sigma_{MS}^{2}");
  
  resolutionPlotInvGeorg->Fit("fit3","R");
  resolutionPlotInvGeorg->Draw("AEP");

  // TPaveStats * ps = (TPaveStats *)fit1->FindObject("stats");
  // ps->SetX1NDC(0.2);
  // ps->SetX2NDC(0.2);
  gStyle->SetStatX(0.5);
  gStyle->SetStatY(0.85);

  ostringstream strfit[2];
  TString ss_fit[2];
  ostringstream strfit_err[2];
  TString ss_fit_err[2];
  for(int i=0; i<2;i++)
    {
      strfit[i] << setprecision(3) << sqrt(fit3->GetParameter(i));
      ss_fit[i]=strfit[i].str();
      strfit_err[i] << setprecision(1) << 0.5*fit3->GetParError(i)/sqrt(fit3->GetParameter(i));
      ss_fit_err[i]=strfit_err[i].str();
    }
  
  TLatex Tl_2;
  Tl_2.SetTextAlign(12);
  Tl_2.SetTextSize(0.04);
  Tl_2.DrawLatexNDC(0.15,0.6,"#sigma_{hit} = ("+ss_fit[0]+" #pm "+ss_fit_err[0]+") #mum");
  Tl_2.DrawLatexNDC(0.15,0.52,"#sigma_{MS} = ("+ss_fit[1]+" #pm "+ss_fit_err[1]+") #mum*GeV");
  
  
  TLegend* leg4 = new TLegend(0.75,0.4,0.85,0.6);
  leg4->SetLineColor(0);
  leg4->SetTextSize(0.03);
  leg4->AddEntry(resolutionPlotInvGeorg,"DUT "+detectorB , "ep");
  leg4->AddEntry(fit3,"#sigma_{hit}^{2}+(#sigma_{MS}/p)^{2}" , "l");
  leg4->Draw();

  
  c4->SaveAs("/home/irene/Documents/CMS_material/DESY_TB/Plot/Res_vs_InvEnergySquared_newFit_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".png");  
  c4->SaveAs("/home/irene/Documents/CMS_material/DESY_TB/Plot/Res_vs_InvEnergySquared_newFit_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".eps");
  c4->SaveAs("/home/irene/Documents/CMS_material/DESY_TB/Plot/Res_vs_InvEnergySquared_newFit_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".pdf");

  

  TCanvas *c3 = new TCanvas("c3", "resolution vs inverse beam energy old", 1500, 900);
  gPad->SetTicks(1,1);

  TGraph* resolutionPlotInvOLD = new TGraphErrors(Angles,BeamEnergyInverse,Resolution,BeamEnergyError,ResolutionError);

  resolutionPlotInvOLD->SetTitle(" ");
  resolutionPlotInvOLD->GetYaxis()->SetTitle("Resolution [#mum]");
  resolutionPlotInvOLD->GetXaxis()->SetTitle("1/Beam Energy [GeV^{-1}]");
  resolutionPlotInvOLD->SetMarkerSize(2.5);
  resolutionPlotInvOLD->SetLineColor(2);
  resolutionPlotInvOLD->SetMarkerColor(2);
  resolutionPlotInvOLD->SetMarkerStyle(20);
  resolutionPlotInvOLD->SetLineWidth(2);
  resolutionPlotInvOLD->SetLineStyle(2);
  resolutionPlotInvOLD->GetYaxis()->SetRangeUser(0.,20.);
  resolutionPlotInvOLD->GetXaxis()->SetLimits(0.,1.6);

  TF1 *fit2 = new TF1("fit2","pol1", 0.15, 1.255); 
  //  TF1 *fit1 = new TF1("fit1","sqrt([0]*[0]+([1]*x)*([1]*x))", 0., 0.8); 
  fit2->SetLineColor(kBlue);
  // fit1->SetParameter(0,3);
  // fit1->SetParameter(1,15);
  // fit1->SetParameterName(0,"#sigma_{hit}");
  // fit1->SetParameterName(1,"#sigma_{MS}");
  
  resolutionPlotInvOLD->Fit("fit2");
  resolutionPlotInvOLD->Draw("AEP");

  // TPaveStats * ps = (TPaveStats *)fit1->FindObject("stats");
  // ps->SetX1NDC(0.2);
  // ps->SetX2NDC(0.2);
  gStyle->SetStatX(0.5);
  gStyle->SetStatY(0.85);
  
  TLegend* leg3 = new TLegend(0.65,0.4,0.75,0.6);
  leg3->SetLineColor(0);
  leg3->SetTextSize(0.03);
  leg3->AddEntry(resolutionPlotInvOLD,"DUT "+detectorB , "ep");
  leg3->AddEntry(fit2,"p_{0}+x*p_{1}" , "l");
  //  leg->AddEntry(fit1,"#sqrt{#sigma_{hit}^{2}+(#sigma_{MS}/p)^{2}}" , "l");
  leg3->Draw();

  
  c3->SaveAs("/home/irene/Documents/CMS_material/DESY_TB/Plot/Res_vs_InvEnergy_oldFit_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".eps");

  */

  
}//resolution 

