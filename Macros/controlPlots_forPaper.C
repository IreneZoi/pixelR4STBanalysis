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
#include "TH1.h"
#include "TArrayD.h"
#include <sstream>
#include "fileHandler.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"


bool print=true;
using namespace std;
#define runNONirr "2743"
#define runPirr2 "2801"
#define runNirr4 "3778" //"3839"
#define irradiations 3

void TDR();
void TDR2(TCanvas * c_all, int period=0, int pos= 11);
Double_t ScaleX(Double_t x);
void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t));

void controlPlots_forPaper(TString name = "preliminary_beamdiv")
{

  TString inputDir="/home/zoiirene/Output/";
  TString outputDir="/home/zoiirene/Output/Plots/";

  TString irr[irradiations];
  irr[0] = "Non-irradiated, 120 V"; // "no irr, 5.6 GeV";
  irr[1] = "#phi_{eq} = 2.1 #times 10^{15} cm^{-2}, proton, 800 V";
  irr[2] = "#phi_{eq} = 3.6 #times 10^{15} cm^{-2}, neutron, 800 V";


  std::map<std::pair<TString, TString>, TH1F *>::iterator it;
  TH1F * h_res;  
  TString Label[1] = {"straightTracksY_isoAandCandB_straightTracksX"};
  TString Hist = "clphB";
  bool dphcut = true;
  //TString pitch[1] = {"25"}; //,"25","25"};
  TString pitch = "25"; //,"25","25"};

  // Non-irr 
  MapTH1 landau_map;
  TString Run[1] = {runNONirr};
  int i_dphcut[1] = {12}; 
  Double_t d_dphcut[1] = {12};
  TString ss_dphcut[1] = {"12"};
  TString inputfileNONirr= "beamdiv_A13C14";
  GetHists(&landau_map,1, 1, 1,dphcut, Run, i_dphcut, Label, Hist, h_res,true,inputfileNONirr);
  
  //proton
  Run[0] = runPirr2;
  i_dphcut[0] = 15; 
  d_dphcut[0] = 15;
  ss_dphcut[0] = "15";
  TString  inputfilePirr2="beamdiv_A12C15";
  GetHists(&landau_map,1, 1, 1,dphcut, Run, i_dphcut, Label, Hist, h_res,true,inputfilePirr2);

  //neutron
  Run[0] = runNirr4;
  i_dphcut[0] = 15; 
  d_dphcut[0] = 15;
  ss_dphcut[0] = "15";
  TString inputfileNirr4="beamdiv_A12C13";
  GetHists(&landau_map,1, 1, 1,dphcut, Run, i_dphcut, Label, Hist, h_res,true,inputfileNirr4);

  TH1F  h_landau[irradiations];
  TH1F  h_landau_orig[irradiations];
  
  float integral[irradiations];
  float integral90[irradiations];
  float bin90[irradiations];
  TString ss_dphcuts[irradiations] = {"12","15","15"};
  TString runs[irradiations] = {runNONirr,runPirr2,runNirr4};
  for(int i=0; i<irradiations; i++)
    {
      auto it2 = landau_map.find(std::make_pair(runs[i]+"_"+ss_dphcuts[i]+"_"+Label[0],Hist));
      if(it2  != landau_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
	  h_landau_orig[i] = *it2->second;
          h_landau[i] = *it2->second;
	  integral[i] = h_landau[i].Integral(0, h_landau[i].GetNbinsX()+1);
	  cout << "integral " << integral[i] << endl;
	  
      	  h_landau_orig[i].Scale(1./integral[i]);
	  integral90[i]= 0.9*h_landau[i].Integral();

          int j = 0;
	  while(integral[i]>integral90[i])
	    {
	      //cout << " while "<< j << endl;
	      integral[i] =  h_landau[i].Integral(0, h_landau[i].GetNbinsX()-j);
	      bin90[i] =  h_landau[i].GetBinCenter( h_landau[i].GetNbinsX()-j);
	      //cout << " integral " << integral[i] << " high " << bin90[i] << endl;
	      j++;
	    }
	  cout << " integral90 " << integral90[i] << endl;
	}
      else cout << "map key " << it2->first.first << " not found " << endl;
    }

  /////// plots!

  
  TCanvas *cFDB2 = new TCanvas("cFDB2", "FDB resolution", 600, 600);
  cFDB2->SetLeftMargin(0.15);  
  cFDB2->SetRightMargin(1.5);  
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTextFont(43);
  gStyle->SetTextSize(10);
  gStyle->SetLegendFont(43);
  gStyle->SetLegendTextSize(15);
  cout << "starting making hists pretty" << endl;
  h_landau_orig[0].SetTitle(" ");
  h_landau_orig[0].GetYaxis()->SetTitle("Normalized number of events");
  h_landau_orig[0].GetXaxis()->SetTitle("Cluster charge [ADC]");
  //  h_landau_orig[0].SetMarkerSize(2.5);
  h_landau_orig[0].SetLineColor(kBlack);
  h_landau_orig[0].SetLineStyle(1);
  h_landau_orig[0].SetLineWidth(2);
  h_landau_orig[0].Rebin(2);
  h_landau_orig[0].Scale(0.5);
  h_landau_orig[0].GetXaxis()->SetLimits(0,1000.);
  h_landau_orig[0].GetXaxis()->SetMaxDigits(3); //SetNoExponent(true);
  h_landau_orig[0].GetYaxis()->SetRangeUser(0.0001,1.);
  h_landau_orig[0].GetXaxis()->SetRange(0.,h_landau_orig[0].GetNbinsX()+1);
  h_landau_orig[0].Draw("histe");
  //  TGaxis::SetMaxDigits(3);
  //  h_landau_orig[1].SetMarkerSize(2.5);
  h_landau_orig[1].GetXaxis()->SetLimits(0,1000.);
  //h_landau_orig[1].GetYaxis()->SetRangeUser(0.,0.08);
  h_landau_orig[1].Rebin(2);
  h_landau_orig[1].Scale(0.5);
  h_landau_orig[1].GetYaxis()->SetRangeUser(0.0001,1.);
  h_landau_orig[1].SetLineColor(kGreen+1);
  h_landau_orig[1].SetMarkerColor(kGreen+1);
  h_landau_orig[1].SetLineStyle(1);
  h_landau_orig[1].SetLineWidth(2);
  h_landau_orig[1].GetXaxis()->SetRange(0.,h_landau_orig[1].GetNbinsX()+1);
  //  h_landau_orig[1].GetXaxis()->SetNoExponent(true);
  h_landau_orig[1].Draw("histesame");
 				      
  // h_landau_orig[2].SetMarkerSize(2.5);
  h_landau_orig[2].GetXaxis()->SetLimits(0,1000.);
  h_landau_orig[2].GetYaxis()->SetRangeUser(0.0001,1.);
  h_landau_orig[2].GetXaxis()->SetRange(0.,h_landau_orig[2].GetNbinsX()+1);
  //h_landau_orig[2].GetYaxis()->SetRangeUser(0.,0.08);
  h_landau_orig[2].SetLineColor(kMagenta);
  h_landau_orig[2].SetMarkerColor(kMagenta);
  h_landau_orig[2].Rebin(2);
  h_landau_orig[2].Scale(0.5);
  h_landau_orig[2].SetLineWidth(2);
  h_landau_orig[2].SetLineStyle(1);
  //  h_landau_orig[2].GetXaxis()->SetNoExponent(true);
  h_landau_orig[2].Draw("histesame");

  TLegend* legFDB2 = new TLegend(0.4,0.7,0.8,0.85);
  legFDB2->SetLineColor(0);
  //legFDB2->SetLegendFont(43);
  for(int i =0; i < irradiations; i++)
    legFDB2->AddEntry(&(h_landau_orig[i]),irr[i],"le");

  legFDB2->Draw();

  //landau90[i]
  int color[irradiations] = {1,417,616};
  TLine *  line2[irradiations];
  TArrow *ar[irradiations];
  //  float height[irradiations] = {;
  
  for(int i =0; i < irradiations; i++){
    //    cout << landau90[i] << endl;
    line2[i]    = new TLine( bin90[i],0.,bin90[i],0.015);
    line2[i]->SetLineColor(color[i]);
    line2[i]->SetLineWidth(2);
    line2[i]->SetLineStyle(1);
    line2[i]->Draw("same");
    ar[i]  = new TArrow(bin90[i],0.001,bin90[i]-50,0.001,0.03,"|>"); //,"<|");
    ar[i]->SetLineColor(color[i]);
    ar[i]->SetFillColor(color[i]);
    ar[i]->SetAngle(30);
    
    ar[i]->Draw();

  }
  cFDB2->SetLogy();
  
  TString  outname = outputDir+"Landau_3irr_bestAngle_"+name;
  cFDB2->SaveAs(outname+".eps");
  cFDB2->SaveAs(outname+".png");
  cFDB2->SaveAs(outname+".pdf");
  cFDB2->SaveAs(outname+".root");
  cFDB2->SaveAs(outname+".C");


  // ####################    residual   ################
  Hist = "dx3_clphABC90evR";
  TH1F * hdx3_clchargeABC90evR;
  // Non-irr 
  MapTH1 rmsph_map;
  Run[0] = runNONirr;
  i_dphcut[0] = 12; 
  d_dphcut[0] = 12;
  ss_dphcut[0] = "12";
  //  GetHists(&rmsph_map,1, 1, 1,dphcut, Run, i_dphcut, Label, Hist, h_res);
  //hdx3_clchargeABC90evR =
  rmsph_map = GetCheckHists(&rmsph_map, 1, 1,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,inputfileNONirr);
  
  //proton

  Run[0] = runPirr2;
  i_dphcut[0] = 15; 
  d_dphcut[0] = 15;
  ss_dphcut[0] = "15";
  //  GetHists(&rmsph_map,1, 1, 1,dphcut, Run, i_dphcut, Label, Hist, h_res);
  //hdx3_clchargeABC90evR =
  rmsph_map = GetCheckHists(&rmsph_map, 1, 1,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,inputfilePirr2);
  
  

  //neutron

  Run[0] = runNirr4;
  i_dphcut[0] = 15; 
  d_dphcut[0] = 15;
  ss_dphcut[0] = "15";
  //rmsph_map = GetCheckHists(&rmsph_map, 1, 1,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,inputfileNirr4);
  rmsph_map = GetCheckHists(&rmsph_map, 1, 1,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,inputfileNirr4);



  TH1F  h_resq[irradiations];
  TH1F  h_resq_orig[irradiations];
  //  TH1F  h_resq[irradiations];
  for(int i=0; i<irradiations; i++)
    {
      auto it2 = rmsph_map.find(std::make_pair(runs[i]+"_"+ss_dphcuts[i]+"_"+pitch,Hist));
      //      if( i==2) it2 = rmsph_map.find(std::make_pair(runs[i]+"_"+pitch,Hist));
      if(it2  != rmsph_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
	  h_resq[i] = *it2->second;
	  h_resq_orig[i] = *it2->second;
	  cout << "integral " << h_resq[i].Integral() << endl;
	  //	  landau90[i]= 0.9*h_resq[i].Integral();
	  h_resq[i].Scale(1./h_resq[i].Integral());
	  h_resq[i].Draw();
	  ScaleAxis(h_resq[i].GetXaxis(), ScaleX);
	  
	}
      else cout << "map key " << it2->first.first << " not found " << endl;
    }


  /////// plots!

  
  TCanvas *cresph = new TCanvas("cresph", "FDB resolution", 600, 600);
  cresph->SetLeftMargin(0.12);  
  cresph->SetRightMargin(-0.1);  
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);
  gStyle->SetTextFont(43);
  gStyle->SetTextSize(10);
  
  cout << "starting making hists pretty" << endl;
  h_resq[0].SetTitle(" ");
  cout << "title fine" << endl;
  
  h_resq[0].GetYaxis()->SetTitle("Normalized number of events");
  h_resq[0].GetXaxis()->SetTitle("#Deltax [#mum]");
  //  h_resq[0].SetMarkerSize(2.5);
  h_resq[0].SetLineColor(kBlack);
  h_resq[0].SetLineStyle(1);
  h_resq[0].SetLineWidth(2);
  //  h_resq[0].GetXaxis()->SetNdivisions(20,5,3);
  h_resq[0].GetXaxis()->SetRangeUser(-0.1,0.1);
  h_resq[0].GetXaxis()->SetMaxDigits(3); //SetNoExponent(true);
  h_resq[0].GetYaxis()->SetRangeUser(0.00001,10.);
  h_resq[0].GetXaxis()->SetRange(0.,h_resq[0].GetNbinsX() + 1);
  h_resq[0].Draw("histe");
  cout << " Xaxis " << h_resq[0].GetXaxis()->GetNdivisions() << endl;
  cout << " RMS " << h_resq[0].GetRMS() << endl;
  cresph->Update();
  TPaveStats *Stats =   (TPaveStats*)h_resq[0].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");

  Stats->SetX1NDC(0.15);
  Stats->SetX2NDC(.35);
  Stats->SetY1NDC(.6);
  Stats->SetY2NDC(.7);
  //  gPad->Modified();
  cout << "stats fine" << endl;

  
  //  TGaxis::SetMaxDigits(3);
  //  h_resq[1].SetMarkerSize(2.5);
  h_resq[1].GetXaxis()->SetRangeUser(-0.1,0.1);
  h_resq[1].GetXaxis()->SetRange(0.,h_resq[1].GetNbinsX() + 1);
  h_resq[1].GetYaxis()->SetRangeUser(0.00001,10.);
  h_resq[1].SetLineColor(kGreen+1);
  h_resq[1].SetLineStyle(2);
  h_resq[1].SetLineWidth(2);
  //  h_resq[1].GetXaxis()->SetNoExponent(true);
  h_resq[1].Draw("histesames");
  cout << "second hist fine" << endl;

  cresph->Update();
  TPaveStats *Stats1 =   (TPaveStats*)h_resq[1].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");

  Stats1->SetX1NDC(0.15);
  Stats1->SetX2NDC(.35);
  Stats1->SetY1NDC(.49);
  Stats1->SetY2NDC(.59);
  Stats1->SetLineColor(kGreen+1);
  Stats1->SetLineStyle(2);
  gPad->Modified();

  
  // h_resq[2].SetMarkerSize(2.5);
  h_resq[2].GetXaxis()->SetRangeUser(-0.1,0.1);
  h_resq[2].GetYaxis()->SetRangeUser(0.00001,10.);
  h_resq[2].GetXaxis()->SetRange(0.,h_resq[2].GetNbinsX() + 1);
  h_resq[2].SetLineColor(kMagenta);
  h_resq[2].SetLineWidth(2);
  h_resq[2].SetLineStyle(3);
  //  h_resq[2].GetXaxis()->SetNoExponent(true);
  h_resq[2].Draw("histesames");
  cresph->Update();
  cout << "third hist also" << endl;
  
  TPaveStats *Stats2 =   (TPaveStats*)h_resq[2].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");

  Stats2->SetX1NDC(0.15);
  Stats2->SetX2NDC(.35);
  Stats2->SetY1NDC(.38);
  Stats2->SetY2NDC(.48);
  Stats2->SetLineColor(kMagenta);
  Stats2->SetLineStyle(3);
  gPad->Modified();
  
  TLegend* legresq = new TLegend(0.35,0.72,0.8,0.87);
  legresq->SetLineColor(0);
  //  legresq->SetTextSize(0.03);
  for(int i =0; i < irradiations; i++)
    legresq->AddEntry(&(h_resq[i]),irr[i],"le");
  
  legresq->Draw();

  /*
  int color[irradiations] = {632,1,600};
  TLine *  line2[irradiations];
  TArrow *ar[irradiations];
  
  for(int i =0; i < irradiations; i++){
    cout << landau90[i] << endl;
    line2[i]    = new TLine( landau90[i],0.,landau90[i],0.007);
    line2[i]->SetLineColor(color[i]);
    line2[i]->SetLineWidth(2);
    line2[i]->SetLineStyle(1);
    line2[i]->Draw("same");
    ar[i]  = new TArrow(landau90[i],0.006,landau90[i]-6000,0.006,0.03,"|>"); //,"<|");
    ar[i]->SetLineColor(color[i]);
    ar[i]->SetFillColor(color[i]);
    ar[i]->SetAngle(30);
    
    ar[i]->Draw();

  }
  */

  cresph->SetLogy();
  
  outname = outputDir+"ResPHAll_3irr_bestAngle_"+name;
  cresph->SaveAs(outname+".eps");
  cresph->SaveAs(outname+".png");
  cresph->SaveAs(outname+".pdf");
  cresph->SaveAs(outname+".root");
  cresph->SaveAs(outname+".C");






    /////// non log version
  Double_t perc[irradiations];
  Double_t rmserr[irradiations];
  Double_t minX[irradiations];
  Double_t maxX[irradiations];
  Double_t RMS[irradiations];
  for( int i = 0 ; i < irradiations; i ++)  FitTH1(&(h_resq_orig[i]), &(RMS[i]), &(rmserr[i]),runs[i]+"_", "A", "B", "C", "RMSself", &(perc[i]),&(minX[i]),&(maxX[i]));

  
  TCanvas *cresph2 = new TCanvas("cresph2", "FDB resolution", 600, 600);
  cresph2->SetLeftMargin(0.12);  
  cresph2->SetRightMargin(-0.1);  
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTextFont(43);
  gStyle->SetTextSize(10);
  int colors[irradiations]={1,417,616};
    
  Float_t Mean[irradiations];


  TLegend* legB = new TLegend(0.25,0.85,0.4,0.88);
  legB->SetLineColor(0);
  legB->AddEntry(&(h_resq[0]),"#theta = 8.8 deg","");
  
  TLegend* legresq2 = new TLegend(0.13,0.67,0.69,0.85);
  legresq2->SetLineColor(0);

  
  cout << "starting making hists pretty" << endl;
  for( int i = 0 ; i < irradiations; i ++){
    cout << "            " << irr[i] << endl;
    h_resq[i].SetTitle(" ");
    h_resq[i].GetYaxis()->SetTitle("Normalized number of events");
    h_resq[i].GetXaxis()->SetTitle("#Deltax [#mum]");

    h_resq[i].SetMarkerColor(color[i]);
    h_resq[i].SetLineColor(color[i]);
    h_resq[i].SetLineStyle(1+i);
    h_resq[i].SetLineWidth(2);
    h_resq[i].GetXaxis()->SetRangeUser(minX[i],maxX[i]);
    h_resq[i].GetYaxis()->SetRangeUser(0.,0.35);
    
    h_resq[i].Draw("histesames");
    cresph->Update();
    cout << " changed range " <<endl;
    Mean[i]=h_resq[i].GetMean();
    cout << " mean " << Mean[i] << endl;
    cout << " rms "<< RMS[i] << endl;


    std::stringstream stream_m;
    stream_m << std::fixed << std::setprecision(2) << Mean[i];
    std::string s_m = stream_m.str();
    std::stringstream stream_r;
    stream_r << std::fixed << std::setprecision(2) << RMS[i];
    std::string s_r = stream_r.str();
    std::stringstream stream_p;
    stream_p << std::fixed << std::setprecision(1) << perc[i]*100;
    std::string s_p = stream_p.str();
    
    if(i!=0)  gPad->Modified();
    irr[0] = "Non-irradiated";//, 120 V"; // "no irr, 5.6 GeV";
    irr[1] = "#phi_{eq} = 2.1 #times 10^{15} cm^{-2}, proton";//, 800 V";
    irr[2] = "#phi_{eq} = 3.6 #times 10^{15} cm^{-2}, neutron";//, 800 V";

    legresq2->AddEntry(&(h_resq[i]),irr[i],"le"); // +", #mu = "+s_m+" #mum, RMS = "+s_r+" #mum, tracks: "+s_p+" %","le");
    irr[i] = "";
    legresq2->AddEntry(&(h_resq[i]),irr[i]+"#mu = "+s_m+" #mum, RMS = "+s_r+" #mum, tracks: "+s_p+" %","");
    h_resq[i].GetXaxis()->SetRangeUser(-50.,50.);

  }
    
   
  legB->Draw();
  
  legresq2->Draw();


  for (int i = 0 ; i < irradiations; i ++){
    TLine *  linemin = new TLine( minX[i],0.,minX[i],0.2);
    linemin->SetLineColor(color[i]);
    linemin->SetLineWidth(2);
    linemin->SetLineStyle(1+i);
    linemin->Draw("same");

    TArrow *armin = new TArrow(minX[i],0.1,minX[i]+6,0.1,0.03,"|>"); //,"<|");
    armin->SetLineColor(color[i]);
    armin->SetFillColor(color[i]);
    armin->SetAngle(30);
    armin->Draw();

    TLine *  linemax = new TLine( maxX[i],0.,maxX[i],0.2);
    linemax->SetLineColor(color[i]);
    linemax->SetLineWidth(2);
    linemax->SetLineStyle(1+i);
    linemax->Draw("same");

    TArrow *armax = new TArrow(maxX[i],0.1,maxX[i]-6,0.1,0.03,"|>"); //,"<|");
    armax->SetLineColor(color[i]);
    armax->SetFillColor(color[i]);
    armax->SetAngle(30);
    armax->Draw();

  }
  


  
  outname = outputDir+"ResPHAllnoLog_3irr_bestAngle_"+name;
  cresph2->SaveAs(outname+".eps");
  cresph2->SaveAs(outname+".png");
  cresph2->SaveAs(outname+".pdf");
  cresph2->SaveAs(outname+".root");
  cresph2->SaveAs(outname+".C");



  ////// clsize 2 versions
  // version 90 cut
  Hist = "nrowB_clphABC90evR";

  MapTH1 nr_map;
  Run[0] = runNONirr;
  i_dphcut[0] = 12;
  d_dphcut[0] = 12;
  ss_dphcut[0] = "12";
  //  GetHists(&nr_map,1, 1, 1,dphcut, Run, i_dphcut, Label, Hist, h_res);
  //hdx3_clchargeABC90evR =
  nr_map = GetCheckHists(&nr_map, 1, 1,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,inputfileNONirr);

  //proton

  Run[0] = runPirr2;
  i_dphcut[0] = 15;
  d_dphcut[0] = 15;
  ss_dphcut[0] = "15";
  //  GetHists(&nr_map,1, 1, 1,dphcut, Run, i_dphcut, Label, Hist, h_res);
  //hdx3_clchargeABC90evR =
  nr_map = GetCheckHists(&nr_map, 1, 1,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,inputfilePirr2);



  //neutron

  Run[0] = runNirr4;
  i_dphcut[0] = 15;
  d_dphcut[0] = 15;
  ss_dphcut[0] = "15";
  //nr_map = GetCheckHists(&nr_map, 1, 1,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,inputfileNirr4);
  nr_map = GetCheckHists(&nr_map, 1, 1,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,inputfileNirr4);

  
  TH1F  h_NrowB[irradiations];
  for(int i=0; i<irradiations; i++)
    {
      auto it2 = nr_map.find(std::make_pair(runs[i]+"_"+ss_dphcuts[i]+"_"+pitch,Hist));
      //auto it2 = res_map.find(std::make_pair(Run[i]+"_"+ss_dphcut[k]+"_"+pitch,Hist));
      if(it2  != nr_map.end())
	{
	    if(print)               cout << " found map " << endl;
	    if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
	    cout << " Getting clusters " << endl;

	    h_NrowB[i] = * it2->second;
	    //	    h_NrowB[i].Sumw2();
	    cout << h_NrowB[i].Integral() << endl;
	    cout << 1./h_NrowB[i].Integral() << endl;
	    Double_t norm = 1./h_NrowB[i].Integral();
	    if (h_NrowB[i].GetSumw2N() == 0) h_NrowB[i].Sumw2(kTRUE);
	    h_NrowB[i].Scale(norm,"width"); //1./h_NrowB[i].Integral());
	    cout << h_NrowB[i].Integral() << endl;
	    if (h_NrowB[i].GetSumw2N() == 0) h_NrowB[i].Sumw2(kTRUE);

	}

    }


  TCanvas *cnB90 = new TCanvas("cnB90", "FDB resolution", 600, 600);
  cnB90->SetLeftMargin(0.12);  
  cnB90->SetRightMargin(-0.1);  
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);
  gStyle->SetTextFont(43);
  gStyle->SetTextSize(10);
  
  cout << "starting making hists pretty" << endl;
  h_NrowB[0].SetTitle(" ");
  h_NrowB[0].GetYaxis()->SetTitle("Normalized number of events");
  h_NrowB[0].GetXaxis()->SetTitle("Cluster size on central plane [pixels]");
  //  h_NrowB[0].SetMarkerSize(2.5);
  h_NrowB[0].SetLineColor(kBlack);
  h_NrowB[0].SetLineStyle(1);
  h_NrowB[0].SetLineWidth(2);
  //h_NrowB[0].GetXaxis()->SetNdivisions(20,5,3);
  h_NrowB[0].GetXaxis()->SetRangeUser(0,10);
  //h_NrowB[0].GetXaxis()->SetMaxDigits(3); //SetNoExponent(true);
  //h_NrowB[0].GetYaxis()->SetRangeUser(0.,40000);
  h_NrowB[0].GetYaxis()->SetLimits(0.,0.1); //RangeUser(0.,40000);
  //  h_NrowB[0].GetXaxis()->SetRange(0.,h_NrowB[0].GetNbinsX() + 1);
  h_NrowB[0].Draw("histe");
  cout << " Xaxis " << h_NrowB[0].GetXaxis()->GetNdivisions() << endl;
  cout << " RMS " << h_NrowB[0].GetRMS() << endl;
  cnB90->Update();
  //TPaveStats *
  Stats =   (TPaveStats*)h_NrowB[0].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");

  Stats->SetX1NDC(0.65);
  Stats->SetX2NDC(.85);
  Stats->SetY1NDC(.55);
  Stats->SetY2NDC(.65);
  //  gPad->Modified();

  
  //  TGaxis::SetMaxDigits(3);
  //  h_NrowB[1].SetMarkerSize(2.5);
  //h_NrowB[1].GetXaxis()->SetRangeUser(-0.1,0.1);
  //  h_NrowB[1].GetXaxis()->SetRange(0.,h_NrowB[1].GetNbinsX() + 1);
  //h_NrowB[1].GetYaxis()->SetRangeUser(0.00001,10.);
  h_NrowB[1].SetLineColor(kGreen+1);
  h_NrowB[1].SetLineStyle(2);
  h_NrowB[1].SetLineWidth(2);
  //  h_NrowB[1].GetXaxis()->SetNoExponent(true);
  h_NrowB[1].Draw("histesames");
  cnB90->Update();
  //TPaveStats *
  Stats1 =   (TPaveStats*)h_NrowB[1].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");

  Stats1->SetX1NDC(0.65);
  Stats1->SetX2NDC(.85);
  Stats1->SetY1NDC(.44);
  Stats1->SetY2NDC(.54);
  Stats1->SetLineColor(kGreen+1);
  Stats1->SetLineStyle(2);
  gPad->Modified();

  
  // h_NrowB[2].SetMarkerSize(2.5);
  //h_NrowB[2].GetXaxis()->SetRangeUser(-0.1,0.1);
  //h_NrowB[2].GetYaxis()->SetRangeUser(0.00001,10.);
  //h_NrowB[2].GetXaxis()->SetRange(0.,h_NrowB[2].GetNbinsX() + 1);
  h_NrowB[2].SetLineColor(kMagenta);
  h_NrowB[2].SetLineWidth(2);
  h_NrowB[2].SetLineStyle(3);
  //  h_NrowB[2].GetXaxis()->SetNoExponent(true);
  h_NrowB[2].Draw("histesames");
  cnB90->Update();
  //TPaveStats *
  Stats2 =   (TPaveStats*)h_NrowB[2].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");

  Stats2->SetX1NDC(0.65);
  Stats2->SetX2NDC(.85);
  Stats2->SetY1NDC(.33);
  Stats2->SetY2NDC(.43);
  Stats2->SetLineColor(kMagenta);
  Stats2->SetLineStyle(3);
  gPad->Modified();
  
  TLegend* legnB90 = new TLegend(0.35,0.72,0.8,0.87);
  legnB90->SetLineColor(0);
  //  legnB90->SetTextSize(0.03);
  for(int i =0; i < irradiations; i++)
    legnB90->AddEntry(&(h_NrowB[i]),irr[i],"le");
  
  legnB90->Draw();

  /*
  int color[irradiations] = {632,1,600};
  TLine *  line2[irradiations];
  TArrow *ar[irradiations];
  
  for(int i =0; i < irradiations; i++){
    cout << landau90[i] << endl;
    line2[i]    = new TLine( landau90[i],0.,landau90[i],0.007);
    line2[i].SetLineColor(color[i]);
    line2[i].SetLineWidth(2);
    line2[i].SetLineStyle(1);
    line2[i].Draw("same");
    ar[i]  = new TArrow(landau90[i],0.006,landau90[i]-6000,0.006,0.03,"|>"); //,"<|");
    ar[i].SetLineColor(color[i]);
    ar[i].SetFillColor(color[i]);
    ar[i].SetAngle(30);
    
    ar[i].Draw();

  }
  */

 
  
  outname = outputDir+"ClB90PHAll_3irr_bestAngle_"+name;
  cnB90->SaveAs(outname+".eps");
  cnB90->SaveAs(outname+".png");
  cnB90->SaveAs(outname+".pdf");
  cnB90->SaveAs(outname+".root");
  cnB90->SaveAs(outname+".C");





  

  /// no 90 cut

  Hist = "nrowB";

  cout << " Getting " << Hist << endl;

  // GetHists(&res_map,Angles, dphcuts, comparisons,dphcut, Run, i_dphcut, Label, Hist, h_res,true,label);

  MapTH1 nrB_map;
  Run[0] = runNONirr;
  i_dphcut[0] = 12;
  d_dphcut[0] = 12;
  ss_dphcut[0] = "12";
  //  GetHists(&nr_map,1, 1, 1,dphcut, Run, i_dphcut, Label, Hist, h_res);
  //hdx3_clargeABC90evR =
  //  void GetHists(MapTH1 * map,int runs, int dphcuts, int comparisons,  bool dphcut, TString * Run, int * i_dphcut, TString * Label, TString Hist, TH1F * h, bool extraname = false, TString name = " "){
  GetHists(&nrB_map, 1, 1,1,dphcut, Run,i_dphcut, Label,Hist,hdx3_clchargeABC90evR,true,inputfileNONirr);

    
  //proton

  Run[0] = runPirr2;
  i_dphcut[0] = 15;
  d_dphcut[0] = 15;
  ss_dphcut[0] = "15";
  //  GetHists(&nr_map,1, 1, 1,dphcut, Run, i_dphcut, Label, Hist, h_res);
  //hdx3_clchargeABC90evR =
  GetHists(&nrB_map, 1, 1,1,dphcut, Run,i_dphcut,Label, Hist,hdx3_clchargeABC90evR,true,inputfilePirr2);



  //neutron

  Run[0] = runNirr4;
  i_dphcut[0] = 15;
  d_dphcut[0] = 15;
  ss_dphcut[0] = "15";
  //nr_map = GetCheckHists(&nr_map, 1, 1,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,inputfileNirr4);
  GetHists(&nrB_map, 1, 1, 1,dphcut, Run,i_dphcut,Label, Hist,hdx3_clchargeABC90evR,true,inputfileNirr4);

  for(int i=0; i<irradiations; i++)
    {
      auto it2 = nrB_map.find(std::make_pair(runs[i]+"_"+ss_dphcuts[i]+"_"+Label[0],Hist));
      //auto it2 = res_map.find(std::make_pair(Run[i]+"_"+ss_dphcut[k]+"_"+pitch,Hist));
      if(it2  != nrB_map.end())
	{
	    if(print)               cout << " found map " << endl;
	    if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
	    cout << " Getting clusters " << endl;


	    h_NrowB[i] =* it2->second; 

	  }

    }


  TCanvas *cnB = new TCanvas("cnB", "FDB resolution", 600, 600);
  cnB->SetLeftMargin(0.12);  
  cnB->SetRightMargin(-0.1);  
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);
  gStyle->SetTextFont(43);
  gStyle->SetTextSize(10);
  
  cout << "starting making hists pretty" << endl;
  h_NrowB[0].SetTitle(" ");
  h_NrowB[0].GetYaxis()->SetTitle("Normalized number of events");
  h_NrowB[0].GetXaxis()->SetTitle("Cluster size on central plane [pixels]");
  //  h_NrowB[0].SetMarkerSize(2.5);
  h_NrowB[0].SetLineColor(kBlack);
  h_NrowB[0].SetLineStyle(1);
  h_NrowB[0].SetLineWidth(2);
  //h_NrowB[0].GetXaxis()->SetNdivisions(20,5,3);
  h_NrowB[0].GetXaxis()->SetRangeUser(0,10);
  h_NrowB[0].GetXaxis()->SetMaxDigits(3); //SetNoExponent(true);
  h_NrowB[0].GetYaxis()->SetRangeUser(0.,40000.);
  h_NrowB[0].GetXaxis()->SetRange(0.,h_NrowB[0].GetNbinsX() + 1);
  h_NrowB[0].Draw("histe");
  cout << " Xaxis " << h_NrowB[0].GetXaxis()->GetNdivisions() << endl;
  cout << " RMS " << h_NrowB[0].GetRMS() << endl;
  cnB->Update();
  //TPaveStats *
  Stats =   (TPaveStats*)h_NrowB[0].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");

  Stats->SetX1NDC(0.15);
  Stats->SetX2NDC(.35);
  Stats->SetY1NDC(.6);
  Stats->SetY2NDC(.7);
  //  gPad->Modified();

  
  //  TGaxis::SetMaxDigits(3);
  //  h_NrowB[1].SetMarkerSize(2.5);
  h_NrowB[1].GetXaxis()->SetRangeUser(-0.1,0.1);
  h_NrowB[1].GetXaxis()->SetRange(0.,h_NrowB[1].GetNbinsX() + 1);
  h_NrowB[1].GetYaxis()->SetRangeUser(0.00001,10.);
  h_NrowB[1].SetLineColor(kGreen+1);
  h_NrowB[1].SetLineStyle(2);
  h_NrowB[1].SetLineWidth(2);
  //  h_NrowB[1].GetXaxis()->SetNoExponent(true);
  h_NrowB[1].Draw("histesames");
  cnB->Update();
  //TPaveStats *
  Stats1 =   (TPaveStats*)h_NrowB[1].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");

  Stats1->SetX1NDC(0.15);
  Stats1->SetX2NDC(.35);
  Stats1->SetY1NDC(.49);
  Stats1->SetY2NDC(.59);
  Stats1->SetLineColor(kGreen+1);
  Stats1->SetLineStyle(2);
  gPad->Modified();

  
  // h_NrowB[2].SetMarkerSize(2.5);
  h_NrowB[2].GetXaxis()->SetRangeUser(-0.1,0.1);
  h_NrowB[2].GetYaxis()->SetRangeUser(0.00001,10.);
  h_NrowB[2].GetXaxis()->SetRange(0.,h_NrowB[2].GetNbinsX() + 1);
  h_NrowB[2].SetLineColor(kMagenta);
  h_NrowB[2].SetLineWidth(2);
  h_NrowB[2].SetLineStyle(3);
  //  h_NrowB[2].GetXaxis()->SetNoExponent(true);
  h_NrowB[2].Draw("histesames");
  cnB->Update();
  //TPaveStats *
  Stats2 =   (TPaveStats*)h_NrowB[2].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");

  Stats2->SetX1NDC(0.15);
  Stats2->SetX2NDC(.35);
  Stats2->SetY1NDC(.38);
  Stats2->SetY2NDC(.48);
  Stats2->SetLineColor(kMagenta);
  Stats2->SetLineStyle(3);
  gPad->Modified();
  /*
  TLegend* legnB90 = new TLegend(0.35,0.72,0.8,0.87);
  legnB90->SetLineColor(0);
  //  legnB90->SetTextSize(0.03);
  for(int i =0; i < irradiations; i++)
    legnB90->AddEntry(&(h_NrowB[i]),irr[i],"le");
  */
  legnB90->Draw();

  /*
  int color[irradiations] = {632,1,600};
  TLine *  line2[irradiations];
  TArrow *ar[irradiations];
  
  for(int i =0; i < irradiations; i++){
    cout << landau90[i] << endl;
    line2[i]    = new TLine( landau90[i],0.,landau90[i],0.007);
    line2[i].SetLineColor(color[i]);
    line2[i].SetLineWidth(2);
    line2[i].SetLineStyle(1);
    line2[i].Draw("same");
    ar[i]  = new TArrow(landau90[i],0.006,landau90[i]-6000,0.006,0.03,"|>"); //,"<|");
    ar[i].SetLineColor(color[i]);
    ar[i].SetFillColor(color[i]);
    ar[i].SetAngle(30);
    
    ar[i].  Draw();

  }
  */

 
  
  outname = outputDir+"ClB_PHAll_3irr_bestAngle_"+name;
  cnB->SaveAs(outname+".eps");
  cnB->SaveAs(outname+".png");
  cnB->SaveAs(outname+".pdf");
  cnB->SaveAs(outname+".root");
  cnB->SaveAs(outname+".C");






  
}//resolution 



Double_t ScaleX(Double_t x)
{
  Double_t v;
  v = 1000 * x; // "linear scaling" function example
  return v;
}

void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t))
{
  if (!a) return; // just a precaution
  if (a->GetXbins()->GetSize())
    {
      // an axis with variable bins
      // note: bins must remain in increasing order, hence the "Scale"
      // function must be strictly (monotonically) increasing
      TArrayD X(*(a->GetXbins()));
      for(Int_t i = 0; i < X.GetSize(); i++) X[i] = Scale(X[i]);
      a->Set((X.GetSize() - 1), X.GetArray()); // new Xbins
    }
  else
    {
      // an axis with fix bins
      // note: we modify Xmin and Xmax only, hence the "Scale" function
      // must be linear (and Xmax must remain greater than Xmin)
      a->Set( a->GetNbins(),
	      Scale(a->GetXmin()), // new Xmin
	      Scale(a->GetXmax()) ); // new Xmax
    }
  return;
}



void TDR()
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_sqrtS = "";
  bool drawLogo = true;
  int H_ref = 600;
  int W_ref = 600;
  float T = 0.07*H_ref;//0.07
  float B = 0.11*H_ref;//0.12
  float L = 0.12*W_ref;
  float R = 0.01*W_ref;


}

void TDR2(TCanvas * c_all, int period = 0, int pos = 11)
{
  CMS_lumi( c_all, period,pos);
  //  CMS_lumi( c_all, 0, 11);
  c_all->Update();
  c_all->RedrawAxis();
}
