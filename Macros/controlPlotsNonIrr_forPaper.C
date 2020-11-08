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
#include <iomanip>
#include <sstream>
#include "fileHandler.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"


bool print=true;
using namespace std;
#define runNONirr0 "2735"
#define runNONirrB "2743"
#define runNONirrS "2758"

#define angles 3

void TDR();
void TDR2(TCanvas * c_all, int period=0, int pos= 11);
Double_t ScaleX(Double_t x);
void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t));

void controlPlotsNonIrr_forPaper(TString name = "preliminary_dycut")
{

  TString inputDir="/home/zoiirene/Output/";
  TString outputDir="/home/zoiirene/Output/Plots/";
  TString base="Non-irradiated, 120 V, 5.6 GeV";  
  //ang[1] = "#phi_{eq} = 2.1 #times 10^{15} cm^{-2}, proton, 800 V";

  TString ang[angles];
  ang[0] = "0 deg"; 
  ang[1] = "8.8 deg";
  ang[2] = "27.5 deg";


  std::map<std::pair<TString, TString>, TH1F *>::iterator it;
  TH1F * h_res;  
  TString Label[1] = {"straightTracksY_isoAandCandB_straightTracksX"};
  TString Hist = "clphB";
  bool dphcut = true;
  //TString pitch[1] = {"25"}; //,"25","25"};
  TString pitch = "25"; //,"25","25"};

  // Non-irr 
  MapTH1 landau_map;
  TString Runs[angles] = {runNONirr0,runNONirrB,runNONirrS};
  int i_dphcut[1] = {12};
  Double_t d_dphcut[1] = {12};
  TString ss_dphcut[1] = {"12"};
  TString inputfileNONirr= "dycut_A13C14";

  GetHists(&landau_map,angles, 1, 1,dphcut, Runs, i_dphcut, Label, Hist, h_res,true,inputfileNONirr);
  
  int color[angles] = {920,1,803};
      
  TH1F  h_landau[angles];
  TH1F  h_landau_orig[angles];
  
  float integral[angles];
  float integral90[angles];
  float bin90[angles];
  TString ss_dphcuts[angles] = {"12","12","12"};
  TString runs[angles] = {runNONirr0,runNONirrB,runNONirrS};
  for(int i=0; i<angles; i++)
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
  h_landau_orig[0].SetLineColor(color[0]);
  h_landau_orig[0].SetMarkerColor(color[0]);
  h_landau_orig[0].SetLineStyle(1);
  h_landau_orig[0].SetLineWidth(2);
  h_landau_orig[0].Rebin(2);
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
  h_landau_orig[1].GetYaxis()->SetRangeUser(0.0001,1.);
  h_landau_orig[1].SetLineColor(color[1]);
  //  h_landau_orig[1].SetMarkerColor(kGray+3);
  h_landau_orig[1].SetLineStyle(2);
  h_landau_orig[1].SetLineWidth(2);
  h_landau_orig[1].GetXaxis()->SetRange(0.,h_landau_orig[1].GetNbinsX()+1);
  //  h_landau_orig[1].GetXaxis()->SetNoExponent(true);
  h_landau_orig[1].Draw("histesame");
 				      
  // h_landau_orig[2].SetMarkerSize(2.5);
  h_landau_orig[2].GetXaxis()->SetLimits(0,1000.);
  h_landau_orig[2].GetYaxis()->SetRangeUser(0.0001,1.);
  h_landau_orig[2].GetXaxis()->SetRange(0.,h_landau_orig[2].GetNbinsX()+1);
  //h_landau_orig[2].GetYaxis()->SetRangeUser(0.,0.08);
  h_landau_orig[2].SetLineColor(color[2]);
  h_landau_orig[2].SetMarkerColor(color[2]);
  h_landau_orig[2].Rebin(2);
  h_landau_orig[2].SetLineWidth(2);
  h_landau_orig[2].SetLineStyle(3);
  //  h_landau_orig[2].GetXaxis()->SetNoExponent(true);
  h_landau_orig[2].Draw("histesame");

  TLegend* legFDB = new TLegend(0.3,0.8,0.7,0.85);
  legFDB->SetLineColor(0);
  //legFDB2->SetLegendFont(43);
  legFDB->AddEntry(&(h_landau_orig[0]),base,"");
  legFDB->Draw();
  TLegend* legFDB2 = new TLegend(0.4,0.7,0.8,0.8);
  legFDB2->SetLineColor(0);
  for(int i =0; i < angles; i++){

    legFDB2->AddEntry(&(h_landau_orig[i]),ang[i],"le");
  }
  legFDB2->Draw();

  //landau90[i]
  //int color[angles] = {923,1,921}; //417,616};
  TLine *  line2[angles];
  TArrow *ar[angles];
  //  float height[angles] = {;
  
  for(int i =0; i < angles; i++){
    //    cout << landau90[i] << endl;
    line2[i]    = new TLine( bin90[i],0.,bin90[i],0.015);
    line2[i]->SetLineColor(color[i]);
    line2[i]->SetLineWidth(2);
    line2[i]->SetLineStyle(i+1);
    line2[i]->Draw("same");
    ar[i]  = new TArrow(bin90[i],0.001,bin90[i]-50,0.001,0.03,"|>"); //,"<|");
    ar[i]->SetLineColor(color[i]);
    ar[i]->SetFillColor(color[i]);
    ar[i]->SetAngle(30);
    
    ar[i]->Draw();

  }
  cFDB2->SetLogy();
  
  TString  outname = outputDir+"Landau_NonIrr_3Angles_"+name;
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
  rmsph_map = GetCheckHists(&rmsph_map, angles, 1,dphcut, Runs,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,inputfileNONirr);
  
  TH1F  h_resq[angles];   
  //  TH1F  h_resq[angles];
  for(int i=0; i<angles; i++)
    {
      auto it2 = rmsph_map.find(std::make_pair(runs[i]+"_"+ss_dphcuts[i]+"_"+pitch,Hist));
      //    if( i==2) it2 = rmsph_map.find(std::make_pair(runs[i]+"_"+pitch,Hist));
      if(it2  != rmsph_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
	  h_resq[i] = *it2->second;
	  cout << "integral " << h_resq[i].Integral() << endl;
	  //	  landau90[i]= 0.9*h_resq[i].Integral();
	  h_resq[i].Scale(1./h_resq[i].Integral());
	  h_resq[i].Draw();
	  ScaleAxis(h_resq[i].GetXaxis(), ScaleX);
	  
	}
      else cout << "map key " << it2->first.first << " not found " << endl;
    }


  /////// plots!
  TPaveStats *Stats[angles];
  Float_t Mean[angles];
  Float_t RMS[angles];

  TString Mean_ss[angles];
  TString RMS_ss[angles];
  

  TCanvas *cresph = new TCanvas("cresph", "FDB resolution", 600, 600);
  cresph->SetLeftMargin(0.12);  
  cresph->SetRightMargin(-0.1);  
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  //gStyle->SetOptStat(111110);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTextFont(43);
  gStyle->SetTextSize(10);

  TLegend* legresq = new TLegend(0.25,0.73,0.7,0.82);
  legresq->SetLineColor(0);

  
  cout << "starting making hists pretty" << endl;
  for( int i = 0 ; i < angles; i ++){
    cout << "            " << ang[i] << endl;
    h_resq[i].SetTitle(" ");
    h_resq[i].GetYaxis()->SetTitle("Normalized number of events");
    h_resq[i].GetXaxis()->SetTitle("#Deltax [#mum]");

    h_resq[i].SetLineColor(color[i]);
    h_resq[i].SetLineStyle(1+i);
    h_resq[i].SetLineWidth(2);
    //Mean[i]=h_resq[i].GetMean();
    //RMS[i]=h_resq[i].GetRMS();
    
    //cout << " mean " << Mean[i] << endl;
    //cout << " rms "<< RMS[i] << endl;
    //h_resq[i].GetXaxis()->SetRange(0,h_resq[i].GetNbinsX()+1);
    h_resq[i].GetXaxis()->SetRangeUser(-200.,200.);
    h_resq[i].GetYaxis()->SetRangeUser(0.00001,10.);
    //h_resq[i].GetYaxis()->SetRangeUser(0.1,30000.);
    h_resq[i].Draw("histesames");
    cresph->Update();
    cout << " changed range " <<endl;
    Mean[i]=h_resq[i].GetMean();
    RMS[i]=h_resq[i].GetRMS();
    cout << " mean " << Mean[i] << endl;
    cout << " rms "<< RMS[i] << endl;
    /*
    Stats[i] =   (TPaveStats*)h_resq[i].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");

    Stats[i]->SetX1NDC(0.15);
    Stats[i]->SetX2NDC(.35);
    Stats[i]->SetY1NDC(.6-i*0.11);
    Stats[i]->SetY2NDC(.7-i*0.11);
    */

    std::stringstream stream_m;
    stream_m << std::fixed << std::setprecision(2) << Mean[i];
    std::string s_m = stream_m.str();
    std::stringstream stream_r;
    stream_r << std::fixed << std::setprecision(2) << RMS[i];
    std::string s_r = stream_r.str();

    //Mean_ss[i].Form("%f",Mean[i]);
    //    RMS_ss[i].Form("%f",RMS[i]);
    if(i!=0)  gPad->Modified();
    legresq->AddEntry(&(h_resq[i]),ang[i]+", #mu = "+s_m+" #mum, RMS = "+s_r+" #mum","le");
  }

  TLegend* legB = new TLegend(0.2,0.83,0.6,0.88);
  legB->SetLineColor(0);
  //legFDB2->SetLegendFont(43);
  legB->AddEntry(&(h_landau_orig[0]),base,"");
  legB->Draw();
  
  
  //  legresq->SetTextSize(0.03);


  
  legresq->Draw();


  cresph->SetLogy();
  
  outname = outputDir+"ResPHAll_NonIrr_3Angles_"+name;
  cresph->SaveAs(outname+".eps");
  cresph->SaveAs(outname+".png");
  cresph->SaveAs(outname+".pdf");
  cresph->SaveAs(outname+".root");
  cresph->SaveAs(outname+".C");






    /////// non log version

  
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
  TLegend* legresq2 = new TLegend(0.25,0.73,0.7,0.82);//0.3,0.7,0.8,0.8);
  legresq2->SetLineColor(0);
  Double_t minrange[angles]={-0.0406891*1000.,-0.0233376*1000.,-0.0436472*1000};
  Double_t maxrange[angles]={0.0409722*1000.,0.024754*1000.,0.0474025*1000.};
  cout << "starting making hists pretty" << endl;
for( int i = 0 ; i < angles; i ++){
    cout << "            " << ang[i] << endl;
    h_resq[i].GetXaxis()->SetRangeUser(minrange[i],maxrange[i]);
    h_resq[i].GetYaxis()->SetRangeUser(0.,0.35);
    //h_resq[i].GetYaxis()->SetRangeUser(0.1,30000.);
    h_resq[i].Draw("histesames");
    cresph->Update();
    cout << " changed range " <<endl;
    Mean[i]=h_resq[i].GetMean();
    RMS[i]=h_resq[i].GetRMS();
    cout << " mean " << Mean[i] << endl;
    cout << " rms "<< RMS[i] << endl;
    /*
    Stats[i] =   (TPaveStats*)h_resq[i].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");

    Stats[i]->SetX1NDC(0.15);
    Stats[i]->SetX2NDC(.35);
    Stats[i]->SetY1NDC(.6-i*0.11);
    Stats[i]->SetY2NDC(.7-i*0.11);
    */

    
    std::stringstream stream_m;
    stream_m << std::fixed << std::setprecision(2) << Mean[i];
    std::string s_m = stream_m.str();
    std::stringstream stream_r;
    stream_r << std::fixed << std::setprecision(2) << RMS[i];
    std::string s_r = stream_r.str();

    //Mean_ss[i].Form("%f",Mean[i]);
    //    RMS_ss[i].Form("%f",RMS[i]);
    if(i!=0)  gPad->Modified();
    legresq2->AddEntry(&(h_resq[i]),ang[i]+", #mu = "+s_m+" #mum, RMS = "+s_r+" #mum","le");
    h_resq[i].GetXaxis()->SetRangeUser(-50.,50.);

  }

/*
  //  TGaxis::SetMaxDigits(3);
  //  h_resq[1].SetMarkerSize(2.5);
  h_resq[1].GetYaxis()->SetTitle("Normalized number of events");
  h_resq[1].GetXaxis()->SetTitle("#Deltax [#mum]");
  
  h_resq[1].GetXaxis()->SetRangeUser(-50,50);
  //h_resq[1].GetXaxis()->SetRange(0.,h_resq[1].GetNbinsX() + 1);
  h_resq[1].GetYaxis()->SetRangeUser(0.0,0.25);
  //h_resq[1].GetYaxis()->SetRangeUser(0.00001,1.);
  h_resq[1].Draw("histesame");
 				      
  h_resq[2].GetXaxis()->SetRangeUser(-25,25);
  h_resq[2].GetYaxis()->SetRangeUser(0.0,0.35);
  //h_resq[2].GetXaxis()->SetRange(0.,h_resq[2].GetNbinsX() + 1);
  h_resq[2].Draw("histesame");
*/
  legB->Draw();
  //  legresq2->SetTextSize(0.03);
  //for(int i =0; i < angles; i++)
  //  legresq2->AddEntry(&(h_resq[1]),ang[1]+", "+base,"le");
  
  legresq2->Draw();
  
  // minrange to 0.0260115 //final range of the resolution calculation with self consistent rms 

  for (int i = 0 ; i < angles; i ++){
    TLine *  linemin = new TLine( minrange[i],0.,minrange[i],0.2);
    linemin->SetLineColor(color[i]);
    linemin->SetLineWidth(2);
    linemin->SetLineStyle(1+i);
    linemin->Draw("same");

    TArrow *armin = new TArrow(minrange[i],0.1,minrange[i]+6,0.1,0.03,"|>"); //,"<|");
    armin->SetLineColor(color[i]);
    armin->SetFillColor(color[i]);
    armin->SetAngle(30);
    armin->Draw();

    TLine *  linemax = new TLine( maxrange[i],0.,maxrange[i],0.2);
    linemax->SetLineColor(color[i]);
    linemax->SetLineWidth(2);
    linemax->SetLineStyle(1+i);
    linemax->Draw("same");
  
    TArrow *armax = new TArrow(maxrange[i],0.1,maxrange[i]-6,0.1,0.03,"|>"); //,"<|");
    armax->SetLineColor(color[i]);
    armax->SetFillColor(color[i]);
    armax->SetAngle(30);
    armax->Draw();

  }
    //cresph2->SetLogy();
     
  
  outname = outputDir+"ResPHAllnoLog_NonIrr_3Angles_"+name;
  cresph2->SaveAs(outname+".eps");
  cresph2->SaveAs(outname+".png");
  cresph2->SaveAs(outname+".pdf");
  cresph2->SaveAs(outname+".root");
  cresph2->SaveAs(outname+".C");



  ////// clsize 2 versions
  // version 90 cut
  Hist = "nrowB_clphABC90evR";

  MapTH1 nr_map;
  nr_map = GetCheckHists(&nr_map, angles, 1,dphcut, Runs,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,inputfileNONirr);

  TH1F  h_NrowB[angles];
  for(int i=0; i<angles; i++)
    {
      auto it2 = nr_map.find(std::make_pair(runs[i]+"_"+ss_dphcuts[i]+"_"+pitch,Hist));
      if(it2  != nr_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
	  cout << " Getting clusters " << endl;

	  h_NrowB[i] = * it2->second;
	  cout << h_NrowB[i].Integral() << endl;
	  cout << 1./h_NrowB[i].Integral() << endl;
	  Double_t norm = 1./h_NrowB[i].Integral();
	  if (h_NrowB[i].GetSumw2N() == 0) h_NrowB[i].Sumw2(kTRUE);
	  h_NrowB[i].Scale(norm,"width"); //1./h_NrowB[i].Integral());
	  cout << h_NrowB[i].Integral() << endl;
	  if (h_NrowB[i].GetSumw2N() == 0) h_NrowB[i].Sumw2(kTRUE);
	  h_NrowB[i].SetLineColor(color[i]);
	  h_NrowB[i].SetLineStyle(i+1);
	  h_NrowB[i].SetLineWidth(2);
	  h_NrowB[i].SetMarkerColor(color[i]);

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
  h_NrowB[0].GetXaxis()->SetTitle("Rows in cluster on DUT plane [pixels]");
  h_NrowB[0].GetXaxis()->SetRangeUser(0,10);
  h_NrowB[0].GetYaxis()->SetRangeUser(0.,1.); //RangeUser(0.,40000);
  h_NrowB[0].Draw("histe");
  cout << " Xaxis " << h_NrowB[0].GetXaxis()->GetNdivisions() << endl;
  cout << " RMS " << h_NrowB[0].GetRMS() << endl;
  cnB90->Update();
  //TPaveStats *
  Stats[0] =   (TPaveStats*)h_NrowB[0].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");

  Stats[0]->SetX1NDC(0.65);
  Stats[0]->SetX2NDC(.85);
  Stats[0]->SetY1NDC(.55);
  Stats[0]->SetY2NDC(.65);
  h_NrowB[1].GetYaxis()->SetRangeUser(0.,1.);
  h_NrowB[1].Draw("histesames");
  cnB90->Update();
  //TPaveStats *
  Stats[1] =   (TPaveStats*)h_NrowB[1].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");

  Stats[1]->SetX1NDC(0.65);
  Stats[1]->SetX2NDC(.85);
  Stats[1]->SetY1NDC(.44);
  Stats[1]->SetY2NDC(.54);
  Stats[1]->SetLineColor(color[1]);
  Stats[1]->SetLineStyle(2);
  gPad->Modified();



  h_NrowB[2].GetYaxis()->SetRangeUser(0.,1.);
  h_NrowB[2].Draw("histesames");
  cnB90->Update();
  //TPaveStats *
  Stats[2] =   (TPaveStats*)h_NrowB[2].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");
  //  Stats[2] =   (TPaveStats*)h_NrowB[2].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");

  Stats[2]->SetX1NDC(0.65);
  Stats[2]->SetX2NDC(.85);
  Stats[2]->SetY1NDC(.33);
  Stats[2]->SetY2NDC(.43);
  Stats[2]->SetLineColor(color[2]);
  Stats[2]->SetLineStyle(3);
  gPad->Modified();

  TLegend*legB2 = new TLegend(0.45,0.8,0.8,0.85);
  legB2->SetLineColor(0);
  //legFDB2->SetLegendFont(43);
  legB2->AddEntry(&(h_landau_orig[0]),base,"");
  legB2->Draw();

  
  
  TLegend* legnB90 = new TLegend(0.55,0.7,0.8,0.8);
  legnB90->SetLineColor(0);
  //  legnB90->SetTextSize(0.03);
  for(int i =0; i < angles; i++)
    legnB90->AddEntry(&(h_NrowB[i]),ang[i],"le");

  legnB90->Draw();
  outname = outputDir+"ClB90PHAll_NonIrr_3Angles_"+name;
  cnB90->SaveAs(outname+".eps");
  cnB90->SaveAs(outname+".png");
  cnB90->SaveAs(outname+".pdf");
  cnB90->SaveAs(outname+".root");
  cnB90->SaveAs(outname+".C");

  







  
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
