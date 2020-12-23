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

void controlPlotsNonIrr_forPaper(TString name = "preliminary_beamdiv")
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
  TString inputfileNONirr= "beamdiv_A13C14";

  GetHists(&landau_map,angles, 1, 1,dphcut, Runs, i_dphcut, Label, Hist, h_res,true,inputfileNONirr);
  
  int color[angles] = {39,1,803}; // at zero color was 920
      
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
  h_landau_orig[2].Scale(0.5);
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
  TH1F  h_resq_orig[angles];
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
  TPaveStats *Stats[angles];
  Float_t Mean[angles];
  Double_t RMS[angles];

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
  Double_t perc[angles];
  Double_t rmserr[angles];
  Double_t minX[angles]={-0.0406891*1000.,-0.0233376*1000.,-0.0436472*1000};
  Double_t maxX[angles]={0.0409722*1000.,0.024754*1000.,0.0474025*1000.};
  for( int i = 0 ; i < angles; i ++)  FitTH1(&(h_resq_orig[i]), &(RMS[i]), &(rmserr[i]),Runs[i]+"_", "A", "B", "C", "RMSself", &(perc[i]),&(minX[i]),&(maxX[i]));
  
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
  TLegend* legresq2 = new TLegend(0.13,0.73,0.7,0.82);//0.3,0.7,0.8,0.8);
  legresq2->SetLineColor(0);
  cout << "starting making hists pretty" << endl;
for( int i = 0 ; i < angles; i ++){
    cout << "            " << ang[i] << endl;

    cout << " min & max " << minX[i] << " " << maxX[i] << endl;
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
    legresq2->AddEntry(&(h_resq[i]),ang[i]+", #mu = "+s_m+" #mum, RMS = "+s_r+" #mum, tracks: "+s_p+" %","le");
    h_resq[i].GetXaxis()->SetRangeUser(-50.,50.);

  }

  legB->Draw();
  legresq2->Draw();

  for (int i = 0 ; i < angles; i ++){
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
    //cresph2->SetLogy();
     
  
  outname = outputDir+"ResPHAllnoLog_NonIrr_3Angles_"+name;
  cresph2->SaveAs(outname+".eps");
  cresph2->SaveAs(outname+".png");
  cresph2->SaveAs(outname+".pdf");
  cresph2->SaveAs(outname+".root");
  cresph2->SaveAs(outname+".C");

  // ---------------------
  //vertical incidence residual for different cluster size
  // getting resudual from tree
  TFile *f = new TFile(inputDir+"drei-r"+Runs[0]+"_irene_dphcutB"+ss_dphcut[0]+"_"+inputfileNONirr+".root");
  TTree *t1 = (TTree*)f->Get("charge_res");
  TH1F * hdx3treeph = new TH1F("hdx3treeph", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
  t1->Draw("dx3tree>>hdx3treeph","clphAiiitree<387&&clphBiiitree<302&&clphCiiitree<402","goff");
  hdx3treeph = (TH1F*)gDirectory->Get("hdx3treeph");
  hdx3treeph->Draw();
  ScaleAxis(hdx3treeph->GetXaxis(), ScaleX);


  TH1F * hdx1 = new TH1F("hdx1", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
  t1->Draw("dx3tree>>hdx1","clphAiiitree<387&&clphBiiitree<302&&clphCiiitree<402&&nrowBtree==1","goff");
  hdx1 = (TH1F*)gDirectory->Get("hdx1");
  hdx1->Draw();
  ScaleAxis(hdx1->GetXaxis(), ScaleX);
  TH1F * hdx2 = new TH1F("hdx2", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
  t1->Draw("dx3tree>>hdx2","clphAiiitree<387&&clphBiiitree<302&&clphCiiitree<402&&nrowBtree==2","goff");
  hdx2 = (TH1F*)gDirectory->Get("hdx2");
  hdx2->Draw();
  ScaleAxis(hdx2->GetXaxis(), ScaleX);
  TH1F * hdx3 = new TH1F("hdx3", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
  t1->Draw("dx3tree>>hdx3","clphAiiitree<387&&clphBiiitree<302&&clphCiiitree<402&&nrowBtree>=3","goff");
  hdx3 = (TH1F*)gDirectory->Get("hdx3");
  hdx3->Draw();
  ScaleAxis(hdx3->GetXaxis(), ScaleX);


  TCanvas *cresph0 = new TCanvas("cresph0", "FDB resolution", 600, 600);
  cresph0->SetLeftMargin(0.12);
  cresph0->SetRightMargin(-0.1);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTextFont(43);
  gStyle->SetTextSize(10);
  TLegend* legresq0 = new TLegend(0.13,0.73,0.7,0.82);//0.3,0.7,0.8,0.8);
  legresq0->SetLineColor(0);
  cout << "starting making hists pretty" << endl;
  hdx1->Scale(1./hdx3treeph->Integral());
  hdx2->Scale(1./hdx3treeph->Integral());
  hdx3->Scale(1./hdx3treeph->Integral());
  hdx3treeph->Scale(1./hdx3treeph->Integral());


  hdx1->GetXaxis()->SetRangeUser(-50,50);
  hdx2->GetXaxis()->SetRangeUser(-50,50);
  hdx3->GetXaxis()->SetRangeUser(-50,50);
  hdx3treeph->GetXaxis()->SetRangeUser(-50,50);

  hdx1->SetFillStyle(3005);
  hdx1->SetFillColor(8);
  hdx2->SetFillColor(4);
  hdx3->SetFillColor(2);
  hdx2->SetFillStyle(3004);
  hdx3->SetFillStyle(3001);
  THStack *hs = new THStack("hs","");
  hs->Add(hdx3);
  hs->Add(hdx2);
  hs->Add(hdx1);
  h_resq[0].GetYaxis()->SetRangeUser(0,0.25);
  h_resq[0].Draw("histe");
  cresph->Update();

  //hdx3treeph->SetLineStyle(2);
  //hdx3treeph->Draw("histsame");
  hs->Draw("histsame");
  cout << " hdx3treeph " << hdx3treeph->GetMean() << endl;
  TLegend* legB0 = new TLegend(0.15,0.83,0.6,0.88);
  legB0->SetLineColor(0);
  legB0->AddEntry(&(h_landau_orig[0]),base+", "+ang[0],"");
  legB0->Draw();

  legresq0->AddEntry(hdx1,"Cluster size = 1","f");
  legresq0->AddEntry(hdx2,"Cluster size = 2","f");
  legresq0->AddEntry(hdx3,"Cluster size #geq 3","f");
  legresq0->Draw();
  
  outname = outputDir+"ResPHAllnoLog_NonIrr_VerticalCLsize_"+name;
  cresph0->SaveAs(outname+".eps");
  cresph0->SaveAs(outname+".png");
  cresph0->SaveAs(outname+".pdf");
  cresph0->SaveAs(outname+".root");
  cresph0->SaveAs(outname+".C");


  ////// clsize distribution
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
  Stats[0] =   (TPaveStats*)h_NrowB[0].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");
  Stats[0]->SetX1NDC(0.65);
  Stats[0]->SetX2NDC(.85);
  Stats[0]->SetY1NDC(.55);
  Stats[0]->SetY2NDC(.65);

  h_NrowB[1].GetYaxis()->SetRangeUser(0.,1.);
  h_NrowB[1].Draw("histesames");
  cnB90->Update();
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
  Stats[2] =   (TPaveStats*)h_NrowB[2].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");
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




  // add also no threshold
  MapTH1 nr_map0;
  TString thr="1";
  TString  ss_dphcut0[1]={thr};
  TString  ss_dphcuts0[angles]={thr,thr,thr};
  nr_map0 = GetCheckHists(&nr_map0, angles, 1,dphcut, Runs,ss_dphcut0,pitch, Hist,hdx3_clchargeABC90evR,true,inputfileNONirr);

  TH1F  h_NrowB0[angles];
  for(int i=0; i<angles; i++)
    {
      auto it2 = nr_map0.find(std::make_pair(runs[i]+"_"+ss_dphcuts0[i]+"_"+pitch,Hist));
      if(it2  != nr_map0.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
	  cout << " Getting clusters " << endl;

	  h_NrowB0[i] = * it2->second;
	  cout << h_NrowB0[i].Integral() << endl;
	  cout << 1./h_NrowB0[i].Integral() << endl;
	  Double_t norm0 = 1./h_NrowB0[i].Integral();
	  if (h_NrowB0[i].GetSumw2N() == 0) h_NrowB0[i].Sumw2(kTRUE);
	  h_NrowB0[i].Scale(norm0,"width"); //1./h_NrowB[i].Integral());
	  cout << h_NrowB0[i].Integral() << endl;
	  if (h_NrowB0[i].GetSumw2N() == 0) h_NrowB0[i].Sumw2(kTRUE);
	  h_NrowB0[i].SetLineColor(color[i]);
	  h_NrowB0[i].SetLineStyle(i+1);
	  h_NrowB0[i].SetLineWidth(2);
	  h_NrowB0[i].SetMarkerColor(color[i]);

	}

    }
  TCanvas *cnB90_0 = new TCanvas("cnB90_0", "FDB resolution", 600, 600);
  cnB90_0->SetLeftMargin(0.12);
  cnB90_0->SetRightMargin(-0.1);
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
  h_NrowB[0].GetXaxis()->SetRangeUser(0,30);
  h_NrowB[0].GetYaxis()->SetRangeUser(0.,1.); //RangeUser(0.,40000);
  h_NrowB[0].Draw("histe");
  cout << " Xaxis " << h_NrowB[0].GetXaxis()->GetNdivisions() << endl;
  cout << " RMS " << h_NrowB[0].GetRMS() << endl;
  cnB90_0->Update();
  Stats[0] =   (TPaveStats*)h_NrowB[0].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");
  Stats[0]->SetX1NDC(0.65);
  Stats[0]->SetX2NDC(.85);
  Stats[0]->SetY1NDC(.55);
  Stats[0]->SetY2NDC(.65);

  h_NrowB[1].GetYaxis()->SetRangeUser(0.,1.);
  h_NrowB[1].Draw("histesames");
  cnB90_0->Update();
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
  cnB90_0->Update();
  Stats[2] =   (TPaveStats*)h_NrowB[2].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");
  Stats[2]->SetX1NDC(0.65);
  Stats[2]->SetX2NDC(.85);
  Stats[2]->SetY1NDC(.33);
  Stats[2]->SetY2NDC(.43);
  Stats[2]->SetLineColor(color[2]);
  Stats[2]->SetLineStyle(3);
  gPad->Modified();

  ///// now thr0
  TPaveStats *Stats0[angles];
  h_NrowB0[0].SetTitle(" ");
  h_NrowB0[0].GetYaxis()->SetTitle("Normalized number of events");
  h_NrowB0[0].GetXaxis()->SetTitle("Rows in cluster on DUT plane [pixels]");
  h_NrowB0[0].GetXaxis()->SetRangeUser(0,10);
  h_NrowB0[0].GetYaxis()->SetRangeUser(0.,1.); //RangeUser(0.,40000);
  h_NrowB0[0].Draw("*histesames");
  cout << " Xaxis " << h_NrowB[0].GetXaxis()->GetNdivisions() << endl;
  cout << " RMS " << h_NrowB[0].GetRMS() << endl;
  cnB90_0->Update();
  Stats0[0] =   (TPaveStats*)h_NrowB0[0].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");
  Stats0[0]->SetX1NDC(0.45);
  Stats0[0]->SetX2NDC(.65);
  Stats0[0]->SetY1NDC(.55);
  Stats0[0]->SetY2NDC(.65);

  h_NrowB0[1].GetYaxis()->SetRangeUser(0.,1.);
  h_NrowB0[1].Draw("*histesames");
  cnB90_0->Update();
  Stats0[1] =   (TPaveStats*)h_NrowB0[1].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");
  Stats0[1]->SetX1NDC(0.45);
  Stats0[1]->SetX2NDC(.65);
  Stats0[1]->SetY1NDC(.44);
  Stats0[1]->SetY2NDC(.54);
  Stats0[1]->SetLineColor(color[1]);
  Stats0[1]->SetLineStyle(2);
  gPad->Modified();

  h_NrowB0[2].GetYaxis()->SetRangeUser(0.,1.);
  h_NrowB0[2].Draw("*histesames");
  cnB90_0->Update();
  Stats0[2] =   (TPaveStats*)h_NrowB0[2].FindObject("stats");     //GetListOfFunctions()->FindObject("stats");
  Stats0[2]->SetX1NDC(0.45);
  Stats0[2]->SetX2NDC(.65);
  Stats0[2]->SetY1NDC(.33);
  Stats0[2]->SetY2NDC(.43);
  Stats0[2]->SetLineColor(color[2]);
  Stats0[2]->SetLineStyle(3);
  gPad->Modified();

  legB2->Draw();
  legnB90->Draw();

  outname = outputDir+"ClB90PHAll_NonIrr_THR"+thr+"_3Angles_"+name;
  cnB90_0->SaveAs(outname+".eps");
  cnB90_0->SaveAs(outname+".png");
  cnB90_0->SaveAs(outname+".pdf");
  cnB90_0->SaveAs(outname+".root");
  cnB90_0->SaveAs(outname+".C");


  // ---------------------
  //pixel charge
  // getting resudual from tree
  TCanvas *cpx = new TCanvas("cpx", "FDB resolution", 600, 600);
  cpx->SetLeftMargin(0.12);
  cpx->SetRightMargin(-0.1);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTextFont(43);
  gStyle->SetTextSize(10);
  TFile *f1 = new TFile(inputDir+"drei-r"+Runs[1]+"_irene_dphcutB"+ss_dphcut[0]+"_"+inputfileNONirr+".root");
  TTree *t = (TTree*)f1->Get("charge_res");
  TH1F * hpxphB = new TH1F( "hpxphB", "B pixel PH;pixel PH [ADC];B pixels on tracks", 250, 0, 500 );
  t->Draw("pxphBtree>>hpxphB","1","goff");
  hpxphB = (TH1F*)gDirectory->Get("hpxphB");
  TH1F * hpxphB90 = new TH1F( "hpxphB90", "B pixel PH;pixel PH [ADC];B pixels on tracks", 250, 0, 500 );
  t->Draw("pxphBtree>>hpxphB90","clphAiiitree<397&&clphBiiitree<317&&clphCiiitree<417","goff");
  hpxphB90 = (TH1F*)gDirectory->Get("hpxphB90");

  TLegend* legresq90 = new TLegend(0.2,0.73,0.7,0.82);//0.3,0.7,0.8,0.8);
  legresq90->SetLineColor(0);
  cout << "starting making hists pretty" << endl;
  hpxphB->GetYaxis()->SetTitle("Pixels above threshold");
  hpxphB->GetXaxis()->SetTitle("Pixel Pulse Height [ADC]");
  hpxphB->GetYaxis()->SetTitleFont(43);
  hpxphB->GetYaxis()->SetTitleSize(20);
  hpxphB->GetYaxis()->SetTitleOffset(1.8);
  hpxphB->GetXaxis()->SetTitleFont(43);
  hpxphB->GetXaxis()->SetTitleSize(20);
  hpxphB->GetYaxis()->SetLabelFont(43);
  hpxphB->GetYaxis()->SetLabelSize(20);
  hpxphB->GetXaxis()->SetLabelFont(43);
  hpxphB->GetXaxis()->SetLabelSize(20);
  hpxphB->Rebin(4);
  hpxphB90->Rebin(4);
  hpxphB->SetLineWidth(2);
  hpxphB90->SetLineWidth(2);
  hpxphB->Draw("histe");
  hpxphB90->SetLineStyle(2);
  hpxphB90->Draw("histesame");
  cresph->Update();

  //hdx3treeph->SetLineStyle(2);
  //hdx3treeph->Draw("histsame");
  TLegend* legB90 = new TLegend(0.2,0.83,0.6,0.88);
  legB90->SetLineColor(0);
  legB90->AddEntry(&(h_landau_orig[0]),base+", "+ang[1],"");
  legB90->Draw();

  legresq90->AddEntry(hpxphB,"All","le");
  legresq90->AddEntry(hpxphB90,"90% events with lowest cluster charge","le");
  legresq90->Draw();
  
  outname = outputDir+"PxPH_NonIrr_BestAngle_"+name;
  cpx->SaveAs(outname+".eps");
  cpx->SaveAs(outname+".png");
  cpx->SaveAs(outname+".pdf");
  cpx->SaveAs(outname+".root");
  cpx->SaveAs(outname+".C");



  
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
