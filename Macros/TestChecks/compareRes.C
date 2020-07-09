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

#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "fileHandler.h"

bool print=true;
bool debug=false;
using namespace std;
#define Events 90000
#define Cuts 2
#define Hists 5
#define runs 1
#define dphcuts 1

void TDR();
void TDR2(TCanvas * c_all, int period=0, int pos= 11);
typedef std::map<int, std::pair<double, double>> MapEvtRes; // evt number, pair (dxyCA, dxy)
void mapfilling(MapEvtRes * map, int evt, double dxyCA,  double dxy);


		
void compareRes(TString name = "PIrrABC_1clusterClosest"){
  
  TString inputDir="/home/zoiirene/Output/TextFiles/";
  TString outputDir="/home/zoiirene/Output/Plots/";
  
  TString filenames[Cuts];
  TString interval[Cuts]    ={"1clusterABC","closest_test2"};
  for(int i = 0; i < Cuts; i++){
    filenames[i] = interval[i]+".txt";
  }
  
  TString irr =  " Non-irradiated, 5.6 GeV"; // "no irr, 5.6 GeV";


  MapEvtRes cluster1,closest;  
  TString filename;

  TH1I * h_event[Cuts];
  TH1I * h_eventdiff[Cuts];// = new TH1I("ev","ev",Events+2,200,90200);  
  TH1D * h_dxyCA[Cuts];
  TH1D * h_dxy[Cuts];
  
  int evt;
  double dxyCA,dxy;
  string dummy;
  for(int i =0; i< Cuts; i++){
    h_event[i] =  new TH1I("ev","ev",Events+2,200,90200);
    h_eventdiff[i] =  new TH1I("ev","ev",Events+2,200,90200);
    h_dxyCA[i]  =  new TH1D("dxyCA","dxyCA",100,0.,200.);
    h_dxy[i]  =  new TH1D("dxy","dxy",100,0.,200.);

    filename      = inputDir+filenames[i];
    cout << filename << endl;
    ifstream stream(filename);
      
    std::string line;
    if(!stream.is_open()){
      cout << " File " << filename[i] << " not opened" << endl;
    }else{

      while(stream >> dummy >> evt >> dummy >> dummy >> dummy >> dxyCA >> dummy >> dxy){

	if(debug)	 cout << "evt " << evt << " dxyCA " << dxyCA << " dxy " << dxy << endl;
	if(i==0) mapfilling(&(cluster1),evt, dxyCA*1000, dxy*1000);
	if(i==1) mapfilling(&(closest), evt, dxyCA*1000, dxy*1000);
      } //while
    }//file scanning
  }//cuts 


   
  int j = 0;
  int k =0;
  for(int i =0; i<= Events; i++){

    int evtnum = i+200;
    if(debug) cout << " Event " << evtnum << endl;
    auto  it = cluster1.find(evtnum);
    if(it !=cluster1.end()){
      if(debug) cout << "cluster1 key " << it->first <<endl;
      
      if(debug)      cout <<" content " << it->second.first << " " << it->second.second << endl;
      h_event[0]->Fill(evtnum,1);
      h_dxyCA[0]->Fill(it->second.first);
      h_dxy[0]->Fill(it->second.second);
      
    }// scan cluster 1
    
    auto  it2 = closest.find(evtnum);
    if(it2 !=closest.end()){
      if(debug)      cout << "closest key " << it2->first <<endl;
      
      if(debug)      cout <<" content " << it2->second.first << " " << it2->second.second << endl;
      h_event[1]->Fill(evtnum,1);
      h_dxyCA[1]->Fill(it2->second.first);
      h_dxy[1]->Fill(it2->second.second);
      
      
    }//scan closest
    if(it !=cluster1.end() && it2 ==closest.end()){
      h_eventdiff[0]->Fill(evtnum,1);
      if(debug) cout << " event " << evtnum << " only in 1cluster" << endl;
      k++;
    }
      if(it ==cluster1.end() && it2 !=closest.end())	{
	  h_eventdiff[1]->Fill(evtnum,1);
	  if(debug) cout << " event "<< evtnum << " only in closest" << endl;
	  j++;
      }
  }// events


  cout << " events only in  1cluster " << k << endl;
  cout << " events only in  closest " << j << endl;
    


 
  /////// plots!

  
  TCanvas *cFDB2 = new TCanvas("cFDB2", "FDB resolution", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  for(int i =0; i< Cuts; i++){
  
    h_event[i]->SetTitle(" ");
    h_event[i]->GetYaxis()->SetTitle("found");
    h_event[i]->GetXaxis()->SetTitle("event");
    //h_event[i]->SetMarkerSize(2.5);
    
    h_event[i]->SetLineColor(i+1);

    h_event[i]->SetLineStyle(i+1);
      
  
    h_event[i]->GetXaxis()->SetLimits(200,90200);
    h_event[i]->GetYaxis()->SetRangeUser(0.,2.);
    
    if(i==0)    h_event[i]->Draw("hist");
    else h_event[i]->Draw("histsame");
  }


  //  if(logar) cFDB2->SetLogy();   
  TLegend* legFDB2 = new TLegend(0.2,0.7,0.7,0.85);
  legFDB2->SetLineColor(0);
  legFDB2->SetTextSize(0.035);
  legFDB2->SetNColumns(2);
  for(int i =0; i < Cuts; i++){
    if(i==0) legFDB2->AddEntry(h_event[i],"1cluster","l");
    else legFDB2->AddEntry(h_event[i],"closest","l");
  }
  legFDB2->Draw();

  
  TString  outname = outputDir+"events_bestAngle_"+name;
  
  cFDB2->SaveAs(outname+".eps");
  cFDB2->SaveAs(outname+".png");
  cFDB2->SaveAs(outname+".pdf");
  cFDB2->SaveAs(outname+".root");
  cFDB2->SaveAs(outname+".C");



  TCanvas *evdiff[Cuts];
  for(int i =0; i < Cuts; i++){
    evdiff[i]      = new TCanvas("evdiff_"+interval[i], "FDB resolution", 1500, 900);
    gPad->SetTicks(1,1);
    gROOT->SetStyle("Plain");
    gStyle->SetPadGridX(0);
    gStyle->SetPadGridY(0);
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    h_eventdiff[i]->SetTitle(" ");
    h_eventdiff[i]->GetYaxis()->SetTitle("found");
    h_eventdiff[i]->GetXaxis()->SetTitle("event");
    //h_eventdiff->SetMarkerSize(2.5);
    
    h_eventdiff[i]->SetLineColor(kBlue);
    
    h_eventdiff[i]->GetXaxis()->SetLimits(200,90200);
    h_eventdiff[i]->GetYaxis()->SetRangeUser(0.,2.);
    
    h_eventdiff[i]->Draw("hist");
    
    outname = outputDir+"eventfoundonly_"+interval[i]+"_bestAngle_"+name;
    
    evdiff[i]->SaveAs(outname+".eps");
    evdiff[i]->SaveAs(outname+".png");
    evdiff[i]->SaveAs(outname+".pdf");
    evdiff[i]->SaveAs(outname+".root");
    evdiff[i]->SaveAs(outname+".C");
  }


  TCanvas *cdxyCA = new TCanvas("dxyCA", "FDB resolution", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  for(int i =0; i< Cuts; i++){

    h_dxyCA[i]->SetTitle(" ");
    h_dxyCA[i]->GetYaxis()->SetTitle("found");
    h_dxyCA[i]->GetXaxis()->SetTitle("dxyCA [#mum]");
    //h_dxyCA[i]->SetMarkerSize(2.5);

    h_dxyCA[i]->SetLineColor(i+1);

    h_dxyCA[i]->SetLineStyle(i+1);


    h_dxyCA[i]->GetXaxis()->SetLimits(0.,200.);
    h_dxyCA[i]->GetYaxis()->SetRangeUser(0.,4000.);

    if(i==0)    h_dxyCA[i]->Draw("hist");
    else h_dxyCA[i]->Draw("histsame");
  }


  //  if(logar) dxyCA->SetLogy();
  TLegend* legdxyCA = new TLegend(0.2,0.7,0.7,0.85);
  legdxyCA->SetLineColor(0);
  legdxyCA->SetTextSize(0.035);
  legdxyCA->SetNColumns(2);
  for(int i =0; i < Cuts; i++){
    if(i==0) legdxyCA->AddEntry(h_dxyCA[i],"1cluster","l");
    else legdxyCA->AddEntry(h_dxyCA[i],"closest","l");
  }
  legdxyCA->Draw();


   outname = outputDir+"dxyCA_bestAngle_"+name;

  cdxyCA->SaveAs(outname+".eps");
  cdxyCA->SaveAs(outname+".png");
  cdxyCA->SaveAs(outname+".pdf");
  cdxyCA->SaveAs(outname+".root");
  cdxyCA->SaveAs(outname+".C");






  TCanvas *cdxy = new TCanvas("cdxy", "FDB resolution", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  for(int i =0; i< Cuts; i++){

    h_dxy[i]->SetTitle(" ");
    h_dxy[i]->GetYaxis()->SetTitle("found");
    h_dxy[i]->GetXaxis()->SetTitle("dxy [#mum]");
    //h_dxy[i]->SetMarkerSize(2.5);

    h_dxy[i]->SetLineColor(i+1);

    h_dxy[i]->SetLineStyle(i+1);

    h_dxy[i]->GetXaxis()->SetLimits(0.,200.);
    h_dxy[i]->GetYaxis()->SetRangeUser(0.,10000.);
    

    if(i==0)    h_dxy[i]->Draw("hist");
    else h_dxy[i]->Draw("histsame");
  }


  //  if(logar) dxy->SetLogy();
  TLegend* legdxy = new TLegend(0.2,0.7,0.7,0.85);
  legdxy->SetLineColor(0);
  legdxy->SetTextSize(0.035);
  legdxy->SetNColumns(2);
  for(int i =0; i < Cuts; i++){
    if(i==0) legdxy->AddEntry(h_dxy[i],"1cluster","l");
    else legdxy->AddEntry(h_dxy[i],"closest","l");
  }
  legdxy->Draw();


   outname = outputDir+"dxy_bestAngle_"+name;

  cdxy->SaveAs(outname+".eps");
  cdxy->SaveAs(outname+".png");
  cdxy->SaveAs(outname+".pdf");
  cdxy->SaveAs(outname+".root");
  cdxy->SaveAs(outname+".C");
  

  ////// compare hists after cuts

  MapTH1 res_map[Cuts];
  int i_dphcut[dphcuts] = {22};
  TString ss_dphcut[dphcuts] = {"22"};
  bool dphcut=true;

  TString Run[runs]={"2743"};
  TString pitch = "25";
  TString Hist[Hists] = {"dx3_clphABC90evR99","dxyCA_ph99","dxy_ph99","hdx3tree2","dx3_clphABC90evR"};
  for(int i = 0; i < Cuts; i++){
    //filenames[i] = "drei-r2743_irene_dphcutB22_"+interval[i]+".root";
    std::map<std::pair<TString, TString>, TH1F *>::iterator it;
    TH1F * hdx3_clchargeABC90evR; // = new TH1I("dx3_clchargeABC90evR ", "triplet dx_clchargeABC90evR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 90% events in Landau (only high tail)
    for(int j =0; j < Hists; j++){

      res_map[i] =    GetCheckHists(&(res_map[i]), runs, dphcuts,dphcut, Run,ss_dphcut,pitch, Hist[j],hdx3_clchargeABC90evR,true,interval[i]);
    }
  
    for ( it = res_map[i].begin(); it != res_map[i].end(); it++ )       {
      cout << it->first.first << " " << it->first.second << " " <<  it->second->GetEntries() << endl;
    }
  }

  float xminl[Hists] ={-0.2,0.,0.,-1,-0.5};
  float xmin[Hists] ={-0.15,0.,0.,-0.15,-0.15};
  float xmaxl[Hists] ={0.2,0.3,0.15,1,0.5};
  float xmax[Hists] ={0.15,0.3,0.15,0.15,0.15};
  float yminl[Hists] = {0.1,0.1,0.1,0.1,0.1};
  float ymin[Hists] = {0.,0.,0.,0.,0.};
  float ymaxl[Hists] = {10000,10000,20000,10000,10000};
  float ymax[Hists] = {10000,3000,12000,10000,10000};
  for(int j =0; j < Hists; j++){
    TString key = Run[0]+"_"+ss_dphcut[0]+"_"+pitch;
    auto it1 = res_map[0].find(std::make_pair(key,Hist[j]));
    TH1F * h1 = it1->second;
    cout << " h1 " << it1->second->GetEntries() << endl;
    auto it2 = res_map[1].find(std::make_pair(key,Hist[j]));
    TH1F * h2 = it2->second;
    cout << " h2 " << it2->second->GetEntries()<< endl;
    
    DrawTwoHist(h1,h2, interval[0], interval[1],Run[0]+"_"+Hist[j]+"_log",xminl[j],xmaxl[j],yminl[j], ymaxl[j]);    
    DrawTwoHist(h1,h2, interval[0], interval[1],Run[0]+"_"+Hist[j],xmin[j],xmax[j],ymin[j], ymax[j],false);    
  }



}//resolution 




void mapfilling(MapEvtRes * map, int evt, double dxyCA,  double dxy){
  
  auto it = map->find(evt);


  if(it  != map->end())
    {
      if(debug) cout << " it  != map->end() " << endl;
      if(debug)           cout << " found map " << endl;
      if(debug)  cout << "map key " << it->first <<endl;

      if(debug) cout <<" content " << it->second.first << " " << it->second.second << endl;

    }
  else
    {
      if(debug) cout << " it  == map->end() " << endl;
      map->insert(std::make_pair(evt,std::make_pair(dxyCA,dxy)));
      it = map->find(evt);

      if(debug) cout << "map key " << it->first <<endl;

      if(debug) cout <<" content " << it->second.first << " " << it->second.second << endl;
      
    }
  


}//mapfilling







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
