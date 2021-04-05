#include "resolution.h"
typedef std::map<std::pair<TString, TString>, TH1F*> MapTH1;

void DrawMoreTGraphWithError(int comparisons, int *i_dphcut, float * mean, float * rms,float * mean1, float * rms1,float * mean2, float * rms2, TString Hist,TString Run, TString * Label, TString Name, TString yaxistitle, float xmin,float xmax,float ymin, float ymax, TString scan = "thresholds", TString xaxistitle =  "ph cut [ADC]")
{

  TString outputDir = "/home/zoiirene/Output/Plots/";
  TString ss_dphcut[comparisons];
  float f_dphcut[comparisons];
  float errx[comparisons];
  for(int i =0; i< comparisons;i++)
    {
      ss_dphcut[i].Form("%d",i_dphcut[i]);
      f_dphcut[i] = (float) i_dphcut[i];
      errx[i]=0;
    }

  TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
  c2->SetLeftMargin(0.15);
  c2->SetRightMargin(0.05);
  c2->SetTopMargin(0.05);

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
  gStyle->SetLegendTextSize(24);

  TGraphErrors * gr = new TGraphErrors(comparisons, f_dphcut, mean,errx,rms);//, errx, rms);
  TGraphErrors * gr1 = new TGraphErrors(comparisons, f_dphcut, mean1,errx,rms1);//, errx, rms);
  TGraphErrors * gr2 = new TGraphErrors(comparisons, f_dphcut, mean2,errx,rms2);//, errx, rms);

  gr->SetLineColor(kAzure);
  gr->SetLineWidth(2);
  gr->SetMarkerColor(kAzure);
  gr->SetMarkerSize(1.5);
  gr->SetMarkerStyle(21);
  gr->GetXaxis()->SetRangeUser(xmin,xmax);
  gr->SetMinimum(ymin);
  gr->SetMaximum(ymax);

  gr2->SetLineColor(kBlue-8);
  gr2->SetLineWidth(2);
  gr2->SetLineStyle(3);
  gr2->SetMarkerColor(kBlue-8);
  gr2->SetMarkerSize(1.5);
  gr2->SetMarkerStyle(25);

  gr1->SetLineColor(kBlue+4);
  gr1->SetLineWidth(2);
  gr1->SetLineStyle(2);
  gr1->SetMarkerColor(kBlue+4);
  gr1->SetMarkerSize(1.5);
  gr1->SetMarkerStyle(31);


  gr->GetYaxis()->SetMaxDigits(1);
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetTitleOffset(0.9);
  gr->GetYaxis()->SetTitleOffset(1.4);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetLabelSize(0.05);

  gr->GetXaxis()->SetTitle(xaxistitle);
  gr->GetYaxis()->SetTitle(yaxistitle);
  gr->GetXaxis()->SetLimits(xmin,xmax);
  gr->SetMinimum(ymin);
  gr->SetMaximum(ymax);

  gr->Draw("ALPE");
  gr1->Draw("LPEsame");
  gr2->Draw("LPEsame");
  TLegend* leg2;
  if (Run  == "2801" && Hist == "dx3vsx"){
    leg2 = new TLegend(0.2,0.65,0.4,0.8);
  } else if (Run  == "2801" && (Hist == "alignfA" || Hist == "alignxC" )){
    leg2 = new TLegend(0.7,0.65,0.9,0.8);
  }else{
    leg2 = new TLegend(0.7,0.15,0.9,0.3);
  }

  leg2->SetLineColor(0);
  leg2->AddEntry(gr,Label[0], "pl");
  leg2->AddEntry(gr1,Label[1], "pl");
  leg2->AddEntry(gr2,Label[2], "pl");

  leg2->Draw("same");
  TString name;
  name = outputDir+"compare_"+scan+"_"+Run+"_"+Hist+"_"+Name;
  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".pdf");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".root");
}


void DrawTGraphWithError(int comparisons, int *i_dphcut, float * mean, float * rms, TString Hist,TString Run, TString Label, TString Name, TString yaxistitle, float xmin,float xmax,float ymin, float ymax, TString scan = "thresholds", TString xaxistitle =  "ph cut [ADC]")
{

  TString outputDir = "/home/zoiirene/Output/Plots/";
  TString ss_dphcut[comparisons];
  float f_dphcut[comparisons];
  float errx[comparisons];
  for(int i =0; i< comparisons;i++)
    {
      ss_dphcut[i].Form("%d",i_dphcut[i]);
      f_dphcut[i] = (float) i_dphcut[i];
      errx[i]=0;
    }

  TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
  c2->SetLeftMargin(0.15);
  c2->SetRightMargin(0.05);
  c2->SetTopMargin(0.05);

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
  gStyle->SetLegendTextSize(24);

  TGraphErrors * gr = new TGraphErrors(comparisons, f_dphcut, mean,errx,rms);//, errx, rms);
  gr->SetLineColor(kAzure);
  gr->SetLineWidth(2);
  gr->SetMarkerColor(kAzure);
  gr->SetMarkerSize(1.5);
  gr->SetMarkerStyle(21);
  gr->GetXaxis()->SetRangeUser(xmin,xmax);
  gr->SetMinimum(ymin);
  gr->SetMaximum(ymax);


  gr->GetYaxis()->SetMaxDigits(1);
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetTitleOffset(0.9);
  gr->GetYaxis()->SetTitleOffset(1.4);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetLabelSize(0.05);

  gr->GetXaxis()->SetTitle(xaxistitle);
  gr->GetYaxis()->SetTitle(yaxistitle);
  gr->GetXaxis()->SetLimits(xmin,xmax);
  gr->SetMinimum(ymin);
  gr->SetMaximum(ymax);

  gr->Draw("ALPE");


  TString name;
  name = outputDir+"compare_"+scan+"_"+Run+"_"+Label+"_"+Hist+"_"+Name;
  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".pdf");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".root");
}


void DrawTGraphWithErrorDouble(int comparisons, double * i_dphcut, double * mean, double * rms, TString Hist,TString Run, TString Label, TString Name, TString yaxistitle, float xmin,float xmax,float ymin, float ymax, TString scan = "thresholds", TString xaxistitle =  "ph cut [ADC]")
{

  TString outputDir = "/home/zoiirene/Output/Plots/";
  TString ss_dphcut[comparisons];

  double errx[comparisons];
  for(int i =0; i< comparisons;i++)
    {
      ss_dphcut[i].Form("%f",i_dphcut[i]);
      
      errx[i]=0;
    }
  
  TCanvas *c2 = new TCanvas("c2", "c2", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  //  gStyle->SetOptStat();
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);

  TGraphErrors * gr = new TGraphErrors(comparisons, i_dphcut, mean,errx,rms);//, errx, rms);
  gr->SetLineColor(kAzure);
  gr->SetLineWidth(2);
  gr->SetMarkerColor(kAzure);
  gr->SetMarkerSize(1.5);
  gr->SetMarkerStyle(21);
  gr->GetXaxis()->SetRangeUser(xmin,xmax);
  gr->SetMinimum(ymin);
  gr->SetMaximum(ymax);


  
  //  gr->SetTitle("Option ACP example");
  gr->GetXaxis()->SetTitle(xaxistitle);
  gr->GetYaxis()->SetTitle(yaxistitle);
  //gr->SetFillStyle(3003);
  //gr->SetFillColor(kRed-8);


  gr->GetXaxis()->SetRangeUser(xmin,xmax);
  gr->SetMinimum(ymin);
  gr->SetMaximum(ymax);
  
  gr->Draw("ALPE");


  TString name;
  name = outputDir+"compare_"+scan+"_"+Run+"_"+Label+"_"+Hist+"_"+Name;
  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".pdf");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".root");
}



void DrawTGraph(int comparisons, int * i_dphcut, int * entries, TString Hist,TString Run, TString Label, TString Name, int xmin,int xmax,int ymin, int ymax)
{

  TString outputDir = "/home/zoiirene/Output/Plots/";
  
  TCanvas *c2 = new TCanvas("c2", "c2", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  //  gStyle->SetOptStat();
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);

  TGraph * gr = new TGraph(comparisons, i_dphcut, entries);
  gr->SetLineColor(2);
  gr->SetLineWidth(2);
  gr->SetMarkerColor(2);
  gr->SetMarkerSize(1.5);
  gr->SetMarkerStyle(21);
  //  gr->SetTitle("Option ACP example");
  gr->GetXaxis()->SetTitle("ph cut [ADC]");
  gr->GetXaxis()->SetRangeUser(xmin,xmax);
  gr->GetYaxis()->SetTitle("# tracks");
  gr->SetMinimum(ymin);
  gr->SetMaximum(ymax);
  gr->Draw("ALP");
  
  
  TString name;
  name = outputDir+"compare_thresholds_"+Run+"_"+Label+"_"+Hist+"_"+Name;
  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".pdf");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".root");
}

void DrawTGraphError(int comparisons, int * i_dphcut, float * mean, float * rms, TString Hist,TString Run, TString Label, TString Name, TString yaxistitle, float xmin,float xmax,float ymin, float ymax, TString scan = "thresholds", TString xaxistitle =  "ph cut [ADC]")
{

  TString outputDir = "/home/zoiirene/Output/Plots/";
  TString ss_dphcut[comparisons];
  float f_dphcut[comparisons];
  float errx[comparisons];
  for(int i =0; i< comparisons;i++)
    {
      ss_dphcut[i].Form("%d",i_dphcut[i]);
      f_dphcut[i] = (float) i_dphcut[i];
      errx[i]=0;
    }
  
  TCanvas *c2 = new TCanvas("c2", "c2", 1500, 900);
  c2->SetLeftMargin(0.1);
  c2->SetTopMargin(0.05);
  c2->SetBottomMargin(0.2);
  
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  //  gStyle->SetOptStat();
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);

  TGraph *  grgreen = new TGraph(2*comparisons);

  TGraph * gr = new TGraph(comparisons, f_dphcut, mean);//, errx, rms);
  gr->SetLineColor(kAzure);
  gr->SetLineWidth(2);
  gr->SetMarkerColor(kAzure);
  gr->SetMarkerSize(1.5);
  gr->SetMarkerStyle(21);

  grgreen->GetXaxis()->SetLabelFont(42);
  grgreen->GetXaxis()->SetLabelSize(0.07);
  grgreen->GetXaxis()->SetTitleSize(0.07);
  grgreen->GetXaxis()->SetTitleOffset(0.9);
  grgreen->GetXaxis()->SetTitleFont(42);

  grgreen->GetYaxis()->SetLabelFont(42);
  grgreen->GetYaxis()->SetLabelSize(0.07);
  grgreen->GetYaxis()->SetTitleSize(0.07);
  grgreen->GetYaxis()->SetTitleOffset(0.6);
  grgreen->GetYaxis()->SetTitleFont(42);
  grgreen->GetYaxis()->SetNoExponent(3);
  

  grgreen->GetXaxis()->SetRangeUser(xmin,xmax);
  grgreen->SetMinimum(ymin);
  grgreen->SetMaximum(ymax);
  
  //  gr->SetTitle("Option ACP example");
  grgreen->GetXaxis()->SetTitle(xaxistitle);
  grgreen->GetYaxis()->SetTitle(yaxistitle);
  //gr->SetFillStyle(3003);
  //gr->SetFillColor(kRed-8);
  for(int i= 0; i<comparisons; i++)
    {
      grgreen->SetPoint(i,f_dphcut[i],mean[i]+rms[i]);
      grgreen->SetPoint(comparisons+i,f_dphcut[comparisons-i-1],mean[comparisons-i-1]-rms[comparisons-i-1]);
    }
  grgreen->SetFillStyle(1001);
  grgreen->SetFillColor(kAzure-9);

  grgreen->Draw("AF");
  gr->GetXaxis()->SetRangeUser(xmin,xmax);
  gr->SetMinimum(ymin);
  gr->SetMaximum(ymax);
  
  gr->Draw("LPsame");


  TString name;
  name = outputDir+"compare_"+scan+"_"+Run+"_"+Label+"_"+Hist+"_"+Name;
  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".pdf");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".root");
}


TH1F* GetMoreCheckHists(MapTH1 * map,int runs, int dphcuts, bool dphcut, TString * Run, TString * ss_dphcut, TString pitch, TString * Hist,int hists , TH1F * h, bool extralabel = false, TString label = " ")
{
  TString inputDir="/home/zoiirene/Output/";
  TString inputfile, Path;

  TFile * file[runs][dphcuts];
  bool print = true;
  TH1F * g;
  TString key;
  if(print) cout << "Getting runs "  << endl;
  for(int i=0; i<runs; i++)
    {
      for(int l=0; l<dphcuts; l++)
	{

	  if(print) cout << "Run: "<< i <<" " << Run[i] << endl;
	  inputfile = "drei-r"+Run[i]+"_irene.root";
	  if(extralabel) 	  inputfile = "drei-r"+Run[i]+"_irene_"+label+".root";
	  if(dphcut)
	    {
	      inputfile = "drei-r"+Run[i]+"_irene_dphcutB"+ss_dphcut[l]+".root";
	      if(extralabel) inputfile = "drei-r"+Run[i]+"_irene_dphcutB"+ss_dphcut[l]+"_"+label+".root";
	    }
	  if(print) cout << "File Name: " << inputfile << endl;
	  Path=inputDir+inputfile;
	  file[i][l] = new TFile(Path);
	  if(print) cout << "File Path: " << Path << endl;


	  if(dphcut)
	    {
	      key = Run[i]+"_"+ss_dphcut[l]+"_"+pitch;
	    }
	  else  key = Run[i]+"_"+pitch;

	  for(int k = 0; k < hists ; k++){
	    cout << "key " << key << " hist " << Hist[k] << endl;
	  
	    auto it = map->find(std::make_pair(key,Hist[k]));
	  
	  
	    if(it  != map->end()){
	      cout << " it  != map->end() " << endl;
	      if(print)		  cout << " found map " << endl;
	      cout << "map key " << it->first.first << " " << it->first.second << endl;
	      
	      if(print) cout << Path << "  entries " << it->second->GetEntries() << endl;
	      
	    }
	    else{
	      cout << " it  == map->end() " << endl;
              cout << "l*runs+i " << l*runs+i << endl;
	      g = (TH1F*)file[i][l]->Get(Hist[k]);
              h = (TH1F *)g->Clone();
	      cout << " Hist " << Hist[k] << " " << h->GetEntries() << endl;
	      map->insert(std::make_pair(std::make_pair(key,Hist[k]),g));
	      it = map->find(std::make_pair(key,Hist[k]));
	      
	      cout << "map key " << it->first.first << " " << it->first.second << endl;
	    }
	      


	  //	      if(print) cout << Path << "  entries " << h->GetEntries() << endl;
	  //file[i][l]->Close();	  
	}//dphcuts
    }//runs
      return h;
    }
}

  //TH1F*
MapTH1   GetCheckHists(MapTH1 * map,int runs, int dphcuts, bool dphcut, TString * Run, TString * ss_dphcut, TString pitch, TString Hist, TH1F * h, bool extralabel = false, TString label = " ")
{
  TString inputDir="/home/zoiirene/Output/";
  TString inputfile, Path;

  TFile * file[runs][dphcuts];
  bool print = true;
  TH1F * g;
  TString key;
  if(print) cout << "Getting runs "  << endl;
  for(int i=0; i<runs; i++)
    {
      for(int l=0; l<dphcuts; l++)
	{

	  if(print) cout << "Run: "<< i <<" " << Run[i] << endl;
	  inputfile = "drei-r"+Run[i]+"_irene.root";
	  if(extralabel) 	  inputfile = "drei-r"+Run[i]+"_irene_"+label+".root";
	  if(dphcut)
	    {
	      inputfile = "drei-r"+Run[i]+"_irene_dphcutB"+ss_dphcut[l]+".root";
	      if(extralabel) inputfile = "drei-r"+Run[i]+"_irene_dphcutB"+ss_dphcut[l]+"_"+label+".root";
	    }
	  if(print) cout << "File Name: " << inputfile << endl;
	  Path=inputDir+inputfile;
	  file[i][l] = new TFile(Path);
	  if(print) cout << "File Path: " << Path << endl;


	  if(dphcut)
	    {
	      key = Run[i]+"_"+ss_dphcut[l]+"_"+pitch;
	    }
	  else  key = Run[i]+"_"+pitch;
	  cout << "key " << key << " hist " << Hist << endl;
	  
	  auto it = map->find(std::make_pair(key,Hist));
	  
	  
	  if(it  != map->end())
	    {
	      cout << " it  != map->end() " << endl;
	      if(print)		  cout << " found map " << endl;
	      cout << "map key " << it->first.first << " " << it->first.second << endl;

	      if(print) cout << Path << "  entries " << it->second->GetEntries() << endl;
	  	  
	    }
	  else
	    {
	      cout << " it  == map->end() " << endl;
              cout << "l*runs+i " << l*runs+i << endl;
	      g = (TH1F*)file[i][l]->Get(Hist);
	      g->SetDirectory(0);
              h = (TH1F *)g->Clone();
	      cout << " Hist " << Hist << " " << h->GetEntries() << endl;
	      map->insert(std::make_pair(std::make_pair(key,Hist),g));
	      it = map->find(std::make_pair(key,Hist));
	      
	      cout << "map key " << it->first.first << " " << it->first.second << endl;
	    }
	      


	  //	      if(print) cout << Path << "  entries " << h->GetEntries() << endl;
	  file[i][l]->Close();	  
	}//dphcuts
    }//runs


 return *map;
}//

void GetPixelavHists(MapTH1 * map,int runs, TString * Run, TString * pitch, TString Hist, TH1F * h)
{
  TString inputDir="/home/zoiirene/Output/";
  TString inputfile, Path;

  TFile * file[runs];
  bool print = true;
  TString key;
  if(print) cout << "Getting runs "  << endl;
  for(int i=0; i<runs; i++)
    {
      
      
      if(print) cout << "Run: "<< i <<" " << Run[i] << endl;
      inputfile = "pixelav-r"+Run[i]+".out.root";
      if(print) cout << "File Name: " << inputfile << endl;
      Path=inputDir+inputfile;
      file[i] = new TFile(Path);
      if(print) cout << "File Path: " << Path << endl;


       key = Run[i]+"_"+pitch[i];
       cout << "key " << key << " hist " << Hist << endl;
	  
       auto it = map->find(std::make_pair(key,Hist));
	  
	  
       if(it  != map->end())
	 {
	   cout << " it  != map->end() " << endl;
	   if(print)		  cout << " found map " << endl;
	   cout << "map key " << it->first.first << " " << it->first.second << endl;
	   
	   if(print) cout << Path << "  entries " << it->second->GetEntries() << endl;
	  	  
	 }
       else
	 {
	   cout << " it  == map->end() " << endl;
	   
	   h = (TH1F*)file[i]->Get(Hist);
	   h->SetDirectory(0);
	   cout << " Hist " << Hist << " " << h->GetEntries() << endl;
	   map->insert(std::make_pair(std::make_pair(key,Hist),h));
	   it = map->find(std::make_pair(key,Hist));
	   
	   cout << "map key " << it->first.first << " " << it->first.second << endl;
	 }
       


	  //	      if(print) cout << Path << "  entries " << h->GetEntries() << endl;
	

    }//runs
    
}//pixelav hists

void GetHists(MapTH1 * map,int runs, int dphcuts, int comparisons,  bool dphcut, TString * Run, int * i_dphcut, TString * Label, TString Hist, TH1F * h, bool extraname = false, TString name = " "){
  cout << "getting files " << endl;
  TString inputDir="/home/zoiirene/Output/";
  TString inputfile, Path;
  TString ss_dphcut[dphcuts];
  for(int i =0; i< dphcuts;i++)
    ss_dphcut[i].Form("%d",i_dphcut[i]);
  
  TFile * file[runs][dphcuts];
  TDirectory * histodir;
  TString key;
  bool print = true;
  if(print) cout << "Getting runs "  << endl;
  for(int i=0; i<runs; i++){
    for(int l=0; l<dphcuts; l++){
      if(print) cout << "Run: "<< i <<" " << Run[i] << endl;
      inputfile = "drei-r"+Run[i]+"_irene.root";
      if(extraname) inputfile = "drei-r"+Run[i]+"_irene_"+name+".root";
      if(dphcut){
	inputfile = "drei-r"+Run[i]+"_irene_dphcutB"+ss_dphcut[l]+".root";
	if(extraname) inputfile = "drei-r"+Run[i]+"_irene_dphcutB"+ss_dphcut[l]+"_"+name+".root";
      }	
      
      if(print) cout << "File Name: " << inputfile << endl;
      Path=inputDir+inputfile;
      file[i][l] = new TFile(Path);
      if(print) cout << "File Path: " << Path << endl;
      


      for(int j=0; j<comparisons; j++){
	key = Run[i]+"_"+Label[j];
	if(dphcut) key = Run[i]+"_"+ss_dphcut[l]+"_"+Label[j];
	if(print) cout << " key " << key << endl;
	auto it = map->find(std::make_pair(key,Hist));
	if(print) cout << "iterator" << endl;
	if(print) cout << "histodir " << Label[j] <<endl;
	
	histodir =(TDirectoryFile*)file[i][l]->Get(Label[j]);
	if(print) cout << "histodir" << endl;
	
	if(it  != map->end()){
	  cout << " it  != map->end() " << endl;
	if(print)           cout << " found map " << endl;
	cout << "map key " << it->first.first << " " << it->first.second << endl;	
		  
      }else
	 {
	   
	   h = (TH1F*)histodir->Get(Hist);

	   map->insert(std::make_pair(std::make_pair(key,Hist),h));
	   it = map->find(std::make_pair(key,Hist));
	   cout << "map key " << it->first.first << " " << it->first.second << endl;
	   
	 }
	      

      
      if(print) cout << Path << " " << Label[j] <<"  entries " << h->GetEntries() << endl;
      
    }//dphcuts
  }//comparisons
}//runs
  cout << " DONE GETTING HISTS!! ********************" << endl;
}//



void DrawThrScanHists(MapTH1 * map,int comparisons, TString Run, int * i_dphcut, TString Label, int hists,TString * Hist,float * xmin,float * xmax,float * ymin, float * ymax)
{

  TString ss_dphcut[comparisons];
  TString f_dphcut[comparisons];
  TString outputDir = "/home/zoiirene/Output/Plots/";
  int entries[hists][comparisons];
  float mean[hists][comparisons];
  float rms[hists][comparisons];
  float rmserror[hists][comparisons];
  double rms_d;
  double rmserr;
  TString name;
  name = outputDir+"compare_thresholds_"+Run+"_"+Label;//+"_"+Hist;
  ofstream myfile;
  myfile.open (name+".txt");
  myfile << "\\begin{adjustbox}{angle=90}\n";
  myfile << "\\begin{tabular}{|c|";
  for(int k=0; k<hists; k++)
    {
      if(Hist[k] == "dx3" || Hist[k] == "dx3_clchargeAC90evR" || Hist[k] == "dx3_clchargeABC90evR" )       myfile << "c|c|c|";
      else myfile << "c|c|";
    }
  myfile << "}\t\n";

  myfile << "\\hline \n";
  myfile << "thr &"; 

  
  for(int k=0; k<hists; k++)
    {
      if(Hist[k] == "dx3" || Hist[k] == "dx3_clchargeAC90evR" || Hist[k] == "dx3_clchargeABC90evR" )      myfile<< "  \\multicolumn{3}{ |c| }{" << Hist[k] << "}";
      else myfile<< "  \\multicolumn{2}{ |c| }{" << Hist[k] << "}";
      if(k!=hists-1)  myfile <<  "&";
    }
    myfile << "\\\\ \n";

    myfile << "\\hline \n";

  myfile << " &"; 
  for(int k=0; k<hists; k++)
    {
      if(Hist[k] == "dx3" || Hist[k] == "dx3_clchargeAC90evR" || Hist[k] == "dx3_clchargeABC90evR" )    myfile <<  " entries&  mean &rms ";
      else myfile <<  "   mean &rms ";
      if(k!=hists-1)  myfile <<  "&";
    }
  myfile << "\\\\ \n";
  myfile << "\\hline \n";

  for(int i =0; i< comparisons;i++)
    {
      ss_dphcut[i].Form("%d",i_dphcut[i]);
      f_dphcut[i] = (float) i_dphcut[i];
    }
  for(int k=0; k<hists; k++)
    { 
      TCanvas *c2 = new TCanvas("c2", "c2", 1500, 900);
      gPad->SetTicks(1,1);
      gROOT->SetStyle("Plain");
      gStyle->SetPadGridX(0);
      gStyle->SetPadGridY(0);
      gStyle->SetPalette(1);
      //      gStyle->SetOptStat();
      gStyle->SetOptStat(1110);
      gStyle->SetOptTitle(0);

  //Hist[0] = "dx3"; (entries, mean, rms)
  //Hist[1] = "clchargeB"; (entries, bin center at maximum)
  //Hist[2] = "clphB"; (entries, bin center at maximum)
  //Hist[3] = "clsizeB"; (entries, mean)
  //Hist[4] = "nrowvsxmB3";( mean on y, rms on y)
  //Hist[5] = "dx3_clchargeAC90evR";(entries, mean, rms)     
  //Hist[6] = "dx3_clchargeABC90evR";(entries, mean, rms)

  TString yaxistitle = "center of dx3 [mm] ";
  if(Hist[k] == "clsizeB") yaxistitle = " mean cluster size ";
  if(Hist[k] == "clphB") yaxistitle = " mpv [ADC] ";
  if(Hist[k] == "clchargeB") yaxistitle = " mpv [e-] ";
  
  
  int axis = 1;
  int color;
  for(int i = 0; i< comparisons; i++)
    {
      color = i+1;
      if(i+1 ==5) color = 95;
      if(i+1 ==10) color = 29;

      TString finding;
      finding = Run+"_"+ss_dphcut[i]+"_"+Label;
      cout << finding << endl;
      auto it3 = map->find(std::make_pair(finding,Hist[k]));

      if(it3  != map->end())
	{
	  if(Hist[k] == "nrowvsxmB3")
	    {
	      yaxistitle = "nrow";
	      axis = 2;
	    }

	  entries[k][i] = it3->second->GetEntries();
	  mean[k][i] = it3->second->GetMean(axis);
	  rms[k][i] = it3->second->GetRMS(axis);

	  if(Hist[k] == "clchargeB" || Hist[k] == "clphB")
	    mean[k][i] = it3->second->GetBinCenter(it3->second->GetMaximumBin());

	  
	  if(Hist[k] == "dx3"||   Hist[k] == "dx3_clchargeAC90evR" ||   Hist[k] == "dx3_clchargeABC90evR"    )
	    {
	      cout << " RMS 95" << endl;
	      cout << " COMMENTED FitTH1 " << endl;
	      //FitTH1(it3->second, &(rms_d), &(rmserr), "dphcut"+ss_dphcut[i], Run, Label, Hist[k],"RMS");
	      cout << " success" << endl;
	      rms[k][i] = (float) rms_d;
	      rmserror[k][i] = (float) rmserr;

	    }
	  
	  cout << entries[k][i] << " " << mean[k][i] << " " << rms[k][i] << endl;


	  it3->second->SetLineWidth(2);
	  it3->second->SetLineStyle(color);
	  it3->second->SetLineColor(color);
	  it3->second->GetXaxis()->SetRangeUser(xmin[k],xmax[k]);
	  it3->second->GetYaxis()->SetRangeUser(ymin[k],ymax[k]);
	  it3->second->Draw("histsames");
	  c2->Update();

	  /*
	  TPaveStats *Stats =   (TPaveStats*)it3->second->GetListOfFunctions()->FindObject("stats");
	      
	  Stats->SetX1NDC(0.55);
	  Stats->SetX2NDC(.65);
	  Stats->SetY1NDC(.3+0.1*i);
	  Stats->SetY2NDC(.4+0.1*i);
	  
	  if(i>5)
	    {
	      Stats->SetX1NDC(0.65);
	      Stats->SetX2NDC(.75);
	      Stats->SetY1NDC(.3+0.1*(i-6));
	      Stats->SetY2NDC(.4+0.1*(i-6));
	    }
	  if(i>11)
	    {
	      Stats->SetX1NDC(0.75);
	      Stats->SetX2NDC(.85);
	      Stats->SetY1NDC(.3+0.1*(i-12));
	      Stats->SetY2NDC(.4+0.1*(i-12));
	    }

	  if(Hist[k] == "nrowvsxmB3")
	    {
	      Stats->SetX1NDC(0.15);
	      Stats->SetX2NDC(.25);
	      Stats->SetY1NDC(.15+0.1*i);
	      Stats->SetY2NDC(.25+0.1*i);

	      if(i>1)
		{
		  Stats->SetX1NDC(0.25);
		  Stats->SetX2NDC(.35);
		  Stats->SetY1NDC(.15+0.1*(i-2));
		  Stats->SetY2NDC(.25+0.1*(i-2));
		}
	      if(i>3)
		{
		  Stats->SetX1NDC(0.35);
		  Stats->SetX2NDC(.45);
		  Stats->SetY1NDC(.15+0.1*(i-4));
		  Stats->SetY2NDC(.25+0.1*(i-4));
		}
	      if(i>5)
		{
		  Stats->SetX1NDC(0.45);
		  Stats->SetX2NDC(.55);
		  Stats->SetY1NDC(.15+0.1*(i-6));
		  Stats->SetY2NDC(.25+0.1*(i-6));
		}
	      if(i>7)
		{
		  Stats->SetX1NDC(0.55);
		  Stats->SetX2NDC(.65);
		  Stats->SetY1NDC(.15+0.1*(i-8));
		  Stats->SetY2NDC(.25+0.1*(i-8));
		}
	      if(i>9)
		{
		  Stats->SetX1NDC(0.65);
		  Stats->SetX2NDC(.75);
		  Stats->SetY1NDC(.15+0.1*(i-10));
		  Stats->SetY2NDC(.25+0.1*(i-10));
		}


	    }  
	  Stats->SetTextColor(color);
	  Stats->SetTextSize(0.02);
	  */	  
	  gPad->Modified();
	  //      c2->Update();
	}
      else cout << "not found! " <<endl;
    


    }


  
  name += "_"+Hist[k]; 
  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".pdf");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".root");

  DrawTGraph(comparisons,i_dphcut,entries[k],Hist[k],Run,Label,"entries", i_dphcut[0],i_dphcut[comparisons-1],0,60000);
  if(Hist[k] == "nrowvsxmB3") xmax[k]=4;
  DrawTGraphError(comparisons,i_dphcut, mean[k], rms[k],Hist[k],Run,Label,"meanandRMS" , yaxistitle,  i_dphcut[0],i_dphcut[comparisons-1],xmin[k],xmax[k]);
  if(Hist[k] == "dx3" ||   Hist[k] == "dx3_clchargeAC90evR" ||   Hist[k] == "dx3_clchargeABC90evR"     )
    {
      yaxistitle = "RMS 95% [#mum]";
      double xminr,xmaxr;
      if(Hist[k] == "dx3")
	{
	  xminr=0.005;//0.0038;
	  xmaxr=0.007;//0.0045;
	}
      if( Hist[k] == "dx3_clchargeAC90evR")
	{
	  xminr=0.004;////0.0035;
	  xmaxr=0.006;//0.0038;
	}
      if( Hist[k] == "dx3_clchargeABC90evR"     )
	{
	  xminr=0.004;//0.0032;//0.004
	  xmaxr=0.005;//0.0035;//0.006
	}
      DrawTGraphWithError(comparisons, i_dphcut, rms[k], rmserror[k], Hist[k],Run, Label, "RMS95_error", yaxistitle,  i_dphcut[0],i_dphcut[comparisons-1],xminr,xmaxr);
      DrawTGraphError(comparisons,i_dphcut, rms[k], rmserror[k],Hist[k],Run,Label,"RMS95" , yaxistitle,  i_dphcut[0],i_dphcut[comparisons-1],0.,7.);
    }
    }
  
  
     for(int i=0;i<comparisons;i++)
	{
	  myfile << i_dphcut[i] ;

	  for(int k=0; k<hists; k++)
	    {
	      if(Hist[k] == "dx3" || Hist[k] == "dx3_clchargeAC90evR" || Hist[k] == "dx3_clchargeABC90evR" )   	      myfile <<  "& " << entries[k][i] << "& " << setprecision(2) <<mean[k][i] << "& " << setprecision(2) << rms[k][i] ;
	      else myfile <<  "& " << setprecision(2) <<mean[k][i] << "& " << setprecision(2) << rms[k][i] ;
	    }
	  myfile << "\\\\ \n";
	  myfile << "\\hline \n";
	}
     myfile<< "\\end{tabular} \n";
     myfile << "\\end{adjustbox} \n";

    myfile.close();

}




void DrawHists(MapTH1 * map,int comparisons, TString Run, bool dphcut, TString ss_dphcut, TString * Label, TString Hist,float xmin,float xmax,float ymin, float ymax)
{

  TString outputDir = "/home/zoiirene/Output/Plots/";
  
  TCanvas *c2 = new TCanvas("c2", "c2", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  //  gStyle->SetOptStat();
  gStyle->SetOptStat(1110);

  gStyle->SetOptTitle(0);

  TLegend* leg2 = new TLegend(0.15,0.55,0.45,0.75);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.02);
  int color;
  for(int i = 0; i< comparisons; i++)
    {
      color = i+1;
      if(i+1 ==5) color = 95;
      if(i+1 ==10) color = 29;

      TString finding;
      if(dphcut) finding = Run+"_"+ss_dphcut+"_"+Label[i];
      else finding = Run+"_"+Label[i];
      cout << finding << endl;
      auto it3 = map->find(std::make_pair(finding,Hist));
      if(it3  != map->end())
	{

	  it3->second->GetXaxis()->SetLabelFont(42);
	  it3->second->GetXaxis()->SetLabelSize(0.04);
	  it3->second->GetXaxis()->SetTitleSize(0.05);
	  it3->second->GetXaxis()->SetTitleOffset(0.8);
	  it3->second->GetXaxis()->SetTitleFont(42);

	  it3->second->GetYaxis()->SetLabelFont(42);
	  it3->second->GetYaxis()->SetLabelSize(0.04);
	  it3->second->GetYaxis()->SetTitleSize(0.05);
	  it3->second->GetYaxis()->SetTitleOffset(1.);
	  it3->second->GetYaxis()->SetTitleFont(42);
	  it3->second->GetYaxis()->SetNoExponent(3);	  
	  //	  SetExponentOffset()
	  it3->second->SetLineWidth(2);
	  it3->second->SetLineStyle(color);
	  it3->second->SetLineColor(color);
	  it3->second->GetXaxis()->SetRangeUser(xmin,xmax);
	  it3->second->GetYaxis()->SetRangeUser(ymin,ymax);
	  //	  leg2->AddEntry(it3->second,Label[i], "l");
	  it3->second->Draw("histsames");
	  c2->Update();
	  TPaveStats *Stats =   (TPaveStats*)it3->second->GetListOfFunctions()->FindObject("stats");
	      
	  Stats->SetX1NDC(0.55);
	  Stats->SetX2NDC(.75);
	  Stats->SetY1NDC(.3+0.1*i);
	  Stats->SetY2NDC(.4+0.1*i);
	  
	  if(i>3)
	    {
	      Stats->SetX1NDC(0.75);
	      Stats->SetX2NDC(.95);
	      Stats->SetY1NDC(.3+0.1*(i-3));
	      Stats->SetY2NDC(.4+0.1*(i-3));
	    }

	  if(Hist == "nrowvsxmB3")
	    {
	      Stats->SetX1NDC(0.15);
	      Stats->SetX2NDC(.35);
	      Stats->SetY1NDC(0.15+0.1*i);
	      Stats->SetY2NDC(0.25+0.1*i);

	      if(i>3)
		{
		  Stats->SetX1NDC(0.35);
		  Stats->SetX2NDC(.55);
		  Stats->SetY1NDC(0.15+0.1*(i-4));
		  Stats->SetY2NDC(0.25+0.1*(i-4));
		}
	      if(i>7)
		{
		  Stats->SetX1NDC(0.55);
		  Stats->SetX2NDC(.75);
		  Stats->SetY1NDC(0.15+0.1*(i-8));
		  Stats->SetY2NDC(0.25+0.1*(i-8));
		}
	      

	    }  
	  Stats->SetTextColor(color);
	  Stats->SetTextSize(0.02);
	  gPad->Modified();
	  //      c2->Update();
	}
      else cout << "not found! " <<endl;
    }



  //  if(Hist == "clchargeB") gStyle->SetOptStat(0);

  //  leg2->Draw();

  TString name;
  if(dphcut)  name = outputDir+"compare_"+Run+"_"+ss_dphcut+"_"+Hist;
  else name = outputDir+"compare_"+Run+"_"+Hist;
  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".pdf");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".root");
  //gStyle->SetOptStat(1);
}//DrawHists


void DrawTwoHist(TH1F * h1,TH1F * h2, TString Label1, TString Label2, TString Name,float xmin,float xmax,float ymin, float ymax, bool log = true)
{

  TString outputDir = "/home/zoiirene/Output/Plots/";
  
  TCanvas *c2 = new TCanvas("c2", "c2", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  //gStyle->SetOptStat();
  gStyle->SetOptStat(111111);
  gStyle->SetOptTitle(0);

  
  h1->SetLineWidth(2);
  h1->SetLineColor(kBlue);
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->GetYaxis()->SetRangeUser(ymin,ymax);
	  //	  leg2->AddEntry(h,Label[i], "l");
  h1->Draw("hist");
  c2->Update();
  TPaveStats *Stats =   (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
  Stats->SetX1NDC(0.55);
  Stats->SetX2NDC(.75);
  Stats->SetY1NDC(.4+0.15);
  Stats->SetY2NDC(.5+0.15);
  Stats->SetTextColor(kBlue);
  Stats->SetTextSize(0.02);
  gPad->Modified();
  

  
  h2->SetLineWidth(2);
  h2->SetLineColor(kRed);
  h2->SetLineStyle(2);
  h2->GetXaxis()->SetRangeUser(xmin,xmax);
  h2->GetYaxis()->SetRangeUser(ymin,ymax);
  
  h2->Draw("histsames");
  c2->Update();
  TPaveStats *Stats2 =   (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
  Stats2->SetX1NDC(0.55);
  Stats2->SetX2NDC(.75);
  Stats2->SetY1NDC(.4+0.3);
  Stats2->SetY2NDC(.5+0.3);
  Stats2->SetTextColor(kRed);
  Stats2->SetTextSize(0.02);
  gPad->Modified();
  
  
  TLegend* leg2 = new TLegend(0.7,0.3,0.85,0.45);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.02);
  leg2->AddEntry(h1,Label1, "l");
  leg2->AddEntry(h2,Label2, "l");
  leg2->Draw("same");
  c2->Update();
  if(log) c2->SetLogy();
  TString name = outputDir+"compare_"+Label1+"_"+Label2+"_"+"_"+Name;
  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".pdf");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".root");
}


void DrawHist(TH1F * h, TString Run, TString ss_dphcut, TString Label, TString Hist,float xmin,float xmax,float ymin, float ymax)
{

  TString outputDir = "/home/zoiirene/Output/Plots/";
  
  TCanvas *c2 = new TCanvas("c2", "c2", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  //gStyle->SetOptStat();
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);


  h->SetLineWidth(2);
  h->SetLineColor(kBlue);
  h->GetXaxis()->SetRangeUser(xmin,xmax);
  h->GetYaxis()->SetRangeUser(ymin,ymax);
	  //	  leg2->AddEntry(h,Label[i], "l");
  h->Draw("hist");
  if(Hist == "clchargeB") gStyle->SetOptStat(0);
  c2->Update();
  TString name = outputDir+"compare_"+Run+"_"+Label+"_"+ss_dphcut+"_"+Hist;
  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".pdf");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".root");
}


void DrawCheckHists(MapTH1 * map, int runs, TString * Run, bool dphcut, TString ss_dphcut, TString pitch, TString Hist,float xmin,float xmax,float ymin, float ymax)
{

  TString outputDir = "/home/zoiirene/Output/Plots/";
  
  TCanvas *c2 = new TCanvas("c2", "c2", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);

  int color;
  for(int i = 0; i< runs; i++)
    {
      color = i+1;
      if(i+1 ==5) color = 95;
      if(i+1 ==10) color = 29;

      TString finding;
      if(dphcut) finding = Run[i]+"_"+ss_dphcut+"_"+pitch;
      else finding = Run[i]+"_"+pitch;
      cout << finding << endl;
      auto it3 = map->find(std::make_pair(finding,Hist));
      if(it3  != map->end())
	{

	  it3->second->SetLineWidth(2);
	  it3->second->SetLineStyle(color);
	  it3->second->SetLineColor(color);
	  it3->second->GetXaxis()->SetRangeUser(xmin,xmax);
	  it3->second->GetYaxis()->SetRangeUser(ymin,ymax);
	  it3->second->Draw("histsames");
	  c2->Update();
	  TPaveStats *Stats =   (TPaveStats*)it3->second->GetListOfFunctions()->FindObject("stats");
	  Stats->SetX1NDC(0.55);
	  Stats->SetX2NDC(.75);
	  Stats->SetY1NDC(.3+0.1*i);
	  Stats->SetY2NDC(.4+0.1*i);
	  
	  if(i>5)
	    {
	      Stats->SetX1NDC(0.75);
	      Stats->SetX2NDC(.95);
	      Stats->SetY1NDC(.3+0.1*(i-5));
	      Stats->SetY2NDC(.4+0.1*(i-5));
	    }



	  Stats->SetTextColor(color);
	  Stats->SetTextSize(0.02);
	  gPad->Modified();

	}
      else cout << "not found! " <<endl;
    }


  TString name;
  name = outputDir+"compare_alignment_"+pitch+"_"+Hist;
  
  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".pdf");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".root");
}


void DrawCheckLegend(MapTH1 * map,int runs, TString * Run, TString ss_dphcut, bool dphcut, TString pitch, TString Hist)
{
  TString outputDir = "/home/zoiirene/Output/Plots/";
  
  TCanvas *c2 = new TCanvas("c2", "c2", 1500, 900);
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);

  TLegend* leg2 = new TLegend(0.15,0.25,0.85,0.85);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.04);
  int color;
  for(int i = 0; i< runs; i++)
    {
      color = i+1;
      if(i+1 ==5) color = 95;
      if(i+1 ==10) color = 29;

      TString finding;
      if(dphcut) finding = Run[i]+"_"+ss_dphcut+"_"+pitch;
      else finding = Run[i]+"_"+pitch;
      cout << finding << endl;
      auto it3 = map->find(std::make_pair(finding,Hist));
      if(it3  != map->end())
	{
	  it3->second->SetLineColor(color);
	  leg2->AddEntry(it3->second,finding, "l");
	  c2->Update();
	}
      else cout << "not found! " <<endl;
    }

  leg2->Draw();

  TString name;
  name = outputDir+"compare_alignment_"+pitch+"_"+Hist+"_legend";

  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".pdf");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".root");
}



void DrawLegend(MapTH1 * map,int comparisons, TString Run, TString ss_dphcut, bool dphcut, TString * Label, TString Hist)
{
  TString outputDir = "/home/zoiirene/Output/Plots/";
  
  TCanvas *c2 = new TCanvas("c2", "c2", 1500, 900);
  //gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  //gStyle->SetPadGridX(0);
  //gStyle->SetPadGridY(0);
  //gStyle->SetPalette(1);
  //gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);

  TLegend* leg2 = new TLegend(0.15,0.25,0.85,0.85);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.04);
  int color;
  for(int i = 0; i< comparisons; i++)
    {
      color = i+1;
      if(i+1 ==5) color = 95;
      if(i+1 ==10) color = 29;

      TString finding;
      if(dphcut) finding = Run+"_"+ss_dphcut+"_"+Label[i];
      else finding = Run+"_"+Label[i];
      cout << finding << endl;
      auto it3 = map->find(std::make_pair(finding,Hist));
      if(it3  != map->end())
	{

	  //	  it3->second->SetLineWidth(2);
	  //it3->second->SetLineStyle(color);
	  it3->second->SetLineColor(color);
	  //it3->second->GetXaxis()->SetRangeUser(xmin,xmax);
	  //it3->second->GetYaxis()->SetRangeUser(ymin,ymax);
	  leg2->AddEntry(it3->second,Label[i], "l");
	  //it3->second->Draw("histsames");
	  c2->Update();
	      

	}
      else cout << "not found! " <<endl;
    }





  leg2->Draw();
  TString name;
  if(dphcut)  name = outputDir+"compare_"+Run+"_"+ss_dphcut+"_"+Hist+"_legend";
  else name = outputDir+"compare_"+Run+"_"+Hist+"_legend";

  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".pdf");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".root");
}

void DrawThresholdLegend(MapTH1 * map,int comparisons, TString Run, TString * ss_dphcut, TString  Label, TString Hist)
{
  TString outputDir = "/home/zoiirene/Output/Plots/";
  
  TCanvas *c2 = new TCanvas("c2", "c2", 1500, 900);
   gROOT->SetStyle("Plain");
   gStyle->SetOptTitle(0);

  TLegend* leg2 = new TLegend(0.15,0.25,0.85,0.85);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.04);
  int color;
  for(int i = 0; i< comparisons; i++)
    {
      color = i+1;
      if(i+1 ==5) color = 95;
      if(i+1 ==10) color = 29;

      TString finding;
      finding = Run+"_"+ss_dphcut[i]+"_"+Label;
      cout << finding << endl;
      auto it3 = map->find(std::make_pair(finding,Hist));
      if(it3  != map->end())
	{
	  it3->second->SetLineColor(color);
	  leg2->AddEntry(it3->second,ss_dphcut[i], "l");
	  c2->Update();
	      

	}
      else cout << "not found! " <<endl;
    }





  leg2->Draw();
  TString name;
  name = outputDir+"compare_threshold_"+Run+"_"+Label+"_"+Hist+"_legend";


  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".pdf");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".root");
}
