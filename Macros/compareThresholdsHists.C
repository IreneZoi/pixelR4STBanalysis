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

#define comparisons 1
#define runs 1
#define dphcuts 12
#define hists 5
bool print=true;
using namespace std;
TString outputDir = "/home/zoiirene/Output/Plots/";



void compareThresholdsHists()//TString hist1 = "clszBiii",TString hist2 = "clszBiii",TString run1 = "2743", TString label1 = "p-stop", TString run2 = "2832", TString label2 = "p-spray", bool norm = true)
{

  TString inputDir="/home/zoiirene/Output/";
  TString inputfile, Path;
 
  TFile * file[runs];
  TDirectory * histodir;
  TString Run[runs];
  Run[0] = "2801";
  TString Hist[hists];
  Hist[0] = "dx3";
  Hist[1] = "clchargeB";
  Hist[2] = "clphB";
  Hist[3] = "clsizeB";
  Hist[4] = "nrowvsxmB3";
  //  Hist[5] = "dx3_clchargeAC90evR";
  //Hist[6] = "dx3_clchargeABC90evR";
  int i_dphcut[dphcuts];
  i_dphcut[0] =  9 ;
  i_dphcut[1] =  10 ;
  i_dphcut[2] =  11 ;
  i_dphcut[3] =  12 ;
  i_dphcut[4] =  13 ;
  i_dphcut[5] =  14 ;
  i_dphcut[6] =  15 ;
  i_dphcut[7] =  18 ;
  i_dphcut[8] =  20 ;
  i_dphcut[9] =  22 ;
  i_dphcut[10] =  25 ;
  i_dphcut[11] =  27 ;

  double sigma[dphcuts];
  double sigmaerr[dphcuts];

  TString ss_dphcut[dphcuts];
  for(int i =0; i< dphcuts;i++)
    ss_dphcut[i].Form("%d",i_dphcut[i]);

  TString Label[comparisons];
  /*
  Label[0] = "raw";
  Label[1] = "dxCAcut";
  */  
  //Label[0] = "beforeCorrections";
  //Label[1] = "nocuts";
  //Label[2] = "straightTracksY";
  //Label[3] = "straightTracksY_isoAandC";
  //Label[4] = "straightTracksY_isoAandCandB";
  //Label[5] = "straightTracksY_isoAandCandB_straightTracksX";
  Label[0] = "straightTracksY_isoAandCandB_straightTracksX";
  //Label[6] = "hitsOnTrack";
  /*
  Label[7] = "straightTracksY_isoAandCandB_chargerAandC";
  Label[8] = "straightTracksY_isoAandCandB_chargerAandCandB";
  Label[9] = "straightTracksY_isoAandCandB_chargeAandC";
  Label[10] = "straightTracksY_isoAandCandB_chargeAandCandB";
  */  
bool dphcut = true;

  MapTH1 m_hists;
  std::map<std::pair<TString, TString>, TH1F *>::iterator it;
  TH1F * h[hists];
  TString outputDir = "/home/zoiirene/Output/Plots/";

  TString name;

  for(int k=0; k<hists; k++)
    {
      GetHists(&m_hists,runs, dphcuts, comparisons,  dphcut, Run, i_dphcut, Label,Hist[k], h[k]);

      for ( it = m_hists.begin(); it != m_hists.end(); it++ )
	{
	  it->second->GetEntries();
	}
    }

  if(print) cout << " ***************** GOT ALL HISTS ********************** " << endl;
  ofstream myfile[hists];


  for(int k=0; k<hists; k++)
    {
      if(k ==0 || k == 5 || k ==6 )
	{

	  name = outputDir+"compare_thresholdsResolution_"+Run[0]+"_"+Label[0]+"_"+Hist[k];//+"_"+Hist;
	  myfile[k].open (name+".txt");
	  
	  myfile[k] << "\\begin{tabular}{|c|c|c|c|}\t\n";
	  myfile[k] << "\\hline \n";
	  myfile[k] << "dphcut &  \\multicolumn{3}{ |c| }{" << Hist[k] << "}\\\\ \n";
	  myfile[k] << "\\hline \n";
	  
	  myfile[k] << " & entries&  sigma [$\\mu m$]& sigma err [$\\mu m$]\\\\ \n";
	  myfile[k] << "\\hline \n";
	}
    }

  
  for(int i=0; i<runs; i++)
    {
      for(int l=0; l<dphcuts; l++)
	{
	  
	  for(int j=0; j<comparisons; j++)
	    {
	      
	      for(int k=0; k<hists; k++)
		{
		  if(k ==0 || k == 5 || k ==6 )
		    {
		      myfile[k] << i_dphcut[l] ;
		      
		      auto it2 = m_hists.find(std::make_pair(Run[i]+"_"+ss_dphcut[l]+"_"+Label[j],Hist[k]));
		      if(it2  != m_hists.end())
			{
			  if(print)               cout << " found map " << endl;
			  if(print)		      cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;  

			  FitTH1(it2->second, &(sigma[k]), &(sigmaerr[k]), "dphcut"+ss_dphcut[l], Run[i], Label[j], Hist[k], "RMS");
			  myfile[k] <<  "& " <<  setprecision(5) <<  it2->second->GetEntries() << "& " << setprecision(3) <<sigma[k] << "& " << setprecision(1) << sigmaerr[k];
			    if(print) cout << " ***************** done RMS ********************** " << endl;
			}
		      myfile[k] << "\\\\ \n";
		      myfile[k] << "\\hline \n";

		    }  
		      
		}//hists
	      
	    }//dphcuts
	}//comparisons
    }//runs
  
  for(int k=0; k<hists; k++)
    {
      if(k ==0 || k == 5 || k ==6)
	{
	  
	  myfile[k]<< "\\end{tabular} \n";
	  
	  
	  myfile[k].close();
	}
    }
  


  if(print) cout << " ***************** wrote tab file ********************** " << endl;
  
  if(print) cout << "Plotting"  << endl;
  float xmin[hists]={-0.05,0,0,0,0};//,-0.05,-0.05};
  float xmax[hists]={0.05,50,600,10,25};//,0.05,0.05};
  float ymin[hists]={0,0,0,0,0};//,0,0};
  float ymax[hists]={12000,4000,2500,50000,4.};//,10000,10000};

  DrawThrScanHists(&m_hists,dphcuts, Run[0], i_dphcut, Label[0], hists,Hist,xmin,xmax,ymin,ymax);
  DrawThresholdLegend(&m_hists,dphcuts, Run[0], ss_dphcut, Label[0], Hist[0]);
  
    
  if(print) cout << " ***************** DONEEEE ********************** " << endl;
 
}//resolution 

  //  Run[0] = "2775";2773
/*
i_dphcut[0] =  0 ;
  i_dphcut[1] =  9 ;
  i_dphcut[2] =  10 ;
  i_dphcut[3] =  11 ;
  i_dphcut[4] =  12 ;
  i_dphcut[5] =  13 ;
  i_dphcut[6] =  14 ;
  i_dphcut[7] =  15 ;
  i_dphcut[8] =  18 ;
  i_dphcut[9] =  20 ;
  i_dphcut[10] =  22 ;
  i_dphcut[11] =  25 ;
  i_dphcut[12] =  27 ;
  i_dphcut[13] =  30 ;
  i_dphcut[14] =  32 ;
  i_dphcut[15] =  35 ;
  i_dphcut[16] =  37 ;
  i_dphcut[17] =  40 ;
*/

/*

  Run[0] = "2743";//2832
i_dphcut[0] =  0 ;
  i_dphcut[1] =  6 ;
  i_dphcut[2] =  9 ;
  i_dphcut[3] =  10 ;
  i_dphcut[4] =  11 ;
  i_dphcut[5] =  12 ;
  i_dphcut[6] =  13 ;
  i_dphcut[7] =  14 ;
  i_dphcut[8] =  15 ;
  i_dphcut[9] =  18 ;
  i_dphcut[10] =  20 ;
  i_dphcut[11] =  22 ;
  i_dphcut[12] =  25 ;
  i_dphcut[13] =  27 ;
  i_dphcut[14] =  30 ;
  i_dphcut[15] =  32 ;
  i_dphcut[16] =  35 ;
  i_dphcut[17] =  37 ;
  i_dphcut[18] =  40 ;
*/
