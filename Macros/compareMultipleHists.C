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
#define dphcuts 1
#define hists 7
bool print=true;
using namespace std;
TString outputDir = "/home/zoiirene/Output/Plots/";



void compareMultipleHists()//TString hist1 = "clszBiii",TString hist2 = "clszBiii",TString run1 = "2743", TString label1 = "p-stop", TString run2 = "2832", TString label2 = "p-spray", bool norm = true)
{

  TString inputDir="/home/zoiirene/Output/";
  TString inputfile, Path;
 
  TFile * file[runs];
  TDirectory * histodir;
  //TH1F * h[runs*dphcuts][comparisons][hists];
  TString Run[runs];
  Run[0] = "3778";
  TString Hist[hists];
  Hist[0] = "dx3";
  Hist[1] = "clchargeB";
  Hist[2] = "clphB";
  Hist[3] = "clsizeB";
  Hist[4] = "nrowvsxmB3";
  Hist[5] = "dx3_clchargeAC90evR";
  Hist[6] = "dx3_clchargeABC90evR";
  int i_dphcut[dphcuts];
  i_dphcut[0] = 15;
  TString ss_dphcut[dphcuts];
  ss_dphcut[0] = "15";
  TString Label[comparisons];
  Label[0] = "straightTracksY_isoAandCandB_straightTracksX";
  //Label[1] = "dxCAcut";
  /*  Label[0] = "beforeCorrections";
  Label[1] = "nocuts";
  Label[2] = "straightTracksY";
  Label[3] = "straightTracksY_isoAandC";
  Label[4] = "straightTracksY_isoAandCandB";
  Label[5] = "straightTracksY_isoAandCandB_straightTracksX";
  Label[6] = "hitsOnTrack";

  Label[7] = "straightTracksY_isoAandCandB_chargerAandC";
  Label[8] = "straightTracksY_isoAandCandB_chargerAandCandB";
  Label[9] = "straightTracksY_isoAandCandB_chargeAandC";
  Label[10] = "straightTracksY_isoAandCandB_chargeAandCandB";
  */  
bool dphcut = false;

  MapTH1 m_hists;
  std::map<std::pair<TString, TString>, TH1F *>::iterator it;
  TH1F * h[hists];
  for(int k=0; k<hists; k++)
    {
      GetHists(&m_hists,runs, dphcuts, comparisons,  dphcut, Run, i_dphcut, Label,Hist[k], h[k]);
      for ( it = m_hists.begin(); it != m_hists.end(); it++ )
	{
	  it->second->GetEntries();
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
		  auto it2 = m_hists.find(std::make_pair(Run[i]+"_"+ss_dphcut[l]+"_"+Label[j],Hist[k]));
		  if(it2  != m_hists.end())
		    {
		      if(print)               cout << " found map " << endl;
		      if(print)		      cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;  
		    }


		}//hists
	    }//dphcuts
	}//comparisons
    }//runs
  
  /*
  if(norm)
    {
      float normaliz1 = h[0]->Integral();
      float normaliz2 = h[1]->Integral();
      if(print) cout << "normaliz 1 "  << normaliz1 << "  2  " << normaliz2 << " ratio 2/1 " << normaliz2/normaliz1 << endl;
      h[0]->Scale(normaliz2/normaliz1);
    }
  */
  float xmin[hists]={-0.1,0,0,0,0,-0.1,-0.1};
  float xmax[hists]={0.1,40,600,10,25,0.1,0.1};
  float ymin[hists]={0,0,0,0,0,0,0};
  float ymax[hists]={8000,8000,4000,40000,4.,8000,8000};

  for(int i=0; i<runs; i++)
    {
      for(int l=0; l<dphcuts; l++)
	{

	  for(int k=0; k<hists; k++)
	    {
	      DrawHists(&m_hists,comparisons, Run[i], dphcut, ss_dphcut[l], Label, Hist[k],xmin[k],xmax[k],ymin[k],ymax[k]);
	      //	      DrawLegend(&m_hists,comparisons, Run[i], ss_dphcut[l], dphcut,Label, Hist[k]);
    
	    }
	}
    }
  if(print) cout << "Plotting"  << endl;
 
}//resolution 

