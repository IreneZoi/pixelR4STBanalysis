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
#define hists 5
bool print=true;
using namespace std;
TString outputDir = "/home/zoiirene/Output/Plots/";



void compare2743Alignments()//TString hist1 = "clszBiii",TString hist2 = "clszBiii",TString run1 = "2743", TString label1 = "p-stop", TString run2 = "2832", TString label2 = "p-spray", bool norm = true)
{

  TString inputDir="/home/zoiirene/Output/";
  TString inputfile, Path;
 
  TFile * file[runs];
  TDirectory * histodir;
  TString Run[runs];
  Run[0] = "2743";
  TString Hist[hists];
  Hist[0] = "dx3";
  Hist[1] = "clchargeB";
  Hist[2] = "clphB";
  Hist[3] = "clsizeB";
  Hist[4] = "nrowvsxmB3";
  //  Hist[5] = "dx3_clchargeAC90evR";
  //Hist[6] = "dx3_clchargeABC90evR";
  TString Label[comparisons];
  Label[0] = "straightTracksY_isoAandCandB_straightTracksX";
  bool dphcut = false;

  MapTH1 m_hists_irene;
  MapTH1 m_hists_irene0;
  MapTH1 m_hists_daniel;

  int dph[dphcuts] = {12};  

  std::map<std::pair<TString, TString>, TH1F *>::iterator it;
  TH1F * h[hists];
  for(int k=0; k<hists; k++)
    {
      GetHists(&m_hists_irene,runs, dphcuts, comparisons,  dphcut, Run, dph, Label,Hist[k], h[k]);
      GetHists(&m_hists_daniel,runs, dphcuts, comparisons,  dphcut, Run, dph, Label,Hist[k], h[k],true, "alignDaniel");
      GetHists(&m_hists_irene0,runs, dphcuts, comparisons,  dphcut, Run, dph, Label,Hist[k], h[k],true, "iteration0");
    }

  double rms;
  double rmserr;
  for(int i=0; i<runs; i++)
    {
      for(int l=0; l<dphcuts; l++)
	{
	  
	  for(int j=0; j<comparisons; j++)
	  {

	    for(int k=0; k<hists; k++)
	  	{
		  if(k==0|| k==5 || k==6)
		    {
		      auto it2 = m_hists_irene.find(std::make_pair(Run[i]+"_"+Label[j],Hist[k]));
		      if(it2  != m_hists_irene.end())
			{
			  if(print)               cout << " found map " << endl;
			  if(print)		      cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
			  FitTH1(it2->second, &rms,&rmserr, "irene", Run[0], Label[0], Hist[k],"RMS");
			}
		      auto it3 = m_hists_daniel.find(std::make_pair(Run[i]+"_"+Label[j],Hist[k]));
		      if(it3  != m_hists_daniel.end())
			{
			  if(print)               cout << " found map " << endl;
			  if(print)		      cout << "map key " << it3->first.first << " " << it3->first.second << " " << it3->second->GetEntries() << endl;
			  FitTH1(it3->second, &rms,&rmserr, "daniel", Run[0], Label[0], Hist[k],"RMS");
			}
		      
		      auto it4 = m_hists_irene0.find(std::make_pair(Run[i]+"_"+Label[j],Hist[k]));
		      if(it4  != m_hists_irene0.end())
			{
			  if(print)               cout << " found map " << endl;
			  if(print)		      cout << "map key " << it4->first.first << " " << it4->first.second << " " << it4->second->GetEntries() << endl;
			  FitTH1(it4->second, &rms,&rmserr, "irene0", Run[0], Label[0], Hist[k],"RMS");
			}
		      
		    }

		}//hists
	    }//dphcuts
	}//comparisons
    }//runs


 
}//resolution 

