void drawps2d( int chipNr , const char* path = "../x/" ){

    // needs to match scanhold2d comand in cmd_analyzer.cpp
  const int bVpr = 400;
  const int eVpr = 801;
  const int sVpr = 100;
  
  const int bVsh = 400;
  const int eVsh = 901;
  const int sVsh = 100;
  
  const int szVpr = 1 + ( eVpr - bVpr ) / sVpr;
  const int szVsh = 1 + ( eVsh - bVsh ) / sVsh;
  
    // book canvas
  TCanvas* cCnstVpr[szVsh];
  TCanvas* cCnstVsh[szVpr];
  
    // calculate colors
  int np = szVpr;
  if( szVpr < szVsh )
    np = szVsh;
  int nc = gStyle->GetNumberOfColors();
    
    // style
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetPalette(kRainBow);
  
    //loop over VgPr and VgSh settings
  for( int ipr = 0; ipr < szVpr; ++ipr ){
    for( int ish = 0; ish < szVsh; ++ish ){
        
      int vpr = bVpr + ipr * sVpr;
      int vsh = bVsh + ish * sVsh;
      
        // open root files
      cout << "Trying to open: '" << Form( "%sschld2d/c%i_scanhold_pr%i_sh%i.root" , path, chipNr, vpr, vsh )  << "'" << endl;
      TFile* f = new TFile( Form( "%sschld2d/c%i_scanhold_pr%i_sh%i.root" , path, chipNr, vpr, vsh ) );
        // check for success
      if( f->IsOpen() )
        cout << "Succeed !" << endl;
      else
        return 0;
      
        // get analogue current
      TProfile* iavsev = new TProfile();
      iavsev = (TProfile*) f->Get( "iavsev" );
      int iana = (int)iavsev->GetMean(2);
        
        //set name for the pulse shape histogram
      char *histname = new char[100];
      sprintf( histname, "VgPr %i, VgSh %i, Vana %i" , vpr, vsh, iana );
      
        // get pulse shape
      TProfile* phvshld = new TProfile();
      phvshld = (TProfile*) f->Get( "phvshld" );
      phvshld->SetTitle( histname );
      phvshld->SetMinimum( -100 );
      phvshld->SetMaximum(  700 );
      
        // draw for constant VgPr
      if( ish == 0 ){
        cCnstVpr[ipr] = new TCanvas( Form( "cc%iPr%i" , chipNr, vpr), Form( "cc%iPr%i" , chipNr, vpr), 650, 500);
        phvshld->SetMarkerStyle( 20 + ish );
        phvshld->SetMarkerColor(gStyle->GetColorPalette( (float)nc / np* ish ));
        phvshld->SetLineWidth(2);
        phvshld->SetLineColor(gStyle->GetColorPalette( (float)nc / np* ish ));
        phvshld->DrawCopy();
      }
      else{
        cCnstVpr[ipr]->cd();
        phvshld->SetMarkerStyle( 20 + ish );
        phvshld->SetMarkerColor(gStyle->GetColorPalette( (float)nc / np* ish ));
        phvshld->SetLineWidth(2);
        phvshld->SetLineColor(gStyle->GetColorPalette( (float)nc / np* ish ));
        phvshld->DrawCopy("SAME");
      }
        // draw for constant VgSh
      if( ipr == 0 ){
        cCnstVsh[ish] = new TCanvas( Form( "cc%iSh%i" , chipNr, vsh), Form( "cc%iSh%i" , chipNr, vsh), 650, 500);
        phvshld->SetMarkerStyle( 20 + ipr );
        phvshld->SetMarkerColor(gStyle->GetColorPalette( (float)nc / np* ipr ));
        phvshld->SetLineWidth(2);
        phvshld->SetLineColor(gStyle->GetColorPalette( (float)nc / np* ipr ));
        phvshld->DrawCopy();
      }
      else{
        cCnstVsh[ish]->cd();
        phvshld->SetMarkerStyle( 20 + ipr );
        phvshld->SetMarkerColor(gStyle->GetColorPalette( (float)nc / np* ipr ));
        phvshld->SetLineWidth(2);
        phvshld->SetLineColor(gStyle->GetColorPalette( (float)nc / np* ipr ));
        phvshld->DrawCopy("SAME");
      }
      
    }
  } // VgPr and VgSh
  
    // re-iterate to create legends
  for( int ipr = 0; ipr < szVpr; ++ipr ){
    int vpr = bVpr + ipr * sVpr;
    cCnstVpr[ipr]->cd();
    cCnstVpr[ipr]->BuildLegend(0.60,0.60,0.90,0.90, Form(" Constant VgPr %i", vpr ) );
  }
  for( int ish = 0; ish < szVsh; ++ish ){
    int vsh = bVsh + ish * sVsh;
    cCnstVsh[ish]->cd();
    cCnstVsh[ish]->BuildLegend(0.60,0.60,0.90,0.90, Form(" Constant VgSh %i", vsh ) );
  }
  
}
