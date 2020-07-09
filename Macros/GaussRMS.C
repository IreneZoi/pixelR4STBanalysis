void GaussRMS(){


  // histogram
  TH1D* hFull = new TH1D("hFull","hFull", 500, -10, 10 );


  // random generator
  gRandom = new TRandom3(0);
  //gRandom->SetSeed(1337); // reproducible


  // fill histogram with 2D gaussian random numbers
  Double_t x = 0;
  for( Int_t i = 0; i < 30000; i++ ){

    x = gRandom->Gaus( 0, 1);
    hFull->Fill( x );

  }


  ///////////////////////////////
  // Self consistent RMS in N RMS
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>

  TH1D* hSubR = (TH1D*)hFull->Clone();
  Double_t N = 4;
  Double_t NFull = hFull->Integral();
  Double_t NSubR = NFull;
  Double_t prev_rms = 0;
  Double_t this_rms = 0;
  Double_t x1 = 0;
  Double_t x2 = 0;
  Int_t cnt = 0;

  do{

    cnt++;
    prev_rms = this_rms;
    this_rms = hSubR->GetRMS();
    x1 = hSubR->GetMean() - this_rms * N;
    x2 = hSubR->GetMean() + this_rms * N;
    NSubR = hSubR->Integral();

    hSubR->GetXaxis()->SetRangeUser(x1,x2);

    cout << "Iterations   " << cnt << endl
	 << "  RMS        " << this_rms << endl
	 << "  Percentage " << NSubR / NFull << endl
	 << "  New Range  " << x1 << " to " << x2 << endl << endl;


  } while( this_rms != prev_rms );

  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<
  // Self consistent RMS in n RMS
  ///////////////////////////////


  // draw
  //TCanvas* cFull = new TCanvas( "cFull", "cFull",   0,   0, 600, 600 );
  //hFull->Draw("");


  return;

} // end main
