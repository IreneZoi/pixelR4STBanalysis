

// Daniel Pitzl, Jul 2012
// fit generalized error + p0 at mid
// .x fitep0.C("h012")

//----------------------------------------------------------------------
// generalized error function:
Double_t ep0Fit( Double_t *x, Double_t *par )
{
  static int nn = 0;
  nn++;
  static double dx = 0.1;
  static double b1 = 0;
  if( nn == 1 ) b1 = x[0];
  if( nn == 2 ) dx = x[0] - b1;

  // Mean and width:

  double bet = par[2]; // exponent

  double pk = 0;

  if( bet > 0 ) {

    double xm  = par[0];
    double sig = par[1];

    double t = ( x[0] - xm ) / sig;

    double u = pow( fabs( t/sqrt(2.0) ), bet );

    // Normalization needs Gamma function:

    double aa = dx / sig / sqrt(8.0) * bet / TMath::Gamma( 1/bet );

    pk = par[3] * aa * exp( -u );
  }

  return pk + par[4];
}

//----------------------------------------------------------------------
void fitep0( string hs, double x1=1, double x9=0 )
{
  TH1 *h = (TH1*)gDirectory->Get( hs.c_str() );

  if( h == NULL ) {
    cout << hs << " does not exist\n";
    return 1;
  }

  h->SetMarkerStyle(21);
  h->SetMarkerSize(0.8);
  h->SetStats(1);
  gStyle->SetOptFit(101);

  gROOT->ForceStyle();

  int nb = h->GetNbinsX();

  double xl = h->GetBinCenter(1); // first
  double xr = h->GetBinCenter(nb); // last

  if( x9 < x1 ) {
    x1 = xl;
    x9 = xr;
  }

  // protect against out-of-bounds:

  if( x1 < xl ) x1 = xl; // left
  if( x1 > xr ) x1 = xl; // left

  if( x9 > xr ) x9 = xr; // right
  if( x9 < xl ) x9 = xr; // right

  cout << hs << ": " << x1 << " - " << x9 << endl;

  // find mid:

  double xmid = 0.5*(x1+x9);
  double imid = h->FindBin(xmid);
  double ymid = h->GetBinContent(imid);

  int i1 = h->FindBin(x1);
  int i9 = h->FindBin(x9);
  double y1 = h->GetBinContent(i1);
  double y9 = h->GetBinContent(i9);

  double bg = 0.5*(y1+y9);
  double slp = (y9-y1)/(x9-x1);

  double dx = h->GetBinWidth(imid);
  double sm = 2*dx;
  double aa = 2.5 * ( ymid - bg ) * sm / dx; // Gaussian normalization

  // create a TF1 with the range from x1 to x9 and 5 parameters

  TF1 *ep0Fcn = new TF1( "ep0Fcn", ep0Fit, x1, x9, 5 );

  ep0Fcn->SetParName( 0, "mean" );
  ep0Fcn->SetParName( 1, "sigma" );
  ep0Fcn->SetParName( 2, "pow" );
  ep0Fcn->SetParName( 3, "area" );
  ep0Fcn->SetParName( 4, "BG" );
 
  // set start values for some parameters:

  cout << hs << ": fit from " << x1 << " to " << x9 << endl;
  cout << "start mid " << ymid << " at " << xmid << endl;
  cout << "start sm " << sm << endl;
  cout << "start BG " << bg << endl;
  cout << "start area " << aa << endl;

  ep0Fcn->SetParameter( 0, xmid ); // peak position
  ep0Fcn->SetParameter( 1, sm ); // width
  ep0Fcn->SetParameter( 2, 3.3 ); // pow
  ep0Fcn->SetParameter( 3, aa ); // area
  ep0Fcn->SetParameter( 4, bg );

  ep0Fcn->SetNpx(500);
  ep0Fcn->SetLineWidth(4);
  ep0Fcn->SetLineColor(kMagenta);
  ep0Fcn->SetLineColor(kGreen);

  h->Fit( "ep0Fcn", "R", "ep" );
  // h->Fit("ep0Fcn","V+","ep");

  h->Draw("histepsame");  // data again on top

  xmid = ep0Fcn->GetParameter(0);
  bg = ep0Fcn->GetParameter(4);
  cout << "height " << ep0Fcn->Eval(xmid) - bg << endl;

}
