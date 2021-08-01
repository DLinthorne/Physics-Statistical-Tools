// make the data globally available
TGraphErrors *data;

using namespace std;
// calculate the chiSquare quantity assuming a linear function
// with the parameters alpha and beta
double chi2_linear(double alpha, double beta) {
  double chi2 = 0;
  for (int i=0;i<data->GetN();++i)
    chi2 += pow( data->GetY()[i] - alpha - beta*data->GetX()[i],2)/pow(data->GetEY()[i],2);
  return chi2;
}

// Minuit function used by fitter
// Minuit minimize, so this function should return either the chi^2 or -2 ln L
// Error will be given by +/- 1 of this quantity
void minuitChi2(int &nDim, double* gout, double& result, double *par, int flg) {
  double alpha = par[0], beta = par[1];
  result = chi2_linear(alpha,beta);
}

void root_fit_example() {
  // true parameters
  double a = 1.0, b=0.5, stdDev = 1.0;
  
  // 1. Generate the data
  gRandom->SetSeed(12345);
  data = new TGraphErrors();
  for (int i=1;i<=10;++i) {
    data->SetPoint(i-1,i,gRandom->Gaus(a+b*i,stdDev));
    data->SetPointError(i-1,0,stdDev);
  }
  data->SetMarkerStyle(20);
  //data->Draw("AP");
  //gPad->Print("toy_data.pdf");
  
  TFitter *fitter = new TFitter(10);
  fitter->SetFCN(minuitChi2);
  
  // Set the fit starting values
  // Declare the parameter:
  // index, name, startign value, error, low and high range
  double err = 0.0001, low_val = 0.0, high_val = 0.0;
  fitter->SetParameter(0,"alpha",1,err,low_val,high_val);
  fitter->SetParameter(1,"beta",1,err,low_val,high_val);
  
  double arg = 0;
  fitter->ExecuteCommand("MIGRAD",&arg,0);

  double alphaHat = fitter->GetParameter(0), err_alphaHat = fitter->GetParError(0);
  double betaHat  = fitter->GetParameter(1), err_betaHat  = fitter->GetParError(1);
  
  printf("alpha: %.3f +/- %.3f\n",alphaHat,err_alphaHat);
  printf("beta:  %.3f +/- %.3f\n",betaHat,err_betaHat);
  
  printf("%.3f\n",chi2_linear(alphaHat,betaHat));
  printf("%.3f\n",chi2_linear(alphaHat - err_alphaHat,betaHat + err_betaHat));

  TF1  *LS = new TF1("chi", "chi2_linear([0],x)", -1, 4);
  LS -> SetParameter(0,alphaHat + err_alphaHat);
  printf("%.3f\n", LS -> GetMinimum());
  printf("%.3f\n", LS -> GetMinimumX());

  //LS -> Draw();

  fitter->GetMinuit()->SetErrorDef(5.99);
  TGraph *conf_95 = (TGraph*)fitter->GetMinuit()->Contour(120,0,1);
  conf_95->SetFillColor(kGray);
  fitter->GetMinuit()->SetErrorDef(2.3);
  TGraph *conf_68 = (TGraph*)fitter->GetMinuit()->Contour(80,0,1);
  conf_68->SetFillColor(kAzure);

  conf_95->Draw("af"); // A: draw with axis and F: fill
  conf_68->Draw("f"); // only fill "f" (draw on top, i.e. no axis)

  TMarker *m = new TMarker();
  m->SetMarkerStyle(29);
  m->SetMarkerSize(1.4);
  m->DrawMarker( a,b);

  leg = new TLegend(0.1,0.7,0.48,0.9);
  //leg->SetHeader("The Legend Title"); // option "C" allows to center the header
  leg->AddEntry(conf_68,"65% Confidence Interval");
  leg->AddEntry(conf_95,"95% Confidence Interval");
  leg->AddEntry(m,"(#alpha, #beta) = (1.0, 5.0)");

  leg-> Draw();
}
