// make the data globally available
TGraphErrors *data;

double height[6] = {1000,828,800,600,300,0};
double dist[6] = {1500,1340,1328,1172,800,0};

using namespace std;
double linear(int n, double alpha){

	return alpha*height[n];
}
double quad(int n, double alpha, double beta){

	return alpha*height[n] + beta*pow(height[n],2);
}
double power(int n, double alpha, double beta){

	return alpha*pow(height[n], beta);
}
double sqroot(int n, double alpha){
	return alpha*pow(height[n], .5);
}
// calculate the chiSquare quantity assuming a linear function
// with the parameters alpha and beta
double chi2_linear(double alpha) {
  double chi2 = 0;
  for (int i=0;i< 6;i++)
    chi2 += pow( (dist[i] - linear(i, alpha))/15,2);
  return chi2;
}

double chi2_quad(double alpha, double beta) {
  double chi2 = 0;
  for (int i=0;i< 6;i++)
    chi2 += pow( (dist[i] - quad(i, alpha, beta))/15,2);
  return chi2;
}

double chi2_power(double alpha, double beta) {
  double chi2 = 0;
  for (int i=0;i< 6;i++)
    chi2 += pow( (dist[i] - power(i, alpha, beta))/15,2);
  return chi2;
}
double chi2_sqr(double alpha) {
  double chi2 = 0;
  for (int i=0;i< 6;i++)
    chi2 += pow( (dist[i] - sqroot(i, alpha))/15,2);
  return chi2;
}
// Minuit function used by fitter
// Minuit minimize, so this function should return either the chi^2 or -2 ln L
// Error will be given by +/- 1 of this quantity
void minuitChi2lin(int &nDim, double* gout, double& result, double *par, int flg) {
  double alpha = par[0];
  result = chi2_linear(alpha);
}
void minuitChi2(int &nDim, double* gout, double& result, double *par, int flg) {
  double alpha = par[0], beta = par[1];
  result = chi2_quad(alpha,beta);
}
void minuitChi2pow(int &nDim, double* gout, double& result, double *par, int flg) {
  double alpha = par[0], beta = par[1];
  result = chi2_power(alpha,beta);
}
void minuitChi2sqr(int &nDim, double* gout, double& result, double *par, int flg) {
  double alpha = par[0];
  result = chi2_sqr(alpha);
}
void ques1() {
  // true parameters

  TFitter *qfitter = new TFitter(6);
  qfitter->SetFCN(minuitChi2);
  
  double err = 0.0001, low_val = 0.0, high_val = 0.0;
  qfitter->SetParameter(0,"alpha",1,err,low_val,high_val);
  qfitter->SetParameter(1,"beta",1,err,low_val,high_val);
  double arg = 0;
  qfitter->ExecuteCommand("MIGRAD",&arg,0);

  TFitter *lfitter = new TFitter(6);
  lfitter->SetFCN(minuitChi2lin);
  lfitter->SetParameter(0,"alpha",1,err,low_val,high_val);
  lfitter->ExecuteCommand("MIGRAD",&arg,0);

  TFitter *pfitter = new TFitter(6);
  pfitter->SetFCN(minuitChi2pow);
  pfitter->SetParameter(0,"alpha",1,err,low_val,high_val);
  pfitter->SetParameter(1,"beta",1,err,low_val,high_val);
  pfitter->ExecuteCommand("MIGRAD",&arg,0);
  
  TFitter *sfitter = new TFitter(6);
  sfitter->SetFCN(minuitChi2sqr);
  sfitter->SetParameter(0,"alpha",1,err,low_val,high_val);
  sfitter->ExecuteCommand("MIGRAD",&arg,0);
  
  double qalphahat = qfitter->GetParameter(0), qerr_alphaHat = qfitter->GetParError(0);
  double qbetaHat  = qfitter->GetParameter(1), qerr_betaHat  = qfitter->GetParError(1);
  double lalphahat = lfitter->GetParameter(0), lerr_alphaHat = lfitter->GetParError(0);
  double palphahat = pfitter->GetParameter(0), perr_alphaHat = pfitter->GetParError(0);
  double pbetaHat  = pfitter->GetParameter(1), perr_betaHat  = pfitter->GetParError(1);
  double salphahat = sfitter->GetParameter(0), serr_alphaHat = sfitter->GetParError(0);

  
  printf("%.3f\n",chi2_quad(qalphahat,qbetaHat));
  printf("%.3f\n",chi2_linear(lalphahat));
  printf("%.3f\n",chi2_power(palphahat,pbetaHat));
  printf("%.3f\n",chi2_sqr(salphahat));

  cout << TMath::Prob(chi2_quad(qalphahat,qbetaHat), 6 -2) << endl;
  cout << TMath::Prob(chi2_linear(lalphahat), 6 -1) << endl;
  cout << TMath::Prob(chi2_power(palphahat,pbetaHat), 6 -2) << endl;
  cout << TMath::Prob(chi2_sqr(salphahat), 6 -1) << endl;



  }