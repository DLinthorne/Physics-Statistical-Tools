
using namespace std;

TGraphErrors *data;
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

void q1(){

	// true parameters
	double a = 1.0, b=0.5, stdDev = 1.0;
	// 1. Generate the data
	gRandom->SetSeed(12345);
	data = new TGraphErrors();
	for (int x=1;x<=10;++x) {
		// Sample each datapoint at x=1, 2, ... 10
		data->SetPoint(x-1,x,gRandom->Gaus(a+b*x,stdDev));
		data->SetPointError(x-1,0,stdDev);
	}
	data-> Draw();
}