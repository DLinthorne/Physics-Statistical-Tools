// Include the general helper functions
// saves a lot of typing!
#include "Utils.h"

double binomialProb(int k, double p, int N) {
  if (k>N) return 0;
  return TMath::Binomial(N,k)*pow(p,k)*pow(1.0-p,N-k);
}

TF1 *drawExp(double m, int col=kRed) {
  TF1 *exp = new TF1("","exp(-x/[0])/[0]",0,100);
  // C++ usually start
  exp->SetParameter(0,m);
  // Let's set the line color
  exp->SetLineColor(col);
  exp->Draw("same");

  // calculate the mean and variance by numerical
  // evaluation from integrating the 0-100 interval
  printf("Exponential with m = %.2f has ",m);
  printf("mean = %.2f, variance = %.2f\n",
         exp->Mean(0,100),exp->Variance(0,100));
  return exp;
}

TF1 *drawGauss(double mu, double sigma, int col=kRed) {
  TF1 *gaus = new TF1("","TMath::Gaus(x,[0],[1],1)",-100,100);
  gaus->SetLineColor(col);
  gaus->SetParameters(mu,sigma);
  gaus->Draw("same");
  return gaus;
}

TF1 *drawLogNormal(double mu, double sigma, int col=kRed) {
  TF1 *logNorm = new TF1("","pow(TMath::TwoPi()*[1]*[1],-0.5)*1/x*exp(-pow(log(x)-[0],2)/(2*[1]*[1]))",0,100);
  logNorm->SetLineColor(col);
  logNorm->SetParameters(mu,sigma);
  logNorm->Draw("same");
  return logNorm;
}

TF1 *drawChiQuare(int ndof, int col=kRed, double y_text=0.88) {
  // a chi2 distribution is a Gamma distribution with
  // shape parameter n/2 and normalization parameter =2
  TF1 *chi2 = new TF1("","TMath::GammaDist(x,[0]/2,0,2.0)",0,100);
  chi2->SetLineColor(col); chi2->SetParameter(0,ndof);
  chi2->Draw("same");
  drawText(0.5,y_text,Form("Chi-Square, #it{n}_{dof} = %i",ndof),col,0.05);
  return chi2;
}

void drawPDFs() {
  setStyle();
  
  // create a canvas and split it into 3x2 = 6 parts
  TCanvas *can = new TCanvas();
  can->Divide(3,2);
  
  /************
   * Binomial discrete probability function
   */

  // what N and p should we use?
  vector<int> Ns = {5,10,20,20,20,20};
  vector<double> Ps = {0.5,0.5,0.5,0.1,0.2,0.6};
  for (int ican=1;ican<=6;++ican) {
    // get the current N and p
    int N = Ns[ican-1];
    double p = Ps[ican-1];

    can->cd(ican);
    // defined in Utils.h
    drawAxis(-0.5,0,19.5,0.35,"k");//,"f(k|p,n)");
    for (int i=0;i<20;++i)
      drawVertLine(i,0,binomialProb(i,p,N),kBlue,5);
    drawText(0.4,0.88,"Binomial",kBlue,0.08);
    drawText(0.4,0.8,Form("p = %.2f, N = %i",p,N),kBlue,0.08);
  }
  can->Print("binomial.pdf");
  
  /************
   * Poisson
   */
  vector<double> poisMean = {2,4,12};
  for (int ican=1;ican<=3;++ican) {
    can->cd(ican);
    double nu = poisMean[ican-1];
    drawAxis(-0.5,0,19.5,0.35,"k");
    for (int i=0;i<20;++i)
      drawVertLine(i,0,TMath::Poisson(i,nu),kRed,5);
    drawText(0.4,0.88,Form("Poisson, #nu = %.0f",nu),kRed,0.08);
  }
  can->Print("poisson.pdf");

  /************
   * Exponential p.d.f.
   */
  can->cd();
  drawAxis(0,0,5,1,"x","f(x;m)");
  drawExp(1.0); drawText(0.5,0.88,"m = 1.0",kRed,0.05);
  drawExp(2.0,kBlue); drawText(0.5,0.83,"m = 2.0",kBlue,0.05);
  drawExp(3.0,kGreen+1); drawText(0.5,0.78,"m = 3.0",kGreen+1,0.05);
  can->Print("exponential.pdf");
  
  /************
   * Normal distribution
   */
  drawAxis(-5,0,5,0.5,"x","f(x ; \\mu,\\sigma)");
  drawGauss(0,1,kRed);
  drawText(0.6,0.88,"Gaussian, #mu = 0, #sigma = 1",kRed,0.05);
  drawGauss(-1,1,kBlue);
  drawText(0.6,0.83,"Gaussian, #mu = -1, #sigma = 1",kBlue,0.05);
  drawGauss(-2,2,kGreen+1);
  drawText(0.6,0.78,"Gaussian, #mu = -2, #sigma = 2",kGreen+1,0.05);
  can->Print("gauss.pdf");
  
  // let's copy and paste some code from above, and chagne it a bit
  drawAxis(-0.5,0,19.5,0.35,"x");
  for (int i=0;i<20;++i)
    drawVertLine(i,0,TMath::Poisson(i,5),kRed,5);
  drawGauss(5,sqrt(5.0),kBlue);
  drawText(0.6,0.88,"Poisson, #nu = 5",kRed,0.05);
  drawText(0.6,0.83,"Gaussian, #mu = 5, #sigma = #sqrt{5}",kBlue,0.05);
  can->Print("pois_vs_gauss1.pdf");

  drawAxis(-0.5,0,19.5,0.35,"x");
  for (int i=0;i<20;++i)
    drawVertLine(i,0,TMath::Poisson(i,12),kRed,5);
  drawGauss(12,sqrt(12.0),kBlue);
  drawText(0.6,0.88,"Poisson, #nu = 12",kRed,0.05);
  drawText(0.6,0.83,"Gaussian, #mu = 12, #sigma = #sqrt{12}",kBlue,0.05);
  can->Print("pois_vs_gauss12.pdf");
  
  /************
   * Log-Normal 
   */
  drawAxis(0,0,5,1,"x","f(x ; \\mu,\\sigma)");
  drawLogNormal(0,1,kRed);
  drawText(0.5,0.88,"Log-Normal, #mu = 0, #sigma = 1",kRed,0.05);
  drawLogNormal(1,0.5,kBlue);
  drawText(0.5,0.83,"Log-Normal, #mu = 1, #sigma = 0.5",kBlue,0.05);
  drawLogNormal(1,0.2,kGreen+1);
  drawText(0.5,0.78,"Log-Normal, #mu = 1, #sigma = 0.2",kGreen+1,0.05);

  can->Print("logNorm.pdf");
  
  /************
   * Chi-square
   */
  drawAxis(0,0,15,0.5,"x","f(x ; n)");
  drawChiQuare(1,kRed);
  drawChiQuare(2,kBlue,0.83);
  drawChiQuare(5,kGreen+1,0.78);
  drawChiQuare(8,kOrange+2,0.73);
  can->Print("chiSquare.pdf");
  
}