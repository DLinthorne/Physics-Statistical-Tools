#include <iostream>
#include <vector>
#include <cmath>

#include <TCanvas.h>
#include <TString.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TH2F.h>
#include <TF1.h>
#include <TError.h>

using namespace std;

// define a global variable
TString g_pdf("Central_Limit_Theorem.pdf");

void fatal(TString msg) { 
  cout << endl << "FATAL" << endl << "   " << msg << endl << endl;
  abort();
}
// draw text to the screen
void drawText(double x, double y, TString txt, int col=kBlack,
	      double size=0.04) {
  static TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextFont(42); tex->SetTextSize(size);
  tex->SetTextColor(col); tex->DrawLatex(x,y,txt);
}

/*
 * This is where the bulk of the work gets done!
 * See below
 */
void centralLimitTheorem(int nvars, int ntrials);

int main() {
  gErrorIgnoreLevel = 2001;
  int seed = 123456;
  cout << "Using random seed " << seed << endl;
  gRandom->SetSeed(seed);


  // create a new Canvas and open the pdf for writing
  new TCanvas();
  gPad->Print(g_pdf+"[");

  // How many trials should we do ?
  int ntrials = 100000;

  // How many random variables shoudl add ?
  for (int nvar : {1,2,10} ) {

    // produce results for the current
    // number of varaibles and trials
    centralLimitTheorem(nvar,ntrials);

  }

  // close the pdf for writing
  gPad->Print(g_pdf+"]");
  cout << "Produced the file " << g_pdf << endl;
}


/*****
 * IMPLEMENTAION
 */
void centralLimitTheorem(int nvars, int ntrials) {
  cout << endl << "Will add up " << nvars << " random variables " << ntrials << " times" << endl;

  /*
   *  1. Draw the sum of nvar random variables
   */
  
  TH1F *x_hist = new TH1F("",";x",100,-1,4);
  TH1F *y_hist = new TH1F("","",101,-0.5,0.5+nvars);
  TString yTitle="Sum of "; yTitle+=nvars; yTitle+=" variables";
  y_hist->SetXTitle(yTitle);
  // 
  TH1F *y_pull = new TH1F("",";(y - #mu_{y})/#sigma_{y}",100,-5,5);
  vector<double> y_values;

  // loop over the trials
  for (int i_trial=0;i_trial<ntrials;++i_trial) {
    // what is the sum for the given trial?
    double y = 0;
    for (int ivar=0;ivar<nvars;++ivar) {
      double x = gRandom->Uniform(0,1);
      //double x = gRandom->Exp(0.6);
      y += x;
      x_hist->Fill(x);
    }
    y_hist -> Fill(y);
    y_values.push_back(y);
  }

  double p1 = y_hist->GetMean();
  double p2 = y_hist->GetRMS();
  double p0 = 1;

  for (int i=0; i<ntrials; i++){
    y_pull -> Fill((y_values[i] - p1)/p2);
    cout << (y_values[i] - p1)/p2 << endl;	
  }


  double mu_x = x_hist->GetMean(), sigma_x = x_hist->GetRMS();
  cout << "mu_x = " << mu_x << ", sigma_x = " << sigma_x << endl;
  cout << "mu_y = " << y_hist->GetMean() << ", sigma_y = " << y_hist->GetRMS() << endl;
    

  // distribution of x
  x_hist->Draw();
  gPad->Print(g_pdf);

  // distribution of y = sum(x)
  //y_hist->Draw();

  y_pull->Draw();

  drawText(0.5,0.8,Form("y = #sum_{i=0}^{%i} x_{i}",nvars));
  //  drawText(0.5,0.7,Form("Sum of %i random variables",nvars));
  gPad->Print(g_pdf);
  
  /*
   *  2. Draw (y - E(y))/STD(y)
   */ 

  //TF1 gaus("","[p0]*TMath::Gaus(x,[p1],[p2],1)",0,nvars); 
  // what is the bin width?
  //double binWidth = y_hist->GetBinWidth(1);
  //gaus.SetParameters(p0,p1,p2);
  //gaus.Draw("same");


}


