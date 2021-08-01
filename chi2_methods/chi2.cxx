#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <TCanvas.h>
#include <TString.h>
#include <TH2F.h>
#include <TMath.h>
#include <TFile.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TF1.h>


using namespace std;

// if something goes wrong, this method will abort the program
void fatal(TString message) {
  cout << endl << "FATAL" << endl << message << endl << endl;
  abort();
}

// access a histogram from a ROOT file
TH1F *getHisto(TFile *f, TString hname) {
  TH1F *h = (TH1F*)f->Get(hname);
  if (h==nullptr) fatal("Cannot access histo "+hname);
  return h;
}

// draw text to the screen
void drawText(double x, double y, TString txt, int col=kBlack,
              double size=0.04) {
  static TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextFont(42); tex->SetTextSize(size);
  tex->SetTextColor(col); tex->DrawLatex(x,y,txt);
}

// calcualte the chi21
double chi2f(TH1F *h_data, TH1F *h_pred) {
  int Nbins=h_data->GetNbinsX();
  double chi2=0;
  for (int bin=1;bin<=Nbins;++bin) {
    double n_obs = h_data->GetBinContent(bin);
    double n_pre = h_pred->GetBinContent(bin);

    chi2 += (n_obs - n_pre)*(n_obs - n_pre)/(n_pre);
  }
  return chi2;
}

int chi2() {
  // open the file
  TFile *f = TFile::Open("Chi2_input.root");
  if (f==nullptr) fatal("cannot open file");
  
  TH1F *h_theory1 = getHisto(f,"theory1");
  TH1F *h_theory2 = getHisto(f,"theory2");
  TH1F *h_data = getHisto(f,"data");
  int Nbins = h_data->GetNbinsX();

  // print out what is in the histogram
  for (int bin=1;bin<=h_data->GetNbinsX();++bin) {
    double Nobs = h_data->GetBinContent(bin);
    printf("Bin %2i has %5.1f observations in data\n",bin,Nobs);
    //cout << "Bin = " << bin << " has " << Nobs << " observations in data." << endl;
  }
 

  // draw the histograms
  h_theory1->SetStats(0);
  h_theory1->Draw();
  h_data->Draw("same");


  double chi2_data_theory1 = chi2f(h_data,h_theory1);
  double p_test_theory1 = TMath::Prob(chi2_data_theory1, Nbins );
  text = Form("#chi^{2}(data,theory1) = %.2f, #it{n}_{dof} = %i",chi2_data_theory1,Nbins);
  
  drawText( 0.4, 0.85,text);

  TString pdf("chi2.pdf");
  gPad->Print(pdf+"["); 
  gPad->Print(pdf);


  h_theory2->SetStats(0);
  h_theory2->Draw();
  h_theory1->Draw("same");
  h_data->Draw("same");


  double chi2_data_theory2 = chi2f(h_data,h_theory2);
  double p_test_theory2 = TMath::Prob(chi2_data_theory2, Nbins);
  TString text1 = Form("#chi^{2}(data,theory2) = %.2f, #it{n}_{dof} = %i",chi2_data_theory2,Nbins);
  TString text2 = Form("P-value (data,theory1) = %.4f",p_test_theory1);
  TString text3 = Form("P-value (data,theory2) = %.4f",p_test_theory2);

  drawText( 0.4, 0.85,text);
  drawText( 0.4, 0.75,text2);
  drawText( 0.4, 0.65,text1);
  drawText( 0.4, 0.55,text3);
  cout << "p_test theory 1 = " << p_test_theory1 << endl;
  cout << "p_test theory 2 = " << p_test_theory2 << endl;


  gPad->Print(pdf);

  
  // PART 2
  // construct pseudo experiments
  gRandom->SetSeed(123456); 
  int Ntrials=10000;
  // creat an empty histogram with pseudo dataa
  TH1F *h_pseudodata = (TH1F*)h_data->Clone();
  TH1F *h_chi2_theory2 = new TH1F("chi1","chi1", 200, 0, 50);
  TH1F *h_p_value_1 = new TH1F("pvalue","pvalue", 200, 0, 1);
  TH1F *h_p_value_2 = new TH1F("pvalue2","pvalue2", 200, 0, 1);


  for (int trial=0;trial<Ntrials;++trial) {

    // Let's construct pseudodata from the null-hypthesis, theory1
    for (int bin=1;bin<=Nbins;++bin) {

      double mean = h_theory2->GetBinContent(bin);
      double nobs = gRandom->Poisson(mean);

      h_pseudodata->SetBinContent(bin,nobs);

    }

    double value = chi2f(h_pseudodata, h_theory1);
    double value_2 = chi2f(h_pseudodata, h_theory2);


    double p_v_1 = TMath::Prob(value, Nbins );
    double p_v_2 = TMath::Prob(value_2, Nbins );

    h_chi2_theory2 -> Fill(value_2);
    h_p_value_1 -> Fill(p_v_1);
    h_p_value_2 -> Fill(p_v_2);


  }

  h_p_value_1->SetStats(0);
  h_p_value_2->SetStats(0);
  h_chi2_theory2->SetStats(0);

  //double binWidth = h_chi2_theory1->GetBinWidth(1);
  double area = h_chi2_theory2->Integral("width");

  h_chi2_theory2->Scale(1.0/area);
  h_chi2_theory2 -> Draw();

  cout << "integral " << h_chi2_theory2->Integral("width") << endl;
  TF1 *fa2 = new TF1("fa2","ROOT::Math::chisquared_pdf(x,20)",0,50);
  fa2 -> Draw("same");

 // h_p_value_1 -> Draw();
  //h_p_value_2 -> Draw("same");

  leg = new TLegend(0.1,0.7,0.48,0.9);
  //leg->SetHeader("The Legend Title"); // option "C" allows to center the header
  leg->AddEntry(h_chi2_theory2,"#chi^{2} Poisson Pseudo-data");
  leg->AddEntry("fa2","#chi^{2} Guassian; NDof = 20");
  //leg->AddEntry(h_p_value_2,"P-value distro: Theory 2");
  //leg -> AddEntry(h_theory1," Theory 1");
  //leg -> AddEntry(h_theory2," Theory 2");
  //leg -> AddEntry(h_data,"Data");

  leg->Draw();


 

  gPad->Print(pdf);
  
  gPad->Print(pdf+"]");
  return 0;
}
