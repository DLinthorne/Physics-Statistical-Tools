#include "Utils.h"

void fatal(TString msg) {
printf("Fatal:\n\n %s\n\n",msg.Data());
abort();
}

void Q4_2() {


TH1D * h = new TH1D("", ";Energy(GeV);Density", 1000,170,555);
TH1D * ht = new TH1D("", ";t;Density", 1000,0,555);



TFile *file = TFile::Open("particleX_search.root");
if (file==nullptr)
  fatal("Could not open file");

TTree *tree = (TTree*)file->Get("ee_tree");
if (tree==nullptr)
  fatal("Failed accessing tree");

float el_pt = 0,  el_eta = 0, el_phi = 0, el_m = 0;
float pos_pt = 0, pos_eta = 0, pos_phi = 0, pos_m = 0;
tree->SetBranchAddress("electron_pt",&el_pt);
tree->SetBranchAddress("electron_eta",&el_eta);
tree->SetBranchAddress("electron_phi",&el_phi);
tree->SetBranchAddress("electron_m",&el_m);
tree->SetBranchAddress("positron_pt",&pos_pt);
tree->SetBranchAddress("positron_eta",&pos_eta);
tree->SetBranchAddress("positron_phi",&pos_phi);
tree->SetBranchAddress("positron_m",&pos_m);


double_t pex, pey, pez, Ee, ppx, ppy, ppz, Ep, Etot, Mx, pepp;

  TVector3 p_e, p_p, tot_p;
  TLorentzVector Lp_e, Lp_p, Ltot_p;


  int Nevents = tree->GetEntries();
  printf("Will loop over %i events\n",Nevents);
  for (int ievent = 0; ievent<Nevents;++ievent) {

  tree->GetEntry(ievent);




  //Calculating Momentum of Electron
  pex = el_pt*cos(el_phi);
  pey = el_pt*sin(el_phi);
  pez = el_pt*sinh(el_eta);

  //Calculating Momentum of Positron
  ppx = pos_pt*cos(pos_phi);
  ppy = pos_pt*sin(pos_phi);
  ppz = pos_pt*sinh(pos_eta);


  //Momentum three vector of Electron
  p_e.SetXYZ(pex,pey,pez);
  Ee = p_e.Mag();

  //Momentum three vector of Positron
  p_p.SetXYZ(ppx,ppy,ppz);
  Ep= p_p.Mag();


  //Didn't use these. Just found them online and thought they might be helpful.
  Lp_p.SetXYZM(ppx,ppy,ppz,Ep);
  Lp_e.SetXYZM(pex,pey,pez,Ee);


  //Dot product of the two Momentum vectors
  pepp = (pex*ppx + pey*ppy + pez*ppz);


  //Calculating Particle X mass
  Mx = pow(el_m,2) + pow(pos_m,2)  +  2*(Ee*Ep - pepp);
  Mx = sqrt(Mx);



  h->Fill(Mx);

}


  data1 = new TGraph();




  TF1 *f2 = new TF1("f2", "gaus");
  //double_t P2 = f2->GetProb();


  TF1 *f1 = new TF1("f2", "[0]*pow(x,[1])");
  f1->SetParameter(0,1);
  f1->SetParameter(1,-3);




double_t sig, bkg, t;
for (int i=175;i<=545;){


  h->Fit(f1,"r","",i,i+5);
  h->Fit(f2,"r","",i,i+5);



  bkg = f1->Integral(i,i+5);
  sig = f2->Integral(i,i+5);


  t = bkg/sig;

  data1->SetPoint((i-175)/5,i,t);

  i = i+5;
}


//c2 = new TCanvas("c2","",200,10,700,500);
//c2->SetLogy();

gPad->SetLogy(1);


data1->SetMarkerStyle(20);
data1->GetYaxis()->SetRangeUser(0.1,10);
data1->GetHistogram()->SetTitle("Test Statistics for Particle Discovery");
data1->GetYaxis()->SetTitle("Background/Signal");
data1->GetXaxis()->SetTitle("T (GeV)");
data1->Draw("ALP");




}
