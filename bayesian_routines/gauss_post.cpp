double GaussPosterior(double x) {
 // p(nu|n) = L(n|nu)*pi(nu) / \int L(n|nu') * pi(nu') dnu'
 // Flat prior:
 // p(nu|n) = L(n|nu) / \int L(n|nu') dnu'
 static TF1 *pois = new TF1("","ROOT::Math::gaussian_pdf(x, [p0], [p1])",-100,100);
 pois->SetParameters(3,2,14.1);

 // integrate from 0 -> expectation value + 5 sigma
 return ROOT::Math::gaussian_pdf(x,3.2, 14.1)/pois->Integral(-100,100.0);
}
void ques2() {
 int n=5;
 int Nsteps=1000;
 double max=50;
 // build posterior. 200 ppoints between 0 and 20
 //TGraph *posterior = new TGraph();
 //for (int step=0;step<=Nsteps;++step) {
 //double x=max*step/Nsteps;
 //posterior->SetPoint(step,x,GaussPosterior(x));
 //}
 //TH1F *axis = new TH1F("",";#nu;probability density",1,0,20);
 //axis->SetStats(0); axis->SetMaximum(0.25); axis->Draw();

TF1 *gaus = new TF1("Guass", "GaussPosterior(x)",0,40);
gaus-> Draw();

 //posterior->SetLineColor(kBlue);
 //posterior->Draw();
 //TF1 *gaus = new TF1("","ROOT::Math::gaussian_pdf(x,[p0],[p1])",0,100);
 //gaus->SetParameters(sqrt(n),n);
 //gaus->SetLineColor(kRed);
 //gaus->Draw("same");

 //TLatex *tex = new TLatex();
 //tex->SetTextFont(42); tex->SetTextSize(0.04); tex->SetNDC();
 //tex->DrawLatex(0.4,0.85,Form("Baysian posterior, Possion with #it{n}_{obs} = %i",n));
 //tex->DrawLatex(0.4,0.80,Form("ML #hat{#nu} given #it{n}_{obs} = %i, assume Gauss. uncert.",n));
 //gPad->Print("bayes_vs_ml.pdf");
}