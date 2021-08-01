
TH1F *m_ee;
TTree *tree;
TGraphErrors *data;
double m_0 = 200;

double guass_A(double m){

	return .23*exp(-0.5*(pow((m - m_0)/3.0,2)));
}
double guass_B(double m){

	return exp(-0.5*(pow((m - m_0)/8.0,2)));
}


void ques4(){

	TFile *file = TFile::Open("particleX_search.root");

	tree = (TTree*)file->Get("ee_tree");
	

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

	int Nevents = tree->GetEntries();
	printf("Will loop over %i events\n",Nevents);

	double px_e, py_e, pz_e, px_p
		 , py_p, pz_p, m_2, p_e_2, p_p_2, p_ep;

	m_ee = new TH1F("invar mass","mass",375,175,550);
	TH1F * testh = new TH1F("tes","tes",200,0,50);
	TH1F * hyp_A = new TH1F("hyp a","mass",200,170,560);



	for (int ievent = 0; ievent<Nevents;++ievent) {
  		tree->GetEntry(ievent);
  		if (ievent<10) {
    	
    	printf("Event %5i:\n",ievent);
    	printf("  The electron has (pT,eta,phi,m)=(%.1f GeV,%5.2f,%5.2f,%4.1f GeV)\n",el_pt,el_eta,el_phi,el_m);
  		}

  		px_e = el_pt*cos(el_phi);
  		px_p = pos_pt*cos(pos_phi);
    	
    	py_e = el_pt*sin(el_phi); 
    	py_p = pos_pt*sin(pos_phi);
    	
    	pz_e = el_pt*sinh(el_eta);
    	pz_p = pos_pt*sinh(pos_eta);
    	
    	p_e_2 = px_e*px_e + py_e*py_e + pz_e*pz_e;
    	p_p_2 = px_p*px_p + py_p*py_p + pz_p*pz_p;
    	p_ep = px_e*px_p + py_e *py_p + pz_e*pz_p;

    	m_2 = el_m*el_m + pos_m*pos_m + 2*(sqrt(p_p_2)*sqrt(p_e_2) - p_ep );

  		m_ee -> Fill(sqrt(m_2));
 
	} 

	double area_a = m_ee->Integral("width");

  	m_ee->Scale(1.0/area_a);

	//m_ee -> Fit("f_bkg");
  	data1 = new TGraph();
  	data2 = new TGraph();

  	TF1 *f2_a = new TF1("f2", "[0]*guass_B(x) + [1]*pow(x,[2])", 175,550);

	//TF1 *f2_a = new TF1("f2", "[0]*guass_B(x) + 1.08533*pow(x,-3.08)", 175, 550);
  	f2_a->SetParameters(1.2,1,-3);
  	//f2_b->SetParameter(0,1.2);

  	TF1 *f1 = new TF1("f2", "[0]*pow(x,[1])", 175, 550);

  	f1->SetParameter(0,1);
  	f1->SetParameter(1,-3);
  	m_ee->Fit(f1,"r","",175,550);


	double_t sig, sig_b, bkg, t, t_b;
	for (int i=200;i<=550;){

		m_0 = i;
  		
    //m_ee->Fit(f1,"r","", 175, 550);
 		m_ee->Fit(f2_a,"r","",175, 550);
 		//m_ee->Fit(f2_b,"r","",175,550);

      double p;
      t = -2*TMath::Log(f1->Eval(i)/f2_a->Eval(i));

 		 p = 1 - ROOT::Math::chisquared_cdf(t,2);


  		//t = -2*TMath::Log(bkg/sig);
  		//t = bkg/sig;


  		//t_b = bkg/sig_b; 

  		data1->SetPoint((i-200)/5,i,p);
  		//data2->SetPoint((i-200)/5,i,t_b);

  		i = i+5;
  	}
gPad->SetLogy(1);
data1->SetMarkerStyle(20);
data1->GetYaxis()->SetRangeUser(0.1,10);
data1->GetHistogram()->SetTitle("Test Statistics for Particle Discovery");
data1->GetYaxis()->SetTitle("Background/Signal");
data1->GetXaxis()->SetTitle("T (GeV)");
//testh -> Draw();
data1->Draw("ALP");
//data2->Draw("* SAME");


}