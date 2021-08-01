////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//	Computational Take home Assignament
//		Construction of Fisher discriminant
//	
//		by Dylan Linthorne
//////////////////////////////////////////////////////////// 
////////////////////////////////////////////////////////////s

//Standard Namespace
using namespace std;

#include <TString.h>
#include <TFile.h>
#include <TTree.h>


void fatal(TString msg) {
	printf("Fatal:\n\n %s\n\n",msg.Data());
	abort();
}

// Fills a Tvector of means from the distros
TVectorD Fill(double values[], int n){

	TVectorD vector(n);

	for (int i = 0; i<n; i++){
		vector(i) = (double) values[i];
	}
	return vector;
}
// The signifigance level function 
double alpha(TH1F* hist, int bin, int max_x, double width ){

	return hist -> Integral(bin, max_x)*width;
}
// The beta function used to calculate the power
double beta(TH1F* hist, int bin, double width){

	return hist -> Integral(0, bin)*width;
}
//function for the signal purity
double purity(double alp, double bet){

	return (1-bet)/((1-bet) + (alp)*10);
}

void project(){

	//opening file
	TFile *file = TFile::Open("dylanlinthorne_training_sets.root");

	if (file==nullptr) fatal("Could not open file");
	// now let’s access the signal and background training data
	TTree *signal_tree = (TTree*)file->Get("signal");
	if (signal_tree==nullptr) fatal("Failed accessing signal tree");
	TTree *bkg_tree = (TTree*)file->Get("background");
	if (bkg_tree==nullptr)	fatal("Failed accessing background tree");

	// define tree vairables and extraction
	double s_a = 0, s_b = 0, s_c = 0;
	signal_tree->SetBranchAddress("a",&s_a);
	signal_tree->SetBranchAddress("b",&s_b);
	signal_tree->SetBranchAddress("c",&s_c);
	
	double b_a = 0, b_b = 0, b_c = 0;
	bkg_tree->SetBranchAddress("a",&b_a);
	bkg_tree->SetBranchAddress("b",&b_b);
	bkg_tree->SetBranchAddress("c",&b_c);
	
	//intitialize histograms 
	TH1F * h_sig_a = new TH1F("signal_a","Signal a",100,-40,40);
	TH1F * h_sig_b = new TH1F("signal_b","Signal b",100,-40,40);
	TH1F * h_sig_c = new TH1F("signal_c","Signal c",100,-40,40);
	
	TH2F * th_sig_ab = new TH2F("2sig_ab","",100,-40,40,100,-40,40);
	TH2F * th_sig_ac = new TH2F("2sig_ac","",100,-40,40,100,-40,40);
	TH2F * th_sig_bc = new TH2F("2sig_bc","",100,-40,40,100,-40,40);

	TH2F * th_bkg_ab = new TH2F("2bkg_ab","",100,-40,40,100,-40,40);
	TH2F * th_bkg_ac = new TH2F("2bkg_ac","",100,-40,40,100,-40,40);
	TH2F * th_bkg_bc = new TH2F("2bkg_bc","",100,-40,40,100,-40,40);

	TH1F * h_bkg_a = new TH1F("bkg_a","bkg a",100,-40,40);
	TH1F * h_bkg_b = new TH1F("bkg_b","bkg b",100,-40,40);
	TH1F * h_bkg_c = new TH1F("bkg_c","bkg c",100,-40,40);

	// number of events for both trees
	int N_signal_events = signal_tree->GetEntries();
	int N_bkg_events = bkg_tree->GetEntries();

	for (int i=0;i<N_signal_events; ++i) {

	// the next line fills the variables a, b and c
	// with observations i
		signal_tree->GetEntry(i);
		bkg_tree->GetEntry(i);

		h_sig_a -> Fill(s_a);
		h_sig_b -> Fill(s_b);
		h_sig_c -> Fill(s_c);

		h_bkg_a -> Fill(b_a);
		h_bkg_b -> Fill(b_b);
		h_bkg_c -> Fill(b_c);

		th_sig_ab -> Fill(s_a,s_b);
		th_sig_ac -> Fill(s_a,s_c);
		th_sig_bc -> Fill(s_b,s_c);

		th_bkg_ab -> Fill(b_a,b_b);
		th_bkg_ac -> Fill(b_a,b_c);
		th_bkg_bc -> Fill(b_b,b_c);


	}

	// create standard array of  means
	double s_mean[3] = { h_sig_a->GetMean(), h_sig_b->GetMean(), h_sig_c->
		GetMean()};
	double b_mean[3] = { h_bkg_a->GetMean(), h_bkg_b->GetMean(), h_bkg_c->
		GetMean()};

	// initalize covariance matrices
	TMatrixD s_V(3,3), b_V(3,3), tot_V(3,3);

	// diagonals of variances
	s_V(0,0) = h_sig_a ->GetStdDev()*h_sig_a ->GetStdDev();
	s_V(1,1) = h_sig_b ->GetStdDev()*h_sig_b ->GetStdDev();
	s_V(2,2) = h_sig_c ->GetStdDev()*h_sig_c ->GetStdDev();

	// Off diagonals using extracted corrleation coeffs
	s_V(0,1) = s_V(1,0) = th_sig_ab -> GetCorrelationFactor()*h_sig_a ->GetStdDev()*
							h_sig_b ->GetStdDev();
	s_V(0,2) = s_V(2,0) = th_sig_ac -> GetCorrelationFactor()*h_sig_a ->GetStdDev()*
							h_sig_c ->GetStdDev();
	s_V(1,2) = s_V(2,1) = th_sig_bc -> GetCorrelationFactor()*h_sig_b ->GetStdDev()*
							h_sig_c ->GetStdDev();
	
	b_V(0,0) = h_bkg_a ->GetStdDev()*h_bkg_a ->GetStdDev();
	b_V(1,1) = h_bkg_b ->GetStdDev()*h_bkg_b ->GetStdDev();
	b_V(2,2) = h_bkg_c ->GetStdDev()*h_bkg_c ->GetStdDev();

	b_V(0,1) = b_V(1,0) = th_bkg_ab -> GetCorrelationFactor()*h_bkg_a ->GetStdDev()*
							h_bkg_b ->GetStdDev();
	b_V(0,2) = b_V(2,0) = th_bkg_ac -> GetCorrelationFactor()*h_bkg_a ->GetStdDev()*
							h_bkg_c ->GetStdDev();
	b_V(1,2) = b_V(2,1) = th_bkg_bc -> GetCorrelationFactor()*h_bkg_b ->GetStdDev()*
							h_bkg_c ->GetStdDev();

	
	// add the covarianve matrices
	tot_V = s_V + b_V;

	s_V.Print();
	b_V.Print();
	TVectorD s_vec = Fill(s_mean, 3);
	TVectorD b_vec = Fill(b_mean, 3);

	TVectorD fisher = s_vec - b_vec;

	// create fisher matrix	
	fisher = tot_V.Invert()*fisher;


	double Nbin = 100; 
	double width = (20 - (-5))/Nbin; 

	TH1F * h_fisher_sig = new TH1F("signal_fisher","Signal a",Nbin,-5,20);
	TH1F * h_fisher_bkg = new TH1F("bkg_fisher","bkg a",Nbin,-5,20);
	for (int j=0;j<N_signal_events; ++j) {
	// the next line fills the variables a, b and c
	// with observations i
		signal_tree->GetEntry(j);
		bkg_tree->GetEntry(j);

		h_fisher_sig -> Fill(fisher[0]*s_a + fisher[1]*s_b +fisher[2]*s_c);
		h_fisher_bkg -> Fill(fisher[0]*b_a + fisher[1]*b_b +fisher[2]*b_c);


	}

	// normalize the histograms
	double area_sig = h_fisher_sig->Integral("width");
	double area_bkg = h_fisher_bkg->Integral("width");

  	h_fisher_sig->Scale(1.0/area_sig);
  	h_fisher_bkg->Scale(1.0/area_bkg);
	
	h_fisher_sig -> SetStats(0);
	h_fisher_bkg -> SetStats(0);

	h_fisher_sig -> Draw();
	h_fisher_bkg ->SetLineColor(kRed);
	h_fisher_bkg -> Draw("same");

	leg = new TLegend(0.1,0.7,0.48,0.9);
	leg->AddEntry(h_fisher_sig,"Signal");
	leg->AddEntry(h_fisher_bkg,"Background");

	leg->Draw();

	TLine *l = new TLine(8,0,8,10);
	l->Draw();



	// algorithm which gets the cut bin with optimal alpha
  	double cut_bin = 0;
  	double alp = 1;
  	for( int n = 0; n < Nbin; n++){
  		
  		if (alp < 0.05)  continue;

  		alp = alpha(h_fisher_bkg, n, Nbin, width);
  		cut_bin = n;

  	}
  	
  	// get beta with degined alpha and t_cut
  	double bet =  beta(h_fisher_sig, cut_bin, width );

  	
  	//open mystery data
	TH1F * h_myst = new TH1F("mystery","mystery",Nbin,-5,20);

  	TFile *file_2 = TFile::Open("dylanlinthorne_mystery_data.root");

  	TTree *mtree = (TTree*)file_2->Get("data");

  	// set branch addresses
  	double m_a = 0, m_b = 0, m_c = 0;
	mtree->SetBranchAddress("a",&m_a);
	mtree->SetBranchAddress("b",&m_b);
	mtree->SetBranchAddress("c",&m_c);
	
	// create mystery histogram in terms of fisher statistic
	int N_events = mtree->GetEntries();	
	for (int i = 0; i < N_events; i++){

		mtree -> GetEntry(i);
		h_myst -> Fill(fisher[0]*m_a + fisher[1]*m_b +fisher[2]*m_c);
	}
	h_myst -> SetStats(0);
	double area_msy = h_myst->Integral("width");

  	h_myst->Scale(1.0/area_msy);


	h_myst -> Draw();

	// fit the mystery data using the guassian fits of the hypothesis
	TF1 *mys = new TF1("f"
		,"([0])*.343*exp(-0.5*(pow((x - 10.21)/1.35,2))) + ([1])*.2917*exp(-0.5*(pow((x - 4.828)/1.87,2)))",-5, 25);
	mys -> SetParameters(.2,.2);
	
	h_myst -> Fit(mys,"r","", -5,25);

	leg = new TLegend(0.1,0.7,0.48,0.9);
	leg->AddEntry(h_myst,"Mystery Data");
	leg->AddEntry(mys,"n(s) f(H(1)) + n(b) f(H(0))");

	leg->Draw();


}
