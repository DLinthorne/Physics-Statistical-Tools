

using namespace std;

double mean[3] = {6.2, 3.1, 7.2};
double stdv[3] = {4.0, 6.0, 7.0};

TVector3 T(mean);

double weight(int i, TMatrixD V){

	double w = 0;
	for (int j = 0; j < 3; j++){

		w += V[i][j];
	}
	return w;
}
double Tweight(int j, TMatrixD V){
double w = 0;
	for (int i = 0; i < 3; i++){

		w += V[i][j];
	}
	return w;

}

void resonance(){

	TString func="[p0]*ROOT::Math::breitwigner_pdf(x,[p1],[p2])";
	func+="+[p3]*exp(-[p4]*x/100)";
	TF1 *s_plus_b = new TF1("",func,70,110);
	s_plus_b->SetParameters(1.0,2.0,90,2.0,4.0);
	//s_plus_b->Draw();

	TH1F * h_chi = new TH1F("","",200,0,40);
	TH1F * h_true = new TH1F("","",200,80,100);
 	TH1F * h_meas = new TH1F("","",200,80,100);
 	TH1F * h_sum  = new TH1F("","",100, -50,50);
 	TH1F * h_uni  = new TH1F("","",100, -4,4);

 	int N=500000;
 	for (int meas = 0; meas<N; ++meas) {
 		double value = s_plus_b->GetRandom(70,110);
 		double detect = value + gRandom->Gaus(0,2.0);
 		h_true->Fill(value);
 		h_meas->Fill(detect);

 		double s = gRandom -> Gaus(0,2 );

 		double x_1 = gRandom->Gaus(-1 2);
 		double x_2 = gRandom->Gaus(-1,2);
 		double x_3 = gRandom->Gaus(-1,2);
 		h_sum -> Fill(4*x_1 - 3*x_2 + x_3);

 		double x_4 = gRandom->Gaus(6.2,4);
 		double x_5 = gRandom->Gaus(3.1,6);
 		double x_6 = gRandom->Gaus(7.2,7);


 		double y_1 = gRandom->Uniform(-1,1);
 		double y_2 = gRandom->Uniform(-1,1);
 		h_uni -> Fill(y_1 + y_2);

 		double chi = .5*((x_1 +1)*(x_1 +1)+ (x_2 +1)*(x_2 +1) +(x_3 +1)*(x_3 +1));
 		h_chi -> Fill((x_4 + x_5 + x_6)/3 + gRandom->Gaus(0,2)/3);

 	}
 	h_true->SetStats(0);
 	h_meas->SetStats(0);

 	//h_true->Draw();
 	//h_meas->SetLineColor(kRed);
 	//h_meas->Draw("same");

 	//leg = new TLegend(0.1,0.7,0.48,0.9);
	//leg->AddEntry(h_true,"Breit-Wigner");
	//leg->AddEntry(h_meas,"Breit-Wigner + detector (#sigma = 2GeV)");

	//leg->Draw();

 	//h_chi -> Draw();
 	h_sum -> Draw("same");

 	// completely correlated
 	TMatrixD V_sym(3,3);
 	TMatrixD V_stat(3,3);
 	TMatrixD V_tot(3,3);
 	TMatrixD V_tot_org(3,3);
 	for (int i = 0; i < 3; i++){

 		for(int j = 0; j < 3; j++){
 			
 			V_stat[i][j] = 4.0;

 			if (i==j) {
 				V_sym[i][j] = 1.0*stdv[i]*stdv[j];
 				V_tot[i][j] = 1.0*stdv[i]*stdv[j] + 4.0;
 			}
 			if (i!=j) V_tot[i][j] = 4.0 ;

 		}
 	}
 	V_tot.Print();
 	V_tot_org = V_tot;
 	V_tot.Invert();
 	V_tot.Print();
 	double V_sum = 0;
 	for (int f = 0; f < 3; f++){
 		for(int t = 0; t < 3; t++){
 			 V_sum += V_tot[f][t];

 		}
 	}
 	V_tot_org.Print();
 	cout << V_sum << endl;

 	//cout << weight(1, V_tot)/V_sum << endl;
 	double lamda = 0;
 	double l_stdv= 0;
 	for (k = 0; k < 3;k++){

 		lamda += (weight(k, V_tot)*mean[k])/V_sum;
 	}



 	double err = 0;
 	for(int i = 0; i < 3; i++){
 		for(int j = 0; j < 3; j++){
 			err += (weight(i,V_tot)/V_sum)*V_tot_org[i][j]*(Tweight(j,V_tot)/V_sum);

 		}
 	} 
 	cout << lamda << " +/- " << sqrt(err) << endl;

 	

}