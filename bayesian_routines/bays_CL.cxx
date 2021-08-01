

using namespace std;

double g_esti(int n, double tau, double tau_e){

	double coeff = pow(n,n)/TMath::Factorial(n - 1);
	double ratio = pow(tau_e,n-1)/pow(tau,n);

	return coeff*ratio*exp((-n*tau_e)/tau);
}
double u_inv(int n, double tau, double a){

	return (tau/(2*n))*ROOT::MathMore::chisquared_quantile(1-a, 2*n);
}
double v_inv(int n, double tau, double b){

	return (tau/(2*n))*ROOT::MathMore::chisquared_quantile(b, 2*n);
}

void ques3(){

	TH1F * h_esti_a = new TH1F("","",100,0,7);
	TH1F * h_esti_b = new TH1F("","",100,0,7);
	
	int n = 5;
	double tau_a = 2.5; 
	double tau_b = 1.0;


	double x_a[n], x_b[n];
	int Ntrials=10000;
	for(int trial = 0; trial < Ntrials; trial++){
		double sum_a = 0;
		double sum_b = 0; 

		for(int i = 0; i<n; i++){

			x_a[i] = gRandom-> Exp(tau_a);
			x_b[i] = gRandom-> Exp(tau_b);

			sum_a += x_a[i];
			sum_b += x_b[i];
		}

		h_esti_a -> Fill(sum_a/n);
		h_esti_b -> Fill(sum_b/n);
	}

	double area_a = h_esti_a->Integral("width");
	double area_b = h_esti_b->Integral("width");

  	h_esti_a->Scale(1.0/area_a);
  	h_esti_b->Scale(1.0/area_b);

	h_esti_b -> Draw();
	h_esti_a -> Draw("same");

	TF1 *g_a = new TF1("","g_esti([0],[1],x)",0,7);
	TF1 *g_b = new TF1("","g_esti([0],[1],x)",0,7);

	g_a -> SetParameters(n, tau_a);
	g_b -> SetParameters(n, tau_b);

	g_a -> Draw("same");
	g_b -> Draw("same");

	TF1 *u_95 = new TF1("","u_inv([0],x,[1])",0,5);
	u_95 -> SetParameters(5,.025);
	TF1 *v_95 = new TF1("","v_inv([0],x,[1])",0,5);
	v_95 -> SetParameters(5,.025);

	TF1 *u_68 = new TF1("","u_inv([0],x,[1])",0,5);
	u_68 -> SetParameters(5,.1585);
	TF1 *v_68 = new TF1("","v_inv([0],x,[1])",0,5);
	v_68 -> SetParameters(5,.1585);


	//u_68 -> Draw();
	//v_68 -> Draw("same");
	//u_95 -> Draw("same");
	//v_95 -> Draw("same");

	leg = new TLegend(0.1,0.7,0.48,0.9);
  	//leg->SetHeader("The Legend Title"); // option "C" allows to center the header
  	//leg->AddEntry(u_68,"65% Confidence Interval U");
  	//leg->AddEntry(u_95,"95% Confidence Interval U");
	//leg->AddEntry(v_68,"65% Confidence Interval V");
  	//leg->AddEntry(v_95,"95% Confidence Interval V");


  	leg-> AddEntry(h_esti_a,"Sim Distro for #tau = 2.5 ");
	leg-> AddEntry(h_esti_b,"Sim Distro for #tau = 1.0");
	leg-> AddEntry(g_a,"Expected Distro for #tau = 2.5 ");
	leg-> AddEntry(g_b,"Expected Distro for #tau = 2.5");

	leg -> Draw();

	cout << u_68->GetX(1.2)  << endl;
	cout << v_68->GetX(1.2) << endl;

}