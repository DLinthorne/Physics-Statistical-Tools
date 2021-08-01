#include "Utils.h"

void nsphere() {

	int Ntrial = 100;	
	int Ndarts = 1000;
	
	vector<double> n_two;
	vector<double> n_three;
	vector<double> n_four;
	vector<double> n_five;

	TH1F *vol = new TH1F("vol","Volume of n-Sphere",200,0,1);

	for( int itrial = 0; itrial < Ntrial; ++itrial){

		int xyhits = 0;
		int xyzhits = 0;
		int xyzwhits = 0;
		int xyzwshits = 0;

		for (int idart = 0; idart < Ndarts; ++idart) {

			double x = gRandom->Uniform(-1,1);
			double y = gRandom->Uniform(-1,1);
			double z = gRandom->Uniform(-1,1);
			double w = gRandom->Uniform(-1,1);
			double s = gRandom->Uniform(-1,1);


			bool xyhit = x*x + y*y < 1;
			if (xyhit) xyhits++;

			bool xyzhit = x*x + y*y + z*z < 1;
			if (xyzhit) xyzhits++;

			bool xyzwhit = x*x + y*y + z*z + w*w < 1;
			if (xyzwhit) xyzwhits++;
	
			bool xyzwshit = x*x + y*y + z*z + w*w + s*s< 1;
			if (xyzwshit) xyzwshits++;

		}
		
		n_two.push_back(xyhits);
		n_three.push_back(xyzhits);
		n_four.push_back(xyzwhits);
		n_five.push_back(xyzwshits);
		
		double value = xyhits;
		double dart  = Ndarts;	
		//cout << value << endl;
		vol -> Fill(value/dart);

	}
	
	double avg_two = accumulate( n_two.begin(), n_two.end(), 0.0)/n_two.size();
	double avg_three = accumulate( n_three.begin(), n_three.end(), 0.0)/n_three.size(); 
	double avg_four = accumulate( n_four.begin(), n_four.end(), 0.0)/n_four.size(); 
	double avg_five = accumulate( n_five.begin(), n_five.end(), 0.0)/n_five.size(); 

	double sum_two = 0;	
	double sum_three = 0;	
	double sum_four = 0;	
	double sum_five = 0;	


	for (double num:n_two) {
		sum_two += (num - avg_two)*(num - avg_two);
	}
	for (double jnum:n_three) {
		sum_three += (jnum - avg_three)*(jnum - avg_three);
	}
	for (double knum:n_four) {
		sum_four += (knum - avg_four)*(knum - avg_four);
	}
	for (double jnum:n_five) {
		sum_five += (jnum - avg_five)*(jnum - avg_five);
	}
	//cout << "dem = 4: " << " eff         = " << avg_five/Ndarts << endl;
	cout << "dem : " << " eff std     = " << (sqrt((avg_five/Ndarts)*(1.- avg_five/Ndarts)/Ndarts)) << endl;
	//cout << "dem = 4: " << " Simulated V = " << 2.0*2.0*2.0*2.0*2.0*avg_five/Ndarts << endl;
	//cout << "         " << " Standard    = " << sqrt(sum_five/n_five.size())/sqrt(Ntrial)<< endl;	
	//out << "         " << " Analytic  V = " << (8./15.)*(3.14*3.14) << endl;
	//cout << "         " << " Analytic  V = " <<  << endl;
	//vol -> Draw();
}
