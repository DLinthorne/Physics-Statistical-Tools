// This exercise is about calculating the arithmetic mean,
// the variance and the standard deviation of two datasets
// First, let's create methods for these
double mean(vector<double> data);
double variance(vector<double> data);
double standardDeviation(vector<double> data);
double covariance(vector<double> data_h, vector<double> data_j);
double correlation(vector<double> data_h, vector<double> data_j);

void VarianceAndMean() {
//  // 1. Read in the data from on of the file and put it in a vector
	vector<double> data_H;
	vector<double> data_J;

	ifstream infile_H;
	ifstream infile_J;

        infile_H.open("H_pT.txt");
  	infile_J.open("jet1_pT.txt");
	
        double value;
        if (infile_H.is_open())
		while (infile_H>>value)
        		data_H.push_back(value);
	double nvalue;
        if (infile_J.is_open())
		while (infile_J>>nvalue)
        		data_J.push_back(nvalue);
	
	infile_H.close();
	infile_J.close();

	int Ndata_H = data_H.size();
	int Ndata_J = data_J.size();

        cout << "Read " << Ndata_H << " numbers from the file!" << endl; 
	
	//Create a canvas object where the histogram will be drawn
	TCanvas can;

	TH1F *h = new TH1F("h","Higgs boson transverse momentum",50,0,200);
	TH1F *j = new TH1F("j","Parton Jet transverse momentum",50,0,200); 
	// increase the range of Pt to include all event data (no overflow)
	TH2F *hj = new TH2F("hj","Tansverse momentum of Higgs vs parton jet",50,0,400,50,0,400);

	for (double num:data_H) {
 		h->Fill(num);
 	}
	for (double jnum:data_J) {
		j->Fill(jnum);
	}
	
	int i = 0;
	while(i < Ndata_H){
	
		hj->Fill(data_H[i],data_J[i]);
		i++;
	}
			

	//h->Draw();
	hj->Draw("SAME");
	//can.Print("GG_H.pdf");
	//can.Print("Jet_H.pdf");
	can.Print("both.pdf");	
	
// 2. Calculate the mean and the variance
cout << "The mean of the dataset is " << mean(data_J) << endl;
cout << "The variance of the dataset is " << variance(data_H) << endl;
cout << "The standard deviation is " << standardDeviation(data_H) << endl;
cout << "The covariance of the datasets is " << covariance(data_H,data_J) << endl;
cout << "The correlation coeff of the datasets is " 
		<< correlation(data_H,data_J) << endl;
cout << "The correlation coeff from the histogram is " << hj -> GetCorrelationFactor() << endl;
}
double mean(vector<double> data) {
	int N = data.size();

	double sum = 0.0;
 	for ( double num : data ) {
 		sum = sum + num;
 	} 
	return sum/N;
	
}
double variance(vector<double> data) {
	int N = data.size();

        double sum = 0.0;
        for ( double num : data ) {
                sum = sum + (num - mean(data))*(num - mean(data));
        }
        return sum/N;
}
double standardDeviation(vector<double> data) {
	return sqrt(variance(data));
}
double covariance(vector<double> data_h, vector<double> data_j) {
	
	int N = data_h.size();

        double sum = 0.0;
        for( int i = 0; i < N; i++){
	
	sum = sum + (data_h[i] - mean(data_h))*(data_j[i] - mean(data_j));

	}
        return sum/N;

}
double correlation(vector<double> data_h, vector<double> data_j){

	return covariance(data_h,data_j)/(standardDeviation(data_h)*standardDeviation(data_j));

}
