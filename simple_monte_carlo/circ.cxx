#include "Utils.h"

void circ() {

	TCanvas *can = new TCanvas("can","can",600,600);
	drawAxis(-1,-1,1,1,"#it{x}","#it{y}");

	// how many darts should we throw?
	int Ndarts = 10000;
	int Nhits = 0;

	for (int idart=0;idart<Ndarts;++idart) {
		// randomly sample our coordinates
		double x = gRandom->Uniform(-1,1);
		double y = gRandom->Uniform(-1,1);


		// is it a hit?
		bool hit = x*x + y*y < 1;
		if (hit) Nhits++;

		// draw our dart!
		static TMarker *dart = new TMarker();
		dart->SetMarkerStyle(27);
		
		if (hit) dart->SetMarkerColor(kRed);
		else dart->SetMarkerColor(kBlue);

		dart->DrawMarker(x,y);
	}
	
	drawText(0.4,0.96,Form("%i hits out of %i",Nhits,Ndarts));
	drawText(0.4,0.92,Form("%.1f%% #rightarrow #it{A} = %.4f",
								100.0*Nhits/Ndarts,4.0*Nhits/Ndarts));
	can->Print(Form("pi_%idarts.pdf",Ndarts));
}
