// A bunch of helper functions
// Dag Gillberg, Sept 2016

void setStyle() {
  double tsize=0.045;
  gStyle->SetTextSize(tsize);
  gStyle->SetLabelSize(tsize,"x"); gStyle->SetTitleSize(tsize,"x");
  gStyle->SetLabelSize(tsize,"y"); gStyle->SetTitleSize(tsize,"y");
  
  gStyle->SetPadLeftMargin(0.1); gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadBottomMargin(0.12); gStyle->SetPadTopMargin(0.04);
  gStyle->SetOptTitle(0); gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetHistMinimumZero();
}

// draw text to the screen
void drawText(double x, double y, TString txt, int col=kBlack,
              double size=0.04) {
  static TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextFont(42); tex->SetTextSize(size);
  tex->SetTextColor(col); tex->DrawLatex(x,y,txt);
}

// Draw a line to the screen
void drawLine(double x1, double y1, double x2, double y2,
              int col=kBlue, int width=2) {
  static TLine *line = new TLine();
  line->SetLineColor(col); line->SetLineWidth(width);
  line->DrawLine(x1,y1,x2,y2);
}

// draw a vertical line
void drawVertLine(double x, double y1, double y2,
                  int col=kBlue, int width=2) {
  drawLine(x,y1,x,y2,col,width);
}

// Draw an empty histogram with axis
TH1F *drawAxis(double x1, double y1, double x2, double y2,
              TString xtitle="", TString ytitle="") {
  TH1F *axis = new TH1F("","",1,x1,x2);
  axis->SetXTitle(xtitle); axis->SetYTitle(ytitle);
  axis->GetYaxis()->SetRangeUser(y1,y2);
  axis->SetStats(0);
  axis->Draw();
  return axis;
}
