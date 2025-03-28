void pull1D()
{
//=========Macro generated from canvas: cpull1D/
//=========  (Sat Mar 23 20:54:10 2024) by ROOT version 6.24/06
   TCanvas *cpull1D = new TCanvas("cpull1D", "",0,0,500,500);
   cpull1D->Range(-10,-24.9375,10,224.4375);
   cpull1D->SetFillColor(0);
   cpull1D->SetBorderMode(0);
   cpull1D->SetBorderSize(2);
   cpull1D->SetFrameBorderMode(0);
   cpull1D->SetFrameBorderMode(0);
   
   TH1D *hpull1D__3 = new TH1D("hpull1D__3","",40,-8,8);
   hpull1D__3->SetBinContent(6,1);
   hpull1D__3->SetBinContent(7,1);
   hpull1D__3->SetBinContent(8,2);
   hpull1D__3->SetBinContent(9,2);
   hpull1D__3->SetBinContent(10,4);
   hpull1D__3->SetBinContent(11,8);
   hpull1D__3->SetBinContent(12,7);
   hpull1D__3->SetBinContent(13,18);
   hpull1D__3->SetBinContent(14,28);
   hpull1D__3->SetBinContent(15,48);
   hpull1D__3->SetBinContent(16,80);
   hpull1D__3->SetBinContent(17,108);
   hpull1D__3->SetBinContent(18,119);
   hpull1D__3->SetBinContent(19,188);
   hpull1D__3->SetBinContent(20,190);
   hpull1D__3->SetBinContent(21,176);
   hpull1D__3->SetBinContent(22,169);
   hpull1D__3->SetBinContent(23,158);
   hpull1D__3->SetBinContent(24,93);
   hpull1D__3->SetBinContent(25,57);
   hpull1D__3->SetBinContent(26,60);
   hpull1D__3->SetBinContent(27,37);
   hpull1D__3->SetBinContent(28,13);
   hpull1D__3->SetBinContent(29,16);
   hpull1D__3->SetBinContent(30,5);
   hpull1D__3->SetBinContent(31,5);
   hpull1D__3->SetBinContent(32,3);
   hpull1D__3->SetBinContent(33,1);
   hpull1D__3->SetBinContent(34,3);
   hpull1D__3->SetEntries(1600);
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *ptstats_LaTex = ptstats->AddText("hpull1D");
   ptstats_LaTex->SetTextSize(0.0368);
   ptstats_LaTex = ptstats->AddText("Entries = 1600   ");
   ptstats_LaTex = ptstats->AddText("Mean  = 0.01461");
   ptstats_LaTex = ptstats->AddText("Std Dev   =   1.46");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   hpull1D__3->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(hpull1D__3);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   hpull1D__3->SetLineColor(ci);
   hpull1D__3->GetXaxis()->SetLabelFont(42);
   hpull1D__3->GetXaxis()->SetTitleOffset(1);
   hpull1D__3->GetXaxis()->SetTitleFont(42);
   hpull1D__3->GetYaxis()->SetLabelFont(42);
   hpull1D__3->GetYaxis()->SetTitleFont(42);
   hpull1D__3->GetZaxis()->SetLabelFont(42);
   hpull1D__3->GetZaxis()->SetTitleOffset(1);
   hpull1D__3->GetZaxis()->SetTitleFont(42);
   hpull1D__3->Draw("");
   cpull1D->Modified();
   cpull1D->cd();
   cpull1D->SetSelected(cpull1D);
}
