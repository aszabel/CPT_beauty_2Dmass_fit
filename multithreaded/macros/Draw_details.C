#include "D_M_fit_shape.h"

using namespace cpt_b0_analysis;
   double minx = 1800.;
   double maxx = 1940.;
   double miny = 2700.;
   double maxy = 8300.;
   const int nvar_md = 7;
   const int nvar_mb = 7;
   const int ncontr = 6;
   const int nbins = 30;

void Draw_pull(TCanvas *c, TH1D *hist, TF1 *func[], TF1 *tf1_sum);

double chebyshev(const double *x, const double *par){
  double m_rec = x[0];
  double m_corr = x[1];
  double a1 = par[0];
  double a2 = par[1];
  double cheb = 1.0 + a1*m_rec + a2*(2.0*m_rec*m_rec-1.0);
  double int_cheb = (1.0-a2)*maxx + 0.5*a1*maxx*maxx+2./3.*a2*maxx*maxx*maxx- (1.0-a2)*minx-0.5*a1*minx*minx-2./3.*a2*minx*minx*minx;
  cheb/=int_cheb;
  return cheb;
}

void Draw_details(int sign){

	TH2D *hist2D = new TH2D("h2D", "", nbins, minx, maxx, nbins, miny, maxy);
	TH1D *histMD = new TH1D("hMD", "", nbins, minx, maxx);
	TH1D *histMB = new TH1D("hMB", "", nbins, miny, maxy);
       TChain ch("BlindedTree");
   ch.Add("/home/szabelskia/LHCb/data2016/tree_missPT_D_M_MagDown25102023_nomassDmuCut/selected_data2016MagDown.root");
     double D_M, mu_PT, mu_P, mu_eta, K_PT, B_M, missPT;
   bool charge;


   std::vector<std::pair<double, double>> vect_2D;
   ch.SetBranchAddress("B_M", &B_M);
   ch.SetBranchAddress("missPT", &missPT);
   ch.SetBranchAddress("D_M", &D_M);
   ch.SetBranchAddress("mu_PT", &mu_PT);
   ch.SetBranchAddress("mu_P", &mu_P);
   ch.SetBranchAddress("mu_eta", &mu_eta);
   ch.SetBranchAddress("K_PT", &K_PT);
   ch.SetBranchAddress("truecharge", &charge);


   for (int i=0; i<ch.GetEntries(); ++i){
   //for (int i=0; i<5.0e5; ++i){
      ch.GetEntry(i);
      if (mu_PT<500 || mu_P<5000 || mu_eta<2 || mu_eta>4.5) continue;
      double B_MMcorr = B_M +2.0*missPT;
      if (B_MMcorr<miny || B_MMcorr>maxy) continue;
      if (D_M<minx || D_M>maxx) continue;
      if (int(charge) == sign){
              hist2D->Fill(D_M, B_MMcorr);
	      histMD->Fill(D_M);
	      histMB->Fill(B_MMcorr);
      }
   }
        double  nevents = double(hist2D->Integral());
 
   	double res[(nvar_md+nvar_mb)*ncontr+ncontr-1];
	ifstream input(Form("results_%d.txt", sign));
	double x, dx;
	int i=0;
	while (input>>x>>dx){res[i] = x; i++;}
	input.close();


	md_fit md_shape(minx, maxx);
        mb_2misspt_fit mb_shape(miny, maxy);
        TF2 *func2D[ncontr];
	TF1 *funcMD[ncontr], *funcMB[ncontr];
        double frac[ncontr];
        double *pa = &res[ncontr*(nvar_md+nvar_mb)];
                frac[0] = 1-abs(pa[0]);
                frac[1] = abs(pa[0])*(1.0-abs(pa[1]));
                frac[2] = abs(pa[0])*abs(pa[1])*(1.0-abs(pa[2]));
                frac[3] = abs(pa[0])*abs(pa[1])*abs(pa[2])*(1.0-abs(pa[3]));
                frac[4] = abs(pa[0])*abs(pa[1])*abs(pa[2])*abs(pa[3])*(1.0-abs(pa[4]));
                frac[5] = abs(pa[0])*abs(pa[1])*abs(pa[2])*abs(pa[3])*abs(pa[4]);

                double param[ncontr*(nvar_md+nvar_mb)];
                for (int i=0; i<ncontr; i++){
                        for (int ivar=0; ivar<nvar_md; ivar++){
                                param[i*nvar_md+ivar] = res[i*nvar_md+ivar];
                                if (ivar==1 && i!=2 && i!=4)param[i*nvar_md+ivar] = res[1];
                        }
                        for (int ivar=0; ivar<nvar_mb; ivar++){
                                param[ncontr*nvar_md+i*nvar_mb+ivar] = res[ncontr*nvar_md+i*nvar_mb+ivar];
                        }
                }


        for(int i=0; i<ncontr; i++){
                cout << i << "  " << frac [i] << endl;
        }
                for (int i=0; i<ncontr; i++){
                        auto wrap = [&md_shape, &mb_shape, nevents, param, i, &frac] (double *x, double *par)->double{
                                double bin_widthx = (maxx-minx)/double(nbins);
                                double bin_widthy = (maxy-miny)/double(nbins);
                                double fval_md, fval_mb;
                                if (i==4) fval_md = chebyshev(x, &param[i*nvar_md]);
                                else fval_md = md_shape.func_full(x, &param[i*nvar_md]);
                                fval_mb = mb_shape.func_full(&x[1], &param[ncontr*nvar_md+i*nvar_mb]);
                                return nevents*bin_widthx*bin_widthy*frac[i]*fval_md*fval_mb;
                        };
                        func2D[i] = new TF2(Form("tf2_%d", i), wrap, minx, maxx, miny, maxy, 0);
			
			auto wrap_md = [&md_shape, nevents, param, i, &frac](double *x, double *par)->double{
                                double bin_widthx = (maxx-minx)/double(nbins);
				double fval_md;
                                if (i==4) fval_md = chebyshev(x, &param[i*nvar_md]);
                                else fval_md = md_shape.func_full(x, &param[i*nvar_md]);
				return nevents*bin_widthx*frac[i]*fval_md;
                        };
			funcMD[i] = new TF1(Form("tf1_mD_%d", i), wrap_md, minx, maxx, 0);
                  
			auto wrap_mb = [&mb_shape, nevents, param, i, &frac](double *x, double *par)->double{
                                double bin_widthy = (maxy-miny)/double(nbins);
				double fval_mb;
                                fval_mb = mb_shape.func_full(x, &param[ncontr*nvar_md+i*nvar_mb]);
				return nevents*bin_widthy*frac[i]*fval_mb;
                        };
			funcMB[i] = new TF1(Form("tf1_mb_%d", i), wrap_mb, miny, maxy, 0);
                }
		auto func_sumMD = [funcMD](double *x, double *par)->double{
                        double sum = 0.0;
                        for (int i=0; i<ncontr; i++){
                                sum+= funcMD[i]->Eval(x[0]);
                        }
                        return sum;
                };
                TF1 *tf1_sumMD = new TF1("tf1_sumMD", func_sumMD, minx, maxx, 0);
               auto func_sumMB = [funcMB](double *x, double *par)->double{
                        double sum = 0.0;
                        for (int i=0; i<ncontr; i++){
                                sum+= funcMB[i]->Eval(x[0]);
                        }
                        return sum;
                };
                TF1 *tf1_sumMB = new TF1("tf1_sumMB", func_sumMB, miny, maxy, 0);
               auto func_sum2D = [func2D](double *x, double *par)->double{
                        double sum = 0.0;
                        for (int i=0; i<ncontr; i++){
                                sum+= func2D[i]->Eval(x[0], x[1]);
                        }
                        return sum;
                };
                TF2 *tf2_sum2D = new TF2("tf2_sum2D", func_sum2D, minx, maxx, miny, maxy, 0);

		TCanvas *c2D = new TCanvas("c2D", "", 700, 700);
		       TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
                pad1->Draw();
                pad1->cd();
		hist2D->DrawClone("ep");
		tf2_sum2D->DrawClone("surf same");

        		TH2D histpull2D(*hist2D);
        		histpull2D.SetMinimum();
			TH1D *hpull1D = new TH1D("hpull1D", "", 50, -15., 15.);
        		for (int binx=1; binx<=hist2D->GetNbinsX(); binx++){
                		for (int biny=1; biny<=hist2D->GetNbinsY(); biny++){
                        		double err = hist2D->GetBinError(binx, biny);
                        		if (err==0.0) continue;
                        		double diff  = hist2D->GetBinContent(binx, biny)-tf2_sum2D->Eval(hist2D->GetXaxis()->GetBinCenter(binx), hist2D->GetYaxis()->GetBinCenter(biny));
                       			 histpull2D.SetBinContent(binx, biny, diff/err);
					 hpull1D->Fill(diff/err);
                		}
        		}
			c2D->cd();
			      TPad *pad2 = new TPad("pad2", "", 0.0, 0.0, 1.0, 0.3);
        pad2->Draw();
        pad2->cd();
        histpull2D.SetStats(kFALSE);
        histpull2D.SetFillColor(kBlue);
        histpull2D.GetXaxis()->SetLabelSize(0.15);
        histpull2D.GetXaxis()->SetNdivisions(4);
        histpull2D.GetYaxis()->SetLabelSize(0.15);
        histpull2D.GetYaxis()->SetNdivisions(4);
        histpull2D.GetZaxis()->SetLabelSize(0.15);
        histpull2D.GetZaxis()->SetNdivisions(4);
        histpull2D.DrawClone("surf");

	       double nrot = 100.;
        for (int i=0; i<nrot; ++i){
                c2D->cd();
                pad1->cd();
                pad1->SetPhi(30+i*360./double(nrot));
                pad1->SetTheta(30+i*180./double(nrot));
                pad1->Modified();
                pad1->Update();
                c2D->cd();
                pad2->cd();
                pad2->SetPhi(30+i*360./double(nrot));
                pad2->SetTheta(30+i*180./double(nrot));
                pad2->Modified();
                pad2->Update();
                c2D->Modified();
                c2D->Update();
                c2D->Print("fit_2missPT.gif+50");
        }

        c2D->Print("fit_2missPT.gif++");


		
		TCanvas *cpull1D = new TCanvas("cpull1D", "", 500, 500);
		hpull1D->SetFillColor(kYellow);
		hpull1D->Draw("hist");



 
		TCanvas *c = new TCanvas("c", "", 500, 500);
		Draw_pull(c, histMD, funcMD, tf1_sumMD);
		TCanvas *c_mb = new TCanvas("c_mb", "", 500, 500);
		Draw_pull(c_mb, histMB, funcMB, tf1_sumMB);
   }

void Draw_pull(TCanvas *c, TH1D *hist, TF1 *func[], TF1 *tf1_sum){    		
	c->cd();
        TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
                pad1->SetLogy();
                pad1->Draw();
                pad1->cd();

        hist->SetStats(kFALSE);
        hist->SetMinimum(1.0);
        hist->DrawClone("ep");
        tf1_sum->SetLineColor(kBlack);
        tf1_sum->DrawClone("same");
        int color[ncontr]={kRed, kBlue, kMagenta, kOrange, kCyan, kViolet};
        for (int i=0; i<ncontr; i++){
                func[i]->SetLineColor(color[i]);
                func[i]->DrawClone("same");
        }
	TString name[ncontr] = {"signal", "BuDmunu", "BsDsMunu", "B02DpDsm", "sidebands", "Bu2D0Dsm"};
        TLegend *leg = new TLegend(0.8, 0.5, 1.0, 1.0);
        leg->AddEntry(hist, "data", "p");
        leg->AddEntry(tf1_sum, "sum", "l");
        for(int i=0; i<ncontr; i++) leg->AddEntry(func[i], name[i].Data(), "l");
        leg->DrawClone();

        //hist->Sumw2();
        TH1D histpull1D(*hist);
        histpull1D.SetMinimum();
        for (int bin=1; bin<=hist->GetNbinsX(); bin++){
                double err = hist->GetBinError(bin);
                if (err==0.0) continue;
                double diff  = hist->GetBinContent(bin)-tf1_sum->Eval(hist->GetBinCenter(bin));
                histpull1D.SetBinContent(bin, diff/err);
        }
       c->cd();
        TPad *pad2 = new TPad("pad2", "", 0.0, 0.0, 1.0, 0.3);
        pad2->Draw();
        pad2->cd();
        histpull1D.SetStats(kFALSE);
        histpull1D.SetFillColor(kBlue);
        histpull1D.GetYaxis()->SetLabelSize(0.1);
        histpull1D.GetXaxis()->SetLabelSize(0.1);
        histpull1D.DrawClone("hist");
}
