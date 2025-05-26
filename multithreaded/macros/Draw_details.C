#include "D_M_fit_shape.h"
#include "M_B_2missPT_fit.h"
#include "ChebyshevPDF.h"


using namespace cpt_b0_analysis;
   const int nbins = 40;

void Draw_pull(TCanvas *c, TH1D *hist, TF1 *func[], TF1 *tf1_sum);



using namespace cpt_b0_analysis;
void Draw_details(std::string config_file){

        // Load config
        std::cout<<config_file<<endl;
        if (Config::load(config_file)){
                std::cerr<< " Bad config file! " << std::endl;
                return;
        }

	TH2D *hist2D = new TH2D("h2D", "", nbins, Config::minDM, Config::maxDM, nbins, Config::minBMcorr, Config::maxBMcorr);
	TH1D *histMD = new TH1D("hMD", "", nbins, Config::minDM, Config::maxDM);
	TH1D *histMB = new TH1D("hMB", "", nbins, Config::minBMcorr, Config::maxBMcorr);
   TChain ch(Config::chainName.c_str());
   ch.Add(Config::input_files[0].c_str());
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

   
   int Nentries = Config::nentries;
   if (Nentries == -1)
	Nentries = ch.GetEntries();

   for (int i=0; i<Nentries; ++i){
      ch.GetEntry(i);
      if (mu_PT<Config::muPTmin || mu_P<Config::muPmin || mu_eta< Config::eta_min || mu_eta>Config::eta_max) continue;
      double B_MMcorr = B_M +2.0*missPT;
      if (B_MMcorr<Config::minBMcorr || B_MMcorr> Config::maxBMcorr) continue;
      if (D_M<Config::minDM || D_M>Config::maxDM) continue;
      if (int(charge) == Config::sign){
              hist2D->Fill(D_M, B_MMcorr);
	      histMD->Fill(D_M);
	      histMB->Fill(B_MMcorr);
      }
   }
        double  nevents = double(hist2D->Integral());

   	double res[(Config::nvar_md+Config::nvar_mb)*Config::ncontr+Config::ncontr];
        TString resname = Form("../toy_res/toy_%d/best_results_unbinned/fit2D_best/results_%d_3.txt",Config::randSeed, Config::sign);
        cout << resname << endl;
	ifstream input(resname.Data());
	double x, dx;
	int i=0;
	cout << " Hereeeeeeeeeeeee\n";
	input>>x>>dx;
	cout << x << endl; 
	while (input>>x>>dx){
		res[i] = x;
		cout << res[i] << endl;
	 	i++;
	}
	input.close();


        TF2 *func2D[Config::ncontr];
	TF1 *funcMD[Config::ncontr], *funcMB[Config::ncontr];
        double frac[Config::ncontr];
        double *pa = &res[Config::ncontr*(Config::nvar_md+Config::nvar_mb)];

		double sum_frac=0.0;
               for (int i=0; i < Config::ncontr-1; i++){
                        frac[i] = abs(pa[i]);
                        sum_frac+=frac[i];
                }
                frac[Config::ncontr-1] = abs(1.0-sum_frac);

		const int nall = Config::ncontr*(Config::nvar_md+Config::nvar_mb);

                double param[nall];
                for (int i=0; i<Config::ncontr; i++){
                        for (int ivar=0; ivar<Config::nvar_md; ivar++){
                                param[i*Config::nvar_md+ivar] = res[i*Config::nvar_md+ivar];
                                if (ivar==1 && i!=2 && i!=4)param[i*Config::nvar_md+ivar] = res[1];
                        }
                        for (int ivar=0; ivar<Config::nvar_mb; ivar++){
                                param[Config::ncontr*Config::nvar_md+i*Config::nvar_mb+ivar] = res[Config::ncontr*Config::nvar_md+i*Config::nvar_mb+ivar];
                        }
                }


        for(int i=0; i<Config::ncontr; i++){
                cout << i << "  " << frac [i] << endl;
        }


        auto D_PDFs = Config::getVectorPDFs("Dmass");
	std::vector<cpt_b0_analysis::PDFInterface *> D_PDFs_get = {};
        if (int(D_PDFs.size()) != Config::ncontr){
                std::cout<< " NO D_PDFs \n";
                return 1;
        }
        for (int i=0; i<Config::ncontr; ++i){
                //auto pdf = D_PDFs[i].get();
		D_PDFs_get.push_back( D_PDFs[i].get());
                if (!D_PDFs_get[i]){
                        std::cout<< "Nullptr passed as pdf\n";
                        return 1;
                }
        }

        auto B_PDFs = Config::getVectorPDFs("Bmass");
	std::vector<cpt_b0_analysis::PDFInterface *> B_PDFs_get={};
        if (int(B_PDFs.size()) != Config::ncontr){
                std::cout<< " NO B_PDFs \n";
                return 1;
        }
        for (int i=0; i<Config::ncontr; ++i){
                //auto pdf = B_PDFs[i].get();
		B_PDFs_get.push_back(B_PDFs[i].get());
                if (!B_PDFs_get[i]){
                        std::cout<< "Nullptr passed as pdf\n";
                        return 1;
                }
        }


        for (int i = 0; i < Config::ncontr; i++)
                {
                        D_PDFs_get[i]->CalcIntegral(&param[i * Config::nvar_md], Config::minDM, Config::maxDM);
                        B_PDFs_get[i]->CalcIntegral(&param[Config::ncontr * Config::nvar_md + i * Config::nvar_mb], Config::minBMcorr, Config::maxBMcorr);
                }

                double bin_widthx = (Config::maxDM-Config::minDM)/double(nbins);
                double bin_widthy = (Config::maxBMcorr-Config::minBMcorr)/double(nbins);
                for (int i=0; i<Config::ncontr; i++){
                        auto wrap = [&D_PDFs_get, &B_PDFs_get, nevents, &param, i, &frac, bin_widthx, bin_widthy] (double *x, double *par)->double{
                                double fval_md, fval_mb;
                                fval_md = D_PDFs_get[i]->EvalPDF(x, &param[i*Config::nvar_md]);
                                fval_mb = B_PDFs_get[i]->EvalPDF(&x[1], &param[Config::ncontr*Config::nvar_md+i*Config::nvar_mb]);
                                return nevents*bin_widthx*bin_widthy*frac[i]*fval_md*fval_mb;
                        };
                        func2D[i] = new TF2(Form("tf2_%d", i), wrap, Config::minDM, Config::maxDM, Config::minBMcorr, Config::maxBMcorr, 0);
			
			auto wrap_md = [&D_PDFs_get, nevents, &param, i, &frac, bin_widthx](double *x, double *par)->double{
				double fval_md = D_PDFs_get[i]->EvalPDF(x, &param[i*Config::nvar_md]);
				return nevents*bin_widthx*frac[i]*fval_md;
                        };
			funcMD[i] = new TF1(Form("tf1_mD_%d", i), wrap_md, Config::minDM, Config::maxDM, 0);
                  
			auto wrap_mb = [&B_PDFs_get, nevents, &param, i, &frac, bin_widthy](double *x, double *par)->double{
				double fval_mb = B_PDFs_get[i]->EvalPDF(&x[0], &param[Config::ncontr*Config::nvar_md+i*Config::nvar_mb]);;
				return nevents*bin_widthy*frac[i]*fval_mb;
                        };
			funcMB[i] = new TF1(Form("tf1_mb_%d", i), wrap_mb, Config::minBMcorr, Config::maxBMcorr, 0);
                }
		auto func_sumMD = [&funcMD](double *x, double *par)->double{
                        double sum = 0.0;
                        for (int i=0; i<Config::ncontr; i++){
                                sum+= funcMD[i]->Eval(x[0]);
                        }
                        return sum;
                };
                TF1 *tf1_sumMD = new TF1("tf1_sumMD", func_sumMD, Config::minDM, Config::maxDM, 0);
               auto func_sumMB = [&funcMB](double *x, double *par)->double{
                        double sum = 0.0;
                        for (int i=0; i<Config::ncontr; i++){
                                sum+= funcMB[i]->Eval(x[0]);
                        }
                        return sum;
                };
                TF1 *tf1_sumMB = new TF1("tf1_sumMB", func_sumMB, Config::minBMcorr, Config::maxBMcorr, 0);
               auto func_sum2D = [&func2D](double *x, double *par)->double{
                        double sum = 0.0;
                        for (int i=0; i<Config::ncontr; i++){
                                sum+= func2D[i]->Eval(x[0], x[1]);
                        }
                        return sum;
                };
                TF2 *tf2_sum2D = new TF2("tf2_sum2D", func_sum2D, Config::minDM, Config::maxDM, Config::minBMcorr, Config::maxBMcorr, 0);

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
        
        c2D->SaveAs("../results/fit2D.pdf");
        c2D->SaveAs("../results/fit2D.C");
/*
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
*/

		
		TCanvas *cpull1D = new TCanvas("cpull1D", "", 500, 500);
		hpull1D->SetFillColor(kYellow);
		hpull1D->Draw("hist");
        
		cpull1D->SaveAs("../results/pull1D.pdf");

		

 
		TCanvas *c = new TCanvas("c", "", 500, 500);
		Draw_pull(c, histMD, funcMD, tf1_sumMD);
		
        	c->SaveAs("../results/fitMD.pdf");
        	c->SaveAs("../results/fitMD.C");
		
		TCanvas *c_mb = new TCanvas("c_mb", "", 500, 500);
		Draw_pull(c_mb, histMB, funcMB, tf1_sumMB);
        	c_mb->SaveAs("../results/fitMB.pdf");
        	c_mb->SaveAs("../results/fitMB.C");
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
        int color[]={kRed, kBlue, kMagenta, kOrange, kCyan, kViolet};
        for (int i=0; i<Config::ncontr; i++){
                func[i]->SetLineColor(color[i]);
                func[i]->DrawClone("same");
        }
	auto name = Config::contrName;
        TLegend *leg = new TLegend(0.8, 0.5, 1.0, 1.0);
        leg->AddEntry(hist, "data", "p");
        leg->AddEntry(tf1_sum, "sum", "l");
        for(int i=0; i<Config::ncontr; i++) leg->AddEntry(func[i], name[i].c_str(), "l");
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
