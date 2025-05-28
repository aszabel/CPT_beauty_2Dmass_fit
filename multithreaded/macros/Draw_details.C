#include "D_M_fit_shape.h"
#include "M_B_2missPT_fit.h"
#include "ChebyshevPDF.h"
#include <cmath>


using namespace cpt_b0_analysis;
   const int nbins = 40;

void Draw_pull(TCanvas *c, TH1D *hist, TF1 *func[], TF1 *tf1_sum);

double LikelihoodRatioTest(TH1* h_data, TH1* h_mc);


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

	TH1D *h_tau_signal_sweighted = new TH1D("h_tau_sw", "time signal MC vs sweighted", nbins, Config::tMin, Config::tMax);
	TH1D *h_MC_signal = new TH1D("h_tau_mc", "time signal MC vs sweighted", nbins, Config::tMin, Config::tMax);

        std::string fileWeightsName = "Tree_sWeights";
        std::string TreeName = "Tree_sWeights";
        std::vector<double> *vec_Tau_ptr = nullptr;
        std::vector<double> *vec_sWeights_ptr = nullptr;

	std::string sign_name = "";
	if (Config::sign == 1)
		sign_name = "plus";
	else
		sign_name = "minus";
        TFile weightFile_plus(Form("../toy_res/toy_%d/%s_%s.root", Config::randSeed, fileWeightsName.c_str(), sign_name.c_str()), "readonly");
        TTree *tree_sWeights = (TTree*)weightFile_plus.Get(TreeName.c_str());
        tree_sWeights -> SetBranchAddress("vec_Tau", &vec_Tau_ptr);
        tree_sWeights -> SetBranchAddress("vec_sWeights", &vec_sWeights_ptr);

        tree_sWeights -> GetEntry(0);

        auto *vec_Tau = vec_Tau_ptr->data();
        auto *vec_sWeights = vec_sWeights_ptr->data();



        for (int e=0; e<int(vec_Tau_ptr->size()); e++){
                h_tau_signal_sweighted->Fill(vec_Tau[e], vec_sWeights[e]);
        }
   TChain chMC(Config::chainName.c_str());
   chMC.Add("/mnt/home/share/lhcb/CPT_beauty/MC2016/selected/selectedMagDown_B2Dmunu_signal_taustrip_nomassDmuCut/selected_MCsignal25102023.root");

   double tau;

   chMC.SetBranchAddress("B_M", &B_M);
   chMC.SetBranchAddress("missPT", &missPT);
   chMC.SetBranchAddress("D_M", &D_M);
   chMC.SetBranchAddress("mu_PT", &mu_PT);
   chMC.SetBranchAddress("mu_P", &mu_P);
   chMC.SetBranchAddress("mu_eta", &mu_eta);
   chMC.SetBranchAddress("K_PT", &K_PT);
   chMC.SetBranchAddress("truecharge", &charge);
   chMC.SetBranchAddress("Tau", &tau);

   

   for (int i=0; i<chMC.GetEntries(); ++i){
      chMC.GetEntry(i);
      if (mu_PT<Config::muPTmin || mu_P<Config::muPmin || mu_eta< Config::eta_min || mu_eta>Config::eta_max) continue;
      double B_MMcorr = B_M +2.0*missPT;
      if (B_MMcorr<Config::minBMcorr || B_MMcorr> Config::maxBMcorr) continue;
      if (D_M<Config::minDM || D_M>Config::maxDM) continue;
      if (int(charge) == Config::sign){
              h_MC_signal->Fill(tau);
      }
   }
        h_MC_signal->Scale(h_tau_signal_sweighted->Integral()/h_MC_signal->Integral());

		TCanvas *c_mc = new TCanvas("c_mc", "", 500, 500);
	        TPad *mc_pad1 = new TPad("mc_pad1", "", 0.0, 0.3, 1.0, 1.0);
        	        mc_pad1->SetLogy();
                	mc_pad1->Draw();
               		mc_pad1->cd();

        		h_tau_signal_sweighted->SetStats(kFALSE);
        		h_tau_signal_sweighted->SetMinimum(1.0);
        		h_tau_signal_sweighted->DrawClone("ep");
			h_MC_signal->DrawClone("same hist");
	

			//h_tau_signal_sweighted->Scale(1./h_tau_signal_sweighted->Integral());
			//h_MC_signal->Scale(1./h_MC_signal->Integral());
			//double p_value = h_tau_signal_sweighted->KolmogorovTest(h_MC_signal, "X");  
			double p_value = h_MC_signal->KolmogorovTest(h_tau_signal_sweighted, "X");  

			cout << p_value << " <=====  p-value K-S \n";


			p_value =  h_tau_signal_sweighted->Chi2Test(h_MC_signal, "WW");
			cout << p_value << " <=====  p-value chi2 test \n";


                        p_value = LikelihoodRatioTest(h_tau_signal_sweighted, h_MC_signal);
	
		        TH1D histpull1D(*h_tau_signal_sweighted);
        		histpull1D.SetMinimum();
       			 for (int bin=1; bin<=h_tau_signal_sweighted->GetNbinsX(); bin++){
                		double err = h_tau_signal_sweighted->GetBinError(bin);
				double err2 = h_MC_signal->GetBinError(bin);
				err = TMath::Sqrt(err*err+err2*err2);
                		if (err==0.0) continue;
                		double diff  = h_tau_signal_sweighted->GetBinContent(bin)-h_MC_signal->GetBinContent(bin);
                		histpull1D.SetBinContent(bin, diff/err);
       			 }
			
        	c_mc->cd();
      	  	TPad *mc_pad2 = new TPad("mc_pad2", "", 0.0, 0.0, 1.0, 0.3);
        	mc_pad2->Draw();
       		 mc_pad2->cd();
       		 histpull1D.SetStats(kFALSE);
        	histpull1D.SetFillColor(kBlue);
        	histpull1D.GetYaxis()->SetLabelSize(0.1);
        	histpull1D.GetXaxis()->SetLabelSize(0.1);
        	histpull1D.DrawClone("hist");	

        	c_mc->SaveAs("../results/swtau_signalMC.pdf");
        	c_mc->SaveAs("../results/swtau_signalMC.C");



   	double res[(Config::nvar_md+Config::nvar_mb)*Config::ncontr+Config::ncontr];
        TString resname = Form("../toy_res/toy_%d/best_results_unbinned/fit2D_best/results_%d_3.txt",Config::randSeed, Config::sign);
        cout << resname << endl;
	ifstream input(resname.Data());
	double x, dx;
	int i=0;
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

        		TH2D *histpull2D = (TH2D*)hist2D->Clone("histpull2D");
        		histpull2D->SetMinimum();
			TH1D *hpull1D = new TH1D("hpull1D", "", 50, -15., 15.);
        		for (int binx=1; binx<=hist2D->GetNbinsX(); binx++){
                		for (int biny=1; biny<=hist2D->GetNbinsY(); biny++){
                        		double err = hist2D->GetBinError(binx, biny);
                        		if (err==0.0) continue;
                        		double diff  = hist2D->GetBinContent(binx, biny)-tf2_sum2D->Eval(hist2D->GetXaxis()->GetBinCenter(binx), hist2D->GetYaxis()->GetBinCenter(biny));
                       			 histpull2D->SetBinContent(binx, biny, diff/err);
					 hpull1D->Fill(diff/err);
                		}
        		}
			c2D->cd();
			      TPad *pad2 = new TPad("pad2", "", 0.0, 0.0, 1.0, 0.3);
        pad2->Draw();
        pad2->cd();
        histpull2D->SetStats(kFALSE);
        histpull2D->SetFillColor(kBlue);
        histpull2D->GetXaxis()->SetLabelSize(0.15);
        histpull2D->GetXaxis()->SetNdivisions(4);
        histpull2D->GetYaxis()->SetLabelSize(0.15);
        histpull2D->GetYaxis()->SetNdivisions(4);
        histpull2D->GetZaxis()->SetLabelSize(0.15);
        histpull2D->GetZaxis()->SetNdivisions(4);
        histpull2D->DrawClone("surf");
        
	c2D->cd();
	c2D->Update();
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
		cpull1D->SaveAs("../results/pull1D.C");

		

 
		TCanvas *c = new TCanvas("c", "", 500, 500);
		Draw_pull(c, histMD, funcMD, tf1_sumMD);
		
        	c->SaveAs("../results/fitMD.pdf");
        	c->SaveAs("../results/fitMD.C");
		
		TCanvas *c_mb = new TCanvas("c_mb", "", 500, 500);
		Draw_pull(c_mb, histMB, funcMB, tf1_sumMB);
        	c_mb->SaveAs("../results/fitMB.pdf");
        	c_mb->SaveAs("../results/fitMB.C");

		TFile *histout = new TFile ("../results/output_hists.root", "recreate");
		hpull1D->Write();
		histMD->Write();
		histMB->Write();
		for (int i = 0; i<Config::ncontr; ++i){
			funcMD[i]->Write();
			funcMB[i]->Write();
		}
		histpull2D->Write();
		hist2D->Write();
                tf2_sum2D->Write();		
                tf1_sumMD->Write();		
                tf1_sumMB->Write();		
		h_tau_signal_sweighted->Write();
		h_MC_signal->Write();
		histout->Close();
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

double LikelihoodRatioTest(TH1* h_data, TH1* h_mc) {
    if (h_data->GetNbinsX() != h_mc->GetNbinsX()) {
        std::cerr << "Histograms must have the same binning!" << std::endl;
        return -1;
    }

    double scale = h_data->Integral() / h_mc->Integral();
    TH1* h_mc_scaled = (TH1*) h_mc->Clone("h_mc_scaled");
    h_mc_scaled->Scale(scale);

    double logL_data = 0.0;
    double logL_mc = 0.0;

    for (int i = 1; i <= h_data->GetNbinsX(); ++i) {
        double n_obs = h_data->GetBinContent(i);
        double n_exp = h_mc_scaled->GetBinContent(i);

        if (n_exp <= 0) continue;

        if (n_obs > 0)
            logL_mc += n_obs * std::log(n_exp) - n_exp - std::lgamma(n_obs + 1);
        else
            logL_mc += -n_exp;

        logL_data += n_obs > 0 ? n_obs * std::log(n_obs) - n_obs - std::lgamma(n_obs + 1)
                               : -n_obs;
    }

    double LR = -2.0 * (logL_mc - logL_data);

    std::cout << "Likelihood Ratio: -2ln(L_MC / L_Data) = " << LR << std::endl;

    delete h_mc_scaled;
    return LR;
}

