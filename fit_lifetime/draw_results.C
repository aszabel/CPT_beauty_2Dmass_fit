#include "pdfs_cpt_model.h"
#include "decay_parameters.h"

using namespace cpt_b0_analysis;
TH1D* make_log_bins(TH1D *hAsym);
void draw_results(){
	TFile *file = new TFile("trees/tree_hists_0.root", "readonly");
	TH1D *hist_B = dynamic_cast<TH1D *>(file->Get("hist_B")); 
	TH1D *hist_Bbar = dynamic_cast<TH1D *>(file->Get("hist_Bbar")); 
	TH1D *hist_rB = dynamic_cast<TH1D *>(file->Get("hist_rB")); 
	TH1D *hist_rBbar = dynamic_cast<TH1D *>(file->Get("hist_rBbar")); 
	TH1D *hist_res = dynamic_cast<TH1D *>(file->Get("hist_res")); 
	TTree *outtree = dynamic_cast<TTree *>(file->Get("outtree"));

	TH1D *logHistB = dynamic_cast<TH1D *>(file->Get("logHistB")); 
	TH1D *logHistBbar = dynamic_cast<TH1D *>(file->Get("logHistBbar")); 
	TH1D *logHist_resB = dynamic_cast<TH1D *>(file->Get("logHist_resB")); 
	TH1D *logHist_resBbar = dynamic_cast<TH1D *>(file->Get("logHist_resBbar")); 
	TH1D *logAsym = dynamic_cast<TH1D *>(file->Get("logAsym")); 
	TH1D *logAsym_res = dynamic_cast<TH1D *>(file->Get("logAsym_res")); 


/*
	//TH1D *logHistB = dynamic_cast<TH1D *>(make_log_bins(hist_B));
	//TH1D *logHistBbar = dynamic_cast<TH1D *>(make_log_bins(hist_Bbar));

	TH1D *hAsym = dynamic_cast<TH1D *>(logHistB->GetAsymmetry(logHistBbar));
	hAsym->SetName("hAsym");
	//hAsym->GetAsymmetry(hist_Bbar);		
	double output[10];
	outtree->SetBranchAddress("rez", &output[0]);
        outtree->SetBranchAddress("imz", &output[1]);
        outtree->SetBranchAddress("Gamma", &output[2]);
        outtree->SetBranchAddress("eta", &output[3]);
        outtree->SetBranchAddress("t_shift", &output[4]);
        outtree->SetBranchAddress("alpha", &output[5]);
        outtree->SetBranchAddress("beta", &output[6]);
        outtree->SetBranchAddress("t_shift_anit", &output[7]);
        outtree->SetBranchAddress("alpha_anti", &output[8]);
        outtree->SetBranchAddress("beta_anti", &output[9]);

	outtree->GetEntry(0);

	int Nrbins = hist_rB->GetNbinsX();
	int Nbins = hist_B->GetNbinsX();
        std::vector<double> r_bin_center, r_bin_cont, r_bin_center_bar, r_bin_cont_bar;
        for (int r_bin=1; r_bin<=Nrbins; ++r_bin){
                r_bin_center.push_back(hist_rB->GetBinCenter(r_bin));
                r_bin_cont.push_back(hist_rB->GetBinContent(r_bin));
                r_bin_center_bar.push_back(hist_rBbar->GetBinCenter(r_bin));
                r_bin_cont_bar.push_back(hist_rBbar->GetBinContent(r_bin));
        }

	pdfs PDFs;
        auto func_conv = [&PDFs, r_bin_center, r_bin_cont, Nrbins] (double x, const double *par, bool bbar) -> double{
                double value =0.0;
                double Tau = x;
                for (int r_bin=1; r_bin<=Nrbins; ++r_bin){
                        double r_rec = r_bin_center[r_bin-1];
                        double n_r = r_bin_cont[r_bin-1];

                        double time_c =Tau-r_rec; //time of mean exp value in bin
                        if (time_c<Params::tMin_gen || time_c> Params::tMax_gen ) continue;
                        if (bbar)
                                value += PDFs.time_acceptance(&time_c, &par[4])*PDFs.func1D_f(&time_c, par)*n_r;
                        else
                                value += PDFs.time_acceptance(&time_c, &par[7])*PDFs.func1D_fbar(&time_c, par)*n_r;
                }
        return value;
        };

				TH1D *hist_resB = new TH1D("hist_resB", "", Nbins, Params::tMin, Params::tMax);
				double norm = 0.0;
                                for (int bin=1; bin<=Nbins; bin++){
                                        double t = hist_B->GetBinCenter(bin);
                                        double expect = func_conv(t, output, true);
					hist_resB->SetBinContent(bin, expect);
					hist_resB->SetBinError(bin, 0.0);
                                        norm+=expect;
                                }
				hist_resB->Scale(hist_B->Integral()/norm);

				cout << output[2] << " GGGGGGGAaaaaammmmaaaaaa!\n\n";
				
				TH1D *hist_resBbar = new TH1D("hist_resBbar", "", Nbins, Params::tMin, Params::tMax);
				double normbar = 0.0;
                                for (int bin=1; bin<=Nbins; bin++){
                                        double t = hist_Bbar->GetBinCenter(bin);
                                        double expect = func_conv(t, output, false);
					hist_resBbar->SetBinContent(bin, expect);
					hist_resBbar->SetBinError(bin, 0.0);
                                        normbar+=expect;
                                }
				hist_resBbar->Scale(hist_Bbar->Integral()/norm);
  	auto func_norm = [func_conv, norm, hist_B] (double *x, double *par)->double{
               return func_conv(x[0], par, true)/norm*hist_B->Integral();
        };
        auto func_norm_bar = [func_conv, normbar, hist_Bbar] (double *x, double *par)->double{
                return func_conv(x[0], par, false)/normbar*hist_Bbar->Integral();
        };

	TF1 *funcB = new TF1("funcB", func_norm, Params::tMin, Params::tMax, 10);
	funcB->SetParameters(output);



	TH1D *logHist_resB = dynamic_cast<TH1D *>(make_log_bins(hist_resB));
*/
	//TH1D *pull_B = new TH1D("pull_B", "", Nbins, Params::tMin, Params::tMax);
	TH1D *pull_B =dynamic_cast<TH1D *> (logHist_resB->Clone("pullB")) ;
	pull_B->Reset();
	for (int bin=1; bin<=pull_B->GetNbinsX(); bin++){
		//double time = hist_B->GetBinCenter(bin);
		//double expect = funcB->Eval(time);
		//logHist_resB->SetBinContent(bin, logHist_resB->GetBinWidth(bin)*logHist_resB->GetBinContent(bin));
		double expect = logHist_resB->GetBinContent(bin);
		double cont  = logHistB->GetBinContent(bin);
		double err  = logHistB->GetBinError(bin);
		if (err!=0)
			pull_B->SetBinContent(bin, (expect-cont)/err);
	}

//	TF1 *funcBbar = new TF1 ("funcBbar", func_norm_bar, Params::tMin, Params::tMax, 10);
//	funcBbar->SetParameters(output);
//	TH1D *logHist_resBbar = dynamic_cast<TH1D *>(make_log_bins(hist_resBbar));
	
	TH1D *pull_Bbar =dynamic_cast<TH1D *> (logHist_resBbar->Clone("pullBbar")) ;
	pull_Bbar->Reset();
	for (int bin=1; bin<=pull_Bbar->GetNbinsX(); bin++){
		//double time = hist_Bbar->GetBinCenter(bin);
		//double expect = funcBbar->Eval(time);
		double expect = logHist_resBbar->GetBinContent(bin);
		double cont  = logHistBbar->GetBinContent(bin);
		double err  = logHistBbar->GetBinError(bin);
		if (err!=0)
			pull_Bbar->SetBinContent(bin, (expect-cont)/err);
	}


	logHistB->SetStats(kFALSE);
	pull_B->SetStats(kFALSE);
	logHistBbar->SetStats(kFALSE);
	pull_Bbar->SetStats(kFALSE);
	logAsym->SetStats(kFALSE);

	TCanvas *c = new TCanvas("c", "", 1000, 500);
	c->Divide(2,1);
	c->cd(1);
	gPad->Draw();
	       TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
                pad1->SetLogx();
                pad1->SetLogy();
                pad1->Draw();
                pad1->cd();

		logHistB->GetXaxis()->SetMoreLogLabels();
		logHistB->Draw();
		logHist_resB->SetLineWidth(1);
		logHist_resB->SetFillColor(0);
		logHist_resB->SetLineColor(kRed);
		logHist_resB->Draw("same l");
	c->cd(1);
	TPad *pad2 = new TPad("pad2", "", 0.0, 0.0, 1.0, 0.3);
        	pad2->Draw();
        	pad2->cd();
		pad2->SetLogx();
		pull_B->GetXaxis()->SetMoreLogLabels();
		pull_B->SetFillColor(kBlue);
		pull_B->Draw("hist");
 
	//gPad->Update();
	c->cd(2);
	gPad->Draw();
	//cout << hist_Bbar->Integral() << endl;
	TPad *pad1b = new TPad("pad1b", "", 0.0, 0.3, 1.0, 1.0);
                pad1b->SetLogx();
                pad1b->SetLogy();

                pad1b->Draw();
                pad1b->cd();
		
		logHistBbar->GetXaxis()->SetMoreLogLabels();
		logHistBbar->Draw();
		//funcBbar->Draw("same");
		logHist_resBbar->SetLineColor(kRed);
		logHist_resBbar->Draw("same l");

		c->cd(2);
	       TPad *pad2b = new TPad("pad2b", "", 0.0, 0.0, 1.0, 0.3);
        	pad2b->Draw();
        	pad2b->cd();
		pad2b->SetLogx();
		pull_Bbar->GetXaxis()->SetMoreLogLabels();
		pull_Bbar->SetFillColor(kBlue);
		pull_Bbar->Draw("hist");
 
	c->Update();
    


    // Cleanup

	TCanvas *c1 =  new TCanvas("c1", "", 500, 500); 
		gPad->SetLogx();
		logAsym->GetXaxis()->SetMoreLogLabels();
		logAsym->Draw();
		logAsym->SetMinimum(-0.15);
		logAsym->SetMaximum(0.15);
		logAsym_res->SetLineColor(kRed);
		logAsym_res->Draw("same");
		c1->Update();

	//file->Close();

}

TH1D* make_log_bins(TH1D *hAsym){
    // Define the number of logarithmic bins and the range
    int nLogBins = 50;      // Number of logarithmic bins
    double xMin = hAsym->GetXaxis()->GetXmin();  // Get the original xMin
    double xMax = hAsym->GetXaxis()->GetXmax();  // Get the original xMax

    // Create an array for the logarithmic bin edges
    double *logBins = new double[nLogBins + 1];

    // Logarithmic binning: Calculate the bin edges
    for (int i = 0; i <= nLogBins; ++i) {
        logBins[i] = TMath::Power(10, (i * (TMath::Log10(xMax) - TMath::Log10(xMin)) / nLogBins) + TMath::Log10(xMin));
    }

    // Create the new histogram with logarithmic bins
    TH1D *logHist = new TH1D("logHist", "Rebinned Histogram with Logarithmic Bins", nLogBins, logBins);

    // Loop through the bins of the original histogram
    for (int i = 1; i <= hAsym->GetNbinsX(); ++i) {
        double binCenter = hAsym->GetBinCenter(i);
        double binContent = hAsym->GetBinContent(i);
        double binError = hAsym->GetBinError(i);

        // Find the corresponding bin in the new logarithmic histogram
        int logBin = logHist->FindBin(binCenter);

        // Add the content and error from the original histogram to the new logarithmic histogram
        logHist->AddBinContent(logBin, binContent);
        logHist->SetBinError(logBin, sqrt(pow(logHist->GetBinError(logBin), 2) + pow(binError, 2)));
    }
    return logHist;
}


