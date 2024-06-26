#include "D_M_fit_shape.h"

using namespace cpt_b0_analysis;
void D_M_fit_Every(int choice, int sign, bool isUnbinned=true){

	std::string minName = "Minuit2";
 	std::string algoName = "";
	ROOT::Math::Minimizer* min =
	ROOT::Math::Factory::CreateMinimizer(minName, algoName);

	// set tolerance , etc...
   	min->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
   	min->SetMaxIterations(10000);  // for GSL
   	min->SetTolerance(2.5e0);
	if(choice==2) min->SetTolerance(2.5e2);
   	min->SetPrintLevel(2);

   	//min->SetStrategy(2);
   	//min->SetPrecision(0.00001);
   	const int nvar = 7;
	const int ncontr = 6;

   double minx = 1800.;
   double maxx = 1940.;
   double miny = 2700.;
   double maxy = 8300.;
	

   	// create funciton wrapper for minmizer
   	// a IMultiGenFunction type

   	TChain ch("BlindedTree");
	switch (choice){
		case 0:
   		ch.Add("/home/szabelskia/LHCb/MC2016/nomassDmuCut/selected_MCsignal25102023.root");
		break;
		case 1:
		ch.Add("/home/szabelskia/LHCb/MC2016/nomassDmuCut/selected_BuDmunu.root");
		break;
		case 2:
		ch.Add("/home/szabelskia/LHCb/MC2016/nomassDmuCut/selected_MCBsDsMunu25102023.root");
		break;
		case 3:
		ch.Add("/home/szabelskia/LHCb/MC2016/nomassDmuCut/selected_B02DpDsm.root");
		break;
		case 4:
		cout << " No MC for combinatorial " << endl;
		return;
		case 5:
		ch.Add("/home/szabelskia/LHCb/MC2016/nomassDmuCut/selected_MCBu2D0Dsm.root");
		break;
		default:
		cout << " Unknown channel" << endl;
		return;
	}
	double B_M, missPT, D_M, mu_PT, mu_P, mu_eta, K_PT;
	bool charge;
	ch.SetBranchAddress("B_M", &B_M);
   	ch.SetBranchAddress("missPT", &missPT);
   	ch.SetBranchAddress("D_M", &D_M);
   	ch.SetBranchAddress("mu_PT", &mu_PT);
   	ch.SetBranchAddress("mu_P", &mu_P);
   	ch.SetBranchAddress("mu_eta", &mu_eta);
   	ch.SetBranchAddress("K_PT", &K_PT);
   	ch.SetBranchAddress("truecharge", &charge);

	std::vector<double> vect_Dmass;
	int nbins = 40;
	TH1D *hist = new TH1D("hist", "", nbins, minx, maxx);
   	for (int i=0; i<ch.GetEntries(); ++i){
      		ch.GetEntry(i);
     		if (mu_PT<500 || mu_P<5000 || mu_eta<2 || mu_eta>4.5) continue;
      		double B_MMcorr = sqrt(B_M * B_M) +2.0*TMath::Abs(missPT);
      		if (B_MMcorr<miny || B_MMcorr>maxy) continue;
      		if (D_M<minx || D_M>maxx) continue;
      		if (int(charge) == sign){
        	      	vect_Dmass.push_back(D_M);
              		hist->Fill(D_M);
      		}
   	}
	double  nevents = double(vect_Dmass.size());
	md_fit md_shape(minx, maxx);

	auto fchi2 = [&md_shape, vect_Dmass](const double *par)->double{
		double chi2 = 0.0;
		for (auto dmass: vect_Dmass){
			double likelihood = md_shape.func_full(&dmass, par);
			chi2 -= 2.0*log(likelihood);
		}
		return chi2;
	};
	auto fchi2binned = [&md_shape, hist, minx, maxx, nevents, nbins](const double *par)->double{
		double chi2 = 0;
		double bin_width = (maxx-minx)/double(nbins);
		   for (int bin=1; bin<=nbins; bin++){
			   double bincenter = hist->GetBinCenter(bin);
			   double bincont = hist->GetBinContent(bin);
			   double err = hist->GetBinError(bin);
			   double estim = bin_width*double(nevents)*md_shape.func_full(&bincenter, par);
			   if(err!=0.0) chi2+=(bincont-estim)*(bincont-estim)/err/err;
		   }
		   return chi2;
	};
			   
	TString varname[nvar] = {"sigma", "mean", "sigmaCB", "f12", "alpha", "n", "alpha_h"};

        double step = 0.01;
	double initial[ncontr] [nvar] = {
		8.7673, 1870.23, 5.88232, 0.443328, 1.75869, 1.22507, 2.12662,
		8.7673, 1870.23, 5.88232, 0.443328, 1.75869, 1.22507, 2.12662,
		17.7149, 1880.59, 12.893, 2.92653e-07, 0.326914, 4.1634, 1.68152,
		8.7673, 1870.23, 5.88232, 0.443328, 1.75869, 1.22507, 2.12662,
		8.7673, 1870.23, 5.88232, 0.443328, 1.75869, 1.22507, 2.12662,
		8.7673, 1870.23, 5.88232, 0.443328, 1.75869, 1.22507, 2.12662
	};

	ROOT::Math::Functor f_binned(fchi2binned, nvar);
	ROOT::Math::Functor f(fchi2, nvar);
   	if (isUnbinned) min->SetFunction(f);
	else min->SetFunction(f_binned);

	for (int ivar=0; ivar<nvar; ivar++){
		min->SetVariable(ivar, varname[ivar].Data(), initial[choice][ivar], step);
		min->SetVariableLowerLimit(ivar, 0.0);
	}
	
   	double CL_normal = ROOT::Math::normal_cdf(1) -  ROOT::Math::normal_cdf(-1);

	min -> SetErrorDef(ROOT::Math::chisquared_quantile(CL_normal, min->NFree()));
	min->Minimize();
	double err_up, err_down;
	min->GetMinosError(1, err_up, err_down);
	cout << "Minos mean " << err_up << "  " << err_down << endl;

	for (int i=0; i<nvar; i++){
		cout << min->X()[i] << ", ";
	}
	cout << endl;

	TCanvas *c = new TCanvas("c", "", 500, 500);
	   	TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
		pad1->SetLogy();
   	 	pad1->Draw();
   	 	pad1->cd();

	auto funcDraw = [&md_shape, min, nevents, minx, maxx, nbins](double *x, double *par)->double{
		double bin_width = (maxx-minx)/double(nbins);
		return bin_width*double(nevents)*md_shape.func_full(x, par);
	};
	TF1 *tf1 = new TF1("tf1", funcDraw, minx, maxx, nvar);
	tf1->SetParameters(min->X());
	hist->DrawClone("ep");
	tf1->DrawClone("same");
	      
	hist->Sumw2();
    	TH1D histpull1D(*hist);
    	for (int bin=1; bin<=hist->GetNbinsX(); bin++){
        	double err = hist->GetBinError(bin);
        	if (err==0) continue;
        	double diff  = hist->GetBinContent(bin)-tf1->Eval(hist->GetBinCenter(bin));
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
	gSystem->Exec("mkdir -p D_M_figures D_M_results");
	TString name[6] = {"signal", "BuDmunu", "BsDsMunu", "B02DpDsm" , "sidebands", "Bu2D0Dsm"};
	c->SaveAs(Form("D_M_figures/D_M_%s_%d.pdf", name[choice].Data(), sign));
	
	ofstream outfile(Form("D_M_results/res_%s_%d.txt", name[choice].Data(), sign));
   	for (int i =0; i< nvar; i++){
        	outfile << min->X()[i] << "  " << min->Errors()[i] << endl;
   	}
   	outfile.close();

}
      

