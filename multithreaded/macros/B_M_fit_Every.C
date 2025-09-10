#include "M_B_2missPT_fit.h"
#include "config.h"

using namespace cpt_b0_analysis;
void B_M_fit_Every(std::string config_file){

	std::string minName = "Minuit2";
 	std::string algoName = "";
	

	// Load config
	std::cout<<config_file<<endl;
	if (Config::load(config_file)){
		std::cerr<< " Bad config file! " << std::endl;
		return;
	}


	TString binning = "";
	if (Config::binned)
		binning = "binned";
	else 	
		binning = "unbinned";

	TString path_results = Form("%s_%s", Config::MC_directory_MB.c_str(), binning.Data());

	gSystem->Exec(Form("mkdir -p %s %s_figures", path_results.Data(), path_results.Data()));
	int fit_id = 0;
	for (auto& choice: Config::int_choose_fits)
	{
	std::cout<<"##### Running fit: "<<Config::Fits[fit_id]<<endl;

	// set tolerance , etc...
	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
   	min->SetMaxFunctionCalls(Config::functionCalls); // for Minuit/Minuit2
   	min->SetMaxIterations(10000);  // for GSL
   	min->SetTolerance(Config::tolerance[fit_id]);
   	min->SetPrintLevel(Config::printLevel);
	TRandom rand;
	if (Config::randSeed >-1)
		rand.SetSeed(Config::randSeed);

   	//min->SetStrategy(2);
   	//min->SetPrecision(0.00001);
  	const int nvar = Config::nvar_mb;
	const int ncontr = Config::ncontr;

   	// create funciton wrapper for minmizer
   	// a IMultiGenFunction type

	std::cout<<"Load data: "<<Config::input_files[fit_id]<<std::endl;
   	TChain ch(Config::chainName.c_str());
	ch.Add(Config::input_files[fit_id].c_str());
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

	std::vector<double> vect_Bmass;
	int nbins = 100;
	TH1D *hist = new TH1D("hist", "", nbins, Config::minBMcorr, Config::maxBMcorr);

	int nentries = Config::nentries;
	if (nentries < 0) nentries = ch.GetEntries(); 
	if (nentries > ch.GetEntries()){
		std::cerr << "The value of 'nentries' exceeds the number of events in the file." << std::endl;
                return;
        }
	std::cout<<"Events: "<<nentries<<std::endl;
	double B_MMcorr;
   	for (int i=0; i<nentries; ++i){
      		ch.GetEntry(i);
			if (mu_PT < Config::muPTmin || mu_P < Config::muPmin || mu_eta < Config::eta_min || mu_eta > Config::eta_max)
				continue;	
      		B_MMcorr = sqrt(B_M * B_M) +2.0*TMath::Abs(missPT);
			if (B_MMcorr < Config::minBMcorr || B_MMcorr > Config::maxBMcorr)
				continue;
			if (D_M < Config::minDM || D_M > Config::maxDM)
				continue;
      		if (int(charge) == Config::sign){
        	      	vect_Bmass.push_back(B_M+2.0*missPT);
              		hist->Fill(B_M+2.0*missPT);
      		}
   	}
	double  nevents = double(vect_Bmass.size());
	const auto& B_PDFs = Config::getVectorPDFs("Bmass");


	auto fchi2 = [&B_PDFs, vect_Bmass, choice](const double *par)->double{
		double chi2 = 0.0;
		B_PDFs[choice]->CalcIntegral(par, Config::minBMcorr, Config::maxBMcorr);
		for (auto bmass: vect_Bmass){
			double likelihood = B_PDFs[choice]->EvalPDF(&bmass, par);
			if (likelihood>=1.0|| likelihood <=0.0)
			{
				continue;
			}
			chi2 -= 2.0*log(likelihood);
		}
		return chi2;
	};
	auto fchi2binned = [&B_PDFs, hist, nevents, nbins, choice](const double *par)->double{
		double chi2 = 0;
		B_PDFs[choice]->CalcIntegral(par, Config::minBMcorr, Config::maxBMcorr);
		double bin_width = (Config::maxBMcorr-Config::minBMcorr)/double(nbins);
		   for (int bin=1; bin<=nbins; bin++){
			   double bincenter = hist->GetBinCenter(bin);
			   double bincont = hist->GetBinContent(bin);
			   double err = hist->GetBinError(bin);
			   double estim = bin_width*double(nevents)*B_PDFs[choice]->EvalPDF(&bincenter, par);
			   if(err!=0.0) chi2+=(bincont-estim)*(bincont-estim)/err/err;
		   }
		   return chi2/double(nbins);
	};
			   
	double step = 0.1;

	ROOT::Math::Functor f_binned(fchi2binned, nvar);
	ROOT::Math::Functor f(fchi2, nvar);
   	if (Config::binned) min->SetFunction(f_binned);
	else min->SetFunction(f);

	for (int ivar=0; ivar<nvar; ivar++){
		/*if (ivar== 2 || ivar == 3 || ivar == 5 || ivar == 7) 
			step = 0.01; 
		else
			step = 10.0;*/
		step = abs(0.01*Config::init_values[choice][ivar])+0.01;
		min->SetVariable(ivar, Config::varname_mb[ivar].c_str(), Config::init_values[choice][ivar], step);
		//min -> FixVariable(ivar);		
// TODO from config
	}
	min->SetVariableLimits(min->VariableIndex("f12"), -1.0, 1.0);
	min->SetVariableLimits(min->VariableIndex("f12_gaus"), -1.0, 1.0);
	//list of fixed variables form config
	for (const auto& fix: Config::fixVect){
		if (min->VariableIndex(fix) >= 0) {
			min->FixVariable(min->VariableIndex(fix));
			std::cout << fix << "  " << min->VariableIndex(fix) << std::endl;
		}
	}
	
   	double CL_normal = ROOT::Math::normal_cdf(1) -  ROOT::Math::normal_cdf(-1);

	min -> SetErrorDef(ROOT::Math::chisquared_quantile(CL_normal, min->NFree()));
	min->Minimize();
	//double err_up, err_down;
	//min->GetMinosError(1, err_up, err_down);
	//cout << "Minos mean " << err_up << "  " << err_down << endl;

	for (int i=0; i<nvar; i++){
		cout << min->X()[i] << ", ";
	}
	cout << endl;

	TCanvas *c = new TCanvas("c", "", 500, 500);
	   	TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
		pad1->SetLogy();
   	 	pad1->Draw();
   	 	pad1->cd();

	auto funcDraw = [&B_PDFs, min, nevents, nbins, choice](double *x, double *par)->double{
		B_PDFs[choice]->CalcIntegral(par, Config::minBMcorr, Config::maxBMcorr);
		double bin_width = (Config::maxBMcorr-Config::minBMcorr)/double(nbins);
		//return double(nevents)*B_PDFs[choice]->EvalPDF(x, par);
		return bin_width*double(nevents)*B_PDFs[choice]->EvalPDF(x, par);
	};
	TF1 *tf1 = new TF1("tf1", funcDraw, Config::minBMcorr, Config::maxBMcorr, nvar);
	tf1->SetParameters(min->X());
	hist->Draw("ep");
	tf1->DrawClone("same");

    std::vector<int> colors = {4,6,7,8,9,30,40,41,38,42,46,28,39};
	for (int i=0; i<B_PDFs[choice]->getComponentCount(); i++) {
    	auto funcDrawComponent = [&B_PDFs, min, nevents, nbins, choice, i](double *x, double *par)->double{
    		double bin_width = (Config::maxBMcorr-Config::minBMcorr)/double(nbins);
    		//return double(nevents)*B_PDFs[choice]->EvalPDF(x, par);
    		return bin_width*double(nevents)*B_PDFs[choice]->EvalPDF(x, par, i);
    	};
        TF1 *tfc = new TF1((std::string("tfc_") + std::to_string(i)).c_str(), funcDrawComponent, Config::minBMcorr, Config::maxBMcorr, nvar);
        tfc->SetParameters(min->X());
        tfc->SetLineColor(colors[i]);
        tfc->DrawClone("same");
    }
	      
	hist->Sumw2();
    	TH1D histpull1D(*hist);
	histpull1D.SetName("histpull1D");
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
	c->SaveAs(Form("%s_figures/B_M_%s_%d.pdf", path_results.Data(), Config::contrName[choice].c_str(), Config::sign));
	
	ofstream outfile(Form("%s/res_%s_%d.txt", path_results.Data(),Config::contrName[choice].c_str(), Config::sign));
   	for (int i =0; i< nvar; i++){
        	outfile << min->X()[i] << "  " << min->Errors()[i] << endl;
   	}
   	outfile.close();

	fit_id++;
	if (hist)
		delete hist;
	if(c)
		delete c;
	}
}
      

