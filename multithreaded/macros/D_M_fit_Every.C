#include "D_M_fit_shape.h"
#include "config.h"

using namespace cpt_b0_analysis;
void D_M_fit_Every(std::string config_file){

// Load config
std::cout<<config_file<<endl;
if (Config::load(config_file)){
    std::cerr<< " Bad config file! " << std::endl;
    return;
}

std::string minName = "Minuit2";
std::string algoName = "";

TString binning = "";
if (Config::binned)
    binning = "binned";
else
    binning = "unbinned";


TString path_results = Form("%s_%s", Config::MC_directory_MD.c_str(), binning.Data());

gSystem->Exec(Form("mkdir -p %s %s_figures", path_results.Data(), path_results.Data()));

int fit_id = 0;
std::cout << Config::int_choose_fits.size() << " size fits " << std::endl;
for (auto& choice: Config::int_choose_fits){
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
	const int nvar = Config::nvar_md;
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

	std::vector<double> vect_Dmass;
	int nbins = 100;
	TH1D *hist = new TH1D("hist", "", nbins, Config::minDM, Config::maxDM);

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
			vect_Dmass.push_back(D_M);
			hist->Fill(D_M);
		}
	}
	double  nevents = double(vect_Dmass.size());
	const auto& D_PDFs = Config::getVectorPDFs("Dmass");

	double step = 0.1;
	for (int ivar=0; ivar<nvar; ivar++){
		step = abs(0.1*Config::init_values[choice][ivar])+0.01;
		min->SetVariable(ivar, (Config::Fits[fit_id] + std::string("_") + Config::varname_md[ivar]).c_str(), Config::init_values[choice][ivar], step);
	}

	//Set Limits on variables
	for (auto it=Config::varLimitsMap.begin(); it!=Config::varLimitsMap.end(); ++it)
	{
		if (it->first.rfind(Config::Fits[fit_id], 0) != 0) continue;
		int index_var = min->VariableIndex(it->first);
		if (index_var == -1 ){
			std::cerr<< "Error in limiting parameters: param " << it->first << " not found.\n";
			return 1;
		}
		auto pairlims = it->second;
		min->SetVariableLimits(index_var, pairlims.first, pairlims.second);
	}

	//list of fixed variables form config
	for (const auto& fix: Config::fixVect){
		if (min->VariableIndex(fix) >= 0) {
			min->FixVariable(min->VariableIndex(fix));
			std::cout << fix << "  " << min->VariableIndex(fix) << std::endl;
		}
	}

// TODO add variable replace

	int ndim = min->NFree();
	auto fchi2 = [&D_PDFs, vect_Dmass, choice](const double *par)->double{
		double chi2 = 0.0;
		D_PDFs[choice]->CalcIntegral(par, Config::minDM, Config::maxDM);
		for (auto dmass: vect_Dmass){
			double likelihood = D_PDFs[choice]->EvalPDF(&dmass, par);
			chi2 -= 2.0*log(likelihood);
		}
		return chi2;
	};
	auto fchi2binned = [&D_PDFs, hist, nevents, nbins, choice, nvar, ndim](const double *par)->double{
		double chi2 = 0;
		D_PDFs[choice]->CalcIntegral(par, Config::minDM, Config::maxDM);
		double bin_width = (Config::maxDM-Config::minDM)/double(nbins);
		   for (int bin=1; bin<=nbins; bin++){
			   double bincenter = hist->GetBinCenter(bin);
			   double bincont = hist->GetBinContent(bin);
			   double err = hist->GetBinError(bin);

			   /// WHY ???
			   double param[nvar];
			   for (int ipar=0; ipar <nvar; ipar++)
				   param[ipar] = par[ipar];
			   //param[2] = par[0];
			   //param[6] = par[4];
			   double estim = bin_width*double(nevents)*D_PDFs[choice]->EvalPDF(&bincenter, param);
			   if(err!=0.0) chi2+=(bincont-estim)*(bincont-estim)/err/err;
		   }
		   return chi2/(nbins-ndim);
	};


	ROOT::Math::Functor f_binned(fchi2binned, nvar);
	ROOT::Math::Functor f(fchi2, nvar);
	if (Config::binned) min->SetFunction(f_binned);
	else min->SetFunction(f);

	double CL_normal = ROOT::Math::normal_cdf(1) -  ROOT::Math::normal_cdf(-1);

	min -> SetErrorDef(ROOT::Math::chisquared_quantile(CL_normal, min->NFree()));
	std::cout<<"Minimize ..."<<std::endl;
	min->Minimize();
	min->Hesse();
	double err_up, err_down;
	//min->GetMinosError(1, err_up, err_down);
	//cout << "Minos mean " << err_up << "  " << err_down << endl;

	// Print correlation matrix
	cout<<"Correlation matrix:"<<endl;
	size_t headerWidths[nvar];
	size_t max_header = 0;
	for (int i=0; i<nvar; i++){
		const auto& name = Config::varname_md[i];
		headerWidths[i] = name.size() > 7 ? name.size() : 7;
		if (max_header < name.size()) max_header = name.size();
	}
	for (int i=0; i<nvar; i++){
		if (i == 0) {
			cout<<std::setw(max_header) << " " <<std::setw(0);
			for (int j=0; j<nvar; j++){
				cout<<" | "<<std::setw(headerWidths[j])<<Config::varname_md[j];
			}
			cout<<endl;
		}
		cout<<std::setw(max_header)<<Config::varname_md[i];
		for (int j=0; j<nvar; j++){
			cout<<std::setw(0)<<" | "<<std::setw(headerWidths[j])<<std::setprecision(3)<<min->Correlation(i, j);
		}
		cout<<endl;
	}

	for (int i=0; i<nvar; i++){
		cout << min->X()[i] << ", ";
	}
	cout << endl;

	TCanvas *c = new TCanvas("c", "", 500, 500);
		TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
		pad1->SetLogy();
		pad1->Draw();
		pad1->cd();

	auto funcDraw = [&D_PDFs, min, nevents, nbins, choice, nvar](double *x, double *par)->double{
		D_PDFs[choice]->CalcIntegral(par, Config::minDM, Config::maxDM);
		double bin_width = (Config::maxDM-Config::minDM)/double(nbins);
			   double param[nvar];
			   for (int ipar=0; ipar <nvar; ipar++)
				   param[ipar] = par[ipar];
			   //param[2] = par[0];
			   //param[6] = par[4];
		return bin_width*double(nevents)*D_PDFs[choice]->EvalPDF(x, param);
	};
	TF1 *tf1 = new TF1("tf1", funcDraw, Config::minDM, Config::maxDM, nvar);
	tf1->SetParameters(min->X());
	cout << tf1->Integral(Config::minDM, Config::maxDM)<< " <<<<=====\n";
	hist->DrawClone("ep");
	tf1->DrawClone("same");

	std::vector<int> colors = {4,6,7,8,9,30,40,41,38,42,46,28,39};
	for (int i=0; i<D_PDFs[choice]->getComponentCount(); i++) {
		auto funcDrawComponent = [&D_PDFs, min, nevents, nbins, choice, i](double *x, double *par)->double{
			double bin_width = (Config::maxDM-Config::minDM)/double(nbins);
			//return double(nevents)*B_PDFs[choice]->EvalPDF(x, par);
			return bin_width*double(nevents)*D_PDFs[choice]->EvalPDF(x, par, i);
		};
		TF1 *tfc = new TF1((std::string("tfc_") + std::to_string(i)).c_str(), funcDrawComponent, Config::minDM, Config::maxDM, nvar);
		tfc->SetParameters(min->X());
		tfc->SetLineColor(colors[i]);
		tfc->DrawClone("same");
	}

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
	c->SaveAs(Form("%s_figures/D_M_%s_%d.pdf", path_results.Data(), Config::contrName[choice].c_str(), Config::sign));
	
	ofstream outfile(Form("%s/res_%s_%d.txt", path_results.Data(), Config::contrName[choice].c_str(), Config::sign));
	outfile << min->Status() << endl;
	outfile << min->MinValue() << endl;
	for (int i =0; i< nvar; i++){
		outfile << min->X()[i] << "  " << min->Errors()[i] << endl;
	}
	outfile.close();

	fit_id++;
}
}


