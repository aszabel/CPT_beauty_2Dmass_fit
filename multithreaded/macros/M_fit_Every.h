#include "config.h"

using namespace cpt_b0_analysis;
void M_fit_Every(
	const TString& path_results,
	const int& nvar,
	const double& minM,
	const double& maxM,
	const TString& fit_type,
	const std::vector<std::string>& varnames,
	const int& nbins=100
){

std::string minName = "Minuit2";
std::string algoName = "";

TString binning = "";
if (Config::binned)
	binning = "binned";
else
	binning = "unbinned";

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

	std::vector<double> vect_mass;
	TH1D *hist = new TH1D("hist", "", nbins, minM, maxM);

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
			switch (Config::int_mass_variable){
				case 0: // D_M
					vect_mass.push_back(D_M);
					hist->Fill(D_M);
					break;
				case 1: // B_M
					vect_mass.push_back(B_M);
					hist->Fill(B_M);
					break;
				case 2: // B_Mcorr
					vect_mass.push_back(B_M+2.0*missPT);
					hist->Fill(B_M+2.0*missPT);
					break;
				case 3: // B_MMcorr
					vect_mass.push_back(B_MMcorr);
					hist->Fill(B_MMcorr);
					break;
				case 4: // missPT
					vect_mass.push_back(missPT);
					hist->Fill(missPT);
					break;
				default:
					std::cerr << "Unsupported mass variable " << Config::mass_variable << std::endl;
					return;
			}
		}
	}
	double  nevents = double(vect_mass.size());
	const auto& PDFs = Config::getVectorPDFs(fit_type.Data());

	double step = 0.1;
	for (int ivar=0; ivar<nvar; ivar++){
		/*if (ivar== 2 || ivar == 3 || ivar == 5 || ivar == 7) 
			step = 0.01; 
		else
			step = 10.0;*/
		//TODO step is to small ?? Fit tends to keep to initial values
		//step = abs(0.01*Config::init_values[choice][ivar])+0.01;
		step = abs(0.1*Config::init_values[choice][ivar])+0.01;
		//TODO name variables after fit name so that we can fix, set limits per contribution
		min->SetVariable(ivar, (Config::Fits[fit_id] + std::string("_") + varnames[ivar]).c_str(), Config::init_values[choice][ivar], step);
		//min -> FixVariable(ivar);		
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

	std::vector<std::pair<int, int>> replaceIndexVect = {};
	for (const auto& rep_var: Config::replace_var){
		if (rep_var.first.rfind(Config::Fits[fit_id], 0) != 0) continue;
		int index_replaced = min->VariableIndex(rep_var.first);
		int index_substitute = min->VariableIndex(rep_var.second);
		if (index_replaced == -1 ){
			std::cerr<< "Error in substituting parameters param " << rep_var.first << " not found.\n";
			return 1;
		}
		if (index_substitute == -1 ){
			std::cerr<< "Error in substituting parameters param " << rep_var.second << " not found.\n";
			return 1;
		}
		replaceIndexVect.push_back(std::make_pair(index_replaced, index_substitute));
	}

	int ndim = min->NFree();
	auto fchi2 = [&PDFs, nvar, vect_mass, minM, maxM, choice, replaceIndexVect](const double *par)->double {
		double param[nvar];
		for (int ivar = 0; ivar < nvar; ivar++)
		{
			param[ivar] = par[ivar];
		}
		for (const auto& irep_var: replaceIndexVect){
			param[irep_var.first] = par[irep_var.second];
		}

		double chi2 = 0.0;
		PDFs[choice]->CalcIntegral(param, minM, maxM);
		for (auto bmass: vect_mass){
			double likelihood = PDFs[choice]->EvalPDF(&bmass, param);
			if (likelihood>=1.0|| likelihood <=0.0)
			{
				continue;
			}
			chi2 -= 2.0*log(likelihood);
		}
		return chi2;
	};
	auto fchi2binned = [&PDFs, nvar, hist, nevents, nbins, minM, maxM, choice, ndim, replaceIndexVect](const double *par)->double{
		double param[nvar];
		for (int ivar = 0; ivar < nvar; ivar++)
		{
			param[ivar] = par[ivar];
		}
		for (const auto& irep_var: replaceIndexVect){
			param[irep_var.first] = par[irep_var.second];
		}

		double chi2 = 0;
		PDFs[choice]->CalcIntegral(param, minM, maxM);
		double bin_width = (maxM-minM)/double(nbins);
		for (int bin=1; bin<=nbins; bin++){
			double bincenter = hist->GetBinCenter(bin);
			double bincont = hist->GetBinContent(bin);
			double err = hist->GetBinError(bin);
			double estim = bin_width*double(nevents)*PDFs[choice]->EvalPDF(&bincenter, param);
			if(err!=0.0) chi2+=(bincont-estim)*(bincont-estim)/err/err;
		}
		return chi2/double(nbins-ndim);
	};


	ROOT::Math::Functor f_binned(fchi2binned, nvar);
	ROOT::Math::Functor f(fchi2, nvar);
	if (Config::binned) min->SetFunction(f_binned);
	else min->SetFunction(f);
	
	double CL_normal = ROOT::Math::normal_cdf(1) -  ROOT::Math::normal_cdf(-1);

	min -> SetErrorDef(ROOT::Math::chisquared_quantile(CL_normal, min->NFree()));
	min->Minimize();
	min->Hesse();
	//double err_up, err_down;
	//min->GetMinosError(1, err_up, err_down);
	//cout << "Minos mean " << err_up << "  " << err_down << endl;

/*
Contour

TGraph g(n); 
minimizer->Contour(ipar, jpar, n, g.GetX(), g.GetY() ); 
g.Draw("AC");
*/

	// Print correlation matrix
	cout<<"Correlation matrix:"<<endl;
	size_t headerWidths[nvar];
	size_t max_header = 0;
	for (int i=0; i<nvar; i++){
		const auto& name = varnames[i];
		headerWidths[i] = name.size() > 7 ? name.size() : 7;
		if (max_header < name.size()) max_header = name.size();
	}
	for (int i=0; i<nvar; i++){
		if (i == 0) {
			cout<<std::setw(max_header) << " " <<std::setw(0);
			for (int j=0; j<nvar; j++){
				cout<<" | "<<std::setw(headerWidths[j])<<varnames[j];
			}
			cout<<endl;
		}
		cout<<std::setw(max_header)<<varnames[i];
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

	auto funcDraw = [&PDFs, nvar, min, nevents, nbins, minM, maxM, choice, replaceIndexVect](double *x, double *par)->double
	{
		double param[nvar];
		for (int ivar = 0; ivar < nvar; ivar++)
		{
			param[ivar] = par[ivar];
		}
		for (const auto& irep_var: replaceIndexVect){
			param[irep_var.first] = par[irep_var.second];
		}

		PDFs[choice]->CalcIntegral(param, minM, maxM);
		double bin_width = (maxM-minM)/double(nbins);
		//return double(nevents)*PDFs[choice]->EvalPDF(x, par);
		return bin_width*double(nevents)*PDFs[choice]->EvalPDF(x, param);
	};
	TF1 *tf1 = new TF1("tf1", funcDraw, minM, maxM, nvar);
	tf1->SetParameters(min->X());
	hist->Draw("ep");
	tf1->DrawClone("same");

	std::vector<int> colors = {4,6,7,8,9,30,40,41,38,42,46,28,39};
	for (int i=0; i<PDFs[choice]->getComponentCount(); i++) {
		auto funcDrawComponent = [&PDFs, nvar, min, nevents, nbins, minM, maxM, choice, i, replaceIndexVect](double *x, double *par)->double{
			double param[nvar];
			for (int ivar = 0; ivar < nvar; ivar++)
			{
				param[ivar] = par[ivar];
			}
			for (const auto& irep_var: replaceIndexVect){
				param[irep_var.first] = par[irep_var.second];
			}

			double bin_width = (maxM-minM)/double(nbins);
			//return double(nevents)*PDFs[choice]->EvalPDF(x, par);
			return bin_width*double(nevents)*PDFs[choice]->EvalPDF(x, param, i);
		};
		TF1 *tfc = new TF1((std::string("tfc_") + std::to_string(i)).c_str(), funcDrawComponent, minM, maxM, nvar);
		tfc->SetParameters(min->X());
		tfc->SetLineColor(colors[i]);
		tfc->DrawClone("same");
	}

	hist->Sumw2();
	TH1D histpull1D(*hist);
	histpull1D.SetName("histpull1D");
	for (int bin=1; bin<=hist->GetNbinsX(); bin++) {
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
	c->SaveAs(Form("%s_figures/%s_%s_%d.pdf", path_results.Data(), fit_type.Data(), Config::contrName[choice].c_str(), Config::sign));

	ofstream outfile(Form("%s/res_%s_%d.txt", path_results.Data(),Config::contrName[choice].c_str(), Config::sign));
	outfile << min->Status() << endl;
	outfile << min->MinValue() << endl;
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

// vim: tabstop=4 softtabstop=0 noexpandtab shiftwidth=4
