#include "D_M_fit_shape.h"

using namespace cpt_b0_analysis;
   double minx = 1800.;
   double maxx = 1940.;
   double miny = 2700.;
   double maxy = 8300.;
   const int nvar_md = 7;
   const int nvar_mb = 7;
   const int ncontr = 6;
   const int nbins = 40;

   bool uncorr = false;
   TH2D *hist;

void Draw(const double *res, double nevents, TString name[], int sign);
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

void Dfit_2D(int sign){

        std::string minName = "Minuit2";
        std::string algoName = "";
        ROOT::Math::Minimizer* min =
        ROOT::Math::Factory::CreateMinimizer(minName, algoName);

        // set tolerance , etc...
        min->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
        min->SetMaxIterations(10000);  // for GSL
        min->SetTolerance(2.50e5);
        //min->SetTolerance(1.0e-3);
        min->SetPrintLevel(2);

        //min->SetStrategy(2);
        //min->SetPrecision(0.00001);




        // create funciton wrapper for minmizer
        // a IMultiGenFunction type

        TChain ch("BlindedTree");
   //ch.Add("/mnt/home/share/lhcb/CPT_beauty/data2016/selected/selected_data2016MagDown.root");
   ch.Add("/home/szabelskia/LHCb/data2016/tree_missPT_D_M_MagDown25102023_nomassDmuCut/selected_data2016MagDown.root");
     double D_M, mu_PT, mu_P, mu_eta, K_PT, B_M, missPT;
   bool charge;


   std::vector<std::pair<double, double>> vect_2D;
   hist = new TH2D("hist", "", nbins, minx, maxx, nbins, miny, maxy);  
   ch.SetBranchAddress("B_M", &B_M);
   ch.SetBranchAddress("missPT", &missPT);
   ch.SetBranchAddress("D_M", &D_M);
   ch.SetBranchAddress("mu_PT", &mu_PT);
   ch.SetBranchAddress("mu_P", &mu_P);
   ch.SetBranchAddress("mu_eta", &mu_eta);
   ch.SetBranchAddress("K_PT", &K_PT);
   ch.SetBranchAddress("truecharge", &charge);


   //for (int i=0; i<ch.GetEntries(); ++i){
   for (int i=0; i<5.0e3; ++i){
      ch.GetEntry(i);
      if (mu_PT<500 || mu_P<5000 || mu_eta<2 || mu_eta>4.5) continue;
      double B_MMcorr = B_M +2.0*missPT;
      if (B_MMcorr<miny || B_MMcorr>maxy) continue;
      if (D_M<minx || D_M>maxx) continue;
      if (int(charge) == sign){
              vect_2D.push_back(std::make_pair(D_M, B_MMcorr));
              hist->Fill(D_M, B_MMcorr);
      }
   }
        //double  nevents = double(vect_2D.size());
        double  nevents = double(hist->Integral());
        md_fit md_shape(minx, maxx);
	mb_2misspt_fit mb_shape(miny, maxy);

	double xx[ncontr][nvar_md];
   	double dxx[ncontr][nvar_md];
   	TString name[ncontr] = {"signal", "BuDmunu", "BsDsMunu", "B02DpDsm", "sidebands", "Bu2D0Dsm"};
   	for (int ifile=0; ifile<ncontr; ifile++){
		if (ifile==4) continue;
        	ifstream infile(Form("D_M_results/res_%s_%d.txt", name[ifile].Data(), sign));
           	int k=0;
           	while(infile>>xx[ifile][k]>> dxx[ifile][k]){
			cout << xx[ifile][k] << endl;
			xx[4][k]=1.0;
			dxx[4][k]=0.001;
			
                	k++;
           	}
   	}
	xx[4][0] = 1.88518e-08;
	dxx[4][0] = 2.97034e-07;
	xx[4][1] = 1.31834e-07;
	dxx[4][1] = 7.13695e-11;

	double xx_mcorr[ncontr][nvar_mb];
	double dxx_mcorr[ncontr][nvar_mb];
	for (int ifile=0; ifile<ncontr; ifile++){
                ifstream infile_mb(Form("B_M_results/res_%s_%d.txt", name[ifile].Data(), sign));
                int k=0;
                while(infile_mb>>xx_mcorr[ifile][k]>> dxx_mcorr[ifile][k]){
                        k++;
                }
        }


        auto fchi2 = [&md_shape, &mb_shape, vect_2D, xx, dxx, xx_mcorr, dxx_mcorr, min](const double *par)->double{
               	const double *pa = &par[ncontr*(nvar_md+nvar_mb)];
		double frac [ncontr];
		frac[0] = 1-abs(pa[0]);
		frac[1] = abs(pa[0])*(1.0-abs(pa[1]));
    		frac[2] = abs(pa[0])*abs(pa[1])*(1.0-abs(pa[2]));
    		frac[3] = abs(pa[0])*abs(pa[1])*abs(pa[2])*(1.0-abs(pa[3]));
   		frac[4] = abs(pa[0])*abs(pa[1])*abs(pa[2])*abs(pa[3])*(1.0-abs(pa[4]));
    		frac[5] = abs(pa[0])*abs(pa[1])*abs(pa[2])*abs(pa[3])*abs(pa[4]);
		double param[ncontr*(nvar_md+nvar_mb)];
		for (int i=0; i<ncontr; i++){
			for (int ivar=0; ivar<nvar_md; ivar++){
				param[i*nvar_md+ivar] = par[i*nvar_md+ivar];
				if (ivar==1 && i!=2 && i!=4)param[i*nvar_md+ivar] = par[1];
			}
		}
	       for (int i=0; i<ncontr; i++){
                        for (int ivar=0; ivar<nvar_mb; ivar++){
				int ipar = ncontr*nvar_md+i*nvar_mb+ivar;
                                param[ipar] = par[ipar];
                        }
                }




 		double chi2 = 0.0;
		int count = 0;
                for (auto &[mdass, mcorr]: vect_2D){
			if (++count>5.0e4 && !uncorr) break;
			double likelihood = 0.0;
			for (int i=0; i<ncontr; i++){
				double md_like, mb_like;
				if (i==4) md_like = chebyshev(&mdass, &param[i*nvar_md]);
				else md_like = md_shape.func_full(&mdass, &param[i*nvar_md]);
				mb_like = mb_shape.func_full(&mcorr, &param[ncontr*nvar_md+i*nvar_mb]);
				likelihood+= md_like*mb_like*frac[i];
			}
                        if(likelihood>0.0) chi2 -= 2.0*log(likelihood);
		}
		for (int i=0; i<ncontr; i++){
			if(i==4) continue;
			for (int ivar=0; ivar<nvar_md; ivar++){
				//if (ivar==1) continue;
				if(dxx[i][ivar]!=0)chi2+= (xx[i][ivar]-param[i*nvar_md+ivar])*(xx[i][ivar]-param[i*nvar_md+ivar])/(dxx[i][ivar]*dxx[i][ivar]);   // use the results of MC fits
                	}
		}
		for (int i=0; i<ncontr; i++){
                        for (int ivar=0; ivar<nvar_mb; ivar++){
				//cout << min->VariableName(ncontr*nvar_md+i*nvar_mb+ivar) << "  " << xx_mcorr[i][ivar] << "  " << dxx_mcorr[i][ivar] << endl;
                        	if(dxx_mcorr[i][ivar]!=0)chi2+= (xx_mcorr[i][ivar]-param[ncontr*nvar_md+i*nvar_mb+ivar])*(xx_mcorr[i][ivar]-param[ncontr*nvar_md+i*nvar_mb+ivar])/(dxx_mcorr[i][ivar]*dxx_mcorr[i][ivar]);   // use the results of MC fits
                        }
                }
		
                return chi2;
        };


	TString varname_md[nvar_md] = {"sigma", "mean", "sigmaCB", "f12", "alpha", "n", "alpha_h"};

        for (int i=0; i<ncontr; i++){
		for (int ivar=0; ivar<nvar_md; ivar++){
			min->SetVariable(i*nvar_md+ivar, Form("%s_%s", name[i].Data(), varname_md[ivar].Data()), xx[i][ivar], dxx[i][ivar]+1.0e-11);
			//min->FixVariable(i*nvar_md+ivar);
			min->SetVariableLowerLimit(i*nvar_md+ivar, 1e-13);
			if(ivar==3)min->SetVariableLimits(i*nvar_md+ivar, -1.0, 1.0);

			cout << ivar << endl;
			if ( ivar==6) min->FixVariable(i*nvar_md+ivar);
			if (ivar==5) min->FixVariable(i*nvar_md+ivar);
			if (i==4 && ivar!=0 && ivar!=1) min->FixVariable(i*nvar_md+ivar); // combinatorial
			//if (i==4) min->FixVariable(i*nvar_md+ivar); // combinatorial
			//if(ivar!=1&&ivar!=0)min->FixVariable(i*nvar_md+ivar);
			if(i!=0 && ivar==1) min->FixVariable(i*nvar_md+ivar);
		}
	}


	
	
	
	TString varname_mb[nvar_mb] = {"mean_rc", "sigma_rc", "f12_gaus", "mean1_gaus", "sigma1_gaus", "mean2_gaus", "sigma2_gaus"};

        for (int i=0; i<ncontr; i++){
                for (int ivar=0; ivar<nvar_mb; ivar++){
                        min->SetVariable(ncontr*nvar_md+i*nvar_mb+ivar, Form("%s_%s", name[i].Data(), varname_mb[ivar].Data()), xx_mcorr[i][ivar], dxx_mcorr[i][ivar]+1.0e-11);
			//min->FixVariable(ncontr*nvar_md+i*nvar_mb+ivar);
                        min->SetVariableLowerLimit(ncontr*nvar_md+i*nvar_mb+ivar, 1e-13);
                        if(ivar==2)min->SetVariableLimits(ncontr*nvar_md+i*nvar_mb+ivar, -1.0, 1.0);

                }
        }



      	//      for (auto i: vect_fix_mb) min->FixVariable(ncontr*nvar_md+i);
        double frac_init[ncontr-1] = {0.468636, 0.736568, 0.973902, 0.949009, -0.00476566};
	double frac_init_plus[ncontr-1] = {0.578309, 0.816153, 0.519927, 0.627151, 0.0151929}; 
        //double frac_init[ncontr-1] = {0.657945, 0.901491, 0.968305, 0.527389, 0.0034851};
        for (int i=0; i<ncontr-1; i++){

                if (sign==0) min->SetVariable((nvar_md+nvar_mb)*ncontr+i, Form("par_frac%d", i), frac_init[i], 0.01);
		else min->SetVariable((nvar_md+nvar_mb)*ncontr+i, Form("par_frac%d", i), frac_init_plus[i], 0.01);

                min->SetVariableLimits((nvar_md+nvar_mb)*ncontr+i, -1.0, 1.0);
                //min->FixVariable((nvar_md+nvar_mb)*ncontr+i);
        }
	ROOT::Math::Functor f(fchi2, (nvar_md+nvar_mb)*ncontr+ncontr-1);
	min->SetFunction(f);
	double CL_normal = ROOT::Math::normal_cdf(1) -  ROOT::Math::normal_cdf(-1);

        min->SetErrorDef(ROOT::Math::chisquared_quantile(CL_normal, min->NFree()));

	double x_fix[ncontr*(nvar_md+nvar_mb)+ncontr-1];
	double dx_fix[ncontr*(nvar_md+nvar_mb)+ncontr-1];
	ifstream input_fix(Form("results_fix_%d.txt", sign));
	int k=0;
	while (input_fix >> x_fix[k] >> dx_fix[k]) k++; 
	input_fix.close();
	uncorr = true;
	if (k==ncontr*(nvar_md+nvar_mb)+ncontr-1) cout << "TUUUUUUU \n\n";

	
 	for(int ivar=0; ivar<ncontr*(nvar_md+nvar_mb)+ncontr-1; ivar++){ 
	//	cout << x_fix[ivar] << "  " << ivar << endl;
		//min->FixVariable(ivar);
		//min->SetVariableValue(ivar, x_fix[ivar]);
	}
        min->Minimize();
	bool over[ncontr*(nvar_md+nvar_mb)];
	for (int i=0; i<ncontr*(nvar_md+nvar_mb); i++){
		over[i]=false;
		for (int j=0; j<ncontr-1; j++){
			double cor = min->Correlation(i, ncontr*(nvar_md+nvar_mb)+j);
			if (abs(cor) < 0.02) {
				cout << min->VariableName(i) << "  " << cor << endl;
				
			}else over[i]= true;
		}
	}
	std::vector<int> vect_fix = {5, 8, 12, 15, 19, 22, 26, 29, 30, 31, 32, 33, 34, 36, 40, 57, 62};
	for (int i=0; i<ncontr*(nvar_md+nvar_mb)+ncontr-1; i++) min->SetVariableValue(i, x_fix[i]);
	uncorr = true;	
	for (auto i: vect_fix) if (uncorr) min->FixVariable(i);

	//min->Minimize();

	ofstream results(Form("results_%d.txt", sign));
	for (int i=0; i<(nvar_md+nvar_mb)*ncontr+ncontr-1; i++){
		results<< min->X()[i]<< "  " << min->Errors()[i] << endl;
	}
	results.close();
	Draw(min->X(), nevents, name, sign);
}
void Draw(const double *res, double nevents, TString name[], int sign){        
        md_fit md_shape(minx, maxx);
	mb_2misspt_fit mb_shape(miny, maxy);
	TF2 *func[ncontr];
	double frac[ncontr];
        const double *pa = &res[ncontr*(nvar_md+nvar_mb)];
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

	double integ = 0.0;
		for (int i=0; i<ncontr; i++){
			auto wrap = [&md_shape, &mb_shape, nevents, param, i, &frac] (double *x, double *par)->double{
				double bin_widthx = (maxx-minx)/double(nbins);
				double bin_widthy = (maxy-miny)/double(nbins);
				double fval_md, fval_mb;
				if (i==4) fval_md = chebyshev(x, &param[i*nvar_md]);
				else fval_md = md_shape.func_full(x, &param[i*nvar_md]);
				fval_mb	= mb_shape.func_full(&x[1], &param[ncontr*nvar_md+i*nvar_mb]);
				return nevents*bin_widthx*bin_widthy*frac[i]*fval_md*fval_mb;
			};
			func[i] = new TF2(Form("tf1_%d", i), wrap, minx, maxx, miny, maxy, 0);
			double bin_widthx = (maxx-minx)/double(nbins);
			double bin_widthy = (maxy-miny)/double(nbins);
			integ+=func[i]->Integral(minx, maxx, miny, maxy)/bin_widthx/bin_widthy;
			cout << func[i]->Integral(minx, maxx, miny, maxy) << " int " << i << endl;
		}
		cout << integ << "  " << nevents << endl;

		auto func_sum = [func](double *x, double *par)->double{
			double sum = 0.0;
			for (int i=0; i<ncontr; i++){
				sum+= func[i]->Eval(x[0], x[1]);
			}
			return sum;
		};
		TF2 *tf2_sum = new TF2("tf2_sum", func_sum, minx, maxx, miny, maxy, 0);
		cout  << tf2_sum->Integral(minx, maxx) << endl;
	TCanvas *c = new TCanvas("c", "", 500, 500);

	hist->SetStats(kFALSE);
	hist->SetMinimum(1.0);
	hist->DrawClone("ep");
	tf2_sum->SetLineColor(kRed);
	tf2_sum->DrawClone("same surf");

	TCanvas *cpull  = new TCanvas("cpull", "", 500, 500);
        TH2D histpull2D(*hist);
	histpull2D.SetMinimum();
        for (int binx=1; binx<=hist->GetNbinsX(); binx++){
        	for (int biny=1; biny<=hist->GetNbinsY(); biny++){
                	double err = hist->GetBinError(binx, biny);
                	if (err==0.0) continue;
                	double diff  = hist->GetBinContent(binx, biny)-tf2_sum->Eval(hist->GetXaxis()->GetBinCenter(binx), hist->GetYaxis()->GetBinCenter(biny));
                	histpull2D.SetBinContent(binx, biny, diff/err);
		}
        }
        histpull2D.SetStats(kFALSE);
        histpull2D.SetFillColor(kBlue);
        //histpull2D.GetYaxis()->SetLabelSize(0.1);
        //histpull2D.GetXaxis()->SetLabelSize(0.1);
        histpull2D.DrawClone("surf");
}
