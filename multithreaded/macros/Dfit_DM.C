#include "D_M_fit_shape.h"

using namespace cpt_b0_analysis;
   double minx = 1800.;
   double maxx = 1940.;
   double miny = 2700.;
   double maxy = 8300.;
   const int nvar = 7;
   const int ncontr = 6;
   const int nbins = 100;

   TH1D *hist;

void Draw(const double *res, double nevents, TString name[], int sign);
double chebyshev(double *x, const double *par){
  double m_rec = x[0];
  double m_corr = x[1];
  double a1 = par[0];
  double a2 = par[1];
  double cheb = 1.0 + a1*m_rec + a2*(2.0*m_rec*m_rec-1.0);
  double int_cheb = (1.0-a2)*maxx + 0.5*a1*maxx*maxx+2./3.*a2*maxx*maxx*maxx- (1.0-a2)*minx-0.5*a1*minx*minx-2./3.*a2*minx*minx*minx;
  cheb/=int_cheb;
  return cheb;
}

void Dfit_DM(int sign){

        std::string minName = "Minuit2";
        std::string algoName = "";
        ROOT::Math::Minimizer* min =
        ROOT::Math::Factory::CreateMinimizer(minName, algoName);

        // set tolerance , etc...
        min->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
        min->SetMaxIterations(10000);  // for GSL
        min->SetTolerance(1.0e4);
        //min->SetTolerance(1.0e-3);
        min->SetPrintLevel(2);

        //min->SetStrategy(2);
        //min->SetPrecision(0.00001);




        // create funciton wrapper for minmizer
        // a IMultiGenFunction type

        TChain ch("BlindedTree");
   ch.Add("/mnt/home/share/lhcb/CPT_beauty/data2016/selected/selected_data2016MagDown.root");
     double D_M, mu_PT, mu_P, mu_eta, K_PT, B_M, missPT;
   bool charge;


   std::vector<double> vect_Dmass;
   hist = new TH1D("hist", "", nbins, minx, maxx);  
   ch.SetBranchAddress("B_M", &B_M);
   ch.SetBranchAddress("missPT", &missPT);
   ch.SetBranchAddress("D_M", &D_M);
   ch.SetBranchAddress("mu_PT", &mu_PT);
   ch.SetBranchAddress("mu_P", &mu_P);
   ch.SetBranchAddress("mu_eta", &mu_eta);
   ch.SetBranchAddress("K_PT", &K_PT);
   ch.SetBranchAddress("truecharge", &charge);


   for (int i=0; i<ch.GetEntries(); ++i){
   //for (int i=0; i<5.0e4; ++i){
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
        //double  nevents = double(vect_Dmass.size());
        double  nevents = double(hist->Integral());
        md_fit md_shape(minx, maxx);

	double xx[ncontr][nvar];
   	double dxx[ncontr][nvar];
   	TString name[ncontr] = {"signal", "BuDmunu", "BsDsMunu", "B02DpDsm", "sidebands", "Bu2D0Dsm"};
   	for (int ifile=0; ifile<ncontr; ifile++){
		if (ifile==4) continue;
        	ifstream infile(Form("D_M_results/res_%s_%d.txt", name[ifile].Data(), sign));
           	int k=0;
           	while(infile>>xx[ifile][k]>> dxx[ifile][k]){
			xx[4][k]=1.0;
			dxx[4][k]=0.001;
			
                	k++;
           	}
   	}
	xx[4][0] = -0.00102368;
	dxx[4][0] = 2.97034e-07;
	xx[4][1] = 1.31834e-07;
	dxx[4][1] = 7.13695e-11;


        auto fchi2 = [&md_shape, vect_Dmass, xx, dxx](const double *par)->double{
               	const double *pa = &par[ncontr*nvar];
		double frac [ncontr];
		frac[0] = 1-abs(pa[0]);
		frac[1] = abs(pa[0])*(1.0-abs(pa[1]));
    		frac[2] = abs(pa[0])*abs(pa[1])*(1.0-abs(pa[2]));
    		frac[3] = abs(pa[0])*abs(pa[1])*abs(pa[2])*(1.0-abs(pa[3]));
   		frac[4] = abs(pa[0])*abs(pa[1])*abs(pa[2])*abs(pa[3])*(1.0-abs(pa[4]));
    		frac[5] = abs(pa[0])*abs(pa[1])*abs(pa[2])*abs(pa[3])*abs(pa[4]);
		double param[ncontr*nvar];
		for (int i=0; i<ncontr; i++){
			for (int ivar=0; ivar<nvar; ivar++){
				param[i*nvar+ivar] = par[i*nvar+ivar];
				if (ivar==1 && i!=2 && i!=4)param[i*nvar+ivar] = par[1];
			}
		}

 		double chi2 = 0.0;
                for (auto dmass: vect_Dmass){
				double likelihood = 0.0;
			for (int i=0; i<ncontr; i++){
				if (i==4) likelihood+= frac[i]*chebyshev(&dmass, &param[i*nvar]);
				else {
					likelihood += frac[i]*md_shape.func_full(&dmass, &param[i*nvar]);
				}
			}
                        if(likelihood>0.0) chi2 -= 2.0*log(likelihood);
			for (int i=0; i<ncontr; i++){
				if(i==4) continue;
				for (int ivar=0; ivar<nvar; ivar++){
					if (ivar==1) continue;
					if(dxx[i][ivar]!=0)chi2+= (xx[i][ivar]-param[i*nvar+ivar])*(xx[i][ivar]-param[i*nvar+ivar])/(dxx[i][ivar]*dxx[i][ivar]);   // use the results of MC fits
                		}
			}
		}
                return chi2;
        };


	TString varname[nvar] = {"sigma", "mean", "sigmaCB", "f12", "alpha", "n", "alpha_h"};

        for (int i=0; i<ncontr; i++){
		for (int ivar=0; ivar<nvar; ivar++){
			min->SetVariable(i*nvar+ivar, Form("%s_%s", name[i].Data(), varname[ivar].Data()), xx[i][ivar], dxx[i][ivar]+1.0e-11);
			min->SetVariableLowerLimit(i*nvar+ivar, 1e-13);
			if(ivar==3)min->SetVariableLimits(i*nvar+ivar, -1.0, 1.0);

			if (ivar==5) min->FixVariable(i*nvar+ivar);
			if (i==4 && ivar!=0 && ivar!=1) min->FixVariable(i*nvar+ivar); // combinatorial
			if (i==4) min->FixVariable(i*nvar+ivar); // combinatorial
			//if(ivar!=1&&ivar!=0)min->FixVariable(i*nvar+ivar);
			if(i!=0 && ivar==1) min->FixVariable(i*nvar+ivar);
		}
	}

	//std::vector<int> vect_fix= {0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41};
	std::vector<int> vect_fix= {0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41};
	for (auto i: vect_fix) min->FixVariable(i);
	//double frac_init[ncontr-1] = {0.53491, 0.792216, 0.863499, 0.770001, 0.253301};  
	double frac_init[ncontr-1] = {0.657945, 0.901491, -0.968305, 0.527389, 0.0034851};  
	for (int i=0; i<ncontr-1; i++){
		min->SetVariable(nvar*ncontr+i, Form("par_frac%d", i), frac_init[i], 0.01);
		min->SetVariableLimits(nvar*ncontr+i, -1.0, 1.0);
		//min->FixVariable(nvar*ncontr+i);
	}
	ROOT::Math::Functor f(fchi2, nvar*ncontr+ncontr-1);
	min->SetFunction(f);
	double CL_normal = ROOT::Math::normal_cdf(1) -  ROOT::Math::normal_cdf(-1);

        min->SetErrorDef(ROOT::Math::chisquared_quantile(CL_normal, min->NFree()));
        min->Minimize();
	bool over[ncontr*nvar];
	for (int i=0; i<ncontr*nvar; i++){
		over[i]=false;
		for (int j=0; j<ncontr-1; j++){
			double cor = min->Correlation(i, ncontr*nvar+j);
			if (abs(cor) < 0.05) {
				cout << min->VariableName(i) << "  " << cor << endl;
			}else over[i]= true;
		}
	}
	for (int i=0; i<ncontr*nvar; i++) if(!over[i]) cout<< i << ", ";
        cout << endl;	
	Draw(min->X(), nevents, name, sign);
}
void Draw(const double *res, double nevents, TString name[], int sign){        
        md_fit md_shape(minx, maxx);
	TF1 *func[ncontr];
	double frac[ncontr];
        const double *pa = &res[ncontr*nvar];
		frac[0] = 1-abs(pa[0]);
		frac[1] = abs(pa[0])*(1.0-abs(pa[1]));
    		frac[2] = abs(pa[0])*abs(pa[1])*(1.0-abs(pa[2]));
    		frac[3] = abs(pa[0])*abs(pa[1])*abs(pa[2])*(1.0-abs(pa[3]));
   		frac[4] = abs(pa[0])*abs(pa[1])*abs(pa[2])*abs(pa[3])*(1.0-abs(pa[4]));
    		frac[5] = abs(pa[0])*abs(pa[1])*abs(pa[2])*abs(pa[3])*abs(pa[4]);

		double param[ncontr*nvar];
		for (int i=0; i<ncontr; i++){
			for (int ivar=0; ivar<nvar; ivar++){
				param[i*nvar+ivar] = res[i*nvar+ivar];
				if (ivar==1 && i!=2 && i!=4)param[i*nvar+ivar] = res[1];
			}
		}


        for(int i=0; i<ncontr; i++){
		cout << i << "  " << frac [i] << endl;
	}

	double integ = 0.0;
		for (int i=0; i<ncontr; i++){
			auto wrap = [&md_shape, nevents, param, i, &frac] (double *x, double *par)->double{
				double bin_width = (maxx-minx)/double(nbins);
				double fval;
				if (i==4) fval = chebyshev(x, &param[i*nvar]);
				else fval = md_shape.func_full(x, &param[i*nvar]);	
				return nevents*bin_width*frac[i]*fval;
			};
			func[i] = new TF1(Form("tf1_%d", i), wrap, minx, maxx, 0);
			double bin_width = (maxx-minx)/double(nbins);
			integ+=func[i]->Integral(minx, maxx)/bin_width;
			cout << func[i]->Integral(minx, maxx) << " int " << i << endl;
		}
		cout << integ << "  " << nevents << endl;

		auto func_sum = [func](double *x, double *par)->double{
			double sum = 0.0;
			for (int i=0; i<ncontr; i++){
				sum+= func[i]->Eval(x[0]);
			}
			return sum;
		};
		TF1 *tf1_sum = new TF1("tf1_sum", func_sum, minx, maxx, 0);
		cout  << tf1_sum->Integral(minx, maxx) << endl;
	TCanvas *c = new TCanvas("c", "", 500, 500);
	TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
		pad1->SetLogy();
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
	c->SaveAs(Form("DM_fit_MagDown_%d.pdf", sign));
}
