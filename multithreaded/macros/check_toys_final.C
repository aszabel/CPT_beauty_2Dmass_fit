//void check_fits(int &best, int itoy){
void check_toys_final(string results){
	double min_chi2 = 1.0e45;
	double chi2;
	int status;
	int best = -1;
	const int nfits = 100;
		const int nfrac = 6;
	for (int i=1; i<=nfits; i++){
	//ifstream resin(Form("results/fit2D_%d/results_1_3.txt", i));
	
		double res[(7+7)*(nfrac)+2*nfrac];
		double dres[(7+7)*(nfrac)+2*nfrac];
	ifstream resin(Form("%s/fit2D_%d/results_0_3.txt", results.c_str(), i));
		resin>>status>>chi2;

		int j=0;
		double  x, dx;
  	        while (resin>>x>>dx){
                        res[j] = abs(x);
                        dres[j] = dx;
                        j++;
                }
  
		resin.close();
		
		int n0 = (nfrac) * (7 + 7);
		if ((status==0||status==1) && chi2<min_chi2 && res[n0+1]<0.2 && res[n0+2] <0.2 && res[n0+3]<0.2&& res[n0+5]<0.2){// && i!=10){
			min_chi2 = chi2;
			best = i;
		}
	}
	cout <<best << endl; //"  " << min_chi2 << endl;
		/*ifstream input(Form("%s/fit2D_%d/results_1_3.txt", results.c_str(), best));
                input>>status >> chi2;
		int i=0;
        	while (input>>x>>dx){
                 	res[i] = abs(x);
			dres[i] = dx;  
			i++;
        	}
        	input.close();

		
		cout << "======================================\n";
		for (int ifrac=0; ifrac<nfrac; ifrac++){
			cout << res[n0+ifrac] << " frac" << ifrac << " " << dres[n0+ifrac] << endl << endl;
			
		}*/

}
/*
void check_toys(){

	gStyle->SetOptFit(1);
	const int nfrac = 6;
		double res[(7+7)*(nfrac)+nfrac];
		double dres[(7+7)*(nfrac)+nfrac];
	TCanvas *c[nfrac];
        for (int ifrac=0; ifrac<nfrac; ifrac++){
                c[ifrac] = new TCanvas(Form("c%d", ifrac), "", 500, 500);
	}  
	TH1D *pull[nfrac];
	double gen_frac[nfrac] = {0.50, 0.06, 0.027, 0.06, 0.24, 0.113}; 
	const int ntoys = 100;
	for (int ifrac=0; ifrac<nfrac; ifrac++){
		pull[ifrac] = new TH1D( Form("pull%d", ifrac), "", 100, -0.1, 0.1);
	}
	for (int itoy=1; itoy<=ntoys; itoy++){
		int best=0;
		check_fits(best, itoy);
		ifstream input(Form("toy_res/toy_%d/results/fit2D_%d/results_1_0.txt", itoy, best));
		int status;
		double chi2, x, dx;
                input>>status >> chi2;
		int i=0;
        	while (input>>x>>dx){
                 	res[i] = abs(x);
			dres[i] = dx;  
			i++;
        	}
        	input.close();

		int n0 = (nfrac) * (7 + 7);
		
		cout << "======================================\n";
		for (int ifrac=0; ifrac<nfrac; ifrac++){
			pull[ifrac]->Fill((res[n0+ifrac]-gen_frac[ifrac]));
			cout << res[n0+ifrac] << " frac" << ifrac << " " << dres[n0+ifrac] << endl << endl;
			
		}
	}
	
	for (int ifrac=0; ifrac<nfrac; ifrac++){
		c[ifrac]->cd();
		pull[ifrac]->Draw();
		pull[ifrac]->Fit("gaus");
		c[ifrac]->SaveAs(Form("frac%d.pdf", ifrac));
	}
	cout << " TUUUUUUUUUUUUUUUUUU\n";
}
*/
