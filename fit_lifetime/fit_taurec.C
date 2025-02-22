void fit_taurec(){
	TFile *file_r = new TFile("/home/szabelskia/LHCb/MC2016/studies/k_factor/newMC_r_factor.root", "readonly");
        auto tree_r = (TTree*)file_r->Get("BlindedTree");
        bool charge;
        double truetau, tau;
        tree_r->SetBranchAddress("truecharge", &charge);
        tree_r->SetBranchAddress("Tau", &tau);
        tree_r->SetBranchAddress("TrueTau", &truetau);
 
	TH1D *histtau = new TH1D("histtau", "", 20, 1.0, 10.0);
	TH1D *histcount= new TH1D("histtau", "", 20, 1.0, 10.0);
    	for (int i=0; i<tree_r->GetEntries(); i++){
                tree_r->GetEntry(i);
                if (tau<1.0 || tau>10.0)
                        continue;
                if (charge){
			histtau->Fill(tau, tau-truetau);
			histcount->Fill(tau);
		}
			
	}
	histtau->Divide(histcount);
	TF1 *linear = new TF1("linear", "[0]*x*x+[1]*x+[2]", 0, 10);
    	linear->SetParameters(1.0, 0.0, 0.0); 

    // Fit the histogram with the custom function
    histtau->Fit(linear, "R");

	cout << histtau->Integral() << endl;
	TCanvas *c = new TCanvas("c", "", 500, 500);
	histtau->Draw("ep");
	linear->Draw("same");
	c->Update();
}
 
