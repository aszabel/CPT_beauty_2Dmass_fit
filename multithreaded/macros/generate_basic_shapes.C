#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"
#include "TRandom.h"

void generate_basic_shapes() {
	double D_M;
    double B_M;
    double mu_PT = 1.0;
    double mu_P = 1.0;
    double mu_eta = 1.0;
    double K_PT = 1.0;
    double missPT = 0.0;
    double Tau = 1.0;
	bool charge = 0;

    TFile *file = file = TFile::Open("gauss.root","RECREATE");

    TTree* tree = new TTree("BlindedTree","");
    tree->Branch("B_M", &B_M, "B_M/D");
    tree->Branch("missPT", &missPT, "missPT/D");
    tree->Branch("D_M", &D_M, "D_M/D");
    tree->Branch("mu_PT", &mu_PT, "mu_PT/D");
    tree->Branch("mu_P", &mu_P, "mu_P/D");
    tree->Branch("mu_eta", &mu_eta, "mu_eta/D");
    tree->Branch("K_PT", &K_PT, "K_PT/D");
    tree->Branch("Tau", &Tau, "Tau/D");
    tree->Branch("truecharge", &charge, "truecharge/O");

    int type = 1;

    if (type == 0) {
        for(int i=0; i<100000; i++) {
            B_M = gRandom->Gaus(100,1);
            D_M = B_M;
            charge = 0;
            tree->Fill();
            charge = 1;
            tree->Fill();
        }
    } else if(type == 1) {
        for(int i=0; i<100000; i++) {
            // Signal
            B_M = gRandom->Gaus(200,10);
            D_M = gRandom->Gaus(100,1);
            charge = 0;
            tree->Fill();
            charge = 1;
            tree->Fill();
        }
        for(int i=0; i<50000; i++) {
            // BG
            B_M = gRandom->Uniform() * 2 * 5 * 10 + 200 - 5*10;
            D_M = gRandom->Uniform() * 2 * 5 * 1 + 100 - 5*1;
            charge = 0;
            tree->Fill();
            charge = 1;
            tree->Fill();
        }
    }
    tree->Write();
    file->Close();
}