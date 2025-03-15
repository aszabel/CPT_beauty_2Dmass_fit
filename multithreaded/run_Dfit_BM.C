void run_Dfit_BM(int sign){
	TStopwatch t;
   	t.Start();
   	gROOT->ProcessLine(".L M_B_2missPT_fit.cpp+");
   	gROOT->ProcessLine(Form(".x Dfit_BM.C(%d)", sign));
	t.Stop();
 	t.Print();
}
