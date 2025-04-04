void run_Dfit_BM(int sign){
	TStopwatch t;
   	t.Start();
	gInterpreter->AddIncludePath("../common/include");
   	gROOT->ProcessLine(".L ../common/src/M_B_2missPT_fit.cpp+");
   	gROOT->ProcessLine(Form(".x Dfit_BM.C(%d)", sign));
	t.Stop();
 	t.Print();
}
