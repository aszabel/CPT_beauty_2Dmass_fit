void run_B_M_Every(int icontr, int sign, bool isUnbinned=true){
	TStopwatch t;
	t.Start();
	gInterpreter->AddIncludePath("../../common/include");
	gROOT->ProcessLine(".L ../../common/src/M_B_2missPT_fit.cpp+");
	gROOT->ProcessLine(Form(".x ../macros/B_M_fit_Every.C(%d, %d, %d)", icontr, sign, isUnbinned));
	t.Stop();
	t.Print();
}
