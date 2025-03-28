void run_Dfit_DM(int sign){
	TStopwatch t;
   	t.Start();
	gInterpreter->AddIncludePath("../common/include");
   	gROOT->ProcessLine(".L ../common/src/D_M_fit_shape.cpp+");
   	gROOT->ProcessLine(Form(".x Dfit_DM.C(%d)", sign));
	t.Stop();
 	t.Print();
}
