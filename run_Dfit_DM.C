void run_Dfit_DM(int sign){
	TStopwatch t;
   	t.Start();
   	gROOT->ProcessLine(".L D_M_fit_shape.cpp+");
   	gROOT->ProcessLine(Form(".x Dfit_DM.C(%d)", sign));
	t.Stop();
 	t.Print();
}
