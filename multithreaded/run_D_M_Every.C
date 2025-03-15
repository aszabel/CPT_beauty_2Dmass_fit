void run_D_M_Every(int icontr, int sign, bool isUnbinned=true){
	if (icontr==4) return;
	TStopwatch t;
   	t.Start();
   	gROOT->ProcessLine(".L D_M_fit_shape.cpp+");
   	gROOT->ProcessLine(Form(".x D_M_fit_Every.C(%d, %d, %d)", icontr, sign, isUnbinned));
	t.Stop();
 	t.Print();
}
