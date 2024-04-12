void runDraw(int sign){
	TStopwatch t;
   	t.Start();
   	gROOT->ProcessLine(".L D_M_fit_shape.cpp+");
	gROOT->ProcessLine(".L M_B_2missPT_fit.cpp+");
   	gROOT->ProcessLine(Form(".x Draw_details.C(%d)", sign));
	t.Stop();
 	t.Print();
}
