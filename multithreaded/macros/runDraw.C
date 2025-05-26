void runDraw(string config_file){
	TStopwatch t;
   	t.Start();
	gInterpreter->AddIncludePath("../include");
	gInterpreter->AddIncludePath("../../common/include");

        gROOT->ProcessLine(".L ../../common/src/M_B_2missPT_fit.cpp+");
   	gROOT->ProcessLine(".L ../../common/src/ChebyshevPDF.cpp+");
   	gROOT->ProcessLine(".L ../../common/src/D_M_fit_shape.cpp+");
   	gROOT->ProcessLine(".L ../src/config.cpp+");
   	gROOT->ProcessLine(Form(".x ../macros/Draw_details.C(\"%s\")", config_file.c_str()));
	t.Stop();
 	t.Print();
}
