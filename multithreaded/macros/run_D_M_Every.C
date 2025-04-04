void run_D_M_Every(string config_file){
	TStopwatch t;
   	t.Start();
	gInterpreter->AddIncludePath("../include");
	gInterpreter->AddIncludePath("../../common/include");
   	gROOT->ProcessLine(".L ../src/config.cpp+");
   	gROOT->ProcessLine(".L ../../common/src/ChebyshevPDF.cpp+");
   	gROOT->ProcessLine(".L ../../common/src/D_M_fit_shape.cpp+");
   	gROOT->ProcessLine(Form(".x ../macros/D_M_fit_Every.C(\"%s\")", config_file.c_str()));
	t.Stop();
 	t.Print();
}
