void run_B_M_Every(string config_file){
	TStopwatch t;
	t.Start();
	gInterpreter->AddIncludePath("../include");
	gInterpreter->AddIncludePath("../../common/include");
   	gROOT->ProcessLine(".L ../src/config.cpp+");
	gROOT->ProcessLine(".L ../../common/src/M_B_2missPT_fit.cpp+");
   	gROOT->ProcessLine(Form(".x ../macros/B_M_fit_Every.C(\"%s\")", config_file.c_str()));
	t.Stop();
	t.Print();
}
