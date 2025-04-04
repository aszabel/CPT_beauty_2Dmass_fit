void runDrawRes(){
	gInterpreter->AddIncludePath("../include");
	gROOT->ProcessLine(".L ../src/pdfs_cpt_model.cpp+");
	gROOT->ProcessLine(".x draw_results.C");
}
