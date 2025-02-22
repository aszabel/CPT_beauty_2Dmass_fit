void runDrawRes(){
	gROOT->ProcessLine(".L pdfs_cpt_model.cpp+");
	gROOT->ProcessLine(".x draw_results.C");
}
