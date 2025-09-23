void runDraw(string config_file) {
	TStopwatch t;
	t.Start();
	gSystem->AddIncludePath(
		"-I/mnt/opt/spack/0.20/opt/spack/linux-centos7-ivybridge/gcc-12.3.0/"
		"boost-1.82.0-bntx2bbdfbbpoq52tjj5zvl5ddbngidf/include");
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
