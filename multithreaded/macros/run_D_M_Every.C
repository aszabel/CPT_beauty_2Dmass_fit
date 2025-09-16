void run_D_M_Every(string config_file, string install_path="../.."){
	TStopwatch t;
	t.Start();
	gSystem->AddIncludePath("-I/mnt/opt/spack/0.20/opt/spack/linux-centos7-ivybridge/gcc-12.3.0/boost-1.82.0-bntx2bbdfbbpoq52tjj5zvl5ddbngidf/include");
	gInterpreter->AddIncludePath(Form("%s/multithreaded/include", install_path.c_str()));
	gInterpreter->AddIncludePath(Form("%s/common/include", install_path.c_str()));
	gROOT->ProcessLine(Form(".L %s/common/src/M_B_2missPT_fit.cpp+", install_path.c_str()));
	gROOT->ProcessLine(Form(".L %s/common/src/D_M_fit_shape.cpp+", install_path.c_str()));
	gROOT->ProcessLine(Form(".L %s/common/src/BasicShapes.cpp+", install_path.c_str()));
	gROOT->ProcessLine(Form(".L %s/common/src/ChebyshevPDF.cpp+", install_path.c_str()));
	gROOT->ProcessLine(Form(".L %s/multithreaded/src/config.cpp+", install_path.c_str()));
	gROOT->ProcessLine(Form(".x %s/multithreaded/macros/D_M_fit_Every.C(\"%s\")", install_path.c_str(), config_file.c_str()));
	t.Stop();
	t.Print();
}
