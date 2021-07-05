{
   // gSystem->Load("../build/libPrt.so");
   
  gROOT->ProcessLine(".L ../src/PrtLutNode.cxx+");  
  gROOT->ProcessLine(".L ../src/PrtHit.cxx+");
  gROOT->ProcessLine(".L ../src/PrtEvent.cxx+");
  gROOT->ProcessLine(".L ../src/PrtRun.cxx+");
  gROOT->ProcessLine(".L ../src/PrtTools.cxx+");
}
