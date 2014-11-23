#include <iostream>
#include <fstream>


void writeInfo(TString path, TString info){
  ofstream myfile;
  myfile.open (path+"/"+"readme");
  myfile << info+"\n";
  myfile.close();
}

void save(TPad *c= NULL, TString dir="rdata", TString name="", TString info="", Int_t flag=0){
  if(flag==-1) return;
  TString path = "";
  if(flag == 0){
    TDatime *time = new TDatime();
    TString stime = Form("%d.%d.%d", time->GetDay(),time->GetMonth(),time->GetYear()); 
    Int_t btime = gSystem->mkdir(dir+"/"+stime);
    for(Int_t i=0; i<1000; i++){
      path = dir+"/"+stime+"/"+Form("arid-%d",i);
      if(gSystem->mkdir(path)==0) break;
    }
    writeInfo(path, info);
  }else{
    path = dir;
  }
  
  if(c) {
    c->Modified();
    c->Update();
    c->Print(path+"/"+name+".png");
    c->Print(path+"/"+name+".pdf");
    c->Print(path+"/"+name+".root");
  }
}

TString createDir(TString dir="rdata", TString info = "", Int_t flag=0){
  if(flag==-1) return "";
  TDatime *time = new TDatime();
  TString stime = Form("%d.%d.%d", time->GetDay(),time->GetMonth(),time->GetYear()); 
  TString path = "";
  Int_t btime = gSystem->mkdir(dir+"/"+stime);
  for(Int_t i=0; i<1000; i++){
    path = dir+"/"+stime+"/"+Form("arid-%d",i);
    if(gSystem->mkdir(path)==0) break;
  }
  writeInfo(path, info);
  return path;
}
