
void lutmean(TString baseFile = "../data/lut"){  
  TString inFile =baseFile+".root";
  TString outFile =baseFile+"_avr.root";

  TFile* f = new TFile(inFile);
  TTree *t=(TTree *) f->Get("prtlut") ;
  TClonesArray* fLut = new TClonesArray("PrtLutNode");
  t->SetBranchAddress("LUT",&fLut);
  t->GetEntry(0);
 
 
  TFile *fFileNew = TFile::Open(outFile, "RECREATE");
  TClonesArray *fLutNew;

  TTree *fTreeNew = new TTree("prtlut","Look-up table for DIRC. Averaged");

  fLutNew = new TClonesArray("PrtLutNode");
  fTreeNew->Branch("LUT",&fLutNew,256000,2); 

  Int_t Nnodes = 100000;
  TClonesArray &fLutaNew = *fLutNew;
  for (Long64_t n=0; n<Nnodes; n++) {
    new((fLutaNew)[n]) PrtLutNode(-1);
  }

  // TCanvas* c = new TCanvas("c","c",0,0,800,1200); c->Divide(1,2);
  // TH1F * histNode = new TH1F("LutNode","Node vs Multiplicity",30000,0,150000);
  // TH1F * hTime = new TH1F("hTime","Time",5000,0,10);
  // TH1F * hDir = new TH1F("hDir","X component",1000,-1,1);


  std::vector<TVector3> vArray[100];
  std::vector<Double_t> tArray[100];
  std::vector<Double_t> pArray;
  
  TVector3 dir, dir2, sum;
  Double_t angle, minangle,pathid,time,sumt;
  PrtLutNode *node;
  
  for (Int_t inode=0; inode<fLut->GetEntriesFast(); inode++){
    if(inode%1000==0) cout<<"Node # "<<inode<<endl;
    node= (PrtLutNode*) fLut->At(inode);
    //histNode->Fill(node->GetNodeId(),node->Entries());
    Int_t size = node->Entries();
    if(size<1) continue;
    for(int i=0; i<size; i++){
      dir = node->GetEntry(i);
      time = node->GetTime(i);
      pathid = node->GetPathId(i);
      
      // hDir->Fill(dir.X());
      // hTime->Fill(time);

      bool newid = true;
      for(int j=0; j<pArray.size(); j++){
	if(pathid == pArray[j]){
	  vArray[j].push_back(dir);
	  tArray[j].push_back(time);
	  newid= false;
	}
      }
      if(newid) {
	vArray[pArray.size()].push_back(dir);
	tArray[pArray.size()].push_back(time);
	pArray.push_back(pathid);
      }
    }
  
    for(int j=0; j<pArray.size(); j++){
      sum = TVector3(0,0,0);
      sumt=0;
      for(int v=0; v<vArray[j].size(); v++) {
	sum += vArray[j][v]; 
	sumt += tArray[j][v]; 

	// hDir->Fill(vArray[j][v].X());
	// hTime->Fill(tArray[j][v]);
      }
      
      // c->cd(1);
      // hTime->Draw();
      // c->cd(2);
      // hDir->Draw();
      // c->Update();  
      // c->WaitPrimitive();
      // hDir->Reset();
      // hTime->Reset();

      //if(vArray[j].size()<5) continue;
      sum *= 1/(Double_t)vArray[j].size();
      sumt *= 1./(Double_t)tArray[j].size();
      std::cout<<"inode "<<inode << " "<<node->GetDetectorId() <<std::endl;
      
      ((PrtLutNode*)(fLutNew->At(inode)))->AddEntry(node->GetDetectorId(), sum,j,j,sumt, node->GetDigiPos(),node->GetDigiPos()); 
    }
    for(int i=0; i<100; i++) {vArray[i].clear();  tArray[i].clear();}
    pArray.clear();
  }

  fTreeNew->Fill();
  fTreeNew->Write();
  fFileNew->Write();

}
