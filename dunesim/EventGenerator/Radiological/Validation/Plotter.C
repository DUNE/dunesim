//macro to plot the distributions of decay0 on VD geometry

void Plotter(TString inputBKG= "ValAr39GenInLAr"){
  TString inFile = Form("BackgroundValidation_%s.root",inputBKG.Data());
  
  TFile *fInput = new TFile(inFile.Data());
  if(fInput->IsOpen())
    Int_t a=0;
  else{
    cout<<"fail to open file. Skipping..."<<endl;
  }
  
  Double_t Energy, Momentum, StartX, StartY, StartZ;
  Int_t    PDG, Event;

  TTree *tSNSim = (TTree*)fInput->Get(Form("%s/BackgroundValidation",inputBKG.Data()));
  tSNSim->SetBranchAddress("Energy"        ,&Energy       );
  tSNSim->SetBranchAddress("Momentum"      ,&Momentum     );
  tSNSim->SetBranchAddress("StartX"        ,&StartX       );
  tSNSim->SetBranchAddress("StartY"        ,&StartY       );
  tSNSim->SetBranchAddress("StartZ"        ,&StartZ       );
  tSNSim->SetBranchAddress("PDG"           ,&PDG          );
  tSNSim->SetBranchAddress("Event"         ,&Event        );

  TH1D *hX   = new TH1D("hX"  ,"X Position",100, -500, 500 ); hX->SetLineWidth(3);
  TH1D *hY   = new TH1D("hY"  ,"Y Position",100, -900, 900 ); hY->SetLineWidth(3);
  TH1D *hZ   = new TH1D("hZ"  ,"Z Position",100, -200, 1200); hZ->SetLineWidth(3);
  TH1D *hKEe = new TH1D("hKEe","Electron"  ,500, 0   , 10  ); hKEe->SetLineWidth(3);
  TH1D *hKEg = new TH1D("hKEg","Gamma"     ,500, 0   , 10  ); hKEg->SetLineWidth(3); hKEg->SetLineColor(kBlack);
  TH1D *hKEa = new TH1D("hKEa","Alpha"     ,500, 0   , 10  ); hKEa->SetLineWidth(3); hKEa->SetLineColor(kRed);
  TH1D *hKEn = new TH1D("hKEn","Neutron"   ,500, 0   , 10  ); hKEn->SetLineWidth(3); hKEn->SetLineColor(kMagenta);

  TH2D *hYZ   = new TH2D("hYZ"  ,"",100, -200, 1200, 100, -900, 900);
  TH2D *hYX   = new TH2D("hYX"  ,"",100, -500, 500,  100, -900, 900);
  TH2D *hXZ   = new TH2D("hXZ"  ,"",100, -200, 1200, 100, -500, 500);

  Int_t nevt = tSNSim->GetEntries();

  Double_t maxEe = 0, maxEa = 0, maxEg = 0;

  for(int i=0;i<nevt;++i){
    tSNSim->GetEntry(i);
    //cout<<Energy<<" "<<Momentum<<" "<<sqrt(Energy*Energy-Momentum*Momentum)*1e3<<endl;

    hX->Fill(StartX);
    hY->Fill(StartY);
    hZ->Fill(StartZ);
    hYZ->Fill(StartZ,StartY);
    hYX->Fill(StartX,StartY);
    hXZ->Fill(StartZ,StartX);
    
    Double_t KE = Energy*1e3 - TMath::Sqrt(Energy*Energy - Momentum*Momentum)*1e3;
    if(PDG == 11){
      hKEe->Fill(KE);
      if(KE>maxEe) maxEe=KE;
    }
    else if(PDG == 22){
      hKEg->Fill(KE);
      if(KE>maxEg) maxEg=KE;
    }
    else if(PDG == 1000020040){
      hKEa->Fill(KE);
      if(KE>maxEa) maxEa=KE;
    }
    else if(PDG == 2112){
      hKEn->Fill(KE);
    }
    else
      cout<<"PDG "<< PDG <<" not known!"<<endl;
  }

  gStyle->SetOptStat("");
  gStyle->SetPadTickX(kTRUE);
  gStyle->SetPadTickY(kTRUE);
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  gStyle->SetNumberContours(99);

  hKEe->Scale(1./Event);
  hKEg->Scale(1./Event);
  hKEa->Scale(1./Event);
  hKEn->Scale(1./Event);
  hX->Scale(1./Event);
  hY->Scale(1./Event);
  hZ->Scale(1./Event);
  hYZ->Scale(1./Event);
  hYX->Scale(1./Event);
  hXZ->Scale(1./Event);

  TString outPDF = Form("Plots_%s.pdf",inputBKG.Data());

  TCanvas c("c","");
  c.Print(Form("%s[",outPDF.Data()));
  Double_t KEs[] = {maxEe, maxEg, maxEa}; Double_t maxE = *std::max_element(std::begin(KEs), std::end(KEs));
  hKEn->GetXaxis()->SetRangeUser(0,maxE*1.1); hKEn->Draw("hist"); hKEn->GetXaxis()->SetTitle("Kinectic Energy (MeV)"); hKEn->GetYaxis()->SetTitle("Decays/Event");
  //hKEa->GetXaxis()->SetRangeUser(0,maxE*1.1); hKEa->Draw("hist"); hKEa->GetXaxis()->SetTitle("Kinectic Energy (MeV)"); hKEa->GetYaxis()->SetTitle("Decays/Event");
  hKEa->GetXaxis()->SetRangeUser(0,maxEg*1.1); hKEa->Draw("hist same");
  hKEg->GetXaxis()->SetRangeUser(0,maxEg*1.1); hKEg->Draw("hist same");
  hKEe->GetXaxis()->SetRangeUser(0,maxEe*1.1); hKEe->Draw("hist same");
  gPad->BuildLegend(); hKEa->SetTitle("");
  c.Print(outPDF.Data());

  hX->Draw("hist"); hX->GetXaxis()->SetTitle("X (cm)"); hX->GetYaxis()->SetTitle("Decays/Event"); c.Print(outPDF.Data());
  hY->Draw("hist"); hY->GetXaxis()->SetTitle("Y (cm)"); hY->GetYaxis()->SetTitle("Decays/Event"); c.Print(outPDF.Data());
  hZ->Draw("hist"); hZ->GetXaxis()->SetTitle("Z (cm)"); hZ->GetYaxis()->SetTitle("Decays/Event"); c.Print(outPDF.Data());

  hYX->Draw("colz");  hYX->GetXaxis()->SetTitle("X (cm)"); hYX->GetYaxis()->SetTitle("Y (cm)"); c.Print(outPDF.Data());
  hYZ->Draw("colz");  hYZ->GetXaxis()->SetTitle("Z (cm)"); hYZ->GetYaxis()->SetTitle("Y (cm)"); c.Print(outPDF.Data());

  hXZ->Draw("colz");  hXZ->GetXaxis()->SetTitle("Z (cm)"); hXZ->GetYaxis()->SetTitle("X (cm)"); c.Print(outPDF.Data());

  c.Print(Form("%s]",outPDF.Data()));

  TFile f1(Form("Plots_%s.root",inputBKG.Data()),"recreate");
  hX->Write("hX");
  hY->Write("hY");
  hZ->Write("hZ");
  hKEe->Write("hKEe");
  hKEg->Write("hKEg");
  hKEa->Write("hKEa");
  hYX->Write("hYX");
  hYX->Write("hYZ");
  hXZ->Write("hXZ");

  
}
