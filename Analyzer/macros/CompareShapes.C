{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111111);
  gStyle->SetPalette(55);
  TFile *a = new TFile("CandidatePlots.root","READ");
  TFile *b = new TFile("BHPlots.root","READ");
  TFile *c = new TFile("WPlots.root","READ");
  TFile *d = new TFile("SpikePlots.root","READ");
  TFile *e = new TFile("FakePlots.root","READ");
  TFile *f = new TFile("ZPlots.root","READ");
  

  TH1F *CanTime = (TH1F*)a->Get("Timing");
  TH1F *BHTime = (TH1F*)b->Get("Timing");
  TH1F *WTime = (TH1F*)c->Get("Timing");
  TH1F *SpikeTime = (TH1F*)d->Get("Timing");


  TH1F *CanWidth = (TH1F*)a->Get("SigIetaIeta");
  TH1F *BHWidth = (TH1F*)b->Get("SigIetaIeta");
  TH1F *WWidth = (TH1F*)c->Get("SigIetaIeta");
  TH1F *SpikeWidth = (TH1F*)d->Get("SigIetaIeta");
  TH1F *FakeWidth = (TH1F*)e->Get("SigIetaIeta");


  CanTime->Sumw2();
  BHTime->Sumw2();
  WTime->Sumw2();
  SpikeTime->Sumw2();
  FakeWidth->Sumw2();

  TH1F *data = CanTime->Clone();               //data histogram
  TH1F *mc0 = SpikeTime->Clone();              // first MC histogram
  TH1F *mc1 = BHTime->Clone();                 // second MC histogram
  TH1F *mc2 = WTime->Clone();    // third MC histogram
  
  TH1F *mc2_2 = WWidth->Clone();
  TH1F *mc3 = FakeWidth->Clone();
  
  Float_t NCan = CanTime->GetEntries();
  Float_t NBH = BHTime->GetEntries();
  Float_t NW = WTime->GetEntries();
  Float_t NSp = SpikeTime->GetEntries();
  Float_t NFake = FakeWidth->GetEntries();

  // CanTime->Scale(1./NCan);
  // BHTime->Scale(1./NBH);
  //WTime->Scale(1./NW);
  //SpikeTime->Scale(1./NSp);
  
  //CanTime->Rebin(2);
  //BHTime->Rebin(2);
  //WTime->Rebin(2);
  //SpikeTime->Rebin(2);
  
  TCanvas *can1 = new TCanvas("can1","Timing");
  WTime->SetLineColor(2);
  WTime->SetLineWidth(2);
  BHTime->SetLineColor(4);
  BHTime->SetLineWidth(2);
  SpikeTime->SetLineColor(6);
  SpikeTime->SetLineWidth(2);
  CanTime->SetLineWidth(2);
  WTime->Draw();
  WTime->GetXaxis()->SetRangeUser(-25,25);
  CanTime->Draw("sames");
  BHTime->Draw("sames");
  SpikeTime->Draw("sames");

  TCanvas *can11 = new TCanvas("can11","quad");
  can11->Divide(2,2);
  can11->cd(1);
  CanTime->Draw("histe");
  CanTime->GetXaxis()->SetRangeUser(-25,25);
  can11->cd(2);
  WTime->Draw("histe");
  can11->cd(3);
  BHTime->Draw("histe");
  BHTime->GetXaxis()->SetRangeUser(-25,25);
  can11->cd(4);
  SpikeTime->Draw("histe");
  SpikeTime->GetXaxis()->SetRangeUser(-25,25);

  
  TCanvas *can2= new TCanvas("can2","FractionFit");

  // retrieve histograms
  TObjArray *mc = new TObjArray(3);        // MC histograms are put in this array
  mc->Add(mc0);
  mc->Add(mc1);
  mc->Add(mc2);
  TFractionFitter* fit = new TFractionFitter(data, mc); // initialise
  //fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
  //fit->SetRangeX(1,15);                    // use only the first 15 bins in the fit
  Int_t status = fit->Fit();               // perform the fit
  cout << "fit status: " << status << endl;
  if (status == 0) {                       // check on fit status
    TH1F* result = (TH1F*) fit->GetPlot();
    data->Draw("Ep");
    result->Draw("same");
  
  
    Double_t FracSpikes=0;
    Double_t ErrSpikes=0;
    Double_t FracBH=0;
    Double_t ErrBH=0;
    Double_t FracPrompt=0;
    Double_t ErrPrompt=0;
    fit->GetResult(0, FracSpikes, ErrSpikes);
    fit->GetResult(1, FracBH, ErrBH);
    fit->GetResult(2, FracPrompt, ErrPrompt);

    cout << "FracSpikes: " << FracSpikes << " Est: " << NCan*FracSpikes << endl;
    cout << "FracBH: " << FracBH << " Est: " << NCan*FracBH << endl;
    cout << "FracPrompt: " << FracPrompt << " Est: " << NCan*FracPrompt << endl;
    
    TH1F *sp = mc0->Clone();
    sp->Scale(NCan*FracSpikes/NSp);
    sp->SetLineColor(4);
    sp->Draw("samehist");
    TH1F *bh = mc1->Clone();
    bh->Scale(NCan*FracBH/NBH);
    bh->SetLineColor(6);
    bh->Draw("samehist");
    TH1F *pr = mc2->Clone();
    pr->Scale(NCan*FracPrompt/NW);
    pr->SetLineColor(2);
    pr->Draw("samehist");
    
    TCanvas *can3 = new TCanvas("can","Stack Canvas");
    THStack *stacker = new THStack("stacker","Stack of Time Profiles");
    stacker->Add(sp);
    stacker->Add(bh);
    stacker->Add(pr);
    data->Draw("Ep");
    stacker->Draw("samehist");



    TCanvas *can4 = new TCanvas("can4","Width Canvas"); 
    CanWidth->SetLineWidth(2);
    CanWidth->Draw();
    BHWidth->Scale(NCan*FracBH/NBH);
    BHWidth->SetLineColor(4);
    BHWidth->SetLineWidth(2);
    WWidth->SetLineColor(2);
    WWidth->SetLineWidth(2);
    SpikeWidth->SetLineColor(6);
    SpikeWidth->SetLineWidth(2);
    WWidth->Scale(NCan*FracPrompt/NW);
    SpikeWidth->Scale(NCan*FracSpikes/NSp);

    THStack *stacker2 = new THStack("stacker2","Stack of Widths");
    stacker2->Add(SpikeWidth);
    stacker2->Add(BHWidth);
    //stacker2->Add(WWidth);
    stacker2->Draw("samehist");

    TH1F *CanSub = CanWidth->Clone();
    CanSub->Add(SpikeWidth,-1);
    CanSub->Add(BHWidth, -1);
    TCanvas *can5 = new TCanvas("can5","subtracted");
    CanSub->Draw();
    cout << "Subtracted Entries: " << CanSub->Integral() << endl;
    for (int ik=0;ik<CanSub->GetNbinsX()+1;++ik){
      Float_t Bincon = CanSub->GetBinContent(ik);
      if (Bincon<0){
	CanSub->SetBinContent(ik,0);
      }
    }
    Float_t NCanSub = CanSub->Integral();
    cout << "After zeroing negative bins: " << NCanSub << endl;
    TCanvas *can6= new TCanvas("can6","FractionFit2");

  // retrieve histograms
    TObjArray *mc_2 = new TObjArray(2);        // MC histograms are put in this array
    mc_2->Add(mc2_2);
    mc_2->Add(mc3);
    TFractionFitter* fit2 = new TFractionFitter(CanSub, mc_2); // initialise
    //    fit2->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
    //    fit2->Constrain(0,0.0,1.0);
   //fit->SetRangeX(1,15);                    // use only the first 15 bins in the fit
    Int_t status2 = fit2->Fit();               // perform the fit
    cout << "fit status2: " << status2 << endl;
    if (status2 == 0) {                       // check on fit status
      TH1F* result2 = (TH1F*) fit2->GetPlot();
      CanSub->Draw("Ep");
      result2->Draw("same");

      Float_t FracFake=0;
      Float_t ErrFake=0;
      Float_t FracReal=0;
      Float_t ErrReal=0;
      fit2->GetResult(0, FracReal, ErrReal);
      fit2->GetResult(1, FracFake, ErrFake);
      
      cout << "FracFake: " << FracFake << " Est: " << NCanSub*.417 << endl;
      cout << "FracReal: " << FracReal << " Est: " << NCanSub*.583 << endl;
    
      TCanvas *canShow = new TCanvas("canShow","show dist");
      CanSub->Draw();
      TH1F *fp = mc3->Clone();
      fp->Scale(NCanSub*.417/NFake);
      fp->SetLineColor(4);
      fp->Draw("samehist");
      TH1F *rp = mc2_2->Clone();
      rp->Scale(NCanSub*.583/NW);
      rp->SetLineColor(2);
      rp->Draw("samehist");
    
      TCanvas *can7 = new TCanvas("can7","Stack2 Canvas");
      THStack *stacker3 = new THStack("stacker3","Stack of Width Profiles");
      stacker3->Add(fp);
      stacker3->Add(rp);
      
      CanSub->Draw("Ep");
      stacker3->Draw("samehist");

    }
  }

}
