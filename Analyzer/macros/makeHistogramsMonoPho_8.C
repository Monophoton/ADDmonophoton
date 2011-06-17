#define makeHistogramsMonoPho_8_cxx
#include "makeHistogramsMonoPho_8.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
/*
 upgrade from
 _backup1: gonna add sigieie, and timing and phi histos before and after cut
 _backup2: gonna add inv. mass and MT plots
 _backup3: gonna add spike/nonspike Et, time distros
 _backup4: gonna take halo out of monophoton sample
 _backup5: gonna examine spike phi, eta, run range
 _backup6: gonna examine cosmics
 _backup7: gonna examine Zee and Wenu to get number of W faking monophoton, and plot MET
 makeHistogramsW: add bookkeeping histos
 makeHistogramsW_2: update selections to match with Andrew (19 April 2011)
 makeHistogramsW_5: update photon id to include sig ieie
 */

void makeHistogramsMonoPho_8::Loop(char subset, string outName, bool applyTriggerCut, bool applyCosMuCut, bool applyScraping, bool onlyOneVertex)
{
    //   In a ROOT session, you can do:
    //      Root > .L makeHistogramsMonoPho_8.C
    //      Root > makeHistogramsMonoPho_8 t
    //      Root > t.GetEntry(12); // Fill t data members with entry number 12
    //      Root > t.Show();       // Show values of entry 12
    //      Root > t.Show(16);     // Read and show values of entry 16
    //      Root > t.Loop();       // Loop on all entries
    //
    
    //     This is the loop skeleton where:
    //    jentry is the global entry number in the chain
    //    ientry is the entry number in the current Tree
    //  Note that the argument to GetEntry must be:
    //    jentry for TChain::GetEntry
    //    ientry for TTree::GetEntry and TBranch::GetEntry
    //
    //       To read only selected branches, Insert statements like:
    // METHOD1:
    //    fChain->SetBranchStatus("*",0);  // disable all branches
    //    fChain->SetBranchStatus("branchname",1);  // activate branchname
    // METHOD2: replace line
    //    fChain->GetEntry(jentry);       //read all branches
    //by  b_branchname->GetEntry(ientry); //read only this branch
    cout<<"subset="<<subset<<endl;
    cout<<"outName="<<outName<<endl;
    cout<<"applyCosMuCut="<<applyCosMuCut<<endl;
    cout<<"applyTriggerCut="<<applyTriggerCut<<endl;
    
    if (fChain == 0) return;
    
    TFile *outputFile = new TFile(outName.data(),"RECREATE");
    
    //backup1
    TH1F *h_pho_n = new TH1F("h_pho_n","Number of Photons per Event; Number of Photons; Events",20,0,20);
    TH1F *h_pho_E = new TH1F("h_pho_E","Photons; E; Events",1000,0,500);
    TH1F *h_pho_pt = new TH1F("h_pho_pt","Photons; p_{T}; Events",1000,0,500);
    TH1F *h_pho_phi = new TH1F("h_pho_phi","Photons; #phi; Events",1000,0,6.3);
    TH1F *h_pho_eta = new TH1F("h_pho_eta","Photons; #eta; Events",1000,-3,3);
    
    //backup2
    TH1F *h_pho_sieie = new TH1F("h_pho_sieie","Photons; #sigma_{i#etai#eta}; Events",1000,0,0.03);
    TH1F *h_pho_sipip = new TH1F("h_pho_sipip","Photons; #sigma_{i#phii#phi}; Events",1000,0,0.03);
    TH1F *h_pho_tSeed = new TH1F("h_pho_tSeed","Photons; t_{seed}; Events",1000,-20,20);
    
    //backup3
    TH1F *h_tcmet_pt = new TH1F("h_pfmet_pt","PFMet; p_{T}; Events",5000,0,500);
    TH1F *h_tcmet_phi = new TH1F("h_pfmet_phi","PFMet; #phi; Events",1000,0,6.3);
    
    //backup7
    TH1F *h_sigmaR = new TH1F("h_sigmaR"," Events;#sigmaR;",10000,-1,1000);
    TH1F *h_deltaD = new TH1F("h_deltaD"," Events;#deltaD;",10000,-1,1000);
    //sigR<20 -> 99.3% of signal retained, 93% bg suppressed
    
    //makeHistogramsW
    TH1F *h_massT = new TH1F("h_massT","Transverse Mass of Pho Candidates;M_{T} GeV;",4000,0,400);
    TH1F *h_cands = new TH1F("h_cands","cands passing req.; electrons; Events",10,0,10);
    
    //new
    TH2F *h_nJet_pT = new TH2F("h_nJet_pT","Jets; p_{T Cut}; N_{jets}",500,0,500,100,0,100);
    TH2F *h_nJet_pT_eta25 = new TH2F("h_nJet_pT_eta25","Jets |#eta|<2.5; p_{T Cut}; N_{jets}",500,0,500,100,0,100);
    TH2F *h_nTrk_pT = new TH2F("h_nTrk_pT","Tracks; p_{T Cut}; N_{tracks}",500,0,500,200,0,200);
    TH2F *h_nTrk_pT_eta20 = new TH2F("h_nTrk_pT_eta20","Tracks |#eta|<2.0; p_{T Cut}; N_{tracks}",500,0,500,200,0,200);
    TH2F *h_nTrk_pT_eta25 = new TH2F("h_nTrk_pT_eta25","Tracks |#eta|<2.5; p_{T Cut}; N_{tracks}",500,0,500,200,0,200);
    TH2F *h_recoPhoPix_pT_barrel = new TH2F("h_recoPhoPix_pT_barrel","recoPho w/ pixel seed match;p_{T Cut};N_{recoPho w/ pix}",500,0,500,100,0,100);
    
    TH1F *h_mu_isGlobal = new TH1F("h_mu_isGlobal","Muons (Global);N_{#mu}",100,0,100);
    TH1F *h_mu_isTracker = new TH1F("h_mu_isTracker","Muons (Tracker);N_{#mu}",100,0,100);
    TH1F *h_mu_isSA = new TH1F("h_mu_isSA","Muons (Stand Alone);N_{#mu}",100,0,100);
    
    TH1F *h_cosmu_isGlobal = new TH1F("h_cosmu_isGlobal","Cosmic Muons (Global);N_{#mu}",100,0,100);
    TH1F *h_cosmu_isTracker = new TH1F("h_cosmu_isTracker","Cosmic Muons (Tracker);N_{#mu}",100,0,100);
    TH1F *h_cosmu_isSA = new TH1F("h_cosmu_isSA","Cosmic Muons (Stand Alone);N_{#mu}",100,0,100);
    TH1F *h_e2e9 = new TH1F("h_e2e9","E2/E9;E2/E9",300,-1,2);
    
    //newer
    TH1F *h_dPhi_PhoMet = new TH1F("h_dPhi_PhoMet","d#phi(phos, PFMet);d#phi",1000,0,4);
    TH1F *h_sumEt = new TH1F("h_sumEt","sum Et;sum Et",1000,0,1000);
    TH1F *h_leadJet_leadPho_deltaR = new TH1F("h_leadJet_leadPho_deltaR","d#phi(lead Jet,lead Pho);d#phi",1000,0,8);
    TH1F *h_dPhi_lJetMet = new TH1F("h_dPhi_lJetMet","d#phi(lead Jet, PFMet);d#phi",1000,0,4);
    
    //more shit
    TH1F *h_roundness = new TH1F("h_roundness","Roundness;Roundness",1000,0,1);
    TH1F *h_LITCD = new TH1F("h_LITCD","Photons; LITCD; Events",1000,-5,5);
    TH2F *h_seed_iphi_ieta = new TH2F("h_seed_iphi_ieta","iphi v ieta",364,-2,362,174,-87,87);
    TH1F *h_lJet_pt = new TH1F("h_lJet_pt","Leading Jet; p_{T}; Events",1000,0,500);
    TH1F *h_VertexZ = new TH1F("h_VertexZ","Vertex Z (pre-cut at #pm25); Z cm; Events",1000,-50,50);
    TH1F *h_trkIso = new TH1F("h_trkIso","Track Iso;trkSumPtHollowConeDR04",1000,0,10);
    TH1F *h_ecalIso = new TH1F("h_ecalIso","Ecal Iso;ecalRecHitSumEtConeDR04",1000,0,10);
    TH1F *h_hcalIso = new TH1F("h_hcalIso","Hcal Iso;hcalTowerSumEtConeDR04",1000,0,10);
    TH1F *h_HoE = new TH1F("h_HoE","Had/Em;H/E",1000,0,0.5);
    
    //Andrew request 5-16-2011
    TH1F *n_ifVeto_none = new TH1F("n_ifVeto_none","Number of Original Cands",1,0,1);
    TH1F *n_ifVeto_tpho = new TH1F("n_ifVeto_tpho","Cands left after veto tight photon",1,0,1);
    TH1F *n_ifVeto_lpho = new TH1F("n_ifVeto_lpho","Cands left after veto loose photon",1,0,1);
    TH1F *n_ifVeto_em = new TH1F("n_ifVeto_em","Cands left after veto em",1,0,1);
    TH1F *n_ifVeto_em20 = new TH1F("n_ifVeto_em20","Cands left after veto em p_{T}>20",1,0,1);
    TH1F *n_ifVeto_em30 = new TH1F("n_ifVeto_em30","Cands left after veto em p_{T}>30",1,0,1);
    TH1F *n_ifVeto_tpho_barrel = new TH1F("n_ifVeto_tpho_barrel","Cands left after veto tight barrel photon",1,0,1);
    TH1F *n_ifVeto_lpho_barrel = new TH1F("n_ifVeto_lpho_barrel","Cands left after veto loose barrel photon",1,0,1);
    TH1F *n_ifVeto_em_barrel = new TH1F("n_ifVeto_em_barrel","Cands left after veto barrel em",1,0,1);
    TH1F *n_ifVeto_em20_barrel = new TH1F("n_ifVeto_em20_barrel","Cands left after veto barrel em p_{T}>20",1,0,1);
    TH1F *n_ifVeto_em30_barrel = new TH1F("n_ifVeto_em30_barrel","Cands left after veto barrel em p_{T}>30",1,0,1);
    
    TH1F *n_ifVeto_cosmu = new TH1F("n_ifVeto_cosmu","Cands left after veto cosmic #mu",1,0,1);
    TH1F *n_ifVeto_mu = new TH1F("n_ifVeto_mu","Cands left after veto #mu",1,0,1);
    TH1F *n_ifVeto_cosmu_or_mu = new TH1F("n_ifVeto_cosmu_or_mu","Cands left after veto #mu or cosmic #mu",1,0,1);
    TH1F *n_ifVeto_muWithVertex = new TH1F("n_ifVeto_muWithVertex","Cands left after veto #mu with matching vertex",1,0,1);
    
    TH1F *n_ifVeto_trk5 = new TH1F("n_ifVeto_trk5","Cands left after veto tracks p_{T}>5",1,0,1);
    TH1F *n_ifVeto_trk20 = new TH1F("n_ifVeto_trk20","Cands left after veto tracks p_{T}>20",1,0,1);
    TH1F *n_ifVeto_trk5_eta20 = new TH1F("n_ifVeto_trk5_eta20","Cands left after veto tracks p_{T}>5 within |#eta|<2.0",1,0,1);
    TH1F *n_ifVeto_trk20_eta20 = new TH1F("n_ifVeto_trk20_eta20","Cands left after veto tracks p_{T}>20 within |#eta|<2.0",1,0,1);
    TH1F *n_ifVeto_trk5_eta25 = new TH1F("n_ifVeto_trk5_eta25","Cands left after veto tracks p_{T}>5 within |#eta|<2.5",1,0,1);
    TH1F *n_ifVeto_trk20_eta25 = new TH1F("n_ifVeto_trk20_eta25","Cands left after veto tracks p_{T}>20 within |#eta|<2.5",1,0,1);
    TH1F *n_ifVeto_trk10 = new TH1F("n_ifVeto_trk10","Cands left after veto tracks p_{T}>10",1,0,1);
    TH1F *n_ifVeto_trk10_eta20 = new TH1F("n_ifVeto_trk10_eta20","Cands left after veto tracks p_{T}>10 within |#eta|<2.0",1,0,1);
    TH1F *n_ifVeto_trk10_eta25 = new TH1F("n_ifVeto_trk10_eta25","Cands left after veto tracks p_{T}>10 within |#eta|<2.5",1,0,1);

    
    TH1F *n_ifVeto_jet50 = new TH1F("n_ifVeto_jet50","Cands left after veto jets p_{T}>50",1,0,1);
    TH1F *n_ifVeto_jet50_eta25 = new TH1F("n_ifVeto_jet50_eta25","Cands left after veto jets p_{T}>50 within |#eta|<2.5",1,0,1);
    TH1F *n_ifVeto_jet30 = new TH1F("n_ifVeto_jet30","Cands left after veto jets p_{T}>30",1,0,1);
    TH1F *n_ifVeto_jet30_eta25 = new TH1F("n_ifVeto_jet30_eta25","Cands left after veto jets p_{T}>30 within |#eta|<2.5",1,0,1);
    
    TH1F *n_ifVeto_trk20_or_jet50 = new TH1F("n_ifVeto_trk20_or_jet50","Cands left after veto (tracks p_{T}>20 or jets p_{T}>50)",1,0,1);
    TH1F *n_ifVeto_trk20_eta20_or_jet50_eta25 = new TH1F("n_ifVeto_trk20_eta20_or_jet50_eta25","Cands left after veto (tracks p_{T}>20 (within |#eta|<2.0) or jets p_{T}>50 (within |#eta|<2.5))",1,0,1);
    TH1F *n_ifVeto_trk20_eta25_or_jet50_eta25 = new TH1F("n_ifVeto_trk20_eta25_or_jet50_eta25","Cands left after veto (tracks p_{T}>20 (within |#eta|<2.5) or jets p_{T}>50 (within |#eta|<2.5))",1,0,1);
    
    TH1F *n_ifVeto_trk20_or_jet30 = new TH1F("n_ifVeto_trk20_or_jet30","Cands left after veto (tracks p_{T}>20 or jets p_{T}>30)",1,0,1);
    TH1F *n_ifVeto_trk10_or_jet30 = new TH1F("n_ifVeto_trk10_or_jet30","Cands left after veto (tracks p_{T}>10 or jets p_{T}>30)",1,0,1);
    TH1F *n_ifVeto_trk20_eta25_or_jet30_eta25 = new TH1F("n_ifVeto_trk20_eta25_or_jet30_eta25","Cands left after veto (tracks p_{T}>20 (within |#eta|<2.5) or jets p_{T}>30 (within |#eta|<2.5))",1,0,1);
    TH1F *n_ifVeto_trk10_eta25_or_jet30_eta25 = new TH1F("n_ifVeto_trk10_eta25_or_jet30_eta25","Cands left after veto (tracks p_{T}>10 (within |#eta|<2.5) or jets p_{T}>30 (within |#eta|<2.5))",1,0,1);

    
    TH1F *n_ifVeto_em20_or_mu = new TH1F("n_ifVeto_em_or_mu","Cands left after veto (EM p_{T}>20 or muons)",1,0,1);
    TH1F *n_ifVeto_em20_or_cosmu = new TH1F("n_ifVeto_em_or_cosmu","Cands left after veto (EM p_{T}>20 or cosmic muons)",1,0,1);
    TH1F *n_ifVeto_em20_or_mu_or_cosmu = new TH1F("n_ifVeto_em_or_mu_or_cosmu","Cands left after veto (EM p_{T}>20 or #mu or cosmic #mu)",1,0,1);
    
    TH1F *n_ifVeto_em20_or_trk20 = new TH1F("n_ifVeto_em20_or_trk20","Cands left after veto (EM p_{T}>20 or tracks p_{T}>20)",1,0,1);
    TH1F *n_ifVeto_em20_or_jet50 = new TH1F("n_ifVeto_em20_or_jet50","Cands left after veto (EM p_{T}>20 or jets p_{T}>50)",1,0,1);
    TH1F *n_ifVeto_em20_or_trk20_or_jet50 = new TH1F("n_ifVeto_em20_or_trk20_or_jet50","Cands left after veto (EM p_{T}>20 or tracks p_{T}>20 or jets p_{T}>50)",1,0,1);
    
    TH1F *n_ifVeto_em20_or_trk20_or_jet50_or_mu_or_cosmu = new TH1F("n_ifVeto_em20_or_trk20_or_jet50_or_mu_or_cosmu","Cands left after veto (EM p_{T}>20 or tracks p_{T}>20 or jets p_{T}>50 or #mu or cosmic #mu)",1,0,1);
    TH1F *n_ifVeto_em20_barrel_or_trk20_eta25_or_jet50_eta25_or_mu_or_cosmu = new TH1F("n_ifVeto_em20_barrel_or_trk20_eta25_or_jet50_eta25_or_mu_or_cosmu","Cands left after veto (barrel EM p_{T}>20 or tracks p_{T}>20 (within |#eta|<2.5) or jets p_{T}>50 (within |#eta|<2.5) or #mu or cosmic #mu)",1,0,1);
    TH1F *n_ifVeto_em_or_trk10_or_jet30 = new TH1F("n_ifVeto_em_or_trk10_or_jet30","Cands left after veto ( EM p_{T}>10 or tracks p_{T}>10 or jets p_{T}>30)",1,0,1);
    TH1F *n_ifVeto_em_or_trk10_eta25_or_jet30_eta25 = new TH1F("n_ifVeto_em_or_trk10_eta25_or_jet30_eta25","Cands left after veto ( EM p_{T}>10 or tracks p_{T}>10 (within |#eta|<2.5) or jets p_{T}>30 (within |#eta|<2.5))",1,0,1);
    TH1F *n_ifVeto_em_or_trk10_or_jet30_or_mu_or_cosmu = new TH1F("n_ifVeto_em_or_trk10_or_jet30_or_mu_or_cosmu","Cands left after veto ( EM p_{T}>10 or tracks p_{T}>10 or jets p_{T}>30 or #mu or cosmic #mu)",1,0,1);
    TH1F *n_ifVeto_em_or_trk10_eta25_or_jet30_eta25_or_mu_or_cosmu = new TH1F("n_ifVeto_em_or_trk10_eta25_or_jet30_eta25_or_mu_or_cosmu","Cands left after veto ( EM p_{T}>10 or tracks p_{T}>10 (within |#eta|<2.5) or jets p_{T}>30 (within |#eta|<2.5) or #mu or cosmic #mu)",1,0,1);

    TH1F *h_vertexDistance_pho_pho = new TH1F("h_vertexDistance_pho_pho","Distance between candidate pho and other pho; cm; Events",1000,0,100);
    TH1F *h_vertexDistance_pho_trk = new TH1F("h_vertexDistance_pho_trk","Distance between candidate pho and other Track; cm; Events",1000,0,100);
    TH1F *h_vertexDistance_pho_jet = new TH1F("h_vertexDistance_pho_jet","Distance between candidate pho and other Jet; cm; Events",1000,0,100);
    TH1F *h_vertexDistance_pho_mu  = new TH1F("h_vertexDistance_pho_mu" ,"Distance between candidate pho and other #mu; cm; Events",1000,0,100);
    
    TH1F *h_vertexDistance_pho_pho_noR = new TH1F("h_vertexDistance_pho_pho_noR","Distance between candidate pho and other pho ignoring #deltaR; cm; Events",1000,0,5);
    TH1F *h_vertexDistance_pho_trk_noR = new TH1F("h_vertexDistance_pho_trk_noR","Distance between candidate pho and other Track ignoring #deltaR; cm; Events",1000,0,5);
    TH1F *h_vertexDistance_pho_jet_noR = new TH1F("h_vertexDistance_pho_jet_noR","Distance between candidate pho and other Jet ignoring #deltaR; cm; Events",1000,0,5);
    TH1F *h_vertexDistance_pho_mu_noR  = new TH1F("h_vertexDistance_pho_mu_noR" ,"Distance between candidate pho and other #mu ignoring #deltaR; cm; Events",1000,0,5);
    
    
    Long64_t nentries = fChain->GetEntriesFast();
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        
        //for (int i=0; i<(*triggernames).size();i++)
        //{
           // cout << (*triggernames)[i] <<endl;
        //}
        
        //select events with the triggers I want
        vector <string> myDesiredTrigger;
        myDesiredTrigger.push_back("HLT_Photon75_CaloIdVL_IsoL_v1");
        myDesiredTrigger.push_back("HLT_Photon75_CaloIdVL_IsoL_v2");
        myDesiredTrigger.push_back("HLT_Photon75_CaloIdVL_IsoL_v3");
        if(applyTriggerCut && false==triggerPassed(ntriggers,myDesiredTrigger,*triggernames,*ifTriggerpassed)) continue;
        //cout<<"event = "<<event<<endl;

        
        //select nonscraping events
	if(applyScraping && Scraping_isScrapingEvent) continue;
        
        for(int v=0; v<Vertex_n && v<10; v++){
            h_VertexZ->Fill(Vertex_z[v]);
        }
        //select events that have a good vertex
        if(!eventHasGoodVertex(Vertex_n,Vertex_ndof,Vertex_z,Vertex_d0,Vertex_isFake,onlyOneVertex)) continue;
        //remove events with cosmic muon objects
        if(applyCosMuCut && eventHasMu(CosmicMuon_n, CosmicMuon_OuterTrack_isNonnull, CosmicMuon_isStandAloneMuon)) continue;
        
        
        //ID photons
        vector<int> indexTightBarrelPho;
	
	//get photon vectors
	// these are IDs for EGM-10-006
	for(int x=0;x<Photon_n && x<200;x++){
            //remove bad photons
            if(!passesPhotonID(Photon_ncrys[x], Photon_pt[x],Photon_sc_eta[x],Photon_HoE[x],
                Photon_trkSumPtHollowConeDR04[x],Photon_ecalRecHitSumEtConeDR04[x],Photon_hcalTowerSumEtConeDR04[x],
                Photon_SigmaIetaIeta[x])) continue;
            //require LICTD small
            if(fabs(computeLICTD(Photon_ncrys[x], Photon_timing_xtal[x], Photon_energy_xtal[x]))>5.) continue;
            //require pixel seed veto
            if(Photon_hasPixelSeed[x]) continue;
            //remove spikes
            if(Photon_SigmaIetaIeta[x]<0.001) continue;
            if(Photon_SigmaIphiIphi[x]<0.001) continue;
            //if(GetE2overE9(Photon_ncrys[x], Photon_ieta_xtalEB[x], Photon_iphi_xtalEB[x], Photon_energy_xtal[x])>0.95) continue;
            //remove halos
            //if(isCSCHalo(Photon_sc_phi[x], CSCseg_n, CSCseg_x, CSCseg_y, CSCseg_time, false)) continue;
            //remove out of time bg
            if(fabs(Photon_timing_xtal[x][0])>3.) continue;
            
            
            indexTightBarrelPho.push_back(x);
            
        }// get photon IDs
        
        int isCandidate=0; 
        int leadPhoIdx=-1;
        float leadPhoPt=-1;
        float leadMETphi=-1;
        float leadMET=-1;      
        //tight barrel phos
        for(unsigned int index=0; index<indexTightBarrelPho.size(); index++){
            int i = indexTightBarrelPho[index];
            /*
            float seedT = Photon_timing_xtal[i][GetSeed(Photon_ncrys[i], Photon_energy_xtal[i])];
            float MET = newMetPt(PFMetPt[0], PFMetPhi[0], Photon_pt[i], Photon_phi[i],seedT);
            float METphi = newMetPhi(PFMetPt[0], PFMetPhi[0], Photon_pt[i], Photon_phi[i],seedT);
            */
            float METphi=PFMetPhi[0];
            float MET=PFMetPt[0];            
            
            if(MET>80.){
                float m_T = massT(Photon_phi[i],Photon_et[i],METphi,MET);
                h_massT->Fill(m_T);
                isCandidate++;
                h_sigmaR->Fill( sigmaRCosmic(Photon_sc_x[i],  Photon_sc_y[i], Photon_sc_z[i],CosmicMuon_OuterTrack_InnerPoint_x,CosmicMuon_OuterTrack_InnerPoint_y,CosmicMuon_OuterTrack_InnerPoint_z, CosmicMuon_OuterTrack_InnerPoint_px,CosmicMuon_OuterTrack_InnerPoint_py,CosmicMuon_OuterTrack_InnerPoint_pz, CosmicMuon_n) );
                h_deltaD->Fill( deltaRCosmic(Photon_sc_x[i],  Photon_sc_y[i], Photon_sc_z[i],CosmicMuon_OuterTrack_InnerPoint_x,CosmicMuon_OuterTrack_InnerPoint_y,CosmicMuon_OuterTrack_InnerPoint_z, CosmicMuon_OuterTrack_InnerPoint_px,CosmicMuon_OuterTrack_InnerPoint_py,CosmicMuon_OuterTrack_InnerPoint_pz, CosmicMuon_n) );
                h_pho_E->Fill(Photon_E[i]);
                h_pho_pt->Fill(Photon_pt[i]);
                h_pho_phi->Fill(Photon_phi[i]);
                h_pho_eta->Fill(Photon_eta[i]);
                h_pho_sieie->Fill(Photon_SigmaIetaIeta[i]);
                h_pho_sipip->Fill(Photon_SigmaIphiIphi[i]);
                h_pho_tSeed->Fill(Photon_timing_xtal[i][0]);//crystals are energy sorted :D
                h_e2e9->Fill(GetE2overE9(Photon_ncrys[i], Photon_ieta_xtalEB[i], Photon_iphi_xtalEB[i], Photon_energy_xtal[i]));
                h_roundness->Fill(Photon_Roundness[i]);
                
                float iLICTD = computeLICTD(Photon_ncrys[i], Photon_timing_xtal[i], Photon_energy_xtal[i]);
                h_LITCD->Fill(iLICTD);
                
                h_seed_iphi_ieta->Fill(Photon_iphi_xtalEB[i][0],Photon_ieta_xtalEB[i][0]);
                h_trkIso->Fill(Photon_trkSumPtHollowConeDR04[i]);
                h_ecalIso->Fill(Photon_ecalRecHitSumEtConeDR04[i]);
                h_hcalIso->Fill(Photon_hcalTowerSumEtConeDR04[i]);
                h_HoE->Fill(Photon_HoE[i]);
                            
                if(Photon_pt[i]>leadPhoPt){
                    leadPhoPt=Photon_pt[i];
                    leadPhoIdx=i;
                    leadMETphi=METphi;
                    leadMET = MET;
                }
                
            }//MetPt req.
            
        }//tight barrel phos
        h_cands->Fill(isCandidate);
        if(isCandidate>0){
            h_pho_n->Fill(Photon_n);
            h_sumEt->Fill(PFMetSumEt[0]);
        }
        if(leadPhoIdx>-1){
            h_tcmet_pt->Fill(leadMET);
            h_tcmet_phi->Fill(leadMETphi);
            float dPhi_PhoMet = deltaR(0,Photon_pt[leadPhoIdx],0,leadMETphi);
            h_dPhi_PhoMet->Fill(dPhi_PhoMet);
        }

        
        if(isCandidate>0){
            
            for(int pt=0;pt<500;pt++){
                int nJetsAboveCut=0;
                int nJetsAboveCut_eta25=0;
                int nTrksAboveCut=0;
                int nTrksAboveCut_eta20=0;
                int nTrksAboveCut_eta25=0;
                int nPhosAboveCut=0;
                
                for(int jet=0; jet<pfJet_n && jet<200;jet++){
                    float dR=deltaR(pfJet_eta[jet],pfJet_phi[jet],Photon_eta[leadPhoIdx],Photon_phi[leadPhoIdx]);
                    if(dR<0.3) continue; //ah! this jet lines up with the photon!
                    if(pfJet_pt[jet]>=pt) nJetsAboveCut++;
                    if(fabs(pfJet_eta[jet])<2.5 && pfJet_pt[jet]>=pt) nJetsAboveCut_eta25++;
                }//jets
                h_nJet_pT->Fill(0.5+pt,nJetsAboveCut);
                h_nJet_pT_eta25->Fill(0.5+pt,nJetsAboveCut_eta25);
                
                for(int trk=0; trk<Track_n && trk<400;trk++){
                    float dR=deltaR(Track_eta[trk],Track_phi[trk],Photon_eta[leadPhoIdx],Photon_phi[leadPhoIdx]);
                    if(dR<0.3) continue; //ah! this track lines up with the photon!
                    if(Track_pt[trk]>=pt) nTrksAboveCut++;
                    if(fabs(Track_eta[trk])<2.0 && Track_pt[trk]>pt) nTrksAboveCut_eta20++;
                    if(fabs(Track_eta[trk])<2.5 && Track_pt[trk]>pt) nTrksAboveCut_eta25++;
                }//jets
                h_nTrk_pT->Fill(0.5+pt,nTrksAboveCut);
                h_nTrk_pT_eta20->Fill(0.5+pt,nTrksAboveCut_eta20);
                h_nTrk_pT_eta25->Fill(0.5+pt,nTrksAboveCut_eta25);
                
                for(unsigned int pho=0; pho<indexTightBarrelPho.size();pho++){
                    int i = indexTightBarrelPho[pho];
                    if(i==leadPhoIdx) continue; //ah! this photon is a candidate object!
                    if(Photon_hasPixelSeed[i]==true &&Photon_pt[i]>=pt) nPhosAboveCut++;
                }//phos
                h_recoPhoPix_pT_barrel->Fill(0.5+pt,nPhosAboveCut);
                
            }//pt cut
            
            float leadJetPt=-1;
            int leadJetIdx=-1;
            for(int jet=0; jet<pfJet_n && jet<200;jet++){
                float d = sqrt( (Photon_vx[leadPhoIdx]-pfJet_vz[jet])*(Photon_vx[leadPhoIdx]-pfJet_vz[jet])
                    + (Photon_vy[leadPhoIdx]-pfJet_vy[jet])*(Photon_vy[leadPhoIdx]-pfJet_vy[jet])
                    + (Photon_vz[leadPhoIdx]-pfJet_vz[jet])*(Photon_vz[leadPhoIdx]-pfJet_vz[jet])
                    );
                h_vertexDistance_pho_jet_noR->Fill(d);
                float dR=deltaR(pfJet_eta[jet],pfJet_phi[jet],Photon_eta[leadPhoIdx],Photon_phi[leadPhoIdx]);
                if(dR<0.3) continue; //ah! this jet lines up with the photon!
                h_vertexDistance_pho_jet->Fill(d);
                if(pfJet_pt[jet]>leadJetPt){
                    leadJetPt=pfJet_pt[jet];
                    leadJetIdx = jet;
                }
            }
            if(leadJetIdx>-1 && leadPhoIdx>-1){
                float dR_lJetlPho=deltaR(pfJet_eta[leadJetIdx],pfJet_phi[leadJetIdx],Photon_eta[leadPhoIdx],Photon_phi[leadPhoIdx]);
                h_leadJet_leadPho_deltaR->Fill(dR_lJetlPho);
                    
                float dPhi_lJetMet=deltaR(0,pfJet_phi[leadJetIdx],0,leadMETphi);
                h_dPhi_lJetMet->Fill(dPhi_lJetMet);
                
                h_lJet_pt->Fill(pfJet_pt[leadJetIdx]);
            }
            
            int muG=0;
            int muT=0;
            int muSA=0;
            
            for(int m=0;m<Muon_n && m<200;m++){
                if(Muon_isGlobalMuon[m]==true) muG++;
                if(Muon_isTrackerMuon[m]==true) muT++;
                if(Muon_isStandAloneMuon[m]==true) muSA++;
                if(!Muon_isStandAloneMuon[m]) continue;
                float d = sqrt( (Photon_vx[leadPhoIdx]-Muon_vz[m])*(Photon_vx[leadPhoIdx]-Muon_vz[m])
                    + (Photon_vy[leadPhoIdx]-Muon_vy[m])*(Photon_vy[leadPhoIdx]-Muon_vy[m])
                    + (Photon_vz[leadPhoIdx]-Muon_vz[m])*(Photon_vz[leadPhoIdx]-Muon_vz[m])
                    );
                h_vertexDistance_pho_mu_noR->Fill(d);
                float dR=deltaR(Muon_eta[m],Muon_phi[m],Photon_eta[leadPhoIdx],Photon_phi[leadPhoIdx]);
                if(dR<0.3) continue; //ah! this mu lines up with the photon!
                h_vertexDistance_pho_mu->Fill(d);

            }//muons
            h_mu_isGlobal->Fill(muG);
            h_mu_isTracker->Fill(muT);
            h_mu_isSA->Fill(muSA);
            
            int cosmuG=0;
            int cosmuT=0;
            int cosmuSA=0;
            
            for(int m=0;m<CosmicMuon_n && m<200;m++){
                if(CosmicMuon_isGlobalMuon[m]==true) cosmuG++;
                if(CosmicMuon_isTrackerMuon[m]==true) cosmuT++;
                if(CosmicMuon_isStandAloneMuon[m]==true) cosmuSA++;
            }//cosmuons
            h_cosmu_isGlobal->Fill(cosmuG);
            h_cosmu_isTracker->Fill(cosmuT);
            h_cosmu_isSA->Fill(cosmuSA);
            
            for(int w=0;w<Photon_n && w<200; w++){
                if(w==leadPhoIdx) continue;
                float d = sqrt( (Photon_vx[leadPhoIdx]-Photon_vz[w])*(Photon_vx[leadPhoIdx]-Photon_vz[w])
                    + (Photon_vy[leadPhoIdx]-Photon_vy[w])*(Photon_vy[leadPhoIdx]-Photon_vy[w])
                    + (Photon_vz[leadPhoIdx]-Photon_vz[w])*(Photon_vz[leadPhoIdx]-Photon_vz[w])
                    );
                h_vertexDistance_pho_pho_noR->Fill(d);
                float dR=deltaR(Photon_eta[w],Photon_phi[w],Photon_eta[leadPhoIdx],Photon_phi[leadPhoIdx]);
                if(dR<0.3) continue; //ah! this mu lines up with the photon!
                h_vertexDistance_pho_pho->Fill(d);
            }//phos
            
            for(int trk=0; trk<Track_n && trk<400;trk++){
                float d = sqrt( (Photon_vx[leadPhoIdx]-Track_vz[trk])*(Photon_vx[leadPhoIdx]-Track_vz[trk])
                    + (Photon_vy[leadPhoIdx]-Track_vy[trk])*(Photon_vy[leadPhoIdx]-Track_vy[trk])
                    + (Photon_vz[leadPhoIdx]-Track_vz[trk])*(Photon_vz[leadPhoIdx]-Track_vz[trk])
                    );
                h_vertexDistance_pho_trk_noR->Fill(d);
                float dR=deltaR(Track_eta[trk],Track_phi[trk],Photon_eta[leadPhoIdx],Photon_phi[leadPhoIdx]);
                if(dR<0.3) continue; //ah! this track lines up with the photon!
                h_vertexDistance_pho_trk->Fill(d);
            }//trks
            
            n_ifVeto_none->Fill(0);
            bool hasOtherTightPho = hasOtherPho(0,leadPhoIdx,Photon_n,Photon_et,Photon_HoE,Photon_trkSumPtHollowConeDR04,Photon_ecalRecHitSumEtConeDR04,Photon_hcalTowerSumEtConeDR04,Photon_SigmaIetaIeta,Photon_sc_eta,Photon_vx,Photon_vy,Photon_vz);
            bool hasOtherLoosePho = hasOtherPho(1,leadPhoIdx,Photon_n,Photon_et,Photon_HoE,Photon_trkSumPtHollowConeDR04,Photon_ecalRecHitSumEtConeDR04,Photon_hcalTowerSumEtConeDR04,Photon_SigmaIetaIeta,Photon_sc_eta,Photon_vx,Photon_vy,Photon_vz);
            bool hasOtherEM = hasOtherPho(2,leadPhoIdx,Photon_n,Photon_et,Photon_HoE,Photon_trkSumPtHollowConeDR04,Photon_ecalRecHitSumEtConeDR04,Photon_hcalTowerSumEtConeDR04,Photon_SigmaIetaIeta,Photon_sc_eta,Photon_vx,Photon_vy,Photon_vz);
            bool hasOtherEM_20 = hasOtherPho(2,leadPhoIdx,Photon_n,Photon_et,Photon_HoE,Photon_trkSumPtHollowConeDR04,Photon_ecalRecHitSumEtConeDR04,Photon_hcalTowerSumEtConeDR04,Photon_SigmaIetaIeta,Photon_sc_eta,Photon_vx,Photon_vy,Photon_vz,20.);
            bool hasOtherEM_30 = hasOtherPho(2,leadPhoIdx,Photon_n,Photon_et,Photon_HoE,Photon_trkSumPtHollowConeDR04,Photon_ecalRecHitSumEtConeDR04,Photon_hcalTowerSumEtConeDR04,Photon_SigmaIetaIeta,Photon_sc_eta,Photon_vx,Photon_vy,Photon_vz,30.);
            bool hasOtherTightPho_barrel = hasOtherPho(0,leadPhoIdx,Photon_n,Photon_et,Photon_HoE,Photon_trkSumPtHollowConeDR04,Photon_ecalRecHitSumEtConeDR04,Photon_hcalTowerSumEtConeDR04,Photon_SigmaIetaIeta,Photon_sc_eta,Photon_vx,Photon_vy,Photon_vz,true,1.4442);
            bool hasOtherLoosePho_barrel = hasOtherPho(1,leadPhoIdx,Photon_n,Photon_et,Photon_HoE,Photon_trkSumPtHollowConeDR04,Photon_ecalRecHitSumEtConeDR04,Photon_hcalTowerSumEtConeDR04,Photon_SigmaIetaIeta,Photon_sc_eta,Photon_vx,Photon_vy,Photon_vz,true,1.4442);
            bool hasOtherEM_barrel = hasOtherPho(2,leadPhoIdx,Photon_n,Photon_et,Photon_HoE,Photon_trkSumPtHollowConeDR04,Photon_ecalRecHitSumEtConeDR04,Photon_hcalTowerSumEtConeDR04,Photon_SigmaIetaIeta,Photon_sc_eta,Photon_vx,Photon_vy,Photon_vz,true,1.4442);
            bool hasOtherEM_20_barrel = hasOtherPho(2,leadPhoIdx,Photon_n,Photon_et,Photon_HoE,Photon_trkSumPtHollowConeDR04,Photon_ecalRecHitSumEtConeDR04,Photon_hcalTowerSumEtConeDR04,Photon_SigmaIetaIeta,Photon_sc_eta,Photon_vx,Photon_vy,Photon_vz,20.,true,1.4442);
            bool hasOtherEM_30_barrel = hasOtherPho(2,leadPhoIdx,Photon_n,Photon_et,Photon_HoE,Photon_trkSumPtHollowConeDR04,Photon_ecalRecHitSumEtConeDR04,Photon_hcalTowerSumEtConeDR04,Photon_SigmaIetaIeta,Photon_sc_eta,Photon_vx,Photon_vy,Photon_vz,30.,true,1.4442);
            
            if(!hasOtherTightPho) n_ifVeto_tpho->Fill(0);
            if(!hasOtherLoosePho) n_ifVeto_lpho->Fill(0);
            if(!hasOtherEM) n_ifVeto_em->Fill(0);
            if(!hasOtherEM_20) n_ifVeto_em20->Fill(0);
            if(!hasOtherEM_30) n_ifVeto_em30->Fill(0);
            if(!hasOtherTightPho_barrel) n_ifVeto_tpho_barrel->Fill(0);
            if(!hasOtherLoosePho_barrel) n_ifVeto_lpho_barrel->Fill(0);
            if(!hasOtherEM_barrel) n_ifVeto_em_barrel->Fill(0);
            if(!hasOtherEM_20_barrel) n_ifVeto_em20_barrel->Fill(0);
            if(!hasOtherEM_30_barrel) n_ifVeto_em30_barrel->Fill(0);
            
            bool hasCosMuon = eventHasMu(CosmicMuon_n, CosmicMuon_OuterTrack_isNonnull, CosmicMuon_isStandAloneMuon);
            bool hasMuon = eventHasMu(Muon_n, Muon_OuterTrack_isNonnull, Muon_isStandAloneMuon);
            bool hasMuonWSameVertex = phoHasMu(Muon_n, Muon_OuterTrack_isNonnull, Muon_isStandAloneMuon, Muon_pt,
             Muon_eta, Muon_phi, Muon_vx, Muon_vy, Muon_vz,Photon_phi[leadPhoIdx], Photon_eta[leadPhoIdx],
             Photon_vx[leadPhoIdx],Photon_vy[leadPhoIdx],Photon_vz[leadPhoIdx],-10);
            
            if(!hasCosMuon) n_ifVeto_cosmu->Fill(0);
            if(!hasMuon) n_ifVeto_mu->Fill(0);
            if(!(hasCosMuon||hasMuon)) n_ifVeto_cosmu_or_mu->Fill(0);
            if(!hasMuonWSameVertex) n_ifVeto_muWithVertex->Fill(0);
            
            bool hasTrack5 = hasObj(Photon_phi[leadPhoIdx], Photon_eta[leadPhoIdx],Photon_vx[leadPhoIdx],Photon_vy[leadPhoIdx],Photon_vz[leadPhoIdx],Track_n,Track_pt,Track_eta,Track_phi,Track_vx,Track_vy,Track_vz,5.);
            bool hasTrack10 = hasObj(Photon_phi[leadPhoIdx], Photon_eta[leadPhoIdx],Photon_vx[leadPhoIdx],Photon_vy[leadPhoIdx],Photon_vz[leadPhoIdx],Track_n,Track_pt,Track_eta,Track_phi,Track_vx,Track_vy,Track_vz,10.);
            bool hasTrack20 = hasObj(Photon_phi[leadPhoIdx], Photon_eta[leadPhoIdx],Photon_vx[leadPhoIdx],Photon_vy[leadPhoIdx],Photon_vz[leadPhoIdx],Track_n,Track_pt,Track_eta,Track_phi,Track_vx,Track_vy,Track_vz,20.);
            bool hasTrack5_eta20 = hasObj(Photon_phi[leadPhoIdx], Photon_eta[leadPhoIdx],Photon_vx[leadPhoIdx],Photon_vy[leadPhoIdx],Photon_vz[leadPhoIdx],Track_n,Track_pt,Track_eta,Track_phi,Track_vx,Track_vy,Track_vz,5.,true,2.0);
            bool hasTrack10_eta20 = hasObj(Photon_phi[leadPhoIdx], Photon_eta[leadPhoIdx],Photon_vx[leadPhoIdx],Photon_vy[leadPhoIdx],Photon_vz[leadPhoIdx],Track_n,Track_pt,Track_eta,Track_phi,Track_vx,Track_vy,Track_vz,10.,true,2.0);
            bool hasTrack20_eta20 = hasObj(Photon_phi[leadPhoIdx], Photon_eta[leadPhoIdx],Photon_vx[leadPhoIdx],Photon_vy[leadPhoIdx],Photon_vz[leadPhoIdx],Track_n,Track_pt,Track_eta,Track_phi,Track_vx,Track_vy,Track_vz,20.,true,2.0);
            bool hasTrack5_eta25 = hasObj(Photon_phi[leadPhoIdx], Photon_eta[leadPhoIdx],Photon_vx[leadPhoIdx],Photon_vy[leadPhoIdx],Photon_vz[leadPhoIdx],Track_n,Track_pt,Track_eta,Track_phi,Track_vx,Track_vy,Track_vz,5.,true,2.5);
            bool hasTrack10_eta25 = hasObj(Photon_phi[leadPhoIdx], Photon_eta[leadPhoIdx],Photon_vx[leadPhoIdx],Photon_vy[leadPhoIdx],Photon_vz[leadPhoIdx],Track_n,Track_pt,Track_eta,Track_phi,Track_vx,Track_vy,Track_vz,10.,true,2.5);
            bool hasTrack20_eta25 = hasObj(Photon_phi[leadPhoIdx], Photon_eta[leadPhoIdx],Photon_vx[leadPhoIdx],Photon_vy[leadPhoIdx],Photon_vz[leadPhoIdx],Track_n,Track_pt,Track_eta,Track_phi,Track_vx,Track_vy,Track_vz,20.,true,2.5);
            
            if(!hasTrack5) n_ifVeto_trk5->Fill(0);
            if(!hasTrack20) n_ifVeto_trk20->Fill(0);
            if(!hasTrack10) n_ifVeto_trk10->Fill(0);
            if(!hasTrack5_eta20) n_ifVeto_trk5_eta20->Fill(0);
            if(!hasTrack10_eta20) n_ifVeto_trk10_eta20->Fill(0);
            if(!hasTrack20_eta20) n_ifVeto_trk20_eta20->Fill(0);
            if(!hasTrack5_eta25) n_ifVeto_trk5_eta25->Fill(0);
            if(!hasTrack10_eta25) n_ifVeto_trk10_eta25->Fill(0);
            if(!hasTrack20_eta25) n_ifVeto_trk20_eta25->Fill(0);    
            
            bool hasJet50 = hasObj(Photon_phi[leadPhoIdx], Photon_eta[leadPhoIdx],Photon_vx[leadPhoIdx],Photon_vy[leadPhoIdx],Photon_vz[leadPhoIdx],pfJet_n,pfJet_pt,pfJet_eta,pfJet_phi,pfJet_vx,pfJet_vy,pfJet_vz,50.);
            bool hasJet50_eta25 = hasObj(Photon_phi[leadPhoIdx], Photon_eta[leadPhoIdx],Photon_vx[leadPhoIdx],Photon_vy[leadPhoIdx],Photon_vz[leadPhoIdx],pfJet_n,pfJet_pt,pfJet_eta,pfJet_phi,pfJet_vx,pfJet_vy,pfJet_vz,50.,true,2.5);
            bool hasJet30 = hasObj(Photon_phi[leadPhoIdx], Photon_eta[leadPhoIdx],Photon_vx[leadPhoIdx],Photon_vy[leadPhoIdx],Photon_vz[leadPhoIdx],pfJet_n,pfJet_pt,pfJet_eta,pfJet_phi,pfJet_vx,pfJet_vy,pfJet_vz,30.);
            bool hasJet30_eta25 = hasObj(Photon_phi[leadPhoIdx], Photon_eta[leadPhoIdx],Photon_vx[leadPhoIdx],Photon_vy[leadPhoIdx],Photon_vz[leadPhoIdx],pfJet_n,pfJet_pt,pfJet_eta,pfJet_phi,pfJet_vx,pfJet_vy,pfJet_vz,30.,true,2.5);
            
            
            if(!hasJet50) n_ifVeto_jet50->Fill(0);
            if(!hasJet50_eta25) n_ifVeto_jet50_eta25->Fill(0);
            if(!hasJet30) n_ifVeto_jet30->Fill(0);
            if(!hasJet30_eta25) n_ifVeto_jet30_eta25->Fill(0);
             
                   
            if( !(hasTrack20 || hasJet50) ) n_ifVeto_trk20_or_jet50->Fill(0);
            if( !(hasTrack20_eta20 || hasJet50_eta25) ) n_ifVeto_trk20_eta20_or_jet50_eta25->Fill(0);
            if( !(hasTrack20_eta25 || hasJet50_eta25) ) n_ifVeto_trk20_eta25_or_jet50_eta25->Fill(0);
            
            if( !(hasTrack20 || hasJet30) ) n_ifVeto_trk20_or_jet30->Fill(0);
            if( !(hasTrack10 || hasJet30) ) n_ifVeto_trk10_or_jet30->Fill(0);
            if( !(hasTrack20_eta25 || hasJet30_eta25) ) n_ifVeto_trk20_eta25_or_jet30_eta25->Fill(0);
            if( !(hasTrack10_eta25 || hasJet30_eta25) ) n_ifVeto_trk10_eta25_or_jet30_eta25->Fill(0);
            
            if(! (hasOtherEM_20 || hasMuon) ) n_ifVeto_em20_or_mu->Fill(0);
            if(! (hasOtherEM_20 || hasCosMuon) ) n_ifVeto_em20_or_cosmu->Fill(0);
            if(! (hasOtherEM_20 || hasMuon || hasCosMuon) ) n_ifVeto_em20_or_mu_or_cosmu->Fill(0);
            
            if(! (hasOtherEM_20 || hasTrack20)) n_ifVeto_em20_or_trk20->Fill(0);
            if(! (hasOtherEM_20 || hasJet50)) n_ifVeto_em20_or_jet50->Fill(0);
            if(! (hasOtherEM_20 || hasTrack20 || hasJet50)) n_ifVeto_em20_or_trk20_or_jet50->Fill(0);
            
            if(! (hasOtherEM_20 || hasTrack20 || hasJet50 || hasMuon || hasCosMuon))
                n_ifVeto_em20_or_trk20_or_jet50_or_mu_or_cosmu->Fill(0);
            if(! (hasOtherEM_20_barrel || hasTrack20_eta25 || hasJet50_eta25 || hasMuon || hasCosMuon))
                n_ifVeto_em20_barrel_or_trk20_eta25_or_jet50_eta25_or_mu_or_cosmu->Fill(0);
            
            if(! (hasOtherEM || hasTrack10 || hasJet30))
                n_ifVeto_em_or_trk10_or_jet30->Fill(0);
            if(! (hasOtherEM || hasTrack10 || hasJet30 || hasMuon || hasCosMuon))
                n_ifVeto_em_or_trk10_or_jet30_or_mu_or_cosmu->Fill(0);
            if(! (hasOtherEM || hasTrack10_eta25 || hasJet30_eta25))
                n_ifVeto_em_or_trk10_eta25_or_jet30_eta25->Fill(0);
            if(! (hasOtherEM || hasTrack10_eta25 || hasJet30_eta25 || hasMuon || hasCosMuon))
                n_ifVeto_em_or_trk10_eta25_or_jet30_eta25_or_mu_or_cosmu->Fill(0);
                
            
            
        }//if cands>0

    }//loop over entries jentry
    cout<<endl;
    cout<<"check for zombies"<<endl;
    if( outputFile->IsZombie() ){
        cout << "Error opening "<<"output file !" << endl;
        exit(1);
    }else{
        cout<<"about to write histos"<<endl;
        outputFile->Write();
        outputFile->Close();
    }
}
