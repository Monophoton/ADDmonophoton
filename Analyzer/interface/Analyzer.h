#include "TFile.h"
#include "TBrowser.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TTree.h"
#include "ADDmonophoton/Analyzer/interface/CrystalInfo.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
//#include "FWCore/Framework/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include <string>

class Analyzer : public edm::EDAnalyzer {
 public:
  explicit Analyzer(const edm::ParameterSet&);
  ~Analyzer();
  
  
 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  
  // ----------member data ---------------------------
  edm::ESHandle<CaloTopology> theCaloTopo_;
  int nevents;
  bool is_signal_event, is_Z_event, is_W_event;
  bool is_Znunu_event, is_Zelec_event, is_Zmu_event, is_Ztau_event ;
  bool is_Welec_event, is_Wmu_event, is_Wtau_event ;
  bool is_SingleHardPhoton_event;
  bool is_diphoton_event;
  bool is_isr_photon_event;
  float gen_pthat;
  
  int n_signal_events,n_Z_events,n_W_events;
  int n_Zelec_events, n_Zmu_events, n_Ztau_events, n_Znunu_events; 
  int n_Welec_events, n_Wmu_events, n_Wtau_events;
  int n_diphoton_events, n_SingleHardPhoton_events;
  
  unsigned int RunNumber, EventNumber, LumiNumber, BXNumber;
  unsigned int totalIntensityBeam1, totalIntensityBeam2;
  float avgInsDelLumi, avgInsDelLumiErr, avgInsRecLumi, avgInsRecLumiErr;
  int ngenphotons;
  int nhardphotons;
  int Photon_n;
  int Vertex_n;
  int Muon_n;
  int CosmicMuon_n;
  int Tau_n;
  int Electron_n;
  int Track_n;
  int Jet_n;
  int HERecHit_subset_n;
  //HLT
  
  edm::TriggerNames triggerNames_;  // TriggerNames class
  std::vector<std::string>  hlNames_;  // name of each HLT algorithm
  
  bool HLT_MET50_event;
  bool HLT_MET75_event;
  bool HLT_Photon15_event;
  bool HLT_Photon25_event;
  bool HLT_DoubleEle10_event;
  bool HLT_DoubleMu3_event;
  bool HLT_Photon20_event;
  bool HLT_Photon20_Cleaned_event;
  bool HLT_Photon30_event;
  bool HLT_Photon30_L1R_8E29_event;
  bool HLT_Photon30_L1R_1E31_event;
  bool HLT_Photon30_Cleaned_event;
  bool HLT_Photon30_Isol_EBOnly_Cleaned_event;
  bool HLT_Photon35_Isol_Cleaned_event;
  bool HLT_Photon50_Cleaned_event;
  bool HLT_Photon70_NoHE_Cleaned_event;
  bool HLT_Photon100_NoHE_Cleaned_event;
  bool HLT_DoublePhoton17_L1R_event;
  bool HLT_DoublePhoton5_CEP_L1R_event;
  bool HLT_Photon100_NoHE_Cleaned_L1R_v1_event;
  bool HLT_Photon10_Cleaned_L1R_event;
  bool HLT_Photon15_Cleaned_L1R_event;
  bool HLT_Photon17_SC17HE_L1R_v1_event;
  bool HLT_Photon20_NoHE_L1R_event;
  bool HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1_event;
  bool HLT_Photon35_Isol_Cleaned_L1R_v1_event;
  bool HLT_Photon50_Cleaned_L1R_v1_event;
  bool HLT_Photon50_NoHE_L1R_event; 
  bool HLT_Photon70_NoHE_Cleaned_L1R_v1_event;
 
  std::map<std::string,int> HLT_chosen;
  std::map<std::string,int> L1_chosen;
  
  std::vector<std::string> JET_CORR;
  
  edm::InputTag eleLabel_;
  edm::InputTag muoLabel_;
  edm::InputTag jetLabel_;
  edm::InputTag tauLabel_;
  edm::InputTag metLabel_;
  edm::InputTag PFmetLabel_;
  edm::InputTag TCmetLabel_;
  edm::InputTag phoLabel_;
  edm::InputTag rechitBLabel_;
  edm::InputTag rechitELabel_;
  edm::InputTag hlTriggerResults_;  // Input tag for TriggerResults
  edm::InputTag Tracks_;
  edm::InputTag Vertices_;
  std::string outFile_;
  bool rungenParticleCandidates_;
  bool runphotons_;
  bool runmet_;
  bool rungenmet_;
  bool runPFmet_;
  bool runTCmet_;
  bool runelectrons_;
  bool runmuons_;
  bool runcosmicmuons_;
  bool runjets_;
  bool runtaus_;
  bool runHLT_;
  bool runL1_;
  bool runscraping_;
  bool runtracks_;
  bool runrechit_;
  bool runHErechit_;
  bool runvertex_;
  bool init_;
  
  //output root file
  TFile *f;
  
  // output root tree
  TTree* myEvent;
  //variables to be filled in the tree
  
  //vertex variables
  float vx[10];
  float vy[10];
  float vz[10];
  float chi2[10];
  float vtracksize[10];
  float vndof[10];
  bool  v_isFake[10];
  float v_d0[10];
  
  //scraping variables
  bool  Scraping_isScrapingEvent;
  int   Scraping_numOfTracks;
  float Scraping_fractionOfGoodTracks;
  
  //track variables
  float trk_pt[200];
  float trk_px[200];
  float trk_py[200];
  float trk_pz[200];
  float trk_eta[200];
  float trk_phi[200];
  
  //jet variables
  float jet_pt[100];
  float jet_px[100];
  float jet_py[100];
  float jet_pz[100];
  float jet_eta[100];
  float jet_phi[100];
  float jet_emEnergyFraction[100];
  float jet_energyFractionHadronic[100];
  int   jet_hitsInN90[100];
  int   jet_n90Hits[100];
  float jet_fHPD[100];
  float jet_fRBX[100];
  float jet_RHF[100];
  int   jet_nTowers[100];
  
  //electron variables
  float electron_pt[100];
  float electron_px[100];
  float electron_py[100];
  float electron_pz[100];
  float electron_energy[100];
  float electron_charge[100];
  float electron_eta[100];
  float electron_phi[100];
  float electron_trkIso[100];
  
  //muon variables
  float muon_pt[100];
  float muon_px[100];
  float muon_py[100];
  float muon_pz[100];
  float muon_energy[100];
  float muon_charge[100];
  float muon_eta[100];
  float muon_phi[100];
  bool  muon_isGlobalMuon[100];
  bool  muon_isTrackerMuon[100];
  bool  muon_isStandAloneMuon[100];
  bool  muon_InnerTrack_isNonnull[100];
  bool  muon_OuterTrack_isNonnull[100];
  
  float muon_OuterTrack_InnerPoint_x[100];
  float muon_OuterTrack_InnerPoint_y[100];
  float muon_OuterTrack_InnerPoint_z[100];
  float muon_OuterTrack_InnerPoint_px[100];
  float muon_OuterTrack_InnerPoint_py[100];
  float muon_OuterTrack_InnerPoint_pz[100];
  float muon_OuterTrack_OuterPoint_x[100];
  float muon_OuterTrack_OuterPoint_y[100];
  float muon_OuterTrack_OuterPoint_z[100];
  float muon_OuterTrack_OuterPoint_px[100];
  float muon_OuterTrack_OuterPoint_py[100];
  float muon_OuterTrack_OuterPoint_pz[100];
  float muon_InnerTrack_InnerPoint_x[100];
  float muon_InnerTrack_InnerPoint_y[100];
  float muon_InnerTrack_InnerPoint_z[100];
  float muon_InnerTrack_InnerPoint_px[100];
  float muon_InnerTrack_InnerPoint_py[100];
  float muon_InnerTrack_InnerPoint_pz[100];
  float muon_InnerTrack_OuterPoint_x[100];
  float muon_InnerTrack_OuterPoint_y[100];
  float muon_InnerTrack_OuterPoint_z[100];
  float muon_InnerTrack_OuterPoint_px[100];
  float muon_InnerTrack_OuterPoint_py[100];
  float muon_InnerTrack_OuterPoint_pz[100];
  
  //cosmicmuon variables
  float cosmicmuon_pt[100];
  float cosmicmuon_px[100];
  float cosmicmuon_py[100];
  float cosmicmuon_pz[100];
  float cosmicmuon_energy[100];
  float cosmicmuon_charge[100];
  float cosmicmuon_eta[100];
  float cosmicmuon_phi[100];
  bool  cosmicmuon_isGlobalMuon[100];
  bool  cosmicmuon_isTrackerMuon[100];
  bool  cosmicmuon_isStandAloneMuon[100];
  bool  cosmicmuon_InnerTrack_isNonnull[100];
  bool  cosmicmuon_OuterTrack_isNonnull[100];
  
  float cosmicmuon_OuterTrack_InnerPoint_x[100];
  float cosmicmuon_OuterTrack_InnerPoint_y[100];
  float cosmicmuon_OuterTrack_InnerPoint_z[100];
  float cosmicmuon_OuterTrack_InnerPoint_px[100];
  float cosmicmuon_OuterTrack_InnerPoint_py[100];
  float cosmicmuon_OuterTrack_InnerPoint_pz[100];
  float cosmicmuon_OuterTrack_OuterPoint_x[100];
  float cosmicmuon_OuterTrack_OuterPoint_y[100];
  float cosmicmuon_OuterTrack_OuterPoint_z[100];
  float cosmicmuon_OuterTrack_OuterPoint_px[100];
  float cosmicmuon_OuterTrack_OuterPoint_py[100];
  float cosmicmuon_OuterTrack_OuterPoint_pz[100];
  float cosmicmuon_InnerTrack_InnerPoint_x[100];
  float cosmicmuon_InnerTrack_InnerPoint_y[100];
  float cosmicmuon_InnerTrack_InnerPoint_z[100];
  float cosmicmuon_InnerTrack_InnerPoint_px[100];
  float cosmicmuon_InnerTrack_InnerPoint_py[100];
  float cosmicmuon_InnerTrack_InnerPoint_pz[100];
  float cosmicmuon_InnerTrack_OuterPoint_x[100];
  float cosmicmuon_InnerTrack_OuterPoint_y[100];
  float cosmicmuon_InnerTrack_OuterPoint_z[100];
  float cosmicmuon_InnerTrack_OuterPoint_px[100];
  float cosmicmuon_InnerTrack_OuterPoint_py[100];
  float cosmicmuon_InnerTrack_OuterPoint_pz[100];
  
  //tau variables
  float tau_pt[100];
  float tau_px[100];
  float tau_py[100];
  float tau_pz[100];
  float tau_energy[100];
  float tau_charge[100];
  float tau_eta[100];
  float tau_phi[100];
  
  //gen level variables
  float gen_pho_pt[1000];
  float gen_pho_px[1000];
  float gen_pho_py[1000];
  float gen_pho_pz[1000];
  float gen_pho_phi[1000];
  float gen_pho_eta[1000];
  float gen_pho_E[1000];
  int   gen_pho_status[1000];
  int   gen_pho_motherID[1000];
  int   gen_pho_motherStatus[1000];
  float gen_pho_motherPt[1000];
  float gen_pho_motherEta[1000];
  float gen_pho_motherPhi[1000];
  int   gen_pho_GrandmotherID[1000];
  int   gen_pho_GrandmotherStatus[1000];
  float gen_pho_GrandmotherPt[1000];
  float gen_pho_GrandmotherEta[1000];
  float gen_pho_GrandmotherPhi[1000];
  
  float gen_Hpho_pt[2];
  float gen_Hpho_px[2];
  float gen_Hpho_py[2];
  float gen_Hpho_pz[2];
  float gen_Hpho_phi[2];
  float gen_Hpho_eta[2];
  float gen_Hpho_E[2];
  
  //gen level variables for Graviton
  float gen_graviton_pt;
  float gen_graviton_px;
  float gen_graviton_py;
  float gen_graviton_pz;
  float gen_graviton_phi;
  float gen_graviton_eta;
  float gen_graviton_E;
  
  //gen level variables for Zdaughter
  float gen_Zdaughter_pt[2];
  float gen_Zdaughter_px[2];
  float gen_Zdaughter_py[2];
  float gen_Zdaughter_pz[2];
  float gen_Zdaughter_phi[2];
  float gen_Zdaughter_eta[2];
  float gen_Zdaughter_E[2];
  int gen_Zdaughter_charge[2];
  int gen_Zdaughter_ID[2];
  
  //gen level variables for Z
  float gen_Zboson_pt;
  float gen_Zboson_px;
  float gen_Zboson_py;
  float gen_Zboson_pz;
  float gen_Zboson_phi;
  float gen_Zboson_eta;
  float gen_Zboson_E;
  
  //gen level variables for Wdaughter
  float gen_Wdaughter_pt[2];
  float gen_Wdaughter_px[2];
  float gen_Wdaughter_py[2];
  float gen_Wdaughter_pz[2];
  float gen_Wdaughter_phi[2];
  float gen_Wdaughter_eta[2];
  float gen_Wdaughter_E[2];
  int gen_Wdaughter_charge[2];
  int gen_Wdaughter_ID[2];
  
  //gen level variables for W
  float gen_Wboson_pt;
  float gen_Wboson_px;
  float gen_Wboson_py;
  float gen_Wboson_pz;
  float gen_Wboson_phi;
  float gen_Wboson_eta;
  float gen_Wboson_E;
  int gen_Wboson_charge;
  int gen_Wboson_ID;
  
  //gen level variables for  mu/tau daughter
  float gen_Muon_ID[3];
  float gen_Muon_Status[3];
  float gen_Muon_Pt[3];
  float gen_MuonDaughter_pt[3];
  float gen_MuonDaughter_px[3];
  float gen_MuonDaughter_py[3];
  float gen_MuonDaughter_pz[3];
  float gen_MuonDaughter_phi[3];
  float gen_MuonDaughter_eta[3];
  float gen_MuonDaughter_E[3];
  int gen_MuonDaughter_charge[3];
  int gen_MuonDaughter_status[3];
  int gen_MuonDaughter_ID[3];
  
  float gen_tau_ID[3];
  float gen_tau_Status[3];
  float gen_tau_Pt[3];
  float gen_tauDaughter_pt[3];
  float gen_tauDaughter_px[3];
  float gen_tauDaughter_py[3];
  float gen_tauDaughter_pz[3];
  float gen_tauDaughter_phi[3];
  float gen_tauDaughter_eta[3];
  float gen_tauDaughter_E[3];
  int gen_tauDaughter_charge[3];
  int gen_tauDaughter_status[3];
  int gen_tauDaughter_ID[3];
  
  float pho_E[100];
  float pho_pt[100];
  float pho_eta[100];
  float pho_phi[100];
  float pho_px[100];
  float pho_py[100];
  float pho_pz[100];
  float pho_r9[100];
  bool  pho_isEB[100];
  bool  pho_isEE[100];
  bool  pho_isEBGap[100];
  bool  pho_isEEGap[100];
  bool  pho_isEBEEGap[100];
  float pho_e1x5[100];
  float pho_e2x5[100];
  float pho_e3x3[100];
  float pho_e5x5[100];
  float pho_r1x5[100];
  float pho_r2x5[100];
  float pho_SigmaEtaEta[100];
  float pho_SigmaIetaIeta[100];
  float pho_roundness[100];
  float pho_angle[100];
  float pho_maxEnergyXtal[100];
  float pho_theta[100];
  float pho_et[100];
  float pho_swissCross[100];
  bool    pho_isConverted[100];
  
  //isolation variables
  float pho_ecalRecHitSumEtConeDR03[100];
  float pho_hcalTowerSumEtConeDR03[100];
  float pho_hcalDepth1TowerSumEtConeDR03[100];
  float pho_hcalDepth2TowerSumEtConeDR03[100];
  float pho_trkSumPtSolidConeDR03[100];
  float pho_trkSumPtHollowConeDR03[100];
  int   pho_nTrkSolidConeDR03[100];
  int   pho_nTrkHollowConeDR03[100];
  float pho_ecalRecHitSumEtConeDR04[100];
  float pho_hcalTowerSumEtConeDR04[100];
  float pho_hcalDepth1TowerSumEtConeDR04[100];
  float pho_hcalDepth2TowerSumEtConeDR04[100];
  float pho_trkSumPtSolidConeDR04[100];
  float pho_trkSumPtHollowConeDR04[100];
  int   pho_nTrkSolidConeDR04[100];
  int   pho_nTrkHollowConeDR04[100];
  float pho_HoE[100];
  bool  pho_hasPixelSeed[100];   
  
  //SC variables
  float pho_sc_energy[100];
  int   pho_size[100];
  float pho_sc_eta[100];
  float pho_sc_phi[100];
  float pho_sc_etaWidth[100];
  float pho_sc_phiWidth[100];
  float pho_sc_et[100];
  
  //gen matched photon 
  float matchpho_E[100];
  float matchpho_pt[100];
  float matchpho_eta[100];
  float matchpho_phi[100];
  float matchpho_px[100];
  float matchpho_py[100];
  float matchpho_pz[100];
  bool  ismatchedpho[100];
  
  //converted photon variabes
  unsigned int pho_nTracks[100];
  float pho_pairInvariantMass[100];
  float pho_pairCotThetaSeparation[100];
  float pho_pairMomentum_x[100];
  float pho_pairMomentum_y[100];
  float pho_pairMomentum_z[100];
  float pho_EoverP[100];
  float pho_vertex_x[100];
  float pho_vertex_y[100];
  float pho_vertex_z[100];
  float pho_zOfPrimaryVertex[100];
  float pho_distOfMinimumApproach[100];
  float pho_dPhiTracksAtVtx[100];      
  float pho_dPhiTracksAtEcal[100];     
  float pho_dEtaTracksAtEcal[100];     
  
  //rechit information
  int ncrysPhoton[100];
  float pho_timing_xtal[100][100];
  float pho_timingavg_xtal[100];
  float pho_energy_xtal[100][100];
  int   pho_ieta_xtalEB[100][100];
  int   pho_iphi_xtalEB[100][100];
  float pho_rookFraction[100];
  float pho_s9[100];
  
  //HErechit information
  unsigned int HERecHit_subset_detid[10000];
  float HERecHit_subset_energy[10000];
  float HERecHit_subset_time[10000];
  float HERecHit_subset_depth[10000];
  float HERecHit_subset_phi[10000];
  float HERecHit_subset_eta[10000];
  float HERecHit_subset_x[10000];
  float HERecHit_subset_y[10000];
  float HERecHit_subset_z[10000];
  
  //calomet variables
  float CaloMetSig;
  float CaloMetCorr;
  float CaloMetEt;
  float CaloMetEx;
  float CaloMetEy;
  float CaloMetPhi;
  float CaloEtFractionHadronic;
  float CaloEmEtFraction;
  float CaloHadEtInHB;
  float CaloHadEtInHO;
  float CaloHadEtInHE;
  float CaloHadEtInHF;
  float CaloEmEtInEB;
  float CaloEmEtInEE;
  float CaloEmEtInHF;
  float CaloMetEz;
  float CaloMaxEtInEmTowers;
  float CaloSumEt;
  float CaloMaxEtInHadTowers;
  float Delta_phi;
  
  //PFMet variables
  float PFMetPt;
  float PFMetPx;
  float PFMetPy;
  float PFMetPhi;
  float PFMetSumEt;
  float Delta_phiPF;
  
  //TCMet variables
  float TCMetPt;
  float TCMetPx;
  float TCMetPy;
  float TCMetPhi;
  float TCMetSumEt;
  float Delta_phiTC;
  
  //genMet variables
  float genMetPt;
  float genMetPx;
  float genMetPy;
  float genMetPhi;
  float genMetSumEt;
  float Delta_phiGEN;
  
};

