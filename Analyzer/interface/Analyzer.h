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
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
//#include "RecoLocalCalo/EcalRecAlgos/interface/EcalCleaningAlgo.h"
#include <string>
#include <map>



class Analyzer : public edm::EDAnalyzer {
 public:
  explicit Analyzer(const edm::ParameterSet&);
  ~Analyzer();
  
  
 private:
  virtual void beginJob() ;
  virtual void beginRun(const edm::Run& , const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

 double GetE2OverE9(const DetId id, const EcalRecHitCollection & recHits); 
 double  e4e1(const DetId& id, const EcalRecHitCollection & rhs); 
 double  Gete6e2(const DetId& id, const EcalRecHitCollection & rhs); 

  HLTConfigProvider hltConfig_;
  std::vector<std::string> photon_triggers_in_run;
  std::vector<std::string> met_triggers_in_run;
  std::vector<std::string> cosmic_triggers_in_run;
  std::vector<std::string> halo_triggers_in_run;
  std::vector<std::string> all_triggers;
  std::vector<int> all_triggerprescales;
  std::vector<bool> all_ifTriggerpassed;
  std::vector<std::string> *triggernames;
  std::vector<int> *triggerprescales;
  std::vector<bool> *ifTriggerpassed;
  int ntriggers;

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
  int ucPhoton_n;
  int CaloTower_n;
  int Vertex_n;
  int Muon_n;
  int CosmicMuon_n;
  int Tau_n;
  int Electron_n;
  int Track_n;
  int Jet_n;
  int pfJet_n;
  int HERecHit_subset_n;
  int ucHERecHit_subset_n;
  int CSCseg_n;
  int RPChit_n;
  int genJet_n;
  //HLT

 //To get e6e2
// EcalCleaningAlgo * ecalCleaningTool_;

 std::vector<std::string>  hlNames_;           
  //HLT      
  TString module_type[50];
  float trobjpt[100][100][10];
  float trobjeta[100][100][10];
  float trobjphi[100][100][10];


  //maximum size to be filled in tree branches
  size_t MaxN;

  edm::TriggerNames triggerNames_;  // TriggerNames class
 
  std::map<std::string,int> HLT_chosen;
  std::map<std::string,int> L1_chosen;
  
  std::vector<std::string> JET_CORR;
  
  edm::InputTag eleLabel_;
  edm::InputTag muoLabel_;
  edm::InputTag cosMuoLabel_;
  edm::InputTag jetLabel_;
  edm::InputTag pfjetLabel_;
  edm::InputTag genjetLabel_;
  edm::InputTag tauLabel_;
  edm::InputTag metLabel_;
  edm::InputTag PFmetLabel_;
  edm::InputTag TCmetLabel_;
  edm::InputTag phoLabel_;
  edm::InputTag ucphoLabel_;
  edm::InputTag caloTowerLabel_;
  edm::InputTag cscLabel_;
  edm::InputTag rpcLabel_;
  edm::InputTag rechitBLabel_;
  edm::InputTag rechitELabel_;
  edm::InputTag hcalrechitLabel_;
  edm::InputTag hlTriggerResults_;  // Input tag for TriggerResults
  edm::InputTag triggerEventTag_; 
  edm::InputTag Tracks_;
  edm::InputTag BeamHaloSummaryLabel_;
  edm::InputTag Vertices_;
  edm::InputTag pileupLabel_;
  edm::InputTag rhoLabel_;
  edm::InputTag sigmaLabel_;

  std::string outFile_;
  std::string hltlabel_;

  bool rungenParticleCandidates_;
  bool runphotons_;
  bool runucphotons_;
  bool runmet_;
  bool rungenmet_;
  bool runPFmet_;
  bool runTCmet_;
  bool runelectrons_;
  bool runmuons_;
  bool runcosmicmuons_;
  bool runjets_;
  bool runpfjets_;
  bool rungenjets_;
  bool runtaus_;
  bool runDetailTauInfo_;
  bool runHLT_;
  bool runL1_;
  bool runscraping_;
  bool runtracks_;
  bool runrechit_;
  bool runHErechit_;
  bool runvertex_;
  bool runCSCseg_;
  bool runBeamHaloSummary_;      
  bool runRPChit_;
  bool runPileUp_;
  bool runcaloTower_;
  bool isAOD_;
  bool debug_;
  bool init_;
  
  //output root file
  TFile *f;
  
  // output root tree
  TTree* myEvent;
  //variables to be filled in the tree
  
  //vertex variables
  float vx[200];
  float vy[200];
  float vz[200];
  float chi2[200];
  int vtracksize[200];
  int vndof[200];
  bool  v_isFake[200];
  float v_d0[200];
  
  //scraping variables
  bool  Scraping_isScrapingEvent;
  int   Scraping_numOfTracks;
  float Scraping_fractionOfGoodTracks;
  
  //track variables
  float trk_pt[400];
  float trk_px[400];
  float trk_py[400];
  float trk_pz[400];
  float trk_vx[400];
  float trk_vy[400];
  float trk_vz[400];
  float trk_eta[400];
  float trk_phi[400];
  
  //jet variables
  float jet_pt[200];
  float jet_px[200];
  float jet_py[200];
  float jet_E[200];
  float jet_pz[200];
  float jet_vx[200];
  float jet_vy[200];
  float jet_vz[200];
  float jet_eta[200];
  float jet_phi[200];
  float jet_emEnergyFraction[200];
  float jet_energyFractionHadronic[200];
  int   jet_hitsInN90[200];
  int   jet_n90Hits[200];
  float jet_fHPD[200];
  float jet_fRBX[200];
  float jet_RHF[200];
  int   jet_nTowers[200];
  float jet_jecUncer[200];
  float jet_jecCorr[200];

   //some uncorrected jet infor
  float ucjet_pt[200];
  float ucjet_px[200];
  float ucjet_py[200];
  float ucjet_pz[200];
  float ucjet_E[200];
  float ucjet_eta[200];
  float ucjet_phi[200];


  //pfjet variables
  float pfjet_pt[200];
  float pfjet_px[200];
  float pfjet_py[200];
  float pfjet_E[200];
  float pfjet_pz[200];
  float pfjet_vx[200];
  float pfjet_vy[200];
  float pfjet_vz[200];
  float pfjet_eta[200];
  float pfjet_phi[200];
  float pfjet_jecUncer[200];
  float pfjet_jecCorr[200];
  /*
  float pfjet_emEnergyFraction[100];
  float pfjet_energyFractionHadronic[100];
  int   pfjet_hitsInN90[100];
  int   pfjet_n90Hits[100];
  float pfjet_fHPD[100];
  float pfjet_fRBX[100];
  float pfjet_RHF[100];
  int   pfjet_nTowers[100];
 */
 
  //some uncorrected jet infor
  float ucpfjet_pt[200];
  float ucpfjet_px[200];
  float ucpfjet_py[200];
  float ucpfjet_pz[200];
  float ucpfjet_E[200];
  float ucpfjet_eta[200];
  float ucpfjet_phi[200];

  
  //ak5GenJets
  float genjet_pt[200];
  float genjet_px[200];
  float genjet_py[200];
  float genjet_pz[200];
  float genjet_E[200];
  float genjet_eta[200];
  float genjet_phi[200];
 


  //electron variables
  float electron_pt[200];
  float electron_px[200];
  float electron_py[200];
  float electron_pz[200];
  float electron_vx[200];
  float electron_vy[200];
  float electron_vz[200];
  float electron_energy[200];
  float electron_charge[200];
  float electron_eta[200];
  float electron_phi[200];
  float electron_trkIso[200];
  float electron_ecalIso[200];
  float electron_hcalIso[200];
  float electron_HoE[200];
  float electron_SigmaIetaIeta[200];
  float electron_dEtaIn[200];
  float electron_dPhiIn[200];
  float electron_sc_energy[200];
  float electron_sc_eta[200];
  float electron_sc_phi[200];
  
  //muon variables
  float muon_pt[200];
  float muon_px[200];
  float muon_py[200];
  float muon_pz[200];
  float muon_vx[200];
  float muon_vy[200];
  float muon_vz[200];
  float muon_energy[200];
  float muon_charge[200];
  float muon_eta[200];
  float muon_phi[200];
  bool  muon_isGlobalMuon[200];
  bool  muon_isTrackerMuon[200];
  bool  muon_isStandAloneMuon[200];
  bool  muon_InnerTrack_isNonnull[200];
  bool  muon_OuterTrack_isNonnull[200];
  
  float muon_OuterTrack_InnerPoint_x[200];
  float muon_OuterTrack_InnerPoint_y[200];
  float muon_OuterTrack_InnerPoint_z[200];
  float muon_OuterTrack_InnerPoint_px[200];
  float muon_OuterTrack_InnerPoint_py[200];
  float muon_OuterTrack_InnerPoint_pz[200];
  float muon_OuterTrack_OuterPoint_x[200];
  float muon_OuterTrack_OuterPoint_y[200];
  float muon_OuterTrack_OuterPoint_z[200];
  float muon_OuterTrack_OuterPoint_px[200];
  float muon_OuterTrack_OuterPoint_py[200];
  float muon_OuterTrack_OuterPoint_pz[200];
  float muon_InnerTrack_InnerPoint_x[200];
  float muon_InnerTrack_InnerPoint_y[200];
  float muon_InnerTrack_InnerPoint_z[200];
  float muon_InnerTrack_InnerPoint_px[200];
  float muon_InnerTrack_InnerPoint_py[200];
  float muon_InnerTrack_InnerPoint_pz[200];
  float muon_InnerTrack_OuterPoint_x[200];
  float muon_InnerTrack_OuterPoint_y[200];
  float muon_InnerTrack_OuterPoint_z[200];
  float muon_InnerTrack_OuterPoint_px[200];
  float muon_InnerTrack_OuterPoint_py[200];
  float muon_InnerTrack_OuterPoint_pz[200];
  
  //FOR AOD only
  float muon_OuterPoint_x[200];
  float muon_OuterPoint_y[200];
  float muon_OuterPoint_z[200];
  float muon_InnerPoint_x[200];
  float muon_InnerPoint_y[200];
  float muon_InnerPoint_z[200];

 
  //cosmicmuon variables
  float cosmicmuon_pt[200];
  float cosmicmuon_px[200];
  float cosmicmuon_py[200];
  float cosmicmuon_pz[200];
  float cosmicmuon_energy[200];
  float cosmicmuon_charge[200];
  float cosmicmuon_eta[200];
  float cosmicmuon_phi[200];
  bool  cosmicmuon_isGlobalMuon[200];
  bool  cosmicmuon_isTrackerMuon[200];
  bool  cosmicmuon_isStandAloneMuon[200];
  bool  cosmicmuon_InnerTrack_isNonnull[200];
  bool  cosmicmuon_OuterTrack_isNonnull[200];
  
  float cosmicmuon_OuterTrack_InnerPoint_x[200];
  float cosmicmuon_OuterTrack_InnerPoint_y[200];
  float cosmicmuon_OuterTrack_InnerPoint_z[200];
  float cosmicmuon_OuterTrack_InnerPoint_px[200];
  float cosmicmuon_OuterTrack_InnerPoint_py[200];
  float cosmicmuon_OuterTrack_InnerPoint_pz[200];
  float cosmicmuon_OuterTrack_OuterPoint_x[200];
  float cosmicmuon_OuterTrack_OuterPoint_y[200];
  float cosmicmuon_OuterTrack_OuterPoint_z[200];
  float cosmicmuon_OuterTrack_OuterPoint_px[200];
  float cosmicmuon_OuterTrack_OuterPoint_py[200];
  float cosmicmuon_OuterTrack_OuterPoint_pz[200];
  float cosmicmuon_InnerTrack_InnerPoint_x[200];
  float cosmicmuon_InnerTrack_InnerPoint_y[200];
  float cosmicmuon_InnerTrack_InnerPoint_z[200];
  float cosmicmuon_InnerTrack_InnerPoint_px[200];
  float cosmicmuon_InnerTrack_InnerPoint_py[200];
  float cosmicmuon_InnerTrack_InnerPoint_pz[200];
  float cosmicmuon_InnerTrack_OuterPoint_x[200];
  float cosmicmuon_InnerTrack_OuterPoint_y[200];
  float cosmicmuon_InnerTrack_OuterPoint_z[200];
  float cosmicmuon_InnerTrack_OuterPoint_px[200];
  float cosmicmuon_InnerTrack_OuterPoint_py[200];
  float cosmicmuon_InnerTrack_OuterPoint_pz[200];
 
  //FOR AOD only
  float cosmicmuon_OuterPoint_x[200];
  float cosmicmuon_OuterPoint_y[200];
  float cosmicmuon_OuterPoint_z[200];


 
  //tau variables
  float tau_pt[100];
  float tau_px[100];
  float tau_py[100];
  float tau_pz[100];
  float tau_vx[100];
  float tau_vy[100];
  float tau_vz[100];
  float tau_energy[100];
  float tau_charge[100];
  float tau_eta[100];
  float tau_phi[100];
 
  //detailed infor 
  int nPions[100];                                                                                                                                  
  int nPi0[100];
  int nPhotons[100];
  int oneProng0Pi0[100];
  int oneProng1Pi0[100];
  int oneProng2Pi0[100];
  int threeProng0Pi0[100];
  int threeProng1Pi0[100];
  int tauelectron[100];
  int taumuon[100];
         
  int nthreeProng0Pi0;
  int nthreeProng1Pi0;
  int ntauelectron;
  int ntaumuon;

  int PionPdgId[100][5];
  int PhotonPdgId[100][5];
  int Pi0PdgId[100][5];
        
  float PionPt[100][5];
  float PionEta[100][5];
  float PionPhi[100][5];
  float Pi0Pt[100][5];
  float Pi0Eta[100][5];
  float Pi0Phi[100][5];
  float PhotonPt[100][5];
  float PhotonEta[100][5];
  float PhotonPhi[100][5];
        
  float genHadTauPt[100];
  float genHadTauEta[100];
  float genHadTauPhi[100];
  float oneProng0Pi0Pt[100];
  float oneProng0Pi0Eta[100];
  float oneProng0Pi0Phi[100];

  std::vector <const reco::GenParticle*>  genParticleList ;
  std::vector<std::string> genTauDecayMode1;                                                                                                          
  std::string genTauDecayMode;


 
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
  
  float pho_E[200];
  float pho_pt[200];
  float pho_eta[200];
  float pho_phi[200];
  float pho_px[200];
  float pho_py[200];
  float pho_pz[200];
  float pho_vx[200];
  float pho_vy[200];
  float pho_vz[200];
  float pho_r9[200];
  bool  pho_isEB[200];
  bool  pho_isEE[200];
  bool  pho_isEBGap[200];
  bool  pho_isEEGap[200];
  bool  pho_isEBEEGap[200];
  float pho_e1x5[200];
  float pho_e2x5[200];
  float pho_e3x3[200];
  float pho_e5x5[200];
  float pho_r1x5[200];
  float pho_r2x5[200];
  float pho_SigmaEtaEta[200];
  float pho_SigmaIetaIeta[200];
  float pho_SigmaEtaPhi[200];
  float pho_SigmaIetaIphi[200];
  float pho_SigmaPhiPhi[200];
  float pho_SigmaIphiIphi[200];
  float pho_roundness[200];
  float pho_angle[200];
  float pho_maxEnergyXtal[200];
  float pho_theta[200];
  float pho_et[200];
  float pho_swissCross[200];
  float pho_e6e2[200];
  float pho_e4e1[200];
  bool  pho_isConverted[200];
  bool  pho_hasConvTrk[200]; 
 
  //isolation variables
  float pho_ecalRecHitSumEtConeDR03[200];
  float pho_hcalTowerSumEtConeDR03[200];
  float pho_hcalDepth1TowerSumEtConeDR03[200];
  float pho_hcalDepth2TowerSumEtConeDR03[200];
  float pho_trkSumPtSolidConeDR03[200];
  float pho_trkSumPtHollowConeDR03[200];
  int   pho_nTrkSolidConeDR03[200];
  int   pho_nTrkHollowConeDR03[200];
  float pho_ecalRecHitSumEtConeDR04[200];
  float pho_hcalTowerSumEtConeDR04[200];
  float pho_hcalDepth1TowerSumEtConeDR04[200];
  float pho_hcalDepth2TowerSumEtConeDR04[200];
  float pho_trkSumPtSolidConeDR04[200];
  float pho_trkSumPtHollowConeDR04[200];
  int   pho_nTrkSolidConeDR04[200];
  int   pho_nTrkHollowConeDR04[200];
  float pho_HoE[200];
  bool  pho_hasPixelSeed[200];   
  
  //SC variables
  float pho_sc_energy[200];
  int   pho_size[200];
  float pho_sc_eta[200];
  float pho_sc_phi[200];
  float pho_sc_etaWidth[200];
  float pho_sc_phiWidth[200];
  float pho_sc_et[200];
  float pho_sc_x[200];
  float pho_sc_y[200];
  float pho_sc_z[200];
  //gen matched photon 
  float matchpho_E[200];
  float matchpho_pt[200];
  float matchpho_eta[200];
  float matchpho_phi[200];
  float matchpho_px[200];
  float matchpho_py[200];
  float matchpho_pz[200];
  bool  ismatchedpho[200];
  
  //converted photon variabes
  unsigned int pho_nTracks[200];
  float pho_pairInvariantMass[200];
  float pho_pairCotThetaSeparation[200];
  float pho_pairMomentum_x[200];
  float pho_pairMomentum_y[200];
  float pho_pairMomentum_z[200];
  float pho_EoverP[200];
  float pho_conv_vx[200];
  float pho_conv_vy[200];
  float pho_conv_vz[200];
  float pho_zOfPrimaryVertex[200];
  float pho_distOfMinimumApproach[200];
  float pho_dPhiTracksAtVtx[200];      
  float pho_dPhiTracksAtEcal[200];     
  float pho_dEtaTracksAtEcal[200];     
  
  //rechit information
  int ncrysPhoton[200];
  float pho_timing_xtal[200][100]; 
  float pho_timeError_xtal[200][100]; 
  float pho_timingavg_xtal[200];
  float pho_energy_xtal[200][100];
  int   pho_ieta_xtalEB[200][100];
  int   pho_iphi_xtalEB[200][100];
  int   pho_recoFlag_xtalEB[200][100];
  float pho_rookFraction[200];
  float pho_s9[200];
  float pho_e2e9[200];
  
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
  
  //CSCseg info
  float CSCseg_time[10000];
  float CSCseg_x[10000];
  float CSCseg_y[10000];
  float CSCseg_z[10000];
  float CSCseg_phi[10000];
  float CSCseg_DirectionX[10000];
  float CSCseg_DirectionY[10000];
  float CSCseg_DirectionZ[10000];
 
 //BeamHaloSummary                                                                                                              
  bool isBeamHaloIDTightPass;
  bool isBeamHaloIDLoosePass;
         
  bool isBeamHaloEcalLoosePass;
  bool isBeamHaloHcalLoosePass;
  bool isBeamHaloCSCLoosePass;
  bool isBeamHaloGlobalLoosePass;
         
  bool isBeamHaloEcalTightPass;
  bool isBeamHaloHcalTightPass;
  bool isBeamHaloCSCTightPass;
  bool isBeamHaloGlobalTightPass;

  bool isSmellsLikeHalo_Tag;
  bool isLooseHalo_Tag;
  bool isTightHalo_Tag;
  bool isExtremeTightHalo_Tag;

 
  //RPChit info



  float RPChit_x[10000];
  float RPChit_y[10000];
  float RPChit_z[10000];
  int RPChit_BunchX[10000];
  
  //calomet variables
  float CaloMetSig;
  float CaloMetCorr;
  float CaloMetPt[6];
  float CaloMetPx[6];
  float CaloMetPy[6];
  float CaloMetPhi[6];
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
  float CaloMetSumEt[6];
  float CaloMaxEtInHadTowers;
  float Delta_phi;
  
  //PFMet variables
  float PFMetPt[6];
  float PFMetPx[6];
  float PFMetPy[6];
  float PFMetPhi[6];
  float PFMetSumEt[6];
  float Delta_phiPF;
  
  //TCMet variables
  float TCMetPt[6];
  float TCMetPx[6];
  float TCMetPy[6];
  float TCMetPhi[6];
  float TCMetSumEt[6];
  float Delta_phiTC;
  
  //genMet variables
  float genMetPt;
  float genMetPx;
  float genMetPy;
  float genMetPhi;
  float genMetSumEt;
  float Delta_phiGEN;

  //EBrecHit variables;

  int EBRecHit_size;
  float EBRecHit_eta[10000];
  float EBRecHit_phi[10000];
  int EBRecHit_ieta[10000];
  int EBRecHit_iphi[10000];
  float EBRecHit_e[10000];
  float EBRecHit_et[10000];
  int EBRecHit_flag[10000];
  float EBRecHit_time[10000];

  

  //EErecHit variables;
  int EERecHit_size;
  float EERecHit_eta[10000];
  float EERecHit_phi[10000];
  int EERecHit_ieta[10000];
  int EERecHit_iphi[10000];
  float EERecHit_e[10000];
  float EERecHit_et[10000];
  int EERecHit_flag[10000];
  float EERecHit_time[10000];

 //pileup info
 int npuVertices;
 int ootnpuVertices;




//Uncleand Photon variables
  float ucpho_E[200];
  float ucpho_pt[200];
  float ucpho_eta[200];
  float ucpho_phi[200];
  float ucpho_px[200];
  float ucpho_py[200];
  float ucpho_pz[200];
  float ucpho_vx[200];
  float ucpho_vy[200];
  float ucpho_vz[200];
  float ucpho_r9[200];
  bool  ucpho_isEB[200];
  bool  ucpho_isEE[200];
  bool  ucpho_isEBGap[200];
  bool  ucpho_isEEGap[200];
  bool  ucpho_isEBEEGap[200];
  float ucpho_e1x5[200];
  float ucpho_e2x5[200];
  float ucpho_e3x3[200];
  float ucpho_e5x5[200];
  float ucpho_r1x5[200];
  float ucpho_r2x5[200];
  float ucpho_SigmaEtaEta[200];
  float ucpho_SigmaIetaIeta[200];
  float ucpho_SigmaEtaPhi[200];
  float ucpho_SigmaIetaIphi[200];
  float ucpho_SigmaPhiPhi[200];
  float ucpho_SigmaIphiIphi[200];
  float ucpho_roundness[200];
  float ucpho_angle[200];
  float ucpho_maxEnergyXtal[200];
  float ucpho_theta[200];
  float ucpho_et[200];
  float ucpho_swissCross[200];
  float ucpho_e6e2[200];
  float ucpho_e4e1[200];
  bool  ucpho_isConverted[200];
  bool  ucpho_hasConvTrk[200]; 
             
  //isolation variables
  float ucpho_ecalRecHitSumEtConeDR03[200];
  float ucpho_hcalTowerSumEtConeDR03[200];
  float ucpho_hcalDepth1TowerSumEtConeDR03[200];
  float ucpho_hcalDepth2TowerSumEtConeDR03[200];
  float ucpho_trkSumPtSolidConeDR03[200];
  float ucpho_trkSumPtHollowConeDR03[200];
  int   ucpho_nTrkSolidConeDR03[200];
  int   ucpho_nTrkHollowConeDR03[200];
  float ucpho_ecalRecHitSumEtConeDR04[200];
  float ucpho_hcalTowerSumEtConeDR04[200];
  float ucpho_hcalDepth1TowerSumEtConeDR04[200];
  float ucpho_hcalDepth2TowerSumEtConeDR04[200];
  float ucpho_trkSumPtSolidConeDR04[200];
  float ucpho_trkSumPtHollowConeDR04[200];
  int   ucpho_nTrkSolidConeDR04[200];
  int   ucpho_nTrkHollowConeDR04[200];
  float ucpho_HoE[200];
  bool  ucpho_hasPixelSeed[200];   
             
  //SC variables
  float ucpho_sc_energy[200];
  int   ucpho_size[200];
  float ucpho_sc_eta[200];
  float ucpho_sc_phi[200];
  float ucpho_sc_etaWidth[200];
  float ucpho_sc_phiWidth[200];
  float ucpho_sc_et[200];
  float ucpho_sc_x[200];
  float ucpho_sc_y[200];
  float ucpho_sc_z[200];
  //gen matched photon 
  float matchucpho_E[200];
  float matchucpho_pt[200];
  float matchucpho_eta[200];
  float matchucpho_phi[200];
  float matchucpho_px[200];
  float matchucpho_py[200];
  float matchucpho_pz[200];
  bool  ismatcheducpho[200];

  //converted photon variabes
  unsigned int ucpho_nTracks[200];
  float ucpho_pairInvariantMass[200];
  float ucpho_pairCotThetaSeparation[200];
  float ucpho_pairMomentum_x[200];
  float ucpho_pairMomentum_y[200];
  float ucpho_pairMomentum_z[200];
  float ucpho_EoverP[200];
  float ucpho_conv_vx[200];
  float ucpho_conv_vy[200];
  float ucpho_conv_vz[200];
  float ucpho_zOfPrimaryVertex[200];
  float ucpho_distOfMinimumApproach[200];
  float ucpho_dPhiTracksAtVtx[200];      
  float ucpho_dPhiTracksAtEcal[200];     
  float ucpho_dEtaTracksAtEcal[200];     
             
  //rechit information
  int uc_ncrysPhoton[200];
  float ucpho_timing_xtal[200][100];
  float ucpho_timeError_xtal[200][100];
  float ucpho_timingavg_xtal[200];
  float ucpho_energy_xtal[200][100];
  int   ucpho_ieta_xtalEB[200][100];
  int   ucpho_iphi_xtalEB[200][100];
  int   ucpho_recoFlag_xtalEB[200][100];
  float ucpho_rookFraction[200];
  float ucpho_s9[200];
  float ucpho_e2e9[200];
             
  //HErechit information
  unsigned int ucHERecHit_subset_detid[10000];
  float ucHERecHit_subset_energy[10000];
  float ucHERecHit_subset_time[10000];
  float ucHERecHit_subset_depth[10000];
  float ucHERecHit_subset_phi[10000];
  float ucHERecHit_subset_eta[10000];
  float ucHERecHit_subset_x[10000];
  float ucHERecHit_subset_y[10000];
  float ucHERecHit_subset_z[10000];

  //CaloTwoers
  float caloTower_eta[5000];                      
  float caloTower_phi[5000];                      
  float caloTower_E[5000];                   
  float caloTower_Et[5000];
  float caloTower_emEnergy[5000];
  float caloTower_hadEnergy[5000];
  float caloTower_p[5000];
  float caloTower_EMEt[5000];
  float caloTower_HadEt[5000];
  float caloTower_HadPhi[5000];
  float caloTower_HadEta[5000];
  float caloTower_EMPhi[5000];
  float caloTower_EMEta[5000];
  float caloTower_HadX[5000];
  float caloTower_HadY[5000];
  float caloTower_HadZ[5000];
  float caloTower_HE_E[5000];
  float caloTower_HB_E[5000];
  float caloTower_EMTime[5000];
  float caloTower_HadTime[5000];
  //float caloTower_recoFlag[5000];


  float rho;
  float sigma;
 
};

