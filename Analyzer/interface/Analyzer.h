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
  static double rookFractionBarrelCalculator( const reco::SuperCluster & ,const EcalRecHitCollection &);
  
  
  // ----------member data ---------------------------
  edm::ESHandle<CaloTopology> theCaloTopo_;
  int nevents;
  int is_signal_event, is_Z_event, is_W_event;
  int is_Znunu_event, is_Zelec_event, is_Zmu_event, is_Ztau_event ;
  int is_Welec_event, is_Wmu_event, is_Wtau_event ;
  int is_SingleHardPhoton_event;
  int is_diphoton_event;
 
  int n_signal_events,n_Z_events,n_W_events;
  int n_Zelec_events, n_Zmu_events, n_Ztau_events, n_Znunu_events; 
  int n_Welec_events, n_Wmu_events, n_Wtau_events;
  int n_diphoton_events, n_SingleHardPhoton_events;

  int RunNumber,EventNumber;
  int ngenphotons;
  int nhardphotons;
  int Photon_n;
  int Vertex_n;
  int Muon_n;
  int Tau_n;
  int Electron_n;
  int Track_n;
  int Jet_n;
  //HLT

  edm::TriggerNames triggerNames_;  // TriggerNames class
  std::vector<std::string>  hlNames_;  // name of each HLT algorithm
 
  int HLT_MET50_event;
  int HLT_MET75_event;
  int HLT_Photon15_event;
  int HLT_Photon25_event;
  int HLT_DoubleEle10_event;
  int HLT_DoubleMu3_event;

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
  bool runjets_;
  bool runtaus_;
  bool runHLT_;
  bool runL1_;
  bool runtracks_;
  bool runrechit_;
  bool runvertex_;
  bool init_;

  //output root file
  TFile *f;
  
  // output root tree
  TTree* myEvent;
  //variables to be filled in the tree

  //vertex variables
  double vx[10];
  double vy[10];
  double vz[10];
  double chi2[10];
  double vtracksize[10];
  double vndof[10];
  double v_isFake[10];
  double v_d0[10];

  //track variables
  double trk_pt[200];
  double trk_px[200];
  double trk_py[200];
  double trk_pz[200];
  double trk_eta[200];
  double trk_phi[200];

  //jet variables
  double jet_pt[100];
  double jet_px[100];
  double jet_py[100];
  double jet_pz[100];
  double jet_eta[100];
  double jet_phi[100];
  double jet_emEnergyFraction[100];
  double jet_energyFractionHadronic[100];

  //electron variables
  double electron_pt[100];
  double electron_px[100];
  double electron_py[100];
  double electron_pz[100];
  double electron_energy[100];
  double electron_charge[100];
  double electron_eta[100];
  double electron_phi[100];
  double electron_trkIso[100];

  //muon variables
  double muon_pt[100];
  double muon_px[100];
  double muon_py[100];
  double muon_pz[100];
  double muon_energy[100];
  double muon_charge[100];
  double muon_eta[100];
  double muon_phi[100];

  //tau variables
  double tau_pt[100];
  double tau_px[100];
  double tau_py[100];
  double tau_pz[100];
  double tau_energy[100];
  double tau_charge[100];
  double tau_eta[100];
  double tau_phi[100];

  //gen level variables
  double gen_pho_pt[1000];
  double gen_pho_px[1000];
  double gen_pho_py[1000];
  double gen_pho_pz[1000];
  double gen_pho_phi[1000];
  double gen_pho_eta[1000];
  double gen_pho_E[1000];
  int    gen_pho_status[1000];
  int gen_pho_motherID[1000];
  int gen_pho_motherStatus[1000];
  double gen_pho_motherPt[1000];
  double gen_pho_motherEta[1000];
  double gen_pho_motherPhi[1000];
  int gen_pho_GrandmotherID[1000];
  int gen_pho_GrandmotherStatus[1000];
  double gen_pho_GrandmotherPt[1000];
  double gen_pho_GrandmotherEta[1000];
  double gen_pho_GrandmotherPhi[1000];

  double gen_Hpho_pt[2];
  double gen_Hpho_px[2];
  double gen_Hpho_py[2];
  double gen_Hpho_pz[2];
  double gen_Hpho_phi[2];
  double gen_Hpho_eta[2];
  double gen_Hpho_E[2];

  //gen level variables for Graviton
  double gen_graviton_pt;
  double gen_graviton_px;
  double gen_graviton_py;
  double gen_graviton_pz;
  double gen_graviton_phi;
  double gen_graviton_eta;
  double gen_graviton_E;
  
  //gen level variables for Zdaughter
  double gen_Zdaughter_pt[2];
  double gen_Zdaughter_px[2];
  double gen_Zdaughter_py[2];
  double gen_Zdaughter_pz[2];
  double gen_Zdaughter_phi[2];
  double gen_Zdaughter_eta[2];
  double gen_Zdaughter_E[2];
  int gen_Zdaughter_charge[2];
  int gen_Zdaughter_ID[2];

  //gen level variables for Z
  double gen_Zboson_pt;
  double gen_Zboson_px;
  double gen_Zboson_py;
  double gen_Zboson_pz;
  double gen_Zboson_phi;
  double gen_Zboson_eta;
  double gen_Zboson_E;

  //gen level variables for Wdaughter
  double gen_Wdaughter_pt[2];
  double gen_Wdaughter_px[2];
  double gen_Wdaughter_py[2];
  double gen_Wdaughter_pz[2];
  double gen_Wdaughter_phi[2];
  double gen_Wdaughter_eta[2];
  double gen_Wdaughter_E[2];
  int gen_Wdaughter_charge[2];
  int gen_Wdaughter_ID[2];

  //gen level variables for W
  double gen_Wboson_pt;
  double gen_Wboson_px;
  double gen_Wboson_py;
  double gen_Wboson_pz;
  double gen_Wboson_phi;
  double gen_Wboson_eta;
  double gen_Wboson_E;
  int gen_Wboson_charge;
  int gen_Wboson_ID;

  //gen level variables for  mu/tau daughter
  double gen_Muon_ID[3];
  double gen_Muon_Status[3];
  double gen_Muon_Pt[3];
  double gen_MuonDaughter_pt[3];
  double gen_MuonDaughter_px[3];
  double gen_MuonDaughter_py[3];
  double gen_MuonDaughter_pz[3];
  double gen_MuonDaughter_phi[3];
  double gen_MuonDaughter_eta[3];
  double gen_MuonDaughter_E[3];
  int gen_MuonDaughter_charge[3];
  int gen_MuonDaughter_status[3];
  int gen_MuonDaughter_ID[3];
 
  double gen_tau_ID[3];
  double gen_tau_Status[3];
  double gen_tau_Pt[3];
  double gen_tauDaughter_pt[3];
  double gen_tauDaughter_px[3];
  double gen_tauDaughter_py[3];
  double gen_tauDaughter_pz[3];
  double gen_tauDaughter_phi[3];
  double gen_tauDaughter_eta[3];
  double gen_tauDaughter_E[3];
  int gen_tauDaughter_charge[3];
  int gen_tauDaughter_status[3];
  int gen_tauDaughter_ID[3];

  double pho_E[100];
  double pho_pt[100];
  double pho_eta[100];
  double pho_phi[100];
  double pho_px[100];
  double pho_py[100];
  double pho_pz[100];
  double pho_r9[100];
  double pho_isEB[100];
  double pho_isEE[100];
  double pho_isEBGap[100];
  double pho_isEEGap[100];
  double pho_isEBEEGap[100];
  double pho_e1x5[100];
  double pho_e2x5[100];
  double pho_e3x3[100];
  double pho_e5x5[100];
  double pho_r1x5[100];
  double pho_r2x5[100];
  double pho_SigmaEtaEta[100];
  double pho_SigmaIetaIeta[100];
  double pho_roundness[100];
  double pho_angle[100];
  double pho_maxEnergyXtal[100];
  double pho_theta[100];
  double pho_et[100];
  double pho_swissCross[100];
  int    pho_isConverted[100];
  
  //isolation variables
  double pho_ecalRecHitSumEtConeDR03[100];
  double pho_hcalTowerSumEtConeDR03[100];
  double pho_hcalDepth1TowerSumEtConeDR03[100];
  double pho_hcalDepth2TowerSumEtConeDR03[100];
  double pho_trkSumPtSolidConeDR03[100];
  double pho_trkSumPtHollowConeDR03[100];
  int   pho_nTrkSolidConeDR03[100];
  int   pho_nTrkHollowConeDR03[100];
  double pho_ecalRecHitSumEtConeDR04[100];
  double pho_hcalTowerSumEtConeDR04[100];
  double pho_hcalDepth1TowerSumEtConeDR04[100];
  double pho_hcalDepth2TowerSumEtConeDR04[100];
  double pho_trkSumPtSolidConeDR04[100];
  double pho_trkSumPtHollowConeDR04[100];
  int   pho_nTrkSolidConeDR04[100];
  int   pho_nTrkHollowConeDR04[100];
  double   pho_HoE[100];
  int pho_hasPixelSeed[100];   

  //SC variables
  double pho_sc_energy[100];
  int    pho_size[100];
  double pho_sc_eta[100];
  double pho_sc_phi[100];
  double pho_sc_etaWidth[100];
  double pho_sc_phiWidth[100];
  double pho_sc_et[100];

  //gen matched photon 
  double matchpho_E[100];
  double matchpho_pt[100];
  double matchpho_eta[100];
  double matchpho_phi[100];
  double matchpho_px[100];
  double matchpho_py[100];
  double matchpho_pz[100];
  int    ismatchedpho[100];
  
  //converted photon variabes
  unsigned int pho_nTracks[100];
  double pho_pairInvariantMass[100];
  double pho_pairCotThetaSeparation[100];
  double pho_pairMomentum_x[100];
  double pho_pairMomentum_y[100];
  double pho_pairMomentum_z[100];
  double pho_EoverP[100];
  double pho_vertex_x[100];
  double pho_vertex_y[100];
  double pho_vertex_z[100];
  double pho_zOfPrimaryVertex[100];
  double pho_distOfMinimumApproach[100];
  double pho_dPhiTracksAtVtx[100];      
  double pho_dPhiTracksAtEcal[100];     
  double pho_dEtaTracksAtEcal[100];     
  
  //rechit information
  int ncrysPhoton[100];
  double pho_timing_xtal[100][100];
  double pho_timingavg_xtal[100];
  double pho_energy_xtal[100][100];
  int    pho_ieta_xtalEB[100][100];
  int    pho_iphi_xtalEB[100][100];
  double pho_rookFraction[100];
  double pho_s9[100];

  //calomet variables
  double CaloMetSig;
  double CaloMetCorr;
  double CaloMetEt;
  double CaloMetEx;
  double CaloMetEy;
  double CaloMetPhi;
  double CaloEtFractionHadronic;
  double CaloEmEtFraction;
  double CaloHadEtInHB;
  double CaloHadEtInHO;
  double CaloHadEtInHE;
  double CaloHadEtInHF;
  double CaloEmEtInEB;
  double CaloEmEtInEE;
  double CaloEmEtInHF;
  double CaloMetEz;
  double CaloMaxEtInEmTowers;
  double CaloSumEt;
  double CaloMaxEtInHadTowers;
  double Delta_phi;

  //PFMet variables
  double PFMetPt;
  double PFMetPx;
  double PFMetPy;
  double PFMetPhi;
  double PFMetSumEt;
  double Delta_phiPF;

  //TCMet variables
  double TCMetPt;
  double TCMetPx;
  double TCMetPy;
  double TCMetPhi;
  double TCMetSumEt;
  double Delta_phiTC;

  //genMet variables
  double genMetPt;
  double genMetPx;
  double genMetPy;
  double genMetPhi;
  double genMetSumEt;
  double Delta_phiGEN;

};

