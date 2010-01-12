#include "TFile.h"
#include "TBrowser.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TTree.h"
#include "/afs/cern.ch/user/s/sandhya/scratch0/CMSSW_3_1_2/src/Analysis/Analyzer/interface/CrystalInfo.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include <string>

class Analyzer : public edm::EDAnalyzer {
 public:
  explicit Analyzer(const edm::ParameterSet&);
  ~Analyzer();


 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
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

  int ngenphotons;
  int nhardphotons;
  int nrecPhotons;
  int ncrysPhoton;
  int nvertices;
  int nmuons;
  int ntaus;
  int nelectrons;
  int nTracks;
  int njets;
  //HLT

  edm::TriggerNames triggerNames_;  // TriggerNames class
  std::vector<std::string>  hlNames_;  // name of each HLT algorithm
 
  int is_HLT_MET50_event;
  int is_HLT_MET75_event;
  int is_HLT_Photon25_event;
  std::map<std::string,int> HLT_chosen;
  std::map<std::string,int> L1_chosen;

  std::vector<std::string> JET_CORR;

  edm::InputTag eleLabel_;
  edm::InputTag muoLabel_;
  edm::InputTag jetLabel_;
  edm::InputTag tauLabel_;
  edm::InputTag metLabel_;
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
  bool runelectrons_;
  bool runmuons_;
  bool runjets_;
  bool runtaus_;
  bool runHLT_;
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
  float vx[10];
  float vy[10];
  float vz[10];
  float chi2[10];

  //track variables
  double trk_pt[100];
  double trk_px[100];
  double trk_py[100];
  double trk_pz[100];
  double trk_eta[100];
  double trk_phi[100];

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
  double gen_pho_pt[250];
  double gen_pho_px[250];
  double gen_pho_py[250];
  double gen_pho_pz[250];
  double gen_pho_phi[250];
  double gen_pho_eta[250];
  double gen_pho_E[250];
  int gen_pho_motherID[250];
  double gen_pho_motherPt[250];
  double gen_pho_motherEta[250];
  double gen_pho_motherPhi[250];

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
  double gen_Lepton_MotherID[3];
  double gen_Lepton_MotherStatus[3];
  double gen_Lepton_pt[3];
  double gen_Lepton_px[3];
  double gen_Lepton_py[3];
  double gen_Lepton_pz[3];
  double gen_Lepton_phi[3];
  double gen_Lepton_eta[3];
  double gen_Lepton_E[3];
  int gen_Lepton_charge[3];
  int gen_Lepton_status[3];
  int gen_Lepton_ID[3];

  double pho_E[100];
  double pho_pt[100];
  double pho_eta[100];
  double pho_phi[100];
  double pho_px[100];
  double pho_py[100];
  double pho_pz[100];
  double pho_r9[100];
  float  pho_e1x5[100];
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
  double pho_theta[100];
  double pho_et[100];
  int    pho_isConverted[100];
  
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
  double   pho_HoE[100];
  

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
  
  /*
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
  */

  //rechit information
  double pho_timing_xtalEB[100][100];
  double pho_timingavg_xtalEB[100];
  double pho_energy_xtalEB[100][100];
  double pho_ieta_xtalEB[100][100];
  double pho_iphi_xtalEB[100][100];

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
};

