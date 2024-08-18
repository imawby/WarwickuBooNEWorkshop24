////////////////////////////////////////////////////////////////////////
// Class:       VisualiseSlice
// Plugin Type: analyzer (art v3_01_02)
// File:        VisualiseSlice_module.cc
//
// Generated at Mon Sep  4 08:24:35 2023 by Isobel Mawby using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

// ROOT
#include "TTree.h"

// Services
#include "art/Framework/Services/Optional/TFileService.h"

// Hyperon
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleG4Truth.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleGeneratorTruth.h"

// larpandora
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

// lardataobj
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"

// nusimdata
#include "nusimdata/SimulationBase/MCParticle.h"

// searchingfornues
#include "ubana/searchingfornues/Selection/CommonDefs/PIDFuncs.h"
#include "ubana/searchingfornues/Selection/CommonDefs/LLR_PID.h"
#include "ubana/searchingfornues/Selection/CommonDefs/LLRPID_proton_muon_lookup.h"

namespace pandora {
  class VisualiseSlice;
}

double DEFAULT_INT = -999;
double DEFAULT_DOUBLE = -999.0;
double DEFAULT_BOOL = false;

class pandora::VisualiseSlice : public art::EDAnalyzer {
public:
  explicit VisualiseSlice(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  VisualiseSlice(VisualiseSlice const&) = delete;
  VisualiseSlice(VisualiseSlice&&) = delete;
  VisualiseSlice& operator=(VisualiseSlice const&) = delete;
  VisualiseSlice& operator=(VisualiseSlice&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  enum PandoraView {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W};

  void Reset();
  void ResetPfo();
  void InitialiseTrees(); 
  void FillPandoraMaps(art::Event const& evt);
  void FillMCParticleMaps(art::Event const& evt);
  bool IsEM(const art::Ptr<simb::MCParticle> &mcParticle);
  int GetLeadEMTrackID(const art::Ptr<simb::MCParticle> &mcParticle);
  void FillMCSliceInfo(art::Event const& evt);

  void FillEventInformation(art::Event const& evt);
  void FillEventVisualisationInfo(art::Event const& evt);
  void FillEventNuVertexInfo(art::Event const& evt);
  void TrueToRecoSCE(const double t, TVector3 &inputPosition);
  void FillPFPInformation(art::Event const& evt);
  void FillPFPMatchInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
  void CollectHitsFromClusters(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle, 
      std::vector<art::Ptr<recob::Hit>> &hits);
  void FillPFPVisualisationInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
  PandoraView GetPandora2DView(const art::Ptr<recob::Hit> &hit);
  TVector3 ObtainPandoraHitPosition(const art::Event &evt, const art::Ptr<recob::Hit> hit, 
      const PandoraView hitType);
  TVector3 ProjectIntoPandoraView(const TVector3 &inputPosition3D, const PandoraView pandoraView);
  float YZToU(const float yCoord, const float zCoord);
  float YZToV(const float yCoord, const float zCoord);
  float YZToW(const float yCoord, const float zCoord);
  void FillPFPHitInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
  void FillPFPNuVertexInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
  void FillPFPMiscInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
  void FillPFPShowerInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
  void FillPFPEnergyInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
  double GetMedianValue(const std::vector<float> &inputVector);
  void FillPIDInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
private:

    bool m_applySCECorrections;

    // Product labels
    std::string m_GeneratorModuleLabel;    
    std::string m_MCParticleModuleLabel;
    std::string m_HitModuleLabel;
    std::string m_FlashMatchModuleLabel;
    std::string m_BacktrackModuleLabel;
    std::string m_TrackModuleLabel;
    std::string m_ShowerModuleLabel;
    std::string m_CalorimetryModuleLabel;
    std::string m_PIDModuleLabel;

    //////////////////////////////
    // BDT trees
    //////////////////////////////
    TTree * m_tree;
    //////////////////////////////
    // BDT tree variables
    //////////////////////////////
    // Event vars
    int m_event;
    int m_run;
    int m_subrun;
    // True event vars
    double m_sliceCompleteness;
    double m_slicePurity;
    double m_trueNuVertexX;
    double m_trueNuVertexY;
    double m_trueNuVertexZ;
    double m_uTrueNuVertex_wire;
    double m_vTrueNuVertex_wire;
    double m_wTrueNuVertex_wire;
    // Event hit info
    std::vector<double> m_allUHits_wire;
    std::vector<double> m_allUHits_drift;
    std::vector<int> m_allUHits_owner;
    std::vector<bool> m_allUHits_isInSlice;
    std::vector<double> m_allVHits_wire;
    std::vector<double> m_allVHits_drift;
    std::vector<int> m_allVHits_owner;
    std::vector<bool> m_allVHits_isInSlice;
    std::vector<double> m_allWHits_wire;
    std::vector<double> m_allWHits_drift;
    std::vector<int> m_allWHits_owner;
    std::vector<bool> m_allWHits_isInSlice;
    // Reco event vars
    int m_trueNuSliceID;
    int m_flashMatchNuSliceID;
    double m_recoNuVertexX;
    double m_recoNuVertexY;
    double m_recoNuVertexZ;
    double m_uRecoNuVertex_wire;
    double m_vRecoNuVertex_wire;
    double m_wRecoNuVertex_wire;    
    // TODO: Topological score?
    // True PFP vars
    std::vector<int> m_truePDG_out;
    std::vector<double> m_completeness_out;
    std::vector<double> m_purity_out;
    // Reco PFP vars
    std::vector<double> m_recoVertexX_out;
    std::vector<double> m_recoVertexY_out;
    std::vector<double> m_recoVertexZ_out;
    std::vector<double> m_uRecoVertex_wire_out;
    std::vector<double> m_uRecoVertex_drift_out;
    std::vector<double> m_vRecoVertex_wire_out;
    std::vector<double> m_vRecoVertex_drift_out;
    std::vector<double> m_wRecoVertex_wire_out;
    std::vector<double> m_wRecoVertex_drift_out;
    std::vector<std::vector<double>> m_uHits_wire_out;
    std::vector<std::vector<double>> m_uHits_drift_out;
    std::vector<std::vector<double>> m_vHits_wire_out;
    std::vector<std::vector<double>> m_vHits_drift_out;
    std::vector<std::vector<double>> m_wHits_wire_out;
    std::vector<std::vector<double>> m_wHits_drift_out;
    std::vector<int> m_nSpacePoints_out;
    std::vector<std::vector<double>> m_spacePointsX_out;
    std::vector<std::vector<double>> m_spacePointsY_out;
    std::vector<std::vector<double>> m_spacePointsZ_out;
    std::vector<int> m_generation_out;
    std::vector<int> m_pandoraPFPCode_out;
    std::vector<double> m_trackShowerScore_out;
    std::vector<double> m_nuVertexSeparation_out;
    std::vector<double> m_dca_out;
    std::vector<double> m_nuVertexChargeDistribution_out;
    std::vector<double> m_trackParentSeparation_out;
    std::vector<double> m_showerOpeningAngle_out;
    std::vector<double> m_showerLength_out;
    std::vector<int> m_nHits2D_out;
    std::vector<int> m_nHits3D_out;
    std::vector<int> m_nHitsU_out;
    std::vector<int> m_nHitsV_out;
    std::vector<int> m_nHitsW_out;
    std::vector<double> m_totalEnergy_out;
    std::vector<double> m_initialdEdx_out;
    std::vector<double> m_chiPIDProton_out;
    std::vector<double> m_chiPIDMuon_out;
    std::vector<double> m_chiPIDPion_out;
    std::vector<double> m_braggPIDProton_out;
    std::vector<double> m_braggPIDMuon_out;
    std::vector<double> m_braggPIDPion_out;
    std::vector<double> m_LLRPIDReduced_out;

    //////////////////////////////
    // Analyzer Variables - THIS IS NOT SUPER CLEAN :'(
    //////////////////////////////
    // True PFP vars
    int m_truePDG;
    double m_completeness;
    double m_purity;
    // Reco PFP vars
    double m_recoVertexX;
    double m_recoVertexY;
    double m_recoVertexZ;
    double m_uRecoVertex_wire;
    double m_uRecoVertex_drift;
    double m_vRecoVertex_wire;
    double m_vRecoVertex_drift;
    double m_wRecoVertex_wire;
    double m_wRecoVertex_drift;
    std::vector<double> m_uHits_wire;
    std::vector<double> m_uHits_drift;
    std::vector<double> m_vHits_wire;
    std::vector<double> m_vHits_drift;
    std::vector<double> m_wHits_wire;
    std::vector<double> m_wHits_drift;
    int m_nSpacePoints;
    std::vector<double> m_spacePointsX;
    std::vector<double> m_spacePointsY;
    std::vector<double> m_spacePointsZ;
    int m_generation;
    int m_pandoraPFPCode;
    double m_trackShowerScore;
    double m_nuVertexSeparation;
    double m_dca;
    double m_nuVertexChargeDistribution;
    double m_trackParentSeparation;
    double m_showerOpeningAngle;
    double m_showerLength;
    int m_nHits2D;
    int m_nHits3D;
    int m_nHitsU;
    int m_nHitsV;
    int m_nHitsW;
    double m_totalEnergy;
    double m_initialdEdx;
    double m_chiPIDProton;
    double m_chiPIDMuon;
    double m_chiPIDPion;
    double m_braggPIDProton;
    double m_braggPIDMuon;
    double m_braggPIDPion;
    double m_LLRPIDReduced;

    // FlashMatchNuSlice key
    int m_flashMatchNuSliceKey;

    // Linking hits and true owners
    std::map<int, int> m_hitToTrackID;
    std::map<int, std::vector<int>> m_trackIDToHits;

    // Linking TrackID -> MCParticle
    lar_pandora::MCParticleMap m_mcParticleMap;

    // Linking Self() -> PFParticle
    lar_pandora::PFParticleMap m_pfpMap;

    // LLR_PID
    searchingfornues::LLRPID m_LLRPIDCalculator;
    searchingfornues::ProtonMuonLookUpParameters m_ProtonMuonParams;

    // Get print statements?
    bool m_debug;
};

pandora::VisualiseSlice::VisualiseSlice(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset},
    m_applySCECorrections(pset.get<bool>("ApplySCECorrections")),
    m_GeneratorModuleLabel(pset.get<std::string>("GeneratorModuleLabel")),    
    m_MCParticleModuleLabel(pset.get<std::string>("MCParticleModuleLabel")),
    m_HitModuleLabel(pset.get<std::string>("HitModuleLabel")),
    m_FlashMatchModuleLabel(pset.get<std::string>("FlashMatchModuleLabel")),
    m_BacktrackModuleLabel(pset.get<std::string>("BacktrackModuleLabel")),
    m_TrackModuleLabel(pset.get<std::string>("TrackModuleLabel")),
    m_ShowerModuleLabel(pset.get<std::string>("ShowerModuleLabel")),
    m_CalorimetryModuleLabel(pset.get<std::string>(m_applySCECorrections ? "CalorimetryModuleLabelSCE" : "CalorimetryModuleLabel")),
    m_PIDModuleLabel(pset.get<std::string>(m_applySCECorrections ? "PIDModuleLabelSCE" : "PIDModuleLabel")),
    m_debug(pset.get<bool>("Debug"))
{
    Reset();
    InitialiseTrees();

    m_LLRPIDCalculator.set_dedx_binning(0, m_ProtonMuonParams.dedx_edges_pl_0);
    m_LLRPIDCalculator.set_par_binning(0, m_ProtonMuonParams.parameters_edges_pl_0);
    m_LLRPIDCalculator.set_lookup_tables(0, m_ProtonMuonParams.dedx_pdf_pl_0);

    m_LLRPIDCalculator.set_dedx_binning(1, m_ProtonMuonParams.dedx_edges_pl_1);
    m_LLRPIDCalculator.set_par_binning(1, m_ProtonMuonParams.parameters_edges_pl_1);
    m_LLRPIDCalculator.set_lookup_tables(1, m_ProtonMuonParams.dedx_pdf_pl_1);

    m_LLRPIDCalculator.set_dedx_binning(2, m_ProtonMuonParams.dedx_edges_pl_2);
    m_LLRPIDCalculator.set_par_binning(2, m_ProtonMuonParams.parameters_edges_pl_2);
    m_LLRPIDCalculator.set_lookup_tables(2, m_ProtonMuonParams.dedx_pdf_pl_2);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::Reset() 
{
    //////////////////////////////
    // BDT tree variables
    //////////////////////////////
    // Event vars
    m_event = DEFAULT_INT;
    m_run = DEFAULT_INT;
    m_subrun = DEFAULT_INT;
    // True event vars
    m_sliceCompleteness = DEFAULT_DOUBLE;
    m_slicePurity = DEFAULT_DOUBLE;
    m_trueNuVertexX = DEFAULT_DOUBLE;
    m_trueNuVertexY = DEFAULT_DOUBLE;
    m_trueNuVertexZ = DEFAULT_DOUBLE;
    m_uTrueNuVertex_wire = DEFAULT_DOUBLE;
    m_vTrueNuVertex_wire = DEFAULT_DOUBLE;
    m_wTrueNuVertex_wire = DEFAULT_DOUBLE;
    // Event hit info
    m_allUHits_wire.clear();
    m_allUHits_drift.clear();
    m_allUHits_owner.clear();
    m_allUHits_isInSlice.clear();
    m_allVHits_wire.clear();
    m_allVHits_drift.clear();
    m_allVHits_owner.clear();
    m_allVHits_isInSlice.clear();
    m_allWHits_wire.clear();
    m_allWHits_drift.clear();
    m_allWHits_owner.clear();
    m_allWHits_isInSlice.clear();
    // Reco event vars
    m_trueNuSliceID = DEFAULT_INT;
    m_flashMatchNuSliceID = DEFAULT_INT;
    m_recoNuVertexX = DEFAULT_DOUBLE;
    m_recoNuVertexY = DEFAULT_DOUBLE;
    m_recoNuVertexZ = DEFAULT_DOUBLE;
    m_uRecoNuVertex_wire = DEFAULT_DOUBLE;
    m_vRecoNuVertex_wire = DEFAULT_DOUBLE;
    m_wRecoNuVertex_wire = DEFAULT_DOUBLE;
    // TODO: Topological score?
    // True PFP vars
    m_truePDG_out.clear();
    m_completeness_out.clear();
    m_purity_out.clear();
    // Reco PFP vars
    m_recoVertexX_out.clear();
    m_recoVertexY_out.clear();
    m_recoVertexZ_out.clear();
    m_uRecoVertex_wire_out.clear();
    m_uRecoVertex_drift_out.clear();
    m_vRecoVertex_wire_out.clear();
    m_vRecoVertex_drift_out.clear();
    m_wRecoVertex_wire_out.clear();
    m_wRecoVertex_drift_out.clear();
    m_uHits_wire_out.clear();
    m_uHits_drift_out.clear();
    m_vHits_wire_out.clear();
    m_vHits_drift_out.clear();
    m_wHits_wire_out.clear();
    m_wHits_drift_out.clear();
    m_nSpacePoints_out.clear();
    m_spacePointsX_out.clear();
    m_spacePointsY_out.clear();
    m_spacePointsZ_out.clear();
    m_generation_out.clear();
    m_pandoraPFPCode_out.clear();
    m_trackShowerScore_out.clear();
    m_nuVertexSeparation_out.clear();
    m_dca_out.clear();
    m_nuVertexChargeDistribution_out.clear();
    m_trackParentSeparation_out.clear();
    m_showerOpeningAngle_out.clear();
    m_showerLength_out.clear();
    m_nHits2D_out.clear();
    m_nHits3D_out.clear();
    m_nHitsU_out.clear();
    m_nHitsV_out.clear();
    m_nHitsW_out.clear();
    m_totalEnergy_out.clear();
    m_initialdEdx_out.clear();
    m_chiPIDProton_out.clear();
    m_chiPIDMuon_out.clear();
    m_chiPIDPion_out.clear();
    m_braggPIDProton_out.clear();
    m_braggPIDMuon_out.clear();
    m_braggPIDPion_out.clear();
    m_LLRPIDReduced_out.clear();

    //////////////////////////////
    // Analyzer Variables - THIS IS NOT SUPER CLEAN :'(
    //////////////////////////////
    ResetPfo();

    // FlashMatchNuSlice key
    m_flashMatchNuSliceKey = DEFAULT_INT;
    m_hitToTrackID.clear();
    m_trackIDToHits.clear();
    m_mcParticleMap.clear();
    m_pfpMap.clear();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::ResetPfo() 
{
    m_truePDG = DEFAULT_INT;

    m_recoVertexX = DEFAULT_DOUBLE;
    m_recoVertexY = DEFAULT_DOUBLE;
    m_recoVertexZ = DEFAULT_DOUBLE;
    m_uRecoVertex_wire = DEFAULT_DOUBLE;
    m_uRecoVertex_drift = DEFAULT_DOUBLE;
    m_vRecoVertex_wire = DEFAULT_DOUBLE;
    m_vRecoVertex_drift = DEFAULT_DOUBLE;
    m_wRecoVertex_wire = DEFAULT_DOUBLE;
    m_wRecoVertex_drift = DEFAULT_DOUBLE;
    m_uHits_wire.clear();
    m_uHits_drift.clear();
    m_vHits_wire.clear();
    m_vHits_drift.clear();
    m_wHits_wire.clear();
    m_wHits_drift.clear();
    m_nSpacePoints = DEFAULT_INT;
    m_spacePointsX.clear();
    m_spacePointsY.clear();
    m_spacePointsZ.clear();

    m_completeness = DEFAULT_DOUBLE;
    m_purity = DEFAULT_DOUBLE;
    m_generation = DEFAULT_INT;
    m_pandoraPFPCode = DEFAULT_INT;
    m_nHits2D = DEFAULT_INT;
    m_nHits3D = DEFAULT_INT;
    m_nHitsU = DEFAULT_INT;
    m_nHitsV = DEFAULT_INT;
    m_nHitsW = DEFAULT_INT;
    m_nuVertexSeparation = DEFAULT_DOUBLE;
    m_dca = DEFAULT_DOUBLE;
    m_nuVertexChargeDistribution = DEFAULT_DOUBLE;
    m_totalEnergy = DEFAULT_DOUBLE;
    m_initialdEdx = DEFAULT_DOUBLE;
    m_trackShowerScore = DEFAULT_DOUBLE;
    m_showerOpeningAngle = DEFAULT_DOUBLE;
    m_showerLength = DEFAULT_DOUBLE;
    m_trackParentSeparation = DEFAULT_DOUBLE;
    m_chiPIDProton = DEFAULT_DOUBLE;
    m_chiPIDMuon = DEFAULT_DOUBLE;
    m_chiPIDPion = DEFAULT_DOUBLE;
    m_braggPIDProton = DEFAULT_DOUBLE;
    m_braggPIDMuon = DEFAULT_DOUBLE;
    m_braggPIDPion = DEFAULT_DOUBLE;
    m_LLRPIDReduced = DEFAULT_DOUBLE;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::InitialiseTrees()
{
    art::ServiceHandle<art::TFileService> tfs;

    m_tree = tfs->make<TTree>("VisualisationTree", "VisualisationTree");

    //////////////////////////////
    // BDT tree variables
    //////////////////////////////
    // Event vars
    m_tree->Branch("Event", &m_event);
    m_tree->Branch("Run", &m_run);
    m_tree->Branch("Subrun", &m_subrun);
    // True event vars
    m_tree->Branch("SliceCompleteness", &m_sliceCompleteness);
    m_tree->Branch("SlicePurity", &m_slicePurity);
    m_tree->Branch("TrueNuVertexX", &m_trueNuVertexX);
    m_tree->Branch("TrueNuVertexY", &m_trueNuVertexY);
    m_tree->Branch("TrueNuVertexZ", &m_trueNuVertexZ);    
    m_tree->Branch("UTrueNuVertex_wire", &m_uTrueNuVertex_wire);
    m_tree->Branch("VTrueNuVertex_wire", &m_vTrueNuVertex_wire);
    m_tree->Branch("WTrueNuVertex_wire", &m_wTrueNuVertex_wire);    
    // Event hit info
    m_tree->Branch("AllUHits_wire", &m_allUHits_wire);
    m_tree->Branch("AllUHits_drift", &m_allUHits_drift);
    m_tree->Branch("AllUHits_owner", &m_allUHits_owner);
    m_tree->Branch("AllUHits_isInSlice", &m_allUHits_isInSlice);
    m_tree->Branch("AllVHits_wire", &m_allVHits_wire);
    m_tree->Branch("AllVHits_drift", &m_allVHits_drift);
    m_tree->Branch("AllVHits_owner", &m_allVHits_owner);
    m_tree->Branch("AllVHits_isInSlice", &m_allVHits_isInSlice);
    m_tree->Branch("AllWHits_wire", &m_allWHits_wire);
    m_tree->Branch("AllWHits_drift", &m_allWHits_drift);
    m_tree->Branch("AllWHits_owner", &m_allWHits_owner);
    m_tree->Branch("AllWHits_isInSlice", &m_allWHits_isInSlice);
    // Reco event vars
    m_tree->Branch("TrueNuSliceID", &m_trueNuSliceID);
    m_tree->Branch("FlashMatchNuSliceID", &m_flashMatchNuSliceID);
    m_tree->Branch("RecoNuVertexX", &m_recoNuVertexX);
    m_tree->Branch("RecoNuVertexY", &m_recoNuVertexY);
    m_tree->Branch("RecoNuVertexZ", &m_recoNuVertexZ);
    m_tree->Branch("URecoNuVertex_wire", &m_uRecoNuVertex_wire);
    m_tree->Branch("VRecoNuVertex_wire", &m_vRecoNuVertex_wire);
    m_tree->Branch("WRecoNuVertex_wire", &m_wRecoNuVertex_wire);        
    // TODO: Topological score?
    // True PFP vars
    m_tree->Branch("TruePDG", &m_truePDG_out);
    m_tree->Branch("Completeness", &m_completeness_out);
    m_tree->Branch("Purity", &m_purity_out);
    // Reco PFP vars
    m_tree->Branch("RecoVertexX", &m_recoVertexX_out);
    m_tree->Branch("RecoVertexY", &m_recoVertexY_out);
    m_tree->Branch("RecoVertexZ", &m_recoVertexZ_out);
    m_tree->Branch("URecoVertex_wire", &m_uRecoVertex_wire_out);
    m_tree->Branch("URecoVertex_drift", &m_uRecoVertex_drift_out);
    m_tree->Branch("VRecoVertex_wire", &m_vRecoVertex_wire_out);
    m_tree->Branch("VRecoVertex_drift", &m_vRecoVertex_drift_out);
    m_tree->Branch("WRecoVertex_wire", &m_wRecoVertex_wire_out);
    m_tree->Branch("WRecoVertex_drift", &m_wRecoVertex_drift_out);
    m_tree->Branch("UHits_wire", &m_uHits_wire_out);
    m_tree->Branch("UHits_drift", &m_uHits_drift_out);
    m_tree->Branch("VHits_wire", &m_vHits_wire_out);
    m_tree->Branch("VHits_drift", &m_vHits_drift_out);
    m_tree->Branch("WHits_wire", &m_wHits_wire_out);
    m_tree->Branch("WHits_drift", &m_wHits_drift_out);
    m_tree->Branch("NSpacePoints", &m_nSpacePoints_out);
    m_tree->Branch("SpacePointsX", &m_spacePointsX_out);
    m_tree->Branch("SpacePointsY", &m_spacePointsY_out);
    m_tree->Branch("SpacePointsZ", &m_spacePointsZ_out);
    m_tree->Branch("Generation", &m_generation_out);
    m_tree->Branch("PandoraPFPCode", &m_pandoraPFPCode_out);
    m_tree->Branch("NHits3D", &m_nHits3D_out);
    m_tree->Branch("NHits2D", &m_nHits2D_out);
    m_tree->Branch("NHitsU", &m_nHitsU_out);
    m_tree->Branch("NHitsV", &m_nHitsV_out);
    m_tree->Branch("NHitsW", &m_nHitsW_out);
    m_tree->Branch("NuVertexSeparation", &m_nuVertexSeparation_out);
    m_tree->Branch("DCA", &m_dca_out);
    m_tree->Branch("NuVertexChargeDistribution", &m_nuVertexChargeDistribution_out);
    m_tree->Branch("TotalEnergy", &m_totalEnergy_out);
    m_tree->Branch("InitialdEdx", &m_initialdEdx_out);
    m_tree->Branch("TrackShowerScore", &m_trackShowerScore_out);
    m_tree->Branch("ShowerOpeningAngle", &m_showerOpeningAngle_out);
    m_tree->Branch("ShowerLength", &m_showerLength_out);
    m_tree->Branch("TrackParentSeparation", &m_trackParentSeparation_out);
    m_tree->Branch("ChiPIDProton", &m_chiPIDProton_out);
    m_tree->Branch("ChiPIDMuon", &m_chiPIDMuon_out);
    m_tree->Branch("ChiPIDPion", &m_chiPIDPion_out);
    m_tree->Branch("BraggPIDProton", &m_braggPIDProton_out);
    m_tree->Branch("BraggPIDMuon", &m_braggPIDMuon_out);
    m_tree->Branch("BraggPIDPion", &m_braggPIDPion_out);
    m_tree->Branch("LLRPIDReduced", &m_LLRPIDReduced_out);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::analyze(art::Event const& evt)
{
    Reset();

    m_event = evt.event();
    m_run = evt.run();
    m_subrun = evt.subRun();

    if (m_debug) std::cout << "Filling Pandora Maps..." << std::endl;
    FillPandoraMaps(evt);
    if (m_debug) std::cout << "Filling MCParticle Maps..." << std::endl;
    FillMCParticleMaps(evt);
    if (m_debug) std::cout << "Filling MC Slice Info..." << std::endl;
    FillMCSliceInfo(evt);

    // i.e. we have found a slice...
    if (m_flashMatchNuSliceID >= 0)
    {
        if (m_debug) std::cout << "Filling Event Info..." << std::endl;
        FillEventInformation(evt);
        if (m_debug) std::cout << "Filling PFParticle Information..." << std::endl;
        FillPFPInformation(evt);
    }

    m_tree->Fill();
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillPandoraMaps(art::Event const& evt)
{
    // MCParticle map
    art::Handle<std::vector<simb::MCParticle>> mcParticleHandle;
    std::vector<art::Ptr<simb::MCParticle>> mcParticleVector;

    if (!evt.getByLabel(m_MCParticleModuleLabel, mcParticleHandle))
        throw cet::exception("VisualiseSlice") << "No MCParticle Data Products Found!" << std::endl;

    art::fill_ptr_vector(mcParticleVector, mcParticleHandle);
    lar_pandora::LArPandoraHelper::BuildMCParticleMap(mcParticleVector, m_mcParticleMap);

    // PFParticle map
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    std::vector<art::Ptr<recob::PFParticle>> pfpVector;

    if (!evt.getByLabel(m_FlashMatchModuleLabel, pfpHandle))
        throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

    art::fill_ptr_vector(pfpVector, pfpHandle);

    lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfpVector, m_pfpMap);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillMCParticleMaps(art::Event const& evt)
{
    // Get event hits
    art::Handle<std::vector<recob::Hit>> hitHandle;
    std::vector<art::Ptr<recob::Hit>> hitVector;

    if (!evt.getByLabel(m_HitModuleLabel, hitHandle))
        throw cet::exception("VisualiseSlice") << "No Hit Data Products Found!" << std::endl;

    art::fill_ptr_vector(hitVector, hitHandle);

    // Get backtracker info
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> assocMCPart = art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitHandle, evt, m_BacktrackModuleLabel);

    // Truth match
    for (unsigned int hitIndex = 0; hitIndex < hitVector.size(); hitIndex++)
    {
        const art::Ptr<recob::Hit> &hit = hitVector[hitIndex];
        const std::vector<art::Ptr<simb::MCParticle>> &matchedMCParticles = assocMCPart.at(hit.key());
        auto matchedDatas = assocMCPart.data(hit.key());

        for (unsigned int mcParticleIndex = 0; mcParticleIndex < matchedMCParticles.size(); ++mcParticleIndex)
        {
            const art::Ptr<simb::MCParticle> &matchedMCParticle = matchedMCParticles.at(mcParticleIndex);
            auto matchedData = matchedDatas.at(mcParticleIndex);

            if (matchedData->isMaxIDE != 1)
                continue;

            // Get hit view
            const int trackID = IsEM(matchedMCParticle) ? GetLeadEMTrackID(matchedMCParticle) : matchedMCParticle->TrackId();

            m_hitToTrackID[hit.key()] = trackID;
            m_trackIDToHits[trackID].push_back(hit.key());

        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

bool pandora::VisualiseSlice::IsEM(const art::Ptr<simb::MCParticle> &mcParticle)
{
    return ((std::abs(mcParticle->PdgCode()) == 11) || (mcParticle->PdgCode() == 22));
}

///////////////////////////////////////////////////////////////////////////////////////////

// If its an EM particle, we have to move up the EM hierarchy
int pandora::VisualiseSlice::GetLeadEMTrackID(const art::Ptr<simb::MCParticle> &mcParticle)
{
    int trackID = mcParticle->TrackId();
    art::Ptr<simb::MCParticle> motherMCParticle = mcParticle;

    do
    {
        trackID = motherMCParticle->TrackId();
        const int motherID = motherMCParticle->Mother();

        if (m_mcParticleMap.find(motherID) == m_mcParticleMap.end())
            break;

        motherMCParticle = m_mcParticleMap.at(motherID);
    }
    while (IsEM(motherMCParticle));

    return trackID;
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillMCSliceInfo(art::Event const& evt)
{
    // Get slice information
    art::Handle<std::vector<recob::Slice>> sliceHandle;
    std::vector<art::Ptr<recob::Slice>> sliceVector;

    if (!evt.getByLabel(m_FlashMatchModuleLabel, sliceHandle))
        throw cet::exception("PhotonBDTModuleLabel") << "No Slice Data Products Found!" << std::endl;

    art::FindManyP<recob::Hit> hitAssoc = art::FindManyP<recob::Hit>(sliceHandle, evt, m_FlashMatchModuleLabel);
    art::fill_ptr_vector(sliceVector, sliceHandle);

    /////////////////////////////////////////////
    // Find the true nu slice ID 
    /////////////////////////////////////////////
    int highestHitNumber(-1);
    std::map<int, int> sliceSignalHitMap;
    int totalTrueHits(0);

    for (art::Ptr<recob::Slice> &slice : sliceVector)
    {
        sliceSignalHitMap[slice->ID()] = 0;

        const std::vector<art::Ptr<recob::Hit>> &sliceHits(hitAssoc.at(slice.key()));

        for (const art::Ptr<recob::Hit> &sliceHit : sliceHits)
        {
            if (m_hitToTrackID.find(sliceHit.key()) == m_hitToTrackID.end())
                continue;

            ++sliceSignalHitMap[slice->ID()];
            ++totalTrueHits;
        }

        if ((sliceSignalHitMap[slice->ID()] > highestHitNumber) && (sliceSignalHitMap[slice->ID()] > 0))
        {
            highestHitNumber = sliceSignalHitMap[slice->ID()];
            m_trueNuSliceID = slice->ID();
        }
    }

    /////////////////////////////////////////////
    // Find the flash match nu slice ID 
    /////////////////////////////////////////////
    art::Handle<std::vector<recob::PFParticle>> flashMatchPFPHandle;
    std::vector<art::Ptr<recob::PFParticle>> flashMatchPFPVector;

    if (!evt.getByLabel(m_FlashMatchModuleLabel, flashMatchPFPHandle))
        throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;

    art::fill_ptr_vector(flashMatchPFPVector, flashMatchPFPHandle);

    std::vector<art::Ptr<recob::PFParticle>> neutrinoPFPs;
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(flashMatchPFPVector, neutrinoPFPs);

    if (neutrinoPFPs.size() > 1)
    {
        throw cet::exception("SigmaRecoAnalyser") << "Too many neutrinos found!" << std::endl;
    }
    else if (neutrinoPFPs.size() == 1)
    {
        art::FindManyP<recob::Slice> flashMatchSliceAssoc = art::FindManyP<recob::Slice>(flashMatchPFPHandle, evt, m_FlashMatchModuleLabel);
        const std::vector<art::Ptr<recob::Slice>> &flashMatchSlices = flashMatchSliceAssoc.at(neutrinoPFPs[0].key());

        if (!flashMatchSlices.empty())
        {
            const art::Ptr<recob::Slice> &flashMatchSlice = flashMatchSlices[0];

            m_flashMatchNuSliceID = flashMatchSlice->ID();
            m_flashMatchNuSliceKey = flashMatchSlice.key();

            // Fill completeness and purity
            const std::vector<art::Ptr<recob::Hit>> &sliceHits(hitAssoc.at(flashMatchSlice.key()));

            const int nSliceHits = sliceHits.size();
            const int nSliceTrueHits = sliceSignalHitMap.find(flashMatchSlice->ID()) == sliceSignalHitMap.end() ? 0 : sliceSignalHitMap.at(flashMatchSlice->ID());

            m_sliceCompleteness = totalTrueHits == 0 ? 0.0 : static_cast<float>(nSliceTrueHits) / static_cast<float>(totalTrueHits);
            m_slicePurity = nSliceHits == 0 ? 0.0 : static_cast<float>(nSliceTrueHits) / static_cast<float>(nSliceHits);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillEventInformation(art::Event const& evt)
{
    if (m_debug) std::cout << "Filling Event Visualisation Info..." << std::endl;
    FillEventVisualisationInfo(evt);

    if (m_debug) std::cout << "Filling Event Info..." << std::endl;
    FillEventNuVertexInfo(evt);
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillEventVisualisationInfo(art::Event const& evt)
{
    // Get event hits
    art::Handle<std::vector<recob::Hit>> hitHandle;
    std::vector<art::Ptr<recob::Hit>> hitVector;

    if (!evt.getByLabel(m_HitModuleLabel, hitHandle))
        throw cet::exception("VisualiseSlice") << "No Hit Data Products Found!" << std::endl;

    art::fill_ptr_vector(hitVector, hitHandle);

    // Get hit assoc to slice
    art::FindManyP<recob::Slice> sliceAssn = art::FindManyP<recob::Slice>(hitHandle, evt, m_FlashMatchModuleLabel);

    // Loop through event hits
    for (const art::Ptr<recob::Hit> &hit : hitVector)
    {
        PandoraView pandoraView = GetPandora2DView(hit);

        std::vector<double> &allHits_wire = (pandoraView == TPC_VIEW_U) ? m_allUHits_wire : (pandoraView == TPC_VIEW_V) ? m_allVHits_wire : m_allWHits_wire;
        std::vector<double> &allHits_drift = (pandoraView == TPC_VIEW_U) ? m_allUHits_drift : (pandoraView == TPC_VIEW_V) ? m_allVHits_drift : m_allWHits_drift;
        std::vector<int> &allHits_owner = (pandoraView == TPC_VIEW_U) ? m_allUHits_owner : (pandoraView == TPC_VIEW_V) ? m_allVHits_owner : m_allWHits_owner;
        std::vector<bool> &allHits_isInSlice = (pandoraView == TPC_VIEW_U) ? m_allUHits_isInSlice : (pandoraView == TPC_VIEW_V) ? m_allVHits_isInSlice : m_allWHits_isInSlice;

        // Store who owns it
        if (m_hitToTrackID.find(hit.key()) == m_hitToTrackID.end())
        {
            allHits_owner.push_back(DEFAULT_DOUBLE);
        }
        else
        {
            const int ownerTrackID = m_hitToTrackID.at(hit.key());

            allHits_owner.push_back((m_mcParticleMap.find(ownerTrackID) == m_mcParticleMap.end()) ? DEFAULT_DOUBLE : 
                                    m_mcParticleMap.at(ownerTrackID)->PdgCode());
        }

        // Store position
        TVector3 pandoraPosition = ObtainPandoraHitPosition(evt, hit, pandoraView);

        allHits_wire.push_back(pandoraPosition.Z());
        allHits_drift.push_back(pandoraPosition.X());

        // Store is in slice
        std::vector<art::Ptr<recob::Slice>> sliceOwner = sliceAssn.at(hit.key());

        if (sliceOwner.empty())
            allHits_isInSlice.push_back(false);
        else
            allHits_isInSlice.push_back(sliceOwner[0]->ID() == m_flashMatchNuSliceID);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillEventNuVertexInfo(art::Event const& evt)
{
    ///////////////////////
    // True
    ///////////////////////
    art::Handle<std::vector<simb::MCTruth>> mcTruthHandle;
    std::vector<art::Ptr<simb::MCTruth>> mcTruthVector;

    if (!evt.getByLabel(m_GeneratorModuleLabel, mcTruthHandle))
        throw cet::exception("PhotonBDTModuleLabel") << "No MCTruth Data Products Found!" << std::endl;

    art::fill_ptr_vector(mcTruthVector, mcTruthHandle);

    if (!mcTruthVector.empty())
    {
       const simb::MCNeutrino &mcNeutrino = mcTruthVector[0]->GetNeutrino();
       TVector3 trueNuVertex = TVector3(mcNeutrino.Nu().Vx(), mcNeutrino.Nu().Vy(), mcNeutrino.Nu().Vz());

       TrueToRecoSCE(mcNeutrino.Nu().T(), trueNuVertex);
       
       m_trueNuVertexX = trueNuVertex.X();
       m_trueNuVertexY = trueNuVertex.Y();
       m_trueNuVertexZ = trueNuVertex.Z();       

       const TVector3 uProjection = ProjectIntoPandoraView(trueNuVertex, TPC_VIEW_U);
       const TVector3 vProjection = ProjectIntoPandoraView(trueNuVertex, TPC_VIEW_V);
       const TVector3 wProjection = ProjectIntoPandoraView(trueNuVertex, TPC_VIEW_W);       

       m_uTrueNuVertex_wire = uProjection.Z();
       m_vTrueNuVertex_wire = vProjection.Z();
       m_wTrueNuVertex_wire = wProjection.Z();       
    }

    ///////////////////////
    // Reco
    ///////////////////////
    
    // Get slice information
    art::Handle<std::vector<recob::Slice>> sliceHandle;
    std::vector<art::Ptr<recob::Slice>> sliceVector;

    if (!evt.getByLabel(m_FlashMatchModuleLabel, sliceHandle))
        throw cet::exception("PhotonBDTModuleLabel") << "No Slice Data Products Found!" << std::endl;

    art::FindManyP<recob::PFParticle> pfpAssoc = art::FindManyP<recob::PFParticle>(sliceHandle, evt, m_FlashMatchModuleLabel);
    art::fill_ptr_vector(sliceVector, sliceHandle);

    // Get PFP handle
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_FlashMatchModuleLabel, pfpHandle))
        throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

    for (art::Ptr<recob::Slice> &slice : sliceVector)
    {
        if (slice->ID() == m_flashMatchNuSliceID)
        {
            const std::vector<art::Ptr<recob::PFParticle>> &slicePFPs(pfpAssoc.at(slice.key()));

            std::vector<art::Ptr<recob::PFParticle>> neutrinoPFPs;
            lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(slicePFPs, neutrinoPFPs);

            if (neutrinoPFPs.size() != 1)
                return;

            art::FindManyP<recob::Vertex> vertexAssoc = art::FindManyP<recob::Vertex>(pfpHandle, evt, m_FlashMatchModuleLabel);

            const std::vector<art::Ptr<recob::Vertex>> &nuVertex(vertexAssoc.at(neutrinoPFPs.at(0).key()));

            if (nuVertex.empty())
                return;

            const TVector3 recoNuVertex = TVector3(nuVertex[0]->position().X(), nuVertex[0]->position().Y(), nuVertex[0]->position().Z());

            m_recoNuVertexX = recoNuVertex.X();
            m_recoNuVertexY = recoNuVertex.Y();
            m_recoNuVertexZ = recoNuVertex.Z();       

            const TVector3 uProjection = ProjectIntoPandoraView(recoNuVertex, TPC_VIEW_U);
            const TVector3 vProjection = ProjectIntoPandoraView(recoNuVertex, TPC_VIEW_V);
            const TVector3 wProjection = ProjectIntoPandoraView(recoNuVertex, TPC_VIEW_W);       

            m_uRecoNuVertex_wire = uProjection.Z();
            m_vRecoNuVertex_wire = vProjection.Z();
            m_wRecoNuVertex_wire = wProjection.Z();       

            break;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::TrueToRecoSCE(const double t, TVector3 &inputPosition)
{
    auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    //auto const &detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    //auto const &detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
    
    auto offset = SCE->GetPosOffsets(geo::Point_t(inputPosition.X(), inputPosition.Y(), inputPosition.Z()));
    double newX = inputPosition.X() - offset.X();
    double newY = inputPosition.Y() + offset.Y();
    double newZ = inputPosition.Z() + offset.Z();    

    //double g4Ticks = detClocks->TPCG4Time2Tick(t) + detProperties->GetXTicksOffset(0, 0, 0) - detProperties->TriggerOffset();    
    //float _xtimeoffset = detProperties->ConvertTicksToX(g4Ticks, 0, 0, 0);

    //newX += _xtimeoffset;
    //newX += 0.6;
    
    inputPosition = TVector3(newX, newY, newZ);
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillPFPInformation(art::Event const& evt)
{
    // Get slice information
    art::Handle<std::vector<recob::Slice>> sliceHandle;
    std::vector<art::Ptr<recob::Slice>> sliceVector;

    if (!evt.getByLabel(m_FlashMatchModuleLabel, sliceHandle))
        throw cet::exception("PhotonBDTModuleLabel") << "No Slice Data Products Found!" << std::endl;

    art::fill_ptr_vector(sliceVector, sliceHandle);

    art::FindManyP<recob::PFParticle> pfpAssoc = art::FindManyP<recob::PFParticle>(sliceHandle, evt, m_FlashMatchModuleLabel);

    for (const art::Ptr<recob::Slice> &slice : sliceVector)
    {
        if (slice->ID() != m_flashMatchNuSliceID)
            continue;

        if (m_debug) std::cout << "Found Slice, Filling PFParticle Information..." << std::endl;

        const std::vector<art::Ptr<recob::PFParticle>> &pfps(pfpAssoc.at(slice.key()));

        for (const art::Ptr<recob::PFParticle> &pfp : pfps)
        {
            // Is the PFP in the neutrino hierarchy??
            const art::Ptr<recob::PFParticle> parentPFP = lar_pandora::LArPandoraHelper::GetParentPFParticle(m_pfpMap, pfp);

            if (!lar_pandora::LArPandoraHelper::IsNeutrino(parentPFP))
                continue;

            if (m_debug) std::cout << "New PFParticle!" << std::endl;

            ResetPfo();

            if (m_debug) std::cout << "Fill PFP Match Information..." << std::endl;
            FillPFPMatchInfo(evt, pfp);
            if (m_debug) std::cout << "Fill PFP Visualisation Information..." << std::endl;
            FillPFPVisualisationInfo(evt, pfp);
            if (m_debug) std::cout << "Fill PFP Hit Information..." << std::endl;
            FillPFPHitInfo(evt, pfp);
            if (m_debug) std::cout << "Fill Nu Vertex Information..." << std::endl;
            FillPFPNuVertexInfo(evt, pfp);
            if (m_debug) std::cout << "Fill PFP Misc Information..." << std::endl;
            FillPFPMiscInfo(evt, pfp);
            if (m_debug) std::cout << "Fill PFP Shower Information..." << std::endl;
            FillPFPShowerInfo(evt, pfp);
            if (m_debug) std::cout << "Fill PFP Energy Information..." << std::endl;
            FillPFPEnergyInfo(evt, pfp);
            if (m_debug) std::cout << "Fill PID Information..." << std::endl;
            FillPIDInfo(evt, pfp);
            if (m_debug) std::cout << "Fill Tree..." << std::endl;

            // IKR i hate this too
            m_truePDG_out.push_back(m_truePDG);
            m_completeness_out.push_back(m_completeness);
            m_purity_out.push_back(m_purity);
            m_recoVertexX_out.push_back(m_recoVertexX);
            m_recoVertexY_out.push_back(m_recoVertexY);
            m_recoVertexZ_out.push_back(m_recoVertexZ);
            m_uRecoVertex_wire_out.push_back(m_uRecoVertex_wire);
            m_uRecoVertex_drift_out.push_back(m_uRecoVertex_drift);
            m_vRecoVertex_wire_out.push_back(m_vRecoVertex_wire);
            m_vRecoVertex_drift_out.push_back(m_vRecoVertex_drift);
            m_wRecoVertex_wire_out.push_back(m_wRecoVertex_wire);
            m_wRecoVertex_drift_out.push_back(m_wRecoVertex_drift);
            m_uHits_wire_out.push_back(m_uHits_wire);
            m_uHits_drift_out.push_back(m_uHits_drift);
            m_vHits_wire_out.push_back(m_vHits_wire);
            m_vHits_drift_out.push_back(m_vHits_drift);
            m_wHits_wire_out.push_back(m_wHits_wire);
            m_wHits_drift_out.push_back(m_wHits_drift);
            m_nSpacePoints_out.push_back(m_nSpacePoints);
            m_spacePointsX_out.push_back(m_spacePointsX);
            m_spacePointsY_out.push_back(m_spacePointsY);
            m_spacePointsZ_out.push_back(m_spacePointsZ);
            m_generation_out.push_back(m_generation);
            m_pandoraPFPCode_out.push_back(m_pandoraPFPCode);
            m_nHits3D_out.push_back(m_nHits3D);
            m_nHits2D_out.push_back(m_nHits2D);
            m_nHitsU_out.push_back(m_nHitsU);
            m_nHitsV_out.push_back(m_nHitsV);
            m_nHitsW_out.push_back(m_nHitsW);
            m_nuVertexSeparation_out.push_back(m_nuVertexSeparation);
            m_dca_out.push_back(m_dca);
            m_nuVertexChargeDistribution_out.push_back(m_nuVertexChargeDistribution);
            m_totalEnergy_out.push_back(m_totalEnergy);
            m_initialdEdx_out.push_back(m_initialdEdx);
            m_trackShowerScore_out.push_back(m_trackShowerScore);
            m_showerOpeningAngle_out.push_back(m_showerOpeningAngle);
            m_showerLength_out.push_back(m_showerLength);
            m_trackParentSeparation_out.push_back(m_trackParentSeparation);
            m_chiPIDProton_out.push_back(m_chiPIDProton);
            m_chiPIDMuon_out.push_back(m_chiPIDMuon);
            m_chiPIDPion_out.push_back(m_chiPIDPion);
            m_braggPIDProton_out.push_back(m_braggPIDProton);
            m_braggPIDMuon_out.push_back(m_braggPIDMuon);
            m_braggPIDPion_out.push_back(m_braggPIDPion);
            m_LLRPIDReduced_out.push_back(m_LLRPIDReduced);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillPFPMatchInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    std::vector<art::Ptr<recob::Hit>> pfpHits;
    CollectHitsFromClusters(evt, pfp, pfpHits);

    std::map<int, int> countingMap;
    for (const art::Ptr<recob::Hit> pfpHit : pfpHits)
    {
        if (m_hitToTrackID.find(pfpHit.key()) == m_hitToTrackID.end())
            continue;

        const int trackID(m_hitToTrackID.at(pfpHit.key()));

        if (countingMap.find(trackID) == countingMap.end())
            countingMap[trackID] = 1;
        else
            ++countingMap[trackID];
    }

    if (countingMap.empty())
        return;

    int maxHits = -1;
    int matchedTrackID(-1);

    for (auto &entry : countingMap)
    {
        if ((entry.second > maxHits) || ((entry.second == maxHits) && (entry.first > matchedTrackID)))
        {
            maxHits = entry.second;
            matchedTrackID = entry.first;
        }
    }

    if (m_mcParticleMap.find(matchedTrackID) == m_mcParticleMap.end())
        return;

    const art::Ptr<simb::MCParticle> &matchedMCParticle(m_mcParticleMap.at(matchedTrackID)); 

    m_truePDG = matchedMCParticle->PdgCode();

    const int nTrueHits = m_trackIDToHits.at(matchedTrackID).size();

    m_completeness = static_cast<double>(maxHits) / static_cast<double>(nTrueHits);
    m_purity = static_cast<double>(maxHits) / pfpHits.size();
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::CollectHitsFromClusters(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle, 
   std::vector<art::Ptr<recob::Hit>> &hits)
{
   art::Handle<std::vector<recob::PFParticle>> pfparticleHandle;

   if (!evt.getByLabel(m_FlashMatchModuleLabel, pfparticleHandle))
       throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

   art::Handle<std::vector<recob::Cluster>> clusterHandle;

   if (!evt.getByLabel(m_FlashMatchModuleLabel, clusterHandle)) 
       throw cet::exception("VisualiseSlice") << "No Cluster Data Products Found!" << std::endl;

   art::FindManyP<recob::Cluster> pfparticleClustersAssoc = art::FindManyP<recob::Cluster>(pfparticleHandle, evt, m_FlashMatchModuleLabel);
   art::FindManyP<recob::Hit> clusterHitAssoc = art::FindManyP<recob::Hit>(clusterHandle, evt, m_FlashMatchModuleLabel);

   std::vector<art::Ptr<recob::Cluster>> clusters = pfparticleClustersAssoc.at(pfparticle.key());

   for (const art::Ptr<recob::Cluster> cluster : clusters)
   {
       std::vector<art::Ptr<recob::Hit>> clusterHits = clusterHitAssoc.at(cluster.key());
       hits.insert(hits.end(), clusterHits.begin(), clusterHits.end());
   }
}


///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillPFPVisualisationInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    // PFP vertex
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_FlashMatchModuleLabel, pfpHandle))
        throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<recob::Vertex> vertexAssoc = art::FindManyP<recob::Vertex>(pfpHandle, evt, m_FlashMatchModuleLabel);

    const std::vector<art::Ptr<recob::Vertex>> &pfpVertex(vertexAssoc.at(pfp.key()));

    if (!pfpVertex.empty())
    {
        const TVector3 vertexPosition = TVector3(pfpVertex[0]->position().X(), pfpVertex[0]->position().Y(), pfpVertex[0]->position().Z());

        m_recoVertexX = vertexPosition.X();
        m_recoVertexY = vertexPosition.Y();
        m_recoVertexZ = vertexPosition.Z();

        TVector3 uProjection = ProjectIntoPandoraView(vertexPosition, TPC_VIEW_U);
        TVector3 vProjection = ProjectIntoPandoraView(vertexPosition, TPC_VIEW_V);
        TVector3 wProjection = ProjectIntoPandoraView(vertexPosition, TPC_VIEW_W);

        m_uRecoVertex_wire = uProjection.Z();
        m_uRecoVertex_drift = uProjection.X();
        m_vRecoVertex_wire = vProjection.Z();
        m_vRecoVertex_drift = vProjection.X();
        m_wRecoVertex_wire = wProjection.Z();
        m_wRecoVertex_drift = wProjection.X();
    }

    // Spacepoints
    art::FindManyP<recob::SpacePoint> spacePointAssoc = art::FindManyP<recob::SpacePoint>(pfpHandle, evt, m_FlashMatchModuleLabel);
    const std::vector<art::Ptr<recob::SpacePoint>> &pfpSpacePoints = spacePointAssoc.at(pfp.key());

    m_nSpacePoints = pfpSpacePoints.size();

    for (const art::Ptr<recob::SpacePoint> &spacePoint : pfpSpacePoints)
    {
        m_spacePointsX.push_back(spacePoint->XYZ()[0]);
        m_spacePointsY.push_back(spacePoint->XYZ()[1]);
        m_spacePointsZ.push_back(spacePoint->XYZ()[2]);
    }

    // Hits
    std::vector<art::Ptr<recob::Hit>> pfpHits;
    CollectHitsFromClusters(evt, pfp, pfpHits);

    for (const art::Ptr<recob::Hit> &hit : pfpHits)
    {
        PandoraView pandoraView = GetPandora2DView(hit);
        TVector3 pandoraPosition = ObtainPandoraHitPosition(evt, hit, pandoraView);

        if (pandoraView == TPC_VIEW_W)
        {
            m_wHits_wire.push_back(pandoraPosition.Z());
            m_wHits_drift.push_back(pandoraPosition.X());
        }
        else if (pandoraView == TPC_VIEW_U)
        {
            m_uHits_wire.push_back(pandoraPosition.Z());
            m_uHits_drift.push_back(pandoraPosition.X());
        }
        else if (pandoraView == TPC_VIEW_V)
        {
            m_vHits_wire.push_back(pandoraPosition.Z());
            m_vHits_drift.push_back(pandoraPosition.X());
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

pandora::VisualiseSlice::PandoraView pandora::VisualiseSlice::GetPandora2DView(const art::Ptr<recob::Hit> &hit)
{
    const geo::WireID hitWireID(hit->WireID());
    const geo::View_t hitView(hit->View());
    const geo::View_t thisPandoraView(lar_pandora::LArPandoraGeometry::GetGlobalView(hitWireID.Cryostat, hitWireID.TPC, hitView));

    if (thisPandoraView == geo::kW || thisPandoraView == geo::kY)
        return TPC_VIEW_W;
    else if (thisPandoraView == geo::kU)
        return TPC_VIEW_U;
    else if (thisPandoraView == geo::kV)
        return TPC_VIEW_V;
    else
        throw cet::exception("ivysaur::GridManager") << "wire view not recognised";
}

///////////////////////////////////////////////////////////////////////////////////////////

TVector3 pandora::VisualiseSlice::ObtainPandoraHitPosition(const art::Event &evt, const art::Ptr<recob::Hit> hit, 
    const PandoraView hitType)
{
    art::ServiceHandle<geo::Geometry> theGeometry;
    auto const* theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // Get basic hit properties (view, time, charge)
    const geo::WireID hitWireID(hit->WireID());
    const double hitTime(hit->PeakTime());

    // Get hit X coordinate and, if using a single global drift volume, remove any out-of-time hits here
    const double xCoord = theDetector->ConvertTicksToX(hitTime, hitWireID.Plane, hitWireID.TPC, hitWireID.Cryostat);

    // Get hit Y and Z coordinates, based on central position of wire
    TVector3 xyz = theGeometry->Cryostat(hitWireID.Cryostat).TPC(hitWireID.TPC).Plane(hitWireID.Plane).Wire(hitWireID.Wire).GetCenter();

    return TVector3(xCoord, 0.f, hitType == TPC_VIEW_U ? YZToU(xyz.Y(), xyz.Z()) : hitType == TPC_VIEW_V ? YZToV(xyz.Y(), xyz.Z()) : YZToW(xyz.Y(), xyz.Z()));
}

/////////////////////////////////////////////////////////////

TVector3 pandora::VisualiseSlice::ProjectIntoPandoraView(const TVector3 &inputPosition3D, const PandoraView pandoraView)
{
    const float xCoord = inputPosition3D.X();
    const float yCoord = inputPosition3D.Y();
    const float zCoord = inputPosition3D.Z();

    return TVector3(xCoord, 0.f, pandoraView == TPC_VIEW_U ? YZToU(yCoord, zCoord) : pandoraView == TPC_VIEW_V ? YZToV(yCoord, zCoord) : YZToW(yCoord, zCoord));
}

/////////////////////////////////////////////////////////////

float pandora::VisualiseSlice::YZToU(const float yCoord, const float zCoord)
{
    const float m_uWireAngle = 1.04719758034;

    return (zCoord * std::cos(m_uWireAngle)) - (yCoord * std::sin(m_uWireAngle));
}

/////////////////////////////////////////////////////////////

float pandora::VisualiseSlice::YZToV(const float yCoord, const float zCoord)
{
    const float m_vWireAngle = -1.04719758034;

    return (zCoord * std::cos(m_vWireAngle)) - (yCoord * std::sin(m_vWireAngle));
}

/////////////////////////////////////////////////////////////

float pandora::VisualiseSlice::YZToW(const float yCoord, const float zCoord)
{
    const float m_wWireAngle = 0.0;

    return (zCoord * std::cos(m_wWireAngle)) - (yCoord * std::sin(m_wWireAngle));
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillPFPHitInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_FlashMatchModuleLabel, pfpHandle))
        throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<recob::SpacePoint> spacePointAssoc = art::FindManyP<recob::SpacePoint>(pfpHandle, evt, m_FlashMatchModuleLabel);
    const std::vector<art::Ptr<recob::SpacePoint>> &pfpSpacePoints = spacePointAssoc.at(pfp.key());

    m_nHits3D = pfpSpacePoints.size();

    std::vector<art::Ptr<recob::Hit>> pfpHits;
    CollectHitsFromClusters(evt, pfp, pfpHits);

    int uCount(0), vCount(0), wCount(0);

    for (const art::Ptr<recob::Hit> &hit : pfpHits)
    {
        // Get hit view
        const geo::WireID hit_WireID(hit->WireID());
        const geo::View_t hit_View(hit->View());
        const geo::View_t pandoraView(lar_pandora::LArPandoraGeometry::GetGlobalView(hit_WireID.Cryostat, hit_WireID.TPC, hit_View));

        if (pandoraView == geo::kW || pandoraView == geo::kY)
            ++wCount;
        else if (pandoraView == geo::kU)
            ++uCount;
        else if (pandoraView == geo::kV)
            ++vCount;
    }

    m_nHits2D = pfpHits.size();
    m_nHitsU = uCount;
    m_nHitsV = vCount;
    m_nHitsW = wCount;
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillPFPNuVertexInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    ////////////////////////
    // NuVertex Separation
    ////////////////////////
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_FlashMatchModuleLabel, pfpHandle))
        throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<recob::Vertex> vertexAssoc = art::FindManyP<recob::Vertex>(pfpHandle, evt, m_FlashMatchModuleLabel);

    const std::vector<art::Ptr<recob::Vertex>> &pfpVertex(vertexAssoc.at(pfp.key()));

    if (pfpVertex.empty())
    {
        std::cout << "pfp has no vertex!?" << std::endl;
        return;
    }

    const double dX(m_recoNuVertexX - pfpVertex.at(0)->position().X());
    const double dY(m_recoNuVertexY - pfpVertex.at(0)->position().Y());
    const double dZ(m_recoNuVertexZ - pfpVertex.at(0)->position().Z());

    m_nuVertexSeparation = std::sqrt((dX * dX) + (dY * dY) + (dZ * dZ));

    ////////////////////////
    // NuVertex Charge Distribution (in 3D)
    ////////////////////////
    const geo::Vector_t geoNuVertex(m_recoNuVertexX, m_recoNuVertexY, m_recoNuVertexZ);
    const geo::Vector_t geoPFPVertex(pfpVertex.at(0)->position().X(), pfpVertex.at(0)->position().Y(), pfpVertex.at(0)->position().Z());
    const geo::Vector_t nuVertexAxis((geoPFPVertex - geoNuVertex).Unit());

    // Need SpacePoints
    art::FindManyP<recob::SpacePoint> spacePointAssoc = art::FindManyP<recob::SpacePoint>(pfpHandle, evt, m_FlashMatchModuleLabel);
    const std::vector<art::Ptr<recob::SpacePoint>> &spacePoints(spacePointAssoc.at(pfp.key()));

    if (!spacePoints.empty())
    {
        // Need charge info
        art::Handle<std::vector<recob::SpacePoint>> spacePointHandle;

        if (!evt.getByLabel(m_FlashMatchModuleLabel, spacePointHandle))
            throw cet::exception("PhotonBDTNtuple") << "No SpacePoint Data Products Found!" << std::endl;

        art::FindManyP<recob::Hit> hitAssoc = art::FindManyP<recob::Hit>(spacePointHandle, evt, m_FlashMatchModuleLabel);

        // Nu vertex charge distribution
        double totalCharge(0.f);
        double nuVertexChargeDistribution = 0.f;

        for (const art::Ptr<recob::SpacePoint> &spacePoint : spacePoints)
        {
            const std::vector<art::Ptr<recob::Hit>> &hit(hitAssoc.at(spacePoint.key()));

            if (hit.empty())
                continue;

            const double hitCharge(std::fabs(hit.at(0)->Integral()));
            totalCharge += hitCharge;

            const geo::Vector_t geoSpacePoint(spacePoint->XYZ()[0], spacePoint->XYZ()[1], spacePoint->XYZ()[2]);
            const geo::Vector_t displacement(geoSpacePoint - geoNuVertex);
            const double thisT(std::sqrt(nuVertexAxis.Cross(displacement).Mag2()));

            nuVertexChargeDistribution += (thisT * hitCharge);
        }

        if (totalCharge > std::numeric_limits<double>::epsilon())
            nuVertexChargeDistribution /= totalCharge;

        m_nuVertexChargeDistribution = nuVertexChargeDistribution;
    }

    ////////////////////////
    // Distance of closest approach
    ////////////////////////
    art::FindManyP<recob::Shower> showerAssoc = art::FindManyP<recob::Shower>(pfpHandle, evt, m_ShowerModuleLabel);
    const std::vector<art::Ptr<recob::Shower>> &showers = showerAssoc.at(pfp.key());

    if (!showers.empty())
    {
        const art::Ptr<recob::Shower> &shower = showers.at(0);

        const TVector3 nuVertexPosition(m_recoNuVertexX, m_recoNuVertexY, m_recoNuVertexZ);
        const double alpha = std::fabs((shower->ShowerStart() - nuVertexPosition).Dot(shower->Direction()));
        const TVector3 r = shower->ShowerStart() - (alpha * shower->Direction());
        m_dca = (r - nuVertexPosition).Mag();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillPFPMiscInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    ////////////////////////
    // Generation
    ////////////////////////
    m_generation = lar_pandora::LArPandoraHelper::GetGeneration(m_pfpMap, pfp);
    m_pandoraPFPCode = pfp->PdgCode();

    ////////////////////////
    // TrackShower Score
    ////////////////////////
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_FlashMatchModuleLabel, pfpHandle))
        throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssn = art::FindManyP<larpandoraobj::PFParticleMetadata>(pfpHandle, evt, m_FlashMatchModuleLabel);
    std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadata = metadataAssn.at(pfp.key());

    if (!pfpMetadata.empty() && (pfpMetadata[0]->GetPropertiesMap().find("TrackScore") != pfpMetadata[0]->GetPropertiesMap().end()))
        m_trackShowerScore = pfpMetadata[0]->GetPropertiesMap().at("TrackScore");

    ////////////////////////
    // Parent Separation
    ////////////////////////
    if (m_generation > 2)
    {
        const int parentID(pfp->Parent());
        art::FindManyP<recob::SpacePoint> spacePointAssoc = art::FindManyP<recob::SpacePoint>(pfpHandle, evt, m_FlashMatchModuleLabel);

        if (m_pfpMap.find(parentID) != m_pfpMap.end())
        {
            const art::Ptr<recob::PFParticle> parentPFP(m_pfpMap.at(parentID));
            const std::vector<art::Ptr<recob::SpacePoint>> &pfpSpacePoints(spacePointAssoc.at(pfp.key()));
            const std::vector<art::Ptr<recob::SpacePoint>> &parentSpacePoints(spacePointAssoc.at(parentPFP.key()));

            double closestDistanceSq(std::numeric_limits<double>::max());

            for (const art::Ptr<recob::SpacePoint> &pfpSpacePoint : pfpSpacePoints)
            {
                for (const art::Ptr<recob::SpacePoint> &parentSpacePoint : parentSpacePoints)
                {
                    const double dX = pfpSpacePoint->XYZ()[0] - parentSpacePoint->XYZ()[0];
                    const double dY = pfpSpacePoint->XYZ()[1] - parentSpacePoint->XYZ()[1];
                    const double dZ = pfpSpacePoint->XYZ()[2] - parentSpacePoint->XYZ()[2];
                    const double thisSep((dX * dX) + (dY * dY) + (dZ * dZ));

                    if (thisSep < closestDistanceSq)
                        closestDistanceSq = thisSep;
                }
            }

            m_trackParentSeparation = std::sqrt(closestDistanceSq);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillPFPShowerInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_FlashMatchModuleLabel, pfpHandle))
        throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<recob::Shower> showerAssoc = art::FindManyP<recob::Shower>(pfpHandle, evt, m_ShowerModuleLabel);
    const std::vector<art::Ptr<recob::Shower>> &shower = showerAssoc.at(pfp.key());

    if (shower.empty())
        return;

    m_showerOpeningAngle = shower.at(0)->OpenAngle() * 180.0 / 3.14;
    m_showerLength = shower.at(0)->Length();
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillPFPEnergyInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_FlashMatchModuleLabel, pfpHandle))
        throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<recob::Track> trackAssoc = art::FindManyP<recob::Track>(pfpHandle, evt, m_TrackModuleLabel);

    // Get the track
    const std::vector<art::Ptr<recob::Track>> &track = trackAssoc.at(pfp.key());

    if (track.empty())
        return;

    art::Handle<std::vector<recob::Track>> trackHandle;

    if (!evt.getByLabel(m_TrackModuleLabel, trackHandle))
        throw cet::exception("VisualiseSlice") << "No Track Data Products Found!" << std::endl;

    art::FindManyP<anab::Calorimetry> caloAssoc = art::FindManyP<anab::Calorimetry>(trackHandle, evt, m_CalorimetryModuleLabel);

    const std::vector<art::Ptr<anab::Calorimetry>> &trackCalo = caloAssoc.at(track.at(0).key());

    if (trackCalo.empty())
        return;

    // Get the correct plane (collection) because it's not always the same                                                                                                        
    bool foundCorrectPlane = false;
    size_t index = 0;

    for (size_t i = 0; i < trackCalo.size(); ++i)
    {
        if (trackCalo.at(i)->PlaneID().Plane == 2)
        {
            foundCorrectPlane = true;
            index = i;
            break;
        }
    }

    if (!foundCorrectPlane)
        return;

    std::vector<float> calibrated_dEdX = trackCalo.at(index)->dEdx();
    const std::vector<float> residualRange = trackCalo.at(index)->ResidualRange();
    const std::vector<geo::Point_t> theXYZPositions = trackCalo.at(index)->XYZ();

    float minRange = std::numeric_limits<float>::max();
    float maxRange = -std::numeric_limits<float>::max();

    for (size_t i = 0; i < calibrated_dEdX.size(); ++i)
    {
        // Sometimes I think the dEdx fails? sometimes, everybody cries...
        // I think you can pick this out with a magnitude of zero
        const geo::Point_t xyzPosition = theXYZPositions[i];
        double magnitude = sqrt((xyzPosition.X() * xyzPosition.X()) + (xyzPosition.Y() * xyzPosition.Y()) + (xyzPosition.Z() * xyzPosition.Z()));

         if (magnitude < std::numeric_limits<double>::epsilon())
             continue;

         if (residualRange[i] < minRange)
             minRange = residualRange[i];

         if (residualRange[i] > maxRange)
             maxRange = residualRange[i];
    }

    // Median energy in first 4cm
    std::vector<float> initialSegmentdEdx;

    for (size_t i = 0; i < calibrated_dEdX.size(); ++i)
    {
        if (((maxRange - residualRange[i]) > 0) && ((maxRange - residualRange[i]) < 4.0))
            initialSegmentdEdx.push_back(calibrated_dEdX[i]);
    }

    m_initialdEdx = GetMedianValue(initialSegmentdEdx);


    // SUM THE ENERGY OF THE HITS!!!!! - I have no idea how uboone do this..
    m_totalEnergy = trackCalo.at(index)->KineticEnergy();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double pandora::VisualiseSlice::GetMedianValue(const std::vector<float> &inputVector)
{
    if (inputVector.empty())
        return -999.0;

    // if odd
    if (inputVector.size() % 2 != 0)
    {
        const int index = std::floor(static_cast<double>(inputVector.size()) / 2.0);
        return inputVector.at(index);
    }

    // if even
    const int firstIndex = inputVector.size() / static_cast<int>(2);
    const int secondIndex = firstIndex - 1;

    return (inputVector.at(firstIndex) + inputVector.at(secondIndex)) / 2.0;
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillPIDInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_FlashMatchModuleLabel, pfpHandle))
        throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<recob::Track> trackAssoc = art::FindManyP<recob::Track>(pfpHandle, evt, m_TrackModuleLabel);

    // Get the track
    const std::vector<art::Ptr<recob::Track>> &track = trackAssoc.at(pfp.key());

    if (track.empty())
        return;

    // Now get the PID
    art::Handle<std::vector<recob::Track>> trackHandle;

    if (!evt.getByLabel(m_TrackModuleLabel, trackHandle))
        throw cet::exception("VisualiseSlice") << "No Track Data Products Found!" << std::endl;

    art::FindManyP<anab::ParticleID> pidAssoc = art::FindManyP<anab::ParticleID>(trackHandle, evt, m_PIDModuleLabel);

    const std::vector<art::Ptr<anab::ParticleID>> &pids = pidAssoc.at(track.at(0).key());

    if (pids.empty())
        return;

    const art::Ptr<anab::ParticleID> &pid = pids.at(0);

    // Look at the chi2 PID for induction plane
    m_chiPIDProton = searchingfornues::PID(pid, "Chi2", anab::kGOF, anab::kForward, 2212, 2);
    m_chiPIDMuon = searchingfornues::PID(pid, "Chi2", anab::kGOF, anab::kForward, 13, 2);
    m_chiPIDPion = searchingfornues::PID(pid, "Chi2", anab::kGOF, anab::kForward, 211, 2);

    // Look at Bragg peak for induction plane
    m_braggPIDProton = std::max(searchingfornues::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kForward, 2212, 2),
                                searchingfornues::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 2212, 2));
    m_braggPIDMuon = std::max(searchingfornues::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kForward, 13, 2),
                              searchingfornues::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 13, 2));
    m_braggPIDPion = std::max(searchingfornues::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kForward, 211, 2),
                              searchingfornues::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 211, 2));

    // PID LLR calculator
    art::FindManyP<anab::Calorimetry> caloAssoc = art::FindManyP<anab::Calorimetry>(trackHandle, evt, m_CalorimetryModuleLabel);
    const std::vector<art::Ptr<anab::Calorimetry>> &trackCalo = caloAssoc.at(track.at(0).key());

    if (trackCalo.empty())
        return;

    float LLRPID = 0;

    for (const art::Ptr<anab::Calorimetry> &calo : trackCalo)
    {
        if (calo->ResidualRange().size() == 0) continue;

        auto const &plane = calo->PlaneID().Plane;
        auto const &dEdxValues = calo->dEdx();
        auto const &resRange = calo->ResidualRange();
        auto const &pitch = calo->TrkPitchVec();
        std::vector<std::vector<float>> paramValues;
        paramValues.push_back(resRange);
        paramValues.push_back(pitch);

        float caloEnergy = 0;

        for (unsigned int i = 0; i < dEdxValues.size(); i++)
            caloEnergy += dEdxValues[i] * pitch[i];

        const float thisLLRPID = m_LLRPIDCalculator.LLR_many_hits_one_plane(dEdxValues, paramValues, plane);

        LLRPID += thisLLRPID;
    }

    m_LLRPIDReduced = atan(LLRPID / 100.f) * 2.f / 3.14159266;  
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::beginJob()
{
    if (m_debug) std::cout << "Starting..." << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::endJob()
{
    if (m_debug) std::cout << "Finished..." << std::endl;
}

DEFINE_ART_MODULE(pandora::VisualiseSlice)
