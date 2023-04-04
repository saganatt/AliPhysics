#ifndef ALIANALYSISTASKTRACKINGEFFPID
#define ALIANALYSISTASKTRACKINGEFFPID

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskTrackingEffPID
// AliAnalysisTaskSE to compute tracking and PID efficiencies for 
//  different particle species
//
// Authors:
//          M. Puccio
//          F. Prino
//          
//*************************************************************************


#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"
#include "AliPID.h"
#include "AliEventCuts.h"
#include <THnSparse.h>

class TList;

class AliAnalysisTaskTrackingEffPID : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskTrackingEffPID();
  virtual ~AliAnalysisTaskTrackingEffPID();

  static bool  HasTOF(AliVTrack *t);
  bool ConvertAndSelectAODTrack(AliAODTrack* aTrack, const AliESDVertex vESD, Double_t magField);

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);

  void SetTrackCuts(AliESDtrackCuts* cuts){
    if(fTrackCuts) delete fTrackCuts;
    fTrackCuts = new AliESDtrackCuts(*cuts);
    float c1, c2, c3, c4, c5;
    fTrackCuts->GetMaxCovDiagonalElements(c1, c2, c3, c4, c5);
    float pmin, pmax, ptmin, ptmax, pxmin, pxmax, pymin, pymax, pzmin, pzmax;
    float etamin, etamax, rapmin, rapmax;
    fTrackCuts->GetPRange(pmin, pmax);
    fTrackCuts->GetPtRange(ptmin, ptmax);
    fTrackCuts->GetPxRange(pxmin, pxmax);
    fTrackCuts->GetPyRange(pymin, pymax);
    fTrackCuts->GetPzRange(pzmin, pzmax);
    fTrackCuts->GetEtaRange(etamin, etamax);
    fTrackCuts->GetRapRange(rapmin, rapmax);
    std::cout << "Setting new track cuts " << std::endl;
    std::cout << " filter mask: " << fTrackCuts->GetFilterMask() << std::endl;
    std::cout << " selected: " << fTrackCuts->Selected() << std::endl;
    std::cout << " fCutMinNClusterTPC: " << fTrackCuts->GetMinNClusterTPC() << " ITS: " << fTrackCuts->GetMinNClustersITS() << std::endl;
    std::cout << " fCutMinNCrossedRowsTPC: " << fTrackCuts->GetMinNCrossedRowsTPC() << std::endl;
    std::cout << " min ratio crossed rows over findable clusters TPC: " << fTrackCuts->GetMinRatioCrossedRowsOverFindableClustersTPC() << std::endl;
    //std::cout << " fCutMaxPtDepNClustersTPC: " << fTrackCuts->GetMaxPtDepNClustersTPC() << std::endl;
    std::cout << " fCutMinLengthActiveVolumeTPC: " << fTrackCuts->GetMinLengthActiveVolumeTPC() << std::endl;
    //std::cout << " fDeadZoneWidth: " << fTrackCuts->GetDeadZoneWidth() << std::endl;
    //std::cout << " fCutGeoNcrNclLength: " << fTrackCuts->GetGeoNcrNclLength() << std::endl;
    //std::cout << " fCutGeoNcrNclGeom1Pt: " << fTrackCuts->GetGeoNcrNclGeom1Pt() << std::endl;
    //std::cout << " fCutGeoNcrNclFractionNcr: " << fTrackCuts->GetGeoNcrNclFractionNcr() << std::endl;
    //std::cout << " fCutGeoNcrNclFractionNcl: " << fTrackCuts->GetGeoNcrNclFractionNcl() << std::endl;
    //std::cout << " fCutOutDistortedRegionTPC: " << fTrackCuts->GetOutDistortedRegionTPC() << std::endl;
    std::cout << " fCutMaxChi2PerClusterTPC: " << fTrackCuts->GetMaxChi2PerClusterTPC() << std::endl;
    std::cout << " fCutMaxChi2PerClusterITS: " << fTrackCuts->GetMaxChi2PerClusterITS() << std::endl;
    std::cout << " fCutMaxChi2TPCConstrainedVsGlobal: " << fTrackCuts->GetMaxChi2TPCConstrainedGlobal() << std::endl;
    std::cout << " fCutMaxChi2TPCConstrainedVsGlobalVertexType: " << fTrackCuts->GetMaxChi2TPCConstrainedGlobalVertexType() << std::endl;
    std::cout << " fCutMaxMissingITSPoints: " << fTrackCuts->GetMaxNOfMissingITSPoints() << std::endl;
    //std::cout << " fCutClusterRequirementITS: " << fTrackCuts->GetClusterRequirementITS(kSPD) << " " << fTrackCuts->GetClusterRequirementITS(kSDD) << " " << fTrackCuts->GetClusterRequirementITS(kSSD) << std::endl;
    std::cout << " fCutMax cov diag elements: " << c1 << " " << c2 << " " << c3 << " " << c4 << " " << c5 << std::endl;
    std::cout << " fCutMaxRel1PtUncertainty: " << fTrackCuts->GetMaxRel1PtUncertainty() << std::endl;
    std::cout << " fCutMaxRel1PtUncertaintyPtDep: " << fTrackCuts->GetMaxRel1PtUncertaintyPtDep() << std::endl;
    std::cout << " fCutAcceptKinkDaughters: " << fTrackCuts->GetAcceptKinkDaughters() << std::endl;
    std::cout << " fCutAcceptSharedTPCClusters: " << fTrackCuts->GetAcceptSharedTPCClusters() << std::endl;
    std::cout << " fCutMaxFractionSharedTPCClusters: " << fTrackCuts->GetMaxFractionSharedTPCClusters() << std::endl;
    std::cout << " fCutRequireTPCRefit: " << fTrackCuts->GetRequireTPCRefit() << std::endl;
    std::cout << " fCutRequireTPCStandAlone: " << fTrackCuts->GetRequireTPCStandAlone() << std::endl;
    std::cout << " fCutRequireITSRefit: " << fTrackCuts->GetRequireITSRefit() << std::endl;
    //std::cout << " fCutRequireITSPid: " << fTrackCuts->GetRequireITSPid() << std::endl;
    std::cout << " fCutRequireITSStandAlone: " << fTrackCuts->GetRequireITSStandAlone() << std::endl;
    std::cout << " fCutRequireITSpureSA: " << fTrackCuts->GetRequireITSpureSA() << std::endl;
    std::cout << " fCutNsigmaToVertex: " << fTrackCuts->GetMaxNsigmaToVertex() << std::endl;
    std::cout << " fCutSigmaToVertexRequired: " << fTrackCuts->GetRequireSigmaToVertex() << std::endl;
    std::cout << " fCutMaxDCAToVertexXY: " << fTrackCuts->GetMaxDCAToVertexXY() << std::endl;
    std::cout << " fCutMaxDCAToVertexZ: " << fTrackCuts->GetMaxDCAToVertexZ() << std::endl;
    std::cout << " fCutDCAToVertex2D: " << fTrackCuts->GetDCAToVertex2D() << std::endl;
    std::cout << " fCutMinDCAToVertexXY: " << fTrackCuts->GetMinDCAToVertexXY() << std::endl;
    std::cout << " fCutMinDCAToVertexZ: " << fTrackCuts->GetMinDCAToVertexZ() << std::endl;
    std::cout << " fCutMaxDCAToVertexXYPtDep: " << fTrackCuts->GetMaxDCAToVertexXYPtDep() << std::endl;
    std::cout << " fCutMaxDCAToVertexZPtDep: " << fTrackCuts->GetMaxDCAToVertexZPtDep() << std::endl;
    std::cout << " fCutMinDCAToVertexXYPtDep: " << fTrackCuts->GetMinDCAToVertexXYPtDep() << std::endl;
    std::cout << " fCutMinDCAToVertexZPtDep: " << fTrackCuts->GetMinDCAToVertexZPtDep() << std::endl;
    std::cout << " fPMin, fPMax: " << pmin << ", " << pmax << std::endl;
    std::cout << " fPtMin, fPtMax: " << ptmin << ", " << ptmax << std::endl;
    std::cout << " fPxMin, fPxMax: " << pxmin << ", " << pxmax << std::endl;
    std::cout << " fPyMin, fPyMax: " << pymin << ", " << pymax << std::endl;
    std::cout << " fPzMin, fPzMax: " << pzmin << ", " << pzmax << std::endl;
    std::cout << " fEtaMin, fEtaMax: " << etamin << ", " << etamax << std::endl;
    std::cout << " fRapMin, fRapMax: " << rapmin << ", " << rapmax << std::endl;
    std::cout << " fFlagCutTOFdistance: " << fTrackCuts->GetFlagCutTOFdistance() << std::endl;
    std::cout << " fCutTOFdistance: " << fTrackCuts->GetCutTOFdistance() << std::endl;
    //std::cout << " fCutRequireTOFout: " << fTrackCuts->GetRequireTOFout() << std::endl;
  }
  AliESDtrackCuts* GetTrackCuts() const{
    return fTrackCuts;
  }
  void SetPrimarySelectionOption(Int_t opt){
    fPrimarySelectionOpt=opt;
  }
  void SetTrackletMultiplicityEstimator(){ fMultEstimator=0;}
  void SetVertexContribMultiplicityEstimator(){ fMultEstimator=1;}
  void SetTracksMultiplicityEstimator(){ fMultEstimator=2;}
  void SetTPCTracksMultiplicityEstimator(){ fMultEstimator=3;}
  void SetTPCClustersMultiplicityEstimator(){ fMultEstimator=4;}

  void SetCollisionSystem(TString collsy){
    collsy.ToLower();
    collsy.ReplaceAll("-","");
    if(collsy=="pbpb" || collsy=="xexe") fIsAA=kTRUE;
    else fIsAA=kFALSE;
  }
  void SetFilterBitCutForAODTracks(int fb){
    fFilterBit=fb;
  }
  void UseOnlyFilterBitCutForAODTracks(){
    fUseTrackCutsForAOD=kFALSE;    
  }
  void UseTrackCutObjectForAODTracks(){
    fUseTrackCutsForAOD=kTRUE;    
  }

  void SetUseGeneratedKine(bool flag) {fUseGeneratedKine=flag;}
  void RejectEventsWithSameBunchPileup(bool flag) {fRejectEventsWithSameBunchPileup=flag;}
  void RejectGeneratedParticlesFromPileup(bool flag) {fRejectPileupParticles=flag;}
  void RejectTracksOfPileupParticles(bool flag) {fRejectTracksOfPileupPart=flag;}
  
  void SetTriggerMask(ULong64_t mask){
    fEventCut.OverrideAutomaticTriggerSelection(mask);
  }
  void SetUseSPDPileup(bool multDep, int nContrCut=5, double dzCut=0.8){
    fEventCut.OverridePileUpCuts(nContrCut,dzCut,3.,2.,5.,kTRUE);
    fEventCut.fUseMultiplicityDependentPileUpCuts=multDep;
  }
  void SetUseMVPileup(bool flag) {fEventCut.fPileUpCutMV=flag;}

  void SetUseImpactParameter(bool flag) {fUseImpPar=flag;}
  void SetUseLocalTrackDensity(bool flag, double deltaRcut=0.2, double maxNtracks=-1.) {
    fUseLocDen=flag; fDeltaRcut=deltaRcut; fMaxTracksInCone=maxNtracks;
  }
  void SetPtHardRange(double pmin, double pmax){
    fSelectPtHardRange=kTRUE; fMinPtHard=pmin; fMaxPtHard=pmax;
  }
  void SelectedGeneratorName(TString name){
    fGenerToKeep=name.Data(); fSelectOnGenerator=kTRUE;}
  void ExcludedGeneratorName(TString name){
    fGenerToExclude=name.Data(); fSelectOnGenerator=kTRUE;}
  void KeepOnlyInjectedParticles(bool opt){fKeepOnlyInjected=opt;}
  void KeepOnlyUnderlyingEventParticles(bool opt){fKeepOnlyUE=opt;}
  TString GetGenerator(int label, TList *lh);
  bool IsInjectedParticle(int lab, TList *lh);
  double GetLocalTrackDens(TNtuple* trEtaPhiMap, double eta, double phi) const;
  AliEventCuts  fEventCut;


private:
  AliAnalysisTaskTrackingEffPID (const AliAnalysisTaskTrackingEffPID &source);
  AliAnalysisTaskTrackingEffPID &operator=(const AliAnalysisTaskTrackingEffPID &source);

  bool fUseTrackCutsForAOD;       /// flag to switch off/on fTrackCuts for AOD
  bool fUseGeneratedKine;         /// flag to use the generated pt, eta phi
  bool fRejectEventsWithSameBunchPileup; /// flag to reject events with generated same-bunch pileup
  bool fRejectPileupParticles;    /// flag to reject from generated particles those from pileup
  bool fRejectTracksOfPileupPart; /// flag to reject from reco particles those from pileup
  int  fPrimarySelectionOpt;      /// 0=no selection, 1=IsPhysicalPrimary, 2= cut on the origin
  int  fMultEstimator;            /// multiplicity estimator: 0=trackelts, 1=ITS+TPCtracks, 2=primary vertex contributors
  bool fIsAA;                     /// flag to control collision system
  int  fFilterBit;                /// filter-bit selection for AOD tracks
  AliESDtrackCuts* fTrackCuts;                            /// cut object
  bool fSelectOnGenerator;       /// flag to select events with generator name
  TString fGenerToKeep;          /// generator name to analyse
  TString fGenerToExclude;       /// generator name to exclude
  bool fKeepOnlyInjected;        /// flag to keep only injected particles
  bool fKeepOnlyUE;              /// flag to keep only underlying event
  bool fUseImpPar;               /// flag to enable plots vs. impact parameter
  bool fUseLocDen;               /// flag to enable plots vs. local track density
  double fDeltaRcut;             /// radius cut to count local track density
  double fMaxTracksInCone;       /// upper limit for track density axis
  bool fSelectPtHardRange;       /// flag to enable the cut on pthard
  double fMinPtHard;             /// min pthard
  double fMaxPtHard;             /// max pthard
  TList* fOutputList;                                     //!<! Output list
  TList* fListCuts;                                       //!<! Output with cuts
  TH1F*  fHistNEvents;                                    //!<!  histo with N of events  
  TH1D*  fHistNParticles;                                 //!<!  histo with N of particles
  TH1D*  fHistNTracks;                                    //!<!  histo with N of tracks
  TH2D*  fHistPileupTagAOD;                               //!<!  control plot for AOD pileup
  TH1D*  hHistXsecVsPtHard;                               //!<!  control plot
  THnSparseF* fGenerated[AliPID::kSPECIESC][2];           //!<! Generated particles (pt, eta, phi, mult, zvert)
  THnSparseF* fGeneratedEvSel[AliPID::kSPECIESC][2];      //!<! Generated particles after event selection
  THnSparseF* fReconstructed[AliPID::kSPECIESC][2];       //!<! Reconstructed particles (pt, eta, phi, mult, zvert)
  THnSparseF* fReconstructedTOF[AliPID::kSPECIESC][2];    //!<! Reconstructed particles after PID (pt, eta, phi, mult, zvert)
  THnSparseF* fReconstructedPID[AliPID::kSPECIESC][2];    //!<! Reconstructed particles after PID (pt, eta, phi, mult, zvert)
  TH1D* fRawPt;
  TH1D* fRawEta;
  TH1D* fRawPhi;

  bool collisionProcessed;


  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskTrackingEffPID, 12);
  /// \endcond
};


#endif 
