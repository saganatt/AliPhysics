/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, proviyaded that the above copyright notice appears in all *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purapose. It is         *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                                                                       //
//          Analysis for event-by-event particle ratio studies           //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TROOT.h"
#include "TGrid.h"
#include "TSystem.h"
#include "TCutG.h"
#include "TH1F.h"
#include "THn.h"
#include "THnSparse.h"
#include "TList.h"
#include "TMath.h"
#include "TMatrixF.h"
#include "TVectorF.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TRandom.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include "AliPDG.h"
#include "AliMathBase.h"
#include "AliESDFMD.h"
#include "AliFMDFloatMap.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliTPCdEdxInfo.h"
#include "AliESDv0KineCuts.h"
#include "AliKFVertex.h"
#include "AliLumiTools.h"
#include "AliKFParticle.h"
#include "AliCollisionGeometry.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenHepMCEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEposEventHeader.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "AliPID.h"
#include "AliESDtrackCuts.h"
#include "AliESDv0Cuts.h"
#include "AliESDv0.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliCentrality.h"
#include "AliESDUtils.h"
#include "AliMultiplicity.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliTPCParam.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliKFParticle.h"
#include "AliAnalysisTaskFilteredTree.h"
#include "AliAnalysisTaskEbyeIterPID.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"
#include "AliRunLoader.h"
#include "AliEventCuts.h"
#include "AliAnalysisUtils.h"
#include "AliESDtools.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <bitset>
using namespace std;
using std::cout;
using std::setw;



ClassImp(AliAnalysisTaskEbyeIterPID)

const char* AliAnalysisTaskEbyeIterPID::fEventInfo_centEstStr[] = {"V0M","CL0","CL1"};
// {"V0M","V0A","V0C","FMD","TRK","TKL","CL0","CL1","V0MvsFMD","ZNA","TKLvsV0M","ZEMvsZDC","V0A123","V0A0","V0S", "MB", "Ref", "V0av"};
// const char* AliAnalysisTaskEbyeIterPID::fEventInfo_centEstStr[] = {"V0M","CL0","CL1","TRK","TKL","V0MvsFMD","TKLvsV0M","ZEMvsZDC","kV0A"
// ,"V0C","ZNA","ZNC","ZPA","ZPC","CND","FMD","NPA","V0A0","V0A123","V0A23","V0C01","V0S","V0MEq","V0AEq","V0CEq","SPDClusters","SPDTracklets"};

#define USE_STREAMER 1


// -----------------------------------------------------------------------
//                            Constructor and Destructor
// -----------------------------------------------------------------------
//________________________________________________________________________
AliAnalysisTaskEbyeIterPID::AliAnalysisTaskEbyeIterPID()
: AliAnalysisTaskSE("TaskEbyeRatios"), fEventCuts(0), fPIDResponse(0),fESD(0), fListHist(0),
fESDtrackCuts(0),
fESDtrackCuts_Bit96(0),
fESDtrackCuts_Bit96_spd(0),
fESDtrackCuts_Bit96_sdd(0),
fESDtrackCuts_Bit128(0),
fESDtrackCuts_Bit768(0),
fESDtrackCutsLoose(0),
fESDtrackCutsV0(0),
fESDtrackCutsCleanSamp(0),
fPIDCombined(0x0),
fTPCdEdxInfo(0x0),
fMCStack(0x0),
fV0OpenCuts(0x0),
fV0StrongCuts(0x0),
fK0sPionCuts(0x0),
fLambdaProtonCuts(0x0),
fLambdaPionCuts(0x0),
fGammaElectronCuts(0x0),
fVertex(0x0),
fESDtool(nullptr),
fArmPodTree(0x0),
fTreeSRedirector(0x0),
fTreeMCFull(0x0),
fTreeMCgen(0x0),
fTreeDnchDeta(0x0),
fTreeMC(0x0),
fTreedEdxCheck(0x0),
fTreeCuts(0x0),
fTreeMCFullAcc(0x0),
fTreeResonance(0x0),
fTreeMCgenMoms(0x0),
fTreeEvents(0x0),
fTreeDScaled(0x0),
fTreeMCEffCorr(0x0),
fRandom(0),
fPeriodName(""),
fYear(0),
fPassIndex(0),
fPileUpBit(0),
fHistCent(0),
fHistPhi(0),
fHistGenMult(0),
fHistRapDistFullAccPr(0),
fHistRapDistFullAccAPr(0),
fHistInvK0s(0),
fHistInvLambda(0),
fHistInvAntiLambda(0),
fHistInvPhoton(0),
fHistPhiTPCcounterA(0),
fHistPhiTPCcounterC(0),
fHistPhiTPCcounterAITS(0),
fHistPhiTPCcounterCITS(0),
fHistPhiITScounterA(0),
fHistPhiITScounterC(0),
fHndEdx(),
fHnExpected(),
fH2MissCl(),
fChunkName(""),
fTrackCutBits(0),
fSystClass(0),
fEtaDown(0),
fEtaUp(0),
fNEtaBins(0),
fPercentageOfEvents(0),
fRunOnGrid(kFALSE),
fMCtrue(kFALSE),
fEventInfo(kFALSE),
fWeakAndMaterial(kFALSE),
fEffMatrix(kFALSE),
fDEdxCheck(kFALSE),
fIncludeITS(kTRUE),
fFillTracks(kFALSE),
fFillOnlyHists(kFALSE),
fFillEffLookUpTable(kFALSE),
fFillHigherMomentsMCclosure(kFALSE),
fFillArmPodTree(kTRUE),
fRunFastSimulation(kFALSE),
fRunFastHighMomentCal(kFALSE),
fFillDistributions(kFALSE),
fFillTreeMC(kFALSE),
fDefaultTrackCuts(kFALSE),
fDefaultEventCuts(kFALSE),
fFillNudynFastGen(kFALSE),
fCorrectForMissCl(0),
fUsePtCut(1),
fTrackOriginOnlyPrimary(0),
fRapidityType(0),
fSisterCheck(0),
fFillDnchDeta(kFALSE),
fIncludeTOF(kFALSE),
fUseThnSparse(kFALSE),
fUseCouts(kFALSE),
fNSettings(22),
fNMomBins(0),
fMomDown(0),
fMomUp(0),
fDEdxBinWidth(0),
fDEdxUp(0),
fDEdxDown(0),
fDEdxCleanUp(0),
fArmPodTPCSignal(0),
fArmPodptot(0),
fArmPodEta(0),
fArmPodCentrality(0),
fQt(0),
fAlfa(0),
fNSigmasElTOF(0),
fNSigmasPiTOF(0),
fNSigmasKaTOF(0),
fNSigmasPrTOF(0),
fNSigmasDeTOF(0),
fDEdxEl(0),
fDEdxKa(0),
fDEdxPi(0),
fDEdxPr(0),
fDEdxDe(0),
fSigmaEl(0),
fSigmaKa(0),
fSigmaPi(0),
fSigmaPr(0),
fSigmaDe(0),
fNSigmasElTPC(0),
fNSigmasPiTPC(0),
fNSigmasKaTPC(0),
fNSigmasPrTPC(0),
fNSigmasDeTPC(0),
fTPCSignalMC(0),
fPtotMC(0),
fPtotMCtruth(0),
fPtMC(0),
fEtaMC(0),
fSignMC(0),
fPxMC(0),
fPyMC(0),
fPzMC(0),
fElMC(0),
fPiMC(0),
fKaMC(0),
fPrMC(0),
fDeMC(0),
fMuMC(0),
fLaMC(0),
fMCImpactParameter(0),
fNHardScatters(0),
fNProjectileParticipants(0),
fNTargetParticipants(0),
fNNColl(0),
fNNwColl(0),
fNwNColl(0),
fNwNwColl(0),
fElMCgen(0),
fPiMCgen(0),
fKaMCgen(0),
fPrMCgen(0),
fDeMCgen(0),
fMuMCgen(0),
fLaMCgen(0),
fBaMCgen(0),
fPx(0),
fPy(0),
fPz(0),
fPtot(0),
fPVertex(0),
fPt(0),
fY(0),
fMultiplicity(0),
fMultiplicityMC(0),
fCentrality(0),
fCentImpBin(0),
fVz(0),
fEventGID(0),
fEventGIDMC(0),
fEventCountInFile(0),
fEvent(0),
fEventMC(0),
fEventMCgen(0),
fTPCSignal(0),
fEta(0),
fNContributors(0),
fTheta(0),
fPhi(0),
fSign(0),
fTPCShared(0),
fNcl(0),
fNclCorr(0),
fNResBins(0),
fNBarBins(0),
fNEtaWinBinsMC(-100),
fNMomBinsMC(-100),
fNCentBinsMC(-100),
fGenprotonBins(-100),
fNResModeMC(2),
fNCentbinsData(10),
fMissingCl(0.),
fTPCMult(0),
fEventMult(0),
fTimeStamp(0),
fIntRate(0),
fRunNo(0),
fBField(0),
fBeamType(0),
fTrackProbElTPC(0),
fTrackProbPiTPC(0),
fTrackProbKaTPC(0),
fTrackProbPrTPC(0),
fTrackProbDeTPC(0),
fTrackProbElTOF(0),
fTrackProbPiTOF(0),
fTrackProbKaTOF(0),
fTrackProbPrTOF(0),
fTrackProbDeTOF(0),
fTrackTPCCrossedRows(0),
fTrackChi2TPC(0),
fTrackChi2TPCcorr(0),
fTrackDCAxy(0),
fTrackDCAz(0),
fTrackLengthInActiveZone(0),
fTrackTPCSignalN(0),
fTrackIsFirstITSlayer(0),
fTrackIsSecondITSlayer(0),
fTrackNewITScut(0),
fTrackRequireITSRefit(0),
fIsITSpixel01(0),
fNITSclusters(0),
fPrimRestriction(0),
fTPCvZ(0),
fSPDvZ(0),
fCleanPionsFromK0(0),
fCleanPion0FromK0(0),
fCleanPion1FromK0(0),
fCleanPion0FromLambda(0),
fCleanPion1FromLambda(0),
fCleanProton0FromLambda(0),
fCleanProton1FromLambda(0),
fHasTrack0FirstITSlayer(0),
fHasTrack1FirstITSlayer(0),
fHasV0FirstITSlayer(0),
fSystCentEstimatetor(0),
fSystCrossedRows(0),
fSystDCAxy(0),
fSystChi2(0),
fSystVz(0),
fetaDownArr(),
fetaUpArr(),
fcentDownArr(),
fcentUpArr(),
fpDownArr(),
fpUpArr(),
fxCentBins(),
fResonances(),
fBaryons(),
fHistPosEffMatrixRec(0),
fHistNegEffMatrixRec(0),
fHistPosEffMatrixGen(0),
fHistNegEffMatrixGen(0),
fHistPosEffMatrixScanRec(0),
fHistNegEffMatrixScanRec(0),
fHistPosEffMatrixScanGen(0),
fHistNegEffMatrixScanGen(0),
fHistEmptyEvent(0),
fHistCentrality(0),
fHistCentralityImpPar(0),
fHistImpParam(0),
fHistVertex(0),
fHistdEdxTPC(0),
fHistArmPod(0),
fEventInfo_PhiTPCdcarA(0),
fEventInfo_PhiTPCdcarC(0),
fEventInfo_CacheTrackCounters(0),
fEventInfo_CacheTrackdEdxRatio(0),
fEventInfo_CacheTrackNcl(0),
fEventInfo_CacheTrackChi2(0),
fEventInfo_CacheTrackMatchEff(0),
fEventInfo_CentralityEstimates(0),
fEventInfo_LumiGraph(0),
fEventInfo_HisTPCVertexA(0),
fEventInfo_HisTPCVertexC(0),
fEventInfo_HisTPCVertexACut(0),
fEventInfo_HisTPCVertexCCut(0),
fEventInfo_HisTPCVertex(0),
fEventInfo_CacheTrackTPCCountersZ(0),
fPileUpTightnessCut1(0),
fPileUpTightnessCut2(0),
fPileUpTightnessCut3(0),
fPileUpTightnessCut4(0)
{
  // default Constructor
  /* fast compilation test
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/AliAnalysisTaskEbyeIterPID.cxx++
  .L /lustre/nyx/alice/users/marsland/train/trunk/marsland_EbyeRatios/AliAnalysisTaskEbyeIterPID.cxx++
  */
}

//________________________________________________________________________
AliAnalysisTaskEbyeIterPID::AliAnalysisTaskEbyeIterPID(const char *name)
: AliAnalysisTaskSE(name), fEventCuts(0), fPIDResponse(0), fESD(0), fListHist(0),
fESDtrackCuts(0),
fESDtrackCuts_Bit96(0),
fESDtrackCuts_Bit96_spd(0),
fESDtrackCuts_Bit96_sdd(0),
fESDtrackCuts_Bit128(0),
fESDtrackCuts_Bit768(0),
fESDtrackCutsLoose(0),
fESDtrackCutsV0(0),
fESDtrackCutsCleanSamp(0),
fPIDCombined(0x0),
fTPCdEdxInfo(0x0),
fMCStack(0x0),
fV0OpenCuts(0x0),
fV0StrongCuts(0x0),
fK0sPionCuts(0x0),
fLambdaProtonCuts(0x0),
fLambdaPionCuts(0x0),
fGammaElectronCuts(0x0),
fVertex(0x0),
fESDtool(nullptr),
fArmPodTree(0x0),
fTreeSRedirector(0x0),
fTreeMCFull(0x0),
fTreeMCgen(0x0),
fTreeDnchDeta(0x0),
fTreeMC(0x0),
fTreedEdxCheck(0x0),
fTreeCuts(0x0),
fTreeMCFullAcc(0x0),
fTreeResonance(0x0),
fTreeMCgenMoms(0x0),
fTreeEvents(0x0),
fTreeDScaled(0x0),
fTreeMCEffCorr(0x0),
fRandom(0),
fPeriodName(""),
fYear(0),
fPassIndex(0),
fPileUpBit(0),
fHistCent(0),
fHistPhi(0),
fHistGenMult(0),
fHistRapDistFullAccPr(0),
fHistRapDistFullAccAPr(0),
fHistInvK0s(0),
fHistInvLambda(0),
fHistInvAntiLambda(0),
fHistInvPhoton(0),
fHistPhiTPCcounterA(0),
fHistPhiTPCcounterC(0),
fHistPhiTPCcounterAITS(0),
fHistPhiTPCcounterCITS(0),
fHistPhiITScounterA(0),
fHistPhiITScounterC(0),
fHndEdx(),
fHnExpected(),
fH2MissCl(),
fChunkName(""),
fTrackCutBits(0),
fSystClass(0),
fEtaDown(0),
fEtaUp(0),
fNEtaBins(0),
fPercentageOfEvents(0),
fRunOnGrid(kFALSE),
fMCtrue(kFALSE),
fEventInfo(kFALSE),
fWeakAndMaterial(kFALSE),
fEffMatrix(kFALSE),
fDEdxCheck(kFALSE),
fIncludeITS(kTRUE),
fFillTracks(kFALSE),
fFillOnlyHists(kFALSE),
fFillEffLookUpTable(kFALSE),
fFillHigherMomentsMCclosure(kFALSE),
fFillArmPodTree(kTRUE),
fRunFastSimulation(kFALSE),
fRunFastHighMomentCal(kFALSE),
fFillDistributions(kFALSE),
fFillTreeMC(kFALSE),
fDefaultTrackCuts(kFALSE),
fDefaultEventCuts(kFALSE),
fFillNudynFastGen(kFALSE),
fCorrectForMissCl(0),
fUsePtCut(1),
fTrackOriginOnlyPrimary(0),
fRapidityType(0),
fSisterCheck(0),
fFillDnchDeta(kFALSE),
fIncludeTOF(kFALSE),
fUseThnSparse(kFALSE),
fUseCouts(kFALSE),
fNSettings(22),
fNMomBins(0),
fMomDown(0),
fMomUp(0),
fDEdxBinWidth(0),
fDEdxUp(0),
fDEdxDown(0),
fDEdxCleanUp(0),
fArmPodTPCSignal(0),
fArmPodptot(0),
fArmPodEta(0),
fArmPodCentrality(0),
fQt(0),
fAlfa(0),
fNSigmasElTOF(0),
fNSigmasPiTOF(0),
fNSigmasKaTOF(0),
fNSigmasPrTOF(0),
fNSigmasDeTOF(0),
fDEdxEl(0),
fDEdxKa(0),
fDEdxPi(0),
fDEdxPr(0),
fDEdxDe(0),
fSigmaEl(0),
fSigmaKa(0),
fSigmaPi(0),
fSigmaPr(0),
fSigmaDe(0),
fNSigmasElTPC(0),
fNSigmasPiTPC(0),
fNSigmasKaTPC(0),
fNSigmasPrTPC(0),
fNSigmasDeTPC(0),
fTPCSignalMC(0),
fPtotMC(0),
fPtotMCtruth(0),
fPtMC(0),
fEtaMC(0),
fSignMC(0),
fPxMC(0),
fPyMC(0),
fPzMC(0),
fElMC(0),
fPiMC(0),
fKaMC(0),
fPrMC(0),
fDeMC(0),
fMuMC(0),
fLaMC(0),
fMCImpactParameter(0),
fNHardScatters(0),
fNProjectileParticipants(0),
fNTargetParticipants(0),
fNNColl(0),
fNNwColl(0),
fNwNColl(0),
fNwNwColl(0),
fElMCgen(0),
fPiMCgen(0),
fKaMCgen(0),
fPrMCgen(0),
fDeMCgen(0),
fMuMCgen(0),
fLaMCgen(0),
fBaMCgen(0),
fPx(0),
fPy(0),
fPz(0),
fPtot(0),
fPVertex(0),
fPt(0),
fY(0),
fMultiplicity(0),
fMultiplicityMC(0),
fCentrality(0),
fCentImpBin(0),
fVz(0),
fEventGID(0),
fEventGIDMC(0),
fEventCountInFile(0),
fEvent(0),
fEventMC(0),
fEventMCgen(0),
fTPCSignal(0),
fEta(0),
fNContributors(0),
fTheta(0),
fPhi(0),
fSign(0),
fTPCShared(0),
fNcl(0),
fNclCorr(0),
fNResBins(0),
fNBarBins(0),
fNEtaWinBinsMC(-100),
fNMomBinsMC(-100),
fNCentBinsMC(-100),
fGenprotonBins(-100),
fNResModeMC(2),
fNCentbinsData(10),
fMissingCl(0.),
fTPCMult(0),
fEventMult(0),
fTimeStamp(0),
fIntRate(0),
fRunNo(0),
fBField(0),
fBeamType(0),
fTrackProbElTPC(0),
fTrackProbPiTPC(0),
fTrackProbKaTPC(0),
fTrackProbPrTPC(0),
fTrackProbDeTPC(0),
fTrackProbElTOF(0),
fTrackProbPiTOF(0),
fTrackProbKaTOF(0),
fTrackProbPrTOF(0),
fTrackProbDeTOF(0),
fTrackTPCCrossedRows(0),
fTrackChi2TPC(0),
fTrackChi2TPCcorr(0),
fTrackDCAxy(0),
fTrackDCAz(0),
fTrackLengthInActiveZone(0),
fTrackTPCSignalN(0),
fTrackIsFirstITSlayer(0),
fTrackIsSecondITSlayer(0),
fTrackNewITScut(0),
fTrackRequireITSRefit(0),
fIsITSpixel01(0),
fNITSclusters(0),
fPrimRestriction(0),
fTPCvZ(0),
fSPDvZ(0),
fCleanPionsFromK0(0),
fCleanPion0FromK0(0),
fCleanPion1FromK0(0),
fCleanPion0FromLambda(0),
fCleanPion1FromLambda(0),
fCleanProton0FromLambda(0),
fCleanProton1FromLambda(0),
fHasTrack0FirstITSlayer(0),
fHasTrack1FirstITSlayer(0),
fHasV0FirstITSlayer(0),
fSystCentEstimatetor(0),
fSystCrossedRows(0),
fSystDCAxy(0),
fSystChi2(0),
fSystVz(0),
fetaDownArr(),
fetaUpArr(),
fcentDownArr(),
fcentUpArr(),
fpDownArr(),
fpUpArr(),
fxCentBins(),
fResonances(),
fBaryons(),
fHistPosEffMatrixRec(0),
fHistNegEffMatrixRec(0),
fHistPosEffMatrixGen(0),
fHistNegEffMatrixGen(0),
fHistPosEffMatrixScanRec(0),
fHistNegEffMatrixScanRec(0),
fHistPosEffMatrixScanGen(0),
fHistNegEffMatrixScanGen(0),
fHistEmptyEvent(0),
fHistCentrality(0),
fHistCentralityImpPar(0),
fHistImpParam(0),
fHistVertex(0),
fHistdEdxTPC(0),
fHistArmPod(0),
fEventInfo_PhiTPCdcarA(0),
fEventInfo_PhiTPCdcarC(0),
fEventInfo_CacheTrackCounters(0),
fEventInfo_CacheTrackdEdxRatio(0),
fEventInfo_CacheTrackNcl(0),
fEventInfo_CacheTrackChi2(0),
fEventInfo_CacheTrackMatchEff(0),
fEventInfo_CentralityEstimates(0),
fEventInfo_LumiGraph(0),
fEventInfo_HisTPCVertexA(0),
fEventInfo_HisTPCVertexC(0),
fEventInfo_HisTPCVertexACut(0),
fEventInfo_HisTPCVertexCCut(0),
fEventInfo_HisTPCVertex(0),
fEventInfo_CacheTrackTPCCountersZ(0),
fPileUpTightnessCut1(0),
fPileUpTightnessCut2(0),
fPileUpTightnessCut3(0),
fPileUpTightnessCut4(0)
{
  //
  //         standard constructur which should be used
  //
  /* fast compilation test
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  gSystem->AddIncludePath("-I/lustre/nyx/alice/users/marsland/alicehub/sw/ubuntu1404_x86-64/AliPhysics/master-1/include");
  gSystem->AddIncludePath("-I/lustre/nyx/alice/users/marsland/alicehub/AliPhysics/PWGPP");
  .L /home/marsland/Desktop/RUN_ON_GRID/Ebye/code/AliAnalysisTaskEbyeIterPID.cxx++
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/AliAnalysisTaskEbyeIterPID.cxx++
  .L /lustre/nyx/alice/users/marsland/train/trunk/marsland_EbyeRatios/AliAnalysisTaskEbyeIterPID.cxx++
  */
  std::cout << " Info::marsland:===================================================================================="<< std::endl;
  std::cout << " Info::marsland:===================================================================================="<< std::endl;
  std::cout << " Info::marsland:===================================================================================="<< std::endl;
  std::cout << " Info::marsland:***************** CONSTRUCTOR CALLED: AliAnalysisTaskEbyeIterPID  *****************"<< std::endl;
  std::cout << " Info::marsland:===================================================================================="<< std::endl;
  std::cout << " Info::marsland:===================================================================================="<< std::endl;
  std::cout << " Info::marsland:===================================================================================="<< std::endl;
  // ==========================================
  //
  // ==========================================
  // Initialize arrays
  for (Int_t ires=0; ires<2; ires++){
    for (Int_t imom=0; imom<4; imom++){
      for (Int_t icent=0; icent<10; icent++){
        for (Int_t ieta=0; ieta<8; ieta++){
          fNetPiFirstMoments[ires][imom][icent][ieta]=0.;
          fNetKaFirstMoments[ires][imom][icent][ieta]=0.;
          fNetPrFirstMoments[ires][imom][icent][ieta]=0.;
          fNetLaFirstMoments[ires][imom][icent][ieta]=0.;
          fNetChFirstMoments[ires][imom][icent][ieta]=0.;
        }
      }
    }
  }
  // Initialize arrays
  for (Int_t imom=0; imom<4; imom++){
    for (Int_t icent=0; icent<10; icent++){
      for (Int_t ieta=0; ieta<8; ieta++){
        fNetPiFirstMomentsGen[imom][icent][ieta]=0.;
        fNetKaFirstMomentsGen[imom][icent][ieta]=0.;
        fNetPrFirstMomentsGen[imom][icent][ieta]=0.;
        fNetPiFirstMomentsRec[imom][icent][ieta]=0.;
        fNetKaFirstMomentsRec[imom][icent][ieta]=0.;
        fNetPrFirstMomentsRec[imom][icent][ieta]=0.;
      }
    }
  }
  for (Int_t imom=0; imom<4; imom++){
    for (Int_t icent=0; icent<10; icent++){
      for (Int_t ieta=0; ieta<8; ieta++){
        fCrossPiFirstMomentsGen[imom][icent][ieta]=0.;
        fCrossKaFirstMomentsGen[imom][icent][ieta]=0.;
        fCrossPrFirstMomentsGen[imom][icent][ieta]=0.;
        fCrossPiFirstMomentsRec[imom][icent][ieta]=0.;
        fCrossKaFirstMomentsRec[imom][icent][ieta]=0.;
        fCrossPrFirstMomentsRec[imom][icent][ieta]=0.;
      }
    }
  }
  // Initialize arrays
  for (Int_t isign=0; isign<2; isign++){
    for (Int_t imom=0; imom<4; imom++){
      for (Int_t icent=0; icent<10; icent++){
        for (Int_t ieta=0; ieta<8; ieta++){
          fPiFirstMomentsGen[isign][imom][icent][ieta]=0.;
          fKaFirstMomentsGen[isign][imom][icent][ieta]=0.;
          fPrFirstMomentsGen[isign][imom][icent][ieta]=0.;
          fPiFirstMomentsRec[isign][imom][icent][ieta]=0.;
          fKaFirstMomentsRec[isign][imom][icent][ieta]=0.;
          fPrFirstMomentsRec[isign][imom][icent][ieta]=0.;
        }
      }
    }
  }
  // Initialize miaaing cluter map
  // const Int_t nParticles=4;
  // const Int_t nCentBins=9;
  // for (Int_t ipart=0; ipart<nParticles; ipart++){
  //   for (Int_t icent=0; icent<nCentBins; icent++){
  //     fH2MissCl[ipart][icent] = 0;
  //   }
  // }
  //
  // ==========================================
  //
  // Define outputs
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class());
  DefineOutput(5, TTree::Class());
  DefineOutput(6, TTree::Class());
  DefineOutput(7, TTree::Class());
  DefineOutput(8, TTree::Class());
  DefineOutput(9, TTree::Class());
  DefineOutput(10, TTree::Class());
  DefineOutput(11, TTree::Class());
  DefineOutput(12, TTree::Class());
  DefineOutput(13, TTree::Class());
  DefineOutput(14, TTree::Class());
  // ==========================================

}
//________________________________________________________________________
AliAnalysisTaskEbyeIterPID::~AliAnalysisTaskEbyeIterPID()
{

  //
  // Destructor
  //
  std::cout << " Info::marsland: ===== In the Destructor ===== " << std::endl;
  if (fHistPosEffMatrixRec) delete fHistPosEffMatrixRec;
  if (fHistNegEffMatrixRec) delete fHistNegEffMatrixRec;
  if (fHistPosEffMatrixGen) delete fHistPosEffMatrixGen;
  if (fHistNegEffMatrixGen) delete fHistNegEffMatrixGen;
  if (fHistPosEffMatrixScanRec) delete fHistPosEffMatrixScanRec;
  if (fHistNegEffMatrixScanRec) delete fHistNegEffMatrixScanRec;
  if (fHistPosEffMatrixScanGen) delete fHistPosEffMatrixScanGen;
  if (fHistNegEffMatrixScanGen) delete fHistNegEffMatrixScanGen;
  if (fHistdEdxTPC)         delete fHistdEdxTPC;
  if (fHistEmptyEvent)      delete fHistEmptyEvent;
  if (fHistCentrality)      delete fHistCentrality;
  if (fHistCentralityImpPar)delete fHistCentralityImpPar;
  if (fHistImpParam)        delete fHistImpParam;
  if (fHistVertex)          delete fHistVertex;
  if (fHistArmPod)          delete fHistArmPod;
  if (fHistCent)            delete fHistCent;
  if (fHistPhi)             delete fHistPhi;
  if (fHistGenMult)         delete fHistGenMult;
  if (fHistRapDistFullAccPr)   delete fHistRapDistFullAccPr;
  if (fHistRapDistFullAccAPr)  delete fHistRapDistFullAccAPr;
  if (fHistInvK0s)          delete fHistInvK0s;
  if (fHistInvLambda)       delete fHistInvLambda;
  if (fHistInvAntiLambda)   delete fHistInvAntiLambda;
  if (fHistInvPhoton)       delete fHistInvPhoton;
  if (fHndEdx)              delete fHndEdx;
  for (Int_t i=0;i<fNSettings;i++) { if (fHnExpected[i])   delete fHnExpected[i]; }
  //
  // Marians histograms
  if (fHistPhiTPCcounterA)  delete fHistPhiTPCcounterA;
  if (fHistPhiTPCcounterC)  delete fHistPhiTPCcounterC;
  if (fHistPhiTPCcounterAITS)  delete fHistPhiTPCcounterAITS;
  if (fHistPhiTPCcounterCITS)  delete fHistPhiTPCcounterCITS;
  if (fHistPhiITScounterA)  delete fHistPhiITScounterA;
  if (fHistPhiITScounterC)  delete fHistPhiITScounterC;
  if (fEventInfo_HisTPCVertexA)  delete fEventInfo_HisTPCVertexA;
  if (fEventInfo_HisTPCVertexC)  delete fEventInfo_HisTPCVertexC;
  if (fEventInfo_HisTPCVertex)  delete fEventInfo_HisTPCVertex;
  if (fEventInfo_HisTPCVertexACut)  delete fEventInfo_HisTPCVertexACut;
  if (fEventInfo_HisTPCVertexCCut)  delete fEventInfo_HisTPCVertexCCut;
  if (fEventInfo_PhiTPCdcarA)  delete fEventInfo_PhiTPCdcarA;
  if (fEventInfo_PhiTPCdcarC)  delete fEventInfo_PhiTPCdcarC;
  if (fEventInfo_CacheTrackCounters)  delete fEventInfo_CacheTrackCounters;
  if (fEventInfo_CacheTrackdEdxRatio)  delete fEventInfo_CacheTrackdEdxRatio;
  if (fEventInfo_CacheTrackNcl)  delete fEventInfo_CacheTrackNcl;
  if (fEventInfo_CacheTrackChi2)  delete fEventInfo_CacheTrackChi2;
  if (fEventInfo_CacheTrackMatchEff)  delete fEventInfo_CacheTrackMatchEff;
  if (fEventInfo_CentralityEstimates)  delete fEventInfo_CentralityEstimates;
  if (fEventInfo_CacheTrackTPCCountersZ)  delete fEventInfo_CacheTrackTPCCountersZ;
  if (fPIDCombined) delete fPIDCombined;
  if (fESDtrackCuts)          delete fESDtrackCuts;
  if (fESDtrackCuts_Bit96)    delete fESDtrackCuts_Bit96;
  if (fESDtrackCuts_Bit96_spd)    delete fESDtrackCuts_Bit96_spd;
  if (fESDtrackCuts_Bit96_sdd)    delete fESDtrackCuts_Bit96_sdd;
  if (fESDtrackCuts_Bit128)   delete fESDtrackCuts_Bit128;
  if (fESDtrackCuts_Bit768)   delete fESDtrackCuts_Bit768;
  if (fESDtrackCutsLoose)     delete fESDtrackCutsLoose;
  if (fESDtrackCutsV0)        delete fESDtrackCutsV0;
  if (fTreeSRedirector)       delete fTreeSRedirector;
  if (fESDtrackCutsCleanSamp) delete fESDtrackCutsCleanSamp;
  if (fPileUpTightnessCut4) delete fPileUpTightnessCut4;
  if (fPileUpTightnessCut3) delete fPileUpTightnessCut3;
  if (fPileUpTightnessCut2) delete fPileUpTightnessCut2;
  if (fPileUpTightnessCut1) delete fPileUpTightnessCut1;




}
//
// ---------------------------------------------------------------------------------
//                                     Functions
// ---------------------------------------------------------------------------------
//
void AliAnalysisTaskEbyeIterPID::Initialize()
{
  //
  // updating parameters in case of changes (standard cuts and the eta window)
  //
  std::cout << " Info::marsland: ===== In the Initialize ===== " << std::endl;
  if (fRunFastSimulation)    { std::cout << " Info::marsland: !!! We are running fast simulation return !!! " << std::endl; return; }
  if (fRunFastHighMomentCal) { std::cout << " Info::marsland: !!! We are running fast high moment calculation return !!! " << std::endl; return; }
  // fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
  //
  // ------------------------------------------------
  //
  // tight DCA cut used by Emil
  fESDtrackCuts_Bit96_spd = new AliESDtrackCuts("Bit96_spd",""); fESDtrackCuts_Bit96_spd->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
  fESDtrackCuts_Bit96_sdd = new AliESDtrackCuts("Bit96_sdd",""); fESDtrackCuts_Bit96_sdd->SetClusterRequirementITS(AliESDtrackCuts::kSDD, AliESDtrackCuts::kFirst);
  fESDtrackCuts_Bit96     = new AliESDtrackCuts("fESDtrackCuts_Bit96","");
  fESDtrackCuts_Bit96->SetEtaRange(-100.,100.);
  fESDtrackCuts_Bit96->SetPtRange(0.,100000.);
  fESDtrackCuts_Bit96->SetMinNClustersTPC(70); // ???? should be --> fESDtrackCuts_Bit96->SetMinNCrossedRowsTPC(70);
  fESDtrackCuts_Bit96->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fESDtrackCuts_Bit96->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts_Bit96->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts_Bit96->SetRequireITSRefit(kTRUE);
  fESDtrackCuts_Bit96->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  fESDtrackCuts_Bit96->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
  fESDtrackCuts_Bit96->SetMaxChi2TPCConstrainedGlobal(36);
  fESDtrackCuts_Bit96->SetMaxDCAToVertexZ(2);
  fESDtrackCuts_Bit96->SetDCAToVertex2D(kFALSE);
  fESDtrackCuts_Bit96->SetRequireSigmaToVertex(kFALSE);
  fESDtrackCuts_Bit96->SetMaxChi2PerClusterITS(36);
  if ( (fYear==2015&&fPassIndex==2) || (fYear==2018&&fPassIndex==3) ){
    fESDtrackCuts_Bit96->SetMaxChi2PerClusterTPC(2.5);
  } else {
    fESDtrackCuts_Bit96->SetMaxChi2PerClusterTPC(4);
  }
  //
  // TPC only tracks
  fESDtrackCuts_Bit128 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fESDtrackCuts_Bit128->SetName("Bit128");
  fESDtrackCuts_Bit128->SetEtaRange(-100.,100.);
  fESDtrackCuts_Bit128->SetPtRange(0.,100000.);
  fESDtrackCuts_Bit128->SetMinNClustersTPC(70);
  fESDtrackCuts_Bit128->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts_Bit128->SetMaxDCAToVertexZ(3.2);
  fESDtrackCuts_Bit128->SetMaxDCAToVertexXY(2.4);
  fESDtrackCuts_Bit128->SetDCAToVertex2D(kTRUE);
  if ( (fYear==2015&&fPassIndex==2) || (fYear==2018&&fPassIndex==3) ){
    fESDtrackCuts_Bit128->SetMaxChi2PerClusterTPC(2.5);
  } else {
    fESDtrackCuts_Bit128->SetMaxChi2PerClusterTPC(4);
  }
  //
  // Hybrid tracks
  fESDtrackCuts_Bit768 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
  fESDtrackCuts_Bit768->SetName("Bit768");
  fESDtrackCuts_Bit768->SetEtaRange(-100.,100.);
  fESDtrackCuts_Bit768->SetPtRange(0.,100000.);
  fESDtrackCuts_Bit768->SetMaxDCAToVertexXY(2.4);
  fESDtrackCuts_Bit768->SetMaxDCAToVertexZ(3.2);
  fESDtrackCuts_Bit768->SetDCAToVertex2D(kTRUE);
  fESDtrackCuts_Bit768->SetMaxChi2TPCConstrainedGlobal(36);
  fESDtrackCuts_Bit768->SetMaxFractionSharedTPCClusters(0.4);
  fESDtrackCuts_Bit768->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
  fESDtrackCuts_Bit768->SetRequireITSRefit(kTRUE);
  if ( (fYear==2015&&fPassIndex==2) || (fYear==2018&&fPassIndex==3) ){
    fESDtrackCuts_Bit768->SetMaxChi2PerClusterTPC(2.5);
  } else {
    fESDtrackCuts_Bit768->SetMaxChi2PerClusterTPC(4);
  }
  //
  // for the systematic check fill all tracks and tag them with cutbit but for MC do not
  //
  fESDtrackCuts = new AliESDtrackCuts("esdTrackCuts","");
  fESDtrackCuts->SetEtaRange(-100.,100.);
  fESDtrackCuts->SetPtRange(0.,100000.);
  fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
  fESDtrackCuts->SetMaxChi2PerClusterITS(36);
  fESDtrackCuts->SetMaxFractionSharedTPCClusters(0.4);    // ?? FROM MARIAN
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetRequireITSRefit(kTRUE);
  fESDtrackCuts->SetMinNCrossedRowsTPC(80);
  fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0208+0.04/pt^1.01");
  fESDtrackCuts->SetMaxDCAToVertexXY(2.4);   // hybrid cuts  TODO
  fESDtrackCuts->SetMaxDCAToVertexZ(3.2);    // hybrid cuts  TODO
  fESDtrackCuts->SetRequireSigmaToVertex(kFALSE);
  fESDtrackCuts->SetDCAToVertex2D(kTRUE);  // fESDtrackCuts->SetDCAToVertex2D(kFALSE);    TODO
  if ( (fYear==2015&&fPassIndex==2) || (fYear==2018&&fPassIndex==3) ){
    fESDtrackCuts->SetMaxChi2PerClusterTPC(2.5);
  } else {
    fESDtrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if (fIncludeITS) {
    // require ITS pixels  -->  Reason for the empty events and structure in phi
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  }
  //
  // very loose cuts --> cuts will be tightened using the bipmap
  fESDtrackCutsLoose = new AliESDtrackCuts("esdTrackCutsLoose","");
  fESDtrackCutsLoose->SetEtaRange(-100.,100.);
  fESDtrackCutsLoose->SetPtRange(0.1,100000.);
  fESDtrackCutsLoose->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fESDtrackCutsLoose->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCutsLoose->SetMaxFractionSharedTPCClusters(0.4);
  fESDtrackCutsLoose->SetMinNClustersTPC(50);
  fESDtrackCutsLoose->SetMinNCrossedRowsTPC(50);
  fESDtrackCutsLoose->SetMaxDCAToVertexXY(10);   // hybrid cuts  TODO
  fESDtrackCutsLoose->SetMaxDCAToVertexZ(10);    // hybrid cuts  TODO
  //
  // track cuts to be used for v0s
  fESDtrackCutsCleanSamp = new AliESDtrackCuts("AliESDtrackCutsV0","");
  fESDtrackCutsCleanSamp -> SetEtaRange(-1.5,1.5);
  fESDtrackCutsCleanSamp -> SetPtRange(0.,1e10);
  fESDtrackCutsCleanSamp -> SetMinNCrossedRowsTPC(80);
  fESDtrackCutsCleanSamp -> SetRequireTPCRefit(kTRUE);
  fESDtrackCutsCleanSamp -> SetAcceptKinkDaughters(kFALSE);
  fESDtrackCutsCleanSamp -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fESDtrackCutsCleanSamp -> SetMaxChi2PerClusterITS(36);
  fESDtrackCutsCleanSamp -> SetMaxFractionSharedTPCClusters(0.4);
  // ------------------------------------------------
  // V0 selection
  // ------------------------------------------------
  fESDtrackCutsV0 = new AliESDv0Cuts("AliESDCutsV0","");
  fESDtrackCutsV0 ->SetMaxDcaV0Daughters(1.0);
  // ------------------------------------------------
  //
  // Special selection for clean samples
  fV0OpenCuts   = new AliESDv0KineCuts();
  fV0StrongCuts = new AliESDv0KineCuts();
  SetSpecialV0Cuts(fV0OpenCuts);
  SetSpecialV0Cuts(fV0StrongCuts);

  fPileUpTightnessCut4 = new AliEventCuts(kFALSE);
  fPileUpTightnessCut3 = new AliEventCuts(kFALSE);
  fPileUpTightnessCut2 = new AliEventCuts(kFALSE);
  fPileUpTightnessCut1 = new AliEventCuts(kFALSE);

  fPileUpTightnessCut4->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,4);
  fPileUpTightnessCut3->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,3);
  fPileUpTightnessCut2->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,2);
  fPileUpTightnessCut1->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,1);
  //
  //
  std::cout << " Info::marsland: ===================================================== " << std::endl;
  std::cout << " Info::marsland: =============== Summary of Track Cuts =============== " << std::endl;
  std::cout << " Info::marsland: ===================================================== " << std::endl;
}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::UserCreateOutputObjects()
{
  //
  // Create output histograms, trees of the analysis (called once)
  //
  if (!fDEdxCheck) Initialize();
  std::cout << " Info::marsland: ===== In the UserCreateOutputObjects ===== " << std::endl;
  // ------------  setup PIDCombined  ---------------
  fPIDCombined = new AliPIDCombined("pidCombined","");
  //
  // **********************   Input handler to get the PID object *********************
  if (!(fRunFastSimulation || fRunFastHighMomentCal)) {
    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
    if (!inputHandler)
    AliFatal("Input handler needed");
    else {
      fPIDResponse = inputHandler->GetPIDResponse();       // PID response object
      if (!fPIDResponse) std::cout << " Info::marsland: ======= PIDResponse object was not created ====== " << std::endl;
    }
  }
  //
  // ************************************************************************
  //   OpenFile output --> one can open several files
  // ************************************************************************
  //
  OpenFile(1);
  fTreeSRedirector = new TTreeSRedirector();
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);
  //
  //
  if (fDefaultEventCuts) fEventCuts.AddQAplotsToList(fListHist); /// fList is your output TList
  //
  // ************************************************************************
  //   histogram of splines
  // ************************************************************************
  //
  if(!fEffMatrix && !fMCtrue)
  {
    // 0 --> assumed particle: 0. electron, 1. pion, 2. kaon, 3. proton, 5. deuteron
    // 1 --> sign
    // 2 --> Centrality
    // 3 --> eta
    // 4 --> momentum
    // 5 --> Expected Sigma of a given track
    // 6 --> Expected Mean of a given track
    const Int_t nExpectedbins = 7;
    //                                        0    1,    2,                 3,           4,        5,       6
    Int_t   binsExpected[nExpectedbins]  = {  5,   3,  fNCentbinsData,   fNEtaBins,   fNMomBins,   240,   4000 };  // ????
    Double_t xminExpected[nExpectedbins] = {  0., -2.,   0.,             fEtaDown,    fMomDown,     1.,   20.  };
    Double_t xmaxExpected[nExpectedbins] = {  5.,  2.,  80.,             fEtaUp,      fMomUp,      61.,   1020.};
    TString axisNameExpected[nExpectedbins]   = {"particleType","sign","Centrality"    ,"eta" ,"momentum" ,"ExSigma","ExMean"};
    TString axisTitleExpected[nExpectedbins]  = {"particleType","sign","Centrality [%]","#eta","#it{p} (GeV/#it{c})", "#sigma","#mu"};
    for (Int_t i=0;i<fNSettings;i++){
      if (fDefaultTrackCuts && i==1) break;
      fHnExpected[i] = new THnSparseF(Form("hExpected_%d",i),Form("hExpected_%d",i),nExpectedbins,binsExpected,xminExpected,xmaxExpected);
      std::cout << "Number of centrality bins: " << fxCentBins.size() << std::endl;
      fHnExpected[i]->GetAxis(2)->Set(fNCentbinsData-1,fxCentBins.data());
      for (Int_t iaxis=0; iaxis<nExpectedbins;iaxis++){
        fHnExpected[i]->GetAxis(iaxis)->SetName(axisNameExpected[iaxis]);
        fHnExpected[i]->GetAxis(iaxis)->SetTitle(axisTitleExpected[iaxis]);
      }
      fListHist->Add(fHnExpected[i]);
    }
  }
  //
  // ************************************************************************
  //   dEdx histograms
  // ************************************************************************
  //
  // 0 --> sign
  // 1 --> Centrality
  // 2 --> eta
  // 3 --> momentum
  // 4 --> TPC dEdx
  //                                  0,         1,           2,           3,            4
  const Int_t nhistbins = 5;
  Int_t dEdxnBins = Int_t((fDEdxUp-fDEdxDown)/fDEdxBinWidth);
  TString axisNamedEdx[nhistbins]   = {"Sign" ,"Centrality"     ,"eta"  ,"momentum"  ,"TPC dEdx"};
  TString axisTitledEdx[nhistbins]  = {"Sign" ,"Centrality [%]" ,"#eta" ,"#it{p} (GeV/#it{c})" ,"TPC d#it{E}/d#it{x} Signal (a.u.)"};
  if (fUseThnSparse)
  {
    //
    // inclusive spectra
    Int_t   binsdEdx[nhistbins]  = { 2,  fNCentbinsData, fNEtaBins,   fNMomBins,    dEdxnBins};
    Double_t xmindEdx[nhistbins] = {-2,  0.,             fEtaDown,    fMomDown,     fDEdxDown};
    Double_t xmaxdEdx[nhistbins] = { 2,  80.,            fEtaUp,      fMomUp,       fDEdxUp};
    fHndEdx= new THnSparseF("hdEdx","Inclusive dEdx Spectrum"  ,nhistbins,binsdEdx,xmindEdx,xmaxdEdx);
    fHndEdx->GetAxis(1)->Set(fNCentbinsData-1,fxCentBins.data());
    // Set the branch names
    for (Int_t iaxis=0; iaxis<nhistbins;iaxis++){
      fHndEdx   ->GetAxis(iaxis)->SetName(axisNamedEdx[iaxis]);
      fHndEdx   ->GetAxis(iaxis)->SetTitle(axisTitledEdx[iaxis]);
    }
    fListHist->Add(fHndEdx);
  }
  //
  // ************************************************************************
  //   Efficiency matrix histograms
  // ************************************************************************
  //
  if(fEffMatrix && !fRunOnGrid)
  {
    const Int_t ndim=5;
    Int_t nbins0[ndim]  ={3,8, 40       ,16        ,50  };
    Double_t xmin0[ndim]={0,0, fMomDown ,fEtaDown  ,0.  };
    Double_t xmax0[ndim]={3,80,fMomUp   ,fEtaUp    ,6.25};
    fHistPosEffMatrixRec  =new THnF("fHistPosEffMatrixRec","fHistPosEffMatrixRec",ndim, nbins0,xmin0,xmax0);
    fHistNegEffMatrixRec  =new THnF("fHistNegEffMatrixRec","fHistNegEffMatrixRec",ndim, nbins0,xmin0,xmax0);
    fHistPosEffMatrixGen  =new THnF("fHistPosEffMatrixGen","fHistPosEffMatrixGen",ndim, nbins0,xmin0,xmax0);
    fHistNegEffMatrixGen  =new THnF("fHistNegEffMatrixGen","fHistNegEffMatrixGen",ndim, nbins0,xmin0,xmax0);
    TString axisNameEff[ndim]  = {"particle"      ,"Centrality"     ,"momentum"      ,"eta"  ,"phi"};
    TString axisTitleEff[ndim] = {"particle type" ,"Centrality (%)" ,"#it{p}_{T} (GeV/#it{c})" ,"#eta" ,"#phi"};
    for (Int_t iEff=0; iEff<ndim;iEff++){
      fHistPosEffMatrixRec->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistPosEffMatrixRec->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
      fHistNegEffMatrixRec->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistNegEffMatrixRec->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
      fHistPosEffMatrixGen->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistPosEffMatrixGen->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
      fHistNegEffMatrixGen->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistNegEffMatrixGen->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
    }
    fListHist->Add(fHistPosEffMatrixRec);
    fListHist->Add(fHistNegEffMatrixRec);
    fListHist->Add(fHistPosEffMatrixGen);
    fListHist->Add(fHistNegEffMatrixGen);
    //
    //
    const Int_t ndimScan=6;
    Int_t nbinsScan[ndimScan]   = {2, fNSettings,           3,8, 40       ,16   };
    Double_t xminScan[ndimScan] = {0, 0,                    0,0, fMomDown ,fEtaDown };
    Double_t xmaxScan[ndimScan] = {2, Double_t(fNSettings), 3,80,fMomUp   ,fEtaUp };
    fHistPosEffMatrixScanRec  =new THnF("fHistPosEffMatrixScanRec","fHistPosEffMatrixScanRec",ndimScan, nbinsScan,xminScan,xmaxScan);
    fHistNegEffMatrixScanRec  =new THnF("fHistNegEffMatrixScanRec","fHistNegEffMatrixScanRec",ndimScan, nbinsScan,xminScan,xmaxScan);
    fHistPosEffMatrixScanGen  =new THnF("fHistPosEffMatrixScanGen","fHistPosEffMatrixScanGen",ndimScan, nbinsScan,xminScan,xmaxScan);
    fHistNegEffMatrixScanGen  =new THnF("fHistNegEffMatrixScanGen","fHistNegEffMatrixScanGen",ndimScan, nbinsScan,xminScan,xmaxScan);
    TString axisNameEffScan[ndimScan]  = {"detector", "systematic", "particle"      ,"Centrality"     ,"momentum"      ,"eta"  };
    TString axisTitleEffScan[ndimScan] = {"detector", "setting",    "particle type" ,"Centrality (%)" ,"#it{p}_{T} (GeV/#it{c})" ,"#eta"};
    for (Int_t iEff=0; iEff<ndimScan;iEff++){
      fHistPosEffMatrixScanRec->GetAxis(iEff)->SetName(axisNameEffScan[iEff]);  fHistPosEffMatrixScanRec->GetAxis(iEff)->SetTitle(axisTitleEffScan[iEff]);
      fHistNegEffMatrixScanRec->GetAxis(iEff)->SetName(axisNameEffScan[iEff]);  fHistNegEffMatrixScanRec->GetAxis(iEff)->SetTitle(axisTitleEffScan[iEff]);
      fHistPosEffMatrixScanGen->GetAxis(iEff)->SetName(axisNameEffScan[iEff]);  fHistPosEffMatrixScanGen->GetAxis(iEff)->SetTitle(axisTitleEffScan[iEff]);
      fHistNegEffMatrixScanGen->GetAxis(iEff)->SetName(axisNameEffScan[iEff]);  fHistNegEffMatrixScanGen->GetAxis(iEff)->SetTitle(axisTitleEffScan[iEff]);
    }
    fListHist->Add(fHistPosEffMatrixScanRec);
    fListHist->Add(fHistNegEffMatrixScanRec);
    fListHist->Add(fHistPosEffMatrixScanGen);
    fListHist->Add(fHistNegEffMatrixScanGen);
  }
  //
  // ************************************************************************
  //   Event histograms
  // ************************************************************************
  //
  Int_t nSafeBins = (fRunFastSimulation) ? 1200 : 5;
  fHistEmptyEvent        = new TH1F("hEmptyEvent",           "control histogram to count empty events"    , 10,  0., 10.);
  fHistCentrality        = new TH1F("hCentrality",           "control histogram for centrality"           , 10,  0., 100.);
  fHistCentralityImpPar  = new TH1F("hCentralityImpPar",     "control histogram for centrality imppar"    , 10,  0., 100.);
  fHistImpParam          = new TH1F("hImpParam",             "control histogram for impact parameter"     , 200, 0., 20.);
  fHistVertex            = new TH1F("hVertex",               "control histogram for vertex Z position"    , 200, -20., 20.);
  fHistGenMult           = new TH1F("hGenPrMult",            "generated protons"                          , fGenprotonBins,0., 200.);
  fHistRapDistFullAccPr  = new TH2F("hRapDistFullAccPr",     "rapidity dist of protons in full acceptance"    , nSafeBins, -15, 15., 10, 0., 100.);
  fHistRapDistFullAccAPr = new TH2F("hRapDistFullAccApr",    "rapidity dist of antiprotons in full acceptance", nSafeBins, -15, 15., 10, 0., 100.);
  fListHist->Add(fHistEmptyEvent);
  fListHist->Add(fHistCentrality);
  fListHist->Add(fHistCentralityImpPar);
  fListHist->Add(fHistImpParam);
  fListHist->Add(fHistVertex);
  fListHist->Add(fHistGenMult);
  fListHist->Add(fHistRapDistFullAccPr);
  fListHist->Add(fHistRapDistFullAccAPr);
  fEventInfo_CentralityEstimates  = new TVectorF(3);
  for (Int_t i=0;i<3;i++) (*fEventInfo_CentralityEstimates)[i]=-10.;
  //
  // ************************************************************************
  //   Marians counters
  // ************************************************************************
  //
  if (fEventInfo)
  {
    // vectors
    fEventInfo_PhiTPCdcarA         = new TVectorF(36);
    fEventInfo_PhiTPCdcarC         = new TVectorF(36);
    fEventInfo_CacheTrackCounters  = new TVectorF(20);
    fEventInfo_CacheTrackdEdxRatio = new TVectorF(30);
    fEventInfo_CacheTrackNcl       = new TVectorF(20);
    fEventInfo_CacheTrackChi2      = new TVectorF(20);
    fEventInfo_CacheTrackMatchEff  = new TVectorF(20);
    fEventInfo_CacheTrackTPCCountersZ = new TVectorF(8);
    for (Int_t i=0;i<8;i++) (*fEventInfo_CacheTrackTPCCountersZ)[i]=0.;
    for (Int_t i=0;i<36;i++){
      (*fEventInfo_PhiTPCdcarA)[i]=0.;
      (*fEventInfo_PhiTPCdcarC)[i]=0.;
    }
    for (Int_t i=0;i<20;i++){
      (*fEventInfo_CacheTrackCounters)[i]=0.;
      (*fEventInfo_CacheTrackdEdxRatio)[i]=0.;
      (*fEventInfo_CacheTrackNcl)[i]=0.;
      (*fEventInfo_CacheTrackChi2)[i]=0.;
      (*fEventInfo_CacheTrackMatchEff)[i]=0.;
    }
    // Hists
    fHistPhiTPCcounterA    = new TH1F("hPhiTPCcounterC",       "control histogram to count tracks on the A side in phi ", 36, 0.,18.);
    fHistPhiTPCcounterC    = new TH1F("hPhiTPCcounterA",       "control histogram to count tracks on the C side in phi ", 36, 0.,18.);
    fHistPhiTPCcounterAITS = new TH1F("hPhiTPCcounterAITS",    "control histogram to count tracks on the A side in phi ", 36, 0.,18.);
    fHistPhiTPCcounterCITS = new TH1F("hPhiTPCcounterCITS",    "control histogram to count tracks on the C side in phi ", 36, 0.,18.);
    fHistPhiITScounterA    = new TH1F("hPhiITScounterA",       "control histogram to count tracks on the A side in phi ", 36, 0.,18.);
    fHistPhiITScounterC    = new TH1F("hPhiITScounterC",       "control histogram to count tracks on the C side in phi ", 36, 0.,18.);
    fEventInfo_HisTPCVertexA = new TH1F("hisTPCZA", "hisTPCZA", 1000, -250, 250);
    fEventInfo_HisTPCVertexC = new TH1F("hisTPCZC", "hisTPCZC", 1000, -250, 250);
    fEventInfo_HisTPCVertex = new TH1F("hisTPCZ", "hisTPCZ", 1000, -250, 250);
    fEventInfo_HisTPCVertexACut = new TH1F("hisTPCZACut", "hisTPCZACut", 1000, -250, 250);
    fEventInfo_HisTPCVertexCCut = new TH1F("hisTPCZCCut", "hisTPCZCCut", 1000, -250, 250);
    fEventInfo_HisTPCVertex->SetLineColor(1);
    fEventInfo_HisTPCVertexA->SetLineColor(2);
    fEventInfo_HisTPCVertexC->SetLineColor(4);
    fEventInfo_HisTPCVertexACut->SetLineColor(3);
    fEventInfo_HisTPCVertexCCut->SetLineColor(6);
  }
  //
  // ************************************************************************
  //   Clean sample helper histograms
  // ************************************************************************
  //
  if (fUseCouts && fFillArmPodTree)
  {
    fHistArmPod            = new TH2F("hArmPod",           "Armenteros-Podolanski plot"                     , 100,-1.,1., 110,0.,0.22);
    fHistInvK0s            = new TH1F("fHistInvK0s",       "control histogram for K0s invariant mass"       , 1000, 0.3,  0.70);
    fHistInvLambda         = new TH1F("fHistInvLambda",    "control histogram for lambda invariant mass"    , 1000, 1.07, 1.16);
    fHistInvAntiLambda     = new TH1F("fHistInvAntiLambda","control histogram for antilambda invariant mass", 1000, 1.07, 1.16);
    fHistInvPhoton         = new TH1F("fHistInvPhoton",    "control histogram for photon invariant mass"    , 1000, 0.,   0.05);
    fListHist->Add(fHistInvK0s);
    fListHist->Add(fHistInvLambda);
    fListHist->Add(fHistInvAntiLambda);
    fListHist->Add(fHistInvPhoton);
    fListHist->Add(fHistArmPod);
  }
  //
  // ************************************************************************
  //   Trees
  // ************************************************************************
  //
  fArmPodTree    = ((*fTreeSRedirector)<<"fArmPodTree").GetTree();
  fTreeMCFull    = ((*fTreeSRedirector)<<"mcFull").GetTree();
  fTreeMCgen     = ((*fTreeSRedirector)<<"mcGen").GetTree();
  fTreeDnchDeta  = ((*fTreeSRedirector)<<"dnchdeta").GetTree();
  fTreeMC        = ((*fTreeSRedirector)<<"fTreeMC").GetTree();
  fTreedEdxCheck = ((*fTreeSRedirector)<<"dEdxCheck").GetTree();
  fTreeCuts      = ((*fTreeSRedirector)<<"tracks").GetTree();
  fTreeMCFullAcc = ((*fTreeSRedirector)<<"fullacc").GetTree();
  fTreeResonance = ((*fTreeSRedirector)<<"resonance").GetTree();
  fTreeMCgenMoms = ((*fTreeSRedirector)<<"mcGenMoms").GetTree();
  fTreeEvents    = ((*fTreeSRedirector)<<"eventInfo").GetTree();
  fTreeDScaled   = ((*fTreeSRedirector)<<"dscaled").GetTree();
  fTreeMCEffCorr = ((*fTreeSRedirector)<<"mcMoms").GetTree();
  //
  // ************************************************************************
  //   Send output objects to container
  // ************************************************************************
  //
  PostData(1, fListHist);
  PostData(2, fArmPodTree);
  PostData(3, fTreeMCFull);
  PostData(4, fTreeMCgen);
  PostData(5, fTreeMC);
  PostData(6, fTreedEdxCheck);
  PostData(7, fTreeCuts);
  PostData(8, fTreeDnchDeta);
  PostData(9, fTreeMCFullAcc);
  PostData(10, fTreeResonance);
  PostData(11, fTreeMCgenMoms);
  PostData(12, fTreeEvents);
  PostData(13, fTreeDScaled);
  PostData(14, fTreeMCEffCorr);

  std::cout << " Info::marsland: ===== Out of UserCreateOutputObjects ===== " << std::endl;

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::UserExec(Option_t *)
{
  //
  // main event loop
  //
  if (fRunOnGrid) fUseCouts=kFALSE; // for security
  if (fUseCouts) std::cout << " Info::marsland: ===== In the UserExec ===== " << std::endl;
  //
  // Check Monte Carlo information and other access first:
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) fMCtrue = kFALSE;
  fEventCountInFile++;
  //
  //  get the filename
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { Printf(" Error::marsland: Could not receive input chain"); return; }
  TString tmpChunkname = fChunkName;
  TObjString fileName(chain->GetCurrentFile()->GetName());
  fChunkName = (TString)fileName.GetString();
  if (tmpChunkname != fChunkName) std::cout <<  " Info::marsland: ===== Current chunk name is ===== " << fChunkName << std::endl;
  //
  // ======================================================================
  // ========================== See if MC or Real =========================
  // ======================================================================
  //
  if (eventHandler) fMCEvent = eventHandler->MCEvent();
  AliGenEventHeader* genHeader = 0x0;
  TString genheaderName;
  if (fMCEvent){
    genHeader = fMCEvent->GenEventHeader();
    genheaderName = genHeader->GetName();
    if(!genHeader){ Printf(" Error::marsland: Event generator header not available!!!\n"); return; }
  }
  //
  // If ESDs exist get some event variables
  //
  fCentrality = -5;
  fCentImpBin =-10.;
  AliCentrality    *esdCentrality = 0x0;
  AliMultSelection *MultSelection = 0x0;
  AliMultSelectionTask *MultSelectionTask = 0x0;
  ULong64_t gid=0;
  fESD = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (fESD)
  {
    //
    // Init magnetic filed for golden chi2 cut
    fESD->InitMagneticField();
    //
    // event selection
    fPileUpBit=0;
    if (fDefaultEventCuts && fESD){

      if ( (fPassIndex==3 || fPassIndex==2) && fYear>2013){
        //
        // pileup bit: 0bxxxx, where the leftmost bit is the tightest and the rightmost is the loosest cut
        // OOB pileup cut (for Pb-Pb) based on ITS and TPC clusters: 0-> no cut; 1-> default cut (remove all OOB pileup); 2-> looser cut; 3-> even more looser cut; 4-> very loose cut
        if (fPileUpTightnessCut4->AcceptEvent(fESD)) { fPileUpBit |= 1 << 3; if (fUseCouts) std::cout << "pileupbit: " << std::bitset<4>(fPileUpBit) << std::endl;}
        if (fPileUpTightnessCut3->AcceptEvent(fESD)) { fPileUpBit |= 1 << 2; if (fUseCouts) std::cout << "pileupbit: " << std::bitset<4>(fPileUpBit) << std::endl;}
        if (fPileUpTightnessCut2->AcceptEvent(fESD)) { fPileUpBit |= 1 << 1; if (fUseCouts) std::cout << "pileupbit: " << std::bitset<4>(fPileUpBit) << std::endl;}
        if (fPileUpTightnessCut1->AcceptEvent(fESD)) { fPileUpBit |= 1 << 0; if (fUseCouts) std::cout << "pileupbit: " << std::bitset<4>(fPileUpBit) << std::endl;}
        fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,0); // do not apply any pile cut
        if (!fEventCuts.AcceptEvent(fESD)) {cout<< "pileup event " << endl; return;}
      }
    }
    //
    // calculate TPC mult
    fTPCMult = 0;
    for (Int_t itrack=0;itrack<fESD->GetNumberOfTracks();++itrack){
      AliESDtrack *track = fESD->GetTrack(itrack);
      if (track->IsOn(AliESDtrack::kTPCin)) fTPCMult++;
    }
    //
    //
    esdCentrality = fESD->GetCentrality();
    MultSelection = (AliMultSelection*) fESD-> FindListObject("MultSelection");
    fRunNo = fESD->GetRunNumber();
    //
    static Int_t timeStampCache = -1;
    if (!fMCtrue) {
      fTimeStamp = fESD->GetTimeStampCTPBCCorr();
      const char *ocdb;
      if(!fRunOnGrid) ocdb = Form("local:///cvmfs/alice.cern.ch/calibration/data/%d/OCDB",fYear);
      else ocdb = "raw://";
      //
      // retrieve interaction rate
      if (timeStampCache!=Int_t(fTimeStamp)) timeStampCache=Int_t(fTimeStamp);
      if (!gGrid && timeStampCache>0) {
        AliInfo("Trying to connect to AliEn ...");
        TGrid::Connect("alien://");
      } else {
        fEventInfo_LumiGraph = (TGraph*)AliLumiTools::GetLumiFromCTP(fRunNo,ocdb);
        fIntRate   = fEventInfo_LumiGraph->Eval(fTimeStamp); delete fEventInfo_LumiGraph;
      }
    }
    fEventMult = fESD->GetNumberOfTracks();
    fBField    = fESD->GetMagneticField();
    fEvent     = fESD->GetEventNumberInFile();
    fBeamType  = fESD->GetBeamType();
    //
    // Global id for the event --> which is made with Hashing
    //
    ULong64_t orbitID      = (ULong64_t)fESD->GetOrbitNumber();
    ULong64_t bunchCrossID = (ULong64_t)fESD->GetBunchCrossNumber();
    ULong64_t periodID     = (ULong64_t)fESD->GetPeriodNumber();
    gid = ((periodID << 36) | (orbitID << 12) | bunchCrossID);
    fEventGID = gid;    // uniqe event id for real data
    // fTimeStamp             = fESD->GetTimeStamp();
    // fEventGID              = TMath::Abs(Int_t(TString::Hash(&gid,sizeof(Int_t))));    // uniqe event id for real data
  }
  //
  // Get rid of "E-AliESDpid::GetTPCsignalTunedOnData: Tune On Data requested, but MC event not set. Call SetCurrentMCEvent before!" errors
  if (!fPIDResponse && !(fRunFastSimulation || fRunFastHighMomentCal)) fPIDResponse = ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetESDpid();
  if (fMCEvent && fPIDResponse) fPIDResponse->SetCurrentMCEvent(fMCEvent);
  //
  if(fMCtrue)
  {
    //
    // Add different generator particles to PDG Data Base to avoid problems when reading MC generator particles
    AliPDG::AddParticlesToPdgDataBase();
    //
    // ========================== MC =========================
    //
    fMCStack = fMCEvent->Stack();
    if (!fMCStack) { Printf(" Error::marsland: No MC stack available !!!\n"); return;}
    //
    if (MultSelection) {
      fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
      if (fUseCouts)  std::cout << " Info::marsland: Centralitity is taken from MultSelection " << fCentrality << std::endl;
    } else if (esdCentrality) {
      fCentrality = esdCentrality->GetCentralityPercentile("V0M");
      if (fUseCouts)  std::cout << " Info::marsland: Centralitity is taken from esdCentrality " << fCentrality << std::endl;
    }
    //
    // AliCollisionGeometry* colGeometry = dynamic_cast<AliCollisionGeometry*>(genHeader);
    // std::cout << " aaaa  " <<  ((AliGenEposEventHeader*) genHeader)->ImpactParameter() << std::endl;
    // std::cout << " bbbb  " <<  colGeometry->ImpactParameter() << std::endl;
    //
    //
    // impact parameters to use: 0.0, 3.72, 5.23, 7.31, 8.88, 10.20, 11.38, 12.47, 13.50, 14.5
    // corresponding Centrality:  0     5    10    20    30     40     50    60      70    80
    Double_t impParArr[10] = {0.0, 3.72, 5.23, 7.31, 8.88, 10.20, 11.38, 12.47, 13.50, 14.5};
    AliGenHijingEventHeader *lHIJINGHeader = 0x0;  // event header for HIJING
    AliGenHepMCEventHeader *lHepMCHeader = 0x0;    // event header for EPOS
    //
    // If EPOS based MC
    if ( genheaderName.Contains("EPOSLHC") && fMCEvent ){
      if (genHeader->InheritsFrom(AliGenHepMCEventHeader::Class()))
      lHepMCHeader = (AliGenHepMCEventHeader*)genHeader;
      if (lHepMCHeader ){
        fNHardScatters = lHepMCHeader->Ncoll_hard(); // Number of hard scatterings
        fNProjectileParticipants = lHepMCHeader->Npart_proj(); // Number of projectile participants
        fNTargetParticipants     = lHepMCHeader->Npart_targ(); // Number of target participants
        fNNColl   = lHepMCHeader->Ncoll(); // Number of NN (nucleon-nucleon) collisions
        fNNwColl  = lHepMCHeader->N_Nwounded_collisions(); // Number of N-Nwounded collisions
        fNwNColl  = lHepMCHeader->Nwounded_N_collisions(); // Number of Nwounded-N collisons
        fNwNwColl = lHepMCHeader->Nwounded_Nwounded_collisions();// Number of Nwounded-Nwounded collisions
        fMCImpactParameter = lHepMCHeader->impact_parameter();
        if (fUseCouts)  std::cout << " Info::marsland: EPOS: Centralitity is taken from ImpactParameter = " << fMCImpactParameter << "  "  << ((AliGenEposEventHeader*) genHeader)->GetName() << std::endl;
        fHistImpParam->Fill(fMCImpactParameter);
        if (fMCImpactParameter>=impParArr[0] && fMCImpactParameter<impParArr[1]) fCentImpBin=2.5;
        if (fMCImpactParameter>=impParArr[1] && fMCImpactParameter<impParArr[2]) fCentImpBin=7.5;
        if (fMCImpactParameter>=impParArr[2] && fMCImpactParameter<impParArr[3]) fCentImpBin=15.;
        if (fMCImpactParameter>=impParArr[3] && fMCImpactParameter<impParArr[4]) fCentImpBin=25.;
        if (fMCImpactParameter>=impParArr[4] && fMCImpactParameter<impParArr[5]) fCentImpBin=35.;
        if (fMCImpactParameter>=impParArr[5] && fMCImpactParameter<impParArr[6]) fCentImpBin=45.;
        if (fMCImpactParameter>=impParArr[6] && fMCImpactParameter<impParArr[7]) fCentImpBin=55.;
        if (fMCImpactParameter>=impParArr[7] && fMCImpactParameter<impParArr[8]) fCentImpBin=65.;
        if (fMCImpactParameter>=impParArr[8] && fMCImpactParameter<impParArr[9]) fCentImpBin=75.;
        if (fMCImpactParameter<impParArr[0]  || fMCImpactParameter>impParArr[9]) fCentImpBin=-10.;
        fHistCentralityImpPar->Fill(fCentImpBin);
      }
    } else if (!TMath::IsNaN(((AliGenHijingEventHeader*) genHeader)->ImpactParameter()) && fMCEvent){
      //
      // If HIJING based MC
      lHIJINGHeader = (AliGenHijingEventHeader*) genHeader;
      fNHardScatters = lHIJINGHeader->HardScatters();
      fNProjectileParticipants = lHIJINGHeader->ProjectileParticipants();
      fNTargetParticipants     = lHIJINGHeader->TargetParticipants();
      fNNColl   = lHIJINGHeader->NN();
      fNNwColl  = lHIJINGHeader->NNw();
      fNwNColl  = lHIJINGHeader->NwN();
      fNwNwColl = lHIJINGHeader->NwNw();
      fMCImpactParameter = lHIJINGHeader->ImpactParameter();
      fHistImpParam->Fill(fMCImpactParameter);
      if (fUseCouts)  std::cout << " Info::marsland: HIJING: Centralitity is taken from ImpactParameter = " << fMCImpactParameter << "  "  << ((AliGenHijingEventHeader*) genHeader)->GetName() << std::endl;
      if (fMCImpactParameter>=impParArr[0] && fMCImpactParameter<impParArr[1]) fCentImpBin=2.5;
      if (fMCImpactParameter>=impParArr[1] && fMCImpactParameter<impParArr[2]) fCentImpBin=7.5;
      if (fMCImpactParameter>=impParArr[2] && fMCImpactParameter<impParArr[3]) fCentImpBin=15.;
      if (fMCImpactParameter>=impParArr[3] && fMCImpactParameter<impParArr[4]) fCentImpBin=25.;
      if (fMCImpactParameter>=impParArr[4] && fMCImpactParameter<impParArr[5]) fCentImpBin=35.;
      if (fMCImpactParameter>=impParArr[5] && fMCImpactParameter<impParArr[6]) fCentImpBin=45.;
      if (fMCImpactParameter>=impParArr[6] && fMCImpactParameter<impParArr[7]) fCentImpBin=55.;
      if (fMCImpactParameter>=impParArr[7] && fMCImpactParameter<impParArr[8]) fCentImpBin=65.;
      if (fMCImpactParameter>=impParArr[8] && fMCImpactParameter<impParArr[9]) fCentImpBin=75.;
      if (fMCImpactParameter<impParArr[0]  || fMCImpactParameter>impParArr[9]) fCentImpBin=-10.;
      fHistCentralityImpPar->Fill(fCentImpBin);
    }

    if (fCentrality<0) fCentrality=fCentImpBin;
    //
    // Use file name in Hashing to create unique event ID
    fEventGIDMC  = TMath::Abs(Int_t(TString::Hash(&fEventCountInFile,sizeof(Int_t))));    // uniqe event id for real data
    // fEventGIDMC += TMath::Abs(Int_t(fCentrality)+fEventCountInFile+(1000*TMath::Abs(fMCImpactParameter)));  // ????
    fEventGID    = fEventGIDMC;
    if (fUseCouts) {
      std::cout << " Info::marsland: ========================================================================================== " << std::endl;
      std::cout << " Info::marsland: " << fEventCountInFile << " ----- eventIDMC = " << fEventGIDMC << "   " << fChunkName << std::endl;
      std::cout << " Info::marsland: Centrality = " << fCentrality << " ----- Impact Param = " << fMCImpactParameter << " fCentralityImp = " << fCentImpBin << std::endl;
      std::cout << " Info::marsland: ========================================================================================== " << std::endl;
    }
    //
  }
  //
  if (!(fRunFastSimulation || fRunFastHighMomentCal))
  {
    //
    // ========================== Real =========================
    //
    if (!fESD)          { Printf(" Error::marsland: fESD not available"); return; }
    if (!fESDtrackCuts) { Printf(" Error::marsland: fESDtrackCuts not available"); return; }
    //
    // ------------------------------------------------
    // ------- monitor vertex position =---------------
    // ------------------------------------------------
    //
    Bool_t isVertexOk = kTRUE;
    fVertex = fESD->GetPrimaryVertexTracks();
    const AliESDVertex *vertexSPD = fESD->GetPrimaryVertexTracks();
    const AliESDVertex *vertexTPC = fESD->GetPrimaryVertexTracks();
    if( fVertex->GetNContributors()<1) isVertexOk = kFALSE;
    if( fVertex->GetNContributors()>1) {
      vertexSPD = fESD->GetPrimaryVertexSPD();    // SPD vertex
      vertexTPC = fESD->GetPrimaryVertexTPC();    // TPC vertex
      fTPCvZ = vertexTPC->GetZ();
      fSPDvZ = vertexSPD->GetZ();
      fVz    = fVertex->GetZ();
      TString vertexType = fVertex->GetTitle();    // ??? Put condition Abs(vertex-vertexTPC) as a bool_t into ttree
      if ( vertexType.Contains("vertexer: Z") && (fVertex->GetDispersion() > 0.04 || fVertex->GetZRes() > 0.25) ) isVertexOk = kFALSE; // TODO
    }
    fMultiplicity    = fVertex->GetNContributors();    // fMultiplicity = fESD -> GetNumberOfTracks();
    fNContributors   = fVertex->GetNContributors();
    fMultiplicityMC  = fMultiplicity;
    //
    // ------------------------------------------------
    // ------- event vertex cut along Z ---------------
    // ------------------------------------------------
    //
    // if (fMCtrue && TMath::Abs(fVz) > 15) return;   // For MC put fixed cut
    // if (fDefaultTrackCuts && (TMath::Abs(fVz)>7 || TMath::Abs(fVz)<0.15) ) return;
    // else if (TMath::Abs(fVz)>15) return;
    if (TMath::Abs(fVz)>15) return;
    //
    if (fVertex && isVertexOk) fHistVertex->Fill(fVz);
    else return;
    //
    // ------------------------------------------------
    // ---------- Centrality definition ---------------
    // ------------------------------------------------
    //
    Int_t nEst = sizeof(fEventInfo_centEstStr)/sizeof(char*);
    fEventInfo_CentralityEstimates->Zero();  // matchEff counter
    if (fBeamType.CompareTo("A-A") == 0) { // PbPb
      if (MultSelection) {
        if(fSystCentEstimatetor == -1) fCentrality = MultSelection->GetMultiplicityPercentile("TRK");
        if(fSystCentEstimatetor ==  0) fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
        if(fSystCentEstimatetor ==  1) fCentrality = MultSelection->GetMultiplicityPercentile("CL1");
        for (Int_t i=0;i<nEst;i++) (*fEventInfo_CentralityEstimates)[i]=MultSelection->GetMultiplicityPercentile(fEventInfo_centEstStr[i]);
      } else if (esdCentrality) {
        if(fSystCentEstimatetor == -1) fCentrality = esdCentrality->GetCentralityPercentile("TRK");
        if(fSystCentEstimatetor ==  0) fCentrality = esdCentrality->GetCentralityPercentile("V0M");
        if(fSystCentEstimatetor ==  1) fCentrality = esdCentrality->GetCentralityPercentile("CL1");
        for (Int_t i=0;i<nEst;i++) (*fEventInfo_CentralityEstimates)[i]=esdCentrality->GetCentralityPercentile(fEventInfo_centEstStr[i]);
      } else {
        std::cout << " Info::marsland: Error: There is no cent info " << std::endl;
      }
    }
    //
    if (fUseCouts) {
      std::cout << " Info::marsland: =============================================================================================== " << std::endl;
      std::cout << " Info::marsland: Event counter = " << fEventCountInFile << " - cent =  " << fCentrality << " = gid = " << gid << " = fEventGID = " << fEventGID << std::endl;
      std::cout << " Info::marsland: =============================================================================================== " << std::endl;
    }
  }
  fHistCentrality->Fill(fCentrality);  // count events after physics and vertex selection
  //
  // if (!fESDtool) {
  //   fESDtool = new AliESDtools();
  //   fESDtool->SetStreamer(fTreeSRedirector);
  // }
  // fESDtool->Init(NULL,fESD);
  // fESDtool->CalculateEventVariables();
  // fESDtool->SetMCEvent(fMCEvent);
  // fESDtool->DumpEventVariables();
  //
  // in case small stat is enough
  if (fPercentageOfEvents>0 && (fEventCountInFile%fPercentageOfEvents)==0) return;
  //
  // ======================================================================
  //   Filling part
  // ======================================================================
  //
  if (fDEdxCheck)                          { FillTPCdEdxCheck(); return;}
  if (fRunFastSimulation && fFillDnchDeta) { FillDnchDeta(); return;}
  if (fRunFastHighMomentCal)               { FastGenHigherMoments(); return;}
  if (fRunFastSimulation)                  { FastGen_NetParticles(); return;}
  //
  // Real Data Analysis
  //
  if (!fMCtrue && fFillTracks && fESD){
    FillTPCdEdxReal();
    if (fFillArmPodTree) FillCleanSamples();
    if (fEventInfo) {CalculateEventInfo(); CreateEventInfoTree();}
    if (fUseCouts)  std::cout << " Info::marsland: (Real Data Analysis) End of Filling part = " << fEventCountInFile << std::endl;
    return;
  }
  //
  // full MC analysis
  //
  if (fMCtrue && fEffMatrix && fESD){
    // FillEffMatrix();
    FillMCFull_NetParticles();
    if (fFillArmPodTree) FillCleanSamples();
    if (fEventInfo) {CalculateEventInfo(); CreateEventInfoTree();}
    if (fUseCouts)  std::cout << " Info::marsland: (full MC analysis) End of Filling part = " << fEventCountInFile << std::endl;
    return;
  }
  //
  if (fMCtrue && fFillHigherMomentsMCclosure){
    MCclosureHigherMoments();
    if (fUseCouts)  std::cout << " Info::marsland: (MCclosure Higher Moments) End of Filling part = " << fEventCountInFile << std::endl;
    return;
  }


}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FillTPCdEdxReal()
{
  //
  // Fill dEdx information for the TPC and also clean kaon and protons
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillTPCdEdxReal ===== " << std::endl;
  // --------------------------------------------------------------
  // Get the event
  AliVEvent *event=InputEvent();
  if (CountEmptyEvents(2)<1) return;
  //
  // --------------------------------------------------------------
  //  Main track loop
  // --------------------------------------------------------------
  //
  Int_t tpcClusterMultiplicity   = fESD->GetNumberOfTPCClusters();
  const AliMultiplicity *multObj = fESD->GetMultiplicity();
  Int_t itsNumberOfTracklets   = multObj->GetNumberOfTracklets();
  for (Int_t itrack=0;itrack<event->GetNumberOfTracks();++itrack) {   // Track loop

    fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
    fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
    //
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack *track = fESD->GetTrack(itrack);
    if (!track) continue;
    //
    // --------------------------------------------------------------
    //      Get relevant track info and set cut bits
    // --------------------------------------------------------------
    //
    Bool_t ifDefaultCuts = fESDtrackCuts->AcceptTrack(track);
    Bool_t fBit96_spd    = fESDtrackCuts_Bit96_spd->AcceptTrack(track);
    Bool_t fBit96_sdd    = fESDtrackCuts_Bit96_sdd->AcceptTrack(track);
    Bool_t fBit96_base   = fESDtrackCuts_Bit96->AcceptTrack(track);
    Bool_t fBit128       = fESDtrackCuts_Bit128->AcceptTrack(track);
    Bool_t fBit768       = fESDtrackCuts_Bit768->AcceptTrack(track);
    Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(track);
    if (!track->GetInnerParam()) continue;               // Ask if track is in the TPC
    if (!fESDtrackCutsLoose->AcceptTrack(track))  continue;    // Loose cuts
    if (!(track->GetTPCsignalN()>0)) continue;
    //
    // Get the track variables
    Float_t closestPar[3];
    GetExpecteds(track,closestPar);
    SetCutBitsAndSomeTrackVariables(track,0);
    Int_t tpcNcls = track->GetTPCncls();
    //
    // --------------------------------------------------------------
    //  Some print out
    // --------------------------------------------------------------
    //
    // Tree for the all cut variables
    if (fUseCouts && fEventCountInFile==5 && !fRunOnGrid) {
      std::cout << " Info::marsland: CutBinMap --> " <<fTrackTPCCrossedRows << " " << fTrackChi2TPC << " " <<  fTrackNewITScut  << std::endl;
      PrintNumInBinary(fTrackCutBits);
    }
    //
    // --------------------------------------------------------------
    //   Fill the trees
    // --------------------------------------------------------------
    //
    Int_t nTPCClusters = fESD->GetNumberOfTPCClusters();
    Int_t nITSClusters = 0;
    AliVMultiplicity *multiObj = fESD->GetMultiplicity();
    for(Int_t i=2;i<6;i++) nITSClusters += multiObj->GetNumberOfITSClusters(i);
    //
    // different dca cuts
    // TMath::Abs(fTrackDCAxy)< 0.3
    Bool_t dca11h     = TMath::Abs(fTrackDCAxy)<0.0105+0.0350/TMath::Power(fPt,1.1);    // 10h tuned loose cut
    Bool_t dca10h     = TMath::Abs(fTrackDCAxy)<0.0182+0.0350/TMath::Power(fPt,1.01);    // 10h tuned loose cut
    Bool_t dcaBaseCut = TMath::Abs(fTrackDCAxy)<0.0208+0.04/TMath::Power(fPt,1.01);  // 10h tuned loose cut
    //
    if (fFillTracks && !fFillOnlyHists)
    {
      if(!fTreeSRedirector) return;
      if (fFillDistributions) FillTrackVariables(track);
      (*fTreeSRedirector)<<"tracks"<<
      "defCut="    << ifDefaultCuts <<  // default cuts tuned by hand
      "bit96spd="  << fBit96_spd <<     // cut account for the regions where no ITS-TPC match
      "bit96sdd="  << fBit96_sdd <<     // cut account for the regions where no ITS-TPC match
      "bit96="     << fBit96_base <<    // tight cuts of 2011 tuned data
      "bit128="    << fBit128 <<        // TPC only tracks cuts
      "bit768="    << fBit768 <<        // Hybrid track cuts
      "pixCut="    << ifDCAcutIfNoITSPixel <<    // cut: apply a DCAcut If No ITS Pixel
      "run="       << fRunNo <<                  // run Number
      "bField="    << fBField <<                 // magnetic filed
      "pileupbit=" << fPileUpBit <<              // flag for pileup selection
      "primMult="  << fNContributors <<          //  #prim tracks
      "tpcClMult=" << tpcClusterMultiplicity <<  //  TPC cluster multiplicity
      //
      "dcabase="  << dcaBaseCut <<  //  TPC multiplicity
      "dca10h="   << dca10h <<  //  TPC multiplicity
      "dca11h="   << dca11h <<  //  TPC multiplicity
      //
      "tpcmult="   << fTPCMult <<                //  TPC track multiplicity
      "itsmult="   << itsNumberOfTracklets <<    // ITS multiplicity
      "itsclmult=" << nITSClusters <<    // ITS multiplicity
      "tpcclmult=" << nTPCClusters <<    // ITS multiplicity
      //
      "gid="       << fEventGID             <<  //  global event ID
      "eventtime=" << fTimeStamp            <<  // event timeStamp
      "intrate="   << fIntRate              <<  // interaction rate
      "cutBit="    << fTrackCutBits         <<  //  Systematic Cuts
      "dEdx="      << fTPCSignal            <<  //  dEdx of the track
      "sign="      << fSign                 <<  //  charge
      "ptot="      << fPtot                 <<  //  TPC momentum
      "p="         << fPVertex              <<  //  momentum at vertex
      "pT="        << fPt                   <<  // transverse momentum
      "eta="       << fEta                  <<  //  eta
      "cent="      << fCentrality           <<  //  centrality
      //
      "phi="       << fPhi                  <<  //  phi
      "dcaxy="     << fTrackDCAxy           <<  // dca cut on xy plane
      "dcaz="      << fTrackDCAz            <<  // dca cut along z
      "ncltpc="    << fNcl                  <<  // number of clusters
      "cRows="     << fTrackTPCCrossedRows  <<  // crossed Rows in TPC
      "chi2tpc="   << fTrackChi2TPC         <<  // TPC chi2
      "missCl="    << fMissingCl            <<  // fraction of missing clusters
      //
      "eltpcpid="  << fNSigmasElTPC         <<  // nsigma TPC for electrons
      "pitpcpid="  << fNSigmasPiTPC         <<
      "katpcpid="  << fNSigmasKaTPC         <<
      "prtpcpid="  << fNSigmasPrTPC         <<
      "prtofpid="  << fNSigmasPrTOF         <<
      "\n";
    }
    //
    // --------------------------------------------------------------
    //  Fill the THnSparseF for the Expected values form PID response
    // --------------------------------------------------------------
    //
    // define acceptance of interest
    Bool_t etaAcc  = (fEta >=fEtaDown       && fEta<=fEtaUp);
    Bool_t momAcc  = (fPVertex>=fMomDown    && fPVertex<=fMomUp);
    Bool_t dEdxAcc = (fTPCSignal>=fDEdxDown && fTPCSignal<=fDEdxUp);
    Bool_t fAcceptance = (etaAcc && momAcc && dEdxAcc);
    Bool_t nSigmasElTPCCut = (TMath::Abs(fNSigmasElTPC)<2);
    Bool_t nSigmasPiTPCCut = (TMath::Abs(fNSigmasPiTPC)<2);
    Bool_t nSigmasKaTPCCut = (TMath::Abs(fNSigmasKaTPC)<2);
    Bool_t nSigmasPrTPCCut = (TMath::Abs(fNSigmasPrTPC)<2);
    Bool_t nSigmasDeTPCCut = (TMath::Abs(fNSigmasDeTPC)<2);
    Bool_t nSigmaTPCall = (nSigmasElTPCCut || nSigmasPiTPCCut || nSigmasKaTPCCut || nSigmasPrTPCCut || nSigmasDeTPCCut);
    Bool_t ndEdxTPCall  = (fDEdxEl>20 || fDEdxPi>20 || fDEdxKa>20 || fDEdxPr>20 || fDEdxDe>20);

    if(!fEffMatrix && fAcceptance && !fMCtrue){
      for (Int_t i=0; i<fNSettings; i++) {
        Double_t exMean[5]  = {fDEdxEl,  fDEdxPi,  fDEdxKa,  fDEdxPr,  fDEdxDe};
        Double_t exSigma[5] = {fSigmaEl, fSigmaPi, fSigmaKa, fSigmaPr, fSigmaDe};
        if ( ndEdxTPCall && nSigmaTPCall && ifDefaultCuts )
        {

          if (fDefaultTrackCuts && i==0){
            for (Int_t iPart = 0; iPart< 5; iPart++){
              Double_t weightExpected[7] = {Double_t(iPart),Double_t(fSign),fCentrality,fEta,fPtot, exSigma[iPart], exMean[iPart]};
              fHnExpected[i]->Fill(weightExpected);
            }
            break;
          }

          if (GetSystematicClassIndex(fTrackCutBits,i))
          {
            for (Int_t iPart = 0; iPart< 5; iPart++){
              Double_t weightExpected[7] = {Double_t(iPart),Double_t(fSign),fCentrality,fEta,fPtot, exSigma[iPart], exMean[iPart]};
              fHnExpected[i]->Fill(weightExpected);
            }
          }
        }

      }
    }
    //
    // --------------------------------------------------------------
    //  Fill thnsparse for inclusive data
    // --------------------------------------------------------------
    //
    if (fUseThnSparse){
      Double_t trackdEdx[5] = {Double_t(fSign),fCentrality, fEta,fPtot, fTPCSignal};
      if(fUseThnSparse) fHndEdx->Fill(trackdEdx);
    }

  } // end of track loop

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FillTrackVariables(AliESDtrack *track)
{

  /*

  // ITSTPC standard cuts
  // TPC
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  // ITS
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);

  // Additional cuts
  esdTrackCuts->SetDCAToVertex2D(kTRUE);
  esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
  esdTrackCuts->SetMaxFractionSharedTPCClusters(0.4);

  */


  //
  // Apply downclaing condition
  Double_t qP      = track->Charge()/track->P();
  if ( !(fRandom.Rndm()*(qP*qP) < 0.05) ) return;
  AliTPCdEdxInfo tpcdEdxInfo;
  //
  // Calculate further variables
  Int_t tpcCrossedRows=0, tpcSignalN=0;
  Double_t eta   = track->Eta();
  Double_t tgl   = track->Pz()/track->Pt();
  Double_t phi   = track->Phi()-TMath::Pi(); // ????
  Int_t sign     = track->GetSign();
  Double_t phi85 = track->GetParameterAtRadius(85,5,7);
  ULong64_t flag = track->GetStatus();
  Bool_t isOnITS = track->IsOn(AliESDtrack::kITSrefit);
  Bool_t isOnTRD = track->IsOn(AliESDtrack::kTRDrefit);
  Bool_t isOnTPC = track->IsOn(AliESDtrack::kTPCrefit);
  Int_t nclTPC   = track->GetTPCncls(); if (nclTPC<1) nclTPC=-1;
  Int_t nclITS   = track->GetITSNcls(); if (nclITS<1) nclITS=-1;
  Int_t nclTRD   = track->GetTRDncls(); if (nclTRD<1) nclTRD=-1;
  Double_t chi2TPC = TMath::Sqrt(TMath::Abs(track->GetTPCchi2()/nclTPC));
  Double_t chi2ITS = TMath::Sqrt(TMath::Abs(track->GetITSchi2()));
  Double_t chi2TRD = TMath::Sqrt(TMath::Abs(track->GetTRDchi2()));
  Double_t itsdEdx = track->GetITSsignal();
  Double_t trddEdx = track->GetTRDsignal();
  Double_t tpcdEdx = track->GetTPCsignal();
  Double_t ptot0   = track->GetP();
  Double_t pt      = track->Pt();
  UChar_t itsclmap = track->GetITSClusterMap();
  Float_t pv[2],cov[3];
  track->GetImpactParameters(pv,cov); // p[0]=fD; p[1]=fZ; cov[0]=fCdd; cov[1]=fCdz; cov[2]=fCzz;
  //
  // --------------------------------------------------------------
  //      TPC related observables
  // --------------------------------------------------------------
  //
  Double_t sharedTPCClusters=0.;
  Double_t phiTPC=-100.;
  Double_t ptotTPC=0.;
  Double_t lengthInActiveZone=0.;
  Float_t pTPC[2],covTPC[3];          // p[0]=fdTPC; p[1]=fzTPC; cov[0]=fCddTPC; cov[1]=fCdzTPC; cov[2]=fCzzTPC;
  if (track->GetInnerParam()){
    track->GetTPCdEdxInfo(tpcdEdxInfo);
    phiTPC  = track->GetInnerParam()->GetParameterAtRadius(85,5,7);
    ptotTPC = track->GetInnerParam()->GetP();
    lengthInActiveZone = track->GetLengthInActiveZone(1,3,230, track->GetBz(),0,0);
    track->GetImpactParametersTPC(pTPC,covTPC);
    tpcCrossedRows = track->GetTPCCrossedRows();
    tpcSignalN = track->GetTPCsignalN();
    Float_t closestPar[3];
    GetExpecteds(track,closestPar);
    SetCutBitsAndSomeTrackVariables(track,0);
    if(nclTPC) sharedTPCClusters = static_cast<Double_t>(track->GetTPCnclsS())/static_cast<Double_t>(nclTPC);

  }
  //
  // pt dependent nclusters cut
  TFormula *f1NClustersTPCLinearPtDepSmall = new TFormula("f1NClustersTPCLinearPtDepSmall","60.+30./20.*x");
  TFormula *f1NClustersTPCLinearPtDep      = new TFormula("f1NClustersTPCLinearPtDep"     ,"70.+30./20.*x");
  TFormula *f1NClustersTPCLinearPtDepLArge = new TFormula("f1NClustersTPCLinearPtDepLArge","80.+30./20.*x");
  Int_t nClsPtDepSmall = (Int_t)(f1NClustersTPCLinearPtDepSmall->Eval(track->Pt()));
  Int_t nClsPtDep      = (Int_t)(f1NClustersTPCLinearPtDep     ->Eval(track->Pt()));
  Int_t nClsPtDepLarge = (Int_t)(f1NClustersTPCLinearPtDepLArge->Eval(track->Pt()));
  delete f1NClustersTPCLinearPtDepSmall;
  delete f1NClustersTPCLinearPtDep;
  delete f1NClustersTPCLinearPtDepLArge;
  //
  // --------------------------------------------------------------
  //   Fill downscaled tree
  // --------------------------------------------------------------
  //
  if(!fTreeSRedirector) return;
  (*fTreeSRedirector)<<"dscaled"<<
  "gid="                  << fEventGID             <<  //  global event ID
  "eventtime="            << fTimeStamp            <<
  "intrate="              << fIntRate              <<  // interaction rate
  "cutBit="               << fTrackCutBits         <<  //  Systematic Cuts
  "dEdx="                 << fTPCSignal            <<  //  dEdx of the track
  "sign="                 << fSign                 <<  //  charge
  "ptot="                 << fPtot                 <<  //  TPC momentum
  "p="                    << fPVertex              <<  //  TPC momentum
  "pT="                   << fPt                   <<
  "eta="                  << fEta                  <<  //  eta
  "cent="                 << fCentrality           <<  //  centrality
  "nclTPC="               << nclTPC                <<  //  #ITS clusters
  "dcaxy="                << fTrackDCAxy           <<  // fD pv[0]
  "dcaz="                 << fTrackDCAz            <<  // fZ pv[1]
  //
  "isOnITS="              << isOnITS         <<
  "isOnTRD="              << isOnTRD         <<
  //
  "itsclmap="             << itsclmap              <<  //  vertex Z
  "flag="                 << flag                  <<
  "sharedTPCClusters="    << sharedTPCClusters     <<
  "tpcSignalN="           << tpcSignalN            <<  //  number of cl used in dEdx
  "cRows="                << tpcCrossedRows        <<  //  crossed rows
  "lengthInActiveZone="   << lengthInActiveZone    <<  //  fTrackLengthInActiveZone in TPC
  "phi="                  << phi                   <<  //  ph
  "phiTPC="               << phiTPC                <<
  "phi85="                << phi85                 <<
  "qP="                   << qP                    <<  //  charge/momentu,
  "tgl="                  << tgl                   <<  //  tangent
  "pt="                   << pt                    <<  //  pT
  //
  "vz="                   << fVz                   <<  //  vertex Z
  "centArray.="           << fEventInfo_CentralityEstimates  <<  // track counter
  "fdTPC="                << pTPC[0]               <<
  "fzTPC="                << pTPC[1]               <<
  "fCddTPC="              << covTPC[0]                <<  //  DCAxy
  "fCdzTPC="              << covTPC[1]                <<  //  DCAz
  "fCzzTPC="              << covTPC[2]                <<  //  DCAz
  //
  "fd="                   << pv[0]               <<
  "fz="                   << pv[1]               <<
  "fCdd="                 << cov[0]                <<  //  DCAxy
  "fCdz="                 << cov[1]                <<  //  DCAz
  "fCzz="                 << cov[2]                <<  //  DCAz
  //
  "itsdEdx="              << itsdEdx               <<
  "trddEdx="              << trddEdx               <<
  "nclits="               << nclITS                <<  //  #ITS clusters
  "ncltrd="               << nclTRD                <<  //  #TRD clusters
  "chi2tpc="              << chi2TPC               <<  //  TPC chi2
  "chi2its="              << chi2ITS               <<  //  ITS chi2
  "chi2trd="              << chi2TRD               <<  //  TRD chi2
  //
  "itspixel01="           << fIsITSpixel01         <<
  "primRes="              << fPrimRestriction      <<
  "\n";

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FillMCFull()
{

  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  // Assign subsample index
  Int_t sampleNo = 0;
  Int_t nSubSample = 20;
  sampleNo = Int_t(fEventGID)%nSubSample;
  Int_t runNumber = fESD->GetRunNumber();
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillMCFull ===== " << std::endl;
  //
  // ======================================================================
  // For Marian
  AliTPCdEdxInfo tpcdEdxInfo;
  // ======================================================================
  //
  // Fill dEdx tree for MC closure
  for(Int_t irectrack = 0; irectrack < fESD->GetNumberOfTracks(); irectrack++)
  {
    //
    // initialize the dummy particle id
    fElMC =-100.; fPiMC =-100.; fKaMC =-100.; fPrMC =-100.; fDeMC =-100.; fMuMC =-100.;
    // Esd track
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack *trackReal = fESD->GetTrack(irectrack);
    if (trackReal==NULL) continue;
    // Get generated track info
    Int_t lab = TMath::Abs(trackReal->GetLabel());           // avoid from negatif labels, they include some garbage
    TParticle *trackMC  = fMCStack->Particle(lab);
    Int_t pdg = trackMC->GetPdgCode();
    //
    Float_t dcaprim[2], covprim[3];
    Float_t dcaweak[2], covweak[3];
    Float_t dcamaterial[2], covmaterial[3];
    if(fMCStack->IsSecondaryFromMaterial(lab))  trackReal->GetImpactParameters(dcamaterial, covmaterial);
    if(fMCStack->IsSecondaryFromWeakDecay(lab)) trackReal->GetImpactParameters(dcaweak,     covweak);
    if(fMCStack->IsPhysicalPrimary(lab))        trackReal->GetImpactParameters(dcaprim,     covprim);
    //
    // track cuts
    Bool_t ifDefaultCuts = fESDtrackCuts->AcceptTrack(trackReal);
    Bool_t fBit96_spd = fESDtrackCuts_Bit96_spd->AcceptTrack(trackReal);
    Bool_t fBit96_sdd = fESDtrackCuts_Bit96_sdd->AcceptTrack(trackReal);
    Bool_t fBit96_base = fESDtrackCuts_Bit96->AcceptTrack(trackReal);
    Bool_t fBit128       = fESDtrackCuts_Bit128->AcceptTrack(trackReal);
    Bool_t fBit768       = fESDtrackCuts_Bit768->AcceptTrack(trackReal);
    Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(trackReal);
    if (!fWeakAndMaterial) {
      if (!fMCStack -> IsPhysicalPrimary(lab)) continue;
      if (!trackReal -> GetInnerParam()) continue;
      if (!fESDtrackCutsLoose -> AcceptTrack(trackReal)) continue; // real track cuts
      if (!(trackReal->GetTPCsignalN()>0)) continue;
    }
    //
    // match the track with mc track
    Int_t iPart = -10;
    if (TMath::Abs(pdg) == kPDGel) iPart = 0; // select el
    if (TMath::Abs(pdg) == kPDGpi) iPart = 1; // select pi
    if (TMath::Abs(pdg) == kPDGka) iPart = 2; // select ka
    if (TMath::Abs(pdg) == kPDGpr) iPart = 3; // select pr
    if (TMath::Abs(pdg) == kPDGde) iPart = 4; // select de
    if (TMath::Abs(pdg) == kPDGmu) iPart = 5; // select mu
    //
    if (iPart == -10) continue;
    //
    if (iPart == 0 ) fElMC = trackReal->GetTPCsignal();
    if (iPart == 1 ) fPiMC = trackReal->GetTPCsignal();
    if (iPart == 2 ) fKaMC = trackReal->GetTPCsignal();
    if (iPart == 3 ) fPrMC = trackReal->GetTPCsignal();
    if (iPart == 4 ) fDeMC = trackReal->GetTPCsignal();
    if (iPart == 5 ) fMuMC = trackReal->GetTPCsignal();
    //
    fEtaMC        = trackReal->Eta();
    fPtotMC       = trackReal->GetInnerParam()->GetP();
    fPtMC         = trackReal->Pt();
    fSignMC       = trackReal->GetSign();
    fTPCSignalMC  = trackReal->GetTPCsignal();
    fPxMC         = trackReal->Px();
    fPyMC         = trackReal->Py();
    fPzMC         = trackReal->Pz();
    fMissingCl    = trackReal->GetTPCClusterInfo(3,0,0,159);
    Float_t pMC   = trackReal->P();
    Float_t fPhiMC= trackReal->Phi();
    trackReal->GetTPCdEdxInfo(tpcdEdxInfo);
    Float_t pv[2],cov[3];
    trackReal->GetImpactParameters(pv,cov); // p[0]=fD; p[1]=fZ; cov[0]=fCdd; cov[1]=fCdz; cov[2]=fCzz;
    Int_t itsNcls = trackReal->GetITSNcls();
    //
    if (trackReal -> GetInnerParam())  {
      Float_t closestPar[3];
      GetExpecteds(trackReal,closestPar);
      SetCutBitsAndSomeTrackVariables(trackReal,iPart);
    }
    //
    // --------------------------------------------------------------
    //                        Fill the trees
    // --------------------------------------------------------------
    //
    // Fill MC closure tree
    if(!fTreeSRedirector) return;
    if (fWeakAndMaterial){
      (*fTreeSRedirector)<<"fTreeMC"<<
      "defCut="    << ifDefaultCuts <<  // default cuts
      "bit96spd="  << fBit96_spd <<  // run Number
      "bit96sdd="  << fBit96_sdd <<  // run Number
      "bit96="     << fBit96_base <<  // run Number
      "bit128="    << fBit128 <<  // run Number
      "bit768="    << fBit768 <<  // run Number
      "pixCut="    << ifDCAcutIfNoITSPixel <<  // run Number
      "dEdx="      << fTPCSignalMC <<    // dEdx of mc track
      "gid="       << fEventGID             <<  //  global event ID
      "part="      << iPart <<
      "ptot="      << fPtotMC <<         // mc momentum
      "p="         << pMC <<
      "pT="        << fPtMC <<         // mc momentum
      "eta="       << fEtaMC <<          // mc eta
      "phi="       << fPhiMC <<          // mc eta
      "cent="      << fCentrality <<     // Centrality
      "centimp="   << fCentImpBin <<
      "sign="      << fSignMC <<         // sign
      "fZ="        << fTrackDCAz <<
      "fD="        << fTrackDCAxy <<
      "fCdd="      << cov[0] <<
      "fCdz="      << cov[1] <<
      "fCzz="      << cov[2] <<
      "missCl="    << fMissingCl <<
      "dEdxInfo.=" << &tpcdEdxInfo <<
      "nMult="     << fNContributors <<
      "tpcMult="   << fTPCMult <<
      "itsNcls="   << itsNcls <<
      "dcaxyprim="      << dcaprim[0] <<           // electron dEdx
      "dcaxyweak="      << dcaweak[0] <<           // electron dEdx
      "dcaxymaterial="  << dcamaterial[0] <<           // electron dEdx
      "dcazprim="       << dcaprim[1] <<           // electron dEdx
      "dcazweak="       << dcaweak[1] <<           // electron dEdx
      "dcazmaterial="   << dcamaterial[1] <<           // electron dEdx
      "\n";
    } else {
      (*fTreeSRedirector)<<"fTreeMC"<<
      "isample="   << sampleNo <<        // sample id for subsample method
      "gid="       << fEventGID <<       //  global event ID
      "part="      << iPart <<
      "dEdx="      << fTPCSignalMC <<    // dEdx of mc track
      "cutBit="    << fTrackCutBits <<   //  Systematic Cuts
      "sign="      << fSignMC <<         // sign
      "ptot="      << fPtotMC <<         // TPC momentum
      "p="         << pMC <<             // momentum at vertex
      "pT="        << fPtMC <<           // mc momentum
      "eta="       << fEtaMC <<          // mc eta
      "phi="       << fPhiMC <<          // mc eta
      "cent="      << fCentrality <<     // Centrality
      "centimp="   << fCentImpBin <<
      //
      "dcaxy="     << fTrackDCAxy           <<
      "dcaz="      << fTrackDCAz            <<
      "ncltpc="    << fNcl                  <<  //  centrality
      "ncltpccorr="<< fNclCorr                  <<  //  centrality
      "cRows="     << fTrackTPCCrossedRows  <<
      "chi2tpc="   << fTrackChi2TPC         <<
      "chi2tpccorr=" << fTrackChi2TPCcorr         <<
      // "intrate="   << fIntRate              <<  // interaction rate
      "\n";
    }

  } // ======= end of track loop for MC dEdx filling =======
  //
  // ======================================================================
  // --------------   MC information with ideal PID   ---------------------
  // ======================================================================
  //
  const Int_t nMoments = 11;
  for (Int_t ieta=0; ieta<fNEtaWinBinsMC; ieta++){
    for (Int_t imom=0; imom<fNMomBinsMC; imom++){

      Int_t nTracksrec=0, trCountrec=0;
      Int_t nTracksgen=0, trCountgen=0;
      // vectors to hold moments
      TVectorF genMoments(nMoments);
      TVectorF genMomentsPos(nMoments);
      TVectorF genMomentsNeg(nMoments);
      TVectorF genMomentsCross(nMoments);
      TVectorF recMoments(nMoments);
      TVectorF recMomentsPos(nMoments);
      TVectorF recMomentsNeg(nMoments);
      TVectorF recMomentsCross(nMoments);
      //
      // Moments without resonances
      TVectorF nRrecMoments(nMoments);
      TVectorF nRrecMomentsPos(nMoments);
      TVectorF nRrecMomentsNeg(nMoments);
      TVectorF nRrecMomentsCross(nMoments);
      TVectorF nRgenMoments(nMoments);
      TVectorF nRgenMomentsPos(nMoments);
      TVectorF nRgenMomentsNeg(nMoments);
      TVectorF nRgenMomentsCross(nMoments);
      //
      // initialize counters
      for(Int_t i=0;i<nMoments; i++){
        genMoments[i]=0.;          nRgenMoments[i]=0.;
        genMomentsPos[i]=0.;       nRgenMomentsPos[i]=0.;
        genMomentsNeg[i]=0.;       nRgenMomentsNeg[i]=0.;
        genMomentsCross[i]=0.;     nRgenMomentsCross[i]=0.;
        recMoments[i]=0.;          nRrecMoments[i]=0.;
        recMomentsPos[i]=0.;       nRrecMomentsPos[i]=0.;
        recMomentsNeg[i]=0.;       nRrecMomentsNeg[i]=0.;
        recMomentsCross[i]=0.;     nRrecMomentsCross[i]=0.;
      }
      //
      // -----------------------------------------------------------------------------------------
      // ----------------------------   reconstructed MC particles  ------------------------------
      // -----------------------------------------------------------------------------------------
      //
      for(Int_t irectrack = 0; irectrack < fESD->GetNumberOfTracks(); irectrack++)
      {
        // track loop
        // initialize the dummy particle id
        fElMC =-100.; fPiMC =-100.; fKaMC =-100.; fPrMC =-100.; fDeMC =-100.; fMuMC =-100.; fLaMC =-100.;
        // Esd track
        fTrackCutBits=0;  // reset the bits for the next track
        AliESDtrack *trackReal = fESD->GetTrack(irectrack);
        if (trackReal==NULL) continue;
        // Get generated track info
        Int_t lab = TMath::Abs(trackReal->GetLabel());           // avoid from negatif labels, they include some garbage
        AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(lab);
        if (!fMCStack->IsPhysicalPrimary(lab)) continue;
        //
        // Select real trigger event and reject other pile up vertices
        Bool_t isTPCPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(lab,fMCEvent);
        Bool_t isITSPileup = AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCEvent, "Hijing");
        if (isTPCPileup || isITSPileup) continue;
        //
        TParticle *trackMC  = fMCStack->Particle(lab);
        Int_t pdg = trackMC->GetPdgCode();
        //
        // acceptance cuts
        Double_t ptotMCrec = trackReal->P();
        Double_t etaMCrec  = trackReal->Eta();
        if (etaMCrec<fetaDownArr[ieta] || etaMCrec>fetaUpArr[ieta]) continue;
        if (ptotMCrec<fpDownArr[imom]  || ptotMCrec>fpUpArr[imom]) continue;
        //
        // Track cuts from detector
        Bool_t ifDefaultCuts = fESDtrackCuts->AcceptTrack(trackReal);
        Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(trackReal);
        if (!trackReal -> GetInnerParam()) continue;
        if (!ifDefaultCuts) continue;
        if (!ifDCAcutIfNoITSPixel) continue;
        //
        // Identify particle wrt pdg code
        Int_t iPart = -10;
        if (TMath::Abs(pdg) == kPDGel) { iPart = 0; fElMC = trackReal->GetTPCsignal(); } // select el-
        if (TMath::Abs(pdg) == kPDGpi) { iPart = 1; fPiMC = trackReal->GetTPCsignal(); } // select pi+
        if (TMath::Abs(pdg) == kPDGka) { iPart = 2; fKaMC = trackReal->GetTPCsignal(); } // select ka+
        if (TMath::Abs(pdg) == kPDGpr) { iPart = 3; fPrMC = trackReal->GetTPCsignal(); } // select pr+
        if (TMath::Abs(pdg) == kPDGde) { iPart = 4; fDeMC = trackReal->GetTPCsignal(); } // select de
        if (TMath::Abs(pdg) == kPDGmu) { iPart = 5; fMuMC = trackReal->GetTPCsignal(); } // select mu-
        if (TMath::Abs(pdg) == kPDGla) { iPart = 6; fLaMC = iPart; } // select Lambda
        if (iPart == -10) continue;
        //
        // additional TOF requirement
        fPtotMC = trackReal->GetInnerParam()->GetP();
        if (fIncludeTOF && fPtotMC>0.8){
          Float_t nSigmasPiTOF = fPIDResponse->NumberOfSigmasTOF(trackReal, AliPID::kPion);
          Float_t nSigmasKaTOF = fPIDResponse->NumberOfSigmasTOF(trackReal, AliPID::kKaon);
          Float_t nSigmasPrTOF = fPIDResponse->NumberOfSigmasTOF(trackReal, AliPID::kProton);
          if ( !(
            ((TMath::Abs(nSigmasPiTOF)<=3) && iPart==1) ||
            ((TMath::Abs(nSigmasKaTOF)<=3) && iPart==2) ||
            ((TMath::Abs(nSigmasPrTOF)<=3) && iPart==3)
          ) ) continue;
        }
        //
        // Resonance control
        Bool_t parInterest = (fPiMC>-1||fKaMC>-1||fPrMC>-1||fElMC>-1||fLaMC>-1) ? kTRUE : kFALSE;
        Bool_t acceptRes = CheckIfFromResonance(0,trackMCgen,lab,parInterest,ptotMCrec,etaMCrec,fCentrality,kFALSE);
        //
        // count first moments for given Centrality and momentum window
        nTracksrec++;
        if ( fPiMC>-1 || fKaMC>-1 || fPrMC>-1 || fLaMC>-1) trCountrec++;
        if ( fPiMC>-1 ) recMoments[kPi]++;
        if ( fKaMC>-1 ) recMoments[kKa]++;
        if ( fPrMC>-1 ) recMoments[kPr]++;
        if ( fPiMC>-1 && pdg<0) recMomentsNeg[kPi]++;
        if ( fKaMC>-1 && pdg<0) recMomentsNeg[kKa]++;
        if ( fPrMC>-1 && pdg<0) recMomentsNeg[kPr]++;
        if ( fPiMC>-1 && pdg>0) recMomentsPos[kPi]++;
        if ( fKaMC>-1 && pdg>0) recMomentsPos[kKa]++;
        if ( fPrMC>-1 && pdg>0) recMomentsPos[kPr]++;

        // Lambdas for alice
        if ( fLaMC>-1 ) recMoments[kLa]++;
        if ( fLaMC>-1 && pdg>0) recMomentsPos[kLa]++;
        if ( fLaMC>-1 && pdg<0) recMomentsNeg[kLa]++;

        if ( acceptRes ) {
          if ( fPiMC>-1 ) nRrecMoments[kPi]++;
          if ( fKaMC>-1 ) nRrecMoments[kKa]++;
          if ( fPrMC>-1 ) nRrecMoments[kPr]++;
          if ( fPiMC>-1 && pdg<0) nRrecMomentsNeg[kPi]++;
          if ( fKaMC>-1 && pdg<0) nRrecMomentsNeg[kKa]++;
          if ( fPrMC>-1 && pdg<0) nRrecMomentsNeg[kPr]++;
          if ( fPiMC>-1 && pdg>0) nRrecMomentsPos[kPi]++;
          if ( fKaMC>-1 && pdg>0) nRrecMomentsPos[kKa]++;
          if ( fPrMC>-1 && pdg>0) nRrecMomentsPos[kPr]++;

          // Lambdas for alice
          if ( fLaMC>-1 ) nRrecMoments[kLa]++;
          if ( fLaMC>-1 && pdg>0) nRrecMomentsPos[kLa]++;
          if ( fLaMC>-1 && pdg<0) nRrecMomentsNeg[kLa]++;
        }

      } // ======= end of track loop =======
      //
      // -----------------------------------------------------------------------------------------
      // ----------------------------   MC generated pure MC particles  --------------------------
      // -----------------------------------------------------------------------------------------
      //
      AliMCParticle *trackMCgen;
      Int_t nMCStackTracks = fMCEvent->GetNumberOfTracks();
      for (Int_t iTrack = 0; iTrack < nMCStackTracks; iTrack++)
      { // track loop
        //
        // Select real trigger event and reject other pile up vertices
        Bool_t isTPCPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iTrack,fMCEvent);
        Bool_t isITSPileup = AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCEvent, "Hijing");
        if (isTPCPileup || isITSPileup) continue;
        //
        // initialize the dummy particle id
        fElMCgen =-100.; fPiMCgen =-100.; fKaMCgen =-100.; fPrMCgen =-100.; fDeMCgen =-100.; fMuMCgen =-100.;
        trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
        Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
        //
        // apply primary track and acceptance cuts
        Double_t ptotMCgen = trackMCgen->P();
        Double_t etaMCgen = (fRapidityType==0) ? trackMCgen->Eta() :  trackMCgen->Y();
        if (etaMCgen<fetaDownArr[ieta]  || etaMCgen>fetaUpArr[ieta]) continue;
        if (ptotMCgen<fpDownArr[imom]   || ptotMCgen>fpUpArr[imom]) continue;
        if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
        //
        // select particle of interest
        Int_t iPart = -10;
        if (TMath::Abs(pdg) == kPDGel) {iPart = 0; fElMCgen = iPart;} // select el-
        if (TMath::Abs(pdg) == kPDGpi) {iPart = 1; fPiMCgen = iPart;} // select pi+
        if (TMath::Abs(pdg) == kPDGka) {iPart = 2; fKaMCgen = iPart;} // select ka+
        if (TMath::Abs(pdg) == kPDGpr) {iPart = 3; fPrMCgen = iPart;} // select pr+
        if (TMath::Abs(pdg) == kPDGde) {iPart = 4; fDeMCgen = iPart;} // select de
        if (TMath::Abs(pdg) == kPDGmu) {iPart = 5; fMuMCgen = iPart;} // select mu-
        if (TMath::Abs(pdg) == kPDGla) {iPart = 6; fLaMCgen = iPart;} // select Lambda
        if (iPart == -10) continue;
        //
        // Resonance control
        Bool_t parInterest = (fPiMCgen>-1||fKaMCgen>-1||fPrMCgen>-1||fElMCgen>-1||fLaMCgen>-1) ? kTRUE : kFALSE;
        Bool_t acceptRes = CheckIfFromResonance(1,trackMCgen,iTrack,parInterest,ptotMCgen,etaMCgen,fCentrality,kTRUE);
        //
        // fill the moments
        nTracksgen++;
        if ( fPiMCgen>-1 || fKaMCgen>-1 || fPrMCgen>-1 || fLaMCgen>-1) trCountgen++;
        if ( fPiMCgen>-1 ) genMoments[kPi]++;
        if ( fKaMCgen>-1 ) genMoments[kKa]++;
        if ( fPrMCgen>-1 ) genMoments[kPr]++;
        if ( fPiMCgen>-1 && pdg<0) genMomentsNeg[kPi]++;
        if ( fKaMCgen>-1 && pdg<0) genMomentsNeg[kKa]++;
        if ( fPrMCgen>-1 && pdg<0) genMomentsNeg[kPr]++;
        if ( fPiMCgen>-1 && pdg>0) genMomentsPos[kPi]++;
        if ( fKaMCgen>-1 && pdg>0) genMomentsPos[kKa]++;
        if ( fPrMCgen>-1 && pdg>0) genMomentsPos[kPr]++;
        //
        // Lambdas for alice
        if ( fLaMCgen>-1 ) genMoments[kLa]++;
        if ( fLaMCgen>-1 && pdg>0) genMomentsPos[kLa]++;
        if ( fLaMCgen>-1 && pdg<0) genMomentsNeg[kLa]++;
        //
        // reject resonances
        if ( acceptRes ) {
          // std::cout << pdg << "   " << parName.GetString() << " ---  " << pdgMom << "  " << momName.GetString() << std::endl;
          if ( fPiMCgen>-1 ) nRgenMoments[kPi]++;
          if ( fKaMCgen>-1 ) nRgenMoments[kKa]++;
          if ( fPrMCgen>-1 ) nRgenMoments[kPr]++;
          if ( fPiMCgen>-1 && pdg<0) nRgenMomentsNeg[kPi]++;
          if ( fKaMCgen>-1 && pdg<0) nRgenMomentsNeg[kKa]++;
          if ( fPrMCgen>-1 && pdg<0) nRgenMomentsNeg[kPr]++;
          if ( fPiMCgen>-1 && pdg>0) nRgenMomentsPos[kPi]++;
          if ( fKaMCgen>-1 && pdg>0) nRgenMomentsPos[kKa]++;
          if ( fPrMCgen>-1 && pdg>0) nRgenMomentsPos[kPr]++;
          //
          // Lambdas for alice
          if ( fLaMCgen>-1 ) nRgenMoments[kLa]++;
          if ( fLaMCgen>-1 && pdg>0) nRgenMomentsPos[kLa]++;
          if ( fLaMCgen>-1 && pdg<0) nRgenMomentsNeg[kLa]++;
        }

      } // ======= end of track loop for generated particles =======
      //
      // -----------------------------------------------------------------------------------------
      // --------------------   Calculation of moments on the event level  -----------------------
      // -----------------------------------------------------------------------------------------
      //
      // ************************************************************************
      //   Moments with resonances
      // ************************************************************************
      //
      // generated level moments                                                  // reconstruced level moments
      genMoments[kPiPi]=genMoments[kPi]*genMoments[kPi];                          recMoments[kPiPi]=recMoments[kPi]*recMoments[kPi];
      genMoments[kKaKa]=genMoments[kKa]*genMoments[kKa];                          recMoments[kKaKa]=recMoments[kKa]*recMoments[kKa];
      genMoments[kPrPr]=genMoments[kPr]*genMoments[kPr];                          recMoments[kPrPr]=recMoments[kPr]*recMoments[kPr];
      genMoments[kPiKa]=genMoments[kPi]*genMoments[kKa];                          recMoments[kPiKa]=recMoments[kPi]*recMoments[kKa];
      genMoments[kPiPr]=genMoments[kPi]*genMoments[kPr];                          recMoments[kPiPr]=recMoments[kPi]*recMoments[kPr];
      genMoments[kKaPr]=genMoments[kKa]*genMoments[kPr];                          recMoments[kKaPr]=recMoments[kKa]*recMoments[kPr];
      genMomentsNeg[kPiPi]=genMomentsNeg[kPi]*genMomentsNeg[kPi];                 recMomentsNeg[kPiPi]=recMomentsNeg[kPi]*recMomentsNeg[kPi];
      genMomentsNeg[kKaKa]=genMomentsNeg[kKa]*genMomentsNeg[kKa];                 recMomentsNeg[kKaKa]=recMomentsNeg[kKa]*recMomentsNeg[kKa];
      genMomentsNeg[kPrPr]=genMomentsNeg[kPr]*genMomentsNeg[kPr];                 recMomentsNeg[kPrPr]=recMomentsNeg[kPr]*recMomentsNeg[kPr];
      genMomentsNeg[kPiKa]=genMomentsNeg[kPi]*genMomentsNeg[kKa];                 recMomentsNeg[kPiKa]=recMomentsNeg[kPi]*recMomentsNeg[kKa];
      genMomentsNeg[kPiPr]=genMomentsNeg[kPi]*genMomentsNeg[kPr];                 recMomentsNeg[kPiPr]=recMomentsNeg[kPi]*recMomentsNeg[kPr];
      genMomentsNeg[kKaPr]=genMomentsNeg[kKa]*genMomentsNeg[kPr];                 recMomentsNeg[kKaPr]=recMomentsNeg[kKa]*recMomentsNeg[kPr];
      genMomentsPos[kPiPi]=genMomentsPos[kPi]*genMomentsPos[kPi];                 recMomentsPos[kPiPi]=recMomentsPos[kPi]*recMomentsPos[kPi];
      genMomentsPos[kKaKa]=genMomentsPos[kKa]*genMomentsPos[kKa];                 recMomentsPos[kKaKa]=recMomentsPos[kKa]*recMomentsPos[kKa];
      genMomentsPos[kPrPr]=genMomentsPos[kPr]*genMomentsPos[kPr];                 recMomentsPos[kPrPr]=recMomentsPos[kPr]*recMomentsPos[kPr];
      genMomentsPos[kPiKa]=genMomentsPos[kPi]*genMomentsPos[kKa];                 recMomentsPos[kPiKa]=recMomentsPos[kPi]*recMomentsPos[kKa];
      genMomentsPos[kPiPr]=genMomentsPos[kPi]*genMomentsPos[kPr];                 recMomentsPos[kPiPr]=recMomentsPos[kPi]*recMomentsPos[kPr];
      genMomentsPos[kKaPr]=genMomentsPos[kKa]*genMomentsPos[kPr];                 recMomentsPos[kKaPr]=recMomentsPos[kKa]*recMomentsPos[kPr];
      genMomentsCross[kPiPosPiNeg]=genMomentsPos[kPi]*genMomentsNeg[kPi];         recMomentsCross[kPiPosPiNeg]=recMomentsPos[kPi]*recMomentsNeg[kPi];
      genMomentsCross[kPiPosKaNeg]=genMomentsPos[kPi]*genMomentsNeg[kKa];         recMomentsCross[kPiPosKaNeg]=recMomentsPos[kPi]*recMomentsNeg[kKa];
      genMomentsCross[kPiPosPrNeg]=genMomentsPos[kPi]*genMomentsNeg[kPr];         recMomentsCross[kPiPosPrNeg]=recMomentsPos[kPi]*recMomentsNeg[kPr];
      genMomentsCross[kKaPosPiNeg]=genMomentsPos[kKa]*genMomentsNeg[kPi];         recMomentsCross[kKaPosPiNeg]=recMomentsPos[kKa]*recMomentsNeg[kPi];
      genMomentsCross[kKaPosKaNeg]=genMomentsPos[kKa]*genMomentsNeg[kKa];         recMomentsCross[kKaPosKaNeg]=recMomentsPos[kKa]*recMomentsNeg[kKa];
      genMomentsCross[kKaPosPrNeg]=genMomentsPos[kKa]*genMomentsNeg[kPr];         recMomentsCross[kKaPosPrNeg]=recMomentsPos[kKa]*recMomentsNeg[kPr];
      genMomentsCross[kPrPosPiNeg]=genMomentsPos[kPr]*genMomentsNeg[kPi];         recMomentsCross[kPrPosPiNeg]=recMomentsPos[kPr]*recMomentsNeg[kPi];
      genMomentsCross[kPrPosKaNeg]=genMomentsPos[kPr]*genMomentsNeg[kKa];         recMomentsCross[kPrPosKaNeg]=recMomentsPos[kPr]*recMomentsNeg[kKa];
      genMomentsCross[kPrPosPrNeg]=genMomentsPos[kPr]*genMomentsNeg[kPr];         recMomentsCross[kPrPosPrNeg]=recMomentsPos[kPr]*recMomentsNeg[kPr];
      // net lambda for Alice
      genMoments[kLaLa]=genMoments[kLa]*genMoments[kLa];                          recMoments[kLaLa]=recMoments[kLa]*recMoments[kLa];
      genMomentsNeg[kLaLa]=genMomentsNeg[kLa]*genMomentsNeg[kLa];                 recMomentsNeg[kLaLa]=recMomentsNeg[kLa]*recMomentsNeg[kLa];
      genMomentsPos[kLaLa]=genMomentsPos[kLa]*genMomentsPos[kLa];                 recMomentsPos[kLaLa]=recMomentsPos[kLa]*recMomentsPos[kLa];
      genMomentsCross[kLaPosLaNeg]=genMomentsPos[kLa]*genMomentsNeg[kLa];         recMomentsCross[kLaPosLaNeg]=recMomentsPos[kLa]*recMomentsNeg[kLa];
      //
      // ************************************************************************
      //   Moments without resonances
      // ************************************************************************
      //
      // Generated level                                                           // Generated level
      nRgenMoments[kPiPi]=nRgenMoments[kPi]*nRgenMoments[kPi];                     nRrecMoments[kPiPi]=nRrecMoments[kPi]*nRrecMoments[kPi];
      nRgenMoments[kKaKa]=nRgenMoments[kKa]*nRgenMoments[kKa];                     nRrecMoments[kKaKa]=nRrecMoments[kKa]*nRrecMoments[kKa];
      nRgenMoments[kPrPr]=nRgenMoments[kPr]*nRgenMoments[kPr];                     nRrecMoments[kPrPr]=nRrecMoments[kPr]*nRrecMoments[kPr];
      nRgenMoments[kPiKa]=nRgenMoments[kPi]*nRgenMoments[kKa];                     nRrecMoments[kPiKa]=nRrecMoments[kPi]*nRrecMoments[kKa];
      nRgenMoments[kPiPr]=nRgenMoments[kPi]*nRgenMoments[kPr];                     nRrecMoments[kPiPr]=nRrecMoments[kPi]*nRrecMoments[kPr];
      nRgenMoments[kKaPr]=nRgenMoments[kKa]*nRgenMoments[kPr];                     nRrecMoments[kKaPr]=nRrecMoments[kKa]*nRrecMoments[kPr];
      nRgenMomentsNeg[kPiPi]=nRgenMomentsNeg[kPi]*nRgenMomentsNeg[kPi];            nRrecMomentsNeg[kPiPi]=nRrecMomentsNeg[kPi]*nRrecMomentsNeg[kPi];
      nRgenMomentsNeg[kKaKa]=nRgenMomentsNeg[kKa]*nRgenMomentsNeg[kKa];            nRrecMomentsNeg[kKaKa]=nRrecMomentsNeg[kKa]*nRrecMomentsNeg[kKa];
      nRgenMomentsNeg[kPrPr]=nRgenMomentsNeg[kPr]*nRgenMomentsNeg[kPr];            nRrecMomentsNeg[kPrPr]=nRrecMomentsNeg[kPr]*nRrecMomentsNeg[kPr];
      nRgenMomentsNeg[kPiKa]=nRgenMomentsNeg[kPi]*nRgenMomentsNeg[kKa];            nRrecMomentsNeg[kPiKa]=nRrecMomentsNeg[kPi]*nRrecMomentsNeg[kKa];
      nRgenMomentsNeg[kPiPr]=nRgenMomentsNeg[kPi]*nRgenMomentsNeg[kPr];            nRrecMomentsNeg[kPiPr]=nRrecMomentsNeg[kPi]*nRrecMomentsNeg[kPr];
      nRgenMomentsNeg[kKaPr]=nRgenMomentsNeg[kKa]*nRgenMomentsNeg[kPr];            nRrecMomentsNeg[kKaPr]=nRrecMomentsNeg[kKa]*nRrecMomentsNeg[kPr];
      nRgenMomentsPos[kPiPi]=nRgenMomentsPos[kPi]*nRgenMomentsPos[kPi];            nRrecMomentsPos[kPiPi]=nRrecMomentsPos[kPi]*nRrecMomentsPos[kPi];
      nRgenMomentsPos[kKaKa]=nRgenMomentsPos[kKa]*nRgenMomentsPos[kKa];            nRrecMomentsPos[kKaKa]=nRrecMomentsPos[kKa]*nRrecMomentsPos[kKa];
      nRgenMomentsPos[kPrPr]=nRgenMomentsPos[kPr]*nRgenMomentsPos[kPr];            nRrecMomentsPos[kPrPr]=nRrecMomentsPos[kPr]*nRrecMomentsPos[kPr];
      nRgenMomentsPos[kPiKa]=nRgenMomentsPos[kPi]*nRgenMomentsPos[kKa];            nRrecMomentsPos[kPiKa]=nRrecMomentsPos[kPi]*nRrecMomentsPos[kKa];
      nRgenMomentsPos[kPiPr]=nRgenMomentsPos[kPi]*nRgenMomentsPos[kPr];            nRrecMomentsPos[kPiPr]=nRrecMomentsPos[kPi]*nRrecMomentsPos[kPr];
      nRgenMomentsPos[kKaPr]=nRgenMomentsPos[kKa]*nRgenMomentsPos[kPr];            nRrecMomentsPos[kKaPr]=nRrecMomentsPos[kKa]*nRrecMomentsPos[kPr];
      nRgenMomentsCross[kPiPosPiNeg]=nRgenMomentsPos[kPi]*nRgenMomentsNeg[kPi];    nRrecMomentsCross[kPiPosPiNeg]=nRrecMomentsPos[kPi]*nRrecMomentsNeg[kPi];
      nRgenMomentsCross[kPiPosKaNeg]=nRgenMomentsPos[kPi]*nRgenMomentsNeg[kKa];    nRrecMomentsCross[kPiPosKaNeg]=nRrecMomentsPos[kPi]*nRrecMomentsNeg[kKa];
      nRgenMomentsCross[kPiPosPrNeg]=nRgenMomentsPos[kPi]*nRgenMomentsNeg[kPr];    nRrecMomentsCross[kPiPosPrNeg]=nRrecMomentsPos[kPi]*nRrecMomentsNeg[kPr];
      nRgenMomentsCross[kKaPosPiNeg]=nRgenMomentsPos[kKa]*nRgenMomentsNeg[kPi];    nRrecMomentsCross[kKaPosPiNeg]=nRrecMomentsPos[kKa]*nRrecMomentsNeg[kPi];
      nRgenMomentsCross[kKaPosKaNeg]=nRgenMomentsPos[kKa]*nRgenMomentsNeg[kKa];    nRrecMomentsCross[kKaPosKaNeg]=nRrecMomentsPos[kKa]*nRrecMomentsNeg[kKa];
      nRgenMomentsCross[kKaPosPrNeg]=nRgenMomentsPos[kKa]*nRgenMomentsNeg[kPr];    nRrecMomentsCross[kKaPosPrNeg]=nRrecMomentsPos[kKa]*nRrecMomentsNeg[kPr];
      nRgenMomentsCross[kPrPosPiNeg]=nRgenMomentsPos[kPr]*nRgenMomentsNeg[kPi];    nRrecMomentsCross[kPrPosPiNeg]=nRrecMomentsPos[kPr]*nRrecMomentsNeg[kPi];
      nRgenMomentsCross[kPrPosKaNeg]=nRgenMomentsPos[kPr]*nRgenMomentsNeg[kKa];    nRrecMomentsCross[kPrPosKaNeg]=nRrecMomentsPos[kPr]*nRrecMomentsNeg[kKa];
      nRgenMomentsCross[kPrPosPrNeg]=nRgenMomentsPos[kPr]*nRgenMomentsNeg[kPr];    nRrecMomentsCross[kPrPosPrNeg]=nRrecMomentsPos[kPr]*nRrecMomentsNeg[kPr];
      // net lambda for Alice
      nRgenMoments[kLaLa]=nRgenMoments[kLa]*nRgenMoments[kLa];                     nRrecMoments[kLaLa]=nRrecMoments[kLa]*nRrecMoments[kLa];
      nRgenMomentsNeg[kLaLa]=nRgenMomentsNeg[kLa]*nRgenMomentsNeg[kLa];            nRrecMomentsNeg[kLaLa]=nRrecMomentsNeg[kLa]*nRrecMomentsNeg[kLa];
      nRgenMomentsPos[kLaLa]=nRgenMomentsPos[kLa]*nRgenMomentsPos[kLa];            nRrecMomentsPos[kLaLa]=nRrecMomentsPos[kLa]*nRrecMomentsPos[kLa];
      nRgenMomentsCross[kLaPosLaNeg]=nRgenMomentsPos[kLa]*nRgenMomentsNeg[kLa];    nRrecMomentsCross[kLaPosLaNeg]=nRrecMomentsPos[kLa]*nRrecMomentsNeg[kLa];
      //
      // fill tree which contains moments
      if(!fTreeSRedirector) return;
      (*fTreeSRedirector)<<"mcFull"<<
      "run="          << runNumber <<               // run number
      "isample="      << sampleNo <<                // sample id for subsample method
      "alltrgen"      << nTracksgen <<
      "alltrrec"      << nTracksrec <<
      "trgen="        << trCountgen <<              // number of identified tracks within the given cent and mom range
      "trrec="        << trCountrec <<                 // number od identified tracks within the given cent and mom range
      "vZ="           << fVz <<
      "cent="         << fCentrality <<
      "centimp="      << fCentImpBin <<
      "impPar="       << fMCImpactParameter <<      // impact parameter taken from MC event header
      "pDown="        << fpDownArr[imom] <<         // lower edge of momentum bin
      "pUp="          << fpUpArr[imom] <<           // upper edge of momentum bin
      "etaDown="      << fetaDownArr[ieta] <<       // lower edge of eta bin
      "etaUp="        << fetaUpArr[ieta] <<         // upper edge of eta bin
      //
      "momentGen.="      << &genMoments <<             // second moments for particle+antiparticle
      "momentPosGen.="   << &genMomentsPos <<          // second moment of positive particles
      "momentNegGen.="   << &genMomentsNeg <<          // second moment of negative particles
      "momentCrossGen.=" << &genMomentsCross <<        // second moment of unlikesign particles
      "noResmomentGen.="      << &nRgenMoments <<             // second moments for particle+antiparticle
      "noResmomentPosGen.="   << &nRgenMomentsPos <<          // second moment of positive particles
      "noResmomentNegGen.="   << &nRgenMomentsNeg <<          // second moment of negative particles
      "noResmomentCrossGen.=" << &nRgenMomentsCross <<        // second moment of unlikesign particles
      "momentRec.="      << &recMoments <<             // second moments for particle+antiparticle
      "momentPosRec.="   << &recMomentsPos <<          // second moment of positive particles
      "momentNegRec.="   << &recMomentsNeg <<          // second moment of negative particles
      "momentCrossRec.=" << &recMomentsCross <<        // second moment of unlikesign particles
      "noResmomentRec.="      << &nRrecMoments <<             // second moments for particle+antiparticle
      "noResmomentPosRec.="   << &nRrecMomentsPos <<          // second moment of positive particles
      "noResmomentNegRec.="   << &nRrecMomentsNeg <<          // second moment of negative particles
      "noResmomentCrossRec.=" << &nRrecMomentsCross <<        // second moment of unlikesign particles
      "\n";
      // tree filling
      //
    }// ======= end of momentum loop =======
  } // ======= end of eta loop =======

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FillMCFull_NetParticles()
{

  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  // Assign subsample index
  Int_t sampleNo = 0;
  Int_t nSubSample = 20;
  sampleNo = Int_t(fEventGID)%nSubSample;
  // Int_t runNumber = (Int_t)fESD->GetRunNumber();
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillMCFull_NetParticles ===== " << std::endl;
  //
  // Count charged primary particles in gen level
  Int_t nChGen = 0;
  for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++)
  {
    AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
    if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
    Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
    Int_t sign = (pdg>0) ? 1 : -1;
    if (TMath::Abs(sign)==1) nChGen++;
  }
  //
  // Get V0 amplitudes
  Int_t multV0A = 0, multV0C = 0;
  for (Int_t i=0;i<32;i++) { multV0A += fESD->GetVZEROData()-> GetMultiplicityV0A(i); }
  for (Int_t i=0;i<32;i++) { multV0C += fESD->GetVZEROData()-> GetMultiplicityV0C(i); }
  //
  // fill mcGen for the nudyn
  if (fFillNudynFastGen) FastGen();
  //
  // reconstructed MC particles with dEdx information
  if (fFillTreeMC) FillTreeMC();
  //
  // Full acceptance kinematic information
  if (fFillDistributions) FillGenDistributions();
  //
  // ======================================================================
  // --------------   MC information with ideal PID   ---------------------
  // ======================================================================
  //
  const Int_t nParticles  = 3;
  const Int_t nMoments    = 14;
  const Int_t nPileUpSettings = 2;
  Int_t nOriginType = (fTrackOriginOnlyPrimary!=1) ? 4 : 1;
  //
  // counters with resonances         // counters without resonances
  TVectorF genPos(nParticles);        TVectorF nRgenPos(nParticles);
  TVectorF genNeg(nParticles);        TVectorF nRgenNeg(nParticles);
  TVectorF fMomNetPiGen(nMoments);    TVectorF fNRMomNetPiGen(nMoments);
  TVectorF fMomNetKaGen(nMoments);    TVectorF fNRMomNetKaGen(nMoments);
  TVectorF fMomNetPrGen(nMoments);    TVectorF fNRMomNetPrGen(nMoments);
  TVectorF recPos(nParticles);        TVectorF nRrecPos(nParticles);
  TVectorF recNeg(nParticles);        TVectorF nRrecNeg(nParticles);
  TVectorF fMomNetPiRec(nMoments);    TVectorF fNRMomNetPiRec(nMoments);
  TVectorF fMomNetKaRec(nMoments);    TVectorF fNRMomNetKaRec(nMoments);
  TVectorF fMomNetPrRec(nMoments);    TVectorF fNRMomNetPrRec(nMoments);
  //
  Bool_t bEventVertexZSmall      = (TMath::Abs(fVz)<6 && TMath::Abs(fVz)>0.2);
  Bool_t bCutReference           = (TMath::Abs(fVz)<7 && TMath::Abs(fVz)>0.15);
  Bool_t bEventVertexZLarge      = (TMath::Abs(fVz)<8 && TMath::Abs(fVz)>0.1);
  Bool_t bEventVertexZALICE      = (TMath::Abs(fVz)<10);
  Bool_t bEventVertexZALICETight = (TMath::Abs(fVz)<7);
  //
  // setting scan
  for (Int_t iset=0; iset<fNSettings; iset++){

    Bool_t mainSettings = (iset==0 || iset==11 || iset==15 || iset==21);    // ????
    if (!mainSettings) continue;
    //
    // event Vz cuts
    if      (iset==kCutEventVertexZSmall      && !bEventVertexZSmall)      continue;
    else if (iset==kCutEventVertexZLarge      && !bEventVertexZLarge)      continue;
    else if (iset==kCutEventVertexZALICE      && !bEventVertexZALICE)      continue;
    else if (iset==kCutEventVertexZALICETight && !bEventVertexZALICETight) continue;
    else if ( !bCutReference ) continue;
    //
    // Acceptance scan
    for (Int_t ieta=0; ieta<fNEtaWinBinsMC; ieta++){
      for (Int_t imom=0; imom<fNMomBinsMC; imom++){
        for (Int_t iorig=0; iorig<nOriginType; iorig++){
          for (Int_t ipileup = 0; ipileup < nPileUpSettings; ipileup++)
          {
            //
            // Initialize counters
            Int_t nTracksrec=0, nTracksgen=0, nTracksTPC = 0, nTracksITS = 0;
            Bool_t isTPCPileup = kFALSE, isITSPileup = kFALSE;
            for(Int_t i=0;i<nParticles; i++){
              genPos[i]=0.;   nRgenPos[i]=0.;
              genNeg[i]=0.;   nRgenNeg[i]=0.;
              recPos[i]=0.;   nRrecPos[i]=0.;
              recNeg[i]=0.;   nRrecNeg[i]=0.;
            }
            for(Int_t i=0;i<nMoments; i++){
              fMomNetPiGen[i]=0.;       fNRMomNetPiGen[i]=0.;
              fMomNetKaGen[i]=0.;       fNRMomNetKaGen[i]=0.;
              fMomNetPrGen[i]=0.;       fNRMomNetPrGen[i]=0.;
              fMomNetPiRec[i]=0.;       fNRMomNetPiRec[i]=0.;
              fMomNetKaRec[i]=0.;       fNRMomNetKaRec[i]=0.;
              fMomNetPrRec[i]=0.;       fNRMomNetPrRec[i]=0.;
            }
            //
            // -----------------------------------------------------------------------------------------
            // ----------------------------   MC generated pure MC particles  --------------------------
            // -----------------------------------------------------------------------------------------
            //
            for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++)
            {
              // track loop
              //
              // Select real trigger event and reject other pile up vertices
              isTPCPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iTrack,fMCEvent);
              isITSPileup = AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCEvent, "Hijing");
              Bool_t ispileup = (isTPCPileup || isITSPileup);
              if (ipileup==0 && ispileup) continue;
              //
              // initialize the dummy particle id
              fElMCgen =-100.; fPiMCgen =-100.; fKaMCgen =-100.; fPrMCgen =-100.;
              AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
              if (!trackMCgen) continue;
              //
              // check the origin of the track
              Bool_t bPrim     = fMCStack->IsPhysicalPrimary(iTrack);
              Bool_t bMaterial = fMCStack->IsSecondaryFromMaterial(iTrack);
              Bool_t bWeak     = fMCStack->IsSecondaryFromWeakDecay(iTrack);
              Bool_t bAcceptOrigin = kFALSE;
              if (iorig==0) bAcceptOrigin = bPrim;
              if (iorig==1) bAcceptOrigin = (bPrim || bWeak);
              if (iorig==2) bAcceptOrigin = (bPrim || bMaterial);
              if (iorig==3) bAcceptOrigin = (bPrim || bMaterial || bWeak);
              if (!bAcceptOrigin) continue;
              //
              // select particle of interest
              Int_t iPart = -10;
              Int_t pdg = trackMCgen->Particle()->GetPdgCode();
              if (TMath::Abs(pdg) == kPDGpi) {iPart = 0; fPiMCgen = iPart;} // select pi+
              if (TMath::Abs(pdg) == kPDGka) {iPart = 1; fKaMCgen = iPart;} // select ka+
              if (TMath::Abs(pdg) == kPDGpr) {iPart = 2; fPrMCgen = iPart;} // select pr+
              if (iPart == -10) continue; // perfect PID cut
              //
              Double_t ptotMCgen = 0.;
              if(fUsePtCut==0) ptotMCgen = trackMCgen->P();
              if(fUsePtCut==1) ptotMCgen = trackMCgen->P();
              if(fUsePtCut==2) ptotMCgen = trackMCgen->Pt();
              Float_t phiMCGen  = trackMCgen->Phi();
              Double_t etaMCgen = trackMCgen->Eta();
              Bool_t etaAcc  = (etaMCgen>=fetaDownArr[ieta] && etaMCgen<=fetaUpArr[ieta]);
              Bool_t momAcc  = (ptotMCgen>=fpDownArr[imom]  && ptotMCgen<fpUpArr[imom]);
              Bool_t etaAccMaxWindow = (etaMCgen>=fEtaDown  && etaMCgen<=fEtaUp);
              Bool_t momAccMaxWindow = (ptotMCgen>=fMomDown && ptotMCgen<=fMomUp);
              //
              // Fill eff Matrix
              if (fEffMatrix && !fRunOnGrid && etaAccMaxWindow && momAccMaxWindow && ieta==0 && imom==0 && iorig==0 && ispileup==0){
                //
                // Eff matrix phi etc.
                if (iset==kCutReference){
                  Double_t xxxGen[5]={Float_t(iPart),fCentrality,ptotMCgen,etaMCgen,phiMCGen};
                  if (pdg>0) fHistPosEffMatrixGen->Fill(xxxGen);
                  if (pdg<0) fHistNegEffMatrixGen->Fill(xxxGen);
                }
                //
                // systematic scan
                Double_t xxxGenSystScan[6]={0.,Float_t(iset),Float_t(iPart),fCentrality,ptotMCgen,etaMCgen};
                if (pdg>0) fHistPosEffMatrixScanGen->Fill(xxxGenSystScan);
                if (pdg<0) fHistNegEffMatrixScanGen->Fill(xxxGenSystScan);
                // generated does not know about TOF
                Double_t xxxGenSystScanTOF[6]={1.,Float_t(iset),Float_t(iPart),fCentrality,ptotMCgen,etaMCgen};
                if (pdg>0) fHistPosEffMatrixScanGen->Fill(xxxGenSystScanTOF);
                if (pdg<0) fHistNegEffMatrixScanGen->Fill(xxxGenSystScanTOF);
              }
              //
              // fill the moments
              if (etaAcc && momAcc){
                nTracksgen++;
                if ( fPiMCgen>-1 && pdg<0) genNeg[kPi]++;
                if ( fKaMCgen>-1 && pdg<0) genNeg[kKa]++;
                if ( fPrMCgen>-1 && pdg<0) genNeg[kPr]++;
                if ( fPiMCgen>-1 && pdg>0) genPos[kPi]++;
                if ( fKaMCgen>-1 && pdg>0) genPos[kKa]++;
                if ( fPrMCgen>-1 && pdg>0) genPos[kPr]++;
                //
                // reject resonances
                Bool_t acceptRes = CheckIfFromAnyResonance(trackMCgen,fetaDownArr[ieta],fetaUpArr[ieta],fpDownArr[imom],fpUpArr[imom]);
                if ( acceptRes ) {
                  if ( fPiMCgen>-1 && pdg<0) nRgenNeg[kPi]++;
                  if ( fKaMCgen>-1 && pdg<0) nRgenNeg[kKa]++;
                  if ( fPrMCgen>-1 && pdg<0) nRgenNeg[kPr]++;
                  if ( fPiMCgen>-1 && pdg>0) nRgenPos[kPi]++;
                  if ( fKaMCgen>-1 && pdg>0) nRgenPos[kKa]++;
                  if ( fPrMCgen>-1 && pdg>0) nRgenPos[kPr]++;
                }
              }

            } // ======= end of track loop for generated particles =======
            //
            // -----------------------------------------------------------------------------------------
            // ----------------------------   reconstructed MC particles  ------------------------------
            // -----------------------------------------------------------------------------------------
            //
            for(Int_t irectrack = 0; irectrack < fESD->GetNumberOfTracks(); irectrack++)
            {
              // track loop
              //
              // select pile up
              Bool_t ispileup = (isTPCPileup || isITSPileup);
              if (ipileup==0 && ispileup) continue;
              //
              // initialize the dummy particle id
              fElMC =-100.; fPiMC =-100.; fKaMC =-100.; fPrMC =-100.;
              // Esd track
              fTrackCutBits=0;  // reset the bits for the next track
              AliESDtrack *trackReal = fESD->GetTrack(irectrack);
              if (trackReal==NULL) continue;
              Int_t lab = TMath::Abs(trackReal->GetLabel());           // avoid from negatif labels, they include some garbage
              //
              // check the origin of the track
              Bool_t bPrim     = fMCStack->IsPhysicalPrimary(lab);
              Bool_t bMaterial = fMCStack->IsSecondaryFromMaterial(lab);
              Bool_t bWeak     = fMCStack->IsSecondaryFromWeakDecay(lab);
              Bool_t bAcceptOrigin = kFALSE;
              if (iorig==0) bAcceptOrigin = bPrim;
              if (iorig==1) bAcceptOrigin = (bPrim || bWeak);
              if (iorig==2) bAcceptOrigin = (bPrim || bMaterial);
              if (iorig==3) bAcceptOrigin = (bPrim || bMaterial || bWeak);
              if (!bAcceptOrigin) continue;
              //
              // Identify particle wrt pdg code
              AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(lab); // TParticle *trackMC  = fMCStack->Particle(lab);
              Int_t pdg = trackMCgen->Particle()->GetPdgCode();                     // Int_t pdg = trackMC->GetPdgCode();
              //
              Int_t iPart = -10;
              if (TMath::Abs(pdg) == kPDGpi) { iPart = 0; fPiMC = trackReal->GetTPCsignal(); } // select pi+
              if (TMath::Abs(pdg) == kPDGka) { iPart = 1; fKaMC = trackReal->GetTPCsignal(); } // select ka+
              if (TMath::Abs(pdg) == kPDGpr) { iPart = 2; fPrMC = trackReal->GetTPCsignal(); } // select pr+
              if (iPart == -10) continue; // perfect PID cut
              //
              // apply detector cuts
              Bool_t ifDefaultCuts = fESDtrackCuts->AcceptTrack(trackReal);
              Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(trackReal);
              if (!trackReal-> GetInnerParam()) continue;
              if (!fESDtrackCutsLoose->AcceptTrack(trackReal))  continue;    // Loose Cuts
              if (!(trackReal->GetTPCsignalN()>0)) continue;
              if (!ifDefaultCuts) continue;
              if (!ifDCAcutIfNoITSPixel) continue;
              //
              // acceptance cuts
              Double_t ptotMCrec = 0.;
              if(fUsePtCut==0) ptotMCrec = trackReal->GetInnerParam()->GetP();
              if(fUsePtCut==1) ptotMCrec = trackReal->P();
              if(fUsePtCut==2) ptotMCrec = trackReal->Pt();
              Float_t phiMCRec  = trackReal->Phi();
              Double_t etaMCrec = trackReal->Eta();
              Bool_t etaAcc  = (etaMCrec>=fetaDownArr[ieta] && etaMCrec<fetaUpArr[ieta]);
              Bool_t momAcc  = (ptotMCrec>=fpDownArr[imom]  && ptotMCrec<fpUpArr[imom]);
              Bool_t etaAccMaxWindow = (etaMCrec>=fEtaDown  && etaMCrec<=fEtaUp);
              Bool_t momAccMaxWindow = (ptotMCrec>=fMomDown && ptotMCrec<=fMomUp);
              //
              // Get the cut bit information apply track cuts
              Float_t closestPar[3];
              GetExpecteds(trackReal,closestPar);
              SetCutBitsAndSomeTrackVariables(trackReal,iPart);
              if (!GetSystematicClassIndex(fTrackCutBits,iset)) continue;
              //
              // Fill Efficiency Matrices
              if (fEffMatrix && !fRunOnGrid && etaAccMaxWindow && momAccMaxWindow && ieta==0 && imom==0 && iorig==0 && ispileup==0){
                //
                // Eff matrix with phi etc. for the default setting
                if (iset==kCutReference){
                  Double_t xxxRec[5]={Float_t(iPart),fCentrality,ptotMCrec,etaMCrec,phiMCRec};
                  if (pdg>0) fHistPosEffMatrixRec->Fill(xxxRec);
                  if (pdg<0) fHistNegEffMatrixRec->Fill(xxxRec);
                }
                //
                // TPC eff matrix for all settings
                Double_t xxxRecSystScan[6]={0.,Float_t(iset),Float_t(iPart),fCentrality,ptotMCrec,etaMCrec};
                if (pdg>0) fHistPosEffMatrixScanRec->Fill(xxxRecSystScan);
                if (pdg<0) fHistNegEffMatrixScanRec->Fill(xxxRecSystScan);
                //
                // TOF+TPC eff matrix for all settings
                Bool_t piTOF = (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trackReal, AliPID::kPion))  <=2.5);
                Bool_t kaTOF = (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trackReal, AliPID::kKaon))  <=2.5);
                Bool_t prTOF = (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trackReal, AliPID::kProton))<=2.5);
                if ( (piTOF && iPart==0) ||  (kaTOF && iPart==1) || (prTOF && iPart==2) ) {
                  Double_t xxxRecTOF[6]={1.,Float_t(iset),Float_t(iPart),fCentrality,ptotMCrec,etaMCrec};
                  if (pdg>0) fHistPosEffMatrixScanRec->Fill(xxxRecTOF);
                  if (pdg<0) fHistNegEffMatrixScanRec->Fill(xxxRecTOF);
                }
              }
              //
              // count first moments for given Centrality and momentum window
              if (etaAcc && momAcc){
                //
                if (trackReal->IsOn(AliESDtrack::kTPCrefit)) nTracksTPC++;
                if (trackReal->IsOn(AliESDtrack::kITSrefit)) nTracksITS++;
                nTracksrec++;
                if ( fPiMC>-1 && pdg<0) recNeg[kPi]++;
                if ( fKaMC>-1 && pdg<0) recNeg[kKa]++;
                if ( fPrMC>-1 && pdg<0) recNeg[kPr]++;
                if ( fPiMC>-1 && pdg>0) recPos[kPi]++;
                if ( fKaMC>-1 && pdg>0) recPos[kKa]++;
                if ( fPrMC>-1 && pdg>0) recPos[kPr]++;
                //
                // count first moments for given Centrality and momentum window without resonances
                Bool_t acceptRes = CheckIfFromAnyResonance(trackMCgen,fetaDownArr[ieta],fetaUpArr[ieta],fpDownArr[imom],fpUpArr[imom]);
                if ( acceptRes ) {
                  if ( fPiMC>-1 && pdg<0) nRrecNeg[kPi]++;
                  if ( fKaMC>-1 && pdg<0) nRrecNeg[kKa]++;
                  if ( fPrMC>-1 && pdg<0) nRrecNeg[kPr]++;
                  if ( fPiMC>-1 && pdg>0) nRrecPos[kPi]++;
                  if ( fKaMC>-1 && pdg>0) nRrecPos[kKa]++;
                  if ( fPrMC>-1 && pdg>0) nRrecPos[kPr]++;
                }
              }

            } // ======= end of track loop =======
            //
            // -----------------------------------------------------------------------------------------
            // --------------------   Calculation of moments on the event level  -----------------------
            // -----------------------------------------------------------------------------------------
            //
            // ************************************************************************
            //   Moments with resonances
            // ************************************************************************
            //
            // Generated with reosnances                                                    // Reconstructed with reosnances
            // Net Pions
            fMomNetPiGen[kA]    = genPos[kPi];                                              fMomNetPiRec[kA]    = recPos[kPi];
            fMomNetPiGen[kB]    = genNeg[kPi];                                              fMomNetPiRec[kB]    = recNeg[kPi];
            fMomNetPiGen[kAA]   = genPos[kPi]*genPos[kPi];                                  fMomNetPiRec[kAA]   = recPos[kPi]*recPos[kPi];
            fMomNetPiGen[kBB]   = genNeg[kPi]*genNeg[kPi];                                  fMomNetPiRec[kBB]   = recNeg[kPi]*recNeg[kPi];
            fMomNetPiGen[kAB]   = genPos[kPi]*genNeg[kPi];                                  fMomNetPiRec[kAB]   = recPos[kPi]*recNeg[kPi];
            fMomNetPiGen[kAAA]  = genPos[kPi]*genPos[kPi]*genPos[kPi];                      fMomNetPiRec[kAAA]  = recPos[kPi]*recPos[kPi]*recPos[kPi];
            fMomNetPiGen[kBBB]  = genNeg[kPi]*genNeg[kPi]*genNeg[kPi];                      fMomNetPiRec[kBBB]  = recNeg[kPi]*recNeg[kPi]*recNeg[kPi];
            fMomNetPiGen[kAAB]  = genPos[kPi]*genPos[kPi]*genNeg[kPi];                      fMomNetPiRec[kAAB]  = recPos[kPi]*recPos[kPi]*recNeg[kPi];
            fMomNetPiGen[kBBA]  = genNeg[kPi]*genNeg[kPi]*genPos[kPi];                      fMomNetPiRec[kBBA]  = recNeg[kPi]*recNeg[kPi]*recPos[kPi];
            fMomNetPiGen[kABBB] = genPos[kPi]*genNeg[kPi]*genNeg[kPi]*genNeg[kPi];          fMomNetPiRec[kABBB] = recPos[kPi]*recNeg[kPi]*recNeg[kPi]*recNeg[kPi];
            fMomNetPiGen[kAABB] = genPos[kPi]*genPos[kPi]*genNeg[kPi]*genNeg[kPi];          fMomNetPiRec[kAABB] = recPos[kPi]*recPos[kPi]*recNeg[kPi]*recNeg[kPi];
            fMomNetPiGen[kAAAB] = genPos[kPi]*genPos[kPi]*genPos[kPi]*genNeg[kPi];          fMomNetPiRec[kAAAB] = recPos[kPi]*recPos[kPi]*recPos[kPi]*recNeg[kPi];
            fMomNetPiGen[kAAAA] = genPos[kPi]*genPos[kPi]*genPos[kPi]*genPos[kPi];          fMomNetPiRec[kAAAA] = recPos[kPi]*recPos[kPi]*recPos[kPi]*recPos[kPi];
            fMomNetPiGen[kBBBB] = genNeg[kPi]*genNeg[kPi]*genNeg[kPi]*genNeg[kPi];          fMomNetPiRec[kBBBB] = recNeg[kPi]*recNeg[kPi]*recNeg[kPi]*recNeg[kPi];

            //
            // Net Kaons
            fMomNetKaGen[kA]    = genPos[kKa];                                              fMomNetKaRec[kA]    = recPos[kKa];
            fMomNetKaGen[kB]    = genNeg[kKa];                                              fMomNetKaRec[kB]    = recNeg[kKa];
            fMomNetKaGen[kAA]   = genPos[kKa]*genPos[kKa];                                  fMomNetKaRec[kAA]   = recPos[kKa]*recPos[kKa];
            fMomNetKaGen[kBB]   = genNeg[kKa]*genNeg[kKa];                                  fMomNetKaRec[kBB]   = recNeg[kKa]*recNeg[kKa];
            fMomNetKaGen[kAB]   = genPos[kKa]*genNeg[kKa];                                  fMomNetKaRec[kAB]   = recPos[kKa]*recNeg[kKa];
            fMomNetKaGen[kAAA]  = genPos[kKa]*genPos[kKa]*genPos[kKa];                      fMomNetKaRec[kAAA]  = recPos[kKa]*recPos[kKa]*recPos[kKa];
            fMomNetKaGen[kBBB]  = genNeg[kKa]*genNeg[kKa]*genNeg[kKa];                      fMomNetKaRec[kBBB]  = recNeg[kKa]*recNeg[kKa]*recNeg[kKa];
            fMomNetKaGen[kAAB]  = genPos[kKa]*genPos[kKa]*genNeg[kKa];                      fMomNetKaRec[kAAB]  = recPos[kKa]*recPos[kKa]*recNeg[kKa];
            fMomNetKaGen[kBBA]  = genNeg[kKa]*genNeg[kKa]*genPos[kKa];                      fMomNetKaRec[kBBA]  = recNeg[kKa]*recNeg[kKa]*recPos[kKa];
            fMomNetKaGen[kABBB] = genPos[kKa]*genNeg[kKa]*genNeg[kKa]*genNeg[kKa];          fMomNetKaRec[kABBB] = recPos[kKa]*recNeg[kKa]*recNeg[kKa]*recNeg[kKa];
            fMomNetKaGen[kAABB] = genPos[kKa]*genPos[kKa]*genNeg[kKa]*genNeg[kKa];          fMomNetKaRec[kAABB] = recPos[kKa]*recPos[kKa]*recNeg[kKa]*recNeg[kKa];
            fMomNetKaGen[kAAAB] = genPos[kKa]*genPos[kKa]*genPos[kKa]*genNeg[kKa];          fMomNetKaRec[kAAAB] = recPos[kKa]*recPos[kKa]*recPos[kKa]*recNeg[kKa];
            fMomNetKaGen[kAAAA] = genPos[kKa]*genPos[kKa]*genPos[kKa]*genPos[kKa];          fMomNetKaRec[kAAAA] = recPos[kKa]*recPos[kKa]*recPos[kKa]*recPos[kKa];
            fMomNetKaGen[kBBBB] = genNeg[kKa]*genNeg[kKa]*genNeg[kKa]*genNeg[kKa];          fMomNetKaRec[kBBBB] = recNeg[kKa]*recNeg[kKa]*recNeg[kKa]*recNeg[kKa];

            //
            // Net Protons
            fMomNetPrGen[kA]    = genPos[kPr];                                              fMomNetPrRec[kA]    = recPos[kPr];
            fMomNetPrGen[kB]    = genNeg[kPr];                                              fMomNetPrRec[kB]    = recNeg[kPr];
            fMomNetPrGen[kAA]   = genPos[kPr]*genPos[kPr];                                  fMomNetPrRec[kAA]   = recPos[kPr]*recPos[kPr];
            fMomNetPrGen[kBB]   = genNeg[kPr]*genNeg[kPr];                                  fMomNetPrRec[kBB]   = recNeg[kPr]*recNeg[kPr];
            fMomNetPrGen[kAB]   = genPos[kPr]*genNeg[kPr];                                  fMomNetPrRec[kAB]   = recPos[kPr]*recNeg[kPr];
            fMomNetPrGen[kAAA]  = genPos[kPr]*genPos[kPr]*genPos[kPr];                      fMomNetPrRec[kAAA]  = recPos[kPr]*recPos[kPr]*recPos[kPr];
            fMomNetPrGen[kBBB]  = genNeg[kPr]*genNeg[kPr]*genNeg[kPr];                      fMomNetPrRec[kBBB]  = recNeg[kPr]*recNeg[kPr]*recNeg[kPr];
            fMomNetPrGen[kAAB]  = genPos[kPr]*genPos[kPr]*genNeg[kPr];                      fMomNetPrRec[kAAB]  = recPos[kPr]*recPos[kPr]*recNeg[kPr];
            fMomNetPrGen[kBBA]  = genNeg[kPr]*genNeg[kPr]*genPos[kPr];                      fMomNetPrRec[kBBA]  = recNeg[kPr]*recNeg[kPr]*recPos[kPr];
            fMomNetPrGen[kABBB] = genPos[kPr]*genNeg[kPr]*genNeg[kPr]*genNeg[kPr];          fMomNetPrRec[kABBB] = recPos[kPr]*recNeg[kPr]*recNeg[kPr]*recNeg[kPr];
            fMomNetPrGen[kAABB] = genPos[kPr]*genPos[kPr]*genNeg[kPr]*genNeg[kPr];          fMomNetPrRec[kAABB] = recPos[kPr]*recPos[kPr]*recNeg[kPr]*recNeg[kPr];
            fMomNetPrGen[kAAAB] = genPos[kPr]*genPos[kPr]*genPos[kPr]*genNeg[kPr];          fMomNetPrRec[kAAAB] = recPos[kPr]*recPos[kPr]*recPos[kPr]*recNeg[kPr];
            fMomNetPrGen[kAAAA] = genPos[kPr]*genPos[kPr]*genPos[kPr]*genPos[kPr];          fMomNetPrRec[kAAAA] = recPos[kPr]*recPos[kPr]*recPos[kPr]*recPos[kPr];
            fMomNetPrGen[kBBBB] = genNeg[kPr]*genNeg[kPr]*genNeg[kPr]*genNeg[kPr];          fMomNetPrRec[kBBBB] = recNeg[kPr]*recNeg[kPr]*recNeg[kPr]*recNeg[kPr];
            //
            // ************************************************************************
            //   Moments without resonances
            // ************************************************************************
            //
            // Generated with reosnances
            //
            // Net Pions
            fNRMomNetPiGen[kA]    = nRgenPos[kPi];                                                    fNRMomNetPiRec[kA]   = nRrecPos[kPi];
            fNRMomNetPiGen[kB]    = nRgenNeg[kPi];                                                    fNRMomNetPiRec[kB]   = nRrecNeg[kPi];
            fNRMomNetPiGen[kAA]   = nRgenPos[kPi]*nRgenPos[kPi];                                      fNRMomNetPiRec[kAA]  = nRrecPos[kPi]*nRrecPos[kPi];
            fNRMomNetPiGen[kBB]   = nRgenNeg[kPi]*nRgenNeg[kPi];                                      fNRMomNetPiRec[kBB]  = nRrecNeg[kPi]*nRrecNeg[kPi];
            fNRMomNetPiGen[kAB]   = nRgenPos[kPi]*nRgenNeg[kPi];                                      fNRMomNetPiRec[kAB]  = nRrecPos[kPi]*nRrecNeg[kPi];
            fNRMomNetPiGen[kAAA]  = nRgenPos[kPi]*nRgenPos[kPi]*nRgenPos[kPi];                        fNRMomNetPiRec[kAAA] = nRrecPos[kPi]*nRrecPos[kPi]*nRrecPos[kPi];
            fNRMomNetPiGen[kBBB]  = nRgenNeg[kPi]*nRgenNeg[kPi]*nRgenNeg[kPi];                        fNRMomNetPiRec[kBBB] = nRrecNeg[kPi]*nRrecNeg[kPi]*nRrecNeg[kPi];
            fNRMomNetPiGen[kAAB]  = nRgenPos[kPi]*nRgenPos[kPi]*nRgenNeg[kPi];                        fNRMomNetPiRec[kAAB] = nRrecPos[kPi]*nRrecPos[kPi]*nRrecNeg[kPi];
            fNRMomNetPiGen[kBBA]  = nRgenNeg[kPi]*nRgenNeg[kPi]*nRgenPos[kPi];                        fNRMomNetPiRec[kBBA] = nRrecNeg[kPi]*nRrecNeg[kPi]*nRrecPos[kPi];
            fNRMomNetPiGen[kABBB] = nRgenPos[kPi]*nRgenNeg[kPi]*nRgenNeg[kPi]*nRgenNeg[kPi];          fNRMomNetPiRec[kABBB] = nRrecPos[kPi]*nRrecNeg[kPi]*nRrecNeg[kPi]*nRrecNeg[kPi];
            fNRMomNetPiGen[kAABB] = nRgenPos[kPi]*nRgenPos[kPi]*nRgenNeg[kPi]*nRgenNeg[kPi];          fNRMomNetPiRec[kAABB] = nRrecPos[kPi]*nRrecPos[kPi]*nRrecNeg[kPi]*nRrecNeg[kPi];
            fNRMomNetPiGen[kAAAB] = nRgenPos[kPi]*nRgenPos[kPi]*nRgenPos[kPi]*nRgenNeg[kPi];          fNRMomNetPiRec[kAAAB] = nRrecPos[kPi]*nRrecPos[kPi]*nRrecPos[kPi]*nRrecNeg[kPi];
            fNRMomNetPiGen[kAAAA] = nRgenPos[kPi]*nRgenPos[kPi]*nRgenPos[kPi]*nRgenPos[kPi];          fNRMomNetPiRec[kAAAA] = nRrecPos[kPi]*nRrecPos[kPi]*nRrecPos[kPi]*nRrecPos[kPi];
            fNRMomNetPiGen[kBBBB] = nRgenNeg[kPi]*nRgenNeg[kPi]*nRgenNeg[kPi]*nRgenNeg[kPi];          fNRMomNetPiRec[kBBBB] = nRrecNeg[kPi]*nRrecNeg[kPi]*nRrecNeg[kPi]*nRrecNeg[kPi];

            //
            // Net Kaons
            fNRMomNetKaGen[kA]    = nRgenPos[kKa];                                                    fNRMomNetKaRec[kA]   = nRrecPos[kKa];
            fNRMomNetKaGen[kB]    = nRgenNeg[kKa];                                                    fNRMomNetKaRec[kB]   = nRrecNeg[kKa];
            fNRMomNetKaGen[kAA]   = nRgenPos[kKa]*nRgenPos[kKa];                                      fNRMomNetKaRec[kAA]  = nRrecPos[kKa]*nRrecPos[kKa];
            fNRMomNetKaGen[kBB]   = nRgenNeg[kKa]*nRgenNeg[kKa];                                      fNRMomNetKaRec[kBB]  = nRrecNeg[kKa]*nRrecNeg[kKa];
            fNRMomNetKaGen[kAB]   = nRgenPos[kKa]*nRgenNeg[kKa];                                      fNRMomNetKaRec[kAB]  = nRrecPos[kKa]*nRrecNeg[kKa];
            fNRMomNetKaGen[kAAA]  = nRgenPos[kKa]*nRgenPos[kKa]*nRgenPos[kKa];                        fNRMomNetKaRec[kAAA] = nRrecPos[kKa]*nRrecPos[kKa]*nRrecPos[kKa];
            fNRMomNetKaGen[kBBB]  = nRgenNeg[kKa]*nRgenNeg[kKa]*nRgenNeg[kKa];                        fNRMomNetKaRec[kBBB] = nRrecNeg[kKa]*nRrecNeg[kKa]*nRrecNeg[kKa];
            fNRMomNetKaGen[kAAB]  = nRgenPos[kKa]*nRgenPos[kKa]*nRgenNeg[kKa];                        fNRMomNetKaRec[kAAB] = nRrecPos[kKa]*nRrecPos[kKa]*nRrecNeg[kKa];
            fNRMomNetKaGen[kBBA]  = nRgenNeg[kKa]*nRgenNeg[kKa]*nRgenPos[kKa];                        fNRMomNetKaRec[kBBA] = nRrecNeg[kKa]*nRrecNeg[kKa]*nRrecPos[kKa];
            fNRMomNetKaGen[kABBB] = nRgenPos[kKa]*nRgenNeg[kKa]*nRgenNeg[kKa]*nRgenNeg[kKa];          fNRMomNetKaRec[kABBB] = nRrecPos[kKa]*nRrecNeg[kKa]*nRrecNeg[kKa]*nRrecNeg[kKa];
            fNRMomNetKaGen[kAABB] = nRgenPos[kKa]*nRgenPos[kKa]*nRgenNeg[kKa]*nRgenNeg[kKa];          fNRMomNetKaRec[kAABB] = nRrecPos[kKa]*nRrecPos[kKa]*nRrecNeg[kKa]*nRrecNeg[kKa];
            fNRMomNetKaGen[kAAAB] = nRgenPos[kKa]*nRgenPos[kKa]*nRgenPos[kKa]*nRgenNeg[kKa];          fNRMomNetKaRec[kAAAB] = nRrecPos[kKa]*nRrecPos[kKa]*nRrecPos[kKa]*nRrecNeg[kKa];
            fNRMomNetKaGen[kAAAA] = nRgenPos[kKa]*nRgenPos[kKa]*nRgenPos[kKa]*nRgenPos[kKa];          fNRMomNetKaRec[kAAAA] = nRrecPos[kKa]*nRrecPos[kKa]*nRrecPos[kKa]*nRrecPos[kKa];
            fNRMomNetKaGen[kBBBB] = nRgenNeg[kKa]*nRgenNeg[kKa]*nRgenNeg[kKa]*nRgenNeg[kKa];          fNRMomNetKaRec[kBBBB] = nRrecNeg[kKa]*nRrecNeg[kKa]*nRrecNeg[kKa]*nRrecNeg[kKa];

            //
            // Net Protons
            fNRMomNetPrGen[kA]    = nRgenPos[kPr];                                                    fNRMomNetPrRec[kA]   = nRrecPos[kPr];
            fNRMomNetPrGen[kB]    = nRgenNeg[kPr];                                                    fNRMomNetPrRec[kB]   = nRrecNeg[kPr];
            fNRMomNetPrGen[kAA]   = nRgenPos[kPr]*nRgenPos[kPr];                                      fNRMomNetPrRec[kAA]  = nRrecPos[kPr]*nRrecPos[kPr];
            fNRMomNetPrGen[kBB]   = nRgenNeg[kPr]*nRgenNeg[kPr];                                      fNRMomNetPrRec[kBB]  = nRrecNeg[kPr]*nRrecNeg[kPr];
            fNRMomNetPrGen[kAB]   = nRgenPos[kPr]*nRgenNeg[kPr];                                      fNRMomNetPrRec[kAB]  = nRrecPos[kPr]*nRrecNeg[kPr];
            fNRMomNetPrGen[kAAA]  = nRgenPos[kPr]*nRgenPos[kPr]*nRgenPos[kPr];                        fNRMomNetPrRec[kAAA] = nRrecPos[kPr]*nRrecPos[kPr]*nRrecPos[kPr];
            fNRMomNetPrGen[kBBB]  = nRgenNeg[kPr]*nRgenNeg[kPr]*nRgenNeg[kPr];                        fNRMomNetPrRec[kBBB] = nRrecNeg[kPr]*nRrecNeg[kPr]*nRrecNeg[kPr];
            fNRMomNetPrGen[kAAB]  = nRgenPos[kPr]*nRgenPos[kPr]*nRgenNeg[kPr];                        fNRMomNetPrRec[kAAB] = nRrecPos[kPr]*nRrecPos[kPr]*nRrecNeg[kPr];
            fNRMomNetPrGen[kBBA]  = nRgenNeg[kPr]*nRgenNeg[kPr]*nRgenPos[kPr];                        fNRMomNetPrRec[kBBA] = nRrecNeg[kPr]*nRrecNeg[kPr]*nRrecPos[kPr];
            fNRMomNetPrGen[kABBB] = nRgenPos[kPr]*nRgenNeg[kPr]*nRgenNeg[kPr]*nRgenNeg[kPr];          fNRMomNetPrRec[kABBB] = nRrecPos[kPr]*nRrecNeg[kPr]*nRrecNeg[kPr]*nRrecNeg[kPr];
            fNRMomNetPrGen[kAABB] = nRgenPos[kPr]*nRgenPos[kPr]*nRgenNeg[kPr]*nRgenNeg[kPr];          fNRMomNetPrRec[kAABB] = nRrecPos[kPr]*nRrecPos[kPr]*nRrecNeg[kPr]*nRrecNeg[kPr];
            fNRMomNetPrGen[kAAAB] = nRgenPos[kPr]*nRgenPos[kPr]*nRgenPos[kPr]*nRgenNeg[kPr];          fNRMomNetPrRec[kAAAB] = nRrecPos[kPr]*nRrecPos[kPr]*nRrecPos[kPr]*nRrecNeg[kPr];
            fNRMomNetPrGen[kAAAA] = nRgenPos[kPr]*nRgenPos[kPr]*nRgenPos[kPr]*nRgenPos[kPr];          fNRMomNetPrRec[kAAAA] = nRrecPos[kPr]*nRrecPos[kPr]*nRrecPos[kPr]*nRrecPos[kPr];
            fNRMomNetPrGen[kBBBB] = nRgenNeg[kPr]*nRgenNeg[kPr]*nRgenNeg[kPr]*nRgenNeg[kPr];          fNRMomNetPrRec[kBBBB] = nRrecNeg[kPr]*nRrecNeg[kPr]*nRrecNeg[kPr]*nRrecNeg[kPr];
            //
            // fill tree which contains moments
            if(!fTreeSRedirector) return;
            if (nTracksrec>0){
              (*fTreeSRedirector)<<"mcFull"<<
              "ispileup="     << ipileup <<                  // pileup or not
              "syst="         << iset <<                     // systematic setting index
              "isample="      << sampleNo <<                 // sample id for subsample method
              "orig="         << iorig <<                    // origin type primary, or several combinations
              "nTPC="         << nTracksTPC <<               // number of tracks in the TPC
              "nITS="         << nTracksITS <<               // number of tracks in the ITS
              "vZ="           << fVz <<                      // event vertex z
              "cent="         << fCentrality <<              // centrality from V0
              "centimp="      << fCentImpBin <<              // centraltiy from impact parameter
              "impPar="       << fMCImpactParameter <<       // impact parameter taken from MC event header
              //
              "nhard="        << fNHardScatters <<           // Number of hard scatterings
              "nproj="        << fNProjectileParticipants << // Number of projectiles participants
              "ntarget="      << fNTargetParticipants <<     // Number of target participants
              "nn="           << fNNColl <<                  // Number of N-N collisions
              "nnw="          << fNNwColl <<                 // Number of N-Nwounded collisions
              "nwn="          << fNwNColl <<                 // Number of Nwounded-N collisons
              "nwnw="         << fNwNwColl <<                // Number of Nwounded-Nwounded collisions
              "nch="          << nChGen <<                   // Number of charged particles in 4pi
              "multv0a="      << multV0A <<                  // V0A multiplicity
              "multv0c="      << multV0C <<                  // V0C multiplicity
              //
              "pDown="        << fpDownArr[imom] <<          // lower edge of momentum bin
              "pUp="          << fpUpArr[imom] <<            // upper edge of momentum bin
              "etaDown="      << fetaDownArr[ieta] <<        // lower edge of eta bin
              "etaUp="        << fetaUpArr[ieta] <<          // upper edge of eta bin
              //
              "netPiMomGen.="   << &fMomNetPiGen <<         // momnets up to 3rd order for (net)pions on generated level with resonances
              "netKaMomGen.="   << &fMomNetKaGen <<         // momnets up to 3rd order for (net)kaons on generated level with resonances
              "netPrMomGen.="   << &fMomNetPrGen <<         // momnets up to 3rd order for (net)protons on generated level with resonances
              "netPiMomRec.="   << &fMomNetPiRec <<         // momnets up to 3rd order for (net)pions on reconstruced level with resonances
              "netKaMomRec.="   << &fMomNetKaRec <<         // momnets up to 3rd order for (net)kaons on reconstruced level with resonances
              "netPrMomRec.="   << &fMomNetPrRec <<         // momnets up to 3rd order for (net)protons on reconstruced level with resonances
              //
              "nRnetPiMomGen.="   << &fNRMomNetPiGen <<     // momnets up to 3rd order for (net)pions on generated level without resonances
              "nRnetKaMomGen.="   << &fNRMomNetKaGen <<     // momnets up to 3rd order for (net)kaons on generated level without resonances
              "nRnetPrMomGen.="   << &fNRMomNetPrGen <<     // momnets up to 3rd order for (net)protons on generated level without resonances
              "nRnetPiMomRec.="   << &fNRMomNetPiRec <<     // momnets up to 3rd order for (net)pions on reconstruced level without resonances
              "nRnetKaMomRec.="   << &fNRMomNetKaRec <<     // momnets up to 3rd order for (net)kaons on reconstruced level without resonances
              "nRnetPrMomRec.="   << &fNRMomNetPrRec <<     // momnets up to 3rd order for (net)protons on reconstruced level without resonances
              //
              // "posGen.="     << &genPos <<                  // counters for generated positive particles with resonances
              // "negGen.="     << &genNeg <<                  // counters for generated positive particles with resonances
              // "posRec.="     << &recPos <<                  // counters for reconstruced positive particles with resonances
              // "negRec.="     << &recNeg <<                  // counters for reconstruced positive particles with resonances
              // "nRposGen.="   << &nRgenPos <<                // counters for generated positive particles without resonances
              // "nRnegGen.="   << &nRgenNeg <<                // counters for generated positive particles without resonances
              // "nRposRec.="   << &nRrecPos <<                // counters for reconstruced positive particles without resonances
              // "nRnegRec.="   << &nRrecNeg <<                // counters for reconstruced positive particles without resonances
              "\n";
            }
            //
          }
        } // track origin loop
      } // ======= end of momentum loop =======
    } // ======= end of eta loop =======
  } // settitngs loop

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FastGen_NetParticles()
{
  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  // Assign subsample index
  Int_t sampleNo = 0;
  Int_t nSubSample = 20;
  sampleNo = Int_t(fEventGID)%nSubSample;
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FastGen_NetParticles ===== " << std::endl;
  //
  // ======================================================================
  // --------------   MC information with ideal PID   ---------------------
  // ======================================================================
  //
  const Int_t nParticles = 3;
  const Int_t nMoments   = 14;
  TVectorF genPos(nParticles);
  TVectorF genNeg(nParticles);
  TVectorF fMomNetPiGen(nMoments);
  TVectorF fMomNetKaGen(nMoments);
  TVectorF fMomNetPrGen(nMoments);
  //
  // Acceptance scan
  for (Int_t ieta=0; ieta<fNEtaWinBinsMC; ieta++){
    for (Int_t imom=0; imom<fNMomBinsMC; imom++){
      //
      // Initialize counters
      Int_t nTracksgen=0;
      for(Int_t i=0;i<nParticles; i++){ genPos[i]=0.; genNeg[i]=0.; }
      for(Int_t i=0;i<nMoments; i++)  { fMomNetPiGen[i]=0.; fMomNetKaGen[i]=0.; fMomNetPrGen[i]=0.; }
      //
      // ----------------------------   MC generated pure MC particles  --------------------------
      AliMCParticle *trackMCgen;
      for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++)
      {
        //
        // initialize the dummy particle id
        fElMCgen =-100.; fPiMCgen =-100.; fKaMCgen =-100.; fPrMCgen =-100.;
        trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
        if (!trackMCgen) continue;
        // 
        // Accepty only primary particles
        if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
        //
        // select particle of interest
        Int_t sign = trackMCgen->Charge();
        Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
        Int_t absPDG = TMath::Abs(pdg);
        Int_t iPart = -10;
        if (absPDG == kPDGpi) {iPart = 0; fPiMCgen = iPart;} // select pi+
        if (absPDG == kPDGka) {iPart = 1; fKaMCgen = iPart;} // select ka+
        if (absPDG == kPDGpr) {iPart = 2; fPrMCgen = iPart;} // select pr+
        if (iPart == -10) continue; // perfect PID cut
        //
        // Acceptance selection
        Double_t ptotMCgen = 0.;
        if(fUsePtCut==0) ptotMCgen = trackMCgen->P();
        if(fUsePtCut==1) ptotMCgen = trackMCgen->P();
        if(fUsePtCut==2) ptotMCgen = trackMCgen->Pt();
        Double_t etaMCgen = trackMCgen->Eta();
        Bool_t etaAcc  = (etaMCgen  >= fetaDownArr[ieta] && etaMCgen  <= fetaUpArr[ieta]);
        Bool_t momAcc  = (ptotMCgen >= fpDownArr[imom]   && ptotMCgen <= fpUpArr[imom]);
        Bool_t etaAccMaxWindow = (etaMCgen>=fEtaDown  && etaMCgen<=fEtaUp);
        Bool_t momAccMaxWindow = (ptotMCgen>=fMomDown && ptotMCgen<=fMomUp);
        Bool_t centAcc = (fCentrality>0);
        //
        // fill the moments
        if (etaAcc && momAcc && centAcc){
          nTracksgen++;
          if ( fPiMCgen>-1 && pdg<0) genNeg[kPi]++;
          if ( fKaMCgen>-1 && pdg<0) genNeg[kKa]++;
          if ( fPrMCgen>-1 && pdg<0) genNeg[kPr]++;
          // 
          if ( fPiMCgen>-1 && pdg>0) genPos[kPi]++;
          if ( fKaMCgen>-1 && pdg>0) genPos[kKa]++;
          if ( fPrMCgen>-1 && pdg>0) genPos[kPr]++;
        }

      } // ======= end of track loop for generated particles =======
      //
      // -----------------------------------------------------------------------------------------
      // --------------------   Calculation of moments on the event level  -----------------------
      // -----------------------------------------------------------------------------------------
      //
      // Net Pions
      fMomNetPiGen[kA]    = genPos[kPi];
      fMomNetPiGen[kB]    = genNeg[kPi];
      fMomNetPiGen[kAA]   = genPos[kPi]*genPos[kPi];
      fMomNetPiGen[kBB]   = genNeg[kPi]*genNeg[kPi];
      fMomNetPiGen[kAB]   = genPos[kPi]*genNeg[kPi];
      fMomNetPiGen[kAAA]  = genPos[kPi]*genPos[kPi]*genPos[kPi];
      fMomNetPiGen[kBBB]  = genNeg[kPi]*genNeg[kPi]*genNeg[kPi];
      fMomNetPiGen[kAAB]  = genPos[kPi]*genPos[kPi]*genNeg[kPi];
      fMomNetPiGen[kBBA]  = genNeg[kPi]*genNeg[kPi]*genPos[kPi];
      fMomNetPiGen[kABBB] = genPos[kPi]*genNeg[kPi]*genNeg[kPi]*genNeg[kPi];
      fMomNetPiGen[kAABB] = genPos[kPi]*genPos[kPi]*genNeg[kPi]*genNeg[kPi];
      fMomNetPiGen[kAAAB] = genPos[kPi]*genPos[kPi]*genPos[kPi]*genNeg[kPi];
      fMomNetPiGen[kAAAA] = genPos[kPi]*genPos[kPi]*genPos[kPi]*genPos[kPi];
      fMomNetPiGen[kBBBB] = genNeg[kPi]*genNeg[kPi]*genNeg[kPi]*genNeg[kPi];
      //
      // Net Kaons
      fMomNetKaGen[kA]    = genPos[kKa];
      fMomNetKaGen[kB]    = genNeg[kKa];
      fMomNetKaGen[kAA]   = genPos[kKa]*genPos[kKa];
      fMomNetKaGen[kBB]   = genNeg[kKa]*genNeg[kKa];
      fMomNetKaGen[kAB]   = genPos[kKa]*genNeg[kKa];
      fMomNetKaGen[kAAA]  = genPos[kKa]*genPos[kKa]*genPos[kKa];
      fMomNetKaGen[kBBB]  = genNeg[kKa]*genNeg[kKa]*genNeg[kKa];
      fMomNetKaGen[kAAB]  = genPos[kKa]*genPos[kKa]*genNeg[kKa];
      fMomNetKaGen[kBBA]  = genNeg[kKa]*genNeg[kKa]*genPos[kKa];
      fMomNetKaGen[kABBB] = genPos[kKa]*genNeg[kKa]*genNeg[kKa]*genNeg[kKa];
      fMomNetKaGen[kAABB] = genPos[kKa]*genPos[kKa]*genNeg[kKa]*genNeg[kKa];
      fMomNetKaGen[kAAAB] = genPos[kKa]*genPos[kKa]*genPos[kKa]*genNeg[kKa];
      fMomNetKaGen[kAAAA] = genPos[kKa]*genPos[kKa]*genPos[kKa]*genPos[kKa];
      fMomNetKaGen[kBBBB] = genNeg[kKa]*genNeg[kKa]*genNeg[kKa]*genNeg[kKa];
      //
      // Net Protons
      fMomNetPrGen[kA]    = genPos[kPr];
      fMomNetPrGen[kB]    = genNeg[kPr];
      fMomNetPrGen[kAA]   = genPos[kPr]*genPos[kPr];
      fMomNetPrGen[kBB]   = genNeg[kPr]*genNeg[kPr];
      fMomNetPrGen[kAB]   = genPos[kPr]*genNeg[kPr];
      fMomNetPrGen[kAAA]  = genPos[kPr]*genPos[kPr]*genPos[kPr];
      fMomNetPrGen[kBBB]  = genNeg[kPr]*genNeg[kPr]*genNeg[kPr];
      fMomNetPrGen[kAAB]  = genPos[kPr]*genPos[kPr]*genNeg[kPr];
      fMomNetPrGen[kBBA]  = genNeg[kPr]*genNeg[kPr]*genPos[kPr];
      fMomNetPrGen[kABBB] = genPos[kPr]*genNeg[kPr]*genNeg[kPr]*genNeg[kPr];
      fMomNetPrGen[kAABB] = genPos[kPr]*genPos[kPr]*genNeg[kPr]*genNeg[kPr];
      fMomNetPrGen[kAAAB] = genPos[kPr]*genPos[kPr]*genPos[kPr]*genNeg[kPr];
      fMomNetPrGen[kAAAA] = genPos[kPr]*genPos[kPr]*genPos[kPr]*genPos[kPr];
      fMomNetPrGen[kBBBB] = genNeg[kPr]*genNeg[kPr]*genNeg[kPr]*genNeg[kPr];
      //
      // fill tree which contains moments
      if(!fTreeSRedirector) return;
      if (nTracksgen>0){
        (*fTreeSRedirector)<<"mcGen"<<
        "ngenacc="      << nTracksgen <<
        "isample="      << sampleNo <<                 // sample id for subsample method
        "cent="         << fCentrality <<              // centrality from V0
        "centimp="      << fCentImpBin <<              // centraltiy from impact parameter
        "impPar="       << fMCImpactParameter <<       // impact parameter taken from MC event header
        //
        "nhard="        << fNHardScatters <<           // Number of hard scatterings
        "nproj="        << fNProjectileParticipants << // Number of projectiles participants
        "ntarget="      << fNTargetParticipants <<     // Number of target participants
        "nn="           << fNNColl <<                  // Number of N-N collisions
        "nnw="          << fNNwColl <<                 // Number of N-Nwounded collisions
        "nwn="          << fNwNColl <<                 // Number of Nwounded-N collisons
        "nwnw="         << fNwNwColl <<                // Number of Nwounded-Nwounded collisions
        //
        "pDown="        << fpDownArr[imom] <<          // lower edge of momentum bin
        "pUp="          << fpUpArr[imom] <<            // upper edge of momentum bin
        "etaDown="      << fetaDownArr[ieta] <<        // lower edge of eta bin
        "etaUp="        << fetaUpArr[ieta] <<          // upper edge of eta bin
        //
        "netPiMomGen.="   << &fMomNetPiGen <<         // momnets up to 3rd order for (net)pions on generated level with resonances
        "netKaMomGen.="   << &fMomNetKaGen <<         // momnets up to 3rd order for (net)kaons on generated level with resonances
        "netPrMomGen.="   << &fMomNetPrGen <<         // momnets up to 3rd order for (net)protons on generated level with resonances
        "\n";
      }


    } // ======= end of momentum loop =======
  } // ======= end of eta loop =======


}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FillTreeMC()
{

  Int_t trackOrigin = -10;
  Int_t sampleNo = 0;
  Int_t nSubSample = 20;
  sampleNo = Int_t(fEventGID)%nSubSample;
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillTreeMC ===== " << std::endl;
  //
  // ======================================================================
  // ------   reconstructed MC particles with dEdx information-------------
  // ======================================================================
  //
  Int_t tpcClusterMultiplicity   = fESD->GetNumberOfTPCClusters();
  const AliMultiplicity *multObj = fESD->GetMultiplicity();
  Int_t itsNumberOfTracklets   = multObj->GetNumberOfTracklets();
  for(Int_t irectrack = 0; irectrack < fESD->GetNumberOfTracks(); irectrack++)
  {
    //
    // Esd track
    //
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack *trackReal = fESD->GetTrack(irectrack);
    if (trackReal==NULL) continue;
    //
    // Get generated track info
    Int_t lab = TMath::Abs(trackReal->GetLabel());
    AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(lab);
    Int_t pdg = trackMCgen->Particle()->GetPdgCode();
    //
    Bool_t isTPCPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(lab,fMCEvent);
    Bool_t isITSPileup = AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCEvent, "Hijing");
    //
    // check the origin of the track
    trackOrigin = -10;
    if (fMCStack->IsPhysicalPrimary(lab))        trackOrigin = 0;
    if (fMCStack->IsSecondaryFromMaterial(lab))  trackOrigin = 1;
    if (fMCStack->IsSecondaryFromWeakDecay(lab)) trackOrigin = 2;
    if (trackOrigin<-1) continue;
    //
    // Track cuts from dtector
    Bool_t ifDefaultCuts = fESDtrackCuts->AcceptTrack(trackReal);
    Bool_t fBit96_spd    = fESDtrackCuts_Bit96_spd->AcceptTrack(trackReal);
    Bool_t fBit96_sdd    = fESDtrackCuts_Bit96_sdd->AcceptTrack(trackReal);
    Bool_t fBit96_base   = fESDtrackCuts_Bit96->AcceptTrack(trackReal);
    Bool_t fBit128       = fESDtrackCuts_Bit128->AcceptTrack(trackReal);
    Bool_t fBit768       = fESDtrackCuts_Bit768->AcceptTrack(trackReal);
    Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(trackReal);
    if (!trackReal->GetInnerParam()) continue;             // Ask if track is in the TPC
    if (!fESDtrackCutsLoose->AcceptTrack(trackReal))  continue;    // Loose Cuts
    if (!(trackReal->GetTPCsignalN()>0)) continue;
    //
    // match the track with mc track
    Int_t iPart = -10;
    if (TMath::Abs(pdg) == kPDGel) { iPart = 0; } // select el
    if (TMath::Abs(pdg) == kPDGpi) { iPart = 1; } // select pi
    if (TMath::Abs(pdg) == kPDGka) { iPart = 2; } // select ka
    if (TMath::Abs(pdg) == kPDGpr) { iPart = 3; } // select pr
    if (TMath::Abs(pdg) == kPDGde) { iPart = 4; } // select de
    if (iPart == -10) continue;
    //
    Float_t closestPar[3];
    GetExpecteds(trackReal,closestPar);
    SetCutBitsAndSomeTrackVariables(trackReal,iPart);
    //
    if (trackReal-> GetInnerParam()){
      fPtotMC       = trackReal->GetInnerParam()->GetP();
      fTPCSignalMC  = trackReal->GetTPCsignal();
    }
    fEtaMC        = trackReal->Eta();
    fPtMC         = trackReal->Pt();
    fSignMC       = trackReal->GetSign();
    Float_t pMC   = trackReal->P();
    Float_t fPhiMC= trackReal->Phi();
    fMissingCl    = trackReal->GetTPCClusterInfo(3,0,0,159);
    //
    Int_t nTPCClusters = fESD->GetNumberOfTPCClusters();
    Int_t nITSClusters = 0;
    AliVMultiplicity *multiObj = fESD->GetMultiplicity();
    for(Int_t i=2;i<6;i++) nITSClusters += multiObj->GetNumberOfITSClusters(i);
    //
    // different dca cuts
    // TMath::Abs(fTrackDCAxy)< 0.3
    Bool_t dca11h     = TMath::Abs(fTrackDCAxy)<0.0105+0.0350/TMath::Power(fPtMC,1.1);    // 10h tuned loose cut
    Bool_t dca10h     = TMath::Abs(fTrackDCAxy)<0.0182+0.0350/TMath::Power(fPtMC,1.01);    // 10h tuned loose cut
    Bool_t dcaBaseCut = TMath::Abs(fTrackDCAxy)<0.0208+0.04/TMath::Power(fPtMC,1.01);  // 10h tuned loose cut
    //
    // Fill MC closure tree
    if(!fTreeSRedirector) return;
    (*fTreeSRedirector)<<"fTreeMC"<<
    "defCut="    << ifDefaultCuts <<  // default cut
    "bit96spd="  << fBit96_spd <<  // run Number
    "bit96sdd="  << fBit96_sdd <<  // run Number
    "bit96="     << fBit96_base <<  // run Number
    "bit128="    << fBit128 <<  // run Number
    "bit768="    << fBit768 <<  // run Number
    "pixCut="    << ifDCAcutIfNoITSPixel <<  // run Number
    "run="       << fRunNo <<  // run Number
    "bField="    << fBField <<  // run Number
    "pileupbit=" << fPileUpBit <<
    "primmult="  << fNContributors <<  //  #prim tracks
    "tpcClMult=" << tpcClusterMultiplicity <<  //  TPC multiplicity
    //
    "orig="     << trackOrigin <<   // origin of the track
    "dcabase="  << dcaBaseCut <<  //  TPC multiplicity
    "dca10h="   << dca10h <<  //  TPC multiplicity
    "dca11h="   << dca11h <<  //  TPC multiplicity
    //
    "tpcpileup=" << isTPCPileup <<
    "itspileup=" << isITSPileup <<
    "tpcmult="   << fTPCMult <<  //  TPC multiplicity
    "itsmult="   << itsNumberOfTracklets <<
    "itsclmult="   << nITSClusters <<    // ITS multiplicity
    "tpcclmult="   << nTPCClusters <<    // ITS multiplicity
    //
    "isample="   << sampleNo <<                // sample id for subsample method
    "part="      << iPart <<
    "gid="       << fEventGID <<  //  global event ID
    "dEdx="      << fTPCSignalMC <<    // dEdx of mc track
    "cutBit="    << fTrackCutBits <<  //  Systematic Cuts
    "sign="      << fSignMC <<         // sign
    "ptot="      << fPtotMC <<         // tpc momentum
    "p="         << pMC <<             // vertex momentum
    "pT="        << fPtMC <<           // transverse momentum
    "eta="       << fEtaMC <<          // mc eta
    "phi="       << fPhiMC <<          // mc eta
    "cent="      << fCentrality <<     // Centrality
    "centimp="   << fCentImpBin <<
    "vZ="        << fVz <<
    "nsigmatofka="  << fNSigmasKaTOF         <<  // interaction rate
    "nsigmatofpr="  << fNSigmasPrTOF         <<  // interaction rate
    "dcaxy="     << fTrackDCAxy           <<
    "dcaz="      << fTrackDCAz            <<
    "ncltpc="    << fNcl                  <<  //  centrality
    "ncltpccorr="<< fNclCorr              <<  //  centrality
    "cRows="     << fTrackTPCCrossedRows  <<
    "missCl="    << fMissingCl <<
    "chi2tpc="   << fTrackChi2TPC         <<
    "chi2tpccorr=" << fTrackChi2TPCcorr         <<
    "\n";

  } // ======= end of track loop for MC dEdx filling =======



}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FillGenDistributions()
{

  for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++)
  { // track loop
    //
    Bool_t isTPCPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iTrack,fMCEvent);
    Bool_t isITSPileup = AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCEvent, "Hijing");
    if (isTPCPileup || isITSPileup) continue;
    //
    // initialize the dummy particle id
    fPiMCgen =-100.; fKaMCgen =-100.; fPrMCgen =-100.;
    AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
    //
    // check the origin of the track
    Int_t trackOrigin = -10;
    if (fMCStack->IsPhysicalPrimary(iTrack))        trackOrigin = 0;
    if (fMCStack->IsSecondaryFromMaterial(iTrack))  trackOrigin = 1;
    if (fMCStack->IsSecondaryFromWeakDecay(iTrack)) trackOrigin = 2;
    if (trackOrigin<-1) continue;
    //
    Float_t ptotMCgen = trackMCgen->P();
    Float_t pTMCgen   = trackMCgen->Pt();
    Float_t phiMCGen  = trackMCgen->Phi();
    Float_t etaMCgen  = trackMCgen->Eta();
    Float_t rapMCgen  = trackMCgen->Y();
    //
    // select particle of interest
    Int_t iPart = -10;
    Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
    Int_t sign = (pdg>0) ? 1 : -1;
    if (TMath::Abs(pdg) == kPDGpi) {iPart = 1; fPiMCgen = iPart;} // select pi+
    if (TMath::Abs(pdg) == kPDGka) {iPart = 2; fKaMCgen = iPart;} // select ka+
    if (TMath::Abs(pdg) == kPDGpr) {iPart = 3; fPrMCgen = iPart;} // select pr+
    if (iPart == -10) continue;
    //
    // Resonance control
    Bool_t parInterest = (fPiMCgen>-1||fKaMCgen>-1||fPrMCgen>-1) ? kTRUE : kFALSE;
    Bool_t acceptRes = CheckIfFromResonance(1,trackMCgen,iTrack,parInterest,ptotMCgen,etaMCgen,fCentrality,kTRUE);
    //
    // Fill MC closure tree
    if(!fTreeSRedirector) return;
    (*fTreeSRedirector)<<"mcGen"<<
    "gid="       << fEventGID <<
    "acceptRes=" << acceptRes <<                // sample id for subsample method
    "part="      << iPart <<                // sample id for subsample method
    "origin="    << trackOrigin <<
    "sign="      << sign <<         // sign
    "p="         << ptotMCgen <<             // vertex momentum
    "pT="        << pTMCgen <<           // transverse momentum
    "eta="       << etaMCgen <<          // mc eta
    "rap="       << rapMCgen <<          // mc eta
    "phi="       << phiMCGen <<          // mc eta
    "cent="      << fCentrality <<     // Centrality
    "vZ="        << fVz <<
    "\n";

  } // ======= end of track loop for generated particles to see distributions =======

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FastGen()
{

  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FastGen ===== " << std::endl;
  Int_t sampleNo = 0;
  Int_t nSubSample = 20;
  sampleNo = Int_t(fEventGID)%nSubSample;  // ======================================================================
  // ======================================================================
  //
  // ========= Efficiency Check Eta momentum and Centrality scan ==========
  // vectors to hold moments
  const Int_t nMoments = 15;
  TVectorF genMoments(nMoments);
  TVectorF genMomentsPos(nMoments);
  TVectorF genMomentsNeg(nMoments);
  TVectorF genMomentsCross(nMoments);
  // Moments without resonances
  TVectorF nRgenMoments(nMoments);
  TVectorF nRgenMomentsPos(nMoments);
  TVectorF nRgenMomentsNeg(nMoments);
  TVectorF nRgenMomentsCross(nMoments);
  //
  // Protection for full acceptance case
  Int_t nOriginType = (fTrackOriginOnlyPrimary!=1) ? 4 : 1;
  //
  // Acceptance scan
  for (Int_t ieta=0; ieta<fNEtaWinBinsMC; ieta++){
    for (Int_t imom=0; imom<fNMomBinsMC; imom++){
      for (Int_t iorig=0; iorig<nOriginType; iorig++){
        //
        // initialize counters
        for(Int_t i=0;i<nMoments; i++){
          genMoments[i]=0.;       nRgenMoments[i]=0.;
          genMomentsPos[i]=0.;    nRgenMomentsPos[i]=0.;
          genMomentsNeg[i]=0.;    nRgenMomentsNeg[i]=0.;
          genMomentsCross[i]=0.;  nRgenMomentsCross[i]=0.;
        }
        //
        // enter track loop for counting
        Float_t trCountgen=0;
        AliMCParticle *trackMCgen;
        for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++) {
          // track loop
          //
          // Select real trigger event and reject other pile up vertices
          Bool_t isTPCPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iTrack,fMCEvent);
          Bool_t isITSPileup = AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCEvent, "Hijing");
          if (isTPCPileup || isITSPileup) continue;
          //
          // initialize the dummy particle id
          fElMCgen =-100.; fPiMCgen =-100.; fKaMCgen =-100.; fPrMCgen =-100.; fDeMCgen =-100.; fMuMCgen =-100.; fLaMCgen =-100., fBaMCgen =-100.;
          trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
          if (!trackMCgen) continue;
          //
          // check the origin of the track
          Bool_t bPrim     = fMCStack->IsPhysicalPrimary(iTrack);
          Bool_t bMaterial = fMCStack->IsSecondaryFromMaterial(iTrack);
          Bool_t bWeak     = fMCStack->IsSecondaryFromWeakDecay(iTrack);
          Bool_t bAcceptOrigin = kFALSE;
          if (iorig==0) bAcceptOrigin = bPrim;
          if (iorig==1) bAcceptOrigin = (bPrim || bWeak);
          if (iorig==2) bAcceptOrigin = (bPrim || bMaterial);
          if (iorig==3) bAcceptOrigin = (bPrim || bMaterial || bWeak);
          if (!bAcceptOrigin) continue; // origin cut
          //
          // select sigle particle type
          Int_t sign = trackMCgen->Charge();
          Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
          Int_t absPDG = TMath::Abs(pdg);
          Int_t iPart = -10;
          if (absPDG == kPDGel) {iPart = 0; fElMCgen = iPart;} // select el-
          if (absPDG == kPDGpi) {iPart = 1; fPiMCgen = iPart;} // select pi+
          if (absPDG == kPDGka) {iPart = 2; fKaMCgen = iPart;} // select ka+
          if (absPDG == kPDGpr) {iPart = 3; fPrMCgen = iPart;} // select pr+
          if (absPDG == kPDGde) {iPart = 4; fDeMCgen = iPart;} // select de
          if (absPDG == kPDGmu) {iPart = 5; fMuMCgen = iPart;} // select mu-
          if (absPDG == kPDGla) {iPart = 6; fLaMCgen = iPart;} // select Lambda
          //
          // tag baryons
          for (Int_t ibar=0;ibar<fNBarBins;ibar++){
            if ( fBaryons[ibar] == absPDG ){
              iPart = 7; fBaMCgen = iPart; break;
            }
          }
          //
          Double_t ptotMCgen = 0.;
          if(fUsePtCut==0) ptotMCgen = trackMCgen->P();
          if(fUsePtCut==1) ptotMCgen = trackMCgen->P();
          if(fUsePtCut==2) ptotMCgen = trackMCgen->Pt();
          Double_t etaMCgen = (fRapidityType==0) ? trackMCgen->Eta() :  trackMCgen->Y();
          Bool_t etaAcc  = (etaMCgen>=fetaDownArr[ieta] && etaMCgen<=fetaUpArr[ieta]);
          Bool_t momAcc  = (ptotMCgen>=fpDownArr[imom]  && ptotMCgen<fpUpArr[imom]);
          Bool_t centAcc = (fCentrality>0);
          //
          if (iorig==0 && fPrMCgen>0 && fetaUpArr[ieta]>11.7 && fpUpArr[imom]>100.)
          {
            if (pdg>0) fHistRapDistFullAccPr->Fill(etaMCgen, fCentrality, 1.);
            if (pdg<0) fHistRapDistFullAccAPr->Fill(etaMCgen, fCentrality, 1.);
          }
          //
          // Check particle is from a Resonance
          // Bool_t acceptRes = CheckIfFromAnyResonance(trackMCgen,fetaDownArr[ieta],fetaUpArr[ieta],fpDownArr[imom],fpUpArr[imom]);
          //
          // Resonance control
          Bool_t parInterest = (fPiMCgen>-1||fKaMCgen>-1||fPrMCgen>-1||fElMCgen>-1||fLaMCgen>-1) ? kTRUE : kFALSE;
          // Bool_t acceptRes = CheckIfFromResonance(1,trackMCgen,iTrack,parInterest,ptotMCgen,etaMCgen,fCentrality,kTRUE);
          Bool_t acceptRes = kTRUE;  // TODO
          //
          // count first moments
          if (etaAcc && momAcc && centAcc){
            //
            // count charged particles
            trCountgen++;
            if (sign>0 || sign<0) genMoments[kCh]++;
            if (sign>0) genMomentsPos[kCh]++;
            if (sign<0) genMomentsNeg[kCh]++;
            if ( acceptRes ) {
              if (sign>0 || sign<0) nRgenMoments[kCh]++;
              if (sign>0) nRgenMomentsPos[kCh]++;
              if (sign<0) nRgenMomentsNeg[kCh]++;
            }
            //
            if ( iPart == -10) continue; // perfect PID cut
            //
            if ( fPiMCgen>-1 ) genMoments[kPi]++;
            if ( fKaMCgen>-1 ) genMoments[kKa]++;
            if ( fPrMCgen>-1 ) genMoments[kPr]++;
            //
            if ( fPiMCgen>-1 && pdg<0) genMomentsNeg[kPi]++;
            if ( fKaMCgen>-1 && pdg<0) genMomentsNeg[kKa]++;
            if ( fPrMCgen>-1 && pdg<0) genMomentsNeg[kPr]++;
            //
            if ( fPiMCgen>-1 && pdg>0) genMomentsPos[kPi]++;
            if ( fKaMCgen>-1 && pdg>0) genMomentsPos[kKa]++;
            if ( fPrMCgen>-1 && pdg>0) genMomentsPos[kPr]++;
            // Lambdas for alice
            if ( fLaMCgen>-1 ) genMoments[kLa]++;
            if ( fLaMCgen>-1 && pdg>0) genMomentsPos[kLa]++;
            if ( fLaMCgen>-1 && pdg<0) genMomentsNeg[kLa]++;
            // Net baryons
            if ( fBaMCgen>-1 ) genMoments[kBa]++;
            if ( fBaMCgen>-1 && pdg>0) genMomentsPos[kBa]++;
            if ( fBaMCgen>-1 && pdg<0) genMomentsNeg[kBa]++;
            // reject resonances
            if ( acceptRes ) {
              if ( fPiMCgen>-1 ) nRgenMoments[kPi]++;
              if ( fKaMCgen>-1 ) nRgenMoments[kKa]++;
              if ( fPrMCgen>-1 ) nRgenMoments[kPr]++;

              if ( fPiMCgen>-1 && pdg<0) nRgenMomentsNeg[kPi]++;
              if ( fKaMCgen>-1 && pdg<0) nRgenMomentsNeg[kKa]++;
              if ( fPrMCgen>-1 && pdg<0) nRgenMomentsNeg[kPr]++;

              if ( fPiMCgen>-1 && pdg>0) nRgenMomentsPos[kPi]++;
              if ( fKaMCgen>-1 && pdg>0) nRgenMomentsPos[kKa]++;
              if ( fPrMCgen>-1 && pdg>0) nRgenMomentsPos[kPr]++;
              // Lambdas for alice
              if ( fLaMCgen>-1 ) nRgenMoments[kLa]++;
              if ( fLaMCgen>-1 && pdg>0) nRgenMomentsPos[kLa]++;
              if ( fLaMCgen>-1 && pdg<0) nRgenMomentsNeg[kLa]++;
              // Net baryon
              if ( fBaMCgen>-1 ) nRgenMoments[kBa]++;
              if ( fBaMCgen>-1 && pdg>0) nRgenMomentsPos[kBa]++;
              if ( fBaMCgen>-1 && pdg<0) nRgenMomentsNeg[kBa]++;
            }
          }
        } // ======= end of track loop =======

        // calculate second moments                                              // calculate second moments with resonances
        genMoments[kPiPi]=genMoments[kPi]*genMoments[kPi];                       nRgenMoments[kPiPi]=nRgenMoments[kPi]*nRgenMoments[kPi];
        genMoments[kKaKa]=genMoments[kKa]*genMoments[kKa];                       nRgenMoments[kKaKa]=nRgenMoments[kKa]*nRgenMoments[kKa];
        genMoments[kPrPr]=genMoments[kPr]*genMoments[kPr];                       nRgenMoments[kPrPr]=nRgenMoments[kPr]*nRgenMoments[kPr];
        genMoments[kPiKa]=genMoments[kPi]*genMoments[kKa];                       nRgenMoments[kPiKa]=nRgenMoments[kPi]*nRgenMoments[kKa];
        genMoments[kPiPr]=genMoments[kPi]*genMoments[kPr];                       nRgenMoments[kPiPr]=nRgenMoments[kPi]*nRgenMoments[kPr];
        genMoments[kKaPr]=genMoments[kKa]*genMoments[kPr];                       nRgenMoments[kKaPr]=nRgenMoments[kKa]*nRgenMoments[kPr];
        genMomentsNeg[kPiPi]=genMomentsNeg[kPi]*genMomentsNeg[kPi];              nRgenMomentsNeg[kPiPi]=nRgenMomentsNeg[kPi]*nRgenMomentsNeg[kPi];
        genMomentsNeg[kKaKa]=genMomentsNeg[kKa]*genMomentsNeg[kKa];              nRgenMomentsNeg[kKaKa]=nRgenMomentsNeg[kKa]*nRgenMomentsNeg[kKa];
        genMomentsNeg[kPrPr]=genMomentsNeg[kPr]*genMomentsNeg[kPr];              nRgenMomentsNeg[kPrPr]=nRgenMomentsNeg[kPr]*nRgenMomentsNeg[kPr];
        genMomentsNeg[kPiKa]=genMomentsNeg[kPi]*genMomentsNeg[kKa];              nRgenMomentsNeg[kPiKa]=nRgenMomentsNeg[kPi]*nRgenMomentsNeg[kKa];
        genMomentsNeg[kPiPr]=genMomentsNeg[kPi]*genMomentsNeg[kPr];              nRgenMomentsNeg[kPiPr]=nRgenMomentsNeg[kPi]*nRgenMomentsNeg[kPr];
        genMomentsNeg[kKaPr]=genMomentsNeg[kKa]*genMomentsNeg[kPr];              nRgenMomentsNeg[kKaPr]=nRgenMomentsNeg[kKa]*nRgenMomentsNeg[kPr];
        genMomentsPos[kPiPi]=genMomentsPos[kPi]*genMomentsPos[kPi];              nRgenMomentsPos[kPiPi]=nRgenMomentsPos[kPi]*nRgenMomentsPos[kPi];
        genMomentsPos[kKaKa]=genMomentsPos[kKa]*genMomentsPos[kKa];              nRgenMomentsPos[kKaKa]=nRgenMomentsPos[kKa]*nRgenMomentsPos[kKa];
        genMomentsPos[kPrPr]=genMomentsPos[kPr]*genMomentsPos[kPr];              nRgenMomentsPos[kPrPr]=nRgenMomentsPos[kPr]*nRgenMomentsPos[kPr];
        genMomentsPos[kPiKa]=genMomentsPos[kPi]*genMomentsPos[kKa];              nRgenMomentsPos[kPiKa]=nRgenMomentsPos[kPi]*nRgenMomentsPos[kKa];
        genMomentsPos[kPiPr]=genMomentsPos[kPi]*genMomentsPos[kPr];              nRgenMomentsPos[kPiPr]=nRgenMomentsPos[kPi]*nRgenMomentsPos[kPr];
        genMomentsPos[kKaPr]=genMomentsPos[kKa]*genMomentsPos[kPr];              nRgenMomentsPos[kKaPr]=nRgenMomentsPos[kKa]*nRgenMomentsPos[kPr];
        genMomentsCross[kPiPosPiNeg]=genMomentsPos[kPi]*genMomentsNeg[kPi];      nRgenMomentsCross[kPiPosPiNeg]=nRgenMomentsPos[kPi]*nRgenMomentsNeg[kPi];
        genMomentsCross[kPiPosKaNeg]=genMomentsPos[kPi]*genMomentsNeg[kKa];      nRgenMomentsCross[kPiPosKaNeg]=nRgenMomentsPos[kPi]*nRgenMomentsNeg[kKa];
        genMomentsCross[kPiPosPrNeg]=genMomentsPos[kPi]*genMomentsNeg[kPr];      nRgenMomentsCross[kPiPosPrNeg]=nRgenMomentsPos[kPi]*nRgenMomentsNeg[kPr];
        genMomentsCross[kKaPosPiNeg]=genMomentsPos[kKa]*genMomentsNeg[kPi];      nRgenMomentsCross[kKaPosPiNeg]=nRgenMomentsPos[kKa]*nRgenMomentsNeg[kPi];
        genMomentsCross[kKaPosKaNeg]=genMomentsPos[kKa]*genMomentsNeg[kKa];      nRgenMomentsCross[kKaPosKaNeg]=nRgenMomentsPos[kKa]*nRgenMomentsNeg[kKa];
        genMomentsCross[kKaPosPrNeg]=genMomentsPos[kKa]*genMomentsNeg[kPr];      nRgenMomentsCross[kKaPosPrNeg]=nRgenMomentsPos[kKa]*nRgenMomentsNeg[kPr];
        genMomentsCross[kPrPosPiNeg]=genMomentsPos[kPr]*genMomentsNeg[kPi];      nRgenMomentsCross[kPrPosPiNeg]=nRgenMomentsPos[kPr]*nRgenMomentsNeg[kPi];
        genMomentsCross[kPrPosKaNeg]=genMomentsPos[kPr]*genMomentsNeg[kKa];      nRgenMomentsCross[kPrPosKaNeg]=nRgenMomentsPos[kPr]*nRgenMomentsNeg[kKa];
        genMomentsCross[kPrPosPrNeg]=genMomentsPos[kPr]*genMomentsNeg[kPr];      nRgenMomentsCross[kPrPosPrNeg]=nRgenMomentsPos[kPr]*nRgenMomentsNeg[kPr];
        //
        // net lambda for Alice                                                  // net lambda for Alice
        genMoments[kLaLa]=genMoments[kLa]*genMoments[kLa];                       nRgenMoments[kLaLa]=nRgenMoments[kLa]*nRgenMoments[kLa];
        genMomentsNeg[kLaLa]=genMomentsNeg[kLa]*genMomentsNeg[kLa];              nRgenMomentsNeg[kLaLa]=nRgenMomentsNeg[kLa]*nRgenMomentsNeg[kLa];
        genMomentsPos[kLaLa]=genMomentsPos[kLa]*genMomentsPos[kLa];              nRgenMomentsPos[kLaLa]=nRgenMomentsPos[kLa]*nRgenMomentsPos[kLa];
        genMomentsCross[kLaPosLaNeg]=genMomentsPos[kLa]*genMomentsNeg[kLa];      nRgenMomentsCross[kLaPosLaNeg]=nRgenMomentsPos[kLa]*nRgenMomentsNeg[kLa];

        genMoments[kChCh]=genMoments[kCh]*genMoments[kCh];                       nRgenMoments[kChCh]=nRgenMoments[kCh]*nRgenMoments[kCh];
        genMomentsNeg[kChCh]=genMomentsNeg[kCh]*genMomentsNeg[kCh];              nRgenMomentsNeg[kChCh]=nRgenMomentsNeg[kCh]*nRgenMomentsNeg[kCh];
        genMomentsPos[kChCh]=genMomentsPos[kCh]*genMomentsPos[kCh];              nRgenMomentsPos[kChCh]=nRgenMomentsPos[kCh]*nRgenMomentsPos[kCh];
        genMomentsCross[kChPosChNeg]=genMomentsPos[kCh]*genMomentsNeg[kCh];      nRgenMomentsCross[kChPosChNeg]=nRgenMomentsPos[kCh]*nRgenMomentsNeg[kCh];

        genMoments[kBaBa]=genMoments[kBa]*genMoments[kBa];                       nRgenMoments[kBaBa]=nRgenMoments[kBa]*nRgenMoments[kBa];
        genMomentsNeg[kBaBa]=genMomentsNeg[kBa]*genMomentsNeg[kBa];              nRgenMomentsNeg[kBaBa]=nRgenMomentsNeg[kBa]*nRgenMomentsNeg[kBa];
        genMomentsPos[kBaBa]=genMomentsPos[kBa]*genMomentsPos[kBa];              nRgenMomentsPos[kBaBa]=nRgenMomentsPos[kBa]*nRgenMomentsPos[kBa];
        genMomentsCross[kBaPosBaNeg]=genMomentsPos[kBa]*genMomentsNeg[kBa];      nRgenMomentsCross[kBaPosBaNeg]=nRgenMomentsPos[kBa]*nRgenMomentsNeg[kBa];
        //
        // fill tree which contains moments
        if(!fTreeSRedirector) return;
        // if there is at least one track in an event fill the tree
        if (trCountgen>0){
          (*fTreeSRedirector)<<"mcGen"<<
          "isample="      << sampleNo <<                // sample id for subsample method
          "orig="         << iorig <<
          "vZ="           << fVz <<
          "nhard="        << fNHardScatters <<           // Number of hard scatterings
          "nproj="        << fNProjectileParticipants << // Number of projectiles participants
          "ntarget="      << fNTargetParticipants <<     // Number of target participants
          "nn="           << fNNColl <<                  // Number of N-N collisions
          "nnw="          << fNNwColl <<                 // Number of N-Nwounded collisions
          "nwn="          << fNwNColl <<                 // Number of Nwounded-N collisons
          "nwnw="         << fNwNwColl <<                // Number of Nwounded-Nwounded collisions
          "cent="         << fCentrality <<
          "centimp="      << fCentImpBin <<
          "impPar="       << fMCImpactParameter <<      // impact parameter taken from MC event header
          "pDown="        << fpDownArr[imom] <<         // lower edge of momentum bin
          "pUp="          << fpUpArr[imom] <<           // upper edge of momentum bin
          "etaDown="      << fetaDownArr[ieta] <<       // lower edge of eta bin
          "etaUp="        << fetaUpArr[ieta] <<         // upper edge of eta bin
          "momentGen.="      << &genMoments <<             // second moments for particle+antiparticle
          "momentPosGen.="   << &genMomentsPos <<          // second moment of positive particles
          "momentNegGen.="   << &genMomentsNeg <<          // second moment of negative particles
          "momentCrossGen.=" << &genMomentsCross <<        // second moment of unlikesign particles
          // "noResmomentGen.="      << &nRgenMoments <<             // second moments for particle+antiparticle
          // "noResmomentPosGen.="   << &nRgenMomentsPos <<          // second moment of positive particles
          // "noResmomentNegGen.="   << &nRgenMomentsNeg <<          // second moment of negative particles
          // "noResmomentCrossGen.=" << &nRgenMomentsCross <<        // second moment of unlikesign particles
          "\n";
        } // tree filling
      } // track origin loop
    }// ======= end of momentum loop =======
  } // ======= end of eta loop =======
  //
  // ======================================================================
  //
  // ======================================================================


}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FastGenHigherMoments()
{

  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FastGenHigherMoments ===== " << std::endl;

  Int_t sampleNo = 0;
  // Int_t evtNuminFile = fMCEvent -> GetEventNumberInFile();

  // ======================================================================
  // ======================================================================
  //
  // ========= Efficiency Check Eta momentum and Centrality scan ==========
  const Int_t nHighMoments = 5;
  const Int_t nMoments = 13;
  for (Int_t ieta=0; ieta<fNEtaWinBinsMC; ieta++){
    for (Int_t imom=0; imom<fNMomBinsMC; imom++){
      for (Int_t icent=0; icent<fNCentbinsData-1; icent++){

        // vectors to hold moments
        TVectorF netPi(nHighMoments);
        TVectorF netKa(nHighMoments);
        TVectorF netPr(nHighMoments);
        TVectorF netLa(nHighMoments);
        TVectorF netCh(nHighMoments);
        TVectorF noResnetPi(nHighMoments);
        TVectorF noResnetKa(nHighMoments);
        TVectorF noResnetPr(nHighMoments);
        TVectorF noResnetLa(nHighMoments);
        TVectorF noResnetCh(nHighMoments);

        // initialize counters
        for(Int_t i=0;i<nHighMoments; i++){
          netPi[i]=0.;  noResnetPi[i]=0.;
          netKa[i]=0.;  noResnetKa[i]=0.;
          netPr[i]=0.;  noResnetPr[i]=0.;
          netLa[i]=0.;  noResnetLa[i]=0.;
          netCh[i]=0.;  noResnetCh[i]=0.;
        }

        // vectors to hold moments
        TVectorF genMomentsPos(nMoments);
        TVectorF genMomentsNeg(nMoments);
        TVectorF genMomentsCross(nMoments);
        TVectorF nRgenMomentsPos(nMoments);
        TVectorF nRgenMomentsNeg(nMoments);
        TVectorF nRgenMomentsCross(nMoments);
        // initialize counters
        for(Int_t i=0;i<nMoments; i++){
          genMomentsPos[i]=0.;    nRgenMomentsPos[i]=0.;
          genMomentsNeg[i]=0.;    nRgenMomentsNeg[i]=0.;
          genMomentsCross[i]=0.;  nRgenMomentsCross[i]=0.;
        }

        Double_t centBin = (fcentDownArr[icent]+fcentUpArr[icent])/2.;
        Float_t nTracksgen=0, trCountgen=0;
        AliMCParticle *trackMCgen;
        Int_t nStackTracks = fMCEvent->GetNumberOfTracks();
        // TRACK LOOP
        for (Int_t iTrack = 0; iTrack < nStackTracks; iTrack++) {
          // track loop
          //
          // Select real trigger event and reject other pile up vertices
          Bool_t isTPCPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iTrack,fMCEvent);
          Bool_t isITSPileup = AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCEvent, "Hijing");
          if (isTPCPileup || isITSPileup) continue;
          //
          // initialize the dummy particle id
          fElMCgen =-100.; fPiMCgen =-100.; fKaMCgen =-100.; fPrMCgen =-100.; fDeMCgen =-100.; fMuMCgen =-100.; fLaMCgen =-100.;
          trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
          //
          // apply primary vertex and eta cut
          if ((trackMCgen->Eta()<fetaDownArr[ieta]) || (trackMCgen->Eta()>fetaUpArr[ieta])) continue;
          //
          // iwith or wihout weak decays
          if (fWeakAndMaterial){
            if ( !(fMCStack->IsPhysicalPrimary(iTrack) || fMCStack->IsSecondaryFromWeakDecay(iTrack)) ) continue;
          } else if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
          //
          // get the pdg info for maother and daughter
          Int_t sign = trackMCgen->Particle()->GetPDG()->Charge();
          Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
          Int_t labMom = trackMCgen->GetMother();
          TObjString parName(trackMCgen->Particle()->GetName());
          Int_t pdgMom = 0;
          TObjString momName="xxx";
          if ((labMom>=0) && (labMom < nStackTracks)){
            pdgMom = fMCStack->Particle(labMom)->GetPdgCode();
            momName = fMCStack->Particle(labMom)->GetName();
          }
          //
          // Check if the particle is in the black list of resonances
          Bool_t acceptRes = kTRUE;
          for (Int_t ires=0;ires<fNResBins;ires++){

            if (fResonances[ires].find("xxx")){
              // reject all resonances
              if (!(momName.GetString().Contains(fResonances[ires]))) {acceptRes=kFALSE; break;}
            } else {
              // reject resonances in the array
              if (momName.GetString().Contains(fResonances[ires])) {acceptRes=kFALSE; break;}
            }
          }

          Int_t iPart = -10;
          if (TMath::Abs(pdg) == kPDGel) {iPart = 0; fElMCgen = iPart;} // select el-
          if (TMath::Abs(pdg) == kPDGpi) {iPart = 1; fPiMCgen = iPart;} // select pi+
          if (TMath::Abs(pdg) == kPDGka) {iPart = 2; fKaMCgen = iPart;} // select ka+
          if (TMath::Abs(pdg) == kPDGpr) {iPart = 3; fPrMCgen = iPart;} // select pr+
          if (TMath::Abs(pdg) == kPDGde) {iPart = 4; fDeMCgen = iPart;} // select de
          if (TMath::Abs(pdg) == kPDGmu) {iPart = 5; fMuMCgen = iPart;} // select mu-
          if (TMath::Abs(pdg) == kPDGla) {iPart = 6; fLaMCgen = iPart;} // select Lambda

          // dump resonance info
          if(fEventCountInFile==2 && !fRunOnGrid) {
            Bool_t parInterest = (fPiMCgen>-1||fKaMCgen>-1||fPrMCgen>-1||fElMCgen>-1||fLaMCgen>-1) ? kTRUE : kFALSE;
            if(!fTreeSRedirector) return;
            (*fTreeSRedirector)<<"resonance"<<
            "acceptRes="   << acceptRes <<
            "parInterest=" << parInterest <<          // only pi, ka, and proton
            "centBin="     << centBin <<                 // cent bin
            "pDown="       << fpDownArr[imom] <<         // lower edge of momentum bin
            "etaDown="     << fetaDownArr[ieta] <<       // lower edge of eta bin
            "pdg="         << pdg      <<         // pdg of prim particle
            "lab="         << iTrack   <<         // index of prim particle
            "pdgMom="      << pdgMom   <<         // pdg of mother
            "labMom="      << labMom   <<         // index of mother
            "parName.="    << &parName <<         //  full path - file name with ESD
            "momName.="    << &momName <<         //  full path - file name with ESD
            "\n";
          }

          // count first moments
          Float_t ptotMCgen = trackMCgen->P();
          if ((fCentrality>=fcentDownArr[icent])
          &&(fCentrality<fcentUpArr[icent])
          &&(ptotMCgen>=fpDownArr[imom])
          &&(ptotMCgen<=fpUpArr[imom]))
          {
            nTracksgen++;
            //
            // count charged particles
            if (sign>0) genMomentsPos[kCh]++;
            if (sign<0) genMomentsNeg[kCh]++;
            if ( acceptRes ) {
              if (sign>0) nRgenMomentsPos[kCh]++;
              if (sign<0) nRgenMomentsNeg[kCh]++;
            }
            // Count identified particles
            if (iPart == -10) continue;
            if ( fPiMCgen>-1 || fKaMCgen>-1 || fPrMCgen>-1 || fLaMCgen>-1) trCountgen++;
            //
            if ( fPiMCgen>-1 && pdg<0) genMomentsNeg[kPi]++;
            if ( fKaMCgen>-1 && pdg<0) genMomentsNeg[kKa]++;
            if ( fPrMCgen>-1 && pdg<0) genMomentsNeg[kPr]++;
            //
            if ( fPiMCgen>-1 && pdg>0) genMomentsPos[kPi]++;
            if ( fKaMCgen>-1 && pdg>0) genMomentsPos[kKa]++;
            if ( fPrMCgen>-1 && pdg>0) genMomentsPos[kPr]++;
            // Lambdas for alice
            if ( fLaMCgen>-1 && pdg>0) genMomentsPos[kLa]++;
            if ( fLaMCgen>-1 && pdg<0) genMomentsNeg[kLa]++;
            // reject resonances
            if ( acceptRes ) {
              if ( fPiMCgen>-1 && pdg<0) nRgenMomentsNeg[kPi]++;
              if ( fKaMCgen>-1 && pdg<0) nRgenMomentsNeg[kKa]++;
              if ( fPrMCgen>-1 && pdg<0) nRgenMomentsNeg[kPr]++;
              //
              if ( fPiMCgen>-1 && pdg>0) nRgenMomentsPos[kPi]++;
              if ( fKaMCgen>-1 && pdg>0) nRgenMomentsPos[kKa]++;
              if ( fPrMCgen>-1 && pdg>0) nRgenMomentsPos[kPr]++;
              // Lambdas for alice
              if ( fLaMCgen>-1 && pdg>0) nRgenMomentsPos[kLa]++;
              if ( fLaMCgen>-1 && pdg<0) nRgenMomentsNeg[kLa]++;
            }
          }
        } // ======= end of track loop =======

        //       moments from Lookup table
        //         std::cout << " Info::marsland: ====================== " << imom << "  " << icent << "  " << ieta << " ====================== " << std::endl;
        //         std::cout << fNetPiFirstMoments[0][imom][icent][ieta] << std::endl;
        //         std::cout << fNetKaFirstMoments[0][imom][icent][ieta] << std::endl;
        //         std::cout << fNetPrFirstMoments[0][imom][icent][ieta] << std::endl;
        //         std::cout << fNetLaFirstMoments[0][imom][icent][ieta] << std::endl;
        //         std::cout << fNetChFirstMoments[0][imom][icent][ieta] << std::endl;
        //
        //         std::cout << fNetPiFirstMoments[1][imom][icent][ieta] << std::endl;
        //         std::cout << fNetKaFirstMoments[1][imom][icent][ieta] << std::endl;
        //         std::cout << fNetPrFirstMoments[1][imom][icent][ieta] << std::endl;
        //         std::cout << fNetLaFirstMoments[1][imom][icent][ieta] << std::endl;
        //         std::cout << fNetChFirstMoments[1][imom][icent][ieta] << std::endl;
        //         std::cout << " Info::marsland: ============================================================================================= " << std::endl;

        // net lambda for Alice and
        genMomentsNeg[kLaLa]=genMomentsNeg[kLa]*genMomentsNeg[kLa];
        genMomentsPos[kLaLa]=genMomentsPos[kLa]*genMomentsPos[kLa];
        genMomentsCross[kLaPosLaNeg]=genMomentsPos[kLa]*genMomentsNeg[kLa];

        genMomentsNeg[kChCh]=genMomentsNeg[kCh]*genMomentsNeg[kCh];
        genMomentsPos[kChCh]=genMomentsPos[kCh]*genMomentsPos[kCh];
        genMomentsCross[kChPosChNeg]=genMomentsPos[kCh]*genMomentsNeg[kCh];

        netPi[0]=(genMomentsPos[kPi]-genMomentsNeg[kPi]);
        netKa[0]=(genMomentsPos[kKa]-genMomentsNeg[kKa]);
        netPr[0]=(genMomentsPos[kPr]-genMomentsNeg[kPr]);
        netLa[0]=(genMomentsPos[kLa]-genMomentsNeg[kLa]);
        netCh[0]=(genMomentsPos[kCh]-genMomentsNeg[kCh]);

        netPi[1]=netPi[0]-fNetPiFirstMoments[0][imom][icent][ieta];
        netKa[1]=netKa[0]-fNetKaFirstMoments[0][imom][icent][ieta];
        netPr[1]=netPr[0]-fNetPrFirstMoments[0][imom][icent][ieta];
        netLa[1]=netLa[0]-fNetLaFirstMoments[0][imom][icent][ieta];
        netCh[1]=netCh[0]-fNetChFirstMoments[0][imom][icent][ieta];

        netPi[2]=netPi[1]*netPi[1];
        netKa[2]=netKa[1]*netKa[1];
        netPr[2]=netPr[1]*netPr[1];
        netLa[2]=netLa[1]*netLa[1];
        netCh[2]=netCh[1]*netCh[1];

        netPi[3]=netPi[1]*netPi[1]*netPi[1];
        netKa[3]=netKa[1]*netKa[1]*netKa[1];
        netPr[3]=netPr[1]*netPr[1]*netPr[1];
        netLa[3]=netLa[1]*netLa[1]*netLa[1];
        netCh[3]=netCh[1]*netCh[1]*netCh[1];

        netPi[4]=netPi[1]*netPi[1]*netPi[1]*netPi[1];
        netKa[4]=netKa[1]*netKa[1]*netKa[1]*netKa[1];
        netPr[4]=netPr[1]*netPr[1]*netPr[1]*netPr[1];
        netLa[4]=netLa[1]*netLa[1]*netLa[1]*netLa[1];
        netCh[4]=netCh[1]*netCh[1]*netCh[1]*netCh[1];


        // Moments without resonances
        nRgenMomentsNeg[kLaLa]=nRgenMomentsNeg[kLa]*nRgenMomentsNeg[kLa];
        nRgenMomentsPos[kLaLa]=nRgenMomentsPos[kLa]*nRgenMomentsPos[kLa];
        nRgenMomentsCross[kLaPosLaNeg]=nRgenMomentsPos[kLa]*nRgenMomentsNeg[kLa];

        nRgenMomentsNeg[kChCh]=nRgenMomentsNeg[kCh]*nRgenMomentsNeg[kCh];
        nRgenMomentsPos[kChCh]=nRgenMomentsPos[kCh]*nRgenMomentsPos[kCh];
        nRgenMomentsCross[kChPosChNeg]=nRgenMomentsPos[kCh]*nRgenMomentsNeg[kCh];

        noResnetPi[0]=(nRgenMomentsPos[kPi]-nRgenMomentsNeg[kPi]);
        noResnetKa[0]=(nRgenMomentsPos[kKa]-nRgenMomentsNeg[kKa]);
        noResnetPr[0]=(nRgenMomentsPos[kPr]-nRgenMomentsNeg[kPr]);
        noResnetLa[0]=(nRgenMomentsPos[kLa]-nRgenMomentsNeg[kLa]);
        noResnetCh[0]=(nRgenMomentsPos[kCh]-nRgenMomentsNeg[kCh]);

        noResnetPi[1]=noResnetPi[0]-fNetPiFirstMoments[1][imom][icent][imom];
        noResnetKa[1]=noResnetKa[0]-fNetKaFirstMoments[1][imom][icent][imom];
        noResnetPr[1]=noResnetPr[0]-fNetPrFirstMoments[1][imom][icent][imom];
        noResnetLa[1]=noResnetLa[0]-fNetLaFirstMoments[1][imom][icent][imom];
        noResnetCh[1]=noResnetCh[0]-fNetChFirstMoments[1][imom][icent][imom];

        noResnetPi[2]=noResnetPi[1]*noResnetPi[1];
        noResnetKa[2]=noResnetKa[1]*noResnetKa[1];
        noResnetPr[2]=noResnetPr[1]*noResnetPr[1];
        noResnetLa[2]=noResnetLa[1]*noResnetLa[1];
        noResnetCh[2]=noResnetCh[1]*noResnetCh[1];

        noResnetPi[3]=noResnetPi[1]*noResnetPi[1]*noResnetPi[1];
        noResnetKa[3]=noResnetKa[1]*noResnetKa[1]*noResnetKa[1];
        noResnetPr[3]=noResnetPr[1]*noResnetPr[1]*noResnetPr[1];
        noResnetLa[3]=noResnetLa[1]*noResnetLa[1]*noResnetLa[1];
        noResnetCh[3]=noResnetCh[1]*noResnetCh[1]*noResnetCh[1];

        noResnetPi[4]=noResnetPi[1]*noResnetPi[1]*noResnetPi[1]*noResnetPi[1];
        noResnetKa[4]=noResnetKa[1]*noResnetKa[1]*noResnetKa[1]*noResnetKa[1];
        noResnetPr[4]=noResnetPr[1]*noResnetPr[1]*noResnetPr[1]*noResnetPr[1];
        noResnetLa[4]=noResnetLa[1]*noResnetLa[1]*noResnetLa[1]*noResnetLa[1];
        noResnetCh[4]=noResnetCh[1]*noResnetCh[1]*noResnetCh[1]*noResnetCh[1];

        // fill tree which contains moments
        if(!fTreeSRedirector) return;
        // if there is at least one track in an event fill the tree
        if ( trCountgen>0 ){
          (*fTreeSRedirector)<<"mcGenMoms"<<
          "trCount="      << nTracksgen <<              // number of identified tracks within the given cent and mom range
          "isample="      << sampleNo <<                // sample id for subsample method
          "centDown="     << fcentDownArr[icent] <<     // lower edge of cent bin
          "centUp="       << fcentUpArr[icent] <<       // upper edge of cent bin
          "centBin="      << centBin <<                 // cent bin
          "vZ="           << fVz <<
          "impPar="       << fMCImpactParameter <<      // impact parameter taken from MC event header
          "pDown="        << fpDownArr[imom] <<         // lower edge of momentum bin
          "pUp="          << fpUpArr[imom] <<           // upper edge of momentum bin
          "etaDown="      << fetaDownArr[ieta] <<       // lower edge of eta bin
          "etaUp="        << fetaUpArr[ieta] <<         // upper edge of eta bin
          "momPos.="      << &genMomentsPos <<          // second moment of positive particles
          "momNeg.="      << &genMomentsNeg <<          // second moment of negative particles
          "noResmomPos.=" << &nRgenMomentsPos <<     // second moment of positive particles
          "noResmomNeg.=" << &nRgenMomentsNeg <<     // second moment of negative particles
          "netPi.="       << &netPi <<                  // second moments for particle+antiparticle
          "netKa.="       << &netKa <<                  // second moment of positive particles
          "netPr.="       << &netPr <<                  // second moment of negative particles
          "netLa.="       << &netLa <<                  // second moment of unlikesign particles
          "netCh.="       << &netCh <<                  // second moment of unlikesign particles
          "noResnetPi.="  << &noResnetPi <<             // second moments for particle+antiparticle
          "noResnetKa.="  << &noResnetKa <<             // second moment of positive particles
          "noResnetPr.="  << &noResnetPr <<             // second moment of negative particles
          "noResnetLa.="  << &noResnetLa <<             // second moment of unlikesign particles
          "noResnetCh.="  << &noResnetCh <<             // second moment of unlikesign particles
          "\n";
        } // tree filling

      } // ======= end of Centrality loop =======
    }// ======= end of momentum loop =======
  } // ======= end of eta loop =======
  //
  // ======================================================================
  //
  // ======================================================================

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::MCclosureHigherMoments()
{

  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  Int_t sampleNo = 0;
  Int_t nSubSample = 20;
  sampleNo = Int_t(fEventGID)%nSubSample;
  Int_t evtNuminFile = fMCEvent -> GetEventNumberInFile();
  if (fUseCouts) {
    std::cout << " Info::marsland: ===== In the MCclosureHigherMoments ===== sampleNo = " << sampleNo;
    std::cout << " gid = " << fEventGID << " evtNuminFile " << evtNuminFile << std::endl;
  }
  //
  // ======================================================================
  // ======================================================================
  //
  // ========= Efficiency Check Eta momentum and Centrality scan ==========
  const Int_t nHighMoments = 5;
  const Int_t nMoments = 3;
  for (Int_t ieta=0; ieta<fNEtaWinBinsMC; ieta++){
    for (Int_t imom=0; imom<fNMomBinsMC; imom++){
      for (Int_t icent=0; icent<fNCentbinsData-1; icent++){

        Double_t centBin = (fcentDownArr[icent]+fcentUpArr[icent])/2.;
        Float_t nTracksgen=0, nTracksrec=0, trCountgen=0, trCountrec=0;
        // vectors to hold moments
        TVectorF netPiGen(nHighMoments);   TVectorF piPosGen(nHighMoments);   TVectorF piNegGen(nHighMoments);
        TVectorF netKaGen(nHighMoments);   TVectorF kaPosGen(nHighMoments);   TVectorF kaNegGen(nHighMoments);
        TVectorF netPrGen(nHighMoments);   TVectorF prPosGen(nHighMoments);   TVectorF prNegGen(nHighMoments);
        TVectorF netPiRec(nHighMoments);   TVectorF piPosRec(nHighMoments);   TVectorF piNegRec(nHighMoments);
        TVectorF netKaRec(nHighMoments);   TVectorF kaPosRec(nHighMoments);   TVectorF kaNegRec(nHighMoments);
        TVectorF netPrRec(nHighMoments);   TVectorF prPosRec(nHighMoments);   TVectorF prNegRec(nHighMoments);
        for(Int_t i=0;i<nHighMoments; i++){
          netPiGen[i]=0.;  piPosGen[i]=0.;  piNegGen[i]=0.;
          netKaGen[i]=0.;  kaPosGen[i]=0.;  kaNegGen[i]=0.;
          netPrGen[i]=0.;  prPosGen[i]=0.;  prNegGen[i]=0.;
          netPiRec[i]=0.;  piPosRec[i]=0.;  piNegRec[i]=0.;
          netKaRec[i]=0.;  kaPosRec[i]=0.;  kaNegRec[i]=0.;
          netPrRec[i]=0.;  prPosRec[i]=0.;  prNegRec[i]=0.;
        }
        // vectors to hold moments
        TVectorF genMomentsPos(nMoments);     TVectorF recMomentsPos(nMoments);
        TVectorF genMomentsNeg(nMoments);     TVectorF recMomentsNeg(nMoments);
        TVectorF genMomentsCross(nMoments);   TVectorF recMomentsCross(nMoments);
        for(Int_t i=0;i<nMoments; i++){
          genMomentsPos[i]=0.;    recMomentsPos[i]=0.;
          genMomentsNeg[i]=0.;    recMomentsNeg[i]=0.;
          genMomentsCross[i]=0.;  recMomentsCross[i]=0.;
        }
        //
        // ************************************************************************
        //   Constructed track counters
        // ************************************************************************
        //
        for(Int_t iTrack = 0; iTrack < fESD->GetNumberOfTracks(); iTrack++)
        { // track loop
          //
          // initialize the dummy particle id
          Int_t piMCrec =-100., kaMCrec =-100., prMCrec =-100.;
          //
          // Esd track
          fTrackCutBits=0;  // reset the bits for the next track
          AliESDtrack *trackReal = fESD->GetTrack(iTrack);
          if (trackReal==NULL) continue;
          // Get generated track info
          Int_t lab = TMath::Abs(trackReal->GetLabel());
          if (!fMCStack->IsPhysicalPrimary(lab))continue;
          TParticle *trackMC  = fMCStack->Particle(lab);
          Int_t pdg = trackMC->GetPdgCode();
          //
          // acceptance cuts
          Double_t ptotMCrec = trackReal->P();
          Double_t etaMCrec  = trackReal->Eta();
          if (etaMCrec<fetaDownArr[ieta] || etaMCrec>fetaUpArr[ieta]) continue;
          if (ptotMCrec<fpDownArr[imom]  || ptotMCrec>fpUpArr[imom]) continue;
          //
          // detector cuts
          Bool_t ifDefaultCuts = fESDtrackCuts->AcceptTrack(trackReal);
          Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(trackReal);
          if (!trackReal -> GetInnerParam()) continue;
          if (!ifDefaultCuts) continue;  // real track cuts
          if (!ifDCAcutIfNoITSPixel) continue;  // TODO
          //
          Int_t iPart = -10;
          if (TMath::Abs(pdg) == kPDGpi) {iPart = 1; piMCrec = iPart;} // select pi+
          if (TMath::Abs(pdg) == kPDGka) {iPart = 2; kaMCrec = iPart;} // select ka+
          if (TMath::Abs(pdg) == kPDGpr) {iPart = 3; prMCrec = iPart;} // select pr+
          if (iPart == -10) continue;
          // count first moments
          if ( (fCentrality>=fcentDownArr[icent]) && (fCentrality<fcentUpArr[icent]) )
          {
            nTracksrec++;
            // Count identified particles
            if ( piMCrec>-1 || kaMCrec>-1 || prMCrec>-1) trCountrec++;
            //
            if ( piMCrec>-1 && pdg<0) recMomentsNeg[kPi]++;
            if ( kaMCrec>-1 && pdg<0) recMomentsNeg[kKa]++;
            if ( prMCrec>-1 && pdg<0) recMomentsNeg[kPr]++;
            //
            if ( piMCrec>-1 && pdg>0) recMomentsPos[kPi]++;
            if ( kaMCrec>-1 && pdg>0) recMomentsPos[kKa]++;
            if ( prMCrec>-1 && pdg>0) recMomentsPos[kPr]++;
          }

        } // ======= end of rec track loop =======
        //
        // ************************************************************************
        //   Generated track counters
        // ************************************************************************
        //
        for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++)
        { // track loop
          //
          // Select real trigger event and reject other pile up vertices
          Bool_t isTPCPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iTrack,fMCEvent);
          Bool_t isITSPileup = AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCEvent, "Hijing");
          if (isTPCPileup || isITSPileup) continue;
          //
          // initialize the dummy particle id
          Int_t piMCgen =-100., kaMCgen =-100., prMCgen =-100.;
          AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
          Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
          //
          // apply primary track and acceptance cuts
          Double_t ptotMCgen = trackMCgen->P();
          Double_t etaMCgen = (fRapidityType==0) ? trackMCgen->Eta() :  trackMCgen->Y();
          if (etaMCgen<fetaDownArr[ieta]  || etaMCgen>fetaUpArr[ieta]) continue;
          if (ptotMCgen<fpDownArr[imom]   || ptotMCgen>fpUpArr[imom]) continue;
          if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
          //
          // select particle of interest
          Int_t iPart = -10;
          if (TMath::Abs(pdg) == kPDGpi) {iPart = 1; piMCgen = iPart;} // select pi+
          if (TMath::Abs(pdg) == kPDGka) {iPart = 2; kaMCgen = iPart;} // select ka+
          if (TMath::Abs(pdg) == kPDGpr) {iPart = 3; prMCgen = iPart;} // select pr+
          if (iPart == -10) continue;
          // count first moments
          if ( (fCentrality>=fcentDownArr[icent]) && (fCentrality<fcentUpArr[icent]) )
          {
            nTracksgen++;
            // Count identified particles
            if ( piMCgen>-1 || kaMCgen>-1 || prMCgen>-1) trCountgen++;
            //
            if ( piMCgen>-1 && pdg<0) genMomentsNeg[kPi]++;
            if ( kaMCgen>-1 && pdg<0) genMomentsNeg[kKa]++;
            if ( prMCgen>-1 && pdg<0) genMomentsNeg[kPr]++;
            //
            if ( piMCgen>-1 && pdg>0) genMomentsPos[kPi]++;
            if ( kaMCgen>-1 && pdg>0) genMomentsPos[kKa]++;
            if ( prMCgen>-1 && pdg>0) genMomentsPos[kPr]++;
          }
        } // ======= end of track loop =======
        //
        // net particle cumulants  ---- GEN                                // net particle cumulants  ---- REC
        netPiGen[0]=(genMomentsPos[kPi]-genMomentsNeg[kPi]);               netPiRec[0]=(recMomentsPos[kPi]-recMomentsNeg[kPi]);
        netKaGen[0]=(genMomentsPos[kKa]-genMomentsNeg[kKa]);               netKaRec[0]=(recMomentsPos[kKa]-recMomentsNeg[kKa]);
        netPrGen[0]=(genMomentsPos[kPr]-genMomentsNeg[kPr]);               netPrRec[0]=(recMomentsPos[kPr]-recMomentsNeg[kPr]);
        netPiGen[1]=netPiGen[0]-fNetPiFirstMomentsGen[imom][icent][ieta];  netPiRec[1]=netPiRec[0]-fNetPiFirstMomentsRec[imom][icent][ieta];
        netKaGen[1]=netKaGen[0]-fNetKaFirstMomentsGen[imom][icent][ieta];  netKaRec[1]=netKaRec[0]-fNetKaFirstMomentsRec[imom][icent][ieta];
        netPrGen[1]=netPrGen[0]-fNetPrFirstMomentsGen[imom][icent][ieta];  netPrRec[1]=netPrRec[0]-fNetPrFirstMomentsRec[imom][icent][ieta];
        netPiGen[2]=netPiGen[1]*netPiGen[1];                               netPiRec[2]=netPiRec[1]*netPiRec[1];
        netKaGen[2]=netKaGen[1]*netKaGen[1];                               netKaRec[2]=netKaRec[1]*netKaRec[1];
        netPrGen[2]=netPrGen[1]*netPrGen[1];                               netPrRec[2]=netPrRec[1]*netPrRec[1];
        netPiGen[3]=netPiGen[1]*netPiGen[1]*netPiGen[1];                   netPiRec[3]=netPiRec[1]*netPiRec[1]*netPiRec[1];
        netKaGen[3]=netKaGen[1]*netKaGen[1]*netKaGen[1];                   netKaRec[3]=netKaRec[1]*netKaRec[1]*netKaRec[1];
        netPrGen[3]=netPrGen[1]*netPrGen[1]*netPrRec[1];                   netPrRec[3]=netPrRec[1]*netPrRec[1]*netPrRec[1];
        netPiGen[4]=netPiGen[1]*netPiGen[1]*netPiGen[1]*netPiGen[1];       netPiRec[4]=netPiRec[1]*netPiRec[1]*netPiRec[1]*netPiRec[1];
        netKaGen[4]=netKaGen[1]*netKaGen[1]*netKaGen[1]*netKaGen[1];       netKaRec[4]=netKaRec[1]*netKaRec[1]*netKaRec[1]*netKaRec[1];
        netPrGen[4]=netPrGen[1]*netPrGen[1]*netPrGen[1]*netPrGen[1];       netPrRec[4]=netPrRec[1]*netPrRec[1]*netPrRec[1]*netPrRec[1];
        //
        // particle cumulants  ---- GEN                                    // particle cumulants  ---- REC
        piPosGen[0]=genMomentsPos[kPi];                                    piPosRec[0]=recMomentsPos[kPi];
        kaPosGen[0]=genMomentsPos[kKa];                                    kaPosRec[0]=recMomentsPos[kKa];
        prPosGen[0]=genMomentsPos[kPr];                                    prPosRec[0]=recMomentsPos[kPr];
        piPosGen[1]=piPosGen[0]-fPiFirstMomentsGen[0][imom][icent][ieta];  piPosRec[1]=piPosRec[0]-fPiFirstMomentsRec[0][imom][icent][ieta];
        kaPosGen[1]=kaPosGen[0]-fKaFirstMomentsGen[0][imom][icent][ieta];  kaPosRec[1]=kaPosRec[0]-fKaFirstMomentsRec[0][imom][icent][ieta];
        prPosGen[1]=prPosGen[0]-fPrFirstMomentsGen[0][imom][icent][ieta];  prPosRec[1]=prPosRec[0]-fPrFirstMomentsRec[0][imom][icent][ieta];
        piPosGen[2]=piPosGen[1]*piPosGen[1];                               piPosRec[2]=piPosRec[1]*piPosRec[1];
        kaPosGen[2]=kaPosGen[1]*kaPosGen[1];                               kaPosRec[2]=kaPosRec[1]*kaPosRec[1];
        prPosGen[2]=prPosGen[1]*prPosGen[1];                               prPosRec[2]=prPosRec[1]*prPosRec[1];
        piPosGen[3]=piPosGen[1]*piPosGen[1]*piPosGen[1];                   piPosRec[3]=piPosRec[1]*piPosRec[1]*piPosRec[1];
        kaPosGen[3]=kaPosGen[1]*kaPosGen[1]*kaPosGen[1];                   kaPosRec[3]=kaPosRec[1]*kaPosRec[1]*kaPosRec[1];
        prPosGen[3]=prPosGen[1]*prPosGen[1]*prPosGen[1];                   prPosRec[3]=prPosRec[1]*prPosRec[1]*prPosRec[1];
        //
        // Anti particle cumulants  ---- GEN                               // Anti particle cumulants  ---- REC
        piNegGen[0]=genMomentsNeg[kPi];                                    piNegRec[0]=recMomentsNeg[kPi];
        kaNegGen[0]=genMomentsNeg[kKa];                                    kaNegRec[0]=recMomentsNeg[kKa];
        prNegGen[0]=genMomentsNeg[kPr];                                    prNegRec[0]=recMomentsNeg[kPr];
        piNegGen[1]=piNegGen[0]-fPiFirstMomentsGen[1][imom][icent][ieta];  piNegRec[1]=piNegRec[0]-fPiFirstMomentsRec[1][imom][icent][ieta];
        kaNegGen[1]=kaNegGen[0]-fKaFirstMomentsGen[1][imom][icent][ieta];  kaNegRec[1]=kaNegRec[0]-fKaFirstMomentsRec[1][imom][icent][ieta];
        prNegGen[1]=prNegGen[0]-fPrFirstMomentsGen[1][imom][icent][ieta];  prNegRec[1]=prNegRec[0]-fPrFirstMomentsRec[1][imom][icent][ieta];
        piNegGen[2]=piNegGen[1]*piNegGen[1];                               piNegRec[2]=piNegRec[1]*piNegRec[1];
        kaNegGen[2]=kaNegGen[1]*kaNegGen[1];                               kaNegRec[2]=kaNegRec[1]*kaNegRec[1];
        prNegGen[2]=prNegGen[1]*prNegGen[1];                               prNegRec[2]=prNegRec[1]*prNegRec[1];
        piNegGen[3]=piNegGen[1]*piNegGen[1]*piNegGen[1];                   piNegRec[3]=piNegRec[1]*piNegRec[1]*piNegRec[1];
        kaNegGen[3]=kaNegGen[1]*kaNegGen[1]*kaNegGen[1];                   kaNegRec[3]=kaNegRec[1]*kaNegRec[1]*kaNegRec[1];
        prNegGen[3]=prNegGen[1]*prNegGen[1]*prNegGen[1];                   prNegRec[3]=prNegRec[1]*prNegRec[1]*prNegRec[1];
        //
        // fill tree which contains moments
        if(!fTreeSRedirector) return;
        (*fTreeSRedirector)<<"mcMoms"<<
        "gid="          << fEventGID <<
        "isample="      << sampleNo <<                // sample id for subsample method
        "impPar="       << fMCImpactParameter <<      // impact parameter taken from MC event header
        "centDown="     << fcentDownArr[icent] <<     // lower edge of cent bin
        "centUp="       << fcentUpArr[icent] <<       // upper edge of cent bin
        "centBin="      << centBin <<                 // cent bin
        "pDown="        << fpDownArr[imom] <<         // lower edge of momentum bin
        "pUp="          << fpUpArr[imom] <<           // upper edge of momentum bin
        "etaDown="      << fetaDownArr[ieta] <<       // lower edge of eta bin
        "etaUp="        << fetaUpArr[ieta] <<         // upper edge of eta bin
        //
        "trCountGen="   << nTracksgen <<              // number of identified tracks within the given cent and mom range
        "trCountRec="   << nTracksrec <<              // number of identified tracks within the given cent and mom range
        //
        "momPosGen.="   << &genMomentsPos <<          // second moment of positive particles
        "momNegGen.="   << &genMomentsNeg <<          // second moment of negative particles
        "netPiGen.="    << &netPiGen <<               // second moments for particle+antiparticle
        "netKaGen.="    << &netKaGen <<               // second moment of positive particles
        "netPrGen.="    << &netPrGen <<               // second moment of negative particles
        //
        "piPosGen.="    << &piPosGen <<               // second moments for particle+antiparticle
        "kaPosGen.="    << &kaPosGen <<               // second moment of positive particles
        "prPosGen.="    << &prPosGen <<               // second moment of negative particles
        "piNegGen.="    << &piNegGen <<               // second moments for particle+antiparticle
        "kaNegGen.="    << &kaNegGen <<               // second moment of positive particles
        "prNegGen.="    << &prNegGen <<               // second moment of negative particles
        //
        "momPosRec.="   << &recMomentsPos <<          // second moment of positive particles
        "momNegRec.="   << &recMomentsNeg <<          // second moment of negative particles
        "netPiRec.="    << &netPiRec <<               // second moments for particle+antiparticle
        "netKaRec.="    << &netKaRec <<               // second moment of positive particles
        "netPrRec.="    << &netPrRec <<               // second moment of negative particles
        //
        "piPosRec.="    << &piPosRec <<               // second moments for particle+antiparticle
        "kaPosRec.="    << &kaPosRec <<               // second moment of positive particles
        "prPosRec.="    << &prPosRec <<               // second moment of negative particles
        "piNegRec.="    << &piNegRec <<               // second moments for particle+antiparticle
        "kaNegRec.="    << &kaNegRec <<               // second moment of positive particles
        "prNegRec.="    << &prNegRec <<               // second moment of negative particles
        "\n";

      } // ======= end of Centrality loop =======
    }// ======= end of momentum loop =======
  } // ======= end of eta loop =======
  //
  // ======================================================================
  //
  // ======================================================================

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FillEffMatrix()
{

  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  Int_t trackOrigin = -10;
  Int_t sampleNo = 0;
  Int_t nSubSample = 20;
  sampleNo = Int_t(fEventGID)%nSubSample;
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillEffMatrix ===== " << std::endl;
  //
  // -----------------------------------------------------------------------------------------
  // ----------------------------   reconstructed MC particles  ------------------------------
  // -----------------------------------------------------------------------------------------
  //
  //
  Bool_t bEventVertexZSmall      = (TMath::Abs(fVz)<6 && TMath::Abs(fVz)>0.2);
  Bool_t bCutReference           = (TMath::Abs(fVz)<7 && TMath::Abs(fVz)>0.15);
  Bool_t bEventVertexZLarge      = (TMath::Abs(fVz)<8 && TMath::Abs(fVz)>0.1);
  Bool_t bEventVertexZALICE      = (TMath::Abs(fVz)<10);
  Bool_t bEventVertexZALICETight = (TMath::Abs(fVz)<7);
  //
  // loop over tracks
  for (Int_t iset=0; iset<fNSettings; iset++) {
    //
    // event Vz cuts
    if      (iset==kCutEventVertexZSmall      && !bEventVertexZSmall)      continue;
    else if (iset==kCutEventVertexZLarge      && !bEventVertexZLarge)      continue;
    else if (iset==kCutEventVertexZALICE      && !bEventVertexZALICE)      continue;
    else if (iset==kCutEventVertexZALICETight && !bEventVertexZALICETight) continue;
    else if ( !bCutReference ) continue;
    //
    for(Int_t irectrack = 0; irectrack < fESD->GetNumberOfTracks(); irectrack++)
    { // track loop
      //
      fTrackCutBits=0;  // reset the bits for the next track
      AliESDtrack *trackReal = fESD->GetTrack(irectrack);
      if (trackReal==NULL) continue;
      Int_t lab = TMath::Abs(trackReal->GetLabel()); // avoid from negatif labels, they include some garbage

      //
      Float_t etaMCrec = (fRapidityType==0) ? trackReal->Eta() :  trackReal->Y();
      if ((etaMCrec<fEtaDown) || (etaMCrec>fEtaUp)) continue;   // eta [-0.8,0.8]
      if (!fMCStack->IsPhysicalPrimary(lab)) continue;   // MC primary track check
      //
      Bool_t ifDefaultCuts = fESDtrackCuts->AcceptTrack(trackReal);
      Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(trackReal);
      if (!trackReal->GetInnerParam()) continue;   // If track in TPC
      if (!ifDefaultCuts) continue;   // apply esdtrack cuts
      if (!ifDCAcutIfNoITSPixel) continue;   // apply esdtrack cuts
      //
      // access generated track info
      AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(lab);
      Int_t pdg = trackMCgen->Particle()->GetPdgCode();  //  Int_t pdg = trackMC->GetPdgCode();   TODO
      //
      // get track info
      Float_t fPtRec = 0.;
      if(fUsePtCut==1) fPtRec = trackReal->P();
      if(fUsePtCut==2) fPtRec = trackReal->Pt();
      if ((fPtRec<fMomDown) || (fPtRec>fMomUp)) continue;
      //
      Int_t fPartID  = -10;
      if (TMath::Abs(pdg) == kPDGpi) fPartID=0; // select pi
      if (TMath::Abs(pdg) == kPDGka) fPartID=1; // select ka
      if (TMath::Abs(pdg) == kPDGpr) fPartID=2; // select pr
      if (fPartID == -10) continue;
      //
      // Loop over all track settings
      Float_t closestPar[3];
      GetExpecteds(trackReal,closestPar);
      SetCutBitsAndSomeTrackVariables(trackReal,fPartID);
      //
      if (GetSystematicClassIndex(fTrackCutBits,iset)){
        //
        // Fill TPC eff matrix
        Double_t xxxRec[6]={0.,Float_t(iset),Float_t(fPartID),fCentrality,fPtRec,etaMCrec};
        if (pdg>0) fHistPosEffMatrixScanRec->Fill(xxxRec);
        if (pdg<0) fHistNegEffMatrixScanRec->Fill(xxxRec);
        //
        // Fill TOF eff matrix
        Bool_t piTOF = (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trackReal, AliPID::kPion))  <=2.5);
        Bool_t kaTOF = (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trackReal, AliPID::kKaon))  <=2.5);
        Bool_t prTOF = (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trackReal, AliPID::kProton))<=2.5);
        if ( (piTOF && fPartID==0) ||  (kaTOF && fPartID==1) || (prTOF && fPartID==2) ) {
          Double_t xxxRecTOF[6]={1.,Float_t(iset),Float_t(fPartID),fCentrality,fPtRec,etaMCrec};
          if (pdg>0) fHistPosEffMatrixScanRec->Fill(xxxRecTOF);
          if (pdg<0) fHistNegEffMatrixScanRec->Fill(xxxRecTOF);
        }
      }

    } // ======= end of rec track loop =======
    //
    // -----------------------------------------------------------------------------------------
    // ----------------------------   MC generated pure MC particles  --------------------------
    // -----------------------------------------------------------------------------------------
    //
    // only event Vz cuts can be applied on generated level
    if ( iset==kCutReference || iset==kCutEventVertexZSmall || iset==kCutEventVertexZLarge || iset==kCutEventVertexZALICE || iset==kCutEventVertexZALICETight ) {
      //
      for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++)
      { // track loop
        //
        // Select real trigger event and reject other pile up vertices
        Bool_t isTPCPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iTrack,fMCEvent);
        Bool_t isITSPileup = AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCEvent, "Hijing");
        if (isTPCPileup || isITSPileup) continue;
        //
        AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
        if (!trackMCgen) continue;
        //
        // apply primary track and acceptance cuts
        if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
        //
        // Aplly eta acceptance
        Float_t etaMCgen = (fRapidityType==0) ? trackMCgen->Eta() :  trackMCgen->Y();
        if ((etaMCgen<fEtaDown) || (etaMCgen>fEtaUp))  continue;
        //
        // Aplly momentum acceptance
        Float_t ptotMCgen =0.;
        if(fUsePtCut==1) ptotMCgen = trackMCgen->P();
        if(fUsePtCut==2) ptotMCgen = trackMCgen->Pt();
        if ((ptotMCgen<fMomDown) || (ptotMCgen>fMomUp)) continue;
        //
        // Efficiency matices for individual particles
        Int_t fPartID  = -10;
        Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
        if (TMath::Abs(pdg) == kPDGpi) fPartID=0; // select pi
        if (TMath::Abs(pdg) == kPDGka) fPartID=1; // select ka
        if (TMath::Abs(pdg) == kPDGpr) fPartID=2; // select pr
        if (fPartID == -10) continue;
        //
        Double_t xxxGen[6]={0.,Float_t(iset),Float_t(fPartID),fCentrality,ptotMCgen,etaMCgen};
        if (pdg>0) fHistPosEffMatrixScanGen->Fill(xxxGen);
        if (pdg<0) fHistNegEffMatrixScanGen->Fill(xxxGen);
        // generated does not know about TOF
        Double_t xxxGenTOF[6]={1.,Float_t(iset),Float_t(fPartID),fCentrality,ptotMCgen,etaMCgen};
        if (pdg>0) fHistPosEffMatrixScanGen->Fill(xxxGenTOF);
        if (pdg<0) fHistNegEffMatrixScanGen->Fill(xxxGenTOF);

      } // ======= end of gen track loop =======
    }
    //
    if (fDefaultTrackCuts && iset==kCutReference) break; // for default track cuts run only once

  }


  /*
  //
  // -----------------------------------------------------------------------------------------
  // ----------------------------   MC generated pure MC particles  --------------------------
  // -----------------------------------------------------------------------------------------
  //
  for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++)
  { // track loop
  AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
  if (!trackMCgen) {fTrackCutBits=0; continue;}
  //
  // check the origin of the track
  trackOrigin = -10;
  if (fMCStack->IsPhysicalPrimary(iTrack))        trackOrigin = 0;
  if (fMCStack->IsSecondaryFromMaterial(iTrack))  trackOrigin = 1;
  if (fMCStack->IsSecondaryFromWeakDecay(iTrack)) trackOrigin = 2;
  if (trackOrigin<-1) {fTrackCutBits=0; continue;}
  //
  // Aplly eta acceptance
  if ((trackMCgen->Eta()<-1) || (trackMCgen->Eta()>1)) {fTrackCutBits=0; continue;}
  if ((trackMCgen->P()<0.1)  || (trackMCgen->P()>4))   {fTrackCutBits=0; continue;}
  //
  // Efficiency matices for individual particles
  Int_t fPartID  = -10;
  Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
  if (TMath::Abs(pdg) == kPDGpi) fPartID=0; // select pi
  if (TMath::Abs(pdg) == kPDGka) fPartID=1; // select ka
  if (TMath::Abs(pdg) == kPDGpr) fPartID=2; // select pr
  if (fPartID == -10) {fTrackCutBits=0; continue;}
  //
  // Get real tracks:
  AliESDtrack *trackReal = NULL;
  trackReal = fESD->GetTrack(iTrack);
  fTPCSignal = 0; fTrackCutBits = 0; fSign = -100; fPtot = 0; fPVertex = 0; fPt = 0;
  fEta = 0; fPhi = 0; fTrackTPCCrossedRows=0; fTrackChi2TPC=0;
  Int_t tpcNcls = 0 ;
  if (trackReal){
  SetCutBitsAndSomeTrackVariables(trackReal);
  if (trackReal->GetInnerParam()){
  Float_t closestPar[3];
  GetExpecteds(trackReal,closestPar);
  tpcNcls = trackReal->GetTPCncls();
}
}
Float_t etagen = trackMCgen->Eta();
Float_t pgen   = trackMCgen->P();
Float_t pTgen  = trackMCgen->Pt();
Float_t phigen = trackMCgen->Phi();
Int_t signgen  = (trackMCgen->Charge()<0) ? -1 : 1;

if(!fTreeSRedirector) return;
(*fTreeSRedirector)<<"mcGen"<<
"cent="      << fCentrality <<     // Centrality
"centimp="   << fCentImpBin <<
"isample="   << sampleNo <<                // sample id for subsample method
"origin="    << trackOrigin <<
"part="      << fPartID <<
"gid="       << fEventGID <<  //  global event ID
"vZ="        << fVz <<
//
"dEdx="      << fTPCSignal <<    // dEdx of mc track
"cutBit="    << fTrackCutBits <<  //  Systematic Cuts
"sign="      << fSign <<         // sign
"ptot="      << fPtot <<         // tpc momentum
"p="         << fPVertex <<             // vertex momentum
"pT="        << fPt <<           // transverse momentum
"eta="       << fEta <<          // mc eta
"phi="       << fPhi <<          // mc eta
"dcaxy="     << fTrackDCAxy <<
"dcaz="      << fTrackDCAz <<
"nclTPC="    << tpcNcls               <<  //  centrality
"cRows="     << fTrackTPCCrossedRows  <<
"chi2tpc="   << fTrackChi2TPC         <<
"nsigmatofka=" << fNSigmasKaTOF         <<  // interaction rate
"nsigmatofpr=" << fNSigmasPrTOF         <<  // interaction rate
//
"etagen="    << etagen <<          // mc eta
"phigen="    << phigen <<          // mc eta
"pgen="      << pgen <<          // mc eta
"pTgen="     << pTgen <<          // mc eta
"signgen="   << signgen <<
"\n";
//
fTrackCutBits=0;  // reset the bits for the next track
} // ======= end of track loop =======
*/

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FillCleanSamples()
{

  // Fill Clean Pions from K0s
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillCleanSamples ===== " << std::endl;
  if (fPIDResponse) {
    fPIDResponse->GetTPCResponse().SetBetheBlochParameters(1.28778e+00/50., 3.13539e+01, TMath::Exp(-3.16327e+01), 1.87901e+00, 6.41583e+00);
  }
  AliKFParticle::SetField(fESD->GetMagneticField());
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};
  Double_t mm[3] = {0,0,0};
  const Double_t cProtonMass  =TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  const Double_t cPionMass    =TDatabasePDG::Instance()->GetParticle(211)->Mass();
  const Double_t cElectronMass=TDatabasePDG::Instance()->GetParticle(11)->Mass();
  //
  // Selection From Ionut
  const AliESDVertex *primaryVertex = fESD->GetPrimaryVertex();
  AliKFVertex primaryVertexKF(*primaryVertex);
  if(fV0OpenCuts) {
    fV0OpenCuts->SetEvent(fESD);
    fV0OpenCuts->SetPrimaryVertex(&primaryVertexKF);
  }
  if(fV0StrongCuts) {
    fV0StrongCuts->SetEvent(fESD);
    fV0StrongCuts->SetPrimaryVertex(&primaryVertexKF);
  }
  //
  //
  TObjArray* listCrossV0 = fESDtrackCutsV0->GetAcceptedV0s(fESD);
  Int_t nGoodV0s         = listCrossV0->GetEntries();
  delete listCrossV0;
  //
  // Loop over V0s
  Int_t pdgV0=0; Int_t pdgP=0; Int_t pdgN=0;
  for(Int_t iV0MI = 0; iV0MI < nGoodV0s; iV0MI++) {
    //
    Int_t v0purity = 0;
    AliESDv0 * fV0s = fESD->GetV0(iV0MI);
    Int_t lOnFlyStatus = 0;
    lOnFlyStatus = fV0s->GetOnFlyStatus();
    if (!lOnFlyStatus) {fTrackCutBits=0; continue;}
    //
    AliESDtrack* trackPosTest = fESD->GetTrack(fV0s->GetPindex());
    AliESDtrack* trackNegTest = fESD->GetTrack(fV0s->GetNindex());
    //
    // ----------------------------------------------------------------------------------------------------------
    //  Selections from ionuts
    // ----------------------------------------------------------------------------------------------------------
    //
    if(trackPosTest->GetSign() == trackNegTest->GetSign()) {fTrackCutBits=0; continue;}
    Bool_t v0ChargesAreCorrect = (trackPosTest->GetSign()==+1 ? kTRUE : kFALSE);
    trackPosTest = (!v0ChargesAreCorrect ? fESD->GetTrack(fV0s->GetNindex()) : trackPosTest);
    trackNegTest = (!v0ChargesAreCorrect ? fESD->GetTrack(fV0s->GetPindex()) : trackNegTest);
    //
    Bool_t goodK0s = kTRUE, goodLambda = kTRUE, goodALambda = kTRUE, goodGamma = kTRUE;
    if(fV0OpenCuts) {
      goodK0s = kFALSE, goodLambda = kFALSE, goodALambda = kFALSE, goodGamma = kFALSE;
      pdgV0=0; pdgP=0; pdgN=0;
      Bool_t processV0 = fV0OpenCuts->ProcessV0(fV0s, pdgV0, pdgP, pdgN);
      if (processV0 && TMath::Abs(pdgV0)==310 &&  TMath::Abs(pdgP)==kPDGpi && TMath::Abs(pdgN)==kPDGpi) {
        goodK0s = kTRUE;
        if(fK0sPionCuts && (!fK0sPionCuts->IsSelected(trackPosTest) || !fK0sPionCuts->IsSelected(trackNegTest))) goodK0s = kFALSE;
      }
      if (processV0 && pdgV0== kPDGla         && (TMath::Abs(pdgP)==kPDGpi || TMath::Abs(pdgP)==kPDGpr) && (TMath::Abs(pdgN)==kPDGpi || TMath::Abs(pdgN)==kPDGpr)) {
        goodLambda = kTRUE;
        if(fLambdaProtonCuts && !fLambdaProtonCuts->IsSelected(trackPosTest)) goodLambda = kFALSE;
        if(fLambdaPionCuts && !fLambdaPionCuts->IsSelected(trackNegTest)) goodLambda = kFALSE;
      }
      if (processV0 && pdgV0==-kPDGla         && (TMath::Abs(pdgP)==kPDGpi || TMath::Abs(pdgP)==kPDGpr) && (TMath::Abs(pdgN)==kPDGpi || TMath::Abs(pdgN)==kPDGpr)) {
        goodALambda = kTRUE;
        if(fLambdaProtonCuts && !fLambdaProtonCuts->IsSelected(trackNegTest)) goodALambda = kFALSE;
        if(fLambdaPionCuts && !fLambdaPionCuts->IsSelected(trackPosTest)) goodALambda = kFALSE;
      }
      if (processV0 && TMath::Abs(pdgV0)==22  &&  TMath::Abs(pdgP)==kPDGel && TMath::Abs(pdgN)==kPDGel) {
        goodGamma = kTRUE;
        if(fGammaElectronCuts && (!fGammaElectronCuts->IsSelected(trackPosTest) || !fGammaElectronCuts->IsSelected(trackNegTest))) goodGamma = kFALSE;
      }
    }
    //
    Bool_t veryGoodK0s = kFALSE, veryGoodLambda = kFALSE, veryGoodALambda = kFALSE, veryGoodGamma = kFALSE;
    if(fV0StrongCuts && (goodK0s || goodLambda || goodALambda || goodGamma)) {
      pdgV0=0; pdgP=0; pdgN=0;
      Bool_t processV0 = fV0StrongCuts->ProcessV0(fV0s, pdgV0, pdgP, pdgN);
      if (processV0 && goodK0s     && TMath::Abs(pdgV0)==310 &&  TMath::Abs(pdgP)==kPDGpi && TMath::Abs(pdgN)==kPDGpi) veryGoodK0s = kTRUE;
      if (processV0 && goodLambda  && pdgV0== kPDGla         && (TMath::Abs(pdgP)==kPDGpi || TMath::Abs(pdgP)==kPDGpr) && (TMath::Abs(pdgN)==kPDGpi || TMath::Abs(pdgN)==kPDGpr)) veryGoodLambda = kTRUE;
      if (processV0 && goodALambda && pdgV0==-kPDGla         && (TMath::Abs(pdgP)==kPDGpi || TMath::Abs(pdgP)==kPDGpr) && (TMath::Abs(pdgN)==kPDGpi || TMath::Abs(pdgN)==kPDGpr)) veryGoodALambda = kTRUE;
      if (processV0 && goodGamma   && TMath::Abs(pdgV0)==22  &&  TMath::Abs(pdgP)==kPDGel && TMath::Abs(pdgN)==kPDGel) veryGoodGamma = kTRUE;
    }
    //
    //
    if( goodK0s || goodLambda || goodALambda || goodGamma ) v0purity = 1;
    if( veryGoodK0s || veryGoodLambda || veryGoodALambda || veryGoodGamma ) v0purity = 2;
    //
    // ----------------------------------------------------------------------------------------------------------
    //  My cuts
    // ----------------------------------------------------------------------------------------------------------
    //
    if (!fESDtrackCutsCleanSamp->AcceptTrack(trackPosTest)) {fTrackCutBits=0; continue;} // To FIX
    if (!fESDtrackCutsCleanSamp->AcceptTrack(trackNegTest)) {fTrackCutBits=0; continue;} // To FIX
    if (!trackPosTest->GetInnerParam()) {fTrackCutBits=0; continue;}
    if (!trackNegTest->GetInnerParam()) {fTrackCutBits=0; continue;}

    if( trackPosTest->GetSign() >0 && trackNegTest->GetSign() <0){
      fV0s->GetNPxPyPz(mn[0],mn[1],mn[2]); //reconstructed cartesian momentum components of negative daughter
      fV0s->GetPPxPyPz(mp[0],mp[1],mp[2]); //reconstructed cartesian momentum components of positive daughter
    }

    if( trackPosTest->GetSign() <0 && trackNegTest->GetSign() >0){
      fV0s->GetPPxPyPz(mn[0],mn[1],mn[2]); //reconstructed cartesian momentum components of negative daughter
      fV0s->GetNPxPyPz(mp[0],mp[1],mp[2]); //reconstructed cartesian momentum components of positive daughter
    }

    fV0s->GetPxPyPz(mm[0],mm[1],mm[2]); //reconstructed cartesian momentum components of mother

    TVector3 vecN(mn[0],mn[1],mn[2]);
    TVector3 vecP(mp[0],mp[1],mp[2]);
    TVector3 vecM(mm[0],mm[1],mm[2]);

    if ((vecP.Mag() * vecM.Mag())<0.00001) {fTrackCutBits=0; continue;}
    if ((vecN.Mag() * vecM.Mag())<0.00001) {fTrackCutBits=0; continue;}
    Double_t thetaP  = acos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
    Double_t thetaN  = acos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));
    if ( ((vecP.Mag())*cos(thetaP)+(vecN.Mag())*cos(thetaN)) <0.00001) {fTrackCutBits=0; continue;}
    fAlfa = ((vecP.Mag())*cos(thetaP)-(vecN.Mag())*cos(thetaN))/((vecP.Mag())*cos(thetaP)+(vecN.Mag())*cos(thetaN));
    fQt   = vecP.Mag()*sin(thetaP);
    if (fUseCouts) fHistArmPod->Fill(fAlfa,fQt);
    // fV0s->ChangeMassHypothesis(22);   // ????
    // fV0s->ChangeMassHypothesis(310); // ????
    //
    // main armentoros podolanki cuts
    if (TMath::Abs(fAlfa)>0.9) {fTrackCutBits=0; continue;}
    if (fQt >0.22) {fTrackCutBits=0; continue;}
    if (fQt >0.02 && fQt<0.12 && TMath::Abs(fAlfa)<0.4) {fTrackCutBits=0; continue;}
    SelectCleanSamplesFromV0s(fV0s,trackPosTest,trackNegTest);
    //
    TLorentzVector posE, negE, photon, posP, negP, posPi, negPi, lambda, antiLambda, kaon, posProton, k0sProton;
    negE.SetXYZM(mn[0],mn[1],mn[2],cElectronMass);
    posE.SetXYZM(mp[0],mp[1],mp[2],cElectronMass);
    negPi.SetXYZM(mn[0],mn[1],mn[2],cPionMass);
    posPi.SetXYZM(mp[0],mp[1],mp[2],cPionMass);
    negP.SetXYZM(mn[0],mn[1],mn[2],cProtonMass);
    posP.SetXYZM(mp[0],mp[1],mp[2],cProtonMass);
    kaon=posPi+negPi;
    photon=posE+negE;
    lambda=posP+negPi;
    antiLambda=posPi+negP;
    //
    //
    Float_t ptotForBetaGamma0 = trackPosTest->GetInnerParam()->GetP();
    Float_t ptotForBetaGamma1 = trackNegTest->GetInnerParam()->GetP();
    Float_t ptotForBetaGammaThr = 0.2;
    Double_t posNTPCSigmaPi = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPosTest, AliPID::kPion));
    Double_t negNTPCSigmaPi = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNegTest, AliPID::kPion));
    Double_t posNTPCSigmaPr = (ptotForBetaGamma0>ptotForBetaGammaThr) ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPosTest, AliPID::kProton)) : 0.;
    Double_t negNTPCSigmaPr = (ptotForBetaGamma1>ptotForBetaGammaThr) ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNegTest, AliPID::kProton)) : 0.;
    Double_t posNTPCSigmaEl = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPosTest, AliPID::kElectron));
    Double_t negNTPCSigmaEl = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNegTest, AliPID::kElectron));
    //
    // --------------------------------------------------------------
    //  Invariant mass cuts
    // --------------------------------------------------------------
    //
    Bool_t isK0sMass        = (kaon.M()>0.485 && kaon.M()<0.51); // (kaon.M()>0.490 && kaon.M()<0.504);
    Bool_t isLambdaMass     = (lambda.M()>1.112 && lambda.M()<1.119); //(lambda.M()>1.113 && lambda.M()<1.118);
    Bool_t isAntiLambdaMass = (antiLambda.M()>1.112 && antiLambda.M()<1.119); // (antiLambda.M()>1.113 && antiLambda.M()<1.118);
    Bool_t isPhotonMass     = (photon.M()<0.005); // (photon.M()<0.005);
    Double_t oneLegSigma = 3.5;
    if (fQt<0.02 && TMath::Abs(fAlfa)<0.5){
      // Cuts concerns Electrons
      //
      // Apply one leg cut for electrons
      if (!(negNTPCSigmaEl<2. || posNTPCSigmaEl<2.)) {fTrackCutBits=0; continue;}
      //
      if (fUseCouts) fHistInvPhoton->Fill(photon.M());
      if (isK0sMass) {fTrackCutBits=0; continue;}
      if (isLambdaMass) {fTrackCutBits=0; continue;}
      if (isAntiLambdaMass) {fTrackCutBits=0; continue;}
      if (!isPhotonMass) {fTrackCutBits=0; continue;}
    } else {
      if (fUseCouts) {
        fHistInvK0s->Fill(kaon.M());
        fHistInvLambda->Fill(lambda.M());
        fHistInvAntiLambda->Fill(antiLambda.M());
      }
      if ( !(isK0sMass || isLambdaMass || isAntiLambdaMass) ) {fTrackCutBits=0; continue;}
      //
      // Apply one leg cut for K0s
      if (fQt>0.11 && (!(negNTPCSigmaPi < oneLegSigma || posNTPCSigmaPi < oneLegSigma))) {fTrackCutBits=0; continue;}
      // //
      // // Apply one leg cut for antilambda
      // if (fQt<1.1 && fAlfa<0) {
      //   if (!(negNTPCSigmaPr < oneLegSigma || posNTPCSigmaPi < oneLegSigma)) {fTrackCutBits=0; continue;}
      //   if (negNTPCSigmaPr > 4) {fTrackCutBits=0; continue;}
      // }
      // //
      // // Apply one leg cut for lambda
      // if (fQt<1.1 && fAlfa>0) {
      //   if (!(negNTPCSigmaPi < oneLegSigma || posNTPCSigmaPr < oneLegSigma)) {fTrackCutBits=0; continue;}
      //   if (posNTPCSigmaPr > 4) {fTrackCutBits=0; continue;}
      // }
    }
    //
    // Set the variables to be filled in the tree
    Bool_t selectPosLegs = (posNTPCSigmaPi<3 || posNTPCSigmaPr<3 || posNTPCSigmaEl<3);
    Bool_t selectNegLegs = (negNTPCSigmaPi<3 || negNTPCSigmaPr<3 || negNTPCSigmaEl<3);
    //
    // positive leg
    UInt_t cutBit0 = 0, cutBit1=0, fTrackCutBits=0;
    cutBit0 = SetCutBitsAndSomeTrackVariables(trackPosTest,0);
    Float_t dEdx0    = trackPosTest->GetTPCsignal();
    Float_t itsdEdx0 = trackPosTest->GetITSsignal();
    Float_t ptot0  = trackPosTest->GetInnerParam()->GetP();
    Float_t p0     = trackPosTest->P();
    Float_t pT0    = trackPosTest->Pt();
    Float_t eta0   = trackPosTest->Eta();
    Int_t sign0    = trackPosTest->Charge();
    Float_t phi0   = trackPosTest->Phi();
    Float_t nSigmasPiTOF0 = fPIDResponse->NumberOfSigmasTOF(trackPosTest, AliPID::kPion,  fPIDResponse->GetTOFResponse().GetTimeZero());
    Float_t nSigmasPrTOF0 = fPIDResponse->NumberOfSigmasTOF(trackPosTest, AliPID::kProton,fPIDResponse->GetTOFResponse().GetTimeZero());
    //
    // negative leg
    cutBit1 = SetCutBitsAndSomeTrackVariables(trackNegTest,0);
    Float_t dEdx1  = trackNegTest->GetTPCsignal();
    Float_t itsdEdx1 = trackNegTest->GetITSsignal();
    Float_t ptot1  = trackNegTest->GetInnerParam()->GetP();
    Float_t p1     = trackNegTest->P();
    Float_t pT1    = trackNegTest->Pt();
    Float_t eta1   = trackNegTest->Eta();
    Int_t sign1    = trackNegTest->Charge();
    Float_t phi1   = trackNegTest->Phi();
    Float_t nSigmasPiTOF1 = fPIDResponse->NumberOfSigmasTOF(trackNegTest, AliPID::kPion,  fPIDResponse->GetTOFResponse().GetTimeZero());
    Float_t nSigmasPrTOF1 = fPIDResponse->NumberOfSigmasTOF(trackNegTest, AliPID::kProton,fPIDResponse->GetTOFResponse().GetTimeZero());
    //
    // --------------------------------------------------------------
    //  Fill Clean Samples tree
    // --------------------------------------------------------------
    //
    if (selectNegLegs || selectPosLegs)
    {
      if (fFillArmPodTree)
      {
        if(!fTreeSRedirector) return;
        (*fTreeSRedirector)<<"fArmPodTree"<<
        "gid="                  << fEventGID             <<  //  global event ID
        "eventtime="            << fTimeStamp            <<
        "intrate="              << fIntRate              <<  // interaction rate
        "piFromK0="             << fCleanPionsFromK0    <<  // K0s cut for pions
        "v0haspixel="           << fHasV0FirstITSlayer  <<  // ITS pixel cout
        "purity="               << v0purity             <<
        "qt="                   << fQt                  <<  // qT
        "alfa="                 << fAlfa                <<  // alpha
        "cent="                 << fCentrality          <<  // centrality
        //
        "cutBit0="               << cutBit0        <<  // cut bits
        "itsdEdx0="              << itsdEdx0       <<  // TPC dEdx
        "dEdx0="                 << dEdx0          <<  // TPC dEdx
        "sign0="                 << sign0          <<
        "ptot0="                 << ptot0          <<  // momentum
        "p0="                    << p0             <<
        "pT0="                   << pT0            <<
        "eta0="                  << eta0           <<  // eta
        "phi0="                  << phi0           <<  // eta
        "nSigmasPiTOF0="         << nSigmasPiTOF0  <<  // TOF nsigma cut for pions
        "nSigmasPrTOF0="         << nSigmasPrTOF0  <<  // TOF nsigma cut for protons
        //
        "cutBit1="               << cutBit1        <<  // cut bits
        "dEdx1="                 << dEdx1          <<  // TPC dEdx
        "itsdEdx1="              << itsdEdx1       <<  // TPC dEdx
        "sign1="                 << sign1          <<
        "ptot1="                 << ptot1          <<  // momentum
        "p1="                    << p1             <<
        "pT1="                   << pT1            <<
        "eta1="                  << eta1           <<  // eta
        "phi1="                  << phi1           <<  // eta
        "nSigmasPiTOF1="         << nSigmasPiTOF1  <<  // TOF nsigma cut for pions
        "nSigmasPrTOF1="         << nSigmasPrTOF1  <<  // TOF nsigma cut for protons
        //
        "\n";

      }
    }
    cutBit0 = 0, cutBit1=0; fTrackCutBits=0;

  } // end of V0 loop

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::GetExpecteds(AliESDtrack *track, Float_t closestPar[3])
{

  //
  // bettaGamma is not well deifned below bg=0.01 --> below 200MeV protons and deuterons
  Float_t ptotForBetaGamma = track->GetInnerParam()->GetP();
  Float_t ptotForBetaGammaThr = 0.2;
  // if (ptotForBetaGamma<ptotForBetaGammaThr) return;
  //
  // --------------------------------------------------------------
  //  Calculates expected sigma and dEdx for a given track and returns colesest expected particle and its index
  // --------------------------------------------------------------
  //
  fNSigmasElTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
  fNSigmasPiTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  fNSigmasKaTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
  //
  fNSigmasElTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron, fPIDResponse->GetTOFResponse().GetTimeZero());
  fNSigmasPiTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion,     fPIDResponse->GetTOFResponse().GetTimeZero());
  fNSigmasKaTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon,     fPIDResponse->GetTOFResponse().GetTimeZero());
  //
  if (ptotForBetaGamma>ptotForBetaGammaThr){
    fNSigmasPrTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    fNSigmasDeTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
    fNSigmasPrTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton,   fPIDResponse->GetTOFResponse().GetTimeZero());
    fNSigmasDeTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron, fPIDResponse->GetTOFResponse().GetTimeZero());
  }

  //
  //
  Int_t nSigmaTmp = (fEventInfo) ? 10000 : 2;
  //
  // Electron Expected mean and sigma within 2nsigmaTPC
  if (TMath::Abs(fNSigmasElTPC)<nSigmaTmp && ptotForBetaGamma>ptotForBetaGammaThr) {
    // fDEdxEl  = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kElectron, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    // fSigmaEl = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kElectron, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    fDEdxEl  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kElectron);
    fSigmaEl = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kElectron);
  }
  //
  // Pion Expected mean and sigma within 2nsigmaTPC
  if (TMath::Abs(fNSigmasPiTPC)<nSigmaTmp && ptotForBetaGamma>ptotForBetaGammaThr) {
    // fDEdxPi  = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kPion,     AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    // fSigmaPi = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kPion,     AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    fDEdxPi  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kPion);
    fSigmaPi = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kPion);
  }
  //
  // Kaon Expected mean and sigma within 2nsigmaTPC
  if (TMath::Abs(fNSigmasKaTPC)<nSigmaTmp && ptotForBetaGamma>ptotForBetaGammaThr) {
    // fDEdxKa  = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kKaon,  AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    // fSigmaKa = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kKaon,   AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    fDEdxKa  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kKaon);
    fSigmaKa = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kKaon);
  }
  //
  // Proton Expected mean and sigma within 2nsigmaTPC
  if (TMath::Abs(fNSigmasPrTPC)<nSigmaTmp && ptotForBetaGamma>ptotForBetaGammaThr) {
    // fDEdxPr  = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kProton, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    // fSigmaPr = fPIDResponse->GetTPCResponse().GetExpectedSigma(track,  AliPID::kProton, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    fDEdxPr  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kProton);
    fSigmaPr = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kProton);
  }
  //
  // Deuteron Expected mean and sigma within 2nsigmaTPC
  if (TMath::Abs(fNSigmasDeTPC)<nSigmaTmp && ptotForBetaGamma>ptotForBetaGammaThr) {
    // fDEdxDe  = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kDeuteron, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    // fSigmaDe = fPIDResponse->GetTPCResponse().GetExpectedSigma(track,  AliPID::kDeuteron, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    fDEdxDe  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kDeuteron);
    fSigmaDe = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kDeuteron);
  }
  //
  // --------------------------------------------------------------
  //  Find Closest dEdx its corresponding particle and mass
  // --------------------------------------------------------------
  //
  Float_t values[] = {fDEdxEl, fDEdxPi, fDEdxKa, fDEdxPr, fDEdxDe};
  Float_t tpcdEdx = track->GetTPCsignal();
  Float_t smallestDiff = TMath::Abs(tpcdEdx - values[0]);
  Int_t closestIndex = 0;
  for (Int_t i = 0; i < 5; i++) {
    Double_t currentDiff = TMath::Abs(tpcdEdx - values[i]);
    if (currentDiff < smallestDiff) {
      smallestDiff = currentDiff;
      closestIndex = i;
    }
  }
  //
  // TF1 f1("f1","AliExternalTrackParam::BetheBlochAleph(x/0.1)",0.1,5)
  Float_t partMass = 0.;
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  if (closestIndex == 0 ) partMass = pdg->GetParticle(kPDGel)->Mass();   // GetParticle("e+")
  if (closestIndex == 1 ) partMass = pdg->GetParticle(kPDGpi)->Mass();  // GetParticle("pi+")
  if (closestIndex == 2 ) partMass = pdg->GetParticle(kPDGka)->Mass();  // GetParticle("K+")
  if (closestIndex == 3 ) partMass = pdg->GetParticle(kPDGpr)->Mass(); // GetParticle("proton")
  if (closestIndex == 4 ) partMass = 2.01410178;                     // pdg->GetParticle(1000010020)->Mass();
  //
  closestPar[0]=values[closestIndex];
  closestPar[1]=closestIndex;
  closestPar[2]=partMass;

}
//________________________________________________________________________
Bool_t AliAnalysisTaskEbyeIterPID::CheckIfFromResonance(Int_t mcType, AliMCParticle *trackMCgen, Int_t trackIndex, Bool_t parInterest, Double_t ptot, Double_t eta, Double_t cent, Bool_t fillTree)
{

  //
  // default is accept resonances
  Bool_t acceptRes = kTRUE;
  //
  TObjString momName="xxx";
  Int_t labMom = trackMCgen->Particle()->GetFirstMother();
  Int_t pdgMom = 0;
  if ( (labMom>=0) && (labMom < fMCEvent->GetNumberOfTracks()) ){
    pdgMom  = fMCStack->Particle(labMom)->GetPdgCode();
    momName = fMCStack->Particle(labMom)->GetName();
  }
  //
  //Check if the particle is in the black list of resonances
  for (Int_t ires=0;ires<fNResBins;ires++){

    if ( fResonances[ires]=="xxx" ){
      if (!(momName.GetString().Contains(fResonances[ires]))) {acceptRes=kFALSE; break;}
    } else {
      if (momName.GetString().Contains(fResonances[ires])) {acceptRes=kFALSE; break;}
    }
  }
  //
  // dump resonance info
  Int_t pdg = trackMCgen->Particle()->GetPdgCode();
  TObjString parName(trackMCgen->Particle()->GetName());
  if(fEventCountInFile==2 && !fRunOnGrid && fillTree) {
    if(!fTreeSRedirector) return kFALSE;
    (*fTreeSRedirector)<<"resonance"<<
    "acceptRes="   << acceptRes   <<
    "mcType="      << mcType       <<        // lower edge of momentum bin
    "ptot="        << ptot       <<          // lower edge of momentum bin
    "eta="         << eta     <<             // lower edge of eta bin
    "cent="        << cent        <<         // cent bin
    "parInterest=" << parInterest <<         // only pi, ka, and proton
    "pdg="         << pdg         <<         // pdg of prim particle
    "lab="         << trackIndex  <<         // index of prim particle
    "pdgMom="      << pdgMom      <<         // pdg of mother
    "labMom="      << labMom      <<         // index of mother
    "parName.="    << &parName    <<         //  full path - file name with ESD
    "momName.="    << &momName    <<         //  full path - file name with ESD
    "\n";
  }

  return acceptRes;

}
//________________________________________________________________________
Bool_t AliAnalysisTaskEbyeIterPID::CheckIfFromAnyResonance(AliMCParticle *trackMCgen, Float_t etaLow, Float_t etaUp, Float_t pDown, Float_t pUp)
{

  Bool_t sisterInAcceptance = kTRUE;
  Bool_t motherInAcceptance = kTRUE;

  Int_t labMom = trackMCgen->GetMother();
  Int_t pdgMom = 0;
  if ((labMom>=0) && (labMom < fMCEvent->GetNumberOfTracks())){
    //
    // Check if the mother is also in the acceptance
    AliMCParticle *momTrack = (AliMCParticle *)fMCEvent->GetTrack(labMom);
    pdgMom = momTrack->Particle()->GetPdgCode();
    if ( pdgMom!=0 ) {
      Int_t    motherSign = momTrack->Charge();
      Double_t motherEta  = (fRapidityType==0) ? momTrack->Eta() :  momTrack->Y();
      Double_t motherMom  = momTrack->P();
      Bool_t etaAcc  = (motherEta>=etaLow && motherEta<etaUp);
      Bool_t momAcc  = (motherMom>=pDown  && motherMom<pUp  );
      if ( !(etaAcc && momAcc) ) motherInAcceptance=kFALSE;
    }
    //
    // Check if the sister is also in the acceptance
    Int_t labSister = momTrack->GetDaughterLast();
    if ((labSister>=0) && (labSister < fMCEvent->GetNumberOfTracks())){
      AliMCParticle *sisterTrack = (AliMCParticle *)fMCEvent->GetTrack(labSister);
      Int_t    sisterSign = sisterTrack->Charge();
      Double_t sisterEta  = (fRapidityType==0) ? sisterTrack->Eta() :  sisterTrack->Y();
      Double_t sisterMom  = sisterTrack->P();
      Bool_t etaAcc  = (sisterEta>=etaLow && sisterEta<etaUp);
      Bool_t momAcc  = (sisterMom>=pDown  && sisterMom<pUp  );
      if ( !(sisterSign!=0 && etaAcc && momAcc) ) sisterInAcceptance=kFALSE;
    }
  }
  // default is accept resonances
  Bool_t acceptRes = kTRUE;
  if ( pdgMom!=0 &&  fSisterCheck==0 ) acceptRes = kFALSE;                          // in anycase if mother exist      reject particle
  if ( pdgMom!=0 &&  motherInAcceptance &&  sisterInAcceptance && fSisterCheck==1 ) acceptRes = kFALSE;   // if sister and mother are in acc reject particle
  if ( pdgMom!=0 &&  motherInAcceptance && !sisterInAcceptance && fSisterCheck==2 ) acceptRes = kFALSE;   // if sister and mother are in acc reject particle
  if ( pdgMom!=0 &&  motherInAcceptance && !sisterInAcceptance && fSisterCheck==3 ) acceptRes = kTRUE;    // if sister and mother are in acc accept particle
  //
  if ( pdgMom!=0 && !motherInAcceptance && !sisterInAcceptance && fSisterCheck==4 ) acceptRes = kFALSE;    // if sister and mother are in acc accept particle
  if ( pdgMom!=0 && !motherInAcceptance &&  sisterInAcceptance && fSisterCheck==5 ) acceptRes = kFALSE;    // if sister and mother are in acc accept particle
  if ( pdgMom!=0 && !motherInAcceptance &&  fSisterCheck==6 ) acceptRes = kFALSE;    // if sister and mother are in acc accept particle
  return acceptRes;

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::WeakAndMaterial()
{
  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the WeakAndMaterial ===== " << std::endl;
  //
  // ======================================================================
  // ======================================================================
  //
  // ========= Efficiency Check Eta momentum and Centrality scan ==========
  AliMCParticle *trackMCgen;
  for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++) {
    // track loop
    //
    // Select real trigger event and reject other pile up vertices
    Bool_t isTPCPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iTrack,fMCEvent);
    Bool_t isITSPileup = AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCEvent, "Hijing");
    if (isTPCPileup || isITSPileup) continue;
    //
    // initialize the dummy particle id
    Int_t fElWeak =-100., fPiWeak =-100., fKaWeak =-100., fPrWeak =-100., fDeWeak =-100., fMuWeak =-100.;
    trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);

    // apply primary vertex and eta cut
    if ( !fMCStack->IsPhysicalPrimary(iTrack) ) continue;

    Int_t pdg   = trackMCgen->Particle()->GetPdgCode();
    Int_t iPart = -10;
    if (TMath::Abs(pdg) == kPDGel) {iPart = 0; fElWeak = iPart;} // select el-
    if (TMath::Abs(pdg) == kPDGpi) {iPart = 1; fPiWeak = iPart;} // select pi+
    if (TMath::Abs(pdg) == kPDGka) {iPart = 2; fKaWeak = iPart;} // select ka+
    if (TMath::Abs(pdg) == kPDGpr) {iPart = 3; fPrWeak = iPart;} // select pr+
    if (TMath::Abs(pdg) == kPDGde) {iPart = 4; fDeWeak = iPart;} // select de
    if (TMath::Abs(pdg) == kPDGmu) {iPart = 5; fMuWeak = iPart;} // select mu-
    if (iPart == -10) continue;
    Float_t fSignWeak = (pdg<0) ? -1:1;
    Float_t fPtotWeak = trackMCgen->P();
    Float_t fPtWeak   = trackMCgen->Pt();
    Float_t fYWeak    = trackMCgen->Y();
    Float_t fEtaWeak  = trackMCgen->Eta();
    Float_t fPhiWeak  = trackMCgen->Phi();

    if(!fTreeSRedirector) return;
    (*fTreeSRedirector)<<"fullacc"<<
    "ptot="      << fPtotWeak <<         // mc momentum
    "pT="        << fPtWeak <<         // mc momentum
    "Y="         << fYWeak <<         // mc momentum
    "eta="       << fEtaWeak <<          // mc eta
    "phi="       << fPhiWeak <<          // mc eta
    "cent="      << fCentrality <<     // Centrality
    "centimp="   << fCentImpBin <<     // Centrality
    "sign="      << fSignWeak <<         // sign
    "el="        << fElWeak <<         // sign
    "pi="        << fPiWeak <<         // sign
    "ka="        << fKaWeak <<         // sign
    "pr="        << fPrWeak <<         // sign
    "de="        << fDeWeak <<         // sign
    "mu="        << fMuWeak <<         // sign
    "\n";

  } // ======= end of track loop =======

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FillDnchDeta()
{

  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillDnchDeta ===== " << std::endl;
  //
  // ======================================================================
  // ======================================================================
  //
  // ========= Efficiency Check Eta momentum and Centrality scan ==========
  const Int_t netabins = 20;
  Double_t etaDownArray[netabins] ={0.};
  Double_t etaUpArray[netabins]   ={0.};
  for (Int_t i=0; i<netabins; i++){
    etaUpArray[i]=0.1*(i+1);
    etaDownArray[i]=etaUpArray[i]*-1.;
  }
  Double_t centDownArray[9]={0., 5.,  10., 20., 30., 40., 50., 60., 70.};
  Double_t centUpArray[9]  ={5., 10., 20., 30., 40., 50., 60., 70., 80.};
  for (Int_t ieta=0; ieta<netabins; ieta++){
    for (Int_t icent=0; icent<9; icent++){

      AliMCParticle *trackMCgen;
      Int_t trCount=0,    elCount=0,    piCount=0,    kaCount=0,    prCount=0;
      for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++) {
        // track loop
        //
        // Select real trigger event and reject other pile up vertices
        Bool_t isTPCPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iTrack,fMCEvent);
        Bool_t isITSPileup = AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCEvent, "Hijing");
        if (isTPCPileup || isITSPileup) continue;
        //
        // initialize the dummy particle id
        trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);

        // apply primary vertex and eta cut
        if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
        // 	if ((trackMCgen->Eta()<fEtaDown) || (trackMCgen->Eta()>fEtaUp)) continue;
        Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
        Double_t etaGen  = trackMCgen->Eta();
        // skip neutral particles
        if ( TMath::Abs(trackMCgen->Charge()) < 0.0001 ) continue;

        Int_t iPart = -10;
        if (TMath::Abs(pdg) == kPDGel) iPart = 0;  // select el-
        if (TMath::Abs(pdg) == kPDGpi) iPart = 1;  // select pi+
        if (TMath::Abs(pdg) == kPDGka) iPart = 2;  // select ka+
        if (TMath::Abs(pdg) == kPDGpr) iPart = 3;  // select pr+

        // count first moments
        if ((fCentrality>=centDownArray[icent])
        &&(fCentrality<centUpArray[icent])
        &&(etaGen>=etaDownArray[ieta])
        &&(etaGen<=etaUpArray[ieta]))
        {
          trCount++;
          if ( iPart==0   ) elCount++;
          if ( iPart==1   ) piCount++;
          if ( iPart==2   ) kaCount++;
          if ( iPart==3   ) prCount++;
        }

      } // ======= end of track loop =======

      // fill tree which contains moments
      Double_t etaBin  = (TMath::Abs(etaDownArray[ieta])+etaUpArray[ieta]);
      Double_t centBin = (centDownArray[icent]+centUpArray[icent])/2.;
      if(!fTreeSRedirector) return;
      if ( trCount>0 ){
        (*fTreeSRedirector)<<"dnchdeta"<<
        "event="    << fEventCountInFile <<
        "centbin="  << centBin <<                 // cent bin
        "etabin="   << etaBin <<                  // eta bin
        "imppar="   << fMCImpactParameter <<      // impact parameter taken from MC event header
        "cent="     << fCentrality <<             // impact parameter taken from MC event header
        "trcount="  << trCount <<                 // number of identified tracks within the given cent and mom range
        "el="       << elCount <<                 // first moment of pions
        "pi="       << piCount <<                 // first moment of pions
        "ka="       << kaCount <<                 // first moment of kaons
        "pr="       << prCount <<                 // first moment of protons
        "\n";
      } // tree filling

    } // ======= end of Centrality loop =======
  } // ======= end of eta loop =======
  //
  // ======================================================================
  //
  // ======================================================================

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::SelectCleanSamplesFromV0s(AliESDv0 *v0, AliESDtrack *track0, AliESDtrack *track1)
{
  //
  // SetAliases and Metadata for the V0 trees
  //
  AliKFParticle kfparticle; //
  AliAnalysisTaskFilteredTree filteredV0;
  // Int_t type=filteredV0.GetKFParticle(v0,fESD,kfparticle);
  filteredV0.GetKFParticle(v0,fESD,kfparticle);
  //
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  Double_t massLambda = pdg->GetParticle("Lambda0")->Mass();
  Double_t massK0 = pdg->GetParticle("K0")->Mass();
  // const Double_t massProton  =pdg->GetParticle(kPDGpr)->Mass();
  // const Double_t massPion    =pdg->GetParticle(kPDGpi)->Mass();

  const Double_t livetimeK0=2.684341668932;  // livetime in cm (surpisely missing info in PDG - see root forum)
  const Double_t livetimeLambda=7.8875395;  // livetime in cm (missing info in PDG - see root forum)
  fHasTrack0FirstITSlayer = track0->HasPointOnITSLayer(0);
  fHasTrack1FirstITSlayer = track1->HasPointOnITSLayer(0);
  fHasV0FirstITSlayer = (fHasTrack0FirstITSlayer||fHasTrack1FirstITSlayer);
  //
  //
  Double_t v0Rr = v0->GetRr();  // rec position of the vertex CKBrev
  Double_t v0P  = v0->P();      // TMath::Sqrt(Px()*Px()+Py()*Py()+Pz()*Pz())

  //   tree->SetAlias("livetimeLikeK0",TString::Format("exp(-v0.fRr/(sqrt((v0.P()/%f)^2+1)*%f))",massK0, livetimeK0));
  //   tree->SetAlias("livetimeLikeLambda",TString::Format("exp(-v0.fRr/(sqrt((v0.P()/%f)^2+1)*%f))",massLambda,livetimeLambda));
  //   tree->SetAlias("livetimeLikeGamma","v0.fRr/80");
  //   tree->SetAlias("livetimeLikeBkg","v0.fRr/80");
  Double_t livetimeLikeK0 = TMath::Exp(-v0Rr/(TMath::Sqrt((v0P/massK0)*(v0P/massK0)+1)*livetimeK0));
  Double_t livetimeLikeLambda = TMath::Exp(-v0Rr/(TMath::Sqrt((v0P/massLambda)*(v0P/massLambda)+1)*livetimeLambda));
  Double_t livetimeLikeGamma = v0Rr/80.;
  // Double_t livetimeLikeBkg   = v0Rr/80.;

  // delta of mass
  Double_t K0Delta = v0->GetEffMass(2,2)-massK0;        //   tree->SetAlias("K0Delta","(v0.GetEffMass(2,2)-massK0)");
  Double_t LDelta  = v0->GetEffMass(4,2)-massLambda;    //   tree->SetAlias("LDelta","(v0.GetEffMass(4,2)-massLambda)");
  Double_t ALDelta = v0->GetEffMass(2,4)-massLambda;    //   tree->SetAlias("ALDelta","(v0.GetEffMass(2,4)-massLambda)");
  Double_t EDelta  = v0->GetEffMass(0,0);               //   tree->SetAlias("EDelta","(v0.GetEffMass(0,0))");

  // pull of the mass
  if (v0->GetKFInfo(2,2,1)==0. || v0->GetKFInfo(4,2,1)==0. || v0->GetKFInfo(2,4,1)==0. || v0->GetKFInfo(0,0,1)==0.) return;
  Double_t K0Pull = (v0->GetEffMass(2,2)-massK0)/v0->GetKFInfo(2,2,1);        //   tree->SetAlias("K0Pull","(v0.GetEffMass(2,2)-massK0)/v0.GetKFInfo(2,2,1)");
  Double_t LPull  = (v0->GetEffMass(4,2)-massLambda)/v0->GetKFInfo(4,2,1);    //   tree->SetAlias("LPull","(v0.GetEffMass(4,2)-massLambda)/v0.GetKFInfo(4,2,1)");
  Double_t ALPull = (v0->GetEffMass(2,4)-massLambda)/v0->GetKFInfo(2,4,1);    //   tree->SetAlias("ALPull","(v0.GetEffMass(2,4)-massLambda)/v0.GetKFInfo(2,4,1)");
  Double_t EPull  = EDelta/v0->GetKFInfo(0,0,1);                              //   tree->SetAlias("EPull","EDelta/v0.GetKFInfo(0,0,1)");

  // effective pull of the mass - (empirical values from fits)
  //   tree->SetAlias("K0PullEff","K0Delta/sqrt((3.63321e-03)**2+(5.68795e-04*v0.Pt())**2)");
  //   tree->SetAlias("LPullEff","LDelta/sqrt((1.5e-03)**2+(1.8e-04*v0.Pt())**2)");
  //   tree->SetAlias("ALPullEff","ALDelta/sqrt((1.5e-03)**2+(1.8e-04*v0.Pt())**2)");
  //   tree->SetAlias("EPullEff","v0.GetEffMass(0,0)/sqrt((5e-03)**2+(1.e-04*v0.Pt())**2)");
  // Double_t K0PullEff = K0Delta/TMath::Sqrt((3.63321e-03)*(3.63321e-03)+(5.68795e-04*v0->Pt())*(5.68795e-04*v0->Pt()));
  // Double_t LPullEff  = LDelta/TMath::Sqrt((1.5e-03)*(1.5e-03)+(1.8e-04*v0->Pt())*(1.8e-04*v0->Pt()));
  // Double_t ALPullEff = ALDelta/TMath::Sqrt((1.5e-03)*(1.5e-03)+(1.8e-04*v0->Pt())*(1.8e-04*v0->Pt()));
  // Double_t EPullEff  = v0->GetEffMass(0,0)/TMath::Sqrt((5e-03)*(5e-03)+(1.e-04*v0->Pt())*(1.e-04*v0->Pt()));

  //
  //    tree->SetAlias("dEdx0DProton","AliMathBase::BetheBlochAleph(track0.fIp.P()/massProton)");
  //    tree->SetAlias("dEdx1DProton","AliMathBase::BetheBlochAleph(track1.fIp.P()/massProton)");
  //    tree->SetAlias("dEdx0DPion","AliMathBase::BetheBlochAleph(track0.fIp.P()/massPion)");
  //    tree->SetAlias("dEdx1DPion","AliMathBase::BetheBlochAleph(track1.fIp.P()/massPion)");
  // Double_t dEdx0DProton = AliTPCParam::BetheBlochAleph(track0->GetInnerParam()->GetP()/massProton);
  // Double_t dEdx1DProton = AliTPCParam::BetheBlochAleph(track1->GetInnerParam()->GetP()/massProton);
  // Double_t dEdx0DPion   = AliTPCParam::BetheBlochAleph(track0->GetInnerParam()->GetP()/massPion);
  // Double_t dEdx1DPion   = AliTPCParam::BetheBlochAleph(track1->GetInnerParam()->GetP()/massPion);

  //   tree->SetAlias("K0Like0","exp(-K0Pull^2)*livetimeLikeK0");
  //   tree->SetAlias("LLike0","exp(-LPull^2)*livetimeLikeLambda");
  //   tree->SetAlias("ALLike0","exp(-ALPull^2)*livetimeLikeLambda");
  //   tree->SetAlias("ELike0","exp(-abs(EPull)*0.2)*livetimeLikeGamma");
  //   tree->SetAlias("V0Like","exp(-acos(v0.fPointAngle)*v0.fRr/0.36)*exp(-sqrt(kf.GetChi2())/0.5)");
  Double_t K0Like0 = TMath::Exp(-K0Pull*K0Pull)*livetimeLikeK0;
  Double_t LLike0  = TMath::Exp(-LPull*LPull)*livetimeLikeLambda;
  Double_t ALLike0 = TMath::Exp(-ALPull*ALPull)*livetimeLikeLambda;
  Double_t ELike0  = TMath::Exp(-abs(EPull)*0.2)*livetimeLikeGamma;
  Double_t V0Like  = TMath::Exp(-TMath::ACos(v0->GetV0CosineOfPointingAngle())*v0Rr/0.36)*TMath::Exp(-TMath::Sqrt(kfparticle.GetChi2())/0.5);


  //   tree->SetAlias("BkgLike","0.000005*ntracks");  // backround coeefecint  to be fitted - depends on other cuts
  Int_t ntracks = fESD->GetNumberOfTracks();
  Double_t BkgLike = 0.000005*ntracks;    // backround coeefecint  to be fitted - depends on other cuts

  //   tree->SetAlias("ELike","(V0Like*ELike0)/(V0Like*(K0Like0+LLike0+ALLike0+ELike0)+BkgLike)");
  // Double_t ELike = (V0Like*ELike0)/(V0Like*(K0Like0+LLike0+ALLike0+ELike0)+BkgLike);
  //   tree->SetAlias("K0Like","K0Like0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike)");
  // Double_t K0Like = K0Like0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike);
  //   tree->SetAlias("LLike","LLike0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike)");
  Double_t LLike = LLike0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike);
  //   tree->SetAlias("ALLike","ALLike0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike)");
  Double_t ALLike = ALLike0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike);

  Double_t tr0NTPCSigmaPi = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track0, AliPID::kPion));
  Double_t tr1NTPCSigmaPi = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track1, AliPID::kPion));

  Float_t ptotForBetaGamma0 = track0->GetInnerParam()->GetP();
  Float_t ptotForBetaGamma1 = track1->GetInnerParam()->GetP();
  Float_t ptotForBetaGammaThr = 0.2;
  Double_t tr0NTPCSigmaPr = (ptotForBetaGamma0>ptotForBetaGammaThr) ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track0, AliPID::kProton)) : 0.;
  Double_t tr1NTPCSigmaPr = (ptotForBetaGamma1>ptotForBetaGammaThr) ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track1, AliPID::kProton)) : 0.;

  //   treeV0->SetAlias("cleanPion0FromK0","K0Like0>0.05&&K0Like0>(LLike0+ALLike0+ELike0)*3&&abs(K0Delta)<0.006&&V0Like>0.1&&abs(track1.fTPCsignal/dEdx1DPion-50)<8&&v0.PtArmV0()>0.06");
  fCleanPion0FromK0 = (K0Like0>0.05) && (K0Like0>(LLike0+ALLike0+ELike0)*3) && (TMath::Abs(K0Delta)<0.006) && (V0Like>0.1) && tr1NTPCSigmaPi<2 && (v0->PtArmV0()>0.06);
  //   treeV0->SetAlias("cleanPion1FromK0","K0Like0>0.05&&K0Like0>(LLike0+ALLike0+ELike0)*3&&abs(K0Delta)<0.006&&V0Like>0.1&&abs(track0.fTPCsignal/dEdx0DPion-50)<8&&v0.PtArmV0()>0.06");
  fCleanPion1FromK0 = (K0Like0>0.05) && (K0Like0>(LLike0+ALLike0+ELike0)*3) && (TMath::Abs(K0Delta)<0.006) && (V0Like>0.1) && tr0NTPCSigmaPi<2 && (v0->PtArmV0()>0.06);
  //   treeV0->SetAlias("cleanPion0FromLambda","ALLike>0.05&&ALLike0>(K0Like0+LLike0+ELike0)*3&&abs(ALDelta)<0.006&&V0Like>0.1&&abs(track1.fTPCsignal/dEdx1DProton-50)<8");
  fCleanPion0FromLambda = (ALLike>0.05) && (ALLike0>(K0Like0+LLike0+ELike0)*3) && (TMath::Abs(ALDelta)<0.006) && (V0Like>0.1) && tr1NTPCSigmaPr<2;
  //   treeV0->SetAlias("cleanPion1FromLambda","LLike>0.05&&LLike0>(K0Like0+ALLike0+ELike0)*3&&abs(LDelta)<0.006&&V0Like>0.1&&abs(track0.fTPCsignal/dEdx0DProton-50)<8");
  fCleanPion1FromLambda = (LLike>0.05) && (LLike0>(K0Like0+ALLike0+ELike0)*3) && (TMath::Abs(LDelta)<0.006) && (V0Like>0.1) && tr0NTPCSigmaPr<2;
  //   treeV0->SetAlias("cleanProton0FromLambda","LLike>0.05&&LLike0>(K0Like0+ALLike0+ELike0)*3&&abs(LDelta)<0.001&&V0Like>0.1&&abs(track1.fTPCsignal/dEdx1DPion-50)<8");
  fCleanProton0FromLambda = (LLike>0.05) && (LLike0>(K0Like0+ALLike0+ELike0)*3) && (TMath::Abs(LDelta)<0.006) && (V0Like>0.1) && tr1NTPCSigmaPi<2;
  //   treeV0->SetAlias("cleanProton1FromLambda","ALLike>0.05&&ALLike0>(K0Like0+LLike0+ELike0)*3&&abs(ALDelta)<0.001&&V0Like>0.1&&abs(track0.fTPCsignal/dEdx0DPion-50)<8");
  fCleanProton1FromLambda = (ALLike>0.05) && (ALLike0>(K0Like0+LLike0+ELike0)*3) && (TMath::Abs(ALDelta)<0.006) && (V0Like>0.1) && tr0NTPCSigmaPi<2;

  fCleanPionsFromK0 =  (fCleanPion0FromK0 || fCleanPion1FromK0);

  //    //   treeV0->SetAlias("cleanPion0FromK0","K0Like0>0.05&&K0Like0>(LLike0+ALLike0+ELike0)*3&&abs(K0Delta)<0.006&&V0Like>0.1&&abs(track1.fTPCsignal/dEdx1DPion-50)<8&&v0.PtArmV0()>0.06");
  //   fCleanPion0FromK0 = (K0Like0>0.05) && (K0Like0>(LLike0+ALLike0+ELike0)*3) && (TMath::Abs(K0Delta)<0.006) && (V0Like>0.1) && (TMath::Abs(track1->GetTPCsignal()/dEdx1DPion-50)<8) && (v0->PtArmV0()>0.06);
  //   //   treeV0->SetAlias("cleanPion1FromK0","K0Like0>0.05&&K0Like0>(LLike0+ALLike0+ELike0)*3&&abs(K0Delta)<0.006&&V0Like>0.1&&abs(track0.fTPCsignal/dEdx0DPion-50)<8&&v0.PtArmV0()>0.06");
  //   fCleanPion1FromK0 = (K0Like0>0.05) && (K0Like0>(LLike0+ALLike0+ELike0)*3) && (TMath::Abs(K0Delta)<0.006) && (V0Like>0.1) && (TMath::Abs(track0->GetTPCsignal()/dEdx0DPion-50)<8) && (v0->PtArmV0()>0.06);
  //   //   treeV0->SetAlias("cleanPion0FromLambda","ALLike>0.05&&ALLike0>(K0Like0+LLike0+ELike0)*3&&abs(ALDelta)<0.006&&V0Like>0.1&&abs(track1.fTPCsignal/dEdx1DProton-50)<8");
  //   fCleanPion0FromLambda = (ALLike>0.05) && (ALLike0>(K0Like0+LLike0+ELike0)*3) && (TMath::Abs(ALDelta)<0.006) && (V0Like>0.1) && TMath::Abs(track1->GetTPCsignal()/dEdx1DProton-50)<8;
  //   //   treeV0->SetAlias("cleanPion1FromLambda","LLike>0.05&&LLike0>(K0Like0+ALLike0+ELike0)*3&&abs(LDelta)<0.006&&V0Like>0.1&&abs(track0.fTPCsignal/dEdx0DProton-50)<8");
  //   fCleanPion1FromLambda = (LLike>0.05) && (LLike0>(K0Like0+ALLike0+ELike0)*3) && (TMath::Abs(LDelta)<0.006) && (V0Like>0.1) && TMath::Abs(track0->GetTPCsignal()/dEdx0DProton-50)<8;
  //   //   treeV0->SetAlias("cleanProton0FromLambda","LLike>0.05&&LLike0>(K0Like0+ALLike0+ELike0)*3&&abs(LDelta)<0.001&&V0Like>0.1&&abs(track1.fTPCsignal/dEdx1DPion-50)<8");
  //   fCleanProton0FromLambda = (LLike>0.05) && (LLike0>(K0Like0+ALLike0+ELike0)*3) && (TMath::Abs(LDelta)<0.006) && (V0Like>0.1) && TMath::Abs(track1->GetTPCsignal()/dEdx1DPion-50)<8;
  //   //   treeV0->SetAlias("cleanProton1FromLambda","ALLike>0.05&&ALLike0>(K0Like0+LLike0+ELike0)*3&&abs(ALDelta)<0.001&&V0Like>0.1&&abs(track0.fTPCsignal/dEdx0DPion-50)<8");
  //   fCleanProton1FromLambda = (ALLike>0.05) && (ALLike0>(K0Like0+LLike0+ELike0)*3) && (TMath::Abs(ALDelta)<0.006) && (V0Like>0.1) && TMath::Abs(track0->GetTPCsignal()/dEdx0DPion-50)<8;
  //


  //   // V0 - cuts - for PID
  //   tree->SetAlias("cutDist","sqrt((track0.fIp.fP[0]-track1.fIp.fP[0])**2+(track0.fIp.fP[1]-track1.fIp.fP[1])**2)>3");
  //   tree->SetAlias("cutLong","track0.GetTPCClusterInfo(3,1,0)-5*abs(track0.fP[4])>130&&track1.GetTPCClusterInfo(3,1,0)>130-5*abs(track0.fP[4])");
  //   tree->SetAlias("cutPID","track0.fTPCsignal>0&&track1.fTPCsignal>0");
  //   tree->SetAlias("cutResol","sqrt(track0.fC[14]/track0.fP[4])<0.15&&sqrt(track1.fC[14]/track1.fP[4])<0.15");
  //   tree->SetAlias("cutV0","cutPID&&cutLong&&cutResol");
  //   //
  //   tree->SetAlias("K0PullBkg","min(min(abs(LPull),abs(ALPull)),abs(EPull))+0");
  //   tree->SetAlias("LambdaPullBkg","min(min(abs(K0Pull),abs(ALPull)),abs(EPull)+0)");
  //   tree->SetAlias("ALambdaPullBkg","min(min(abs(K0Pull),abs(LPull)),abs(EPull)+0)");
  //   tree->SetAlias("EPullBkg","min(min(abs(K0Pull),abs(LPull)),abs(ALPull)+0)");
  //   //
  //   tree->SetAlias("K0Selected",      "abs(K0Pull)<3. &&abs(K0PullEff)<3.  && abs(LPull)>3  && abs(ALPull)>3 &&v0.PtArmV0()>0.11");
  //   tree->SetAlias("LambdaSelected",  "abs(LPull)<3.  &&abs(LPullEff)<3.   && abs(K0Pull)>3 && abs(EPull)>3  && abs(EDelta)>0.05");
  //   tree->SetAlias("ALambdaSelected", "abs(ALPull)<3. &&abs(ALPullEff)<3   && abs(K0Pull)>3 && abs(EPull)>3  &&abs(EDelta)>0.05");
  //   tree->SetAlias("GammaSelected", "abs(EPull)<3     && abs(K0Pull)>3 && abs(LPull)>3 && abs(ALPull)>3");
  //   tree->SetAlias("BkgLike","0.000005*ntracks");  // backround coeefecint  to be fitted - depends on other cuts
  //   //
  //   tree->SetAlias("ELike","(V0Like*ELike0)/(V0Like*(K0Like0+LLike0+ALLike0+ELike0)+BkgLike)");
  //   tree->SetAlias("K0Like","K0Like0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike)");
  //   tree->SetAlias("LLike","LLike0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike)");
  //   tree->SetAlias("ALLike","ALLike0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike)");
  //   //
  //   tree->SetAlias("K0PIDPull","(abs(track0.fTPCsignal/dEdx0DPion-50)+abs(track1.fTPCsignal/dEdx1DPion-50))/5.");
  //   tree->SetAlias("mpt","1/v0.Pt()");                 //
  //   tree->SetAlias("tglV0","v0.Pz()/v0.Pt()");                 //
  //   tree->SetAlias("alphaV0","atan2(v0.Py(),v0.Px()+0)");
  //   tree->SetAlias("dalphaV0","alphaV0-((int(36+9*(alphaV0/pi))-36)*pi/9.)");

}
//________________________________________________________________________
Bool_t AliAnalysisTaskEbyeIterPID::ApplyDCAcutIfNoITSPixel(AliESDtrack *track)
{

  //     treeITS->SetAlias("ITS01","(Tracks[].HasPointOnITSLayer(0)||Tracks[].HasPointOnITSLayer(1))");
  //     treeITS->SetAlias("isPrimPtDep","abs(Tracks[].fD)<0.0182+0.0350/(Tracks[].Pt()^1.01)");
  //     treeITS->SetAlias("isPrimPtDep2","abs(Tracks[].fD/2)<0.0182+0.0350/(Tracks[].Pt()^1.01)");
  //     treeITS->SetAlias("isPrim2","sqrt(Tracks[].fD**2/Tracks[].fCdd+Tracks[].fZ**2/Tracks[].fCzz+0)<2");
  //     treeITS->SetAlias("isPrim5","sqrt(Tracks[].fD**2/Tracks[].fCdd+Tracks[].fZ**2/Tracks[].fCzz+0)<5");
  //     treeITS->SetAlias("IsPrimCA","((isPrim2&&Tracks[].fITSncls>2)||(isPrim5&&ITS01))");

  Float_t p[2],cov[3];
  track->GetImpactParameters(p,cov); // p[0]=fD; p[1]=fZ; cov[0]=fCdd; cov[1]=fCdz; cov[2]=fCzz;
  Bool_t isFirstITSlayer  = track->HasPointOnITSLayer(0);
  Bool_t isSecondITSlayer = track->HasPointOnITSLayer(1);

  fIsITSpixel01 = (isFirstITSlayer || isSecondITSlayer);
  fNITSclusters = track->GetNumberOfITSClusters();

  if (!cov[0] || !cov[2]) {
    return kFALSE;
  } else {
    fPrimRestriction = TMath::Sqrt((p[0]*p[0])/cov[0] + (p[1]*p[1])/cov[2]);
    return (fPrimRestriction<2 && fNITSclusters>2) || (fPrimRestriction<5 && fIsITSpixel01);
  }

  // t->SetAlias("less5","TMath::Sqrt((dcaxy*dcaxy)/cov0 + (dcaz*dcaz)/cov2)<5");
  // t->SetAlias("less2","TMath::Sqrt((dcaxy*dcaxy)/cov0 + (dcaz*dcaz)/cov2)<2");
  // t->SetAlias("lay01","its6");
  // t->SetAlias("pixCut","(less2 && itscls>2) || (less5 && lay01)");

}
//________________________________________________________________________
UInt_t AliAnalysisTaskEbyeIterPID::SetCutBitsAndSomeTrackVariables(AliESDtrack *track, Int_t particleType)
{
  //
  // Set some track variables
  //
  //
  // --------------------------------------------------------------
  //  calculate some variables by hand
  // --------------------------------------------------------------
  //
  // Double_t p[3];
  // track->GetPxPyPz(p);
  // Double_t momentum = TMath::Sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
  // Double_t pt       = TMath::Sqrt(p[0]*p[0] + p[1]*p[1]);
  // Double_t mass     = track->GetMass();  // assumed to be pion mass
  // Double_t energy   = TMath::Sqrt(mass*mass + momentum*momentum);
  // Float_t eta = -100.;
  // Float_t rap   = -100.;
  // if((momentum != TMath::Abs(p[2]))&&(momentum != 0)) eta = 0.5*TMath::Log((momentum + p[2])/(momentum - p[2]));
  // if((energy != TMath::Abs(p[2]))&&(energy != 0))     rap = 0.5*TMath::Log((energy + p[2])/(energy - p[2]));
  //
  // --------------------------------------------------------------
  //  some extra getters
  // --------------------------------------------------------------
  //
  Double_t goldenChi2 = 0.;
  fPVertex = track->P();
  fTheta=track->Theta();
  fSign= track->GetSign();
  fPx  =track->Px();
  fPy  =track->Py();
  fPz  =track->Pz();
  fPt  =track->Pt();
  fY   =track->Y();
  fPhi =track->Phi()-TMath::Pi();
  fEta = track->Eta();
  Float_t pv[2],cov[3];
  track->GetImpactParameters(pv,cov); // p[0]=fD; p[1]=fZ; cov[0]=fCdd; cov[1]=fCdz; cov[2]=fCzz;
  fTrackDCAxy = pv[0];
  fTrackDCAz  = pv[1];
  //
  // TPC related quantities
  Bool_t cleanDeTPC = kFALSE;
  if (track->GetInnerParam()){
    fPtot      = track->GetInnerParam()->GetP();
    fTPCSignal = track->GetTPCsignal();
    fTrackTPCSignalN     = track->GetTPCsignalN();
    fTrackTPCCrossedRows = Float_t(track->GetTPCCrossedRows());
    fTPCShared = track->GetTPCnclsS();
    fMissingCl = track->GetTPCClusterInfo(3,0,0,159);
    goldenChi2 = track->GetChi2TPCConstrainedVsGlobal(fVertex);
    fTrackLengthInActiveZone = track->GetLengthInActiveZone(1,3,230, track->GetBz(),0,0);
    //
    Float_t ptotForBetaGamma = track->GetInnerParam()->GetP();
    Float_t ptotForBetaGammaThr = 0.2;
    Double_t nSigmasDeTPC = (ptotForBetaGamma>ptotForBetaGammaThr) ? fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron) : 0.;
    cleanDeTPC = ((TMath::Abs(nSigmasDeTPC)<=2.));
    fNcl       = track->GetTPCncls();
    // fTrackChi2TPC  = (fNcl>0) ? TMath::Sqrt(TMath::Abs(track->GetTPCchi2()/fNcl)) : -1;  // ???
    fTrackChi2TPC  = (fNcl>0) ? TMath::Abs(track->GetTPCchi2()/fNcl) : -1;  // ???
    //
    // correct for the missing clusters
    fNclCorr = fNcl;
    fTrackChi2TPCcorr = fTrackChi2TPC;
    if (fCorrectForMissCl==2){
      if ( (fPtot>0.2 && fPtot<3.2) && TMath::Abs(fEta)<0.8 && fCentrality<5 && particleType<4 ){
        Float_t missCl_binContent = 0.;
        Int_t multBin = 0;
        if (fTPCMult>=6000  && fTPCMult<6500)  multBin=0;
        if (fTPCMult>=6500  && fTPCMult<7000)  multBin=1;
        if (fTPCMult>=7000  && fTPCMult<7500)  multBin=2;
        if (fTPCMult>=7500  && fTPCMult<8000)  multBin=3;
        if (fTPCMult>=8000  && fTPCMult<8500)  multBin=4;
        if (fTPCMult>=8500  && fTPCMult<9000)  multBin=5;
        if (fTPCMult>=9000  && fTPCMult<9500)  multBin=6;
        if (fTPCMult>=9500  && fTPCMult<10000) multBin=7;
        if (fTPCMult>=10000 && fTPCMult<10500) multBin=8;
        if (fTPCMult>=10500 && fTPCMult<11000) multBin=9;
        if (fTPCMult>=11000 && fTPCMult<11500) multBin=10;
        if (fTPCMult>=11500 && fTPCMult<12000) multBin=11;
        if (fTPCMult>=12000 && fTPCMult<12500) multBin=12;
        if (fTPCMult>=12500 && fTPCMult<13000) multBin=13;
        if (fTPCMult>=13000 && fTPCMult<13500) multBin=14;
        if (fTPCMult>=13500 && fTPCMult<14000) multBin=15;
        missCl_binContent = fH2MissCl[particleType][multBin].GetBinContent(fH2MissCl[particleType][multBin].FindBin(fPtot,fEta));
        if(missCl_binContent<30 && missCl_binContent>0) {
          fNclCorr = fNcl-(missCl_binContent/100.)*fNcl;
          fTrackChi2TPCcorr  = (fNclCorr>0) ? TMath::Abs(track->GetTPCchi2()/fNclCorr) : -1;  // ???
        }
      }
    }
    //
    // --------------------------------------------------------------
    //      Bayesian PID part
    // --------------------------------------------------------------
    //
    if (ptotForBetaGamma>ptotForBetaGammaThr) {
      fPIDCombined->SetDefaultTPCPriors();
      Double_t probTPC[AliPID::kSPECIES]={0.};
      Double_t probTOF[AliPID::kSPECIES]={0.};
      // Get TPC probabilities
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
      fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPC);
      fTrackProbPiTPC = probTPC[AliPID::kPion];
      fTrackProbKaTPC = probTPC[AliPID::kKaon];
      fTrackProbPrTPC = probTPC[AliPID::kProton];
      // fTrackProbDeTPC = probTPC[AliPID::kDeuteron];
      // Get TOF probabilities
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
      fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTOF);
      fTrackProbPiTOF = probTOF[AliPID::kPion];
      fTrackProbKaTOF = probTOF[AliPID::kKaon];
      fTrackProbPrTOF = probTOF[AliPID::kProton];
    }
  }
  //
  // its restrictions
  Bool_t isFirstITSlayer  = track->HasPointOnITSLayer(0);
  Bool_t isSecondITSlayer = track->HasPointOnITSLayer(1);
  fIsITSpixel01    = (isFirstITSlayer || isSecondITSlayer);
  fNITSclusters    = track->GetNumberOfITSClusters();
  if (cov[0]>0 && cov[1]>0){
    fPrimRestriction = TMath::Sqrt((pv[0]*pv[0])/cov[0] + (pv[1]*pv[1])/cov[2]);
  }
  //
  fTrackRequireITSRefit  = track->IsOn(AliESDtrack::kITSrefit); // track->IsOn(AliESDtrack::kTPCrefit);
  fTrackIsFirstITSlayer  = track->HasPointOnITSLayer(0);
  fTrackIsSecondITSlayer = track->HasPointOnITSLayer(1);
  fTrackNewITScut        = ApplyDCAcutIfNoITSPixel(track);
  //
  Double_t nclsTRD     = (Float_t)track->GetTRDncls();
  Double_t TOFSignalDx = track->GetTOFsignalDx();
  Double_t TOFSignalDz = track->GetTOFsignalDz();
  //
  //
  Float_t nSigmasKaTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon,    fPIDResponse->GetTOFResponse().GetTimeZero());
  Float_t nSigmasPrTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton,  fPIDResponse->GetTOFResponse().GetTimeZero());
  Float_t nSigmasDeTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron,fPIDResponse->GetTOFResponse().GetTimeZero());
  Bool_t cleanPrTOF = ((TMath::Abs(nSigmasPrTOF)<=2.5));
  Bool_t cleanDeTOF = ((TMath::Abs(nSigmasDeTOF)<=3.0));
  Bool_t cleanKaTOF = ((TMath::Abs(nSigmasKaTOF)<=2.5));
  Bool_t cleanKaTOFTRD = ((TMath::Abs(nSigmasKaTOF)<=1.2) && TOFSignalDz<1. && TOFSignalDx<1. && nclsTRD>100);
  //
  // Systematic settings
  fTrackCutBits=0;
  //
  if (fCorrectForMissCl==1){
    if (fNcl>=70) (fTrackCutBits |= 1 << kNCrossedRowsTPC60);
    if (fNcl>=80) (fTrackCutBits |= 1 << kNCrossedRowsTPC80);
    if (fNcl>=90) (fTrackCutBits |= 1 << kNCrossedRowsTPC100);
  } else if (fCorrectForMissCl==2){
    if (fNclCorr>=70) (fTrackCutBits |= 1 << kNCrossedRowsTPC60);
    if (fNclCorr>=80) (fTrackCutBits |= 1 << kNCrossedRowsTPC80);
    if (fNclCorr>=90) (fTrackCutBits |= 1 << kNCrossedRowsTPC100);
    // cout <<  "ncls  = " << fNclCorr <<  " --- " << fNcl << endl;
  } else {
    if (fTrackTPCCrossedRows>=60)  (fTrackCutBits |= 1 << kNCrossedRowsTPC60);
    if (fTrackTPCCrossedRows>=80)  (fTrackCutBits |= 1 << kNCrossedRowsTPC80);
    if (fTrackTPCCrossedRows>=100) (fTrackCutBits |= 1 << kNCrossedRowsTPC100);
  }
  //
  // Special treatment of the 2018 pass3 and 2015 pass2 data
  if ( (fYear==2015&&fPassIndex==2) || (fYear==2018&&fPassIndex==3) ){
    if (fTrackChi2TPC<2.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCSmall);
    if (fTrackChi2TPC<2.5) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPC);
    if (fTrackChi2TPC<3.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCLarge);
  } else {
    //
    // correction for missing clusters
    if (fCorrectForMissCl==2){
      if (fTrackChi2TPCcorr<3.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCSmall);   // ????
      if (fTrackChi2TPCcorr<4.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPC);        // ????
      if (fTrackChi2TPCcorr<5.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCLarge);   // ????
    } else {
      if (fTrackChi2TPC<3.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCSmall);   // ????
      if (fTrackChi2TPC<4.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPC);        // ????
      if (fTrackChi2TPC<5.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCLarge);   // ????
    }
  }

  Bool_t dca11h     = TMath::Abs(fTrackDCAxy)<0.0105+0.0350/TMath::Power(fPt,1.1);    // 10h tuned loose cut
  Bool_t dca10h     = TMath::Abs(fTrackDCAxy)<0.0182+0.0350/TMath::Power(fPt,1.01);    // 10h tuned loose cut
  Bool_t dcaBaseCut = TMath::Abs(fTrackDCAxy)<0.0208+0.04/TMath::Power(fPt,1.01);  // 10h tuned loose cut
  Bool_t dcaLoose   = TMath::Abs(fTrackDCAxy)<0.4;  // 10h tuned loose cut
  if (dca10h)     (fTrackCutBits |= 1 << kMaxDCAToVertexXYPtDepSmall);
  if (dcaBaseCut) (fTrackCutBits |= 1 << kMaxDCAToVertexXYPtDep);
  if (dcaLoose)   (fTrackCutBits |= 1 << kMaxDCAToVertexXYPtDepLarge);
  //
  if (TMath::Abs(fTrackDCAz)<2.0) (fTrackCutBits |= 1 << kVertexZSmall);
  if (TMath::Abs(fTrackDCAz)<3.0) (fTrackCutBits |= 1 << kVertexZ);
  if (TMath::Abs(fTrackDCAz)<4.0) (fTrackCutBits |= 1 << kVertexZLarge);
  //
  if (TMath::Abs(fVz)<6 && TMath::Abs(fVz)>0.2 ) (fTrackCutBits |= 1 << kEventVertexZSmall);
  if (TMath::Abs(fVz)<7 && TMath::Abs(fVz)>0.15) (fTrackCutBits |= 1 << kEventVertexZ);
  if (TMath::Abs(fVz)<8 && TMath::Abs(fVz)>0.1 ) (fTrackCutBits |= 1 << kEventVertexZLarge);
  if (TMath::Abs(fVz)<10 ) (fTrackCutBits |= 1 << kEventVertexZALICE);
  if (TMath::Abs(fVz)<7  ) (fTrackCutBits |= 1 << kEventVertexZALICETight);
  //
  if (fTrackRequireITSRefit)                           (fTrackCutBits |= 1 << kRequireITSRefit);
  if (fTrackIsFirstITSlayer || fTrackIsSecondITSlayer) (fTrackCutBits |= 1 << kPixelRequirementITS);
  if (fTrackNewITScut)                                 (fTrackCutBits |= 1 << kNewITSCut);
  //
  // dangerous cuts
  if (fTrackLengthInActiveZone>=90)  (fTrackCutBits |= 1 << kActiveZoneSmall);
  if (fTrackLengthInActiveZone>=100) (fTrackCutBits |= 1 << kActiveZone);
  if (fTrackLengthInActiveZone>=110) (fTrackCutBits |= 1 << kActiveZoneLarge);
  //
  if (fTrackTPCSignalN>=60) (fTrackCutBits |= 1 << kTPCSignalNSmall);
  if (fTrackTPCSignalN>=70) (fTrackCutBits |= 1 << kTPCSignalN);
  if (fTrackTPCSignalN>=80) (fTrackCutBits |= 1 << kTPCSignalNLarge);
  //
  // --------------------------------------------------------------------
  //                    Clean sample selections
  // --------------------------------------------------------------------
  //
  // 2.5 nsigma TOF protons and kaons for amplitude estimation
  if (cleanPrTOF) (fTrackCutBits |= 1 << kCleanPrTOF);
  if (cleanKaTOF) (fTrackCutBits |= 1 << kCleanKaTOF);
  //
  // Clean Kaons protons and deuterons
  if (cleanKaTOFTRD)        (fTrackCutBits |= 1 << kCleanKaTOFTRD);
  if (fTrackProbKaTOF>=0.8) (fTrackCutBits |= 1 << kTrackProbKaTOF);
  if (fTrackProbPrTOF>=0.8) (fTrackCutBits |= 1 << kTrackProbPrTOF);
  if (cleanDeTOF && cleanDeTPC) (fTrackCutBits |= 1 << kCleanDeTOF);
  //
  return fTrackCutBits;

}
//________________________________________________________________________
Bool_t AliAnalysisTaskEbyeIterPID::GetSystematicClassIndex(UInt_t cut,Int_t syst)
{
  /*
  syst:
  0 -->  Reference
  1 -->  CRows60
  2 -->  CRows100
  3 -->  Chi2TPCDmall
  4 -->  Chi2TPCLarge
  5 -->  DCAXYSmall
  6 -->  DCAXYLarge
  7 -->  VZSmall
  8 -->  VZLarge
  9 -->  EventVertexZSmall
  10 --> EventVertexZLarge
  11 --> RequireITSRefit
  12 --> NewITSCut
  13 --> PixelRequirementITS
  // extra settings
  14 --> kTPCSignalN,
  15 --> kActiveZone,
  16 --> kTPCSignalN+kActiveZone,
  17 --> kTPCSignalNSmall+kActiveZoneSmall,
  18 --> kTPCSignalNLarge,kActiveZoneLarge,
  //
  19 --> kTPCSignalNLarge,kActiveZoneLarge,
  20 --> kTPCSignalNLarge,kActiveZoneLarge,
  */

  const Int_t fnCutBins=10;
  Int_t fCutArr[fnCutBins]={0};

  switch(syst) {

    case kCutReference:   // 0 -->  Reference
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutCrossedRowsTPC60:  // 1 -->  CRows60
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC60,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutCrossedRowsTPC100:  // 2 -->  CRows100
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC100, kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutMaxChi2PerClusterTPCSmall:   // 3 -->  Chi2TPCSmall
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPCSmall, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutMaxChi2PerClusterTPCLarge:   // 4 -->  Chi2TPCLarge
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPCLarge, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutMaxDCAToVertexXYPtDepSmall:   // 5 -->  kMaxDCAToVertexXYPtDepSmall
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDepSmall, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutMaxDCAToVertexXYPtDepLarge:   // 6 -->  kMaxDCAToVertexXYPtDepLarge
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDepLarge, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutVertexZSmall:   // 7 -->  kVertexZSmall
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZSmall, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutVertexZLarge:   // 8 -->  kVertexZLarge
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZLarge, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutEventVertexZSmall:  // 9 -->  kEventVertexZSmall
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZSmall, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutEventVertexZLarge:  // 10 -->  kEventVertexZLarge
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZLarge, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutRequireITSRefit:  // 11 -->  no kRequireITSRefit
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, 1, 1, 1, 1 ,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutPixelRequirementITS:  // 12 --> no  kNewITSCut
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, 1, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutNewITSCut:  // 13 -->  kPixelRequirementITS
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, 1, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    // ----------------------------------------------------------------------------------------------------
    //                                              Dangerous cuts
    // ----------------------------------------------------------------------------------------------------
    //
    case kCutTPCSignalN:  // 14 -->  kTPCSignalN
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,kTPCSignalN,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutActiveZone:  // 15 -->  kActiveZone
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,kActiveZone};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutTPCSignalNActiveZone:  // 16 -->  kTPCSignalN + kActiveZone
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,kTPCSignalN,kActiveZone};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutTPCSignalNSmallActiveZoneSmall:  // 17 -->  tightest cut kTPCSignalNSmall + kActiveZoneSmall
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,kTPCSignalNSmall,kActiveZoneSmall};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutTPCSignalNLargeActiveZoneLarge:  // 18 -->  kTPCSignalNLarge + kActiveZoneLarge
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,kTPCSignalNLarge,kActiveZoneLarge};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    // ----------------------------------------------------------------------------------------------------
    //                                              Extra vZ cuts
    // ----------------------------------------------------------------------------------------------------
    //
    //
    case kCutEventVertexZALICE:  // 19 -->  kEventVertexZALICE
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZALICE, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutEventVertexZALICETight:  // 20 -->  kEventVertexZALICETight
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZALICETight, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;

    default:
    {
      Int_t fCutArrTmp[fnCutBins] = {0};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }

  }
  //
  //  Apply conditions
  for (Int_t i=0;i<fnCutBins;i++){
    if( ((cut >> fCutArr[i]) & 1) == 0 )
    {
      return kFALSE;
    }
  }
  return kTRUE;

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::FillTPCdEdxCheck()
{
  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillTPCdEdxCheck ===== " << std::endl;
  AliVEvent *event=InputEvent();
  //
  // Main track loop
  //
  Int_t mult=event->GetNumberOfTracks();
  for (Int_t itrack=0;itrack<event->GetNumberOfTracks();++itrack) {   // Track loop

    // get the track object
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack *track = fESD->GetTrack(itrack);

    // get track info
    if (!track->GetInnerParam()) continue;            // check if track in TPC
    fTPCSignal = track->GetTPCsignal();
    fPtot      = track->GetInnerParam()->GetP();
    fEta       = track->Eta();

    // Track cuts
    if (fTPCSignal>400)                continue;
    if (fPtot>2.)                      continue;
    if (fEta<fEtaDown && fEta>fEtaUp)  continue;
    if (track->GetTPCNcls()<80)        continue;

    // Fill the tree
    if(!fTreeSRedirector) return;
    (*fTreeSRedirector)<<"dEdxCheck"<<
    "dEdx="     << fTPCSignal <<    // dEdx of mc track
    "ptot="     << fPtot <<         // mc momentum
    "eta="      << fEta <<          // mc eta
    "mult="     << mult <<          // multiplicity
    "\n";
  }

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::SetSpecialV0Cuts(AliESDv0KineCuts* cuts)
{

  cuts->SetMode(0, 0);    // cuts->SetMode(mode, type); mode 0: purely kinematical selection   type 0: pp 1:PbPb not yet ready
  //
  // leg cuts
  cuts->SetNTPCclusters(50);
  cuts->SetTPCrefit(kTRUE);
  cuts->SetTPCchi2perCls(4.0);
  cuts->SetTPCclusterratio(0.6);
  cuts->SetNoKinks(kTRUE);
  //
  // gamma cuts
  cuts->SetGammaCutChi2NDF(10.0);
  Float_t cosPoint[2] = {0.0, 0.02};
  cuts->SetGammaCutCosPoint(cosPoint);
  Float_t cutDCA[2] = {0.0, 0.25};
  cuts->SetGammaCutDCA(cutDCA);
  Float_t vtxR[2] = {3.0, 90.0};
  cuts->SetGammaCutVertexR(vtxR);
  Float_t psiPairCut[2]={0.0,0.05};
  cuts->SetGammaCutPsiPair(psiPairCut);
  cuts->SetGammaCutInvMass(0.05);
  // K0s cuts
  cuts->SetK0CutChi2NDF(10.0);
  Float_t cosPointK0s[2] = {0.0, 0.02};
  cuts->SetK0CutCosPoint(cosPointK0s);
  Float_t cutDCAK0s[2] = {0.0, 0.2};
  cuts->SetK0CutDCA(cutDCAK0s);
  Float_t vtxRK0s[2] = {2.0, 30.0};
  cuts->SetK0CutVertexR(vtxRK0s);
  Float_t k0sInvMass[2] = {0.486, 0.508};
  cuts->SetK0CutInvMass(k0sInvMass);
  // Lambda and anti-Lambda cuts
  cuts->SetLambdaCutChi2NDF(10.0);
  Float_t cosPointLambda[2] = {0.0, 0.02};
  cuts->SetLambdaCutCosPoint(cosPointLambda);
  Float_t cutDCALambda[2] = {0.0, 0.2};
  cuts->SetLambdaCutDCA(cutDCALambda);
  Float_t vtxRLambda[2] = {2.0, 40.0};
  cuts->SetLambdaCutVertexR(vtxRLambda);
  Float_t lambdaInvMass[2] = {1.11, 1.12};
  cuts->SetLambdaCutInvMass(lambdaInvMass);

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::BinLogAxis(TH1 *h)
{
  //
  // Method for the correct logarithmic binning of histograms
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the BinLogAxis ===== " << std::endl;
  TAxis *axis       = h->GetXaxis();
  Int_t bins        = axis->GetNbins();

  Double_t from     = axis->GetXmin();
  Double_t to       = axis->GetXmax();
  std::vector<double>  newBins;
  newBins.resize(bins + 1);

  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);

  for (Int_t i = 1; i <= bins; i++) {
    newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins.data());

}
//________________________________________________________________________
Int_t AliAnalysisTaskEbyeIterPID::CountEmptyEvents(Int_t counterBin)
{

  //
  // count Empty Events
  //
  Int_t emptyCount=0;
  for (Int_t itrack=0;itrack<fESD->GetNumberOfTracks();++itrack) {   // Track loop
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack *track = fESD->GetTrack(itrack);
    if (!track->GetInnerParam()) continue;
    Float_t momtrack = track->GetInnerParam()->GetP();
    if (momtrack<fMomDown || momtrack>fMomUp) continue;
    if (!fESDtrackCuts->AcceptTrack(track)) continue;
    if (track->GetTPCsignalN()<60) continue;
    if (track->GetTPCsignal()>0) emptyCount++;
  }
  //
  // check if the event is empty
  if (emptyCount<1) {
    fHistEmptyEvent->Fill(counterBin);
    std::cout << " Info::marsland: Empty event in " << fChunkName << std::endl;
  }
  if (fUseCouts) std::cout << " Info::marsland: ====== EVENT IS COOL GO AHEAD ======= " << std::endl;
  return emptyCount;

}
//
void AliAnalysisTaskEbyeIterPID::PrintNumInBinary(UInt_t num)
{
  TString bin="";
  Int_t numberOfBits = sizeof(UInt_t)*8;
  for (Int_t i=numberOfBits-1; i>=0; i--) {
    Bool_t isBitSet = (num & (1<<i));
    if (isBitSet) {
      bin+="1";
    } else {
      bin+="0";
    }
  }
  std::cout << "Info::marsland: fTrackCutBits = " << bin << std::endl;
}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  std::cout << " Info::marsland: ===== In the Terminate ===== " << std::endl;

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::CreateEventInfoTree()
{

  if (fUseCouts) std::cout << " Info::marsland: ===== In the CreateEventInfoTree ===== " << std::endl;

  Int_t tpcClusterMultiplicity   = fESD->GetNumberOfTPCClusters();
  const AliMultiplicity *multObj = fESD->GetMultiplicity();
  Int_t itsNumberOfTracklets   = multObj->GetNumberOfTracklets();

  TVectorF phiCountA(36);
  TVectorF phiCountC(36);
  TVectorF phiCountAITS(36);
  TVectorF phiCountCITS(36);
  TVectorF phiCountAITSonly(36);
  TVectorF phiCountCITSonly(36);
  TVectorF tzeroMult(24);  for (Int_t i=1;i<24;i++) tzeroMult[i] = 0.;
  TVectorF vzeroMult(64);  for (Int_t i=1;i<64;i++) vzeroMult[i] = 0.;
  TVectorF itsClustersPerLayer(6); for (Int_t i=1;i<6;i++) itsClustersPerLayer[i] = 0.;
  //
  for (Int_t i=1;i<37;i++){
    phiCountA[i-1] = fHistPhiTPCcounterA->GetBinContent(i);
    phiCountC[i-1] = fHistPhiTPCcounterC->GetBinContent(i);
    phiCountAITS[i-1] = fHistPhiTPCcounterAITS->GetBinContent(i);
    phiCountCITS[i-1] = fHistPhiTPCcounterCITS->GetBinContent(i);
    phiCountAITSonly[i-1] = fHistPhiITScounterA->GetBinContent(i);
    phiCountCITSonly[i-1] = fHistPhiITScounterC->GetBinContent(i);
  }
  //
  // Additional counters for ITS TPC V0 and T0
  // AliESDFMD* esdFMD = fESD->GetFMDData();
  // AliFMDFloatMap fmdMult = esdFMD->MultiplicityMap();
  const AliESDTZERO *esdTzero = fESD->GetESDTZERO();
  const Double32_t *t0amp=esdTzero->GetT0amplitude();
  ULong64_t triggerMask = fESD->GetTriggerMask();
  Short_t   eventMult = fESD->GetNumberOfTracks();
  Int_t eventMultESD = fESD->GetNumberOfESDTracks();
  Int_t tpcTrackBeforeClean =fESD->GetNTPCTrackBeforeClean();

  for (Int_t i=0;i<24;i++) { tzeroMult[i] = t0amp[i]; }
  for (Int_t i=0;i<64;i++) { vzeroMult[i] = fESD->GetVZEROData()-> GetMultiplicity(i); }
  for (Int_t i=0;i<6;i++)  { itsClustersPerLayer[i] = multObj->GetNumberOfITSClusters(i); }

  if(!fTreeSRedirector) return;
  DumpDownScaledTree();
  (*fTreeSRedirector)<<"eventInfo"<<
  "run="                  << fRunNo                 <<  // run Number
  "pileupbit="            << fPileUpBit             <<
  "bField="               << fBField                <<  // run Number
  "gid="                  << fEventGID              <<  // global event ID
  "timestamp="            << fTimeStamp             <<  // timestamp
  "triggerMask="          << triggerMask            <<  //trigger mask
  "intrate="              << fIntRate               <<  // interaction rate
  "centV0M="              << fCentrality            <<  // centrality
  "centrality.="          << fEventInfo_CentralityEstimates   <<  // track counter
  "vz="                   << fVz                    <<  // vertex Z
  "tpcvz="                << fTPCvZ                 <<
  "spdvz="                << fSPDvZ                 <<
  "tpcMult="              << fTPCMult               <<  //  TPC multiplicity
  "eventMult="            << fEventMult             <<  //  event multiplicity
  "eventMultESD="         << eventMultESD           <<  //  event multiplicity ESD
  "primMult="             << fNContributors         <<  //  #prim tracks
  "tpcClusterMult="       << tpcClusterMultiplicity <<  // tpc cluster multiplicity
  "tpcTrackBeforeClean="  << tpcTrackBeforeClean <<   // tpc track before cleaning
  "itsTracklets="         << itsNumberOfTracklets   <<  // number of ITS tracklets
  //
  // "fmdMult.="             << &fmdMult               <<  // T0 multiplicity
  "tZeroMult.="           << &tzeroMult             <<  // T0 multiplicity
  "vZeroMult.="           << &vzeroMult             <<  // V0 multiplicity
  "itsClustersPerLayer.=" << &itsClustersPerLayer   <<  // its clusters per layer
  "trackCounters.="       << fEventInfo_CacheTrackCounters    <<  // track counter
  "trackdEdxRatio.="      << fEventInfo_CacheTrackdEdxRatio   <<  // dEdx conter
  "trackNcl.="            << fEventInfo_CacheTrackNcl         <<  // nCluster counter
  "trackChi2.="           << fEventInfo_CacheTrackChi2        <<  // Chi2 counter
  "trackMatchEff.="       << fEventInfo_CacheTrackMatchEff    <<  // Chi2 counter
  "trackTPCCountersZ.="   << fEventInfo_CacheTrackTPCCountersZ    <<  // Chi2 counter
  "hisTPCVertexA.="       << fEventInfo_HisTPCVertexA         <<  // Chi2 counter
  "hisTPCVertexC.="       << fEventInfo_HisTPCVertexC         <<  // Chi2 counter
  "hisTPCVertex.="        << fEventInfo_HisTPCVertex          <<  // Chi2 counter
  "hisTPCVertexACut.="    << fEventInfo_HisTPCVertexACut      <<  // Chi2 counter
  "hisTPCVertexCCut.="    << fEventInfo_HisTPCVertexCCut      <<  // Chi2 counter
  "phiTPCdcarA.="         << fEventInfo_PhiTPCdcarA           <<  // track counter
  "phiTPCdcarC.="         << fEventInfo_PhiTPCdcarC           <<  // dEdx conter
  "phiCountA.="           << &phiCountA             <<  // TPC track count on A side
  "phiCountC.="           << &phiCountC             <<  // TPC track count on C side
  "phiCountAITS.="        << &phiCountAITS          <<  // track count fitted ITS on A side
  "phiCountCITS.="        << &phiCountCITS          <<  // track count fitted ITS on C side
  "phiCountAITSOnly.="    << &phiCountAITSonly      <<  // track count only ITS on A side
  "phiCountCITSOnly.="    << &phiCountCITSonly      <<  // track count only ITS on C side
  "\n";

}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::CalculateEventInfo()
{

  if (fUseCouts) std::cout << " Info::marsland: ===== In the CalculateEventInfo ===== " << std::endl;

  AliVEvent *event=InputEvent();
  CacheTPCEventInformation();
  //
  //
  const Int_t kNclTPCcut=60;
  const Int_t kDCACut=5;  // 5 cm primary cut
  const Int_t kMindEdxClustersRegion=15;
  const Float_t kTglCut=1.5;
  const Float_t kPtCut=0.100;
  const Float_t kDCAtpcNULL = -10000;
  // const Float_t kNTrackletCut=1.5;
  //
  fEventInfo_PhiTPCdcarA->Zero();
  fEventInfo_PhiTPCdcarC->Zero();
  fEventInfo_CacheTrackCounters->Zero();   // track counter
  fEventInfo_CacheTrackdEdxRatio->Zero(); // dedx info counter
  fEventInfo_CacheTrackNcl->Zero();       // ncl counter
  fEventInfo_CacheTrackChi2->Zero();      // chi2 counter
  fEventInfo_CacheTrackMatchEff->Zero();  // matchEff counter
  //
  if (fHistPhiTPCcounterA)    fHistPhiTPCcounterA->Reset();
  if (fHistPhiTPCcounterC)    fHistPhiTPCcounterC->Reset();
  if (fHistPhiTPCcounterAITS) fHistPhiTPCcounterAITS->Reset();
  if (fHistPhiTPCcounterCITS) fHistPhiTPCcounterCITS->Reset();
  if (fHistPhiITScounterA)    fHistPhiITScounterA->Reset();
  if (fHistPhiITScounterC)    fHistPhiITScounterC->Reset();
  //
  //
  Int_t nNumberOfTracks = event->GetNumberOfTracks();
  Float_t tpcDCAarrPhiA[36][nNumberOfTracks];
  Float_t tpcDCAarrPhiC[36][nNumberOfTracks];
  for (Int_t i=0;i<36;i++){
    for (Int_t j=0;j<nNumberOfTracks;j++){
      tpcDCAarrPhiA[i][j]=kDCAtpcNULL;
      tpcDCAarrPhiC[i][j]=kDCAtpcNULL;
    }
  }
  //
  // --------------------------------------------------------------
  //      Track LOOP
  // --------------------------------------------------------------
  //
  AliTPCdEdxInfo tpcdEdxInfo;
  TRandom r;
  for (Int_t itrack=0;itrack<nNumberOfTracks;++itrack)
  {

    //
    Double_t eta=-100., phiTPC=0.,sectorNumber=0., tpcdEdx=0., ptotTPC=0.;
    //
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack *track = fESD->GetTrack(itrack);
    if (track == NULL) continue;
    //
    Double_t tgl        = track->Pz()/track->Pt();
    Double_t phiGlobal  = track->Phi()-TMath::Pi(); // ????
    // Int_t sign          = track->GetSign();
    Double_t phi        = track->GetParameterAtRadius(85,5,7);
    Double_t sectorNumbertmp = (9*phi/TMath::Pi()+18*(phi<0));
    eta = track->Eta();
    if (TMath::Abs(eta)>0.9) continue;
    Bool_t isOnITS = track->IsOn(AliESDtrack::kITSrefit);
    Bool_t isOnTPC = track->IsOn(AliESDtrack::kTPCrefit);
    // Bool_t isOnTRD = track->IsOn(AliESDtrack::kTRDrefit);
    //
    // --------------------------------------------------------------
    //      TPC track information
    // --------------------------------------------------------------
    //
    if (track->GetInnerParam()) {
      tpcdEdx = track->GetTPCsignal();
      ptotTPC = track->GetInnerParam()->GetP();
      phiTPC  = track->GetInnerParam()->GetParameterAtRadius(85,5,7);
      sectorNumber = (9*phiTPC/TMath::Pi()+18*(phiTPC<0));
    }
    //
    // --------------------------------------------------------------
    //      Count only ITS tracks
    // --------------------------------------------------------------
    //
    if ( isOnITS && !isOnTPC ) {
      if (TMath::Abs(phi)>1e-10){
        if (tgl>0) fHistPhiITScounterA->Fill(sectorNumbertmp);
        if (tgl<0) fHistPhiITScounterC->Fill(sectorNumbertmp);
      }
    }
    //
    if (!track->GetInnerParam()) continue;  // ????
    if (track->IsOn(AliVTrack::kTPCout)==kFALSE)  continue;  // ????
    (*fEventInfo_CacheTrackCounters)[4]++;      // all TPC track with out flag
    // TPC track counters with DCAZ
    for (Int_t izCut=1; izCut<4; izCut++){
      Float_t impactParam[2];
      track->GetImpactParameters(impactParam[0],impactParam[1]);
      if (TMath::Abs(impactParam[0])>kDCACut) continue;
      if (TMath::Abs(track->GetInnerParam()->GetParameter()[1])<10.*(izCut+1.)) (*fEventInfo_CacheTrackTPCCountersZ)[izCut]++;
      if (TMath::Abs(impactParam[1])<10.*(izCut+1.)) (*fEventInfo_CacheTrackTPCCountersZ)[izCut+4]++;
    }
    //
    //
    Float_t dcaRPhi, dcaZ;
    track->GetImpactParameters(dcaRPhi, dcaZ);
    Int_t nclTPC    = track->GetTPCncls(); if (nclTPC<1) nclTPC=-1;
    Int_t nclITS    = track->GetITSNcls(); if (nclITS<1) nclITS=-1;
    Int_t nclTRD    = track->GetTRDncls(); if (nclTRD<1) nclTRD=-1;
    // Int_t nclTOF    = track->IsOn(AliVTrack::kTOFout);
    Double_t chi2TPC = TMath::Sqrt(TMath::Abs(track->GetTPCchi2()/nclTPC));
    Double_t chi2ITS = TMath::Sqrt(TMath::Abs(track->GetITSchi2()));
    Double_t chi2TRD = TMath::Sqrt(TMath::Abs(track->GetTRDchi2()));
    // Double_t qP      = track->Charge()/track->P();
    // Double_t ptot0   = track->GetP();
    //
    // --------------------------------------------------------------
    //      Some track selections
    // --------------------------------------------------------------
    //
    if (nclTPC<kNclTPCcut) continue;
    if (TMath::Abs(tgl)>kTglCut) continue;
    if (track->Pt()<kPtCut) continue;
    if (TMath::Abs(dcaRPhi)>kDCACut || TMath::Abs(dcaZ)>kDCACut) continue;
    // if ( !( isOnITS||isOnTRD ) ) continue;   // ????
    //
    // --------------------------------------------------------------
    //      Fill TPC dca information for a given phi bin for each track
    // --------------------------------------------------------------
    //
    Int_t phiBin = fHistPhi->FindBin(phi)-1;
    Float_t pTPC[2],covTPC[3];          // p[0]=fdTPC; p[1]=fzTPC; cov[0]=fCddTPC; cov[1]=fCdzTPC; cov[2]=fCzzTPC;
    track->GetImpactParametersTPC(pTPC,covTPC);
    if (tgl>0) tpcDCAarrPhiA[phiBin][itrack]=pTPC[0];
    if (tgl<0) tpcDCAarrPhiC[phiBin][itrack]=pTPC[0];
    //
    // --------------------------------------------------------------
    //      TPC Phi counter
    // --------------------------------------------------------------
    //
    (*fEventInfo_CacheTrackCounters)[5]++;
    if (TMath::Abs(phiTPC)>1e-10){
      if (tgl>0) fHistPhiTPCcounterA->Fill(sectorNumber);
      if (tgl<0) fHistPhiTPCcounterC->Fill(sectorNumber);
      if(isOnITS){
        if (tgl>0) fHistPhiTPCcounterAITS->Fill(sectorNumber);
        if (tgl<0) fHistPhiTPCcounterCITS->Fill(sectorNumber);
      }
    }
    //
    // --------------------------------------------------------------
    //      track counter after pile up  ????
    // --------------------------------------------------------------
    //
    Bool_t pileUpCut=  ( (nclITS>2) || (nclTRD>40) );
    if (pileUpCut==kFALSE) continue;
    if (TMath::Min(chi2TPC,100.)<0) continue;
    (*fEventInfo_CacheTrackCounters)[1]++;
    //
    Bool_t itsOK=track->IsOn(AliVTrack::kITSout) && nclITS>2  && chi2ITS>0;
    Bool_t trdOK=track->IsOn(AliVTrack::kTRDout) && nclTRD>35 && chi2TRD>0;
    Bool_t tofOK=track->IsOn(AliVTrack::kTOFout);
    //
    // --------------------------------------------------------------
    //      number of clusters cut
    // --------------------------------------------------------------
    //
    (*fEventInfo_CacheTrackNcl)[4]+=track->GetTPCncls(0, 63);
    (*fEventInfo_CacheTrackNcl)[5]+=track->GetTPCncls(64, 127);
    (*fEventInfo_CacheTrackNcl)[6]+=track->GetTPCncls(128, 159);
    (*fEventInfo_CacheTrackNcl)[1] += nclTPC;
    (*fEventInfo_CacheTrackChi2)[1]+= (chi2TPC>0) ? TMath::Sqrt(chi2TPC):2;   // sometimes negative chi2?

    if (itsOK && track->GetTPCdEdxInfo(tpcdEdxInfo)){

      Bool_t isOK=(tpcdEdxInfo.GetNumberOfCrossedRows(0)>kMindEdxClustersRegion);
      isOK&=(tpcdEdxInfo.GetNumberOfCrossedRows(1)>kMindEdxClustersRegion);
      isOK&=(tpcdEdxInfo.GetNumberOfCrossedRows(2)>kMindEdxClustersRegion);
      isOK&=((tpcdEdxInfo.GetSignalMax(0)>0) && (tpcdEdxInfo.GetSignalMax(1)>0) && (tpcdEdxInfo.GetSignalMax(2)>0));
      isOK&=((tpcdEdxInfo.GetSignalTot(0)>0) && (tpcdEdxInfo.GetSignalTot(1)>0) && (tpcdEdxInfo.GetSignalTot(2)>0));
      isOK&=(itsOK||trdOK);      // stronger pile-up cut requiring ITS or TRD

      if (isOK) {
        (*fEventInfo_CacheTrackCounters)[6]+=1;         // Counter with accepted TPC dEdx info
        (*fEventInfo_CacheTrackdEdxRatio)[0]+=TMath::Log(tpcdEdxInfo.GetSignalMax(3));
        (*fEventInfo_CacheTrackdEdxRatio)[1]+=TMath::Log(tpcdEdxInfo.GetSignalTot(3));
        (*fEventInfo_CacheTrackdEdxRatio)[2]+=TMath::Log(tpcdEdxInfo.GetSignalMax(0)/tpcdEdxInfo.GetSignalTot(0));
        (*fEventInfo_CacheTrackdEdxRatio)[3]+=TMath::Log(tpcdEdxInfo.GetSignalMax(1)/tpcdEdxInfo.GetSignalTot(1));
        (*fEventInfo_CacheTrackdEdxRatio)[4]+=TMath::Log(tpcdEdxInfo.GetSignalMax(2)/tpcdEdxInfo.GetSignalTot(2));
        (*fEventInfo_CacheTrackdEdxRatio)[5]+=TMath::Log(tpcdEdxInfo.GetSignalMax(3)/tpcdEdxInfo.GetSignalTot(3));
        (*fEventInfo_CacheTrackdEdxRatio)[6]+=TMath::Log(tpcdEdxInfo.GetSignalMax(1)/tpcdEdxInfo.GetSignalMax(0));
        (*fEventInfo_CacheTrackdEdxRatio)[7]+=TMath::Log(tpcdEdxInfo.GetSignalMax(1)/tpcdEdxInfo.GetSignalMax(2));
        (*fEventInfo_CacheTrackdEdxRatio)[8]+=TMath::Log(tpcdEdxInfo.GetSignalTot(1)/tpcdEdxInfo.GetSignalTot(0));
        (*fEventInfo_CacheTrackdEdxRatio)[9]+=TMath::Log(tpcdEdxInfo.GetSignalTot(1)/tpcdEdxInfo.GetSignalTot(2));
        //
        // --------------------------------------------------------------
        //      dEdx counter wrt splines and Bethe bloch
        // --------------------------------------------------------------
        //
        Float_t closestPar[3];    // closestPar[0] --> closest spline, Int_t(closestPar[1]) --> particle index,  closestPar[2] --> corresponding particle mass
        GetExpecteds(track,closestPar);
        // std::cout << " aaaaaaaa  " << tpcdEdx << "     "  << closestPar[0] << "    " << closestPar[1] << "      " << closestPar[2]  << "   " << ptotTPC << std::endl;
        (*fEventInfo_CacheTrackdEdxRatio)[10]+=TMath::Log(tpcdEdx/closestPar[0]);    // ???
        (*fEventInfo_CacheTrackdEdxRatio)[11]+=TMath::Log((tpcdEdxInfo.GetSignalMax(0)/50.)/AliExternalTrackParam::BetheBlochAleph(ptotTPC/closestPar[2]));
        (*fEventInfo_CacheTrackdEdxRatio)[12]+=TMath::Log((tpcdEdxInfo.GetSignalMax(1)/50.)/AliExternalTrackParam::BetheBlochAleph(ptotTPC/closestPar[2]));
        (*fEventInfo_CacheTrackdEdxRatio)[13]+=TMath::Log((tpcdEdxInfo.GetSignalMax(2)/50.)/AliExternalTrackParam::BetheBlochAleph(ptotTPC/closestPar[2]));
        (*fEventInfo_CacheTrackdEdxRatio)[14]+=TMath::Log((tpcdEdxInfo.GetSignalMax(3)/50.)/AliExternalTrackParam::BetheBlochAleph(ptotTPC/closestPar[2]));
        (*fEventInfo_CacheTrackdEdxRatio)[15]+=TMath::Log((tpcdEdxInfo.GetSignalTot(0)/50.)/AliExternalTrackParam::BetheBlochAleph(ptotTPC/closestPar[2]));
        (*fEventInfo_CacheTrackdEdxRatio)[16]+=TMath::Log((tpcdEdxInfo.GetSignalTot(1)/50.)/AliExternalTrackParam::BetheBlochAleph(ptotTPC/closestPar[2]));
        (*fEventInfo_CacheTrackdEdxRatio)[17]+=TMath::Log((tpcdEdxInfo.GetSignalTot(2)/50.)/AliExternalTrackParam::BetheBlochAleph(ptotTPC/closestPar[2]));
        (*fEventInfo_CacheTrackdEdxRatio)[18]+=TMath::Log((tpcdEdxInfo.GetSignalTot(3)/50.)/AliExternalTrackParam::BetheBlochAleph(ptotTPC/closestPar[2]));
        (*fEventInfo_CacheTrackdEdxRatio)[19]+=TMath::Log((tpcdEdxInfo.GetSignalMax(0)/50.));
        (*fEventInfo_CacheTrackdEdxRatio)[20]+=TMath::Log((tpcdEdxInfo.GetSignalMax(1)/50.));
        (*fEventInfo_CacheTrackdEdxRatio)[21]+=TMath::Log((tpcdEdxInfo.GetSignalMax(2)/50.));
        (*fEventInfo_CacheTrackdEdxRatio)[22]+=TMath::Log((tpcdEdxInfo.GetSignalMax(3)/50.));
        (*fEventInfo_CacheTrackdEdxRatio)[23]+=TMath::Log((tpcdEdxInfo.GetSignalTot(0)/50.));
        (*fEventInfo_CacheTrackdEdxRatio)[24]+=TMath::Log((tpcdEdxInfo.GetSignalTot(1)/50.));
        (*fEventInfo_CacheTrackdEdxRatio)[25]+=TMath::Log((tpcdEdxInfo.GetSignalTot(2)/50.));
        (*fEventInfo_CacheTrackdEdxRatio)[26]+=TMath::Log((tpcdEdxInfo.GetSignalTot(3)/50.));

      }
    }

    if (itsOK) {  // ITS
      (*fEventInfo_CacheTrackCounters)[0]++;
      (*fEventInfo_CacheTrackNcl)[0] += nclITS;
      (*fEventInfo_CacheTrackChi2)[0] += TMath::Min(TMath::Sqrt(chi2ITS),10.); // cutoff chi2 10
      (*fEventInfo_CacheTrackMatchEff)[2]+=trdOK;
      (*fEventInfo_CacheTrackMatchEff)[3]+=tofOK;
      (*fEventInfo_CacheTrackChi2)[4]+= (fTrackChi2TPC>0) ? TMath::Sqrt(fTrackChi2TPC):2; // TPC chi2 in case prolongation to ITS
      // long tracks properties
      if (nclITS>4){
        (*fEventInfo_CacheTrackCounters)[7]++;
        (*fEventInfo_CacheTrackNcl)[7] += nclITS;
        (*fEventInfo_CacheTrackChi2)[7]+=TMath::Min(TMath::Sqrt(chi2ITS),10.);
      }
    }
    if (trdOK) {// TRD    ///TODO - why chi2TRD could be smaller than 0?
      (*fEventInfo_CacheTrackCounters)[2]++;
      (*fEventInfo_CacheTrackNcl)[2] += nclTRD;
      (*fEventInfo_CacheTrackChi2)[2] += TMath::Sqrt(chi2TRD);
      (*fEventInfo_CacheTrackMatchEff)[0]+=itsOK;
      (*fEventInfo_CacheTrackChi2)[5]+= (fTrackChi2TPC>0) ? TMath::Sqrt(fTrackChi2TPC):2; // TPC chi2 in case prolongation to TRD
      if (nclTRD>80){
        (*fEventInfo_CacheTrackCounters)[8]++;
        (*fEventInfo_CacheTrackNcl)[8] += nclTRD;
        (*fEventInfo_CacheTrackChi2)[8]+=TMath::Min(TMath::Sqrt(chi2TRD),10.);
      }
    }
    if (tofOK) {  // TOF
      (*fEventInfo_CacheTrackCounters)[3]++;
      (*fEventInfo_CacheTrackdEdxRatio)[3] += 1;   // dummy for the moment
      (*fEventInfo_CacheTrackChi2)[3]+= 1;   //
    }
  } // end of track LOOP
  //
  // ======================================================================
  //  calculate event averages
  // ======================================================================
  //
  for (Int_t i=0; i<9; i++) if ((*fEventInfo_CacheTrackCounters)[i]>0) (*fEventInfo_CacheTrackdEdxRatio)[i]/=(*fEventInfo_CacheTrackCounters)[i];
  for (Int_t i=0; i<4; i++) if ((*fEventInfo_CacheTrackCounters)[i]>0) (*fEventInfo_CacheTrackChi2)[i]/=(*fEventInfo_CacheTrackCounters)[i];

  for (Int_t i=4; i<7; i++)  if ((*fEventInfo_CacheTrackCounters)[1]>0) (*fEventInfo_CacheTrackdEdxRatio)[i]/=(*fEventInfo_CacheTrackCounters)[1];
  //
  if ((*fEventInfo_CacheTrackCounters)[6]>0){
    for (Int_t i=0; i<27; i++)   (*fEventInfo_CacheTrackdEdxRatio)[i]/=(*fEventInfo_CacheTrackCounters)[6];
  }
  //
  // conditional matching efficiency and chi2
  if ((*fEventInfo_CacheTrackCounters)[0]>0){
    (*fEventInfo_CacheTrackMatchEff)[2]/=(*fEventInfo_CacheTrackCounters)[0];  // TRD if ITS
    (*fEventInfo_CacheTrackMatchEff)[3]/=(*fEventInfo_CacheTrackCounters)[0];  // TOF if ITS
    (*fEventInfo_CacheTrackChi2)[4]/=(*fEventInfo_CacheTrackCounters)[0];
  }
  if ((*fEventInfo_CacheTrackCounters)[2]>0) {
    (*fEventInfo_CacheTrackMatchEff)[0]/=(*fEventInfo_CacheTrackCounters)[2];
    (*fEventInfo_CacheTrackChi2)[5]/=(*fEventInfo_CacheTrackCounters)[2];
  } //ITS if TRD
  (*fEventInfo_CacheTrackCounters)[9]=event->GetNumberOfTracks();  // original number of ESDtracks
  //
  //
  for (Int_t iphi=0; iphi<36; iphi++){

    // count nonzero entries
    Int_t countNonZerosA=0;
    Int_t countNonZerosC=0;
    for (Int_t itrack=0;itrack<nNumberOfTracks;itrack++) {
      if (tpcDCAarrPhiA[iphi][itrack]>kDCAtpcNULL) countNonZerosA++;
      if (tpcDCAarrPhiC[iphi][itrack]>kDCAtpcNULL) countNonZerosC++;
    }
    //
    // create arrays from nonzero entries for A side
    Float_t tmpA[countNonZerosA];
    Int_t j=0;
    for (Int_t itrack=0;itrack<nNumberOfTracks;itrack++) {
      if (tpcDCAarrPhiA[iphi][itrack]>kDCAtpcNULL) { tmpA[j]=tpcDCAarrPhiA[iphi][itrack]; j++; }
    }
    //
    // create arrays from nonzero entries for C side
    Float_t tmpC[countNonZerosC];
    Int_t k=0;
    for (Int_t itrack=0;itrack<nNumberOfTracks;itrack++) {
      if (tpcDCAarrPhiC[iphi][itrack]>kDCAtpcNULL) { tmpC[k]=tpcDCAarrPhiC[iphi][itrack]; k++; }
    }
    (*fEventInfo_PhiTPCdcarA)[iphi]=TMath::Median(countNonZerosA,tmpA);
    (*fEventInfo_PhiTPCdcarC)[iphi]=TMath::Median(countNonZerosC,tmpC);

  }

}
//________________________________________________________________________
Int_t AliAnalysisTaskEbyeIterPID::CacheTPCEventInformation()
{

  AliVEvent *event=InputEvent();
  const Int_t kNCRCut=80;
  const Double_t kDCACut=5;
  const Float_t knTrackletCut=1.5;
  // FILL DCA histograms
  fEventInfo_HisTPCVertexA->Reset();
  fEventInfo_HisTPCVertexC->Reset();
  fEventInfo_HisTPCVertexACut->Reset();
  fEventInfo_HisTPCVertexCCut->Reset();
  fEventInfo_HisTPCVertex->Reset();

  Int_t nTracks=event->GetNumberOfTracks();
  Int_t selected=0;
  for (Int_t iTrack=0; iTrack<nTracks; iTrack++){
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack * pTrack = fESD->GetTrack(iTrack);
    if (pTrack==nullptr) continue;
    if (pTrack->IsOn(AliESDtrack::kTPCin)==0) continue;
    if (pTrack->GetTPCClusterInfo(3,1)<kNCRCut) continue;
    Float_t dcaxy,dcaz;
    pTrack->GetImpactParameters(dcaxy,dcaz);
    if (TMath::Abs(dcaxy)>kDCACut) continue;
    pTrack->SetESDEvent(fESD);
    selected++;
    if ((pTrack->GetNumberOfTRDClusters()/20.+pTrack->GetNumberOfITSClusters())>knTrackletCut){
      fEventInfo_HisTPCVertex->Fill(pTrack->GetTPCInnerParam()->GetZ());
      if (pTrack->GetTgl()>0) fEventInfo_HisTPCVertexACut->Fill(pTrack->GetTPCInnerParam()->GetZ());
      if (pTrack->GetTgl()<0) fEventInfo_HisTPCVertexCCut->Fill(pTrack->GetTPCInnerParam()->GetZ());
    }else{
      if (pTrack->GetTgl()>0) fEventInfo_HisTPCVertexA->Fill(pTrack->GetTPCInnerParam()->GetZ());
      if (pTrack->GetTgl()<0) fEventInfo_HisTPCVertexC->Fill(pTrack->GetTPCInnerParam()->GetZ());
    }
  }
  if (fUseCouts) printf(" Info::marsland: ===== In the CacheTPCEventInformation:  %d\n",selected);
  return selected;
}
//________________________________________________________________________
void AliAnalysisTaskEbyeIterPID::DumpDownScaledTree()
{

  if (fUseCouts) std::cout << " Info::marsland: ===== In the DumpDownScaledTree ===== " << std::endl;

  AliVEvent *event=InputEvent();
  const Int_t kNclTPCcut=60;
  // const Int_t kDCACut=5;  // 5 cm primary cut
  // const Int_t kMindEdxClustersRegion=15;
  const Float_t kTglCut=1.5;
  const Float_t kPtCut=0.100;
  //
  // --------------------------------------------------------------
  //  Counter only for the TPC track multiplicity
  // --------------------------------------------------------------
  //
  // Go into track loop
  AliTPCdEdxInfo tpcdEdxInfo;
  TRandom r;
  Int_t nNumberOfTracks = event->GetNumberOfTracks();
  for (Int_t itrack=0;itrack<nNumberOfTracks;++itrack)
  {
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack *track = fESD->GetTrack(itrack);
    if (!track) continue;
    if (!fESDtrackCutsLoose->AcceptTrack(track))  continue;    // Loose cuts
    if (!(track->GetTPCsignalN()>0)) continue;
    Int_t tpcCrossedRows=0, tpcSignalN=0;
    Double_t eta=-100.;
    // Double_t sectorNumber=0.;
    Double_t tgl  = track->Pz()/track->Pt();
    Double_t phi  = track->Phi()-TMath::Pi(); // ????
    Double_t pt   = track->Pt();
    Int_t sign    = track->GetSign();
    Double_t phi85 = track->GetParameterAtRadius(85,5,7);
    // Double_t sectorNumbertmp = (9*phi85/TMath::Pi()+18*(phi85<0));
    ULong64_t flag = track->GetStatus();
    eta = track->Eta();
    // Bool_t isOnITS = track->IsOn(AliESDtrack::kITSrefit);
    // Bool_t isOnTRD = track->IsOn(AliESDtrack::kTRDrefit);
    // Bool_t isOnTPC = track->IsOn(AliESDtrack::kTPCrefit);
    Int_t nclTPC   = track->GetTPCncls(); if (nclTPC<1) nclTPC=-1;
    Int_t nclITS   = track->GetITSNcls(); if (nclITS<1) nclITS=-1;
    Int_t nclTRD   = track->GetTRDncls(); if (nclTRD<1) nclTRD=-1;
    // Int_t nclTOF   = track->IsOn(AliVTrack::kTOFout);
    Double_t chi2TPC = TMath::Sqrt(TMath::Abs(track->GetTPCchi2()/nclTPC));
    Double_t chi2ITS = TMath::Sqrt(TMath::Abs(track->GetITSchi2()));
    Double_t chi2TRD = TMath::Sqrt(TMath::Abs(track->GetTRDchi2()));
    Double_t itsdEdx = track->GetITSsignal();
    Double_t trddEdx = track->GetTRDsignal();
    Double_t tpcdEdx = track->GetTPCsignal();
    Double_t ptot0   = track->GetP();
    Double_t qP      = track->Charge()/track->P();
    Float_t dcaRPhi, dcaZ;
    track->GetImpactParameters(dcaRPhi, dcaZ);
    //
    // --------------------------------------------------------------
    //      Preparation for downscaled tree
    // --------------------------------------------------------------
    //
    Double_t phiTPC=-100.;
    Double_t ptotTPC=0.;
    Double_t lengthInActiveZone=0.;
    Float_t pTPC[2],covTPC[3];          // p[0]=fdTPC; p[1]=fzTPC; cov[0]=fCddTPC; cov[1]=fCdzTPC; cov[2]=fCzzTPC;
    if (track->GetInnerParam()){
      track->GetTPCdEdxInfo(tpcdEdxInfo);
      phiTPC  = track->GetInnerParam()->GetParameterAtRadius(85,5,7);
      ptotTPC = track->GetInnerParam()->GetP();
      lengthInActiveZone = track->GetLengthInActiveZone(1,3,230, track->GetBz(),0,0);
      track->GetImpactParametersTPC(pTPC,covTPC);
      tpcCrossedRows = track->GetTPCCrossedRows();
      tpcSignalN = track->GetTPCsignalN();
      Float_t closestPar[3];
      GetExpecteds(track,closestPar);
    }
    UChar_t itsclmap = track->GetITSClusterMap();
    Float_t pv[2],cov[3];
    track->GetImpactParameters(pv,cov); // p[0]=fD; p[1]=fZ; cov[0]=fCdd; cov[1]=fCdz; cov[2]=fCzz;
    // High dEdx its conditon
    TDatabasePDG *pdg = TDatabasePDG::Instance();
    Double_t mproton = pdg->GetParticle(2212)->Mass();  // Double_t mproton = 9.3827199e-01;
    Bool_t itsHighDeDx = kFALSE;
    if (ptot0>0.1 && ptot0<100 && itsdEdx>0){
      itsHighDeDx = ( (TMath::Log(itsdEdx/AliExternalTrackParam::BetheBlochSolid(ptot0/mproton))>11.5) && (nclITS>4) );
    }
    //
    // --------------------------------------------------------------
    //   Fill downscaled tree
    // --------------------------------------------------------------
    //
    if (nclTPC<kNclTPCcut) continue;
    if (TMath::Abs(tgl)>kTglCut) continue;
    if (track->Pt()<kPtCut) continue;
    if ( (r.Rndm()*(qP*qP) < 0.05) || itsHighDeDx )
    {
      if(!fTreeSRedirector) return;
      (*fTreeSRedirector)<<"dscaled"<<
      "centV0M="              << fCentrality           <<  //  centrality
      "cent.="                << fEventInfo_CentralityEstimates  <<  // track counter
      "gid="                  << fEventGID             <<  //  global event ID
      "intrate="              << fIntRate              <<  // interaction rate
      "timestamp="            << fTimeStamp            <<  // timestamp
      "tpcMult="              << fTPCMult              <<  //  TPC multiplicity
      "primMult="             << fNContributors        <<  //  #prim tracks
      "run="                  << fRunNo                <<  //  run number
      "vz="                   << fVz                   <<  //  vertex Z
      "tpcvz="                << fTPCvZ                <<  //  TPC event vertex
      "spdvz="                << fSPDvZ                <<
      //
      "expel="                << fDEdxEl               <<
      "exppi="                << fDEdxPi               <<
      "expka="                << fDEdxKa               <<
      "exppr="                << fDEdxPr               <<
      "expde="                << fDEdxDe               <<
      "tofnsigmael="          << fNSigmasElTOF         <<
      "tofnsigmapi="          << fNSigmasPiTOF         <<
      "tofnsigmaka="          << fNSigmasKaTOF         <<
      "tofnsigmapr="          << fNSigmasPrTOF         <<
      "tofnsigmade="          << fNSigmasDeTOF         <<
      //
      "flag="                 << flag                  <<
      "eventmult="            << fEventMult            <<  //  event multiplicity
      "nclTPC="               << nclTPC                <<  //  number of cl used in dEdx
      "tpcSignalN="           << tpcSignalN            <<  //  number of cl used in dEdx
      "cRows="                << tpcCrossedRows        <<  //  crossed rows
      "lengthInActiveZone="   << lengthInActiveZone    <<  //  fTrackLengthInActiveZone in TPC
      "ptotTPC="              << ptotTPC               <<  //  TPC momentum
      "ptot="                 << ptot0                 <<  //  TPC momentum
      "eta="                  << eta                   <<  //  eta
      "phi="                  << phi                   <<  //  ph
      "phiTPC="               << phiTPC                <<
      "phi85="                << phi85                 <<
      "sign="                 << sign                  <<  //  charge
      "tpcdEdx="              << tpcdEdx               <<  //  dEdx of the track
      "itsdEdx="              << itsdEdx               <<
      "trddEdx="              << trddEdx               <<
      "fdTPC="                << pTPC[0]               <<
      "fzTPC="                << pTPC[1]               <<
      "fD="                   << pv[0]                 <<  //  DCAxy
      "fZ="                   << pv[1]                 <<  //  DCAz
      "dEdxInfo.="            << &tpcdEdxInfo          <<  //  TPC dEdx info
      "ncltpc="               << nclTPC                <<  //  #ITS clusters
      "nclits="               << nclITS                <<  //  #ITS clusters
      "ncltrd="               << nclTRD                <<  //  #TRD clusters
      "chi2tpc="              << chi2TPC               <<  //  TPC chi2
      "chi2its="              << chi2ITS               <<  //  ITS chi2
      "chi2trd="              << chi2TRD               <<  //  TRD chi2
      "qP="                   << qP                    <<  //  charge/momentu,
      "tgl="                  << tgl                   <<  //  tangent
      "pt="                   << pt                    <<  //  pT
      "itsclmap="             << itsclmap              <<  //  vertex Z
      "\n";
    }

  } // end of track LOOP

}
