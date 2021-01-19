//
// Original Author:
//         Created:  Fri December 7, 2017 by  N.Z. Jomhari
//                   with contributions from  A. Geiser
//                                            A. Anuar
//                   and pieces similar to previous analysis examples
// $Id$
// ..
//
// ***************************************************************************
//  Higgs-to-four-lepton analysis example at research level                  *
// ***************************************************************************
//                                                                           *
// Built upon the DEMO setup provided by CMS open data access team,          *
// expanded/upgraded to contain a pedagocigal analysis example for the       *
// approximate reproducttion of the Higgs to four lepton mass spectrum       *
// published in CMS-HIG-12-028                                               *
// Phys.Lett. B716 (2012) 30-61,  arXiv:1207.7235                            *
//                                                                           *
// This research level example is a strongly simplified reimplementation of  *
// parts of the original CMS Higgs to four lepton analysis                   *
//                                                                           *
// The published reference plot which is being approximated in this example  *
// (in addition to many auxiliary plots) is                                  *
// https://inspirehep.net/record/1124338/files/H4l_mass_v3.png               *
// Other Higgs final states (e.g. Higgs to two photons), which were also     *
// part of the same CMS paper and strongly contributed to the Higgs          *
// discovery, are not covered by this example.                               *
//                                                                           *
// The example addresses users who feel they have at least some minimal      *
// understanding of the content of the CMS paper and of the meaning of this  *
// reference plot, which can be reached via (separate) educational exercises.*
// A Root ntuple containing the kinematic information about the event        *
// candidates, which could be used for educational purposes, is part of the  *
// Root output file.                                                         *
//                                                                           *
// The analysis code provided here recodes the spirit of the original        *
// analysis and (approximately) recodes many of the original cuts on         *
// original data objects, but does not provide the original analysis code    *
// itself. Also, for the sake of simplicity, it skips some of the more       *
// advanced analysis methods of the original paper, and does not use any     *
// corrections beyond those already implicit in the objects used, and no     *
// systematic uncertainties.                                                 *
// Another reason why the published spectrum is only reproduced very         *
// approximately is that the data sets only partially overlap, and that the  *
// legacy software version and corresponding calibrations differ from those  *
// of the original paper.                                                    *
// Nevertheless, the example provides a qualitative insight about how the    *
// original result was obtained.                                             *
//                                                                           *
// In addition to the documented core results, the resulting  Root files     *
// also contain many undocumented plots which grew as a side product from    *
// setting up this example and earlier examples.                             *
// And it contains an ntuple with candidate four-vectors as mentioned above. *
// ***************************************************************************

// ***************************************************************************
// Analysis strategy                                                         *
// ***************************************************************************
//
// The analysis strategy is the following: Get the histograms for the 4mu    *
// and 2mu2e final states from the DoubleMu datasets and those for 4e final  *
// state from the DoubleElectron dataset. This avoids double counting due to *
// trigger overlaps.                                                         *
// The code itself is agnostic with respect to the input data set used, and  *
// the appropriate histograms have to selected at the subsequent root        *
// analysis step.                                                            *
// ***************************************************************************

// system include files
#include <memory>
#include <vector>
#include <algorithm>
#include <utility>

// user include files, general
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//------ EXTRA HEADER FILES--------------------//
#include "math.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"

// for Root histogramming
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TLorentzVector.h"

// for tracking information
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

// for vertex information
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// for muon information
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

//for beamspot information
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//for electron informaton
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

// for particle flow information
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

// class declaration
class HiggsDemoAnalyzerGit : public edm::EDAnalyzer
{
public:
	explicit HiggsDemoAnalyzerGit(const edm::ParameterSet &);
	~HiggsDemoAnalyzerGit();

private:
	virtual void beginJob();
	virtual void analyze(const edm::Event &, const edm::EventSetup &);
	virtual void endJob();
	bool providesGoodLumisection(const edm::Event &iEvent);

	// Declare Root histograms or tree
	// For a description of their content see below

	TH1D *h_globalmu_size;

	TTree *ttree;

	UInt_t nRun, nEvt, nLumi;

	Bool_t Muon_isPFMuon[9999];
	Bool_t Muon_isPFIsolationValid[9999];
	Bool_t Muon_hasglobalTrack[9999];
	Float_t Muon_p[9999];
	Float_t Muon_px[9999];
	Float_t Muon_py[9999];
	Float_t Muon_pz[9999];
	Float_t Muon_pt[9999];
	Float_t Muon_eta[9999];
	Float_t Muon_phi[9999];
	Float_t Muon_globalTrack_chi2[9999];
	Float_t Muon_globalTrack_ndof[9999];
	Float_t Muon_globalTrack_normalizedChi2[9999];
	Float_t Muon_pfRelIso04_all[9999];
	Float_t Muon_dxy[9999];
	Float_t Muon_dxyErr[9999];
	Float_t Muon_dz[9999];
	Float_t Muon_dzErr[9999];
	Int_t Muon_charge[9999];
	UInt_t nMuon;

	Bool_t Electron_isEB[9999];
	Bool_t Electron_isEE[9999];
	Float_t Electron_p[9999];
	Float_t Electron_px[9999];
	Float_t Electron_py[9999];
	Float_t Electron_pz[9999];
	Float_t Electron_pt[9999];
	Float_t Electron_eta[9999];
	Float_t Electron_sc_eta[9999];
	Float_t Electron_phi[9999];
	Float_t Electron_misshits[9999];
	Float_t Electron_pfRelIso_all[9999];
	Float_t Electron_dxy[9999];
	Float_t Electron_dxyErr[9999];
	Float_t Electron_dz[9999];
	Float_t Electron_dzErr[9999];
	Int_t Electron_charge[9999];
	UInt_t nElectron;
};
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

HiggsDemoAnalyzerGit::HiggsDemoAnalyzerGit(const edm::ParameterSet &iConfig)
{

	// *****************************************************************
	// This is the main analysis routine
	// The goal is to approximately reproduce the Higgs-to-four-lepton
	// mass spectrum from HIG-12-028
	// *****************************************************************

	// now do what ever initialization is needed
	edm::Service<TFileService> fs;

	// ************************************
	// book histograms and set axis labels
	// (called once for initialization)
	// ************************************

	// Global Muon (GM) size
	h_globalmu_size = fs->make<TH1D>("NGMuons", "GMuon size", 10, 0., 10.);
	h_globalmu_size->GetXaxis()->SetTitle("Number of GMuons");
	h_globalmu_size->GetYaxis()->SetTitle("Number of Events");
}

HiggsDemoAnalyzerGit::~HiggsDemoAnalyzerGit()
{
	//do anything here that needs to be done at destruction time
	// (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------//
void HiggsDemoAnalyzerGit::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{

	// **********************************************
	// here each relevant event will get analyzed
	// **********************************************

	nRun = iEvent.run();
	nEvt = (iEvent.id()).event(); // iEvent: no class named event()
	nLumi = iEvent.luminosityBlock();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example", pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif

	// Event is to be analyzed

	edm::LogInfo("Demo")
		<< "Starting to analyze \n"
		<< "Event number: " << (iEvent.id()).event()
		<< ", Run number: " << iEvent.run()
		<< ", Lumisection: " << iEvent.luminosityBlock();

	//------------------Load (relevant) Event information------------------------//
	// INFO: Getting Data From an Event
	// https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookChapter4#GetData
	// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMGetDataFromEvent#get_ByLabel
	// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAodDataTable
	// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideRecoDataTable

	// INFO: globalMuons
	// NB: note that when using keyword "globalMuons" getByLabel-function returns
	// reco::TrackCollection
	edm::Handle<reco::TrackCollection> tracks;
	iEvent.getByLabel("generalTracks", tracks);

	edm::Handle<reco::TrackCollection> gmuons;
	iEvent.getByLabel("globalMuons", gmuons);

	edm::Handle<reco::MuonCollection> muons;
	iEvent.getByLabel("muons", muons);

	edm::Handle<reco::BeamSpot> beamSpotHandle;
	iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);

	edm::Handle<reco::VertexCollection> primvtxHandle;
	iEvent.getByLabel("offlinePrimaryVertices", primvtxHandle);

	edm::Handle<reco::GsfElectronCollection> electrons;
	iEvent.getByLabel("gsfElectrons", electrons);

	reco::BeamSpot beamSpot;
	if (beamSpotHandle.isValid())
	{
		beamSpot = *beamSpotHandle;
	}
	else
	{
		edm::LogInfo("Demo")
			<< "No beam spot available from EventSetup \n";
	}

	reco::VertexCollection primvtx;
	if (primvtxHandle.isValid())
	{
		primvtx = *primvtxHandle;
	}
	else
	{
		edm::LogInfo("Demo")
			<< "No primary vertex available from EventSetup \n";
	}

	////////////////////////////////////////////////////////////////////////////////
	///////////////////////// Reco Muon Collection Start ///////////////////////////
	////////////////////////////////////////////////////////////////////////////////

	//******************************************************************************
	// This muon collection is being used here for the Higgs->4 lepton analysis    *
	//******************************************************************************

	// Loop over muons size and select good muons

	int im = 0;
	for (unsigned u = 0; u < muons->size(); u++)
	{
		const reco::Muon &itMuon = (*muons)[u];

		// math::XYZPoint point(beamSpot.position());
		math::XYZPoint point(primvtx[0].position());

		Muon_isPFMuon[im] = itMuon.isPFMuon();
		Muon_isPFIsolationValid[im] = itMuon.isPFIsolationValid();
		Muon_hasglobalTrack[im] = (itMuon.globalTrack()).isNonnull();
		Muon_p[im] = itMuon.p();
		Muon_px[im] = itMuon.px();
		Muon_py[im] = itMuon.py();
		Muon_pz[im] = itMuon.pz();
		Muon_pt[im] = itMuon.pt();
		Muon_eta[im] = itMuon.eta();
		Muon_phi[im] = itMuon.phi();
		Muon_pfRelIso04_all[im] = ((itMuon.pfIsolationR04()).sumChargedHadronPt +
								   (itMuon.pfIsolationR04()).sumNeutralHadronEt +
								   (itMuon.pfIsolationR04()).sumPhotonEt) /
								  itMuon.pt();
		Muon_charge[im] = itMuon.charge();

		if (itMuon.isPFMuon() && itMuon.isPFIsolationValid() && (itMuon.globalTrack()).isNonnull())
		{
			Muon_globalTrack_chi2[im] = (itMuon.globalTrack())->chi2();
			Muon_globalTrack_ndof[im] = (itMuon.globalTrack())->ndof();
			Muon_globalTrack_normalizedChi2[im] = (itMuon.globalTrack())->normalizedChi2();
			Muon_dxy[im] = itMuon.globalTrack()->dxy(point);
			Muon_dxyErr[im] = itMuon.globalTrack()->d0Error();
			Muon_dz[im] = itMuon.globalTrack()->dz(point);
			Muon_dzErr[im] = itMuon.globalTrack()->dzError();
		}
		else
		{
			Muon_globalTrack_chi2[im] = -99;
			Muon_globalTrack_ndof[im] = -99;
			Muon_globalTrack_normalizedChi2[im] = -99;
			Muon_dxy[im] = -99;
			Muon_dxyErr[im] = -99;
			Muon_dz[im] = -99;
			Muon_dzErr[im] = -99;
		}
		im++;
	}
	nMuon = im;

	int in = 0;
	for (unsigned te = 0; te < electrons->size(); te++)
	{
		const reco::GsfElectron &iElectron = (*electrons)[te];

		// math::XYZPoint point(beamSpot.position());
		math::XYZPoint point(primvtx[0].position());

		if (iElectron.passingPflowPreselection())
		{
			Electron_isEB[in] = iElectron.isEB();
			Electron_isEE[in] = iElectron.isEE();
			Electron_p[in] = iElectron.p();
			Electron_px[in] = iElectron.px();
			Electron_py[in] = iElectron.py();
			Electron_pz[in] = iElectron.pz();
			Electron_pt[in] = iElectron.pt();
			Electron_eta[in] = iElectron.eta();
			Electron_sc_eta[in] = (iElectron.superCluster())->eta();
			Electron_phi[in] = iElectron.phi();
			Electron_misshits[in] = ((iElectron.gsfTrack())->trackerExpectedHitsInner()).numberOfHits();
			Electron_pfRelIso_all[in] = ((iElectron.pfIsolationVariables()).chargedHadronIso +
										 (iElectron.pfIsolationVariables()).neutralHadronIso +
										 (iElectron.pfIsolationVariables()).photonIso) /
										iElectron.pt();
			Electron_charge[in] = iElectron.charge();
			Electron_dxy[in] = iElectron.gsfTrack()->dxy(point);
			Electron_dxyErr[in] = iElectron.gsfTrack()->d0Error();
			Electron_dz[in] = iElectron.gsfTrack()->dz(point);
			Electron_dzErr[in] = iElectron.gsfTrack()->dzError();
			in++;
		}
	}
	nElectron = in;
	ttree->Fill();

} // HiggsDemoAnalyzerGit::analyze ends

// ------ method called once each job just before starting event loop ---------//

void HiggsDemoAnalyzerGit::beginJob()
{

	// *******************************************************
	// book the ntuple for the surviving 4 lepton candidates *
	// in the mass range 70 < m4l < 181 GeV                  *
	// *******************************************************
	ttree = new TTree("ttree", "ttree");
	ttree->Branch("nRun", &nRun, "nRun/I");
	ttree->Branch("nEvt", &nEvt, "nEvt/I");
	ttree->Branch("nLumi", &nLumi, "nLumi/I");
	ttree->Branch("nMuon", &nMuon, "nMuon/I");
	ttree->Branch("Muon_isPFMuon", &Muon_isPFMuon, "Muon_isPFMuon[nMuon]/O");
	ttree->Branch("Muon_isPFIsolationValid", &Muon_isPFIsolationValid, "Muon_isPFIsolationValid[nMuon]/O");
	ttree->Branch("Muon_hasglobalTrack", &Muon_hasglobalTrack, "Muon_hasglobalTrack[nMuon]/O");
	ttree->Branch("Muon_p", &Muon_p, "Muon_p[nMuon]/F");
	ttree->Branch("Muon_px", &Muon_px, "Muon_px[nMuon]/F");
	ttree->Branch("Muon_py", &Muon_py, "Muon_py[nMuon]/F");
	ttree->Branch("Muon_pz", &Muon_pz, "Muon_pz[nMuon]/F");
	ttree->Branch("Muon_pt", &Muon_pt, "Muon_pt[nMuon]/F");
	ttree->Branch("Muon_eta", &Muon_eta, "Muon_eta[nMuon]/F");
	ttree->Branch("Muon_phi", &Muon_phi, "Muon_phi[nMuon]/F");
	ttree->Branch("Muon_globalTrack_chi2", &Muon_globalTrack_chi2, "Muon_globalTrack_chi2[nMuon]/F");
	ttree->Branch("Muon_globalTrack_ndof", &Muon_globalTrack_ndof, "Muon_globalTrack_ndof[nMuon]/F");
	ttree->Branch("Muon_globalTrack_normalizedChi2", &Muon_globalTrack_normalizedChi2, "Muon_globalTrack_normalizedChi2[nMuon]/F");
	ttree->Branch("Muon_pfRelIso04_all", &Muon_pfRelIso04_all, "Muon_pfRelIso04_all[nMuon]/F");
	ttree->Branch("Muon_dxy", &Muon_dxy, "Muon_dxy[nMuon]/F");
	ttree->Branch("Muon_dxyErr", &Muon_dxyErr, "Muon_dxyErr[nMuon]/F");
	ttree->Branch("Muon_dz", &Muon_dz, "Muon_dz[nMuon]/F");
	ttree->Branch("Muon_dzErr", &Muon_dzErr, "Muon_dzErr[nMuon]/F");
	ttree->Branch("Muon_charge", &Muon_charge, "Muon_charge[nMuon]/I");

	ttree->Branch("nElectron", &nElectron, "nElectron/I");
	ttree->Branch("Electron_isEB", &Electron_isEB, "Electron_isEB[nElectron]/O");
	ttree->Branch("Electron_isEE", &Electron_isEE, "Electron_isEE[nElectron]/O");
	ttree->Branch("Electron_p", &Electron_p, "Electron_p[nElectron]/F");
	ttree->Branch("Electron_px", &Electron_px, "Electron_px[nElectron]/F");
	ttree->Branch("Electron_py", &Electron_py, "Electron_py[nElectron]/F");
	ttree->Branch("Electron_pz", &Electron_pz, "Electron_pz[nElectron]/F");
	ttree->Branch("Electron_pt", &Electron_pt, "Electron_pt[nElectron]/F");
	ttree->Branch("Electron_eta", &Electron_eta, "Electron_eta[nElectron]/F");
	ttree->Branch("Electron_sc_eta", &Electron_sc_eta, "Electron_sc_eta[nElectron]/F");
	ttree->Branch("Electron_phi", &Electron_phi, "Electron_phi[nElectron]/F");
	ttree->Branch("Electron_misshits", &Electron_misshits, "Electron_misshits[nElectron]/F");
	ttree->Branch("Electron_pfRelIso_all", &Electron_pfRelIso_all, "Electron_pfRelIso_all[nElectron]/F");
	ttree->Branch("Electron_dxy", &Electron_dxy, "Electron_dxy[nElectron]/F");
	ttree->Branch("Electron_dxyErr", &Electron_dxyErr, "Electron_dxyErr[nElectron]/F");
	ttree->Branch("Electron_dz", &Electron_dz, "Electron_dz[nElectron]/F");
	ttree->Branch("Electron_dzErr", &Electron_dzErr, "Electron_dzErr[nElectron]/F");
	ttree->Branch("Electron_charge", &Electron_charge, "Electron_charge[nElectron]/I");
}

// ------------ method called once each job just after ending the event loop  ------------
void HiggsDemoAnalyzerGit::endJob()
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiggsDemoAnalyzerGit);
