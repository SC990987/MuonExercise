// -*- C++ -*-
//
// Package:    CMSDASExercises/MuonExercise1
// Class:      MuonExercise1
// 
/**\class MuonExercise1 MuonExercise1.cc CMSDASExercises/MuonExercise1/plugins/MuonExercise1.cc

 Description: Short Muon exercise for CMSDAS

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Norbert Neumeister
//         Created:  Thu, 14 Dec 2017 09:31:13 GMT
//
//

// system include files
#include <memory>
#include <iomanip>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "Math/GenVector/VectorUtil.h"
#include <iostream>
#include <fstream>
//
// class declaration
//

class MuonExercise1 : public edm::one::EDAnalyzer<edm::one::SharedResources> {

  public:

    explicit MuonExercise1(const edm::ParameterSet&);
    ~MuonExercise1();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:

    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------
     
    // TFileService
    edm::Service<TFileService> fs;
    edm::EDGetTokenT<reco::Muon> muonCollToken;
    edm::EDGetTokenT<reco::GenParticle> genCollToken;
};

//
// constructors and destructor
//
MuonExercise1::MuonExercise1(const edm::ParameterSet& iConfig) {

  usesResource("TFileService");

  edm::InputTag muonTag("muons");
  
  edm::InputTag genPartTag("genParticles");

  muonCollToken = consumes<reco::Muon>(muonTag);
  
  genCollToken = consumes<reco::GenParticle>(genPartTag);
 
}


MuonExercise1::~MuonExercise1() {
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
ofstream MyFile("filename.txt");
// ------------ method called for each event  ------------
void MuonExercise1::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;


  size_t ngen(0);
  //
  // RECO Muons
  //
  edm::Handle<reco::Muon> muonColl;
  iEvent.getByToken(muonCollToken, muonColl);
  //
  // GEN Muons
  //  
  edm::Handle<reco::GenParticle> genColl;
  iEvent.getByToken(genCollToken, genColl);
  for (auto&& genPart : *(genColl.product())) {
    if ( genPart.status() == 1 && std::abs(genPart.pdgId()) == 13 && fabs(genPart.eta()) < 2.4 && genPart.pt() > 1.5 ) ngen++;
  }
 
//
// Reconstruction efficiency code 
//

for (auto&& genPart1: *(genColl.product())) { 
	for(auto&& genPart2: *(genColl.product())){
		int nrec = 0;
		//
		// This part of the code just applies filters to the generated muons, i.e. selects certain eta and pt and filters events with 
		// Two generate muons
		//

		if ( genPart1.status() == 1 && std::abs(genPart1.pdgId()) == 13 && fabs(genPart1.eta()) < 2.4 && genPart1.pt() > 1.5 && ngen == 2 ){
		if ( genPart2.status() == 1 && std::abs(genPart2.pdgId()) == 13 && fabs(genPart2.eta()) < 2.4 && genPart2.pt() > 1.5 && ngen == 2 ){
			//
			// INvariant Mass Calc
			//
			 
			/*float inv = 2.0*(genPart1.pt())*(genPart2.pt())*(coshf(genPart1.eta() - genPart2.eta()) - cos(genPart1.phi() - genPart2.phi()));
                	float inv1 = sqrt(inv);
			if (inv1 >= 0.0001){
                	h_inv->Fill(inv1);}*/			
			//
			// This part of the code calculates the deltaR between the two generated muons
			// An if statement is used then to filter out any deltaR of 0: When the deltaR is calculate between a pair of identical
			// gen muons
			//

			float drgen = deltaR(genPart1,genPart2);
				if (drgen > 0.00001){
					
					//
					// Performing Gen Matching to compute efficiency
					//

					for (auto&& muon: *(muonColl.product())){
						if (muon.isGlobalMuon()){ //This selects just muons from the global muon collection
							float dr2 = deltaR(muon,genPart1); // computes the delta R between the reco and gen muon
							float dr3 = deltaR(muon,genPart2);// computes the delta R between the reco and gen muon
								if (dr2 < 0.3 && dr3 > 0.3) ++nrec; // Matches the reco and gen muons
								if (dr2 > 0.3 && dr3 < 0.3) ++nrec; // Matches the reco and gen muons
								if (dr2 < 0.3 && dr3 < 0.3) ++nrec; // Matches the reco and gen muons
									}
										}
					cout  << drgen << ", " << nrec << endl;
					MyFile  << drgen << ", " << nrec << endl;
						}
					}
				}
			}
		}



}//end of analyzer


// ------------ method called once each job just before starting event loop  ------------
void MuonExercise1::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void MuonExercise1::endJob() {
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuonExercise1::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonExercise1);
//MyFile.close();

