// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
// new includes
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <cmath>
#include "DataFormats/Math/interface/Vector3D.h"
#include <vector>
#include "TLorentzVector.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

//////////Written by Ethan Cannaert [2019] Using Stephen Mrenna's thrust algorithm [2007]////////////////////////////////

class Thrust
{	
	typedef math::XYZVector Vector;

	public:
		
		Thrust(std::vector<reco::LeafCandidate> cands, int nconstituents);
		double thrust(void);
		int ncyc(void);
		Vector axis(void);

	private:
		bool isConverged1 = false;
		bool isConverged2 = false;
		bool isConverged3 = false;
		bool isConverged4 = false;

		int ncyc_ = 0;
		double thrust_ = 0.;
		Vector axis_;

		Vector det_init_dir(std::vector<reco::LeafCandidate> cands, int nconstituents, int num);
		Vector calc_Tnext(std::vector<reco::LeafCandidate> cands, int nconstituents, Vector Tvec);

		double det_eps(Vector Tvec, double px, double py, double pz);
		bool check_convergence(Vector Tvec, Vector Tnext);
		double calc_thrust(std::vector<reco::LeafCandidate> cands, int nconstituents, Vector axis_);
		double calc_ptot(double px, double py, double pz);

};

