////////////////////////////HELP////////////////////////////////
//////////////Analyzes the diquark -> chi chi hadronic decay//////////////
////////////////Last updated Nov 9 2020/////////////////////////////////


// system include files
#include <fastjet/JetDefinition.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/tools/Filter.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ActiveAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>

#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
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
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/PtEtaPhiMass.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
//#include "PhysicsTools/CandUtils/interface/Thrust.h"
#include <TTree.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <cmath>
#include "TLorentzVector.h"
#include "TVector3.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include  "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include <algorithm>   

#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include <DataFormats/Math/interface/deltaR.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <string>
using namespace reco;
typedef math::XYZTLorentzVector LorentzVector;
typedef math::XYZVector Vector;

class diquarkAnalyzerBR : public edm::EDAnalyzer 
{
public:
   explicit diquarkAnalyzerBR(const edm::ParameterSet&);
private:
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   double calc_mag(double px,double py, double pz);
   bool isgoodjet(const float pt, const float fatjet_sd_mass, const float eta, const float NHF,const float NEMF, const size_t NumConst,const float CHF,const int CHM, const float MUF, const float CEMF, const int NumNeutralParticles);
   double calc_mpp_beta(pat::Jet* iM);
   //edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartToken_; 
   edm::EDGetTokenT<std::vector<pat::Jet>> fatJetToken_;
   edm::EDGetTokenT<edm::TriggerResults> triggerBits_;

   TTree * tree;

  // std::vector<std::string> triggs = {"HLT_AK8PFJet800"};
   std::vector<std::string> triggs = {"HLT_PFHT1050"};

   int eventnum = 0;
   int nfatjets = 0;
   int raw_nfatjets;
   int nfatjets_minus1;
   double jet_pt[100], jet_eta[100], jet_mass[100], jet_dr[100], raw_jet_mass[100],raw_jet_pt[100],raw_jet_phi[100];
   double jet_beta[100], beta_T[100], AK4_mass_20[100],AK4_mass_30[100],AK4_mass_50[100],AK4_mass_70[100],AK4_mass_100[100],AK4_mass_150[100];
   double b_W_deltar,h_t_deltar;
   int jet_ndaughters[100], jet_nAK4[100],jet_nAK4_20[100],jet_nAK4_30[100],jet_nAK4_50[100],jet_nAK4_70[100],jet_nAK4_100[100],jet_nAK4_150[100];
};


diquarkAnalyzerBR::diquarkAnalyzerBR(const edm::ParameterSet& iConfig)
{
   //genPartToken_ = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genPartCollection"));
   fatJetToken_ =    consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("fatJetCollection"));
   triggerBits_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"));

   edm::Service<TFileService> fs;      

   tree = fs->make<TTree>("tree", "tree");

   tree->Branch("nfatjets", &nfatjets, "nfatjets/I");

   tree->Branch("jet_ndaughters", jet_ndaughters, "jet_ndaughters[nfatjets]/I");
   //tree->Branch("nfatjets_minus1",&nfatjets_minus1, "nfatjets_minus1/I");
   tree->Branch("jet_nAK4", jet_nAK4, "jet_nAK4[nfatjets]/I");
   tree->Branch("jet_nAK4_20", jet_nAK4_20, "jet_nAK4_20[nfatjets]/I");
   tree->Branch("jet_nAK4_30", jet_nAK4_30, "jet_nAK4_30[nfatjets]/I");
   tree->Branch("jet_nAK4_50", jet_nAK4_50, "jet_nAK4_50[nfatjets]/I");
   tree->Branch("jet_nAK4_70", jet_nAK4_70, "jet_nAK4_70[nfatjets]/I");
   tree->Branch("jet_nAK4_100", jet_nAK4_100, "jet_nAK4_100[nfatjets]/I");
   tree->Branch("jet_nAK4_150", jet_nAK4_150, "jet_nAK4_150[nfatjets]/I");

   tree->Branch("beta_T", beta_T, "beta_T[nfatjets]/D");

   tree->Branch("AK4_mass_20", AK4_mass_20, "AK4_mass_20[nfatjets]/D");
   tree->Branch("AK4_mass_30", AK4_mass_30, "AK4_mass_30[nfatjets]/D");
   tree->Branch("AK4_mass_50", AK4_mass_50, "AK4_mass_50[nfatjets]/D");
   tree->Branch("AK4_mass_70", AK4_mass_70, "AK4_mass_70[nfatjets]/D");
   tree->Branch("AK4_mass_100", AK4_mass_100, "AK4_mass_100[nfatjets]/D");
   tree->Branch("AK4_mass_150", AK4_mass_150, "AK4_mass_150[nfatjets]/D");

   tree->Branch("jet_beta", jet_beta, "jet_beta[nfatjets]/D");

   tree->Branch("jet_pt", jet_pt, "jet_pt[nfatjets]/D");
   tree->Branch("jet_eta", jet_eta, "jet_eta[nfatjets]/D");
   tree->Branch("jet_mass", jet_mass, "jet_mass[nfatjets]/D");
   tree->Branch("raw_nfatjets",&raw_nfatjets,"raw_nfatjets/I");
   tree->Branch("raw_jet_mass",raw_jet_mass,"raw_jet_mass[raw_nfatjets]/D");
   tree->Branch("raw_jet_pt",raw_jet_pt,"raw_jet_pt[raw_nfatjets]/D");
   tree->Branch("raw_jet_phi",raw_jet_phi,"raw_jet_phi[raw_nfatjets]/D");


}


double diquarkAnalyzerBR::calc_mag(double px,double py, double pz)
{
   return sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
}


bool diquarkAnalyzerBR::isgoodjet(const float pt, const float fatjet_sd_mass, const float eta, const float NHF,const float NEMF, const size_t NumConst,const float CHF,const int CHM, const float MUF, const float CEMF, const int NumNeutralParticles)
{
   if( (abs(eta) > 2.4) || (pt < 500.) ||(fatjet_sd_mass < 50.0)) return false;

   if ((NHF>0.9) || (NEMF>0.9) || (NumConst<1) || (CHF<0.) || (CHM<0) || (MUF > 0.8) || (CEMF > 0.8)) 
      {
         return false;
      }
   else{ return true;}
}


void diquarkAnalyzerBR::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{




////////////////Triggers//////////////////////////
   
   edm::Handle<edm::TriggerResults> triggerBits;
   iEvent.getByToken(triggerBits_, triggerBits);
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

   bool pass = false;
   //std::string trigname= "HLT_AK8PFJet550_v";
   std::string trigname= "HLT_PFHT1050_v";
   for (unsigned int i = 0; i < triggerBits->size(); ++i) 
   {
      const std::string name = names.triggerName(i);
      const bool accept = triggerBits->accept(i);
      if ((name.find(trigname) != std::string::npos) &&(accept)) pass =true;
   }
   
   if(!pass) return;


   nfatjets = 0;
   raw_nfatjets = 0;

////////////////Jets//////////////////////////////////////
   edm::Handle<std::vector<pat::Jet> > fatJets;
   iEvent.getByToken(fatJetToken_, fatJets);


   for(auto iJet = fatJets->begin(); iJet != fatJets->end(); iJet++)         ////////Over AK8 Jets
   { 
      if(!(iJet->isPFJet())) return;

      if(sqrt(pow(iJet->energy(),2) - pow(iJet->p(),2)) < 30000)
         {
            raw_jet_pt[raw_nfatjets]   = iJet->pt();
            raw_jet_phi[raw_nfatjets]  = iJet->phi();
            raw_jet_mass[raw_nfatjets] = sqrt(pow(iJet->energy(),2) - pow(iJet->p(),2));
            raw_nfatjets++;
         }
      if(!isgoodjet(iJet->pt(),iJet->userFloat("ak8PFJetsPuppiSoftDropMass"), iJet->eta(),iJet->neutralHadronEnergyFraction(), iJet->neutralEmEnergyFraction(),iJet->numberOfDaughters(),iJet->chargedHadronEnergyFraction(),iJet->chargedMultiplicity(),iJet->muonEnergyFraction(),iJet->chargedEmEnergyFraction(),iJet->neutralMultiplicity())) continue;



      jet_pt[nfatjets]         = iJet->pt();
      jet_eta[nfatjets]        = iJet->eta();
      jet_mass[nfatjets]       = sqrt(pow(iJet->energy(),2) - pow(iJet->p(),2));
      jet_ndaughters[nfatjets] = iJet->numberOfDaughters();


///////////////////////////////mpp frame boost calculation////////////////////
      Vector jet_axis(iJet->px()/iJet->p(),iJet->py()/iJet->p(),iJet->pz()/iJet->p());
      TLorentzVector jet_particles(iJet->px(),iJet->py(),iJet->pz(),iJet->energy());
      double min_pp = 999999999.;
      double min_boost = 0.;

      for(int iii=0;iii<10000;iii++)
      {
         TLorentzVector jet_particles_ = jet_particles;
         double beta_cand = iii/10000.;
         jet_particles_.Boost(-beta_cand*jet_axis.X(),-beta_cand*jet_axis.Y(),-beta_cand*jet_axis.Z());
         if(abs( ( jet_particles_.Px()*iJet->px()+jet_particles_.Py()*iJet->py() +jet_particles_.Pz()*iJet->py() )/iJet->p() ) < min_pp) 
            {
               min_boost = beta_cand; 
               min_pp = abs( ( jet_particles_.Px()*iJet->px()+jet_particles_.Py()*iJet->py() +jet_particles_.Pz()*iJet->py() )/iJet->p() ) ;
            }
      }
      double beta_mag = min_boost;

      jet_beta[nfatjets] = beta_mag;
      beta_T[nfatjets] = beta_mag*sin(iJet->theta());
      ///Recluster//
      //double beta_mag = iJet->p()/iJet->energy();
      Vector beta(beta_mag*iJet->px()/iJet->p(),beta_mag*iJet->py()/iJet->p(),beta_mag*iJet->pz()/iJet->p());
/////////////////////////////////////////////////////////////////////////////////////////////////


      std::vector<fastjet::PseudoJet> cands_;

      for (unsigned int i=0; i<iJet->numberOfDaughters();i++)
      {
         const reco::Candidate* iJ = iJet->daughter(i);
         const pat::PackedCandidate* cand_begin = (pat::PackedCandidate*) iJ;
         double puppiweight = cand_begin->puppiWeight();
         TLorentzVector w(puppiweight*iJ->px(),puppiweight*iJ->py(),puppiweight*iJ->pz(),puppiweight*iJ->energy());
         w.Boost(-beta.X(),-beta.Y(),-beta.Z());
         cands_.push_back(fastjet::PseudoJet(w.Px(),w.Py(),w.Pz(),w.E()));

      }

      double R = 0.4;
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
      fastjet::ClusterSequence cs_jet(cands_, jet_def); 
      std::vector<fastjet::PseudoJet> jetsFJ_jet = fastjet::sorted_by_E(cs_jet.inclusive_jets(5.0));
      
      jet_nAK4[nfatjets]    = 0;
      jet_nAK4_20[nfatjets] = 0;
      jet_nAK4_30[nfatjets] = 0;
      jet_nAK4_50[nfatjets] = 0;
      jet_nAK4_70[nfatjets] = 0;
      jet_nAK4_100[nfatjets] = 0;
      jet_nAK4_150[nfatjets] = 0;

      double jet_px_20 =0, jet_py_20 =0,jet_pz_20 =0,jet_E_20 =0;
      double jet_px_30 =0, jet_py_30 =0,jet_pz_30 =0,jet_E_30 =0;
      double jet_px_50 =0, jet_py_50 =0,jet_pz_50 =0,jet_E_50 =0;
      double jet_px_70 =0, jet_py_70 =0,jet_pz_70 =0,jet_E_70 =0;
      double jet_px_100 =0, jet_py_100 =0,jet_pz_100 =0,jet_E_100 =0;
      double jet_px_150 =0, jet_py_150 =0,jet_pz_150 =0,jet_E_150 =0;

      for (auto i=jetsFJ_jet.begin(); i<jetsFJ_jet.end(); i++)                             
      {

         if (i->E() > 20.)
         {
            jet_nAK4_20[nfatjets]++;
            jet_px_20+= i->px();jet_py_20+= i->py();jet_pz_20+= i->pz();jet_E_20+= i->E();
         }
         if(i->E() > 30.) 
         {
            jet_nAK4_30[nfatjets]++;
            jet_px_30+= i->px();jet_py_30+= i->py();jet_pz_30+= i->pz();jet_E_30+= i->E();
         }
         if (i->E() > 50.)
         {
            jet_nAK4_50[nfatjets]++;
            jet_px_50+= i->px();jet_py_50+= i->py();jet_pz_50+= i->pz();jet_E_50+= i->E();
         }
         if (i->E() > 70.)
         {
            jet_nAK4_70[nfatjets]++;
            jet_px_70+= i->px();jet_py_70+= i->py();jet_pz_70+= i->pz();jet_E_70+= i->E();
         } 
         if (i->E() > 100.)
         {
            jet_nAK4_100[nfatjets]++;
            jet_px_100+= i->px();jet_py_100+= i->py();jet_pz_100+= i->pz();jet_E_100+= i->E();
         }    
         if (i->E() > 150.)
         {
            jet_nAK4_150[nfatjets]++;
            jet_px_150+= i->px();jet_py_150+= i->py();jet_pz_150+= i->pz();jet_E_150+= i->E();
         }       
         jet_nAK4[nfatjets]++;
         AK4_mass_20[nfatjets] = sqrt(pow(jet_E_20,2)-pow(jet_px_20,2) - pow(jet_py_20,2)-pow(jet_pz_20,2));
         AK4_mass_30[nfatjets] = sqrt(pow(jet_E_30,2)-pow(jet_px_30,2) - pow(jet_py_30,2)-pow(jet_pz_30,2));
         AK4_mass_50[nfatjets] = sqrt(pow(jet_E_50,2)-pow(jet_px_50,2) - pow(jet_py_50,2)-pow(jet_pz_50,2));
         AK4_mass_70[nfatjets] = sqrt(pow(jet_E_70,2)-pow(jet_px_70,2) - pow(jet_py_70,2)-pow(jet_pz_70,2));
         AK4_mass_100[nfatjets] = sqrt(pow(jet_E_100,2)-pow(jet_px_100,2) - pow(jet_py_100,2)-pow(jet_pz_100,2));
         AK4_mass_150[nfatjets] = sqrt(pow(jet_E_150,2)-pow(jet_px_150,2) - pow(jet_py_150,2)-pow(jet_pz_150,2));

      }      

      nfatjets++;
   }
   //nfatjets_minus1 = nfatjets -1;
   eventnum++;
   tree->Fill();
}   
DEFINE_FWK_MODULE(diquarkAnalyzerBR);




