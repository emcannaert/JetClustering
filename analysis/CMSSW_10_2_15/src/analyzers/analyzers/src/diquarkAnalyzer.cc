////////////////////////////HELP////////////////////////////////
//////////////Analyzes the diquark -> chi chi hadronic decay//////////////
////////////////Last updated Nov 9 2020/////////////////////////////////


// system include files
#include <fastjet/JetDefinition.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/tools/Filter.hh>
#include <fastjet/ClusterSequence.hh>
//#include <fastjet/ActiveAreaSpec.hh>
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

class diquarkAnalyzer : public edm::EDAnalyzer 
{
public:
   explicit diquarkAnalyzer(const edm::ParameterSet&);
private:
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   double calc_mag(double px,double py, double pz);
   bool isgoodjet(const float pt, const float fatjet_sd_mass, const float eta, const float NHF,const float NEMF, const size_t NumConst,const float CHF,const int CHM, const float MUF, const float CEMF, const int NumNeutralParticles);
   const reco::Candidate* parse_chain(const reco::Candidate* cand);
   edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartToken_; 
   edm::EDGetTokenT<std::vector<pat::Jet>> fatJetToken_;
   edm::EDGetTokenT<edm::TriggerResults> triggerBits_;

   TTree * tree;
   std::vector<std::string> triggs = {"HLT_AK8PFJet800"};

   int eventnum = 0;
   int nhadevents = 0;
   int nfatjets = 0;
   int raw_nfatjets;
   int nfatjets_minus1;
   double jet_pt[100], jet_eta[100], jet_mass[100], jet_dr[100], raw_jet_mass[100],raw_jet_pt[100],raw_jet_phi[100];
   double jet_beta[100], beta_T[100], AK4_mass_20[100],AK4_mass_30[100],AK4_mass_50[100],AK4_mass_70[100],AK4_mass_100[100],AK4_mass_150[100];
   double b_W_deltar,h_t_deltar;
   int nChi = 0;
   int nSuu = 0;
   int nq = 0;
   double tot_jet_mass,decay_inv_mass, chi_inv_mass;
   double tot_HT = 0.;


   int jet_ndaughters[100], jet_nAK4[100],jet_nAK4_20[100],jet_nAK4_30[100],jet_nAK4_50[100],jet_nAK4_70[100],jet_nAK4_100[100],jet_nAK4_150[100];
};


diquarkAnalyzer::diquarkAnalyzer(const edm::ParameterSet& iConfig)
{
   genPartToken_ = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genPartCollection"));
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

   tree->Branch("b_W_deltar",&b_W_deltar,"b_W_deltar/D");
   tree->Branch("h_t_deltar",&h_t_deltar,"h_t_deltar/D");

   tree->Branch("nChi",&nChi,"nChi/I");
   tree->Branch("nSuu",&nSuu,"nSuu/I");
   tree->Branch("nq",&nq,"nq/I");
   tree->Branch("decay_inv_mass", &decay_inv_mass, "decay_inv_mass/D");
   tree->Branch("chi_inv_mass", &chi_inv_mass, "chi_inv_mass/D");

   tree->Branch("tot_HT",&tot_HT,"tot_HT/D");
   tree->Branch("tot_jet_mass",&tot_jet_mass, "tot_jet_mass/D" );
   //tree->Branch("jet_dr", jet_dr, "jet_dr[nfatjets_minus1]/D");

}


double diquarkAnalyzer::calc_mag(double px,double py, double pz)
{
   return sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
}
const reco::Candidate* diquarkAnalyzer::parse_chain(const reco::Candidate* cand)
{  
   for (unsigned int iii=0; iii<cand->numberOfDaughters(); iii++)
   {
      if(cand->daughter(iii)->pdgId() == cand->pdgId()) return parse_chain(cand->daughter(iii));
   }
   return cand;
}

bool diquarkAnalyzer::isgoodjet(const float pt, const float fatjet_sd_mass, const float eta, const float NHF,const float NEMF, const size_t NumConst,const float CHF,const int CHM, const float MUF, const float CEMF, const int NumNeutralParticles)
{
   if( (abs(eta) > 2.4) || (pt < 500.) ||(fatjet_sd_mass < 50.0)) return false;

   if ((NHF>0.9) || (NEMF>0.9) || (NumConst<1) || (CHF<0.) || (CHM<0) || (MUF > 0.8) || (CEMF > 0.8)) 
      {
         return false;
      }
   else{ return true;}

}

void diquarkAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   eventnum++;


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

////////////Gen Particles//////////////////////////////////
/*
   edm::Handle<std::vector<reco::GenParticle>> genParticles;
   iEvent.getByToken(genPartToken_, genParticles);

   int suu_pdgid = 9936661;
   int chi_pdgid = 9936662;
   int nH = 0;
   nq = 0;
   int ntop = 0;
   int nb   = 0;
   int nW   = 0;
   nChi = 0;
   nSuu = 0;
   double top_eta = 0;
   double top_phi = 0;
   double higgs_eta = 0;
   double higgs_phi = 0;
   double b_eta = 0;
   double b_phi = 0;
   double W_eta = 0;
   double W_phi = 0;

   double decay_px =0, decay_py=0,decay_pz = 0, decay_E= 0;
   double chi_px =0 ,chi_py=0,chi_pz=0,chi_E=0;
std::cout << "//////////////////////////////////" << std::endl;

//   double b_W_deltar[100],h_t_deltar[100], min_qqb_deltar[100], min_bbqqb_deltar[100];
//want full hadronic case
   for (auto iG = genParticles->begin(); iG != genParticles->end(); iG++) 
   {
      //if((iG->mass() > 150.) && (iG->isLastCopy())) std::cout << "Mass is " << iG->mass() << " pdgID is " << iG->pdgId() << " Parent is " << iG->mother()->pdgId() <<  std::endl;  
      if((iG->mass() > 75.) && (iG->isLastCopy())) std::cout << "Mass is " << iG->mass() << " pdgID is " << iG->pdgId() <<  " coming from "  << iG->mother()->pdgId() <<  std::endl;  

      if ((abs(iG->pdgId()) == 24) && ((abs(iG->mother()->pdgId()) == chi_pdgid) || (abs(iG->mother()->pdgId()) == 2)) ) 
      {
         decay_px+=iG->px();
         decay_py+=iG->py();
         decay_pz+=iG->pz();
         decay_E +=iG->energy();

         const reco::Candidate* W_final = parse_chain(iG->clone());
         W_eta = W_final->eta();
         W_phi = W_final->phi();


         for (unsigned int iii=0; iii<W_final->numberOfDaughters(); iii++)
         {
            const reco::Candidate* W_daughter = W_final->daughter(iii);
            if (abs(W_daughter->pdgId())<6) nq++;

         }
         nW++;
      }
      else if ( (abs(iG->pdgId()) == 5) &&  ((abs(iG->mother()->pdgId()) == chi_pdgid) || (abs(iG->mother()->pdgId()) == 2))  )
      {
         b_eta = iG->eta();
         b_phi = iG->phi();

         decay_px+=iG->px();
         decay_py+=iG->py();
         decay_pz+=iG->pz();
         decay_E +=iG->energy();

         nq++;
         //std::cout << "Found the first 'guaranteed' quark directly from the chi " << std::endl;
      } 

      else if ( (abs(iG->pdgId()) == 6) && ((abs(iG->mother()->pdgId()) == chi_pdgid)|| (abs(iG->mother()->pdgId()) == 2)) ) 
      {
         const reco::Candidate* t_final = parse_chain(iG->clone());
         top_eta = t_final->eta();
         top_phi = t_final->phi();


         decay_px+=iG->px();
         decay_py+=iG->py();
         decay_pz+=iG->pz();
         decay_E +=iG->energy();

         for (unsigned int iii=0; iii<t_final->numberOfDaughters(); iii++)
         {
            const reco::Candidate* t_daughter = t_final->daughter(iii);
            if (abs(t_daughter->pdgId())==24) 
            {
               const reco::Candidate* W_final = parse_chain(t_daughter->clone());
               for (unsigned int jjj=0; jjj<W_final->numberOfDaughters(); jjj++)
               {
                  const reco::Candidate* W_daughter = W_final->daughter(jjj);
                  if(abs(W_daughter->pdgId()) < 6) nq++;
               }
               nW++;
            }
            else if(abs(t_daughter->pdgId())==5) 
            {
               //std::cout << "Found the second 'guaranteed' quark from the t " << std::endl;
               nq++;
            }
         }
         ntop++;
      }

      else if ( (abs(iG->pdgId()) == 25) && ((abs(iG->mother()->pdgId()) == chi_pdgid) || (abs(iG->mother()->pdgId()) == 2)) ) 
      {
         const reco::Candidate* h_final = parse_chain(iG->clone());


         decay_px+=iG->px();
         decay_py+=iG->py();
         decay_pz+=iG->pz();
         decay_E +=iG->energy();

         higgs_eta = h_final->eta();
         higgs_phi = h_final->phi();

         for (unsigned int iii=0; iii<h_final->numberOfDaughters(); iii++)
         {
            const reco::Candidate* h_daughter = h_final->daughter(iii);
            if (abs(h_daughter->pdgId())<6) nq++;
         }
         nH++;
      }
      else if ((abs(iG->pdgId()) == chi_pdgid) && (iG->isLastCopy()))
      {
         chi_px+=iG->px();
         chi_py+=iG->py();
         chi_pz+=iG->pz();
         chi_E +=iG->energy();
         nChi++;
      } 
      else if ((abs(iG->pdgId()) == suu_pdgid) && (iG->isLastCopy())) nSuu++;
   }

   decay_inv_mass= sqrt(pow(decay_E ,2)-pow(decay_px,2)-pow(decay_py,2)-pow(decay_pz,2));
   chi_inv_mass   =sqrt(pow(chi_E ,2)-pow(chi_px,2)-pow(chi_py,2)-pow(chi_pz,2));
   std::cout << "Chi inv mass is " << decay_inv_mass << std::endl;
   bool no_higgs     = false;
   bool single_higgs = false;
   bool double_higgs = false;
   if     ((nH == 0) && (nW == 2) && (nq == 6)  && (ntop == 0)) no_higgs     = true;
   if((nH == 1) && (nW == 1) && (nq == 8)  && (ntop == 1)) single_higgs = true;
   else if((nH == 2) && (nW == 0) && (nq == 10) && (ntop == 2)) double_higgs = true;

   if         (no_higgs) std::cout << "No Higgs event"     << std::endl;
   else if(single_higgs) std::cout << "Single Higgs event" << std::endl;
   else if(double_higgs) std::cout << "Double Higgs event" << std::endl;
   else{ std::cout << "nHiggs/nW/nq/ntop: " << nH << "/" << nW<< "/"<< nq<< "/"<< ntop << std::endl;}

   //std::cout << "Number of quarks is " << nq << std::endl;
   std::cout << "nSuu/nChi/nHiggs/nW/nq/ntop: "<< nSuu << "/" << nChi << "/" << nH << "/" << nW<< "/"<< nq<< "/"<< ntop << std::endl;
   //if(nq < 8) return;
   nhadevents++;

   h_t_deltar = sqrt(pow(higgs_eta-top_eta,2)+pow(higgs_phi-top_phi,2));
   b_W_deltar = sqrt(pow(b_eta-W_eta,2)+pow(b_phi-W_phi,2));
*/
   //if(!(single_higgs)) return;
   
   edm::Handle<std::vector<pat::Jet> > fatJets;
   iEvent.getByToken(fatJetToken_, fatJets);

   std::vector<fastjet::PseudoJet> allAK8_part_;

   for(auto iJet = fatJets->begin(); iJet != fatJets->end(); iJet++)         ////////Over AK8 Jets
   {
      for (unsigned int i=0; i<iJet->numberOfDaughters();i++)
      {
         const reco::Candidate* iJ = iJet->daughter(i);
         const pat::PackedCandidate* cand_begin = (pat::PackedCandidate*) iJ;
         double puppiweight = cand_begin->puppiWeight();
         allAK8_part_.push_back(fastjet::PseudoJet(puppiweight*iJ->px(),puppiweight*iJ->py(),puppiweight*iJ->pz(),puppiweight*iJ->energy()));
      }
   }
   double R = 1.2;
   fastjet::JetDefinition huge_jet_def(fastjet::antikt_algorithm, R);
   fastjet::ClusterSequence cs_huge_jet(allAK8_part_, huge_jet_def); 
   std::vector<fastjet::PseudoJet> jetsFJ_hugejet = fastjet::sorted_by_E(cs_huge_jet.inclusive_jets(10.0));
   int nhugeJet = 0;
   for (auto iJet=jetsFJ_hugejet.begin(); iJet<jetsFJ_hugejet.end(); iJet++)                             
   {
      if(iJet->m()>400.0)nhugeJet++;
      std::cout << "Huge jet number " << nhugeJet << ", mass is " << iJet->m() << std::endl;
   }
   std::cout << "//////////////////////////////" << std::endl;
   nfatjets = 0;
   raw_nfatjets = 0;
   tot_HT = 0;
////////////////Jets//////////////////////////////////////

   double tot_jet_px =0;
   double tot_jet_py =0;
   double tot_jet_pz =0;
   double tot_jet_E =0;

   for(auto iJet = fatJets->begin(); iJet != fatJets->end(); iJet++)         ////////Over AK8 Jets
   { 
      if( (abs(iJet->eta())<2.5)&& (iJet->pt()>30.)) tot_HT+=abs(iJet->pt());
      if(sqrt(pow(iJet->energy(),2) - pow(iJet->p(),2)) < 30000)
         {
            raw_jet_pt[raw_nfatjets]   = iJet->pt();
            raw_jet_phi[raw_nfatjets]  = iJet->phi();
            raw_jet_mass[raw_nfatjets] = sqrt(pow(iJet->energy(),2) - pow(iJet->p(),2));
            raw_nfatjets++;
         }
      if(!(iJet->isPFJet())) continue;

      
      if(!isgoodjet(iJet->pt(),iJet->userFloat("ak8PFJetsPuppiSoftDropMass"), iJet->eta(),iJet->neutralHadronEnergyFraction(), iJet->neutralEmEnergyFraction(),iJet->numberOfDaughters(),iJet->chargedHadronEnergyFraction(),iJet->chargedMultiplicity(),iJet->muonEnergyFraction(),iJet->chargedEmEnergyFraction(),iJet->neutralMultiplicity())) continue;
      tot_jet_px+= iJet->px();
      tot_jet_py+= iJet->py();
      tot_jet_pz+= iJet->pz();
      tot_jet_E+= iJet->energy();



      jet_pt[nfatjets]         = iJet->pt();
      jet_eta[nfatjets]        = iJet->eta();
      jet_mass[nfatjets]       = sqrt(pow(iJet->energy(),2) - pow(iJet->p(),2));
      jet_ndaughters[nfatjets] = iJet->numberOfDaughters();





      ///////////////////////////////////calc beta
      Vector jet_axis(iJet->px()/iJet->p(),iJet->py()/iJet->p(),iJet->pz()/iJet->p());
      TLorentzVector jet_particles(iJet->px(),iJet->py(),iJet->pz(),iJet->energy());
      double min_pp = 99999999999999.;
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

      /////////////////////////////////////////////
      double beta_mag = min_boost;

      jet_beta[nfatjets] = beta_mag;
      beta_T[nfatjets] = beta_mag*sin(iJet->theta());
      ///Recluster//
      //double beta_mag = iJet->p()/iJet->energy();
      Vector beta(beta_mag*iJet->px()/iJet->p(),beta_mag*iJet->py()/iJet->p(),beta_mag*iJet->pz()/iJet->p());

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
   tot_jet_mass = sqrt(pow(tot_jet_E,2)-pow(tot_jet_px,2)-pow(tot_jet_py,2)-pow(tot_jet_pz,2));
   //std::cout << "Tot jet mass is " << tot_jet_mass << std::endl;
   tree->Fill();
}   
DEFINE_FWK_MODULE(diquarkAnalyzer);




