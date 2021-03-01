


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
#include "DataFormats/PatCandidates/interface/Muon.h"
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

class ttbarAnalyzer : public edm::EDAnalyzer 
{
public:
   explicit ttbarAnalyzer(const edm::ParameterSet&);
private:
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   double calc_mag(double px,double py, double pz);
   bool isgoodjet(pat::Jet* iM);
   bool isgoodAK4(pat::Jet* iM);
   bool isgoodmuon(pat::Muon* iM);
   double calc_mpp_beta(pat::Jet* iM);
   const reco::Candidate* parse_chain(const reco::Candidate* cand);
   //edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartToken_; 

   edm::EDGetTokenT<std::vector<pat::Jet>> PUPPI_AK4Token_;
   edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
   edm::EDGetTokenT<std::vector<pat::Jet>> fatJetToken_;
   edm::EDGetTokenT<edm::TriggerResults> triggerBits_;

   TTree * tree;

   int eventnum = 0;
   int nfatjets = 0;
   int raw_nfatjets;
   int nfatjets_minus1;
   int nmuons = 0;
   int nevents = 0;
   int nPass_iso = 0;
   int nPass_hpt = 0;
   int btagged_AK4 = 0;
   int hpt_muons = 0;
   double jet_pt[100], jet_eta[100], jet_mass[100], jet_dr[100], raw_jet_mass[100],raw_jet_pt[100],raw_jet_phi[100], muon_pt[100],bjet_pt[100];
   double jet_beta[100], beta_T[100], AK4_mass_20[100],AK4_mass_30[100],AK4_mass_50[100],AK4_mass_70[100],AK4_mass_100[100],AK4_mass_150[100];
   double b_W_deltar,h_t_deltar;
   int jet_ndaughters[100], jet_nAK4[100],jet_nAK4_20[100],jet_nAK4_30[100],jet_nAK4_50[100],jet_nAK4_70[100],jet_nAK4_100[100],jet_nAK4_150[100];
};


ttbarAnalyzer::ttbarAnalyzer(const edm::ParameterSet& iConfig)
{
   //genPartToken_ = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genPartCollection"));

   muonToken_ =    consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonCollection"));
   PUPPI_AK4Token_ =  consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("PUPPI_AK4Collection"));
   fatJetToken_ =    consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("fatJetCollection"));
   triggerBits_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"));

   edm::Service<TFileService> fs;      

   tree = fs->make<TTree>("tree", "tree");

   tree->Branch("nfatjets", &nfatjets, "nfatjets/I");

   tree->Branch("nevents", &nevents, "nevents/I");
   tree->Branch("nPass_iso", &nPass_iso, "nPass_iso/I");
   tree->Branch("nPass_hpt", &nPass_hpt, "nPass_hpt/I");

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

   tree->Branch("btagged_AK4", &btagged_AK4, "btagged_AK4/I");
   tree->Branch("hpt_muons", &hpt_muons, "hpt_muons/I");


   tree->Branch("muon_pt", muon_pt, "muon_pt[hpt_muons]/D");
   tree->Branch("bjet_pt", bjet_pt, "bjet_pt[btagged_AK4]/D");

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

   //tree->Branch("jet_dr", jet_dr, "jet_dr[nfatjets_minus1]/D");

}


double ttbarAnalyzer::calc_mag(double px,double py, double pz)
{
   return sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
}
const reco::Candidate* ttbarAnalyzer::parse_chain(const reco::Candidate* cand)
{  
   for (unsigned int iii=0; iii<cand->numberOfDaughters(); iii++)
   {
      if(cand->daughter(iii)->pdgId() == cand->pdgId()) return parse_chain(cand->daughter(iii));
   }
   return cand;
}

bool ttbarAnalyzer::isgoodjet(pat::Jet* iM)
{
   //double fatjet_inv_mass = sqrt(pow(iM->energy(),2) - pow(iM->p(),2)); //
   double fatjet_sd_mass  = iM->userFloat("ak8PFJetsPuppiSoftDropMass");

   if( (iM->isPFJet() == false)|| (abs(iM->eta()) > 2.4) || (iM->pt() < 300.) ||(fatjet_sd_mass < 50.0)) return false;
  // if( (iM->isPFJet() == false) || (abs(iM->eta()) > 2.4) || (iM->pt() < 250.)) return false;

   const float NHF = iM->neutralHadronEnergyFraction();
   const float NEMF = iM->neutralEmEnergyFraction();
   const size_t NumConst = iM->numberOfDaughters();
   const float CHF = iM->chargedHadronEnergyFraction();
   const int CHM = iM->chargedMultiplicity();
   if ((NHF>0.9) || (NEMF>0.9) || (NumConst<1) || (CHF<0.) || (CHM<0)) 
      {
         return false;
      }
   else{ return true;}
}
bool ttbarAnalyzer::isgoodAK4(pat::Jet* iM)
{
   //double fatjet_inv_mass = sqrt(pow(iM->energy(),2) - pow(iM->p(),2)); //

   if( (iM->isPFJet() == false)|| (abs(iM->eta()) > 2.4) || (iM->pt() < 50.)) return false;
  // if( (iM->isPFJet() == false) || (abs(iM->eta()) > 2.4) || (iM->pt() < 250.)) return false;

   const float NHF = iM->neutralHadronEnergyFraction();
   const float NEMF = iM->neutralEmEnergyFraction();
   const size_t NumConst = iM->numberOfDaughters();
   const float CHF = iM->chargedHadronEnergyFraction();
   const int CHM = iM->chargedMultiplicity();
   if ((NHF>0.9) || (NEMF>0.9) || (NumConst<1) || (CHF<0.) || (CHM<0)) 
      {
         return false;
      }
   else{ return true;}
}
bool ttbarAnalyzer::isgoodmuon(pat::Muon* iM)
{
  
  if(iM->isMediumMuon()) return true;

  return false;
}
double ttbarAnalyzer::calc_mpp_beta(pat::Jet* iM)
{
   Vector jet_axis(iM->px()/iM->p(),iM->py()/iM->p(),iM->pz()/iM->p());
   TLorentzVector jet_particles(iM->px(),iM->py(),iM->pz(),iM->energy());
   //double beta_cm = iM->p()/iM->energy();
   double min_pp = 99999999999999.;
   double min_boost = 0.;

   for(int iii=0;iii<100000;iii++)
   {
      TLorentzVector jet_particles_ = jet_particles;
      double beta_cand = iii/100000.;
      jet_particles_.Boost(-beta_cand*jet_axis.X(),-beta_cand*jet_axis.Y(),-beta_cand*jet_axis.Z());
      if(abs( ( jet_particles_.Px()*iM->px()+jet_particles_.Py()*iM->py() +jet_particles_.Pz()*iM->py() )/iM->p() ) < min_pp) 
         {
            min_boost = beta_cand; 
            min_pp = abs( ( jet_particles_.Px()*iM->px()+jet_particles_.Py()*iM->py() +jet_particles_.Pz()*iM->py() )/iM->p() ) ;
         }
   }
   //std::cout << "Min parallel momentum is " << min_pp << " beta is " << min_boost << " cm beta is " << beta_cm << std::endl;
   return min_boost;
}

void ttbarAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
  nevents++;

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

  /*
  for (unsigned int i = 0; i < triggerBits->size(); ++i) {
      const std::string name = names.triggerName(i);
      const bool accept = triggerBits->accept(i);
      //const int prescale = triggerPrescales->getPrescaleForIndex(i);
      std::cout << " name " << name<< std::endl;
   }

   return;

   */
////////////////Triggers//////////////////////////
   
   bool pass = false;
   //std::string trigname= "HLT_AK8PFJet550_v";
   std::string trigname  = "HLT_Mu50";
   std::string trigname2 = "HLT_IsoMu27_v";
   for (unsigned int i = 0; i < triggerBits->size(); ++i) 
   {
      const std::string name = names.triggerName(i);
      const bool accept = triggerBits->accept(i);
      if (  ((name.find(trigname) != std::string::npos) &&(accept)) )
        {
          pass =true;
          nPass_hpt++;
        }
      if ( ((name.find(trigname2) != std::string::npos) &&(accept)) )
      {
        pass = true;
        nPass_iso++;
      }
   }
   
   if(!pass) return;

////////////////////Muons////////////////////////////////////////////
 edm::Handle<std::vector<pat::Muon>> muons;
 iEvent.getByToken(muonToken_, muons);
 hpt_muons = 0;
 for(auto iMuon = muons->begin(); iMuon != muons->end(); iMuon++) 
 {
    if(!(isgoodmuon(iMuon->clone())) || (iMuon->pt()<50.))continue;
    muon_pt[hpt_muons] = iMuon->pt(); 
    hpt_muons++;
 }

 //std::cout << "Number of high pt (>50 GeV) muons: " << hpt_muons << std::endl;
 if(hpt_muons!=1)return;

   nfatjets = 0;
   raw_nfatjets = 0;

////////////////AK4 Jets//////////////////////////////////////
   edm::Handle<std::vector<pat::Jet> > PUPPI_AK4Jets;
   iEvent.getByToken(PUPPI_AK4Token_, PUPPI_AK4Jets);
   btagged_AK4 = 0;
   for(auto iJet = PUPPI_AK4Jets->begin(); iJet != PUPPI_AK4Jets->end(); iJet++)
   {
      double bdisc = iJet->bDiscriminator("pfDeepCSVJetTags:probb") + iJet->bDiscriminator("pfDeepCSVJetTags:probbb");
      if(   !(isgoodAK4(iJet->clone())) || (bdisc<0.4941)  ) continue;
      bjet_pt[btagged_AK4] = iJet->pt();
      btagged_AK4++;
   }         

   //std::cout << "Number of b jets: "  << btagged_AK4 << std::endl;
   if (btagged_AK4!=1)return;



////////////////Jets//////////////////////////////////////
   edm::Handle<std::vector<pat::Jet> > fatJets;
   iEvent.getByToken(fatJetToken_, fatJets);


   for(auto iJet = fatJets->begin(); iJet != fatJets->end(); iJet++)         ////////Over AK8 Jets
   { 
      

      if(nfatjets>2)continue;
      raw_jet_pt[raw_nfatjets]   = iJet->pt();
      raw_jet_phi[raw_nfatjets]  = iJet->phi();
      raw_jet_mass[raw_nfatjets] = sqrt(pow(iJet->energy(),2) - pow(iJet->p(),2));
      raw_nfatjets++;

      //if(!isgoodjet(iJet->clone())) continue;
      jet_pt[nfatjets]         = iJet->pt();
      jet_eta[nfatjets]        = iJet->eta();
      jet_mass[nfatjets]       = sqrt(pow(iJet->energy(),2) - pow(iJet->p(),2));
      jet_ndaughters[nfatjets] = iJet->numberOfDaughters();
      double beta_mag = calc_mpp_beta(iJet->clone());
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
      std::vector<fastjet::PseudoJet> jetsFJ_jet = fastjet::sorted_by_E(cs_jet.inclusive_jets(0.0));
      
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
   if(nfatjets < 1) return;
   eventnum++;
   tree->Fill();
}   
DEFINE_FWK_MODULE(ttbarAnalyzer);




