////////////////////////////HELP////////////////////////////////
//////////////Uses new clustering algorithm to capture heavy resonance jet substructure//////////////
////////////////Last updated Feb 23 2021 ////////////////////////////////////////////////////////////


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
#include "Thrust.h"
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

class clusteringAnalyzer : public edm::EDAnalyzer 
{
public:
   explicit clusteringAnalyzer(const edm::ParameterSet&);
private:
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   double calc_mag(double px,double py, double pz);
   double calcMPP(TLorentzVector superJetTLV ); 
   bool isgoodjet(const float eta, const float NHF,const float NEMF, const size_t NumConst,const float CHF,const int CHM, const float MUF, const float CEMF);
   const reco::Candidate* parse_chain(const reco::Candidate* cand);
   edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartToken_; 
   edm::EDGetTokenT<std::vector<pat::Jet>> fatJetToken_;
   edm::EDGetTokenT<edm::TriggerResults> triggerBits_;

   TTree * tree;
   std::vector<std::string> triggs = {"HLT_AK8PFJet800"};

   int eventnum = 0;
   int nfatjets = 0;
   int raw_nfatjets;
   int tot_nAK4_50,tot_nAK4_70;
   int SJ_nAK4_50[100],SJ_nAK4_70[100];
   double jet_pt[100], jet_eta[100], jet_mass[100], jet_dr[100], raw_jet_mass[100],raw_jet_pt[100],raw_jet_phi[100];
   double jet_beta[100], beta_T[100], AK4_mass_20[100],AK4_mass_30[100],AK4_mass_50[100],AK4_mass_70[100],AK4_mass_100[100],AK4_mass_150[100];
   double SJ_mass_50[100], SJ_mass_70[100],superJet_mass[100],SJ_AK4_50_mass[100],SJ_AK4_70_mass[100];
   double tot_jet_mass,decay_inv_mass, chi_inv_mass;
   int nSuperJets;
   int nhadevents = 0;
   int nBadClusters =0;
   int jet_ndaughters[100], jet_nAK4[100],jet_nAK4_20[100],jet_nAK4_30[100],jet_nAK4_50[100],jet_nAK4_70[100],jet_nAK4_100[100],jet_nAK4_150[100];
};


clusteringAnalyzer::clusteringAnalyzer(const edm::ParameterSet& iConfig)
{
   genPartToken_ = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genPartCollection"));
   fatJetToken_ =    consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("fatJetCollection"));
   triggerBits_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"));

   edm::Service<TFileService> fs;      

   tree = fs->make<TTree>("tree", "tree");

   tree->Branch("nfatjets", &nfatjets, "nfatjets/I");
   tree->Branch("nSuperJets", &nSuperJets, "nSuperJets/I");
   tree->Branch("tot_nAK4_50", &tot_nAK4_50, "tot_nAK4_50/I");             //total #AK4 jets (E>50 GeV) for BOTH superjets
   tree->Branch("tot_nAK4_70", &tot_nAK4_70, "tot_nAK4_70/I");

   tree->Branch("jet_pt", jet_pt, "jet_pt[nfatjets]/D");
   tree->Branch("jet_eta", jet_eta, "jet_eta[nfatjets]/D");
   tree->Branch("jet_mass", jet_mass, "jet_mass[nfatjets]/D");

   tree->Branch("SJ_nAK4_50", SJ_nAK4_50, "SJ_nAK4_50[nSuperJets]/I");
   tree->Branch("SJ_nAK4_70", SJ_nAK4_70, "SJ_nAK4_70[nSuperJets]/I");
   tree->Branch("SJ_mass_50", SJ_mass_50, "SJ_mass_50[nSuperJets]/D");
   tree->Branch("SJ_mass_70", SJ_mass_70, "SJ_mass_70[nSuperJets]/D");
   tree->Branch("superJet_mass", superJet_mass, "superJet_mass[nSuperJets]/D");
   tree->Branch("SJ_AK4_50_mass", SJ_AK4_50_mass, "SJ_AK4_50_mass[tot_nAK4_50]/D");    //mass of individual reclustered AK4 jets
   tree->Branch("SJ_AK4_70_mass", SJ_AK4_70_mass, "SJ_AK4_70_mass[tot_nAK4_70]/D");

}


double clusteringAnalyzer::calc_mag(double px,double py, double pz)
{
   return sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
}
const reco::Candidate* clusteringAnalyzer::parse_chain(const reco::Candidate* cand)
{  
   for (unsigned int iii=0; iii<cand->numberOfDaughters(); iii++)
   {
      if(cand->daughter(iii)->pdgId() == cand->pdgId()) return parse_chain(cand->daughter(iii));
   }
   return cand;
}

bool clusteringAnalyzer::isgoodjet(const float eta, const float NHF,const float NEMF, const size_t NumConst,const float CHF,const int CHM, const float MUF, const float CEMF)
{
   if( (abs(eta) > 2.4)) return false;

   if ((NHF>0.9) || (NEMF>0.9) || (NumConst<1) || (CHF<0.) || (CHM<0) || (MUF > 0.8) || (CEMF > 0.8)) 
      {
         return false;
      }
   else{ return true;}

}

double clusteringAnalyzer::calcMPP(TLorentzVector superJetTLV ) 
{
   Vector jet_axis(superJetTLV.Px()/superJetTLV.P(),superJetTLV.Py()/superJetTLV.P(),superJetTLV.Pz()/superJetTLV.P());
   double min_pp = 99999999999999.;
   double min_boost = 0.;

   for(int iii=0;iii<10000;iii++)
   {
      TLorentzVector superJetTLV_ = superJetTLV;
      double beta_cand = iii/10000.;
      superJetTLV_.Boost(-beta_cand*jet_axis.X(),-beta_cand*jet_axis.Y(),-beta_cand*jet_axis.Z());
      if(abs( ( superJetTLV_.Px()*superJetTLV.Px()+superJetTLV_.Py()*superJetTLV.Py() +superJetTLV_.Pz()*superJetTLV.Py() )/superJetTLV.P() ) < min_pp) 
         {
            min_boost = beta_cand; 
            min_pp = abs( ( superJetTLV_.Px()*superJetTLV.Px()+superJetTLV_.Py()*superJetTLV.Py() +superJetTLV_.Pz()*superJetTLV.Py() )/superJetTLV.P() ) ;
         }
      }
      return min_boost;
}

void clusteringAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   eventnum++;
////////////////Triggers//////////////////////////
   
   edm::Handle<edm::TriggerResults> triggerBits;
   iEvent.getByToken(triggerBits_, triggerBits);
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

   bool pass = false;
   std::string trigname= "HLT_PFHT1050_v";
   for (unsigned int i = 0; i < triggerBits->size(); ++i) 
   {
      const std::string name = names.triggerName(i);
      const bool accept = triggerBits->accept(i);
      if ((name.find(trigname) != std::string::npos) &&(accept)) pass =true;
   }
   if(!pass) return;
//////////////////////////////////////////////////////////////////////



////////////Gen Particles//////////////////////////////////
//////////////////////////////////////////////////////////
   edm::Handle<std::vector<reco::GenParticle>> genParticles;
   iEvent.getByToken(genPartToken_, genParticles);

   int suu_pdgid = 9936661;
   int chi_pdgid = 9936662;
   int nSuu = 0;
   int nChi = 0;
   int nq   = 0;
   int nW   = 0;
   int ntop = 0;
   int nH   = 0;
   TLorentzVector chi1(0,0,0,0);
   TLorentzVector chi2(0,0,0,0);
   TLorentzVector top(0,0,0,0);
   TLorentzVector W_Suu(0,0,0,0);
   TLorentzVector b_Suu(0,0,0,0);
   TLorentzVector higgs(0,0,0,0);
//   double b_W_deltar[100],h_t_deltar[100], min_qqb_deltar[100], min_bbqqb_deltar[100];
//want full hadronic case
   for (auto iG = genParticles->begin(); iG != genParticles->end(); iG++) 
   {
      if ((abs(iG->pdgId()) == 24) && ((abs(iG->mother()->pdgId()) == chi_pdgid)) ) 
      {
         const reco::Candidate* W_final = parse_chain(iG->clone());
         for (unsigned int iii=0; iii<W_final->numberOfDaughters(); iii++)
         {
            const reco::Candidate* W_daughter = W_final->daughter(iii);
            if (abs(W_daughter->pdgId())<6) nq++;
         }
         W_Suu.SetPxPyPzE(W_final->px(),W_final->py(),W_final->pz(),W_final->energy());
         nW++;
      }
      else if ( (abs(iG->pdgId()) == 5) &&  ((abs(iG->mother()->pdgId()) == chi_pdgid))  )
      {
         b_Suu.SetPxPyPzE(iG->px(),iG->py(),iG->pz(),iG->energy());
         nq++;
      } 

      else if ( (abs(iG->pdgId()) == 6) && ((abs(iG->mother()->pdgId()) == chi_pdgid)) ) 
      {
         const reco::Candidate* t_final = parse_chain(iG->clone());
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
               nq++;
            }
         }
         top.SetPxPyPzE(iG->px(),iG->py(),iG->pz(),iG->energy());
         ntop++;
      }

      else if ( (abs(iG->pdgId()) == 25) && ((abs(iG->mother()->pdgId()) == chi_pdgid)) ) 
      {
         const reco::Candidate* h_final = parse_chain(iG->clone());
         for (unsigned int iii=0; iii<h_final->numberOfDaughters(); iii++)
         {
            const reco::Candidate* h_daughter = h_final->daughter(iii);
            if (abs(h_daughter->pdgId())<6) nq++;
         }

         higgs.SetPxPyPzE(iG->px(),iG->py(),iG->pz(),iG->energy());
         nH++;
      }
      else if ((abs(iG->pdgId()) == chi_pdgid) && (iG->isLastCopy()))
      {
         if      (nChi == 0) chi1.SetPxPyPzE(iG->px(),iG->py(),iG->pz(),iG->energy());
         else if (nChi == 1) chi2.SetPxPyPzE(iG->px(),iG->py(),iG->pz(),iG->energy());
         nChi++;
      } 
      else if ((abs(iG->pdgId()) == suu_pdgid) && (iG->isLastCopy())) nSuu++;
   
      //else if  ((abs(iG->pdgId()) == chi_pdgid) && (iG->isLastCopy())) nChi++;
   }
   //std::cout << "nH/nt/nW/nq" << nH << "/" << ntop<< "/" << nW<< "/" <<nq << std::endl;
   if( !((nH == 1) && (nW == 2) && (nq > 7)  && (ntop == 1))) return;
   nhadevents++;




   std::vector<reco::LeafCandidate> candsBoosted;   
   std::vector<reco::LeafCandidate> candsUnboosted;
   std::vector<fastjet::PseudoJet> superJetOne;        //jets for dot prduct #cos(theta) > 0
   std::vector<fastjet::PseudoJet> superJetTwo;       //jets for dot product #cos(theta) < 0

   edm::Handle<std::vector<pat::Jet> > fatJets;
   iEvent.getByToken(fatJetToken_, fatJets);


   //Get AK8 jet info in order to get to COM frame and get a vector of all the (good) AK8 jet particles 
   nfatjets = 0;
   double tot_jet_px=0,tot_jet_py=0,tot_jet_pz=0, tot_jet_E=0;
   for(auto iJet = fatJets->begin(); iJet != fatJets->end(); iJet++)         ////////Over AK8 Jets
   {
      if((sqrt(pow(iJet->mass(),2)+pow(iJet->pt(),2)) < 100.) || (!(iJet->isPFJet())) || (!isgoodjet(iJet->eta(),iJet->neutralHadronEnergyFraction(), iJet->neutralEmEnergyFraction(),iJet->numberOfDaughters(),iJet->chargedHadronEnergyFraction(),iJet->chargedMultiplicity(),iJet->muonEnergyFraction(),iJet->chargedEmEnergyFraction())) || (iJet->userFloat("ak8PFJetsPuppiSoftDropMass") < 15.) ) continue;
      jet_pt[nfatjets] = iJet->pt();
      jet_eta[nfatjets] = iJet->eta();
      jet_mass[nfatjets] = iJet->mass();
      for (unsigned int iii=0; iii<iJet->numberOfDaughters();iii++)   ///////get all jet particles
      {
         const reco::Candidate* iJ = iJet->daughter(iii);
         const pat::PackedCandidate* candJetbegin = (pat::PackedCandidate*) iJ;
         double puppiweight = candJetbegin->puppiWeight();
         candsUnboosted.push_back(LeafCandidate(iJet->daughter(iii)->charge(), Particle::LorentzVector(puppiweight*candJetbegin->px(), puppiweight*candJetbegin->py(), puppiweight*candJetbegin->pz(), puppiweight*candJetbegin->energy())));
      }
      tot_jet_px+=iJet->px();tot_jet_py+=iJet->py();tot_jet_pz+=iJet->pz();tot_jet_E+=iJet->energy();
      nfatjets++;
   }
   if(nfatjets<2)return;


   ///////////calculate COM of all jets
   double tot_jet_p = sqrt(pow(tot_jet_px,2)+pow(tot_jet_py,2)+pow(tot_jet_pz,2));            ///is this right?///
   double tot_jet_beta = tot_jet_p/tot_jet_E;
   TVector3 totJetBeta(tot_jet_beta*tot_jet_px/tot_jet_p,tot_jet_beta*tot_jet_py/tot_jet_p,tot_jet_beta*tot_jet_pz/tot_jet_p);

   //////////boost all jet particles into COM frame
   int nSuperJetConst = 0;
   for(auto iC = candsUnboosted.begin();iC != candsUnboosted.end(); iC++)
   {
      TLorentzVector iC_(iC->px(),iC->py(),iC->pz(),iC->energy());
      iC_.Boost(-totJetBeta.X(),-totJetBeta.Y(),-totJetBeta.Z());
      candsBoosted.push_back(LeafCandidate(iC->charge(), Particle::LorentzVector(iC_.Px(), iC_.Py(), iC_.Pz(), iC_.E())));
      nSuperJetConst++;
   }
   Thrust thrust_(candsBoosted, nSuperJetConst);                             //thrust axis in COM frame
   Vector thrustAxis = thrust_.axis();
   TVector3 thrust_vector(thrustAxis.X(),thrustAxis.Y(),thrustAxis.Z());

   ////////check out chi's and see if they are in the same superjet
   chi1.Boost(-totJetBeta.X(),-totJetBeta.Y(),-totJetBeta.Z());
   chi2.Boost(-totJetBeta.X(),-totJetBeta.Y(),-totJetBeta.Z());
   top.Boost(-totJetBeta.X(),-totJetBeta.Y(),-totJetBeta.Z());
   W_Suu.Boost(-totJetBeta.X(),-totJetBeta.Y(),-totJetBeta.Z());
   b_Suu.Boost(-totJetBeta.X(),-totJetBeta.Y(),-totJetBeta.Z());
   higgs.Boost(-totJetBeta.X(),-totJetBeta.Y(),-totJetBeta.Z());

   TVector3 chi1_vec = chi1.Vect();
   TVector3 chi2_vec = chi2.Vect();
   TVector3 top_vec = top.Vect();
   TVector3 W_Suu_vec = W_Suu.Vect();
   TVector3 b_Suu_vec = b_Suu.Vect();
   TVector3 higgs_vec = higgs.Vect();

   double cosAngleChi1 = cos(chi1_vec.Angle(thrust_vector));
   double cosAngleChi2 = cos(chi2_vec.Angle(thrust_vector));
   double cosAngleTop = cos(top_vec.Angle(thrust_vector));
   double cosAngleW_Suu = cos(W_Suu_vec.Angle(thrust_vector));
   double cosAngleB_Suu = cos(b_Suu_vec.Angle(thrust_vector));
   double cosAngleHiggs = cos(higgs_vec.Angle(thrust_vector));


   if((cosAngleChi1*cosAngleChi2 > 0))
   {
      nBadClusters++;
      //std::cout << "Gen Chis are in same SuperJet - " << nBadClusters << " out of " << eventnum << " events." << std::endl;
   }

   for(auto iJet = fatJets->begin(); iJet != fatJets->end(); iJet++)         //get vector of all sorted superjet particles
   {
      if((sqrt(pow(iJet->mass(),2)+pow(iJet->pt(),2)) < 100.) || (!(iJet->isPFJet())) || (!isgoodjet(iJet->eta(),iJet->neutralHadronEnergyFraction(), iJet->neutralEmEnergyFraction(),iJet->numberOfDaughters(),iJet->chargedHadronEnergyFraction(),iJet->chargedMultiplicity(),iJet->muonEnergyFraction(),iJet->chargedEmEnergyFraction())) || (iJet->userFloat("ak8PFJetsPuppiSoftDropMass") < 15.) ) continue;
      TLorentzVector candJet(iJet->px(),iJet->py(),iJet->pz(),iJet->energy());
      candJet.Boost(-totJetBeta.X(),-totJetBeta.Y(),-totJetBeta.Z());         //boost jet into COM frame

      TVector3 candJet_vec = candJet.Vect();
      double cosAngle = cos(candJet_vec.Angle(thrust_vector));
      //double dotProduct = candJet.X()*thrustAxis.X() + candJet.Y()*thrustAxis.Y() + candJet.Z()*thrustAxis.Z();           //look at dot product of jet in COM frame and thrust axis
      //double magProduct = sqrt(pow(candJet.X(),2)+pow(candJet.Y(),2)+pow(candJet.Z(),2))*sqrt(pow(thrustAxis.X(),2)+pow(thrustAxis.Y(),2)+pow(thrustAxis.Z(),2));
      //double cosTheta   = dotProduct/magProduct;
      

      //sort jet particles into either SuperJet 1 or SuperJet 2
      for (unsigned int i=0; i<iJet->numberOfDaughters();i++)   
      {
         const reco::Candidate* iJ = iJet->daughter(i);
         const pat::PackedCandidate* candJetbegin = (pat::PackedCandidate*) iJ;
         double puppiweight = candJetbegin->puppiWeight();

         if     (cosAngle>0)superJetOne.push_back(fastjet::PseudoJet(puppiweight*iJ->px(),puppiweight*iJ->py(),puppiweight*iJ->pz(),puppiweight*iJ->energy()));
         else if(cosAngle<0)superJetTwo.push_back(fastjet::PseudoJet(puppiweight*iJ->px(),puppiweight*iJ->py(),puppiweight*iJ->pz(),puppiweight*iJ->energy()));
      }
   }
   double superJetpx,superJetpy,superJetpz,superJetE;
   std::vector<std::vector<fastjet::PseudoJet>> superJets;
   superJets.push_back(superJetOne);
   superJets.push_back(superJetTwo);

   nSuperJets = 0;
   tot_nAK4_70 = 0; tot_nAK4_50 = 0;

   //Find MPP frame and boost particles into this frame
   for(auto iSJ = superJets.begin();iSJ!= superJets.end();iSJ++)
   {
      superJetpx=0;superJetpy=0;superJetpz=0;superJetE=0;
      for(auto iP = iSJ->begin();iP != iSJ->end();iP++)
      {
         superJetpx+=iP->px();
         superJetpy+=iP->py();
         superJetpz+=iP->pz();
         superJetE +=iP->E();
      }
      superJet_mass[nSuperJets] = sqrt(pow(superJetE,2)-pow(superJetpx,2)-pow(superJetpy,2)-pow(superJetpz,2));


      TLorentzVector superJetTLV(superJetpx,superJetpy,superJetpz,superJetE);    //Lorentz vector representing jet axis -> now minimize the parallel momentum
      double betaMag = calcMPP(superJetTLV);
      std::vector<fastjet::PseudoJet> boostedSuperJetPart;

      //boost particles in SuperJet to MPP frame
      for(auto iP = iSJ->begin();iP != iSJ->end();iP++)
      {
         TLorentzVector iP_(iP->px(),iP->py(),iP->pz(),iP->E());
         iP_.Boost(-betaMag*superJetTLV.Px()/superJetTLV.P(),-betaMag*superJetTLV.Py()/superJetTLV.P(),-betaMag*superJetTLV.Pz()/superJetTLV.P());
         boostedSuperJetPart.push_back(fastjet::PseudoJet(iP_.Px(),iP_.Py(),iP_.Pz(),iP_.E()));
      }

      ///reclustering SuperJet that is now boosted into the MPP frame
      double R = 0.6;
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
      fastjet::ClusterSequence cs_jet(boostedSuperJetPart, jet_def); 
      std::vector<fastjet::PseudoJet> jetsFJ_jet = fastjet::sorted_by_E(cs_jet.inclusive_jets(5.0));

      double SJ_50_px = 0, SJ_50_py=0,SJ_50_pz=0,SJ_50_E=0;
      double SJ_70_px = 0, SJ_70_py=0,SJ_70_pz=0,SJ_70_E=0;

      SJ_nAK4_50[nSuperJets] = 0;  SJ_nAK4_70[nSuperJets] =0;
      for (auto iPJ=jetsFJ_jet.begin(); iPJ<jetsFJ_jet.end(); iPJ++)                             
      {
         if(iPJ->E()>50.)
         {
            SJ_AK4_50_mass[tot_nAK4_50] = iPJ->m();
            SJ_50_px+=iPJ->px();SJ_50_py+=iPJ->py();SJ_50_pz+=iPJ->pz();SJ_50_E+=iPJ->E();
            SJ_nAK4_50[nSuperJets]++;
            tot_nAK4_50++;
         }
         if(iPJ->E()> 70.)
         {
            SJ_AK4_70_mass[tot_nAK4_70] = iPJ->m();
            SJ_70_px+=iPJ->px();SJ_70_py+=iPJ->py();SJ_70_pz+=iPJ->pz();SJ_70_E+=iPJ->E();
            SJ_nAK4_70[nSuperJets]++;
            tot_nAK4_70++;
         }
      }
      SJ_mass_50[nSuperJets]= sqrt(pow(SJ_50_E,2)-pow(SJ_50_px,2)-pow(SJ_50_py,2)-pow(SJ_50_pz,2));
      SJ_mass_70[nSuperJets]= sqrt(pow(SJ_70_E,2)-pow(SJ_70_px,2)-pow(SJ_70_py,2)-pow(SJ_70_pz,2));
      boostedSuperJetPart.clear();   //shouldn't be needed, just in case

      nSuperJets++; 
   }
   if((cosAngleChi1*cosAngleChi2 > 0))
   {
      std::cout << "Bad event: SJ1 Mass / SJ2 Mass  - " << superJet_mass[0] << "/" << superJet_mass[1];
      std::cout << "    Objects in SuperJet 1: [";
      if(cosAngleChi1<0)  std::cout << "Chi1 "; 
      if(cosAngleChi2<0)  std::cout << "Chi2 "; 
      if(cosAngleTop<0)   std::cout << "Top "; 
      if(cosAngleW_Suu<0) std::cout << "W_Suu "; 
      if(cosAngleB_Suu<0) std::cout << "b_Suu "; 
      if(cosAngleHiggs<0) std::cout << "Higgs "; 
      std::cout << "].  Objects in SuperJet 2: [";
      if(cosAngleChi1>0)  std::cout << "Chi1 "; 
      if(cosAngleChi2>0)  std::cout << "Chi2 "; 
      if(cosAngleTop>0)   std::cout << "Top "; 
      if(cosAngleW_Suu>0) std::cout << "W_Suu "; 
      if(cosAngleB_Suu>0) std::cout << "b_Suu "; 
      if(cosAngleHiggs>0) std::cout << "Higgs ";
      std::cout << " ]";

      std::cout << "   cos chi1/chi2 w/ TA: " <<  cos(chi1_vec.Angle(thrust_vector)) << "/"<< cos(chi2_vec.Angle(thrust_vector)) << " ";
      std::cout << "   chi deltaR: " << sqrt(pow(chi1.Phi()-chi2.Phi(),2)+pow(chi1.Eta()-chi2.Eta(),2)) << std::endl;
      std::cout << " Masses of jets being considered: ";
      for(auto iJet = fatJets->begin(); iJet != fatJets->end(); iJet++)         ////////Over AK8 Jets
      {
         if((sqrt(pow(iJet->mass(),2)+pow(iJet->pt(),2)) < 100.) || (!(iJet->isPFJet())) || (!isgoodjet(iJet->eta(),iJet->neutralHadronEnergyFraction(), iJet->neutralEmEnergyFraction(),iJet->numberOfDaughters(),iJet->chargedHadronEnergyFraction(),iJet->chargedMultiplicity(),iJet->muonEnergyFraction(),iJet->chargedEmEnergyFraction())) || (iJet->userFloat("ak8PFJetsPuppiSoftDropMass") < 15.) ) continue;
         std::cout << iJet->mass() << " ";
      }
      std::cout << std::endl;
   }
   superJetOne.clear();
   superJetTwo.clear();

/*
   nfatjets = 0;

////////////////Jets//////////////////////////////////////



   for(auto iJet = fatJets->begin(); iJet != fatJets->end(); iJet++)         ////////Over AK8 Jets
   { 
      if(!(iJet->isPFJet())) continue;
      if(!isgoodjet(iJet->pt(),iJet->userFloat("ak8PFJetsPuppiSoftDropMass"), iJet->eta(),iJet->neutralHadronEnergyFraction(), iJet->neutralEmEnergyFraction(),iJet->numberOfDaughters(),iJet->chargedHadronEnergyFraction(),iJet->chargedMultiplicity(),iJet->muonEnergyFraction(),iJet->chargedEmEnergyFraction(),iJet->neutralMultiplicity())) continue;


      jet_pt[nfatjets]         = iJet->pt();
      jet_eta[nfatjets]        = iJet->eta();
      jet_mass[nfatjets]       = sqrt(pow(iJet->energy(),2) - pow(iJet->p(),2));




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
     
      nfatjets++;
   }
   */
   tree->Fill();
}   
DEFINE_FWK_MODULE(clusteringAnalyzer);




