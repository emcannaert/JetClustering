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
#include "Thrust.h"
typedef math::XYZVector Vector;


//////////Written by Ethan Cannaert [2019] Using Stephen Mrenna's thrust algorithm [2007]////////////////////////////////

Thrust::Thrust(std::vector<reco::LeafCandidate> cands, int nconstituents)
{	
	Vector Tnext1, Tnext2,Tnext3,Tnext4;
	Vector Tvec1 = det_init_dir(cands, nconstituents, 1);
	Vector Tvec2 = det_init_dir(cands, nconstituents, 2);
	Vector Tvec3 = det_init_dir(cands, nconstituents, 3);
	Vector Tvec4 = det_init_dir(cands, nconstituents, 4);


	while (!((isConverged1) && (isConverged2) && (isConverged3) && (isConverged4)))
	{
	Tnext1.SetCoordinates(0,0,0);
	Tnext2.SetCoordinates(0,0,0);
	Tnext3.SetCoordinates(0,0,0);
	Tnext4.SetCoordinates(0,0,0);

	Tnext1 = calc_Tnext(cands, nconstituents, Tvec1);			//calc Next thrust axis vector
	Tnext2 = calc_Tnext(cands, nconstituents, Tvec2);
	Tnext3 = calc_Tnext(cands, nconstituents, Tvec3);
	Tnext4 = calc_Tnext(cands, nconstituents, Tvec4);

	isConverged1 = check_convergence(Tvec1, Tnext1);			//check to see if this process has converged
	isConverged2 = check_convergence(Tvec2, Tnext2);
	isConverged3 = check_convergence(Tvec3, Tnext3);
	isConverged4 = check_convergence(Tvec4, Tnext4);

	Tvec1.SetCoordinates(Tnext1.X(),Tnext1.Y(),Tnext1.Z());		//set thrust axis vector to be most recent iteration (Tnext)
	Tvec2.SetCoordinates(Tnext2.X(),Tnext2.Y(),Tnext2.Z());
	Tvec3.SetCoordinates(Tnext3.X(),Tnext3.Y(),Tnext3.Z());
	Tvec4.SetCoordinates(Tnext4.X(),Tnext4.Y(),Tnext4.Z());

	ncyc_++;
	if (ncyc_ > 15) break;
	}

///////////////////////////////////////////////////////////////////////
	double max_thrust    = calc_thrust(cands, nconstituents, Tvec1);
	double second_thrust = calc_thrust(cands, nconstituents, Tvec2);
	double third_thrust  = calc_thrust(cands, nconstituents, Tvec3);
	double fourth_thrust = calc_thrust(cands, nconstituents, Tvec4);

	Vector max_axis;
	double max_thrust_val = 0;

	if ((max_thrust >= second_thrust) && (max_thrust >= third_thrust) && (max_thrust>= fourth_thrust))
	{
		max_thrust_val = max_thrust;
		max_axis.SetCoordinates(Tvec1.X(), Tvec1.Y(), Tvec1.Z());
	}
	else if( ( second_thrust >= third_thrust) && (second_thrust >= fourth_thrust))
	{
		max_thrust_val = second_thrust;
		max_axis.SetCoordinates(Tvec2.X(),Tvec2.Y(),Tvec2.Z());
	}
	else if ((third_thrust >= fourth_thrust))
	{
		max_thrust_val = third_thrust;
		max_axis.SetCoordinates(Tvec3.X(),Tvec3.Y(),Tvec3.Z());
	}
	else 
	{
		max_thrust_val = fourth_thrust;
		max_axis.SetCoordinates(Tvec4.X(), Tvec4.Y(),Tvec4.Z());
	}
////////////////////////////////////////////////////////////////////////

	// calc thrust wiith newly found thrust axis
	axis_.SetCoordinates(max_axis.X(),max_axis.Y(),max_axis.Z());
	thrust_ = max_thrust_val;
}
double Thrust::calc_thrust(std::vector<reco::LeafCandidate> cands, int nconstituents, Vector axis_)
{
   double xcomp = axis_.X();
   double ycomp = axis_.Y();
   double zcomp = axis_.Z();
   double numerator   = 0.;
   double denominator = 0.;
   for (int i =0; i <nconstituents; i++) 
   {  
      denominator += calc_ptot(cands[i].px(), cands[i].py(),cands[i].pz());
      numerator   +=   abs(xcomp*cands[i].px() + ycomp*cands[i].py() + zcomp*cands[i].pz());
   }

   return numerator/denominator;
}


//num is the variable that determines which candidate particle to use as starting point
Vector Thrust::det_init_dir(std::vector<reco::LeafCandidate> cands, int nconstituents, int num)
{

	double max_p    = 0.;
	double second_p = 0.;
	double third_p  = 0.;
	double fourth_p = 0.;

	int max_index     = -999;
	int second_index  = -999;
	int third_index   = -999;
	int fourth_index  = -999;
	//0 particle case
	if (nconstituents <1) return Vector(0,0,0);

	for(int i =0; i<nconstituents; i++)
	{
		if(calc_ptot(cands[i].px(), cands[i].py(), cands[i].pz()) > max_p) 
		{	
			max_index = i;
			max_p = calc_ptot(cands[i].px(), cands[i].py(), cands[i].pz());
		}
		else if ((calc_ptot(cands[i].px(), cands[i].py(), cands[i].pz()) < max_p) && (calc_ptot(cands[i].px(), cands[i].py(), cands[i].pz()) > second_p))
		{
			second_index = i;
			second_p = calc_ptot(cands[i].px(), cands[i].py(), cands[i].pz());
		} 
		else if ((calc_ptot(cands[i].px(), cands[i].py(), cands[i].pz()) < second_p) && (calc_ptot(cands[i].px(), cands[i].py(), cands[i].pz()) > third_p))
		{
			third_index = i;
			third_p = calc_ptot(cands[i].px(), cands[i].py(), cands[i].pz());
		} 
		else if ((calc_ptot(cands[i].px(), cands[i].py(), cands[i].pz()) < third_p) && (calc_ptot(cands[i].px(), cands[i].py(), cands[i].pz()) > fourth_p))
		{
			fourth_index = i;
			fourth_p = calc_ptot(cands[i].px(), cands[i].py(), cands[i].pz());
		}  
	}

	//prevents dividing by by 0 if max_p is 0
	if((max_index == -999) || (second_index == -999) || (third_index = -999) || (fourth_index==-999)) return Vector(0,0,0);
	

	if (num == 1)
	{
	reco::LeafCandidate maxp_cand = cands[max_index];
	return Vector(maxp_cand.px()/max_p, maxp_cand.py()/max_p, maxp_cand.pz()/max_p);
	}
	else if (num == 2)
	{
	reco::LeafCandidate second_cand = cands[second_index];
	return Vector(second_cand.px()/second_p, second_cand.py()/second_p, second_cand.pz()/second_p);
	}
	else if (num == 3)
	{
	reco::LeafCandidate third_cand = cands[third_index];
	return Vector(third_cand.px()/third_p, third_cand.py()/third_p, third_cand.pz()/third_p);
	}
	else if (num == 4)
	{
	reco::LeafCandidate fourth_cand = cands[fourth_index];
	return Vector(fourth_cand.px()/fourth_p, fourth_cand.py()/fourth_p, fourth_cand.pz()/fourth_p);
	}
	else {return Vector(0,0,0);}
}

//determines sign of epsilon value used in thrust axis iterated function, used in calc_Tnext
double Thrust::det_eps(Vector Tvec, double px, double py, double pz)
{

	if ((Tvec.X()*px + Tvec.Y()*py + Tvec.Z()*pz) > 0. ) return 1.;
	else {return -1.;}
}

// calculates total momentum
double Thrust::calc_ptot(double px, double py, double pz)
{
	return sqrt(pow(px,2) + pow(py,2) + pow(pz,2));
}

bool Thrust::check_convergence(Vector Tvec, Vector Tnext)
{

	double conv_val = 0.001;
	Vector Tdiff(Tnext.X() - Tvec.X(), Tnext.Y() - Tvec.Y(), Tnext.Z() - Tvec.Z());

	//Check to see if ratio of lengths of new thrust and previous thrust vectors is less than convergence value (0.01)
	if((calc_ptot(Tdiff.X(),Tdiff.Y(),Tdiff.Z())/calc_ptot(Tvec.X(),Tvec.Y(),Tvec.Z())) < conv_val) return true;
	return false;		

}
Vector Thrust::calc_Tnext(std::vector<reco::LeafCandidate> cands, int nconstituents, Vector Tvec)
{
	double Tnextx = 0.;
	double Tnexty = 0.;
	double Tnextz = 0.;
	for (int i =0; i<nconstituents; i++) 
	{
		double epsilon = det_eps(Tvec, cands[i].px(),cands[i].py(), cands[i].pz());
		Tnextx += epsilon*cands[i].px();
		Tnexty += epsilon*cands[i].py();
		Tnextz += epsilon*cands[i].pz();
	}
	double Tnextl = calc_ptot(Tnextx, Tnexty, Tnextz);
	Vector Tnext(Tnextx/Tnextl, Tnexty/Tnextl, Tnextz/Tnextl);
	return Tnext;
}


double Thrust::thrust()
{
	return thrust_;
}

Vector Thrust::axis()
{
	return axis_;
}

int Thrust::ncyc()
{
	return ncyc_;
}
