//------------------------------------------
// Example code to generate an PS for antinu_mu + proton -> neutron + muon,
// evaluate the LLewellyn-Smith DSigma/DQ2 diff. cross-section in the rest
// frame of the (free) nucleon and integrate the diff. cross section
// ----------------------------------------
//
// To run at CERN lxplus:
// 1. Env. setup: source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_96b x86_64-centos7-gcc62-opt
// 2. Compile and link:
//  g++ `gsl-config --cflags` `root-config --cflags` -c generatePS_AntiNuMu_Proton.cpp
//  g++ `gsl-config --libs` `root-config  --evelibs --glibs` generatePS_AntiNuMu_Proton.o -o generatePS_AntiNuMu_Proton
//
// 3. Run: ./generatePS_AntiNuMu_Proton Enu
// Where Enu is the neutrino energy in GeV
// E.g. ./generatePS_AntiNuMu_Proton 0.6
//-----------------------------------------

#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <vector>

#include <TGenPhaseSpace.h>
#include <TRandom2.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TParticlePDG.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_dilog.h>


#include <iostream>
#include <cmath>


// Constants
double au = 0.931494; // GeV
double mAr_A = 39.948; // Atomic mass number, A
unsigned int mAr_Z = 18; // Atomic number, Z
double mAr = mAr_A*au; // GeV
double mC = 12.0*au;
double mProton = 1.007276466879*au; // GeV
double mNeutron = 1.00866491588*au; // GeV
double mMuon = 0.105; // GeV
double mElectron = 0.511/1000.0; // GeV
double mTau = 1.77686; // GeV
double mPion = 0.1396; // pi+, for pi0 it is 0.135 GeV

// Unit conversion: we calculate in natural units, need to convert to SI ?
double perGeV2_to_b = 0.3894e-03; // 1 GeV^{-2} = 0.3894e-03 barn
double barn_to_cm2 = 1.0e-24; // 1 barn = 1e-24 cm2

// https://github.com/GENIE-MC/Generator/blob/master/config/G18_01b/G18_01b_02_11a/CommonParam.xml
double QEL_Mv = 0.840;
double QEL_Ma = 0.961242;
double QEL_Fa0 = -1.2670;
double CabibboAngle = 0.227780466682;

//https://github.com/GENIE-MC/Generator/blob/master/src/Framework/Conventions/Constants.h
double kGF    = 1.16639E-5;            // Fermi const from b-decay, in GeV^-2
double kGF2   = TMath::Power(kGF,2);
double fCos8c2 = TMath::Power(TMath::Cos(CabibboAngle), 2);
double fXSecScale = 1.0; // for tuning?
double kPi    = TMath::Pi();

double SafetyFactor = 1.6;

// Additional structure holding E0 value to be passed into GSL callbacks.
struct BoundParms {
  double Enu;
  double El;
  double Costheta;
  double Phi;
};



// Nucleon at rest cross section as a function of the
// - neutrino energy
// - fixed outgoing lepton Energy
// - angle between neutrino and outgoing lepton in the Lab frame
// - phi angle of the lepton
// - Q2p is just a pointer to return the calculated Q2
double DSigma_dQ2_XSec(double Q2, double Enu, double El, double Costheta, double Phi){

    // We  assume to numerically evaluate the LL-Smith cross section
  // In the nucleon at rest frame to simplfiy

  TLorentzVector neutrino(0,0, Enu, Enu) ; // GeV, fixed
  TLorentzVector nucleon (0, 0, 0, mProton); // target at rest in the Lab frame

  // we assume El to be fixed
  double lepMass = mMuon;// or mElectron or mMuon or mTau
  double outLeptonEnergy = El; // GeV
  double outMomentum = TMath::Sqrt(outLeptonEnergy*outLeptonEnergy - lepMass*lepMass);
  TVector3 lepton3Mom(0., 0., outMomentum);
  lepton3Mom.SetTheta( TMath::ACos(Costheta) );
  lepton3Mom.SetPhi( Phi );
  TLorentzVector leptonMom(lepton3Mom, outLeptonEnergy);
 
  //-----------------------------------------------------------------------------
  // First we need access to all of the particles in the interaction
  // The particles were stored in the lab frame
  double mNucleon = mProton; // GeV
  double mNi = mNucleon; // Initial struck nucleon mass
  double mNf = mNeutron; // Final nucleon mass
  //  TLorentzVector inNucleonMom = nucleon; // this is the 4mom at the rest frame
  //  TLorentzVector outNucleonMom = outNucleon;
  // TLorentzVector neutrinoMom = neutrino;
  // Ordinary 4-momentum transfer
  //  TLorentzVector qP4 = neutrinoMom - leptonMom;
  

  // Return a differential cross section of zero if we're below threshold (and
  // therefore need to sample a new event)
  // Mandelstam s for the probe/hit nucleon system
  TLorentzVector k4 = neutrino + nucleon;
  double s = k4.Dot(k4);
  if ( std::sqrt(s) < lepMass + mNf ){
    std::cout << "WARNING: Below threshold: s: " << s << " < lepMass + mProton: " << lepMass + mNf << " Returning xsec zero!" << std::endl; 
    return 0.;
  }

  double E  = neutrino.E();
  double E2 = TMath::Power(E,2);
  double ml = lepMass;
  double M  = mNucleon;

  //-----------------------------------------------------------------------------
  // Calculate the QEL form factors
  // ref: https://github.com/GENIE-MC/Generator/blob/master/src/Physics/QuasiElastic/XSection/LwlynSmithFF.cxx
  // F1V ---------------
  double fMv2 = TMath::Power(QEL_Mv,2);
  double q2     = -Q2;
  double Mnucl  = mNucleon;
  double Mnucl2 = TMath::Power(Mnucl, 2);
  double t = q2/(4*Mnucl2);
  // elastic form factors
  double gep = 1. / TMath::Power(1-q2/fMv2, 2);
  double gen = 0;
  double gmp = mProton / TMath::Power(1-q2/fMv2, 2);
  double gmn = mNeutron / TMath::Power(1-q2/fMv2, 2);
  double gve = gep - gen;
  double gvm = gmp - gmn;
  double F1V   = (gve - t*gvm) / (1-t);
  
  // xiF2V -------------
  double xiF2V = (gvm-gve) / (1-t);
  
  // Ax FF ------------
  double FA    = QEL_Fa0 / TMath::Power(1-q2/(QEL_Ma*QEL_Ma), 2);
  
  // P FF -----------
  double Fp    = 2. * mNi*mNi * FA/(mPion*mPion-q2);// using the hit nucleon mass
  
  
  bool is_neutrino = false;
  int sign = (is_neutrino) ? -1 : 1;
  
  // Calculate auxiliary parameters
  double ml2     = TMath::Power(ml,    2);
  double M2      = TMath::Power(M,     2);
  double M4      = TMath::Power(M2,    2);
  double FA2     = TMath::Power(FA,    2);
  double Fp2     = TMath::Power(Fp,    2);
  double F1V2    = TMath::Power(F1V,   2);
  double xiF2V2  = TMath::Power(xiF2V, 2);
  double Gfactor = M2*kGF2*fCos8c2 / (8*kPi*E2);
  double s_u     = 4*E*M + q2 - ml2;
  double q2_M2   = q2/M2;
    

    // Compute free nucleon differential cross section
  double A = (0.25*(ml2-q2)/M2) * (
	      (4-q2_M2)*FA2 - (4+q2_M2)*F1V2 - q2_M2*xiF2V2*(1+0.25*q2_M2)
              -4*q2_M2*F1V*xiF2V - (ml2/M2)*(
               (F1V2+xiF2V2+2*F1V*xiF2V)+(FA2+4*Fp2+4*FA*Fp)+(q2_M2-4)*Fp2));
  double B = -1 * q2_M2 * FA*(F1V+xiF2V);
  double C = 0.25*(FA2 + F1V2 - 0.25*q2_M2*xiF2V2);

  double xsec = Gfactor * (A + sign*B*s_u/M2 + C*s_u*s_u/M4);

  // Check the Q2 limits
  double Q2_min = -1;
  double Q2_max = -1;
  double W = mProton;
  double W2  = TMath::Power(W,  2.);
  double s_   = M2 + 2*M*Enu;

  double auxC = 0.5*(s-M2)/s;
  double aux1 = s_ + ml2 - W2;
  double aux2 = aux1*aux1 - 4*s*ml2;

  Q2_max = -ml2 + auxC * (aux1 + aux2); 
  Q2_min = -ml2 + auxC * (aux1 - aux2); 

  Q2_max = TMath::Max(0., Q2_max);
  Q2_min = TMath::Max(0., Q2_min);

  if (Q2 < Q2_min || Q2 > Q2_max){
    //    std::cout << "Q2 limit issue (Q2min, Q2max, Q2): " << Q2_min << ", " << Q2_max << ", " << Q2 << std::endl;
    return 0.;
  }
  // Apply given scaling factor
  xsec *= fXSecScale;
  
  // Number of scattering centers in the target
  //  xsec *= mAr_Z; // For proton target only 1 
  
  
  //    std::cout << "xsec: " << xsec << " GeV^{-2} at Q2: " << Q2 << " GeV^2" << std::endl;
  // std::cout << "xsec: " << xsec * perGeV2_to_b * barn_to_cm2 << " cm^2/GeV^2 at Q2: " << Q2 << " GeV^2" << std::endl;
 

  return xsec;

}

static double _DSigma_dQ2_XSec(double Q2, void * parms_) {
    //BoundParms * parms = (BoundParms*) parms_;  // or, equivalently, in C++ style
    BoundParms * parms = reinterpret_cast<BoundParms*>(parms_);
    // Forward invocation to target method
    return DSigma_dQ2_XSec( Q2, parms->Enu, parms->El, parms->Costheta, parms->Phi );
}

double TotalCrossSection(double Enu, double El, double Costheta, double Phi){
  double sigmaTot;

  gsl_error_handler_t * old_handler=gsl_set_error_handler_off();

  gsl_integration_workspace* w1 = gsl_integration_workspace_alloc (1000);
  double result1, error1;

  //Quick fix
  TLorentzVector neutrino(0,0,Enu, Enu);
  TLorentzVector nucleon(0, 0, 0, mProton);
  TLorentzVector k4 = neutrino + nucleon;
  double s = k4.Dot(k4);
  
  // Check the Q2 limits
  double Q2_min = -1;
  double Q2_max = -1;
  double W = mNeutron;
  double W2  = TMath::Power(W,  2.);

  double ml = mMuon;
  double M  = mProton;
  double M2      = TMath::Power(M,     2);
  double ml2     = TMath::Power(ml,    2);
  
  double s_   = M2 + 2*M*Enu;
  
  double auxC = 0.5*(s-M2)/s;
  double aux1 = s_ + ml2 - W2;
  double aux2 = aux1*aux1 - 4*s*ml2;

  Q2_max = -ml2 + auxC * (aux1 + aux2); 
  Q2_min = -ml2 + auxC * (aux1 - aux2); 

  Q2_max = TMath::Max(0., Q2_max);
  Q2_min = TMath::Max(0., Q2_min);

  std::cout << "Integrating DSigma_DQ2 from Q2min: " << Q2_min << " GeV2 to " << Q2_max << " GeV2" << std::endl;
  
  gsl_function F1;
  BoundParms parms = {Enu, El, Costheta, Phi};
  F1.function = _DSigma_dQ2_XSec;
  F1.params = &parms;

  //gsl_integration_qags (&F1, Xmin1, Xmax1, 0, 1e-7, 1000, w1, &result1, &error1);
  double relerr=1.0e-7;   //initial error tolerance (relative error)
  int status=1;
  while(status) {
    status=gsl_integration_qags (&F1, Q2_min, Q2_max, 0, relerr, 1000, w1, &result1, &error1);
    relerr *= 1.2;
    if(status) std::cout << "Increased tolerance=" << relerr << std::endl;
  }
  //if integration routine returns error code, integration is repeated
  //using increased error tolerance, message is printed out
  gsl_set_error_handler(old_handler); //reset error handler (might be unneccessary.)

  double IntDsDx = result1;
  gsl_integration_workspace_free (w1);

  sigmaTot= IntDsDx;
  return sigmaTot;
  
}



int  main(int argc, char * argv[]){

  TDatabasePDG * fPDGdb = new TDatabasePDG();
  TRandom2 * ran = new TRandom2();
  ran->SetSeed(0);

  double Eneutrino = atof(argv[1]);

  
  // The LL-Smith cross section (with nucleon at Rest) needs some variables "fixed"
  // - neutrino energy (the beam)
  // - outgoing lepton energy
  // - angle between them (cos theta and phi)
  //
  // These are randomly generated by the TGenPhaseSpace according to its rules.
  // The initial state 4-mom is generated in the Lab frame (beam and nucleon at rest)
  // And itt gives back the final state 4-momenta in the Lab frame.
  //
  // Then we will use the Neumann acceptance-rejection method to decide which
  // throws of of 4-momentum in the final state from TGenPhaseSpace to accept.

  
  double mBeam = Eneutrino; // GeV neutrino


  TH1F *hPt_Muon = new TH1F("hPt_muon","MuonPt", 200,0.0,2.0);
  TH1F *hPz_Muon = new TH1F("hPz_muon","MuonPz", 200,-2.0,2.1);
  TH1F *hE_Muon = new TH1F("hE_muon","MuonE", 200,0,2.1);
  TH1F *hPhi_Muon = new TH1F("hPhi_muon","MuonPhi", 200,0,kPi);
  TH1F *hQ2 = new TH1F("hQ2","Q2", 200,0,2);
  TH1F *hEnuReco = new TH1F("hEnuReco","Enu Reco", 200,0,4);

  Int_t Nevents = 100;
  Int_t Nacc = 0.0;
  for (Int_t nEv=0;nEv<Nevents;nEv++) {
  

    double mNi = mProton;
    double mNf = mNeutron;

    // Initial nucleon 4-momentum (nucleus rest, Lab frame) and removal energy from Bodek Ritchie
    TLorentzVector pNi(0.0, 0.0, 0.0, mNi);

    //--------------------------------------------------
    // NOTE: We have various reference frames to consider:
    // 1) Lab frame:
    // - in which the neutrino 4-mom is p(0, 0, Enu, Enu)
    // - and the nucleus is at rest with a 4 mom pN(0, 0, 0, MN)
    // - but in which the nucleons are in Fermi motion
    //
    // 2) Nucleon at rest frame: (in which the Llewellyn-Smith
    // CC QEL diff. cross-section is calculated, i.e. where the nucleon is 'free')
    // - in this frame the nucleon 4-mom is p(0, 0, 0, Mn) as it is at rest
    // - but the neutrino will have an energy distribution
    //
    // 3) The beam + nucleon COM frame:
    // - In this frame the outgoing nucleon and outgoing lepton
    // are back to back
    // - So, if we fix the final state lepton kinematics, then
    // this constrains the outgoing nucleon kinematics in this frame.

    
    // Neutrino beam in the Lab frame
    TLorentzVector beam(0.0,0.0, mBeam, mBeam);//GeV (beam "particle")



    // Get beta for a Lorentz boost from the Lab frame to the Nucleon at rest frame
    // where the target nucleon is at rest, but the beam energy has some distribution
    TVector3 beta_lab_to_nRest = (pNi).BoostVector();// pNi is the nucleon momentum in the Lab frame
    TVector3 beta_nRest_to_lab = -beta_lab_to_nRest;

    // Also get the Lorentz boost beta for the beam+target COM frame
    // in which the outgoing lepton and outgoing nucleon
    // are back-to-back (to be used later)
    TVector3 beta_lab_to_beamtargetCOM = (pNi+beam).BoostVector();
    TVector3 beta_beamtargetCOM_to_lab = -beta_lab_to_beamtargetCOM;

    // Let's boost from the Lab frame to the Nucleon at rest frame
    // where we can evaluate the LL-Smith cross section
    TLorentzVector beam_nRest = TLorentzVector( beam );
    beam_nRest.Boost( beta_lab_to_nRest );
    TLorentzVector nucleon_nRest = TLorentzVector( pNi );
    nucleon_nRest.Boost( -beta_lab_to_nRest );
    
    double Enu_Lab = beam.E();
    double Enu = beam_nRest.E();

    /* // Check */
    /* std::cout << "nu Lab: " << std::endl;  */
    /* beam.Print();  */
    /* std::cout << ", nu nRest: " << std::endl;  */
    /* beam_nRest.Print();  */
    /* std::cout << "Nucleon in Lab: "  << std::endl; */
    /* pNi.Print();  */
    /* std::cout << "Nucleon in nRest: " << std::endl;  */
    /* nucleon_nRest.Print();  */
    /* std::cout << "-----------------------------------------" << std::endl; */


    
    // Generate the random PS in frame in which the nucleon is at rest
    TLorentzVector target = TLorentzVector(nucleon_nRest);

    // The beam neutrino is already boosted into this frame
    TLorentzVector W = beam_nRest + target;

    double lepMass = mMuon;// mMuon or mElectron or mTau
    Double_t masses[2] = { mNeutron , lepMass};

    // Generate the Decay in then nucleon at rest frame
    TGenPhaseSpace event;
    bool allowed = event.SetDecay(W, 2, masses);

    Double_t weight = event.Generate();
    //    std::cout << "weight: " << weight << ", allowed: " << allowed << std::endl;
    TLorentzVector *pProton = event.GetDecay(0);
    TLorentzVector *pMuon  = event.GetDecay(1);
    double El = pMuon->E(); // lepton energy at the nucleon at rest frame 

    //------------------------------ 
    // Before moving forward, let's check the outgoing nucleon momentum
    // and apply Pauli blocking. If outgoing nucleon is below Fermi momentum
    // we skip this event generation

    // The Fermi momentum in the Lab frame
    //    double kF = FGMFermiMomentum(mAr_A, mAr_Z);

    // We need to evaluate for this the outgoing nucleon momentum in the Lab frame
    // We have the outgoing lepton in the nucleon-at-rest frame:
    // We know that the final-state nucleon will have an equal and opposite 3-momentum
    // in the beam+target COM frame and will be on the mass shell.
    // But for this we have to boost the system from the nucleon-at-rest frame
    // to the beam+target COM frame in 2 steps (nucleon at rest -> Lab -> beam+targetCOM)
    TLorentzVector outlepton_COM = TLorentzVector( *pMuon ); // pMuon is at the nucleon at rest frame
    // first boost it to the nucleus at rest( Lab) frame
    outlepton_COM.Boost( beta_nRest_to_lab );
    // then boost it to the beam+target COM frame
    outlepton_COM.Boost( beta_lab_to_beamtargetCOM );

    
    // In this beam+target COM frame the outgoing nucleon and lepton are back-to-back
    TLorentzVector outNucleon_COM(-1*outlepton_COM.Px(),-1*outlepton_COM.Py(),-1*outlepton_COM.Pz(),
			      TMath::Sqrt(outlepton_COM.P()*outlepton_COM.P() + mNf*mNf));

    // Now boost the outgoing nucleon into the Lab frame
    TLorentzVector outNucleon_Lab = TLorentzVector(outNucleon_COM);
    outNucleon_Lab.Boost(beta_beamtargetCOM_to_lab);


    // Apply Pauli blocking
    /* if(outNucleon_Lab.P() < kF){ */
    /*   std::cout << "Fermi Momentum: " << kF << ", Nucleon momentum: " << outNucleon_Lab.P() << std::endl;     */
    /*   continue; */
    /* } */

    //------------------------------ 
    // Evaluate the diff. xsec dsigma/dQ2 at the random throw
    // in the target nucleon at rest frame
    TVector3 beam3V = beam_nRest.Vect();
    TVector3 lepton3V = pMuon->Vect();
    double Costheta = TMath::Cos(lepton3V.Angle(beam3V));

    // To get the phi vector around the random neutrino (beam, whatever) direction,
    // we rotate the lepton 3-mom so that the old z-axis vector to point along the neutrino direction
    TVector3 zvec(0., 0., 1.);
    TVector3 rot_axis = ( zvec.Cross(beam3V) ).Unit(); // the rotation axis
    double rot_angle = beam3V.Angle(zvec); // amount of rotation
    // Rotate the lepton3V
    lepton3V.Rotate(rot_angle, rot_axis);
    // measure the phi now for the rotated lepton momentum vector
    double phi = lepton3V.Phi();

    // finally evaluate the Xsec in the nucleon at rest frame
    
    TLorentzVector neutrino(0,0, Enu, Enu) ; // GeV, fixed
    TLorentzVector qP4 = neutrino - *pMuon;
    double Q2 = -1. * qP4.Mag2();

    double xsec = DSigma_dQ2_XSec(Q2, Enu, El, Costheta, phi);
    //    double R = 1.0;
    //    R = NuclQELXSecSuppression(mAr_A, -Q2, mNi,  0.5);
    //  std::cout << "Q2: " <<  Q2 << ", Supression: " << R << std::endl;
    //    xsec = R*xsec;
    
    //--------------------------------------------------
    // Now, we we apply the Rejection Sampling algorithm

    // Now we need to find the Maximum diff. cross section that covers this 
    // function
    const double acceptable_fraction_of_safety_factor = 0.2;
    const int max_n_layers = 100;
    const int N_theta = 10;
    const int N_phi = 10;
    double costh_range_min = -1.; 
    double costh_range_max = 1.;
    double phi_range_min = 0.;
    double phi_range_max = 2.*TMath::Pi();
    double phi_at_xsec_max = -1.;
    double costh_at_xsec_max = 0.;
    double this_nuc_xsec_max = -1;
    
    double xsec_max = 0;

    for (int ilayer = 0 ; ilayer < max_n_layers ; ilayer++) {
      double last_layer_xsec_max = this_nuc_xsec_max;
      double costh_increment = (costh_range_max-costh_range_min) / N_theta;
      double phi_increment   = (phi_range_max-phi_range_min) / N_phi;
    
      // Now scan through angles coarsely with Enu and El fixed
      for (int itheta = 0; itheta < N_theta; itheta++){
	double costh_ = costh_range_min + itheta * costh_increment;
	for (int iphi = 0; iphi < N_phi; iphi++) { // Scan around phi
	  double phi_ = phi_range_min + iphi * phi_increment;

	  TLorentzVector neutrino_(0,0, Enu, Enu) ; // GeV, fixed
	  double outMomentum_ = TMath::Sqrt(El*El - mMuon*mMuon);
	  TVector3 lepton3Mom_(0., 0., outMomentum_);
	  lepton3Mom_.SetTheta( TMath::ACos(costh_) );
	  lepton3Mom_.SetPhi( phi_ );
	  TLorentzVector lepton_(lepton3Mom_, El);
	  TLorentzVector qP4_ = neutrino_ - lepton_;
	  double Q2_ = -1. * qP4_.Mag2();	  
	  double xs = DSigma_dQ2_XSec(Q2_, Enu, El, costh_, phi_); 
	  //	R= 1;
	  //	R = NuclQELXSecSuppression(mAr_A, -Q2_, mNi,  0.5);
	  //	xs = R*xs;
	  
	  if (xs > this_nuc_xsec_max){
	    phi_at_xsec_max = phi_;
	    costh_at_xsec_max = costh_;
	    this_nuc_xsec_max = xs;
	  }
	  
	} // Done with phi scan
      }// Done with angles coarsely
      
      // Calculate the range for the next layer
      costh_range_min = costh_at_xsec_max - costh_increment;
      costh_range_max = costh_at_xsec_max + costh_increment;
      phi_range_min = phi_at_xsec_max - phi_increment;
      phi_range_max = phi_at_xsec_max + phi_increment;

      double improvement_factor = this_nuc_xsec_max/last_layer_xsec_max;
      if (ilayer && (improvement_factor-1) < acceptable_fraction_of_safety_factor * (SafetyFactor-1)) {
	break;
      }
      
    }
    if (this_nuc_xsec_max > xsec_max) {
      xsec_max = this_nuc_xsec_max;
    }

    // std::cout << " -------------------- best estimate for xsec_max = " << xsec_max*perGeV2_to_b * barn_to_cm2 << std::endl;

    
    // apply safety factor
    xsec_max *= SafetyFactor;

    // Now we apply the acceptance-reject method criteria
    double t = xsec_max * ran->Rndm();

    bool  accept = (t < xsec);
    

    if(accept){
      Nacc += 1;
      // std::cout << "\r " << nEv <<   std::flush;
      std::cout << "Diff. xsec: " << xsec << " GeV^{-2} " << ", " << xsec*perGeV2_to_b * barn_to_cm2 << " cm2 / GeV2" << std::endl;
      double total_xs = TotalCrossSection(Enu, El,  Costheta, phi);
      std::cout << "Total xsec: " << total_xs * perGeV2_to_b * barn_to_cm2 << " cm^2" << std::endl;
      break;

      // OK, Now we need to boost back the particle(s) into the Lab frame,
      // where the target nucleon is not at rest, but the beam and Nucleus is
      TLorentzVector pMuon_Lab = TLorentzVector(*pMuon);
      pMuon_Lab.Boost(beta_nRest_to_lab);

      // Need to recalculate Q2 in the Lab frame
      TLorentzVector qP4 = beam - (pMuon_Lab);
      double Q2Lab = -1. * qP4.Mag2();      
      
      double pt_muon = TMath::Sqrt(pMuon_Lab.Px()*pMuon_Lab.Px() + pMuon_Lab.Py()*pMuon_Lab.Py());
      double pz_muon = pMuon_Lab.Pz();
      double p_muon = TMath::Sqrt(pMuon_Lab.Px()*pMuon_Lab.Px() + pMuon_Lab.Py()*pMuon_Lab.Py() + pz_muon*pz_muon);
      double Theta_Muon_Nu = (pMuon_Lab.Vect()).Angle(beam.Vect());
      double phi_muon = pMuon_Lab.Phi();
      // if we want weighted events we need to multiply with the TGenPS weight
      double w_ = 1.0;// weight;
      hPt_Muon->Fill(pt_muon, w_); 
      hPz_Muon->Fill(pz_muon, w_);
      hPhi_Muon->Fill(phi_muon, w_);
      hE_Muon->Fill(pMuon_Lab.E(), w_);
      hQ2->Fill(Q2Lab, w_);


      // Try to reconstruct the neutrino energy
      double Eb = 0.0;
      double EnuReco =0.5*( mProton*mProton - TMath::Power((mNeutron - Eb), 2) + lepMass*lepMass + 2*(mNeutron - Eb)*pMuon_Lab.E())/ (mNeutron - Eb - pMuon_Lab.E()+ p_muon*TMath::Cos(Theta_Muon_Nu));
      hEnuReco->Fill(EnuReco, w_);
    }

  }
  std::cout << std::endl;

  // TCanvas* c = new TCanvas("c", "c", 1);
  // c->SetTitle("TGenPhaseSpace QEL CC nu on Argon");
  // c->Divide(2,2);
  // c->cd(1);
  // hPt_Muon->Draw("hist");
  // c->cd(2);
  // hPz_Muon->Draw("hist");
  // c->cd(3);
  // hQ2->Draw("hist");
  // c->cd(4);
  // hE_Muon->Draw("hist");

  // TCanvas* c1 = new TCanvas("c1", "c1", 1);
  // hEnuReco->Draw("hist");


  
  return 0;
}

