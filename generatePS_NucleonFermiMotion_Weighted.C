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

#include <iostream>

TDatabasePDG * fPDGdb = new TDatabasePDG();
TRandom2 * ran = new TRandom2();


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
double QEL_Ma = 0.990;
double QEL_Fa0 = -1.2670;
double CabibboAngle = 0.227780466682;

//https://github.com/GENIE-MC/Generator/blob/master/src/Framework/Conventions/Constants.h
double kGF    = 1.16639E-5;            // Fermi const from b-decay, in GeV^-2
double kGF2   = TMath::Power(kGF,2);
double fCos8c2 = TMath::Power(TMath::Cos(CabibboAngle), 2);
double fXSecScale = 1.0; // for tuning?
double kPi    = TMath::Pi();

double SafetyFactor = 1.6;



double BindEnergyPerNucleonParametrization(double A,unsigned int Z){
  
  // Compute the average binding energy per nucleon (in GeV)
  double x = TMath::Power(A,1/3.0) / Z;
  
  return (0.05042772591+x*(-0.11377355795+x*(0.15159890400-0.08825307197*x)));

}


double FermiMomentumForIsoscalarNucleonParametrization(double A){
  // Compute Fermi momentum for isoscalar nucleon (in GeV)
  double x = 1.0 / A;
  return (0.27+x*(-1.12887857491+x*(9.72670908033-39.53095724456*x)));
}

double FmArea(double alpha, double beta, double kf, double pmax)
{
// Adapted from NeuGEN's fm_area() used in r_factor()

  double kf3 = TMath::Power(kf,3.);
  double sum = 4.*kPi* (alpha * kf3/3. + beta*(1./kf - 1./pmax));
  return sum;
}

double FmI2(double alpha, double beta,
               double a, double b,  double kFi, double /*kFf*/, double /*q*/)
{
// Adapted from NeuGEN's fm_integral2() used in r_factor()

  double integral2 = 0;

  if(kFi < a) {
     integral2 = beta * (1./a - 1./b);

  } else if(kFi > b) {
     double a3 = TMath::Power(a,3);
     double b3 = TMath::Power(b,3);
     integral2 = alpha/3. * (b3 - a3);

  } else {
     double a3   = TMath::Power(a,3);
     double kFi3 = TMath::Power(kFi,3);

     integral2 = alpha/3. * (kFi3 - a3) + beta * (1./kFi - 1./b);
  }
  return integral2;
}


double FmI1(double alpha, double beta, 
                       double a, double b,  double kFi, double kFf, double q)
{
// Adapted from NeuGEN's fm_integral1() used in r_factor()

  double f=0;

  double q2   = TMath::Power(q,  2);
  double a2   = TMath::Power(a,  2);
  double b2   = TMath::Power(b,  2);
  double kFi2 = TMath::Power(kFi,2);
  double kFf2 = TMath::Power(kFf,2);

  if(kFi < a) {
     double lg = TMath::Log(b/a);

     f = -beta * (kFf2-q2)/(4.*q) * (1./a2 - 1./b2) + beta/(2.*q) * lg;

  } else if (kFi > b) {
     double a4 = TMath::Power(a2,2);
     double b4 = TMath::Power(b2,2);

     f = - (kFf2-q2) * alpha/(4.*q) * (b2-a2) + alpha/(8.*q) * (b4-a4);

  } else {
     double a4   = TMath::Power(a2,2);
     double kFi4 = TMath::Power(kFi2,2);
     double lg   = TMath::Log(b/kFi);

     f = - alpha*(kFf2-q2)/(4.*q)*(kFi2-a2) + alpha/(8.*q)*(kFi4 - a4)
         - beta*(kFf2-q2)/(4.*q)*(1./kFi2 - 1./b2) + beta/(2.*q)*lg;
  }

  double integral2 = FmI2(alpha,beta,a,b,kFi,kFf,q);
  double integral1 = integral2 + f;
  
  return integral1;
}


// Compute nuclear suppression factor - from GENIE
// q2: momentum transfer
// Mn: hit nucleon mass
// pmax: Fermi momentum cutoff
double NuclQELXSecSuppression(int A, int Z, double q2, double Mn, double pmax = 0.5){
  double R = 1;

  double kFi, kFf;
  https://github.com/GENIE-MC/Generator/blob/master/config/FermiMomentumTables.xml
  //      <!-- Ar40  --> <kf nucleus_pdgc="1000180400"> <p> 0.242 </p> <n> 0.259 </n> </kf> 
  //  kFi = FermiMomentumForIsoscalarNucleonParametrization(mAr_A);
  //  kFf = FermiMomentumForIsoscalarNucleonParametrization(mAr_A);
  kFi = 0.259;// neutron
  kFf = 0.242;// proton
  
  // LFG
  /* double hbarc = kLightSpeed*kPlankConstant/genie::units::fermi; */
  
  /* // kFi */
  /* double numNuci = A-Z; // N */
  /* kFi = TMath::Power(3*kPi2*numNuci* */
  /* 		     genie::utils::nuclear::Density(radius,A),1.0/3.0) *hbarc; */
  /* // kFi */
  /* bool is_p_f = pdg::IsProton(final_nucleon_pdgc); */
  /* double numNucf = (is_p_f) ? (double)tgt->Z():(double)tgt->N(); */
  /* kFf = TMath::Power(3*kPi2*numNucf* */
  /* 		     genie::utils::nuclear::Density(radius,A),1.0/3.0) *hbarc; */
  

  // Compute magnitude of the 3-momentum transfer to the nucleon
  double Mn2    = Mn*Mn;
  double magq2  = q2 * (0.25*q2/Mn2 - 1.);
  double q      = TMath::Sqrt(TMath::Max(0.,magq2));

  double kfa   = kFi * 2./kPi;
  double kfa2  = TMath::Power(kfa,2);
  double kFi4  = TMath::Power(kFi,4);
  double rkf   = 1./(1. - kFi/4.);
  double alpha = 1. - 6.*kfa2;
  double beta  = 2. * rkf * kfa2 * kFi4;

  double fm_area = FmArea(alpha,beta,kFi,pmax);

  if (q <= kFf) {

     double p1   = kFf - q;
     double p2   = kFf + q;
     double fmi1 = FmI1(alpha,beta,p1,p2,  kFi,kFf,q);
     double fmi2 = FmI2(alpha,beta,p2,pmax,kFi,kFf,q);

     R = 2*kPi * (fmi1 + 2*fmi2) / fm_area;

  } else if (q > kFf && q <= (pmax-kFf)) {

     double p1    = q - kFf;
     double p2    = q + kFf;
     double fmi1  = FmI1(alpha,beta, p1, p2,   kFi, kFf,q);
     double fmi2a = FmI2(alpha,beta, 0., p1,   kFi, kFf,q);
     double fmi2b = FmI2(alpha,beta, p2, pmax, kFi, kFf,q);

     R = 2*kPi * (fmi1 + 2*(fmi2a+fmi2b)) / fm_area;

  } else if (q > (pmax-kFf) && q <= (pmax+kFf)) {

     double p1   = q - kFf;
     double fmi2 = FmI2(alpha,beta, 0.,p1,  kFi,kFf,q);
     double fmi1 = FmI1(alpha,beta, p1,pmax,kFi,kFf,q);

     R = 2*kPi * (2.*fmi2 + fmi1) / fm_area;

  } else if (q > (pmax+kFf)) {
     R = 1.; 
  } else {
    std::cerr << "Illegal input q = " << q << std::endl;
     exit(1);
  }
  
  return R;
  

}



// Nucleon at rest cross section as a function of the
// - neutrino energy
// - fixed outgoing lepton Energy
// - angle between neutrino and outgoing lepton in the Lab frame
// - phi angle of the lepton
// - Q2p is just a pointer to return the calculated Q2
double DSigma_dQ2_XSec(double Enu, double El, double Costheta, double Phi, double * Q2p){

  // We  assume to numerically evaluate the LL-Smith cross section
  // In the nucleon at rest frame to simplfiy

  TLorentzVector neutrino(0,0, Enu, Enu) ; // GeV, fixed
  TLorentzVector nucleon (0, 0, 0, mNeutron); // target at rest in the Lab frame
  // we assume El to be fixed
  double lepMass = mMuon;// or mElectron or mMuon or mTau
  double outLeptonEnergy = El; // GeV
  double outMomentum = TMath::Sqrt(outLeptonEnergy*outLeptonEnergy - lepMass*lepMass);
  TVector3 lepton3Mom(0., 0., outMomentum);
  lepton3Mom.SetTheta( TMath::ACos(Costheta) );
  lepton3Mom.SetPhi( Phi );
  TLorentzVector leptonMom(lepton3Mom, outLeptonEnergy);
  //  TLorentzVector outNucleon;

 
  //-----------------------------------------------------------------------------
  // First we need access to all of the particles in the interaction
  // The particles were stored in the lab frame
  double mNucleon = mNeutron; // GeV
  double mNi = mNeutron; // Initial struck nucleon mass
  TLorentzVector inNucleonMom = nucleon; // this is the 4mom at the rest frame
  //  TLorentzVector outNucleonMom = outNucleon;
  TLorentzVector neutrinoMom = neutrino;
  // Ordinary 4-momentum transfer
  TLorentzVector qP4 = neutrinoMom - leptonMom;
  

  // Return a differential cross section of zero if we're below threshold (and
  // therefore need to sample a new event)
  // Mandelstam s for the probe/hit nucleon system
  TLorentzVector k4 = neutrino + nucleon;
  double s = k4.Dot(k4);
  if ( std::sqrt(s) < lepMass + mProton ){
    std::cout << "WARNING: Below threshold: s: " << s << " < lepMass + mProton: " << lepMass + mProton << " Returning xsec zero!" << std::endl; 
    return 0.;
  }
  double Q2 = -1. * qP4.Mag2();
  
  double E  = nucleon.E();
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
  double Fp    = 2. * mProton*mProton * FA/(mPion*mPion-q2);// using the hit nucleon mass
  
  
  bool is_neutrino = true;
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
  double aux1 = s + ml2 - W2;
  double aux2 = aux1*aux1 - 4*s*ml2;

  Q2_max = -ml2 + auxC * (aux1 + aux2); 
  Q2_min = -ml2 + auxC * (aux1 - aux2); 

  Q2_max = TMath::Max(0., Q2_max);
  Q2_min = TMath::Max(0., Q2_min);

  if (Q2 < Q2_min || Q2 > Q2_max) return 0.;
  
  // Apply given scaling factor
  xsec *= fXSecScale;
  
  // Number of scattering centers in the target
  xsec *= mAr_Z; 
  
  
  //    std::cout << "xsec: " << xsec << " GeV^{-2} at Q2: " << Q2 << " GeV^2" << std::endl;
  //  std::cout << "xsec: " << xsec * perGeV2_to_b * barn_to_cm2 << " cm^2/GeV^2 at Q2: " << Q2 << " GeV^2" << std::endl;
  
  // Set Q2
  *Q2p = Q2;

  return xsec;

}
/*********************************
c
c *** E. Iacopini 13/Febbraio/1997
c
c --- Determina la distribuzione in impulso dei nucleoni a causa
    c --- del moto di Fermi, secondo la parametrizzazione di
c --- S. Fantoni e O. Benhar.
c --- L'impulso x deve essere dato in GeV/c e la distribuzione e'
c --- normalizzata ad area unitaria.
c
c     y   =  impulso in fm-1
c     as  =  densita' di probabilita' da onda S
c     ap  =  densita' di probabilita' da onda P
c     at  =  coda della densita' di probabilita'
c     
************************************/

double FB(double x){
  double zs=0.800;
  double zp=0.907;
  double bs=1.700;
  double bp=1.770;
  double alpha=1.5;
  double beta=0.8;
  double c1=2.823397;
  double c2=7.225905;
  double c3=0.00861524;

  double y=x/0.197;

  double as=c1*TMath::Exp(-(bs*y)*(bs*y));
  double ap=c2*((bp*y)*(bp*y))*TMath::Exp(-(bp*y)*(bp*y));
  double at=c3*TMath::Power(y,beta)*TMath::Exp(-alpha*(y-2.));

  double rr=(3.14159265/4.)*((as+ap+at)*(y*y))/0.197;

  // --- rinormalizziamo in modo che la distribuzione in dx (GeV/c)
  // --- abbia area unitaria
  
  double fb=rr/1.01691371;

  return fb;
}

/***********************************************
*CMZ :  3.00/00 17/03/92  11.22.25  by  A. Rubbia
*-- Author :
c --- It evaluates the Fermi momentum of a nucleon,
c --- according to the parametrization given by Fantoni/Benhar
c
****** E. Iacopini     13/feb/1997
c
c --- We extract according to the fact that the distribution
c --- between 0 and 0.4 Gev/c has a maximum value not exceeding 
c --- am1=6.4, and its integral (between 0 and 0.4) is  0.91993,
c --- whereas between 0.4 and 2.5 Gev/c the maximum of the
c --- distribution does not exceed am2=0.25 ...
	    c --- (above 2.5 Gev/c, we can assume that the distribution 
		   c --- value is zero !). 
***********************************************/
TLorentzVector PFermiMotion(){
  TRandom2 ran(0);

  // units [GeV/c]
  double r=0.;
  double cint=0.91993;
  double pp1=0.4;
  double am1=6.4;
  double pp2=2.1;
  double am2=0.25;
  double cutoff = 0.251;

 label100:
  
  // generates random number in the range ]0,1[
  double a=ran.Rndm();
  double x = 0.0;
  if(a < cint){
  label10:
    // first zone
    x=ran.Rndm()*pp1;
    double c=FB(x);
    double y=am1*ran.Rndm();
    if(y > c)
      goto label10;
  }else{
  label20:
    //   second zone
    x=pp1+ran.Rndm()*pp2;
    double c=FB(x);
    double y=am2*ran.Rndm();
    if(y > c)
      goto label20;
      
  }
  double pf=x;
  double COSTHETA = 1.0-2.0*ran.Rndm();
  double SINTHETA = TMath::Sqrt(1.0-COSTHETA*COSTHETA);
  double PHI = 2.*3.14159265*ran.Rndm();
  // Now find the components of the momentum...
  double px = pf*SINTHETA*TMath::Cos(PHI);
  double py = pf*SINTHETA*TMath::Sin(PHI);
  double pz = pf*COSTHETA;
  
  // apply cut
  if(TMath::Sqrt(px*px+py*py+pz*pz) > cutoff){
    //    std::cout << "FERMI: parametrization returned too high value, ignoring..." << " px,py,pz : " << px << ", " << py << ", " << pz << std::endl;
    goto label100;
  }
  
  TLorentzVector nucl;
  TParticlePDG * prot = fPDGdb->GetParticle(2212);
  TParticlePDG * neut = fPDGdb->GetParticle(2112);  
  TParticlePDG * nucleon = neut;
  nucl.SetPx(px);
  nucl.SetPy(py);
  nucl.SetPz(pz);
  nucl.SetE(TMath::Sqrt(nucleon->Mass()*nucleon->Mass() + nucl.Px()*nucl.Px() + nucl.Py()*nucl.Py() + nucl.Pz()*nucl.Pz()));

  return nucl;
  
}


// Setup BR Prob distro
TH1D * setupBodekRitchieProbDistro(){

  // Default value 4.0 from original paper by A. Bodek and J. L. Ritchie. Phys. Rev. D 23, 1070
  double fPMax = 4.0;
  double fPCutOff = 0.5;
  double KF = FermiMomentumForIsoscalarNucleonParametrization(mAr_A);

  double a  = 2.0;
  double C  = 4. * kPi * TMath::Power(KF,3) / 3.;
  double R  = 1. / (1.- KF/fPMax);
   //-- create the probability distribution
  int npbins = (int) (1000*fPMax);
  TH1D * prob = new TH1D("", "", npbins, 0, fPMax);

  double dp = fPMax / (npbins-1);
  double iC = (C>0) ? 1./C : 0.;
  double kfa_pi_2 = TMath::Power(KF*a/kPi,2);
  
  for(int i = 0; i < npbins; i++) {
    double p  = i * dp;
    double p2 = TMath::Power(p,2);
    
    // calculate |phi(p)|^2
    double phi2 = 0;
    if (p <= KF)
      phi2 = iC * (1. - 6.*kfa_pi_2);
    else if ( p > KF && p < fPCutOff)
      phi2 = iC * (2*R*kfa_pi_2*TMath::Power(KF/p,4.));
    
    // calculate probability density : dProbability/dp
    double dP_dp = 4*kPi * p2 * phi2;
    prob->Fill(p, dP_dp);
  }

  //-- normalize the probability distribution
  prob->Scale( 1.0 / prob->Integral("width") );

  return prob;

}
  

// From Genie Physics/NuclearState/FGMBodekRitchie class
TLorentzVector generateNucleonBodekRitchie(TH1D * prob, double * fRemovalE){


  double fCurrRemovalEnergy = 0;
  TVector3 fCurrMomentum(0,0,0);

  //-- set fermi momentum vector
  //
  if ( ! prob ) {
    std::cerr << "Null nucleon momentum probability distribution" << std::endl;
    exit(1);
  }
  double p = prob->GetRandom();
  //  LOG("BodekRitchie", pINFO) << "|p,nucleon| = " << p;



  double costheta = -1. + 2. * ran->Rndm();
  double sintheta = TMath::Sqrt(1.-costheta*costheta);
  double fi       = 2 * kPi * ran->Rndm();
  double cosfi    = TMath::Cos(fi);
  double sinfi    = TMath::Sin(fi);

  double px = p*sintheta*cosfi;
  double py = p*sintheta*sinfi;
  double pz = p*costheta;

  fCurrMomentum.SetXYZ(px,py,pz);

  //-- set removal energy
  //
  fCurrRemovalEnergy = BindEnergyPerNucleonParametrization(mAr_A, mAr_Z);
 
  *fRemovalE = fCurrRemovalEnergy;

  // Or set Manually for Argon?
  //  https://github.com/GENIE-MC/Generator/blob/master/config/G18_01b/G18_01b_02_11a/CommonParam.xml
    /* <param type="double" name="RFG-NucRemovalE@Pdg=1000030060">  0.0170 </param> <!-- Li6    --> */
    /* <param type="double" name="RFG-NucRemovalE@Pdg=1000060120">  0.0250 </param> <!-- C12    --> */
    /* <param type="double" name="RFG-NucRemovalE@Pdg=1000080160">  0.0270 </param> <!-- O16    --> */
    /* <param type="double" name="RFG-NucRemovalE@Pdg=1000120240">  0.0320 </param> <!-- Mg24   --> */
    /* <param type="double" name="RFG-NucRemovalE@Pdg=1000180400">  0.0295 </param> <!-- Ar40   --> */
    /* <param type="double" name="RFG-NucRemovalE@Pdg=1000200400">  0.0280 </param> <!-- Ca40   --> */
    /* <param type="double" name="RFG-NucRemovalE@Pdg=1000220480">  0.0300 </param> <!-- Ti48  --> */
    /* <param type="double" name="RFG-NucRemovalE@Pdg=1000260560">  0.0360 </param> <!-- Fe56   --> */
    /* <param type="double" name="RFG-NucRemovalE@Pdg=1000280580">  0.0360 </param> <!-- Ni58   --> */
    /* <param type="double" name="RFG-NucRemovalE@Pdg=1000822080">  0.0440 </param> <!-- Pb208  --> */

  // *fRemovalE = 0.0295;
  
  TLorentzVector nucl;
  TParticlePDG * prot = fPDGdb->GetParticle(2212);
  TParticlePDG * neut = fPDGdb->GetParticle(2112);  
  TParticlePDG * nucleon = neut;

  nucl.SetPx(fCurrMomentum.X());
  nucl.SetPy(fCurrMomentum.Y());
  nucl.SetPz(fCurrMomentum.Z());
  nucl.SetE(TMath::Sqrt(nucleon->Mass()*nucleon->Mass() + nucl.Px()*nucl.Px() + nucl.Py()*nucl.Py() + nucl.Pz()*nucl.Pz()));
  
  return nucl;
  
}



void generatePS_NucleonFermiMotion_Weighted(){

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

  TH1D * probDistBR = setupBodekRitchieProbDistro();
  
  
  double mBeam = 2.0; // GeV neutrino


  TH1F *hPt_Muon = new TH1F("hPt_muon","MuonPt", 200,0.0,2.0);
  TH1F *hPz_Muon = new TH1F("hPz_muon","MuonPz", 200,-2.0,2.1);
  TH1F *hE_Muon = new TH1F("hE_muon","MuonE", 200,0,2.1);
  TH1F *hPhi_Muon = new TH1F("hPhi_muon","MuonPhi", 200,0,kPi);
  TH1F *hQ2 = new TH1F("hQ2","Q2", 200,0,2);
  TH1F *hEnuReco = new TH1F("hEnuReco","Enu Reco", 200,0,4);

  Int_t Nevents = 1e+4;
  for (Int_t nEv=0;nEv<Nevents;nEv++) {
  
    // Here we get first the final-state nucleon 4-mom
    // after fermi motion


    double fRemovalE = 0;
    double * pfRemovalE = &fRemovalE;
    // Initial nucleon 4-momentum (nucleus rest frame) and removal energy from Bodek Ritchie
    TLorentzVector pNi = generateNucleonBodekRitchie(probDistBR, pfRemovalE);
    TVector3 p3Ni = pNi.Vect();
    
    //---------------------------------------------
    // Add Removal energy (from Genie) for the initial state nucleon
    // Look up the (on-shell) mass of the initial nucleon
    double mNi = mNeutron;
    double mNf = mProton;
    // Initial nucleus mass
    double Mi = mAr;
     // Final nucleus mass 
    double Mf = 0.; 

    double Eb = fRemovalE;
    // This equation is the definition that we assume
    // here for the "removal energy" (Eb) returned by the
    // nuclear model. It matches GENIE's convention for
    // the Bodek/Ritchie Fermi gas model.
    Mf = Mi + Eb - mNi;
    // The (lab-frame) off-shell initial nucleon energy is the difference
    // between the lab frame total energies of the initial and remnant nuclei
    double ENi = Mi - std::sqrt( Mf*Mf + p3Ni.Mag2() );
    // Update the initial nucleon lab-frame 4-momentum in the interaction with
    // its current components
    pNi.SetVect( p3Ni );
    pNi.SetE( ENi );

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
    TVector3 beta_nRest_to_lab = (pNi).BoostVector();// pNi is the nucleon momentum from Bodek-Ritchie model
    TVector3 beta_lab_to_nRest = -beta_nRest_to_lab;

    // Also get the Lorentz boost beta for the beam+target COM frame
    // in which the outgoing lepton and outgoing nucleon
    // are back-to-back (to be used later)
    TVector3 beta_beamtargetCOM_to_lab = (pNi+beam).BoostVector();
    TVector3 beta_lab_to_beamtargetCOM = -beta_beamtargetCOM_to_lab;

    // Let's boost from the Lab frame to the Nucleon at rest frame
    // where we can evaluate the LL-Smith cross section
    TLorentzVector beam_nRest = TLorentzVector( beam );
    beam_nRest.Boost( beta_lab_to_nRest );
    TLorentzVector nucleon_nRest = TLorentzVector( pNi );
    nucleon_nRest.Boost( beta_lab_to_nRest );
    
    double Enu_Lab = beam.E();
    double Enu = beam_nRest.E();

    // Check
    /* std::cout << "nu Lab: " << std::endl; */
    /* beam.Print(); */
    /* std::cout << ", nu COM: " << std::endl; */
    /* beamCOM.Print(); */
    /* std::cout << "Nucleon in Lab: "  << std::endl; */
    /* pNf.Print(); */
    /* std::cout << "Nucleon in COM: " << std::endl; */
    /* nucleonCOM.Print(); */
    /* std::cout << "-----------------------------------------" << std::endl; */


    
    // Generate the random PS in frame in which the nucleon is at rest
    TLorentzVector target(0.0, 0.0, 0.0, mNeutron);// GeV

    // The beam neutrino is already boosted into this frame
    TLorentzVector W = beam_nRest + target;

    double lepMass = mMuon;// mMuon or mElectron or mTau
    Double_t masses[2] = { mProton , lepMass};

    // Generate the Decay in then nucleon at rest frame
    TGenPhaseSpace event;
    //  event.SetDecay(W, 2, masses, "Fermi");
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

    // We need the Fermi momentum
    double kF = FermiMomentumForIsoscalarNucleonParametrization(mAr_A);

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
    TLorentzVector outNucleon(-1*outlepton_COM.Px(),-1*outlepton_COM.Py(),-1*outlepton_COM.Pz(),
			      TMath::Sqrt(outlepton_COM.P()*outlepton_COM.P() + mNf*mNf));

    // Now boost the outgoing nucleon into the Lab frame
    outNucleon.Boost(beta_beamtargetCOM_to_lab);

    // Apply Pauli blocking
    if(outNucleon.P() < kF){
      //std::cout << "Fermi Momentum: " << kF << ", Nucleon momentum: " << outNucleon.P() << std::endl;    
      continue;
    }

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
    double Q2 = 0;// dummy
    double * Q2p = &Q2; // dummy
    double xsec = DSigma_dQ2_XSec(Enu, El, Costheta, phi, Q2p);   // Q2p is dummy
    double R = NuclQELXSecSuppression(mAr_A, mAr_Z, -Q2, mNi,  0.5);
    //  std::cout << "Q2: " <<  Q2 << ", Supression: " << R << std::endl;
    xsec = R*xsec;
    
    //--------------------------------------------------
    // Now, we we apply the Rejection Sampling algorithm

    // Now we need to find the Maximum diff. cross section that covers this 
    // function
    const int N_theta = 100;
    const int N_phi = 100;
    double costh_range_min = -1.; 
    double costh_range_max = 1;
    double phi_range_min = 0.;
    double phi_range_max = 2*TMath::Pi();
    double costh_increment = (costh_range_max-costh_range_min) / N_theta;
    double phi_increment   = (phi_range_max-phi_range_min) / N_phi;

    double phi_at_xsec_max = -1;
    double costh_at_xsec_max = 0;
    
    double xsec_max = 0;
    
    // Now scan through angles coarsely with Enu and El fixed
    for (int itheta = 0; itheta < N_theta; itheta++){
      double costh_ = costh_range_min + itheta * costh_increment;
      for (int iphi = 0; iphi < N_phi; iphi++) { // Scan around phi
	double phi_ = phi_range_min + iphi * phi_increment;
	double Q2_ = 0;// dummy
	double * Q2p_ = &Q2_; // dummy
	double xs = DSigma_dQ2_XSec(Enu, El, costh_, phi_, Q2p_);   // Q2p is dummy
	R = NuclQELXSecSuppression(mAr_A, mAr_Z, -Q2_, mNi,  0.5);
	//  std::cout << "Q2: " <<  Q2 << ", Supression: " << R << std::endl;
	xs = R*xs;
	
	if (xs > xsec_max){
	  phi_at_xsec_max = phi_;
	  costh_at_xsec_max = costh_;
	  xsec_max = xs;
	}
	
      } // Done with phi scan
    }// Done with angles coarsely
      

    // std::cout << " -------------------- best estimate for xsec_max = " << xsec_max*perGeV2_to_b * barn_to_cm2 << std::endl;

    xsec_max *= SafetyFactor;

    
    // Now we apply the acceptance-reject method criteria
    double t = xsec_max * ran->Rndm();

    bool  accept = (t < xsec);
    

    if(accept){
      std::cout << "\r " << nEv << ", RemovalE: " << Eb <<  std::flush;
      
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
      double w_ = 1.0 ;//weight;
      hPt_Muon->Fill(pt_muon, w_); 
      hPz_Muon->Fill(pz_muon, w_);
      hPhi_Muon->Fill(phi_muon, w_);
      hE_Muon->Fill(pMuon_Lab.E(), w_);
      hQ2->Fill(Q2Lab, w_);


      // Try to reconstruct the neutrino energy
      double Eb = 0;// BindEnergyPerNucleonParametrization(mAr_A, mAr_Z);//0.0295; 
      double EnuReco =0.5*( mProton*mProton - TMath::Power((mNeutron - Eb), 2) + lepMass*lepMass + 2*(mNeutron - Eb)*pMuon_Lab.E())/ (mNeutron - Eb - pMuon_Lab.E()+ p_muon*TMath::Cos(Theta_Muon_Nu));
      hEnuReco->Fill(EnuReco, w_);
    }

  }
  std::cout << std::endl;

  TCanvas* c = new TCanvas("c", "c", 1);
  c->SetTitle("TGenPhaseSpace QEL CC nu on Argon");
  c->Divide(2,2);
  c->cd(1);
  hPt_Muon->Draw("hist");
  c->cd(2);
  hPz_Muon->Draw("hist");
  c->cd(3);
  hQ2->Draw("hist");
  c->cd(4);
  hE_Muon->Draw("hist");

  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  hEnuReco->Draw("hist");
  

}

