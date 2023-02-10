#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <vector>

#include <TRandom2.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TParticlePDG.h>

TDatabasePDG * fPDGdb = new TDatabasePDG();

// Constants
double au = 0.931494; // GeV
double mAr = 39.948*au;
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


// Nucleon at rest cross section as a function of the
// - neutrino energy
// - fixed outgoing lepton Energy
// - angle between neutrino and outgoing lepton in the Lab frame
// - phi angle of the lepton
// - Q2p is just a pointer to return the calculated Q2
double XSec(double Enu, double El, double Costheta, double Phi, double * Q2p){

  // We  assume to numerically evaluate the LL-Smith cross section
  // In the nucleon at rest frame to simplfiy

  TLorentzVector neutrino(0,0, Enu, Enu) ; // GeV, fixed
  TLorentzVector nucleon (0, 0, 0, mNeutron); // at rest in the Lab frame
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
  TLorentzVector inNucleonMom = nucleon; // this is the 4mom after FermiMotion
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


    // Apply given scaling factor
    xsec *= fXSecScale;

    // Number of scattering centers in the target
    xsec *= 18; // for argon

    
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

/* double GetFermiMomentum(){ */

/*    kF = TMath::Power(3 * kPi2 * numNuc * */
/*       genie::utils::nuclear::Density(radius, A), 1.0/3.0) * hbarc; */
/*   } */

/* } */


void generatePS_NucleonFermiMotion_StandAlone(){

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

  TRandom2 * ran = new TRandom2();
  
  double mBeam = 2.0; // GeV neutrino


  TH1F *hPt_Muon = new TH1F("hPt_muon","MuonPt", 50,0.0,2.0);
  TH1F *hPz_Muon = new TH1F("hPz_muon","MuonPz", 50,-2.0,2.1);
  TH1F *hE_Muon = new TH1F("hE_muon","MuonE", 50,0,2.1);
  TH1F *hPhi_Muon = new TH1F("hPhi_muon","MuonPhi", 50,0,kPi);
  TH1F *hQ2 = new TH1F("hQ2","Q2", 50,0,2);
  TH1F *hEnuReco = new TH1F("hEnuReco","Enu Reco", 50,0,4);
  TH1F * hxsec = new TH1F("hxsec", "xsec", 50,  0, 50);
  TH1F * hmaxxsec = new TH1F("hmaxxsec", "maxxsec", 50, 0, 50);

  
  Int_t Nevents = 1e+04;
  for (Int_t n=0;n<Nevents;n++) {
  
    // Here we get first the nucleon 4-mom
    TLorentzVector pFermi = PFermiMotion();

    // Now we need to boost the neutrino into the
    // rest frame of the nucleon where the LL-Smith is evaluated

    // Beam in the Lab frame
    TLorentzVector beam(0.0,0.0, mBeam, mBeam);//GeV (beam "particle")

    // Get beta for a Lorentz boost from the lab frame to the Nucleon COM frame
    // where the target Nucleon is at rest.
    TVector3 beta_COM_to_lab = (pFermi).BoostVector();
    TVector3 beta_lab_to_COM = -beta_COM_to_lab;

    TLorentzVector beamCOM = TLorentzVector( beam );
    beamCOM.Boost( beta_lab_to_COM );
    TLorentzVector nucleonCOM = TLorentzVector( pFermi );
    nucleonCOM.Boost( beta_lab_to_COM );
    
    double Enu_Lab = beam.E();
    double Enu = beamCOM.E();

    // Check
    /* std::cout << "nu Lab: " << std::endl; */
    /* beam.Print(); */
    /* std::cout << ", nu COM: " << std::endl; */
    /* beamCOM.Print(); */
    /* std::cout << "Nucleon in Lab: "  << std::endl; */
    /* pFermi.Print(); */
    /* std::cout << "Nucleon in COM: " << std::endl; */
    /* nucleonCOM.Print(); */
    /* std::cout << "-----------------------------------------" << std::endl; */
    
    // In the COM frame the target is at rest
    TLorentzVector target(0.0, 0.0, 0.0, mNeutron);// GeV


    // Now we have to
    // - boost to the COM frame of the beam+Target
    // - throw random directions there, assign the 4 -momentum
    // - boost back to the Nucleon at rest frame 
    
    TVector3 beta_beamTargetCOM_to_Nucleonlab = (beamCOM + target).BoostVector();
    TVector3 beta_Nucleonlab_to_beamTargetCOM = -beta_beamTargetCOM_to_Nucleonlab;

    // We have to generate uniform random final states in the (beam+target) COM frame
    // Pick a uniform direction 
    double costheta = ran->Rndm()*2.0 - 1.0; // [-1,1]
    double phi = ran->Rndm()*2*kPi; // [0, 2Pi]
    double lepMass = mMuon;// mMuon or mElectron or mTau outgoing lepton
    double mNf = mProton; // Mass of final nucleon

    TLorentzVector k4 = beamCOM + target;
    double s = k4.Dot(k4);// Mandelstam s is invariant
    double outLeptonEnergy = ( s - mNf*mNf + lepMass*lepMass ) / (2 * std::sqrt(s));
    double outMomentum = TMath::Sqrt(outLeptonEnergy*outLeptonEnergy - lepMass*lepMass);

    // Assign the 3-momentum for the lepton
    TVector3 lepton3Mom(0., 0., outMomentum);
    lepton3Mom.SetTheta( TMath::ACos(costheta) );
    lepton3Mom.SetPhi( phi );

    // Then rotate the lepton 3-momentum so that the old +z direction now
    // points along the COM frame velocity (beta)
    TVector3 zvec(0., 0., 1.);
    TVector3 rot = ( zvec.Cross(beta_beamTargetCOM_to_Nucleonlab ) ).Unit();
    double angle = beta_beamTargetCOM_to_Nucleonlab.Angle( zvec );

    lepton3Mom.Rotate(angle, rot);

    TLorentzVector lepton(lepton3Mom, outLeptonEnergy);

    // Boost the 4-momenta for both particles into the Nucleon at rest frame
    lepton.Boost(beta_beamTargetCOM_to_Nucleonlab);
    
    TLorentzVector *pMuon  = new TLorentzVector(lepton);
    
    double El = pMuon->E(); // lepton energy

    // Evaluate the diff. xsec dsigma/dQ2 in the Nucleon at rest frame for this random throw
    TVector3 beam3V = beamCOM.Vect();
    TVector3 lepton3V = pMuon->Vect();
    double Costheta = TMath::Cos(lepton3V.Angle(beam3V));

    // To get the phi vector around the neutrino (beam, whatever) direction,
    // we rotate the lepton 3-mom so that the old z-axis vector to point along the neutrino direction
    TVector3 rot_axis = ( zvec.Cross(beam3V) ).Unit(); // the rotation axis
    double rot_angle = beam3V.Angle(zvec); // amount of rotation 
    // Rotate the lepton3V
    lepton3V.Rotate(rot_angle, rot_axis);
    // measure the phi now for the rotated lepton momentum vector
    phi = lepton3V.Phi();

    double Q2 = 0;// dummy
    double * Q2p = &Q2; // dummy
    double xsec = XSec(Enu, El, Costheta, phi, Q2p);   // Q2p is dummy
    

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
    
    // Now scan through angles coarsely
    for (int itheta = 0; itheta < N_theta; itheta++){
      double costh_ = costh_range_min + itheta * costh_increment;
      for (int iphi = 0; iphi < N_phi; iphi++) { // Scan around phi
	double phi_ = phi_range_min + iphi * phi_increment;
	double Q2_ = 0;// dummy
	double * Q2p_ = &Q2_; // dummy
	double xs = XSec(Enu, pMuon->E(), costh_, phi_, Q2p_);   // Q2p is dummy
	
	if (xs > xsec_max){
	  phi_at_xsec_max = phi_;
	  costh_at_xsec_max = costh_;
	  xsec_max = xs;
	}
	
      } // Done with phi scan
    }// Done with angles coarsely
      

    // std::cout << " -------------------- best estimate for xsec_max = " << xsec_max*perGeV2_to_b * barn_to_cm2 << std::endl;

    xsec_max *= SafetyFactor;

    
    // Now we apply the acceptance-reject method
    double t = xsec_max * ran->Rndm();

    bool  accept = (t < xsec);
    
    if(accept){

      // Now we need to boost back the particle into the Lab frame.
      pMuon->Boost(beta_COM_to_lab);

      // Need to recalculate Q2 in the Lab frame
      TLorentzVector qP4 = beam - (*pMuon);
      double Q2Lab = -1. * qP4.Mag2();      
      
      double pt_muon = TMath::Sqrt(pMuon->Px()*pMuon->Px() + pMuon->Py()*pMuon->Py());
      double pz_muon = pMuon->Pz();
      double p_muon = TMath::Sqrt(pMuon->Px()*pMuon->Px() + pMuon->Py()*pMuon->Py() + pz_muon*pz_muon);
      double Theta_Muon_Nu = (pMuon->Vect()).Angle(beam.Vect());
      double phi_muon = pMuon->Phi();
      // if we want weighted events we need to multiply with the (weight) xsec
      double w_ = 1.0;//xsec;
      hPt_Muon->Fill(pt_muon,w_); 
      hPz_Muon->Fill(pz_muon,w_);
      hPhi_Muon->Fill(phi_muon,w_);
      hE_Muon->Fill(pMuon->E(),w_);
      hQ2->Fill(Q2Lab, w_);
      hxsec->Fill(xsec *perGeV2_to_b * barn_to_cm2*1e+38);
      hmaxxsec->Fill(xsec_max *perGeV2_to_b * barn_to_cm2*1e+38);
      // std::cout << xsec *perGeV2_to_b * barn_to_cm2*1e+38 << ", " << xsec_max *perGeV2_to_b * barn_to_cm2*1e+38 << std::endl;
      
      // Try to reconstruct the neutrino energy?
      double Eb = 0; // Eb should be included later...
      double EnuReco =0.5*( mProton*mProton - TMath::Power((mNeutron - Eb), 2) + lepMass*lepMass + 2*(mNeutron - Eb)*pMuon->E())/ (mNeutron - Eb - pMuon->E()+ p_muon*TMath::Cos(Theta_Muon_Nu));
      hEnuReco->Fill(EnuReco);
    }
  }

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
  c1->Divide(2,2);
  c1->cd(1);
  hEnuReco->Draw("hist");
  c1->cd(2);
  hxsec->Draw("hist");
  c1->cd(3);
  hmaxxsec->Draw("hist");
}

