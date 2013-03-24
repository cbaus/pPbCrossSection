#include "CastorTreeVariables.h"

#include <cmath>
#include <iostream>
#include <cstdlib>

// #define PI 3.14159265 // this is in cmath M_PI

using namespace std;

ClassImp(CentralityData)
ClassImp(RecHit)
ClassImp(RecHitEB)
ClassImp(RecHitEE)
ClassImp(RecHitHBHE)
ClassImp(RecHitHF)
ClassImp(TowerHF)
ClassImp(RecHitZDC)
ClassImp(DigiZDC)
ClassImp(DigiCASTOR)
ClassImp(RecHitCASTOR)
ClassImp(TrackData)
ClassImp(GenParticle)

ClassImp(AnalysisEvent)



CentralityData::CentralityData()
: cenBin(0),
  cenNpartic(0),
  cenNpartSigma(0),
  cenNcoll(0),
  cenNcollSigma(0),
  cenBpar(0),
  cenBsigma(0)
{
}

RecHit::RecHit()
: Energy(0),
  Phi(0),
  Eta(0)
{
}

double
RecHit::GetEt() const
{
  return Energy / std::cosh(Eta);
}


RecHitEB::RecHitEB()
: Iphi(0),
  Ieta(0),
  Time(0),
  Severity(0),
  SwissCrx(0)
{
}


RecHitEE::RecHitEE()
: Ix(0),
  Iy(0),
  Time(0),
  Severity(0),
  SwissCrx(0)
{
}

RecHitHBHE::RecHitHBHE()
: Iphi(0),
  Ieta(0),
  Idepth(0),
  Saturation(0)
{
}

RecHitHF::RecHitHF()
: Iphi(0),
  Ieta(0),
  Idepth(0),
  Saturation(0),
  PMThitLS(0),
  PMThitTime(0)
{
}

TowerHF::TowerHF()
: Energy(0), //// l - s ~ e(0), l + s ~ h + e(0), h ~ 2 s ? // 900
  Pt(0),
  Phi(0),
  Eta(0),
  IetaAbs(0)
{
}

RecHitZDC::RecHitZDC()
: Isection(0), //if(zside==-1) Detector=ZDCM(0), if(section==1) SubDetector=EM(0), 2=HAD
  Ichannel(0), //add saturation bits etc
  Saturation(0)
{
}

DigiZDC::DigiZDC()
: Iside(0),
  Isection(0),
  Ichannel(0),
  CapIdTs0(0)
{
  for (int i=0; i<10; ++i) {
    FC[i] = 0;
    ADC[i] = 0;
  }
}


DigiCASTOR::DigiCASTOR()
: Sector(0),
  Module(0),
  CapIdTs0(0)
{
  for (int i=0; i<10; ++i) {
    FC[i] = 0;
    ADC[i] = 0;
  }
}

RecHitCASTOR::RecHitCASTOR()
: Time(0),
  //  Xpos(0),
  //  Ypos(0),
  //  Zpos(0),
  Sector(0),
  Module(0)
{
}

TrackData::TrackData()
: Pt(0),
  Eta(0)
{
}

GenParticle::GenParticle()
: Id(0),
  Status(0),
  Px(0),
  Py(0),
  Pz(0),
  fProp(0)
{
}

const ParticleInfo::Properties&
GenParticle::GetProperties() const
{
  // don't use if (!fProp). Not always set to 0 when reading tree
  fProp = &ParticleInfo::GetInstance().Get(Id);
  return *fProp;
}

double
GenParticle::GetCosTheta() const
{
  return Pz/GetMomentum();
}

double
GenParticle::GetTheta() const
{
  return acos(GetCosTheta());
}


double
GenParticle::GetEta() const
{
  const double tt = tan(GetTheta()/2);
  return (tt<0 ? 1 : -1 ) * log(abs(tt));
}

double
GenParticle::GetRapidity() const
{
  const double E = GetEnergy();
  return 0.5 * log((E+Pz)/(E-Pz));
}

double
GenParticle::GetCharge() const
{
  return GetProperties().GetCharge();
}

double
GenParticle::GetMass() const
{
  return GetProperties().GetMass();
}

double
GenParticle::GetMomentum() const
{
  return sqrt(Px*Px + Py*Py + Pz*Pz);
}

double
GenParticle::GetEnergy() const
{
  const double P = GetMomentum();
  const double M = GetMass();
  return sqrt(P*P + M*M);
}

double
GenParticle::GetPhi() const
{
  const double Phi = atan2(Px,Py);
  return Phi;
}

double
GenParticle::GetTransverseEnergy() const
{
  const double E = GetEnergy();
  const double Eta = abs(GetEta());
  return E / std::cosh(Eta);
}

double
GenParticle::GetTransverseMomentum() const
{
  return sqrt(Px*Px + Py*Py);
}



AnalysisEvent::AnalysisEvent()
: runNb (0),
  evtNb (0),
  bxNb (0),
  orbitNb (0),
  lumiNb (0),
  // centrality
  cenHf(0),
  cenHft(0),
  cenHftp(0),
  cenHftm(0),
  cenEb(0),
  cenEe(0),
  cenEep(0),
  cenEem(0),
  cenZdc(0),
  cenZdcm(0),
  cenZdcp(0),
  cenNpix(0),
  cenEt(0),
  // event plane
  hasEventPlane(0),
  evtPlaneTracksMidEta(0),
  evtPlaneHF(0),
  hasDetData(0),
  hbheNoiseFilter(0),
  nbCasCalPresamples(0), // for the digis
  nVertex(0),
  vertexX(0),
  vertexY(0),
  vertexZ(0),
  isMC(0),
  genImpactParameter(0),
  genNpart(0)
{
  for (int i=0; i<200; ++i) {
    hltAccept[i] = 0;
    hltPrescl[i] = 0;
  }

  for (int i=0; i<128; ++i) {
    algDecision[i] = 0;
    algDecisionB4mask[i] = 0;
    algPrescale[i] = 0;
  }

  for (int i=0; i<64; ++i) {
    ttDecision[i] = 0;
    ttDecisionB4mask[i] = 0;
    ttPrescale[i] = 0;
  }

  /*
    algResultBeforeMask[128](0),
    algResultAfterMask[128](0),
    algTriggerMask[128](0),
  */
}

bool
AnalysisEvent::HasCentrality(const std::string& id) const
{
  for (vector<CentralityData>::const_iterator iC = Centrality.begin(); iC != Centrality.end(); ++iC) {
    if (iC->id == id)
      return true;
  }
  return false;
}


const CentralityData&
AnalysisEvent::GetCentrality(const std::string& id) const
{
  for (vector<CentralityData>::const_iterator iC = Centrality.begin(); iC != Centrality.end(); ++iC) {
    if (iC->id == id)
      return *iC;
  }
  cout << "Did not find centrality data: " << id << ". Use HasCentrality() next time!" << endl;
  exit(1);
}



