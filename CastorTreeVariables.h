#ifndef _include_CastorTreeVariables_h_
#define _include_CastorTreeVariables_h_

#include <TObject.h> // for ROOT macros

#include <map>
#include <vector>
#include <string>

#include "ParticleInfo.h"

namespace ForwardRecord {
  const unsigned int CSectors = 16;
  const unsigned int CModules = 14;

  const unsigned int nMaxEEMRecHits = 14406;
  const unsigned int nMaxHBHERecHits = 5183;
  const unsigned int nMaxEBMRecHits = 60555;
  const unsigned int nMaxHFMRecHits = 1728;
  // const int nMaxHFRecHits = 900;
  const unsigned int nMaxZDCRecHits = 18;
  const unsigned int nMaxCASTORRecHits = CSectors*CModules;
} // end namespace ForwardRecord


// namespace kit {

/* ***********************************************************
   class for Centrality Data
 */
class CentralityData {

 public:
  CentralityData();

  std::string id; // centrality id
  int cenBin;
  float cenNpartic;
  float cenNpartSigma;
  float cenNcoll;
  float cenNcollSigma;
  float cenBpar;
  float cenBsigma;

  ClassDef(CentralityData, 3);
};

/* ***********************************************************
   base class for all RecHit data
 */
class RecHit {
 public:
  RecHit();

  float Energy;
  float Phi;
  float Eta;

  double GetEt() const;

  ClassDef(RecHit, 1);
};


/* ***********************************************************
   class for EB RecHit data
 */
class RecHitEB : public RecHit {
 public:
  RecHitEB();

  int    Iphi;// do we need this?
  int    Ieta;// do we need this?
  float Time; // do we need this?
  int    Severity;// do we need this?
  float SwissCrx;// do we need this?

  ClassDef(RecHitEB, 1);
};

/* ***********************************************************
   class for EE RecHit data
 */
class RecHitEE : public RecHit {
 public:
  RecHitEE();

  int Ix;// do we need this?
  int Iy;// do we need this?
  float Time;// do we need this?
  int Severity;// do we need this?
  float SwissCrx;// do we need this?

  ClassDef(RecHitEE, 2);
};

/* ***********************************************************
   class for HBHE RecHit data
 */
class RecHitHBHE : public RecHit {
 public:
  RecHitHBHE();

  int Iphi;// do we need this?
  int Ieta;// do we need this?
  int Idepth;// do we need this?
  int Saturation;// do we need this?

  ClassDef(RecHitHBHE, 1);
};


/* ***********************************************************
   class for HF RecHit data
 */
class RecHitHF : public RecHit {
 public:
  RecHitHF();

  int Iphi;// do we need this?
  int Ieta;// do we need this?
  int Idepth;// do we need this?
  int Saturation;// do we need this?
  int PMThitLS;// do we need this?
  int PMThitTime;// do we need this?

  ClassDef(RecHitHF, 1);
};

/* ***********************************************************
   class for HF RecHit data

   in general 36 wedges x 13 eta-rings;
   but need only HF- => 18 wedges;
   each ring of one wedge = either one (for 2 rings closest to beam pipe)
   or two towers => in HF- expect 18 x (2 + 11 x 2) = 432 towers
   1st eta-ring overlaps with HE ? => exclude it later ?
 */
class TowerHF {
 public:
  TowerHF();

  float Energy; //// l - s ~ e; l + s ~ h + e; h ~ 2 s ? // 900
  float Pt;
  float Phi;
  float Eta;
  int IetaAbs;// do we need this?

  ClassDef(TowerHF, 1);
};

/* ***********************************************************
   class for ZDC RecHit data
 */
class RecHitZDC : public RecHit {
 public:
  RecHitZDC();

  // do we need these?
  int Isection; //if(zside==-1) Detector=ZDCM; if(section==1) SubDetector=EM; 2=HAD
  int Ichannel; //add saturation bits etc
  int Saturation;

  ClassDef(RecHitZDC, 1);
};

/* ***********************************************************
   class for ZDC Digi data
 */
class DigiZDC {
 public:
  DigiZDC();

  // do we need zdc digis in general?
  float FC[10];
  float ADC[10];
  int Iside;
  int Isection;
  int Ichannel;
  int CapIdTs0;

  ClassDef(DigiZDC, 1);
};


/* ***********************************************************
   class for CASTOR Digi data
 */
class DigiCASTOR {
 public:
  DigiCASTOR();

  // do we need castor digis in general?
  float FC[10];
  float ADC[10];
  int Sector;
  int Module;
  int CapIdTs0;

  ClassDef(DigiCASTOR, 1);
};

/* ***********************************************************
   class for CASTOR RecHit data
 */
class RecHitCASTOR : public RecHit {
 public:
  RecHitCASTOR();

  int GetModuleId() const { return Module-1; }
  int GetSectorId() const { return Sector-1; }

  float Time; // do we need this?
  //  float Xpos; // do we need this? If we need put it in static detector description!
  //  float Ypos; // do we need this?
  //  float Zpos; // do we need this?
  int Sector; // do we need this? // use single index instead !
  int Module; // do we need this?

  ClassDef(RecHitCASTOR, 2);
};


/* ***********************************************************
   class to store reco track data
 */
class TrackData {
 public:
  TrackData();

  // do we need tracks? (no...)
  float Pt;
  float Eta;

  ClassDef(TrackData, 1);
};


/* ***********************************************************
   class to store gen particles
 */
class GenParticle {
 public:
  GenParticle();

  // simplify content:
  int Id;
  int Status;
  float Px;
  float Py;
  float Pz;

  mutable ParticleInfo::Properties const * fProp; //! for caching !
  const ParticleInfo::Properties& GetProperties() const;

  double GetCosTheta() const;
  double GetTheta() const;
  double GetEta() const;
  double GetRapidity() const;
  double GetCharge() const;
  double GetMass() const;
  double GetMomentum() const;
  double GetEnergy() const;
  double GetPhi() const;
  double GetTransverseEnergy() const;
  double GetTransverseMomentum() const;

  ClassDef(GenParticle, 1);
};





/* ***********************************************************
   main analysis event class
 */
class AnalysisEvent {

 public:
  AnalysisEvent();

  // event variables
  unsigned int runNb ;
  unsigned int evtNb ;
  unsigned int bxNb ;
  unsigned int orbitNb ;
  unsigned int lumiNb ;

  // centralities
  float cenHf;
  float cenHft;
  float cenHftp;
  float cenHftm;
  float cenEb;
  float cenEe;
  float cenEep;
  float cenEem;
  float cenZdc;
  float cenZdcm;
  float cenZdcp;
  float cenNpix;
  float cenEt;
  std::vector<CentralityData> Centrality; // map to keep different centralities, if needed
  const CentralityData& GetCentrality(const std::string& id) const;
  bool HasCentrality(const std::string& id) const;

  // event plane
  bool hasEventPlane;
  float evtPlaneTracksMidEta;
  float evtPlaneHF;

  // rechits and towser + digis
  bool hasDetData; // flag do indicate if there is any detector level data stored
  std::vector<RecHitEB> EB;
  std::vector<RecHitEE> EE;
  int hbheNoiseFilter;
  std::vector<RecHitHBHE> HBHE;
  std::vector<RecHitHF> HF;
  std::vector<TowerHF> HFtowers;
  std::vector<RecHitZDC> ZDC;
  std::vector<DigiZDC> ZDCdigis;
  int nbCasCalPresamples; // for the digis
  std::vector<DigiCASTOR> CASTORdigis;
  std::vector<RecHitCASTOR> CASTOR;

  // track list
  std::vector<TrackData> Tracks;

  // vertex location
  int nVertex;
  float vertexX;
  float vertexY;
  float vertexZ;

  // --------- Monte Carlo information ------------
  bool   isMC;
  float  genImpactParameter;
  int    genNpart;
  int    genProcessID;
  std::vector<GenParticle> GEN;

  // ------------ Trigger information ---------------
  int hltAccept[200];
  int hltPrescl[200];
  //
  /*
    int algResultBeforeMask[128];
    int algResultAfterMask[128];
    int algTriggerMask[128];
  */
  int algDecision[128];
  int algDecisionB4mask[128];
  int algPrescale[128];
  //
  int ttDecision[64];
  int ttDecisionB4mask[64];
  int ttPrescale[64];

  // --------- Luminosity ------------
  float instLuminosity;
  float instLuminosityCorr;
  
  ClassDef(AnalysisEvent, 8);
};

//}


#endif
