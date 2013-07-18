#ifndef _include_PARTICLEINFO__H__
#define _include_PARTICLEINFO__H__

#include <string>
#include <map>
#include <vector>
#include <set>


/************************ParticleInfo************************/
// This is a SINGLETON
//
class ParticleInfo {

 public:
 
  ParticleInfo();
  ~ParticleInfo();
  ParticleInfo(const ParticleInfo&);
  ParticleInfo& operator=(const ParticleInfo&);

  static
    ParticleInfo&
    GetInstance()
#ifdef __MAKECINT__
    ;
#else
    {
      static ParticleInfo instance;
      return instance;
    }
#endif
  
    // the class the holds one particle data
    class Properties {
    public:
      Properties() : fMass(0), fLifetime(0), fCharge(0), fName("-") {}
      Properties(const double m, const double t, const int z, const std::string& n)
      : fMass(m), fLifetime(t), fCharge(z), fName(n) {}
      bool IsCharged() const { return fCharge != 0; }
      double GetMass() const { return fMass; }
      double GetLifetime() const { return fLifetime; }
      int GetCharge() const { return fCharge; }
      const std::string& GetName() const { return fName; }
    private:
      double fMass;
      double fLifetime;
      int fCharge;
      std::string fName;
    };

    // fast access:
    const Properties& Get(int id) const;
    
    // deprecated, since slower:
    double GetMass(int id) const { return Get(id).GetMass(); }
    double GetLifetime(int id) const { return Get(id).GetLifetime(); }
    int GetCharge(int id) const { return Get(id).GetCharge(); }
    bool IsCharged(int id) const { return Get(id).IsCharged(); }
    const std::string& GetName(int id) const { return Get(id).GetName(); }
    
    bool IsBaryon    (int id) const { return fBaryons.count(id); }
    bool IsAntiBaryon(int id) const { return fBaryons.count(-id); }
    bool IsMeson     (int id) const { return fMesons.count(id); }
    bool IsAntiMeson (int id) const { return fMesons.count(-id); }
    
    bool AddCollection(const std::set<int>& theSet, int number = -1);
    bool IsInCollection(int id, int number = -1) const;

    void SetErrorIfNotFound(bool error) { fErrorIfNotFound = error; }    


  private:

    bool fErrorIfNotFound;

    std::map<int, Properties> fParticle;
    std::vector< std::set<int> > fSpecialCollection;
    std::set<int> fBaryons;
    std::set<int> fMesons;
    
  };
  /************************ParticleInfo************************/

#endif // #ifndef PARTICLEINFO__H__
