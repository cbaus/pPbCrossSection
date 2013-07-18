#include <ParticleInfo.h>

#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <cstdlib>

using namespace std;


ParticleInfo::ParticleInfo()
  : fErrorIfNotFound(true)
{
  fstream theFile;
  theFile.open ("./PythiaParticleTable.dat", std::ios::in);
  if (!theFile.is_open ())
    throw runtime_error ("cannot open pythia file");

  string line,element;
  while(getline ( theFile, line) )
    {
      int id;
      string name;
      string antiname;
      int charge;
      double mass;
      double lifetime;
      stringstream ssline(line);
      //KF,KC,particle,antiparticle,chg,col,anti,mass,width,w-cut,lifetime,decay

      getline(ssline, element,','); stringstream sselement(element); sselement >> id;
      getline(ssline, element,',');
      getline(ssline, element,','); name = element;
      getline(ssline, element,','); antiname = element;
      getline(ssline, element,','); stringstream sselement2(element); sselement2 >> charge; charge /= 3;
      getline(ssline, element,',');
      getline(ssline, element,',');
      getline(ssline, element,','); stringstream sselement3(element);sselement3 >> mass;
      getline(ssline, element,',');
      getline(ssline, element,',');
      getline(ssline, element,','); stringstream sselement4(element); sselement4 >> lifetime;
      getline(ssline, element,',');

      if (fParticle.count(id) != 0) {
	ostringstream ss;
	ss << "Particle with id=" << id << " occured multiple times!";
	throw runtime_error(ss.str());
      }

      fParticle[id] = Properties(mass, lifetime, charge, name);

      if (antiname != "")
	fParticle[-id] = Properties(mass, lifetime, -charge, antiname);

    }
  theFile.close ();

  fBaryons.clear();
  fBaryons.insert(2212);
  fBaryons.insert(2112);
  fBaryons.insert(2224);
  fBaryons.insert(2214);
  fBaryons.insert(2114);
  fBaryons.insert(1114);
  fBaryons.insert(3122);
  fBaryons.insert(3222);
  fBaryons.insert(3212);
  fBaryons.insert(3112);
  fBaryons.insert(3224);
  fBaryons.insert(3214);
  fBaryons.insert(3114);
  fBaryons.insert(3322);
  fBaryons.insert(3312);
  fBaryons.insert(3324);
  fBaryons.insert(3314);
  fBaryons.insert(3334);
  fBaryons.insert(4122);
  fBaryons.insert(4222);
  fBaryons.insert(4212);
  fBaryons.insert(4112);
  fBaryons.insert(4224);
  fBaryons.insert(4214);
  fBaryons.insert(4114);
  fBaryons.insert(4232);
  fBaryons.insert(4132);
  fBaryons.insert(4322);
  fBaryons.insert(4312);
  fBaryons.insert(4324);
  fBaryons.insert(4314);
  fBaryons.insert(4332);
  fBaryons.insert(4334);
  fBaryons.insert(4412);
  fBaryons.insert(4422);
  fBaryons.insert(4414);
  fBaryons.insert(4424);
  fBaryons.insert(4432);
  fBaryons.insert(4434);
  fBaryons.insert(4444);
  fBaryons.insert(5122);
  fBaryons.insert(5112);
  fBaryons.insert(5212);
  fBaryons.insert(5222);
  fBaryons.insert(5114);
  fBaryons.insert(5214);
  fBaryons.insert(5224);
  fBaryons.insert(5132);
  fBaryons.insert(5232);
  fBaryons.insert(5312);
  fBaryons.insert(5322);
  fBaryons.insert(5314);
  fBaryons.insert(5324);
  fBaryons.insert(5332);
  fBaryons.insert(5334);
  fBaryons.insert(5142);
  fBaryons.insert(5242);
  fBaryons.insert(5412);
  fBaryons.insert(5422);
  fBaryons.insert(5414);
  fBaryons.insert(5424);
  fBaryons.insert(5342);
  fBaryons.insert(5432);
  fBaryons.insert(5434);
  fBaryons.insert(5442);
  fBaryons.insert(5444);
  fBaryons.insert(5512);
  fBaryons.insert(5522);
  fBaryons.insert(5514);
  fBaryons.insert(5524);
  fBaryons.insert(5532);
  fBaryons.insert(5534);
  fBaryons.insert(5542);
  fBaryons.insert(5544);
  fBaryons.insert(5554);

  fMesons.insert(111);
  fMesons.insert(211);
  fMesons.insert(9000111);
  fMesons.insert(9000211);
  fMesons.insert(100111);
  fMesons.insert(100211);
  fMesons.insert(10111);
  fMesons.insert(10211);
  fMesons.insert(9010111);
  fMesons.insert(9010211);
  fMesons.insert(113);
  fMesons.insert(213);
  fMesons.insert(10113);
  fMesons.insert(10213);
  fMesons.insert(20113);
  fMesons.insert(20213);
  fMesons.insert(9000113);
  fMesons.insert(9000213);
  fMesons.insert(100113);
  fMesons.insert(100213);
  fMesons.insert(9010113);
  fMesons.insert(9010213);
  fMesons.insert(9020113);
  fMesons.insert(9020213);
  fMesons.insert(30113);
  fMesons.insert(30213);
  fMesons.insert(9030113);
  fMesons.insert(9030213);
  fMesons.insert(9040113);
  fMesons.insert(9040213);
  fMesons.insert(115);
  fMesons.insert(215);
  fMesons.insert(10115);
  fMesons.insert(10215);
  fMesons.insert(9000115);
  fMesons.insert(9000215);
  fMesons.insert(9010115);
  fMesons.insert(9010215);
  fMesons.insert(117);
  fMesons.insert(217);
  fMesons.insert(9000117);
  fMesons.insert(9000217);
  fMesons.insert(9010117);
  fMesons.insert(9010217);
  fMesons.insert(119);
  fMesons.insert(219);
  fMesons.insert(221);
  fMesons.insert(331);
  fMesons.insert(9000221);
  fMesons.insert(9010221);
  fMesons.insert(100221);
  fMesons.insert(10221);
  fMesons.insert(9020221);
  fMesons.insert(100331);
  fMesons.insert(9030221);
  fMesons.insert(10331);
  fMesons.insert(9040221);
  fMesons.insert(9050221);
  fMesons.insert(9060221);
  fMesons.insert(9070221);
  fMesons.insert(9080221);
  fMesons.insert(223);
  fMesons.insert(333);
  fMesons.insert(10223);
  fMesons.insert(20223);
  fMesons.insert(10333);
  fMesons.insert(20333);
  fMesons.insert(100223);
  fMesons.insert(9000223);
  fMesons.insert(9010223);
  fMesons.insert(30223);
  fMesons.insert(100333);
  fMesons.insert(225);
  fMesons.insert(9000225);
  fMesons.insert(335);
  fMesons.insert(9010225);
  fMesons.insert(9020225);
  fMesons.insert(10225);
  fMesons.insert(9030225);
  fMesons.insert(10335);
  fMesons.insert(9040225);
  fMesons.insert(9050225);
  fMesons.insert(9060225);
  fMesons.insert(9070225);
  fMesons.insert(9080225);
  fMesons.insert(9090225);
  fMesons.insert(227);
  fMesons.insert(337);
  fMesons.insert(229);
  fMesons.insert(9000229);
  fMesons.insert(9010229);
  fMesons.insert(130);
  fMesons.insert(310);
  fMesons.insert(311);
  fMesons.insert(321);
  fMesons.insert(9000311);
  fMesons.insert(9000321);
  fMesons.insert(10311);
  fMesons.insert(10321);
  fMesons.insert(100311);
  fMesons.insert(100321);
  fMesons.insert(9010311);
  fMesons.insert(9010321);
  fMesons.insert(9020311);
  fMesons.insert(9020321);
  fMesons.insert(313);
  fMesons.insert(323);
  fMesons.insert(10313);
  fMesons.insert(10323);
  fMesons.insert(20313);
  fMesons.insert(20323);
  fMesons.insert(100313);
  fMesons.insert(100323);
  fMesons.insert(9000313);
  fMesons.insert(9000323);
  fMesons.insert(30313);
  fMesons.insert(30323);
  fMesons.insert(315);
  fMesons.insert(325);
  fMesons.insert(9000315);
  fMesons.insert(9000325);
  fMesons.insert(10315);
  fMesons.insert(10325);
  fMesons.insert(20315);
  fMesons.insert(20325);
  fMesons.insert(9010315);
  fMesons.insert(9010325);
  fMesons.insert(9020315);
  fMesons.insert(9020325);
  fMesons.insert(317);
  fMesons.insert(327);
  fMesons.insert(9010317);
  fMesons.insert(9010327);
  fMesons.insert(319);
  fMesons.insert(329);
  fMesons.insert(9000319);
  fMesons.insert(9000329);
  fMesons.insert(411);
  fMesons.insert(421);
  fMesons.insert(10411);
  fMesons.insert(10421);
  fMesons.insert(413);
  fMesons.insert(423);
  fMesons.insert(10413);
  fMesons.insert(10423);
  fMesons.insert(20413);
  fMesons.insert(20423);
  fMesons.insert(415);
  fMesons.insert(425);
  fMesons.insert(431);
  fMesons.insert(10431);
  fMesons.insert(433);
  fMesons.insert(10433);
  fMesons.insert(20433);
  fMesons.insert(435);
  fMesons.insert(511);
  fMesons.insert(521);
  fMesons.insert(10511);
  fMesons.insert(10521);
  fMesons.insert(513);
  fMesons.insert(523);
  fMesons.insert(10513);
  fMesons.insert(10523);
  fMesons.insert(20513);
  fMesons.insert(20523);
  fMesons.insert(515);
  fMesons.insert(525);
  fMesons.insert(531);
  fMesons.insert(10531);
  fMesons.insert(533);
  fMesons.insert(10533);
  fMesons.insert(20533);
  fMesons.insert(535);
  fMesons.insert(541);
  fMesons.insert(10541);
  fMesons.insert(543);
  fMesons.insert(10543);
  fMesons.insert(20543);
  fMesons.insert(545);
  fMesons.insert(441);
  fMesons.insert(10441);
  fMesons.insert(100441);
  fMesons.insert(443);
  fMesons.insert(10443);
  fMesons.insert(20443);
  fMesons.insert(100443);
  fMesons.insert(30443);
  fMesons.insert(9000443);
  fMesons.insert(9010443);
  fMesons.insert(9020443);
  fMesons.insert(445);
  fMesons.insert(100445);
  fMesons.insert(551);
  fMesons.insert(10551);
  fMesons.insert(100551);
  fMesons.insert(110551);
  fMesons.insert(200551);
  fMesons.insert(210551);
  fMesons.insert(553);
  fMesons.insert(10553);
  fMesons.insert(20553);
  fMesons.insert(30553);
  fMesons.insert(100553);
  fMesons.insert(110553);
  fMesons.insert(120553);
  fMesons.insert(130553);
  fMesons.insert(200553);
  fMesons.insert(210553);
  fMesons.insert(220553);
  fMesons.insert(300553);
  fMesons.insert(9000553);
  fMesons.insert(9010553);
  fMesons.insert(555);
  fMesons.insert(10555);
  fMesons.insert(20555);
  fMesons.insert(100555);
  fMesons.insert(110555);
  fMesons.insert(120555);
  fMesons.insert(200555);
  fMesons.insert(557);
  fMesons.insert(100557);
}

ParticleInfo::~ParticleInfo() {}


const ParticleInfo::Properties&
ParticleInfo::Get(int id) const {
  map<int,Properties>::const_iterator theIter = fParticle.find(id);
  if (fErrorIfNotFound && theIter == fParticle.end ()) {
    cout << "PDG Id not found: " << id << endl;
    /*
      cout << "Dumping table n=" << fParticle.size() << endl;
      for (map<int,Properties>::const_iterator theIter = fParticle.begin(); theIter!=fParticle.end(); ++theIter)
      cout << theIter->first << " " << theIter->second.GetName() << endl;
    */

    return *(new Properties());
  }
  return theIter->second;
}

bool
ParticleInfo::AddCollection(const set<int>& theSet, int number)
{
  if (number == -1)
    {
      fSpecialCollection.push_back (theSet);
      return true;
    }
  else if (number < 0)
    {
      return false;
    }
  else
    {
      if (int(fSpecialCollection.size ()) < number)
	fSpecialCollection.resize (number);
      fSpecialCollection[number] = theSet;
      return true;
    }
}

bool
ParticleInfo::IsInCollection(int id, int number) const
{
  if (number < -1)
    {
      return false;
    }
  else
    {
      set<int> collection;
      if (number == -1)
	collection = fSpecialCollection.back ();
      else if (number < 0 || number > int(fSpecialCollection.size ()))
	return false;
      else
	collection = fSpecialCollection[number];

      return !(collection.find (id) == collection.end ());
    }
}

