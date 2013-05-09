{
const double lumiPerLS_begin = 2. * 1e6;
const double lumiPerLS_end = 0.9 * 1e6;
const double p1_single = 0.958707;
const double p1_double = 0.981198;

double p1 = p1_single;
double lumiPerLS = lumiPerLS_begin;

double pileupcorr = 0;
const double lambda = lumiPerLS*2.1/11246./23.31/296.;
for (int k=2; k<100; k++)
  {
    double add =  p1/(1.-pow(1.-p1,k)) * k * TMath::Poisson(k,lambda);
    pileupcorr += add;
    if(add < 0.0001)
      break;
  }
cout << "l=" <<  lumiPerLS << "   p1=" << p1 <<"    ----->  corr=" << pileupcorr *100. << "%" << endl;  

p1 = p1_single;
lumiPerLS = lumiPerLS_end;

double pileupcorr = 0;
const double lambda = lumiPerLS*2.1/11246./23.31/296.;
for (int k=2; k<100; k++)
  {
    double add =  p1/(1.-pow(1.-p1,k)) * k * TMath::Poisson(k,lambda);
    pileupcorr += add;
    if(add < 0.0001)
      break;
  }
cout << "l=" <<  lumiPerLS << "   p1=" << p1 <<"    ----->  corr=" << pileupcorr *100. << "%" << endl;

p1 = p1_double;
lumiPerLS = lumiPerLS_begin;

double pileupcorr = 0;
const double lambda = lumiPerLS*2.1/11246./23.31/296.;
for (int k=2; k<100; k++)
  {
    double add =  p1/(1.-pow(1.-p1,k)) * k * TMath::Poisson(k,lambda);
    pileupcorr += add;
    if(add < 0.0001)
      break;
  }
cout << "l=" <<  lumiPerLS << "   p1=" << p1 <<"    ----->  corr=" << pileupcorr *100. << "%" << endl;

p1 = p1_double;
lumiPerLS = lumiPerLS_end;

double pileupcorr = 0;
const double lambda = lumiPerLS*2.1/11246./23.31/296.;
for (int k=2; k<100; k++)
  {
    double add =  p1/(1.-pow(1.-p1,k)) * k * TMath::Poisson(k,lambda);
    pileupcorr += add;
    if(add < 0.0001)
      break;
  }
cout << "l=" <<  lumiPerLS << "   p1=" << p1 <<"    ----->  corr=" << pileupcorr *100. << "%" << endl;  
}
