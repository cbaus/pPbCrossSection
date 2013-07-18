{
const double lumiPerLS_begin = 2. * 1e6;
const double lumiPerLS_end = 0.9 * 1e6;
const double eff_acc_double = 0.9363;
const double eff_acc_single = 0.9643;

double lumiPerLS = lumiPerLS_begin;

double lambda = lumiPerLS*2.1/11246./23.31/296.;

 double f_pileup_single = 0;
 double f_pileup_double = 0;
 double nomi_single = 0;
 double denomi_single = 0;
 double nomi_double = 0;
 double denomi_double = 0;
 for (int k=1; k<100; k++)
   {
     double nomi_add_single =  eff_acc_single * k * TMath::Poisson(k,lambda);
     double denomi_add_single =  (1.-pow(1.-eff_acc_single,k)) * 1 * TMath::Poisson(k,lambda);
     double nomi_add_double =  eff_acc_double * k * TMath::Poisson(k,lambda);
     double denomi_add_double =  (1.-pow(1.-eff_acc_double,k)) * 1 * TMath::Poisson(k,lambda);
     nomi_single += nomi_add_single;
     denomi_single += denomi_add_single;
     nomi_double += nomi_add_double;
     denomi_double += denomi_add_double;
     if(nomi_add_single < 0.0001 && nomi_add_double < 0.0001 && denomi_add_single < 0.0001 && denomi_add_double < 0.0001)
       break;
   }
 if(!denomi_single || !denomi_double)
   {
     cerr << "Error in PU calculation" << endl; return;
   }
 f_pileup_single = nomi_single / denomi_single;
 f_pileup_double = nomi_double / denomi_double;
 cout << "l=" <<  lumiPerLS << "   p1=" << eff_acc_single <<"    ----->  corr=" << (f_pileup_single-1)*100. << "%" << "   corr=" << (f_pileup_double-1) *100. << "%" <<  endl;

lumiPerLS = lumiPerLS_end;
lambda = lumiPerLS*2.1/11246./23.31/296.;
 double f_pileup_single = 0;
 double f_pileup_double = 0;
 double nomi_single = 0;
 double denomi_single = 0;
 double nomi_double = 0;
 double denomi_double = 0;
 for (int k=1; k<100; k++)
   {
     double nomi_add_single =  eff_acc_single * k * TMath::Poisson(k,lambda);
     double denomi_add_single =  (1.-pow(1.-eff_acc_single,k)) * 1 * TMath::Poisson(k,lambda);
     double nomi_add_double =  eff_acc_double * k * TMath::Poisson(k,lambda);
     double denomi_add_double =  (1.-pow(1.-eff_acc_double,k)) * 1 * TMath::Poisson(k,lambda);
     nomi_single += nomi_add_single;
     denomi_single += denomi_add_single;
     nomi_double += nomi_add_double;
     denomi_double += denomi_add_double;
     if(nomi_add_single < 0.0001 && nomi_add_double < 0.0001 && denomi_add_single < 0.0001 && denomi_add_double < 0.0001)
       break;
   }
 if(!denomi_single || !denomi_double)
   {
     cerr << "Error in PU calculation" << endl; return;
   }
 f_pileup_single = nomi_single / denomi_single;
 f_pileup_double = nomi_double / denomi_double;
 cout << "l=" <<  lumiPerLS << "   p1=" << eff_acc_single <<"    ----->  corr=" << (f_pileup_single-1) *100. << "%" << "   corr=" << (f_pileup_double-1) *100. << "%" <<  endl;

}
