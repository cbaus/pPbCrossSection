#ifndef __CASTORCORRFACTORPP2015_H__
#define __CASTORCORRFACTORPP2015_H__

namespace ForwardRecord {
  static const unsigned int CSectors = 16;
  static const unsigned int CModules = 14;
  static const double absCasEscaleFactor = 0.020;
} //end namespace ForwardRecord

namespace castor {
  //E-map swap: s5m10<->s5m12, s6m10<->s6m12, s7m10<->s7m12
  //using for Katerina's value's (some correction values seem too big (for me big is all >3), some are zeros)
  const bool channelQuality[ForwardRecord::CSectors][ForwardRecord::CModules] =             // sector
    //mod 1   2     3     4     5      6    7     8     9    10     11    12   13    14
    {{true ,true ,true ,false,true ,true ,false,true ,true ,true ,true ,true ,true ,true }, // 1
     {true ,true ,true ,true ,true ,true ,false,true ,true ,false,true ,true ,true ,true }, // 2
     {true ,true ,true ,true ,true ,true ,false,false,false,true ,false,true ,true ,true }, // 3 //s3m9 ?
     {true ,true ,true ,true ,true ,true ,false,false,false,true ,false,true ,false,true }, // 4
     {true ,false,true ,true ,true ,true ,false,false,false,false,true ,false,true ,true }, // 5
     {true ,true ,true ,true ,true ,true ,false,false,false,false,true ,false,true ,true }, // 6 //s6m9 ?
     {true ,true ,true ,true ,false,true ,false,true ,true ,false,false,false,false,false}, // 7 //s7m10-14katerina?
     {true ,true ,true ,true ,true ,true ,false,false,true ,false,false,false,false,false}, // 8 //s8m10-14katerina?
     {true ,true ,true ,true ,true ,true ,false,true ,true ,true ,true ,true ,false,true }, // 9 //s9m13?
     {true ,true ,true ,true ,true ,true ,false,true ,true ,true ,true ,true ,true ,true }, // 10
     {true ,true ,true ,true ,true ,true ,false,false,true ,true ,false,true ,true ,true }, // 11 //s11m11katerina?(was also there before)
     {true ,true ,true ,true ,true ,true ,false,false,true ,true ,false,true ,true ,true }, // 12 //s12m11katerina?(was also there before)
     {true ,true ,true ,true ,true ,false,false,false,false,true ,false,true ,false,true }, // 13 //s13m11katerina?(was also there before) //m9s13-SNP//m13s13-SNP
     {true ,true ,true ,true ,true ,true ,false,false,true ,true ,false,true ,true ,false}, // 14 //s14m11katerina?(was also there before)
     {true ,true ,true ,true ,true ,true ,false,false,true ,false,true ,true ,true ,true }, // 15
     {true ,true ,true ,true ,true ,false,false,false,true ,true ,true ,true ,true ,true }};// 16 //m8s16-strange-noise-peak(SNP),closer look needed

  // Katerina's values using halo muon data (w/o TOTEM), already scaled by s9m4 (what about factor 2 EM vs HAD ?)
  const double channelGainQE[ForwardRecord::CSectors][ForwardRecord::CModules] =                                                                                // sector
    //mod 1          2          3          4           5          6          7         8          9          10         11          12        13           14
    {{  0.7510,    0.8700,    2.7370,    0.0000,    0.3630,    0.6430,    0.0000,    0.3100,    0.2120,    0.2740,    0.3030,    0.1690,    0.2650,    0.1550}, //1
     {  0.6190,    0.6160,    1.8130,    0.8690,    0.1820,    0.6280,    0.0000,    0.5070,    0.1680,    0.2910,    0.3380,    0.1460,    0.2490,    0.1250}, //2
     {  1.0700,    0.6510,    1.4250,    0.7660,    0.3040,    0.1930,    8.2170,   13.2900,    0.4650,    0.2350,    0.0000,    0.2950,    0.3430,    0.3510}, //3
     {  0.5310,    0.3300,    0.8910,    0.8260,    0.1170,    0.3300,    0.0000,    0.0000,    0.0000,    0.6390,    0.0000,    0.3170,    0.0000,    0.4580}, //4
     {  0.6120,    0.0000,    1.3410,    0.7020,    0.1560,    0.5690,    0.8360,    0.0000,    0.0000,    0.5230,    0.2360,    0.3290,    0.3990,    0.3420}, //5
     {  1.3130,    0.4870,    1.4000,    0.6320,    0.1990,    0.7950,    1.2090,    0.0000,    0.5100,    0.7060,    0.2330,    0.2800,    0.4830,    0.4410}, //6
     {  0.4160,    0.2820,    1.0430,    0.3130,    0.1140,    0.0860,  250.6690,    0.1950,    0.4200,    6.9160,    3.4790,    1.5110,    4.8590,    3.5340}, //7
     {  0.3420,    0.2950,    1.1980,    1.4030,    0.2130,    1.0730,    0.0000,    0.2060,    0.6350,   27.2690,    9.4210,    3.3400,    3.4880,    1.0100}, //8
     {  0.3030,    0.8460,    1.4120,    1.0000,    0.2180,    0.8830,    0.0000,    0.1320,    0.1950,    0.2490,    0.2250,    0.2270,    0.2990,    0.2780}, //9
     {  0.9040,    1.4030,    2.6580,    1.1900,    0.2350,    1.5570,    0.0000,    0.3160,    0.1990,    0.3100,    0.1790,    0.2510,    0.2510,    0.2520}, //10
     {  1.0160,    0.9930,    1.6950,    0.8870,    0.2850,    0.6230,    0.0000,   10.0790,    0.3730,    0.2440,    9.6350,    0.5240,    0.6990,    0.3790}, //11
     {  1.1690,    1.1300,    2.1400,    1.3920,    0.2630,    1.2470,    0.0000,    0.0000,    0.5670,    0.3030,   99.3510,    0.3510,    0.1980,    0.3560}, //12
     {  0.9160,    1.2700,    1.6430,    0.8070,    0.2310,    2.3020,    0.0000,    0.0000,    0.3230,    0.2910,    0.0000,    0.3430,    0.1280,    0.3080}, //13
     {  0.6010,    0.9840,    2.1400,    0.8210,    0.1770,    1.0970,    0.0000,    0.0000,    0.2030,    0.2920,   16.6350,    0.3020,    0.3510,    0.3680}, //14
     {  0.7590,    1.3650,    2.9620,    1.1740,    0.3800,    2.3370,    0.0000,  517.2540,    0.2690,    0.0000,    0.1940,    0.2740,    0.2800,    0.4100}, //15
     {  0.7420,    0.9720,    2.4600,    0.9240,    0.2200,    0.1630,    3.9070,    0.1970,    0.2700,    0.2580,    0.1510,    0.1340,    0.2790,    0.2620}};//16

  //LED data, r238434(B=3.8T) to r238342(B=0T), HV ~ 1800V, average amplitude w/o ZS
  const double corr0Tto38T[ForwardRecord::CSectors][ForwardRecord::CModules] =
     //mod 1          2              3          4           5               6          7               8          9          10              11          12        13           14
{{     1.48512,     1.39079,     1.23936,     1.20485,     1.53704,     0.53953,     0.00404,     0.73703,     0.99343,     0.97474,     1.07208,     1.00473,     1.21792,     1.38684}, // s 1
 {     1.34932,     1.31645,     1.16597,     1.07258,     1.10938,     0.91177,     0.00558,     0.06076,     0.94764,     0.92360,     0.96920,     1.16183,     1.18510,     1.23530}, // s 2
 {     1.19602,     1.46404,     1.26956,     1.27067,     1.55640,     0.28252,     0.00494,    -0.00000,     0.04254,     0.97351,     1.00799,     1.01488,     1.17762,     1.33190}, // s 3
 {     1.28421,     1.54017,     1.14003,     1.21478,     1.62988,     0.52778,     0.00217,    -0.00000,    -0.00000,     0.85758,     1.01701,     1.11001,     1.15505,     1.45989}, // s 4
 {     1.08401,     1.27912,     1.15338,     1.02110,     1.48480,     0.18463,     0.04877,    -0.00000,     0.00035,     0.93851,     0.91677,     1.22820,     1.04381,     1.25080}, // s 5
 {     1.09123,     1.19290,     1.02762,     1.15109,     1.01973,     0.27162,     0.01023,     0.00001,     0.22811,     0.95087,     1.02469,     0.94455,     1.14685,     1.19158}, // s 6
 {     1.35555,     1.28333,     0.80874,     1.08100,     0.98072,     0.39716,     0.00095,     0.29955,     1.01768,     0.94803,     0.88673,     1.02905,     1.02306,     1.40891}, // s 7
 {     1.40568,     1.28872,     1.01763,     1.12056,     1.12694,     0.05847,     0.00022,     0.90227,     1.01943,     1.14345,     1.04596,     0.96388,     1.15198,     1.29044}, // s 8
 {     1.06251,     1.14400,     1.01962,     1.03241,     1.02967,     0.18843,     0.00002,     0.80946,     0.87329,     0.88913,     0.86814,     1.13983,     1.05866,     1.04589}, // s 9
 {     1.01657,     1.16994,     1.03485,     1.06841,     1.12475,     0.21338,     0.00008,     0.37391,     0.92616,     0.95543,     0.95621,     1.01440,     0.86022,     1.37397}, // s10
 {     0.98176,     0.98911,     0.99062,     1.03866,     1.21427,     0.10842,     0.00007,     0.05552,     0.71987,     0.87018,     0.89760,     0.95953,     1.19473,     1.22346}, // s11
 {     0.98305,     0.99199,     0.99141,     1.08809,     1.13973,     0.14923,     0.00018,     0.00005,     0.36019,     0.86800,     0.93984,     1.06136,     1.08128,     1.07615}, // s12
 {     0.99641,     0.99703,     0.71680,     1.00467,     1.24617,     0.03526,     0.00227,    -0.00000,     0.36686,     0.68199,     0.98805,     1.09305,     1.19635,     1.33507}, // s13
 {     1.00974,     1.00226,     0.98903,     1.03647,     1.27846,     0.39654,     0.00116,     0.12479,     0.66424,     1.14252,     1.01872,     1.12014,     1.26828,     0.77948}, // s14
 {     1.17599,     1.08516,     1.01189,    -0.00000,     1.36638,     0.07736,    -0.00000,     0.00148,     0.95000,     0.86826,     0.97357,     1.04934,     1.09820,     1.07579}, // s15
 {     1.05973,     1.13987,     1.05313,     0.83706,     1.31051,     0.33701,     0.00358,     0.74754,     0.97755,     0.87077,     1.01832,     1.10647,     1.15042,     1.08966}};// s16

} //end namespace castor

#endif
