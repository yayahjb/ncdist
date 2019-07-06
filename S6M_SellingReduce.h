#ifndef CS6M_SELLING_REDUCE
  #define CS6M_SELLING_REDUCE

  /* #define CS6M_DEBUG zz*/ 
  #ifdef CS6M_DEBUG_PRINT_COUT
    #undef CS6M_DEBUG_PRINT_COUT
  #endif
  #ifdef CS6M_DEBUG_PRINT_FPRINTF
    #undef CS6M_DEBUG_PRINT_FPRINTF
  #endif
  #ifdef CS6M_DEBUG
    #ifdef __cplusplus
      #include <iostream>
      #define CS6M_DEBUG_PRINT_COUT yes
    #else
      #include <stdio.h>
      #define CS6M_DEBUG_PRINT_FPRINTF yes
    #endif
  #endif

  #ifdef __cplusplus
    #define CS6M_sqrt(x) std::sqrt(x)
  #else
    #define CS6M_sqrt(x) dsqrt(x)
  #endif


 /* #define CS6M_DEBUG_PRINT_COUT yes */
  #define CS6M_S6BC 0
  #define CS6M_S6AC 1
  #define CS6M_S6AB 2
  #define CS6M_S6AD 3
  #define CS6M_S6BD 4
  #define CS6M_S6CD 5

  #define CS6M_D7A2 0
  #define CS6M_D7B2 1
  #define CS6M_D7C2 2
  #define CS6M_D7D2 3
  #define CS6M_D7_A_D_2 4
  #define CS6M_D7_B_D_2 5
  #define CS6M_D7_C_D_2 6

  #define CS6M_G6A2 0
  #define CS6M_G6B2 1
  #define CS6M_G6C2 2
  #define CS6M_G62BC 3
  #define CS6M_G62AC 4
  #define CS6M_G62AB 5

  #define CS6M_CELLA 0
  #define CS6M_CELLB 1
  #define CS6M_CELLC 2
  #define CS6M_CELLALPHA 3
  #define CS6M_CELLBETA 4
  #define CS6M_CELLGAMMA 5

  #define CS6M_V7A 0
  #define CS6M_V7B 1
  #define CS6M_V7C 2
  #define CS6M_V7ASTARINV 3
  #define CS6M_V7BSTARINV 4
  #define CS6M_V7CSTARINV 5
  #define CS6M_V7VOLCROOT 6
  


  #define CS6M_abs(x) ((x)<0?-(x):(x))

  #define CS6M_CountPositive(s6vec,delta)  \
    ( ((s6vec[CS6M_S6BC]>delta)?1:0)  \
    + ((s6vec[CS6M_S6AC]>delta)?1:0)  \
    + ((s6vec[CS6M_S6AB]>delta)?1:0)  \
    + ((s6vec[CS6M_S6AD]>delta)?1:0)  \
    + ((s6vec[CS6M_S6BD]>delta)?1:0)  \
    + ((s6vec[CS6M_S6CD]>delta)?1:0)  )

  #define CS6M_NegativeSumOfScalars(s6vec) \
     (-s6vec[CS6M_S6BC]-s6vec[CS6M_S6AC]-s6vec[CS6M_S6AB]-s6vec[CS6M_S6AD]-s6vec[CS6M_S6BD]-s6vec[CS6M_S6CD]) 


  #define CS6M_comptovec6(a,b,c,d,e,f,vec6) { \
    vec6[0]=a;vec6[1]=b;vec6[2]=c;vec6[3]=d;vec6[4]=e;vec6[5]=f; }

  #define CS6M_CellradtoG6(cellrad,g6vec) { \
    g6vec[CS6M_G6A2] = cellrad[CS6M_CELLA]*cellrad[CS6M_CELLA]; \
    g6vec[CS6M_G6B2] = cellrad[CS6M_CELLB]*cellrad[CS6M_CELLB]; \
    g6vec[CS6M_G6C2] = cellrad[CS6M_CELLC]*cellrad[CS6M_CELLC]; \
    g6vec[CS6M_G62BC] = cellrad[CS6M_CELLB]*cellrad[CS6M_CELLC]*2.*cos(cellrad[CS6M_CELLALPHA]); \
    g6vec[CS6M_G62AC] = cellrad[CS6M_CELLA]*cellrad[CS6M_CELLC]*2.*cos(cellrad[CS6M_CELLBETA]); \
    g6vec[CS6M_G62AB] = cellrad[CS6M_CELLA]*cellrad[CS6M_CELLB]*2.*cos(cellrad[CS6M_CELLGAMMA]); }

  #define CS6M_D7toG6(d7vec,g6vec) { \
    g6vec[CS6M_G6A2] = d7vec[CS6M_D7A2]; \
    g6vec[CS6M_G6B2] = d7vec[CS6M_D7B2]; \
    g6vec[CS6M_G6C2] = d7vec[CS6M_D7C2]; \
    g6vec[CS6M_G62BC] = d7vec[CS6M_D7_A_D_2]-d7vec[CS6M_D7C2]-d7vec[CS6M_D7B2];\
    g6vec[CS6M_G62AC] = d7vec[CS6M_D7_B_D_2]-d7vec[CS6M_D7C2]-d7vec[CS6M_D7A2];\
    g6vec[CS6M_G62AB] = d7vec[CS6M_D7_C_D_2]-d7vec[CS6M_D7A2]-d7vec[CS6M_D7B2];\
  }


  #define CS6M_D7toS6(d7vec,s6vec) { \
    s6vec[CS6M_S6BC] = -(-d7vec[CS6M_D7_A_D_2]+d7vec[CS6M_D7C2]+d7vec[CS6M_D7B2])/2;\
    s6vec[CS6M_S6AC] = -(-d7vec[CS6M_D7_B_D_2]+d7vec[CS6M_D7C2]+d7vec[CS6M_D7A2])/2;\
    s6vec[CS6M_S6AB] = -(-d7vec[CS6M_D7_C_D_2]+d7vec[CS6M_D7B2]+d7vec[CS6M_D7A2])/2;\
    s6vec[CS6M_S6AD] = -(-d7vec[CS6M_D7_A_D_2]+d7vec[CS6M_D7D2]+d7vec[CS6M_D7A2])/2;\
    s6vec[CS6M_S6BD] = -(-d7vec[CS6M_D7_B_D_2]+d7vec[CS6M_D7D2]+d7vec[CS6M_D7B2])/2;\
    s6vec[CS6M_S6CD] = (-d7vec[CS6M_D7_B_D_2]-d7vec[CS6M_D7_A_D_2]+d7vec[CS6M_D7B2]+d7vec[CS6M_D7A2])/2;\
  }

  #define CS6M_G6toCell(g6vec,cell) {                                                     \
    double cosalpha, cosbeta, cosgamma;                                                   \
    double sinalpha, sinbeta, singamma;                                                   \
    cell[CS6M_CELLA]=CS6M_sqrt(fabs(g6vec[CS6M_G6A2]));                                        \
    cell[CS6M_CELLB]=CS6M_sqrt(fabs(g6vec[CS6M_G6B2]));                                        \
    cell[CS6M_CELLC]=CS6M_sqrt(fabs(g6vec[CS6M_G6C2]));                                        \
    cosalpha=g6vec[CS6M_G62BC]/(cell[CS6M_CELLB]*cell[CS6M_CELLC]);                       \
    cosbeta =g6vec[CS6M_G62BC]/(cell[CS6M_CELLB]*cell[CS6M_CELLC]);                       \
    cosgamma=g6vec[CS6M_G62BC]/(cell[CS6M_CELLB]*cell[CS6M_CELLC]);                       \
    sinalpha=CS6M_sqrt(fabs((1.-cosalpha*cosalpha)));                                           \
    sinbeta =CS6M_sqrt(fabs((1.-cosbeta* cosbeta)));                                            \
    singamma=CS6M_sqrt(fabs((1.-cosgamma*cosgamma)));                                           \
    cell[CS6M_CELLALPHA]=atan2(sinalpha,cosalpha);                                        \
    cell[CS6M_CELLBETA] =atan2(sinbeta, cosbeta);                                         \
    cell[CS6M_CELLGAMMA]=atan2(singamma,cosgamma);                                        \
  }

  #define CS6M_G6toCelldeg(g6vec,cell) {                                             \
    double pi=3.1415926535897932384626433832795;                                          \
    double todeg=180./pi;                                                                 \
    double cosalpha, cosbeta, cosgamma;                                                   \
    double sinalpha, sinbeta, singamma;                                                   \
    cell[CS6M_CELLA]=CS6M_sqrt(fabs(g6vec[CS6M_G6A2]));                                        \
    cell[CS6M_CELLB]=CS6M_sqrt(fabs(g6vec[CS6M_G6B2]));                                        \
    cell[CS6M_CELLC]=CS6M_sqrt(fabs(g6vec[CS6M_G6C2]));                                        \
    cosalpha=g6vec[CS6M_G62BC]/(cell[CS6M_CELLB]*cell[CS6M_CELLC]);                       \
    cosbeta =g6vec[CS6M_G62BC]/(cell[CS6M_CELLB]*cell[CS6M_CELLC]);                       \
    cosgamma=g6vec[CS6M_G62BC]/(cell[CS6M_CELLB]*cell[CS6M_CELLC]);                       \
    sinalpha=CS6M_sqrt(fabs((1.-cosalpha*cosalpha));                                           \
    sinbeta =CS6M_sqrt(fabs((1.-cosbeta* cosbeta));                                            \
    singamma=CS6M_sqrt(fabs((1.-cosgamma*cosgamma));                                           \
    cell[CS6M_CELLALPHA]=todeg*atan2(sinalpha,cosalpha);                                  \
    cell[CS6M_CELLBETA] =todeg*atan2(sinbeta, cosbeta);                                   \
    cell[CS6M_CELLGAMMA]=todeg*atan2(singamma,cosgamma);                                  \
  }

  #define CS6M_CelltoG6(cell,g6vec) {                                                     \
    double delta;                                                                         \
    size_t i;                                                                             \
    delta = fabs(cell[CS6M_CELLA]);                                                       \
    if (cell[CS6M_CELLB] > delta) delta = cell[CS6M_CELLB];                               \
    if (cell[CS6M_CELLC] > delta) delta = cell[CS6M_CELLC];                               \
    delta *= 1.e-12;                                                                      \
    g6vec[CS6M_G6A2] = cell[CS6M_CELLA]*cell[CS6M_CELLA];                                 \
    g6vec[CS6M_G6B2] = cell[CS6M_CELLB]*cell[CS6M_CELLB];                                 \
    g6vec[CS6M_G6C2] = cell[CS6M_CELLC]*cell[CS6M_CELLC];                                 \
    g6vec[CS6M_G62BC] = 2.0*cell[CS6M_CELLB]*cell[CS6M_CELLC]*cos(cell[CS6M_CELLALPHA]);  \
    g6vec[CS6M_G62AC] = 2.0*cell[CS6M_CELLA]*cell[CS6M_CELLC]*cos(cell[CS6M_CELLBETA]);   \
    g6vec[CS6M_G62AB] = 2.0*cell[CS6M_CELLA]*cell[CS6M_CELLB]*cos(cell[CS6M_CELLGAMMA]);  \
    for ( i=3; i<6; ++i ) if ( std::fabs(g6vec[i]) < delta ) g6vec[i] = 0.0;              \
  }

  #define CS6M_CelldegtoG6(cell,g6vec) {                                            \
    double delta;                                                                         \
    double pi=3.1415926535897932384626433832795;                                          \
    double torad=pi/180.;                                                                 \
    size_t i;                                                                             \
    delta = fabs(cell[CS6M_CELLA]);                                                       \
    if (cell[CS6M_CELLB] > delta) delta = cell[CS6M_CELLB];                               \
    if (cell[CS6M_CELLC] > delta) delta = cell[CS6M_CELLC];                               \
    delta *= 1.e-12;                                                                      \
    g6vec[CS6M_G6A2] = cell[CS6M_CELLA]*cell[CS6M_CELLA];                                 \
    g6vec[CS6M_G6B2] = cell[CS6M_CELLB]*cell[CS6M_CELLB];                                 \
    g6vec[CS6M_G6C2] = cell[CS6M_CELLC]*cell[CS6M_CELLC];                                 \
    g6vec[CS6M_G62BC] = 2.0*cell[CS6M_CELLB]*cell[CS6M_CELLC]*cos(cell[CS6M_CELLALPHA]*torad);  \
    g6vec[CS6M_G62AC] = 2.0*cell[CS6M_CELLA]*cell[CS6M_CELLC]*cos(cell[CS6M_CELLBETA]*torad);   \
    g6vec[CS6M_G62AB] = 2.0*cell[CS6M_CELLA]*cell[CS6M_CELLB]*cos(cell[CS6M_CELLGAMMA]*torad);  \
    for ( i=3; i<6; ++i ) if ( std::fabs(g6vec[i]) < delta ) g6vec[i] = 0.0;              \
  }



  #define CS6M_G6toD7(g6vec,d7vec) { \
    d7vec[CS6M_D7A2] = g6vec[CS6M_G6A2]; \
    d7vec[CS6M_D7B2] = g6vec[CS6M_G6B2]; \
    d7vec[CS6M_D7C2] = g6vec[CS6M_G6C2]; \
    d7vec[CS6M_D7D2] = g6vec[CS6M_G6A2]+g6vec[CS6M_G6B2]+g6vec[CS6M_G6C2]\
               +g6vec[CS6M_G62BC]+g6vec[CS6M_G62AC]+g6vec[CS6M_G62AB]; \
    d7vec[CS6M_D7_A_D_2] = g6vec[CS6M_G6B2]+g6vec[CS6M_G6C2]+g6vec[CS6M_G62BC]; \
    d7vec[CS6M_D7_B_D_2] = g6vec[CS6M_G6A2]+g6vec[CS6M_G6C2]+g6vec[CS6M_G62AC]; \
    d7vec[CS6M_D7_C_D_2] = g6vec[CS6M_G6A2]+g6vec[CS6M_G6B2]+g6vec[CS6M_G62AB]; \
  }

  #define CS6M_G6toS6(g6vec,s6vec) { \
    s6vec[CS6M_S6BC] = g6vec[CS6M_G62BC]/2.; \
    s6vec[CS6M_S6AC] = g6vec[CS6M_G62AC]/2.; \
    s6vec[CS6M_S6AB] = g6vec[CS6M_G62AB]/2.; \
    s6vec[CS6M_S6AD] = -g6vec[CS6M_G6A2]-s6vec[CS6M_S6AB]-s6vec[CS6M_S6AC]; \
    s6vec[CS6M_S6BD] = -s6vec[CS6M_S6AB]-g6vec[CS6M_G6B2]-s6vec[CS6M_S6BC]; \
    s6vec[CS6M_S6CD] = -s6vec[CS6M_S6AC]-s6vec[CS6M_S6BC]-g6vec[CS6M_G6C2]; \
  }

  #define CS6M_S6toD7(s6vec,d7vec) { \
    d7vec[CS6M_D7A2] = -s6vec[CS6M_S6AD] -s6vec[CS6M_S6AB] -s6vec[CS6M_S6AC];\
    d7vec[CS6M_D7B2] = -s6vec[CS6M_S6BD] -s6vec[CS6M_S6AB] -s6vec[CS6M_S6BC];\
    d7vec[CS6M_D7C2] = -s6vec[CS6M_S6CD] -s6vec[CS6M_S6AC] -s6vec[CS6M_S6BC];\
    d7vec[CS6M_D7D2] = -s6vec[CS6M_S6AD] -s6vec[CS6M_S6BD] -s6vec[CS6M_S6CD];\
    d7vec[CS6M_D7_A_D_2] = -s6vec[CS6M_S6BD] -s6vec[CS6M_S6AB] -s6vec[CS6M_S6CD] -s6vec[CS6M_S6AC];\
    d7vec[CS6M_D7_B_D_2] = -s6vec[CS6M_S6AD] -s6vec[CS6M_S6AB] -s6vec[CS6M_S6CD] -s6vec[CS6M_S6BC];\
    d7vec[CS6M_D7_C_D_2] = -s6vec[CS6M_S6AD] -s6vec[CS6M_S6AC] -s6vec[CS6M_S6BD] -s6vec[CS6M_S6BC];\
  }

  #define CS6M_S6toG6(s6vec,g6vec) { \
    g6vec[CS6M_G6A2] = -s6vec[CS6M_S6AD] -s6vec[CS6M_S6AB] -s6vec[CS6M_S6AC];\
    g6vec[CS6M_G6B2] = -s6vec[CS6M_S6BD] -s6vec[CS6M_S6AB] -s6vec[CS6M_S6BC];\
    g6vec[CS6M_G6C2] = -s6vec[CS6M_S6CD] -s6vec[CS6M_S6AC] -s6vec[CS6M_S6BC];\
    g6vec[CS6M_G62BC] = 2.*s6vec[CS6M_S6BC];\
    g6vec[CS6M_G62AC] = 2.*s6vec[CS6M_S6AC];\
    g6vec[CS6M_G62AB] = 2.*s6vec[CS6M_S6AB];\
  }

  #ifdef CS6M_DEBUG
     #ifndef CS6M_VOLCHECK_PREP
      #define CS6M_VOLCHECK_PREP(xg6unred)  \
        LRL_Cell cellred; \
        LRL_Cell newcellred; \
        G6 g6cellred; \
        G6 g6newcellred; \
        double entryvol,curvol;\
        CS6M_G6toCell(xg6unred,cellred);\
        CS6M_Cellvolume(cellred,entryvol); 
      #define CS6M_VOLCHECK(label,xg6red)       \
        CS6M_G6toCell(xg6red,newcellred);   \
        CS6M_Cellvolume(newcellred,curvol);\
        CS6M_CellradtoG6(cellred,g6cellred);  \
        CS6M_CellradtoG6(newcellred,g6newcellred); \
        if (fabs(entryvol-curvol)> 0. ) { \
          std::cerr<<label<<": Volume change: "<< entryvol << " to "<<curvol<<std::endl; \
          std::cerr<<"cellred: " << cellred[0] <<", "    \
          << cellred[1] <<", "                           \
          << cellred[2] <<", "                           \
          << cellred[3] <<", "                           \
          << cellred[4] <<", "                           \
          << cellred[5] <<std::endl;                     \
          std::cerr<<"g6cellred: " << g6cellred[0] <<", "    \
          << g6cellred[1] <<", "                           \
          << g6cellred[2] <<", "                           \
          << g6cellred[3] <<", "                           \
          << g6cellred[4] <<", "                           \
          << g6cellred[5] <<std::endl;                     \
          std::cerr<<"newcellred: " << newcellred[0] <<", "    \
            << newcellred[1] <<", "                           \
            << newcellred[2] <<", "                           \
            << newcellred[3] <<", "                           \
            << newcellred[4] <<", "                           \
            << newcellred[5] <<std::endl;                     \
          std::cerr<<"g6newcellred: " << g6newcellred[0] <<", "    \
            << g6newcellred[1] <<", "                           \
            << g6newcellred[2] <<", "                           \
            << g6newcellred[3] <<", "                           \
            << g6newcellred[4] <<", "                           \
            << g6newcellred[5] <<std::endl;                     \
           entryvol=curvol;                                     \
         } 
    #endif
  #else
    #ifndef CS6M_VOLCHECK_PREP
      #define  CS6M_VOLCHECK_PREP(dummy)
      #define  CS6M_VOLCHECK(label,dummy)
    #endif
  #endif

  #define CS6M_G6Reduce(g6unred,g6red,reduced) { \
    double s6in[6];                           \
    double s6out[6];                          \
    int notdone;                              \
    int ii;                                   \
    notdone = 1;                              \
    int redpass;                              \
    double temp;                              \
    double delta = 0;                         \
    CS6M_VOLCHECK_PREP(g6unred);              \
    CS6M_G6toS6(g6unred,s6in);                \
    reduced=1;                                \
    CS6M_S6Reduce(s6in,s6out,reduced);        \
    if (reduced) {                            \
       CS6M_S6toG6(s6out,g6red);              \
    } else {                                  \
      for (ii=0; ii < 6; ii++) g6red[ii]=g6unred[ii];\
    }                                         \
    delta=g6unred[2];                         \
    CS6M_VOLCHECK("S6",g6red);                \
    if (g6unred[1] > delta) delta=g6unred[1]; \
    if (g6unred[0] > delta) delta=g6unred[0]; \
    delta = delta * 1.e-12;                   \
    if (delta < 1.e-12) delta=1.e-12;         \
    if (fabs(g6red[CS6M_G62BC])  < delta ) g6red[CS6M_G62BC] = 0;      \
    if (fabs(g6red[CS6M_G62AC])  < delta ) g6red[CS6M_G62AC] = 0;      \
    if (fabs(g6red[CS6M_G62AB])  < delta ) g6red[CS6M_G62AB] = 0;      \
                                                                       \
    redpass=0;                                                         \
    while ( notdone && redpass < 1000) {                               \
    if (g6red[CS6M_G6A2] < 0.                                          \
        || g6red[CS6M_G6B2] < 0.                                       \
        || g6red[CS6M_G6C2] < 0.                                       \
        || g6red[CS6M_G6B2] + g6red[CS6M_G6C2] + g6red[CS6M_G62BC] < 0.\
        || g6red[CS6M_G6A2] + g6red[CS6M_G6C2] + g6red[CS6M_G62AC] < 0.\
        || g6red[CS6M_G6A2] + g6red[CS6M_G6B2] + g6red[CS6M_G62AB] < 0.\
        || g6red[CS6M_G6B2] + g6red[CS6M_G6C2] - g6red[CS6M_G62BC] < 0.\
        || g6red[CS6M_G6A2] + g6red[CS6M_G6C2] - g6red[CS6M_G62AC] < 0.\
        || g6red[CS6M_G6A2] + g6red[CS6M_G6B2] - g6red[CS6M_G62AB] < 0.\
    ) {                                                                \
      notdone=0; break;                                                \
    }                                                                  \
      /* g1 < g2 or g1=g2 && |g4| < |g5| */                            \
      if ( g6red[CS6M_G6A2] > g6red[CS6M_G6B2]+delta                   \
           || ( CS6M_abs(g6red[CS6M_G6A2] - g6red[CS6M_G6B2]) < delta      \
           && CS6M_abs(g6red[CS6M_G62BC]) > CS6M_abs(g6red[CS6M_G62AC])+delta) ) { \
        temp = g6red[CS6M_G6A2];                                 \
        g6red[CS6M_G6A2] = g6red[CS6M_G6B2];                     \
        g6red[CS6M_G6B2] = temp;                                 \
        temp = g6red[CS6M_G62BC];                                \
        g6red[CS6M_G62BC] = g6red[CS6M_G62AC];                   \
        g6red[CS6M_G62AC] = temp;                                \
        /* std::cout<<"bd1 "<<g6red[0]<<" "<<g6red[1]<<" "       \
        <<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]  \
        <<" "<<std::endl; */                                     \
        CS6M_VOLCHECK("swap a b ",g6red);                        \
        redpass++; continue;                                     \
      }                                                          \
      /* g2 < g3 or g2=g3 && |g5| < |g6| */                      \
      if ( g6red[CS6M_G6B2] > g6red[CS6M_G6C2]+delta             \
            ||  (fabs(g6red[CS6M_G6B2] - g6red[CS6M_G6C2]) < delta \
            && CS6M_abs(g6red[CS6M_G62AC]) > CS6M_abs(g6red[CS6M_G62AB])+delta) ) {\
        temp = g6red[CS6M_G6B2];                                 \
        g6red[CS6M_G6B2] = g6red[CS6M_G6C2];                     \
        g6red[CS6M_G6C2] = temp;                                 \
        temp = g6red[CS6M_G62AC];                                \
        g6red[CS6M_G62AC] = g6red[CS6M_G62AB];                   \
        g6red[CS6M_G62AB] = temp;                                \
        CS6M_VOLCHECK("swap b c ",g6red);                        \
        /* std::cout<<"bd2 "<<g6red[0]<<" "<<g6red[1]<<" "       \
        <<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]  \
        <<" "<<std::endl; */                                     \
        redpass++; continue;                                     \
      }                                                          \
      if ( g6red[CS6M_G62BC]*g6red[CS6M_G62AC]*g6red[CS6M_G62AB] > 0. ) { \
        if (g6red[CS6M_G62BC] < 0. || g6red[CS6M_G62AC] < 0.     \
          || g6red[CS6M_G62AC] < 0. ) {                          \
          if ( g6red[CS6M_G62BC] < 0. ) g6red[CS6M_G62BC] = -g6red[CS6M_G62BC]; \
          if ( g6red[CS6M_G62AC] < 0. ) g6red[CS6M_G62AC] = -g6red[CS6M_G62AC]; \
          if ( g6red[CS6M_G62AB] < 0. ) g6red[CS6M_G62AB] = -g6red[CS6M_G62AB]; \
          /* std::cout<<"+++ "<<g6red[0]<<" "<<g6red[1]<<" "     \
          <<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]\
          <<" "<<std::endl; */                                   \
          CS6M_VOLCHECK("--- ",g6red);                           \
          redpass++; continue;                                   \
        }                                                        \
      } else {                                                   \
        if (g6red[CS6M_G62BC] > 0. || g6red[CS6M_G62AC] > 0.     \
          || g6red[CS6M_G62AC] > 0. ) {                          \
          if ( g6red[CS6M_G62BC] > 0. ) g6red[CS6M_G62BC] = -g6red[CS6M_G62BC]; \
          if ( g6red[CS6M_G62AC] > 0. ) g6red[CS6M_G62AC] = -g6red[CS6M_G62AC]; \
          if ( g6red[CS6M_G62AB] > 0. ) g6red[CS6M_G62AB] = -g6red[CS6M_G62AB]; \
          /* std::cout<<"--- "<<g6red[0]<<" "<<g6red[1]<<" "     \
          <<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]\
          <<" "<<std::endl; */                                   \
          CS6M_VOLCHECK("+++ ",g6red);                           \
          redpass++; continue;                                   \
        }                                                        \
      }                                                          \
                                                                 \
      /* 678 boundaries outside */                               \
      if (g6red[CS6M_G62BC] >  g6red[CS6M_G6B2]+delta) {         \
          g6red[CS6M_G6C2] = g6red[CS6M_G6B2]+g6red[CS6M_G6C2]-g6red[CS6M_G62BC];\
          g6red[CS6M_G62BC] = -2.*g6red[CS6M_G6B2]+g6red[CS6M_G62BC];            \
          g6red[CS6M_G62AC] = -g6red[CS6M_G62AC]+g6red[CS6M_G62AB];               \
          /* std::cout<<"bd6+ "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "     \
          <<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/            \
          CS6M_VOLCHECK("6+ ",g6red);                                            \
          redpass++; continue;                                                   \
      }                                                                          \
      if (g6red[CS6M_G62BC] < -g6red[CS6M_G6B2]-delta) {                         \
          g6red[CS6M_G6C2] = g6red[CS6M_G6B2]+g6red[CS6M_G6C2]+g6red[CS6M_G62BC];\
          g6red[CS6M_G62BC] = 2.*g6red[CS6M_G6B2]+g6red[CS6M_G62BC];             \
          g6red[CS6M_G62AC] = g6red[CS6M_G62AC]+g6red[CS6M_G62AB];               \
          /* std::cout<<"bd8+"<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "      \
          <<g6red[3 ]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/           \
          CS6M_VOLCHECK("8+ ",g6red);                                            \
          redpass++; continue;                                                   \
      }                                                                          \
      /* 9AB boundaries outside */                                               \
      if (g6red[CS6M_G62AC] > g6red[CS6M_G6A2]+delta) {                          \
          g6red[CS6M_G6C2] = g6red[CS6M_G6A2]+g6red[CS6M_G6C2]-g6red[CS6M_G62AC];\
          g6red[CS6M_G62AC] = -2.*g6red[CS6M_G6A2]+g6red[CS6M_G62AC];            \
          g6red[CS6M_G62BC] = g6red[CS6M_G62BC]-g6red[CS6M_G62AB];               \
          /* std::cout<<"bd9+ "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "     \
          <<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/            \
          CS6M_VOLCHECK("9+ ",g6red);                                            \
          redpass++; continue;                                                   \
      }                                                                          \
      if (g6red[CS6M_G62AC] < -g6red[CS6M_G6A2]-delta) {                         \
          g6red[CS6M_G6C2] = g6red[CS6M_G6A2]+g6red[CS6M_G6C2]+g6red[CS6M_G62AC];\
          g6red[CS6M_G62AC] = 2.*g6red[CS6M_G6A2]+g6red[CS6M_G62AC];             \
          g6red[CS6M_G62BC] = g6red[CS6M_G62BC]+g6red[CS6M_G62AB];               \
          /* std::cout<<"bdB+ "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "     \
          <<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/            \
          CS6M_VOLCHECK("B+ ",g6red);                                            \
          redpass++; continue;                                                   \
      }                                                                          \
      /*  CDE boundaries outside */                                              \
      if (g6red[CS6M_G62AB] > g6red[CS6M_G6A2]+delta) {                          \
          g6red[CS6M_G6B2] = g6red[CS6M_G6A2]+g6red[CS6M_G6B2]-g6red[CS6M_G62AB];\
          g6red[CS6M_G62BC] = -g6red[CS6M_G62BC]+g6red[CS6M_G62AC];              \
          g6red[CS6M_G62AC] = -g6red[CS6M_G62AC];                                \
          g6red[CS6M_G62AB] = -2.*g6red[CS6M_G6A2]+g6red[CS6M_G62AB];            \
          /* std::cout<<"bdC+ "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "     \
          <<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/            \
          CS6M_VOLCHECK("C+ ",g6red);                                            \
          redpass++; continue;                                                   \
      }                                                                          \
      if (g6red[CS6M_G62AB] < -g6red[CS6M_G6A2]-delta) {                         \
          g6red[CS6M_G6B2] = g6red[CS6M_G6A2]+g6red[CS6M_G6B2]+g6red[CS6M_G62AB];\
          g6red[CS6M_G62BC] =-g6red[CS6M_G62BC]-g6red[CS6M_G62AC];               \
          g6red[CS6M_G62AB] = 2.*g6red[CS6M_G6A2]+g6red[CS6M_G62AB];             \
          g6red[CS6M_G62AC] = -g6red[CS6M_G62AC];                                \
          /* std::cout<<"bdE+ "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "     \
          <<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/            \
          CS6M_VOLCHECK("E+ ",g6red);                                            \
         redpass++; continue;                                                    \
      }                                                                          \
      /* F+ boundary outside */                                                  \
      if (g6red[CS6M_G6A2] + g6red[CS6M_G6B2] + g6red[CS6M_G62BC]                \
          + g6red[CS6M_G62AC] + g6red[CS6M_G62AB] < -delta) {                    \
          g6red[CS6M_G6C2] = g6red[CS6M_G6A2] + g6red[CS6M_G6B2]                 \
          + g6red[CS6M_G6C2] + g6red[CS6M_G62BC] + g6red[CS6M_G62AC] + g6red[CS6M_G62AB]; \
          g6red[CS6M_G62AC] = -2.*g6red[CS6M_G6A2] - g6red[CS6M_G62AC] - g6red[CS6M_G62AB];\
          g6red[CS6M_G62BC] = -2.*g6red[CS6M_G6B2] - g6red[CS6M_G62BC] - g6red[CS6M_G62AB];\
          /* std::cout<<"bdF+ "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "     \
          <<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/            \
          CS6M_VOLCHECK("F+ ",g6red);                                            \
          redpass++; continue;                                                   \
      }                                                                          \
      /* 678 boundaries exact */                                                 \
      if (g6red[CS6M_G62BC] ==  g6red[CS6M_G6B2] && 2.*g6red[CS6M_G62AC] < g6red[CS6M_G62AB]) { \
          g6red[CS6M_G6C2] = g6red[CS6M_G6B2]+g6red[CS6M_G6C2]-g6red[CS6M_G62BC];\
          g6red[CS6M_G62BC] = -2.*g6red[CS6M_G6B2]+g6red[CS6M_G62BC];            \
          g6red[CS6M_G62AC] = g6red[CS6M_G62AC]-g6red[CS6M_G62AB];               \
          /* std::cout<<"bd6 "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "      \
          <<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/            \
          CS6M_VOLCHECK("6= ",g6red);                                            \
          redpass++; continue;                                                   \
      }                                                                          \
      if (fabs(g6red[CS6M_G62BC] + g6red[CS6M_G6B2]) < delta && g6red[CS6M_G62AB] < -delta) {    \
          g6red[CS6M_G6C2] = g6red[CS6M_G6B2]+g6red[CS6M_G6C2]+g6red[CS6M_G62BC];\
          g6red[CS6M_G62BC] = 2.*g6red[CS6M_G6B2]+g6red[CS6M_G62BC];             \
          g6red[CS6M_G62AC] = g6red[CS6M_G62AC]+g6red[CS6M_G62AB];               \
          /* std::cout<<"bd8 "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<    \
          g6red[3 ]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/             \
          CS6M_VOLCHECK("8= ",g6red);                                            \
          redpass++; continue;                                                   \
      }                                                                          \
      /* 9AB boundaries exsct */                                                 \
      if (fabs(g6red[CS6M_G62AC] - g6red[CS6M_G6A2]) < delta && 2.*g6red[CS6M_G62BC] < g6red[CS6M_G62AB]-delta) { \
          g6red[CS6M_G6C2] = g6red[CS6M_G6A2]+g6red[CS6M_G6C2]-g6red[CS6M_G62AC];\
          g6red[CS6M_G62AC] = -2.*g6red[CS6M_G6A2]+g6red[CS6M_G62AC];            \
          g6red[CS6M_G62BC] = g6red[CS6M_G62BC]-g6red[CS6M_G62AB];               \
          /* std::cout<<"bd9 "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "      \
          <<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/            \
          CS6M_VOLCHECK("9= ",g6red);                                            \
          redpass++; continue;                                                   \
      }                                                                          \
      if (g6red[CS6M_G62AC] < -g6red[CS6M_G6A2]-delta && g6red[CS6M_G62AB] < -delta) {     \
          g6red[CS6M_G6C2] = g6red[CS6M_G6A2]+g6red[CS6M_G6C2]+g6red[CS6M_G62AC];\
          g6red[CS6M_G62AC] = 2.*g6red[CS6M_G6A2]+g6red[CS6M_G62AC];             \
          g6red[CS6M_G62BC] = g6red[CS6M_G62BC]+g6red[CS6M_G62AB];               \
          /* std::cout<<"bdB "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "      \
          <<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/            \
          CS6M_VOLCHECK("B= ",g6red);                                            \
          redpass++; continue;                                                   \
      }                                                                          \
      /*  CDE boundaries exact */                                                \
      if (fabs(g6red[CS6M_G62AB] - g6red[CS6M_G6A2]) < delta &&                  \
          2.*g6red[CS6M_G62BC] < g6red[CS6M_G62AC]-delta) {                      \
          g6red[CS6M_G6B2] = g6red[CS6M_G6A2]+g6red[CS6M_G6B2]-g6red[CS6M_G62AB];\
          g6red[CS6M_G62AB] = -2.*g6red[CS6M_G6A2]+g6red[CS6M_G62AB];            \
          g6red[CS6M_G62BC] = -g6red[CS6M_G62BC]+g6red[CS6M_G62AC];              \
          g6red[CS6M_G62AC] = -g6red[CS6M_G62AC];                                \
          /* std::cout<<"bdC "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "      \
          <<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/            \
          CS6M_VOLCHECK("C= ",g6red);                                            \
          redpass++; continue;                                                   \
      }                                                                          \
      if (g6red[CS6M_G62AB] < -g6red[CS6M_G6A2] -delta &&                        \
          g6red[CS6M_G62AB] < -delta) {                                          \
          g6red[CS6M_G6B2] = g6red[CS6M_G6A2]+g6red[CS6M_G6B2]+g6red[CS6M_G62AB];\
          g6red[CS6M_G62AB] = 2.*g6red[CS6M_G6A2]+g6red[CS6M_G62AB];             \
          g6red[CS6M_G62BC] = -g6red[CS6M_G62BC]-g6red[CS6M_G62AC];              \
          g6red[CS6M_G62AC] = -g6red[CS6M_G62AC];                                \
          /* std::cout<<"bdE "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "      \
          <<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/            \
          CS6M_VOLCHECK("E= ",g6red);                                            \
          redpass++; continue;                                                   \
      }                                                                          \
      /* F boundary exact */                                                     \
      if (fabs(g6red[CS6M_G6A2] + g6red[CS6M_G6B2] + g6red[CS6M_G62BC]           \
          + g6red[CS6M_G62AC] + g6red[CS6M_G62AB]) < -delta &&                   \
          2.*g6red[CS6M_G6A2]+2.*g6red[CS6M_G62AC]+g6red[CS6M_G62AB] > delta) {  \
          g6red[CS6M_G6C2] = g6red[CS6M_G6A2] + g6red[CS6M_G6B2]                 \
          + g6red[CS6M_G6C2] + g6red[CS6M_G62BC] + g6red[CS6M_G62AC] + g6red[CS6M_G62AB]; \
          g6red[CS6M_G62AC] = -2.*g6red[CS6M_G6A2] - g6red[CS6M_G62AC] - g6red[CS6M_G62AB];\
          g6red[CS6M_G62BC] = -2.*g6red[CS6M_G6B2] - g6red[CS6M_G62BC] - g6red[CS6M_G62AB];\
          /* std::cout<<"bdF "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "      \
          <<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/            \
          CS6M_VOLCHECK("F= ",g6red);                                            \
          redpass++; continue;                                                   \
      }                                                                          \
      notdone = 0;                                                               \
      reduced = 1;                                                               \
    }                                                                            \
  }

  #define CS6M_Reduce11(din,dout) { \
    dout[0] = -din[0];             \
    dout[1] = din[0] + din[1];     \
    dout[2] = din[0] + din[4];     \
    dout[3] = -din[0] + din[3];    \
    dout[4] = din[0] + din[2];     \
    dout[5] = din[0] + din[5];     \
  }

  #define CS6M_Reduce21(din,dout) { \
    dout[0] = din[1] + din[5];     \
    dout[1] = -din[1];             \
    dout[2] = din[1] + din[2];     \
    dout[3] = din[1] + din[3];     \
    dout[4] = -din[1] + din[4];    \
    dout[5] = din[1] + din[0];     \
  }

  #define CS6M_Reduce31(din,dout) { \
    dout[0] = din[2] + din[0];     \
    dout[1] = din[2] + din[3];     \
    dout[2] = -din[2];             \
    dout[3] = din[2] + din[1];     \
    dout[4] = din[2] + din[4];     \
    dout[5] = -din[2] + din[5];    \
  }

  #define CS6M_Reduce41(din,dout) { \
    dout[0] = -din[3] + din[0];    \
    dout[1] = din[3] + din[2];     \
    dout[2] = din[3] + din[1];     \
    dout[3] = -din[3];             \
    dout[4] = din[3] + din[4];     \
    dout[5] = din[3] + din[5];     \
  }

  #define CS6M_Reduce51(din,dout) { \
    dout[0] = din[4] + din[2];     \
    dout[1] = -din[4] + din[1];    \
    dout[2] = din[4] + din[0];     \
    dout[3] = din[4] + din[3];     \
    dout[4] = -din[4];             \
    dout[5] = din[4] + din[5];     \
  }

  #define CS6M_Reduce61(din,dout) { \
    dout[0] = din[5] + din[1];     \
    dout[1] = din[5] + din[0];     \
    dout[2] = -din[5] + din[2];    \
    dout[3] = din[5] + din[3];     \
    dout[4] = din[5] + din[4];     \
    dout[5] = -din[5];             \
  }

  #define CS6M_S6Reduce(in,out,reduced) {      \
    size_t maxIndex       ;                  \
    int reductionCycleCount;                 \
    int reduction;                           \
    int ii;                                  \
    double maxScalar;                        \
    double delta=0.0;                        \
    double outtemp[6];                       \
    maxIndex=99999;                          \
    reductionCycleCount=0;                   \
    reduction=1;                             \
    for (ii=0; ii<6; ii++) {                 \
      out[ii] = in[ii];                      \
      delta += ((out[ii]>0.)?out[ii]:(-out[ii])); \
    }                                        \
    delta *= 1.e-12;                         \
    while ( (CS6M_CountPositive(out,delta)) != 0) {\
        maxScalar = -1.0E20;                  \
                                              \
        for ( ii = 0; ii < 6; ii++) {         \
           if (out[ii] > maxScalar) {         \
              maxIndex = ii;                  \
              maxScalar = out[ii];            \
           }                                  \
        }                                     \
                                              \
        if (maxIndex < 6) {                   \
          switch (maxIndex) {                       \
            case 0: CS6M_Reduce11(out,outtemp); break;  \
            case 1: CS6M_Reduce21(out,outtemp); break;  \
            case 2: CS6M_Reduce31(out,outtemp); break;  \
            case 3: CS6M_Reduce41(out,outtemp); break;  \
            case 4: CS6M_Reduce51(out,outtemp); break;  \
            case 5: CS6M_Reduce61(out,outtemp); break;  \
          }                                   \
          for (ii=0; ii<6; ii++) out[ii] = outtemp[ii]; \
        }                                               \
        ++reductionCycleCount;                          \
        if (reductionCycleCount > 1000 ) {              \
          reduction=0; break;                           \
        }                                               \
    }                                                   \
    reduced=0;                                          \
    if (CS6M_CountPositive(out,delta) == 0) reduced=1;  \
  }

  

  #define CS6M_D7Reduce(in,out,reduced) {     \
    double s6in[6];                           \
    double s6out[6];                          \
    double temp;                              \
    double delta;                             \
    int ii;                                   \
    CS6M_D7toS6(in,s6in);                     \
    reduced=1;                                \
    CS6M_S6Reduce(s6in,s6out,reduced);        \
    if (reduced) {                            \
       CS6M_S6toD7(s6out,out);                \
    } else {                                  \
      for (ii=0; ii < 7; ii++) out[ii]=in[ii];\
    }                                         \
    delta=out[3];                             \
    if (out[2] > delta) delta=out[2];         \
    if (out[1] > delta) delta=out[1];         \
    if (out[0] > delta) delta=out[0];         \
    delta = delta * 1.e-12;                   \
    if (delta < 1.e-12) delta=1.e-12;         \
    while (out[0] > out[1]+delta || out[1] > out[2]+delta || out[2] > out[3]+delta ) {\
       if (out[2] > out[3]+delta) {           \
         temp = out[2];                       \
         out[2] = out[3];                     \
         out[3]=temp;                         \
         temp = out[4];                       \
         out[4] = out[5];                     \
         out[5] = temp;                       \
         continue;                            \
       }                                      \
       if (out[1] > out[2]+delta) {           \
         temp = out[1];                       \
         out[1] = out[2];                     \
         out[2]=temp;                         \
         temp = out[5];                       \
         out[5] = out[6];                     \
         out[6] = temp;                       \
         continue;                            \
       }                                      \
       if (out[0] > out[1]+delta) {           \
         temp = out[0];                       \
         out[0] = out[1];                     \
         out[1]=temp;                         \
         temp = out[4];                       \
         out[4] = out[5];                     \
         out[5] = temp;                       \
         continue;                            \
       }                                      \
    }                                         \
  }

  #define CS6M_Cellvolume(cell,volume)       {  \
   double cosAlpha=cos(cell[CS6M_CELLALPHA]);  \
   double cosBeta =cos(cell[CS6M_CELLBETA ]);  \
   double cosGamma=cos(cell[CS6M_CELLGAMMA]);  \
   volume=                                    \
     cell[CS6M_CELLA]*cell[CS6M_CELLB]*cell[CS6M_CELLC] \
       * CS6M_sqrt( fabs(                          \
           1.0                                \
           -cosAlpha*cosAlpha                 \
           -cosBeta*cosBeta                   \
           -cosGamma*cosGamma                 \
           +2.*cosAlpha*cosBeta*cosGamma));   \
  }

  #define CS6M_Cellinverse(cell,inversecell) {        \
    const double cosAlpha=cos(cell[CS6M_CELLALPHA]);  \
    const double cosBeta =cos(cell[CS6M_CELLBETA ]);  \
    const double cosGamma=cos(cell[CS6M_CELLGAMMA]);  \
    const double sinAlpha=sin(cell[CS6M_CELLALPHA]);  \
    const double sinBeta =sin(cell[CS6M_CELLBETA ]);  \
    const double sinGamma=sin(cell[CS6M_CELLGAMMA]);  \
                                                      \
    double v;                                         \
    CS6M_Cellvolume(cell,v);                          \
                                                      \
    const double astar = cell[CS6M_CELLB]*cell[CS6M_CELLC]*sinAlpha/v; \
    const double bstar = cell[CS6M_CELLA]*cell[CS6M_CELLC]*sinBeta /v; \
    const double cstar = cell[CS6M_CELLA]*cell[CS6M_CELLB]*sinGamma/v; \
                                                                    \
    const double cosAlphaStar = (cosBeta *cosGamma-cosAlpha)/fabs(sinBeta*sinGamma); \
    const double cosBetaStar  = (cosAlpha*cosGamma-cosBeta )/fabs(sinAlpha*sinGamma);\
    const double cosGammaStar = (cosAlpha*cosBeta -cosGamma)/fabs(sinAlpha*sinBeta); \
                                                                    \
    inversecell[CS6M_CELLA] = astar;                                          \
    inversecell[CS6M_CELLB] = bstar;                                          \
    inversecell[CS6M_CELLC] = cstar;                                          \
    inversecell[CS6M_CELLALPHA] = atan2( fabs(CS6M_sqrt(1.0-cosAlphaStar*cosAlphaStar)), cosAlphaStar);\
    inversecell[CS6M_CELLBETA]  = atan2( fabs(CS6M_sqrt(1.0-cosBetaStar*cosBetaStar)),   cosBetaStar); \
    inversecell[CS6M_CELLGAMMA] = atan2( fabs(CS6M_sqrt(1.0-cosGammaStar*cosGammaStar)), cosGammaStar);\
  }


  #define CS6M_G6toV7(g6vec,v7vec) {                                  \
    double cell[6], inversecell[6];                                   \
    double g6inversecell[6], g6redinversecell[6], redinversecell[6];  \
    double volume;                                                    \
    int reduced;                                                      \
    int ii;                                                           \
    for(ii=0;ii<7;ii++) v7vec[ii]=0;                                  \
    CS6M_G6toCell(g6vec,cell);                                        \
    CS6M_Cellinverse(cell,inversecell);                               \
    CS6M_CelltoG6(inversecell,g6inversecell);                         \
    CS6M_G6Reduce(g6inversecell,g6redinversecell,reduced);            \
    if (reduced) {                                                    \
      CS6M_G6toCell(g6redinversecell,redinversecell);                 \
      CS6M_Cellvolume(cell,volume);                                   \
      v7vec[CS6M_V7A] = cell[CS6M_CELLA];                             \
      v7vec[CS6M_V7B] = cell[CS6M_CELLB];                             \
      v7vec[CS6M_V7C] = cell[CS6M_CELLC];                             \
      v7vec[CS6M_V7ASTARINV] = 1./redinversecell[CS6M_CELLA];         \
      v7vec[CS6M_V7BSTARINV] = 1./redinversecell[CS6M_CELLB];         \
      v7vec[CS6M_V7CSTARINV] = 1./redinversecell[CS6M_CELLC];         \
      v7vec[CS6M_V7VOLCROOT] = pow(volume,1./3.);                     \
    }                                                                 \
  }

  #define CS6M_Dot_Prod(vectorleft,vectorright,dotprod){            \
    int ii;                                                       \
    dotprod=0.;                                                   \
    for (ii=0; ii<6; ii++) dotprod+=vectorleft[ii]*vectorright[ii];\
  }

  #define CS6M_Mat66_Scalar_Mult(scalar,mat66in,mat66out) { \
    int ii;                                               \
    for (ii=0; ii<36; ii++) mat66out=scalar*mat66in;      \
  }


  #ifdef CS6M_DEBUG_PRINT_COUT
    #define CS6M_PRINT_VECTOR6(name,vec) { \
      std::cerr<<name<<": [" << vec[0] << ", " << vec[1] << ", " << vec[2] \
      << ", "<< vec[3] << ", "<< vec[4] << ", "<< vec[5] <<  " ]" << std::endl; \
    }
    #define CS6M_PRINT_MATRIX66(name,mat) { \
      std::cerr<<name<<":" << std::endl \
      << "  [" << mat[0] << ", " << mat[1] << ", " << mat[2] << ", "<< mat[3] << ", "<< mat[4] << ", "<< mat[5] <<  " ]" << std::endl \
      << "  [" << mat[6] << ", " << mat[7] << ", " << mat[8] << ", "<< mat[9] << ", "<< mat[10] << ", "<< mat[11] <<  " ]" << std::endl \
      << "  [" << mat[12] << ", " << mat[13] << ", " << mat[14] << ", "<< mat[15] << ", "<< mat[16] << ", "<< mat[17] <<  " ]" << std::endl \
      << "  [" << mat[18] << ", " << mat[19] << ", " << mat[15] << ", "<< mat[21] << ", "<< mat[22] << ", "<< mat[23] <<  " ]" << std::endl \
      << "  [" << mat[24] << ", " << mat[25] << ", " << mat[26] << ", "<< mat[27] << ", "<< mat[28] << ", "<< mat[29] <<  " ]" << std::endl \
      << "  [" << mat[30] << ", " << mat[31] << ", " << mat[32] << ", "<< mat[33] << ", "<< mat[34] << ", "<< mat[35] <<  " ]" << std::endl; \
    }
  #else
    #ifdef CS6M_DEBUG_PRINT_FPRINTF
      #define CS6M_PRINT_VECTOR6(name,vec) { \
        fprintf(stderr," %s: [%g, %g, %g, %g, %g, %g ]\n",name,vec[0],vec[1],vec[2],vec[3],vec[4],vec[5]);
      }
      #define CS6M_PRINT_MATRIX66(name,mat) { \
        fprintf(stderr," %s: \n [%g, %g, %g, %g, %g, %g ]" \
                          "\n [%g, %g, %g, %g, %g, %g ]" \
                          "\n [%g, %g, %g, %g, %g, %g ]" \
                          "\n [%g, %g, %g, %g, %g, %g ]" \
                          "\n [%g, %g, %g, %g, %g, %g ]" \
                          "\n [%g, %g, %g, %g, %g, %g ]", name, \
                          mat[0],mat[1],mat[2],mat[3],mat[4],mat[5], \
                          mat[6],mat[7],mat[8],mat[9],mat[10],mat[11], \
                          mat[12],mat[13],mat[14],mat[15],mat[16],mat[17], \
                          mat[18],mat[19],mat[20],mat[21],mat[22],mat[23], \
                          mat[24],mat[25],mat[26],mat[27],mat[28],mat[29], \
                          mat[30],mat[31],mat[32],mat[33],mat[34],mat[35]); \
      }
    #else
      #define CS6M_PRINT_VECTOR6(name,vec)  
      #define CS6M_PRINT_MATRIX66(name,mat)  
    #endif
  #endif

  #define CS6M_Mat66_Vector_Mult(mat66in,vectorin,vectorout) { \
    int ii, jj, kk;                                          \
    CS6M_PRINT_MATRIX66(mat66in,mat66in)                     \
    CS6M_PRINT_VECTOR6(vectorin,vectorin)                    \
    for (ii=0; ii<6; ii++) {                                 \
      vectorout[ii]=0.;                                      \
      for (jj=0; jj<6; jj++) {                               \
         kk=ii*6+jj;                                         \
         vectorout[ii]+=mat66in[kk]*vectorin[jj];            \
      }                                                      \
    }                                                        \
    CS6M_PRINT_VECTOR6(vectorout,vectorout)                  \
  }
  #define CS6M_Mat66_Mat66_Mult(mat66inleft,mat66inright,mat66out) { \
    int ii, jj, kk, ll;                                            \
    CS6M_PRINT_MATRIX66(mat66inleft,mat66inleft)                   \
    CS6M_PRINT_MATRIX66(mat66inright,mat66inright)                 \
    for (ii=0; ii<6; ii++){      \
      for (jj=0; jj<6; jj++) {   \
        ll=ii*6+jj;              \
        mat66out[ll]=0;          \
        for (kk=0;kk<6,kk++) {   \
          mat66out[ll]+=mat66inleft[ii*6+kk]*mat66inright[kk*6+jj]; \
        }                        \
      }                          \
    }                            \
    CS6M_PRINT_MATRIX66(mat66out,mat66out)                          \
  }


   static double CS6M_MAT_P[36]={
      1., 0., 0., 0., 0., 0.,
      0., 1., 0., 0., 0., 0.,
      0., 0., 1., 0., 0., 0.,
      0., 0., 0., 1., 0., 0.,
      0., 0., 0., 0., 1., 0.,
      0., 0., 0., 0., 0., 1.
   };
   static double CS6M_MAT_I[36]={
      1., 0., 0., 0., 0., 0.,
      0., 1., 0., 0., 0., 0.,
     .25,.25,.25,.25,.25,.25,
      0., 1., 0., .5, 0., .5,
      1., 0., 0., 0., .5, .5,
      0., 0., 0., 0., 0., 1.
   };
   static double CS6M_MAT_A[36]={
      1., 0., 0., 0., 0., 0.,
      0., 1., 0., 0., 0., 0.,
      0.,.25,.25,.25, 0., 0.,
      0. ,1., 0., .5, 0., 0.,
      0., 0., 0., 0., .5, .5,
      0., 0., 0., 0., 0., .1
   };
   static double CS6M_MAT_B[36]={
      1., 0., 0., 0., 0., 0.,
      0., 1., 0., 0., 0., 0.,
     .25, 0.,.25, 0.,.25, 0.,
      0., 0., 0., .5, 0., .5,
      1., 0., 0., 0., .5, 0.,
      0., 0., 0., 0., 0., 1.
   };
   static double CS6M_MAT_C[36]={
      1., 0., 0., 0., 0., 0.,
     .25,.25, 0., 0., 0.,.25,
      0., 0., 1., 0., 0., 0.,
      0., 0., 0., .5, .5, 0.,
      0., 0., 0., 0., 1., 0.,
      1., 0., 0., 0. ,0.,.5
   };
   static double CS6M_MAT_F[36]={
     .25,.25, 0., 0., 0.,.25,
     .25, 0.,.25, 0.,.25, 0.,
      0.,.25,.25,.25, 0., 0.,
      0., 0., .5,.25,.25,.25,
      0., .5, 0.,.25,.25,.25,
      .5, 0., 0.,.25,.25,.25
   };
   static double CS6M_MAT_H[36]={
      1., 1., 1., 1.,-1.,-1.,
      4., 1., 1., 1., 2., 2.,
      1., 4., 1.,-2.,-1., 2.,
     -4.,-4., 2.,-1., 1.,-5.,
      2.,-4., 2.,-1.,-2., 1.,
     -4., 2., 2., 2., 1., 1
   };
   static double CS6M_MAT_R[36]={
      1.,0.,0.,0.,0.,0.,
      0.,1.,0.,0.,0.,0.,
      0.,0.,1.,0.,0.,0.,
      0.,0.,0.,1.,0.,0.,
      0.,0.,0.,0.,1.,0.,
      0.,0.,0.,0.,0.,1.
   };
   static double CS6M_HEXPERP[36]={
      .6666666666667,-.3333333333333, 0., 0., 0.,.3333333333333,
     -.3333333333333, .6666666666667, 0., 0., 0.,.3333333333333,
      0., 0.,0., 0., 0., 0.,
      0., 0., 0.,-1., 0., 0.,
      0., 0., 0., 0.,-1., 0.,
      .3333333333333, .333333333333, 0., 0., 0., .6666666666667
   };
   static double CS6M_RHMPERP[36]={
      .6666666666667,-.3333333333333,-.3333333333333, 0., 0., 0.,
     -.3333333333333, .6666666666667,-.3333333333333, 0., 0., 0.,
     -.3333333333333,-.3333333333333, .6666666666667, 0., 0., 0.,
      0., 0., 0., .6666666666667,-.3333333333333,-.3333333333333,
      0., 0., 0.,-.3333333333333, .6666666666667,-.3333333333333,
      0., 0., 0.,-.3333333333333,-.3333333333333, .6666666666667 
   };


  #define CS6M_multi_redcells(g6vecin,s6redout,d7redout,g6redout,v7redout) { \
    int reduced;                                                             \
    int ii;                                                                  \
    double s6vecin(6);                                                       \
    double d7vecin(7);                                                       \
    double v7vecin(7);                                                       \
    CS6M_G6toS6(g6vecin,s6vecin);                                            \
    reduced=0;                                                               \
    CS6M_S6Reduce(s6vecin,s6redout,reduced);                                 \
    if (!reduced) for (ii=0;ii<6;ii++) s6redout(ii)=0.;                      \
    CS6M_G6toD7(g6vecin,d7vecin);                                            \
    reduced=0;                                                               \
    CS6M_D7Reduce(d7vecin,d7redout,reduced);                                 \
    if (!reduced) for (ii=0;ii<7;ii++) d7redout(ii)=0.;                      \
    CS6M_G6toV7(g6vecin,v7vecin);                                            \
    reduced=0;                                                               \
    CS6M_V7Reduce(v7vecin,v7redout,reduced);                                 \
    if (!reduced) for (ii=0;ii<7;ii++) v7redout(ii)=0.;                      \
   };


  #define CS6M_LatSymMat66(g6vec,latsym,mat66,g6vecout){\
    int ii;                                           \
    double temp66[36];                                \
    double hexperp[6];                                \
    double rhmperp[6];                                \
    double hexperpnormsq,rhmperpnormsq;               \
    switch (latsym) {                                 \
      case 'P':                                       \
      case 'p':                                       \
        for (ii=0; ii<36; ii++) mat66[ii]= CS6M_MAT_P[ii]; \
        break;                                        \
      case 'I':                                       \
      case 'i':                                       \
        for (ii=0; ii<36; ii++) mat66[ii]= CS6M_MAT_I[ii]; \
        break;                                        \
      case 'A':                                       \
      case 'a':                                       \
        for (ii=0; ii<36; ii++) mat66[ii]= CS6M_MAT_A[ii]; \
        break;                                        \
      case 'B':                                       \
      case 'b':                                       \
        for (ii=0; ii<36; ii++) mat66[ii]= CS6M_MAT_B[ii]; \
        break;                                        \
      case 'C':                                       \
      case 'c':                                       \
        for (ii=0; ii<36; ii++) mat66[ii]= CS6M_MAT_C[ii]; \
        break;                                        \
      case 'F':                                       \
      case 'f':                                       \
         for (ii=0; ii<36; ii++) mat66[ii]= CS6M_MAT_F[ii];\
        break;                                        \
      case 'H':                                       \
      case 'h':                                       \
      case 'R':                                       \
      case 'r':                                       \
        CS6M_Mat66_Vector_Mult(CS6M_HEXPERP,g6vec,hexperp);  \
        CS6M_Mat66_Vector_Mult(CS6M_RHMPERP,g6vec,rhmperp);  \
        CS6M_Dot_Prod(hexperp,hexperp,hexperpnormsq); \
        CS6M_Dot_Prod(rhmperp,rhmperp,rhmperpnormsq); \
        if (hexperpnormsq < rhmperpnormsq) {          \
          for (ii=0; ii<36; ii++) mat66[ii]= CS6M_MAT_H[ii]/9.; \
          break;                                      \
        }                                             \
      default:                                        \
        for (ii=0; ii<36; ii++) mat66[ii]= CS6M_MAT_P[ii];  \
        break;                                        \
    }                                                 \
    CS6M_Mat66_Vector_Mult(mat66,g6vec,g6vecout);     \
    /* std::cout<< "mat66: " << mat66  << std::endl;*/\
    /* std::cout<< "g6vec: " << g6vec  << std::endl;*/\
    /* std::cout<< "g6vecout: " << g6vecout << std::endl; */\
  }

#endif
