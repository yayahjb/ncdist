#ifndef CS6M_SELLING_REDUCE
#define CS6M_SELLING_REDUCE

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

#define CS6M_abs(x) (x<0?-x:x)

#define CS6M_CountPositive(s6vec)  \
    ( ((s6vec[CS6M_S6BC]>0.0)?1:0)  \
    + ((s6vec[CS6M_S6AC]>0.0)?1:0)  \
    + ((s6vec[CS6M_S6AB]>0.0)?1:0)  \
    + ((s6vec[CS6M_S6AD]>0.0)?1:0)  \
    + ((s6vec[CS6M_S6BD]>0.0)?1:0)  \
    + ((s6vec[CS6M_S6CD]>0.0)?1:0)  )

#define CS6M_NegativeSumOfScalars(s6vec) \
     (-s6vec[CS6M_S6BC]-s6vec[CS6M_S6AC]-s6vec[CS6M_S6AB]-s6vec[CS6M_S6AD]-s6vec[CS6M_S6BD]-s6vec[CS6M_S6CD]) \

#define CS6M_D7toS6(d7vec,s6vec) { \
    s6vec[CS6M_S6BC] = -(-d7vec[CS6M_D7_A_D_2]+d7vec[CS6M_D7C2]+d7vec[CS6M_D7B2])/2;\
    s6vec[CS6M_S6AC] = -(-d7vec[CS6M_D7_B_D_2]+d7vec[CS6M_D7C2]+d7vec[CS6M_D7A2])/2;\
    s6vec[CS6M_S6AB] = (-d7vec[CS6M_D7_B_D_2]-d7vec[CS6M_D7_A_D_2]+d7vec[CS6M_D7D2]+d7vec[CS6M_D7C2])/2;\
    s6vec[CS6M_S6AD] = -(-d7vec[CS6M_D7_A_D_2]+d7vec[CS6M_D7D2]+d7vec[CS6M_D7A2])/2;\
    s6vec[CS6M_S6BD] = -(-d7vec[CS6M_D7_B_D_2]+d7vec[CS6M_D7D2]+d7vec[CS6M_D7B2])/2;\
    s6vec[CS6M_S6CD] = (-d7vec[CS6M_D7_B_D_2]-d7vec[CS6M_D7_A_D_2]+d7vec[CS6M_D7B2]+d7vec[CS6M_D7A2])/2;\
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


#define CS6M_G6Reduce(g6unred,g6red,reduced) { \
    double s6in[6];                           \
    double s6out[6];                          \
    int notdone;                              \
    int ii;                                   \
    notdone = 1;                              \
    int redpass;                              \
    double temp;                              \
    reduced=0;                                \
    CS6M_G6toS6(g6unred,s6in);                \
    reduced=1;                                \
    CS6M_S6Reduce(s6in,s6out,reduced);        \
    if (reduced) {                            \
       CS6M_S6toG6(s6out,g6red);              \
    } else {                                  \
      for (ii=0; ii < 6; ii++) g6red[ii]=g6unred[ii];\
    }                                         \
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
      if ( g6red[CS6M_G6A2] > g6red[CS6M_G6B2]                         \
           || ( g6red[CS6M_G6A2] == g6red[CS6M_G6B2]                   \
           && CS6M_abs(g6red[CS6M_G62BC]) > CS6M_abs(g6red[CS6M_G62AC])) ) { \
        temp = g6red[CS6M_G6A2];                                 \
        g6red[CS6M_G6A2] = g6red[CS6M_G6B2];                     \
        g6red[CS6M_G6B2] = temp;                                 \
        temp = g6red[CS6M_G62BC];                                \
        g6red[CS6M_G62BC] = g6red[CS6M_G62AC];                   \
        g6red[CS6M_G62AC] = temp;                                \
        /*std::cout<<"bd1 "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl; */\
        redpass++; continue;                                     \
      }                                                          \
      if ( g6red[CS6M_G6B2] > g6red[CS6M_G6C2]                   \
            ||  (g6red[CS6M_G6B2] == g6red[CS6M_G6C2]            \
            && CS6M_abs(g6red[CS6M_G62AC]) > CS6M_abs(g6red[CS6M_G62AB])) ) {\
        temp = g6red[CS6M_G6B2];                                 \
        g6red[CS6M_G6B2] = g6red[CS6M_G6C2];                     \
        g6red[CS6M_G6C2] = temp;                                 \
        temp = g6red[CS6M_G62AC];                                \
        g6red[CS6M_G62AC] = g6red[CS6M_G62AB];                   \
        g6red[CS6M_G62AB] = temp;                                \
        /*std::cout<<"bd1 "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/\
        redpass++; continue;                                     \
      }                                                          \
      if ( g6red[CS6M_G62BC]*g6red[CS6M_G62AC]*g6red[CS6M_G62AB] > 0. ) { \
        if (g6red[CS6M_G62BC] < 0. || g6red[CS6M_G62AC] < 0. || g6red[CS6M_G62AC] < 0. ) { \
          if ( g6red[CS6M_G62BC] < 0. ) g6red[CS6M_G62BC] = -g6red[CS6M_G62BC]; \
          if ( g6red[CS6M_G62AC] < 0. ) g6red[CS6M_G62AC] = -g6red[CS6M_G62AC]; \
          if ( g6red[CS6M_G62AB] < 0. ) g6red[CS6M_G62AB] = -g6red[CS6M_G62AB]; \
          /*std::cout<<"+++ "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/\
          redpass++; continue;                                                  \
        }                                                                       \
      } else {                                                                  \
        if (g6red[CS6M_G62BC] > 0. || g6red[CS6M_G62AC] > 0. || g6red[CS6M_G62AC] > 0. ) { \
          if ( g6red[CS6M_G62BC] > 0. ) g6red[CS6M_G62BC] = -g6red[CS6M_G62BC];   \
          if ( g6red[CS6M_G62AC] > 0. ) g6red[CS6M_G62AC] = -g6red[CS6M_G62AC];   \
          if ( g6red[CS6M_G62AB] > 0. ) g6red[CS6M_G62AB] = -g6red[CS6M_G62AB];   \
          /*std::cout<<"--- "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/\
          redpass++; continue;                                                   \
        }                                                                        \
      }                                                                          \
                                                                                 \
      /* 678 boundaries outside */                                               \
      if (g6red[CS6M_G62BC] >  g6red[CS6M_G6B2]) {                               \
          g6red[CS6M_G6C2] = g6red[CS6M_G6B2]+g6red[CS6M_G6C2]-g6red[CS6M_G62BC];\
          g6red[CS6M_G62BC] = -2.*g6red[CS6M_G6B2]+g6red[CS6M_G62BC];            \
          g6red[CS6M_G62AC] = g6red[CS6M_G62AC]-g6red[CS6M_G62AB];               \
          /*std::cout<<"bd6+ "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/\
          redpass++; continue;                                                   \
      }                                                                          \
      if (g6red[CS6M_G62BC] < -g6red[CS6M_G6B2]) {                               \
          g6red[CS6M_G6C2] = g6red[CS6M_G6B2]+g6red[CS6M_G6C2]+g6red[CS6M_G62BC];\
          g6red[CS6M_G62BC] = 2.*g6red[CS6M_G6B2]+g6red[CS6M_G62BC];             \
          g6red[CS6M_G62AC] = g6red[CS6M_G62AC]+g6red[CS6M_G62AB];               \
          /*std::cout<<"bd8+"<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<g6red[3 ]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/\
          redpass++; continue;                                                   \
      }                                                                          \
      /* 9AB boundaries outside */                                               \
      if (g6red[CS6M_G62AC] > g6red[CS6M_G6A2]) {                                \
          g6red[CS6M_G6C2] = g6red[CS6M_G6A2]+g6red[CS6M_G6C2]-g6red[CS6M_G62AC];\
          g6red[CS6M_G62AC] = -2.*g6red[CS6M_G6A2]+g6red[CS6M_G62AC];            \
          g6red[CS6M_G62BC] = g6red[CS6M_G62BC]-g6red[CS6M_G62AB];              \
          /*std::cout<<"bd9+ "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/\
          redpass++; continue;                                                   \
      }                                                                          \
      if (g6red[CS6M_G62AC] < -g6red[CS6M_G6A2]) {                               \
          g6red[CS6M_G6C2] = g6red[CS6M_G6A2]+g6red[CS6M_G6C2]+g6red[CS6M_G62AC];\
          g6red[CS6M_G62AC] = 2.*g6red[CS6M_G6A2]+g6red[CS6M_G62AC];             \
          g6red[CS6M_G62BC] = g6red[CS6M_G62BC]+g6red[CS6M_G62AB];               \
          /*std::cout<<"bdB+ "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/\
          redpass++; continue;                                                   \
      }                                                                          \
      /*  CDE boundaries outside */                                              \
      if (g6red[CS6M_G62AB] > g6red[CS6M_G6A2]) {                                \
          g6red[CS6M_G6B2] = g6red[CS6M_G6A2]+g6red[CS6M_G6B2]-g6red[CS6M_G62AB];\
          g6red[CS6M_G62AB] = -2.*g6red[CS6M_G6A2]+g6red[CS6M_G62AB];            \
          g6red[CS6M_G62BC] = g6red[CS6M_G62BC]-g6red[CS6M_G62AB];               \
          /*std::cout<<"bdC+ "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/\
          redpass++; continue;                                                   \
      }                                                                          \
      if (g6red[CS6M_G62AB] < -g6red[CS6M_G6A2]) {                               \
          g6red[CS6M_G6C2] = g6red[CS6M_G6A2]+g6red[CS6M_G6B2]+g6red[CS6M_G62AB];\
          g6red[CS6M_G62AC] = 2.*g6red[CS6M_G6A2]+g6red[CS6M_G62AB];             \
          g6red[CS6M_G62BC] = g6red[CS6M_G62BC]+g6red[CS6M_G62AB];               \
          /*std::cout<<"bdE+ "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/\
          redpass++; continue;                                                   \
      }                                                                          \
      /* F+ boundary outside */                                                  \
      if (g6red[CS6M_G6A2] + g6red[CS6M_G6B2] + g6red[CS6M_G62BC]                \
          + g6red[CS6M_G62AC] + g6red[CS6M_G62AB] < 0.) {                        \
          g6red[CS6M_G6C2] = g6red[CS6M_G6A2] + g6red[CS6M_G6B2]                 \
          + g6red[CS6M_G6C2] + g6red[CS6M_G62BC] + g6red[CS6M_G62AC] + g6red[CS6M_G62AB]; \
          g6red[CS6M_G62AC] = 2.*g6red[CS6M_G6A2] + g6red[CS6M_G62AC] + g6red[CS6M_G62AB];     \
          g6red[CS6M_G62BC] = 2.*g6red[CS6M_G6B2] + g6red[CS6M_G62BC] + g6red[CS6M_G62AB];     \
          /*std::cout<<"bdF+ "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/\
          redpass++; continue;                                                   \
      }                                                                          \
      /* 678 boundaries exact */                                                 \
      if (g6red[CS6M_G62BC] ==  g6red[CS6M_G6B2] && 2.*g6red[CS6M_G62AC] < g6red[CS6M_G62AB]) { \
          g6red[CS6M_G6C2] = g6red[CS6M_G6B2]+g6red[CS6M_G6C2]-g6red[CS6M_G62BC];\
          g6red[CS6M_G62BC] = -2.*g6red[CS6M_G6B2]+g6red[CS6M_G62BC];            \
          g6red[CS6M_G62AC] = g6red[CS6M_G62AC]-g6red[CS6M_G62AB];               \
          /*std::cout<<"bd6 "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/\
          redpass++; continue;                                                   \
      }                                                                          \
      if (g6red[CS6M_G62BC] == -g6red[CS6M_G6B2] && g6red[CS6M_G62AB] < 0.) {    \
          g6red[CS6M_G6C2] = g6red[CS6M_G6B2]+g6red[CS6M_G6C2]+g6red[CS6M_G62BC];\
          g6red[CS6M_G62BC] = 2.*g6red[CS6M_G6B2]+g6red[CS6M_G62BC];             \
          g6red[CS6M_G62AC] = g6red[CS6M_G62AC]+g6red[CS6M_G62AB];               \
          /*std::cout<<"bd8 "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<g6red[3 ]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/\
          redpass++; continue;                                                   \
      }                                                                          \
      /* 9AB boundaries exsct */                                                 \
      if (g6red[CS6M_G62AC] == g6red[CS6M_G6A2] && 2.*g6red[CS6M_G62BC] < g6red[CS6M_G62AB]) { \
          g6red[CS6M_G6C2] = g6red[CS6M_G6A2]+g6red[CS6M_G6C2]-g6red[CS6M_G62AC];\
          g6red[CS6M_G62AC] = -2.*g6red[CS6M_G6A2]+g6red[CS6M_G62AC];            \
          g6red[CS6M_G62BC] = g6red[CS6M_G62BC]-g6red[CS6M_G62AB];              \
          /*std::cout<<"bd9 "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/\
          redpass++; continue;                                                   \
      }                                                                          \
      if (g6red[CS6M_G62AC] < -g6red[CS6M_G6A2] && g6red[CS6M_G62AB] < 0.) {     \
          g6red[CS6M_G6C2] = g6red[CS6M_G6A2]+g6red[CS6M_G6C2]+g6red[CS6M_G62AC];\
          g6red[CS6M_G62AC] = 2.*g6red[CS6M_G6A2]+g6red[CS6M_G62AC];             \
          g6red[CS6M_G62BC] = g6red[CS6M_G62BC]+g6red[CS6M_G62AB];               \
          /*std::cout<<"bdB "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/\
          redpass++; continue;                                                   \
      }                                                                          \
      /*  CDE boundaries exact */                                                \
      if (g6red[CS6M_G62AB] == g6red[CS6M_G6A2] && 2.*g6red[CS6M_G62BC] < g6red[CS6M_G62AC]) {  \
          g6red[CS6M_G6B2] = g6red[CS6M_G6A2]+g6red[CS6M_G6B2]-g6red[CS6M_G62AB];\
          g6red[CS6M_G62AB] = -2.*g6red[CS6M_G6A2]+g6red[CS6M_G62AB];            \
          g6red[CS6M_G62BC] = g6red[CS6M_G62BC]-g6red[CS6M_G62AB];               \
          /*std::cout<<"bdC "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/\
          redpass++; continue;                                                   \
      }                                                                          \
      if (g6red[CS6M_G62AB] == -g6red[CS6M_G6A2] &&  g6red[CS6M_G62AB] < 0.) {   \
          g6red[CS6M_G6C2] = g6red[CS6M_G6A2]+g6red[CS6M_G6B2]+g6red[CS6M_G62AB];\
          g6red[CS6M_G62AC] = 2.*g6red[CS6M_G6A2]+g6red[CS6M_G62AB];             \
          g6red[CS6M_G62BC] = g6red[CS6M_G62BC]+g6red[CS6M_G62AB];               \
          /*std::cout<<"bdE "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/\
          redpass++; continue;                                                   \
      }                                                                          \
      /* F boundary exact */                                                     \
      if (g6red[CS6M_G6A2] + g6red[CS6M_G6B2] + g6red[CS6M_G62BC]                \
          + g6red[CS6M_G62AC] + g6red[CS6M_G62AB] == 0. && 2.*g6red[CS6M_G6A2]+2.*g6red[CS6M_G62AC]+g6red[CS6M_G62AB]> 0.) { \
          g6red[CS6M_G6C2] = g6red[CS6M_G6A2] + g6red[CS6M_G6B2]                 \
          + g6red[CS6M_G6C2] + g6red[CS6M_G62BC] + g6red[CS6M_G62AC] + g6red[CS6M_G62AB]; \
          g6red[CS6M_G62AC] = 2.*g6red[CS6M_G6A2] + g6red[CS6M_G62AC] + g6red[CS6M_G62AB];\
          g6red[CS6M_G62BC] = 2.*g6red[CS6M_G6B2] + g6red[CS6M_G62BC] + g6red[CS6M_G62AB];\
          /*std::cout<<"bdF "<<g6red[0]<<" "<<g6red[1]<<" "<<g6red[2]<<" "<<g6red[3]<<" "<<g6red[4]<<" "<<g6red[5]<<" "<<std::endl;*/\
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

#define CS6M_S6Reduce(in,out,reduced) {     \
    unsigned long maxIndex;                  \
    int reductionCycleCount;                 \
    int reduction;                           \
    int ii;                                  \
    double maxScalar;                        \
    double outtemp[6];                       \
    maxIndex=99999;                          \
    reductionCycleCount=0;                   \
    reduction=1;                             \
    for (ii=0; ii<6; ii++) out[ii] = in[ii]; \
                                             \
    while ( (CS6M_CountPositive(out)) != 0) {\
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
    if (CS6M_CountPositive(out) == 0) reduced=1;        \
}

  

#define CS6M_D7Reduce(in,out,reduced) {     \
    double s6in[6];                           \
    double s6out[6];                          \
    double temp;                              \
    int ii;                                   \
    CS6M_D7toS6(in,s6in);                     \
    reduced=1;                                \
    CS6M_S6Reduce(s6in,s6out,reduced);        \
    if (reduced) {                            \
       CS6M_S6toD7(s6out,out);                \
    } else {                                  \
      for (ii=0; ii < 7; ii++) out[ii]=in[ii];\
    }                                         \
    while (out[0] > out[1] || out[1] > out[2] || out[2] > out[3] ) {\
       if (out[2] > out[3]) {                 \
         temp = out[2];                       \
         out[2] = out[3];                     \
         out[3]=temp;                         \
         temp = out[4];                       \
         out[4] = out[5];                     \
         out[5] = temp;                       \
         continue;                            \
       }                                      \
       if (out[1] > out[2]) {                 \
         temp = out[1];                       \
         out[1] = out[2];                     \
         out[2]=temp;                         \
         temp = out[5];                       \
         out[5] = out[6];                     \
         out[6] = temp;                       \
         continue;                            \
       }                                      \
       if (out[0] > out[1]) {                 \
         temp = out[0];                       \
         out[0] = out[1];                     \
         out[1]=temp;                         \
         temp = out[4];                       \
         out[4] = out[5];                     \
         out[5] = temp;                       \
         continue;                            \
       }                                      \
    }                                        \
}
#endif
