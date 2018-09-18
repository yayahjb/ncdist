#ifndef CS6M_SELLING_REDUCE
#define CS6M_SELLING_REDUCE

#define CS6M_CountPositive(s6vec)  \
    ( ((s6vec[0]>0.0)?1:0)  \
    + ((s6vec[1]>0.0)?1:0)  \
    + ((s6vec[2]>0.0)?1:0)  \
    + ((s6vec[3]>0.0)?1:0)  \
    + ((s6vec[4]>0.0)?1:0)  \
    + ((s6vec[5]>0.0)?1:0)  )

#define CS6M_NegativeSumOfScalars(s6vec) \
     (-s6vec[0]-s6vec[1]-s6vec[2]-s6vec[3]-s6vec[4]-s6vec[5]) \

#define CS6M_S6toD7(s6vec,d7vec) { \
    d7vec[0] = -s6vec[3] -s6vec[2] -s6vec[1];\
    d7vec[1] = -s6vec[4] -s6vec[2] -s6vec[0];\
    d7vec[2] = -s6vec[5] -s6vec[1] -s6vec[0];\
    d7vec[3] = -s6vec[3] -s6vec[4] -s6vec[5];\
    d7vec[4] = -s6vec[4] -s6vec[2] -s6vec[5] -s6vec[1];\
    d7vec[5] = -s6vec[3] -s6vec[2] -s6vec[5] -s6vec[0];\
    d7vec[6] = -s6vec[3] -s6vec[1] -s6vec[4] -s6vec[0];\
}

#define CS6M_G6toD7(g6vec,d7vec) { \
    d7vec[0] = g6vec[0]; \
    d7vec[1] = g6vec[1]; \
    d7vec[2] = g6vec[2]; \
    d7vec[3] = g6vec[0]+g6vec[1]+g6vec[2]+g6vec[3]+g6vec[4]+g6vec[5]; \
    d7vec[4] = g6vec[1]+g6vec[2]+g6vec[3]; \
    d7vec[5] = g6vec[0]+g6vec[2]+g6vec[4]; \
    d7vec[6] = g6vec[0]+g6vec[1]+g6vec[5]; \
}

#define CS6M_D7toS6(d7vec,s6vec) { \
    s6vec[0] = -(-d7vec[4]+d7vec[2]+d7vec[1])/2;\
    s6vec[1] = -(-d7vec[5]+d7vec[2]+d7vec[0])/2;\
    s6vec[2] = (-d7vec[5]-d7vec[4]+d7vec[3]+d7vec[2])/2;\
    s6vec[3] = -(-d7vec[4]+d7vec[3]+d7vec[0])/2;\
    s6vec[4] = -(-d7vec[5]+d7vec[3]+d7vec[1])/2;\
    s6vec[5] = (-d7vec[5]-d7vec[4]+d7vec[1]+d7vec[0])/2;\
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
        if (reductionCycleCount > 1000 ) { \
          reduction=0; break;                           \
        }                                               \
    }                                                   \
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
    if (0) while (out[0] > out[1] || out[1] > out[2] || out[2] > out[3] ) {\
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
