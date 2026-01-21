#ifndef MACROS_HPP
#define MACROS_HPP

#define FOR_II for (int ii = 0; ii < nn; ++ii)

#define FOR_IP for (int ip = 0; ip < np; ++ip)

#define FOR_INO for (int ino = 0; ino < nno; ++ino)
#define FOR_INO_A for (int ino = 0; ino < nno_a; ++ino)
#define FOR_INO_G_ONLY for (int ino = nno_a; ino < nno; ++ino)

#define FOR_IFA for (int ifa = 0; ifa < nfa; ++ifa)
#define FOR_IBF for (int ibf = 0; ibf < nbf; ++ibf)
#define FOR_IEF for (int ief = 0; ief < nef; ++ief)

#define FOR_IED for (int ied = 0; ied < ned; ++ied)

// simple surface
#define FOR_ISP for (int isp = 0; isp < nsp; ++isp)
#define FOR_IST for (int ist = 0; ist < nst; ++ist)

// tets...
#define FOR_ITE for (int ite = 0; ite < nte; ++ite)

// get rid of this one eventually...
#define FOR_BOUNDARY_IFA for (int ifa = nfa_i; ifa < nfa; ++ifa)

#define FOR_INTERNAL_IFA for (int ifa = 0; ifa < nfa_i; ++ifa)
#define FOR_INTERPROC_IFA for (int ifa = nfa_i; ifa < nfa; ++ifa)

#define FOR_INTERNAL_IEF for (int ief = 0; ief < nef_i; ++ief)
#define FOR_INTERPROC_IEF for (int ief = nef_i; ief < nef; ++ief)

#define FOR_IZONE(VEC) for (int izone = 0, nzone = VEC.size(); izone < nzone; ++izone)

#define FOR_ICV for (int icv = 0; icv < ncv; ++icv)
#define FOR_ICV_REVERSE for (int icv = ncv-1; icv >= 0; --icv)
#define FOR_ICV_G for (int icv = 0; icv < ncv_g; ++icv)
#define FOR_ICV_G2 for (int icv = 0; icv < ncv_g2; ++icv)
#define FOR_ICV_G_ONLY for (int icv = ncv; icv < ncv_g; ++icv)
#define FOR_GHOST_ICV for (int icv = ncv; icv < ncv_g; ++icv)

//#define FOR_IP(T) for (unsigned int ip = 0; ip < (T)->size(); ++ip)

#define FOR_CVZONE for (vector<CvZone>::iterator cvzone = cvZoneVec.begin(); cvzone != cvZoneVec.end(); ++cvzone)
#define FOR_ZONE for(typename vector<Zone<StateT,RhsT>* >::iterator zone = bf_zones.begin(); zone != bf_zones.end(); ++zone)
#define FOR_STATS for(vector<Stats*>::iterator it = stats_vec.begin(); it != stats_vec.end(); ++it)

#define FOR_SCALAR_DATA(ptr) for(list<CtiRealScalar>::iterator it = ptr->sdata_vec.begin(); it != ptr->sdata_vec.end(); ++it)
#define FOR_VECTOR_DATA(ptr) for(list<CtiRealVector>::iterator it = ptr->vdata_vec.begin(); it != ptr->vdata_vec.end(); ++it)
#define FOR_IVALUE_DATA(ptr) for(list<CtiIntValue>::iterator it = ptr->ivalue_vec.begin(); it != ptr->ivalue_vec.end(); ++it)
#define FOR_DVALUE_DATA(ptr) for(list<CtiDoubleValue>::iterator it = ptr->dvalue_vec.begin(); it != ptr->dvalue_vec.end(); ++it)

#define IF_RANK0 if (mpi_rank == 0)
#define FOR_RANK for (int rank = 0; rank < mpi_size; ++rank)

#define FOR_IGR for (int igr = 0; igr < ngr; ++igr)

#define FOR_NOF for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof)
#define FOR_FOC for (int foc = faocv_i[icv]; foc != faocv_i[icv+1]; ++foc)

#define FOR_I2 for (int i = 0; i < 2; ++i)
#define FOR_J2 for (int j = 0; j < 2; ++j)
#define FOR_K2 for (int k = 0; k < 2; ++k)
#define FOR_L2 for (int l = 0; l < 2; ++l)
#define FOR_M2 for (int m = 0; m < 2; ++m)

#define FOR_I3 for (int i = 0; i < 3; ++i)
#define FOR_J3 for (int j = 0; j < 3; ++j)
#define FOR_K3 for (int k = 0; k < 3; ++k)
#define FOR_L3 for (int l = 0; l < 3; ++l)
#define FOR_M3 for (int m = 0; m < 3; ++m)
#define FOR_N3 for (int n = 0; n < 3; ++n)

#define FOR_I4 for (int i = 0; i < 4; ++i)
#define FOR_J4 for (int j = 0; j < 4; ++j)
#define FOR_K4 for (int k = 0; k < 4; ++k)

#define FOR_I5 for (int i = 0; i < 5; ++i)
#define FOR_J5 for (int j = 0; j < 5; ++j)
#define FOR_K5 for (int k = 0; k < 5; ++k)

#define FOR_I6 for (int i = 0; i < 6; ++i)
#define FOR_J6 for (int j = 0; j < 6; ++j)

#define FOR_I7 for (int i = 0; i < 7; ++i)

#define FOR_I8 for (int i = 0; i < 8; ++i)

#define FOR_I10 for (int i = 0; i < 10; ++i)
#define FOR_J10 for (int j = 0; j < 10; ++j)

#define FOR_I12 for (int i = 0; i < 12; ++i)

#define FOR_I16 for (int i = 0; i < 16; ++i)

#define FOR_I32 for (int i = 0; i < 32; ++i)

#define LOOP_l_N(SIZE) for (int l = 0 ; l < (SIZE); ++l)

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#define SET3(A,B) { (A)[0] = (B)[0] ; (A)[1] = (B)[1] ; (A)[2] = (B)[2] ; }
#define GET3(A) { (A)[0] , (A)[1] , (A)[2] }

#define MAX3(A,B,C) max(max((A),(B)),(C))
#define MAX4(A,B,C,D) max(max(max((A),(B)),(C)),(D))
#define MAX5(A,B,C,D,E) max(max(max(max((A),(B)),(C)),(D)),(E))
#define MAX6(A,B,C,D,E,F) max(max(max(max(max((A),(B)),(C)),(D)),(E)),(F))
#define MAX7(A,B,C,D,E,F,G) max(max(max(max(max(max((A),(B)),(C)),(D)),(E)),(F)),(G))
#define MAX8(A,B,C,D,E,F,G,H) max(max(max(max(max(max(max((A),(B)),(C)),(D)),(E)),(F)),(G)),(H))
#define MAX9(A,B,C,D,E,F,G,H,I) max(max(max(max(max(max(max(max((A),(B)),(C)),(D)),(E)),(F)),(G)),(H)),(I))
#define MAX10(A,B,C,D,E,F,G,H,I,J) max(max(max(max(max(max(max(max(max((A),(B)),(C)),(D)),(E)),(F)),(G)),(H)),(I)),(J))
#define MAX11(A,B,C,D,E,F,G,H,I,J,K) max(max(max(max(max(max(max(max(max(max((A),(B)),(C)),(D)),(E)),(F)),(G)),(H)),(I)),(J)),(K))
#define MAX12(A,B,C,D,E,F,G,H,I,J,K,L) max(max(max(max(max(max(max(max(max(max(max((A),(B)),(C)),(D)),(E)),(F)),(G)),(H)),(I)),(J)),(K)),(L))

#define MIN3(A,B,C) min(min((A),(B)),(C))
#define MIN4(A,B,C,D) min(min(min((A),(B)),(C)),(D))
#define MIN5(A,B,C,D,E) min(min(min(min((A),(B)),(C)),(D)),(E))
#define MIN6(A,B,C,D,E,F) min(min(min(min(min((A),(B)),(C)),(D)),(E)),(F))
#define MIN7(A,B,C,D,E,F,G) min(min(min(min(min(min((A),(B)),(C)),(D)),(E)),(F)),(G))
#define MIN8(A,B,C,D,E,F,G,H) min(min(min(min(min(min(min((A),(B)),(C)),(D)),(E)),(F)),(G)),(H))
#define MIN9(A,B,C,D,E,F,G,H,I) min(min(min(min(min(min(min(min((A),(B)),(C)),(D)),(E)),(F)),(G)),(H)),(I))
#define MIN10(A,B,C,D,E,F,G,H,I,J) min(min(min(min(min(min(min(min(min((A),(B)),(C)),(D)),(E)),(F)),(G)),(H)),(I)),(J))
#define MIN11(A,B,C,D,E,F,G,H,I,J,K) min(min(min(min(min(min(min(min(min(min((A),(B)),(C)),(D)),(E)),(F)),(G)),(H)),(I)),(J)),(K))
#define MIN12(A,B,C,D,E,F,G,H,I,J,K,L) min(min(min(min(min(min(min(min(min(min(min((A),(B)),(C)),(D)),(E)),(F)),(G)),(H)),(I)),(J)),(K)),(L))

// returns the scalar product of 2 3-vectors...
#define DOT_PRODUCT(A,B) ((A)[0]*(B)[0] + (A)[1]*(B)[1] + (A)[2]*(B)[2])
#define DOT_PRODUCT_2D(A,B) ((A)[0]*(B)[0] + (A)[1]*(B)[1])

// returns magnitude of a vector
#define MAG(A) (sqrt((A)[0]*(A)[0] + (A)[1]*(A)[1] + (A)[2]*(A)[2]))
#define MAG_2D(A) (sqrt((A)[0]*(A)[0] + (A)[1]*(A)[1]))

// use this for initializing a 3-vector...
// e.g. double v[3] = CROSS_PRODUCT(v1,v2);
#define CROSS_PRODUCT(A,B) { (A)[1]*(B)[2]-(A)[2]*(B)[1] , (A)[2]*(B)[0]-(A)[0]*(B)[2] , (A)[0]*(B)[1]-(A)[1]*(B)[0] }
#define CROSS_PRODUCT_0(A,B) ((A)[1]*(B)[2]-(A)[2]*(B)[1])
#define CROSS_PRODUCT_1(A,B) ((A)[2]*(B)[0]-(A)[0]*(B)[2])
#define CROSS_PRODUCT_2(A,B) ((A)[0]*(B)[1]-(A)[1]*(B)[0])
#define CROSS_PRODUCT_MAG(A,B) sqrt(pow((A)[1]*(B)[2]-(A)[2]*(B)[1],2) + pow((A)[2]*(B)[0]-(A)[0]*(B)[2],2) + pow((A)[0]*(B)[1]-(A)[1]*(B)[0],2))
#define CROSS_DOT(A,B,C) ( ((A)[1]*(B)[2]-(A)[2]*(B)[1])*(C)[0] + ((A)[2]*(B)[0]-(A)[0]*(B)[2])*(C)[1] + ((A)[0]*(B)[1]-(A)[1]*(B)[0])*(C)[2] )

#define SET_CROSS_PRODUCT(S,A,B) {					\
    (S)[0] = (A)[1]*(B)[2]-(A)[2]*(B)[1];				\
    (S)[1] = (A)[2]*(B)[0]-(A)[0]*(B)[2];				\
    (S)[2] = (A)[0]*(B)[1]-(A)[1]*(B)[0];  }

#define NORMALIZE(V) {							\
    const double mag = sqrt( (V)[0]*(V)[0] + (V)[1]*(V)[1] + (V)[2]*(V)[2] ); \
    (V)[0] /= mag;							\
    (V)[1] /= mag;							\
    (V)[2] /= mag;  }

// use this for matrix[3,3]*vector[3]...
// e.g. double v[3] = MATRIX_PRODUCT(M,v1);
#define MATRIX_PRODUCT(M,A) { (M)[0][0]*(A)[0]+(M)[0][1]*(A)[1]+(M)[0][2]*(A)[2] , (M)[1][0]*(A)[0]+(M)[1][1]*(A)[1]+(M)[1][2]*(A)[2] , (M)[2][0]*(A)[0]+(M)[2][1]*(A)[1]+(M)[2][2]*(A)[2] }
#define MAT_VEC_MULT(M,A) MATRIX_PRODUCT(M,A)

// use this for matrix[3,3]*matrix[3,3]...
// e.g. double v[3][3] = MATRIX_MULT(M1,M2);
#define MATRIX_MULT(A,B) {\
  { (A)[0][0]*(B)[0][0] + (A)[0][1]*(B)[1][0] + (A)[0][2]*(B)[2][0], (A)[0][0]*(B)[0][1] + (A)[0][1]*(B)[1][1] + (A)[0][2]*(B)[2][1], (A)[0][0]*(B)[0][2] + (A)[0][1]*(B)[1][2] + (A)[0][2]*(B)[2][2] }, \
  { (A)[1][0]*(B)[0][0] + (A)[1][1]*(B)[1][0] + (A)[1][2]*(B)[2][0], (A)[1][0]*(B)[0][1] + (A)[1][1]*(B)[1][1] + (A)[1][2]*(B)[2][1], (A)[1][0]*(B)[0][2] + (A)[1][1]*(B)[1][2] + (A)[1][2]*(B)[2][2] }, \
  { (A)[2][0]*(B)[0][0] + (A)[2][1]*(B)[1][0] + (A)[2][2]*(B)[2][0], (A)[2][0]*(B)[0][1] + (A)[2][1]*(B)[1][1] + (A)[2][2]*(B)[2][1], (A)[2][0]*(B)[0][2] + (A)[2][1]*(B)[1][2] + (A)[2][2]*(B)[2][2] } \
}
#define MAT_MAT_MULT(A,B) MATRIX_MULT(A,B)

// use this for transpose(matrix[3,3])*matrix[3,3]...
// e.g. double v[3][3] = TMAT_MAT_MULT(M1,M2);
#define TMAT_MAT_MULT(A,B) {\
  { (A)[0][0]*(B)[0][0] + (A)[1][0]*(B)[0][1] + (A)[2][0]*(B)[0][2], (A)[0][0]*(B)[1][0] + (A)[1][0]*(B)[1][1] + (A)[2][0]*(B)[1][2], (A)[0][0]*(B)[2][0] + (A)[1][0]*(B)[2][1] + (A)[2][0]*(B)[2][2] }, \
  { (A)[0][1]*(B)[0][0] + (A)[1][1]*(B)[0][1] + (A)[2][1]*(B)[0][2], (A)[0][1]*(B)[1][0] + (A)[1][1]*(B)[1][1] + (A)[2][1]*(B)[1][2], (A)[0][1]*(B)[2][0] + (A)[1][1]*(B)[2][1] + (A)[2][1]*(B)[2][2] }, \
  { (A)[0][2]*(B)[0][0] + (A)[1][2]*(B)[0][1] + (A)[2][2]*(B)[0][2], (A)[0][2]*(B)[1][0] + (A)[1][2]*(B)[1][1] + (A)[2][2]*(B)[1][2], (A)[0][2]*(B)[2][0] + (A)[1][2]*(B)[2][1] + (A)[2][2]*(B)[2][2] } \
}

// use this for matrix[3,3]*transpose(matrix[3,3])...
// e.g. double v[3][3] = MAT_TMAT_MULT(M1,M2);
#define MAT_TMAT_MULT(A,B) {\
  { (A)[0][0]*(B)[0][0] + (A)[0][1]*(B)[0][1] + (A)[0][2]*(B)[0][2], (A)[0][0]*(B)[1][0] + (A)[0][1]*(B)[1][1] + (A)[0][2]*(B)[1][2], (A)[0][0]*(B)[2][0] + (A)[0][1]*(B)[2][1] + (A)[0][2]*(B)[2][2] }, \
  { (A)[1][0]*(B)[0][0] + (A)[1][1]*(B)[0][1] + (A)[1][2]*(B)[0][2], (A)[1][0]*(B)[1][0] + (A)[1][1]*(B)[1][1] + (A)[1][2]*(B)[1][2], (A)[1][0]*(B)[2][0] + (A)[1][1]*(B)[2][1] + (A)[1][2]*(B)[2][2] }, \
  { (A)[2][0]*(B)[0][0] + (A)[2][1]*(B)[0][1] + (A)[2][2]*(B)[0][2], (A)[2][0]*(B)[1][0] + (A)[2][1]*(B)[1][1] + (A)[2][2]*(B)[1][2], (A)[2][0]*(B)[2][0] + (A)[2][1]*(B)[2][1] + (A)[2][2]*(B)[2][2] } \
}

// use this for matrix[3,3]*diagonal[3]*transpose(matrix[3,3])...
// e.g. double v[3][3] = MAT_DIAG_TMAT_MULT(M1,D,M2);
// note that the 0 off-diagonals are not included in D[3]
#define MAT_DIAG_TMAT_MULT(A,D,B) {\
  { (A)[0][0]*(D)[0]*(B)[0][0] + (A)[0][1]*(D)[1]*(B)[0][1] + (A)[0][2]*(D)[2]*(B)[0][2], (A)[0][0]*(D)[0]*(B)[1][0] + (A)[0][1]*(D)[1]*(B)[1][1] + (A)[0][2]*(D)[2]*(B)[1][2], (A)[0][0]*(D)[0]*(B)[2][0] + (A)[0][1]*(D)[1]*(B)[2][1] + (A)[0][2]*(D)[2]*(B)[2][2] }, \
  { (A)[1][0]*(D)[0]*(B)[0][0] + (A)[1][1]*(D)[1]*(B)[0][1] + (A)[1][2]*(D)[2]*(B)[0][2], (A)[1][0]*(D)[0]*(B)[1][0] + (A)[1][1]*(D)[1]*(B)[1][1] + (A)[1][2]*(D)[2]*(B)[1][2], (A)[1][0]*(D)[0]*(B)[2][0] + (A)[1][1]*(D)[1]*(B)[2][1] + (A)[1][2]*(D)[2]*(B)[2][2] }, \
  { (A)[2][0]*(D)[0]*(B)[0][0] + (A)[2][1]*(D)[1]*(B)[0][1] + (A)[2][2]*(D)[2]*(B)[0][2], (A)[2][0]*(D)[0]*(B)[1][0] + (A)[2][1]*(D)[1]*(B)[1][1] + (A)[2][2]*(D)[2]*(B)[1][2], (A)[2][0]*(D)[0]*(B)[2][0] + (A)[2][1]*(D)[1]*(B)[2][1] + (A)[2][2]*(D)[2]*(B)[2][2] } \
}

// determinant of A[3][3]...
#define DETERMINANT(A) ( (A)[0][0]*((A)[1][1]*(A)[2][2] - (A)[2][1]*(A)[1][2]) \
                        -(A)[0][1]*((A)[1][0]*(A)[2][2] - (A)[2][0]*(A)[1][2]) \
                        +(A)[0][2]*((A)[1][0]*(A)[2][1] - (A)[2][0]*(A)[1][1]) )

// determinant of A[4][4]...
#define DETERMINANT_4(A) \
  ((A)[0][3] * (A)[1][2] * (A)[2][1] * (A)[3][0] - (A)[0][2] * (A)[1][3] * (A)[2][1] * (A)[3][0] - \
   (A)[0][3] * (A)[1][1] * (A)[2][2] * (A)[3][0] + (A)[0][1] * (A)[1][3] * (A)[2][2] * (A)[3][0] + \
   (A)[0][2] * (A)[1][1] * (A)[2][3] * (A)[3][0] - (A)[0][1] * (A)[1][2] * (A)[2][3] * (A)[3][0] - \
   (A)[0][3] * (A)[1][2] * (A)[2][0] * (A)[3][1] + (A)[0][2] * (A)[1][3] * (A)[2][0] * (A)[3][1] + \
   (A)[0][3] * (A)[1][0] * (A)[2][2] * (A)[3][1] - (A)[0][0] * (A)[1][3] * (A)[2][2] * (A)[3][1] - \
   (A)[0][2] * (A)[1][0] * (A)[2][3] * (A)[3][1] + (A)[0][0] * (A)[1][2] * (A)[2][3] * (A)[3][1] + \
   (A)[0][3] * (A)[1][1] * (A)[2][0] * (A)[3][2] - (A)[0][1] * (A)[1][3] * (A)[2][0] * (A)[3][2] - \
   (A)[0][3] * (A)[1][0] * (A)[2][1] * (A)[3][2] + (A)[0][0] * (A)[1][3] * (A)[2][1] * (A)[3][2] + \
   (A)[0][1] * (A)[1][0] * (A)[2][3] * (A)[3][2] - (A)[0][0] * (A)[1][1] * (A)[2][3] * (A)[3][2] - \
   (A)[0][2] * (A)[1][1] * (A)[2][0] * (A)[3][3] + (A)[0][1] * (A)[1][2] * (A)[2][0] * (A)[3][3] + \
   (A)[0][2] * (A)[1][0] * (A)[2][1] * (A)[3][3] - (A)[0][0] * (A)[1][2] * (A)[2][1] * (A)[3][3] - \
   (A)[0][1] * (A)[1][0] * (A)[2][2] * (A)[3][3] + (A)[0][0] * (A)[1][1] * (A)[2][2] * (A)[3][3])
// use DIFF to initialize a 3-vector with a difference of 2 3-vectors
#define DIFF(A,B) { (A)[0]-(B)[0], (A)[1]-(B)[1], (A)[2]-(B)[2] }
#define DIFF_2D(A,B) { (A)[0]-(B)[0], (A)[1]-(B)[1] }

// use this for outputing a 3-vector...
#define COUT_VEC(X) "[ " << (X)[0] << " " << (X)[1] << " " << (X)[2] << " ]"
#define COUT_VEC_2D(X) "[ " << (X)[0] << " " << (X)[1] << " ]"

// 6 times the signed tet volume of 4 points represented as 3-vectors: [A[0],A[1],A[2]], etc...
#define SIGNED_TET_VOLUME_6(A,B,C,D) (((B)[0]-(A)[0])*(((C)[1]-(A)[1])*((D)[2]-(A)[2])-((C)[2]-(A)[2])*((D)[1]-(A)[1])) + \
				      ((B)[1]-(A)[1])*(((C)[2]-(A)[2])*((D)[0]-(A)[0])-((C)[0]-(A)[0])*((D)[2]-(A)[2])) + \
				      ((B)[2]-(A)[2])*(((C)[0]-(A)[0])*((D)[1]-(A)[1])-((C)[1]-(A)[1])*((D)[0]-(A)[0])))

// 6 times the tet volume of 3 integer points, cast to int8...
#define SIGNED_TET_VOLUME_6_INT8(A,B,C,D) (int8((B)[0]-(A)[0])*(int8((C)[1]-(A)[1])*int8((D)[2]-(A)[2])-int8((C)[2]-(A)[2])*int8((D)[1]-(A)[1])) + \
					   int8((B)[1]-(A)[1])*(int8((C)[2]-(A)[2])*int8((D)[0]-(A)[0])-int8((C)[0]-(A)[0])*int8((D)[2]-(A)[2])) + \
					   int8((B)[2]-(A)[2])*(int8((C)[0]-(A)[0])*int8((D)[1]-(A)[1])-int8((C)[1]-(A)[1])*int8((D)[0]-(A)[0])))

// twice the tri normal of 3 points...
#define TRI_NORMAL_2(A,B,C) { ((B)[1]-(A)[1])*((C)[2]-(A)[2])-((B)[2]-(A)[2])*((C)[1]-(A)[1]), \
      ((B)[2]-(A)[2])*((C)[0]-(A)[0])-((B)[0]-(A)[0])*((C)[2]-(A)[2]),	\
      ((B)[0]-(A)[0])*((C)[1]-(A)[1])-((B)[1]-(A)[1])*((C)[0]-(A)[0]) }

// twice the tri normal of 3 integer points, cast to int8...
#define TRI_NORMAL_2_INT8(A,B,C) { int8((B)[1]-(A)[1])*int8((C)[2]-(A)[2])-int8((B)[2]-(A)[2])*int8((C)[1]-(A)[1]), \
      int8((B)[2]-(A)[2])*int8((C)[0]-(A)[0])-int8((B)[0]-(A)[0])*int8((C)[2]-(A)[2]), \
      int8((B)[0]-(A)[0])*int8((C)[1]-(A)[1])-int8((B)[1]-(A)[1])*int8((C)[0]-(A)[0]) }

// delta version of the above: i.e. first point in tet is at 0,0,0
#define SIGNED_TET_VOLUME_6_VEC(B,C,D) (((B)[0])*(((C)[1])*((D)[2])-((C)[2])*((D)[1])) + \
					((B)[1])*(((C)[2])*((D)[0])-((C)[0])*((D)[2])) + \
					((B)[2])*(((C)[0])*((D)[1])-((C)[1])*((D)[0])))

#define DIST(A,B) ( sqrt( pow((A)[0]-(B)[0],2) + pow((A)[1]-(B)[1],2) + pow((A)[2]-(B)[2],2) ) )
#define DIST2(A,B) ( pow((A)[0]-(B)[0],2) + pow((A)[1]-(B)[1],2) + pow((A)[2]-(B)[2],2) )

#define DIST_2D(A,B) ( sqrt( pow((A)[0]-(B)[0],2) + pow((A)[1]-(B)[1],2) ) )
#define DIST2_2D(A,B) ( pow((A)[0]-(B)[0],2) + pow((A)[1]-(B)[1],2) )

// returns 1 or -1, or 0 depending on the sign of X...
#define SGN(X) ((X)<0.0?-1.0:((X)>0.0?1.0:0.0))

#define ERRSTART "\n\n\n=======================================================================\n"
#define ERREND "\n=======================================================================\n\n\n"
#define ERREND_NO_RETURN "=======================================================================\n\n\n"
#define ERR(MESSAGE) ERRSTART << MESSAGE << ERREND << std::endl;


#define COUT1(MESSAGE) { if (mpi_rank==0) std::cout << MESSAGE << std::endl; }
#define COUT2(MESSAGE) { if ((mpi_rank==0)&&(cti_verbose)) std::cout << MESSAGE << std::endl; }
#define COUT3(MESSAGE) { if ((mpi_rank==0)&&(cti_verbose)) std::cout << MESSAGE << std::endl; }
#define CWARN(MESSAGE) { if (mpi_rank==0) std::cout <<			\
                          "\n-----------------------------------------------------------------------\nWARNING: " << \
			  MESSAGE << \
			  "\n-----------------------------------------------------------------------\n" << std::endl; }
// use the CERR macro as a collective/synchronized way to die. e.g.
// CERR("this is a problem: " << value);
#define CERR(MESSAGE) { if (mpi_rank==0) std::cerr <<			\
			  "\n\n=======================================================================\nERROR: " << \
			  MESSAGE <<					\
			  "\nFile: " << __FILE__ << ", line: " << __LINE__ << \
			  "\n=======================================================================\n" << std::endl; throw(0); }

#define CERRT(MESSAGE,THROW) { if (mpi_rank==0) std::cerr <<			\
			  "\n\n=======================================================================\nERROR: " << \
			  MESSAGE <<					\
			  "\nFile: " << __FILE__ << ", line: " << __LINE__ << \
			  "\n=======================================================================\n" << std::endl; throw(THROW); }


#define CERR_START(MESSAGE) { if (mpi_rank==0) std::cerr <<		\
				"\n\n=======================================================================\nERROR: " << \
				MESSAGE << std::endl; }

#define CERR_MIDDLE(MESSAGE) { if (mpi_rank==0) std::cerr << MESSAGE << std::endl; }

#define CERR_FINISH(MESSAGE) { if (mpi_rank==0) std::cerr <<		\
				 MESSAGE <<				\
				 "\nFile: " << __FILE__ << ", line: " << __LINE__ << \
				 "\n=======================================================================\n" << std::endl; throw(0); }

// serial version
#define COUT1_S(MESSAGE) { std::cout << MESSAGE << std::endl; }
#define COUT2_S(MESSAGE) { if (cti_verbose) std::cout << MESSAGE << std::endl; }
#define COUT3_S(MESSAGE) { if (cti_verbose) std::cout << MESSAGE << std::endl; }
#define CWARN_S(MESSAGE) { std::cout <<			\
                          "\n-----------------------------------------------------------------------\nWARNING: " << \
			  MESSAGE << \
			  "\n-----------------------------------------------------------------------" << std::endl; }
#define CERR_S(MESSAGE) { std::cerr <<			\
			  "\n\n=======================================================================\nERROR: " << \
			  MESSAGE << \
			  "\nFile: " << __FILE__ << ", line: " << __LINE__ << \
			  "\n=======================================================================\n" << std::endl; throw(0); }

// the following delete checks for null, deletes if it is null, and sets to NULL...
#define DELETE(ARRAY) { if ((ARRAY) != NULL) { delete[] (ARRAY); (ARRAY) = NULL; } }

// similar to the above for MPI_Data_types...
#define MPI_TYPE_FREE(T) { if ((T) != MPI_DATATYPE_NULL) { MPI_Type_free(&(T)); (T) = MPI_DATATYPE_NULL; } }

#define SWAP3(A,B) {					\
const double tmp[3] = { (A)[0], (A)[1], (A)[2] };	\
(A)[0] = (B)[0];					\
(A)[1] = (B)[1];					\
(A)[2] = (B)[2];					\
(B)[0] = tmp[0];					\
(B)[1] = tmp[1];					\
(B)[2] = tmp[2];					\
}

// for compiler suppression of unused parameters/variables
#define UNUSED(expr) (void)(expr)

// for quickly getting an int as a string...
#define SSTR( x ) static_cast< std::ostringstream & >( ( std::ostringstream() << std::dec << x ) ).str()

#endif
