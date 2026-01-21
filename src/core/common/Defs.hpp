#ifndef CTIDEFS_HPP
#define CTIDEFS_HPP

#ifdef INT8_IS_LONG_INT
typedef long int int8;
typedef unsigned long int uint8;
#else
typedef long long int int8;
typedef unsigned long long int uint8;
#endif

typedef unsigned char uint1;
typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned short uint2;

//TODO: get rid of this...
typedef double cti_real;

typedef double double3[3];

// ######################################################################
// io related definitions
// ######################################################################

#define UGP_IO_MAGIC_NUMBER           123581321
#define UGP_IO_VERSION                4

#define CHEM_IO_MAGIC_NUMBER          123581321
#define CHEM_IO_VERSION               2

#define POINTS_IO_VERSION 2
#define POINTS_IO_MAGIC_NUMBER 1235813
#define POINTS_IO_STAMP 2478

#define PART_IO_MAGIC_NUMBER 1235814
#define PART_IO_VERSION 1

#define MOVE_IO_MAGIC_NUMBER 1235815
#define MOVE_IO_VERSION 1 // moving solver io switched to UGP_IO_VERSION after version 1

#define TETS_IO_MAGIC_NUMBER          1235816
#define TETS_IO_VERSION               1

// restart record indexing - do not change these...
#define UGP_IO_NO_FA_CV_COUNTS        10
#define UGP_IO_NO_FA_BF_CV_COUNTS     10001
#define UGP_IO_I0                     11
#define UGP_IO_D0                     12

#define UGP_IO_FA_CHECK               20
#define UGP_IO_FA_CHECK_INT8          2001
#define UGP_IO_NOOFA_I_AND_V          21
#define UGP_IO_NOOFA_I_AND_V_V3_INT8  2002
#define UGP_IO_NOOBF_I_AND_V_INT8     2005
#define UGP_IO_NOOFA_I_AND_V_INT8     2007
#define UGP_IO_EDOFA_V                28
#define UGP_IO_CVOFA                  22
#define UGP_IO_CVOFA_V3_INT8          2003
#define UGP_IO_CVOBF_INT8             2004
#define UGP_IO_CVOFA_INT8             2006
#define UGP_IO_FA_ZONE_HEADER         23
#define UGP_IO_FA_D1                  24
#define UGP_IO_FA_D2                  25
#define UGP_IO_FA_I1                  26
#define UGP_IO_FA_GLOBAL_INDEX        29
#define UGP_IO_FA_ZONE                27
#define UGP_IO_BF_ZONE_HEADER         777
//#define UGP_IO_STOBF_I_AND_V_INT      2008
#define UGP_IO_SBOBF_I_AND_V          2009
#define UGP_IO_SPOBF_I_V_WGT          2010

// extended Voronoi face geometry
#define UGP_IO_CVOEF                  2110
#define UGP_IO_EF_D1                  2111
#define UGP_IO_EF_D2                  2112
#define UGP_IO_EF_I1                  2113

// special VV bbox record for fast point
// location: use with x_vv and r_vv
#define UGP_IO_VV_BBOX                2115

// boundary face geometry
#define UGP_IO_BF_D1                  2120
#define UGP_IO_BF_D2                  2121
#define UGP_IO_BF_I1                  2123
// note the transition to CIFH naming convention associated
// with new CtiRegister. Should push this through everywhere
// eventually...
#define UGP_IO_BF_DN33                2122

#define UGP_IO_NO_CHECK               30
#define UGP_IO_NO_CHECK_INT8          3001
#define UGP_IO_X_NO                   31
#define UGP_IO_NO_D1                  32
#define UGP_IO_NO_D2                  33
#define UGP_IO_NO_I1                  34
#define UGP_IO_NO_GLOBAL_INDEX        35
#define UGP_IO_PBI_PNO                36 // note this is uint8 and only for periodic

#define UGP_IO_CV_CHECK               40
#define UGP_IO_CV_CHECK_INT8          4001
#define UGP_IO_CV_PART                41
#define UGP_IO_CV_D1                  42
#define UGP_IO_CV_D2                  43
#define UGP_IO_CV_I1                  45
#define UGP_IO_CV_ZONE_HEADER         46
#define UGP_IO_CV_ZONE                47
#define UGP_IO_CV_GLOBAL_INDEX        48

// these new counts use the uint8 part of the header to store their size, rather
// than Msb/Lsb 2-int decomposition used originally...
#define UGP_IO_CV_D1_NEW              4042
#define UGP_IO_CV_D2_NEW              4043
#define UGP_IO_CV_I1_NEW              4045

#define UGP_IO_CV_F1                  420
#define UGP_IO_CV_F2                  421

// edges ...
#define UGP_IO_ED_CHECK               100
#define UGP_IO_ED_I1                  101
#define UGP_IO_ED_D1                  102
#define UGP_IO_ED_D2                  103
#define UGP_IO_NOOED                  104
#define UGP_IO_ED_GLOBAL_INDEX        105

//fazone data
#define UGP_IO_FAZONE_NO_I1           110
#define UGP_IO_FAZONE_FA_I1           111
#define UGP_IO_FAZONE_NO_D1           112
#define UGP_IO_FAZONE_FA_D1           113
#define UGP_IO_FAZONE_NO_D2           114
#define UGP_IO_FAZONE_FA_D2           115
#define UGP_IO_BF_CHECK_INT8          1115

#define UGP_IO_DATA                   50
#define UGP_IO_EOF                    51

// for particle i/o...
#define UGP_IO_LPOCV_COUNTS           60
#define UGP_IO_LP_CHECK               61
#define UGP_IO_LP_XP                  62
#define UGP_IO_LP_I0                  63
#define UGP_IO_LP_I1                  64
#define UGP_IO_LP_D0                  65
#define UGP_IO_LP_D1                  66
#define UGP_IO_LP_D2                  67
#define UGP_IO_LP_EOF                 68
#define UGP_IO_LPOCV_I                69

// for ReadWrite class data...
#define UGP_IO_READWRITE              80

// for History/Pedigree info...
#define UGP_IO_TIMESTAMP              90
#define UGP_IO_PARAMS                 91
#define UGP_IO_JOURNAL                92
#define UGP_IO_WIREFRAME              93
#define UGP_IO_HASHIDS                94
#define UGP_IO_LICENSE                95

#define UGP_IO_BIT_TRANSFORM_DATA     4122
#define UGP_IO_PERIODIC_TRANSFORM     4123

// chemtable io
#define CART_CHEM_TABLE_IO_COUNTS     1024
#define CART_CHEM_TABLE_CRD_DATA      2048
#define CART_CHEM_TABLE_TAB_DATA      4096
#define CART_CHEM_TABLE_IO_EOF        8192

#define UGP_IO_CT_COOR                70
#define UGP_IO_CT_DATA                71
#define UGP_IO_CT_KD_HC               72
#define UGP_IO_CT_KD_NODE             73
#define UGP_IO_CT_KD                  74
#define UGP_IO_CT_CART_1D             75
#define UGP_IO_CT_CART_2D             76
#define UGP_IO_CT_CART_3D             77
#define UGP_IO_CT_CART_4D             78

// surface record...
#define UGP_IO_SURFACE               801
#define UGP_IO_SURFACE_PERIODIC_INFO 802

// moving solver defs use 7XXX for now...
#define UGP_IO_MOVING_HEADER        7000

// FlaggedFaceWriter data layouts - now used by FwhSurface class...
#define FFW_DATA_I1                    1
#define FFW_DATA_D1                    2
#define FFW_DATA_D2                    3
#define FFW_DATA_D1D2                  4
#define FFW_DATA_D2D2                  5
#define FFW_DATA_UNKNOWN              -1

// CTI files are made up of multiple records consisting of the following
// Header, then the data, Header, data, etc...
#define UGP_IO_HEADER_NAME_LEN        52 // select this to make header size 256 bytes

// endianness...
#define CTI_BIG_ENDIAN 0
#define CTI_LITTLE_ENDIAN 1

// various datatypes...
//#define VALUE_DATA     0
//#define CV_DATA        1
//#define FA_DATA        2
//#define NO_DATA        3
//#define LP_DATA        4
//#define SIGNED_FA_DATA 5 // normal-aligned data is flipped when required
//#define CV_G_DATA      6 // data that gets reallocated to ghosts when ghost cvs are added...
//#define ZF_B_DATA      7 // zone-face boundary data (all boundary zones)
//#define ED_DATA        8
//#define ZN_B_DATA      9 // zone-node boundary data (all boundary zones)
//#define ZN_DATA        10
//#define ZF_DATA        11

// additionally, operator data can be registered, but not I/O'd yet...
// for now, fazone's only...
//#define ZNOZN_DATA     12

// bits for io...
//#define NOREAD_DATA   16
//#define NOWRITE_DATA  32

class Header {
public:
  char name[UGP_IO_HEADER_NAME_LEN];
  int id;
  int8 skip;
  int idata[16];
  // note: in march 2019 rdata was reduced from 16 to 12 to
  // make space for 4 uint8's. This should not affect any header
  // reading of rdata because the maximum use of rdata to date
  // is 12 (associated with periodic transform matrix and vector)
  double rdata[12];
  uint8 ui8data[4];
  Header() {
    for (int i = 0; i < UGP_IO_HEADER_NAME_LEN; ++i)
      name[i] = 0;
    id = 0;
    skip = 0;
    for (int i = 0; i < 16; ++i) idata[i] = 0;
    for (int i = 0; i < 12; ++i) rdata[i] = 0.0;
    for (int i = 0; i < 4; ++i) ui8data[i] = 0;
  }
};

// used in MiscUtils::calcThresholdDist
// double scalar would be a 1MB buffer
#define DIST_THRESHOLD 131072

// these defines are used as action arguments for halo data updates...
enum UpdateAction {
  UNDEFINED_ACTION,
  REPLACE_DATA,
  ADD_DATA,
  MIN_DATA,
  MAX_DATA,
  REPLACE_TRANSLATE_DATA,
  TRANSLATE_REPLACE_DATA,
  REPLACE_ROTATE_DATA,
  ROTATE_REPLACE_DATA,
  ADD_ROTATE_DATA,
  BITWISE_OR_DATA,
  BITWISE_OR_NO_PERIODIC_DATA,
  MIN_NO_PERIODIC_DATA,
  MAX_NO_PERIODIC_DATA,
  ADD_NO_PERIODIC_DATA,
  ADD_TRANSLATE_DATA,
  SUBTRACT_DATA,
  SUBTRACT_ROTATE_DATA,
  SUBTRACT_TRANSLATE_DATA,
  NO_CHANGE_DATA,
  MIN_TRANSLATE_DATA,
  MAX_TRANSLATE_DATA,
  ADD_GHOST_TO_ACTIVE_DATA_,
  REPLACE_GHOST_DATA_,
  MIN_ABS_DATA
};

// face zone kind's...
// reserve -1 for unknown
#define FA_ZONE_UNKNOWN_UNNAMED  -2
#define FA_ZONE_UNKNOWN          -1
//#define FA_ZONE_PERIODIC_UNKNOWN -2 // used by MshFilter at some point?

#define FA_ZONE_BOUNDARY          1
#define FA_ZONE_PERIODIC_CART     2
#define FA_ZONE_PERIODIC_CYL_X    3
#define FA_ZONE_PERIODIC_CYL_Y    4
#define FA_ZONE_PERIODIC_CYL_Z    5
#define FA_ZONE_PERIODIC_INTERNAL 6 // Warning: this was internal in io_version == 1 files
#define FA_ZONE_INTERNAL         20

// range for all and periodic face zones...
#define FA_ZONE_PERIODIC_FIRST    2
#define FA_ZONE_PERIODIC_LAST     6

// cv zone kinds...
#define CV_ZONE_UNKNOWN_UNNAMED  -2
#define CV_ZONE_UNKNOWN          -1
#define CV_ZONE_FLUID             1

// mpi tags used for data exchanges: these have to be
// ints -- they are passed as tags in MPI_Send/Recv...
#define UPDATE_R1_TAG                   1002
#define UPDATE_R2_TAG                   1003
#define UPDATE_R2_REVERSE_TAG           1004
#define EXCHANGE_INT_TAG                1005
#define EXCHANGE_INT_REVERSE_TAG        1006
#define UPDATE_I1_TAG                   1007
#define UPDATE_I1_REVERSE_TAG           1008
#define UPDATE_I2_TAG                   1009
#define UPDATE_I2_REVERSE_TAG           1010
#define UPDATE_R3_TAG                   1011
#define UPDATE_SYMMETRIC_R3_TAG         1012
#define EXCHANGE_PRCOMM_INTV_TAG1       1013
#define EXCHANGE_PRCOMM_INTV_TAG2       1014
#define EXCHANGE_PRCOMM_DOUBLEV_TAG1    1015
#define EXCHANGE_PRCOMM_DOUBLEV_TAG2    1016
#define EXCHANGE_PRCOMM_NPACKV_TAG      1017
#define EXCHANGE_PRCOMM_NUNPACKV_TAG    1018
#define UPDATE_ADDRESS_TAG              1019
#define UPDATE_R2R3_TAG                 1020
#define UPDATE_R2R1R1R1R1R1_TAG         1021
#define UPDATE_GIL_TAG                  1022

#define EXCHANGE_PRCOMM_INT8V_TAG1      1030
#define EXCHANGE_PRCOMM_INT8V_TAG2      1031

#define UPDATE_D2_TAG                   1050

#define UPDATE_RN_TAG                   1100
#define UPDATE_RN3_TAG                  1101

// a resonable limit for signed int's...
#define TWO_BILLION 2000000000ll
#define ONE_BILLION 1000000000ll

// use BIG_INT and -BIG_INT to safely represent the "largest"
// signed 4-byte integer (2^31 = 2.147 billion)
#define BIG_INT     2000000000

#define CTI_TIMESTAMP_SIZE 24

// types of elements...
// limit this to 8 types max for use in adaptation
#define HEX_TYPE     0
#define PRISM_TYPE   1
#define PYRAMID_TYPE 2
#define OCT_TYPE     3
#define TET_TYPE     4
#define UNKNOWN_TYPE 5

// types of faces...
// limit this to 4 types max for use in adaptation
#define QUAD_TYPE    0
#define TRI_TYPE     1

#define RGB(r,g,b) (unsigned short)(r + (g << 5) + (b << 10))

// icvp will either contain a positive index corresponding to
// the cv it is likely in, or one of the following...
#define LP_RECYCLE    -1 // particle should be recycled
#define LP_EXCHANGE   -2 // particle should be exchanged, then recycled
#define LP_UNKNOWN    -3 // particle in unknown state
#define LP_IN_IB_CV   -4 // particle is located properly

//VTK
#define VTK_TET 10
#define VTK_HEX 12
#define VTK_PRISM 13
#define VTK_PYRAMID 14

//bit masking
#define MASK_60BITS 1152921504606846975ll
#define MASK_55BITS 36028797018963967ll
#define MASK_52BITS 4503599627370495ll
#define MASK_32BITS 4294967295ll
#define MASK_30BITS 1073741823ll
#define MASK_20BITS 1048575
#define MASK_12BITS 4095
#define MASK_6BITS 63
#define MASK_2BITS 3

// status used for diagnostics
#define CTI_STATUS_CHECK  0
#define CTI_STATUS_ACTIVE 1
#define CTI_STATUS_ERROR -1

enum SubzoneLimits {
  SURFACE_SZ_MIN = 0,
  SURFACE_SZ_MAX = 32768,
  OPEN_E_SZ_MIN = 32769,
  OPEN_E_SZ_MAX = 49151,
  DYN_E_SZ_MIN = 49152,
  DYN_E_SZ_MAX = 65536,
};

// stitch stuff
//#define MESH_INSPECTOR_IST 100100100 // something positive and big

// TODO: this is NOT the right place for this: probably in serial utils some day...

template <class T>
bool from_string(T& t, const std::string &s, std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  bool b_fail = (iss >> f >> t).fail();
  std::string remainder;
  iss >> remainder; //
  if (remainder.size() != 0) b_fail = true;
  return !b_fail;
}

#endif
