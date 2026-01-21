#include "SimpleSurface.hpp"

// TODO: put this somewhere else some day. It was taken from
// CI's parallel implementation in stitch...

// --------------------------------------------------------
// marching cube tables

static int mcEdgeTable[256] = {
  0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
  0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
  0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
  0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
  0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
  0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
  0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
  0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
  0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
  0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
  0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
  0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
  0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
  0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
  0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
  0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
  0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
  0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
  0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
  0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
  0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
  0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
  0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
  0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
  0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
  0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
  0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
  0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
  0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
  0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
  0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
  0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

static int mcTriTable[256][16] =
  {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
  {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
  {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
  {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
  {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
  {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
  {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
  {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
  {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
  {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
  {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
  {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
  {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
  {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
  {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
  {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
  {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
  {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
  {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
  {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
  {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
  {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
  {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
  {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
  {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
  {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
  {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
  {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
  {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
  {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
  {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
  {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
  {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
  {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
  {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
  {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
  {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
  {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
  {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
  {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
  {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
  {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
  {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
  {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
  {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
  {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
  {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
  {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
  {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
  {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
  {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
  {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
  {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
  {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
  {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
  {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
  {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
  {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
  {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
  {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
  {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
  {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
  {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
  {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
  {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
  {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
  {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
  {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
  {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
  {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
  {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
  {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
  {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
  {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
  {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
  {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
  {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
  {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
  {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
  {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
  {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
  {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
  {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
  {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
  {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
  {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
  {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
  {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
  {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
  {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
  {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
  {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
  {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
  {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
  {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
  {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
  {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
  {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
  {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
  {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
  {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
  {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
  {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
  {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
  {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
  {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
  {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
  {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
  {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
  {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
  {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
  {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
  {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
  {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
  {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
  {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
  {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
  {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
  {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
  {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
  {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
  {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
  {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
  {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
  {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
  {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
  {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
  {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
  {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
  {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
  {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
  {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
  {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
  {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
  {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
  {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
  {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
  {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
  {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
  {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
  {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
  {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
  {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
  {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
  {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
  {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
  {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
  {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
  {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
  {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
  {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
  {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
  {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
  {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
  {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
  {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
  {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
  {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
  {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
  {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
  {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
  {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
  {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
  {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
  {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
  {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
  {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
  {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
  {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

// --------------------------------------------------------

int SimpleSurface::addLiftedSurface(const double dn) {

  cout << "SimpleSurface::addLiftedSurface: dn: " << dn << "..." << endl;

  // assume the surfaces to lift have been selected in sz_flag...
  st_flag.setLength(nst);
  int nst_s = 0; // nst selected
  for (int ist = 0; ist < nst; ++ist) {
    if (sz_flag[szost[ist]]) {
      st_flag[ist] = nst_s++;
    }
    else {
      st_flag[ist] = -1;
    }
  }

  cout << " > lifting surface around " << nst_s << " out of " << nst << " tris" << endl;

  // if no surfaces to lift, we return -1...
  if (nst_s == 0)
    return -1;

  // reverse lookup: st-of-st-selected...
  int *stosts = new int[nst_s];
  nst_s = 0;
  for (int ist = 0; ist < nst; ++ist) {
    if (sz_flag[szost[ist]]) {
      stosts[nst_s] = ist;
      ++nst_s;
    }
  }

  // put tris in an adt...

  double (*bbmin)[3] = new double[nst_s][3];
  double (*bbmax)[3] = new double[nst_s][3];
  double xstart[3] = { HUGE_VAL, 0.0, 0.0 };
  for (int ist_s = 0; ist_s < nst_s; ++ist_s) {
    const int ist = stosts[ist_s];
    const double * const x0 = xsp[spost[ist][0]];
    const double * const x1 = xsp[spost[ist][1]];
    const double * const x2 = xsp[spost[ist][2]];
    FOR_I3 bbmin[ist_s][i] = min(x0[i],min(x1[i],x2[i]));
    FOR_I3 bbmax[ist_s][i] = max(x0[i],max(x1[i],x2[i]));
    // find the x0-most point...
    if (x0[0] < xstart[0]) FOR_I3 xstart[i] = x0[i];
    if (x1[0] < xstart[0]) FOR_I3 xstart[i] = x1[i];
    if (x2[0] < xstart[0]) FOR_I3 xstart[i] = x2[i];
  }

  Adt<double> adt(nst_s,bbmin,bbmax);
  delete[] bbmin;
  delete[] bbmax;

  // the grid size...
  // can we just use a common fraction of dn?

  const double dx = dn/2.0;

  // locate x0 at the w/s/b corner. This keeps all
  // the indexing positive. And find a starting point...

  double x0[3],x1[3];
  adt.getBbox(x0,x1);
  FOR_I3 x0[i] -= dn*1.1;
  FOR_I3 x1[i] += dn*1.1;
  xstart[0] -= dn;

  // check how big the grid is going to be...
  FOR_I3 {
    const int n = int((x1[i]-x0[i]+2.2*dn)/dx);
    assert(n < 1024); // if you hit this, we need to rethink dx = dn/4 above
  }

  // decide on the ijk...
  int i = int((xstart[0]-x0[0])/dx); assert(i >= 0);
  int j = int((xstart[1]-x0[1])/dx); assert(j >= 0);
  int k = int((xstart[2]-x0[2])/dx); assert(k >= 0);

  vector<int> triVec;
  double x[3];
  x[0] = x0[0] + dx*i;
  x[1] = x0[1] + dx*j;
  x[2] = x0[2] + dx*k;
  //cout << "X1: " << COUT_VEC(x) << endl;
  assert(triVec.empty());
  adt.buildListForClosestPoint(triVec,x);
  double d2 = HUGE_VAL;
  for (vector<int>::const_iterator it=triVec.begin(); it!=triVec.end(); ++it) {
    const int ist = stosts[*it];
    d2 = min(d2,MiscUtils::getPointToTriDist2(x,xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]));
  }
  triVec.clear();

  if (sqrt(d2) < dn) {

    while (sqrt(d2) < dn) {

      // go further away in the -i direction...

      --i; assert(i >= 0);
      x[0] = x0[0] + dx*i;
      x[1] = x0[1] + dx*j;
      x[2] = x0[2] + dx*k;
      assert(triVec.empty());
      adt.buildListForClosestPoint(triVec,x);
      d2 = HUGE_VAL;
      for (vector<int>::const_iterator it=triVec.begin(); it!=triVec.end(); ++it) {
        const int ist = stosts[*it];
        d2 = min(d2,MiscUtils::getPointToTriDist2(x,xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]));
      }
      triVec.clear();

    }

  }
  else {

    assert(sqrt(d2) >= dn);
    while (sqrt(d2) >= dn) {

      // come closer...

      ++i;
      x[0] = x0[0] + dx*i;
      x[1] = x0[1] + dx*j;
      x[2] = x0[2] + dx*k;
      assert(triVec.empty());
      adt.buildListForClosestPoint(triVec,x);
      d2 = HUGE_VAL;
      for (vector<int>::const_iterator it=triVec.begin(); it!=triVec.end(); ++it) {
        const int ist = stosts[*it];
        d2 = min(d2,MiscUtils::getPointToTriDist2(x,xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]));
      }
      triVec.clear();

      //cout << "++i d2: " << sqrt(d2) << " dn: " << dn << endl;

      // put i back because cube is staggered: (cube i has nodes i->i+1)
      --i;

    }

  }

  // the cv we want to push onto the stack is (i,j,k)...
  stack<int> cvStack;
  cvStack.push((k<<20)|(j<<10)|i);
  set<int> cvSet;
  cvSet.insert((k<<20)|(j<<10)|i);
  double cube[2][2][2];
  map<const int,double> sdMap; // signed distance map
  int ned = 0;
  map<const pair<int,int>,int> edMap;
  vector<int> spostVec;
  int vertlist[12];
  vector<double> xspVec;

  while (!cvStack.empty()) {

    const int ijk = cvStack.top(); cvStack.pop();
    const int i = ijk & 1023;
    const int j = (ijk>>10) & 1023;
    const int k = (ijk>>20);

    //cout << " just popped: " << i << " " << j << " " << k << endl;

    assert(i >= 0);
    assert(j >= 0);
    assert(k >= 0);

    // put a value in each corner of the cube...

    for (int di = 0; di <= 1; ++di) {
      for (int dj = 0; dj <= 1; ++dj) {
        for (int dk = 0; dk <= 1; ++dk) {
          const int this_ijk = (((k+dk)<<20)|((j+dj)<<10)|(i+di));
          map<const int,double>::iterator iter = sdMap.find(this_ijk);
          if (iter == sdMap.end()) {
            x[0] = x0[0]+dx*(i+di);
            x[1] = x0[1]+dx*(j+dj);
            x[2] = x0[2]+dx*(k+dk);
            assert(triVec.empty());
            adt.buildListForClosestPoint(triVec,x);
            double d2 = HUGE_VAL;
            for (vector<int>::const_iterator it=triVec.begin(); it!=triVec.end(); ++it) {
              const int ist = stosts[*it];
              d2 = min(d2,MiscUtils::getPointToTriDist2(x,xsp[spost[ist][0]],xsp[spost[ist][1]],xsp[spost[ist][2]]));
            }
            triVec.clear();
            sdMap[this_ijk] = cube[di][dj][dk] = sqrt(d2) - dn;
          }
          else {
            cube[di][dj][dk] = iter->second;
          }
        }
      }
    }

    // at this point, the 8 cube vertex values are set, so apply the marching cubes algo to build triangles...

    int cube_index = 0;
    if (cube[0][0][0] < 0.0) cube_index |= 1;
    if (cube[1][0][0] < 0.0) cube_index |= 2;
    if (cube[1][1][0] < 0.0) cube_index |= 4;
    if (cube[0][1][0] < 0.0) cube_index |= 8;
    if (cube[0][0][1] < 0.0) cube_index |= 16;
    if (cube[1][0][1] < 0.0) cube_index |= 32;
    if (cube[1][1][1] < 0.0) cube_index |= 64;
    if (cube[0][1][1] < 0.0) cube_index |= 128;

    // cube is entirely in/out of the surface
    // this should not be the case...

    assert(mcEdgeTable[cube_index] != 0);

    // find the vertices where the surface intersects the cube

    for (int i = 0; i < 12; ++i) vertlist[i] = -1;

    if (mcEdgeTable[cube_index] & 1) {
      const pair<int,int> ijk_pair = pair<int,int>(((k)<<20)|((j)<<10)|(i),((k)<<20)|((j)<<10)|(i+1));
      map<const pair<int,int>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[0] = ned++;
        xspVec.push_back((cube[1][0][0]*(x0[0]+dx*(i))-cube[0][0][0]*(x0[0]+dx*(i+1)))/(cube[1][0][0]-cube[0][0][0]));
        xspVec.push_back(x0[1]+dx*(j));
        xspVec.push_back(x0[2]+dx*(k));
      }
      else {
        vertlist[0] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 2) {
      const pair<int,int> ijk_pair = pair<int,int>(((k)<<20)|((j)<<10)|(i+1),((k)<<20)|((j+1)<<10)|(i+1));
      map<const pair<int,int>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[1] = ned++;
        xspVec.push_back(x0[0]+dx*(i+1));
        xspVec.push_back((cube[1][1][0]*(x0[1]+dx*(j))-cube[1][0][0]*(x0[1]+dx*(j+1)))/(cube[1][1][0]-cube[1][0][0]));
        xspVec.push_back(x0[2]+dx*(k));
      }
      else {
        vertlist[1] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 4) {
      const pair<int,int> ijk_pair = pair<int,int>(((k)<<20)|((j+1)<<10)|(i),((k)<<20)|((j+1)<<10)|(i+1));
      map<const pair<int,int>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[2] = ned++;
        xspVec.push_back((cube[1][1][0]*(x0[0]+dx*(i))-cube[0][1][0]*(x0[0]+dx*(i+1)))/(cube[1][1][0]-cube[0][1][0]));
        xspVec.push_back(x0[1]+dx*(j+1));
        xspVec.push_back(x0[2]+dx*(k));
      }
      else {
        vertlist[2] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 8) {
      const pair<int,int> ijk_pair = pair<int,int>(((k)<<20)|((j)<<10)|(i),((k)<<20)|((j+1)<<10)|(i));
      map<const pair<int,int>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[3] = ned++;
        xspVec.push_back(x0[0]+dx*(i));
        xspVec.push_back((cube[0][1][0]*(x0[1]+dx*(j))-cube[0][0][0]*(x0[1]+dx*(j+1)))/(cube[0][1][0]-cube[0][0][0]));
        xspVec.push_back(x0[2]+dx*(k));
      }
      else {
        vertlist[3] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 16) {
      const pair<int,int> ijk_pair = pair<int,int>(((k+1)<<20)|((j)<<10)|(i),((k+1)<<20)|((j)<<10)|(i+1));
      map<const pair<int,int>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[4] = ned++;
        xspVec.push_back((cube[1][0][1]*(x0[0]+dx*(i))-cube[0][0][1]*(x0[0]+dx*(i+1)))/(cube[1][0][1]-cube[0][0][1]));
        xspVec.push_back(x0[1]+dx*(j));
        xspVec.push_back(x0[2]+dx*(k+1));
      }
      else {
        vertlist[4] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 32) {
      const pair<int,int> ijk_pair = pair<int,int>(((k+1)<<20)|((j)<<10)|(i+1),((k+1)<<20)|((j+1)<<10)|(i+1));
      map<const pair<int,int>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[5] = ned++;
        xspVec.push_back(x0[0]+dx*(i+1));
        xspVec.push_back((cube[1][1][1]*(x0[1]+dx*(j))-cube[1][0][1]*(x0[1]+dx*(j+1)))/(cube[1][1][1]-cube[1][0][1]));
        xspVec.push_back(x0[2]+dx*(k+1));
      }
      else {
        vertlist[5] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 64) {
      const pair<int,int> ijk_pair = pair<int,int>(((k+1)<<20)|((j+1)<<10)|(i),((k+1)<<20)|((j+1)<<10)|(i+1));
      map<const pair<int,int>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[6] = ned++;
        xspVec.push_back((cube[1][1][1]*(x0[0]+dx*(i))-cube[0][1][1]*(x0[0]+dx*(i+1)))/(cube[1][1][1]-cube[0][1][1]));
        xspVec.push_back(x0[1]+dx*(j+1));
        xspVec.push_back(x0[2]+dx*(k+1));
      }
      else {
        vertlist[6] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 128) {
      const pair<int,int> ijk_pair = pair<int,int>(((k+1)<<20)|((j)<<10)|(i),((k+1)<<20)|((j+1)<<10)|(i));
      map<const pair<int,int>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[7] = ned++;
        xspVec.push_back(x0[0]+dx*(i));
        xspVec.push_back((cube[0][1][1]*(x0[1]+dx*(j))-cube[0][0][1]*(x0[1]+dx*(j+1)))/(cube[0][1][1]-cube[0][0][1]));
        xspVec.push_back(x0[2]+dx*(k+1));
      }
      else {
        vertlist[7] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 256) {
      const pair<int,int> ijk_pair = pair<int,int>(((k)<<20)|((j)<<10)|(i),((k+1)<<20)|((j)<<10)|(i));
      map<const pair<int,int>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[8] = ned++;
        xspVec.push_back(x0[0]+dx*(i));
        xspVec.push_back(x0[1]+dx*(j));
        xspVec.push_back((cube[0][0][1]*(x0[2]+dx*(k))-cube[0][0][0]*(x0[2]+dx*(k+1)))/(cube[0][0][1]-cube[0][0][0]));
      }
      else {
        vertlist[8] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 512) {
      const pair<int,int> ijk_pair = pair<int,int>(((k)<<20)|((j)<<10)|(i+1),((k+1)<<20)|((j)<<10)|(i+1));
      map<const pair<int,int>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[9] = ned++;
        xspVec.push_back(x0[0]+dx*(i+1));
        xspVec.push_back(x0[1]+dx*(j));
        xspVec.push_back((cube[1][0][1]*(x0[2]+dx*(k))-cube[1][0][0]*(x0[2]+dx*(k+1)))/(cube[1][0][1]-cube[1][0][0]));
      }
      else {
        vertlist[9] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 1024) {
      const pair<int,int> ijk_pair = pair<int,int>(((k)<<20)|((j+1)<<10)|(i+1),((k+1)<<20)|((j+1)<<10)|(i+1));
      map<const pair<int,int>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[10] = ned++;
        xspVec.push_back(x0[0]+dx*(i+1));
        xspVec.push_back(x0[1]+dx*(j+1));
        xspVec.push_back((cube[1][1][1]*(x0[2]+dx*(k))-cube[1][1][0]*(x0[2]+dx*(k+1)))/(cube[1][1][1]-cube[1][1][0]));
      }
      else {
        vertlist[10] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 2048) {
      const pair<int,int> ijk_pair = pair<int,int>(((k)<<20)|((j+1)<<10)|(i),((k+1)<<20)|((j+1)<<10)|(i));
      map<const pair<int,int>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[11] = ned++;
        xspVec.push_back(x0[0]+dx*(i));
        xspVec.push_back(x0[1]+dx*(j+1));
        xspVec.push_back((cube[0][1][1]*(x0[2]+dx*(k))-cube[0][1][0]*(x0[2]+dx*(k+1)))/(cube[0][1][1]-cube[0][1][0]));
      }
      else {
        vertlist[11] = iter->second;
      }
    }

    // add triangles...

    for (int l = 0; mcTriTable[cube_index][l] != -1; l += 3) {
      spostVec.push_back(vertlist[mcTriTable[cube_index][l  ]]);
      spostVec.push_back(vertlist[mcTriTable[cube_index][l+1]]);
      spostVec.push_back(vertlist[mcTriTable[cube_index][l+2]]);
    }

    // now travel across any face that has both positive and negative
    // cvs...

    int count = 0;
    for (int dj = 0; dj <= 1; ++dj) {
      for (int dk = 0; dk <= 1; ++dk) {
        if (cube[0][dj][dk] >= 0.0)
          ++count;
      }
    }
    //cout << "x0: " << count << endl;
    if ((count >= 1)&&(count <= 3)) {
      const int this_ijk = ((k)<<20)|((j)<<10)|(i-1);
      if (cvSet.find(this_ijk) == cvSet.end()) {
        cvSet.insert(this_ijk);
        cvStack.push(this_ijk);
      }
    }

    count = 0;
    for (int dj = 0; dj <= 1; ++dj) {
      for (int dk = 0; dk <= 1; ++dk) {
        if (cube[1][dj][dk] >= 0.0)
          ++count;
      }
    }
    //cout << "x1: " << count << endl;
    if ((count >= 1)&&(count <= 3)) {
      const int this_ijk = ((k)<<20)|((j)<<10)|(i+1);
      if (cvSet.find(this_ijk) == cvSet.end()) {
        cvSet.insert(this_ijk);
        cvStack.push(this_ijk);
      }
    }

    count = 0;
    for (int di = 0; di <= 1; ++di) {
      for (int dk = 0; dk <= 1; ++dk) {
        if (cube[di][0][dk] >= 0.0)
          ++count;
      }
    }
    //cout << "y0: " << count << endl;
    if ((count >= 1)&&(count <= 3)) {
      const int this_ijk = ((k)<<20)|((j-1)<<10)|(i);
      if (cvSet.find(this_ijk) == cvSet.end()) {
        cvSet.insert(this_ijk);
        cvStack.push(this_ijk);
      }
    }

    count = 0;
    for (int di = 0; di <= 1; ++di) {
      for (int dk = 0; dk <= 1; ++dk) {
        if (cube[di][1][dk] >= 0.0)
          ++count;
      }
    }
    //cout << "y1: " << count << endl;
    if ((count >= 1)&&(count <= 3)) {
      const int this_ijk = ((k)<<20)|((j+1)<<10)|(i);
      if (cvSet.find(this_ijk) == cvSet.end()) {
        cvSet.insert(this_ijk);
        cvStack.push(this_ijk);
      }
    }

    count = 0;
    for (int di = 0; di <= 1; ++di) {
      for (int dj = 0; dj <= 1; ++dj) {
        if (cube[di][dj][0] >= 0.0)
          ++count;
      }
    }
    //cout << "z0: " << count << endl;
    if ((count >= 1)&&(count <= 3)) {
      const int this_ijk = ((k-1)<<20)|((j)<<10)|(i);
      if (cvSet.find(this_ijk) == cvSet.end()) {
        cvSet.insert(this_ijk);
        cvStack.push(this_ijk);
      }
    }

    count = 0;
    for (int di = 0; di <= 1; ++di) {
      for (int dj = 0; dj <= 1; ++dj) {
        if (cube[di][dj][1] >= 0.0)
          ++count;
      }
    }
    //cout << "z1: " << count << endl;
    if ((count >= 1)&&(count <= 3)) {
      const int this_ijk = ((k+1)<<20)|((j)<<10)|(i);
      if (cvSet.find(this_ijk) == cvSet.end()) {
        cvSet.insert(this_ijk);
        cvStack.push(this_ijk);
      }
    }

  }

  delete[] stosts;

  cout << " > out of stack loop. ned: " << ned << " spostVec.size(): " << spostVec.size() << endl;

  // use the edMap to build a bunch of new nodes...

  edMap.clear();
  assert(int(xspVec.size()) == ned*3);

  // add the new points...
  int nsp_orig = nsp; nsp += ned;
  growNspData(nsp,nsp_orig);

  for (int ied = 0; ied < ned; ++ied) {
    FOR_I3 xsp[nsp_orig+ied][i] = xspVec[ied*3+i];
  }
  xspVec.clear();

  // now add the lifted surface to zone "LIFTED"...
  string zone_name = "LIFTED";
  int izone,nzn;
  for (izone = 0,nzn=zoneVec.size(); izone < nzn; ++izone) {
    if (zoneVec[izone].getName() == zone_name)
      break;
  }
  if (izone == int(zoneVec.size())) {
    zoneVec.push_back(SurfaceZone(zone_name));
  }

  // and add the new tris...
  assert(spostVec.size()%3 == 0);
  const int nst_inc = spostVec.size()/3;
  int nst_orig = nst; nst += nst_inc;
  growNstData(nst,nst_orig);

  for (int ist = 0; ist < nst_inc; ++ist) {
    FOR_I3 spost[nst_orig+ist][i] = spostVec[ist*3+i]+nsp_orig;
    znost[nst_orig+ist] = izone;
  }

  return 0;

}

// this routine is useful for orienting a surface with disjoint parts
// it is required prior to using the VoxelGrid retriangulation when the
// normals are incorrect (CAD seems to do this a lot)

int SimpleSurface::alignNormalsUsingBloatedSurface(const double delta,const bool seed_with_bloat) {

  cout << "SimpleSurface::alignNormalsUsingBloatedSurface: delta: " << delta << ", nst: " << nst << "..." << endl;

  if (seed_with_bloat)
    ensureTeost(); // require teost on original tris for seeding algo

  // put tris in an adt...

  double (*bbmin)[3] = new double[nst][3];
  double (*bbmax)[3] = new double[nst][3];
  double xstart[3] = { HUGE_VAL, 0.0, 0.0 };
  for (int ist = 0; ist < nst; ++ist) {
    const double * const x0 = xsp[spost[ist][0]];
    const double * const x1 = xsp[spost[ist][1]];
    const double * const x2 = xsp[spost[ist][2]];
    FOR_I3 bbmin[ist][i] = min(x0[i],min(x1[i],x2[i]));
    FOR_I3 bbmax[ist][i] = max(x0[i],max(x1[i],x2[i]));
    // find the x0-most point...
    if (x0[0] < xstart[0]) FOR_I3 xstart[i] = x0[i];
    if (x1[0] < xstart[0]) FOR_I3 xstart[i] = x1[i];
    if (x2[0] < xstart[0]) FOR_I3 xstart[i] = x2[i];
  }

  Adt<double> *adt = new Adt<double>(nst,bbmin,bbmax);
  delete[] bbmin; bbmin = NULL;
  delete[] bbmax; bbmax = NULL;

  // the grid size...

  const double dx = delta;

  // locate x0 at the w/s/b corner. This keeps all
  // the indexing positive. And find a starting point...

  double x0[3],x1[3];
  adt->getBbox(x0,x1);
  FOR_I3 x0[i] -= delta*1.1;
  FOR_I3 x1[i] += delta*1.1;
  xstart[0] -= delta;

  // check how big the grid is going to be...
  cout << " > voxel grid size: ";
  FOR_I3 {
    const int n = int((x1[i]-x0[i]+2.2*delta)/dx);
    assert(n < 2097152); // each dimension only gets 2^21 bits
    cout << n << " ";
  }
  cout << endl;

  // decide on the ijk...
  int i = int((xstart[0]-x0[0])/dx); assert(i >= 0);
  int j = int((xstart[1]-x0[1])/dx); assert(j >= 0);
  int k = int((xstart[2]-x0[2])/dx); assert(k >= 0);

  vector<int> triVec;
  double x[3];
  x[0] = x0[0] + dx*i;
  x[1] = x0[1] + dx*j;
  x[2] = x0[2] + dx*k;
  //cout << "X1: " << COUT_VEC(x) << endl;
  assert(triVec.empty());
  adt->buildListForClosestPoint(triVec,x);
  double d2 = HUGE_VAL;
  for (vector<int>::const_iterator it=triVec.begin(); it!=triVec.end(); ++it) {
    d2 = min(d2,MiscUtils::getPointToTriDist2(x,xsp[spost[*it][0]],xsp[spost[*it][1]],xsp[spost[*it][2]]));
  }
  triVec.clear();

  if (sqrt(d2) < delta) {

    while (sqrt(d2) < delta) {

      // go further away in the -i direction...

      --i; assert(i >= 0);
      x[0] = x0[0] + dx*i;
      x[1] = x0[1] + dx*j;
      x[2] = x0[2] + dx*k;
      assert(triVec.empty());
      adt->buildListForClosestPoint(triVec,x);
      d2 = HUGE_VAL;
      for (vector<int>::const_iterator it=triVec.begin(); it!=triVec.end(); ++it) {
        d2 = min(d2,MiscUtils::getPointToTriDist2(x,xsp[spost[*it][0]],xsp[spost[*it][1]],xsp[spost[*it][2]]));
      }
      triVec.clear();

    }

  }
  else {

    assert(sqrt(d2) >= delta);
    while (sqrt(d2) >= delta) {

      // come closer...

      ++i;
      x[0] = x0[0] + dx*i;
      x[1] = x0[1] + dx*j;
      x[2] = x0[2] + dx*k;
      assert(triVec.empty());
      adt->buildListForClosestPoint(triVec,x);
      d2 = HUGE_VAL;
      for (vector<int>::const_iterator it=triVec.begin(); it!=triVec.end(); ++it) {
        d2 = min(d2,MiscUtils::getPointToTriDist2(x,xsp[spost[*it][0]],xsp[spost[*it][1]],xsp[spost[*it][2]]));
      }
      triVec.clear();

      //cout << "++i d2: " << sqrt(d2) << " delta: " << delta << endl;

      // put i back because cube is staggered: (cube i has nodes i->i+1)
      --i;

    }

  }

  // the cv we want to push onto the stack is (i,j,k)...
  stack<uint8> cvStack;
  cvStack.push((uint8(k)<<42)|(uint8(j)<<21)|uint8(i));
  set<uint8> cvSet;
  cvSet.insert((uint8(k)<<42)|(uint8(j)<<21)|uint8(i));
  double cube[2][2][2];
  map<const uint8,double> sdMap; // signed distance map
  int ned = 0;
  map<const pair<uint8,uint8>,int> edMap;
  vector<int> spostVec;
  int vertlist[12];
  vector<double> xspVec;

  while (!cvStack.empty()) {

    const uint8 ijk = cvStack.top(); cvStack.pop();
    const int i = int(ijk&2097151);
    const int j = int((ijk>>21)&2097151);
    const int k = int(ijk>>42);

    //cout << " just popped: " << i << " " << j << " " << k << endl;

    assert(i >= 0);
    assert(j >= 0);
    assert(k >= 0);

    // put a value in each corner of the cube...

    for (int di = 0; di <= 1; ++di) {
      for (int dj = 0; dj <= 1; ++dj) {
        for (int dk = 0; dk <= 1; ++dk) {
          const uint8 this_ijk = ((uint8(k+dk)<<42)|(uint8(j+dj)<<21)|uint8(i+di));
          map<const uint8,double>::iterator iter = sdMap.find(this_ijk);
          if (iter == sdMap.end()) {
            x[0] = x0[0]+dx*(i+di);
            x[1] = x0[1]+dx*(j+dj);
            x[2] = x0[2]+dx*(k+dk);
            assert(triVec.empty());
            adt->buildListForClosestPoint(triVec,x);
            double d2 = HUGE_VAL;
            for (vector<int>::const_iterator it=triVec.begin(); it!=triVec.end(); ++it) {
              d2 = min(d2,MiscUtils::getPointToTriDist2(x,xsp[spost[*it][0]],xsp[spost[*it][1]],xsp[spost[*it][2]]));
            }
            triVec.clear();
            sdMap[this_ijk] = cube[di][dj][dk] = sqrt(d2) - delta;
          }
          else {
            cube[di][dj][dk] = iter->second;
          }
        }
      }
    }

    // at this point, the 8 cube vertex values are set, so apply the marching cubes algo to build triangles...

    int cube_index = 0;
    if (cube[0][0][0] < 0.0) cube_index |= 1;
    if (cube[1][0][0] < 0.0) cube_index |= 2;
    if (cube[1][1][0] < 0.0) cube_index |= 4;
    if (cube[0][1][0] < 0.0) cube_index |= 8;
    if (cube[0][0][1] < 0.0) cube_index |= 16;
    if (cube[1][0][1] < 0.0) cube_index |= 32;
    if (cube[1][1][1] < 0.0) cube_index |= 64;
    if (cube[0][1][1] < 0.0) cube_index |= 128;

    // cube is entirely in/out of the surface
    // this should not be the case...

    assert(mcEdgeTable[cube_index] != 0);

    // find the vertices where the surface intersects the cube

    for (int i = 0; i < 12; ++i) vertlist[i] = -1;

    if (mcEdgeTable[cube_index] & 1) {
      const pair<uint8,uint8> ijk_pair = pair<uint8,uint8>((uint8(k)<<42)|(uint8(j)<<21)|uint8(i),(uint8(k)<<42)|(uint8(j)<<21)|uint8(i+1));
      map<const pair<uint8,uint8>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[0] = ned++;
        xspVec.push_back((cube[1][0][0]*(x0[0]+dx*(i))-cube[0][0][0]*(x0[0]+dx*(i+1)))/(cube[1][0][0]-cube[0][0][0]));
        xspVec.push_back(x0[1]+dx*(j));
        xspVec.push_back(x0[2]+dx*(k));
      }
      else {
        vertlist[0] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 2) {
      const pair<uint8,uint8> ijk_pair = pair<uint8,uint8>((uint8(k)<<42)|(uint8(j)<<21)|uint8(i+1),(uint8(k)<<42)|(uint8(j+1)<<21)|uint8(i+1));
      map<const pair<uint8,uint8>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[1] = ned++;
        xspVec.push_back(x0[0]+dx*(i+1));
        xspVec.push_back((cube[1][1][0]*(x0[1]+dx*(j))-cube[1][0][0]*(x0[1]+dx*(j+1)))/(cube[1][1][0]-cube[1][0][0]));
        xspVec.push_back(x0[2]+dx*(k));
      }
      else {
        vertlist[1] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 4) {
      const pair<uint8,uint8> ijk_pair = pair<uint8,uint8>((uint8(k)<<42)|(uint8(j+1)<<21)|uint8(i),(uint8(k)<<42)|(uint8(j+1)<<21)|uint8(i+1));
      map<const pair<uint8,uint8>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[2] = ned++;
        xspVec.push_back((cube[1][1][0]*(x0[0]+dx*(i))-cube[0][1][0]*(x0[0]+dx*(i+1)))/(cube[1][1][0]-cube[0][1][0]));
        xspVec.push_back(x0[1]+dx*(j+1));
        xspVec.push_back(x0[2]+dx*(k));
      }
      else {
        vertlist[2] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 8) {
      const pair<uint8,uint8> ijk_pair = pair<uint8,uint8>((uint8(k)<<42)|(uint8(j)<<21)|uint8(i),(uint8(k)<<42)|(uint8(j+1)<<21)|uint8(i));
      map<const pair<uint8,uint8>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[3] = ned++;
        xspVec.push_back(x0[0]+dx*(i));
        xspVec.push_back((cube[0][1][0]*(x0[1]+dx*(j))-cube[0][0][0]*(x0[1]+dx*(j+1)))/(cube[0][1][0]-cube[0][0][0]));
        xspVec.push_back(x0[2]+dx*(k));
      }
      else {
        vertlist[3] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 16) {
      const pair<uint8,uint8> ijk_pair = pair<uint8,uint8>((uint8(k+1)<<42)|(uint8(j)<<21)|uint8(i),(uint8(k+1)<<42)|(uint8(j)<<21)|uint8(i+1));
      map<const pair<uint8,uint8>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[4] = ned++;
        xspVec.push_back((cube[1][0][1]*(x0[0]+dx*(i))-cube[0][0][1]*(x0[0]+dx*(i+1)))/(cube[1][0][1]-cube[0][0][1]));
        xspVec.push_back(x0[1]+dx*(j));
        xspVec.push_back(x0[2]+dx*(k+1));
      }
      else {
        vertlist[4] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 32) {
      const pair<uint8,uint8> ijk_pair = pair<uint8,uint8>((uint8(k+1)<<42)|(uint8(j)<<21)|uint8(i+1),(uint8(k+1)<<42)|(uint8(j+1)<<21)|uint8(i+1));
      map<const pair<uint8,uint8>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[5] = ned++;
        xspVec.push_back(x0[0]+dx*(i+1));
        xspVec.push_back((cube[1][1][1]*(x0[1]+dx*(j))-cube[1][0][1]*(x0[1]+dx*(j+1)))/(cube[1][1][1]-cube[1][0][1]));
        xspVec.push_back(x0[2]+dx*(k+1));
      }
      else {
        vertlist[5] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 64) {
      const pair<uint8,uint8> ijk_pair = pair<uint8,uint8>((uint8(k+1)<<42)|(uint8(j+1)<<21)|uint8(i),(uint8(k+1)<<42)|(uint8(j+1)<<21)|uint8(i+1));
      map<const pair<uint8,uint8>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[6] = ned++;
        xspVec.push_back((cube[1][1][1]*(x0[0]+dx*(i))-cube[0][1][1]*(x0[0]+dx*(i+1)))/(cube[1][1][1]-cube[0][1][1]));
        xspVec.push_back(x0[1]+dx*(j+1));
        xspVec.push_back(x0[2]+dx*(k+1));
      }
      else {
        vertlist[6] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 128) {
      const pair<uint8,uint8> ijk_pair = pair<uint8,uint8>((uint8(k+1)<<42)|(uint8(j)<<21)|uint8(i),(uint8(k+1)<<42)|(uint8(j+1)<<21)|uint8(i));
      map<const pair<uint8,uint8>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[7] = ned++;
        xspVec.push_back(x0[0]+dx*(i));
        xspVec.push_back((cube[0][1][1]*(x0[1]+dx*(j))-cube[0][0][1]*(x0[1]+dx*(j+1)))/(cube[0][1][1]-cube[0][0][1]));
        xspVec.push_back(x0[2]+dx*(k+1));
      }
      else {
        vertlist[7] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 256) {
      const pair<uint8,uint8> ijk_pair = pair<uint8,uint8>((uint8(k)<<42)|(uint8(j)<<21)|uint8(i),(uint8(k+1)<<42)|(uint8(j)<<21)|uint8(i));
      map<const pair<uint8,uint8>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[8] = ned++;
        xspVec.push_back(x0[0]+dx*(i));
        xspVec.push_back(x0[1]+dx*(j));
        xspVec.push_back((cube[0][0][1]*(x0[2]+dx*(k))-cube[0][0][0]*(x0[2]+dx*(k+1)))/(cube[0][0][1]-cube[0][0][0]));
      }
      else {
        vertlist[8] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 512) {
      const pair<uint8,uint8> ijk_pair = pair<uint8,uint8>((uint8(k)<<42)|(uint8(j)<<21)|uint8(i+1),(uint8(k+1)<<42)|(uint8(j)<<21)|uint8(i+1));
      map<const pair<uint8,uint8>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[9] = ned++;
        xspVec.push_back(x0[0]+dx*(i+1));
        xspVec.push_back(x0[1]+dx*(j));
        xspVec.push_back((cube[1][0][1]*(x0[2]+dx*(k))-cube[1][0][0]*(x0[2]+dx*(k+1)))/(cube[1][0][1]-cube[1][0][0]));
      }
      else {
        vertlist[9] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 1024) {
      const pair<uint8,uint8> ijk_pair = pair<uint8,uint8>((uint8(k)<<42)|(uint8(j+1)<<21)|uint8(i+1),(uint8(k+1)<<42)|(uint8(j+1)<<21)|uint8(i+1));
      map<const pair<uint8,uint8>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[10] = ned++;
        xspVec.push_back(x0[0]+dx*(i+1));
        xspVec.push_back(x0[1]+dx*(j+1));
        xspVec.push_back((cube[1][1][1]*(x0[2]+dx*(k))-cube[1][1][0]*(x0[2]+dx*(k+1)))/(cube[1][1][1]-cube[1][1][0]));
      }
      else {
        vertlist[10] = iter->second;
      }
    }
    if (mcEdgeTable[cube_index] & 2048) {
      const pair<uint8,uint8> ijk_pair = pair<uint8,uint8>((uint8(k)<<42)|(uint8(j+1)<<21)|uint8(i),(uint8(k+1)<<42)|(uint8(j+1)<<21)|uint8(i));
      map<const pair<uint8,uint8>,int>::iterator iter = edMap.find(ijk_pair);
      if (iter == edMap.end()) {
        edMap[ijk_pair] = vertlist[11] = ned++;
        xspVec.push_back(x0[0]+dx*(i));
        xspVec.push_back(x0[1]+dx*(j+1));
        xspVec.push_back((cube[0][1][1]*(x0[2]+dx*(k))-cube[0][1][0]*(x0[2]+dx*(k+1)))/(cube[0][1][1]-cube[0][1][0]));
      }
      else {
        vertlist[11] = iter->second;
      }
    }

    // add triangles...

    for (int l = 0; mcTriTable[cube_index][l] != -1; l += 3) {
      spostVec.push_back(vertlist[mcTriTable[cube_index][l  ]]);
      spostVec.push_back(vertlist[mcTriTable[cube_index][l+1]]);
      spostVec.push_back(vertlist[mcTriTable[cube_index][l+2]]);
    }

    // now travel across any face that has both positive and negative
    // cvs...

    int count = 0;
    for (int dj = 0; dj <= 1; ++dj) {
      for (int dk = 0; dk <= 1; ++dk) {
        if (cube[0][dj][dk] >= 0.0)
          ++count;
      }
    }
    //cout << "x0: " << count << endl;
    if ((count >= 1)&&(count <= 3)) {
      const uint8 this_ijk = ((uint8(k)<<42)|(uint8(j)<<21)|uint8(i-1));
      if (cvSet.find(this_ijk) == cvSet.end()) {
        cvSet.insert(this_ijk);
        cvStack.push(this_ijk);
      }
    }

    count = 0;
    for (int dj = 0; dj <= 1; ++dj) {
      for (int dk = 0; dk <= 1; ++dk) {
        if (cube[1][dj][dk] >= 0.0)
          ++count;
      }
    }
    //cout << "x1: " << count << endl;
    if ((count >= 1)&&(count <= 3)) {
      const uint8 this_ijk = ((uint8(k)<<42)|(uint8(j)<<21)|uint8(i+1));
      if (cvSet.find(this_ijk) == cvSet.end()) {
        cvSet.insert(this_ijk);
        cvStack.push(this_ijk);
      }
    }

    count = 0;
    for (int di = 0; di <= 1; ++di) {
      for (int dk = 0; dk <= 1; ++dk) {
        if (cube[di][0][dk] >= 0.0)
          ++count;
      }
    }
    //cout << "y0: " << count << endl;
    if ((count >= 1)&&(count <= 3)) {
      const uint8 this_ijk = ((uint8(k)<<42)|(uint8(j-1)<<21)|uint8(i));
      if (cvSet.find(this_ijk) == cvSet.end()) {
        cvSet.insert(this_ijk);
        cvStack.push(this_ijk);
      }
    }

    count = 0;
    for (int di = 0; di <= 1; ++di) {
      for (int dk = 0; dk <= 1; ++dk) {
        if (cube[di][1][dk] >= 0.0)
          ++count;
      }
    }
    //cout << "y1: " << count << endl;
    if ((count >= 1)&&(count <= 3)) {
      const uint8 this_ijk = ((uint8(k)<<42)|(uint8(j+1)<<21)|uint8(i));
      if (cvSet.find(this_ijk) == cvSet.end()) {
        cvSet.insert(this_ijk);
        cvStack.push(this_ijk);
      }
    }

    count = 0;
    for (int di = 0; di <= 1; ++di) {
      for (int dj = 0; dj <= 1; ++dj) {
        if (cube[di][dj][0] >= 0.0)
          ++count;
      }
    }
    //cout << "z0: " << count << endl;
    if ((count >= 1)&&(count <= 3)) {
      const uint8 this_ijk = ((uint8(k-1)<<42)|(uint8(j)<<21)|uint8(i));
      if (cvSet.find(this_ijk) == cvSet.end()) {
        cvSet.insert(this_ijk);
        cvStack.push(this_ijk);
      }
    }

    count = 0;
    for (int di = 0; di <= 1; ++di) {
      for (int dj = 0; dj <= 1; ++dj) {
        if (cube[di][dj][1] >= 0.0)
          ++count;
      }
    }
    //cout << "z1: " << count << endl;
    if ((count >= 1)&&(count <= 3)) {
      const uint8 this_ijk = ((uint8(k+1)<<42)|(uint8(j)<<21)|uint8(i));
      if (cvSet.find(this_ijk) == cvSet.end()) {
        cvSet.insert(this_ijk);
        cvStack.push(this_ijk);
      }
    }

  }
  delete adt; adt = NULL;

  cout << " > out of stack loop. ned: " << ned << " spostVec.size(): " << spostVec.size() << endl;

  // use the edMap to build a bunch of new nodes...

  edMap.clear();
  assert(int(xspVec.size()) == ned*3);

  // add the new points...
  int nsp_orig = nsp; nsp += ned;
  growNspData(nsp,nsp_orig);

  for (int ied = 0; ied < ned; ++ied) {
    FOR_I3 xsp[nsp_orig+ied][i] = xspVec[ied*3+i];
  }
  xspVec.clear();

  // now add the lifted surface to zone "LIFTED"...
  string zone_name = "BLOATED";
  int izone,nzn;
  for (izone = 0,nzn=zoneVec.size(); izone < nzn; ++izone) {
    if (zoneVec[izone].getName() == zone_name)
      break;
  }
  if (izone == int(zoneVec.size())) {
    zoneVec.push_back(SurfaceZone(zone_name));
  }

  // and add the new tris...
  assert(spostVec.size()%3 == 0);
  const int nst_inc = spostVec.size()/3;
  int nst_orig = nst; nst += nst_inc;
  growNstData(nst,nst_orig);

  for (int ist = 0; ist < nst_inc; ++ist) {
    FOR_I3 spost[nst_orig+ist][i] = spostVec[ist*3+i]+nsp_orig;
    znost[nst_orig+ist] = izone;
  }

  // put tris in an adt...

  assert(bbmin == NULL); bbmin = new double[nst_inc][3];
  assert(bbmax == NULL); bbmax = new double[nst_inc][3];
  for (int ist = nst_orig; ist < nst; ++ist) {
    const double * const x0 = xsp[spost[ist][0]];
    const double * const x1 = xsp[spost[ist][1]];
    const double * const x2 = xsp[spost[ist][2]];
    FOR_I3 bbmin[ist-nst_orig][i] = min(x0[i],min(x1[i],x2[i]));
    FOR_I3 bbmax[ist-nst_orig][i] = max(x0[i],max(x1[i],x2[i]));
  }
  assert(adt == NULL); adt = new Adt<double>(nst_inc,bbmin,bbmax);
  delete[] bbmin;
  delete[] bbmax;

  // search through original tris, find closest new tri, if normals differ flip...

  int skip_count = 0;
  int flip_count = 0;
  int * st_flag = new int[nst_orig]; // -1 means not visited yet, 0 means visited and orientation correct, 1 means visited and orientation needs flipping
  int ist_seed = -1;
  for (int ist = 0; ist < nst_orig; ++ist) {
    const double * const x0 = xsp[spost[ist][0]];
    const double * const x1 = xsp[spost[ist][1]];
    const double * const x2 = xsp[spost[ist][2]];
    double x_st[3]; FOR_I3 x_st[i] = (x0[i]+x1[i]+x2[i])/3.0;
    const double n_st[3] = TRI_NORMAL_2(x0,x1,x2);
    assert(triVec.empty());
    //adt->buildListForClosestPoint(triVec,x_st);
    adt->buildListForSphere(triVec,x_st,2.0*delta);
    if (!triVec.empty()) {
      double d2 = HUGE_VAL;
      int ist2 = -1;
      for (vector<int>::const_iterator it=triVec.begin(); it!=triVec.end(); ++it) {
        const int this_ist2 = *it+nst_orig; assert((this_ist2 >= nst_orig)&&(this_ist2 < nst));
        const double * const this_x02 = xsp[spost[this_ist2][0]];
        const double * const this_x12 = xsp[spost[this_ist2][1]];
        const double * const this_x22 = xsp[spost[this_ist2][2]];
        double this_x_st2[3]; FOR_I3 this_x_st2[i] = (this_x02[i]+this_x12[i]+this_x22[i])/3.0;
        const double dx[3] = DIFF(x_st,this_x_st2);
        const double this_d2 = fabs(DOT_PRODUCT(n_st,dx));
        if (this_d2 < d2) {
          d2 = this_d2;
          ist2 = this_ist2;
        }
      }
      assert(ist2 != -1);
      triVec.clear();
      const double * const x02 = xsp[spost[ist2][0]];
      const double * const x12 = xsp[spost[ist2][1]];
      const double * const x22 = xsp[spost[ist2][2]];
      const double n_st2[3] = TRI_NORMAL_2(x02,x12,x22);
      // note that lifted surface is inside out, so we flip when they are of the same sign
      if (DOT_PRODUCT(n_st,n_st2) > 0.0) {
        const int tmp = spost[ist][0];
        spost[ist][0] = spost[ist][1];
        spost[ist][1] = tmp;
        st_flag[ist] = 1;
        ++flip_count;
      }
      st_flag[ist] = 0;
      if (ist_seed == -1)
        ist_seed = ist;
    }
    else {
      st_flag[ist] = -1;
      ++skip_count;
    }
  }
  delete adt;

  cout << " > flipped tris: " << flip_count << ", skipped tris: " << skip_count << " out of " << nst_orig << endl;

  // if you skip some tris, you probably want to follow up with the traditional normal alignment using the touched
  // tris as seed points

  if ((seed_with_bloat)&&(skip_count > 0)) {
    assert((ist_seed >= 0)&&(ist_seed < nst_orig));

    cout << " > using unskipped tris to seed align normals walking algorithm..." << endl;

    // put a flag in each tri indicating if it has been visited yet, and whether
    // it should be flipped...

    stack<int> stStack;
    stStack.push(ist_seed);

    int disjoint_count = 0;
    int ist0 = 0;
    bool done = false;
    while (!done) {

      while (!stStack.empty()) {
        const int ist = stStack.top(); stStack.pop();
        assert(st_flag[ist] != -1); // should be set
        // visit our nbrs...
        FOR_I3 {
          int ist_nbr,i_nbr,orient_nbr;
          if (getTriNbrData(ist_nbr,i_nbr,orient_nbr,ist,i)) {
            if (st_flag[ist_nbr] == -1) {
              // this nbr has not been visited.
              if (orient_nbr == 1) {
                // nbr is flipped wrt me...
                st_flag[ist_nbr] = 1 - st_flag[ist];
              }
              else {
                assert(orient_nbr == 0);
                st_flag[ist_nbr] = st_flag[ist];
              }
              //st_flag[ist_nbr] = orient_nbr + st_flag[ist] - 2*orient_nbr*st_flag[ist];
              assert((st_flag[ist_nbr] == 0)||(st_flag[ist_nbr] == 1));
              stStack.push(ist_nbr);
            }
            else {
              // should check if multiedge here - the below assert could be caused by that...?
              //assert(st_flag[ist_nbr] == orient_nbr + st_flag[ist] - 2*orient_nbr*st_flag[ist]); // surface is mobius / folded
            }
          }
        }
      }

      ++disjoint_count;
      done = true;
      for (int ist = ist0; ist < nst_orig; ++ist) {
        if (st_flag[ist] == -1) {
          // surface must be multiply-connected. Insert another seed and continue...
          st_flag[ist] = 0; // 0 means leave in current orientation
          stStack.push(ist);
          ist0 = ist;
          done = false;
          break;
        }
        if (st_flag[ist] == 1) {
          ++flip_count;
          // flip this tri by exchanging first 2 nodes...
          const int tmp = spost[ist][0];
          spost[ist][0] = spost[ist][1];
          spost[ist][1] = tmp;
        }
      }

    }

    cout << " > flipped " << flip_count << " of " << nst_orig << " tris in " << disjoint_count << " walkable regions." << endl;


  }
  delete[] st_flag;

  clearTeost();
  clearDynamicEdgeGroups();

  return 0;
}
