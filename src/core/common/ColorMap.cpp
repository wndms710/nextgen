#include "ColorMap.hpp"

//static ColorMap assignements to appear only once in the objects list
const float ColorMap_HOT_METAL::local_r[] = { 2, 127, 255, 255, 255};
const float ColorMap_HOT_METAL::local_g[] = { 0, 0, 0, 255, 255};
const float ColorMap_HOT_METAL::local_b[] = { 0, 0, 0, 0, 255};

const float ColorMap_BLUE_RED_RAINBOW::local_r[] = { 0, 0, 0, 1, 1};
const float ColorMap_BLUE_RED_RAINBOW::local_g[] = { 0, 1, 1, 1, 0};
const float ColorMap_BLUE_RED_RAINBOW::local_b[] = { 1, 1, 0, 0, 0};

const float ColorMap_COOL_WARM::local_r[] = { 0, 0.75, 0.75};
const float ColorMap_COOL_WARM::local_g[] = { 0, 0.75, 0.1};
const float ColorMap_COOL_WARM::local_b[] = { 0.75, 0.75, 0};

const float ColorMap_RAINBOW_DESATURATED::local_r[] = { 0.278431, 0,        0, 0,        1, 1,        0.419608, 0.878431 };
const float ColorMap_RAINBOW_DESATURATED::local_g[] = { 0.278431, 0,        1, 0.501961, 1, 0.380392, 0,        0.301961 };
const float ColorMap_RAINBOW_DESATURATED::local_b[] = { 0.858824, 0.360784, 1, 0,        0, 0,        0,        0.301961 };

map<const string,ColorMap_Base*> ColorMap::colormapsMap;

