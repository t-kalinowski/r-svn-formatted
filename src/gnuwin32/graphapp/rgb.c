/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998--1999  Guido Masarotto
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 *
 * File: rgb.c --
 * Guido Masarotto (July, 1998)
 *
 */

/*
   This file is an add-on  to GraphApp, a cross-platform C graphics library.
 */

#include "ga.h"
#include <string.h>

#define RGBCOLORS 657

/* not static and NULL terminated => we can use it in list boxes..*/
char *ColorName[] = {"AliceBlue",
                     "AntiqueWhite",
                     "AntiqueWhite1",
                     "AntiqueWhite2",
                     "AntiqueWhite3",
                     "AntiqueWhite4",
                     "aquamarine",
                     "aquamarine1",
                     "aquamarine2",
                     "aquamarine3",
                     "aquamarine4",
                     "azure",
                     "azure1",
                     "azure2",
                     "azure3",
                     "azure4",
                     "beige",
                     "bisque",
                     "bisque1",
                     "bisque2",
                     "bisque3",
                     "bisque4",
                     "black",
                     "BlanchedAlmond",
                     "blue",
                     "blue1",
                     "blue2",
                     "blue3",
                     "blue4",
                     "BlueViolet",
                     "brown",
                     "brown1",
                     "brown2",
                     "brown3",
                     "brown4",
                     "burlywood",
                     "burlywood1",
                     "burlywood2",
                     "burlywood3",
                     "burlywood4",
                     "CadetBlue",
                     "CadetBlue1",
                     "CadetBlue2",
                     "CadetBlue3",
                     "CadetBlue4",
                     "chartreuse",
                     "chartreuse1",
                     "chartreuse2",
                     "chartreuse3",
                     "chartreuse4",
                     "chocolate",
                     "chocolate1",
                     "chocolate2",
                     "chocolate3",
                     "chocolate4",
                     "coral",
                     "coral1",
                     "coral2",
                     "coral3",
                     "coral4",
                     "CornflowerBlue",
                     "cornsilk",
                     "cornsilk1",
                     "cornsilk2",
                     "cornsilk3",
                     "cornsilk4",
                     "cyan",
                     "cyan1",
                     "cyan2",
                     "cyan3",
                     "cyan4",
                     "DarkBlue",
                     "DarkCyan",
                     "DarkGoldenrod",
                     "DarkGoldenrod1",
                     "DarkGoldenrod2",
                     "DarkGoldenrod3",
                     "DarkGoldenrod4",
                     "DarkGray",
                     "DarkGreen",
                     "DarkGrey",
                     "DarkKhaki",
                     "DarkMagenta",
                     "DarkOliveGreen",
                     "DarkOliveGreen1",
                     "DarkOliveGreen2",
                     "DarkOliveGreen3",
                     "DarkOliveGreen4",
                     "DarkOrange",
                     "DarkOrange1",
                     "DarkOrange2",
                     "DarkOrange3",
                     "DarkOrange4",
                     "DarkOrchid",
                     "DarkOrchid1",
                     "DarkOrchid2",
                     "DarkOrchid3",
                     "DarkOrchid4",
                     "DarkRed",
                     "DarkSalmon",
                     "DarkSeaGreen",
                     "DarkSeaGreen1",
                     "DarkSeaGreen2",
                     "DarkSeaGreen3",
                     "DarkSeaGreen4",
                     "DarkSlateBlue",
                     "DarkSlateGray",
                     "DarkSlateGray1",
                     "DarkSlateGray2",
                     "DarkSlateGray3",
                     "DarkSlateGray4",
                     "DarkSlateGrey",
                     "DarkTurquoise",
                     "DarkViolet",
                     "DeepPink",
                     "DeepPink1",
                     "DeepPink2",
                     "DeepPink3",
                     "DeepPink4",
                     "DeepSkyBlue",
                     "DeepSkyBlue1",
                     "DeepSkyBlue2",
                     "DeepSkyBlue3",
                     "DeepSkyBlue4",
                     "DimGray",
                     "DimGrey",
                     "DodgerBlue",
                     "DodgerBlue1",
                     "DodgerBlue2",
                     "DodgerBlue3",
                     "DodgerBlue4",
                     "firebrick",
                     "firebrick1",
                     "firebrick2",
                     "firebrick3",
                     "firebrick4",
                     "FloralWhite",
                     "ForestGreen",
                     "gainsboro",
                     "GhostWhite",
                     "gold",
                     "gold1",
                     "gold2",
                     "gold3",
                     "gold4",
                     "goldenrod",
                     "goldenrod1",
                     "goldenrod2",
                     "goldenrod3",
                     "goldenrod4",
                     "gray",
                     "gray0",
                     "gray1",
                     "gray10",
                     "gray100",
                     "gray11",
                     "gray12",
                     "gray13",
                     "gray14",
                     "gray15",
                     "gray16",
                     "gray17",
                     "gray18",
                     "gray19",
                     "gray2",
                     "gray20",
                     "gray21",
                     "gray22",
                     "gray23",
                     "gray24",
                     "gray25",
                     "gray26",
                     "gray27",
                     "gray28",
                     "gray29",
                     "gray3",
                     "gray30",
                     "gray31",
                     "gray32",
                     "gray33",
                     "gray34",
                     "gray35",
                     "gray36",
                     "gray37",
                     "gray38",
                     "gray39",
                     "gray4",
                     "gray40",
                     "gray41",
                     "gray42",
                     "gray43",
                     "gray44",
                     "gray45",
                     "gray46",
                     "gray47",
                     "gray48",
                     "gray49",
                     "gray5",
                     "gray50",
                     "gray51",
                     "gray52",
                     "gray53",
                     "gray54",
                     "gray55",
                     "gray56",
                     "gray57",
                     "gray58",
                     "gray59",
                     "gray6",
                     "gray60",
                     "gray61",
                     "gray62",
                     "gray63",
                     "gray64",
                     "gray65",
                     "gray66",
                     "gray67",
                     "gray68",
                     "gray69",
                     "gray7",
                     "gray70",
                     "gray71",
                     "gray72",
                     "gray73",
                     "gray74",
                     "gray75",
                     "gray76",
                     "gray77",
                     "gray78",
                     "gray79",
                     "gray8",
                     "gray80",
                     "gray81",
                     "gray82",
                     "gray83",
                     "gray84",
                     "gray85",
                     "gray86",
                     "gray87",
                     "gray88",
                     "gray89",
                     "gray9",
                     "gray90",
                     "gray91",
                     "gray92",
                     "gray93",
                     "gray94",
                     "gray95",
                     "gray96",
                     "gray97",
                     "gray98",
                     "gray99",
                     "green",
                     "green1",
                     "green2",
                     "green3",
                     "green4",
                     "GreenYellow",
                     "grey",
                     "grey0",
                     "grey1",
                     "grey10",
                     "grey100",
                     "grey11",
                     "grey12",
                     "grey13",
                     "grey14",
                     "grey15",
                     "grey16",
                     "grey17",
                     "grey18",
                     "grey19",
                     "grey2",
                     "grey20",
                     "grey21",
                     "grey22",
                     "grey23",
                     "grey24",
                     "grey25",
                     "grey26",
                     "grey27",
                     "grey28",
                     "grey29",
                     "grey3",
                     "grey30",
                     "grey31",
                     "grey32",
                     "grey33",
                     "grey34",
                     "grey35",
                     "grey36",
                     "grey37",
                     "grey38",
                     "grey39",
                     "grey4",
                     "grey40",
                     "grey41",
                     "grey42",
                     "grey43",
                     "grey44",
                     "grey45",
                     "grey46",
                     "grey47",
                     "grey48",
                     "grey49",
                     "grey5",
                     "grey50",
                     "grey51",
                     "grey52",
                     "grey53",
                     "grey54",
                     "grey55",
                     "grey56",
                     "grey57",
                     "grey58",
                     "grey59",
                     "grey6",
                     "grey60",
                     "grey61",
                     "grey62",
                     "grey63",
                     "grey64",
                     "grey65",
                     "grey66",
                     "grey67",
                     "grey68",
                     "grey69",
                     "grey7",
                     "grey70",
                     "grey71",
                     "grey72",
                     "grey73",
                     "grey74",
                     "grey75",
                     "grey76",
                     "grey77",
                     "grey78",
                     "grey79",
                     "grey8",
                     "grey80",
                     "grey81",
                     "grey82",
                     "grey83",
                     "grey84",
                     "grey85",
                     "grey86",
                     "grey87",
                     "grey88",
                     "grey89",
                     "grey9",
                     "grey90",
                     "grey91",
                     "grey92",
                     "grey93",
                     "grey94",
                     "grey95",
                     "grey96",
                     "grey97",
                     "grey98",
                     "grey99",
                     "honeydew",
                     "honeydew1",
                     "honeydew2",
                     "honeydew3",
                     "honeydew4",
                     "HotPink",
                     "HotPink1",
                     "HotPink2",
                     "HotPink3",
                     "HotPink4",
                     "IndianRed",
                     "IndianRed1",
                     "IndianRed2",
                     "IndianRed3",
                     "IndianRed4",
                     "ivory",
                     "ivory1",
                     "ivory2",
                     "ivory3",
                     "ivory4",
                     "khaki",
                     "khaki1",
                     "khaki2",
                     "khaki3",
                     "khaki4",
                     "lavender",
                     "LavenderBlush",
                     "LavenderBlush1",
                     "LavenderBlush2",
                     "LavenderBlush3",
                     "LavenderBlush4",
                     "LawnGreen",
                     "LemonChiffon",
                     "LemonChiffon1",
                     "LemonChiffon2",
                     "LemonChiffon3",
                     "LemonChiffon4",
                     "LightBlue",
                     "LightBlue1",
                     "LightBlue2",
                     "LightBlue3",
                     "LightBlue4",
                     "LightCoral",
                     "LightCyan",
                     "LightCyan1",
                     "LightCyan2",
                     "LightCyan3",
                     "LightCyan4",
                     "LightGoldenrod",
                     "LightGoldenrod1",
                     "LightGoldenrod2",
                     "LightGoldenrod3",
                     "LightGoldenrod4",
                     "LightGoldenrodYellow",
                     "LightGray",
                     "LightGreen",
                     "LightGrey",
                     "LightPink",
                     "LightPink1",
                     "LightPink2",
                     "LightPink3",
                     "LightPink4",
                     "LightSalmon",
                     "LightSalmon1",
                     "LightSalmon2",
                     "LightSalmon3",
                     "LightSalmon4",
                     "LightSeaGreen",
                     "LightSkyBlue",
                     "LightSkyBlue1",
                     "LightSkyBlue2",
                     "LightSkyBlue3",
                     "LightSkyBlue4",
                     "LightSlateBlue",
                     "LightSlateGray",
                     "LightSlateGrey",
                     "LightSteelBlue",
                     "LightSteelBlue1",
                     "LightSteelBlue2",
                     "LightSteelBlue3",
                     "LightSteelBlue4",
                     "LightYellow",
                     "LightYellow1",
                     "LightYellow2",
                     "LightYellow3",
                     "LightYellow4",
                     "LimeGreen",
                     "linen",
                     "magenta",
                     "magenta1",
                     "magenta2",
                     "magenta3",
                     "magenta4",
                     "maroon",
                     "maroon1",
                     "maroon2",
                     "maroon3",
                     "maroon4",
                     "MediumAquamarine",
                     "MediumBlue",
                     "MediumOrchid",
                     "MediumOrchid1",
                     "MediumOrchid2",
                     "MediumOrchid3",
                     "MediumOrchid4",
                     "MediumPurple",
                     "MediumPurple1",
                     "MediumPurple2",
                     "MediumPurple3",
                     "MediumPurple4",
                     "MediumSeaGreen",
                     "MediumSlateBlue",
                     "MediumSpringGreen",
                     "MediumTurquoise",
                     "MediumVioletRed",
                     "MidnightBlue",
                     "MintCream",
                     "MistyRose",
                     "MistyRose1",
                     "MistyRose2",
                     "MistyRose3",
                     "MistyRose4",
                     "moccasin",
                     "NavajoWhite",
                     "NavajoWhite1",
                     "NavajoWhite2",
                     "NavajoWhite3",
                     "NavajoWhite4",
                     "navy",
                     "NavyBlue",
                     "OldLace",
                     "OliveDrab",
                     "OliveDrab1",
                     "OliveDrab2",
                     "OliveDrab3",
                     "OliveDrab4",
                     "orange",
                     "orange1",
                     "orange2",
                     "orange3",
                     "orange4",
                     "OrangeRed",
                     "OrangeRed1",
                     "OrangeRed2",
                     "OrangeRed3",
                     "OrangeRed4",
                     "orchid",
                     "orchid1",
                     "orchid2",
                     "orchid3",
                     "orchid4",
                     "PaleGoldenrod",
                     "PaleGreen",
                     "PaleGreen1",
                     "PaleGreen2",
                     "PaleGreen3",
                     "PaleGreen4",
                     "PaleTurquoise",
                     "PaleTurquoise1",
                     "PaleTurquoise2",
                     "PaleTurquoise3",
                     "PaleTurquoise4",
                     "PaleVioletRed",
                     "PaleVioletRed1",
                     "PaleVioletRed2",
                     "PaleVioletRed3",
                     "PaleVioletRed4",
                     "PapayaWhip",
                     "PeachPuff",
                     "PeachPuff1",
                     "PeachPuff2",
                     "PeachPuff3",
                     "PeachPuff4",
                     "peru",
                     "pink",
                     "pink1",
                     "pink2",
                     "pink3",
                     "pink4",
                     "plum",
                     "plum1",
                     "plum2",
                     "plum3",
                     "plum4",
                     "PowderBlue",
                     "purple",
                     "purple1",
                     "purple2",
                     "purple3",
                     "purple4",
                     "red",
                     "red1",
                     "red2",
                     "red3",
                     "red4",
                     "RosyBrown",
                     "RosyBrown1",
                     "RosyBrown2",
                     "RosyBrown3",
                     "RosyBrown4",
                     "RoyalBlue",
                     "RoyalBlue1",
                     "RoyalBlue2",
                     "RoyalBlue3",
                     "RoyalBlue4",
                     "SaddleBrown",
                     "salmon",
                     "salmon1",
                     "salmon2",
                     "salmon3",
                     "salmon4",
                     "SandyBrown",
                     "SeaGreen",
                     "SeaGreen1",
                     "SeaGreen2",
                     "SeaGreen3",
                     "SeaGreen4",
                     "seashell",
                     "seashell1",
                     "seashell2",
                     "seashell3",
                     "seashell4",
                     "sienna",
                     "sienna1",
                     "sienna2",
                     "sienna3",
                     "sienna4",
                     "SkyBlue",
                     "SkyBlue1",
                     "SkyBlue2",
                     "SkyBlue3",
                     "SkyBlue4",
                     "SlateBlue",
                     "SlateBlue1",
                     "SlateBlue2",
                     "SlateBlue3",
                     "SlateBlue4",
                     "SlateGray",
                     "SlateGray1",
                     "SlateGray2",
                     "SlateGray3",
                     "SlateGray4",
                     "SlateGrey",
                     "snow",
                     "snow1",
                     "snow2",
                     "snow3",
                     "snow4",
                     "SpringGreen",
                     "SpringGreen1",
                     "SpringGreen2",
                     "SpringGreen3",
                     "SpringGreen4",
                     "SteelBlue",
                     "SteelBlue1",
                     "SteelBlue2",
                     "SteelBlue3",
                     "SteelBlue4",
                     "tan",
                     "tan1",
                     "tan2",
                     "tan3",
                     "tan4",
                     "thistle",
                     "thistle1",
                     "thistle2",
                     "thistle3",
                     "thistle4",
                     "tomato",
                     "tomato1",
                     "tomato2",
                     "tomato3",
                     "tomato4",
                     "turquoise",
                     "turquoise1",
                     "turquoise2",
                     "turquoise3",
                     "turquoise4",
                     "violet",
                     "VioletRed",
                     "VioletRed1",
                     "VioletRed2",
                     "VioletRed3",
                     "VioletRed4",
                     "wheat",
                     "wheat1",
                     "wheat2",
                     "wheat3",
                     "wheat4",
                     "white",
                     "WhiteSmoke",
                     "yellow",
                     "yellow1",
                     "yellow2",
                     "yellow3",
                     "yellow4",
                     "YellowGreen",
                     NULL};

static int RgbValue[RGBCOLORS][3] = {
    {240, 248, 255}, {250, 235, 215}, {255, 239, 219}, {238, 223, 204}, {205, 192, 176}, {139, 131, 120},
    {127, 255, 212}, {127, 255, 212}, {118, 238, 198}, {102, 205, 170}, {69, 139, 116},  {240, 255, 255},
    {240, 255, 255}, {224, 238, 238}, {193, 205, 205}, {131, 139, 139}, {245, 245, 220}, {255, 228, 196},
    {255, 228, 196}, {238, 213, 183}, {205, 183, 158}, {139, 125, 107}, {0, 0, 0},       {255, 235, 205},
    {0, 0, 255},     {0, 0, 255},     {0, 0, 238},     {0, 0, 205},     {0, 0, 139},     {138, 43, 226},
    {165, 42, 42},   {255, 64, 64},   {238, 59, 59},   {205, 51, 51},   {139, 35, 35},   {222, 184, 135},
    {255, 211, 155}, {238, 197, 145}, {205, 170, 125}, {139, 115, 85},  {95, 158, 160},  {152, 245, 255},
    {142, 229, 238}, {122, 197, 205}, {83, 134, 139},  {127, 255, 0},   {127, 255, 0},   {118, 238, 0},
    {102, 205, 0},   {69, 139, 0},    {210, 105, 30},  {255, 127, 36},  {238, 118, 33},  {205, 102, 29},
    {139, 69, 19},   {255, 127, 80},  {255, 114, 86},  {238, 106, 80},  {205, 91, 69},   {139, 62, 47},
    {100, 149, 237}, {255, 248, 220}, {255, 248, 220}, {238, 232, 205}, {205, 200, 177}, {139, 136, 120},
    {0, 255, 255},   {0, 255, 255},   {0, 238, 238},   {0, 205, 205},   {0, 139, 139},   {0, 0, 139},
    {0, 139, 139},   {184, 134, 11},  {255, 185, 15},  {238, 173, 14},  {205, 149, 12},  {139, 101, 8},
    {169, 169, 169}, {0, 100, 0},     {169, 169, 169}, {189, 183, 107}, {139, 0, 139},   {85, 107, 47},
    {202, 255, 112}, {188, 238, 104}, {162, 205, 90},  {110, 139, 61},  {255, 140, 0},   {255, 127, 0},
    {238, 118, 0},   {205, 102, 0},   {139, 69, 0},    {153, 50, 204},  {191, 62, 255},  {178, 58, 238},
    {154, 50, 205},  {104, 34, 139},  {139, 0, 0},     {233, 150, 122}, {143, 188, 143}, {193, 255, 193},
    {180, 238, 180}, {155, 205, 155}, {105, 139, 105}, {72, 61, 139},   {47, 79, 79},    {151, 255, 255},
    {141, 238, 238}, {121, 205, 205}, {82, 139, 139},  {47, 79, 79},    {0, 206, 209},   {148, 0, 211},
    {255, 20, 147},  {255, 20, 147},  {238, 18, 137},  {205, 16, 118},  {139, 10, 80},   {0, 191, 255},
    {0, 191, 255},   {0, 178, 238},   {0, 154, 205},   {0, 104, 139},   {105, 105, 105}, {105, 105, 105},
    {30, 144, 255},  {30, 144, 255},  {28, 134, 238},  {24, 116, 205},  {16, 78, 139},   {178, 34, 34},
    {255, 48, 48},   {238, 44, 44},   {205, 38, 38},   {139, 26, 26},   {255, 250, 240}, {34, 139, 34},
    {220, 220, 220}, {248, 248, 255}, {255, 215, 0},   {255, 215, 0},   {238, 201, 0},   {205, 173, 0},
    {139, 117, 0},   {218, 165, 32},  {255, 193, 37},  {238, 180, 34},  {205, 155, 29},  {139, 105, 20},
    {190, 190, 190}, {0, 0, 0},       {3, 3, 3},       {26, 26, 26},    {255, 255, 255}, {28, 28, 28},
    {31, 31, 31},    {33, 33, 33},    {36, 36, 36},    {38, 38, 38},    {41, 41, 41},    {43, 43, 43},
    {46, 46, 46},    {48, 48, 48},    {5, 5, 5},       {51, 51, 51},    {54, 54, 54},    {56, 56, 56},
    {59, 59, 59},    {61, 61, 61},    {64, 64, 64},    {66, 66, 66},    {69, 69, 69},    {71, 71, 71},
    {74, 74, 74},    {8, 8, 8},       {77, 77, 77},    {79, 79, 79},    {82, 82, 82},    {84, 84, 84},
    {87, 87, 87},    {89, 89, 89},    {92, 92, 92},    {94, 94, 94},    {97, 97, 97},    {99, 99, 99},
    {10, 10, 10},    {102, 102, 102}, {105, 105, 105}, {107, 107, 107}, {110, 110, 110}, {112, 112, 112},
    {115, 115, 115}, {117, 117, 117}, {120, 120, 120}, {122, 122, 122}, {125, 125, 125}, {13, 13, 13},
    {127, 127, 127}, {130, 130, 130}, {133, 133, 133}, {135, 135, 135}, {138, 138, 138}, {140, 140, 140},
    {143, 143, 143}, {145, 145, 145}, {148, 148, 148}, {150, 150, 150}, {15, 15, 15},    {153, 153, 153},
    {156, 156, 156}, {158, 158, 158}, {161, 161, 161}, {163, 163, 163}, {166, 166, 166}, {168, 168, 168},
    {171, 171, 171}, {173, 173, 173}, {176, 176, 176}, {18, 18, 18},    {179, 179, 179}, {181, 181, 181},
    {184, 184, 184}, {186, 186, 186}, {189, 189, 189}, {191, 191, 191}, {194, 194, 194}, {196, 196, 196},
    {199, 199, 199}, {201, 201, 201}, {20, 20, 20},    {204, 204, 204}, {207, 207, 207}, {209, 209, 209},
    {212, 212, 212}, {214, 214, 214}, {217, 217, 217}, {219, 219, 219}, {222, 222, 222}, {224, 224, 224},
    {227, 227, 227}, {23, 23, 23},    {229, 229, 229}, {232, 232, 232}, {235, 235, 235}, {237, 237, 237},
    {240, 240, 240}, {242, 242, 242}, {245, 245, 245}, {247, 247, 247}, {250, 250, 250}, {252, 252, 252},
    {0, 255, 0},     {0, 255, 0},     {0, 238, 0},     {0, 205, 0},     {0, 139, 0},     {173, 255, 47},
    {190, 190, 190}, {0, 0, 0},       {3, 3, 3},       {26, 26, 26},    {255, 255, 255}, {28, 28, 28},
    {31, 31, 31},    {33, 33, 33},    {36, 36, 36},    {38, 38, 38},    {41, 41, 41},    {43, 43, 43},
    {46, 46, 46},    {48, 48, 48},    {5, 5, 5},       {51, 51, 51},    {54, 54, 54},    {56, 56, 56},
    {59, 59, 59},    {61, 61, 61},    {64, 64, 64},    {66, 66, 66},    {69, 69, 69},    {71, 71, 71},
    {74, 74, 74},    {8, 8, 8},       {77, 77, 77},    {79, 79, 79},    {82, 82, 82},    {84, 84, 84},
    {87, 87, 87},    {89, 89, 89},    {92, 92, 92},    {94, 94, 94},    {97, 97, 97},    {99, 99, 99},
    {10, 10, 10},    {102, 102, 102}, {105, 105, 105}, {107, 107, 107}, {110, 110, 110}, {112, 112, 112},
    {115, 115, 115}, {117, 117, 117}, {120, 120, 120}, {122, 122, 122}, {125, 125, 125}, {13, 13, 13},
    {127, 127, 127}, {130, 130, 130}, {133, 133, 133}, {135, 135, 135}, {138, 138, 138}, {140, 140, 140},
    {143, 143, 143}, {145, 145, 145}, {148, 148, 148}, {150, 150, 150}, {15, 15, 15},    {153, 153, 153},
    {156, 156, 156}, {158, 158, 158}, {161, 161, 161}, {163, 163, 163}, {166, 166, 166}, {168, 168, 168},
    {171, 171, 171}, {173, 173, 173}, {176, 176, 176}, {18, 18, 18},    {179, 179, 179}, {181, 181, 181},
    {184, 184, 184}, {186, 186, 186}, {189, 189, 189}, {191, 191, 191}, {194, 194, 194}, {196, 196, 196},
    {199, 199, 199}, {201, 201, 201}, {20, 20, 20},    {204, 204, 204}, {207, 207, 207}, {209, 209, 209},
    {212, 212, 212}, {214, 214, 214}, {217, 217, 217}, {219, 219, 219}, {222, 222, 222}, {224, 224, 224},
    {227, 227, 227}, {23, 23, 23},    {229, 229, 229}, {232, 232, 232}, {235, 235, 235}, {237, 237, 237},
    {240, 240, 240}, {242, 242, 242}, {245, 245, 245}, {247, 247, 247}, {250, 250, 250}, {252, 252, 252},
    {240, 255, 240}, {240, 255, 240}, {224, 238, 224}, {193, 205, 193}, {131, 139, 131}, {255, 105, 180},
    {255, 110, 180}, {238, 106, 167}, {205, 96, 144},  {139, 58, 98},   {205, 92, 92},   {255, 106, 106},
    {238, 99, 99},   {205, 85, 85},   {139, 58, 58},   {255, 255, 240}, {255, 255, 240}, {238, 238, 224},
    {205, 205, 193}, {139, 139, 131}, {240, 230, 140}, {255, 246, 143}, {238, 230, 133}, {205, 198, 115},
    {139, 134, 78},  {230, 230, 250}, {255, 240, 245}, {255, 240, 245}, {238, 224, 229}, {205, 193, 197},
    {139, 131, 134}, {124, 252, 0},   {255, 250, 205}, {255, 250, 205}, {238, 233, 191}, {205, 201, 165},
    {139, 137, 112}, {173, 216, 230}, {191, 239, 255}, {178, 223, 238}, {154, 192, 205}, {104, 131, 139},
    {240, 128, 128}, {224, 255, 255}, {224, 255, 255}, {209, 238, 238}, {180, 205, 205}, {122, 139, 139},
    {238, 221, 130}, {255, 236, 139}, {238, 220, 130}, {205, 190, 112}, {139, 129, 76},  {250, 250, 210},
    {211, 211, 211}, {144, 238, 144}, {211, 211, 211}, {255, 182, 193}, {255, 174, 185}, {238, 162, 173},
    {205, 140, 149}, {139, 95, 101},  {255, 160, 122}, {255, 160, 122}, {238, 149, 114}, {205, 129, 98},
    {139, 87, 66},   {32, 178, 170},  {135, 206, 250}, {176, 226, 255}, {164, 211, 238}, {141, 182, 205},
    {96, 123, 139},  {132, 112, 255}, {119, 136, 153}, {119, 136, 153}, {176, 196, 222}, {202, 225, 255},
    {188, 210, 238}, {162, 181, 205}, {110, 123, 139}, {255, 255, 224}, {255, 255, 224}, {238, 238, 209},
    {205, 205, 180}, {139, 139, 122}, {50, 205, 50},   {250, 240, 230}, {255, 0, 255},   {255, 0, 255},
    {238, 0, 238},   {205, 0, 205},   {139, 0, 139},   {176, 48, 96},   {255, 52, 179},  {238, 48, 167},
    {205, 41, 144},  {139, 28, 98},   {102, 205, 170}, {0, 0, 205},     {186, 85, 211},  {224, 102, 255},
    {209, 95, 238},  {180, 82, 205},  {122, 55, 139},  {147, 112, 219}, {171, 130, 255}, {159, 121, 238},
    {137, 104, 205}, {93, 71, 139},   {60, 179, 113},  {123, 104, 238}, {0, 250, 154},   {72, 209, 204},
    {199, 21, 133},  {25, 25, 112},   {245, 255, 250}, {255, 228, 225}, {255, 228, 225}, {238, 213, 210},
    {205, 183, 181}, {139, 125, 123}, {255, 228, 181}, {255, 222, 173}, {255, 222, 173}, {238, 207, 161},
    {205, 179, 139}, {139, 121, 94},  {0, 0, 128},     {0, 0, 128},     {253, 245, 230}, {107, 142, 35},
    {192, 255, 62},  {179, 238, 58},  {154, 205, 50},  {105, 139, 34},  {255, 165, 0},   {255, 165, 0},
    {238, 154, 0},   {205, 133, 0},   {139, 90, 0},    {255, 69, 0},    {255, 69, 0},    {238, 64, 0},
    {205, 55, 0},    {139, 37, 0},    {218, 112, 214}, {255, 131, 250}, {238, 122, 233}, {205, 105, 201},
    {139, 71, 137},  {238, 232, 170}, {152, 251, 152}, {154, 255, 154}, {144, 238, 144}, {124, 205, 124},
    {84, 139, 84},   {175, 238, 238}, {187, 255, 255}, {174, 238, 238}, {150, 205, 205}, {102, 139, 139},
    {219, 112, 147}, {255, 130, 171}, {238, 121, 159}, {205, 104, 137}, {139, 71, 93},   {255, 239, 213},
    {255, 218, 185}, {255, 218, 185}, {238, 203, 173}, {205, 175, 149}, {139, 119, 101}, {205, 133, 63},
    {255, 192, 203}, {255, 181, 197}, {238, 169, 184}, {205, 145, 158}, {139, 99, 108},  {221, 160, 221},
    {255, 187, 255}, {238, 174, 238}, {205, 150, 205}, {139, 102, 139}, {176, 224, 230}, {160, 32, 240},
    {155, 48, 255},  {145, 44, 238},  {125, 38, 205},  {85, 26, 139},   {255, 0, 0},     {255, 0, 0},
    {238, 0, 0},     {205, 0, 0},     {139, 0, 0},     {188, 143, 143}, {255, 193, 193}, {238, 180, 180},
    {205, 155, 155}, {139, 105, 105}, {65, 105, 225},  {72, 118, 255},  {67, 110, 238},  {58, 95, 205},
    {39, 64, 139},   {139, 69, 19},   {250, 128, 114}, {255, 140, 105}, {238, 130, 98},  {205, 112, 84},
    {139, 76, 57},   {244, 164, 96},  {46, 139, 87},   {84, 255, 159},  {78, 238, 148},  {67, 205, 128},
    {46, 139, 87},   {255, 245, 238}, {255, 245, 238}, {238, 229, 222}, {205, 197, 191}, {139, 134, 130},
    {160, 82, 45},   {255, 130, 71},  {238, 121, 66},  {205, 104, 57},  {139, 71, 38},   {135, 206, 235},
    {135, 206, 255}, {126, 192, 238}, {108, 166, 205}, {74, 112, 139},  {106, 90, 205},  {131, 111, 255},
    {122, 103, 238}, {105, 89, 205},  {71, 60, 139},   {112, 128, 144}, {198, 226, 255}, {185, 211, 238},
    {159, 182, 205}, {108, 123, 139}, {112, 128, 144}, {255, 250, 250}, {255, 250, 250}, {238, 233, 233},
    {205, 201, 201}, {139, 137, 137}, {0, 255, 127},   {0, 255, 127},   {0, 238, 118},   {0, 205, 102},
    {0, 139, 69},    {70, 130, 180},  {99, 184, 255},  {92, 172, 238},  {79, 148, 205},  {54, 100, 139},
    {210, 180, 140}, {255, 165, 79},  {238, 154, 73},  {205, 133, 63},  {139, 90, 43},   {216, 191, 216},
    {255, 225, 255}, {238, 210, 238}, {205, 181, 205}, {139, 123, 139}, {255, 99, 71},   {255, 99, 71},
    {238, 92, 66},   {205, 79, 57},   {139, 54, 38},   {64, 224, 208},  {0, 245, 255},   {0, 229, 238},
    {0, 197, 205},   {0, 134, 139},   {238, 130, 238}, {208, 32, 144},  {255, 62, 150},  {238, 58, 140},
    {205, 50, 120},  {139, 34, 82},   {245, 222, 179}, {255, 231, 186}, {238, 216, 174}, {205, 186, 150},
    {139, 126, 102}, {255, 255, 255}, {245, 245, 245}, {255, 255, 0},   {255, 255, 0},   {238, 238, 0},
    {205, 205, 0},   {139, 139, 0},   {154, 205, 50}};

/* Return transparent if the color doesn't exist.
 *  Case insensitive comparison?
 */
rgb nametorgb(char *name)
{
    int mid, high, low, cmp;

    low = 0;
    mid = 0; /* for -Wall */
    high = RGBCOLORS - 1;
    while (low <= high)
    {
        mid = (low + high) / 2;
        cmp = strcmpi(name, ColorName[mid]);
        if (cmp < 0)
            high = mid - 1;
        else if (cmp > 0)
            low = mid + 1;
        else
            break;
    }
    if (high < low)
        return Transparent;
    else
        return rgb(RgbValue[mid][0], RgbValue[mid][1], RgbValue[mid][2]);
}

/* return "white" not "gray100" */
char *rgbtoname(rgb in)
{
    int i;

    for (i = 0; i < RGBCOLORS; i++)
        if (in == rgb(RgbValue[i][0], RgbValue[i][1], RgbValue[i][2]) && strcmp(ColorName[i], "gray100") &&
            strcmp(ColorName[i], "grey100"))
            return ColorName[i];
    return "";
}

int rgbtonum(rgb in)
{
    int i;

    for (i = 0; i < RGBCOLORS; i++)
        if (in == rgb(RgbValue[i][0], RgbValue[i][1], RgbValue[i][2]) && strcmp(ColorName[i], "gray100") &&
            strcmp(ColorName[i], "grey100"))
            return i;
    return -1;
}

#include <windows.h>
/* Windows uses 0x00bbggrr ! */
rgb myGetSysColor(int x)
{
    int col = GetSysColor(x);
    return rgb((col)&0xFFUL, (col >> 8) & 0xFFUL, (col >> 16) & 0x00FFUL);
}

rgb dialog_bg()
{
    return myGetSysColor(COLOR_BTNFACE);
}
