/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001-3 Paul Murrell
 *                2003 The R Development Core Team
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
 *  A copy of the GNU General Public License is available via WWW at
 *  http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
 *  writing to the Free Software Foundation, Inc., 59 Temple Place,
 *  Suite 330, Boston, MA  02111-1307  USA.
 */

/* Code to register grid functions with R
 */
#include <R_ext/Rdynload.h>
#include "grid.h"

static const R_CallMethodDef callMethods[] = {{"L_initGrid", (DL_FUNC)&L_initGrid, 1},
                                              {"L_killGrid", (DL_FUNC)&L_killGrid, 0},
                                              {"L_gridDirty", (DL_FUNC)&L_gridDirty, 0},
                                              {"L_currentViewport", (DL_FUNC)&L_currentViewport, 0},
                                              {"L_setviewport", (DL_FUNC)&L_setviewport, 2},
                                              {"L_downviewport", (DL_FUNC)&L_downviewport, 2},
                                              {"L_downvppath", (DL_FUNC)&L_downvppath, 3},
                                              {"L_unsetviewport", (DL_FUNC)&L_unsetviewport, 1},
                                              {"L_upviewport", (DL_FUNC)&L_upviewport, 1},
                                              {"L_getDisplayList", (DL_FUNC)&L_getDisplayList, 0},
                                              {"L_setDisplayList", (DL_FUNC)&L_setDisplayList, 1},
                                              {"L_getDLelt", (DL_FUNC)&L_getDLelt, 1},
                                              {"L_setDLelt", (DL_FUNC)&L_setDLelt, 1},
                                              {"L_getDLindex", (DL_FUNC)&L_getDLindex, 0},
                                              {"L_setDLindex", (DL_FUNC)&L_setDLindex, 1},
                                              {"L_getDLon", (DL_FUNC)&L_getDLon, 0},
                                              {"L_setDLon", (DL_FUNC)&L_setDLon, 1},
                                              {"L_getEngineDLon", (DL_FUNC)&L_getEngineDLon, 0},
                                              {"L_setEngineDLon", (DL_FUNC)&L_setEngineDLon, 1},
                                              {"L_getAsk", (DL_FUNC)&L_getAsk, 0},
                                              {"L_setAsk", (DL_FUNC)&L_setAsk, 1},
                                              {"L_getCurrentGrob", (DL_FUNC)&L_getCurrentGrob, 0},
                                              {"L_setCurrentGrob", (DL_FUNC)&L_setCurrentGrob, 1},
                                              {"L_getEngineRecording", (DL_FUNC)&L_getEngineRecording, 0},
                                              {"L_setEngineRecording", (DL_FUNC)&L_setEngineRecording, 1},
                                              {"L_currentGPar", (DL_FUNC)&L_currentGPar, 0},
                                              {"L_newpagerecording", (DL_FUNC)&L_newpagerecording, 0},
                                              {"L_newpage", (DL_FUNC)&L_newpage, 0},
                                              {"L_initGPar", (DL_FUNC)&L_initGPar, 0},
                                              {"L_initViewportStack", (DL_FUNC)&L_initViewportStack, 0},
                                              {"L_initDisplayList", (DL_FUNC)&L_initDisplayList, 0},
                                              {"L_moveTo", (DL_FUNC)&L_moveTo, 2},
                                              {"L_lineTo", (DL_FUNC)&L_lineTo, 2},
                                              {"L_lines", (DL_FUNC)&L_lines, 2},
                                              {"L_segments", (DL_FUNC)&L_segments, 4},
                                              {"L_arrows", (DL_FUNC)&L_arrows, 12},
                                              {"L_polygon", (DL_FUNC)&L_polygon, 3},
                                              {"L_circle", (DL_FUNC)&L_circle, 3},
                                              {"L_rect", (DL_FUNC)&L_rect, 5},
                                              {"L_text", (DL_FUNC)&L_text, 6},
                                              {"L_points", (DL_FUNC)&L_points, 4},
                                              {"L_pretty", (DL_FUNC)&L_pretty, 1},
                                              {"L_locator", (DL_FUNC)&L_locator, 0},
                                              {"L_convert", (DL_FUNC)&L_convert, 4},
                                              {"L_layoutRegion", (DL_FUNC)&L_layoutRegion, 2},
                                              {NULL, NULL, 0}};

void R_init_grid(DllInfo *dll)
{
    R_useDynamicSymbols(dll, FALSE);
    /* No .C, .Fortran, or .External routines => NULL
     */
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
}
