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

#include "grid.h"

/* This stuff always returns an LViewportLocation in "npc" units
 */

int layoutNRow(SEXP l)
{
    return INTEGER(getListElement(l, "nrow"))[0];
}

int layoutNCol(SEXP l)
{
    return INTEGER(getListElement(l, "ncol"))[0];
}

SEXP layoutWidths(SEXP l)
{
    return getListElement(l, "widths");
}

SEXP layoutHeights(SEXP l)
{
    return getListElement(l, "heights");
}

int layoutRespect(SEXP l)
{
    return INTEGER(getListElement(l, "valid.respect"))[0];
}

int *layoutRespectMat(SEXP l)
{
    return INTEGER(getListElement(l, "respect.mat"));
}

Rboolean relativeUnit(SEXP unit, int index)
{
    return pureNullUnit(unit, index);
}

void findRelWidths(SEXP layout, int *relativeWidths)
{
    int i;
    SEXP widths = layoutWidths(layout);
    for (i = 0; i < layoutNCol(layout); i++)
        relativeWidths[i] = relativeUnit(widths, i);
}

void findRelHeights(SEXP layout, int *relativeHeights)
{
    int i;
    SEXP heights = layoutHeights(layout);
    for (i = 0; i < layoutNRow(layout); i++)
        relativeHeights[i] = relativeUnit(heights, i);
}

void allocateKnownWidths(SEXP layout, int *relativeWidths, double parentWidthCM, double parentHeightCM,
                         LViewportContext parentContext, GEDevDesc *dd, double *npcWidths, double *widthLeftCM)
{
    int i;
    SEXP widths = layoutWidths(layout);
    for (i = 0; i < layoutNCol(layout); i++)
        if (!relativeWidths[i])
        {
            npcWidths[i] =
                transformWidth(widths, i, parentContext, parentContext.fontfamily, parentContext.font,
                               parentContext.fontsize, parentContext.lineheight, parentWidthCM, parentHeightCM, dd);
            *widthLeftCM -= npcWidths[i] * parentWidthCM;
        }
}

void allocateKnownHeights(SEXP layout, int *relativeHeights, double parentWidthCM, double parentHeightCM,
                          LViewportContext parentContext, GEDevDesc *dd, double *npcHeights, double *heightLeftCM)
{
    int i;
    SEXP heights = layoutHeights(layout);
    for (i = 0; i < layoutNRow(layout); i++)
        if (!relativeHeights[i])
        {
            npcHeights[i] =
                transformHeight(heights, i, parentContext, parentContext.fontfamily, parentContext.font,
                                parentContext.fontsize, parentContext.lineheight, parentWidthCM, parentHeightCM, dd);
            *heightLeftCM -= npcHeights[i] * parentHeightCM;
        }
}

int colRespected(int col, SEXP layout)
{
    int i;
    int result = 0;
    int respect = layoutRespect(layout);
    int *respectMat = layoutRespectMat(layout);
    if (respect == 1)
        result = 1;
    else
        for (i = 0; i < layoutNRow(layout); i++)
            if (respectMat[col * layoutNRow(layout) + i] != 0)
                result = 1;
    return result;
}

int rowRespected(int row, SEXP layout)
{
    int i;
    int result = 0;
    int respect = layoutRespect(layout);
    int *respectMat = layoutRespectMat(layout);
    if (respect == 1)
        result = 1;
    else
        for (i = 0; i < layoutNCol(layout); i++)
            if (respectMat[i * layoutNRow(layout) + row] != 0)
                result = 1;
    return result;
}

/*
 * These sum up ALL relative widths and heights (unit = "null")
 * Some effort is made to find all truly null units
 * (e.g., including a grobwidth unit where the grob's width is null)
 */
double totalWidth(SEXP layout, int *relativeWidths, LViewportContext parentContext, GEDevDesc *dd)
{
    int i;
    SEXP widths = layoutWidths(layout);
    double totalWidth = 0;
    /* We are calculating "null" units for a layout */
    L_nullLayoutMode = 1;
    for (i = 0; i < layoutNCol(layout); i++)
        if (relativeWidths[i])
            totalWidth += transformWidth(widths, i, parentContext, parentContext.fontfamily, parentContext.font,
                                         parentContext.fontsize, parentContext.lineheight,
                                         /*
                                          * NOTE: 0, 0, here is ok
                                          * because we are only
                                          * obtaining "null" units
                                          */
                                         0, 0, dd);
    L_nullLayoutMode = 0;
    return totalWidth;
}

double totalHeight(SEXP layout, int *relativeHeights, LViewportContext parentContext, GEDevDesc *dd)
{
    int i;
    SEXP heights = layoutHeights(layout);
    double totalHeight = 0;
    /* We are calculating "null" units for a layout */
    L_nullLayoutMode = 1;
    for (i = 0; i < layoutNRow(layout); i++)
        if (relativeHeights[i])
            totalHeight += transformHeight(heights, i, parentContext, parentContext.fontfamily, parentContext.font,
                                           parentContext.fontsize, parentContext.lineheight,
                                           /*
                                            * NOTE: 0, 0, here is ok
                                            * because we are only
                                            * obtaining "null" units
                                            */
                                           0, 0, dd);
    L_nullLayoutMode = 0;
    return totalHeight;
}

void allocateRespected(SEXP layout, int *relativeWidths, int *relativeHeights, double hmult, double vmult,
                       double *reducedWidthCM, double *reducedHeightCM, LViewportContext parentContext, GEDevDesc *dd,
                       double *npcWidths, double *npcHeights)
{
    int i;
    SEXP widths = layoutWidths(layout);
    SEXP heights = layoutHeights(layout);
    int respect = layoutRespect(layout);
    double sumWidth = totalWidth(layout, relativeWidths, parentContext, dd);
    double sumHeight = totalHeight(layout, relativeHeights, parentContext, dd);
    double denom, mult;
    double tempWidthCM = *reducedWidthCM;
    double tempHeightCM = *reducedHeightCM;
    if (respect > 0)
    {
        /* Determine whether aspect ratio of available space is
         * bigger or smaller than aspect ratio of layout
         */
        if ((tempHeightCM / tempWidthCM) > (sumHeight / sumWidth))
        {
            denom = sumWidth;
            mult = tempWidthCM;
        }
        else
        {
            denom = sumHeight;
            mult = tempHeightCM;
        }
        /* Allocate respected widths
         */
        for (i = 0; i < layoutNCol(layout); i++)
            if (relativeWidths[i])
                if (colRespected(i, layout))
                {
                    /* Build a unit SEXP with a single value and no data
                     */
                    SEXP width;
                    PROTECT(width = unit(pureNullUnitValue(widths, i) / denom * mult, L_CM));
                    npcWidths[i] =
                        transformWidth(width, 0, parentContext, parentContext.fontfamily, parentContext.font,
                                       parentContext.fontsize, parentContext.lineheight, tempWidthCM, tempHeightCM, dd);
                    *reducedWidthCM -= npcWidths[i] * tempWidthCM;
                    npcWidths[i] = npcWidths[i] * hmult;
                    UNPROTECT(1);
                }
        /* Allocate respected heights
         */
        for (i = 0; i < layoutNRow(layout); i++)
            if (relativeHeights[i])
                if (rowRespected(i, layout))
                {
                    SEXP height;
                    PROTECT(height = unit(pureNullUnitValue(heights, i) / denom * mult, L_CM));
                    npcHeights[i] = transformHeight(height, 0, parentContext, parentContext.fontfamily,
                                                    parentContext.font, parentContext.fontsize,
                                                    parentContext.lineheight, tempWidthCM, tempHeightCM, dd);
                    *reducedHeightCM -= npcHeights[i] * tempHeightCM;
                    npcHeights[i] = npcHeights[i] * vmult;
                    UNPROTECT(1);
                }
    }
}

/* These sum up unrespected relative widths and heights (unit = "null")
 */
double totalUnrespectedWidth(SEXP layout, int *relativeWidths, LViewportContext parentContext, GEDevDesc *dd)
{
    int i;
    SEXP widths = layoutWidths(layout);
    double totalWidth = 0;
    /* We are calculating "null" units for a layout */
    L_nullLayoutMode = 1;
    for (i = 0; i < layoutNCol(layout); i++)
        if (relativeWidths[i])
            if (!colRespected(i, layout))
                totalWidth += transformWidth(widths, i, parentContext, parentContext.fontfamily, parentContext.font,
                                             parentContext.fontsize, parentContext.lineheight,
                                             /*
                                              * NOTE: 0, 0, here is ok
                                              * because we are only
                                              * obtaining "null" units
                                              */
                                             0, 0, dd);
    L_nullLayoutMode = 0;
    return totalWidth;
}

double totalUnrespectedHeight(SEXP layout, int *relativeHeights, LViewportContext parentContext, GEDevDesc *dd)
{
    int i;
    SEXP heights = layoutHeights(layout);
    double totalHeight = 0;
    /* We are calculating "null" units for a layout */
    L_nullLayoutMode = 1;
    for (i = 0; i < layoutNRow(layout); i++)
        if (relativeHeights[i])
            if (!rowRespected(i, layout))
                totalHeight += transformHeight(heights, i, parentContext, parentContext.fontfamily, parentContext.font,
                                               parentContext.fontsize, parentContext.lineheight,
                                               /*
                                                * NOTE: 0, 0, here is ok
                                                * because we are only
                                                * obtaining "null" units
                                                */
                                               0, 0, dd);
    L_nullLayoutMode = 0;
    return totalHeight;
}

void allocateRemainingWidth(SEXP layout, int *relativeWidths, double multiplier, LViewportContext parentContext,
                            GEDevDesc *dd, double *npcWidths)
{
    int i;
    SEXP widths = layoutWidths(layout);
    double sumWidth = totalUnrespectedWidth(layout, relativeWidths, parentContext, dd);
    /* We are calculating "null" units for a layout */
    L_nullLayoutMode = 1;
    for (i = 0; i < layoutNCol(layout); i++)
        if (relativeWidths[i])
            if (!colRespected(i, layout))
                npcWidths[i] = multiplier *
                               transformWidth(widths, i, parentContext, parentContext.fontfamily, parentContext.font,
                                              parentContext.fontsize, parentContext.lineheight,
                                              /*
                                               * NOTE: 0, 0, here is ok
                                               * because we are only
                                               * obtaining "null" units
                                               */
                                              0, 0, dd) /
                               sumWidth;
    L_nullLayoutMode = 0;
}

void allocateRemainingHeight(SEXP layout, int *relativeHeights, double multiplier, LViewportContext parentContext,
                             GEDevDesc *dd, double *npcHeights)
{
    int i;
    SEXP heights = layoutHeights(layout);
    double sumHeight = totalUnrespectedHeight(layout, relativeHeights, parentContext, dd);
    /* We are calculating "null" units for a layout */
    L_nullLayoutMode = 1;
    for (i = 0; i < layoutNRow(layout); i++)
        if (relativeHeights[i])
            if (!rowRespected(i, layout))
                npcHeights[i] = multiplier *
                                transformHeight(heights, i, parentContext, parentContext.fontfamily, parentContext.font,
                                                parentContext.fontsize, parentContext.lineheight,
                                                /*
                                                 * NOTE: 0, 0, here is ok
                                                 * because we are only
                                                 * obtaining "null" units
                                                 */
                                                0, 0, dd) /
                                sumHeight;
    L_nullLayoutMode = 0;
}

static double sumDims(double dims[], int from, int to)
{
    int i;
    double s = 0;
    for (i = from; i < to + 1; i++)
        s = s + dims[i];
    return s;
}

static void subRegion(SEXP layout, int minrow, int maxrow, int mincol, int maxcol, double widths[], double heights[],
                      double *left, double *bottom, double *width, double *height)
{
    double totalWidth = sumDims(widths, 0, layoutNCol(layout) - 1);
    double totalHeight = sumDims(heights, 0, layoutNRow(layout) - 1);
    *left = (0.5 - totalWidth / 2) + sumDims(widths, 0, mincol - 1);
    *width = sumDims(widths, mincol, maxcol);
    *bottom = (0.5 - totalHeight / 2) + totalHeight - sumDims(heights, 0, maxrow);
    *height = sumDims(heights, minrow, maxrow);
}

void calcViewportLayout(SEXP viewport, double parentWidthCM, double parentHeightCM, LViewportContext parentContext,
                        GEDevDesc *dd)
{
    int i;
    SEXP currentWidths, currentHeights;
    SEXP layout = viewportLayout(viewport);
    double *npcWidths = (double *)R_alloc(layoutNCol(layout), sizeof(double));
    double *npcHeights = (double *)R_alloc(layoutNRow(layout), sizeof(double));
    int *relativeWidths = (int *)R_alloc(layoutNCol(layout), sizeof(int));
    int *relativeHeights = (int *)R_alloc(layoutNRow(layout), sizeof(int));
    double reducedWidthCM = parentWidthCM;
    double reducedHeightCM = parentHeightCM;
    /* Figure out which rows and cols have relative heights and widths
     */
    findRelWidths(layout, relativeWidths);
    findRelHeights(layout, relativeHeights);
    /* For any width or height which has a unit other than "null"
     * we can immediately figure out its physical size and we can convert to
     * "npc" units.  We do this and return the widthCM and heightCM
     * remaining after these widths and heights have been allocated
     */
    allocateKnownWidths(layout, relativeWidths, parentWidthCM, parentHeightCM, parentContext, dd, npcWidths,
                        &reducedWidthCM);
    allocateKnownHeights(layout, relativeHeights, parentWidthCM, parentHeightCM, parentContext, dd, npcHeights,
                         &reducedHeightCM);
    /* Now allocate respected widths and heights and return
     * widthCM and heightCM remaining
     */
    allocateRespected(layout, relativeWidths, relativeHeights, reducedWidthCM / parentWidthCM,
                      reducedHeightCM / parentHeightCM, &reducedWidthCM, &reducedHeightCM, parentContext, dd, npcWidths,
                      npcHeights);
    /* Now allocate relative widths and heights (unit = "null")
     * in the remaining space
     */
    allocateRemainingWidth(layout, relativeWidths, reducedWidthCM / parentWidthCM, parentContext, dd, npcWidths);
    allocateRemainingHeight(layout, relativeHeights, reducedHeightCM / parentHeightCM, parentContext, dd, npcHeights);
    /* Record the widths and heights in the viewport
     */
    PROTECT(currentWidths = allocVector(REALSXP, layoutNCol(layout)));
    PROTECT(currentHeights = allocVector(REALSXP, layoutNRow(layout)));
    for (i = 0; i < layoutNCol(layout); i++)
    {
        REAL(currentWidths)[i] = npcWidths[i];
    }
    for (i = 0; i < layoutNRow(layout); i++)
    {
        REAL(currentHeights)[i] = npcHeights[i];
    }
    SET_VECTOR_ELT(viewport, PVP_WIDTHS, currentWidths);
    SET_VECTOR_ELT(viewport, PVP_HEIGHTS, currentHeights);
    UNPROTECT(2);
}

void calcViewportLocationFromLayout(SEXP layoutPosRow, SEXP layoutPosCol, SEXP parent, LViewportLocation *vpl)
{
    int minrow, maxrow, mincol, maxcol;
    double x, y, width, height;
    SEXP vpx, vpy, vpwidth, vpheight;
    SEXP layout = viewportLayout(parent);
    /* It is possible for ONE of layoutPosRow and layoutPosCol to
     * be NULL;  this is interpreted as "occupy all rows/cols"
     * NOTE: The " - 1" is there because R is 1-based and C is zero-based
     */
    if (isNull(layoutPosRow))
    {
        minrow = 0;
        maxrow = layoutNRow(layout) - 1;
    }
    else
    {
        minrow = INTEGER(layoutPosRow)[0] - 1;
        maxrow = INTEGER(layoutPosRow)[1] - 1;
    }
    if (isNull(layoutPosCol))
    {
        mincol = 0;
        maxcol = layoutNCol(layout) - 1;
    }
    else
    {
        mincol = INTEGER(layoutPosCol)[0] - 1;
        maxcol = INTEGER(layoutPosCol)[1] - 1;
    }
    /* Put the relevant values into vpl */
    subRegion(viewportLayout(parent), minrow, maxrow, mincol, maxcol, REAL(viewportLayoutWidths(parent)),
              REAL(viewportLayoutHeights(parent)), &x, &y, &width, &height);
    PROTECT(vpx = unit(x, L_NPC));
    vpl->x = vpx;
    PROTECT(vpy = unit(y, L_NPC));
    vpl->y = vpy;
    PROTECT(vpwidth = unit(width, L_NPC));
    vpl->width = vpwidth;
    PROTECT(vpheight = unit(height, L_NPC));
    vpl->height = vpheight;
    vpl->hjust = L_LEFT;
    vpl->vjust = L_BOTTOM;
    /* Question:  Is there any chance that these newly-allocated
     * unit SEXPs will get corrupted after this unprotect ??
     */
    UNPROTECT(4);
}
