/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001-3 Paul Murrell
 *                2003-5 The R Development Core Team
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

int isUnitArithmetic(SEXP ua)
{
    return inherits(ua, "unit.arithmetic");
}

int isUnitList(SEXP ul)
{
    return inherits(ul, "unit.list");
}

/* Function to build a single-value unit SEXP internally.
 * Cannot build units requiring data as yet.
 */
SEXP unit(double value, int unit)
{
    SEXP u, units, classname;
    PROTECT(u = allocVector(REALSXP, 1));
    REAL(u)[0] = value;
    PROTECT(units = allocVector(INTSXP, 1));
    INTEGER(units)[0] = unit;
    /* NOTE that we do not set the "unit" attribute */
    setAttrib(u, install("valid.unit"), units);
    setAttrib(u, install("data"), R_NilValue);
    PROTECT(classname = allocVector(STRSXP, 1));
    SET_STRING_ELT(classname, 0, mkChar("unit"));
    classgets(u, classname);
    UNPROTECT(3);
    return u;
}

/* Accessor functions for unit objects
 */

/*
 * This is an attempt to extract a single numeric value from
 * a unit.  This is ONLY designed for use on "simple" units
 * (i.e., NOT unitLists or unitArithmetics)
 */
double unitValue(SEXP unit, int index)
{
    /* Recycle values if necessary (used in unit arithmetic)
     */
    int n = LENGTH(unit);
    return numeric(unit, index % n);
}

int unitUnit(SEXP unit, int index)
{
    SEXP units = getAttrib(unit, install("valid.unit"));
    /* Recycle units if necessary
     */
    int n = LENGTH(units);
    return INTEGER(units)[index % n];
}

SEXP unitData(SEXP unit, int index)
{
    SEXP result;
    SEXP data = getAttrib(unit, install("data"));
    if (isNull(data))
        result = R_NilValue;
    else
    {
        /* Recycle data if necessary
         */
        int n = LENGTH(data);
        result = VECTOR_ELT(data, index % n);
    }
    return result;
}

/* Accessor functions for unit arithmetic object
 */
char *fName(SEXP ua)
{
    return CHAR(STRING_ELT(getListElement(ua, "fname"), 0));
}

SEXP arg1(SEXP ua)
{
    return getListElement(ua, "arg1");
}

SEXP arg2(SEXP ua)
{
    return getListElement(ua, "arg2");
}

int fNameMatch(SEXP ua, char *aString)
{
    return !strcmp(fName(ua), aString);
}

int addOp(SEXP ua)
{
    return fNameMatch(ua, "+");
}

int minusOp(SEXP ua)
{
    return fNameMatch(ua, "-");
}

int timesOp(SEXP ua)
{
    return fNameMatch(ua, "*");
}

int fOp(SEXP ua)
{
    return addOp(ua) || minusOp(ua) || timesOp(ua);
}

int minFunc(SEXP ua)
{
    return fNameMatch(ua, "min");
}

int maxFunc(SEXP ua)
{
    return fNameMatch(ua, "max");
}

int sumFunc(SEXP ua)
{
    return fNameMatch(ua, "sum");
}

/* Functions in lattice.c should use this to determine the length
 * of a unit/unitArithmetic object rather than just LENGTH.
 */
int unitLength(SEXP u)
{
    int result = 0;
    if (isUnitList(u))
        result = LENGTH(u);
    else if (isUnitArithmetic(u))
        if (fOp(u))
        {
            if (timesOp(u))
            {
                /*
                 * arg1 is always the numeric vector
                 */
                int n1 = LENGTH(arg1(u));
                int n2 = unitLength(arg2(u));
                result = (n1 > n2) ? n1 : n2;
            }
            else
            { /* must be "+" or "-" */
                int n1 = unitLength(arg1(u));
                int n2 = unitLength(arg2(u));
                result = (n1 > n2) ? n1 : n2;
            }
        }
        else            /* must be "min" or "max" or "sum" */
            result = 1; /* unitLength(arg1(u)); */
    else                /* Must be a unit object */
        result = LENGTH(u);
    return result;
}

/**************************
 * Code for handling "null" units
 **************************
 */

/* Global mode indicators:
 * The value returned for a "null" unit depends on ...
 * (i) whether layout is calling for evaluation of a "pure null" unit
 *     (in which case, the value of the "null" unit is returned)
 * (ii) the sort of arithmetic that is being performed
 *      (in which case, an "identity" value is returned)
 */

/*
 * Evaluate a "null" _value_ dependent on the evaluation context
 */
static double evaluateNullUnit(double value, double thisCM, int nullLayoutMode, int nullArithmeticMode)
{
    double result = value;
    if (!nullLayoutMode)
        switch (nullArithmeticMode)
        {
        case L_plain:
        case L_adding:
        case L_subtracting:
        case L_summing:
            result = 0;
            break;
        case L_multiplying:
            result = 0;
            break;
        case L_maximising:
            result = 0;
            break;
        case L_minimising:
            result = thisCM;
            break;
        }
    return result;
}

/*
 * Evaluate a "null" _unit_
 * This is used by layout code to get a single "null" _value_
 * from a pureNullUnit (which may be a unitList or a unitArithmetic)
 *
 * This must ONLY be called on a unit which has passed the
 * pureNullUnit test below.
 */
double pureNullUnitValue(SEXP unit, int index)
{
    double result = 0;
    if (isUnitArithmetic(unit))
    {
        int i;
        if (addOp(unit))
        {
            result = unitValue(arg1(unit), index) + unitValue(arg2(unit), index);
        }
        else if (minusOp(unit))
        {
            result = unitValue(arg1(unit), index) - unitValue(arg2(unit), index);
        }
        else if (timesOp(unit))
        {
            result = REAL(arg1(unit))[index] * unitValue(arg2(unit), index);
        }
        else if (minFunc(unit))
        {
            int n = unitLength(arg1(unit));
            double temp = DBL_MAX;
            result = unitValue(arg1(unit), 0);
            for (i = 1; i < n; i++)
            {
                temp = unitValue(arg1(unit), i);
                if (temp < result)
                    result = temp;
            }
        }
        else if (maxFunc(unit))
        {
            int n = unitLength(arg1(unit));
            double temp = DBL_MIN;
            result = unitValue(arg1(unit), 0);
            for (i = 1; i < n; i++)
            {
                temp = unitValue(arg1(unit), i);
                if (temp > result)
                    result = temp;
            }
        }
        else if (sumFunc(unit))
        {
            int n = unitLength(arg1(unit));
            result = 0.0;
            for (i = 0; i < n; i++)
            {
                result += unitValue(arg1(unit), i);
            }
        }
        else
            error(_("Unimplemented unit function"));
    }
    else if (isUnitList(unit))
    {
        result = unitValue(VECTOR_ELT(unit, index), 0);
    }
    else
        result = unitValue(unit, index);
    return result;
}

int pureNullUnitArithmetic(SEXP unit, int index, GEDevDesc *dd);

int pureNullUnit(SEXP unit, int index, GEDevDesc *dd)
{
    int result;
    if (isUnitArithmetic(unit))
        result = pureNullUnitArithmetic(unit, index, dd);
    else if (isUnitList(unit))
    {
        result = pureNullUnit(VECTOR_ELT(unit, index), 0, dd);
    }
    else
    { /* Just a plain unit */
        /* Special case:  if "grobwidth" or "grobheight" unit
         * and width/height(grob) is pure null
         */
        if (unitUnit(unit, index) == L_GROBWIDTH)
        {
            SEXP grob, width;
            SEXP widthPreFn, widthFn, widthPostFn, findGrobFn;
            SEXP R_fcall0, R_fcall1, R_fcall2, R_fcall3;
            SEXP savedgpar, savedgrob;
            /*
             * The data could be a gPath to a grob
             * In this case, need to find the grob first, and in order
             * to do that correctly, need to call pre/postDraw code
             */
            PROTECT(grob = unitData(unit, index));
            PROTECT(savedgpar = gridStateElement(dd, GSS_GPAR));
            PROTECT(savedgrob = gridStateElement(dd, GSS_CURRGROB));
            PROTECT(widthPreFn = findFun(install("preDraw"), R_gridEvalEnv));
            PROTECT(widthFn = findFun(install("width"), R_gridEvalEnv));
            PROTECT(widthPostFn = findFun(install("postDraw"), R_gridEvalEnv));
            if (inherits(grob, "gPath"))
            {
                if (isNull(savedgrob))
                {
                    PROTECT(findGrobFn = findFun(install("findGrobinDL"), R_gridEvalEnv));
                    PROTECT(R_fcall0 = lang2(findGrobFn, getListElement(grob, "name")));
                    grob = eval(R_fcall0, R_gridEvalEnv);
                }
                else
                {
                    PROTECT(findGrobFn = findFun(install("findGrobinChildren"), R_gridEvalEnv));
                    PROTECT(R_fcall0 =
                                lang3(findGrobFn, getListElement(grob, "name"), getListElement(savedgrob, "children")));
                    grob = eval(R_fcall0, R_gridEvalEnv);
                }
                UNPROTECT(2);
            }
            PROTECT(R_fcall1 = lang2(widthPreFn, grob));
            eval(R_fcall1, R_gridEvalEnv);
            PROTECT(R_fcall2 = lang2(widthFn, grob));
            PROTECT(width = eval(R_fcall2, R_gridEvalEnv));
            result = pureNullUnit(width, 0, dd);
            PROTECT(R_fcall3 = lang2(widthPostFn, grob));
            eval(R_fcall3, R_gridEvalEnv);
            setGridStateElement(dd, GSS_GPAR, savedgpar);
            setGridStateElement(dd, GSS_CURRGROB, savedgrob);
            UNPROTECT(10);
        }
        else if (unitUnit(unit, index) == L_GROBHEIGHT)
        {
            SEXP grob, height;
            SEXP heightPreFn, heightFn, heightPostFn, findGrobFn;
            SEXP R_fcall0, R_fcall1, R_fcall2, R_fcall3;
            SEXP savedgpar, savedgrob;
            /*
             * The data could be a gPath to a grob
             * In this case, need to find the grob first, and in order
             * to do that correctly, need to call pre/postDraw code
             */
            PROTECT(grob = unitData(unit, index));
            PROTECT(savedgpar = gridStateElement(dd, GSS_GPAR));
            PROTECT(savedgrob = gridStateElement(dd, GSS_CURRGROB));
            PROTECT(heightPreFn = findFun(install("preDraw"), R_gridEvalEnv));
            PROTECT(heightFn = findFun(install("height"), R_gridEvalEnv));
            PROTECT(heightPostFn = findFun(install("postDraw"), R_gridEvalEnv));
            if (inherits(grob, "gPath"))
            {
                if (isNull(savedgrob))
                {
                    PROTECT(findGrobFn = findFun(install("findGrobinDL"), R_gridEvalEnv));
                    PROTECT(R_fcall0 = lang2(findGrobFn, getListElement(grob, "name")));
                    grob = eval(R_fcall0, R_gridEvalEnv);
                }
                else
                {
                    PROTECT(findGrobFn = findFun(install("findGrobinChildren"), R_gridEvalEnv));
                    PROTECT(R_fcall0 =
                                lang3(findGrobFn, getListElement(grob, "name"), getListElement(savedgrob, "children")));
                    grob = eval(R_fcall0, R_gridEvalEnv);
                }
                UNPROTECT(2);
            }
            PROTECT(R_fcall1 = lang2(heightPreFn, grob));
            eval(R_fcall1, R_gridEvalEnv);
            PROTECT(R_fcall2 = lang2(heightFn, grob));
            PROTECT(height = eval(R_fcall2, R_gridEvalEnv));
            result = pureNullUnit(height, 0, dd);
            PROTECT(R_fcall3 = lang2(heightPostFn, grob));
            eval(R_fcall3, R_gridEvalEnv);
            setGridStateElement(dd, GSS_GPAR, savedgpar);
            setGridStateElement(dd, GSS_CURRGROB, savedgrob);
            UNPROTECT(10);
        }
        else
            result = unitUnit(unit, index) == L_NULL;
    }
    return result;
}

int pureNullUnitArithmetic(SEXP unit, int index, GEDevDesc *dd)
{
    /*
     * Initialised to shut up compiler
     */
    int result = 0;
    if (addOp(unit) || minusOp(unit))
    {
        result = pureNullUnit(arg1(unit), index, dd) && pureNullUnit(arg2(unit), index, dd);
    }
    else if (timesOp(unit))
    {
        result = pureNullUnit(arg2(unit), index, dd);
    }
    else if (minFunc(unit) || maxFunc(unit) || sumFunc(unit))
    {
        int n = unitLength(arg1(unit));
        int i = 0;
        result = 1;
        while (result && i < n)
        {
            result = result && pureNullUnit(arg1(unit), i, dd);
            i += 1;
        }
    }
    else
        error(_("Unimplemented unit function"));
    return result;
}

/**************************
 * Code for handling "grobwidth" units
 **************************
 */

/* NOTE:  this code calls back to R code to perform
 * set.gpar operations, which will impact on grid state variables
 * BUT that's ok(ish) because we save and restore the relevant state
 * variables in here so that the overall effect is NULL.
 */

double evaluateGrobWidthUnit(SEXP grob, double vpheightCM, double vpwidthCM, int nullLMode, int nullAMode,
                             GEDevDesc *dd)
{
    double vpWidthCM, vpHeightCM;
    double rotationAngle;
    LViewportContext vpc;
    R_GE_gcontext gc;
    LTransform transform;
    SEXP currentvp, currentgp;
    SEXP widthPreFn, widthFn, widthPostFn, findGrobFn;
    SEXP R_fcall0, R_fcall1, R_fcall2, R_fcall3;
    SEXP savedgpar, savedgrob;
    SEXP width;
    double result;
    /*
     * Save the current gpar state and restore it at the end
     */
    PROTECT(savedgpar = gridStateElement(dd, GSS_GPAR));
    /*
     * Save the current grob and restore it at the end
     */
    PROTECT(savedgrob = gridStateElement(dd, GSS_CURRGROB));
    /*
     * Set up for calling R functions
     */
    PROTECT(widthPreFn = findFun(install("preDraw"), R_gridEvalEnv));
    PROTECT(widthFn = findFun(install("width"), R_gridEvalEnv));
    PROTECT(widthPostFn = findFun(install("postDraw"), R_gridEvalEnv));
    /*
     * If grob is actually a gPath, use it to find an actual grob
     */
    if (inherits(grob, "gPath"))
    {
        /*
         * If the current grob is NULL then we are at the top level
         * and we search the display list, otherwise we search the
         * children of the current grob
         *
         * NOTE: assume here that only gPath of depth == 1 are valid
         */
        if (isNull(savedgrob))
        {
            PROTECT(findGrobFn = findFun(install("findGrobinDL"), R_gridEvalEnv));
            PROTECT(R_fcall0 = lang2(findGrobFn, getListElement(grob, "name")));
            grob = eval(R_fcall0, R_gridEvalEnv);
        }
        else
        {
            PROTECT(findGrobFn = findFun(install("findGrobinChildren"), R_gridEvalEnv));
            PROTECT(R_fcall0 = lang3(findGrobFn, getListElement(grob, "name"), getListElement(savedgrob, "children")));
            grob = eval(R_fcall0, R_gridEvalEnv);
        }
        UNPROTECT(2);
    }
    /* Call preDraw(grob)
     */
    PROTECT(R_fcall1 = lang2(widthPreFn, grob));
    eval(R_fcall1, R_gridEvalEnv);
    /*
     * The call to preDraw may have pushed viewports and/or
     * enforced gpar settings, SO we need to re-establish the
     * current viewport and gpar settings before evaluating the
     * width unit.
     *
     * NOTE:  we are really relying on the grid state to be coherent
     * when we do stuff like this (i.e., not to have changed since
     * we started evaluating the unit [other than the changes we may
     * have deliberately made above by calling preDraw]).  In other
     * words we are relying on no other drawing occurring at the
     * same time as we are doing this evaluation.  In other other
     * words, we are relying on there being only ONE process
     * (i.e., NOT multi-threaded).
     */
    currentvp = gridStateElement(dd, GSS_VP);
    currentgp = gridStateElement(dd, GSS_GPAR);
    getViewportTransform(currentvp, dd, &vpWidthCM, &vpHeightCM, transform, &rotationAngle);
    fillViewportContextFromViewport(currentvp, &vpc);
    /* Call width(grob)
     * to get the unit representing the with
     */
    PROTECT(R_fcall2 = lang2(widthFn, grob));
    PROTECT(width = eval(R_fcall2, R_gridEvalEnv));
    /*
     * Transform the width
     * NOTE:  We transform into INCHES so can produce final answer in terms
     * of NPC for original context
     */
    /* Special case for "null" units
     */
    if (pureNullUnit(width, 0, dd))
    {
        result = evaluateNullUnit(pureNullUnitValue(width, 0), vpWidthCM, nullLMode, nullAMode);
    }
    else
    {
        gcontextFromgpar(currentgp, 0, &gc);
        result = transformWidthtoINCHES(width, 0, vpc, &gc, vpWidthCM, vpHeightCM, dd);
    }
    /* Call postDraw(grob)
     */
    PROTECT(R_fcall3 = lang2(widthPostFn, grob));
    eval(R_fcall3, R_gridEvalEnv);
    /*
     * Restore the saved gpar state and grob
     */
    setGridStateElement(dd, GSS_GPAR, savedgpar);
    setGridStateElement(dd, GSS_CURRGROB, savedgrob);
    UNPROTECT(9);
    /* Return the transformed width
     */
    return result;
}

/* See evaluateGrobWidthUnit for detailed comments
 */
double evaluateGrobHeightUnit(SEXP grob, double vpheightCM, double vpwidthCM, int nullLMode, int nullAMode,
                              GEDevDesc *dd)
{
    double vpWidthCM, vpHeightCM;
    double rotationAngle;
    LViewportContext vpc;
    R_GE_gcontext gc;
    LTransform transform;
    SEXP currentvp, currentgp;
    SEXP heightPreFn, heightFn, heightPostFn, findGrobFn;
    SEXP R_fcall0, R_fcall1, R_fcall2, R_fcall3;
    SEXP savedgpar, savedgrob;
    SEXP height;
    double result;
    /*
     * Save the current gpar state and restore it at the end
     */
    PROTECT(savedgpar = gridStateElement(dd, GSS_GPAR));
    /*
     * Save the current grob and restore it at the end
     */
    PROTECT(savedgrob = gridStateElement(dd, GSS_CURRGROB));
    PROTECT(heightPreFn = findFun(install("preDraw"), R_gridEvalEnv));
    PROTECT(heightFn = findFun(install("height"), R_gridEvalEnv));
    PROTECT(heightPostFn = findFun(install("postDraw"), R_gridEvalEnv));
    /*
     * If grob is actually a gPath, use it to find an actual grob
     */
    if (inherits(grob, "gPath"))
    {
        /*
         * If the current grob is NULL then we are at the top level
         * and we search the display list, otherwise we search the
         * children of the current grob
         *
         * NOTE: assume here that only gPath of depth == 1 are valid
         */
        if (isNull(savedgrob))
        {
            PROTECT(findGrobFn = findFun(install("findGrobinDL"), R_gridEvalEnv));
            PROTECT(R_fcall0 = lang2(findGrobFn, getListElement(grob, "name")));
            grob = eval(R_fcall0, R_gridEvalEnv);
        }
        else
        {
            PROTECT(findGrobFn = findFun(install("findGrobinChildren"), R_gridEvalEnv));
            PROTECT(R_fcall0 = lang3(findGrobFn, getListElement(grob, "name"), getListElement(savedgrob, "children")));
            grob = eval(R_fcall0, R_gridEvalEnv);
        }
        UNPROTECT(2);
    }
    /* Call preDraw(grob)
     */
    PROTECT(R_fcall1 = lang2(heightPreFn, grob));
    eval(R_fcall1, R_gridEvalEnv);
    /*
     * The call to preDraw may have pushed viewports and/or
     * enforced gpar settings, SO we need to re-establish the
     * current viewport and gpar settings before evaluating the
     * height unit.
     *
     * NOTE:  we are really relying on the grid state to be coherent
     * when we do stuff like this (i.e., not to have changed since
     * we started evaluating the unit [other than the changes we may
     * have deliberately made above by calling preDraw]).  In other
     * words we are relying on no other drawing occurring at the
     * same time as we are doing this evaluation.  In other other
     * words, we are relying on there being only ONE process
     * (i.e., NOT multi-threaded).
     */
    currentvp = gridStateElement(dd, GSS_VP);
    currentgp = gridStateElement(dd, GSS_GPAR);
    getViewportTransform(currentvp, dd, &vpWidthCM, &vpHeightCM, transform, &rotationAngle);
    fillViewportContextFromViewport(currentvp, &vpc);
    /* Call height(grob)
     * to get the unit representing the with
     */
    PROTECT(R_fcall2 = lang2(heightFn, grob));
    PROTECT(height = eval(R_fcall2, R_gridEvalEnv));
    /*
     * Transform the height
     * NOTE:  We transform into INCHES so can produce final answer in terms
     * of NPC for original context
     */
    /* Special case for "null" units
     */
    if (pureNullUnit(height, 0, dd))
    {
        result = evaluateNullUnit(pureNullUnitValue(height, 0), vpHeightCM, nullLMode, nullAMode);
    }
    else
    {
        gcontextFromgpar(currentgp, 0, &gc);
        result = transformHeighttoINCHES(height, 0, vpc, &gc, vpWidthCM, vpHeightCM, dd);
    }
    /* Call postDraw(grob)
     */
    PROTECT(R_fcall3 = lang2(heightPostFn, grob));
    eval(R_fcall3, R_gridEvalEnv);
    /*
     * Restore the saved gpar state and grob
     */
    setGridStateElement(dd, GSS_GPAR, savedgpar);
    setGridStateElement(dd, GSS_CURRGROB, savedgrob);
    UNPROTECT(9);
    /* Return the transformed height
     */
    return result;
}

/**************************
 * TRANSFORMATIONS
 **************************
 */

/* Map a value from arbitrary units to INCHES */

/*
 * NULL units are a special case
 * If L_nullLayoutMode = 1 then the value returned is a NULL unit value
 * Otherwise it is an INCHES value
 */
double transform(double value, int unit, SEXP data, double scalemin, double scalemax, R_GE_gcontext *gc, double thisCM,
                 double otherCM, int nullLMode, int nullAMode, GEDevDesc *dd)
{
    double result = value;
    switch (unit)
    {
    case L_NPC:
        result = (result * thisCM) / 2.54; /* 2.54 cm per inch */
        break;
    case L_CM:
        result = result / 2.54;
        break;
    case L_INCHES:
        break;
    /* FIXME:  The following two assume that the pointsize specified
     * by the user is actually the pointsize provided by the
     * device.  This is NOT a safe assumption
     * One possibility would be to do a call to GReset(), just so
     * that mapping() gets called, just so that things like
     * xNDCPerLine are up-to-date, THEN call GStrHeight("M")
     * or somesuch.
     */
    case L_CHAR:
    case L_MYCHAR:                                 /* FIXME: Remove this when I can */
        result = (result * gc->ps * gc->cex) / 72; /* 72 points per inch */
        break;
    case L_LINES:
    case L_MYLINES: /* FIXME: Remove this when I can */
        result = (result * gc->ps * gc->cex * gc->lineheight) / 72;
        break;
    case L_SNPC:
        if (thisCM <= otherCM)
            result = (result * thisCM) / 2.54;
        else
            result = (result * otherCM) / 2.54;
        break;
    case L_MM:
        result = (result / 10) / 2.54;
        break;
    /* Maybe an opportunity for some constants below here (!)
     */
    case L_POINTS:
        result = result / 72.27;
        break;
    case L_PICAS:
        result = (result * 12) / 72.27;
        break;
    case L_BIGPOINTS:
        result = result / 72;
        break;
    case L_DIDA:
        result = result / 1157 * 1238 / 72.27;
        break;
    case L_CICERO:
        result = result * 12 / 1157 * 1238 / 72.27;
        break;
    case L_SCALEDPOINTS:
        result = result / 65536 / 72.27;
        break;
    case L_STRINGWIDTH:
    case L_MYSTRINGWIDTH: /* FIXME: Remove this when I can */
        if (isExpression(data))
            result = result * fromDeviceWidth(GEExpressionWidth(VECTOR_ELT(data, 0), gc, dd), GE_INCHES, dd);
        else
            result = result * fromDeviceWidth(GEStrWidth(CHAR(STRING_ELT(data, 0)), gc, dd), GE_INCHES, dd);
        break;
    case L_STRINGHEIGHT:
    case L_MYSTRINGHEIGHT: /* FIXME: Remove this when I can */
        if (isExpression(data))
            result = result * fromDeviceHeight(GEExpressionHeight(VECTOR_ELT(data, 0), gc, dd), GE_INCHES, dd);
        else
            result = result * fromDeviceHeight(GEStrHeight(CHAR(STRING_ELT(data, 0)), gc, dd), GE_INCHES, dd);
        break;
    case L_GROBWIDTH:
        result = value * evaluateGrobWidthUnit(data, otherCM, thisCM, nullLMode, nullAMode, dd);
        break;
    case L_GROBHEIGHT:
        result = value * evaluateGrobHeightUnit(data, thisCM, otherCM, nullLMode, nullAMode, dd);
        break;
    case L_NULL:
        result = evaluateNullUnit(result, thisCM, nullLMode, nullAMode);
        break;
    default:
        error(_("Illegal unit or unit not yet implemented"));
    }
    return result;
}

/* FIXME:  scales are only linear at the moment */
double transformLocation(double location, int unit, SEXP data, double scalemin, double scalemax, R_GE_gcontext *gc,
                         double thisCM, double otherCM, int nullLMode, int nullAMode, GEDevDesc *dd)
{
    double result = location;
    switch (unit)
    {
    case L_NATIVE:
        result = ((result - scalemin) / (scalemax - scalemin)) * thisCM / 2.54;
        break;
    default:
        result = transform(location, unit, data, scalemin, scalemax, gc, thisCM, otherCM, nullLMode, nullAMode, dd);
    }
    return result;
}

double transformXArithmetic(SEXP x, int index, LViewportContext vpc, R_GE_gcontext *gc, double widthCM, double heightCM,
                            int nullLMode, GEDevDesc *dd);

double transformX(SEXP x, int index, LViewportContext vpc, R_GE_gcontext *gc, double widthCM, double heightCM,
                  int nullLMode, int nullAMode, GEDevDesc *dd)
{
    double result;
    int unit;
    SEXP data;
    if (isUnitArithmetic(x))
        result = transformXArithmetic(x, index, vpc, gc, widthCM, heightCM, nullLMode, dd);
    else if (isUnitList(x))
    {
        int n = unitLength(x);
        result = transformX(VECTOR_ELT(x, index % n), 0, vpc, gc, widthCM, heightCM, nullLMode, nullAMode, dd);
    }
    else
    { /* Just a plain unit */
        int nullamode;
        if (nullAMode == 0)
            nullamode = L_plain;
        else
            nullamode = nullAMode;
        result = unitValue(x, index);
        unit = unitUnit(x, index);
        data = unitData(x, index);
        result = transformLocation(result, unit, data, vpc.xscalemin, vpc.xscalemax, gc, widthCM, heightCM, nullLMode,
                                   nullamode, dd);
    }
    return result;
}

double transformYArithmetic(SEXP y, int index, LViewportContext vpc, R_GE_gcontext *gc, double widthCM, double heightCM,
                            int nullLMode, GEDevDesc *dd);

double transformY(SEXP y, int index, LViewportContext vpc, R_GE_gcontext *gc, double widthCM, double heightCM,
                  int nullLMode, int nullAMode, GEDevDesc *dd)
{
    double result;
    int unit;
    SEXP data;
    if (isUnitArithmetic(y))
        result = transformYArithmetic(y, index, vpc, gc, widthCM, heightCM, nullLMode, dd);
    else if (isUnitList(y))
    {
        int n = unitLength(y);
        result = transformY(VECTOR_ELT(y, index % n), 0, vpc, gc, widthCM, heightCM, nullLMode, nullAMode, dd);
    }
    else
    { /* Just a unit object */
        int nullamode;
        if (nullAMode == 0)
            nullamode = L_plain;
        else
            nullamode = nullAMode;
        result = unitValue(y, index);
        unit = unitUnit(y, index);
        data = unitData(y, index);
        result = transformLocation(result, unit, data, vpc.yscalemin, vpc.yscalemax, gc, heightCM, widthCM, nullLMode,
                                   nullamode, dd);
    }
    return result;
}

double transformDimension(double dim, int unit, SEXP data, double scalemin, double scalemax, R_GE_gcontext *gc,
                          double thisCM, double otherCM, int nullLMode, int nullAMode, GEDevDesc *dd)
{
    double result = dim;
    switch (unit)
    {
    case L_NATIVE:
        result = ((dim) / (scalemax - scalemin)) * thisCM / 2.54;
        break;
    default:
        result = transform(dim, unit, data, scalemin, scalemax, gc, thisCM, otherCM, nullLMode, nullAMode, dd);
    }
    return result;
}

double transformWidthArithmetic(SEXP width, int index, LViewportContext vpc, R_GE_gcontext *gc, double widthCM,
                                double heightCM, int nullLMode, GEDevDesc *dd);

double transformWidth(SEXP width, int index, LViewportContext vpc, R_GE_gcontext *gc, double widthCM, double heightCM,
                      int nullLMode, int nullAMode, GEDevDesc *dd)
{
    double result;
    int unit;
    SEXP data;
    if (isUnitArithmetic(width))
        result = transformWidthArithmetic(width, index, vpc, gc, widthCM, heightCM, nullLMode, dd);
    else if (isUnitList(width))
    {
        int n = unitLength(width);
        result = transformWidth(VECTOR_ELT(width, index % n), 0, vpc, gc, widthCM, heightCM, nullLMode, nullAMode, dd);
    }
    else
    { /* Just a unit object */
        int nullamode;
        if (nullAMode == 0)
            nullamode = L_plain;
        else
            nullamode = nullAMode;
        result = unitValue(width, index);
        unit = unitUnit(width, index);
        data = unitData(width, index);
        result = transformDimension(result, unit, data, vpc.xscalemin, vpc.xscalemax, gc, widthCM, heightCM, nullLMode,
                                    nullamode, dd);
    }
    return result;
}

double transformHeightArithmetic(SEXP height, int index, LViewportContext vpc, R_GE_gcontext *gc, double widthCM,
                                 double heightCM, int nullLMode, GEDevDesc *dd);

double transformHeight(SEXP height, int index, LViewportContext vpc, R_GE_gcontext *gc, double widthCM, double heightCM,
                       int nullLMode, int nullAMode, GEDevDesc *dd)
{
    double result;
    int unit;
    SEXP data;
    if (isUnitArithmetic(height))
        result = transformHeightArithmetic(height, index, vpc, gc, widthCM, heightCM, nullLMode, dd);
    else if (isUnitList(height))
    {
        int n = unitLength(height);
        result =
            transformHeight(VECTOR_ELT(height, index % n), 0, vpc, gc, widthCM, heightCM, nullLMode, nullAMode, dd);
    }
    else
    { /* Just a unit object */
        int nullamode;
        if (nullAMode == 0)
            nullamode = L_plain;
        else
            nullamode = nullAMode;
        result = unitValue(height, index);
        unit = unitUnit(height, index);
        data = unitData(height, index);
        result = transformDimension(result, unit, data, vpc.yscalemin, vpc.yscalemax, gc, heightCM, widthCM, nullLMode,
                                    nullamode, dd);
    }
    return result;
}

double transformXArithmetic(SEXP x, int index, LViewportContext vpc, R_GE_gcontext *gc, double widthCM, double heightCM,
                            int nullLMode, GEDevDesc *dd)
{
    int i;
    double result = 0;
    if (addOp(x))
    {
        result = transformX(arg1(x), index, vpc, gc, widthCM, heightCM, nullLMode, L_adding, dd) +
                 transformX(arg2(x), index, vpc, gc, widthCM, heightCM, nullLMode, L_adding, dd);
    }
    else if (minusOp(x))
    {
        result = transformX(arg1(x), index, vpc, gc, widthCM, heightCM, nullLMode, L_subtracting, dd) -
                 transformX(arg2(x), index, vpc, gc, widthCM, heightCM, nullLMode, L_subtracting, dd);
    }
    else if (timesOp(x))
    {
        result = REAL(arg1(x))[index % LENGTH(arg1(x))] *
                 transformX(arg2(x), index, vpc, gc, widthCM, heightCM, nullLMode, L_multiplying, dd);
    }
    else if (minFunc(x))
    {
        int n = unitLength(arg1(x));
        double temp = DBL_MAX;
        result = transformX(arg1(x), 0, vpc, gc, widthCM, heightCM, nullLMode, L_minimising, dd);
        for (i = 1; i < n; i++)
        {
            temp = transformX(arg1(x), i, vpc, gc, widthCM, heightCM, nullLMode, L_minimising, dd);
            if (temp < result)
                result = temp;
        }
    }
    else if (maxFunc(x))
    {
        int n = unitLength(arg1(x));
        double temp = DBL_MIN;
        result = transformX(arg1(x), 0, vpc, gc, widthCM, heightCM, nullLMode, L_maximising, dd);
        for (i = 1; i < n; i++)
        {
            temp = transformX(arg1(x), i, vpc, gc, widthCM, heightCM, nullLMode, L_maximising, dd);
            if (temp > result)
                result = temp;
        }
    }
    else if (sumFunc(x))
    {
        int n = unitLength(arg1(x));
        result = 0.0;
        for (i = 0; i < n; i++)
        {
            result += transformX(arg1(x), i, vpc, gc, widthCM, heightCM, nullLMode, L_summing, dd);
        }
    }
    else
        error(_("Unimplemented unit function"));
    return result;
}

double transformYArithmetic(SEXP y, int index, LViewportContext vpc, R_GE_gcontext *gc, double widthCM, double heightCM,
                            int nullLMode, GEDevDesc *dd)
{
    int i;
    double result = 0;
    if (addOp(y))
    {
        result = transformY(arg1(y), index, vpc, gc, widthCM, heightCM, nullLMode, L_adding, dd) +
                 transformY(arg2(y), index, vpc, gc, widthCM, heightCM, nullLMode, L_adding, dd);
    }
    else if (minusOp(y))
    {
        result = transformY(arg1(y), index, vpc, gc, widthCM, heightCM, nullLMode, L_subtracting, dd) -
                 transformY(arg2(y), index, vpc, gc, widthCM, heightCM, nullLMode, L_subtracting, dd);
    }
    else if (timesOp(y))
    {
        result = REAL(arg1(y))[index % LENGTH(arg1(y))] *
                 transformY(arg2(y), index, vpc, gc, widthCM, heightCM, nullLMode, L_multiplying, dd);
    }
    else if (minFunc(y))
    {
        int n = unitLength(arg1(y));
        double temp = DBL_MAX;
        result = transformY(arg1(y), 0, vpc, gc, widthCM, heightCM, nullLMode, L_minimising, dd);
        for (i = 1; i < n; i++)
        {
            temp = transformY(arg1(y), i, vpc, gc, widthCM, heightCM, nullLMode, L_minimising, dd);
            if (temp < result)
                result = temp;
        }
    }
    else if (maxFunc(y))
    {
        int n = unitLength(arg1(y));
        double temp = DBL_MIN;
        result = transformY(arg1(y), 0, vpc, gc, widthCM, heightCM, nullLMode, L_maximising, dd);
        for (i = 1; i < n; i++)
        {
            temp = transformY(arg1(y), i, vpc, gc, widthCM, heightCM, nullLMode, L_maximising, dd);
            if (temp > result)
                result = temp;
        }
    }
    else if (sumFunc(y))
    {
        int n = unitLength(arg1(y));
        result = 0.0;
        for (i = 0; i < n; i++)
        {
            result += transformY(arg1(y), i, vpc, gc, widthCM, heightCM, nullLMode, L_summing, dd);
        }
    }
    else
        error(_("Unimplemented unit function"));
    return result;
}

double transformWidthArithmetic(SEXP width, int index, LViewportContext vpc, R_GE_gcontext *gc, double widthCM,
                                double heightCM, int nullLMode, GEDevDesc *dd)
{
    int i;
    double result = 0;
    if (addOp(width))
    {
        result = transformWidth(arg1(width), index, vpc, gc, widthCM, heightCM, nullLMode, L_adding, dd) +
                 transformWidth(arg2(width), index, vpc, gc, widthCM, heightCM, nullLMode, L_adding, dd);
    }
    else if (minusOp(width))
    {
        result = transformWidth(arg1(width), index, vpc, gc, widthCM, heightCM, nullLMode, L_subtracting, dd) -
                 transformWidth(arg2(width), index, vpc, gc, widthCM, heightCM, nullLMode, L_subtracting, dd);
    }
    else if (timesOp(width))
    {
        result = REAL(arg1(width))[index % LENGTH(arg1(width))] *
                 transformWidth(arg2(width), index, vpc, gc, widthCM, heightCM, nullLMode, L_multiplying, dd);
    }
    else if (minFunc(width))
    {
        int n = unitLength(arg1(width));
        double temp = DBL_MAX;
        result = transformWidth(arg1(width), 0, vpc, gc, widthCM, heightCM, nullLMode, L_minimising, dd);
        for (i = 1; i < n; i++)
        {
            temp = transformWidth(arg1(width), i, vpc, gc, widthCM, heightCM, nullLMode, L_minimising, dd);
            if (temp < result)
                result = temp;
        }
    }
    else if (maxFunc(width))
    {
        int n = unitLength(arg1(width));
        double temp = DBL_MIN;
        result = transformWidth(arg1(width), 0, vpc, gc, widthCM, heightCM, nullLMode, L_maximising, dd);
        for (i = 1; i < n; i++)
        {
            temp = transformWidth(arg1(width), i, vpc, gc, widthCM, heightCM, nullLMode, L_maximising, dd);
            if (temp > result)
                result = temp;
        }
    }
    else if (sumFunc(width))
    {
        int n = unitLength(arg1(width));
        result = 0.0;
        for (i = 0; i < n; i++)
        {
            result += transformWidth(arg1(width), i, vpc, gc, widthCM, heightCM, nullLMode, L_summing, dd);
        }
    }
    else
        error(_("Unimplemented unit function"));
    return result;
}

double transformHeightArithmetic(SEXP height, int index, LViewportContext vpc, R_GE_gcontext *gc, double widthCM,
                                 double heightCM, int nullLMode, GEDevDesc *dd)
{
    int i;
    double result = 0;
    if (addOp(height))
    {
        result = transformHeight(arg1(height), index, vpc, gc, widthCM, heightCM, nullLMode, L_adding, dd) +
                 transformHeight(arg2(height), index, vpc, gc, widthCM, heightCM, nullLMode, L_adding, dd);
    }
    else if (minusOp(height))
    {
        result = transformHeight(arg1(height), index, vpc, gc, widthCM, heightCM, nullLMode, L_subtracting, dd) -
                 transformHeight(arg2(height), index, vpc, gc, widthCM, heightCM, nullLMode, L_subtracting, dd);
    }
    else if (timesOp(height))
    {
        result = REAL(arg1(height))[index % LENGTH(arg1(height))] *
                 transformHeight(arg2(height), index, vpc, gc, widthCM, heightCM, nullLMode, L_multiplying, dd);
    }
    else if (minFunc(height))
    {
        int n = unitLength(arg1(height));
        double temp = DBL_MAX;
        result = transformHeight(arg1(height), 0, vpc, gc, widthCM, heightCM, nullLMode, L_minimising, dd);
        for (i = 1; i < n; i++)
        {
            temp = transformHeight(arg1(height), i, vpc, gc, widthCM, heightCM, nullLMode, L_minimising, dd);
            if (temp < result)
                result = temp;
        }
    }
    else if (maxFunc(height))
    {
        int n = unitLength(arg1(height));
        double temp = DBL_MIN;
        result = transformHeight(arg1(height), 0, vpc, gc, widthCM, heightCM, nullLMode, L_maximising, dd);
        for (i = 1; i < n; i++)
        {
            temp = transformHeight(arg1(height), i, vpc, gc, widthCM, heightCM, nullLMode, L_maximising, dd);
            if (temp > result)
                result = temp;
        }
    }
    else if (sumFunc(height))
    {
        int n = unitLength(arg1(height));
        result = 0.0;
        for (i = 0; i < n; i++)
        {
            result += transformHeight(arg1(height), i, vpc, gc, widthCM, heightCM, nullLMode, L_summing, dd);
        }
    }
    else
        error(_("Unimplemented unit function"));
    return result;
}

/* Code for transforming a location in INCHES using a transformation matrix.
 * We work in INCHES so that rotations can be incorporated within the
 * transformation matrix (i.e., the units are the same in both x- and
 * y-directions).
 * INCHES rather than CM because the R graphics engine only has INCHES.
 */

/* The original transform[X | Y | Width | Height] functions
 * were written to transform to NPC.  Rather than muck with them,
 * I am just wrappering them to get the new transformation to INCHES
 * In other words, the reason for the apparent inefficiency here
 * is historical.
 */

/* It is even more inefficient-looking now because I ended up mucking
 * with transform() to return INCHES (to fix bug if width/heightCM == 0)
 * and by then there was too much code that called transformXtoINCHES
 * to be bothered changing calls to it
 */

/* The difference between transform*toINCHES and transformLocn/Dimn
 * is that the former are just converting from one coordinate system
 * to INCHES;  the latter are converting from INCHES relative to
 * the parent to INCHES relative to the device.
 */
double transformXtoINCHES(SEXP x, int index, LViewportContext vpc, R_GE_gcontext *gc, double widthCM, double heightCM,
                          GEDevDesc *dd)
{
    return transformX(x, index, vpc, gc, widthCM, heightCM, 0, 0, dd);
}

double transformYtoINCHES(SEXP y, int index, LViewportContext vpc, R_GE_gcontext *gc, double widthCM, double heightCM,
                          GEDevDesc *dd)
{
    return transformY(y, index, vpc, gc, widthCM, heightCM, 0, 0, dd);
}

void transformLocn(SEXP x, SEXP y, int index, LViewportContext vpc, R_GE_gcontext *gc, double widthCM, double heightCM,
                   GEDevDesc *dd, LTransform t, double *xx, double *yy)
{
    LLocation lin, lout;
    /* x and y are unit objects (i.e., values in any old coordinate
     * system) so the first step is to convert them both to CM
     */
    *xx = transformXtoINCHES(x, index, vpc, gc, widthCM, heightCM, dd);
    *yy = transformYtoINCHES(y, index, vpc, gc, widthCM, heightCM, dd);
    location(*xx, *yy, lin);
    trans(lin, t, lout);
    *xx = locationX(lout);
    *yy = locationY(lout);
}

double transformWidthtoINCHES(SEXP w, int index, LViewportContext vpc, R_GE_gcontext *gc, double widthCM,
                              double heightCM, GEDevDesc *dd)
{
    return transformWidth(w, index, vpc, gc, widthCM, heightCM, 0, 0, dd);
}

double transformHeighttoINCHES(SEXP h, int index, LViewportContext vpc, R_GE_gcontext *gc, double widthCM,
                               double heightCM, GEDevDesc *dd)
{
    return transformHeight(h, index, vpc, gc, widthCM, heightCM, 0, 0, dd);
}

void transformDimn(SEXP w, SEXP h, int index, LViewportContext vpc, R_GE_gcontext *gc, double widthCM, double heightCM,
                   GEDevDesc *dd, double rotationAngle, double *ww, double *hh)
{
    LLocation din, dout;
    LTransform r;
    *ww = transformWidthtoINCHES(w, index, vpc, gc, widthCM, heightCM, dd);
    *hh = transformHeighttoINCHES(h, index, vpc, gc, widthCM, heightCM, dd);
    location(*ww, *hh, din);
    rotation(rotationAngle, r);
    trans(din, r, dout);
    *ww = locationX(dout);
    *hh = locationY(dout);
}

/*
 * ****************************
 * Inverse Transformations
 * ****************************
 */

/*
 * Take a value in inches within the viewport and convert to some
 * other coordinate system
 */

double transformFromINCHES(double value, int unit, R_GE_gcontext *gc, double thisCM, double otherCM, GEDevDesc *dd)
{
    /*
     * Convert to NPC
     */
    double result = value;
    switch (unit)
    {
    case L_NPC:
        result = result / (thisCM / 2.54);
        break;
    case L_CM:
        result = result * 2.54;
        break;
    case L_INCHES:
        break;
    /* FIXME:  The following two assume that the pointsize specified
     * by the user is actually the pointsize provided by the
     * device.  This is NOT a safe assumption
     * One possibility would be to do a call to GReset(), just so
     * that mapping() gets called, just so that things like
     * xNDCPerLine are up-to-date, THEN call GStrHeight("M")
     * or somesuch.
     */
    case L_CHAR:
        result = (result * 72) / (gc->ps * gc->cex);
        break;
    case L_LINES:
        result = (result * 72) / (gc->ps * gc->cex * gc->lineheight);
        break;
    case L_MM:
        result = result * 2.54 * 10;
        break;
    /* Maybe an opportunity for some constants below here (!)
     */
    case L_POINTS:
        result = result * 72.27;
        break;
    case L_PICAS:
        result = result / 12 * 72.27;
        break;
    case L_BIGPOINTS:
        result = result * 72;
        break;
    case L_DIDA:
        result = result / 1238 * 1157 * 72.27;
        break;
    case L_CICERO:
        result = result / 1238 * 1157 * 72.27 / 12;
        break;
    case L_SCALEDPOINTS:
        result = result * 65536 * 72.27;
        break;
    /*
     * I'm not sure the remaining ones makes any sense.
     * For simplicity, these are just forbidden for now.
     */
    case L_SNPC:
    case L_MYCHAR:
    case L_MYLINES:
    case L_STRINGWIDTH:
    case L_MYSTRINGWIDTH:
    case L_STRINGHEIGHT:
    case L_MYSTRINGHEIGHT:
    case L_GROBWIDTH:
    case L_GROBHEIGHT:
    case L_NULL:
    default:
        error(_("Illegal unit or unit not yet implemented"));
    }
    return result;
}

/*
 * This corresponds to transform[X|Y]toINCHES() because
 * it works only within the current viewport, BUT
 * it is much simpler because it is supplied with a
 * double value in INCHES (rather than a unit object in
 * an arbitrary unit).
 *
 * For conceptual symmetry, it should probably return a
 * unit object, but it only returns a double value.
 * The construction of a unit object with the appropriate
 * unit must be performed by the calling function (or higher).
 * This is probably easiest done right up in R code.
 */
double transformXYFromINCHES(double location, int unit, double scalemin, double scalemax, R_GE_gcontext *gc,
                             double thisCM, double otherCM, GEDevDesc *dd)
{
    double result = location;
    switch (unit)
    {
    case L_NATIVE:
        result = scalemin + (result / (thisCM / 2.54)) * (scalemax - scalemin);
        break;
    default:
        result = transformFromINCHES(location, unit, gc, thisCM, otherCM, dd);
    }
    return result;
}

double transformWidthHeightFromINCHES(double dimension, int unit, double scalemin, double scalemax, R_GE_gcontext *gc,
                                      double thisCM, double otherCM, GEDevDesc *dd)
{
    double result = dimension;
    switch (unit)
    {
    case L_NATIVE:
        result = (result / (thisCM / 2.54)) * (scalemax - scalemin);
        break;
    default:
        result = transformFromINCHES(dimension, unit, gc, thisCM, otherCM, dd);
    }
    return result;
}

/* Attempt to make validating units faster
 */
typedef struct
{
    char *name;
    int code;
} UnitTab;

/* NOTE this table must match the order in grid.h
 */
static UnitTab UnitTable[] = {{"npc", 0},     {"cm", 1},          {"inches", 2},       {"lines", 3},
                              {"native", 4},  {"null", 5},        {"snpc", 6},         {"mm", 7},
                              {"points", 8},  {"picas", 9},       {"bigpts", 10},      {"dida", 11},
                              {"cicero", 12}, {"scaledpts", 13},  {"strwidth", 14},    {"strheight", 15},

                              {"char", 18},   {"grobwidth", 19},  {"grobheight", 20},  {"mylines", 21},
                              {"mychar", 22}, {"mystrwidth", 23}, {"mystrheight", 24},

                              {NULL, -1}};

int convertUnit(SEXP unit, int index)
{
    int i = 0;
    int result = 0;
    int found = 0;
    while (result >= 0 && !found)
    {
        if (UnitTable[i].name == NULL)
            result = -1;
        else
        {
            found = !strcmp(CHAR(STRING_ELT(unit, index)), UnitTable[i].name);
            if (found)
                result = UnitTable[i].code;
        }
        i += 1;
    }
    if (result < 0)
        error(_("Invalid unit"));
    return result;
}

SEXP validUnits(SEXP units)
{
    int i;
    int n = LENGTH(units);
    SEXP answer = R_NilValue;
    if (n > 0)
    {
        if (isString(units))
        {
            PROTECT(answer = allocVector(INTSXP, n));
            for (i = 0; i < n; i++)
                INTEGER(answer)[i] = convertUnit(units, i);
            UNPROTECT(1);
        }
        else
        {
            error(_("Units must be character"));
        }
    }
    else
    {
        error(_("Units must be of length > 0"));
    }
    return answer;
}
