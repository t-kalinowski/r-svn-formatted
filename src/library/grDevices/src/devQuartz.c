/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2007-8  The R Foundation
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
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 *  Modular Quartz device for Mac OS X
 *
 *  Partially based on code by Byron Ellis
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if HAVE_AQUA

#include <Defn.h>
#include <Rinternals.h>
#define R_USE_PROTOTYPES 1
#include <R_ext/GraphicsEngine.h>
#include <R_ext/QuartzDevice.h>

#include "grDevices.h"
#ifdef SUPPORT_MBCS
#include <wchar.h>
#endif

#include <CoreFoundation/CoreFoundation.h>
#include <Carbon/Carbon.h>

#define QBE_NATIVE 1 /* either R.app or Cocoa or Carbon depending on the OS X version */
#define QBE_COCOA 2  /* internal Cocoa */
#define QBE_CARBON 3 /* internal Carbon */
#define QBE_BITMAP 4 /* bitmap file creating */
#define QBE_PDF 5    /* PDF file creating */

typedef struct moduleTypes_s
{
    const char *type;
    const char *subst;
    int qbe; /* Quartz back-end */
} quartz_module_t;

/* list of internally supported output modules */
const quartz_module_t quartz_modules[] = {{"", 0, QBE_NATIVE},
                                          {"native", 0, QBE_NATIVE},
                                          {"cocoa", 0, QBE_COCOA},
                                          {"carbon", 0, QBE_CARBON},
                                          {"pdf", 0, QBE_PDF},
                                          {"png", "public.png", QBE_BITMAP},
                                          {"jpeg", "public.jpeg", QBE_BITMAP},
                                          {"jpg", "public.jpeg", QBE_BITMAP},
                                          {"jpeg2000", "public.jpeg-2000", QBE_BITMAP},
                                          {"tiff", "public.tiff", QBE_BITMAP},
                                          {"tif", "public.tiff", QBE_BITMAP},
                                          {"gif", "com.compuserve.gif", QBE_BITMAP},
                                          {"psd", "com.adobe.photoshop.image", QBE_BITMAP},
                                          {"bmp", "com.microsoft.bmp", QBE_BITMAP},
                                          {"sgi", "com.sgi.sgi-image", QBE_BITMAP},
                                          {"pict", "com.apple.pict", QBE_BITMAP},
                                          {0, 0, 0}};

/* for compatibility with OS X <10.5 */
#ifndef CGFLOAT_DEFINED
typedef float CGFloat;
#define CGFLOAT_MIN FLT_MIN
#define CGFLOAT_MAX FLT_MAX
#define CGFLOAT_IS_DOUBLE 0
#define CGFLOAT_DEFINED 1
#endif

typedef struct QuartzSpecific_s
{
    double ps;
    double scalex, scaley; /* resolution correction: px/pt ratio */
    double width, height;  /* size (in inches) */
    double tscale;         /* text scale (resolution independent,
                              i.e. it constitutes a text zoom factor */
    int dirty;             /* dirtly flag. Not acted upon by the Quartz
                              core, but QC sets it whenever a drawing
                              operation is performed (see detailed
                              description in R_ext/QuartzDevice.h) */
    int gstate;            /* gstate counter */
    int async;             /* asynchronous drawing (i.e. context was
                              not ready for an operation) */
    int bg;                /* background color */
    int antialias, smooth; /* smoothing flags (only aa makes any sense) */
    int flags;             /* additional QDFLAGs */
    int redraw;            /* redraw flag is set when replaying
                              and inhibits syncs on Mode */
    CGRect clipRect;       /* clipping rectangle */
    NewDevDesc *dev;       /* device structure holding this one */

    void *userInfo; /* pointer to a module-dependent space */

    /* callbacks - except for getCGContext all others are optional */
    CGContextRef (*getCGContext)(QuartzDesc_t dev, void *userInfo);
    int (*locatePoint)(QuartzDesc_t dev, void *userInfo, double *x, double *y);
    void (*close)(QuartzDesc_t dev, void *userInfo);
    void (*newPage)(QuartzDesc_t dev, void *userInfo, int flags);
    void (*state)(QuartzDesc_t dev, void *userInfo, int state);
    void *(*par)(QuartzDesc_t dev, void *userInfo, void *par);
    void (*sync)(QuartzDesc_t dev, void *userInfo);
} QuartzDesc;

/* coordinates:
   - R graphics (positions etc., usually points)
   - real size (e.g. inches)
   - display view (usually pixels)

   bookkeeping:
   - QuartzDevice.width/height:  inches
   - R GE size (.._Size): points
   - physical (on-screen) coordinates : pixels

the current implementation uses points as plotting units (i.e. this is what
Quartz tells R), but the canvas is specified in pixels. The scalex/y factors
specify the conversion factor between pixels and points.
We are *not* using R's scaling facilities, because R doesn't work with
non-square pixels (e.g. circles become ellipses).

FIXME: yes it does -- ipr is a two-element array.
Actually, dp not points are used.
*/

#pragma mark QuartzDevice API (for modules)

/* Update should be called when ps or tscale change.
   Conservatively, it should be called on scale change, too, in case
   we decide to abandon the CTM approach */
static void QuartzDevice_Update(QuartzDesc_t desc);

/* this function must be called after a new context is created.
   it primes the context for drawing */
void QuartzDevice_ResetContext(QuartzDesc_t desc)
{
    QuartzDesc *qd = ((QuartzDesc *)desc);
    qd->gstate = 0;
    qd->dirty = 0;
    if (qd->getCGContext)
    {
        CGContextRef ctx = qd->getCGContext(qd, qd->userInfo);
        if (ctx)
        {
            CGContextSetAllowsAntialiasing(ctx, qd->antialias);
            CGContextSetShouldSmoothFonts(ctx, qd->smooth);
            CGContextScaleCTM(ctx, qd->scalex, qd->scaley);
            CGContextSaveGState(ctx);
            qd->gstate = 1;
        }
    }
}

double QuartzDevice_GetScaledWidth(QuartzDesc_t desc)
{
    QuartzDesc *qd = ((QuartzDesc *)desc);
    return qd->scalex * qd->width * 72.0;
}
double QuartzDevice_GetScaledHeight(QuartzDesc_t desc)
{
    QuartzDesc *qd = ((QuartzDesc *)desc);
    return qd->scaley * qd->height * 72.0;
}
void QuartzDevice_SetScaledSize(QuartzDesc_t desc, double width, double height)
{
    QuartzDesc *qd = ((QuartzDesc *)desc);
    QuartzDevice_SetWidth(desc, width / qd->scalex / 72.0);
    QuartzDevice_SetHeight(desc, height / qd->scaley / 72.0);
}

int QuartzDevice_DevNumber(QuartzDesc_t desc)
{
    return ndevNumber((((QuartzDesc *)desc)->dev));
}

double QuartzDevice_GetWidth(QuartzDesc_t desc)
{
    return ((QuartzDesc *)desc)->width;
}
void QuartzDevice_SetWidth(QuartzDesc_t desc, double width)
{
    ((QuartzDesc *)desc)->width = width;
    ((QuartzDesc *)desc)->dev->right = width * 72.0;
}

double QuartzDevice_GetHeight(QuartzDesc_t desc)
{
    return ((QuartzDesc *)desc)->height;
}
void QuartzDevice_SetHeight(QuartzDesc_t desc, double height)
{
    ((QuartzDesc *)desc)->height = height;
    ((QuartzDesc *)desc)->dev->bottom = height * 72.0;
}

double QuartzDevice_GetXScale(QuartzDesc_t desc)
{
    return ((QuartzDesc *)desc)->scalex;
}
double QuartzDevice_GetYScale(QuartzDesc_t desc)
{
    return ((QuartzDesc *)desc)->scaley;
}
void QuartzDevice_SetScale(QuartzDesc_t desc, double scalex, double scaley)
{
    ((QuartzDesc *)desc)->scalex = scalex;
    ((QuartzDesc *)desc)->scaley = scaley;
    QuartzDevice_Update(desc);
}

double QuartzDevice_GetTextScale(QuartzDesc_t desc)
{
    return ((QuartzDesc *)desc)->tscale;
}

void QuartzDevice_SetTextScale(QuartzDesc_t desc, double scale)
{
    ((QuartzDesc *)desc)->tscale = scale;
    QuartzDevice_Update(desc);
}

double QuartzDevice_GetPointSize(QuartzDesc_t desc)
{
    return ((QuartzDesc *)desc)->ps;
}

void QuartzDevice_SetPointSize(QuartzDesc_t desc, double ps)
{
    ((QuartzDesc *)desc)->ps = ps;
    QuartzDevice_Update(desc);
}

int QuartzDevice_GetDirty(QuartzDesc_t desc)
{
    return ((QuartzDesc *)desc)->dirty;
}
void QuartzDevice_SetDirty(QuartzDesc_t desc, int dirty)
{
    ((QuartzDesc *)desc)->dirty = dirty;
}

int QuartzDevice_GetAntialias(QuartzDesc_t desc)
{
    return ((QuartzDesc *)desc)->antialias;
}
void QuartzDevice_SetAntialias(QuartzDesc_t desc, int aa)
{
    QuartzDesc *qd = (QuartzDesc *)desc;
    qd->antialias = aa;
    if (NULL != qd->getCGContext)
        CGContextSetAllowsAntialiasing(qd->getCGContext(qd, qd->userInfo), aa);
}

void QuartzDevice_Kill(QuartzDesc_t desc)
{
    GEDevDesc *dd = GEGetDevice(ndevNumber(((QuartzDesc *)desc)->dev));
    if (dd)
        GEkillDevice(dd);
}

int QuartzDesc_GetFontSmooth(QuartzDesc_t desc)
{
    return ((QuartzDesc *)desc)->smooth;
}
void QuartzDesc_SetFontSmooth(QuartzDesc_t desc, int fs)
{
    QuartzDesc *qd = (QuartzDesc *)desc;
    qd->smooth = fs;
    if (qd->getCGContext)
        CGContextSetShouldSmoothFonts(qd->getCGContext(qd, qd->userInfo), fs);
}

int QuartzDevice_GetBackground(QuartzDesc_t desc)
{
    return ((QuartzDesc *)desc)->bg;
}

static void QuartzDevice_Update(QuartzDesc_t desc)
{
    QuartzDesc *qd = (QuartzDesc *)desc;
    NewDevDesc *dev = qd->dev;

    /* pre-scaling happens in Quartz (using CTM), so scales should not be
       reflected in R measurements. We tell R to use 72dpi which corresponds
       to plotting in pt coordinates */
    dev->cra[0] = 0.9 * qd->ps * qd->tscale;
    dev->cra[1] = 1.2 * qd->ps * qd->tscale;
    dev->ipr[0] = 1.0 / 72.0;
    dev->ipr[1] = 1.0 / 72.0;
}

void QuartzDevice_ReplayDisplayList(QuartzDesc_t desc)
{
    QuartzDesc *qd = (QuartzDesc *)desc;
    int _dirty = qd->dirty;
    qd->redraw = 1;
    if (qd->dev->displayList != R_NilValue)
        GEplayDisplayList(GEGetDevice(ndevNumber(qd->dev)));
    qd->redraw = 0;
    qd->dirty = _dirty; /* we do NOT change the dirty flag */
}

void *QuartzDevice_GetSnapshot(QuartzDesc_t desc, int last)
{
    QuartzDesc *qd = (QuartzDesc *)desc;
    GEDevDesc *gd = GEGetDevice(ndevNumber(qd->dev));
    SEXP snap;
    if (last)
        snap = qd->dev->savedSnapshot;
    else
        snap = GEcreateSnapshot(gd);
    if (R_NilValue == VECTOR_ELT(snap, 0))
        snap = 0;
    return (snap == R_NilValue) ? 0 : snap;
}

void QuartzDevice_RestoreSnapshot(QuartzDesc_t desc, void *snap)
{
    QuartzDesc *qd = (QuartzDesc *)desc;
    GEDevDesc *gd = GEGetDevice(ndevNumber(qd->dev));
    if (NULL == snap)
        return; /*Aw, hell no!*/
    PROTECT((SEXP)snap);
    if (R_NilValue == VECTOR_ELT(snap, 0))
        warning("Tried to restore an empty snapshot?");
    qd->redraw = 1;
    GEplaySnapshot((SEXP)snap, gd);
    qd->redraw = 0;
    qd->dirty = 0; /* we reset the dirty flag */
    UNPROTECT(1);
}

#if 0
/* FIXME: these are concepts in base graphics, and so should not
   be present in a graphics device */
#include <Rdevices.h> /* for GetDevice */
double QuartzDevice_UserX(QuartzDesc_t desc,double x) { 
    return GConvertX(x,GMapUnits(0),GMapUnits(1),
		     GetDevice(ndevNumber(((QuartzDesc*)desc)->dev))); 
}
double QuartzDevice_UserY(QuartzDesc_t desc,double y) {
    return GConvertX(y,GMapUnits(0),GMapUnits(1),
		     GetDevice(ndevNumber(((QuartzDesc*)desc)->dev)));
}
#endif

#pragma mark RGD API Function Prototypes

static Rboolean RQuartz_Open(NewDevDesc *, QuartzDesc *, char *, double, double, int);
static void RQuartz_Close(NewDevDesc *);
static void RQuartz_Activate(NewDevDesc *);
static void RQuartz_Deactivate(NewDevDesc *);
static void RQuartz_Size(double *, double *, double *, double *, NewDevDesc *);
static void RQuartz_NewPage(R_GE_gcontext *, NewDevDesc *);
static void RQuartz_Clip(double, double, double, double, NewDevDesc *);
static double RQuartz_StrWidth(const char *, R_GE_gcontext *, NewDevDesc *);
static void RQuartz_Text(double, double, const char *, double, double, R_GE_gcontext *, NewDevDesc *);
static void RQuartz_Rect(double, double, double, double, R_GE_gcontext *, NewDevDesc *);
static void RQuartz_Circle(double, double, double, R_GE_gcontext *, NewDevDesc *);
static void RQuartz_Line(double, double, double, double, R_GE_gcontext *, NewDevDesc *);
static void RQuartz_Polyline(int, double *, double *, R_GE_gcontext *, NewDevDesc *);
static void RQuartz_Polygon(int, double *, double *, R_GE_gcontext *, NewDevDesc *);
static Rboolean RQuartz_Locator(double *, double *, NewDevDesc *);
static void RQuartz_Mode(int mode, NewDevDesc *);
static void RQuartz_Hold(NewDevDesc *);
static void RQuartz_MetricInfo(int, R_GE_gcontext *, double *, double *, double *, NewDevDesc *);

#pragma mark Quartz device implementation

void *QuartzDevice_Create(
    void *_dev, double scalex, double scaley, double ps, double width, double height, int bg, int aa, int fs,
    CGContextRef (*getCGContext)(QuartzDesc_t dev, void *userInfo), /* Get the context for this device */
    int (*locatePoint)(QuartzDesc_t dev, void *userInfo, double *x, double *y),
    void (*close)(QuartzDesc_t dev, void *userInfo), void (*newPage)(QuartzDesc_t dev, void *userInfo, int flags),
    void (*state)(QuartzDesc_t dev, void *userInfo, int state),
    void *(*par)(QuartzDesc_t dev, void *userInfo, void *par), void (*sync)(QuartzDesc_t dev, void *userInfo),
    void *userInfo)
{
    NewDevDesc *dev = (NewDevDesc *)_dev;
    dev->displayList = R_NilValue;

    dev->startfill = R_RGB(255, 255, 255);
    dev->startcol = R_RGB(0, 0, 0);
    dev->startps = ps;
    dev->startfont = 1;
    dev->startlty = LTY_SOLID;
    dev->startgamma = 1;

    /* Set up some happy pointers */
    dev->close = RQuartz_Close;
    dev->activate = RQuartz_Activate;
    dev->deactivate = RQuartz_Deactivate;
    dev->size = RQuartz_Size;
    dev->newPage = RQuartz_NewPage;
    dev->clip = RQuartz_Clip;
    dev->strWidth = RQuartz_StrWidth;
    dev->text = RQuartz_Text;
    dev->rect = RQuartz_Rect;
    dev->circle = RQuartz_Circle;
    dev->line = RQuartz_Line;
    dev->polyline = RQuartz_Polyline;
    dev->polygon = RQuartz_Polygon;
    dev->locator = RQuartz_Locator;
    dev->mode = RQuartz_Mode;
    dev->hold = RQuartz_Hold;
    dev->metricInfo = RQuartz_MetricInfo;
    dev->hasTextUTF8 = FALSE;

    dev->left = 0;
    dev->top = 0;

    /* Magic numbers from on high. */
    dev->xCharOffset = 0.4900;
    dev->yCharOffset = 0.3333;
    dev->yLineBias = 0.20; /* This is .2 for PS/PDF devices... */

    dev->canResizePlot = TRUE;
    dev->canChangeFont = TRUE;
    dev->canRotateText = TRUE;
    dev->canResizeText = TRUE;
    dev->canClip = TRUE;
    dev->canHAdj = 2;
    dev->canChangeGamma = TRUE;
    dev->displayListOn = TRUE;

    QuartzDesc *qd = calloc(1, sizeof(QuartzDesc));
    qd->width = width;
    qd->height = height;
    qd->userInfo = userInfo;
    qd->getCGContext = getCGContext;
    qd->locatePoint = locatePoint;
    qd->close = close;
    qd->newPage = newPage;
    qd->state = state;
    qd->sync = sync;
    qd->scalex = scalex;
    qd->scaley = scaley;
    qd->tscale = 1.0;
    qd->ps = ps;
    qd->bg = bg;
    qd->antialias = aa;
    qd->flags = fs;
    qd->gstate = 0;

    dev->deviceSpecific = qd;
    qd->dev = dev;

    QuartzDevice_Update(qd);

    dev->right = width * 72.0;
    dev->bottom = height * 72.0;
    qd->clipRect = CGRectMake(0, 0, dev->right, dev->bottom);

    qd->dirty = 0;
    qd->redraw = 0;
    qd->async = 0;
    return (QuartzDesc_t)qd;
}

/* old OS X versions has different names for some of the CGFont stuff */
#if MAC_OS_X_VERSION_MAX_ALLOWED <= MAC_OS_X_VERSION_10_4
#define CGFontCreateWithFontName CGFontCreateWithName
#define CGFontGetGlyphBBoxes CGFontGetGlyphBoundingBoxes
#define CGFontGetGlyphsForUnichars CGFontGetGlyphsForUnicodes
/* and some missing declarations */
extern CGFontRef CGFontCreateWithName(CFStringRef);
extern bool CGFontGetGlyphAdvances(CGFontRef font, const CGGlyph glyphs[], size_t count, int advances[]);
extern int CGFontGetUnitsPerEm(CGFontRef font);
extern bool CGFontGetGlyphBBoxes(CGFontRef font, const CGGlyph glyphs[], size_t count, CGRect bboxes[]);
#endif

/* These are internal (GlyphsForUnichars didn't used to be... Anyway...) */
extern CGFontRef CGContextGetFont(CGContextRef);
extern void CGFontGetGlyphsForUnichars(CGFontRef, const UniChar[], const CGGlyph[], size_t);

#define DEVDESC NewDevDesc *dd
#define CTXDESC R_GE_gcontext *gc, NewDevDesc *dd

#define DEVSPEC                                                                                                        \
    QuartzDesc *xd = (QuartzDesc *)dd->deviceSpecific;                                                                 \
    CGContextRef ctx = xd->getCGContext(xd, xd->userInfo)
#define DRAWSPEC                                                                                                       \
    QuartzDesc *xd = (QuartzDesc *)dd->deviceSpecific;                                                                 \
    CGContextRef ctx = xd->getCGContext(xd, xd->userInfo);                                                             \
    xd->dirty = 1
#define XD QuartzDesc *xd = (QuartzDesc *)dd->deviceSpecific
#pragma mark Device Implementation

CFStringRef RQuartz_FindFont(int fontface, char *fontfamily)
{
    SEXP ns, env, db, names;
    PROTECT_INDEX index;
    CFStringRef fontName = CFSTR("");
    PROTECT(ns = R_FindNamespace(ScalarString(mkChar("grDevices"))));
    PROTECT_WITH_INDEX(env = findVar(install(".Quartzenv"), ns), &index);
    if (TYPEOF(env) == PROMSXP)
        REPROTECT(env = eval(env, ns), index);
    PROTECT(db = findVar(install(".Quartz.Fonts"), env));
    PROTECT(names = getAttrib(db, R_NamesSymbol));
    if (strlen(fontfamily) > 0)
    {
        int i;
        for (i = 0; i < length(names); i++)
            if (0 == strcmp(fontfamily, CHAR(STRING_ELT(names, i))))
                break;
        if (i < length(names))
            fontName = CFStringCreateWithCString(kCFAllocatorDefault, CHAR(STRING_ELT(VECTOR_ELT(db, i), fontface)),
                                                 kCFStringEncodingUTF8);
    }
    UNPROTECT(4);
    return fontName;
}

CGFontRef RQuartz_Font(CTXDESC)
{
    int fontface = gc->fontface;
    CFMutableStringRef fontName = CFStringCreateMutable(kCFAllocatorDefault, 0);
    if ((gc->fontface == 5) || (strcmp(gc->fontfamily, "symbol") == 0))
        CFStringAppend(fontName, CFSTR("Symbol"));
    else
    {
        CFStringRef font = RQuartz_FindFont(gc->fontface, gc->fontfamily);
        if (CFStringGetLength(font) > 0)
        {
            fontface = 1; /* This is handled by the lookup process */
            CFStringAppend(fontName, font);
        }
        CFRelease(font);
    }
    if (CFStringGetLength(fontName) == 0)
        CFStringAppend(fontName, CFSTR("Arial"));
    if (fontface == 2 || fontface == 4)
    {
        CFStringAppend(fontName, CFSTR(" Bold"));
    }
    if (fontface == 3)
    {
        CFStringAppend(fontName, CFSTR(" Italic"));
    }
    CGFontRef font = CGFontCreateWithFontName(fontName);
    if (font == 0)
    {
        /* Fall back on ATS */
        ATSFontRef tmp = ATSFontFindFromName(fontName, kATSOptionFlagsDefault);
        font = CGFontCreateWithPlatformFont(&tmp);
    }
    if (NULL == font)
    {
        CFShow(fontName);
    }
    CFRelease(fontName);
    return font;
}

#define RQUARTZ_FILL (1)
#define RQUARTZ_STROKE (1 << 1)
#define RQUARTZ_LINE (1 << 2)
#define RQUARTZ_FONT (1 << 3)

void RQuartz_Set(CGContextRef ctx, R_GE_gcontext *gc, int flags)
{
    if (flags & RQUARTZ_FILL)
    {
        int fill = gc->fill;
        CGContextSetRGBFillColor(ctx, R_RED(fill) / 256.0, R_GREEN(fill) / 256.0, R_BLUE(fill) / 256.0,
                                 R_ALPHA(fill) / 256.0);
    }
    if (flags & RQUARTZ_STROKE)
    {
        int stroke = gc->col;
        CGContextSetRGBStrokeColor(ctx, R_RED(stroke) / 256.0, R_GREEN(stroke) / 256.0, R_BLUE(stroke) / 256.0,
                                   R_ALPHA(stroke) / 256.0);
    }
    if (flags & RQUARTZ_LINE)
    {
        CGFloat dashlist[8];
        int i, ndash = 0;
        int lty = gc->lty;
        CGContextSetLineWidth(ctx, gc->lwd);

        float lwd = gc->lwd * 0.75;
        for (i = 0; i < 8 && lty; i++)
        {
            dashlist[ndash++] = (lwd >= 1 ? lwd : 1) * (lty & 15);
            lty >>= 4;
        }
        CGContextSetLineDash(ctx, 0, dashlist, ndash);
        CGLineCap cap = kCGLineCapButt;
        switch (gc->lend)
        {
        case GE_ROUND_CAP:
            cap = kCGLineCapRound;
            break;
        case GE_BUTT_CAP:
            cap = kCGLineCapButt;
            break;
        case GE_SQUARE_CAP:
            cap = kCGLineCapSquare;
            break;
        }
        CGContextSetLineCap(ctx, cap);
        CGLineJoin join = kCGLineJoinRound;
        switch (gc->ljoin)
        {
        case GE_ROUND_JOIN:
            join = kCGLineJoinRound;
            break;
        case GE_MITRE_JOIN:
            join = kCGLineJoinMiter;
            break;
        case GE_BEVEL_JOIN:
            join = kCGLineJoinBevel;
            break;
        }
        CGContextSetLineJoin(ctx, join);
        CGContextSetMiterLimit(ctx, gc->lmitre);
    }
    if (flags & RQUARTZ_FONT)
    {
        CGFontRef font = RQuartz_Font(gc, NULL);
        CGContextSetFont(ctx, font);
        CGContextSetFontSize(ctx, gc->cex * gc->ps);
    }
}

#define SET(X) RQuartz_Set(ctx, gc, (X))
#define NOCTX                                                                                                          \
    {                                                                                                                  \
        xd->async = 1;                                                                                                 \
        return;                                                                                                        \
    }
#define NOCTXR(V)                                                                                                      \
    {                                                                                                                  \
        xd->async = 1;                                                                                                 \
        return (V);                                                                                                    \
    }

static Rboolean RQuartz_Open(DEVDESC, QuartzDesc *xd, char *display, double width, double height, int bg)
{
    return TRUE;
}

static void RQuartz_Close(DEVDESC)
{
    XD;
    if (xd->close)
        xd->close(xd, xd->userInfo);
}

static void RQuartz_Activate(DEVDESC)
{
    XD;
    if (xd->state)
        xd->state(xd, xd->userInfo, 1);
}

static void RQuartz_Deactivate(DEVDESC)
{
    XD;
    if (xd->state)
        xd->state(xd, xd->userInfo, 0);
}

static void RQuartz_Size(double *left, double *right, double *bottom, double *top, DEVDESC)
{
    XD;
    *left = *top = 0;
    *right = QuartzDevice_GetWidth(xd) * 72.0;
    *bottom = QuartzDevice_GetHeight(xd) * 72.0;
}

static void RQuartz_NewPage(CTXDESC)
{
    {
        DRAWSPEC;
        ctx = NULL;
        if (xd->newPage)
            xd->newPage(xd, xd->userInfo, xd->redraw ? QNPF_REDRAW : 0);
    }
    { /* we have to re-fetch the status *after* newPage since it may have changed it */
        DRAWSPEC;
        if (!ctx)
            NOCTX;
        SET(RQUARTZ_FILL);
        {
            CGRect bounds = CGRectMake(0, 0, QuartzDevice_GetWidth(xd) * 72.0, QuartzDevice_GetHeight(xd) * 72.0);
            if (R_ALPHA(xd->bg) == 255 && R_ALPHA(gc->fill) == 255)
                CGContextClearRect(ctx, bounds);
            CGContextFillRect(ctx, bounds);
        }
    }
}

static void RQuartz_Clip(double x0, double x1, double y0, double y1, DEVDESC)
{
    DRAWSPEC;
    if (!ctx)
        NOCTX;
    if (xd->gstate > 0)
    {
        --xd->gstate;
        CGContextRestoreGState(ctx);
    }
    CGContextSaveGState(ctx);
    xd->gstate++;
    if (x1 > x0)
    {
        double t = x1;
        x1 = x0;
        x0 = t;
    }
    if (y1 > y0)
    {
        double t = y1;
        y1 = y0;
        y0 = t;
    }
    xd->clipRect = CGRectMake(x0, y0, x1 - x0, y1 - y0);
    CGContextClipToRect(ctx, xd->clipRect);
}

CFStringRef prepareText(CTXDESC, const char *text, UniChar **buffer, int *free)
{
    CFStringRef str;
    if (gc->fontface == 5 || strcmp(gc->fontfamily, "symbol") == 0)
        str = CFStringCreateWithCString(NULL, text, kCFStringEncodingMacSymbol);
    else
    {
        str = CFStringCreateWithCString(NULL, text, kCFStringEncodingUTF8);
        /* Try fallback Latin1 encoding if UTF8 doesn't work. */
        if (NULL == str)
            CFStringCreateWithCString(NULL, text, kCFStringEncodingISOLatin1);
    }
    /* FIXME: this can fail, e.g. for 0x7f in the symbol font.  Why? */
    *buffer = (UniChar *)CFStringGetCharactersPtr(str);
    if (*buffer == NULL)
    {
        CFIndex length = CFStringGetLength(str);
        *buffer = malloc(length * sizeof(UniChar));
        CFStringGetCharacters(str, CFRangeMake(0, length), *buffer);
        *free = 1;
    }
    return str;
}

static double RQuartz_StrWidth(const char *text, CTXDESC)
{
    DEVSPEC;
    if (!ctx)
        NOCTXR(strlen(text) * 10.0); /* for sanity reasons */
    SET(RQUARTZ_FONT);
    {
        CGFontRef font = CGContextGetFont(ctx);
        float aScale = (gc->cex * gc->ps * xd->tscale) / CGFontGetUnitsPerEm(font);
        UniChar *buffer;
        CGGlyph *glyphs;
        int *advances;
        int Free = 0, len, i;
        CFStringRef str = prepareText(gc, dd, text, &buffer, &Free);
        len = CFStringGetLength(str);
        glyphs = malloc(sizeof(CGGlyph) * len);
        advances = malloc(sizeof(int) * len);
        CGFontGetGlyphsForUnichars(font, buffer, glyphs, len);
        CGFontGetGlyphAdvances(font, glyphs, len, advances);
        {
            float width = 0.0; /* aScale*CGFontGetLeading(CGContextGetFont(ctx)); */
            for (i = 0; i < len; i++)
                width += aScale * advances[i];
            free(advances);
            free(glyphs);
            if (Free)
                free(buffer);
            CFRelease(str);
            return width;
        }
    }
}

static void RQuartz_Text(double x, double y, const char *text, double rot, double hadj, CTXDESC)
{
    DRAWSPEC;
    if (!ctx)
        NOCTX;
    /* A stupid hack because R isn't consistent. */
    int fill = gc->fill;
    gc->fill = gc->col;
    SET(RQUARTZ_FILL | RQUARTZ_STROKE | RQUARTZ_FONT);
    gc->fill = fill;
    CGFontRef font = CGContextGetFont(ctx);
    float aScale = (gc->cex * gc->ps * xd->tscale) / (CGFontGetUnitsPerEm(font));
    UniChar *buffer;
    CGGlyph *glyphs;

    int Free = 0, len, i;
    float width = 0.0;
    CFStringRef str = prepareText(gc, dd, text, &buffer, &Free);
    len = CFStringGetLength(str);
    glyphs = malloc(sizeof(CGGlyph) * len);
    CGFontGetGlyphsForUnichars(font, buffer, glyphs, len);
    int *advances = malloc(sizeof(int) * len);
    CGSize *g_adv = malloc(sizeof(CGSize) * len);

    CGFontGetGlyphAdvances(font, glyphs, len, advances);
    for (i = 0; i < len; i++)
    {
        width += advances[i] * aScale;
        g_adv[i] = CGSizeMake(aScale * advances[i] * cos(-0.0174532925 * rot),
                              aScale * advances[i] * sin(-0.0174532925 * rot));
    }
    free(advances);
    CGContextSetTextMatrix(ctx, CGAffineTransformConcat(CGAffineTransformMakeScale(1.0, -1.0),
                                                        CGAffineTransformMakeRotation(-0.0174532925 * rot)));
    double ax = (width * hadj) * cos(-0.0174532925 * rot);
    double ay = (width * hadj) * sin(-0.0174532925 * rot);
    /*      double h  = CGFontGetXHeight(CGContextGetFont(ctx))*aScale; */
    CGContextSetTextPosition(ctx, x - ax, y - ay);
    /*      Rprintf("%s,%.2f %.2f (%.2f,%.2f)
     * (%d,%f)\n",text,hadj,width,ax,ay,CGFontGetUnitsPerEm(CGContextGetFont(ctx)),CGContextGetFontSize(ctx));       */
    CGContextShowGlyphsWithAdvances(ctx, glyphs, g_adv, len);
    free(glyphs);
    free(g_adv);
    if (Free)
        free(buffer);
    CFRelease(str);
}

static void RQuartz_Rect(double x0, double y0, double x1, double y1, CTXDESC)
{
    DRAWSPEC;
    if (!ctx)
        NOCTX;
    SET(RQUARTZ_FILL | RQUARTZ_STROKE | RQUARTZ_LINE);
    CGContextBeginPath(ctx);
    CGContextAddRect(ctx, CGRectMake(x0, y0, x1 - x0, y1 - y0));
    CGContextDrawPath(ctx, kCGPathFillStroke);
}

static void RQuartz_Circle(double x, double y, double r, CTXDESC)
{
    DRAWSPEC;
    if (!ctx)
        NOCTX;
    SET(RQUARTZ_FILL | RQUARTZ_STROKE | RQUARTZ_LINE);
    double r2 = 2.0 * r;
    CGContextBeginPath(ctx);
    CGContextAddEllipseInRect(ctx, CGRectMake(x - r, y - r, r2, r2));
    CGContextDrawPath(ctx, kCGPathFillStroke);
}

static void RQuartz_Line(double x1, double y1, double x2, double y2, CTXDESC)
{
    DRAWSPEC;
    if (!ctx)
        NOCTX;
    SET(RQUARTZ_STROKE | RQUARTZ_LINE);
    CGContextBeginPath(ctx);
    CGContextMoveToPoint(ctx, x1, y1);
    CGContextAddLineToPoint(ctx, x2, y2);
    CGContextStrokePath(ctx);
}

static void RQuartz_Polyline(int n, double *x, double *y, CTXDESC)
{
    if (n < 2)
        return;
    int i;
    DRAWSPEC;
    if (!ctx)
        NOCTX;
    SET(RQUARTZ_STROKE | RQUARTZ_LINE);
    CGContextBeginPath(ctx);
    CGContextMoveToPoint(ctx, x[0], y[0]);
    for (i = 1; i < n; i++)
        CGContextAddLineToPoint(ctx, x[i], y[i]);
    CGContextStrokePath(ctx);
}

static void RQuartz_Polygon(int n, double *x, double *y, CTXDESC)
{
    if (n < 2)
        return;
    int i;
    DRAWSPEC;
    if (!ctx)
        NOCTX;
    SET(RQUARTZ_FILL | RQUARTZ_STROKE | RQUARTZ_LINE);
    CGContextBeginPath(ctx);
    CGContextMoveToPoint(ctx, x[0], y[0]);
    for (i = 1; i < n; i++)
        CGContextAddLineToPoint(ctx, x[i], y[i]);
    CGContextClosePath(ctx);
    CGContextDrawPath(ctx, kCGPathFillStroke);
}

static void RQuartz_Mode(int mode, DEVDESC)
{
    DEVSPEC;
    if (!ctx)
        NOCTX;
    /* don't do anything in redraw as we can signal the end */
    if (xd->redraw)
        return;
    /* mode=0 -> drawing complete, signal sync */
    if (mode == 0)
    {
        if (xd->sync)
            xd->sync(xd, xd->userInfo);
        else
            CGContextSynchronize(ctx);
    }
}

static void RQuartz_Hold(DEVDESC)
{
}

static void RQuartz_MetricInfo(int c, R_GE_gcontext *gc, double *ascent, double *descent, double *width, NewDevDesc *dd)
{
    DRAWSPEC;
    if (!ctx)
    { /* dummy data if we have no context, for sanity reasons */
        *ascent = 10.0;
        *descent = 2.0;
        *width = 9.0;
        NOCTX;
    }
    SET(RQUARTZ_FONT);
    {
        CGFontRef font = CGContextGetFont(ctx);
        float aScale = (gc->cex * gc->ps * xd->tscale) / CGFontGetUnitsPerEm(font);
        UniChar *buffer, single;
        CGGlyph glyphs[8];
        CFStringRef str = NULL;
        int free_buffer = 0, len;
        if (c >= 0 && c <= ((mbcslocale && gc->fontface != 5) ? 127 : 255))
        {
            char text[2] = {c, 0};
            str = prepareText(gc, dd, text, &buffer, &free_buffer);
            len = CFStringGetLength(str);
            if (len > 7)
                return; /* this is basically impossible,
               but you never know */
        }
        else
        {
            single = (UniChar)((c < 0) ? -c : c);
            buffer = &single;
            len = 1;
        }
        *width = 0.0;
        CGFontGetGlyphsForUnichars(font, buffer, glyphs, len);
        {
            int i;
            int advances[8];
            CGRect bboxes[8];
            CGFontGetGlyphAdvances(font, glyphs, len, advances);
            CGFontGetGlyphBBoxes(font, glyphs, len, bboxes);
            for (i = 0; i < len; i++)
                *width += advances[i] * aScale;
            *ascent = aScale * (bboxes[0].size.height + bboxes[0].origin.y);
            *descent = -aScale * bboxes[0].origin.y;
        }
        if (free_buffer)
            free(buffer);
        if (str)
            CFRelease(str);
    }
}

static Rboolean RQuartz_Locator(double *x, double *y, DEVDESC)
{
    Rboolean res;
    DEVSPEC;
    ctx = NULL;
    if (!xd->locatePoint)
        return FALSE;
    res = xd->locatePoint(xd, xd->userInfo, x, y);
    *x /= xd->scalex;
    *y /= xd->scaley;
    return res;
}

#pragma mark -
#pragma mark R Interface

#include "qdCocoa.h"
#include "qdBitmap.h"
#include "qdPDF.h"
/* disabled for now until we get to test in on 10.3 #include "qdCarbon.h" */

/* current fake */
Rboolean QuartzCarbon_DeviceCreate(NewDevDesc *dd, const char *type, const char *file, double width, double height,
                                   double pointsize, const char *family, Rboolean antialias, Rboolean smooth,
                                   Rboolean autorefresh, int quartzpos, int bg, const char *title, double *dpi)
{
    return FALSE;
}

#define ARG(HOW, WHAT)                                                                                                 \
    HOW(CAR(WHAT));                                                                                                    \
    WHAT = CDR(WHAT)

/* C version of the Quartz call (experimental)
   returns 0 on success, error code on failure */
int Quartz_C(const char *type, const char *file, double width, double height, double ps, const char *family, int aa,
             int fsm, const char *title, int bg, double *dpi, quartz_create_fn_t q_create)
{
    if (!q_create)
        return -3;
    {
        char *vmax = vmaxget();
        R_CheckDeviceAvailable();
        {
            /* FIXME: check this allocation */
            NewDevDesc *dev = calloc(1, sizeof(NewDevDesc));
            dev->displayList = R_NilValue;
            dev->savedSnapshot = R_NilValue;

            if (!dev)
                return -1;

            if (!q_create(dev, type, file, width, height, ps, family, aa, fsm, TRUE, 1, bg, title, dpi))
            {
                vmaxset(vmax);
                free(dev);
                return -2;
            }
            gsetVar(install(".Device"), mkString("quartz"), R_BaseEnv);
            GEDevDesc *dd = GEcreateDevDesc(dev);
            GEaddDevice(dd);
            GEinitDisplayList(dd);
            vmaxset(vmax);
        }
    }
    return 0;
}

/* ARGS: type, file, widht, height, ps, family, antialias, fontsm, title, bg, dpi */
SEXP Quartz(SEXP args)
{
    SEXP tmps, bgs;
    double width, height, ps;
    Rboolean antialias, smooth, autorefresh = TRUE, succ = FALSE;
    int quartzpos, bg, module = 0;
    double mydpi[2], *dpi = 0;
    const char *type;
    const char *file;
    const char *family;
    const char *title;

    char *vmax = vmaxget();
    /* Get function arguments */
    args = CDR(args); /* Skip the call */
    if (TYPEOF(CAR(args)) != STRSXP || LENGTH(CAR(args)) < 1)
        type = "";
    else
        type = CHAR(STRING_ELT(CAR(args), 0));
    args = CDR(args);
    /* we may want to support connections at some point, but not yet ... */
    if (TYPEOF(CAR(args)) != STRSXP || LENGTH(CAR(args)) < 1)
        file = 0;
    else
        file = CHAR(STRING_ELT(CAR(args), 0));
    args = CDR(args);
    width = ARG(asReal, args);
    height = ARG(asReal, args);
    ps = ARG(asReal, args);
    family = CHAR(STRING_ELT(CAR(args), 0));
    args = CDR(args);
    antialias = ARG(asLogical, args);
    smooth = ARG(asLogical, args);
    title = CHAR(STRING_ELT(CAR(args), 0));
    args = CDR(args);
    bgs = CAR(args);
    args = CDR(args);
    /* FIXME: we should process bgs here ... somehow ... */
    tmps = CAR(args);
    args = CDR(args);
    if (!isNull(tmps))
    {
        tmps = coerceVector(tmps, REALSXP);
        if (LENGTH(tmps) > 0)
        {
            dpi = mydpi;
            mydpi[0] = REAL(tmps)[0];
            if (LENGTH(tmps) > 1)
                mydpi[1] = REAL(tmps)[1];
            else
                mydpi[1] = mydpi[0];
        }
    }
    /* just in case someone passed NAs/NaNs */
    if (dpi && (ISNAN(dpi[0]) || ISNAN(dpi[1])))
        dpi = 0;

    if (ISNAN(width) || ISNAN(height) || width <= 0.0 || height <= 0.0)
        error(_("invalid Quartz device size"));

    if (type)
    {
        const quartz_module_t *m = quartz_modules;
        while (m->type)
        {
            if (!strcasecmp(type, m->type))
            {
                module = m->qbe;
                if (m->subst)
                    type = m->subst;
                break;
            }
            m++;
        }
    }
    if (!strncasecmp(type, "bitmap:", 7))
    {
        module = QBE_BITMAP;
        type = type + 7;
    }

    bg = 0xffffffff;
    quartzpos = 1;

    R_CheckDeviceAvailable();
    BEGIN_SUSPEND_INTERRUPTS
    {
        NewDevDesc *dev = calloc(1, sizeof(NewDevDesc));
        dev->displayList = R_NilValue;
        dev->savedSnapshot = R_NilValue;

        if (!dev)
            error(_("Unable to create device description."));

        /* re-routed code has the first shot */
        if (ptr_QuartzDeviceCreate)
            succ = ptr_QuartzDeviceCreate(dev, type, file, width, height, ps, family, antialias, smooth, autorefresh,
                                          quartzpos, bg, title, dpi);

        if (!succ)
        { /* try internal modules next */
            switch (module)
            {
            case QBE_COCOA:
                succ = QuartzCocoa_DeviceCreate(dev, type, file, width, height, ps, family, antialias, smooth,
                                                autorefresh, quartzpos, bg, title, dpi);
                break;
            case QBE_NATIVE:
                /* native is essentially cocoa with carbon fall-back */
                succ = QuartzCocoa_DeviceCreate(dev, type, file, width, height, ps, family, antialias, smooth,
                                                autorefresh, quartzpos, bg, title, dpi);
                if (succ)
                    break;
            case QBE_CARBON:
                succ = QuartzCarbon_DeviceCreate(dev, type, file, width, height, ps, family, antialias, smooth,
                                                 autorefresh, quartzpos, bg, title, dpi);
                break;
            case QBE_PDF:
                succ = QuartzPDF_DeviceCreate(dev, type, file, width, height, ps, family, antialias, smooth,
                                              autorefresh, quartzpos, bg, title, dpi);
                break;
            case QBE_BITMAP:
                succ = QuartzBitmap_DeviceCreate(dev, type, file, width, height, ps, family, antialias, smooth,
                                                 autorefresh, quartzpos, bg, dpi);
                break;
            }
        }

        if (!succ)
        {
            vmaxset(vmax);
            free(dev);
            error(_("Unable to create Quartz device target, given type may not be supported."));
        }
        gsetVar(install(".Device"), mkString("quartz"), R_BaseEnv);
        GEDevDesc *dd = GEcreateDevDesc(dev);
        GEaddDevice(dd);
        GEinitDisplayList(dd);
    }
    END_SUSPEND_INTERRUPTS;
    vmaxset(vmax);
    return R_NilValue;
}

#else
/* --- no AQUA support = no Quartz --- */

#include <Defn.h>
#include <Graphics.h>
#include <Rdevices.h>
#include <Rinternals.h>
#include <Rgraphics.h>
#include <R_ext/GraphicsDevice.h>
#include <R_ext/GraphicsEngine.h>

#include "grDevices.h"

SEXP Quartz(SEXP args)
{
    warning(_("Quartz device is not available on this platform."));
    return R_NilValue;
}

#endif
