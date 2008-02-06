/* The beginning of code which represents an R base graphics system
 * separate from an R graphics engine (separate from R devices)
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <Defn.h>
#include <Graphics.h>
#include <Colors.h>
#include <R_ext/GraphicsBase.h>

int attribute_hidden baseRegisterIndex = -1;

static void restoredpSaved(pGEDevDesc dd)
{
    /* NOTE that not all params should be restored before playing */
    /* the display list (e.g., don't restore the device size) */

    int i, j, nr, nc;

    /* do NOT restore basic device driver properties;  they are */
    /* either meant to be different (e.g., left, right, bottom, top */
    /* changed because of window resize) or never change (e.g., ipr) */

    Rf_dpptr(dd)->state = Rf_dpSavedptr(dd)->state;
    Rf_dpptr(dd)->adj = Rf_dpSavedptr(dd)->adj;
    Rf_dpptr(dd)->ann = Rf_dpSavedptr(dd)->ann;
    Rf_dpptr(dd)->bg = Rf_dpSavedptr(dd)->bg;
    Rf_dpptr(dd)->bty = Rf_dpSavedptr(dd)->bty;
    Rf_dpptr(dd)->cex = Rf_dpSavedptr(dd)->cex;
    Rf_gpptr(dd)->lheight = Rf_dpSavedptr(dd)->lheight;
    Rf_dpptr(dd)->col = Rf_dpSavedptr(dd)->col;
    Rf_dpptr(dd)->crt = Rf_dpSavedptr(dd)->crt;
    Rf_dpptr(dd)->err = Rf_dpSavedptr(dd)->err;
    Rf_dpptr(dd)->fg = Rf_dpSavedptr(dd)->fg;
    Rf_dpptr(dd)->font = Rf_dpSavedptr(dd)->font;
    strncpy(Rf_dpptr(dd)->family, Rf_dpSavedptr(dd)->family, 201);
    Rf_dpptr(dd)->gamma = Rf_dpSavedptr(dd)->gamma;
    Rf_dpptr(dd)->lab[0] = Rf_dpSavedptr(dd)->lab[0];
    Rf_dpptr(dd)->lab[1] = Rf_dpSavedptr(dd)->lab[1];
    Rf_dpptr(dd)->lab[2] = Rf_dpSavedptr(dd)->lab[2];
    Rf_dpptr(dd)->las = Rf_dpSavedptr(dd)->las;
    Rf_dpptr(dd)->lty = Rf_dpSavedptr(dd)->lty;
    Rf_dpptr(dd)->lwd = Rf_dpSavedptr(dd)->lwd;
    Rf_dpptr(dd)->lend = Rf_dpSavedptr(dd)->lend;
    Rf_dpptr(dd)->ljoin = Rf_dpSavedptr(dd)->ljoin;
    Rf_dpptr(dd)->lmitre = Rf_dpSavedptr(dd)->lmitre;
    Rf_dpptr(dd)->mgp[0] = Rf_dpSavedptr(dd)->mgp[0];
    Rf_dpptr(dd)->mgp[1] = Rf_dpSavedptr(dd)->mgp[1];
    Rf_dpptr(dd)->mgp[2] = Rf_dpSavedptr(dd)->mgp[2];
    Rf_dpptr(dd)->mkh = Rf_dpSavedptr(dd)->mkh;
    Rf_dpptr(dd)->pch = Rf_dpSavedptr(dd)->pch;
    Rf_dpptr(dd)->ps = Rf_dpSavedptr(dd)->ps; /*was commented out --why?*/
    Rf_dpptr(dd)->smo = Rf_dpSavedptr(dd)->smo;
    Rf_dpptr(dd)->srt = Rf_dpSavedptr(dd)->srt;
    Rf_dpptr(dd)->tck = Rf_dpSavedptr(dd)->tck;
    Rf_dpptr(dd)->tcl = Rf_dpSavedptr(dd)->tcl;
    Rf_dpptr(dd)->xaxp[0] = Rf_dpSavedptr(dd)->xaxp[0];
    Rf_dpptr(dd)->xaxp[1] = Rf_dpSavedptr(dd)->xaxp[1];
    Rf_dpptr(dd)->xaxp[2] = Rf_dpSavedptr(dd)->xaxp[2];
    Rf_dpptr(dd)->xaxs = Rf_dpSavedptr(dd)->xaxs;
    Rf_dpptr(dd)->xaxt = Rf_dpSavedptr(dd)->xaxt;
    Rf_dpptr(dd)->xpd = Rf_dpSavedptr(dd)->xpd;
    Rf_dpptr(dd)->xlog = Rf_dpSavedptr(dd)->xlog;
    Rf_dpptr(dd)->yaxp[0] = Rf_dpSavedptr(dd)->yaxp[0];
    Rf_dpptr(dd)->yaxp[1] = Rf_dpSavedptr(dd)->yaxp[1];
    Rf_dpptr(dd)->yaxp[2] = Rf_dpSavedptr(dd)->yaxp[2];
    Rf_dpptr(dd)->yaxs = Rf_dpSavedptr(dd)->yaxs;
    Rf_dpptr(dd)->yaxt = Rf_dpSavedptr(dd)->yaxt;
    Rf_dpptr(dd)->ylog = Rf_dpSavedptr(dd)->ylog;
    Rf_dpptr(dd)->cexbase = Rf_dpSavedptr(dd)->cexbase;
    Rf_dpptr(dd)->cexmain = Rf_dpSavedptr(dd)->cexmain;
    Rf_dpptr(dd)->cexlab = Rf_dpSavedptr(dd)->cexlab;
    Rf_dpptr(dd)->cexsub = Rf_dpSavedptr(dd)->cexsub;
    Rf_dpptr(dd)->cexaxis = Rf_dpSavedptr(dd)->cexaxis;
    Rf_dpptr(dd)->fontmain = Rf_dpSavedptr(dd)->fontmain;
    Rf_dpptr(dd)->fontlab = Rf_dpSavedptr(dd)->fontlab;
    Rf_dpptr(dd)->fontsub = Rf_dpSavedptr(dd)->fontsub;
    Rf_dpptr(dd)->fontaxis = Rf_dpSavedptr(dd)->fontaxis;
    Rf_dpptr(dd)->colmain = Rf_dpSavedptr(dd)->colmain;
    Rf_dpptr(dd)->collab = Rf_dpSavedptr(dd)->collab;
    Rf_dpptr(dd)->colsub = Rf_dpSavedptr(dd)->colsub;
    Rf_dpptr(dd)->colaxis = Rf_dpSavedptr(dd)->colaxis;

    /* must restore layout parameters;	the different graphics */
    /* regions and coordinate transformations will be recalculated */
    /* but they need all of the layout information restored for this */
    /* to happen correctly */

    Rf_dpptr(dd)->devmode = Rf_dpSavedptr(dd)->devmode;
    Rf_dpptr(dd)->fig[0] = Rf_dpSavedptr(dd)->fig[0];
    Rf_dpptr(dd)->fig[1] = Rf_dpSavedptr(dd)->fig[1];
    Rf_dpptr(dd)->fig[2] = Rf_dpSavedptr(dd)->fig[2];
    Rf_dpptr(dd)->fig[3] = Rf_dpSavedptr(dd)->fig[3];
    Rf_dpptr(dd)->fin[0] = Rf_dpSavedptr(dd)->fin[0];
    Rf_dpptr(dd)->fin[1] = Rf_dpSavedptr(dd)->fin[1];
    Rf_dpptr(dd)->fUnits = Rf_dpSavedptr(dd)->fUnits;
    Rf_dpptr(dd)->defaultFigure = Rf_dpSavedptr(dd)->defaultFigure;
    Rf_dpptr(dd)->mar[0] = Rf_dpSavedptr(dd)->mar[0];
    Rf_dpptr(dd)->mar[1] = Rf_dpSavedptr(dd)->mar[1];
    Rf_dpptr(dd)->mar[2] = Rf_dpSavedptr(dd)->mar[2];
    Rf_dpptr(dd)->mar[3] = Rf_dpSavedptr(dd)->mar[3];
    Rf_dpptr(dd)->mai[0] = Rf_dpSavedptr(dd)->mai[0];
    Rf_dpptr(dd)->mai[1] = Rf_dpSavedptr(dd)->mai[1];
    Rf_dpptr(dd)->mai[2] = Rf_dpSavedptr(dd)->mai[2];
    Rf_dpptr(dd)->mai[3] = Rf_dpSavedptr(dd)->mai[3];
    Rf_dpptr(dd)->mUnits = Rf_dpSavedptr(dd)->mUnits;
    Rf_dpptr(dd)->mex = Rf_dpSavedptr(dd)->mex;
    nr = Rf_dpptr(dd)->numrows = Rf_dpSavedptr(dd)->numrows;
    nc = Rf_dpptr(dd)->numcols = Rf_dpSavedptr(dd)->numcols;
    Rf_dpptr(dd)->currentFigure = Rf_dpSavedptr(dd)->currentFigure;
    Rf_dpptr(dd)->lastFigure = Rf_dpSavedptr(dd)->lastFigure;
    for (i = 0; i < nr && i < MAX_LAYOUT_ROWS; i++)
    {
        Rf_dpptr(dd)->heights[i] = Rf_dpSavedptr(dd)->heights[i];
        Rf_dpptr(dd)->cmHeights[i] = Rf_dpSavedptr(dd)->cmHeights[i];
    }
    for (j = 0; j < nc && j < MAX_LAYOUT_COLS; j++)
    {
        Rf_dpptr(dd)->widths[j] = Rf_dpSavedptr(dd)->widths[j];
        Rf_dpptr(dd)->cmWidths[j] = Rf_dpSavedptr(dd)->cmWidths[j];
    }
    for (i = 0; i < nr * nc && i < MAX_LAYOUT_CELLS; i++)
    {
        Rf_dpptr(dd)->order[i] = Rf_dpSavedptr(dd)->order[i];
        Rf_dpptr(dd)->respect[i] = Rf_dpSavedptr(dd)->respect[i];
    }
    Rf_dpptr(dd)->rspct = Rf_dpSavedptr(dd)->rspct;
    Rf_dpptr(dd)->layout = Rf_dpSavedptr(dd)->layout;
    Rf_dpptr(dd)->mfind = Rf_dpSavedptr(dd)->mfind;
    Rf_dpptr(dd)->new = Rf_dpSavedptr(dd)->new;
    Rf_dpptr(dd)->oma[0] = Rf_dpSavedptr(dd)->oma[0];
    Rf_dpptr(dd)->oma[1] = Rf_dpSavedptr(dd)->oma[1];
    Rf_dpptr(dd)->oma[2] = Rf_dpSavedptr(dd)->oma[2];
    Rf_dpptr(dd)->oma[3] = Rf_dpSavedptr(dd)->oma[3];
    Rf_dpptr(dd)->omi[0] = Rf_dpSavedptr(dd)->omi[0];
    Rf_dpptr(dd)->omi[1] = Rf_dpSavedptr(dd)->omi[1];
    Rf_dpptr(dd)->omi[2] = Rf_dpSavedptr(dd)->omi[2];
    Rf_dpptr(dd)->omi[3] = Rf_dpSavedptr(dd)->omi[3];
    Rf_dpptr(dd)->omd[0] = Rf_dpSavedptr(dd)->omd[0];
    Rf_dpptr(dd)->omd[1] = Rf_dpSavedptr(dd)->omd[1];
    Rf_dpptr(dd)->omd[2] = Rf_dpSavedptr(dd)->omd[2];
    Rf_dpptr(dd)->omd[3] = Rf_dpSavedptr(dd)->omd[3];
    Rf_dpptr(dd)->oUnits = Rf_dpSavedptr(dd)->oUnits;
    Rf_dpptr(dd)->plt[0] = Rf_dpSavedptr(dd)->plt[0];
    Rf_dpptr(dd)->plt[1] = Rf_dpSavedptr(dd)->plt[1];
    Rf_dpptr(dd)->plt[2] = Rf_dpSavedptr(dd)->plt[2];
    Rf_dpptr(dd)->plt[3] = Rf_dpSavedptr(dd)->plt[3];
    Rf_dpptr(dd)->pin[0] = Rf_dpSavedptr(dd)->pin[0];
    Rf_dpptr(dd)->pin[1] = Rf_dpSavedptr(dd)->pin[1];
    Rf_dpptr(dd)->pUnits = Rf_dpSavedptr(dd)->pUnits;
    Rf_dpptr(dd)->defaultPlot = Rf_dpSavedptr(dd)->defaultPlot;
    Rf_dpptr(dd)->pty = Rf_dpSavedptr(dd)->pty;
    Rf_dpptr(dd)->usr[0] = Rf_dpSavedptr(dd)->usr[0];
    Rf_dpptr(dd)->usr[1] = Rf_dpSavedptr(dd)->usr[1];
    Rf_dpptr(dd)->usr[2] = Rf_dpSavedptr(dd)->usr[2];
    Rf_dpptr(dd)->usr[3] = Rf_dpSavedptr(dd)->usr[3];
    Rf_dpptr(dd)->logusr[0] = Rf_dpSavedptr(dd)->logusr[0];
    Rf_dpptr(dd)->logusr[1] = Rf_dpSavedptr(dd)->logusr[1];
    Rf_dpptr(dd)->logusr[2] = Rf_dpSavedptr(dd)->logusr[2];
    Rf_dpptr(dd)->logusr[3] = Rf_dpSavedptr(dd)->logusr[3];
}

static SEXP baseCallback(GEevent task, pGEDevDesc dd, SEXP data)
{
    GESystemDesc *sd;
    baseSystemState *bss, *bss2;
    SEXP result = R_NilValue;

    switch (task)
    {
    case GE_FinaliseState:
        sd = dd->gesd[baseRegisterIndex];
        free(sd->systemSpecific);
        sd->systemSpecific = NULL;
        break;
    case GE_InitState: {
        pDevDesc dev;
        GPar *ddp;
        sd = dd->gesd[baseRegisterIndex];
        dev = dd->dev;
        bss = sd->systemSpecific = malloc(sizeof(baseSystemState));
        ddp = &(bss->dp);
        GInit(ddp);
        ddp->cra[0] = dev->cra[0];
        ddp->cra[1] = dev->cra[1];
#if 0
	/* Some things are set by the device, so copy them across now.
	 */
	ddp->ipr[0] = dev->ipr[0];
	ddp->ipr[1] = dev->ipr[1];
	ddp->cra[0] = dev->cra[0];
	ddp->cra[1] = dev->cra[1];
	ddp->asp = dev->asp;
	ddp->left = dev->left;
	ddp->right = dev->right;
	ddp->top = dev->top;
	ddp->bottom = dev->bottom;
	ddp->xCharOffset = dev->xCharOffset;
	ddp->yCharOffset = dev->yCharOffset;
	ddp->yLineBias = dev->yLineBias;
	ddp->canResizePlot = dev->canResizePlot;
	ddp->canChangeFont = dev->canChangeFont;
	ddp->canRotateText = dev->canRotateText;
	ddp->canResizeText = dev->canResizeText;
	ddp->canClip = dev->canClip;
	ddp->canHAdj = dev->canHAdj;
#endif
        /* For some things, the device sets the starting value at least.
         */
        ddp->ps = dev->startps;
        ddp->col = ddp->fg = dev->startcol;
        ddp->bg = dev->startfill;
        ddp->font = dev->startfont;
        ddp->lty = dev->startlty;
        ddp->gamma = dev->startgamma;
        /* Initialise the gp settings too: formerly in addDevice. */
        copyGPar(ddp, &(bss->gp));
        GReset(dd);
        /*
         * The device has not yet received any base output
         */
        bss->baseDevice = FALSE;
        break;
    }
    case GE_CopyState: {
        GEDevDesc *curdd = GEcurrentDevice();
        bss = dd->gesd[baseRegisterIndex]->systemSpecific;
        bss2 = curdd->gesd[baseRegisterIndex]->systemSpecific;
        copyGPar(&(bss->dpSaved), &(bss2->dpSaved));
        restoredpSaved(curdd);
        copyGPar(&(bss2->dp), &(bss2->gp));
        GReset(curdd);
        break;
    }
    case GE_SaveState:
        bss = dd->gesd[baseRegisterIndex]->systemSpecific;
        copyGPar(&(bss->dp), &(bss->dpSaved));
        break;
    case GE_RestoreState:
        bss = dd->gesd[baseRegisterIndex]->systemSpecific;
        restoredpSaved(dd);
        copyGPar(&(bss->dp), &(bss->gp));
        GReset(dd);
        break;
    case GE_SaveSnapshotState: {
        SEXP state;
        bss = dd->gesd[baseRegisterIndex]->systemSpecific;
        PROTECT(state = allocVector(INTSXP,
                                    /* Got this formula from devga.c
                                     * Not sure why the "+ 1"
                                     * Rounding up? Yes!
                                     */
                                    1 + sizeof(GPar) / sizeof(int)));
        copyGPar(&(bss->dpSaved), (GPar *)INTEGER(state));
        result = state;
        UNPROTECT(1);
        break;
    }
    case GE_RestoreSnapshotState:
        bss = dd->gesd[baseRegisterIndex]->systemSpecific;
        copyGPar((GPar *)INTEGER(data), &(bss->dpSaved));
        restoredpSaved(dd);
        copyGPar(&(bss->dp), &(bss->gp));
        GReset(dd);
        break;
    case GE_CheckPlot:
        /* Check that the current plotting state is "valid"
         */
        bss = dd->gesd[baseRegisterIndex]->systemSpecific;
        result = ScalarLogical(bss->baseDevice ? (bss->gp.state == 1) && bss->gp.valid : TRUE);
        break;
    case GE_ScalePS: {
        GPar *ddp, *ddpSaved;
        bss = dd->gesd[baseRegisterIndex]->systemSpecific;
        ddp = &(bss->dp);
        ddpSaved = &(bss->dpSaved);
        if (isReal(data) && LENGTH(data) == 1)
        {
            double rf = REAL(data)[0];
            ddp->scale *= rf;
            ddp->cra[0] *= rf;
            ddp->cra[1] *= rf;
            /* Modify the saved settings so effect display list too */
            ddpSaved->scale *= rf;
            ddpSaved->cra[0] *= rf;
            ddpSaved->cra[1] *= rf;
        }
        else
            error(_("Event GE_ScalePS requires a single numeric value"));
        break;
    }
    }
    return result;
}

/* (un)Register the base graphics system with the graphics engine
 */
void attribute_hidden registerBase(void)
{
    GEregisterSystem(baseCallback, &baseRegisterIndex);
}

void attribute_hidden unregisterBase(void)
{
    GEunregisterSystem(baseRegisterIndex);
}

/* FIXME: Make this a macro to avoid function call overhead?
   Inline it if you really think it matters.
 */
attribute_hidden GPar *Rf_gpptr(pGEDevDesc dd)
{
    baseSystemState *bss = dd->gesd[baseRegisterIndex]->systemSpecific;
    return &(bss->gp);
}

attribute_hidden GPar *Rf_dpptr(pGEDevDesc dd)
{
    baseSystemState *bss = dd->gesd[baseRegisterIndex]->systemSpecific;
    return &(bss->dp);
}

attribute_hidden GPar *Rf_dpSavedptr(pGEDevDesc dd)
{
    baseSystemState *bss = dd->gesd[baseRegisterIndex]->systemSpecific;
    return &(bss->dpSaved);
}

attribute_hidden /* used in GNewPlot */
    void
    Rf_setBaseDevice(Rboolean val, pGEDevDesc dd)
{
    baseSystemState *bss = dd->gesd[baseRegisterIndex]->systemSpecific;
    bss->baseDevice = val;
}
