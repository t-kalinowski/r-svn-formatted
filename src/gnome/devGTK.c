#include <gnome.h>

#include "Defn.h"
#include "Graphics.h"
#include "devGTK.h"
#include "terminal.h"

#define CURSOR GDK_CROSSHAIR /* Default cursor */
#define MM_PER_INCH 25.4     /* mm -> inch conversion */

typedef struct
{
    /* R Graphics Parameters */
    /* local device copy so that we can detect */
    /* when parameter changes */

    double cex; /* Character expansion */
    double srt; /* String rotation */

    gint bg; /* Background */

    int fontface; /* Typeface */
    int fontsize; /* Size in points */

    gint lty, lwd; /* line params */

    /* GTK Driver Specific */

    int windowWidth;   /* Window width (pixels) */
    int windowHeight;  /* Window height (pixels) */
    int resize;        /* Window resized */
    GtkWidget *window; /* Graphics Window */
    GtkWidget *drawing;

    GdkGC *wgc;
    GdkColor *gcol_bg;
    GdkRectangle clip;
    GdkCursor *gcursor;

    int usefixed;
    GdkFont *fixedfont;
    GdkFont *font;

} gtkDesc;

static int numGTKDevices = 0;

/* Device driver actions */
static void GTK_Activate(DevDesc *);
static void GTK_Circle(double, double, int, double, int, int, DevDesc *);
static void GTK_Clip(double, double, double, double, DevDesc *);
static void GTK_Close(DevDesc *);
static void GTK_Deactivate(DevDesc *);
static void GTK_Hold(DevDesc *);
static void GTK_Line(double, double, double, double, int, DevDesc *);
static int GTK_Locator(double *, double *, DevDesc *);
static void GTK_Mode(int);
static void GTK_NewPage(DevDesc *);
static int GTK_Open(DevDesc *, gtkDesc *, char *, double, double);
static void GTK_Polygon(int, double *, double *, int, int, int, DevDesc *);
static void GTK_Polyline(int, double *, double *, int, DevDesc *);
static void GTK_Rect(double, double, double, double, int, int, int, DevDesc *);
static void GTK_Resize(DevDesc *);
static double GTK_StrWidth(char *, DevDesc *);
static void GTK_Text(double, double, int, char *, double, double, double, DevDesc *);
static void GTK_MetricInfo(int, double *, double *, double *, DevDesc *);

/* Pixel Dimensions (Inches) */

static double pixelWidth(void)
{
    double width, widthMM;
    width = gdk_screen_width();
    widthMM = gdk_screen_width_mm();
    return ((double)widthMM / (double)width) / MM_PER_INCH;
}

static double pixelHeight(void)
{
    double height, heightMM;
    height = gdk_screen_height();
    heightMM = gdk_screen_height_mm();
    return ((double)heightMM / (double)height) / MM_PER_INCH;
}

/* font stuff */

/* set the r, g, b, and pixel values of gcol to color */
static void SetColor(GdkColor *gcol, int color)
{
    int red, green, blue;

    red = R_RED(color);
    green = R_GREEN(color);
    blue = R_BLUE(color);
    gcol->red = red * 257;
    gcol->green = green * 257;
    gcol->blue = blue * 257;
    gcol->pixel = gdk_rgb_xpixel_from_rgb((red << 16) | (green << 8) | (blue));
}

/* set the line type */
static void SetLineType(DevDesc *dd, int newlty, int newlwd)
{
    static gchar dashlist[8];
    gint i, j;
    gtkDesc *gtkd = (gtkDesc *)dd->deviceSpecific;

    if (newlty != gtkd->lty || newlwd != gtkd->lwd)
    {
        gtkd->lty = newlty;
        gtkd->lwd = newlwd;

        if (newlty == 0)
        {
            if (newlwd <= 1)
                newlwd = 0;

            gdk_gc_set_line_attributes(gtkd->wgc, newlwd, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_ROUND);
        }
        else
        {
            if (newlwd < 1)
                newlwd = 1;

            for (i = 0; (i < 8) && (newlty != 0); i++)
            {
                j = newlty & 15;

                if (j == 0)
                    j = 1;

                j = j * newlwd;

                if (j > 255)
                    j = 255;

                dashlist[i] = j;
                newlty = newlty >> 4;
            }

            /* set dashes */
            gdk_gc_set_dashes(gtkd->wgc, 0, dashlist, i);
            gdk_gc_set_line_attributes(gtkd->wgc, newlwd, GDK_LINE_ON_OFF_DASH, GDK_CAP_BUTT, GDK_JOIN_ROUND);
        }
    }
}

/* signal functions */

static gint realize_event(GtkWidget *widget, gpointer data)
{
    DevDesc *dd;
    gtkDesc *gtkd;

    dd = (DevDesc *)data;
    g_return_val_if_fail(dd != NULL, FALSE);

    gtkd = (gtkDesc *)dd->deviceSpecific;
    g_return_val_if_fail(gtkd != NULL, FALSE);
    g_return_val_if_fail(gtkd->drawing != NULL, FALSE);
    g_return_val_if_fail(GTK_IS_DRAWING_AREA(gtkd->drawing), FALSE);

    /* create gc */
    gtkd->wgc = gdk_gc_new(gtkd->drawing->window);

    /* set the cursor */
    gtkd->gcursor = gdk_cursor_new(GDK_CROSSHAIR);
    gdk_window_set_cursor(gtkd->drawing->window, gtkd->gcursor);

    /* set window bg */
    gdk_window_set_background(gtkd->drawing->window, gtkd->gcol_bg);

    return FALSE;
}

static gint configure_event(GtkWidget *widget, GdkEventConfigure *event, gpointer data)
{
    DevDesc *dd;
    gtkDesc *gtkd;

    dd = (DevDesc *)data;
    g_return_val_if_fail(dd != NULL, FALSE);

    gtkd = (gtkDesc *)dd->deviceSpecific;
    g_return_val_if_fail(gtkd != NULL, FALSE);
    g_return_val_if_fail(gtkd->drawing != NULL, FALSE);
    g_return_val_if_fail(GTK_IS_DRAWING_AREA(gtkd->drawing), FALSE);

    /* check for resize */
    if ((GTK_WIDGET_REALIZED(gtkd->drawing)) &&
        ((gtkd->windowWidth != event->width) || (gtkd->windowHeight != event->height)))
    {
        gtkd->windowWidth = event->width;
        gtkd->windowHeight = event->height;

        gtkd->resize = 1;
    }

    return FALSE;
}

static gint expose_event(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
    DevDesc *dd;
    gtkDesc *gtkd;

    dd = (DevDesc *)data;
    g_return_val_if_fail(dd != NULL, FALSE);

    gtkd = (gtkDesc *)dd->deviceSpecific;
    g_return_val_if_fail(gtkd != NULL, FALSE);
    g_return_val_if_fail(gtkd->drawing != NULL, FALSE);
    g_return_val_if_fail(GTK_IS_DRAWING_AREA(gtkd->drawing), FALSE);

    if (gtkd->resize != 0)
    {
        dd->dp.resize(dd);
    }

    playDisplayList(dd);

    return FALSE;
}

static gint delete_event(GtkWidget *widget, GdkEvent *event, gpointer data)
{
    DevDesc *dd;
    gtkDesc *gtkd;

    dd = (DevDesc *)data;
    g_return_val_if_fail(dd != NULL, FALSE);

    KillDevice(dd);

    return TRUE;
}

static void tb_activate_cb(GtkWidget *widget, gpointer data)
{
    DevDesc *dd;

    dd = (DevDesc *)data;
    g_return_if_fail(dd != NULL);

    selectDevice(deviceNumber(dd));
}

static void tb_close_cb(GtkWidget *widget, gpointer data)
{
    DevDesc *dd;

    dd = (DevDesc *)data;
    g_return_if_fail(dd != NULL);

    KillDevice(dd);
}

static GnomeUIInfo graphics_toolbar[] = {
    {GNOME_APP_UI_ITEM, "Activate", "Make this window the current device", tb_activate_cb, NULL, NULL,
     GNOME_APP_PIXMAP_STOCK, GNOME_STOCK_PIXMAP_JUMP_TO, 0, (GdkModifierType)0, NULL},
    GNOMEUIINFO_SEPARATOR,
    {GNOME_APP_UI_ITEM, "Save As", "Save as a PS file", NULL, NULL, NULL, GNOME_APP_PIXMAP_STOCK,
     GNOME_STOCK_PIXMAP_SAVE_AS, 0, (GdkModifierType)0, NULL},
    {GNOME_APP_UI_ITEM, "Print", "Print graphics", NULL, NULL, NULL, GNOME_APP_PIXMAP_STOCK, GNOME_STOCK_PIXMAP_PRINT,
     0, (GdkModifierType)0, NULL},
    GNOMEUIINFO_SEPARATOR,
    {GNOME_APP_UI_ITEM, "Close", "Close this graphics device", tb_close_cb, NULL, NULL, GNOME_APP_PIXMAP_STOCK,
     GNOME_STOCK_PIXMAP_CLOSE, 0, (GdkModifierType)0, NULL},
    NULL};

/* create window etc */
static int GTK_Open(DevDesc *dd, gtkDesc *gtkd, char *dsp, double w, double h)
{
    GdkColor bg;
    GtkStyle *wstyle;

    gint iw, ih, result;

    /*gdk_rgb_set_install(TRUE); */
    gdk_rgb_set_verbose(TRUE);
    gdk_rgb_init();
    gtk_widget_push_visual(gdk_rgb_get_visual());
    gtk_widget_push_colormap(gdk_rgb_get_cmap());

    /* initialise pointers */
    gtkd->drawing = NULL;
    gtkd->wgc = NULL;
    gtkd->gcursor = NULL;

    /* FIXME: SetBaseFont */

    /* create window etc */
    gtkd->windowWidth = iw = w / pixelWidth();
    gtkd->windowHeight = ih = h / pixelHeight();

    gtkd->window = gnome_app_new("R.gnome.graphics", "R Graphics");

    gtk_window_set_policy(GTK_WINDOW(gtkd->window), TRUE, TRUE, FALSE);
    gtk_widget_realize(gtkd->window);

    /* create toolbar */
    gnome_app_create_toolbar_with_data(GNOME_APP(gtkd->window), graphics_toolbar, (gpointer)dd);

    /* create drawingarea */
    gtkd->drawing = gtk_drawing_area_new();

    /* connect to signal handlers, etc */
    gtk_signal_connect(GTK_OBJECT(gtkd->drawing), "realize", (GtkSignalFunc)realize_event, (gpointer)dd);

    /* drawingarea properties */
    gtk_widget_set_usize(gtkd->drawing, iw, ih);

    /* setup background color */
    gtkd->bg = dd->dp.bg = R_RGB(255, 255, 255);
    gtkd->gcol_bg = gdk_color_copy(&bg);
    SetColor(gtkd->gcol_bg, gtkd->bg);

    /* place and realize the drawing area */
    gnome_app_set_contents(GNOME_APP(gtkd->window), gtkd->drawing);
    gtk_widget_realize(gtkd->drawing);

    /* connect to signal handlers, etc */
    gtk_signal_connect(GTK_OBJECT(gtkd->drawing), "configure_event", (GtkSignalFunc)configure_event, (gpointer)dd);
    gtk_signal_connect(GTK_OBJECT(gtkd->drawing), "expose_event", (GtkSignalFunc)expose_event, (gpointer)dd);
    gtk_signal_connect(GTK_OBJECT(gtkd->window), "delete_event", (GtkSignalFunc)delete_event, (gpointer)dd);

    /* show everything */
    gtk_widget_show_all(gtkd->window);

    /* let other widgets use the default settings */
    gtk_widget_pop_visual();
    gtk_widget_pop_colormap();

    /* initialise line params */
    gtkd->lty = -1;
    gtkd->lwd = -1;

    /* we made it! */
    return 1;
}

static double GTK_StrWidth(char *str, DevDesc *dd)
{
    return 6.0;
}

static void GTK_MetricInfo(int c, double *ascent, double *descent, double *width, DevDesc *dd)
{
}

/* set clipping */
static void GTK_Clip(double x0, double x1, double y0, double y1, DevDesc *dd)
{
    gtkDesc *gtkd = (gtkDesc *)dd->deviceSpecific;

    gtkd->clip.x = MIN(x0, x1);
    gtkd->clip.width = abs(x0 - x1);

    gtkd->clip.y = MIN(y0, y1);
    gtkd->clip.height = abs(y0 - y1);

    gdk_gc_set_clip_rectangle(gtkd->wgc, &gtkd->clip);
}

static void GTK_Resize(DevDesc *dd)
{
    gtkDesc *gtkd = (gtkDesc *)dd->deviceSpecific;

    if (gtkd->resize != 0)
    {
        dd->dp.left = dd->gp.left = 0.0;
        dd->dp.right = dd->gp.right = gtkd->windowWidth;
        dd->dp.bottom = dd->gp.bottom = gtkd->windowHeight;
        dd->dp.top = dd->gp.top = 0.0;
        gtkd->resize = 0;
    }
}

/* clear the drawing area */
static void GTK_NewPage(DevDesc *dd)
{
    gtkDesc *gtkd;

    g_return_if_fail(dd != NULL);

    gtkd = (gtkDesc *)dd->deviceSpecific;
    g_return_if_fail(gtkd != NULL);
    g_return_if_fail(gtkd->drawing != NULL);
    g_return_if_fail(GTK_IS_DRAWING_AREA(gtkd->drawing));

    if (gtkd->bg != dd->dp.bg)
    {
        SetColor(gtkd->gcol_bg, dd->dp.bg);
        gtkd->bg = dd->dp.bg;
        gdk_window_set_background(gtkd->drawing->window, gtkd->gcol_bg);
    }

    gdk_window_clear(gtkd->drawing->window);
}

/* kill off the window etc */
static void GTK_Close(DevDesc *dd)
{
    gint i, j;
    gtkDesc *gtkd = (gtkDesc *)dd->deviceSpecific;

    gtk_widget_destroy(gtkd->window);

    numGTKDevices--;

    free(gtkd);
}

#define title_text_inactive "R graphics device %d"
#define title_text_active "R graphics device %d - Active"

static void GTK_Activate(DevDesc *dd)
{
    gtkDesc *gtkd;
    gint devnum;
    gchar *title_text;

    gtkd = (gtkDesc *)dd->deviceSpecific;
    g_return_if_fail(gtkd != NULL);

    devnum = deviceNumber(dd);
    devnum++;

    title_text = g_strdup_printf(title_text_active, devnum);

    gtk_window_set_title(GTK_WINDOW(gtkd->window), title_text);

    g_free(title_text);
}

static void GTK_Deactivate(DevDesc *dd)
{
    gtkDesc *gtkd;
    gint devnum;
    gchar *title_text;

    gtkd = (gtkDesc *)dd->deviceSpecific;
    g_return_if_fail(gtkd != NULL);

    devnum = deviceNumber(dd);
    devnum++;

    title_text = g_strdup_printf(title_text_active, devnum);

    gtk_window_set_title(GTK_WINDOW(gtkd->window), title_text);

    g_free(title_text);
}

/* drawing stuff */

static void GTK_Rect(double x0, double y0, double x1, double y1, int coords, int bg, int fg, DevDesc *dd)
{
    double tmp;
    gtkDesc *gtkd = (gtkDesc *)dd->deviceSpecific;
    GdkColor gcol_fill, gcol_outline;

    GConvert(&x0, &y0, coords, DEVICE, dd);
    GConvert(&x1, &y1, coords, DEVICE, dd);

    if (x0 > x1)
    {
        tmp = x0;
        x0 = x1;
        x1 = tmp;
    }
    if (y0 > y1)
    {
        tmp = y0;
        y0 = y1;
        y1 = tmp;
    }

    if (bg != NA_INTEGER)
    {
        SetColor(&gcol_fill, bg);
        gdk_gc_set_foreground(gtkd->wgc, &gcol_fill);

        SetLineType(dd, dd->gp.lty, dd->gp.lwd);

        gdk_draw_rectangle(gtkd->drawing->window, gtkd->wgc, TRUE, (gint)x0, (gint)y0, (gint)x1 - (gint)x0,
                           (gint)y1 - (gint)y0);
    }
    if (fg != NA_INTEGER)
    {
        SetColor(&gcol_outline, fg);
        gdk_gc_set_foreground(gtkd->wgc, &gcol_outline);

        SetLineType(dd, dd->gp.lty, dd->gp.lwd);

        gdk_draw_rectangle(gtkd->drawing->window, gtkd->wgc, FALSE, (gint)x0, (gint)y0, (gint)x1 - (gint)x0,
                           (gint)y1 - (gint)y0);
    }
}

static void GTK_Circle(double x, double y, int coords, double r, int col, int border, DevDesc *dd)
{
    GdkColor gcol_fill, gcol_outline;
    gint ix, iy, ir;
    gtkDesc *gtkd = (gtkDesc *)dd->deviceSpecific;

    GConvert(&x, &y, coords, DEVICE, dd);

    ix = x - r;
    iy = y - r;
    ir = 2 * floor(r + 0.5);

    if (col != NA_INTEGER)
    {
        SetColor(&gcol_fill, col);
        gdk_gc_set_foreground(gtkd->wgc, &gcol_fill);

        gdk_draw_arc(gtkd->drawing->window, gtkd->wgc, TRUE, ix, iy, ir, ir, 0, 23040);
    }
    if (border != NA_INTEGER)
    {
        SetColor(&gcol_outline, border);
        gdk_gc_set_foreground(gtkd->wgc, &gcol_outline);

        SetLineType(dd, dd->gp.lty, dd->gp.lwd);

        gdk_draw_arc(gtkd->drawing->window, gtkd->wgc, FALSE, ix, iy, ir, ir, 0, 23040);
    }
}

static void GTK_Line(double x1, double y1, double x2, double y2, int coords, DevDesc *dd)
{
    gtkDesc *gtkd = (gtkDesc *)dd->deviceSpecific;
    GdkColor gcol_fill;
    gint ix1, iy1, ix2, iy2;

    GConvert(&x1, &y1, coords, DEVICE, dd);
    GConvert(&x2, &y2, coords, DEVICE, dd);

    ix1 = (gint)x1;
    iy1 = (gint)y1;
    ix2 = (gint)x2;
    iy2 = (gint)y2;

    if (dd->gp.col != NA_INTEGER)
    {
        SetColor(&gcol_fill, dd->gp.col);
        gdk_gc_set_foreground(gtkd->wgc, &gcol_fill);

        SetLineType(dd, dd->gp.lty, dd->gp.lwd);

        gdk_draw_line(gtkd->drawing->window, gtkd->wgc, ix1, iy1, ix2, iy2);
    }
}

static void GTK_Polyline(int n, double *x, double *y, int coords, DevDesc *dd)
{
    gtkDesc *gtkd = (gtkDesc *)dd->deviceSpecific;
    GdkColor gcol_fill;
    GdkPoint *points;
    double devx, devy;
    int i;

    points = g_new0(GdkPoint, n);

    for (i = 0; i < n; i++)
    {
        devx = x[i];
        devy = y[i];
        GConvert(&devx, &devy, coords, DEVICE, dd);
        points[i].x = (gint16)devx;
        points[i].y = (gint16)devy;
    }

    if (dd->gp.col != NA_INTEGER)
    {
        SetColor(&gcol_fill, dd->gp.col);
        gdk_gc_set_foreground(gtkd->wgc, &gcol_fill);

        SetLineType(dd, dd->gp.lty, dd->gp.lwd);

        gdk_draw_lines(gtkd->drawing->window, gtkd->wgc, points, n);
    }

    g_free(points);
}

static void GTK_Polygon(int n, double *x, double *y, int coords, int bg, int fg, DevDesc *dd)
{
    gtkDesc *gtkd = (gtkDesc *)dd->deviceSpecific;
    GdkColor gcol_fill, gcol_outline;
    GdkPoint *points;
    double devx, devy;
    int i;

    points = g_new0(GdkPoint, n + 1);

    for (i = 0; i < n; i++)
    {
        devx = x[i];
        devy = y[i];
        GConvert(&devx, &devy, coords, DEVICE, dd);
        points[i].x = (gint16)devx;
        points[i].y = (gint16)devy;
    }

    if (bg != NA_INTEGER)
    {
        SetColor(&gcol_fill, bg);
        gdk_gc_set_foreground(gtkd->wgc, &gcol_fill);

        gdk_draw_polygon(gtkd->drawing->window, gtkd->wgc, TRUE, points, n);
    }
    if (fg != NA_INTEGER)
    {
        SetColor(&gcol_outline, fg);
        gdk_gc_set_foreground(gtkd->wgc, &gcol_outline);

        SetLineType(dd, dd->gp.lty, dd->gp.lwd);

        gdk_draw_polygon(gtkd->drawing->window, gtkd->wgc, FALSE, points, n);
    }

    g_free(points);
}

extern double deg2rad; /* in devGTK.c */

static void GTK_Text(double x, double y, int coords, char *str, double xc, double yc, double rot, DevDesc *dd)
{
    GnomeCanvasItem *item;
    gtkDesc *gtkd = (gtkDesc *)dd->deviceSpecific;
    GdkColor gcol_fill;

    g_message("text");

    GConvert(&x, &y, coords, DEVICE, dd);

    SetColor(&gcol_fill, dd->gp.col);
}

static int GTK_Locator(double *x, double *y, DevDesc *dd)
{
    g_message("locator");
}

static void GTK_Mode(gint mode)
{
}

static void GTK_Hold(DevDesc *dd)
{
}

/* Device driver entry point */
int X11DeviceDriver(DevDesc *dd, char *display, double width, double height, double pointsize)
{
    int ps;
    gtkDesc *gtkd;

    if (!(gtkd = (gtkDesc *)malloc(sizeof(gtkDesc))))
        return 0;

    dd->deviceSpecific = (void *)gtkd;

    /* FIXME: font loading */
    ps = pointsize;
    if (ps < 6 || ps > 24)
        ps = 12;
    ps = 2 * (ps / 2);
    gtkd->fontface = -1;
    gtkd->fontsize = -1;
    dd->dp.font = 1;
    dd->dp.ps = ps;

    /* device driver start */
    if (!GTK_Open(dd, gtkd, display, width, height))
    {
        free(gtkd);
        return 0;
    }

    /* setup data structure */
    dd->dp.open = GTK_Open;
    dd->dp.close = GTK_Close;
    dd->dp.activate = GTK_Activate;
    dd->dp.deactivate = GTK_Deactivate;
    dd->dp.resize = GTK_Resize;
    dd->dp.newPage = GTK_NewPage;
    dd->dp.clip = GTK_Clip;
    dd->dp.strWidth = GTK_StrWidth;
    dd->dp.text = GTK_Text;
    dd->dp.rect = GTK_Rect;
    dd->dp.circle = GTK_Circle;
    dd->dp.line = GTK_Line;
    dd->dp.polyline = GTK_Polyline;
    dd->dp.polygon = GTK_Polygon;
    dd->dp.locator = GTK_Locator;
    dd->dp.mode = GTK_Mode;
    dd->dp.hold = GTK_Hold;
    dd->dp.metricInfo = GTK_MetricInfo;

    dd->dp.left = 0;
    dd->dp.right = gtkd->windowWidth;
    dd->dp.bottom = gtkd->windowHeight;
    dd->dp.top = 0;

    /* FIXME: nominal character sizes */
    dd->dp.cra[0] = 10;
    dd->dp.cra[1] = 10;

    /* FIXME: character addressing offsets */
    dd->dp.xCharOffset = 0.4900;
    dd->dp.yCharOffset = 0.3333;
    dd->dp.yLineBias = 0.1;

    /* inches per raster unit */
    dd->dp.ipr[0] = pixelWidth();
    dd->dp.ipr[1] = pixelHeight();

    /* device capabilities */
    dd->dp.canResizePlot = 1;
    dd->dp.canChangeFont = 0;
    dd->dp.canRotateText = 1;
    dd->dp.canResizeText = 1;
    dd->dp.canClip = 0;

    // x11 device description stuff
    gtkd->cex = 1.0;
    gtkd->srt = 0.0;
    gtkd->resize = 0;

    dd->displayListOn = 1;

    // finish
    return 1;
}
