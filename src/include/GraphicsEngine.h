
/* The graphics engine will only accept locations and dimensions 
 * in native device coordinates, but it provides the following functions
 * for converting between a couple of simple alternative coordinate 
 * systems and device coordinates:
 *    DEVICE = native units of the device
 *    NDC = Normalised device coordinates 
 *    INCHES = inches (!)
 *    CM = centimetres (!!)
 */

typedef enum {
 GE_DEVICE	= 0,	/* native device coordinates (rasters) */
 GE_NDC	= 1,	/* normalised device coordinates x=(0,1), y=(0,1) */
 GE_INCHES = 2,
 GE_CM     = 3
} GEUnit;

#define MAX_GRAPHICS_SYSTEMS 24

typedef enum {
    /* In response to this event, the registered graphics system
     * should allocate and initialise the systemSpecific structure
     */
    GE_InitState = 0,
    /* This event gives the registered system a chance to undo
     * anything done in the initialisation.
     */
    GE_FinaliseState = 1,
    /* This is sent by the graphics engine prior to initialising
     * the display list.  It give the graphics system the chance
     * to squirrel away information it will need for redrawing the
     * the display list
     */
    GE_SaveState = 2,
    /* This is sent by the graphics engine prior to replaying the
     * display list.  It gives the graphics system the chance to
     * restore any information it saved on the GE_SaveState event
     */
    GE_RestoreState = 6,
    /* Copy system state information to the current device.
     * This is used when copying graphics from one device to another
     * so all the graphics system needs to do is to copy across
     * the bits required for the display list to draw faithfully
     * on the new device.
     */
    GE_CopyState = 3,
    /* Create a snapshot of the system state that is sufficient
     * for the current "image" to be reproduced
     */
    GE_SaveSnapshotState = 4,
    /* Restore the system state that is saved by GE_SaveSnapshotState
     */
    GE_RestoreSnapshotState = 5,
    /* When replaying the display list, the graphics engine 
     * checks, after each replayed action, that the action 
     * produced valid output.  This is the graphics system's
     * chance to say that the output is crap (in which case the
     * graphics engine will abort the display list replay).
     */
    GE_CheckPlot = 7
} GEevent;

/* The full definition should be ...
 *    typedef SEXP (* GEcallback)(GEvent, *GEDevDesc, SEXP);
 *
 * ... but I could not figure out how to use *GEDevDesc before
 * the definition of GEDevDesc.
 */
typedef SEXP (* GEcallback)();

typedef struct {
    /* An array of information about each graphics system that
     * has registered with the graphics engine.
     * This is used to store graphics state for each graphics
     * system on each device.
     */
    void *systemSpecific;
    /* 
     * An array of function pointers, one per graphics system that
     * has registered with the graphics engine.
     *
     * system_Callback is called when the graphics engine wants
     * to give a graphics system the chance to play with its
     * device-specific information (stored in systemSpecific)
     * There are two parameters:  an "event" to tell the graphics
     * system why the graphics engine has called this function,
     * and the systemSpecific pointer.  The graphics engine
     * has to pass the systemSpecific pointer because only
     * the graphics engine will know what array index to use.
     */
    GEcallback callback;
} GESystemDesc;

typedef struct {
    int newDevStruct;
    NewDevDesc *dev;
    /* Information about a device which has nothing to do with
     * R's concept of a graphics engine.
     */
    GESystemDesc *gesd[MAX_GRAPHICS_SYSTEMS];
} GEDevDesc;

GEDevDesc* GEcreateDevDesc(NewDevDesc* dev);
void GEdestroyDevDesc(GEDevDesc* dd);
void* GEsystemState(GEDevDesc *dd, int index);
void GEregisterWithDevice(GEDevDesc *dd);
int GEregisterSystem(GEcallback callback);
void GEunregisterSystem(int registerIndex);

SEXP GEHandleEvent(GEevent event, NewDevDesc *dev, SEXP data);

double fromDeviceX(double value, GEUnit to, GEDevDesc *dd); 
double toDeviceX(double value, GEUnit from, GEDevDesc *dd);
double fromDeviceY(double value, GEUnit to, GEDevDesc *dd); 
double toDeviceY(double value, GEUnit from, GEDevDesc *dd);
double fromDeviceWidth(double value, GEUnit to, GEDevDesc *dd); 
double toDeviceWidth(double value, GEUnit from, GEDevDesc *dd);
double fromDeviceHeight(double value, GEUnit to, GEDevDesc *dd); 
double toDeviceHeight(double value, GEUnit from, GEDevDesc *dd);

/*
 *	Some Notes on Line Textures
 *
 *	Line textures are stored as an array of 4-bit integers within
 *	a single 32-bit word.  These integers contain the lengths of
 *	lines to be drawn with the pen alternately down and then up.
 *	The device should try to arrange that these values are measured
 *	in points if possible, although pixels is ok on most displays.
 *
 *	If newlty contains a line texture description it is decoded
 *	as follows:
 *
 *		ndash = 0;
 *		for(i=0 ; i<8 && newlty & 15 ; i++) {
 *			dashlist[ndash++] = newlty & 15;
 *			newlty = newlty>>4;
 *		}
 *		dashlist[0] = length of pen-down segment
 *		dashlist[1] = length of pen-up segment
 *		etc
 *
 *	An integer containing a zero terminates the pattern.  Hence
 *	ndash in this code fragment gives the length of the texture
 *	description.  If a description contains an odd number of
 *	elements it is replicated to create a pattern with an
 *	even number of elements.  (If this is a pain, do something
 *	different its not crucial).
 *
 */

/*--- The basic numbered & names line types; Here device-independent:
  e.g. "dashed" == "44",  "dotdash" == "1343"
*/

#define LTY_BLANK	-1
#define LTY_SOLID	0
#define LTY_DASHED	4 + (4<<4)
#define LTY_DOTTED	1 + (3<<4)
#define LTY_DOTDASH	1 + (3<<4) + (4<<8) + (3<<12)
#define LTY_LONGDASH	7 + (3<<4)
#define LTY_TWODASH	2 + (2<<4) + (6<<8) + (2<<12)

void GESetClip(double x1, double y1, double x2, double y2, GEDevDesc *dd);
void GENewPage(int fill, GEDevDesc *dd);
void GELine(double x1, double y1, double x2, double y2, 
	    int col, double gamma, int lty, double lwd,
	    GEDevDesc *dd);
void GEPolyline(int n, double *x, double *y, 
		int col, double gamma, int lty, double lwd,
		GEDevDesc *dd);
void GEPolygon(int n, double *x, double *y, 
	       int col, int fill, double gamma, int lty, double lwd,
	       GEDevDesc *dd);
void GECircle(double x, double y, double radius,
	     int col, int fill, double gamma, int lty, double lwd,
	     GEDevDesc *dd);
void GERect(double x0, double y0, double x1, double y1,
	    int col, int fill, double gamma, int lty, double lwd,
	    GEDevDesc *dd);
void GEText(double x, double y, char *str,
	    double xc, double yc, double rot, 
	    int col, double gamma, int font, double cex, double ps,
	    GEDevDesc *dd);
void GEMode(int mode, GEDevDesc* dd);
void GESymbol(double x, double y, int pch, double size,
	      int col, int fill, double gamma, double lty, double lwd,
	      int font, double cex, double ps,
	      GEDevDesc *dd);
void GEPretty(double *lo, double *up, int *ndiv);
void GEMetricInfo(int c, int font, double cex, double ps,
		  double *ascent, double *descent, double *width,
		  GEDevDesc *dd);
double GEStrWidth(char *str, int font, double cex, double ps, GEDevDesc *dd);
double GEStrHeight(char *str, int font, double cex, double ps, GEDevDesc *dd);

#define	DEG2RAD 0.01745329251994329576

GEDevDesc* GEcurrentDevice();
void GEinitDisplayList(GEDevDesc *dd);
void GEplayDisplayList(GEDevDesc *dd);
void GEcopyDisplayList(int fromDevice);
SEXP GEcreateSnapshot(GEDevDesc *dd);
void GEplaySnapshot(SEXP snapshot, GEDevDesc* dd);
