#include <tcl.h>
#include <windows.h>

void tcltk_init();
#include <R_ext/Rdynload.h>

typedef void (*DL3)();
extern __declspec(dllimport) void (*R_tcldo)();

static void _R_tcldo()
{
    Tcl_ServiceAll();
}

static void (*old_R_tcldo)();

void tcltk_start()
{
    HWND active = GetForegroundWindow(); /* ActiveTCL steals the focus */
    tcltk_init();                        /* won't return on error */
    old_R_tcldo = R_tcldo;
    R_tcldo = &_R_tcldo;
    _R_tcldo();                  /* one call to trigger the focus stealing bug */
    SetForegroundWindow(active); /* and fix it */
}

void tcltk_end()
{
    R_tcldo = old_R_tcldo;
}
