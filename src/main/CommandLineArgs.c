/*
  R : A Computer Language for Statistical Data Analysis
  Copyright (C) 1997-2005   Robert Gentleman, Ross Ihaka
                            and the R Development Core Team

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or (at
  your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307,
  U.S.A.
 */

/* <UTF8> char here is handled as a whole string,
   or the strings are ASCII (as in sizes).
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <Defn.h>
#include <R_ext/RStartup.h>

/* Remove and process common command-line arguments
 *  Formally part of ../unix/sys-common.c.
 */

/*
  This copies the command line arguments to the Rstart
  structure. The memory is obtained from calloc, etc.
  since these are permanent and it is not intended that
  they be modified. This is why they are copied before
  being processed and removed from the list.

  We might store these as a SEXP. I have no strong opinion
  about this.
 */

/* Permanent copy of the command line arguments and the number
   of them passed to the application.
   These are populated via the routine R_set_command_line_arguments().
*/
static int NumCommandLineArgs = 0;
static char **CommandLineArgs = NULL;

void R_set_command_line_arguments(int argc, char **argv)
{
    int i;

    NumCommandLineArgs = argc;
    CommandLineArgs = (char **)calloc(argc, sizeof(char *));

    for (i = 0; i < argc; i++)
        CommandLineArgs[i] = strdup(argv[i]);
}

/*
  The .Internal which returns the command line arguments that are stored
  in global variables.
 */
SEXP do_commandArgs(SEXP call, SEXP op, SEXP args, SEXP env)
{
    int i;
    SEXP vals;

    vals = allocVector(STRSXP, NumCommandLineArgs);
    for (i = 0; i < NumCommandLineArgs; i++)
        SET_STRING_ELT(vals, i, mkChar(CommandLineArgs[i]));
    return (vals);
}

void R_common_command_line(int *pac, char **argv, Rstart Rp)
{
    int ac = *pac, newac = 1; /* argv[0] is process name */
    int ierr;
    long lval;
    R_size_t value;
    char *p, **av = argv, msg[1024];
    Rboolean processing = TRUE;

    R_RestoreHistory = 1;
    while (--ac)
    {
        if (processing && **++av == '-')
        {
            if (!strcmp(*av, "--version"))
            {
                PrintVersion(msg);
                R_ShowMessage(msg);
                exit(0);
            }
            else if (!strcmp(*av, "--args"))
            {
                /* copy this through for further processing */
                argv[newac++] = *av;
                processing = FALSE;
            }
            else if (!strcmp(*av, "--save"))
            {
                Rp->SaveAction = SA_SAVE;
            }
            else if (!strcmp(*av, "--no-save"))
            {
                Rp->SaveAction = SA_NOSAVE;
            }
            else if (!strcmp(*av, "--restore"))
            {
                Rp->RestoreAction = SA_RESTORE;
            }
            else if (!strcmp(*av, "--no-restore"))
            {
                Rp->RestoreAction = SA_NORESTORE;
                R_RestoreHistory = 0;
            }
            else if (!strcmp(*av, "--no-restore-data"))
            {
                Rp->RestoreAction = SA_NORESTORE;
            }
            else if (!strcmp(*av, "--no-restore-history"))
            {
                R_RestoreHistory = 0;
            }
            else if (!strcmp(*av, "--silent") || !strcmp(*av, "--quiet") || !strcmp(*av, "-q"))
            {
                Rp->R_Quiet = TRUE;
            }
            else if (!strcmp(*av, "--vanilla"))
            {
                Rp->SaveAction = SA_NOSAVE;       /* --no-save */
                Rp->RestoreAction = SA_NORESTORE; /* --no-restore */
                Rp->LoadSiteFile = FALSE;         /* --no-site-file */
                Rp->LoadInitFile = FALSE;         /* --no-init-file */
                R_RestoreHistory = 0;             /* --no-restore-history */
                Rp->NoRenviron = TRUE;
            }
            else if (!strcmp(*av, "--no-environ"))
            {
                Rp->NoRenviron = TRUE;
            }
            else if (!strcmp(*av, "--verbose"))
            {
                Rp->R_Verbose = TRUE;
            }
            else if (!strcmp(*av, "--slave") || !strcmp(*av, "-s"))
            {
                Rp->R_Quiet = TRUE;
                Rp->R_Slave = TRUE;
                Rp->SaveAction = SA_NOSAVE;
            }
            else if (!strcmp(*av, "--no-site-file"))
            {
                Rp->LoadSiteFile = FALSE;
            }
            else if (!strcmp(*av, "--no-init-file"))
            {
                Rp->LoadInitFile = FALSE;
            }
            else if (!strcmp(*av, "--debug-init"))
            {
                Rp->DebugInitFile = TRUE;
            }
            else if (!strncmp(*av, "--encoding", 10))
            {
                if (strlen(*av) < 12)
                {
                    ac--;
                    av++;
                    p = *av;
                }
                else
                    p = &(*av)[11];
                if (p == NULL)
                {
                    R_ShowMessage(_("WARNING: no value given for --encoding given\n"));
                }
                else
                {
                    strncpy(R_StdinEnc, p, 30);
                    R_StdinEnc[31] = '\0';
                    break;
                }
            }
            else if (!strcmp(*av, "-save") || !strcmp(*av, "-nosave") || !strcmp(*av, "-restore") ||
                     !strcmp(*av, "-norestore") || !strcmp(*av, "-noreadline") || !strcmp(*av, "-quiet") ||
                     !strcmp(*av, "-nsize") || !strcmp(*av, "-vsize") || !strcmp(*av, "-V") || !strcmp(*av, "-n") ||
                     !strcmp(*av, "-v"))
            {
                snprintf(msg, 1024, _("WARNING: option '%s' no longer supported\n"), *av);
                R_ShowMessage(msg);
            }
            /* mop up --max/min/-n/vsize */
            else if (strncmp(*av + 7, "size", 4) == 0)
            {
                if (strlen(*av) < 13)
                {
                    ac--;
                    av++;
                    p = *av;
                }
                else
                    p = &(*av)[12];
                if (p == NULL)
                {
                    snprintf(msg, 1024, _("WARNING: no value given for '%s'\n"), *av);
                    R_ShowMessage(msg);
                    break;
                }
                value = R_Decode2Long(p, &ierr);
                if (ierr)
                {
                    if (ierr < 0)
                        snprintf(msg, 1024, _("WARNING: '%s' value is invalid: ignored\n"), *av);
                    else
                        sprintf(msg, _("WARNING: %s=%lu'%c': too large and ignored\n"), *av, (unsigned long)value,
                                (ierr == 1) ? 'M' : ((ierr == 2) ? 'K' : 'k'));
                    R_ShowMessage(msg);
                }
                else
                {
                    if (!strncmp(*av, "--min-nsize", 11))
                        Rp->nsize = value;
                    if (!strncmp(*av, "--max-nsize", 11))
                        Rp->max_nsize = value;
                    if (!strncmp(*av, "--min-vsize", 11))
                        Rp->vsize = value;
                    if (!strncmp(*av, "--max-vsize", 11))
                        Rp->max_vsize = value;
                }
            }
            else if (strncmp(*av, "--max-ppsize", 12) == 0)
            {
                if (strlen(*av) < 14)
                {
                    ac--;
                    av++;
                    p = *av;
                }
                else
                    p = &(*av)[13];
                if (p == NULL)
                {
                    R_ShowMessage(_("WARNING: no value given for '--max-ppsize'\n"));
                    break;
                }
                lval = strtol(p, &p, 10);
                if (lval < 0)
                    R_ShowMessage(_("WARNING: '-max-ppsize' value is negative: ignored\n"));
                else if (lval < 10000)
                    R_ShowMessage(_("WARNING: '-max-ppsize' value is too small: ignored\n"));

                else if (lval > 100000)
                    R_ShowMessage(_("WARNING: '-max-ppsize' value is too large: ignored\n"));
                else
                    Rp->ppsize = lval;
            }
#if 0
	    else if(strncmp(*av, "--vsize", 7) == 0) {
		if(strlen(*av) < 9) {
		    ac--; av++; p = *av;
		}
		else
		    p = &(*av)[8];
		if (p == NULL) {
		    R_ShowMessage(_("WARNING: no 'vsize' given\n"));
		    break;
		}
		value = R_Decode2Long(p, &ierr);
		if(ierr) {
		    if(ierr < 0) /* R_common_badargs(); */
			sprintf(msg, _("WARNING: '--vsize' value is invalid: ignored\n"));
		    else
			sprintf(msg, _("WARNING: --vsize=%ld'%c': too large and ignored\n"),
				value,
				(ierr == 1) ? 'M': ((ierr == 2) ? 'K' : 'k'));
		    R_ShowMessage(msg);

		} else
		    Rp->vsize = value;
	    }
	    else if(strncmp(*av, "--nsize", 7) == 0) {
		if(strlen(*av) < 9) {
		    ac--; av++; p = *av;
		}
		else
		    p = &(*av)[8];
		if (p == NULL) {
		    R_ShowMessage(_("WARNING: no 'nsize' given\n"));
		    break;
		}
		value = R_Decode2Long(p, &ierr);
		if(ierr) {
		    if(ierr < 0) /* R_common_badargs(); */
			sprintf(msg, _("WARNING: '--nsize' value is invalid: ignored\n"));
		    else
			sprintf(msg, _("WARNING: --nsize=%lu'%c': too large and ignored\n"),
			    value,
			    (ierr == 1) ? 'M': ((ierr == 2) ? 'K':'k'));
		    R_ShowMessage(msg);
		} else
		    Rp->nsize = value;
	    }
#endif
            else
            { /* unknown -option */
                argv[newac++] = *av;
            }
        }
        else
        {
            argv[newac++] = *av;
        }
    }
    *pac = newac;
    return;
}
