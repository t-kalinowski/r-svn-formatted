#include "Defn.h"
#include "Fileio.h"

#include "gnome-find-dialog.h"
#include "terminal.h"
#include "terminal-find.h"

#include <gnome.h>
#include <sys/types.h>
#include <regex.h>

#define ERRBUF_SIZE 500
#define MSGBUF_SIZE 2000

/* indexes for a find match */
typedef struct _find_selection find_selection;

struct _find_selection
{
    int select_start;
    int select_end;
};

static GnomeFindDialogParams find_params; /* find parameters from dialog box */

static gint find_pos, find_pos_max; /* current position of the find 'cursor', and the max (min is 0) */

static gchar *find_text_cache; /* cache of the console text */

static gboolean line_cache_update;            /* TRUE if line cache needs updating */
static gint line_cache_start, line_cache_end; /* indexes into find_text_cache of the area covered by the line cache */
static GList *find_current_match;             /* List of matches on currently cached line */

regex_t *preg; /* compiled regular expression structure */

/* FIXME: find/find again needs to check if regex button state has changed; force line cache update if it has */
/* FIXME: save/restore find params to dialog box */

void find_update_text_cache(void)
{
    /* update text cache */
    find_text_cache = gtk_editable_get_chars(GTK_EDITABLE(R_gtk_terminal_text), 0, -1);
    find_pos_max = gtk_text_get_length(GTK_TEXT(R_gtk_terminal_text)) - 1;

    /* force line cache update */
    line_cache_update = TRUE;
}

void find_compile_regex(GnomeFindDialog *find_dialog)
{
#ifdef HAVE_REGCOMP
    int cflags;
    int regex_result;
    char errbuf[ERRBUF_SIZE];
    char messagebuf[MSGBUF_SIZE];

    GtkWidget *message_dialog;

    /* FIXME: cache compiled regular expression */

    if (find_params.regular_exp == TRUE)
    {
        /* compile the regular expression */
        cflags = REG_EXTENDED;
        if (find_params.case_sensitive == FALSE)
        {
            cflags = cflags | REG_ICASE;
        }

        preg = (regex_t *)g_malloc(sizeof(regex_t));

        regex_result = regcomp(preg, find_params.find_text, cflags);

        if (regex_result != 0)
        {
            regerror(regex_result, preg, errbuf, ERRBUF_SIZE);
            g_snprintf(messagebuf, MSGBUF_SIZE, "Error compiling regular expression: %s", errbuf);

            message_dialog = gnome_message_box_new(messagebuf, GNOME_MESSAGE_BOX_ERROR, GNOME_STOCK_BUTTON_OK, NULL);
            if (find_dialog != NULL)
            {
                gnome_dialog_set_parent(GNOME_DIALOG(message_dialog), GTK_WINDOW(find_dialog));
            }
            else
            {
                gnome_dialog_set_parent(GNOME_DIALOG(message_dialog), GTK_WINDOW(R_gtk_main_window));
            }
            gnome_dialog_run_and_close(GNOME_DIALOG(message_dialog));

            return;
        }
    }
#endif /* HAVE_REGCOMP */
}

void find_free_regex(void)
{
#ifdef HAVE_REGCOMP
    /* FIXME: this will have to change if compiled regex is cached */
    if (preg != NULL)
    {
        regfree(preg);
        preg = NULL;
    }
#endif /* HAVE_REGCOMP */
}

void find_free_select(gpointer data, gpointer user_data)
{
    if (data != NULL)
    {
        g_free(data);
    }
}

int find_update_line_cache(GnomeFindDialog *find_dialog)
{
    int find_text_len;
    int eflags;
    regmatch_t pmatch[2];
    int regex_result;
    char errbuf[ERRBUF_SIZE];
    char messagebuf[MSGBUF_SIZE];
    int strcmp_result;
    int tmp_find_pos;
    gchar *find_line_cache;
    find_selection *find_select;

    GtkWidget *message_dialog;

    /* initialisation */
    find_line_cache = NULL;
    find_text_len = strlen(find_params.find_text);
    line_cache_update = FALSE;

    /* free matches list if necessary */
    if (find_current_match != NULL)
    {
        find_current_match = g_list_first(find_current_match);
        g_list_foreach(find_current_match, find_free_select, NULL);
        g_list_free(find_current_match);
        find_current_match = NULL;
    }

    /* Move one line */
    switch (find_params.direction)
    {
    case GNOME_FIND_FORWARDS:
        while ((find_pos <= find_pos_max) && (find_text_cache[find_pos] == '\n'))
        {
            find_pos++;
        }
        if (find_pos > find_pos_max)
        {
            return GNOME_FIND_NOMATCH;
        }
        break;

    case GNOME_FIND_BACKWARDS:
        while ((find_pos >= 0) && (find_text_cache[find_pos] == '\n'))
        {
            find_pos--;
        }
        if (find_pos < 0)
        {
            return GNOME_FIND_NOMATCH;
        }
        break;
    }

    line_cache_start = line_cache_end = find_pos;

    while ((line_cache_start > 0) && (find_text_cache[line_cache_start - 1] != '\n'))
    {
        line_cache_start--;
    }
    while ((line_cache_end < find_pos_max) && (find_text_cache[line_cache_end + 1] != '\n'))
    {
        line_cache_end++;
    }

    /* Grab the line */
    find_line_cache = g_strndup(find_text_cache + line_cache_start, line_cache_end - line_cache_start + 1);

    /* Search line */
#ifdef HAVE_REGCOMP
    if (find_params.regular_exp == TRUE)
    {
        /* Regular expression search */
        tmp_find_pos = line_cache_start;
        eflags = REG_EXTENDED;
        if (find_params.case_sensitive == FALSE)
        {
            eflags |= REG_ICASE;
        }
        do
        {
            if (tmp_find_pos != line_cache_start)
            {
                eflags |= REG_NOTBOL;
            }

            /* execute the match */
            regex_result = regexec(preg, find_line_cache, 1, pmatch, eflags);

            if (regex_result == 0)
            {
                /* construct list item if found */
                find_select = g_new(find_selection, 1);
                find_select->select_start = tmp_find_pos + pmatch[0].rm_so;
                find_select->select_end = tmp_find_pos + pmatch[0].rm_eo;
                find_current_match = g_list_append(find_current_match, (gpointer)find_select);

                tmp_find_pos += pmatch[0].rm_so + 1;
            }
            else if (regex_result != REG_NOMATCH)
            {
                /* report regexec errors and terminate the search */
                regerror(regex_result, preg, errbuf, ERRBUF_SIZE);
                g_snprintf(messagebuf, MSGBUF_SIZE, "Error matching regular expression: %s", errbuf);

                message_dialog =
                    gnome_message_box_new(messagebuf, GNOME_MESSAGE_BOX_ERROR, GNOME_STOCK_BUTTON_OK, NULL);
                if (find_dialog != NULL)
                {
                    gnome_dialog_set_parent(GNOME_DIALOG(message_dialog), GTK_WINDOW(find_dialog));
                }
                else
                {
                    gnome_dialog_set_parent(GNOME_DIALOG(message_dialog), GTK_WINDOW(R_gtk_main_window));
                }
                gnome_dialog_run_and_close(GNOME_DIALOG(message_dialog));

                if (find_line_cache != NULL)
                {
                    g_free(find_line_cache);
                }

                return GNOME_FIND_NOMATCH;
            }
        } while ((tmp_find_pos <= line_cache_end) && (regex_result == 0));
    }
    else
    {
#endif /* HAVE_REGCOMP */
        /* Literal search */
        for (tmp_find_pos = line_cache_start; tmp_find_pos <= line_cache_end; tmp_find_pos++)
        {
            if (find_params.case_sensitive == TRUE)
            {
                strcmp_result = strncmp(find_params.find_text, find_text_cache + tmp_find_pos, find_text_len);
            }
            else
            {
                strcmp_result = strncasecmp(find_params.find_text, find_text_cache + tmp_find_pos, find_text_len);
            }

            /* construct list item if found */
            if (strcmp_result == 0)
            {
                find_select = g_new(find_selection, 1);
                find_select->select_start = tmp_find_pos;
                find_select->select_end = tmp_find_pos + find_text_len;
                find_current_match = g_list_append(find_current_match, (gpointer)find_select);
            }
        }
#ifdef HAVE_REGCOMP
    }
#endif /* HAVE_REGCOMP */

    if (find_line_cache != NULL)
    {
        g_free(find_line_cache);
    }

    /* Return result */
    if (find_current_match != NULL)
    {
        return GNOME_FIND_MATCH;
    }

    return GNOME_FIND_NOTFOUND;
}

void find_process_result(GnomeFindDialog *find_dialog, int find_result)
{
    find_selection *find_select;

    GtkWidget *message_dialog;

    if (find_dialog != NULL)
        g_return_if_fail(GNOME_IS_FIND_DIALOG(find_dialog));

    switch (find_result)
    {
    case GNOME_FIND_NOMATCH:
        message_dialog = gnome_message_box_new("Could not find text in console output.", GNOME_MESSAGE_BOX_WARNING,
                                               GNOME_STOCK_BUTTON_OK, NULL);
        if (find_dialog != NULL)
        {
            gnome_dialog_set_parent(GNOME_DIALOG(message_dialog), GTK_WINDOW(find_dialog));
        }
        else
        {
            gnome_dialog_set_parent(GNOME_DIALOG(message_dialog), GTK_WINDOW(R_gtk_main_window));
        }
        gnome_dialog_run_and_close(GNOME_DIALOG(message_dialog));
        if (find_dialog != NULL)
            gnome_dialog_set_default(GNOME_DIALOG(find_dialog), GNOME_FIND_BUTTON_FIND);
        break;

    case GNOME_FIND_MATCH:
        find_select = (find_selection *)find_current_match->data;
        gtk_editable_set_position(GTK_EDITABLE(R_gtk_terminal_text), find_select->select_end);
        gtk_editable_select_region(GTK_EDITABLE(R_gtk_terminal_text), find_select->select_start,
                                   find_select->select_end);
        break;
    }
}

void R_gtk_terminal_find(GnomeFindDialog *find_dialog)
{
    int find_result, select_start, select_end;

    GtkWidget *message_dialog;

    if (find_dialog != NULL)
        g_return_if_fail(GNOME_IS_FIND_DIALOG(find_dialog));

    if (find_params.find_text == NULL)
        return;

    /* update the text cache if necessary */
    if (find_text_cache == NULL)
    {
        find_update_text_cache();
    }

    /* compile the regular expression */
    find_compile_regex(find_dialog);

    /* FIXME: wrapped searches */

    /* do the actual find */
    find_result = GNOME_FIND_NOTFOUND;

    while (find_result == GNOME_FIND_NOTFOUND)
    {
        /*    if ((line_cache_update == TRUE) || (find_current_match == NULL) || (find_pos < line_cache_start) ||
         * (find_pos > line_cache_end)) {*/
        if ((line_cache_update == TRUE) || (find_pos < line_cache_start) || (find_pos > line_cache_end))
        {
            find_result = find_update_line_cache(find_dialog);
        }
        else if (find_params.direction == GNOME_FIND_FORWARDS)
        {
            if ((find_current_match == NULL) || (find_current_match->next == NULL))
            {
                find_pos = line_cache_end + 1;
                find_result = find_update_line_cache(find_dialog);
            }
            else
            {
                find_current_match = find_current_match->next;
                find_result = GNOME_FIND_MATCH;
            }
        }
        else if (find_params.direction == GNOME_FIND_BACKWARDS)
        {
            if ((find_current_match == NULL) || (find_current_match->prev == NULL))
            {
                find_pos = line_cache_start - 1;
                find_result = find_update_line_cache(find_dialog);
            }
            else
            {
                find_current_match = find_current_match->prev;
                find_result = GNOME_FIND_MATCH;
            }
        }
    }

    /* update UI based on result */
    find_process_result(find_dialog, find_result);

    /* free the compiled regex */
    find_free_regex();
}

static void find_console_changed(GtkWidget *widget, gpointer data)
{
    if (find_text_cache != NULL)
    {
        g_free(find_text_cache);
        find_text_cache = NULL;
    }
}

static void find_dialog_cb(GtkWidget *widget, gpointer data)
{
    g_return_if_fail(widget != NULL);
    g_return_if_fail(GNOME_IS_FIND_DIALOG(widget));

    if (find_params.find_text != NULL)
        g_free(find_params.find_text);

    find_params = GNOME_FIND_DIALOG(widget)->params;
    find_params.find_text = gnome_find_dialog_get_find_text(GNOME_FIND_DIALOG(widget));

    switch (find_params.start_pos)
    {
    case GNOME_FIND_TOP:
        find_pos = 0;
        break;

    case GNOME_FIND_CURSOR:
        find_pos = gtk_editable_get_position(GTK_EDITABLE(R_gtk_terminal_text));
        break;

    case GNOME_FIND_BOTTOM:
        find_pos = gtk_text_get_length(GTK_TEXT(R_gtk_terminal_text)) - 1;
        break;
    }

    R_gtk_terminal_find(GNOME_FIND_DIALOG(widget));
}

static void find_again_dialog_cb(GtkWidget *widget, gpointer data)
{
    g_return_if_fail(widget != NULL);
    g_return_if_fail(GNOME_IS_FIND_DIALOG(widget));

    if (find_params.find_text != NULL)
        g_free(find_params.find_text);

    find_params = GNOME_FIND_DIALOG(widget)->params;
    find_params.find_text = gnome_find_dialog_get_find_text(GNOME_FIND_DIALOG(widget));

    R_gtk_terminal_find(GNOME_FIND_DIALOG(widget));
}

void edit_find_cb(GtkWidget *widget, gpointer data)
{
    GtkWidget *find_dialog;

#ifdef HAVE_REGCOMP
    find_dialog = gnome_find_dialog_new("Find text", NULL, TRUE, TRUE, TRUE);
#else
    find_dialog = gnome_find_dialog_new("Find text", NULL, TRUE, TRUE, FALSE);
#endif /* HAVE_REGCOMP */

    gnome_dialog_set_parent(GNOME_DIALOG(find_dialog), GTK_WINDOW(R_gtk_main_window));

    gtk_signal_connect(GTK_OBJECT(find_dialog), "find", (GtkSignalFunc)find_dialog_cb, NULL);
    gtk_signal_connect(GTK_OBJECT(find_dialog), "find_again", (GtkSignalFunc)find_again_dialog_cb, NULL);

    gtk_widget_show(find_dialog);
}

void edit_find_again_cb(GtkWidget *widget, gpointer data)
{
    R_gtk_terminal_find(NULL);
}

void R_gtk_terminal_find_init()
{
    /* Find functionality */
    find_params.find_text = NULL;
    find_pos = 0;
    find_text_cache = NULL;
    line_cache_update = TRUE;
    find_current_match = NULL;
    preg = NULL;

    gtk_signal_connect(GTK_OBJECT(R_gtk_terminal_text), "changed", (GtkSignalFunc)find_console_changed, NULL);
}
