/*
 *  R : A Computer Language for Statistical Data Analysis

 *  Copyright (C) 1996, 1997  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998-2001   Robert Gentleman, Ross Ihaka and the
 *                            R Development Core Team
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
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h> /* for NULL */
#include <limits.h>
#include <stdio.h>
#include <string.h>
#if defined(Macintosh)
#include <types.h>
#else
#include <sys/types.h>
#endif
/* #include <errno.h>*/
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "sock.h"

static int sock_inited = 0;

static int enter_sock(int fd)
{
#ifdef DEBUG
    printf("enter_sock(%d)\n", fd);
#endif
    if (fd == -1)
        return 0;
    else
        return fd;
}

static int close_sock(int fd)
{
    return Sock_close(fd, NULL) == -1 ? 0 : 1;
}

static void check_init(void)
{
    if (!sock_inited)
    {
#ifdef DEBUG
        printf("initing\n");
#endif
        Sock_init();
        sock_inited = 1;
    }
}

void Rsockopen(int *port)
{
    check_init();
    *port = enter_sock(Sock_open(*port, NULL));
}

void Rsocklisten(int *sockp, char **buf, int *len)
{
    check_init();
    *sockp = enter_sock(Sock_listen(*sockp, *buf, *len, NULL));
}

void Rsockconnect(int *port, char **host)
{
    check_init();
#ifdef DEBUG
    printf("connect to %d at %s\n", *port, *host);
#endif
    *port = enter_sock(Sock_connect(*port, *host, NULL));
}

void Rsockclose(int *sockp)
{
    *sockp = close_sock(*sockp);
}

void Rsockread(int *sockp, char **buf, int *maxlen)
{
    check_init();
#ifdef DEBUG
    printf("Reading from %d\n", *sockp);
#endif
    *maxlen = (int)Sock_read(*sockp, *buf, *maxlen, NULL);
}

void Rsockwrite(int *sockp, char **buf, int *start, int *end, int *len)
{
    ssize_t n;
    if (*end > *len)
        *end = *len;
    if (*start < 0)
        *start = 0;
    if (*end < *start)
    {
        *len = -1;
        return;
    }
    check_init();
#ifdef DEBUG
    printf("writing %s to %d", *buf, *sockp);
#endif
    n = Sock_write(*sockp, *buf + *start, *end - *start, NULL);
    *len = (int)n;
}

/* --------- unused ? ---------- */

#ifdef Unix
#include <signal.h>
#include <sys/wait.h>
static void sig_child(int sig)
{
    int stat;
    while (waitpid(-1, &stat, WNOHANG) > 0)
        ;
}

static int sig_fork_inited = 0;

void Rsockfork(int *pidno)
{
    pid_t pid;
    if (!sig_fork_inited)
    {
        struct sigaction sa;
        sa.sa_handler = sig_child;
        sa.sa_flags = 0;
        sigaction(SIGCHLD, &sa, NULL);
        sig_fork_inited = 1;
    }
    pid = fork();
    *pidno = (int)pid;
}
#endif

/* --------- for use in socket connections ---------- */

#if defined(Win32) || defined(HAVE_BSD_NETWORKING)
#include <R_ext/R-ftp-http.h>

#ifdef Win32
#include <winsock.h>
#include <io.h>
#define EWOULDBLOCK WSAEWOULDBLOCK
#define EINPROGRESS WSAEINPROGRESS
#endif

#ifdef HAVE_BSD_NETWORKING
#include <netdb.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#endif

#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif
#ifdef HAVE_ERRNO_H
#include <errno.h>
#endif
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef HAVE_SYS_SELECT_H
#include <sys/select.h>
#endif
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#ifndef Win32
#define closesocket(s) close(s)
#define SOCKET int
#endif

static int socket_errno(void)
{
#ifdef Win32
    return (WSAGetLastError());
#else
    return (errno);
#endif
}

#ifdef Unix
#include <R_ext/eventloop.h>

/* modified from src/unix/sys-std.c  */
static int setSelectMask(InputHandler *handlers, fd_set *readMask)
{
    int maxfd = -1;
    InputHandler *tmp = handlers;
    FD_ZERO(readMask);

    while (tmp)
    {
        if (tmp->fileDescriptor > 0)
        {
            FD_SET(tmp->fileDescriptor, readMask);
            maxfd = maxfd < tmp->fileDescriptor ? tmp->fileDescriptor : maxfd;
        }
        tmp = tmp->next;
    }

    return (maxfd);
}
#endif

static unsigned int timeout = 60;

static int R_SocketWait(int sockfd, int write)
{
    fd_set rfd, wfd;
    struct timeval tv;
    double used = 0.0;

    while (1)
    {
        int maxfd = 0, howmany;
#ifdef Unix
        InputHandler *what;

        if (R_wait_usec > 0)
        {
            R_PolledEvents();
            tv.tv_sec = 0;
            tv.tv_usec = R_wait_usec;
        }
        else
        {
            tv.tv_sec = timeout;
            tv.tv_usec = 0;
        }
#elif defined(Win32)
        tv.tv_sec = 0;
        tv.tv_usec = 2e5;
        R_ProcessEvents();
#else
        tv.tv_sec = timeout;
        tv.tv_usec = 0;
#endif

#ifdef Unix
        maxfd = setSelectMask(R_InputHandlers, &rfd);
#else
        FD_ZERO(&rfd);
#endif
        FD_ZERO(&wfd);
        if (write)
            FD_SET(sockfd, &wfd);
        else
            FD_SET(sockfd, &rfd);
        if (maxfd < sockfd)
            maxfd = sockfd;

        howmany = select(maxfd + 1, &rfd, &wfd, NULL, &tv);

        if (howmany < 0)
        {
            return -1;
        }
        if (howmany == 0)
        {
            used += tv.tv_sec + 1e-6 * tv.tv_usec;
            if (used >= timeout)
                return 1;
            continue;
        }

#ifdef Unix
        if (!FD_ISSET(sockfd, &rfd) || howmany > 1)
        {
            /* was one of the extras */
            what = getSelectedHandler(R_InputHandlers, &rfd);
            if (what != NULL)
                what->handler((void *)NULL);
            continue;
        }
#endif
        /* the socket was ready */
        break;
    }
    return 0;
}

void R_SockTimeout(int delay)
{
    timeout = (unsigned int)delay;
}

int R_SockConnect(int port, char *host)
{
    SOCKET s = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP);
    fd_set wfd, rfd;
    struct timeval tv;
    int status;
    double used = 0.0;
    struct sockaddr_in server;
    struct hostent *hp;

    check_init();

    if (s == -1)
        return -1;

#ifdef Win32
    {
        u_long one = 1;

        status = ioctlsocket(s, FIONBIO, &one) == SOCKET_ERROR ? -1 : 0;
    }
#else
    if ((status = fcntl(s, F_GETFL, 0)) != -1)
    {
#ifdef O_NONBLOCK
        status |= O_NONBLOCK;
#else /* O_NONBLOCK */
#ifdef F_NDELAY
        status |= F_NDELAY;
#endif /* F_NDELAY */
#endif /* !O_NONBLOCK */
        status = fcntl(s, F_SETFL, status);
    }
    if (status < 0)
    {
        closesocket(s);
        return (-1);
    }
#endif

    if (!(hp = gethostbyname(host)))
        return -1;

    memcpy((char *)&server.sin_addr, hp->h_addr_list[0], hp->h_length);
    server.sin_port = htons((short)port);
    server.sin_family = AF_INET;

    if ((connect(s, (struct sockaddr *)&server, sizeof(server)) == -1))
    {

        switch (socket_errno())
        {
        case EINPROGRESS:
        case EWOULDBLOCK:
            break;
        default:
            closesocket(s);
            return (-1);
        }
    }

    while (1)
    {
        int maxfd = 0;
#ifdef Unix
        InputHandler *what;

        if (R_wait_usec > 0)
        {
            R_PolledEvents();
            tv.tv_sec = 0;
            tv.tv_usec = R_wait_usec;
        }
        else
        {
            tv.tv_sec = timeout;
            tv.tv_usec = 0;
        }
#elif defined(Win32)
        tv.tv_sec = 0;
        tv.tv_usec = 2e5;
        R_ProcessEvents();
#else
        tv.tv_sec = timeout;
        tv.tv_usec = 0;
#endif

#ifdef Unix
        maxfd = setSelectMask(R_InputHandlers, &rfd);
#else
        FD_ZERO(&rfd);
#endif
        FD_ZERO(&wfd);
        FD_SET(s, &wfd);
        if (maxfd < s)
            maxfd = s;

        switch (select(maxfd + 1, &rfd, &wfd, NULL, &tv))
        {
        case 0:
            /* Time out */
            used += tv.tv_sec + 1e-6 * tv.tv_usec;
            if (used < timeout)
                continue;
            closesocket(s);
            return (-1);
        case -1:
            /* Ermm.. ?? */
            closesocket(s);
            return (-1);
        }

        if (FD_ISSET(s, &wfd))
        {
            SOCKLEN_T len;
            len = sizeof(status);
            if (getsockopt(s, SOL_SOCKET, SO_ERROR, (char *)&status, &len) < 0)
            {
                /* Solaris error code */
                return (-1);
            }
            if (status)
            {
                closesocket(s);
                errno = status;
                return (-1);
            }
            else
                return (s);
#ifdef Unix
        }
        else
        { /* some other handler needed */
            what = getSelectedHandler(R_InputHandlers, &rfd);
            if (what != NULL)
                what->handler((void *)NULL);
            continue;
#endif
        }
    }
    /* not reached */
    return (-1);
}

int R_SockClose(int sockp)
{
    return closesocket(sockp);
}

int R_SockRead(int sockp, void *buf, int len, int blocking)
{
    int res;

    if (blocking && R_SocketWait(sockp, 0) != 0)
        return 0;
    res = (int)recv(sockp, buf, len, 0);
    return (res >= 0) ? res : -socket_errno();
}

int R_SockOpen(int port)
{
    check_init();
    return Sock_open(port, NULL);
}

int R_SockListen(int sockp, char *buf, int len)
{
    check_init();
    return Sock_listen(sockp, buf, len, NULL);
}

int R_SockWrite(int sockp, const void *buf, int len)
{
    int res;

    /* Rprintf("socket %d writing |%s|\n", sockp, buf); */
    res = (int)send(sockp, buf, len, 0);
    return (res >= 0) ? res : -socket_errno();
}

#endif /* defined(Win32) || defined(HAVE_BSD_NETWORKING) */
