/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1998, 1999  Guido Masarotto and Brian Ripley
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
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Declarations (for devX11.c) which are not only there,
 * ------------ but also, e.g., in  system.c  [event loop]
 */
int isX11DeviceActive(void);
void SaveX11DeviceAsGif(char *filename);
int X11DeviceDriver(DevDesc*, char*, double, double, double);

