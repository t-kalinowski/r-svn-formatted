/*******************************************************************************
 *  RProxy: Connector implementation between application and R language
 *  Copyright (C) 1999--2006 Thomas Baier
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Library General Public
 *  License as published by the Free Software Foundation; either
 *  version 2 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Library General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public
 *  License along with this library; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA 02110-1301, USA.
 *
 ******************************************************************************/

#ifndef _RPROXY_DEV_H_
#define _RPROXY_DEV_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "bdx.h"

#include <Defn.h>
#include <Graphics.h>
#include <Rdevices.h>


struct __tag_graphics_device {
#if 0
  unsigned long version;
#endif
  union {
    SC_GraphicsDevice*  d1;
    /* SC_GraphicsDevice2* __d2; */
  } device;
};

/* 06-08-20 | baier | data structures for graphics devices ("inherits" from
 *                    NewDevDesc) */
#define GFXDEV_ACTIVEX  1
#define GFXDEV_RECORDER 2
typedef struct __tag_graphics_device_description_cb
{
  NewDevDesc* dd;
} R_Proxy_Graphics_CB;
typedef struct __tag_graphics_device_description_recorder
{
  NewDevDesc* dd;
  void* recorder_data;                           /* for recording plots */
} R_Proxy_Graphics_Recorder;

extern struct __tag_graphics_device __graphics_device;

/* 06-05-17 | baier | graphics device */
/* 06-08-20 | baier | DEVDESC macro new */
#if 1
#define HASGFXDEV() (__graphics_device.device.d1 != NULL)

#define CLRGFXDEV() __graphics_device.device.d1 = NULL
#define GFXDEV() (__graphics_device.device.d1)
#define SETGFXDEV(d) __graphics_device.device.d1 = d;
#else
#define HASGFXDEV1() (__graphics_device.version == 1)
#define HASGFXDEV2() (__graphics_device.version == 2)

#define CLRGFXDEV() __graphics_device.version = 0; __graphics_device.device.d1 = NULL
#define GFXDEV1() (HASGFXDEV1() ? __graphics_device.device.d1 : NULL)
#define GFXDEV2() (HASGFXDEV2() ? __graphics_device.device.d2 : NULL)
#define SETGFXDEV1(d) __graphics_device.version = 1; __graphics_device.device.d1 = d;
#endif
#define DEVDESC(g) ((g)->dd)


#ifdef __cplusplus
}
#endif

#endif
