/*
 *  RProxy: Connector implementation between application and R language
 *  Copyright (C) 1999--2001 Thomas Baier
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
 *  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 *  MA 02111-1307, USA
 *
 *  $Id: rproxy_impl.h,v 1.6 2003/09/13 15:14:09 murdoch Exp $
 */

#ifndef _RPROXY_IMPL_H_
#define _RPROXY_IMPL_H_

#include "bdx.h"

/* exported functions for implementation */

/* 00-02-18 | baier | init() now receives parameter-string */
int R_Proxy_init (char const* pParameters);
int R_Proxy_evaluate (char const* pCmd,BDX_Data** pData);
int R_Proxy_evaluate_noreturn (char const* pCmd);
int R_Proxy_get_symbol (char const* pSymbol,BDX_Data** pData);
int R_Proxy_set_symbol (char const* pSymbol,BDX_Data const* pData);
int R_Proxy_term ();

#endif
