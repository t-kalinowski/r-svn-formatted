/*
 *  RProxy: Connector implementation between application and R language
 *  Copyright (C) 1999 Thomas Baier
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
 */

#include <windows.h>
#include <stdio.h>
#include "bdx.h"
#include "SC_proxy.h"
#include "rproxy_impl.h"
#include <assert.h>
#include <stdlib.h>

// magic number for data structure
#define R_BDX_MAGIC 0x8a35df01

typedef enum
{
    ps_none,
    ps_initialized
} R_Proxy_Object_State;

typedef struct _R_Proxy_Object_Impl
{
    SC_Proxy_Object_Vtbl *vtbl;
    R_Proxy_Object_State state;
    int ref_count;
} R_Proxy_Object_Impl;

int SYSCALL R_get_version(R_Proxy_Object_Impl *object, unsigned long *version)
{
    if ((object == NULL) || (version == NULL))
    {
        return SC_PROXY_ERR_INVALIDARG;
    }

    *version = SC_PROXY_INTERFACE_VERSION;

    return SC_PROXY_OK;
}

int SYSCALL R_init(R_Proxy_Object_Impl *object)
{
    int lRc = SC_PROXY_ERR_UNKNOWN;

    if (object == NULL)
    {
        return SC_PROXY_ERR_INVALIDARG;
    }

    if (object->state != ps_none)
    {
        return SC_PROXY_ERR_INITIALIZED;
    }

    lRc = R_Proxy_init();

    if (lRc == SC_PROXY_OK)
    {
        object->state = ps_initialized;
    }
    return lRc;
}

int SYSCALL R_terminate(R_Proxy_Object_Impl *object)
{
    int lRc = SC_PROXY_ERR_UNKNOWN;

    if (object == NULL)
    {
        return SC_PROXY_ERR_INVALIDARG;
    }

    if (object->state != ps_initialized)
    {
        return SC_PROXY_ERR_NOTINITIALIZED;
    }

    lRc = R_Proxy_term();

    if (lRc == SC_PROXY_OK)
    {
        object->state = ps_none;
    }

    return lRc;
}

int SYSCALL R_retain(R_Proxy_Object_Impl *object)
{
    if (object == NULL)
    {
        return SC_PROXY_ERR_INVALIDARG;
    }

    assert(object->ref_count > 0);

    (object->ref_count)++;

    return SC_PROXY_OK;
}

int SYSCALL R_release(R_Proxy_Object_Impl *object)
{
    if (object == NULL)
    {
        return SC_PROXY_ERR_INVALIDARG;
    }

    // reference count must not be 0 here
    assert(object->ref_count > 0);

    (object->ref_count)--;

    if (object->ref_count > 0)
    {
        return SC_PROXY_OK;
    }

    if (object->state != ps_none)
    {
        return SC_PROXY_ERR_INITIALIZED;
    }

    free(object);

    return SC_PROXY_OK;
}

int SYSCALL R_set_symbol(R_Proxy_Object_Impl *object, char const *symbol, BDX_Data *data)
{
    int lRc = 0;

    // check parameters
    if ((object == NULL) || (symbol == NULL) || (strlen(symbol) == 0) || (data == NULL))
    {
        return SC_PROXY_ERR_INVALIDARG;
    }

    if (data->magic != R_BDX_MAGIC)
    {
        return SC_PROXY_ERR_INVALIDFORMAT;
    }

    lRc = R_Proxy_set_symbol(symbol, data);

    return lRc;
}

int SYSCALL R_get_symbol(R_Proxy_Object_Impl *object, char const *symbol, BDX_Data **data)
{
    int lRc = 0;

    // check parameters
    if ((object == NULL) || (symbol == NULL) || (strlen(symbol) == 0) || (data == NULL))
    {
        return SC_PROXY_ERR_INVALIDARG;
    }

    lRc = R_Proxy_get_symbol(symbol, data);

    if (lRc == SC_PROXY_OK)
    {
        (*data)->magic = R_BDX_MAGIC;
    }

    return lRc;
}

int SYSCALL R_evaluate(R_Proxy_Object_Impl *object, char const *command, BDX_Data **data)
{
    if ((object == NULL) || (command == NULL) || (strlen(command) == 0) || (data == NULL))
    {
        return SC_PROXY_ERR_INVALIDARG;
    }

    if (object->state != ps_initialized)
    {
        return SC_PROXY_ERR_NOTINITIALIZED;
    }

    return R_Proxy_evaluate(command, data);
}

int SYSCALL R_evaluate_noreturn(R_Proxy_Object_Impl *object, char const *command)
{
    if ((object == NULL) || (command == NULL) || (strlen(command) == 0))
    {
        return SC_PROXY_ERR_INVALIDARG;
    }

    if (object->state != ps_initialized)
    {
        return SC_PROXY_ERR_NOTINITIALIZED;
    }

    return R_Proxy_evaluate_noreturn(command);
}

int SYSCALL R_query_types(R_Proxy_Object_Impl *object, long *type_mask)
{
    if ((object == NULL) || (type_mask == NULL))
    {
        return SC_PROXY_ERR_INVALIDARG;
    }

    *type_mask = 0;

    return SC_PROXY_OK;
}

int SYSCALL R_query_ops(R_Proxy_Object_Impl *object, long *op_mask)
{
    if ((object == NULL) || (op_mask == NULL))
    {
        return SC_PROXY_ERR_INVALIDARG;
    }

    *op_mask = 0;

    return SC_PROXY_OK;
}

int SYSCALL R_free_data_buffer(R_Proxy_Object_Impl *object, BDX_Data *data)
{
    if ((data == NULL) || (object == NULL))
    {
        return SC_PROXY_ERR_INVALIDARG;
    }

    if (data->magic != R_BDX_MAGIC)
    {
        return SC_PROXY_ERR_INVALIDFORMAT;
    }

    assert(data != NULL);
    assert(data->magic == R_BDX_MAGIC);

    bdx_free(data);
    //  free (data);

    return SC_PROXY_OK;
}

// global object table
SC_Proxy_Object_Vtbl global_proxy_object_vtbl = {(SC_PROXY_GET_VERSION)R_get_version,
                                                 (SC_PROXY_INIT)R_init,
                                                 (SC_PROXY_TERMINATE)R_terminate,
                                                 (SC_PROXY_RETAIN)R_retain,
                                                 (SC_PROXY_RELEASE)R_release,
                                                 (SC_PROXY_SET_SYMBOL)R_set_symbol,
                                                 (SC_PROXY_GET_SYMBOL)R_get_symbol,
                                                 (SC_PROXY_EVALUATE)R_evaluate,
                                                 (SC_PROXY_EVALUATE_NORETURN)R_evaluate_noreturn,
                                                 (SC_PROXY_QUERY_TYPES)R_query_types,
                                                 (SC_PROXY_QUERY_OPS)R_query_ops,
                                                 (SC_PROXY_FREE_DATA_BUFFER)R_free_data_buffer};

int SYSCALL EXPORT SC_Proxy_get_object(SC_Proxy_Object **obj, unsigned long version)
{
    R_Proxy_Object_Impl *proxy_object = NULL;

    if (obj == NULL)
    {
        return SC_PROXY_ERR_INVALIDARG;
    }

    if (version != SC_PROXY_INTERFACE_VERSION)
    {
        return SC_PROXY_ERR_INVALIDINTERPRETERVERSION;
    }

    proxy_object = (R_Proxy_Object_Impl *)malloc(sizeof(R_Proxy_Object_Impl));

    proxy_object->vtbl = &global_proxy_object_vtbl;
    proxy_object->state = ps_none;
    proxy_object->ref_count = 1;

    *obj = (SC_Proxy_Object *)proxy_object;

    return SC_PROXY_OK;
}
