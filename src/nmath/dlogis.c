/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
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

#include "Mathlib.h"

double dlogis(double x, double location, double scale)
{
    double e, f;
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(location) || ISNAN(scale))
        return x + location + scale;
#endif
    if (scale <= 0.0)
    {
        ML_ERROR(ME_DOMAIN);
        return ML_NAN;
    }
    e = exp(-(x - location) / scale);
    f = 1.0 + e;
    return e / (scale * f * f);
}
