/*
 * error.c -
 *
 * Error messages in FRiM.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of FRiM which is licensed under the MIT "Expat" License.
 *
 * Copyright (c) 2005-2020: Éric Thiébaut <https://github.com/emmt>
 *
 *-----------------------------------------------------------------------------
 */

#include "frim.h"

char const*frim_error_message(int const status)
{
    switch (status) {
#define _FRIM_ERROR(ident, mesg) case FRIM_PASTE(FRIM_ERROR_, ident): return mesg;
    _FRIM_ERROR_TABLE
#undef _FRIM_ERROR
    default: return "unknown error status";
    }
}
