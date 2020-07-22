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

const char *frim_error_message(const int status)
{
  switch (status) {
#define _FRIM_ERROR(ident, message) \
  case FRIM_JOIN(FRIM_ERROR_, ident): return message;
#undef _FRIM_ERROR
  default: return "unknown error status";
  }
}
