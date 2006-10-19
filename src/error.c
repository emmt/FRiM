/*
 * error.c --
 *
 *	Error messages in FRIM.
 *
 *-----------------------------------------------------------------------------
 *
 *	Copyright (C) 2005-2006 Eric Thiébaut.
 *
 *	This file is part of FRIM (FRactal Iterative Method).
 *
 *	FRIM is free software; you can redistribute it and/or modify it
 *	under the terms of the GNU General Public License version 2 as
 *	published by the Free Software Foundation.
 *
 *	FRIM is distributed in the hope that it will be useful, but WITHOUT
 *	ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *	or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
 *	License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with FRIM (file "COPYING" in the top source directory); if
 *	not, write to the Free Software Foundation, Inc., 51 Franklin St,
 *	Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *-----------------------------------------------------------------------------
 *
 *	$Id$
 *	$Log$
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

