/*
 * yor_frim.c -
 *
 * Yorick wrapper for FRiM.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of FRiM which is licensed under the MIT "Expat" License.
 *
 * Copyright (c) 2005-2020: Éric Thiébaut <https://github.com/emmt>
 *
 *-----------------------------------------------------------------------------
 */

#include <ydata.h>
#include "frim.h"

void Y_frim_error(int argc)
{
  if (argc != 1) YError("frim_error takes exactly one argument");
  YError(frim_error_message(YGetInteger(sp)));
}

#if 0
static void *get_array(Symbol *s, long *number, int *type);

void Y_op_pcg(int argc)
{
  long n, arg_number, number;
  int arg_type, type;
  Symbol *stack;
  void *p, *q, *r, *x, *z, *rho, *state;

  /* pcg(n,p,q,r,x,z,rho,state) */
  single = 0;
  if (argc != 7) {
    YError("op_pcg takes exactly 7 arguments");
  }
  s = sp - argc + 1;

  if (s->ops != &referenceSym) {
    YError("last argument must be a symbol");
  }
  s = &globTab[s->index];

  p = get_array(stack, &number, &type);
  if (! p || ! number || (type != T_DOUBLE && type != T_FLOAT)) {
  bad_arg:
    YError("P, Q, R, X and Z must be a real arrays with same number of elements");
  }
  q = get_array(stack + 1, &arg_number, &arg_type);
  if (! q || arg_number != number || arg_type != type) {
    goto bad_arg;
  }
  r = get_array(stack + 2, &arg_number, &arg_type);
  if (! r || arg_number != number || arg_type != type) {
    goto bad_arg;
  }
  x = get_array(stack + 3, &arg_number, &arg_type);
  if (! x || arg_number != number || arg_type != type) {
    goto bad_arg;
  }
  z = get_array(stack + 4, &arg_number, &arg_type);
  if (! z || arg_number != number || arg_type != type) {
    goto bad_arg;
  }
}

static void *get_array(Symbol *s, long *number, int *type)
{
  Operand op;
  if (! s->ops) {
    YError("unexpected keyword argument");
  }
  s->ops->FormOperand(s, &op);
  *type = op.ops->typeID;
  *number = op.type.number;
  return op.value;
}
#endif
