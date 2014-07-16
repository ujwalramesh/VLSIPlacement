/**********************************************************************

wn_copy_sparse_matrix(&out,in)
wn_symmetrize_sparse_matrix(mat)

**********************************************************************/
/****************************************************************************

COPYRIGHT NOTICE:

  The source code in this directory is provided free of
  charge to anyone who wants it.  It is in the public domain
  and therefore may be used by anybody for any purpose.  It
  is provided "AS IS" with no warranty of any kind
  whatsoever.  For further details see the README files in
  the wnlib parent directory.

AUTHOR:

  Will Naylor

****************************************************************************/

#include "wnlib.h"
#include "wnary.h"
#include "wnsmat.h"



wn_copy_sparse_matrix(pout,in)

wn_sparse_matrix *pout,in;

{
  int len_i,len_j,i;
  wn_sll el;
  wn_sparse_matrix_entry entry;

  len_i = in->len_i;
  len_j = in->len_j;

  wn_make_sparse_matrix(pout,len_i,len_j);

  for(i=0;i<len_i;++i)
  {
    for(wn_loopinit(),el=(in->i_lists)[i];wn_sllloop(&el,&entry);)
    {
      wn_insert_sparse_matrix_value(*pout,entry->value,i,entry->j);
    }
  }

  wn_sort_sparse_matrix(*pout);
}



/* 
  do this right later
*/
wn_symmetrize_sparse_matrix(mat)

wn_sparse_matrix mat;

{
  wn_sparse_matrix tmp_mat;
  int len_i,i,j;
  wn_sll el;
  wn_sparse_matrix_entry entry;

  wn_gpmake("no_free");

  wn_copy_sparse_matrix(&tmp_mat,mat);

  wn_gppush(mat->group);

  len_i = mat->len_i;

  for(i=0;i<len_i;++i)
  {
    for(wn_loopinit(),el=(tmp_mat->i_lists)[i];wn_sllloop(&el,&entry);)
    {
      j = entry->j;

      if(i != j)
      {
        wn_insert_sparse_matrix_value(mat,entry->value,j,i);
      }
    }
  }

  wn_sort_sparse_matrix(mat);

  wn_gppop();

  wn_gpfree();
}




