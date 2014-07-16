/****************************************************************************

COPYRIGHT NOTICE:

  The source code in this directory is provided free of
  charge to anyone who wants it.  It is in the public domain
  and therefore may be used by anybody for any purpose.  It
  is provided "AS IS" with no warranty of any kind
  whatsoever.  For further details see the README files in
  the wnlib parent directory.

AUTHOR:

  Will Naylor, Bill Chapman

****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#define __NO_STRING_INLINES  /* long strcmp macro in asserts gets c89 upset */
#include <string.h>

#include "wnlib.h"
#include "wnasrt.h"
#include "wnmem.h"
#include "wnmemb.h"
#include "wnsll.h"
#include "wncmp.h"
#include "wnrnd.h"
#include "wnrndd.h"

#include "wnsort.h"
#include "wnrsrl.h"


WN_EXTERN_BEGIN
  typedef int  lo_cmp_func(ptr a, ptr b);
  local int lo_cmp_unsigned_top_3_bits(    ptr a_arg, ptr b_arg);
  local int lo_cmp_unsigned_bottom_29_bits(ptr a_arg, ptr b_arg);
  local int lo_cmp_unsigned(ptr arg_a, ptr arg_b);
WN_EXTERN_END

#define LENGTH(a)	((int) (sizeof(a)/sizeof((a)[0])))

/* for debugging */
local int print_string_list(wn_sll list)
{
  for (  ;  list;  list = list->next)
  {
    printf("%s", list->contents);
  }

  return 0;
} /* print_string_list */


local char string_index(char string[],int index)
{
  return(string[index]);
}


void wn_verify_string_list_is_sorted(wn_sll orig_list, wn_sll sorted_list)
{
  wn_sll slli, sllj;
  int orig_dup_count, sorted_dup_count;

  for (slli = sorted_list ;  slli && slli->next;  slli = slli->next)
  {
    wn_assert(strcmp((char *) slli->contents,
    /**/	       (char *) slli->next->contents) <= 0);
  }

  /* count duplicates */
  for (slli = orig_list ;  slli;  slli = slli->next)
  {
    orig_dup_count = 0;
    for (sllj = orig_list;  sllj;  sllj = sllj->next)
    {
      if (!strcmp((char *) slli->contents, (char *) sllj->contents))
      {
	++orig_dup_count;
      }
    }

    sorted_dup_count = 0;
    for (sllj = sorted_list;  sllj;  sllj = sllj->next)
    {
      if (!strcmp((char *) slli->contents, (char *) sllj->contents))
      {
	++sorted_dup_count;
      }
    }

    wn_assert(orig_dup_count == sorted_dup_count);
  }
} /* wn_verify_string_list_is_sorted */


local void test_string_list_sort(void)
{
  wn_sll list, list2, list3;

  wn_gpmake("no_free");
    wn_gplabel("test_string_list_sort_group");

    list = NULL;

    wn_sllins(&list, (char *) "foo");
    wn_sllins(&list, (char *) "bar");
    wn_sllins(&list, (char *) "tim");
    wn_sllins(&list, (char *) "mouse");
    wn_sllins(&list, (char *) "ant");
    wn_sllins(&list, (char *) "turkey");
    wn_sllins(&list, (char *) "let's");
    wn_sllins(&list, (char *) "make");
    wn_sllins(&list, (char *) "this");
    wn_sllins(&list, (char *) "a");
    wn_sllins(&list, (char *) "bigger");
    wn_sllins(&list, (char *) "test");
    wn_sllins(&list, (char *) "case");
    wn_sllins(&list, (char *) "with");
    wn_sllins(&list, (char *) "a");
    wn_sllins(&list, (char *) "duplicate");

    wn_sllcpy(&list2, list);
    wn_sllcpy(&list3, list);

    wn_sort_sll(&list2,(int (*)(ptr,ptr))(strcmp));
    wn_verify_string_list_is_sorted(list, list2);

    /* degenerate case: hand it a sorted list */
    wn_sort_sll(&list2,(int (*)(ptr,ptr))(strcmp));
    wn_verify_string_list_is_sorted(list, list2);

    /* degenerate case: hand it a reverse sorted list */
    wn_sllrev(&list2);
    wn_sort_sll(&list2,(int (*)(ptr,ptr))(strcmp));
    wn_verify_string_list_is_sorted(list, list2);

    wn_radix_sort_sll(&list3, (char (*)(ptr,int))(string_index),
    /**/			    (int (*)(ptr))(strlen));
    wn_verify_string_list_is_sorted(list, list3);

    /* degenerate case: hand it a sorted list */
    wn_radix_sort_sll(&list3, (char (*)(ptr,int))(string_index),
    /**/			    (int (*)(ptr))(strlen));
    wn_verify_string_list_is_sorted(list, list3);

    /* degenerate case: hand it a reverse sorted list */
    wn_sllrev(&list3);
    wn_radix_sort_sll(&list3, (char (*)(ptr,int))(string_index),
    /**/			    (int (*)(ptr))(strlen));
    wn_verify_string_list_is_sorted(list, list3);

  wn_gpfree();
} /* test_string_list_sort */


local char *dump_list_to_file = NULL;

void wn_verify_punsigned_list_is_sorted(wn_sll orig_list, wn_sll sorted_list)
{
  wn_sll slli, sllj;
  int orig_dup_count, sorted_dup_count;
  int i;	/* count iterations for debugging */

  if (dump_list_to_file)
  {
    FILE *fp;

    fp = fopen(dump_list_to_file, "w");

    for (slli = sorted_list;  slli;  slli = slli->next)
    {
      fprintf(fp, "%d\n", * (int *) slli->contents);
    }

    fclose(fp);
  }

  for (slli = sorted_list, i = 0;  slli && slli->next;  slli = slli->next, ++i)
  {
    wn_assert(* (unsigned *) slli->contents <=
    /**/				* (unsigned *) slli->next->contents);
  }

  /* count duplicates */
  for (slli = orig_list ;  slli;  slli = slli->next)
  {
    orig_dup_count = 0;
    for (sllj = orig_list;  sllj;  sllj = sllj->next)
    {
      if (* (unsigned *) slli->contents  ==  * (unsigned *) sllj->contents)
      {
	++orig_dup_count;
      }
    }

    sorted_dup_count = 0;
    for (sllj = sorted_list;  sllj;  sllj = sllj->next)
    {
      if (* (unsigned *) slli->contents  ==  * (unsigned *) sllj->contents)
      {
	++sorted_dup_count;
      }
    }

    wn_assert(orig_dup_count == sorted_dup_count);
  }
} /* wn_verify_punsigned_list_is_sorted */


#if 1
  /* this should work on all architectures */
  local char unsigned_index(int *pint,int index)
  {
    return wn_radix_sort_int_index_char(*pint, index);
  }
#else
  /*   this could also work, unless someone ports to a little-endian
  ** architecture that wasn't thought of here. */
  local char unsigned_index(int *pint,int index)
  {
#   if (defined (__arm__) && ! defined (__ARMEB__)) || defined (__i386__) || \
	defined (__i860__) || defined (__ns32000__) || defined (__vax__)
      /* This is for little-endian machines */

      return(((unsigned char *)pint)[(sizeof(int)-1)-index]);
#   else
      /* big-endian machine */

      return(((unsigned char *)pint)[index]);
#   endif
  }
#endif


local int unsigned_len(int *pint)
{
  return(sizeof(unsigned int));
}


local int punsignedcmp(unsigned *pi1,unsigned *pi2)
{
  return(*pi1 - *pi2);
}


local void test_radix_punsigned_list_sort(void)
{
  int i, r, *r_copy;
  wn_sll list, list2, list3;

  wn_gpmake("no_free");
    wn_gplabel("test_radix_punsigned_list_sort_group");

    list = NULL;

    for(i=0;i<1000;i++)	/* sort verification is O(n**2), this number should
    **			** not be TOO big. */
    {
      r = wn_random_int_between(0,3000);	/* choose a number low
      **		enough that there will be some duplicates */

      r_copy = (int *) wn_alloc(sizeof(int));
      *r_copy = r;
      wn_sllins(&list,(ptr)r_copy);
    }

    wn_sllcpy(&list2, list);
    wn_sllcpy(&list3, list);

    wn_sort_sll(&list2,(int (*)(ptr,ptr))(punsignedcmp));
    wn_verify_punsigned_list_is_sorted(list, list2);

    /* degenerate case: hand it a sorted list */
    wn_sort_sll(&list2,(int (*)(ptr,ptr))(punsignedcmp));
    wn_verify_punsigned_list_is_sorted(list, list2);

    /* degenerate case: hand it a reverse sorted list */
    wn_sllrev(&list2);
    wn_sort_sll(&list2,(int (*)(ptr,ptr))(punsignedcmp));
    wn_verify_punsigned_list_is_sorted(list, list2);


    wn_radix_sort_sll(&list3, (char (*)(ptr,int))(unsigned_index),
    /**/			    (int (*)(ptr))(unsigned_len));
    wn_verify_punsigned_list_is_sorted(list, list3);

    /* degenerate case: hand it a sorted list */
    wn_radix_sort_sll(&list3, (char (*)(ptr,int))(unsigned_index),
    /**/			    (int (*)(ptr))(unsigned_len));
    wn_verify_punsigned_list_is_sorted(list, list3);

    /* degenerate case: hand it a reverse sorted list */
    wn_sllrev(&list2);
    wn_radix_sort_sll(&list3, (char (*)(ptr,int))(unsigned_index),
    /**/			    (int (*)(ptr))(unsigned_len));
    wn_verify_punsigned_list_is_sorted(list, list3);

  wn_gpfree();
} /* test_radix_punsigned_list_sort */


void wn_verify_pint_list_is_sorted(wn_sll orig_list, wn_sll sorted_list)
{
  wn_sll slli, sllj;
  int orig_dup_count, sorted_dup_count;
  int i;	/* count iterations for debugging */

  if (dump_list_to_file)
  {
    FILE *fp;

    fp = fopen(dump_list_to_file, "w");

    for (slli = sorted_list;  slli;  slli = slli->next)
    {
      fprintf(fp, "%d\n", * (int *) slli->contents);
    }

    fclose(fp);
  }

  for (slli = sorted_list, i = 0;  slli && slli->next;  slli = slli->next, ++i)
  {
    wn_assert(* (int *) slli->contents <= * (int *) slli->next->contents);
  }

  /* count duplicates */
  for (slli = orig_list ;  slli;  slli = slli->next)
  {
    orig_dup_count = 0;
    for (sllj = orig_list;  sllj;  sllj = sllj->next)
    {
      if (* (int *) slli->contents  ==  * (int *) sllj->contents)
      {
	++orig_dup_count;
      }
    }

    sorted_dup_count = 0;
    for (sllj = sorted_list;  sllj;  sllj = sllj->next)
    {
      if (* (int *) slli->contents  ==  * (int *) sllj->contents)
      {
	++sorted_dup_count;
      }
    }

    wn_assert(orig_dup_count == sorted_dup_count);
  }
} /* wn_verify_pint_list_is_sorted */


local char int_index(int *pint,int index)
{
  return wn_radix_sort_int_index_char(*pint, index);
}


local int int_len(int *pint)
{
  return(sizeof(int));
}


local int pintcmp(int *pi1,int *pi2)
{
  return(wn_intcmp(*pi1,*pi2));
}


local void test_radix_pint_list_sort(void)
{
  int i,r,*r_copy;
  wn_sll list,list2,list3;

  wn_gpmake("no_free");
    wn_gplabel("test_radix_pint_list_sort_group");

    list = NULL;

    for(i=0;i<1000;i++)	/* sort verification is O(n**2), this number should
    **			** not be TOO big. */
    {
      r = wn_random_int_between(0,3000) - 1500;	/* choose a number low
      **	enough that there will be some duplicates, make sure there are
      **	some -ve numbers */

      r_copy = (int *) wn_alloc(sizeof(int));
      *r_copy = r;
      wn_sllins(&list,(ptr)r_copy);
    }

    wn_sllcpy(&list2, list);
    wn_sllcpy(&list3, list);

    wn_sort_sll(&list2,(int (*)(ptr,ptr))(pintcmp));
    wn_verify_pint_list_is_sorted(list, list2);

    /* degenerate case: hand it a sorted list */
    wn_sort_sll(&list2,(int (*)(ptr,ptr))(pintcmp));
    wn_verify_pint_list_is_sorted(list, list2);

    /* degenerate case: hand it a reverse sorted list */
    wn_sllrev(&list2);
    wn_sort_sll(&list2,(int (*)(ptr,ptr))(pintcmp));
    wn_verify_pint_list_is_sorted(list, list2);


    wn_radix_sort_sll(&list3, (char (*)(ptr,int))(int_index),
    /**/			    (int (*)(ptr))(int_len));
    wn_verify_pint_list_is_sorted(list, list3);

    /* degenerate case: hand it a sorted list */
    wn_radix_sort_sll(&list3, (char (*)(ptr,int))(int_index),
    /**/			    (int (*)(ptr))(int_len));
    wn_verify_pint_list_is_sorted(list, list3);

    /* degenerate case: hand it a reverse sorted list */
    wn_sllrev(&list2);
    wn_radix_sort_sll(&list3, (char (*)(ptr,int))(int_index),
    /**/			    (int (*)(ptr))(int_len));
    wn_verify_pint_list_is_sorted(list, list3);

  wn_gpfree();
} /* test_radix_pint_list_sort */


/****************************************************************************/

local bool lo_verify_list_is_sorted(
  wn_sll sorted_list,
  lo_cmp_func pcmp_func
)
{
  wn_sll slli;

  if (!sorted_list  ||  !sorted_list->next)
  {
    return TRUE;
  }

  for (slli = sorted_list ;  slli && slli->next;  slli = slli->next)
  {
    if ((*pcmp_func)(slli->contents, slli->next->contents) > 0)
    {
      return FALSE;
    }
  }

  return TRUE;
} /* lo_verify_list_is_sorted */


/* ---------------------------------------------------------------- */

local int lo_cmp_unsigned_top_3_bits(ptr a_arg, ptr b_arg) {
  long unsigned a = (long unsigned) a_arg & ((long unsigned) 0x7 << (29));
  long unsigned b = (long unsigned) b_arg & ((long unsigned) 0x7 << (29));

  return a < b ? -1 : (a == b ? 0 : 1);
} /* lo_cmp_unsigned_top_3_bits */


/* ---------------------------------------------------------------- */

local int lo_cmp_unsigned_bottom_29_bits(ptr a_arg, ptr b_arg) {
  long int a = (long int) a_arg & ((1 << 29) -1);
  long int b = (long int) b_arg & ((1 << 29) -1);

  return a < b ? -1 : (a == b ? 0 : 1);
} /* lo_cmp_unsigned_bottom_29_bits */


/* ---------------------------------------------------------------- */

local int lo_cmp_unsigned(ptr arg_a, ptr arg_b)
{
  long unsigned a = (long unsigned) arg_a;
  long unsigned b = (long unsigned) arg_b;

  return a < b ? -1 : (a == b ? 0 : 1);
} /* lo_cmp_unsigned */


/****************************************************************************/

void lo_test_merge_sort_order(void)
{
  wn_sll list_to_sort;
  int count;
  int test_index;
  int r, i, run;
  static int list_lengths[] = {	10000, 12345, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
  /**/				10001, 15555, 12000, 10101, 13579, 3*1024,
  /**/				1111, 1234, 4321, 9876, 1001, 1000 };

  for (run = 0;  run < 10;  ++run) {
    for (test_index = 0;  test_index < LENGTH(list_lengths);  ++test_index)
    {
      WN_GPBEGIN("no_free");
	wn_gplabel("lo_test_merge_sort_order");

	list_to_sort = NULL;

	/* note if you make the list 2^n long, it will all break down
	** into 2's and you will never exercise the 3's code */
	for (i = list_lengths[test_index];  i > 0;  --i) {
	  r = wn_random_int();

	  wn_sllins(&list_to_sort, (ptr) r);
	}

	count = wn_sllcount(list_to_sort);
	wn_assert(count == list_lengths[test_index]);

	wn_sort_sll(&list_to_sort, &lo_cmp_unsigned_bottom_29_bits);
	wn_assert(wn_sllcount(list_to_sort) == count);
	wn_assert( lo_verify_list_is_sorted(list_to_sort,
	/**/				&lo_cmp_unsigned_bottom_29_bits));

	wn_sort_sll(&list_to_sort, &lo_cmp_unsigned_top_3_bits);
	wn_assert(wn_sllcount(list_to_sort) == count);
	wn_assert( lo_verify_list_is_sorted(list_to_sort,
	/**/				&lo_cmp_unsigned_top_3_bits));

	wn_assert(lo_verify_list_is_sorted(list_to_sort, &lo_cmp_unsigned));
      WN_GPEND();
    } /* for test_index */
  } /* for run */
} /* lo_test_merge_sort_order */
void wn_verify_string_array_is_sorted(char **sorted_array,
/**/				char **orig_array, int array_length)
{
  int i, j;
  int orig_dup_count, sorted_dup_count;

  for (i = 0;  i+1 < array_length;  ++i)
  {
    wn_assert(strcmp(sorted_array[i], sorted_array[i+1]) <= 0);
  }

  /* count duplicates */
  for (i = 0;  i < array_length;  ++i)
  {
    orig_dup_count = 0;
    for (j = 0;  j < array_length;  ++j)
    {
      if (!strcmp(orig_array[i], orig_array[j]))
      {
	++orig_dup_count;
      }
    }

    sorted_dup_count = 0;
    for (j = 0;  j < array_length;  ++j)
    {
      if (!strcmp(orig_array[i], sorted_array[j]))
      {
	++sorted_dup_count;
      }
    }

    wn_assert(orig_dup_count > 0);
    wn_assert(orig_dup_count == sorted_dup_count);
  }
} /* wn_verify_string_array_is_sorted */


void test_string_quicksort()
{
  char *orig_array[] = {"foo","bar","tim","mouse","ant","turkey",
  /**/		"let's","make", "this", "a", "bigger", "test", "case",
  /**/		"with","a", "duplicate"};
  char *sorted_array[LENGTH(orig_array)];

  wn_memcpy((ptr) sorted_array, (ptr) orig_array, sizeof(orig_array));

  wn_sort_ptrarray((ptr *)sorted_array,LENGTH(sorted_array),
  /**/						(int (*)(ptr, ptr))(strcmp));

  wn_verify_string_array_is_sorted(sorted_array, orig_array,
  /**/							LENGTH(orig_array));
} /* test_string_quicksort */


local void lo_test_int_array_sorted(long int *sorted_array,
/**/					long int *orig_array, int length)
{
  int orig_dup_count, sorted_dup_count;
  int i, j;

  for (i = 0;  i < length-1;  ++i)
  {
    wn_assert(sorted_array[i] <= sorted_array[i+1]);
  }

  for (i = 0;  i < length;  ++i)
  {
    orig_dup_count = 0;
    for (j = 0;  j < length;  ++j)
    {
      if (orig_array[i] == orig_array[j])
      {
	++orig_dup_count;
      }
    }

    sorted_dup_count = 0;
    for (j = 0;  j < length;  ++j)
    {
      if (orig_array[i] == sorted_array[j])
      {
	++sorted_dup_count;
      }
    }

    wn_assert(orig_dup_count == sorted_dup_count);
    wn_assert(orig_dup_count >= 1);
    wn_assert(orig_dup_count <= length);
  } /* for i */
} /* lo_test_int_array_sorted */


local void lo_test_int_quicksort()     /* array sort */
{
  long int array[1000], orig_array[1000];
  long int tmp;
  int i, j;

  for (i = 0;  i < LENGTH(array);  ++i)
  {
    orig_array[i] =
    array[i] = wn_random_int_between(0,3000) - 1500;	/* choose a number low
    **		enough that there will be some duplicates, make sure there are
    **		some -ve numbers */
  }

  wn_sort_ptrarray((ptr *)array, LENGTH(array), (int(*)(ptr,ptr))(wn_longcmp));
  lo_test_int_array_sorted(array, orig_array, LENGTH(array));

  /* degenerate case: already sorted data */
  wn_sort_ptrarray((ptr *)array, LENGTH(array), (int(*)(ptr,ptr))(wn_longcmp));
  lo_test_int_array_sorted(array, orig_array, LENGTH(array));

  /* degenerate case: reverse sorted data */
  for (i = 0, j = LENGTH(array)-1;  i < j;  ++i, --j)
  {
    tmp = array[i];
    array[i] = array[j];
    array[j] = tmp;
  }

  /* ---------------------------------------------------------------- */
  /* run the whole test again, only with very large numbers */

  wn_sort_ptrarray((ptr *)array, LENGTH(array), (int(*)(ptr,ptr))(wn_longcmp));
  lo_test_int_array_sorted(array, orig_array, LENGTH(array));

  for (i = 0;  i < LENGTH(array);  ++i)
  {
    orig_array[i] =
    array[i] = wn_random_int();
  }

  wn_sort_ptrarray((ptr *)array, LENGTH(array), (int(*)(ptr,ptr))(wn_longcmp));
  lo_test_int_array_sorted(array, orig_array, LENGTH(array));

  /* degenerate case: already sorted data */
  wn_sort_ptrarray((ptr *)array, LENGTH(array), (int(*)(ptr,ptr))(wn_longcmp));
  lo_test_int_array_sorted(array, orig_array, LENGTH(array));

  /* degenerate case: reverse sorted data */
  for (i = 0, j = LENGTH(array)-1;  i < j;  ++i, --j)
  {
    tmp = array[i];
    array[i] = array[j];
    array[j] = tmp;
  }

  wn_sort_ptrarray((ptr *)array, LENGTH(array), (int(*)(ptr,ptr))(wn_longcmp));
  lo_test_int_array_sorted(array, orig_array, LENGTH(array));
} /* lo_test_int_quicksort */


int main(void)
{
  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  fprintf(stderr,"testing sort stuff...\n");

  wn_gpmake("no_free");
    wn_gplabel("main_group");

    test_string_list_sort();	/* merge & radix list sorts */

    test_radix_punsigned_list_sort();
    test_radix_pint_list_sort();

    test_string_quicksort();
    lo_test_int_quicksort();

    lo_test_merge_sort_order();

  wn_gpfree();

  fprintf(stderr,"  ok!!!!!!\n");  

  return(0);
} /* main */
