/****************************************************************************

COPYRIGHT NOTICE:

  The source code in this directory is provided free of
  charge to anyone who wants it.  It is in the public domain
  and therefore may be used by anybody for any purpose.  It
  is provided "AS IS" with no warranty of any kind
  whatsoever.  For further details see the README files in
  the wnlib parent directory.

****************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "wnlib.h"
#include "wnasrt.h"

#include "wncmp.h"
#include "wneq.h"


local void test_cmp()
{
  fprintf(stderr,"testing cmp...\n");

  wn_assert(wn_intcmp(1,1) == 0);
  wn_assert(wn_intcmp(1,0) > 0);
  wn_assert(wn_intcmp(0,1) < 0);
  wn_assert(wn_intcmp(-7,-7) == 0);
  wn_assert(wn_intcmp(-7,-4) < 0);
  wn_assert(wn_intcmp(-7,-401) > 0);
  wn_assert(wn_intcmp(-18,31) < 0);
  wn_assert(wn_intcmp(2000000000, -200000000) > 0);
  wn_assert(wn_intcmp(~ (int)(((unsigned)1) << 31), 
		      ((int)(((unsigned)1) << 31))) > 0);

  {
    char *p1="foo1",*p2="foo2";
    bool success;

    wn_assert(wn_ptrcmp(p1,p1) == 0);
    wn_assert(wn_ptrcmp(p2,p2) == 0);
    wn_assert(wn_ptrcmp(p1,p2) != 0);

    wn_assert(wn_ptrNULLcmp(&success,NULL,NULL) == 0);
    wn_assert(success);

    wn_assert(wn_ptrNULLcmp(&success,NULL,p1) < 0);
    wn_assert(success);

    wn_assert(wn_ptrNULLcmp(&success,p1,NULL) > 0);
    wn_assert(success);

    wn_assert(wn_ptrNULLcmp(&success,p1,p1) == 0);
    wn_assert(success);

    wn_assert(wn_ptrNULLcmp(&success,p1,p2) == 0);
    wn_assert(!success);
  }

  wn_assert(wn_boolcmp(TRUE,TRUE) == 0);
  wn_assert(wn_boolcmp(FALSE,TRUE) < 0);
  wn_assert(wn_boolcmp(TRUE,FALSE) > 0);
  wn_assert(wn_boolcmp(FALSE,FALSE) == 0);
  wn_assert(wn_boolcmp((bool)1,(bool)2) == 0);
  wn_assert(wn_boolcmp((bool)(-1),FALSE) > 0);

  wn_assert(wn_doublecmp(1.2,1.2) == 0);
  wn_assert(wn_doublecmp(1.3,1.2) > 0);
  wn_assert(wn_doublecmp(1.3,1.4) < 0);
  wn_assert(wn_doublecmp(-1.3,-1.4) > 0);
  wn_assert(wn_doublecmp(-1.3,0.4) < 0);

  {
    double *pd1,*pd2,
	   d1,d2;

    pd1 = &d1;
    pd2 = &d2;

    d1 = -1.3;
    d2 = -1.3;
    wn_assert(wn_pdoublecmp(pd1,pd2) == 0);

    d1 = 1.2;
    d2 = -1.3;
    wn_assert(wn_pdoublecmp(pd1,pd2) > 0);
  }

  wn_assert(wn_numstrcmp("wow1","wow1") == 0);
  wn_assert(wn_numstrcmp("wow2","wow1") > 0);
  wn_assert(wn_numstrcmp("wow1","wow2") < 0);
  wn_assert(wn_numstrcmp("wow10","wow2") > 0);
  wn_assert(wn_numstrcmp("wow2","wow10") < 0);

  fprintf(stderr,"  ok!!!!!!\n");
}


local void test_eq()
{
  fprintf(stderr,"testing eq...\n");

  wn_assert(wn_inteq(1,1));
  wn_assert(!wn_inteq(1,0));
  wn_assert(!wn_inteq(0,1));
  wn_assert(wn_inteq(-7,-7));
  wn_assert(!wn_inteq(-7,-4));
  wn_assert(!wn_inteq(-7,-401));
  wn_assert(!wn_inteq(-18,31));

  {
    char *p1="foo1",*p2="foo2";
    bool success;

    wn_assert(wn_ptreq(p1,p1));
    wn_assert(wn_ptreq(p2,p2));
    wn_assert(!wn_ptreq(p1,p2));

    wn_assert(wn_ptrNULLeq(&success,NULL,NULL));
    wn_assert(success);

    wn_assert(!wn_ptrNULLeq(&success,NULL,p1));
    wn_assert(success);

    wn_assert(!wn_ptrNULLeq(&success,p1,NULL) > 0);
    wn_assert(success);

    wn_assert(wn_ptrNULLeq(&success,p1,p1));
    wn_assert(success);

    wn_assert(wn_ptrNULLeq(&success,p1,p2));
    wn_assert(!success);
  }

  wn_assert(wn_streq("arf", "arf"));
  wn_assert(!wn_streq("arf0", "arf"));
  wn_assert(!wn_streq("arf", "arfa"));
  wn_assert(!wn_streq("arfa", "arfb"));
  wn_assert(!wn_streq("arfameow", "arfbmeow"));

  fprintf(stderr,"  ok!!!!!!\n");
}


int main(void)
{
  test_cmp();
  test_eq();

  return(0);
}
