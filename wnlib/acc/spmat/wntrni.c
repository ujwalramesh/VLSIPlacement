/**********************************************************************

wn_trans_problem_simplex_improve(&code,&objective,&delta,&new_result,
                                 result,cost_mat,
                                 max_time)

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

#include <stdio.h>

#include "wnlib.h"
#include "wnasrt.h"
#include "wnmax.h"
#include "wnsll.h"
#include "wnbll.h"
#include "wnhtab.h"
#include "wnhtbl.h"
#include "wnhash.h"
#include "wnnop.h"
#include "wncpy.h"
#include "wnrnd.h"

#include "wnspmat.h"
#include "wntrn.h"



typedef struct perif_data_struct *perif_data;
typedef struct entry_data_struct *entry_data;
typedef struct random_entry_queue_struct *random_entry_queue;

struct perif_data_struct
{
  wn_bll basic_list;
  wn_sll entry_list;

  bool dual_calculated;
  double dual;

  bool no_capacity_left;
  double capacity,remaining_capacity;

  wn_bll unsaturated_el;
};

struct entry_data_struct
{
  double flow,cost;

  bool is_basic;
  wn_bll i_el,j_el;

  int random_entry_queue_index;

  int i,j;
};

struct random_entry_queue_struct
{
  wn_htab entries_by_int;
  int size;
};


local random_entry_queue swap_basic_queue,simplex_queue;

local wn_htab entries_by_ij;
local wn_sll entry_list;

local int len_i,len_j,
          is_remaining,js_remaining,
	  time_remaining,
	  number_of_basics,
	  number_of_entries,acceptable_max_iterations,
	  number_of_improvements_found,
	  pivot_depth;
local wn_memgp top_group,free_group;
local perif_data *i_perif_array,*j_perif_array;
local double old_objective;
static double total = 0.0;
static int count = 0;
local bool transpose,matrix_very_sparse; 


local void entry_is_main_basic
(
  double flow,
  int i,
  int j
);
local void calculate_j_dual_implications(int j);
local void pivot_j
(
  bool *psuccess,
  entry_data entry,
  int target_i,
  double smallest_flow
);



local void insert_into_random_entry_queue
(
  random_entry_queue queue,
  entry_data entry
)
{
  bool success;

  wn_assert(entry->random_entry_queue_index == -1);

  success = wn_hins(entry,queue->entries_by_int,(ptr)queue->size);
  wn_assert(success);

  entry->random_entry_queue_index = queue->size;

  ++(queue->size);
}


#if 0 /* unused - bchapman 041111 */
local void randomize_perif_array_entry_list
(
  perif_data perif_array[],
  int len
)
{
  int i;

  for(i=0;i<len;++i)
  {
    wn_randomize_sll(&(perif_array[i]->entry_list));
  }
}
#endif /* 0 */


local void make_random_entry_queue(random_entry_queue *pqueue)
{
  wn_gppush(free_group);

  *pqueue = (random_entry_queue)
	 wn_zalloc(sizeof(struct random_entry_queue_struct));

  wn_mklonghtab(&((*pqueue)->entries_by_int));
  (*pqueue)->size = 0;

  wn_gppop();
}


local void delete_from_random_entry_queue
(
  random_entry_queue queue,
  entry_data entry
)
{
  int random_entry_queue_index;
  entry_data moved_entry;
  bool success;

  --(queue->size);
  wn_assert(queue->size >= 0);

  random_entry_queue_index = entry->random_entry_queue_index;
  entry->random_entry_queue_index = -1;

  wn_assert(random_entry_queue_index >= 0);

  success = wn_hdel(queue->entries_by_int,(ptr)random_entry_queue_index);
  wn_assert(success);

  if(random_entry_queue_index < queue->size)
  {
    success = wn_hget((ptr *)&moved_entry,queue->entries_by_int,
		      (ptr)queue->size);
    wn_assert(success);
    success = wn_hdel(queue->entries_by_int,(ptr)queue->size);
    wn_assert(success);
    success = wn_hins(moved_entry,queue->entries_by_int,
		      (ptr)random_entry_queue_index);
    wn_assert(success);

    moved_entry->random_entry_queue_index = random_entry_queue_index;
  }
}


local void get_random_entry
(
  entry_data *pentry,
  random_entry_queue queue
)
{
  int index,size;
  bool success;

  size = queue->size;

  if(size == 0)
  {
    *pentry = NULL;
  }
  else
  {
    index = wn_random_mod_int(size);

    success = wn_hget((ptr *)pentry,queue->entries_by_int,(ptr)index);
    wn_assert(success);
  }
}


struct index_struct
{
  int i,j;
};


local int hash_index_struct(struct index_struct *s)
{
  int hash;

  hash = 0;

  wn_cumhasho(&hash,wn_inthash(s->i));
  wn_cumhasho(&hash,wn_inthash(s->j));

  return(hash);
}


local bool index_structs_equal
(
  struct index_struct *s1,
  struct index_struct *s2
)
{
  return(
          (s1->i == s2->i)
            &&
          (s1->j == s2->j)
        );
}


local void copy_index_struct
(
  struct index_struct **ps,
  struct index_struct *s
)
{
  wn_memacpy((ptr *)ps,(ptr)s,sizeof(struct index_struct));
}


local void make_entries_by_ij(void)
{
  wn_mkhtab(&entries_by_ij,
            (int (*)(ptr))&hash_index_struct,
	    (bool (*)(ptr,ptr))&index_structs_equal,
            (void (*)(ptr *,ptr))&copy_index_struct,
	    (void (*)(ptr))&wn_do_nothing); 
}


local void insert_into_entries_by_ij(entry_data entry)
{
  static struct index_struct key_struct;
  struct index_struct *key = &key_struct;
  bool success;

  key->i = entry->i;
  key->j = entry->j;

  success = wn_hins(entry,entries_by_ij,key);
  wn_assert(success);
}


local bool get_from_entries_by_ij
(
  entry_data *pentry,
  int i,
  int j
)
{
  static struct index_struct key_struct;
  struct index_struct *key = &key_struct;

  key->i = i;
  key->j = j;

  return(wn_hget((ptr *)pentry,entries_by_ij,key));
}


local void install_entry(int i,int j,double cost)
{
  entry_data entry;

  entry = (entry_data)wn_zalloc(sizeof(struct entry_data_struct));

  entry->flow = 0.0;
  entry->cost = cost;

  entry->random_entry_queue_index = -1;

  entry->i = i;
  entry->j = j;

  wn_sllins(&(i_perif_array[i]->entry_list),entry);
  wn_sllins(&(j_perif_array[j]->entry_list),entry);

  insert_into_entries_by_ij(entry);
  wn_sllins(&entry_list,entry);
  insert_into_random_entry_queue(simplex_queue,entry);
  entry->random_entry_queue_index = -1;
}


local void compute_old_objective(void)
{
  int j;
  perif_data perif;
  register wn_sll el;
  register entry_data entry;
  register double flow,accum_old_objective;

  accum_old_objective = 0.0;

  for(j=0;j<len_j;++j)
  {
    perif = j_perif_array[j];

    for(el=perif->entry_list;el != NULL;el=el->next)
    {
      entry = (entry_data)(el->contents);

      flow = entry->flow;

      wn_assert(flow >= 0.0);

      if(flow > 0.0)
      {
        accum_old_objective += flow * entry->cost;
      }
    }
  }

  old_objective = accum_old_objective;
}


local void make_entries
(
  wn_sparse_matrix result,
  wn_sparse_matrix cost_mat
)
{
  int i,j;
  wn_sparse_matrix_entry matrix_entry;
  entry_data entry;
  wn_sll el;
  bool success;
  double flow;

  make_random_entry_queue(&swap_basic_queue);
  make_random_entry_queue(&simplex_queue);

  make_entries_by_ij();
  entry_list = NULL;

  for(j=0;j<len_j;++j)
  {
    for(el=(cost_mat->j_lists)[j];el!=NULL;el=el->next)
    {
      matrix_entry = (wn_sparse_matrix_entry)(el->contents);

      install_entry(matrix_entry->i,matrix_entry->j,matrix_entry->value);
    }
  }

  for(i=0;i<len_i;++i)
  {
    for(el=(result->i_lists)[i];el!=NULL;el=el->next)
    {
      matrix_entry = (wn_sparse_matrix_entry)(el->contents);

      flow = matrix_entry->value;

      wn_assert(flow >= 0.0);

      if(flow > 0.0)
      {
        success = get_from_entries_by_ij(&entry,
					 matrix_entry->i,matrix_entry->j);
        wn_assert(success);

        entry->flow = flow;
      }
    }
  }

  /*
  randomize_perif_array_entry_list(i_perif_array,len_i);
  randomize_perif_array_entry_list(j_perif_array,len_j);
  */
}


local void initialize
(
  wn_sparse_matrix result,
  wn_sparse_matrix cost_mat,
  int max_time
)
{
  perif_data perif;
  int i,j;

  top_group = wn_curgp();
  wn_gpmake("no_free");
  wn_gpmake("general_free");
  free_group = wn_curgp();
  wn_gppop();

  wn_assert(result->len_i == cost_mat->len_i);
  wn_assert(result->len_j == cost_mat->len_j);

  time_remaining = max_time;

  number_of_improvements_found = 0;

  /*
  transpose = (result->len_i < result->len_j);
  */
  transpose = FALSE;  /* no gain from transpose */
  if(transpose)
  {
    wn_transpose_sparse_matrix(result);
    wn_transpose_sparse_matrix(cost_mat);
  }

  len_i = cost_mat->len_i;
  len_j = cost_mat->len_j;

  i_perif_array = (perif_data *)wn_zalloc(len_i*sizeof(perif_data));
  j_perif_array = (perif_data *)wn_zalloc(len_j*sizeof(perif_data));

  for(i=0;i<len_i;i++)
  {
    perif = (perif_data)wn_zalloc(sizeof(struct perif_data_struct));

    perif->capacity = 0.0;
    
    i_perif_array[i] = perif;
  }
  for(j=0;j<len_j;j++)
  {
    perif = (perif_data)wn_zalloc(sizeof(struct perif_data_struct));
    
    perif->capacity = 0.0;

    j_perif_array[j] = perif;
  }

  make_entries(result,cost_mat);

  compute_old_objective();

  total = 0.0;
  count = 0;

  number_of_basics = len_i+len_j-1;
  number_of_entries = wn_sllcount(entry_list);
  acceptable_max_iterations = number_of_entries/15;
  matrix_very_sparse = (acceptable_max_iterations < number_of_basics);
}


local void decrement_time_remaining(int decrement)
{
  wn_assert(decrement >= 0);

  if(time_remaining < WN_IHUGE)
  {
    time_remaining -= decrement;
  }
}


local void compute_result
(
  double *pobjective,
  wn_sparse_matrix *presult
)
{
  double objective;
  entry_data entry;
  wn_sll el;

  wn_gppush(top_group);

  objective = 0.0;
  wn_make_sparse_matrix(presult,len_i,len_j);

  for(el=entry_list;el!=NULL;el=el->next)
  {
    entry = (entry_data)(el->contents);

    if(entry->is_basic && (entry->flow > 0.0))  
	 /* omit 0 flows, even if basic */
    {
       objective += (entry->flow * entry->cost);
       wn_insert_sparse_matrix_value(*presult,
                                     entry->flow,entry->i,entry->j);
    }
  }

  wn_sort_sparse_matrix(*presult);
  *pobjective = objective;

  wn_gppop();
}


local void load_remaining_capacities
(
  register perif_data perif_array[],
  register int len
)
{
  register int i;
  register perif_data perif;

  for(i=0;i<len;++i)
  {
    perif = perif_array[i]; 

    perif->remaining_capacity = perif->capacity;
    perif->no_capacity_left = FALSE;
  }
}


local void compute_capacities(wn_sparse_matrix result)
{
  register int i;
  register wn_sll el;
  register wn_sparse_matrix_entry matrix_entry;
  register double flow;

  for(i=0;i<len_i;++i)
  {
    for(el=(result->i_lists)[i];el != NULL;el=el->next)
    {
      matrix_entry = (wn_sparse_matrix_entry)(el->contents);

      flow = matrix_entry->value;

      if(flow > 0.0)
      {
	i_perif_array[i]->capacity += flow;
	j_perif_array[matrix_entry->j]->capacity += flow;
      }
    }
  }

  load_remaining_capacities(i_perif_array,len_i);
  load_remaining_capacities(j_perif_array,len_j);
}


local void compute_main_basics(wn_sparse_matrix result)
{
  int i;
  wn_sll el;
  wn_sparse_matrix_entry matrix_entry;
  double flow;

  for(i=0;i<len_i;++i)
  {
    for(el=(result->i_lists)[i];el!=NULL;el=el->next)
    {
      matrix_entry = (wn_sparse_matrix_entry)(el->contents);

      flow = matrix_entry->value;

      if(flow > 0.0)
      {
	entry_is_main_basic(flow,matrix_entry->i,matrix_entry->j);
      }
    }
  }
}


local void update_saturations_for_new_basic(int i,int j)
{
  perif_data i_perif,j_perif;

  i_perif = i_perif_array[i];
  j_perif = j_perif_array[j];

  wn_assert(!(i_perif->no_capacity_left)&&!(j_perif->no_capacity_left));

  if(i_perif->remaining_capacity == 0.0)
  {
    i_perif->no_capacity_left = TRUE;
    --is_remaining;
  }
  else if(j_perif->remaining_capacity == 0.0)
  {
    j_perif->no_capacity_left = TRUE;
    --js_remaining;
  }
  else
  {
    printf("i=%d,j=%d\n",i,j);
    printf("i_perif->remaining_capacity=%lg,j_perif->remaining_capacity=%lg\n",
           i_perif->remaining_capacity,j_perif->remaining_capacity);

    wn_assert_notreached();
  }
}


local void insert_into_i_basic_list(entry_data entry)
{
  wn_bll el;

  wn_gppush(free_group);

  wn_mkbll(&el);
  el->contents = (ptr)entry;
  wn_assert(entry->i_el == NULL);
  entry->i_el = el;

  wn_blleins(&(i_perif_array[entry->i]->basic_list),el);

  wn_gppop();
}


local void insert_into_j_basic_list(entry_data entry)
{
  wn_bll el;

  wn_gppush(free_group);

  wn_mkbll(&el);
  el->contents = (ptr)entry;
  wn_assert(entry->j_el == NULL);
  entry->j_el = el;

  wn_blleins(&(j_perif_array[entry->j]->basic_list),el);

  wn_gppop();
}


local void insert_into_basic_lists(entry_data entry)
{
  insert_into_i_basic_list(entry);
  insert_into_j_basic_list(entry);
}


local void note_new_basic(entry_data entry)
{
  if(entry->is_basic)
  {
    return;
  }

  wn_assert(!(entry->is_basic));
  entry->is_basic = TRUE;

  insert_into_random_entry_queue(swap_basic_queue,entry);

  insert_into_basic_lists(entry);
}


local void entry_is_main_basic
(
  double flow,
  int i,
  int j
)
{
  bool success;
  entry_data entry;

  i_perif_array[i]->remaining_capacity -= flow;
  j_perif_array[j]->remaining_capacity -= flow;

  update_saturations_for_new_basic(i,j);

  success = get_from_entries_by_ij(&entry,i,j);
  wn_assert(success);            /* can't have flow if no path */

  note_new_basic(entry);
}



local void compute_degenerate_basics(void)
{
  entry_data entry;
  int i,j;
  perif_data j_perif;
  bool success;

  wn_assert(is_remaining == 0);

  i = 0;

  for(j=0;j<len_j;++j)
  {
    j_perif = j_perif_array[j];

    if(!(j_perif->no_capacity_left))  /* unused capacity */
    {
      success = get_from_entries_by_ij(&entry,i,j);
      if(!(success))
      {
        install_entry(i,j,WN_FHUGE);  /* infinity cost */
        success = get_from_entries_by_ij(&entry,i,j);
	wn_assert(success);
      }

      if(!(entry->is_basic))
      {
        j_perif->no_capacity_left = TRUE;
        --js_remaining;

        note_new_basic(entry);
      }
    }
  }
}


local void compute_basic_entries(wn_sparse_matrix result)
{
  is_remaining = len_i;
  js_remaining = len_j;

  compute_capacities(result);
  compute_main_basics(result);

  if(is_remaining+js_remaining > 1)
  {
    compute_degenerate_basics();
  }

/*
  wn_assert(is_remaining+js_remaining == 1);
*/
}


local void delete_from_i_basic_list(entry_data entry)
{
  wn_bll el;

  wn_gppush(free_group);

  el = entry->i_el;
  wn_assert(el != NULL);
  entry->i_el = NULL;

  wn_blledel(el);

  wn_gppop();
}


local void delete_from_j_basic_list(entry_data entry)
{
  wn_bll el;

  wn_gppush(free_group);

  el = entry->j_el;
  wn_assert(el != NULL);
  entry->j_el = NULL;

  wn_blledel(el);

  wn_gppop();
}


local void delete_from_basic_lists(entry_data entry)
{
  delete_from_i_basic_list(entry);
  delete_from_j_basic_list(entry);
}


local void note_deleted_basic(entry_data entry)
{
  wn_assert(entry->is_basic);
  entry->is_basic = FALSE;

  if(entry->random_entry_queue_index != -1)
  {
    delete_from_random_entry_queue(swap_basic_queue,entry);
  }

  delete_from_basic_lists(entry);
}


local void keep_score
(
  char *type,
  double flow,
  double score
)
{
  total += flow*score;
  ++count;

  /*
  if(count%10 == 0)
  if(flow > 0.1)
  */
  {
    printf("%s -- %d: score=%20.20lg,flow=%lg,score=%20.20lg\n",
	   type,
	   count,old_objective + total,flow,score);
  }
}


#if 0	/* unused, commented out - bchapman 041111 */
local void swap_keep_score
(
  double flow,
  entry_data in1,
  entry_data in2,
  entry_data out1,
  entry_data out2
)
{
  double score;

  score = in1->cost + in2->cost - out1->cost - out2->cost;

  keep_score("swap",flow,score);
}


local void find_good_swap_involving_basic_entry
(
  bool *psuccess,
  entry_data *pnon_basic_in,
  entry_data *pbasic_in,
  entry_data *pout1,
  entry_data *pout2,
  entry_data basic_entry
)
{
  double partial_score1,partial_score2,partial_score3,score;
  perif_data i_perif,j_perif;
  wn_bll i_el,j_el;

  *psuccess = FALSE;

  i_perif = i_perif_array[basic_entry->i];
  j_perif = j_perif_array[basic_entry->j];
  *pbasic_in = basic_entry;
  partial_score1 = basic_entry->cost;

  for(i_el=i_perif->basic_list;i_el!=NULL;i_el=i_el->next)
  {
    *pout1 = (entry_data)(i_el->contents);

    if(*pout1 != basic_entry)
    {
      partial_score2 = partial_score1 - (*pout1)->cost; 

      for(j_el=j_perif->basic_list;j_el!=NULL;j_el=j_el->next)
      {
	*pout2 = (entry_data)(j_el->contents);

        if(*pout2 != basic_entry)
        {
          partial_score3 = partial_score2 - (*pout2)->cost; 

          decrement_time_remaining(1);

          if(get_from_entries_by_ij(pnon_basic_in,
				    (*pout2)->i,(*pout1)->j))
                     /* swap possible */
          {
            score = partial_score3 + (*pnon_basic_in)->cost;
  
            if(score < 0.0)  /* improvement */
            {
              *psuccess = TRUE;
  
              return;
	    }
          }
        } 
      }
    }
  }

  i_perif = i_perif_array[basic_entry->i];
  *pout1 = basic_entry;
  partial_score1 = -(*pout1)->cost;

  for(i_el=i_perif->basic_list;i_el!=NULL;i_el=i_el->next)
  {
    *pbasic_in = (entry_data)(i_el->contents);

    if(*pbasic_in != *pout1)
    {
      partial_score2 = partial_score1 + (*pbasic_in)->cost; 

      j_perif = j_perif_array[(*pbasic_in)->j];

      for(j_el=j_perif->basic_list;j_el!=NULL;j_el=j_el->next)
      {
	*pout2 = (entry_data)(j_el->contents);

        if(*pout2 != *pbasic_in)
        {
          partial_score3 = partial_score2 - (*pout2)->cost; 

          decrement_time_remaining(1);

          if(get_from_entries_by_ij(pnon_basic_in,
				    (*pout2)->i,(*pout1)->j))
                     /* swap possible */
          {
            score = partial_score3 + (*pnon_basic_in)->cost;
  
            if(score < 0.0)  /* improvement */
            {
              *psuccess = TRUE;
  
              return;
	    }
          }
        } 
      }
    }
  }

  j_perif = j_perif_array[basic_entry->j];
  *pout1 = basic_entry;
  partial_score1 = -(*pout1)->cost;

  for(j_el=j_perif->basic_list;j_el!=NULL;j_el=j_el->next)
  {
    *pbasic_in = (entry_data)(j_el->contents);

    if(*pbasic_in != *pout1)
    {
      partial_score2 = partial_score1 + (*pbasic_in)->cost; 

      i_perif = i_perif_array[(*pbasic_in)->i];

      for(i_el=i_perif->basic_list;i_el!=NULL;i_el=i_el->next)
      {
	*pout2 = (entry_data)(i_el->contents);

        if(*pout2 != *pbasic_in)
        {
          partial_score3 = partial_score2 - (*pout2)->cost; 

          decrement_time_remaining(1);

          if(get_from_entries_by_ij(pnon_basic_in,
				    (*pout1)->i,(*pout2)->j))
                     /* swap possible */
          {
            score = partial_score3 + (*pnon_basic_in)->cost;
  
            if(score < 0.0)  /* improvement */
            {
              *psuccess = TRUE;
  
              return;
	    }
          }
        } 
      }
    }
  }
}
#endif /* 0 */


#if 0	/* unused, commented out - bchapman 041111 */
local void perform_swap
(
  entry_data non_basic_in,
  entry_data basic_in,
  entry_data out1,
  entry_data out2
)
{
  double flow;

  flow = wn_min(out1->flow,out2->flow);

  wn_assert(flow >= 0.0);

  if(flow > 0.0)
  {
    ++number_of_improvements_found;

    non_basic_in->flow += flow;
    basic_in->flow += flow;

    out1->flow -= flow;
    out2->flow -= flow;
  }

  note_new_basic(non_basic_in);

  if(out1->flow == 0.0)
  {
    if(out2->flow == 0.0)
    {
      if(out1->cost > out2->cost)  /* toss most expensive */
      {
	note_deleted_basic(out1);
      }
      else
      {
	note_deleted_basic(out2);
      }
    }
    else /* out2->flow > 0.0 */
    {
      note_deleted_basic(out1);
    }
  }
  else /* out1->flow > 0.0 */
  {
    if(out2->flow == 0.0)
    {
      note_deleted_basic(out2);
    }
    else /* out2->flow > 0.0 */
    {
      wn_assert_notreached();
    }
  }

  swap_keep_score(flow,non_basic_in,basic_in,out1,out2);
  /*
  */
}


local void try_to_swap_basic_entry(entry_data basic_entry)
{
  bool success;
  entry_data non_basic_in,basic_in,out1,out2;

  find_good_swap_involving_basic_entry(&success,
				       &non_basic_in,&basic_in,&out1,&out2,
				       basic_entry);

  if(success)
  {
    perform_swap(non_basic_in,basic_in,out1,out2);
  }
  else
  {
    delete_from_random_entry_queue(swap_basic_queue,basic_entry);
  }
}
#endif /* 0 */


#if 0 /* unused - bchapman 041111 */
local void swap_improve_solution(void)
{
  entry_data basic_entry;

  while((time_remaining > 0) || (number_of_improvements_found < 4))
  {
    get_random_entry(&basic_entry,swap_basic_queue);

    if(basic_entry == NULL)  /* no more swaps possible */
    {
      break;
    }

    try_to_swap_basic_entry(basic_entry);
  }
}
#endif /* 0 */


local void dual_not_calculated
(
  register perif_data perif_array[],
  register int len
)
{
  register int i;

  for(i=0;i<len;++i)
  {
    perif_array[i]->dual_calculated = FALSE;
  }
}


local void calculate_i_dual_implications(int i)
{
  wn_bll el;
  entry_data entry;
  double i_dual;
  int j;
  perif_data i_perif,j_perif;

  i_perif = i_perif_array[i];
  i_dual = i_perif->dual;

  for(el=i_perif->basic_list;el!=NULL;el=el->next)
  {
    entry = (entry_data)(el->contents);

    decrement_time_remaining(1);

    j = entry->j;

    j_perif = j_perif_array[j];

    if(!(j_perif->dual_calculated))
    {
      j_perif->dual = entry->cost - i_dual;
      j_perif->dual_calculated = TRUE;
  
      calculate_j_dual_implications(j);
    }
  }
}


local void calculate_j_dual_implications(int j)
{
  wn_bll el;
  entry_data entry;
  double j_dual;
  int i;
  perif_data i_perif,j_perif;

  j_perif = j_perif_array[j];
  j_dual = j_perif->dual;

  for(el=j_perif->basic_list;el!=NULL;el=el->next)
  {
    entry = (entry_data)(el->contents);

    decrement_time_remaining(1);

    i = entry->i;

    i_perif = i_perif_array[i];

    if(!(i_perif->dual_calculated))
    {
      i_perif->dual = entry->cost - j_dual;
      i_perif->dual_calculated = TRUE;

      calculate_i_dual_implications(i);
    }
  }
}


local void calculate_duals(void)
{
  int i,j;
  perif_data perif;

  dual_not_calculated(i_perif_array,len_i);
  dual_not_calculated(j_perif_array,len_j);

  for(i=0;i<len_i;++i)
  {
    perif = i_perif_array[i];

    if(!(perif->dual_calculated))
    {
      perif->dual_calculated = TRUE;
      perif->dual = 0.0;

      calculate_i_dual_implications(i);
    }
  }
  for(j=0;j<len_j;++j)
  {
    perif = j_perif_array[j];

    if(!(perif->dual_calculated))
    {
      perif->dual_calculated = TRUE;
      perif->dual = 0.0;

      calculate_j_dual_implications(j);
    }
  }
}


local void calculate_good_random_pivot_entry(entry_data *pentry)
{
  int i;
  entry_data entry,best_entry;
  double cost,best_cost;

  decrement_time_remaining(number_of_basics);

  best_entry = NULL;
  best_cost = 0.0;

  for(i=0;i<number_of_basics;++i)
  {
    get_random_entry(&entry,simplex_queue);

    /* don't need to explicitly exclude basics because these
       will have cost of 0 */
    if(entry->is_basic)
    {
      continue;
    }

    cost = entry->cost - 
	   i_perif_array[entry->i]->dual - j_perif_array[entry->j]->dual;

    if(cost < best_cost)
    {
      best_cost = cost;
      best_entry = entry;
    }
  }

  *pentry = best_entry;
}


local void calculate_any_acceptable_random_pivot_entry(entry_data *pentry)
{
  int i;
  entry_data entry;
  double cost;

  for(i=0;i<acceptable_max_iterations;++i)
  {
    decrement_time_remaining(1);

    get_random_entry(&entry,simplex_queue);

    /* don't need to explicitly exclude basics because these
       will have cost of 0 */
    if(entry->is_basic)
    {
      continue;
    }

    cost = entry->cost - 
	   i_perif_array[entry->i]->dual - j_perif_array[entry->j]->dual;

    if(cost < 0.0)
    {
      *pentry = entry;

      return;
    }
  }

  *pentry = NULL;  /* nothing found */
}


/*
  if called, this routine is the inner loop.
*/
local void calculate_best_pivot_entry(entry_data *pentry)
{
  int i;
  perif_data i_perif;
  register entry_data entry,best_entry,next_best_entry;
  register double cost,best_cost,next_best_cost,i_dual;
  register wn_sll el;

  decrement_time_remaining(number_of_entries);

  next_best_entry = best_entry = NULL;
  next_best_cost = best_cost = 0.0;

  for(i=0;i<len_i;++i)
  {
    i_perif = i_perif_array[i];

    i_dual = i_perif->dual;

    for(el=i_perif->entry_list;el != NULL;el=el->next)
    {
      entry = (entry_data)(el->contents);

      /* don't need to explicitly exclude basics because these
	 will have cost of 0 */
      if(entry->is_basic)
      {
        continue;
      }

      cost = entry->cost - i_dual - j_perif_array[entry->j]->dual;

      if(cost < next_best_cost)
      {
	if(cost < best_cost)
	{
	  next_best_cost = best_cost;
	  next_best_entry = best_entry;

	  best_cost = cost;
	  best_entry = entry;
	}
	else
	{
	  next_best_cost = cost;
	  next_best_entry = entry;
	}
      }
    }
  }

  if(best_entry == NULL)
  {
    *pentry = NULL;
  }
  else if(next_best_entry == NULL)
  {
    *pentry = best_entry;
  }
  else
  {
    switch(wn_random_bit())  
    {
      case(0):
      {
	*pentry = best_entry;
      } break;
      case(1):
      {
	*pentry = next_best_entry;
      } break;
      default:
	wn_assert_notreached();
    }
  }
}


local void calculate_pivot_entry(entry_data *pentry)
{
  if(matrix_very_sparse)
  {
    calculate_best_pivot_entry(pentry);
  }
  else
  {
    calculate_good_random_pivot_entry(pentry);
    if(*pentry != NULL)
    {
      return;
    }

    calculate_any_acceptable_random_pivot_entry(pentry);
    if(*pentry != NULL)
    {
      return;
    }

    calculate_best_pivot_entry(pentry);  /* slow exhaustive search */
  }
}


local double flow_change,leaving_basic_cost;
local entry_data leaving_basic;

local void pivot_keep_score(entry_data entry)
{
  double cost;

  cost = entry->cost - 
	 i_perif_array[entry->i]->dual - 
	 j_perif_array[entry->j]->dual;

  if(flow_change > 0.1)
  {
    keep_score("simplex",flow_change,cost);
  /*
    printf("        pivot_depth=%d\n",pivot_depth);
  */
  }
}


local void pivot_i
(
  bool *psuccess,
  entry_data entry,
  int target_i,
  double smallest_flow
)
{
  entry_data i_entry;
  wn_bll basic_list,el;
  int i;

  *psuccess = FALSE;

  if(entry->flow < smallest_flow)
  {
    smallest_flow = entry->flow; 
  }

  i = entry->i;

  if(i == target_i)
  {
    flow_change = smallest_flow;

    *psuccess = TRUE;
  }
  else
  {
    basic_list = i_perif_array[i]->basic_list;

    for(el=basic_list;el!=NULL;el=el->next)
    {
      i_entry = (entry_data)(el->contents);

      decrement_time_remaining(1);

      if(i_entry != entry)
      {
        pivot_j(psuccess,i_entry,target_i,smallest_flow);
  
        if(*psuccess)
        {
          break;
        }
      }
    }
  }
  
  if(*psuccess)
  {
    /*
    ++pivot_depth;
    */

    entry->flow -= flow_change;
    wn_assert(entry->flow >= 0.0);

    if(entry->flow == 0.0)
    {
      if(entry->cost > leaving_basic_cost)
      {
	leaving_basic = entry;
	leaving_basic_cost = entry->cost;
      }
    }
  }
}


local void pivot_j
(
  bool *psuccess,
  entry_data entry,
  int target_i,
  double smallest_flow
)
{
  entry_data j_entry;
  wn_bll basic_list,el;

  *psuccess = FALSE;

  basic_list = j_perif_array[entry->j]->basic_list;

  for(el=basic_list;el!=NULL;el=el->next)
  {
    j_entry = (entry_data)(el->contents);

    decrement_time_remaining(1);

    if(j_entry != entry)
    {
      pivot_i(psuccess,j_entry,target_i,smallest_flow);

      if(*psuccess)
      {
        break;
      }
    }
  }

  if(*psuccess)
  {
    ++pivot_depth;

    entry->flow += flow_change;
  }
}


local void pivot_entry(entry_data entry)
{
  bool success;

  leaving_basic = NULL;
  leaving_basic_cost = -WN_FHUGE;

  pivot_depth = 0;

  pivot_j(&success,entry,entry->i,WN_FHUGE);
  wn_assert(success);

  if(flow_change > 0.0)
  {
    ++number_of_improvements_found;
  }

  pivot_keep_score(entry);
  /*
  */

  note_deleted_basic(leaving_basic);

  note_new_basic(entry);
}


local void simplex_improve_solution(int *pcode)
{
  entry_data entry;

  calculate_duals();
  calculate_pivot_entry(&entry);

  if(entry == NULL)  /* optimality!! */
  {
    *pcode = WN_SUCCESS;
  }
  else
  {
    pivot_entry(entry);

    *pcode = WN_SUBOPTIMAL;
  }
}


local void improve_solution(int *pcode)
{
  *pcode = WN_SUBOPTIMAL;

  while((time_remaining > 0) || (number_of_improvements_found < 4))
  {
/*
    swap_improve_solution();

    if(!((time_remaining > 0) || (number_of_improvements_found < 4)))
    {
      break;
    }
*/

    simplex_improve_solution(pcode);

    if(*pcode == WN_SUCCESS)
    {
      break;
    }
  }
}


local void finish(void)
{
  wn_gpfree();
}


void wn_trans_problem_simplex_improve
(
  int *pcode,
  double *pobjective,
  double *pdelta,
  wn_sparse_matrix *pnew_result,
  wn_sparse_matrix result,
  wn_sparse_matrix cost_mat,
  int max_time
)
{
  initialize(result,cost_mat,max_time);

  compute_basic_entries(result);

  improve_solution(pcode);

  compute_result(pobjective,pnew_result);
  *pdelta = *pobjective - old_objective;

  if(transpose)
  {
    wn_transpose_sparse_matrix(*pnew_result);
    wn_transpose_sparse_matrix(result);
    wn_transpose_sparse_matrix(cost_mat);
  }

  finish();
}



