/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1995 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* sre_string.c
 * 
 * my library of extra string functions. Some for portability
 * across UNIXes
 */

#include <stdio.h>
#include <string.h>
#include "squid.h"


#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


#ifdef NOSTR
char  *strstr(char *s, char *subs)
{
  int  i;

  for ( ; *s != 0; s++) {
    if (*s == *subs) {
      for (i = 1; subs[i] != 0 && subs[i] == s[i]; i++) ;
      if (subs[i] == 0) return(s);
      }
    }
  return (NULL);
}
#endif /* NOSTR */


#ifdef NO_STRDUP		/* many systems lack strdup() */
char *
strdup(char *s)
{
  char *new;
  if ((new = (char *) malloc (strlen(s) +1)) == NULL) return NULL;
  strcpy(new, s);
  return new;
}
#endif /* NO_STRDUP */

int
strinsert(char  *s1,            /* string to insert a char into  */
	  char   c,		/* char to insert                */
	  int    pos)		/* position in s1 to insert c at */
{
  char    oldc;
  char   *s;

  for (s = s1 + pos; c; s++)
    {
				/* swap current char for inserted one */
      oldc = *s;		/* pick up current */
      *s   = c;   		/* put down inserted one    */
      c    = oldc;		/* old becomes next to insert */
    }
  *s = '\0';

  return 1;
}


int
strdelete(char *s1,             /* string to delete a char from       */
	  int   pos)		/* position of char to delete 0..n-1  */
{
  char *s;                      

  for (s = s1 + pos; *s; s++)
    *s = *(s + 1);

  return 1;
}

void
s2lower(char *s)
{
  for (; *s != '\0'; s++)
    *s = sre_tolower((int) *s);
}

void
s2upper(char *s)
{
  for (; *s != '\0'; s++)
    *s = sre_toupper((int) *s);
}


void *
MallocOrDie(size_t size)
{
  void *ptr;
  if ((ptr = malloc (size)) == NULL)
    Die("malloc failed");
  return ptr;
}

void *
ReallocOrDie(void *p, size_t size)
{
  void *ptr;
  if ((ptr = realloc(p, size)) == NULL)
    Die("realloc failed");
  return ptr;
}
