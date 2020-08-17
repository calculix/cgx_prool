/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */

#include <extUtil.h>

               
/* liest einen string bis ' ' oder '\n' */
/* uebergibt wordlaenge nit leading blanks */
int sword2( char *string, char *word)
{
  int i,n, flag;

    i=n=0;
    while (string[i] == ' ') i++;
    flag = sscanf(string, "%s", word);
    if ((flag <= 0)||(flag == (char)EOF )) { return(-1); }

    while (word[n] != '\0')
      n++;
    return(n+i);

}

/* liest einen string bis ' ' oder '\n' */
/* uebergibt wordlaenge ohne leading blanks*/
int sword( char *string, char *word)
{
  int n, flag;

    n=0;
    flag = sscanf(string, "%s", word);
    if ((flag <= 0)||(flag == (char)EOF )) { return(-1); }

    while (word[n] != '\0')
      n++;
    return(n);

}


