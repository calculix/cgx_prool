/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */


#include "extUtil.h"


/* liest einen Record bis '\n'; uebergibt Anzahl gelesene Zeichen */
int frecord( FILE *handle1,  char *string)
{
  int i, c;

  for (i=0; i<MAX_LINE_LENGTH-1; i++)
  {
    string[i] = getc(handle1);
    if (string[i] == '\n')
      {
      string[i+1] = '\0';
      return(i);
      }
    if (string[i] == '\r')
      {
      c = getc(handle1);
      if (c != '\n')
        ungetc(c, handle1);

      string[i+1] = '\0';
      return(i);
     }
    else if (string[i] == (char)EOF)
      {
      string[i+1] = '\0';
      return(i);
      }
  }
  string[MAX_LINE_LENGTH-1] = '\0';
  return(MAX_LINE_LENGTH-1);
}


