#include <extUtil.h>


/* summiert keystrokes, wenn key= "RETURN" dann 1, sonst 0 */
/*   record ist dann auszuwerten.                           */
int parser( char gkey, char *keystroke, int *curshft, int commandLineFlag)
{
  int i,j,n;
  char echo;

  if (gkey == ( char )0xff0d)  /* RETURN */
  {
    /* new-line, command ready */
    if(!commandLineFlag) putchar('\n');
    return(1);
  }
  else
  {
    i=strlen(keystroke)+*curshft;
    if (gkey == ( char )0xff08) /* backspace */
    {
      if(!commandLineFlag) {
      for(j=0; j<strlen(keystroke); j++)
      {
        /* go left */
        echo=( char )0xff08;
        putchar(echo);
      }
      for(j=0; j<strlen(keystroke); j++)
      {
        /* overwrite old command */
        echo=' ';
        putchar(echo);
      }
      for(j=0; j<strlen(keystroke); j++)
      {
        /* go left */
        echo=( char )0xff08;
        putchar(echo);
      }
      }

      n=strlen(keystroke);
      for (j=i-1; j<n; j++) keystroke[j]=keystroke[j+1];

      if(!commandLineFlag) {
      printf("%s",keystroke);

      for(j=0; j<-*curshft; j++)
      {
        /* go left */
        echo=( char )0xff08;
        putchar(echo);
      }
      }
    }
    else
    {
      if(!commandLineFlag) {
      for(j=0; j<strlen(keystroke); j++)
      {
        /* go left */
        echo=( char )0xff08;
        putchar(echo);
      }
      }

      n=strlen(keystroke);
      for (j=n; j>=i; j--) keystroke[j+1]=keystroke[j];
      keystroke[i] = gkey;

      if(!commandLineFlag) {
      printf("%s",keystroke);

      for(j=0; j<-*curshft; j++)
      {
        /* go left */
        echo=( char )0xff08;
        putchar(echo);
      }
      }
    }
    fflush(stdout);
  }
    DrawCommandLine(keystroke, strlen(keystroke)+*curshft);
  return(0);
}

