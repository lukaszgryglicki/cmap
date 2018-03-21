#include <stdio.h>
#include <stdlib.h>

void ado(FILE *O, int X, int Y)
{       fprintf(O,"%c!PS-Adobe-2.0 EPSF-2.0\n",'%');
        fprintf(O,"%c%cBoundingBox: 0 0 %d %d\n",'%','%',X,Y);
        fprintf(O,"/M {moveto} bind def\n");
        fprintf(O,"/L {lineto} bind def\n");
        fprintf(O,"/S {stroke} bind def\n");
        fprintf(O,"/s {show newpath} bind def\n");
        fprintf(O,"/C {closepath} bind def\n");
        fprintf(O,"/F {fill} bind def\n");
        fprintf(O,"/o {.1 0 360 arc C S} bind def\n");
        fprintf(O,"/times-Roman findfont 20 scalefont setfont\n");
        fprintf(O,"/W {setlinewidth} bind def\n");
        fprintf(O,"/RGB {setrgbcolor} bind def\n");
}

