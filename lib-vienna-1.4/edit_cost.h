#ifndef EDIT_COST_H
#define EDIT_COST_H

/*   

cost.h   ::  global variables for Edit Costs
             included by treedist.c and stringdist.c

*/
#define PRIVATE static

#define  INF -10000  /* infinity */  

typedef double CostMatrix[21][21];

#ifndef CMFINDER /* skip unused variable warnings */
PRIVATE char   sep = ':';

PRIVATE char  *coding = "Null:A:C:G:U:AA:AC:AG:AU:CA:CC:CG:CU:GA:GC:GG:GU:UA:UC:UG:UU";

PRIVATE CostMatrix *EditCost;  /* will point to UsualCost or ShapiroCost */

PRIVATE CostMatrix  UsualCost= 
{
/*  Null,    A,    C,    G,    U,   AA,   AC,   AG,   AU,   CA,   CC,   CG,   CU,   GA,   GC,   GG,   GU,   UA,   UC,   UG,   UU, */
{   0.00,  -1.57,  -1.80,  -1.89,  -1.39,  -7.46,  -7.06,  -7.52,  -4.38,  -7.41,  -8.43,  -4.44,  -7.12,  -7.31,  -4.33,  -7.88,  -4.86,  -4.30,  -6.39,  -4.92,  -6.12},
{   -1.57,  2.22,  -1.86,  -1.46,  -1.39,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF},
{   -1.80,  -1.86,  1.16,  -2.48,  -1.05,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF},
{   -1.89,  -1.46,  -2.48,  1.03,  -1.74,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF},
{   -1.39,  -1.39,  -1.05,  -1.74,  1.65,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF},
{   -7.46,  INF,  INF,  INF,  INF,  -2.49,  -7.04,  -8.24,  -4.32,  -8.84,  -14.37,  -4.68,  -12.64,  -6.86,  -5.03,  -8.39,  -5.84,  -4.01,  -11.32,  -6.16,  -9.05},
{   -7.06,  INF,  INF,  INF,  INF,  -7.04,  -2.11,  -8.89,  -2.04,  -9.37,  -9.08,  -5.86,  -10.45,  -9.73,  -3.81,  -11.05,  -4.72,  -5.33,  -8.67,  -6.93,  -7.83},
{   -7.52,  INF,  INF,  INF,  INF,  -8.24,  -8.89,  -0.80,  -5.13,  -10.41,  -14.53,  -4.57,  -10.14,  -8.61,  -5.77,  -5.38,  -6.60,  -5.43,  -8.87,  -5.94,  -11.07},
{   -4.38,  INF,  INF,  INF,  INF,  -4.32,  -2.04,  -5.13,  4.49,  -5.56,  -6.71,  1.67,  -5.17,  -5.33,  2.70,  -5.61,  0.59,  1.61,  -4.81,  -0.51,  -2.98},
{   -7.41,  INF,  INF,  INF,  INF,  -8.84,  -9.37,  -10.41,  -5.56,  -5.13,  -10.45,  -3.57,  -8.49,  -7.98,  -5.95,  -11.36,  -7.93,  -2.42,  -7.08,  -5.63,  -8.39},
{   -8.43,  INF,  INF,  INF,  INF,  -14.37,  -9.08,  -14.53,  -6.71,  -10.45,  -3.59,  -5.71,  -5.77,  -12.43,  -3.70,  -12.58,  -7.88,  -6.88,  -7.40,  -8.41,  -5.41},
{   -4.44,  INF,  INF,  INF,  INF,  -4.68,  -5.86,  -4.57,  1.67,  -3.57,  -5.71,  5.36,  -4.96,  -6.00,  2.11,  -4.66,  -0.27,  2.75,  -4.91,  1.32,  -3.67},
{   -7.12,  INF,  INF,  INF,  INF,  -12.64,  -10.45,  -10.14,  -5.17,  -8.49,  -5.77,  -4.96,  -2.28,  -7.71,  -5.84,  -13.69,  -5.61,  -4.72,  -3.83,  -7.36,  -5.21},
{   -7.31,  INF,  INF,  INF,  INF,  -6.86,  -9.73,  -8.61,  -5.33,  -7.98,  -12.43,  -6.00,  -7.71,  -1.05,  -4.88,  -8.67,  -6.10,  -5.85,  -6.63,  -7.55,  -11.54},
{   -4.33,  INF,  INF,  INF,  INF,  -5.03,  -3.81,  -5.77,  2.70,  -5.95,  -3.70,  2.11,  -5.84,  -4.88,  5.62,  -4.13,  1.21,  1.60,  -4.49,  -0.08,  -3.90},
{   -7.88,  INF,  INF,  INF,  INF,  -8.39,  -11.05,  -5.38,  -5.61,  -11.36,  -12.58,  -4.66,  -13.69,  -8.67,  -4.13,  -1.98,  -5.77,  -5.75,  -12.01,  -4.27,  -10.79},
{   -4.86,  INF,  INF,  INF,  INF,  -5.84,  -4.72,  -6.60,  0.59,  -7.93,  -7.88,  -0.27,  -5.61,  -6.10,  1.21,  -5.77,  3.47,  -0.57,  -5.30,  -2.09,  -4.45},
{   -4.30,  INF,  INF,  INF,  INF,  -4.01,  -5.33,  -5.43,  1.61,  -2.42,  -6.88,  2.75,  -4.72,  -5.85,  1.60,  -5.75,  -0.57,  4.97,  -2.98,  1.14,  -3.39},
{   -6.39,  INF,  INF,  INF,  INF,  -11.32,  -8.67,  -8.87,  -4.81,  -7.08,  -7.40,  -4.91,  -3.83,  -6.63,  -4.49,  -12.01,  -5.30,  -2.98,  -3.21,  -4.76,  -5.97},
{   -4.92,  INF,  INF,  INF,  INF,  -6.16,  -6.93,  -5.94,  -0.51,  -5.63,  -8.41,  1.32,  -7.36,  -7.55,  -0.08,  -4.27,  -2.09,  1.14,  -4.76,  3.36,  -4.28},
{   -6.12,  INF,  INF,  INF,  INF,  -9.05,  -7.83,  -11.07,  -2.98,  -8.39,  -5.41,  -3.67,  -5.21,  -11.54,  -3.90,  -10.79,  -4.45,  -3.39,  -5.97,  -4.28,  -0.02}
};
#endif

void encode(int code, char* s);
int decode(char* c);
int PairCodeVienna(char c1, char c2, char* b);

#endif