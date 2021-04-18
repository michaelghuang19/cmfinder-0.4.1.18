#include "optimize.h"

#define GOLD    1.618034
#define GLIMIT  100
#define TINY    1e-20

#define CGOLD   0.3819660
#define ZEPS    1e-10
#define ITMAX   100

static int Minimize_step, mnbrak_step, brent_step;
static double guess_ax, guess_bx, brent_tol;
/*
static int IntMin_step, IntMin_ax, IntMin_bx;
*/
/*
   This function is called when an optimization starts. ax and bx
   should be guesses of values that probably has the minimum between
   them.
*/
void InitMinimize(double ax, double bx, double tol)
{
  Minimize_step = 0;
  mnbrak_step = 0;
  brent_step = 0;

  guess_ax = ax;
  guess_bx = bx;
  brent_tol = tol;
}

static int mnbrak(double *ax, double *bx, double *cx,
		  double *fa, double *fb, double *fc,
		  double *x, double *f);

static int brent(double *ax, double *bx, double *cx,
		 double *x_val, double *f);

/*
   This function is used when doing a function minimization. After
   having called InitMinimize, this function is called with any
   parameters.  Minimize sets *x to the next value to be
   evaluated. From then on, Minimize is called with the returned *x
   and the function value in that point. The return value of Minimize
   is as follows. 1: Minimize again. 0: ended succesfully, *x is the
   sought value. -1: an error occured.
*/
int Minimize(double *x, double *f)
{
  static double ax, bx, cx, fa, fb, fc;
  int value;

  if (Minimize_step == 0) {
    ax = guess_ax;
    bx = guess_bx;
    value = mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, x, f);
    Minimize_step = 1;
    return value;
  }
  if (Minimize_step == 1) {
    value = mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, x, f);
    if (value == -1) {
      Minimize_step = 2;
    }
    else
      return value;
  }
  if (Minimize_step == 2) {
    value = brent(&ax, &bx, &cx, x, f);
    return value;
  }

  /*  printf("%8.4f %8.4f\n%8.4f %8.4f\n%8.4f %8.4f\n", ax, fa, bx, fb, cx, fc);*/

  return 0;
}

/*
  The six first parameters as in Press et al. The last two as
  Minimize, return values as Minimize.
*/
static int mnbrak(double *ax, double *bx, double *cx,
		  double *fa, double *fb, double *fc,
		  double *x, double *f)
{
  static double dum, fu, q, r, u, ulim;

  /*    printf("step %d\n%8.4f %8.4f\n%8.4f %8.4f\n%8.4f %8.4f\n%8.4f %8.4f\n\n", mnbrak_step, *x, f, *ax, *fa, *bx, *fb, *cx, *fc);*/

  if (mnbrak_step == 0) {
    *x = *ax;
    mnbrak_step = 1;
    return 0;
  }
  if (mnbrak_step == 1) {
    *fa = *f;
    *x = *bx;
    mnbrak_step = 2;
    return 0;
  }
  if (mnbrak_step == 2) {
    *fb = *f;
    if (*fb > *fa) {
      dum = *ax;
      *ax = *bx;
      *bx = dum;
      dum = *fb;
      *fb = *fa;
      *fa = dum;
    }
    *cx = *bx + GOLD * (*bx - *ax);
    *x = *cx;
    mnbrak_step = 3;
    return 0;
  }
  if (mnbrak_step == 3) {
    *fc = *f;
label1:
    if (*fb >= *fc) {
      r = (*bx - *ax) * (*fb - *fc);
      q = (*bx - *cx) * (*fb - *fa);
      dum = q - r;
      if (dum < 0)
	dum = -dum;
      if (dum < TINY)
	  dum = TINY;
      if (q - r < 0)
	dum = -dum;
      u = *bx - ((*bx - *cx) * q - (*bx - *ax) * r) / (2 * dum);
      ulim = *bx + GLIMIT*(*cx-*bx);
      if ((*bx-u) * (u-*cx) > 0) {
	*x = u;
	mnbrak_step = 4;
	return 0;
      }
      else if ((*cx-u) * (u-ulim) > 0) {
	*x = u;
	mnbrak_step = 5;
	return 0;
      }
      else if ((u - ulim) * (ulim - *cx) >= 0) {
	u = ulim;
	mnbrak_step = 10;
      }
      else {
	u = *cx + GOLD * (*cx - *bx);
	mnbrak_step = 10;
      }
    }
    else
      return -1;
  }
  if (mnbrak_step == 4) {
    fu = *f;
    if (fu < *fc) {
      *ax = *bx;
      *fa = *fb;
      *bx = u;
      *fb = fu;
      return -1;
    }
    else if (fu > *fb) {
      *cx = u;
      *fc = fu;
      return -1;
    }
    u = *cx + GOLD * (*cx - *bx);
    mnbrak_step = 10;
  }
  if (mnbrak_step == 5) {
    fu = *f;
    if (fu < *fc) {
      *bx = *cx;
      *cx = u;
      u = *cx + GOLD * (*cx - *bx);
      *fb = *fc;
      *fc = fu;
      mnbrak_step = 10;
    }
    else
      mnbrak_step = 12;
  }
  if (mnbrak_step == 10) {
    *x = u;
    mnbrak_step = 11;
    return 0;
  }
  if (mnbrak_step == 11) {
    fu = *f;
    mnbrak_step = 12;
  }
  if (mnbrak_step == 12) {
    *ax = *bx;
    *bx = *cx;
    *cx = u;
    *fa = *fb;
    *fb = *fc;
    *fc = fu;
    goto label1;
  }

  return 1;
}


static int brent(double *ax, double *bx, double *cx,
		 double *x_val, double *f)
{
  static double a, b, d, e, etemp, fu, fv, fw, fx, p, q, r;
  static double tol1, tol2, u, v, w, x, xm;
  static int iter;

  /*  printf("step %d\n%8.4f %8.4f\n%8.4f %8.4f %8.4f\n\n", mnbrak_step, *x_val, *f, *ax, *bx, *cx);*/

  if (brent_step == 0) {
    a = *ax>*cx ? *cx : *ax;
    b = *ax>*cx ? *ax : *cx;
    v = *bx;
    w = v;
    x = v;
    e = 0.;
    *x_val = x;
    brent_step = 1;
    return 0;
  }
  if (brent_step == 1) {
    fx = *f;
    fv = fx;
    fw = fx;
    iter = 1;
  label_loop:
    if (iter > ITMAX)
      return 1;
    xm = 0.5*(a+b);
    tol1 = brent_tol*(x>0 ? x : -x) + ZEPS;
    tol2 = 2 * tol1;
    if ((x>xm ? x-xm : xm-x) <= tol2-.5*(b-a)) {
      *x_val = x;
      *f = fx;
      return -1;
    }
    if ((e>0 ? e : -e) > tol1) {
      r = (x-w) * (fx-fv);
      q = (x-v) * (fx-fw);
      p = (x-v) * q - (x-w) * r;
      q = 2*(q-r);
      if (q > 0)
	p = -p;
      q = q>0 ? q : -q;
      etemp = e;
      e = d;
      if ((p>0 ? p : -p) >= (q*etemp>0 ? .5*q*etemp : -.5*q*etemp) ||
	  p <= q*(a-x) || p >= q*(b-x))
	goto label1;
      d = p/q;
      u = x+d;
      if (u-a < tol2 || b-u < tol2)
	d = xm>x ? tol1 : -tol1;
      goto label2;
    }
  label1:
    if (x >= xm)
      e = a-x;
    else
      e = b-x;
    d = CGOLD * e;
  label2:
    if ((d>0 ? d : -d) >= tol1)
      u = x+d;
    else
      u = x + (d>0 ? tol1 : -tol1);
    *x_val = u;
    brent_step = 2;
    return 0;
  }
  if (brent_step == 2) {
    fu = *f;
    if (fu <= fx) {
      if (u >= x)
	a = x;
      else
	b = x;
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
    }
    else {
      if (u < x)
	a = u;
      else
	b = u;
      if (fu <= fw || w == x) {
	v = w;
	fv = fw;
	w = u;
	fw = fu;
      }
      else if (fu <= fv || v == x || v == w) {
	v = u;
	fv = fu;
      }
    }
    iter++;
    goto label_loop;
  }
  return 1;
}

/*
void InitIntMin(int ax, int bx)
{
  IntMin_step = 0;

  IntMin_ax = ax;
  IntMin_bx = bx;
}

int IntMin(int *x, int *f)
{
  static int ax, bx, cx;
  static double ay, by;

  if (IntMin_step == 0) {
    ax = IntMin_ax;
    bx = IntMin_bx;
    *x = ax;
    IntMin_step = 1;
  }
  else if (IntMin_step == 1) {
    ay = *f;
    *x = bx;
    IntMin_step = 2;
  }
  else if (IntMin_step == 2) {
    by = *f;
    cx = (ax+bx)/2;
    if (cx == ax || cx == bx) {
      if (ay < by) {
	*x = ax;
	*f = ay;
      }
      else {
	*x = bx;
	*f = by;
      }
      return 0;
    }
    *x = cx;
    IntMin_step = 3;
  }
  else if (IntMin_step == 3) {
    by = *f;
    cx = (ax+bx)/2;
    *x = cx;
    IntMin_step = 3;
  }

  return 1;
}
*/
