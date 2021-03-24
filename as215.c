/* as215.f -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mlegev_(doublereal *x, integer *n, doublereal *para, 
	doublereal *vcov, integer *monit, integer *ifault)
{
    /* Initialized data */

    static char acti1[8] = "  NEWTON";
    static doublereal zero = 0.;
    static doublereal half = .5;
    static doublereal one = 1.;
    static doublereal accu = 1e-5;
    static doublereal acca = 1e-5;
    static doublereal accg = 1e-5;
    static doublereal stepu = .5;
    static doublereal stepa = .25;
    static doublereal stepg = .2;
    static integer maxit = 30;
    static char acti2[8] = "  ST.ASC";
    static integer maxev = 50;
    static doublereal srf = .25;
    static integer maxsr = 30;
    static doublereal small = .001;
    static doublereal vlneg = -1e37;
    static char acti3[8] = "  RESETK";
    static char acti4[8] = "  SR.INF";
    static char acti5[8] = "  SR.LIK";
    static char acti6[8] = "  MAX.SR";
    static char acti7[8] = "  MAX.EV";
    static char acti8[8] = "  MAX.IT";
    static char acti9[8] = "  CONVGD";

    /* Format strings */
    static char fmt_6000[] = "(/\002 MAXIMUM-LIKELIHOOD ESTIMATION OF GENERA"
	    "LIZED EXTREME\002,1x,\002VALUE DISTRIBUTION\002//\002 ITER EVA"
	    "L\002,8x,\002XI\002,5x,\002ALPHA\002,9x,\002K  ACTION\002,6x,"
	    "\002LOG-L\002,7x,\002GNORM\002)";
    static char fmt_6010[] = "(1x,i4,i5,3f10.4,1x,a,f11.3,1pd12.2)";
    static char fmt_6020[] = "(\0020\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    double log(doublereal), exp(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal a, d__, e, f, g, h__;
    static integer i__;
    static doublereal p, q, r__, u, y, z__, da, ai, dg, he, gg, an, gi, hh, 
	    pa, qa, ra;
    static integer kk;
    static doublereal se, du, rg, sh, ye, pq, pu, qu, ru, sy, daa, dag, gai, 
	    dgg, dua, dug, she, shh, pqg, duu, sye;
    static integer nsr;
    static doublereal dela, aigi, delg, fold, shhe, delu, gipq, syhe, xmin, 
	    xmax, syye, temp1, temp2;
    static integer neval;
    static doublereal ratio;
    static integer niter;
    static doublereal gnorm;
    static integer itype;

    /* Fortran I/O blocks */
    static cilist io___28 = { 0, 6, 0, fmt_6000, 0 };
    static cilist io___37 = { 0, 6, 0, fmt_6010, 0 };
    static cilist io___38 = { 0, 6, 0, fmt_6010, 0 };
    static cilist io___61 = { 0, 6, 0, fmt_6010, 0 };
    static cilist io___65 = { 0, 6, 0, fmt_6010, 0 };
    static cilist io___66 = { 0, 6, 0, fmt_6010, 0 };
    static cilist io___95 = { 0, 6, 0, fmt_6010, 0 };
    static cilist io___97 = { 0, 6, 0, fmt_6010, 0 };
    static cilist io___98 = { 0, 6, 0, fmt_6010, 0 };
    static cilist io___99 = { 0, 6, 0, fmt_6010, 0 };
    static cilist io___100 = { 0, 6, 0, fmt_6010, 0 };
    static cilist io___101 = { 0, 6, 0, fmt_6010, 0 };
    static cilist io___102 = { 0, 6, 0, fmt_6020, 0 };


/*        ALGORITHM AS215   APPL. STATIST. (1985) VOL. 34, NO. 3 */
/*        Modifications in AS R76 (1989) have been incorporated. */

/*        Additional modifications by J. R. M. Hosking, August 1994: */
/*        Modify steepest-ascent step so that it is invariant to */
/*        rescaling the data.  Improves chance of convergence from */
/*        poor initial values. */

/*        MAXIMUM-LIKELIHOOD ESTIMATION OF GENERALIZED EXTREME-VALUE */
/*        DISTRIBUTION */

    /* Parameter adjustments */
    --x;
    --para;
    --vcov;

    /* Function Body */

/*        ADDU,ACCA,ACCG ARE ACCURACY CRITERIA FOR TESTING CONVERGENCE */
/*        STEPU,STEPA,STEPG ARE MAXIMUM STEPLENGTHS FOR ITERATIONS */
/*        ACCU,ACCA,STEPU,STEPA ARE SCALED BY CURRENT VALUE OF A WHEN */
/*        USED IN PROGRAM */


/*        MAXIT IS MAX. NO. OF ITERATIONS */
/*        MAXEV IS MAX. NO. OF EVALUATIONS OF LIKELIHOOD FUNCTION */
/*        SRF IS STEPLENGTH REDUCTION FACTOR */
/*        MAXSR IS MAX. NO. OF STEPLENGTH REDUCTIONS PERMITTED PER */
/*        ITERATION */


/*        SMALL IS A SMALL NUMBER, USED TO ADJUST THE SHAPE PARAMETER TO */
/*        AVOID AN EXACT ZERO VALUE OR BORDERLINE INFEASIBILITY */
/*        ALNEG IS A LARGE NEGATIVE NUMBER, USED TO INITIALIZE */
/*        LOG-LIKELIHOOD */


/*        FIND MIN AND MAX DATA VALUE */

    for (i__ = 1; i__ <= 6; ++i__) {
/* L10: */
	vcov[i__] = zero;
    }
    *ifault = 1;
    if (*n <= 2) {
	goto L170;
    }
    xmin = x[1];
    xmax = x[1];
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (x[i__] < xmin) {
	    xmin = x[i__];
	}
	if (x[i__] > xmax) {
	    xmax = x[i__];
	}
/* L20: */
    }

/*       INITIALIZATION */
/*       U IS LOCATION PARAMETER */
/*       A IS SCALE PARAMETER */
/*       G IS SHAPE PARAMETER */

    if (*monit > 0) {
	s_wsfe(&io___28);
	e_wsfe();
    }
    *ifault = 0;
    niter = 0;
    neval = 0;
    fold = vlneg;
    u = para[1];
    a = para[2];
    g = para[3];
    if (abs(g) < small) {
	g = small;
    }
    if (a <= zero) {
	a = one;
    }
    an = (doublereal) ((real) (*n));

/*        CHECK WHETHER ALL DATA POINTS LIE WITHIN THE RANGE OF THE GEV */
/*        DISTRIBUTION WITH THE INITIAL PARAMETERS - IF NOT, ADJUST THE */
/*        SHAPE PARAMETER SO AS TO BRING ALL POINTS WITHIN RANGE */

    if (g > zero) {
	goto L30;
    }
    if (xmin >= u) {
	goto L40;
    }
    z__ = a / (xmin - u);
    if (g > z__) {
	goto L40;
    }
    if (*monit > 0) {
	s_wsfe(&io___37);
	do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&neval, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&u, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&a, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&g, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, acti3, (ftnlen)8);
	e_wsfe();
    }
    g = z__ + small;
    if (g >= zero) {
	g = half * z__;
    }
    goto L40;
L30:
    if (xmax <= u) {
	goto L40;
    }
    z__ = a / (xmax - u);
    if (g < z__) {
	goto L40;
    }
    if (*monit > 0) {
	s_wsfe(&io___38);
	do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&neval, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&u, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&a, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&g, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, acti3, (ftnlen)8);
	e_wsfe();
    }
    g = z__ - small;
    if (g <= zero) {
	g = half * z__;
    }

/*        START OF MAIN LOOP */

L40:
    i__1 = maxit;
    for (niter = 1; niter <= i__1; ++niter) {

	nsr = 0;
L50:
	if (neval >= maxev) {
	    goto L150;
	}
	++neval;
	ai = one / a;
	gi = one / g;
	gai = g * ai;
	aigi = ai * gi;
	gg = one - g;

/*       ACCUMULATE SUMS OF QUANTITIES OCCURRING IN LIKELIHOOD */
/*       DERIVATIVES */

/*       IN PRESCOTT AND WALDEN'S NOTATION: */
/*       Z IS 1 - K * (X(I)-U) / A */
/*       Y IS THE REDUCED VARIATE - (1/K) * LOG(Z) */
/*       E IS EXP(-Y) */
/*       H IS EXP(K*Y) */

	sy = zero;
	se = zero;
	sye = zero;
	syye = zero;
	sh = zero;
	she = zero;
	syhe = zero;
	shhe = zero;
	shh = zero;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__ = one - gai * (x[i__] - u);
	    y = -gi * log(z__);
	    e = exp(-y);
	    h__ = one / z__;
	    ye = y * e;
	    he = h__ * e;
	    hh = h__ * h__;
	    sy += y;
	    se += e;
	    sye += ye;
	    syye += y * ye;
	    sh += h__;
	    she += he;
	    syhe += y * he;
	    shhe += hh * e;
	    shh += hh;
/* L60: */
	}

/*       F IS CURRENT VALUE OF LIKELIHOOD FUNCTIONN */

	f = -an * log(a) - gg * sy - se;
	if (f > fold) {
	    goto L90;
	}

/*       LIKELIHOOD HAS NOT INCREASED - REDUCE STEPLENGTH AND TRY AGAIN */

	if (*monit > 0) {
	    s_wsfe(&io___61);
	    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&neval, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&u, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&a, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&g, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, acti5, (ftnlen)8);
	    do_fio(&c__1, (char *)&f, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	if (nsr == maxsr) {
	    goto L80;
	}
L70:
	++nsr;
	u -= delu;
	a -= dela;
	g -= delg;
	delu = srf * delu;
	dela = srf * dela;
	delg = srf * delg;
	u += delu;
	a += dela;
	g += delg;
	if (a > g * (xmin - u) && a > g * (xmax - u) && g != zero) {
	    goto L50;
	}
	if (*monit > 0) {
	    s_wsfe(&io___65);
	    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&neval, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&u, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&a, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&g, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, acti4, (ftnlen)8);
	    e_wsfe();
	}
	if (nsr < maxsr) {
	    goto L70;
	}

/*        MAX. NO. OF STEPLENGTH REDUCTIONS REACHED */
/*        IF CURRENT ITERATION IS NEWTON-RAPHSON, TRY STEEPEST ASCENT */
/*        INSTEAD.  IF CURRENT ITERATION IS STEEPEST ASCENT, GIVE UP. */

L80:
	u -= delu;
	a -= dela;
	g -= delg;
	if (*monit > 0) {
	    s_wsfe(&io___66);
	    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&neval, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&u, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&a, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&g, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, acti6, (ftnlen)8);
	    e_wsfe();
	}
	if (itype == 1) {
	    goto L100;
	}
	*ifault = 4;
	goto L160;

/*        P,Q,R, ARE AS DEFINED IN FLOOD STUDIES REPORT */

L90:
	fold = f;
	p = an - se;
	q = she - gg * sh;
	r__ = an - sy + sye;
	pq = p + q;
	gipq = gi * pq;

/*        FIRST DERIVATIVES OF LOG-LIKELIHOOD */

	du = -ai * q;
	da = -aigi * pq;
	dg = -gi * (r__ - gipq);
	if (*monit > 0) {
	    gnorm = sqrt(du * du + da * da + dg * dg);
	}

/*        DERIVATIVES OF P,Q,R */

	pu = -ai * she;
	pa = gi * pu + aigi * se;
	qu = gg * ai * (shhe + g * shh);
	ru = ai * (sh - she + syhe);
	ra = gi * ru - aigi * (an - se + sye);
	rg = gi * (sy - sye + syye - a * ra);
	qa = ai * q + gi * (pu + qu);
	pqg = gipq + a * (ra - gi * (pa + qa));

/*         MINUS SECOND DERIVATIVE OF LOG-LIKELIHOOD (HESSIAN MATRIX) */

	duu = ai * qu;
	dua = aigi * (pu + qu);
	daa = -aigi * (ai * pq - pa - qa);
	dug = gi * (ru - gi * (pu + qu));
	dag = -aigi * (gipq - pqg);
	dgg = gi * (rg - gi * (pqg + r__ - gipq - gipq));

/*         INVERT HESSIAN MATRIX */

	for (kk = 1; kk <= 3; ++kk) {
	    if (duu <= zero) {
		goto L100;
	    }
	    d__ = one / duu;
	    temp1 = -dua * d__;
	    if (kk > 2) {
		temp1 = -temp1;
	    }
	    temp2 = -dug * d__;
	    if (kk > 1) {
		temp2 = -temp2;
	    }
	    duu = daa + temp1 * dua;
	    dua = dag + temp1 * dug;
	    daa = dgg + temp2 * dug;
	    dug = temp1;
	    dag = temp2;
	    dgg = d__;
/* L95: */
	}

/*        CALCULATE STEPLENGTHS */

	itype = 1;
	if (*monit > 0) {
	    s_wsfe(&io___95);
	    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&neval, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&u, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&a, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&g, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, acti1, (ftnlen)8);
	    do_fio(&c__1, (char *)&f, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	delu = duu * du + dua * da + dug * dg;
	dela = dua * du + daa * da + dag * dg;
	delg = dug * du + dag * da + dgg * dg;
/* Computing MAX */
	d__1 = abs(delu) / (stepu * a), d__2 = abs(dela) / (stepa * a), d__1 =
		 max(d__1,d__2), d__2 = abs(delg) / stepg;
	ratio = max(d__1,d__2);
	if (ratio < one) {
	    goto L110;
	}
	ratio = one / ratio;
	delu *= ratio;
	dela *= ratio;
	delg *= ratio;
	goto L110;

/*        HESSIAN IS NOT POSITIVE DEFINITE - MAKE A LARGE STEP IN THE */
/*        DIRECTION OF STEEPEST ASCENT */

L100:
	itype = 2;
	if (*monit > 0) {
	    s_wsfe(&io___97);
	    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&neval, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&u, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&a, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&g, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, acti2, (ftnlen)8);
	    do_fio(&c__1, (char *)&f, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
/* ---------------------------------------  NEXT 11 LINES ADDED, AUG.94 */
	d__ = abs(vlneg);
	temp1 = d__;
	if (du != zero) {
	    temp1 = stepu / (abs(du) * a);
	}
	temp2 = d__;
	if (da != zero) {
	    temp2 = stepa / (abs(da) * a);
	}
	z__ = d__;
	if (dg != zero) {
	    z__ = stepg / abs(dg);
	}
/* Computing MIN */
	d__1 = min(temp1,temp2);
	ratio = min(d__1,z__);
	delu = ratio * du * a * a;
	dela = ratio * da * a * a;
	delg = ratio * dg;
/* ---------------------------------------  NEXT 11 LINES REMOVED, AUG.94 */
/*      D = ABS(VLNEG) */
/*      TEMP1 = D */
/*      IF (DU .NE. ZERO) TEMP1 = STEPU * A / ABS(DU) */
/*      TEMP2 = D */
/*      IF (DA .NE. ZERO) TEMP2 = STEPA * A / ABS(DA) */
/*      Z = D */
/*      IF (DG .NE. ZERO) Z = STEPG / ABS(DG) */
/*      RATIO = MIN(TEMP1, TEMP2, Z) */
/*      DELU = RATIO * DU */
/*      DELA = RATIO * DA */
/*      DELG = RATIO * DG */
/* ---------------------------------------  END OF DELETED CODE */

/*        ADJUST PARAMETERS */

L110:
	u += delu;
	a += dela;
	g += delg;

/*        TEST FOR FEASIBILITY */

	if (a > g * (xmin - u) && a > g * (xmax - u) && g != zero) {
	    goto L130;
	}
	i__2 = maxsr;
	for (nsr = 1; nsr <= i__2; ++nsr) {
	    if (*monit > 0) {
		s_wsfe(&io___98);
		do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&neval, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&u, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&g, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, acti4, (ftnlen)8);
		e_wsfe();
	    }
	    u -= delu;
	    a -= dela;
	    g -= delg;
	    delu = srf * delu;
	    dela = srf * dela;
	    delg = srf * delg;
	    u += delu;
	    a += dela;
	    g += delg;
	    if (a > g * (xmin - u) && a > g * (xmax - y) && g != zero) {
		goto L140;
	    }
/* L120: */
	}
	goto L80;

/*        TEST FOR CONVERGENCE */

L130:
	if (abs(delu) > accu * a) {
	    goto L140;
	}
	if (abs(dela) > acca * a) {
	    goto L140;
	}
	if (abs(delg) > accg) {
	    goto L140;
	}
	if (*monit > 0) {
	    s_wsfe(&io___99);
	    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&neval, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&u, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&a, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&g, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, acti9, (ftnlen)8);
	    e_wsfe();
	}
	vcov[1] = duu;
	vcov[2] = dua;
	vcov[3] = daa;
	vcov[4] = dug;
	vcov[5] = dag;
	vcov[6] = dgg;
	goto L160;

/*        END OF MAIN LOOP */

L140:
	;
    }

/*        ITERATIONS NOT CONVERGED - SET FAULT FLAG */

    *ifault = 2;
    if (*monit > 0) {
	s_wsfe(&io___100);
	do_fio(&c__1, (char *)&maxit, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&neval, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&u, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&a, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&g, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, acti8, (ftnlen)8);
	e_wsfe();
    }
    goto L160;
L150:
    *ifault = 3;
    if (*monit > 0) {
	s_wsfe(&io___101);
	do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&maxev, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&u, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&a, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&g, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, acti7, (ftnlen)8);
	e_wsfe();
    }

/*        ITERATION FINISHED -COPY RESULTS INTO ARRAY PARA */

L160:
    if (*monit > 0) {
	s_wsfe(&io___102);
	e_wsfe();
    }
    para[1] = u;
    para[2] = a;
    para[3] = g;
    return 0;

L170:
    for (i__ = 1; i__ <= 3; ++i__) {
/* L180: */
	para[i__] = zero;
    }
    return 0;

} /* mlegev_ */

