#ifndef CPLOT_H
#define CPLOT_H
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <iostream>
#include <string>
#include <utility>
#include <cstring>
#define NCURVES 100
class Point {
  public:
    double x, y, w;
    Point();
    Point(double,double,double);
};

class Vector {
  public:
    double x, y;
    Vector();
    Vector(double,double);
    Vector(Point *, Point *);
    Vector(Point *);
};

class Curve {
  public:
    int degree;
    Point **points;
    Curve(int);
};

class NURBS {
  public:
    int num_points;
    double *knots;
    Point **points;
    NURBS(int);
};

class Curvature {
   public:
     double r;
     double k;
     Vector center;
     Vector pt;
     Curvature(double,double,Vector,Vector);
};

class Polynomial {
  public:
    int degree;
    double *coefficients;
    Polynomial(int);
};

class Deform {
  public:
    int m, n;
    double xmin, ymin, xmax, ymax;
    Point ***matrix;
    Deform(int, int, double, double, double, double);
};

Curve* readCurve();
NURBS* readNURBS(int);
void plotBezier(Curve*);
void circle(double,double,double);
void plotControlPolygon(Curve*);
void plotControlPolygon(NURBS*);
void initps(const char*);
Point*** deCasteljau(Curve *, double);
Curve *leftCurve(Point ***, int);
Curve *rightCurve(Point ***, int);
Curve *elevate(Curve *c, int nelev);
Curvature curvature(Curve *, double t);
void plotCenterCurvature(Curve *, int);
void plotOffset(Curve * c, double r);
int nck(int, int);
Curve *explicitBezier(int, double *);
Curve *addRoot(Curve *, double);
Curve *rootsToBernstein(double *, int);
Point *cross(Point *p0, Point *p1);
Point *unweight(Point *p);
bool sameSign(double a, double b);
double nearestHullIntersect(Curve *c);
Curve *findRoot(Curve *c);
Curve *deflate(Curve *c);
double *findRoots(Curve *c);
void fprintCurve(FILE *, Curve *c);
void fprintPoint(FILE *, Point *);
void intersectCurveWithLine(Curve *c, Curve *l, double **paramValues, Point ***points);
NURBS *insertKnot(NURBS *, double);
Curve *extractCurve(NURBS *, int);
Point *lij(int, int, Curve *);
Polynomial *Lij(int, int, Curve *, Curve *);
Polynomial *Lc(Curve *, Curve *);
Polynomial operator*(Polynomial lhs, const Polynomial& rhs);
Polynomial *gt(Curve *, Curve *);
Curve *polynomialToExplicit(Polynomial *);
void intersectCurves(Curve *p, Curve *q, double **paramValues, Point ***points);
Curve **grid(Deform *, int, int);
void plotDeformCurve(Deform *deform, Curve * c);

extern FILE *in, *ps;

#endif
