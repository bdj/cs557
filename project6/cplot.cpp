#include "cplot.h"

using namespace std;

Point::Point(){
  x = 0;
  y = 0;
  w = 1;
}

Point::Point(double x, double y, double w){
  this->x = x;
  this->y = y;
  this->w = w;
}

Point operator+(Point lhs, const Point& rhs) {
  lhs.x += rhs.x;
  lhs.y += rhs.y;
  lhs.w += rhs.w;
  return lhs;
}

Point operator*(double lhs, Point rhs) {
  rhs.x *= lhs;
  rhs.y *= lhs;
  rhs.w *= lhs; 
  return rhs;
}

Polynomial operator+(Polynomial lhs, const Polynomial& rhs) {
  for (int i = 0; i <= lhs.degree; i++) {
    lhs.coefficients[i] += rhs.coefficients[i];
  }
  return lhs;
}

Polynomial operator-(Polynomial lhs, const Polynomial& rhs) {
  for (int i = 0; i <= lhs.degree; i++) {
    lhs.coefficients[i] -= rhs.coefficients[i];
  }
  return lhs;
}

Polynomial operator*(Polynomial lhs, const Polynomial& rhs) {
  Polynomial result(lhs.degree + rhs.degree);
  for (int i = 0; i <= lhs.degree; i++) {
    for (int j = 0; j <= rhs.degree; j++) {
      result.coefficients[i + j] += nck(lhs.degree, i) * lhs.coefficients[i] * 
                       nck(rhs.degree, j) * rhs.coefficients[j];
    }
  }

  for (int i = 0; i <= result.degree; i++) {
    result.coefficients[i] /= nck(result.degree, i);
  }

  return result;
}

Curve::Curve(int n) {
  degree = n;
  points = new Point*[n+1]; 
}

NURBS::NURBS(int n) {
  num_points = n;
  knots = new double[n + 2];
  points = new Point*[n];
}

Polynomial::Polynomial(int n) {
  degree = n;
  coefficients = new double[n + 1];
}

Deform::Deform(int m, int n, double xmin, double ymin, double xmax, double ymax) {
  this->m = m;
  this->n = n;
  this->xmin = xmin;
  this->ymin = ymin;
  this->xmax = xmax;
  this->ymax = ymax;
  matrix = new Point**[m + 1];
  for (int i = 0; i <= m; i++) {
    matrix[i] = new Point*[n + 1];
  }
}

Vector::Vector() {
  x = 0;
  y = 0;
}

Vector::Vector(double x, double y) {
  this->x = x;
  this->y = y;
}

Vector::Vector(Point *a, Point *b) {
  x = b->x / b->w - a->x / a->w;
  y = b->y / b->w - a->y / a->w;
}

Vector::Vector(Point *p) {
  x = p->x / p->w;
  y = p->y / p->w;
}

Curvature::Curvature(double r, double k, Vector center, Vector pt) {
  this->r = r;
  this->k = k;
  this->center = center;
  this->pt = pt;
}

Vector operator+(Vector lhs, const Vector& rhs) {
  lhs.x += rhs.x;
  lhs.y += rhs.y;
  return lhs;
}

Vector operator-(Vector lhs, const Vector& rhs) {
  lhs.x -= rhs.x;
  lhs.y -= rhs.y;
  return lhs;
}

Vector operator*(double lhs, Vector rhs) {
  rhs.x *= lhs;
  rhs.y *= lhs;
  return rhs;
}

double lengthSq(const Vector& v) {
  return v.x * v.x + v.y * v.y;
}

Vector perp(Vector v) {
  double z = v.x;
  v.x = -v.y;
  v.y = z;
  return v;
}

double dot(const Vector& a, const Vector& b) {
  return a.x * b.x + a.y * b.y;
}

double dotP(Point *a, Point *b) {
  return a->x * b->x + a->y * b->y + a->w * b->w;
}

Vector norm(Vector v) {
  double len = sqrt(lengthSq(v));
  v.x /= len;
  v.y /= len;
  return v;
}

Curve *leftCurve(Point ***subPoints, int degree) {
  Curve *result = new Curve(degree);
  for (int j = 0; j <= degree; j++) {
    result->points[j] = subPoints[j][0]; 
  }
  return result;
}

Curve *rightCurve(Point ***subPoints, int degree) {
  Curve *result = new Curve(degree);
  for (int i = 0; i <= degree; i++) {
    result->points[i] = subPoints[degree - i][i]; 
  }
  return result;
}

Point ***deCasteljau(Curve *c, double t) {
  int n = c->degree + 1;
  double u = 1 - t;
  Point** points = new Point*[(n + 1) * n / 2];
  Point*** subPoints = new Point**[n];
  // setup index array
  for (int i = 0, m = 0; i < n; m += (n - i), i++) {
    subPoints[i] = points + m;
  }
  // put initial points in
  for (int i = 0; i < n; i++) {
    subPoints[0][i] = c->points[i];
  }
  // calculate subdivisions
  for (int j = 1; j < n; j++) {
    for (int i = 0; i < n - j; i++) {
      subPoints[j][i] = new Point(u * *(subPoints[j - 1][i]) + t * *(subPoints[j - 1][i + 1])); 
    }
  }
  return subPoints;
}

Curve *elevate(Curve *c, int nelev) {
  int m = c->degree + nelev;
  Curve *curve = c;
  Curve *elevated;
  // elevate nelev times
  for (int n = c->degree; n < m; n++) {
    elevated = new Curve(n + 1);
    elevated->points[0] = curve->points[0];
    elevated->points[n + 1] = curve->points[n];
    for (int i = 1; i <= n; i++) {
      double a = i / (n + 1.0);
      elevated->points[i] = new Point(a * *(curve->points[i - 1]) + (1 - a) * *(curve->points[i]));
    }
    curve = elevated;
  }
  
  return elevated;
}

Curvature curvature(Curve *co, double t) {
  Curve *c;
  int i, j, k;
  if (t == 1.0) {
    i = co->degree;
    j = i - 1;
    k = i - 2;
    c = co;
  } else {
    i = 0;
    j = 1;
    k = 2; 
    c = rightCurve(deCasteljau(co, t), co->degree);
  }
  Vector a(c->points[i], c->points[j]);
  Vector np = perp(norm(a));
  Vector b(c->points[j], c->points[k]);
  double h = dot(b, np);
  double ku = (c->points[i]->w * c->points[k]->w / (c->points[j]->w * c->points[j]->w)) *
              (c->degree - 1) / c->degree *
              h / lengthSq(a);
  double r  = 1 / ku;
  Vector pt = Vector(c->points[i]);
  Vector center = pt + r * np;
  if (ku < 0) {
    ku = -ku;
    r = -r;
  }
  return Curvature(r, ku, center, pt);
}

Curve* readCurve(){
  double x, y, w;
  int d;
  Curve* c;
  fscanf(in,"%d",&d);
  c = new Curve(d);
  for(int i = 0; i <= d; i++){
    fscanf(in, "%lf %lf %lf", &x, &y, &w);	
    c->points[i] = new Point(x*w,y*w,w);
  }	
  return c;
}

NURBS* readNURBS(int m) {
  NURBS* n = new NURBS(m);

  int k;
  for (int i = 0; i < m + 2; i++) {
    fscanf(in, "%d", &k);
    n->knots[i] = k;
  }

  double x, y, w;
  for (int i = 0; i < m; i++) {
    fscanf(in, "%lf %lf %lf", &x, &y, &w);
    n->points[i] = new Point(x * w, y * w, w);
  }
  
  return n;
}

void plotBezier(Curve* c){
  fprintf(ps, "[");
  for(int i = 0; i <= c->degree; i++){
    fprintf(ps, "[%lf %lf %lf]", (c->points)[i]->x, (c->points)[i]->y, (c->points)[i]->w);
  }
  fprintf(ps, "] cplot\n");
}

void circle(double x, double y, double r){
  fprintf(ps, "%lf %lf %lf circ\n", x,y,r);
}

void plotControlPolygon(Curve* c){
  fprintf(ps, "%lf %lf mv  ", (c->points)[0]->x/(c->points)[0]->w, (c->points)[0]->y/(c->points)[0]->w);
  for(int i = 1; i <= c->degree; i++){
    fprintf(ps, "%lf %lf ln  ", (c->points)[i]->x/(c->points)[i]->w, (c->points)[i]->y/(c->points)[i]->w);
  }
  fprintf(ps, "stroke\n");
}

void plotControlPolygon(NURBS* n){
  fprintf(ps, "%lf %lf mv  ", (n->points)[0]->x/(n->points)[0]->w, (n->points)[0]->y/(n->points)[0]->w);
  for(int i = 1; i < n->num_points; i++){
    fprintf(ps, "%lf %lf ln  ", (n->points)[i]->x/(n->points)[i]->w, (n->points)[i]->y/(n->points)[i]->w);
  }
  fprintf(ps, "stroke\n");
}

void plotCenterCurvature(Curve * c, int n) {
  for (int i = 0; i <= n; i++) {
      double t = ((double)i)/n;
      Curvature cv = curvature(c, t);
      fprintf(ps, "%lf %lf mv  ", cv.pt.x, 
                                  cv.pt.y);
      fprintf(ps, "%lf %lf ln  ", cv.center.x, cv.center.y);
      fprintf(ps, "stroke\n");
  }
}

void plotOffset(Curve * c, double r) {
   Vector start = Vector(c->points[0]) + r * perp(norm(Vector(c->points[0], c->points[1])));
   fprintf(ps, "%lf %lf mv  ", start.x, start.y);
   for (int i = 1; i < 100; i++) {
     double t = ((double)i)/100;
     Curve *ct = rightCurve(deCasteljau(c, t), c->degree);
     Vector a(ct->points[0], ct->points[1]);
     Vector np = perp(norm(a));
     Vector off = Vector(ct->points[0]) + r * np; 
     fprintf(ps, "%lf %lf ln  ", off.x, off.y);
   }
   // last point
   Vector a(c->points[c->degree], c->points[c->degree - 1]);
   Vector n = norm(a);
   double z = n.x;
   n.x = n.y;
   n.y = -z;
   Vector off = Vector(c->points[c->degree]) + r * n;
   fprintf(ps, "%lf %lf ln  ", off.x, off.y);
   fprintf(ps, "stroke\n");
}

int nck(int n, int k ) {
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

Curve *explicitBezier(int n, double *fs) {
  Curve *result = new Curve(n);
  for (int i = 0; i <= n; i++) {
    fs[i] /= nck(n, i);
  } 
  
  double **bs = new double*[n];
  bs[0] = fs;
  for (int i = 1; i <= n; i++) {
    bs[i] = new double[n - i + 1];
    for (int j = n - i; j >= 0; j--) {
      bs[i][j] = bs[i - 1][j + 1] + bs[i - 1][j];
    }
  } 

  for (int i = 0; i <= n; i++) {
    result->points[i] = new Point(((double)i) / n, bs[i][0], 1);
  }
  return result;
}

Curve *addRoot(Curve *c, double r) {
  int n = c->degree;
  Curve *result = new Curve(n + 1);
  result->points[0] = new Point(0, -r * c->points[0]->y, 1); 
  for (int i = 1; i <= n; i++) {
    result->points[i] = new Point((double)i / (n + 1), 
        (nck(n, i) * c->points[i]->y * -r + 
         nck(n, i - 1) * c->points[i - 1]->y * (1 - r)) /
        nck(n + 1, i),
        1);
  }
  result->points[n+1] = new Point(1, c->points[n]->y * (1 - r), 1);

  return result;
}

Curve *rootsToBernstein(double *roots, int n) {
  Curve *c = new Curve(1);
  c->points[0] = new Point(0, -roots[0], 1);
  c->points[1] = new Point(1, 1 - roots[0], 1);
  for (int i = 1; i < n; i++) {
    c = addRoot(c, roots[i]);
  }
  return c;
}

Point *cross(Point *p0, Point *p1) {
  return new Point(
      p0->y * p1->w - p0->w * p1->y,
      p0->w * p1->x - p0->x * p1->w,
      p0->x * p1->y - p0->y * p1->x);
}

Point *unweight(Point *p) {
  return new Point(p->x / p->w, p->y / p->w, 1);
}

bool sameSign(double a, double b) {
  return (a > 0 && b > 0) || (a < 0 && b < 0);
}

double nearestHullIntersect(Curve *c) {
  double min_x = INFINITY;
  Point A(0, 1, 0);

  for (int i = 0; i <= c->degree; i++) {
    if (!sameSign(c->points[0]->y, c->points[i]->y))
      min_x = min(min_x, unweight(cross(cross(c->points[0], c->points[i]), &A))->x);
  }
  return min_x;
}

Curve *findRoot(Curve *c) {
  Point *p;
  do {
    p = c->points[0];
    c = rightCurve(deCasteljau(c, (nearestHullIntersect(c) - p->x) / (1 - p->x)), c->degree) ;
  } while (fabs(p->x - c->points[0]->x) > 1e-20 && sameSign(p->y, c->points[0]->y));

  return c;
}

Curve *deflate(Curve *c) {
  int n = c->degree;
  double r = c->points[0]->x;
  Curve *deflated = new Curve(n - 1);
  for (int i = 0; i < n; i++) {
    deflated->points[i] = new Point(i * (1 - r) / (n - 1) + r, n * c->points[i + 1]->y / (i + 1), 1);
  }
  return deflated;
}

double *findRoots(Curve *c) {
  double *roots = new double[c->degree];
  int n = c->degree;
  for (int i = 0; i < n; i++) {
      c = findRoot(c);
      roots[i] = c->points[0]->x;
      c = deflate(c);
  }
  return roots;
}

void fprintCurve(FILE *out, Curve *c) {
  fprintf(out, "<");
  for (int i = 0; i <= c->degree; i++) {
    if (i) fprintf (out, ", ");
    fprintPoint(out, c->points[i]);
  }
  fprintf(out, ">\n");
}
void fprintPoint(FILE *out, Point *p) {
    fprintf(out, "(%.2lf, %.2lf, %.2lf)", p->x, p->y, p->w);
}

void intersectCurveWithLine(Curve *c, Curve *l, double **paramValues, Point ***points) {
  Point *L = cross(l->points[0], l->points[1]);
  Curve *proj = new Curve(c->degree);
  for (int i = 0; i <= c->degree; i++) {
    proj->points[i] = new Point((double)i / c->degree, dotP(c->points[i], L), 1);
  }

  *paramValues = findRoots(proj);

  *points = new Point*[c->degree];

  for (int i = 0; i < c->degree; i++) {
    if (!isnan((*paramValues)[i])) {
      (*points)[i] = unweight(rightCurve(deCasteljau(c, (*paramValues)[i]), c->degree)->points[0]); 
    } else {
      (*points)[i] = NULL;
    }
  }
}

NURBS *insertKnot(NURBS *n, double t) {
  NURBS *result = new NURBS(n->num_points + 1);
  int i;
  for (i = 0; i < n->num_points + 2; i++) {
    if (t < n->knots[i]) {
      break;
    } else {
      result->knots[i] = n->knots[i];
      if (i >= 1) {
        result->points[i - 1] = n->points[i - 1];
      }
    }
  }
  // handle inserting identical to end conditions
  if (i == n->num_points + 2) {
    i -= 3; 
  }
  result->knots[i] = t;
  result->knots[i + 1] = n->knots[i];

  result->points[i - 2] = new Point(
      1/(n->knots[i] - n->knots[i - 3]) *
      ((n->knots[i] - t) * *(n->points[i - 3]) +
       (t - n->knots[i - 3]) * *(n->points[i - 2])));

  result->points[i - 1] = new Point(
      1/(n->knots[i + 1] - n->knots[i - 2]) *
      ((n->knots[i + 1] - t) * *(n->points[i - 2]) +
       (t - n->knots[i - 2]) * *(n->points[i - 1])));
  
  result->points[i] = new Point(
      1/(n->knots[i + 2] - n->knots[i - 1]) *
      ((n->knots[i + 2] - t) * *(n->points[i - 1]) +
       (t - n->knots[i - 1]) * *(n->points[i])));

  for (i += 2; i < result->num_points + 2; i++) {
    result->knots[i] = n->knots[i - 1];
    if (i <= result->num_points) {
      result->points[i - 1] = n->points[i - 2];
    }
  }
  return result;
}

Curve *extractCurve(NURBS *n, int i) {
  int k = 1;
  int j;
  for (j = 1; j < n->num_points + 2; j++) {
    if (n->knots[j - 1] != n->knots[j]) {
      if (i == k) {
        break;
      }
      k++; 
    }
  }

  NURBS *m = insertKnot(n, n->knots[j - 1]);
  m = insertKnot(m, n->knots[j - 1]);
  m = insertKnot(m, n->knots[j]);
  m = insertKnot(m, n->knots[j]);
  Curve *result = new Curve(3);
  result->points[0] = m->points[j - 1];
  result->points[1] = m->points[j];
  result->points[2] = m->points[j + 1];
  result->points[3] = m->points[j + 2];
  return result;
}

Point *lij(int i, int j, Curve *c) {
  double scale = nck(3, i) * nck(3, j);
  return new Point(scale * *cross(c->points[i], c->points[j]));
}

Polynomial *Lij(int i, int j, Curve *P, Curve *Q) {
  Point *l_ij = lij(i, j, P);
  Polynomial *L_ij = new Polynomial(Q->degree);
  for (int i = 0; i <= Q->degree; i++) {
    L_ij->coefficients[i] = dotP(l_ij, Q->points[i]);
  }
  
  return L_ij;
}

Polynomial *Lc(Curve *P, Curve *Q) {
  Point *l_30 = lij(3, 0, P);
  Point *l_21 = lij(2, 1, P);
  Polynomial *L_c = new Polynomial(Q->degree);
  for (int i = 0; i <= Q->degree; i++) {
    L_c->coefficients[i] = (l_30->x + l_21->x) * Q->points[i]->x + 
                           (l_30->y + l_21->y) * Q->points[i]->y + 
                           (l_30->w + l_21->w) * Q->points[i]->w;
  }

  return L_c;
} 

Polynomial *gt(Curve *P, Curve *Q) {
  Polynomial L_32 = *Lij(3, 2, P, Q);
  Polynomial L_31 = *Lij(3, 1, P, Q);
  Polynomial L_30 = *Lij(3, 0, P, Q);
  Polynomial L_20 = *Lij(2, 0, P, Q);
  Polynomial L_10 = *Lij(1, 0, P, Q);
  Polynomial L_c = *Lc(P, Q);

  return new Polynomial(L_32 * (L_c * L_10 - L_20 * L_20) -
                        L_31 * (L_31 * L_10 - L_30 * L_20) +
                        L_30 * (L_31 * L_20 - L_30 * L_c));
}

Curve *polynomialToExplicit(Polynomial *p) {
  int n = p->degree;
  Curve *result = new Curve(n);
  for (int i = 0; i <= n; i++) {
    result->points[i] = new Point((double)i / n, p->coefficients[i], 1);
  }
  return result;
}

void intersectCurves(Curve *p, Curve *q, double **paramValues, Point ***points) {

  Curve *c = polynomialToExplicit(gt(p, q));
  *paramValues = findRoots(c);

  *points = new Point*[c->degree];

  for (int i = 0; i < c->degree; i++) {
    if (!isnan((*paramValues)[i])) {
      (*points)[i] = unweight(rightCurve(deCasteljau(q, (*paramValues)[i]), q->degree)->points[0]); 
    } else {
      (*points)[i] = NULL;
    }
  }
}

Point *horner(Curve* c, double t) {
  double u, bc, tn;
  Point tmp;
  int n = c->degree;

  u = 1 - t;
  bc = 1;
  tn = 1;
  tmp = u * (*c->points[0]);
  for (int i = 1; i < n; i++) {
    tn = tn*t;
    bc = bc*(n-i+1)/i;
    tmp = u*(tmp + tn*bc*(*c->points[i]));
  }
  return new Point(tmp + tn*t*(*c->points[n]));
}

Point *stCoords(Deform *deform, Point *p) {
  return new Point((p->x/p->w - deform->xmin) / (deform->xmax - deform->xmin), 
                   (p->y/p->w - deform->ymin) / (deform->ymax - deform->ymin),
                   1);
}

Curve *hCurve(Deform *deform, int n) {
  Curve *result = new Curve(deform->m);
  for (int i = 0; i <= deform->m; i++) {
    result->points[i] = deform->matrix[i][n];
  }
  return result;
}

Curve *vCurve(Deform *deform, int m) {
  Curve *result = new Curve(deform->n);
  for (int i = 0; i <= deform->n; i++) {
    result->points[i] = deform->matrix[m][i];
  }
  return result;
}

Curve *sCurve(Deform *deform, double s) {
  Curve *result = new Curve(deform->n);
  for (int i = 0; i <= deform->n; i++) {
    result->points[i] = horner(hCurve(deform, i), s);
  }
  return result;
}

Curve *tCurve(Deform *deform, double t) {
  Curve *result = new Curve(deform->m);
  for (int i = 0; i <= deform->m; i++) {
    result->points[i] = horner(vCurve(deform, i), t);
  }
  return result;
}

Curve **grid(Deform *deform, int nx, int ny) {
  Curve **result = new Curve*[nx + ny + 2];
  for (int i = 0; i <= nx; i++) {
    result[i] = sCurve(deform, (double)i / nx);
  }
  for (int i = 0; i <= ny; i++) {
    result[i + nx + 1] = tCurve(deform, (double)i / ny);
  }
  return result;
}

void plotDeformCurve(Deform *deform, Curve * c) {
  Point *st = stCoords(deform, c->points[0]);
  Point *p = horner(sCurve(deform, st->x), st->y);
  fprintf(ps, "%lf %lf mv  ", p->x/p->w, p->y/p->w);
  for (int i = 1; i <= 100; i++) {
    double t = ((double)i)/100;
    st = stCoords(deform, horner(c, t));
    p = horner(sCurve(deform, st->x), st->y);
    fprintf(ps, "%lf %lf ln  ", p->x/p->w, p->y/p->w);
  }
  fprintf(ps, "stroke\n");
}

void initps(const char *psfile){  //Initialize the PostScript file	
  printf("In initps: %s\n",psfile);	
  ps = fopen(psfile, "w");	
  fprintf(ps,"%%!PS-Adobe-2.0 EPSF-1.2\n%%%%BoundingBox: 0 0 432 432\n/bb [0 0 432 432] def\n");
  fprintf(ps,"/window {/wy1 exch def  /wx1 exch def  /wy0 exch def  /wx0 exch def} def\n/viewport {\n"); 
  fprintf(ps,"	/vy1 exch bb 3 get bb 1 get sub mul bb 1 get add def	/vx1 exch bb 2 get bb 0 get sub mul bb 0 get add def\n");
  fprintf(ps,"	/vy0 exch bb 3 get bb 1 get sub mul bb 1 get add def	/vx0 exch bb 2 get bb 0 get sub mul bb 0 get add def\n}def\n");
  fprintf(ps,"0 0 1 1 viewport			0 0 1 1 window\n");
  fprintf(ps,"/wvmap {wy0 sub vy1 vy0 sub wy1 wy0 sub div mul vy0 add exch  wx0 sub vx1 vx0 sub wx1 wx0 sub div mul vx0 add exch} def\n");
  fprintf(ps,"/mv {wvmap moveto} def			/ln {wvmap lineto} def\n");
  fprintf(ps,"/h2c  {dup dup 0 get exch 2 get div exch dup 1 get exch 2 get div} def\n");
  fprintf(ps,"/mvh  {h2c  mv} def		/lnh  {h2c ln} def\n/border {newpath wx0 wy0 mv wx1 wy0 ln wx1 wy1 ln wx0 wy1 ln closepath} def\n");
  fprintf(ps,"%%====================== Functions to draw a Bezier Curve and its control polygon  ====================\n");
  fprintf(ps,"/cplot{5 dict begin	%%  plot a rational Bezier curve with homogeneous control points.\n");
  fprintf(ps,"		/cp exch def\n		cp 0 get  mvh   %% move to first point on curve\n");
  fprintf(ps,"		0 .01 1{  cp exch poft  lnh	}for stroke \nend}  bind def\n");
  fprintf(ps,"/cpplot{5 dict begin %% Draw the control polygon of a Bezier curve\n");
  fprintf(ps,"		/p exch def		p 0 get mvh			1 1 p length 1 sub{p exch get lnh }for\nend}def\n");
  fprintf(ps,"/poft{10 dict begin   %%  evaluate a rational Bezier curve.  Syntax:     curve   t  poft  =>  point\n");
  fprintf(ps,"		/t exch def		/p exch def\n		/deg  p length 1 sub  def		/u 1 t sub def		/bc 1 def		/tn 1 def\n");
  fprintf(ps,"		/tmp p 0 get u scaleh def\n		1  1  deg 1 sub { /i exch def\n");
  fprintf(ps,"			/tn tn t mul def\n			/bc bc deg i sub 1 add mul i div def\n");
  fprintf(ps,"			/tmp  p i get tn bc mul scaleh tmp addh u scaleh def\n		} for\n");
  fprintf(ps,"		p deg get tn t mul scaleh tmp addh\nend}  bind def\n");
  fprintf(ps,"/addh{5 dict begin   %% add two vectors\n");
  fprintf(ps,"	/a exch def		/b exch def 	[a 0 get b 0 get add  a 1 get b 1 get add   a 2 get b 2 get add]  end\n} bind def\n");
  fprintf(ps,"/scaleh{5 dict begin   %% scale a vector :   vec scalar  scaleh\n");
  fprintf(ps,"		/b exch def		/a exch def  [a 0 get b mul   a 1 get b mul  a 2 get b mul]	end\n} bind def\n");
  fprintf(ps,"/circ{5 dict begin /r exch def  /y exch def  /x exch def x r add y mv 3 3 360{dup cos r mul x add exch sin r mul y add ln}for stroke} bind def\n");
  fprintf(ps,"%%=========================================================================\n");
}




