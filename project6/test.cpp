#include "cplot.h"
FILE *in, *ps; 
using namespace std;
void printPoint(Point *p) {
    printf("(%.2lf, %.2lf, %.2lf)", p->x, p->y, p->w);
}
void printCurve(Curve *c) {
  printf("<");
  for (int i = 0; i <= c->degree; i++) {
    if (i) printf (", ");
    printPoint(c->points[i]);
  }
  printf(">\n");
}

void printPolynomial(Polynomial *p) {
  printf("<");
  for (int i = 0; i <= p->degree; i++) {
    if (i) printf (", ");
    printf("%.2lf", p->coefficients[i]);
  }
  printf(">\n");
}

int main(int argc, char **argv) {
  Curve *c = new Curve(1);
  c->points[0] = new Point(0, -1, 1);
  c->points[1] = new Point(1, 0, 1);

  printCurve(addRoot(addRoot(c, 2), 3));
  printf("\n");
  double roots[] = {1, 2, 3};
  printCurve(rootsToBernstein(roots, 3));

  printf("\ncross\n");
  printPoint(cross(new Point(4, 2, 2), new Point(6, 5, 1)));
  printf("\n");

  printPoint(cross(new Point(0, 1, 1), new Point(1, -1, 1)));
  printf("\n");

  printf("\nline intersect:\n");
  printPoint(unweight(cross(new Point(0, 1, 1), new Point(0, 1, 0))));
  printf("\n");

  double roots2[] = {-10, .25, .5, .75, 10};
  Curve *c2 = rootsToBernstein(roots2, 5); 

  printf("nearest hull intersect: %lf\n", nearestHullIntersect(rightCurve(deCasteljau(c2, nearestHullIntersect(c2)), 5)));

  Curve *root = findRoot(rootsToBernstein(roots2, 5));

  Curve *c3 = deflate(root);

  Curve *root2 = findRoot(c3);
  printf("Second root: %lf\n", root2->points[0]->x);

  Curve *c4 = deflate(root2);

  printf("4th root: %lf\n", findRoot(deflate(findRoot(c4)))->points[0]->x);


  Curve *c5 = rootsToBernstein(roots2, 5);
  double *roots3 = findRoots(c5);

  for (int i = 0; i < c5->degree; i++) {
    printf("root[%d] = %lf\n", i, roots3[i]);
  }

  double *ps;
  Point **points;
  intersectCurveWithLine(c5, c, &ps, &points);
  for (int i = 0; i < c5->degree; i++) {
    if (!isnan(ps[i])) {
      printf("%d: %lf, (%lf, %lf)\n", i, ps[i], points[i]->x, points[i]->y);
    }
  }

  NURBS *n = new NURBS(5);
  double knots[] = {0, 0, 0, 1, 3, 3, 3};
  n->knots = knots;
  for (int i = 0; i < n->num_points; i++) {
    n->points[i] = new Point(i, i % 2, 1);
  }

  NURBS *m = insertKnot(n, 3);
  printf("[");
  for (int i = 0; i < m->num_points + 2; i++) {
    if (i) printf(", ");
    printf("%lf", m->knots[i]);
  }
  printf("]\n");

    printf("<");
  for (int i = 0; i < m->num_points; i++) {
    if (i) printf(", ");
    printPoint(m->points[i]);
  }
  printf(">\n");


  Curve *c_imp = new Curve(3);
  c_imp->points[0] = new Point(0, 0, 1);
  c_imp->points[1] = new Point(0, 1, 1);
  c_imp->points[2] = new Point(1, 0, 1);
  c_imp->points[3] = new Point(1, 1, 1);


  Curve *c_imp2 = new Curve(1);
  c_imp2->points[0] = new Point(0, 0, 1);
  c_imp2->points[1] = new Point(1, 1, 1);



  Point *l21 = lij(2, 1, c_imp);
  printPoint(l21);
  printf("\n");

  Polynomial *L21 = Lij(2, 1, c_imp, c_imp2);
  printPolynomial(L21);

  Polynomial *L_c = Lc(c_imp, c_imp2);
  printPolynomial(L_c);

  Polynomial P1(2);
  P1.coefficients[0] = 1;
  P1.coefficients[1] = 2;
  P1.coefficients[2] = 3;
  Polynomial P2(2);
  P2.coefficients[0] = 2;
  P2.coefficients[1] = 3;
  P2.coefficients[2] = 4;

  printPolynomial(new Polynomial(P1 * P2));

  Polynomial *g_t = gt(c_imp, c_imp2);
  printPolynomial(g_t);

  double *roots4 = findRoots(polynomialToExplicit(g_t));

  for (int i = 0; i < g_t->degree; i++) {
    printf("root[%d] = %lf\n", i, roots4[i]);
  }

  return 0;
}
