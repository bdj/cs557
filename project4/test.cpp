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

  return 0;
}
