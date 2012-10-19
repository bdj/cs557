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
  double min_t = INFINITY;
  Point A(0, 1, 0);

  for (int i = 0; i <= c->degree; i++) {
    if (!sameSign(c->points[0]->y, c->points[i]->y))
      min_t = min(min_t, unweight(cross(cross(c->points[0], c->points[i]), &A))->x);
  }

  return min_t;
}

double findRoot(Curve *c) {
  double root = 1, new_root = 0;
  int i = 0; 
  while (fabs(root - new_root) > 1e-14) {
    root = new_root;
    new_root = nearestHullIntersect(c);
    printf("new_root = %lf\n", new_root);
    c = rightCurve(deCasteljau(c, (new_root - root) / (1 - root)), c->degree);
    i++;
  }
  printf("did it %d times, found %lf\n", i, root);
  return root;
}

Curve *deflate(Curve *c, double r) {

  return c;
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

  double roots2[] = {-20, .25, .5, .75, 20};

  printf("nearest hull intersect: %lf\n", nearestHullIntersect(rightCurve(deCasteljau(rootsToBernstein(roots2, 5), nearestHullIntersect(rootsToBernstein(roots2, 5))), 5)));

  findRoot(rootsToBernstein(roots2, 5));

  return 0;
}
