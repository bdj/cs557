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

Curve::Curve(int n) {
  degree = n;
  points = new Point*[n+1]; 
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




