#include "cplot.h"

FILE *in, *ps, *out;

int main(int argc, char *argv[]){

  char     line[80];       // String from input file
  Curve *curves[NCURVES];
  memset(curves,0,sizeof(curves));

  if( argc < 2 ){
    printf("Syntax:\ncplot filename psfilename \n");	
    return 0;
  }	
  in = fopen( argv[1], "r" );			

  if (argc < 3){
    initps("output.eps");
  }
  else initps(argv[2]);
  out = fopen("output.txt", "w");	

  for (;;){ // Process input lines
    if (fscanf(in, "%s", line) == EOF){
      printf("EOF reached.\n");	
      fclose(in);	
      fclose(ps);	
      fclose(out);
      return 0;
    }
    line[4] = 0;	// Ignore everying past 4th character
    int i;
    for(i=0; i < 4; i++)line[i] = tolower(line[i]);	 // Convert to lower case

    if (!strcmp(line, "bord")){ 
      fprintf(ps,"border stroke\n"); 
    } 
    //Plot Bezier
    else if (!strcmp(line, "cplo")){ 
      int nCurve;
      fscanf(in, "%d", &nCurve);
      if(!curves[nCurve]){
        printf("ERROR: Cannot Plot Curve %d. Curve has not been initialized\n", nCurve);
      }
      else{
        plotBezier(curves[nCurve]);
      }
    }
    //Plot Circle
    else if (!strcmp(line, "circ")){
      double x , y, r;
      fscanf(in, "%lf %lf %lf", &x, &y, &r);
      circle(x, y, r);
    }
    // Set color.  Next line: red green blue  (0 to 1)
    else if (!strcmp(line, "colo")){
      double r, g, b;					
      fscanf(in, "%lf %lf %lf", &r, &g, &b);
      fprintf(ps,"%lf %lf %lf setrgbcolor\n", r,g,b);
    }
    // Plot control polygon.  Next line:  curve number
    else if (!strcmp(line, "cppl")){
      int nCurve;
      fscanf(in, "%d", &nCurve);
      if(!curves[nCurve]){
        printf("ERROR: Cannot Plot Curve %d's Polygon. Curve has not been initialized\n", nCurve);
      }
      plotControlPolygon(curves[nCurve]);
    }
    // Exit the program; close files.
    else if (!strcmp(line, "exit")){
      printf("File processed successfully!\n");
      fclose(in);	fclose(ps);  fclose(out); return 0; 
    }
    //Store Bezier							
    else if (!strcmp(line, "stor")){
      int nCurve;
      fscanf(in,"%d",&nCurve); 
      curves[nCurve] = readCurve();			
    }
    // Print text.  Next line = font size, (x,y). Next line  = text
    else if (!strcmp(line, "text")){
      double size, x, y;									
      fscanf(in, "%lf %lf %lf", &size, &x, &y);
      line[0] = 0;
      while(line[0]<33)fscanf(in,"%c",line); // Strip white space
      for(i=1; line[i-1] != 10 && i<80; i++)fscanf(in,"%c",line + i);
      line[i-1] = 0;
      fprintf(ps,"%f %f mv /Times-Roman findfont %f scalefont setfont (%s) show\n", x,y,size,line);
    }
    // Change viewport.  Next line: xmin ymin xmax ymax
    else if (!strcmp(line, "view")){
      double xmin, ymin, xmax, ymax;									   
      fscanf(in, "%lf %lf %lf %lf", &xmin, &ymin, &xmax, &ymax);
      fprintf(ps,"%lf %lf %lf %lf viewport\n", xmin, ymin, xmax, ymax);
    }
    // Change window.  Next line: xmin ymin xmax ymax
    else if (!strcmp(line, "wind")){
      double xmin,ymin,xmax,ymax;
      fscanf(in, "%lf %lf %lf %lf", &xmin, &ymin, &xmax, &ymax);
      fprintf(ps,"%lf %lf %lf %lf window\n", xmin, ymin, xmax, ymax);
    }
    // Change line width.  Next line: width in points
    else if (!strcmp(line, "widt")) {
      double size;									
      fscanf(in, "%lf", &size);
      fprintf(ps, "%lf setlinewidth\n", size);
    }
    // Output all curves
    else if (!strcmp(line, "disp")) {
      for(i=0; i< NCURVES; i++){
        if(curves[i]){
          Curve* c = curves[i];
          printf("Curve %d (Degree %d):\n",i,c->degree);	
          for(int j = 0; j<= c->degree; j++){
            Point* p = c->points[j];
            printf("\t%.3lf,%.3lf,%.3lf\n",p->x,p->y,p->w);
          }
        }		
      }
    }
    // Subdivide curve. Next line: ncurve t nleft nright
    else if (!strcmp(line, "subd")) {
      int ncurve, nleft, nright;
      double t;
      fscanf(in, "%d %lf %d %d", &ncurve, &t, &nleft, &nright);
      Point ***subPoints = deCasteljau(curves[ncurve], t);
      curves[nleft] = leftCurve(subPoints, curves[ncurve]->degree);
      curves[nright] = rightCurve(subPoints, curves[ncurve]->degree);
    }
    // Elevate curve. Next line: ncurve nelev nnewcurve
    else if (!strcmp(line, "elev")) {
      int ncurve, nelev, nnewcurve;
      fscanf(in, "%d %d %d", &ncurve, &nelev, &nnewcurve);
      curves[nnewcurve] = elevate(curves[ncurve], nelev);
    }
    // Plot line from P(t) to its center of curvature for n + 1 values
    else if (!strcmp(line, "ccur")) {
      int ncurve, n;
      fscanf(in, "%d %d", &ncurve, &n);
      plotCenterCurvature(curves[ncurve], n);
    }
    // Print curvature value next to P(t) 
    else if (!strcmp(line, "curv")) {
      int ncurve;
      double t;
      fscanf(in, "%d %lf", &ncurve, &t);
      Curvature cv = curvature(curves[ncurve], t);
      fprintf(ps,"%f %f mv /Times-Roman findfont %f scalefont setfont (%lf) show\n", 
                 cv.pt.x,cv.pt.y,12.0,cv.k);
    }
    // Draw osculating circle
    else if (!strcmp(line, "oscu")) {
      int ncurve;
      double t;
      fscanf(in, "%d %lf", &ncurve, &t);
      Curvature cv = curvature(curves[ncurve], t);
      circle(cv.center.x, cv.center.y, cv.r); 
    }
    // Plot offset curve
    else if (!strcmp(line, "offs")) {
      int ncurve;
      double r;
      fscanf(in, "%d %lf", &ncurve, &r);
      plotOffset(curves[ncurve], r);
    }
    // Store polynomial as explicit bezier
    else if (!strcmp(line, "expl")) {
      int ncurve, n;
      fscanf(in, "%d %d", &ncurve, &n);
      double *fs = new double[n + 1];
      for (int i = 0; i <= n; i++) {
      fscanf(in, "%lf", &fs[i]);
      }
      curves[ncurve] = explicitBezier(n, fs);
    }
    // Create a Bernstein polynomial of degree n whose roots are r1, ... , rn
    else if (!strcmp(line, "r2bp")) {
      int ncurve, n;
      fscanf(in, "%d %d", &ncurve, &n);
      double *fs = new double[n];
      for (int i = 0; i < n; i++) {
        fscanf(in, "%lf", &fs[i]);
      }
      curves[ncurve] = rootsToBernstein(fs, n);
    }
    // Find real roots of explicit polynomial in ncurve
    else if (!strcmp(line, "root")) {
      int ncurve;
      fscanf(in, "%d", &ncurve);
      double *roots = findRoots(curves[ncurve]);

      fprintf(out, "Curve %d: <", ncurve);

      for (int i = 0; i <= curves[ncurve]->degree; i++) {
	if (i) fprintf(out, ", ");
	fprintf(out, "%lf", curves[ncurve]->points[i]->y);
      }

      fprintf(out, ">\nRoots: [");

      for (int i = 0; i < curves[ncurve]->degree; i++) {
	if (!isnan(roots[i])) {
	  if (i) fprintf(out, ", ");
	  fprintf(out, "%lf", roots[i]);
	}
      }

      fprintf(out, "]\n\n");
    }
    else if (!strcmp(line, "intl")) {
      int ncurve, mcurve;
      fscanf(in, "%d %d", &ncurve, &mcurve);
      double *paramValues;
      Point **points;
      intersectCurveWithLine(curves[ncurve], curves[mcurve], &paramValues, &points);
      
      fprintf(out, "Intersection of Curve %d and Line %d: \n", ncurve, mcurve);
      fprintf(out, "Parameter Values: [");
      for (int i = 0; i < curves[ncurve]->degree; i++) {
	if (!isnan(paramValues[i])) {
	  if (i) fprintf(out, ", ");
	  fprintf(out, "%lf", paramValues[i]);
	}
      }

      fprintf(out, "]\nCoordinates: ");
      for (int i = 0; i < curves[ncurve]->degree; i++) {
	if (points[i]) {
	  if (i) fprintf(out, ", ");
	  fprintf(out, "(%lf, %lf)", points[i]->x, points[i]->y);
	  circle(points[i]->x, points[i]->y, 0.1);
	}
      }



    }
    else{ 
      printf("Illegal command: %s\n", line);		
      return 0;
    }
  }//end for(;;)
}// end main()
