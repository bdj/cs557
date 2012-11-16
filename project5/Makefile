

cplot: cplot.o run.o
	g++ cplot.o run.o -o cplot

test: cplot.o test.o
	g++ cplot.o test.o -o cplot-test
	./cplot-test

.cpp.o:
	g++ -c -g -Wall $< -o $@

clean:
	rm -f *.o
	rm -f cplot
	rm -f cplot-test
