#ifndef APPBUILD
#include "mex.h"
#endif

#include <fstream>
#include <string>
#include <iostream>
using namespace std;



//Test if 2 2D triangles are overlap
//Ref: https://github.com/benardp/contours/blob/master/freestyle/view_map/triangle_triangle_intersection.c

#define ORIENT_2D(a, b, c)  ((a[0]-c[0])*(b[1]-c[1])-(a[1]-c[1])*(b[0]-c[0]))

#define INTERSECTION_TEST_VERTEX(P1, Q1, R1, P2, Q2, R2) {\
  if (ORIENT_2D(R2,P2,Q1) >= 0.0f)\
    if (ORIENT_2D(R2,Q2,Q1) <= 0.0f)\
      if (ORIENT_2D(P1,P2,Q1) > 0.0f) {\
	if (ORIENT_2D(P1,Q2,Q1) <= 0.0f) return 1; \
	else return 0;} else {\
	if (ORIENT_2D(P1,P2,R1) >= 0.0f)\
	  if (ORIENT_2D(Q1,R1,P2) >= 0.0f) return 1; \
	  else return 0;\
	else return 0;}\
    else \
      if (ORIENT_2D(P1,Q2,Q1) <= 0.0f)\
	if (ORIENT_2D(R2,Q2,R1) <= 0.0f)\
	  if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) return 1; \
	  else return 0;\
	else return 0;\
      else return 0;\
  else\
    if (ORIENT_2D(R2,P2,R1) >= 0.0f) \
      if (ORIENT_2D(Q1,R1,R2) >= 0.0f)\
	if (ORIENT_2D(P1,P2,R1) >= 0.0f) return 1;\
	else return 0;\
      else \
	if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) {\
	  if (ORIENT_2D(R2,R1,Q2) >= 0.0f) return 1; \
	  else return 0; }\
	else return 0; \
    else  return 0; \
};

#define INTERSECTION_TEST_EDGE(P1, Q1, R1, P2, Q2, R2) { \
  if (ORIENT_2D(R2,P2,Q1) >= 0.0f) {\
    if (ORIENT_2D(P1,P2,Q1) >= 0.0f) { \
        if (ORIENT_2D(P1,Q1,R2) >= 0.0f) return 1; \
        else return 0;} else { \
      if (ORIENT_2D(Q1,R1,P2) >= 0.0f){ \
	if (ORIENT_2D(R1,P1,P2) >= 0.0f) return 1; else return 0;} \
      else return 0; } \
  } else {\
    if (ORIENT_2D(R2,P2,R1) >= 0.0f) {\
      if (ORIENT_2D(P1,P2,R1) >= 0.0f) {\
	if (ORIENT_2D(P1,R1,R2) >= 0.0f) return 1;  \
	else {\
	  if (ORIENT_2D(Q1,R1,R2) >= 0.0f) return 1; else return 0;}}\
      else  return 0; }\
    else return 0; }}



int ccw_tri_tri_intersection_2d(double p1[2], double q1[2], double r1[2], 
				double p2[2], double q2[2], double r2[2]) {
  if ( ORIENT_2D(p2,q2,p1) >= 0.0f ) {
    if ( ORIENT_2D(q2,r2,p1) >= 0.0f ) {
      if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) return 1;
      else INTERSECTION_TEST_EDGE(p1,q1,r1,p2,q2,r2)
    } else {  
      if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) 
	INTERSECTION_TEST_EDGE(p1,q1,r1,r2,p2,q2)
      else INTERSECTION_TEST_VERTEX(p1,q1,r1,p2,q2,r2)}}
  else {
    if ( ORIENT_2D(q2,r2,p1) >= 0.0f ) {
      if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) 
	INTERSECTION_TEST_EDGE(p1,q1,r1,q2,r2,p2)
      else  INTERSECTION_TEST_VERTEX(p1,q1,r1,q2,r2,p2)}
    else INTERSECTION_TEST_VERTEX(p1,q1,r1,r2,p2,q2)}
};


int tri_tri_overlap_test_2d(double p1[2], double q1[2], double r1[2], 
			    double p2[2], double q2[2], double r2[2]) {
  if ( ORIENT_2D(p1,q1,r1) < 0.0f )
    if ( ORIENT_2D(p2,q2,r2) < 0.0f )
      return ccw_tri_tri_intersection_2d(p1,r1,q1,p2,r2,q2);
    else
      return ccw_tri_tri_intersection_2d(p1,r1,q1,p2,q2,r2);
  else
    if ( ORIENT_2D(p2,q2,r2) < 0.0f )
      return ccw_tri_tri_intersection_2d(p1,q1,r1,p2,r2,q2);
    else
      return ccw_tri_tri_intersection_2d(p1,q1,r1,p2,q2,r2);

};


void getCollisionIdx(const double* base, const double* baseArea, const int N1,
                      const double* test, const double* testArea, const int N2, bool* idx) {
  // idx is already mem allocated 
  for(int i=0; i < N1; i++) {
    double p1[2] = {base[i*6+0], base[i*6+3]};
    double q1[2] = {base[i*6+1], base[i*6+4]};
    double r1[2] = {base[i*6+2], base[i*6+5]};
    double totalTestArea = 0;
    for(int j=0; j < N2; j++) {
      double p2[2] = {test[j*6+0], test[j*6+3]};
      double q2[2] = {test[j*6+1], test[j*6+4]};
      double r2[2] = {test[j*6+2], test[j*6+5]};
      int res = tri_tri_overlap_test_2d(p1, q1, r1, p2, q2, r2);
      if(res == 1) {
        totalTestArea += testArea[j];
      }
    }
    idx[i] = baseArea[i] < totalTestArea; 
  }
}


#ifdef APPBUILD

// get length of array in txt file
int getArrayLength(string filename) {
  std::ifstream infile(filename);
  std::string line;
  int N = 0;
  while(std::getline(infile, line)) {
      if(line.length() > 0) N++;
  }
  return N;
}

// read array
double* getArray(string filename, int N) {
  std::ifstream infile(filename);
  std::string line;
  double* a = new double[N];
  float v;
  for(int i=0; i < N; i++) {
      infile >> v;
      a[i] = v;
  }
  return a;
}

// save bool array
void saveBoolArray(string filename, const bool* a, const int N) {
  std::ofstream outfile(filename);
  for(int i=0; i < N; i++) {
    outfile << a[i] << endl;
  }
}

// print array
void printArray(const double* a, int N) {
  for(int i=0; i < 6; i++) cout << a[i] << " ";
  cout << " ... ";
  for(int i=N-6; i < N; i++) cout << a[i] << " ";
  cout << endl;
}

void printArray(const bool* a, int N) {
  for(int i=0; i <N; i++) {
    if(a[i] == false) {
      cout << (i+1) << " " << a[i] << endl;
    }
  }
}

int main (int argc, char *argv[]) {
  // read number of triangles
  int N1 = getArrayLength("data/basePolyArea.txt"); // number of base triangles
  int N2 = getArrayLength("data/testPolyArea.txt"); // number of test triangles
  cout << N1 << " " << N2 << endl;

  double* basePolygons = getArray("data/basePolygons.txt", N1*6);
  double* basePolyArea = getArray("data/basePolyArea.txt", N1);
  printArray(basePolygons, N1*6);
  printArray(basePolyArea, N1);

  double* testPolygons = getArray("data/testPolygons.txt", N2*6);
  double* testPolyArea = getArray("data/testPolyArea.txt", N2);
  printArray(testPolygons, N2*6);
  printArray(testPolyArea, N2);

  bool* idx = new bool[N1];
  getCollisionIdx(basePolygons, basePolyArea, N1, testPolygons, testPolyArea, N2, idx);
  saveBoolArray("data/idx.txt", idx, N1);

  // clean up
  if(basePolygons) delete []basePolygons;
  if(basePolyArea) delete []basePolyArea;
  if(testPolygons) delete []testPolygons;
  if(testPolyArea) delete []testPolyArea;
  if(idx) delete []idx;

  return 0;
}

#else // mex build

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){
  if(nlhs != 1 || nrhs != 4)
	{
		mexPrintf("Incorrect # of input/output parameters\n");
		mexPrintf("Usage: idx = getCollisionIdx(basePolygons, basePolyArea, testPolygons, testPolyArea)");
		return;
	}

  // retrieve input points
	double* basePolygons = mxGetPr(prhs[0]);
	double* basePolyArea = mxGetPr(prhs[1]);
  double* testPolygons = mxGetPr(prhs[2]);
	double* testPolyArea = mxGetPr(prhs[3]);
	const mwSize *d1 = mxGetDimensions(prhs[1]);
	int N1 = (int)d1[0];
  const mwSize *d2 = mxGetDimensions(prhs[3]);
	int N2 = (int)d2[0];

  //associate outputs
	plhs[0] = mxCreateLogicalMatrix(N1, 1);
	bool* outIdx = mxGetLogicals(plhs[0]);

  getCollisionIdx(basePolygons, basePolyArea, N1, testPolygons, testPolyArea, N2, outIdx);
}

#endif