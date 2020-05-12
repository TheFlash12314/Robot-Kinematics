/* Matrix.h
  
  Calculations of matrix values.
  This code is created by TheFLash/YWW.
*/
#ifndef Matrix_h
#define Matrix_h

void invtran(float* Titi, float* Titf);
void tran2pos(float* Ttp, float* Xtp);
void pos2tran(float* Xpt, float* Tpt);
void DH1line(float thetadh, float alfadh, float rdh, float ddh, float* Tdh);
//void MatrixPrint(float* A, int m, int n, String label);
void MatrixCopy(float* A, int n, int m, float* B);
void MatrixMultiply(float* A, float* B, int m, int p, int n, float* C);
void MatrixAdd(float* A, float* B, int m, int n, float* C);
void MatrixSubtract(float* A, float* B, int m, int n, float* C);
void MatrixTranspose(float* A, int m, int n, float* C);
void MatrixScale(float* A, int m, int n, float k);

#endif
