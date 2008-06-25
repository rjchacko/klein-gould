#define L1 16 // "x" dimension
#define L2 16 // "y" dimension
#define L3 16 // "z" dimension
#define L4 64 // "time" dimension
#define L1h (L1/2) // half of the full "x" dimension, useful for even/odd lattice indexing

#define N (L1*L2*L3*L4) // total number of lattice points
#define Nh (N/2) // total number of even/odd lattice points

#define GAUGE_FIXED 1 // gauge chosen so that which most temporal links are identity
#define SPATIAL_SCALING 2.38 // scale gauge links in spatial dimensions
#define TIME_SYMMETRY -1 // if -1, make gauge field antisymmetric in time

#define packedGaugeSiteSize 12 // real numbers per link, using SU(3) reconstruction
#define gaugeSiteSize 18 // real numbers per link
#define spinorSiteSize 24 // real numbers per spinor

#define EVEN_ODD


// ---------- blas_cuda.cu ----------

extern "C" void zeroCuda(float* dst, int cnt);
extern "C" void copyCuda(float* dst, float *src, int len);

extern "C" void axpbyCuda(float a, float *x, float b, float *y, int len);
extern "C" void axpyCuda(float a, float *x, float *y, int len);
extern "C" void xpayCuda(float *x, float a, float *y, int len);
extern "C" void mxpyCuda(float *x, float *y, int len);

extern "C" void axpyZpbxCuda(float a, float *x, float *y, float *z, float b, int len);
extern "C" float axpyNormCuda(float a, float *x, float *y, int len);

extern "C" float sumCuda(float *a, int n);
extern "C" float normCuda(float *a, int n);
extern "C" float reDotProductCuda(float *a, float *b, int n);

extern "C" void blasTest();
extern "C" void axpbyTest();


// ---------- blas_reference.cpp ----------

extern "C" void zero(float* a, int cnt);
extern "C" void copy(float* a, float *b, int len);

extern "C" void ax(float a, float *x, int len);

extern "C" void axpbyCuda(float a, float *x, float b, float *y, int len);
extern "C" void axpy(float a, float *x, float *y, int len);
extern "C" void xpay(float *x, float a, float *y, int len);
extern "C" void mxpy(float *x, float *y, int len);

extern "C" float norm(float *vector, int len);
extern "C" float reDotProduct(float *v1, float *v2, int len);
extern "C" float imDotProduct(float *v1, float *v2, int len);
extern "C" double normD(float *vector, int len);
extern "C" double reDotProductD(float *v1, float *v2, int len);
extern "C" double imDotProductD(float *v1, float *v2, int len);


// ---------- dslash_test.cpp ----------

extern "C" void dslashTest();


// ---------- dslash_cuda.cu ----------

typedef struct CudaPGauge_s *CudaPGauge;
typedef struct CudaPSpinor_s *CudaPSpinor;

typedef struct {
    CudaPSpinor odd;
    CudaPSpinor even;
} CudaFullSpinor;

typedef struct {
    CudaPGauge odd;
    CudaPGauge even;
} CudaFullGauge;


extern "C" CudaPSpinor allocateParitySpinor();

extern "C" CudaFullGauge loadGaugeField(float **gauge);
extern "C" CudaPSpinor loadParitySpinor(float *spinor);
extern "C" CudaFullSpinor loadSpinorField(float *spinor);

extern "C" void freeGaugeField(CudaFullGauge gauge);
extern "C" void freeParitySpinor(CudaPSpinor spinor);
extern "C" void freeSpinorField(CudaFullSpinor spinor);

extern "C" void retrieveParitySpinor(float *res, CudaPSpinor spinor);
extern "C" void retrieveSpinorField(float *res, CudaFullSpinor spinor);

extern "C" void dslashCuda(CudaPSpinor res, CudaFullGauge gauge, CudaPSpinor spinor, int oddBit, int daggerBit);
extern "C" int  dslashCudaSharedBytes();

extern "C" void MatPCCuda(CudaPSpinor outEven, CudaFullGauge gauge, CudaPSpinor inEven, float kappa, CudaPSpinor tmp);
extern "C" void MatPCDagCuda(CudaPSpinor outEven, CudaFullGauge gauge, CudaPSpinor inEven, float kappa, CudaPSpinor tmp);
extern "C" void MatPCDagMatPCCuda(CudaPSpinor outEven, CudaFullGauge gauge, CudaPSpinor inEven, float kappa, CudaPSpinor tmp1, CudaPSpinor tmp2);



// ---------- dslash_reference.cpp ----------

extern "C" void dslashReference(float *res, float **gauge, float *spinorField, 
				int oddBit, int daggerBit);

extern "C" void Mat(float *out, float **gauge, float *in, float kappa);
extern "C" void MatDag(float *out, float **gauge, float *in, float kappa);
extern "C" void MatDagMat(float *out, float **gauge, float *in, float kappa);

extern "C" void MatPC(float *out, float **gauge, float *in, float kappa);
extern "C" void MatPCDag(float *out, float **gauge, float *in, float kappa);
extern "C" void MatPCDagMatPC(float *out, float **gauge, float *in, float kappa);

// ---------- cg_test.cpp ----------

extern "C" void cgTest();


// ---------- cg_cuda.cpp ----------

extern "C" void cgCuda(float *out, float **gauge, float *in, float kappa, float tol);


// ---------- cg_reference.cpp ----------

extern "C" void cgReference(float *out, float **gauge, float *in, float kappa, float tol);


// ---------- qcd.cpp ----------

extern "C" int compareFloats(float *a, float *b, int len, float tol);

extern "C" void   stopwatchStart();
extern "C" double stopwatchReadSeconds();

extern "C" void printSpinor(float *spinor);
extern "C" void printSpinorElement(float *spinor, int X);
extern "C" void printGauge(float *gauge);
extern "C" void printGaugeElement(float *gauge, int X);

extern "C" int fullLatticeIndex(int i, int oddBit);
extern "C" int getOddBit(int X);

extern "C" void applyGaugeFieldScaling(float **gauge);
extern "C" void constructUnitGaugeField(float **gauge);
extern "C" void constructGaugeField(float **gauge);
extern "C" void constructPointSpinorField(float *spinor, int i0, int s0, int c0);
extern "C" void constructSpinorField(float *res);

extern "C" void applyGamma5(float *out, float *in, int len);


// ---------- gauge_read.cpp ----------

extern "C" void readGaugeField(char *filename, float *gauge[], int argc, char *argv[]);

