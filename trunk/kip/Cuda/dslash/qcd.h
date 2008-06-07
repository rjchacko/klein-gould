
#define L1 16 // "x" dimension
#define L2 16 // "y" dimension
#define L3 16 // "z" dimension
#define L4 64 // "time" dimension
#define L1h (L1/2) // half of the full "x" dimension, useful for even/odd lattice indexing

#define N (L1*L2*L3*L4) // total number of lattice points
#define Nh (N/2) // total number of even/odd lattice points

#define gaugeSiteSize 18 // real numbers per link
#define spinorSiteSize 24 // real numbers per spinor

#define EVEN_ODD


// ---------- dslash_test.cpp ----------

extern "C" void dslashTest(int argc, char **argv);


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


extern "C" CudaFullGauge loadGaugeField(float **gauge);
extern "C" CudaFullSpinor loadSpinorField(float *spinor);

extern "C" void freeGaugeField(CudaFullGauge gauge);
extern "C" void freeSpinorField(CudaFullSpinor spinor);

extern "C" void retrieveParitySpinor(float *res, CudaPSpinor spinor);
extern "C" void retrieveSpinorField(float *res, CudaFullSpinor spinor);

extern "C" void compareParitySpinors(float *sp1, CudaPSpinor sp2);

extern "C" void dslashCuda(CudaPSpinor res, CudaFullGauge gauge, CudaPSpinor spinor, int oddBit, int daggerBit);
extern "C" void printCudaDslashInfo();


// ---------- dslash_reference.cpp ----------

extern "C" void constructUnitGaugeField(float **gauge);
extern "C" void constructGaugeField(float **gauge);
extern "C" void constructPointSpinorField(float *spinor, int i0, int s0, int c0);
extern "C" void constructSpinorField(float *res);

extern "C" void dslashReference(float *res, float **gauge, float *spinorField, 
				int oddBit, int daggerBit);

extern "C" void Mat(float *out, float **gauge, float *in, float kappa);
extern "C" void MatDag(float *out, float **gauge, float *in, float kappa);
extern "C" void MatDagMat(float *out, float **gauge, float *in, float kappa);

extern "C" void MatPC(float *out, float **gauge, float *in, float kappa);
extern "C" void MatPCDag(float *out, float **gauge, float *in, float kappa);
extern "C" void MatPCDagMatPC(float *out, float **gauge, float *in, float kappa);

extern "C" void printSpinor(float *spinor);
extern "C" void printSpinorElement(float *spinor, int X);

extern "C" void applyGamma5(float *out, float *in, int len);

extern "C" void zero(float* a, int cnt);
extern "C" void ax(float a, float *x, int len);
extern "C" void axpy(float a, float *x, float *y, int len);
extern "C" void xpay(float *x, float a, float *y, int len);
extern "C" void mxpy(float *x, float *y, int len);


// ---------- cg_test.cpp ----------

extern "C" void cgReferenceTest();


// ---------- cg_reference.cpp ----------

extern "C" float norm(float *vector, int len);
extern "C" float reDotProduct(float *v1, float *v2, int len);
extern "C" float imDotProduct(float *v1, float *v2, int len);
extern "C" double normD(float *vector, int len);
extern "C" double reDotProductD(float *v1, float *v2, int len);
extern "C" double imDotProductD(float *v1, float *v2, int len);
extern "C" float copy(float *vector, int len);

extern "C" void cg_reference(float *out, float **gauge, float *in, float kappa, float tol);


extern "C" void cgReferenceTest();

