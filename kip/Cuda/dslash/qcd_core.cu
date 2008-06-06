#define SHARED_BYTES (BLOCK_DIM*19*sizeof(float))

#define i00_re I0.x
#define i00_im I0.y
#define i01_re I0.z
#define i01_im I0.w
#define i02_re I1.x
#define i02_im I1.y
#define i10_re I1.z
#define i10_im I1.w
#define i11_re I2.x
#define i11_im I2.y
#define i12_re I2.z
#define i12_im I2.w
#define i20_re I3.x
#define i20_im I3.y
#define i21_re I3.z
#define i21_im I3.w
#define i22_re I4.x
#define i22_im I4.y
#define i30_re I4.z
#define i30_im I4.w
#define i31_re I5.x
#define i31_im I5.y
#define i32_re I5.z
#define i32_im I5.w

#define g00_re G0.x
#define g00_im G0.y
#define g01_re G0.z
#define g01_im G0.w
#define g02_re G1.x
#define g02_im G1.y
#define g10_re G1.z
#define g10_im G1.w
#define g11_re G2.x
#define g11_im G2.y
#define g12_re G2.z
#define g12_im G2.w
#define g20_re G3.x
#define g20_im G3.y
#define g21_re G3.z
#define g21_im G3.w
#define g22_re G4.x
#define g22_im G4.y

#define gT00_re (+g00_re)
#define gT00_im (-g00_im)
#define gT01_re (+g10_re)
#define gT01_im (-g10_im)
#define gT02_re (+g20_re)
#define gT02_im (-g20_im)
#define gT10_re (+g01_re)
#define gT10_im (-g01_im)
#define gT11_re (+g11_re)
#define gT11_im (-g11_im)
#define gT12_re (+g21_re)
#define gT12_im (-g21_im)
#define gT20_re (+g02_re)
#define gT20_im (-g02_im)
#define gT21_re (+g12_re)
#define gT21_im (-g12_im)
#define gT22_re (+g22_re)
#define gT22_im (-g22_im)

#define STRIDE 1
#define o00_re s[STRIDE*0]
#define o00_im s[STRIDE*1]
#define o01_re s[STRIDE*2]
#define o01_im s[STRIDE*3]
#define o02_re s[STRIDE*4]
#define o02_im s[STRIDE*5]
#define o10_re s[STRIDE*6]
#define o10_im s[STRIDE*7]
#define o11_re s[STRIDE*8]
#define o11_im s[STRIDE*9]
#define o12_re s[STRIDE*10]
#define o12_im s[STRIDE*11]
#define o20_re s[STRIDE*12]
#define o20_im s[STRIDE*13]
#define o21_re s[STRIDE*14]
#define o21_im s[STRIDE*15]
#define o22_re s[STRIDE*16]
#define o22_im s[STRIDE*17]
#define o30_re s[STRIDE*18]
float o30_im;
float o31_re;
float o31_im;
float o32_re;
float o32_im;


// Performs the complex conjugated accumulation: a += b* c*
#define ACC_CONJ_PROD(a, b, c) \
    a##_re += b##_re * c##_re - b##_im * c##_im, \
    a##_im -= b##_re * c##_im + b##_im * c##_re

#if 1
#define READ_GAUGE_MATRIX(gauge) \
    float4 G0 = tex1Dfetch((gauge), ga_idx + 0*Nh); \
    float4 G1 = tex1Dfetch((gauge), ga_idx + 1*Nh); \
    float4 G2 = tex1Dfetch((gauge), ga_idx + 2*Nh); \
    float4 G3 = tex1Dfetch((gauge), ga_idx + 3*Nh); \
    float4 G4 = tex1Dfetch((gauge), ga_idx + 4*Nh);
#else
#define READ_GAUGE_MATRIX(gauge) \
    float4 G0 = tex1Dfetch((gauge), ga_idx + 0*Nh); \
    float4 G1 = tex1Dfetch((gauge), ga_idx + 1*Nh); \
    float4 G2 = tex1Dfetch((gauge), ga_idx + 2*Nh); \
    float4 G3 = make_float4(0,0,0,0); \
    float4 G4 = make_float4(0,0,0,0); \
    ACC_CONJ_PROD(g20, +g01, +g12); \
    ACC_CONJ_PROD(g20, -g02, +g11); \
    ACC_CONJ_PROD(g21, +g02, +g10); \
    ACC_CONJ_PROD(g21, -g00, +g12); \
    ACC_CONJ_PROD(g22, +g00, +g11); \
    ACC_CONJ_PROD(g22, -g01, +g10);    
#endif

#define READ_SPINOR(spinor) \
    float4 I0 = tex1Dfetch((spinor), sp_idx + 0*Nh); \
    float4 I1 = tex1Dfetch((spinor), sp_idx + 1*Nh); \
    float4 I2 = tex1Dfetch((spinor), sp_idx + 2*Nh); \
    float4 I3 = tex1Dfetch((spinor), sp_idx + 3*Nh); \
    float4 I4 = tex1Dfetch((spinor), sp_idx + 4*Nh); \
    float4 I5 = tex1Dfetch((spinor), sp_idx + 5*Nh);

#ifndef DSLASH_DAGGER
    #define g0(x) g##x
    #define g1(x) gT##x
#else
    #undef g0
    #undef g1
    #define g0(x) gT##x
    #define g1(x) g##x
#endif

int sid = BLOCK_DIM*blockIdx.x + threadIdx.x;
int boundaryCrossings = sid/L1h + sid/(L2*L1h) + sid/(L3*L2*L1h);
int X = 2*sid + (boundaryCrossings + oddBit) % 2;
int x4 = X/(L3*L2*L1);
int x3 = (X/(L2*L1)) % L3;
int x2 = (X/L1) % L2;
int x1 = X % L1;

extern __shared__ float s_data[];
volatile float *s = s_data+19*threadIdx.x;

o00_re = o00_im = 0;
o01_re = o01_im = 0;
o02_re = o02_im = 0;
o10_re = o10_im = 0;
o11_re = o11_im = 0;
o12_re = o12_im = 0;
o20_re = o20_im = 0;
o21_re = o21_im = 0;
o22_re = o22_im = 0;
o30_re = o30_im = 0;
o31_re = o31_im = 0;
o32_re = o32_im = 0;

if(1)
{
    // Projector 0
    // 1 0 0 -i 
    // 0 1 -i 0 
    // 0 i 1 0 
    // i 0 0 1 
    
    int sp_idx = ((x1==L1-1) ? X-(L1-1) : X+1) / 2;
    int ga_idx = sid + (0/2)*Nh*(20/4);
    
    // read spinor from device memory
    READ_SPINOR(spinorTex);
    
    // project spinor into half spinors
    float a0_re = +i00_re+i30_im;
    float a0_im = +i00_im-i30_re;
    float a1_re = +i01_re+i31_im;
    float a1_im = +i01_im-i31_re;
    float a2_re = +i02_re+i32_im;
    float a2_im = +i02_im-i32_re;
    
    float b0_re = +i10_re+i20_im;
    float b0_im = +i10_im-i20_re;
    float b1_re = +i11_re+i21_im;
    float b1_im = +i11_im-i21_re;
    float b2_re = +i12_re+i22_im;
    float b2_im = +i12_im-i22_re;
    
    // read gauge matrix from device memory
    READ_GAUGE_MATRIX(gauge0Tex);
    
    // multiply row 0 by half spinors
    {
        float A_re = + (g0(00_re) * a0_re - g0(00_im) * a0_im) + (g0(01_re) * a1_re - g0(01_im) * a1_im) + (g0(02_re) * a2_re - g0(02_im) * a2_im);
        float A_im = + (g0(00_re) * a0_im + g0(00_im) * a0_re) + (g0(01_re) * a1_im + g0(01_im) * a1_re) + (g0(02_re) * a2_im + g0(02_im) * a2_re);
        float B_re = + (g0(00_re) * b0_re - g0(00_im) * b0_im) + (g0(01_re) * b1_re - g0(01_im) * b1_im) + (g0(02_re) * b2_re - g0(02_im) * b2_im);
        float B_im = + (g0(00_re) * b0_im + g0(00_im) * b0_re) + (g0(01_re) * b1_im + g0(01_im) * b1_re) + (g0(02_re) * b2_im + g0(02_im) * b2_re);
        o00_re += +A_re;
        o00_im += +A_im;
        o10_re += +B_re;
        o10_im += +B_im;
        o20_re += -B_im;
        o20_im += +B_re;
        o30_re += -A_im;
        o30_im += +A_re;
    }
    
    // multiply row 1 by half spinors
    {
        float A_re = + (g0(10_re) * a0_re - g0(10_im) * a0_im) + (g0(11_re) * a1_re - g0(11_im) * a1_im) + (g0(12_re) * a2_re - g0(12_im) * a2_im);
        float A_im = + (g0(10_re) * a0_im + g0(10_im) * a0_re) + (g0(11_re) * a1_im + g0(11_im) * a1_re) + (g0(12_re) * a2_im + g0(12_im) * a2_re);
        float B_re = + (g0(10_re) * b0_re - g0(10_im) * b0_im) + (g0(11_re) * b1_re - g0(11_im) * b1_im) + (g0(12_re) * b2_re - g0(12_im) * b2_im);
        float B_im = + (g0(10_re) * b0_im + g0(10_im) * b0_re) + (g0(11_re) * b1_im + g0(11_im) * b1_re) + (g0(12_re) * b2_im + g0(12_im) * b2_re);
        o01_re += +A_re;
        o01_im += +A_im;
        o11_re += +B_re;
        o11_im += +B_im;
        o21_re += -B_im;
        o21_im += +B_re;
        o31_re += -A_im;
        o31_im += +A_re;
    }
    
    // multiply row 2 by half spinors
    {
        float A_re = + (g0(20_re) * a0_re - g0(20_im) * a0_im) + (g0(21_re) * a1_re - g0(21_im) * a1_im) + (g0(22_re) * a2_re - g0(22_im) * a2_im);
        float A_im = + (g0(20_re) * a0_im + g0(20_im) * a0_re) + (g0(21_re) * a1_im + g0(21_im) * a1_re) + (g0(22_re) * a2_im + g0(22_im) * a2_re);
        float B_re = + (g0(20_re) * b0_re - g0(20_im) * b0_im) + (g0(21_re) * b1_re - g0(21_im) * b1_im) + (g0(22_re) * b2_re - g0(22_im) * b2_im);
        float B_im = + (g0(20_re) * b0_im + g0(20_im) * b0_re) + (g0(21_re) * b1_im + g0(21_im) * b1_re) + (g0(22_re) * b2_im + g0(22_im) * b2_re);
        o02_re += +A_re;
        o02_im += +A_im;
        o12_re += +B_re;
        o12_im += +B_im;
        o22_re += -B_im;
        o22_im += +B_re;
        o32_re += -A_im;
        o32_im += +A_re;
    }
    
}

if(1)
{
    // Projector 1
    // 1 0 0 i 
    // 0 1 i 0 
    // 0 -i 1 0 
    // -i 0 0 1 
    
    int sp_idx = ((x1==0)    ? X+(L1-1) : X-1) / 2;
    int ga_idx = sp_idx + (1/2)*Nh*(20/4);
    
    // read spinor from device memory
    READ_SPINOR(spinorTex);
    
    // project spinor into half spinors
    float a0_re = +i00_re-i30_im;
    float a0_im = +i00_im+i30_re;
    float a1_re = +i01_re-i31_im;
    float a1_im = +i01_im+i31_re;
    float a2_re = +i02_re-i32_im;
    float a2_im = +i02_im+i32_re;
    
    float b0_re = +i10_re-i20_im;
    float b0_im = +i10_im+i20_re;
    float b1_re = +i11_re-i21_im;
    float b1_im = +i11_im+i21_re;
    float b2_re = +i12_re-i22_im;
    float b2_im = +i12_im+i22_re;
    
    // read gauge matrix from device memory
    READ_GAUGE_MATRIX(gauge1Tex);
    
    // multiply row 0 by half spinors
    {
        float A_re = + (g1(00_re) * a0_re - g1(00_im) * a0_im) + (g1(01_re) * a1_re - g1(01_im) * a1_im) + (g1(02_re) * a2_re - g1(02_im) * a2_im);
        float A_im = + (g1(00_re) * a0_im + g1(00_im) * a0_re) + (g1(01_re) * a1_im + g1(01_im) * a1_re) + (g1(02_re) * a2_im + g1(02_im) * a2_re);
        float B_re = + (g1(00_re) * b0_re - g1(00_im) * b0_im) + (g1(01_re) * b1_re - g1(01_im) * b1_im) + (g1(02_re) * b2_re - g1(02_im) * b2_im);
        float B_im = + (g1(00_re) * b0_im + g1(00_im) * b0_re) + (g1(01_re) * b1_im + g1(01_im) * b1_re) + (g1(02_re) * b2_im + g1(02_im) * b2_re);
        o00_re += +A_re;
        o00_im += +A_im;
        o10_re += +B_re;
        o10_im += +B_im;
        o20_re += +B_im;
        o20_im += -B_re;
        o30_re += +A_im;
        o30_im += -A_re;
    }
    
    // multiply row 1 by half spinors
    {
        float A_re = + (g1(10_re) * a0_re - g1(10_im) * a0_im) + (g1(11_re) * a1_re - g1(11_im) * a1_im) + (g1(12_re) * a2_re - g1(12_im) * a2_im);
        float A_im = + (g1(10_re) * a0_im + g1(10_im) * a0_re) + (g1(11_re) * a1_im + g1(11_im) * a1_re) + (g1(12_re) * a2_im + g1(12_im) * a2_re);
        float B_re = + (g1(10_re) * b0_re - g1(10_im) * b0_im) + (g1(11_re) * b1_re - g1(11_im) * b1_im) + (g1(12_re) * b2_re - g1(12_im) * b2_im);
        float B_im = + (g1(10_re) * b0_im + g1(10_im) * b0_re) + (g1(11_re) * b1_im + g1(11_im) * b1_re) + (g1(12_re) * b2_im + g1(12_im) * b2_re);
        o01_re += +A_re;
        o01_im += +A_im;
        o11_re += +B_re;
        o11_im += +B_im;
        o21_re += +B_im;
        o21_im += -B_re;
        o31_re += +A_im;
        o31_im += -A_re;
    }
    
    // multiply row 2 by half spinors
    {
        float A_re = + (g1(20_re) * a0_re - g1(20_im) * a0_im) + (g1(21_re) * a1_re - g1(21_im) * a1_im) + (g1(22_re) * a2_re - g1(22_im) * a2_im);
        float A_im = + (g1(20_re) * a0_im + g1(20_im) * a0_re) + (g1(21_re) * a1_im + g1(21_im) * a1_re) + (g1(22_re) * a2_im + g1(22_im) * a2_re);
        float B_re = + (g1(20_re) * b0_re - g1(20_im) * b0_im) + (g1(21_re) * b1_re - g1(21_im) * b1_im) + (g1(22_re) * b2_re - g1(22_im) * b2_im);
        float B_im = + (g1(20_re) * b0_im + g1(20_im) * b0_re) + (g1(21_re) * b1_im + g1(21_im) * b1_re) + (g1(22_re) * b2_im + g1(22_im) * b2_re);
        o02_re += +A_re;
        o02_im += +A_im;
        o12_re += +B_re;
        o12_im += +B_im;
        o22_re += +B_im;
        o22_im += -B_re;
        o32_re += +A_im;
        o32_im += -A_re;
    }
    
}

if(1)
{
    // Projector 2
    // 1 0 0 1 
    // 0 1 -1 0 
    // 0 -1 1 0 
    // 1 0 0 1 
    
    int sp_idx = ((x2==L2-1) ? X-(L2-1)*L1 : X+L1) / 2;
    int ga_idx = sid + (2/2)*Nh*(20/4);
    
    // read spinor from device memory
    READ_SPINOR(spinorTex);
    
    // project spinor into half spinors
    float a0_re = +i00_re+i30_re;
    float a0_im = +i00_im+i30_im;
    float a1_re = +i01_re+i31_re;
    float a1_im = +i01_im+i31_im;
    float a2_re = +i02_re+i32_re;
    float a2_im = +i02_im+i32_im;
    
    float b0_re = +i10_re-i20_re;
    float b0_im = +i10_im-i20_im;
    float b1_re = +i11_re-i21_re;
    float b1_im = +i11_im-i21_im;
    float b2_re = +i12_re-i22_re;
    float b2_im = +i12_im-i22_im;
    
    // read gauge matrix from device memory
    READ_GAUGE_MATRIX(gauge0Tex);
    
    // multiply row 0 by half spinors
    {
        float A_re = + (g0(00_re) * a0_re - g0(00_im) * a0_im) + (g0(01_re) * a1_re - g0(01_im) * a1_im) + (g0(02_re) * a2_re - g0(02_im) * a2_im);
        float A_im = + (g0(00_re) * a0_im + g0(00_im) * a0_re) + (g0(01_re) * a1_im + g0(01_im) * a1_re) + (g0(02_re) * a2_im + g0(02_im) * a2_re);
        float B_re = + (g0(00_re) * b0_re - g0(00_im) * b0_im) + (g0(01_re) * b1_re - g0(01_im) * b1_im) + (g0(02_re) * b2_re - g0(02_im) * b2_im);
        float B_im = + (g0(00_re) * b0_im + g0(00_im) * b0_re) + (g0(01_re) * b1_im + g0(01_im) * b1_re) + (g0(02_re) * b2_im + g0(02_im) * b2_re);
        o00_re += +A_re;
        o00_im += +A_im;
        o10_re += +B_re;
        o10_im += +B_im;
        o20_re += -B_re;
        o20_im += -B_im;
        o30_re += +A_re;
        o30_im += +A_im;
    }
    
    // multiply row 1 by half spinors
    {
        float A_re = + (g0(10_re) * a0_re - g0(10_im) * a0_im) + (g0(11_re) * a1_re - g0(11_im) * a1_im) + (g0(12_re) * a2_re - g0(12_im) * a2_im);
        float A_im = + (g0(10_re) * a0_im + g0(10_im) * a0_re) + (g0(11_re) * a1_im + g0(11_im) * a1_re) + (g0(12_re) * a2_im + g0(12_im) * a2_re);
        float B_re = + (g0(10_re) * b0_re - g0(10_im) * b0_im) + (g0(11_re) * b1_re - g0(11_im) * b1_im) + (g0(12_re) * b2_re - g0(12_im) * b2_im);
        float B_im = + (g0(10_re) * b0_im + g0(10_im) * b0_re) + (g0(11_re) * b1_im + g0(11_im) * b1_re) + (g0(12_re) * b2_im + g0(12_im) * b2_re);
        o01_re += +A_re;
        o01_im += +A_im;
        o11_re += +B_re;
        o11_im += +B_im;
        o21_re += -B_re;
        o21_im += -B_im;
        o31_re += +A_re;
        o31_im += +A_im;
    }
    
    // multiply row 2 by half spinors
    {
        float A_re = + (g0(20_re) * a0_re - g0(20_im) * a0_im) + (g0(21_re) * a1_re - g0(21_im) * a1_im) + (g0(22_re) * a2_re - g0(22_im) * a2_im);
        float A_im = + (g0(20_re) * a0_im + g0(20_im) * a0_re) + (g0(21_re) * a1_im + g0(21_im) * a1_re) + (g0(22_re) * a2_im + g0(22_im) * a2_re);
        float B_re = + (g0(20_re) * b0_re - g0(20_im) * b0_im) + (g0(21_re) * b1_re - g0(21_im) * b1_im) + (g0(22_re) * b2_re - g0(22_im) * b2_im);
        float B_im = + (g0(20_re) * b0_im + g0(20_im) * b0_re) + (g0(21_re) * b1_im + g0(21_im) * b1_re) + (g0(22_re) * b2_im + g0(22_im) * b2_re);
        o02_re += +A_re;
        o02_im += +A_im;
        o12_re += +B_re;
        o12_im += +B_im;
        o22_re += -B_re;
        o22_im += -B_im;
        o32_re += +A_re;
        o32_im += +A_im;
    }
    
}

if(1)
{
    // Projector 3
    // 1 0 0 -1 
    // 0 1 1 0 
    // 0 1 1 0 
    // -1 0 0 1 
    
    int sp_idx = ((x2==0)    ? X+(L2-1)*L1 : X-L1) / 2;
    int ga_idx = sp_idx + (3/2)*Nh*(20/4);
    
    // read spinor from device memory
    READ_SPINOR(spinorTex);
    
    // project spinor into half spinors
    float a0_re = +i00_re-i30_re;
    float a0_im = +i00_im-i30_im;
    float a1_re = +i01_re-i31_re;
    float a1_im = +i01_im-i31_im;
    float a2_re = +i02_re-i32_re;
    float a2_im = +i02_im-i32_im;
    
    float b0_re = +i10_re+i20_re;
    float b0_im = +i10_im+i20_im;
    float b1_re = +i11_re+i21_re;
    float b1_im = +i11_im+i21_im;
    float b2_re = +i12_re+i22_re;
    float b2_im = +i12_im+i22_im;
    
    // read gauge matrix from device memory
    READ_GAUGE_MATRIX(gauge1Tex);
    
    // multiply row 0 by half spinors
    {
        float A_re = + (g1(00_re) * a0_re - g1(00_im) * a0_im) + (g1(01_re) * a1_re - g1(01_im) * a1_im) + (g1(02_re) * a2_re - g1(02_im) * a2_im);
        float A_im = + (g1(00_re) * a0_im + g1(00_im) * a0_re) + (g1(01_re) * a1_im + g1(01_im) * a1_re) + (g1(02_re) * a2_im + g1(02_im) * a2_re);
        float B_re = + (g1(00_re) * b0_re - g1(00_im) * b0_im) + (g1(01_re) * b1_re - g1(01_im) * b1_im) + (g1(02_re) * b2_re - g1(02_im) * b2_im);
        float B_im = + (g1(00_re) * b0_im + g1(00_im) * b0_re) + (g1(01_re) * b1_im + g1(01_im) * b1_re) + (g1(02_re) * b2_im + g1(02_im) * b2_re);
        o00_re += +A_re;
        o00_im += +A_im;
        o10_re += +B_re;
        o10_im += +B_im;
        o20_re += +B_re;
        o20_im += +B_im;
        o30_re += -A_re;
        o30_im += -A_im;
    }
    
    // multiply row 1 by half spinors
    {
        float A_re = + (g1(10_re) * a0_re - g1(10_im) * a0_im) + (g1(11_re) * a1_re - g1(11_im) * a1_im) + (g1(12_re) * a2_re - g1(12_im) * a2_im);
        float A_im = + (g1(10_re) * a0_im + g1(10_im) * a0_re) + (g1(11_re) * a1_im + g1(11_im) * a1_re) + (g1(12_re) * a2_im + g1(12_im) * a2_re);
        float B_re = + (g1(10_re) * b0_re - g1(10_im) * b0_im) + (g1(11_re) * b1_re - g1(11_im) * b1_im) + (g1(12_re) * b2_re - g1(12_im) * b2_im);
        float B_im = + (g1(10_re) * b0_im + g1(10_im) * b0_re) + (g1(11_re) * b1_im + g1(11_im) * b1_re) + (g1(12_re) * b2_im + g1(12_im) * b2_re);
        o01_re += +A_re;
        o01_im += +A_im;
        o11_re += +B_re;
        o11_im += +B_im;
        o21_re += +B_re;
        o21_im += +B_im;
        o31_re += -A_re;
        o31_im += -A_im;
    }
    
    // multiply row 2 by half spinors
    {
        float A_re = + (g1(20_re) * a0_re - g1(20_im) * a0_im) + (g1(21_re) * a1_re - g1(21_im) * a1_im) + (g1(22_re) * a2_re - g1(22_im) * a2_im);
        float A_im = + (g1(20_re) * a0_im + g1(20_im) * a0_re) + (g1(21_re) * a1_im + g1(21_im) * a1_re) + (g1(22_re) * a2_im + g1(22_im) * a2_re);
        float B_re = + (g1(20_re) * b0_re - g1(20_im) * b0_im) + (g1(21_re) * b1_re - g1(21_im) * b1_im) + (g1(22_re) * b2_re - g1(22_im) * b2_im);
        float B_im = + (g1(20_re) * b0_im + g1(20_im) * b0_re) + (g1(21_re) * b1_im + g1(21_im) * b1_re) + (g1(22_re) * b2_im + g1(22_im) * b2_re);
        o02_re += +A_re;
        o02_im += +A_im;
        o12_re += +B_re;
        o12_im += +B_im;
        o22_re += +B_re;
        o22_im += +B_im;
        o32_re += -A_re;
        o32_im += -A_im;
    }
    
}

if(1)
{
    // Projector 4
    // 1 0 -i 0 
    // 0 1 0 i 
    // i 0 1 0 
    // 0 -i 0 1 
    
    int sp_idx = ((x3==L3-1) ? X-(L3-1)*L2*L1 : X+L2*L1) / 2;
    int ga_idx = sid + (4/2)*Nh*(20/4);
    
    // read spinor from device memory
    READ_SPINOR(spinorTex);
    
    // project spinor into half spinors
    float a0_re = +i00_re+i20_im;
    float a0_im = +i00_im-i20_re;
    float a1_re = +i01_re+i21_im;
    float a1_im = +i01_im-i21_re;
    float a2_re = +i02_re+i22_im;
    float a2_im = +i02_im-i22_re;
    
    float b0_re = +i10_re-i30_im;
    float b0_im = +i10_im+i30_re;
    float b1_re = +i11_re-i31_im;
    float b1_im = +i11_im+i31_re;
    float b2_re = +i12_re-i32_im;
    float b2_im = +i12_im+i32_re;
    
    // read gauge matrix from device memory
    READ_GAUGE_MATRIX(gauge0Tex);
    
    // multiply row 0 by half spinors
    {
        float A_re = + (g0(00_re) * a0_re - g0(00_im) * a0_im) + (g0(01_re) * a1_re - g0(01_im) * a1_im) + (g0(02_re) * a2_re - g0(02_im) * a2_im);
        float A_im = + (g0(00_re) * a0_im + g0(00_im) * a0_re) + (g0(01_re) * a1_im + g0(01_im) * a1_re) + (g0(02_re) * a2_im + g0(02_im) * a2_re);
        float B_re = + (g0(00_re) * b0_re - g0(00_im) * b0_im) + (g0(01_re) * b1_re - g0(01_im) * b1_im) + (g0(02_re) * b2_re - g0(02_im) * b2_im);
        float B_im = + (g0(00_re) * b0_im + g0(00_im) * b0_re) + (g0(01_re) * b1_im + g0(01_im) * b1_re) + (g0(02_re) * b2_im + g0(02_im) * b2_re);
        o00_re += +A_re;
        o00_im += +A_im;
        o10_re += +B_re;
        o10_im += +B_im;
        o20_re += -A_im;
        o20_im += +A_re;
        o30_re += +B_im;
        o30_im += -B_re;
    }
    
    // multiply row 1 by half spinors
    {
        float A_re = + (g0(10_re) * a0_re - g0(10_im) * a0_im) + (g0(11_re) * a1_re - g0(11_im) * a1_im) + (g0(12_re) * a2_re - g0(12_im) * a2_im);
        float A_im = + (g0(10_re) * a0_im + g0(10_im) * a0_re) + (g0(11_re) * a1_im + g0(11_im) * a1_re) + (g0(12_re) * a2_im + g0(12_im) * a2_re);
        float B_re = + (g0(10_re) * b0_re - g0(10_im) * b0_im) + (g0(11_re) * b1_re - g0(11_im) * b1_im) + (g0(12_re) * b2_re - g0(12_im) * b2_im);
        float B_im = + (g0(10_re) * b0_im + g0(10_im) * b0_re) + (g0(11_re) * b1_im + g0(11_im) * b1_re) + (g0(12_re) * b2_im + g0(12_im) * b2_re);
        o01_re += +A_re;
        o01_im += +A_im;
        o11_re += +B_re;
        o11_im += +B_im;
        o21_re += -A_im;
        o21_im += +A_re;
        o31_re += +B_im;
        o31_im += -B_re;
    }
    
    // multiply row 2 by half spinors
    {
        float A_re = + (g0(20_re) * a0_re - g0(20_im) * a0_im) + (g0(21_re) * a1_re - g0(21_im) * a1_im) + (g0(22_re) * a2_re - g0(22_im) * a2_im);
        float A_im = + (g0(20_re) * a0_im + g0(20_im) * a0_re) + (g0(21_re) * a1_im + g0(21_im) * a1_re) + (g0(22_re) * a2_im + g0(22_im) * a2_re);
        float B_re = + (g0(20_re) * b0_re - g0(20_im) * b0_im) + (g0(21_re) * b1_re - g0(21_im) * b1_im) + (g0(22_re) * b2_re - g0(22_im) * b2_im);
        float B_im = + (g0(20_re) * b0_im + g0(20_im) * b0_re) + (g0(21_re) * b1_im + g0(21_im) * b1_re) + (g0(22_re) * b2_im + g0(22_im) * b2_re);
        o02_re += +A_re;
        o02_im += +A_im;
        o12_re += +B_re;
        o12_im += +B_im;
        o22_re += -A_im;
        o22_im += +A_re;
        o32_re += +B_im;
        o32_im += -B_re;
    }
    
}

if(1)
{
    // Projector 5
    // 1 0 i 0 
    // 0 1 0 -i 
    // -i 0 1 0 
    // 0 i 0 1 
    
    int sp_idx = ((x3==0)    ? X+(L3-1)*L2*L1 : X-L2*L1) / 2;
    int ga_idx = sp_idx + (5/2)*Nh*(20/4);
    
    // read spinor from device memory
    READ_SPINOR(spinorTex);
    
    // project spinor into half spinors
    float a0_re = +i00_re-i20_im;
    float a0_im = +i00_im+i20_re;
    float a1_re = +i01_re-i21_im;
    float a1_im = +i01_im+i21_re;
    float a2_re = +i02_re-i22_im;
    float a2_im = +i02_im+i22_re;
    
    float b0_re = +i10_re+i30_im;
    float b0_im = +i10_im-i30_re;
    float b1_re = +i11_re+i31_im;
    float b1_im = +i11_im-i31_re;
    float b2_re = +i12_re+i32_im;
    float b2_im = +i12_im-i32_re;
    
    // read gauge matrix from device memory
    READ_GAUGE_MATRIX(gauge1Tex);
    
    // multiply row 0 by half spinors
    {
        float A_re = + (g1(00_re) * a0_re - g1(00_im) * a0_im) + (g1(01_re) * a1_re - g1(01_im) * a1_im) + (g1(02_re) * a2_re - g1(02_im) * a2_im);
        float A_im = + (g1(00_re) * a0_im + g1(00_im) * a0_re) + (g1(01_re) * a1_im + g1(01_im) * a1_re) + (g1(02_re) * a2_im + g1(02_im) * a2_re);
        float B_re = + (g1(00_re) * b0_re - g1(00_im) * b0_im) + (g1(01_re) * b1_re - g1(01_im) * b1_im) + (g1(02_re) * b2_re - g1(02_im) * b2_im);
        float B_im = + (g1(00_re) * b0_im + g1(00_im) * b0_re) + (g1(01_re) * b1_im + g1(01_im) * b1_re) + (g1(02_re) * b2_im + g1(02_im) * b2_re);
        o00_re += +A_re;
        o00_im += +A_im;
        o10_re += +B_re;
        o10_im += +B_im;
        o20_re += +A_im;
        o20_im += -A_re;
        o30_re += -B_im;
        o30_im += +B_re;
    }
    
    // multiply row 1 by half spinors
    {
        float A_re = + (g1(10_re) * a0_re - g1(10_im) * a0_im) + (g1(11_re) * a1_re - g1(11_im) * a1_im) + (g1(12_re) * a2_re - g1(12_im) * a2_im);
        float A_im = + (g1(10_re) * a0_im + g1(10_im) * a0_re) + (g1(11_re) * a1_im + g1(11_im) * a1_re) + (g1(12_re) * a2_im + g1(12_im) * a2_re);
        float B_re = + (g1(10_re) * b0_re - g1(10_im) * b0_im) + (g1(11_re) * b1_re - g1(11_im) * b1_im) + (g1(12_re) * b2_re - g1(12_im) * b2_im);
        float B_im = + (g1(10_re) * b0_im + g1(10_im) * b0_re) + (g1(11_re) * b1_im + g1(11_im) * b1_re) + (g1(12_re) * b2_im + g1(12_im) * b2_re);
        o01_re += +A_re;
        o01_im += +A_im;
        o11_re += +B_re;
        o11_im += +B_im;
        o21_re += +A_im;
        o21_im += -A_re;
        o31_re += -B_im;
        o31_im += +B_re;
    }
    
    // multiply row 2 by half spinors
    {
        float A_re = + (g1(20_re) * a0_re - g1(20_im) * a0_im) + (g1(21_re) * a1_re - g1(21_im) * a1_im) + (g1(22_re) * a2_re - g1(22_im) * a2_im);
        float A_im = + (g1(20_re) * a0_im + g1(20_im) * a0_re) + (g1(21_re) * a1_im + g1(21_im) * a1_re) + (g1(22_re) * a2_im + g1(22_im) * a2_re);
        float B_re = + (g1(20_re) * b0_re - g1(20_im) * b0_im) + (g1(21_re) * b1_re - g1(21_im) * b1_im) + (g1(22_re) * b2_re - g1(22_im) * b2_im);
        float B_im = + (g1(20_re) * b0_im + g1(20_im) * b0_re) + (g1(21_re) * b1_im + g1(21_im) * b1_re) + (g1(22_re) * b2_im + g1(22_im) * b2_re);
        o02_re += +A_re;
        o02_im += +A_im;
        o12_re += +B_re;
        o12_im += +B_im;
        o22_re += +A_im;
        o22_im += -A_re;
        o32_re += -B_im;
        o32_im += +B_re;
    }
    
}

if(1)
{
    // Projector 6
    // 1 0 -1 0 
    // 0 1 0 -1 
    // -1 0 1 0 
    // 0 -1 0 1 
    
    int sp_idx = ((x4==L4-1) ? X-(L4-1)*L3*L2*L1 : X+L3*L2*L1) / 2;
    int ga_idx = sid + (6/2)*Nh*(20/4);
    
    // read spinor from device memory
    READ_SPINOR(spinorTex);
    
    // project spinor into half spinors
    float a0_re = +i00_re-i20_re;
    float a0_im = +i00_im-i20_im;
    float a1_re = +i01_re-i21_re;
    float a1_im = +i01_im-i21_im;
    float a2_re = +i02_re-i22_re;
    float a2_im = +i02_im-i22_im;
    
    float b0_re = +i10_re-i30_re;
    float b0_im = +i10_im-i30_im;
    float b1_re = +i11_re-i31_re;
    float b1_im = +i11_im-i31_im;
    float b2_re = +i12_re-i32_re;
    float b2_im = +i12_im-i32_im;
    
    // read gauge matrix from device memory
    READ_GAUGE_MATRIX(gauge0Tex);
    
    // multiply row 0 by half spinors
    {
        float A_re = + (g0(00_re) * a0_re - g0(00_im) * a0_im) + (g0(01_re) * a1_re - g0(01_im) * a1_im) + (g0(02_re) * a2_re - g0(02_im) * a2_im);
        float A_im = + (g0(00_re) * a0_im + g0(00_im) * a0_re) + (g0(01_re) * a1_im + g0(01_im) * a1_re) + (g0(02_re) * a2_im + g0(02_im) * a2_re);
        float B_re = + (g0(00_re) * b0_re - g0(00_im) * b0_im) + (g0(01_re) * b1_re - g0(01_im) * b1_im) + (g0(02_re) * b2_re - g0(02_im) * b2_im);
        float B_im = + (g0(00_re) * b0_im + g0(00_im) * b0_re) + (g0(01_re) * b1_im + g0(01_im) * b1_re) + (g0(02_re) * b2_im + g0(02_im) * b2_re);
        o00_re += +A_re;
        o00_im += +A_im;
        o10_re += +B_re;
        o10_im += +B_im;
        o20_re += -A_re;
        o20_im += -A_im;
        o30_re += -B_re;
        o30_im += -B_im;
    }
    
    // multiply row 1 by half spinors
    {
        float A_re = + (g0(10_re) * a0_re - g0(10_im) * a0_im) + (g0(11_re) * a1_re - g0(11_im) * a1_im) + (g0(12_re) * a2_re - g0(12_im) * a2_im);
        float A_im = + (g0(10_re) * a0_im + g0(10_im) * a0_re) + (g0(11_re) * a1_im + g0(11_im) * a1_re) + (g0(12_re) * a2_im + g0(12_im) * a2_re);
        float B_re = + (g0(10_re) * b0_re - g0(10_im) * b0_im) + (g0(11_re) * b1_re - g0(11_im) * b1_im) + (g0(12_re) * b2_re - g0(12_im) * b2_im);
        float B_im = + (g0(10_re) * b0_im + g0(10_im) * b0_re) + (g0(11_re) * b1_im + g0(11_im) * b1_re) + (g0(12_re) * b2_im + g0(12_im) * b2_re);
        o01_re += +A_re;
        o01_im += +A_im;
        o11_re += +B_re;
        o11_im += +B_im;
        o21_re += -A_re;
        o21_im += -A_im;
        o31_re += -B_re;
        o31_im += -B_im;
    }
    
    // multiply row 2 by half spinors
    {
        float A_re = + (g0(20_re) * a0_re - g0(20_im) * a0_im) + (g0(21_re) * a1_re - g0(21_im) * a1_im) + (g0(22_re) * a2_re - g0(22_im) * a2_im);
        float A_im = + (g0(20_re) * a0_im + g0(20_im) * a0_re) + (g0(21_re) * a1_im + g0(21_im) * a1_re) + (g0(22_re) * a2_im + g0(22_im) * a2_re);
        float B_re = + (g0(20_re) * b0_re - g0(20_im) * b0_im) + (g0(21_re) * b1_re - g0(21_im) * b1_im) + (g0(22_re) * b2_re - g0(22_im) * b2_im);
        float B_im = + (g0(20_re) * b0_im + g0(20_im) * b0_re) + (g0(21_re) * b1_im + g0(21_im) * b1_re) + (g0(22_re) * b2_im + g0(22_im) * b2_re);
        o02_re += +A_re;
        o02_im += +A_im;
        o12_re += +B_re;
        o12_im += +B_im;
        o22_re += -A_re;
        o22_im += -A_im;
        o32_re += -B_re;
        o32_im += -B_im;
    }
    
}

if(1)
{
    // Projector 7
    // 1 0 1 0 
    // 0 1 0 1 
    // 1 0 1 0 
    // 0 1 0 1 
    
    int sp_idx = ((x4==0)    ? X+(L4-1)*L3*L2*L1 : X-L3*L2*L1) / 2;
    int ga_idx = sp_idx + (7/2)*Nh*(20/4);
    
    // read spinor from device memory
    READ_SPINOR(spinorTex);
    
    // project spinor into half spinors
    float a0_re = +i00_re+i20_re;
    float a0_im = +i00_im+i20_im;
    float a1_re = +i01_re+i21_re;
    float a1_im = +i01_im+i21_im;
    float a2_re = +i02_re+i22_re;
    float a2_im = +i02_im+i22_im;
    
    float b0_re = +i10_re+i30_re;
    float b0_im = +i10_im+i30_im;
    float b1_re = +i11_re+i31_re;
    float b1_im = +i11_im+i31_im;
    float b2_re = +i12_re+i32_re;
    float b2_im = +i12_im+i32_im;
    
    // read gauge matrix from device memory
    READ_GAUGE_MATRIX(gauge1Tex);
    
    // multiply row 0 by half spinors
    {
        float A_re = + (g1(00_re) * a0_re - g1(00_im) * a0_im) + (g1(01_re) * a1_re - g1(01_im) * a1_im) + (g1(02_re) * a2_re - g1(02_im) * a2_im);
        float A_im = + (g1(00_re) * a0_im + g1(00_im) * a0_re) + (g1(01_re) * a1_im + g1(01_im) * a1_re) + (g1(02_re) * a2_im + g1(02_im) * a2_re);
        float B_re = + (g1(00_re) * b0_re - g1(00_im) * b0_im) + (g1(01_re) * b1_re - g1(01_im) * b1_im) + (g1(02_re) * b2_re - g1(02_im) * b2_im);
        float B_im = + (g1(00_re) * b0_im + g1(00_im) * b0_re) + (g1(01_re) * b1_im + g1(01_im) * b1_re) + (g1(02_re) * b2_im + g1(02_im) * b2_re);
        o00_re += +A_re;
        o00_im += +A_im;
        o10_re += +B_re;
        o10_im += +B_im;
        o20_re += +A_re;
        o20_im += +A_im;
        o30_re += +B_re;
        o30_im += +B_im;
    }
    
    // multiply row 1 by half spinors
    {
        float A_re = + (g1(10_re) * a0_re - g1(10_im) * a0_im) + (g1(11_re) * a1_re - g1(11_im) * a1_im) + (g1(12_re) * a2_re - g1(12_im) * a2_im);
        float A_im = + (g1(10_re) * a0_im + g1(10_im) * a0_re) + (g1(11_re) * a1_im + g1(11_im) * a1_re) + (g1(12_re) * a2_im + g1(12_im) * a2_re);
        float B_re = + (g1(10_re) * b0_re - g1(10_im) * b0_im) + (g1(11_re) * b1_re - g1(11_im) * b1_im) + (g1(12_re) * b2_re - g1(12_im) * b2_im);
        float B_im = + (g1(10_re) * b0_im + g1(10_im) * b0_re) + (g1(11_re) * b1_im + g1(11_im) * b1_re) + (g1(12_re) * b2_im + g1(12_im) * b2_re);
        o01_re += +A_re;
        o01_im += +A_im;
        o11_re += +B_re;
        o11_im += +B_im;
        o21_re += +A_re;
        o21_im += +A_im;
        o31_re += +B_re;
        o31_im += +B_im;
    }
    
    // multiply row 2 by half spinors
    {
        float A_re = + (g1(20_re) * a0_re - g1(20_im) * a0_im) + (g1(21_re) * a1_re - g1(21_im) * a1_im) + (g1(22_re) * a2_re - g1(22_im) * a2_im);
        float A_im = + (g1(20_re) * a0_im + g1(20_im) * a0_re) + (g1(21_re) * a1_im + g1(21_im) * a1_re) + (g1(22_re) * a2_im + g1(22_im) * a2_re);
        float B_re = + (g1(20_re) * b0_re - g1(20_im) * b0_im) + (g1(21_re) * b1_re - g1(21_im) * b1_im) + (g1(22_re) * b2_re - g1(22_im) * b2_im);
        float B_im = + (g1(20_re) * b0_im + g1(20_im) * b0_re) + (g1(21_re) * b1_im + g1(21_im) * b1_re) + (g1(22_re) * b2_im + g1(22_im) * b2_re);
        o02_re += +A_re;
        o02_im += +A_im;
        o12_re += +B_re;
        o12_im += +B_im;
        o22_re += +A_re;
        o22_im += +A_im;
        o32_re += +B_re;
        o32_im += +B_im;
    }
    
}


//g_out[0*Nh+sid] = make_float4(o00_re, o00_im, o01_re, o01_im);
//g_out[1*Nh+sid] = make_float4(o02_re, o02_im, o10_re, o10_im);
//g_out[2*Nh+sid] = make_float4(o11_re, o11_im, o12_re, o12_im);
//g_out[3*Nh+sid] = make_float4(o20_re, o20_im, o21_re, o21_im);
//g_out[4*Nh+sid] = make_float4(o22_re, o22_im, o30_re, o30_im);
//g_out[5*Nh+sid] = make_float4(o31_re, o31_im, o32_re, o32_im);

((float*)g_out)[0*Nh+sid] = o00_re;
((float*)g_out)[1*Nh+sid] = o00_im;
((float*)g_out)[2*Nh+sid] = o01_re;
((float*)g_out)[3*Nh+sid] = o01_im;
((float*)g_out)[4*Nh+sid] = o02_re;
((float*)g_out)[5*Nh+sid] = o02_im;
((float*)g_out)[6*Nh+sid] = o10_re;
((float*)g_out)[7*Nh+sid] = o10_im;
((float*)g_out)[8*Nh+sid] = o11_re;
((float*)g_out)[9*Nh+sid] = o11_im;
((float*)g_out)[10*Nh+sid] = o12_re;
((float*)g_out)[11*Nh+sid] = o12_im;
((float*)g_out)[12*Nh+sid] = o20_re;
((float*)g_out)[13*Nh+sid] = o20_im;
((float*)g_out)[14*Nh+sid] = o21_re;
((float*)g_out)[15*Nh+sid] = o21_im;
((float*)g_out)[16*Nh+sid] = o22_re;
((float*)g_out)[17*Nh+sid] = o22_im;
((float*)g_out)[18*Nh+sid] = o30_re;
((float*)g_out)[19*Nh+sid] = o30_im;
((float*)g_out)[20*Nh+sid] = o31_re;
((float*)g_out)[21*Nh+sid] = o31_im;
((float*)g_out)[22*Nh+sid] = o32_re;
((float*)g_out)[23*Nh+sid] = o32_im;

