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

const int sid = BLOCK_DIM*blockIdx.x + threadIdx.x;

extern __shared__ float s_data[];
volatile float *s = s_data+19*threadIdx.x;

int x1 = sid % L1;
int x2 = (sid/L1) % L2;
int x3 = (sid/(L2*L1)) % L3;
int x4 = sid/(L3*L2*L1);

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
    // 1 0 0 i 
    // 0 1 i 0 
    // 0 -i 1 0 
    // -i 0 0 1 
    
    int sp_idx = sid;
    int ga_idx = sp_idx + (0/2)*L*(20/4);
    
    // read spinor from device memory
    float4 I0 = tex1Dfetch(spinorTex, sp_idx + 0*L);
    float4 I1 = tex1Dfetch(spinorTex, sp_idx + 1*L);
    float4 I2 = tex1Dfetch(spinorTex, sp_idx + 2*L);
    float4 I3 = tex1Dfetch(spinorTex, sp_idx + 3*L);
    float4 I4 = tex1Dfetch(spinorTex, sp_idx + 4*L);
    float4 I5 = tex1Dfetch(spinorTex, sp_idx + 5*L);
    
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
    float4 G0 = tex1Dfetch(gaugeTex, ga_idx + 0*L);
    float4 G1 = tex1Dfetch(gaugeTex, ga_idx + 1*L);
    float4 G2 = tex1Dfetch(gaugeTex, ga_idx + 2*L);
    float4 G3 = tex1Dfetch(gaugeTex, ga_idx + 3*L);
    float4 G4 = tex1Dfetch(gaugeTex, ga_idx + 4*L);
    
    // multiply row 0 by half spinors
    {
        float A_re = + (g00_re * a0_re - g00_im * a0_im) + (g01_re * a1_re - g01_im * a1_im) + (g02_re * a2_re - g02_im * a2_im);
        float A_im = + (g00_re * a0_im + g00_im * a0_re) + (g01_re * a1_im + g01_im * a1_re) + (g02_re * a2_im + g02_im * a2_re);
        float B_re = + (g00_re * b0_re - g00_im * b0_im) + (g01_re * b1_re - g01_im * b1_im) + (g02_re * b2_re - g02_im * b2_im);
        float B_im = + (g00_re * b0_im + g00_im * b0_re) + (g01_re * b1_im + g01_im * b1_re) + (g02_re * b2_im + g02_im * b2_re);
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
        float A_re = + (g10_re * a0_re - g10_im * a0_im) + (g11_re * a1_re - g11_im * a1_im) + (g12_re * a2_re - g12_im * a2_im);
        float A_im = + (g10_re * a0_im + g10_im * a0_re) + (g11_re * a1_im + g11_im * a1_re) + (g12_re * a2_im + g12_im * a2_re);
        float B_re = + (g10_re * b0_re - g10_im * b0_im) + (g11_re * b1_re - g11_im * b1_im) + (g12_re * b2_re - g12_im * b2_im);
        float B_im = + (g10_re * b0_im + g10_im * b0_re) + (g11_re * b1_im + g11_im * b1_re) + (g12_re * b2_im + g12_im * b2_re);
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
        float A_re = + (g20_re * a0_re - g20_im * a0_im) + (g21_re * a1_re - g21_im * a1_im) + (g22_re * a2_re - g22_im * a2_im);
        float A_im = + (g20_re * a0_im + g20_im * a0_re) + (g21_re * a1_im + g21_im * a1_re) + (g22_re * a2_im + g22_im * a2_re);
        float B_re = + (g20_re * b0_re - g20_im * b0_im) + (g21_re * b1_re - g21_im * b1_im) + (g22_re * b2_re - g22_im * b2_im);
        float B_im = + (g20_re * b0_im + g20_im * b0_re) + (g21_re * b1_im + g21_im * b1_re) + (g22_re * b2_im + g22_im * b2_re);
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
    // Projector 1
    // 1 0 0 -i 
    // 0 1 -i 0 
    // 0 i 1 0 
    // i 0 0 1 
    
    int sp_idx = (x1==L1-1) ? sid-(L1-1) : sid+1;
    int ga_idx = sp_idx + (1/2)*L*(20/4);
    
    // read spinor from device memory
    float4 I0 = tex1Dfetch(spinorTex, sp_idx + 0*L);
    float4 I1 = tex1Dfetch(spinorTex, sp_idx + 1*L);
    float4 I2 = tex1Dfetch(spinorTex, sp_idx + 2*L);
    float4 I3 = tex1Dfetch(spinorTex, sp_idx + 3*L);
    float4 I4 = tex1Dfetch(spinorTex, sp_idx + 4*L);
    float4 I5 = tex1Dfetch(spinorTex, sp_idx + 5*L);
    
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
    float4 G0 = tex1Dfetch(gaugeTex, ga_idx + 0*L);
    float4 G1 = tex1Dfetch(gaugeTex, ga_idx + 1*L);
    float4 G2 = tex1Dfetch(gaugeTex, ga_idx + 2*L);
    float4 G3 = tex1Dfetch(gaugeTex, ga_idx + 3*L);
    float4 G4 = tex1Dfetch(gaugeTex, ga_idx + 4*L);
    
    // multiply row 0 by half spinors
    {
        float A_re = + (gT00_re * a0_re - gT00_im * a0_im) + (gT01_re * a1_re - gT01_im * a1_im) + (gT02_re * a2_re - gT02_im * a2_im);
        float A_im = + (gT00_re * a0_im + gT00_im * a0_re) + (gT01_re * a1_im + gT01_im * a1_re) + (gT02_re * a2_im + gT02_im * a2_re);
        float B_re = + (gT00_re * b0_re - gT00_im * b0_im) + (gT01_re * b1_re - gT01_im * b1_im) + (gT02_re * b2_re - gT02_im * b2_im);
        float B_im = + (gT00_re * b0_im + gT00_im * b0_re) + (gT01_re * b1_im + gT01_im * b1_re) + (gT02_re * b2_im + gT02_im * b2_re);
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
        float A_re = + (gT10_re * a0_re - gT10_im * a0_im) + (gT11_re * a1_re - gT11_im * a1_im) + (gT12_re * a2_re - gT12_im * a2_im);
        float A_im = + (gT10_re * a0_im + gT10_im * a0_re) + (gT11_re * a1_im + gT11_im * a1_re) + (gT12_re * a2_im + gT12_im * a2_re);
        float B_re = + (gT10_re * b0_re - gT10_im * b0_im) + (gT11_re * b1_re - gT11_im * b1_im) + (gT12_re * b2_re - gT12_im * b2_im);
        float B_im = + (gT10_re * b0_im + gT10_im * b0_re) + (gT11_re * b1_im + gT11_im * b1_re) + (gT12_re * b2_im + gT12_im * b2_re);
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
        float A_re = + (gT20_re * a0_re - gT20_im * a0_im) + (gT21_re * a1_re - gT21_im * a1_im) + (gT22_re * a2_re - gT22_im * a2_im);
        float A_im = + (gT20_re * a0_im + gT20_im * a0_re) + (gT21_re * a1_im + gT21_im * a1_re) + (gT22_re * a2_im + gT22_im * a2_re);
        float B_re = + (gT20_re * b0_re - gT20_im * b0_im) + (gT21_re * b1_re - gT21_im * b1_im) + (gT22_re * b2_re - gT22_im * b2_im);
        float B_im = + (gT20_re * b0_im + gT20_im * b0_re) + (gT21_re * b1_im + gT21_im * b1_re) + (gT22_re * b2_im + gT22_im * b2_re);
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
    // Projector 2
    // 1 0 0 -1 
    // 0 1 1 0 
    // 0 1 1 0 
    // -1 0 0 1 
    
    int sp_idx = (x2==0) ? (sid+(L2-1)*L1) : (sid-L1);
    int ga_idx = sp_idx + (2/2)*L*(20/4);
    
    // read spinor from device memory
    float4 I0 = tex1Dfetch(spinorTex, sp_idx + 0*L);
    float4 I1 = tex1Dfetch(spinorTex, sp_idx + 1*L);
    float4 I2 = tex1Dfetch(spinorTex, sp_idx + 2*L);
    float4 I3 = tex1Dfetch(spinorTex, sp_idx + 3*L);
    float4 I4 = tex1Dfetch(spinorTex, sp_idx + 4*L);
    float4 I5 = tex1Dfetch(spinorTex, sp_idx + 5*L);
    
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
    float4 G0 = tex1Dfetch(gaugeTex, ga_idx + 0*L);
    float4 G1 = tex1Dfetch(gaugeTex, ga_idx + 1*L);
    float4 G2 = tex1Dfetch(gaugeTex, ga_idx + 2*L);
    float4 G3 = tex1Dfetch(gaugeTex, ga_idx + 3*L);
    float4 G4 = tex1Dfetch(gaugeTex, ga_idx + 4*L);
    
    // multiply row 0 by half spinors
    {
        float A_re = + (g00_re * a0_re - g00_im * a0_im) + (g01_re * a1_re - g01_im * a1_im) + (g02_re * a2_re - g02_im * a2_im);
        float A_im = + (g00_re * a0_im + g00_im * a0_re) + (g01_re * a1_im + g01_im * a1_re) + (g02_re * a2_im + g02_im * a2_re);
        float B_re = + (g00_re * b0_re - g00_im * b0_im) + (g01_re * b1_re - g01_im * b1_im) + (g02_re * b2_re - g02_im * b2_im);
        float B_im = + (g00_re * b0_im + g00_im * b0_re) + (g01_re * b1_im + g01_im * b1_re) + (g02_re * b2_im + g02_im * b2_re);
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
        float A_re = + (g10_re * a0_re - g10_im * a0_im) + (g11_re * a1_re - g11_im * a1_im) + (g12_re * a2_re - g12_im * a2_im);
        float A_im = + (g10_re * a0_im + g10_im * a0_re) + (g11_re * a1_im + g11_im * a1_re) + (g12_re * a2_im + g12_im * a2_re);
        float B_re = + (g10_re * b0_re - g10_im * b0_im) + (g11_re * b1_re - g11_im * b1_im) + (g12_re * b2_re - g12_im * b2_im);
        float B_im = + (g10_re * b0_im + g10_im * b0_re) + (g11_re * b1_im + g11_im * b1_re) + (g12_re * b2_im + g12_im * b2_re);
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
        float A_re = + (g20_re * a0_re - g20_im * a0_im) + (g21_re * a1_re - g21_im * a1_im) + (g22_re * a2_re - g22_im * a2_im);
        float A_im = + (g20_re * a0_im + g20_im * a0_re) + (g21_re * a1_im + g21_im * a1_re) + (g22_re * a2_im + g22_im * a2_re);
        float B_re = + (g20_re * b0_re - g20_im * b0_im) + (g21_re * b1_re - g21_im * b1_im) + (g22_re * b2_re - g22_im * b2_im);
        float B_im = + (g20_re * b0_im + g20_im * b0_re) + (g21_re * b1_im + g21_im * b1_re) + (g22_re * b2_im + g22_im * b2_re);
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
    // Projector 3
    // 1 0 0 1 
    // 0 1 -1 0 
    // 0 -1 1 0 
    // 1 0 0 1 
    
    int sp_idx = (x2==L2-1) ? (sid-(L2-1)*L1) : (sid+L1);
    int ga_idx = sp_idx + (3/2)*L*(20/4);
    
    // read spinor from device memory
    float4 I0 = tex1Dfetch(spinorTex, sp_idx + 0*L);
    float4 I1 = tex1Dfetch(spinorTex, sp_idx + 1*L);
    float4 I2 = tex1Dfetch(spinorTex, sp_idx + 2*L);
    float4 I3 = tex1Dfetch(spinorTex, sp_idx + 3*L);
    float4 I4 = tex1Dfetch(spinorTex, sp_idx + 4*L);
    float4 I5 = tex1Dfetch(spinorTex, sp_idx + 5*L);
    
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
    float4 G0 = tex1Dfetch(gaugeTex, ga_idx + 0*L);
    float4 G1 = tex1Dfetch(gaugeTex, ga_idx + 1*L);
    float4 G2 = tex1Dfetch(gaugeTex, ga_idx + 2*L);
    float4 G3 = tex1Dfetch(gaugeTex, ga_idx + 3*L);
    float4 G4 = tex1Dfetch(gaugeTex, ga_idx + 4*L);
    
    // multiply row 0 by half spinors
    {
        float A_re = + (gT00_re * a0_re - gT00_im * a0_im) + (gT01_re * a1_re - gT01_im * a1_im) + (gT02_re * a2_re - gT02_im * a2_im);
        float A_im = + (gT00_re * a0_im + gT00_im * a0_re) + (gT01_re * a1_im + gT01_im * a1_re) + (gT02_re * a2_im + gT02_im * a2_re);
        float B_re = + (gT00_re * b0_re - gT00_im * b0_im) + (gT01_re * b1_re - gT01_im * b1_im) + (gT02_re * b2_re - gT02_im * b2_im);
        float B_im = + (gT00_re * b0_im + gT00_im * b0_re) + (gT01_re * b1_im + gT01_im * b1_re) + (gT02_re * b2_im + gT02_im * b2_re);
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
        float A_re = + (gT10_re * a0_re - gT10_im * a0_im) + (gT11_re * a1_re - gT11_im * a1_im) + (gT12_re * a2_re - gT12_im * a2_im);
        float A_im = + (gT10_re * a0_im + gT10_im * a0_re) + (gT11_re * a1_im + gT11_im * a1_re) + (gT12_re * a2_im + gT12_im * a2_re);
        float B_re = + (gT10_re * b0_re - gT10_im * b0_im) + (gT11_re * b1_re - gT11_im * b1_im) + (gT12_re * b2_re - gT12_im * b2_im);
        float B_im = + (gT10_re * b0_im + gT10_im * b0_re) + (gT11_re * b1_im + gT11_im * b1_re) + (gT12_re * b2_im + gT12_im * b2_re);
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
        float A_re = + (gT20_re * a0_re - gT20_im * a0_im) + (gT21_re * a1_re - gT21_im * a1_im) + (gT22_re * a2_re - gT22_im * a2_im);
        float A_im = + (gT20_re * a0_im + gT20_im * a0_re) + (gT21_re * a1_im + gT21_im * a1_re) + (gT22_re * a2_im + gT22_im * a2_re);
        float B_re = + (gT20_re * b0_re - gT20_im * b0_im) + (gT21_re * b1_re - gT21_im * b1_im) + (gT22_re * b2_re - gT22_im * b2_im);
        float B_im = + (gT20_re * b0_im + gT20_im * b0_re) + (gT21_re * b1_im + gT21_im * b1_re) + (gT22_re * b2_im + gT22_im * b2_re);
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
    // Projector 4
    // 1 0 i 0 
    // 0 1 0 -i 
    // -i 0 1 0 
    // 0 i 0 1 
    
    int sp_idx = (x3==0) ? (sid+(L3-1)*L2*L1) : (sid-L2*L1);
    int ga_idx = sp_idx + (4/2)*L*(20/4);
    
    // read spinor from device memory
    float4 I0 = tex1Dfetch(spinorTex, sp_idx + 0*L);
    float4 I1 = tex1Dfetch(spinorTex, sp_idx + 1*L);
    float4 I2 = tex1Dfetch(spinorTex, sp_idx + 2*L);
    float4 I3 = tex1Dfetch(spinorTex, sp_idx + 3*L);
    float4 I4 = tex1Dfetch(spinorTex, sp_idx + 4*L);
    float4 I5 = tex1Dfetch(spinorTex, sp_idx + 5*L);
    
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
    float4 G0 = tex1Dfetch(gaugeTex, ga_idx + 0*L);
    float4 G1 = tex1Dfetch(gaugeTex, ga_idx + 1*L);
    float4 G2 = tex1Dfetch(gaugeTex, ga_idx + 2*L);
    float4 G3 = tex1Dfetch(gaugeTex, ga_idx + 3*L);
    float4 G4 = tex1Dfetch(gaugeTex, ga_idx + 4*L);
    
    // multiply row 0 by half spinors
    {
        float A_re = + (g00_re * a0_re - g00_im * a0_im) + (g01_re * a1_re - g01_im * a1_im) + (g02_re * a2_re - g02_im * a2_im);
        float A_im = + (g00_re * a0_im + g00_im * a0_re) + (g01_re * a1_im + g01_im * a1_re) + (g02_re * a2_im + g02_im * a2_re);
        float B_re = + (g00_re * b0_re - g00_im * b0_im) + (g01_re * b1_re - g01_im * b1_im) + (g02_re * b2_re - g02_im * b2_im);
        float B_im = + (g00_re * b0_im + g00_im * b0_re) + (g01_re * b1_im + g01_im * b1_re) + (g02_re * b2_im + g02_im * b2_re);
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
        float A_re = + (g10_re * a0_re - g10_im * a0_im) + (g11_re * a1_re - g11_im * a1_im) + (g12_re * a2_re - g12_im * a2_im);
        float A_im = + (g10_re * a0_im + g10_im * a0_re) + (g11_re * a1_im + g11_im * a1_re) + (g12_re * a2_im + g12_im * a2_re);
        float B_re = + (g10_re * b0_re - g10_im * b0_im) + (g11_re * b1_re - g11_im * b1_im) + (g12_re * b2_re - g12_im * b2_im);
        float B_im = + (g10_re * b0_im + g10_im * b0_re) + (g11_re * b1_im + g11_im * b1_re) + (g12_re * b2_im + g12_im * b2_re);
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
        float A_re = + (g20_re * a0_re - g20_im * a0_im) + (g21_re * a1_re - g21_im * a1_im) + (g22_re * a2_re - g22_im * a2_im);
        float A_im = + (g20_re * a0_im + g20_im * a0_re) + (g21_re * a1_im + g21_im * a1_re) + (g22_re * a2_im + g22_im * a2_re);
        float B_re = + (g20_re * b0_re - g20_im * b0_im) + (g21_re * b1_re - g21_im * b1_im) + (g22_re * b2_re - g22_im * b2_im);
        float B_im = + (g20_re * b0_im + g20_im * b0_re) + (g21_re * b1_im + g21_im * b1_re) + (g22_re * b2_im + g22_im * b2_re);
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
    // Projector 5
    // 1 0 -i 0 
    // 0 1 0 i 
    // i 0 1 0 
    // 0 -i 0 1 
    
    int sp_idx = (x3==L3-1) ? (sid-(L3-1)*L2*L1) : (sid+L2*L1);
    int ga_idx = sp_idx + (5/2)*L*(20/4);
    
    // read spinor from device memory
    float4 I0 = tex1Dfetch(spinorTex, sp_idx + 0*L);
    float4 I1 = tex1Dfetch(spinorTex, sp_idx + 1*L);
    float4 I2 = tex1Dfetch(spinorTex, sp_idx + 2*L);
    float4 I3 = tex1Dfetch(spinorTex, sp_idx + 3*L);
    float4 I4 = tex1Dfetch(spinorTex, sp_idx + 4*L);
    float4 I5 = tex1Dfetch(spinorTex, sp_idx + 5*L);
    
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
    float4 G0 = tex1Dfetch(gaugeTex, ga_idx + 0*L);
    float4 G1 = tex1Dfetch(gaugeTex, ga_idx + 1*L);
    float4 G2 = tex1Dfetch(gaugeTex, ga_idx + 2*L);
    float4 G3 = tex1Dfetch(gaugeTex, ga_idx + 3*L);
    float4 G4 = tex1Dfetch(gaugeTex, ga_idx + 4*L);
    
    // multiply row 0 by half spinors
    {
        float A_re = + (gT00_re * a0_re - gT00_im * a0_im) + (gT01_re * a1_re - gT01_im * a1_im) + (gT02_re * a2_re - gT02_im * a2_im);
        float A_im = + (gT00_re * a0_im + gT00_im * a0_re) + (gT01_re * a1_im + gT01_im * a1_re) + (gT02_re * a2_im + gT02_im * a2_re);
        float B_re = + (gT00_re * b0_re - gT00_im * b0_im) + (gT01_re * b1_re - gT01_im * b1_im) + (gT02_re * b2_re - gT02_im * b2_im);
        float B_im = + (gT00_re * b0_im + gT00_im * b0_re) + (gT01_re * b1_im + gT01_im * b1_re) + (gT02_re * b2_im + gT02_im * b2_re);
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
        float A_re = + (gT10_re * a0_re - gT10_im * a0_im) + (gT11_re * a1_re - gT11_im * a1_im) + (gT12_re * a2_re - gT12_im * a2_im);
        float A_im = + (gT10_re * a0_im + gT10_im * a0_re) + (gT11_re * a1_im + gT11_im * a1_re) + (gT12_re * a2_im + gT12_im * a2_re);
        float B_re = + (gT10_re * b0_re - gT10_im * b0_im) + (gT11_re * b1_re - gT11_im * b1_im) + (gT12_re * b2_re - gT12_im * b2_im);
        float B_im = + (gT10_re * b0_im + gT10_im * b0_re) + (gT11_re * b1_im + gT11_im * b1_re) + (gT12_re * b2_im + gT12_im * b2_re);
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
        float A_re = + (gT20_re * a0_re - gT20_im * a0_im) + (gT21_re * a1_re - gT21_im * a1_im) + (gT22_re * a2_re - gT22_im * a2_im);
        float A_im = + (gT20_re * a0_im + gT20_im * a0_re) + (gT21_re * a1_im + gT21_im * a1_re) + (gT22_re * a2_im + gT22_im * a2_re);
        float B_re = + (gT20_re * b0_re - gT20_im * b0_im) + (gT21_re * b1_re - gT21_im * b1_im) + (gT22_re * b2_re - gT22_im * b2_im);
        float B_im = + (gT20_re * b0_im + gT20_im * b0_re) + (gT21_re * b1_im + gT21_im * b1_re) + (gT22_re * b2_im + gT22_im * b2_re);
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
    // Projector 6
    // 1 0 1 0 
    // 0 1 0 1 
    // 1 0 1 0 
    // 0 1 0 1 
    
    int sp_idx = (x4==0) ? (sid+(L4-1)*L3*L2*L1) : (sid-L3*L2*L1);
    int ga_idx = sp_idx + (6/2)*L*(20/4);
    
    // read spinor from device memory
    float4 I0 = tex1Dfetch(spinorTex, sp_idx + 0*L);
    float4 I1 = tex1Dfetch(spinorTex, sp_idx + 1*L);
    float4 I2 = tex1Dfetch(spinorTex, sp_idx + 2*L);
    float4 I3 = tex1Dfetch(spinorTex, sp_idx + 3*L);
    float4 I4 = tex1Dfetch(spinorTex, sp_idx + 4*L);
    float4 I5 = tex1Dfetch(spinorTex, sp_idx + 5*L);
    
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
    float4 G0 = tex1Dfetch(gaugeTex, ga_idx + 0*L);
    float4 G1 = tex1Dfetch(gaugeTex, ga_idx + 1*L);
    float4 G2 = tex1Dfetch(gaugeTex, ga_idx + 2*L);
    float4 G3 = tex1Dfetch(gaugeTex, ga_idx + 3*L);
    float4 G4 = tex1Dfetch(gaugeTex, ga_idx + 4*L);
    
    // multiply row 0 by half spinors
    {
        float A_re = + (g00_re * a0_re - g00_im * a0_im) + (g01_re * a1_re - g01_im * a1_im) + (g02_re * a2_re - g02_im * a2_im);
        float A_im = + (g00_re * a0_im + g00_im * a0_re) + (g01_re * a1_im + g01_im * a1_re) + (g02_re * a2_im + g02_im * a2_re);
        float B_re = + (g00_re * b0_re - g00_im * b0_im) + (g01_re * b1_re - g01_im * b1_im) + (g02_re * b2_re - g02_im * b2_im);
        float B_im = + (g00_re * b0_im + g00_im * b0_re) + (g01_re * b1_im + g01_im * b1_re) + (g02_re * b2_im + g02_im * b2_re);
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
        float A_re = + (g10_re * a0_re - g10_im * a0_im) + (g11_re * a1_re - g11_im * a1_im) + (g12_re * a2_re - g12_im * a2_im);
        float A_im = + (g10_re * a0_im + g10_im * a0_re) + (g11_re * a1_im + g11_im * a1_re) + (g12_re * a2_im + g12_im * a2_re);
        float B_re = + (g10_re * b0_re - g10_im * b0_im) + (g11_re * b1_re - g11_im * b1_im) + (g12_re * b2_re - g12_im * b2_im);
        float B_im = + (g10_re * b0_im + g10_im * b0_re) + (g11_re * b1_im + g11_im * b1_re) + (g12_re * b2_im + g12_im * b2_re);
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
        float A_re = + (g20_re * a0_re - g20_im * a0_im) + (g21_re * a1_re - g21_im * a1_im) + (g22_re * a2_re - g22_im * a2_im);
        float A_im = + (g20_re * a0_im + g20_im * a0_re) + (g21_re * a1_im + g21_im * a1_re) + (g22_re * a2_im + g22_im * a2_re);
        float B_re = + (g20_re * b0_re - g20_im * b0_im) + (g21_re * b1_re - g21_im * b1_im) + (g22_re * b2_re - g22_im * b2_im);
        float B_im = + (g20_re * b0_im + g20_im * b0_re) + (g21_re * b1_im + g21_im * b1_re) + (g22_re * b2_im + g22_im * b2_re);
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

if(1)
{
    // Projector 7
    // 1 0 -1 0 
    // 0 1 0 -1 
    // -1 0 1 0 
    // 0 -1 0 1 
    
    int sp_idx = (x4==L4-1) ? (sid-(L4-1)*L3*L2*L1) : (sid+L3*L2*L1);
    int ga_idx = sp_idx + (7/2)*L*(20/4);
    
    // read spinor from device memory
    float4 I0 = tex1Dfetch(spinorTex, sp_idx + 0*L);
    float4 I1 = tex1Dfetch(spinorTex, sp_idx + 1*L);
    float4 I2 = tex1Dfetch(spinorTex, sp_idx + 2*L);
    float4 I3 = tex1Dfetch(spinorTex, sp_idx + 3*L);
    float4 I4 = tex1Dfetch(spinorTex, sp_idx + 4*L);
    float4 I5 = tex1Dfetch(spinorTex, sp_idx + 5*L);
    
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
    float4 G0 = tex1Dfetch(gaugeTex, ga_idx + 0*L);
    float4 G1 = tex1Dfetch(gaugeTex, ga_idx + 1*L);
    float4 G2 = tex1Dfetch(gaugeTex, ga_idx + 2*L);
    float4 G3 = tex1Dfetch(gaugeTex, ga_idx + 3*L);
    float4 G4 = tex1Dfetch(gaugeTex, ga_idx + 4*L);
    
    // multiply row 0 by half spinors
    {
        float A_re = + (gT00_re * a0_re - gT00_im * a0_im) + (gT01_re * a1_re - gT01_im * a1_im) + (gT02_re * a2_re - gT02_im * a2_im);
        float A_im = + (gT00_re * a0_im + gT00_im * a0_re) + (gT01_re * a1_im + gT01_im * a1_re) + (gT02_re * a2_im + gT02_im * a2_re);
        float B_re = + (gT00_re * b0_re - gT00_im * b0_im) + (gT01_re * b1_re - gT01_im * b1_im) + (gT02_re * b2_re - gT02_im * b2_im);
        float B_im = + (gT00_re * b0_im + gT00_im * b0_re) + (gT01_re * b1_im + gT01_im * b1_re) + (gT02_re * b2_im + gT02_im * b2_re);
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
        float A_re = + (gT10_re * a0_re - gT10_im * a0_im) + (gT11_re * a1_re - gT11_im * a1_im) + (gT12_re * a2_re - gT12_im * a2_im);
        float A_im = + (gT10_re * a0_im + gT10_im * a0_re) + (gT11_re * a1_im + gT11_im * a1_re) + (gT12_re * a2_im + gT12_im * a2_re);
        float B_re = + (gT10_re * b0_re - gT10_im * b0_im) + (gT11_re * b1_re - gT11_im * b1_im) + (gT12_re * b2_re - gT12_im * b2_im);
        float B_im = + (gT10_re * b0_im + gT10_im * b0_re) + (gT11_re * b1_im + gT11_im * b1_re) + (gT12_re * b2_im + gT12_im * b2_re);
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
        float A_re = + (gT20_re * a0_re - gT20_im * a0_im) + (gT21_re * a1_re - gT21_im * a1_im) + (gT22_re * a2_re - gT22_im * a2_im);
        float A_im = + (gT20_re * a0_im + gT20_im * a0_re) + (gT21_re * a1_im + gT21_im * a1_re) + (gT22_re * a2_im + gT22_im * a2_re);
        float B_re = + (gT20_re * b0_re - gT20_im * b0_im) + (gT21_re * b1_re - gT21_im * b1_im) + (gT22_re * b2_re - gT22_im * b2_im);
        float B_im = + (gT20_re * b0_im + gT20_im * b0_re) + (gT21_re * b1_im + gT21_im * b1_re) + (gT22_re * b2_im + gT22_im * b2_re);
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


//g_out[0*L+sid] = make_float4(o00_re, o00_im, o01_re, o01_im);
//g_out[1*L+sid] = make_float4(o02_re, o02_im, o10_re, o10_im);
//g_out[2*L+sid] = make_float4(o11_re, o11_im, o12_re, o12_im);
//g_out[3*L+sid] = make_float4(o20_re, o20_im, o21_re, o21_im);
//g_out[4*L+sid] = make_float4(o22_re, o22_im, o30_re, o30_im);
//g_out[5*L+sid] = make_float4(o31_re, o31_im, o32_re, o32_im);

((float*)g_out)[0*L+sid] = o00_re;
((float*)g_out)[1*L+sid] = o00_im;
((float*)g_out)[2*L+sid] = o01_re;
((float*)g_out)[3*L+sid] = o01_im;
((float*)g_out)[4*L+sid] = o02_re;
((float*)g_out)[5*L+sid] = o02_im;
((float*)g_out)[6*L+sid] = o10_re;
((float*)g_out)[7*L+sid] = o10_im;
((float*)g_out)[8*L+sid] = o11_re;
((float*)g_out)[9*L+sid] = o11_im;
((float*)g_out)[10*L+sid] = o12_re;
((float*)g_out)[11*L+sid] = o12_im;
((float*)g_out)[12*L+sid] = o20_re;
((float*)g_out)[13*L+sid] = o20_im;
((float*)g_out)[14*L+sid] = o21_re;
((float*)g_out)[15*L+sid] = o21_im;
((float*)g_out)[16*L+sid] = o22_re;
((float*)g_out)[17*L+sid] = o22_im;
((float*)g_out)[18*L+sid] = o30_re;
((float*)g_out)[19*L+sid] = o30_im;
((float*)g_out)[20*L+sid] = o31_re;
((float*)g_out)[21*L+sid] = o31_im;
((float*)g_out)[22*L+sid] = o32_re;
((float*)g_out)[23*L+sid] = o32_im;

