package dslash

import dslash.DiracSpinor._


class CodeGen(sharedFloats: Int, dagger: Boolean) {
  def block(code: String) = {
    "{\n"+code.lines.map{"    "+_+"\n"}.reduceLeft[String](_+_)+"}\n"
  }
  
  def sign(x:Double) = {
    x match {case 1 => "+"; case -1 => "-"}
  }

  def nthFloat4(n:Int) = {
    (n/4) + "." + Array("x","y","z","w")(n%4)
  }

  def in_re(s:Int, c:Int) = "i"+s+c+"_re"
  def in_im(s:Int, c:Int) = "i"+s+c+"_im"
  def g_re(d:Int, m:Int, n:Int) = (if (d%2==0) "g" else "gT")+m+n+"_re"
  def g_im(d:Int, m:Int, n:Int) = (if (d%2==0) "g" else "gT")+m+n+"_im"
  def out_re(s:Int, c:Int) = "o"+s+c+"_re"
  def out_im(s:Int, c:Int) = "o"+s+c+"_im"
  def h1_re(h:Int, c:Int) = Array("a","b")(h)+c+"_re"
  def h1_im(h:Int, c:Int) = Array("a","b")(h)+c+"_im"
  def h2_re(h:Int) = Array("A","B")(h)+"_re"
  def h2_im(h:Int) = Array("A","B")(h)+"_im"
  
  
  def prolog() = {
    val str = new StringBuilder()
    
    str.append("#define SHARED_FLOATS_PER_THREAD "+sharedFloats+"\n")
    str.append("#define SHARED_BYTES (BLOCK_DIM*SHARED_FLOATS_PER_THREAD*sizeof(float))\n\n")
    
    for (s <- 0 until 4; c <- 0 until 3; i = 3*s+c) {
      str.append("#define "+in_re(s,c)+" I"+nthFloat4(2*i+0)+"\n")
      str.append("#define "+in_im(s,c)+" I"+nthFloat4(2*i+1)+"\n")
    }
    str.append("\n")
    for (m <- 0 until 3; n <- 0 until 3; i = 3*m+n) {
      str.append("#define "+g_re(0,m,n)+" G"+nthFloat4(2*i+0)+"\n")
      str.append("#define "+g_im(0,m,n)+" G"+nthFloat4(2*i+1)+"\n")
    }
    str.append("\n")
    for (m <- 0 until 3; n <- 0 until 3; i = 3*m+n) {
      str.append("#define "+g_re(1,m,n)+" (+"+g_re(0,n,m)+")\n")
      str.append("#define "+g_im(1,m,n)+" (-"+g_im(0,n,m)+")\n")
    }
    str.append("\n")
    
    for (s <- 0 until 4; c <- 0 until 3; i = 3*s+c) {
      if (2*i < sharedFloats)
        str.append("#define "+out_re(s,c)+" s["+(2*i+0)+"]\n")
      else
        str.append("volatile float "+out_re(s,c)+";\n")
      if (2*i+1 < sharedFloats)
        str.append("#define "+out_im(s,c)+" s["+(2*i+1)+"]\n")
      else
        str.append("volatile float "+out_im(s,c)+";\n")        
    }
    str.append("\n")
    
    str.append(
"""
// Performs the complex conjugated accumulation: a += b* c*
#define ACC_CONJ_PROD(a, b, c) \
    a##_re += b##_re * c##_re - b##_im * c##_im, \
    a##_im -= b##_re * c##_im + b##_im * c##_re

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

#define READ_SPINOR(spinor) \
    float4 I0 = tex1Dfetch((spinor), sp_idx + 0*Nh); \
    float4 I1 = tex1Dfetch((spinor), sp_idx + 1*Nh); \
    float4 I2 = tex1Dfetch((spinor), sp_idx + 2*Nh); \
    float4 I3 = tex1Dfetch((spinor), sp_idx + 3*Nh); \
    float4 I4 = tex1Dfetch((spinor), sp_idx + 4*Nh); \
    float4 I5 = tex1Dfetch((spinor), sp_idx + 5*Nh);

int sid = BLOCK_DIM*blockIdx.x + threadIdx.x;
int boundaryCrossings = sid/L1h + sid/(L2*L1h) + sid/(L3*L2*L1h);
int X = 2*sid + (boundaryCrossings + oddBit) % 2;
int x4 = X/(L3*L2*L1);
int x3 = (X/(L2*L1)) % L3;
int x2 = (X/L1) % L2;
int x1 = X % L1;

"""	)

    str.append("extern __shared__ float s_data[];\n")
    str.append("volatile float *s = s_data+SHARED_FLOATS_PER_THREAD*threadIdx.x;\n\n")

    for (s <- 0 until 4; c <- 0 until 3) {
      str.append(out_re(s,c) + " = " + out_im(s,c)+" = 0;\n")
    }
    str.append("\n")
    
    str.toString()
  }
  
  def epilog() = {
    val str = new StringBuilder()
    str.append(
"""
#ifdef DSLASH_XPAY
    float4 accum0 = tex1Dfetch(accumTex, sid + 0*Nh);
    float4 accum1 = tex1Dfetch(accumTex, sid + 1*Nh);
    float4 accum2 = tex1Dfetch(accumTex, sid + 2*Nh);
    float4 accum3 = tex1Dfetch(accumTex, sid + 3*Nh);
    float4 accum4 = tex1Dfetch(accumTex, sid + 4*Nh);
    float4 accum5 = tex1Dfetch(accumTex, sid + 5*Nh);
""")
    
    for (s <- 0 until 4; c <- 0 until 3; i = 3*s+c) {
      str.append("    "+out_re(s,c) +" = a*"+out_re(s,c)+" + accum"+nthFloat4(2*i+0)+";\n")
      str.append("    "+out_im(s,c) +" = a*"+out_im(s,c)+" + accum"+nthFloat4(2*i+1)+";\n")
    }
    str.append("#endif\n\n")
    
    str.append(
"""
#ifdef WRITE_FLOAT4
// this code exhibits a hardware bug in our C870 card
g_out[0*Nh+sid] = make_float4(o00_re, o00_im, o01_re, o01_im);
g_out[1*Nh+sid] = make_float4(o02_re, o02_im, o10_re, o10_im);
g_out[2*Nh+sid] = make_float4(o11_re, o11_im, o12_re, o12_im);
g_out[3*Nh+sid] = make_float4(o20_re, o20_im, o21_re, o21_im);
g_out[4*Nh+sid] = make_float4(o22_re, o22_im, o30_re, o30_im);
g_out[5*Nh+sid] = make_float4(o31_re, o31_im, o32_re, o32_im);
#endif

#ifdef WRITE_FLOAT1_SMEM
int t = threadIdx.x;
int B = BLOCK_DIM;
int b = blockIdx.x;
int f = SHARED_FLOATS_PER_THREAD;
__syncthreads();
for (int i = 0; i < 6; i++) // spinor indices
    for (int c = 0; c < 4; c++) // components of float4
        ((float*)g_out)[i*(Nh*4) + b*(B*4) + c*(B) + t] = s_data[(c*B/4 + t/4)*(f) + i*(4) + t%4];
#endif

#ifdef WRITE_FLOAT1_STAGGERED
// the alternative to writing float4's directly: almost as fast, a lot more confusing
int t = threadIdx.x;
int B = BLOCK_DIM;
int b = blockIdx.x;
int f = SHARED_FLOATS_PER_THREAD;
__syncthreads();
for (int i = 0; i < 4; i++) // spinor indices
    for (int c = 0; c < 4; c++) // components of float4
        ((float*)g_out)[i*(Nh*4) + b*(B*4) + c*(B) + t] = s_data[(c*B/4 + t/4)*(f) + i*(4) + t%4];
__syncthreads();
s[0] = o22_re;
s[1] = o22_im;
s[2] = o30_re;
s[3] = o30_im;
s[4] = o31_re;
s[5] = o31_im;
s[6] = o32_re;
s[7] = o32_im;
__syncthreads();
for (int i = 0; i < 2; i++)
    for (int c = 0; c < 4; c++)
        ((float*)g_out)[(i+4)*(Nh*4) + b*(B*4) + c*(B) + t] = s_data[(c*B/4 + t/4)*(f) + i*(4) + t%4];
#endif

""")
  }
  
  
  val projectors = Array(
      id-gamma1, id+gamma1,
      id-gamma2, id+gamma2,
      id-gamma3, id+gamma3,
      id-gamma4, id+gamma4
  )
  
  def gen(dir: Int):String = {
    val proj_idx = if (!dagger) dir else dir + (1 - 2*(dir%2))
    val proj = projectors(proj_idx)
    
    // if row(i) = (j, c), then the i'th row of the projector can be represented
    // as a multiple of the j'th row: i = c j
    def row(i:Int) = {
      require(i == 2 || i == 3)
      (proj(i,0), proj(i,1)) match {
      case (Complex(0,0), c) => (1, c)
      case (c, Complex(0,0)) => (0, c)
      }
    }
    
    val str = new StringBuilder()
    
    val projStr = "P"+(dir/2)+Array("-","+")(proj_idx%2)
    str.append("// Projector "+projStr+"\n")
    proj.toString.lines.foreach{l => str.append("// "+l+"\n")}
    str.append("\n")
    
    str.append ( (dir match {
    case 0 => "int sp_idx = ((x1==L1-1) ? X-(L1-1) : X+1) / 2;\n"
    case 1 => "int sp_idx = ((x1==0)    ? X+(L1-1) : X-1) / 2;\n"
    case 2 => "int sp_idx = ((x2==L2-1) ? X-(L2-1)*L1 : X+L1) / 2;\n"
    case 3 => "int sp_idx = ((x2==0)    ? X+(L2-1)*L1 : X-L1) / 2;\n"
    case 4 => "int sp_idx = ((x3==L3-1) ? X-(L3-1)*L2*L1 : X+L2*L1) / 2;\n"
    case 5 => "int sp_idx = ((x3==0)    ? X+(L3-1)*L2*L1 : X-L2*L1) / 2;\n"
    case 6 => "int sp_idx = ((x4==L4-1) ? X-(L4-1)*L3*L2*L1 : X+L3*L2*L1) / 2;\n"
    case 7 => "int sp_idx = ((x4==0)    ? X+(L4-1)*L3*L2*L1 : X-L3*L2*L1) / 2;\n"
    }))
    val baseIdx = if (dir % 2 == 0) "sid" else "sp_idx"
    str.append("int ga_idx = "+baseIdx+" + ("+dir+"/2)*Nh*3;\n\n")
    
    str.append("// read spinor from device memory\n")
    str.append("READ_SPINOR(spinorTex);\n\n")
    
    str.append("// project spinor into half spinors\n")
    for (h <- 0 until 2) {
      for (c <- 0 until 3) {
        val strRe = new StringBuilder()
        val strIm = new StringBuilder()
        for (s <- 0 until 4) {
          proj(h, s) match {
          case Complex(0, 0) => ()
          case Complex(re, 0) => {
            strRe.append(sign(re)+in_re(s,c))
            strIm.append(sign(re)+in_im(s,c))
          }
          case Complex(0, im) => 
            strRe.append(sign(-im)+in_im(s,c))
            strIm.append(sign(im)+in_re(s,c))
          }
        }
        str.append("float "+h1_re(h,c)+ " = "+strRe+";\n")
        str.append("float "+h1_im(h,c)+ " = "+strIm+";\n")
      }
      str.append("\n")
    }
    
    str.append("// read gauge matrix from device memory\n")
    str.append("READ_GAUGE_MATRIX(gauge"+(dir%2)+"Tex);\n\n")
    
    for (m <- 0 until 3) {
      str.append("// multiply row "+m+" by half spinors\n")
      val mstr = new StringBuilder
      
      for (h <- 0 until 2) {
        val re = new StringBuilder("float "+h2_re(h)+" =")
        val im = new StringBuilder("float "+h2_im(h)+" =")
        for (c <- 0 until 3) {
          re.append(" + ("+g_re(dir,m,c)+" * "+h1_re(h,c)+" - "+g_im(dir,m,c)+" * "+h1_im(h,c)+")")
          im.append(" + ("+g_re(dir,m,c)+" * "+h1_im(h,c)+" + "+g_im(dir,m,c)+" * "+h1_re(h,c)+")")
        }
        mstr.append(re.toString+";\n")
        mstr.append(im.toString+";\n")
      }
      
      mstr.append(out_re(0, m) + " += +" + h2_re(0) + ";\n")
      mstr.append(out_im(0, m) + " += +" + h2_im(0) + ";\n")
      mstr.append(out_re(1, m) + " += +" + h2_re(1) + ";\n")
      mstr.append(out_im(1, m) + " += +" + h2_im(1) + ";\n")
    
      for (s <- 2 until 4) {
        row(s) match {
        case (h, Complex(re, 0)) => {
          mstr.append(out_re(s, m) + " += " + sign(re)+h2_re(h)+";\n")
          mstr.append(out_im(s, m) + " += " + sign(re)+h2_im(h)+";\n")
        }
        case (h, Complex(0, im)) => {
          mstr.append(out_re(s, m) + " += " + sign(-im)+h2_im(h)+";\n")
          mstr.append(out_im(s, m) + " += " + sign(+im)+h2_re(h)+";\n")
        }
        }
      }
      str.append(block(mstr.toString) + "\n")
    }

    (if (dir == 0) "if(1)\n" else "if(1)\n") + block(str.toString)+"\n"
  }
  
  def generate() = {
    prolog() + gen(0) + gen(1) + gen(2) + gen(3) + gen(4) + gen(5) + gen(6) + gen(7) + epilog()
  }
  
  def projectorsAsCString() = {
    val str = new StringBuilder
    str.append("{\n")
    for (d <- 0 until 8) {
      str.append("{\n")
      for (i <- 0 until 4) {
        str.append("  {")
        for (j <- 0 until 4) {
          val c = projectors(d)(i, j)
          str.append("{" + c.re.toInt +"," + c.im.toInt +"}")
          if (j < 3)
            str.append(", ")
        }
        str.append("}")
        if (i < 3)
          str.append(",")
        str.append("\n")
      }
      str.append("}")
      if (d < 7)
        str.append(",")
      str.append("\n")
    }
    str.append("}\n")
    
    str.toString
  }
}


object Dslash {
  def main(args: Array[String]) {
    val codeGen = new CodeGen(19, false)
    Console.println("// *** CUDA DSLASH ***\n\n" + codeGen.generate())
  }
}

object DslashDagger {
  def main(args: Array[String]) {
    val codeGen = new CodeGen(19, true)
    Console.println("// *** CUDA DSLASH DAGGER ***\n\n" + codeGen.generate())
  }
}
