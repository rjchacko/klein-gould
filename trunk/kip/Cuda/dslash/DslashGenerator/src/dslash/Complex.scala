package dslash


object Complex {
  val I = new Complex(0, 1) 
  implicit def num2complex[T<%double](x: T): Complex = new Complex(x, 0)
}

case class Complex(val re: double, val im: double) { 
  def unary_- : Complex =
    new Complex(-re, -im)
  def + (that: Complex): Complex =
    new Complex(this.re+ that.re, this.im + that.im)
  def - (that: Complex): Complex =
    new Complex(this.re-that.re, this.im -that.im)
  def * (that: Complex): Complex =
    new Complex(this.re* that.re-this.im * that.im, this.re * that.im+ this.im* that.re)
  def / (that: Complex): Complex = {
    val denom= that.re* that.re+ that.im* that.im 
    new Complex((this.re* that.re+ this.im* that.im) / denom,
                (this.im* that.re-this.re* that.im) / denom)
  }
  def === (that: Complex) : Boolean =
    re == that.re && im == that.im
  def !== (that: Complex) : Boolean =
    ! (this === that)
  
  override def toString = {
    def fltToString(a: Double) = (if (a == a.toInt) a.toInt else a).toString 
    def imToString(a: Double) = a match {
    case 0 	=> "0i"
    case -1 => "-i"
    case 1 	=> "i"
    case _ 	=> fltToString(a)+"i"
    }
    
    if (re == 0 && im == 0)
      "0"
    else if (re == 0)
      imToString(im)
    else if (im == 0)
      fltToString(re)
    else {
      val im_str = if (im < 0) "-"+imToString(-im) else "+"+imToString(im)
      fltToString(re)+im_str
    }
  }
}
