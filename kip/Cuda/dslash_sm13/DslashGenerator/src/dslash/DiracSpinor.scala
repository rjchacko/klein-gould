package dslash

import dslash.Complex._

object DiracSpinor {
  val id = new DiracSpinor(
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
  )
  
  val gamma1 = new DiracSpinor (
     0,  0, 0, I,
     0,  0, I, 0,
     0, -I, 0, 0,
    -I,  0, 0, 0
  )
  
  val gamma2 = new DiracSpinor (
     0, 0, 0, -1,
     0, 0, 1,  0,
     0, 1, 0,  0,
    -1, 0, 0,  0
  )
  
  val gamma3 = new DiracSpinor (
     0, 0, I,  0,
     0, 0, 0, -I,
    -I, 0, 0,  0,
     0, I, 0,  0
  )

  val gamma4 = new DiracSpinor (
    0, 0, 1, 0,
    0, 0, 0, 1,
    1, 0, 0, 0,
    0, 1, 0, 0
  )
}

class DiracSpinor(val elems: Array[Complex]) {
  def this(elems: Complex*) = this(elems.toArray)
  
  def apply(i:Int, j:Int) = elems(i*4 + j)
  
  def +(that: DiracSpinor) = {
    new DiracSpinor(elems.zip(that.elems).map{case (c1, c2) => c1+c2}) 
  }
  
  def -(that: DiracSpinor) = {
    new DiracSpinor(elems.zip(that.elems).map{case (c1, c2) => c1-c2}) 
  }
  
  override def toString = {
    val str = new StringBuilder
    for (i <- 0 until 4) {
      for (j <- 0 until 4) {
        str.append(apply(i, j) + " ")
      }
      str.append("\n")
    }
    str.toString
  }

}

