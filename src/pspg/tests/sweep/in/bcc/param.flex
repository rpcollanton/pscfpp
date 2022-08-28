System{
  Mixture{
    nMonomer  2
    monomers  1.0  
              1.0 
    nPolymer  1
    Polymer{
      type    linear
      nBlock  2
      blocks  0   0.125
              1   0.875
      phi     1.0
    }
    ds   0.01
  }
  Interaction{
    chi  0   0    0.0
         1   0    41.0
         1   1    0.0
  }
  Domain{
    unitCell     cubic  1.7593562142
    mesh         32      32    32
    groupName    I_m_-3_m
  }
  AmIterator{
    maxItr      1000
    epsilon     1.0e-7
    maxHist     40
    isFlexible   1
  }
}
