#!/usr/bin/env python

import h5py
import numpy as np

def umat(U,lorb, uvalue, jvalue,cub)
    l=2
    nlm=2*l+1
    lmaxw=4
    lmmax=(2*lmaxw+2)**2
    c = np.array((2*lmaxw+1,lmmax,lmmax));
    us = np.zeros((nlm,nlm,nlm,nlm))
    u = np.zeros((nlm,nlm,nlm,nlm))
    fk = np.array((l+1))
    yor = np.array((7,7))
    yoi = np.array((7,7))
    am1,am2,am3,am4,amz = 0.0
    amz = 0.0+0.0j
    pi = np.pi
    #----> fk - Slater integrals
    rcl_init(l,uvalue,jvalue,fk)
    print 'FK=',fk
    sgaunt(c,nlm)
    for i in range(0,2*lmaxw+2):
      for j in range(1,lmmax+1):
        print "Clebsh:", c[i,j,:]
    kmax=2*l
    for m1 in range(-l,l+1):
      nm1=m1+l
      lm1=l*(l+1)+m1+1
      for m2 in range(-l,l+1):
        nm2=m2+l
        lm2=l*(l+1)+m2+1
        for m3 in range(-l,l+1):
          nm3=m3+l
          lm3=l*(l+1)+m3+1
          for m4 in range(-l,l+1):
            nm4=m4+l
            lm4=l*(l+1)+m4
            uk=0.0
            for k in range(0,kmax+2,2):
              uq=0.0
              for mk in range(-k,k+1):
                if mk != m1-m3: continue
                cgk1=C(k/2,lm1,lm3)
                if cgk1 == 0.0: continue
                if mk != m4-m2: continue
                cgk2=C(k/2,lm4,lm2)
                if cgk2 == 0.0: continue
                uq=uq+cgk1*cgk2
              if uq == 0.0: continue
              nfk=k/2+1
              uk=uk+uq*fk(nfk)*4*pi/(2*k+1)
              print "add Uk", nm1,nm2,nm3,nm4, nfk
            us[nm1,nm2,nm3,nm4]=uk
    avu=0.0
    avj=0.0
    for i in range(-l,l+1):
      li=i+l
      for j in range(-l,l+1):
        lj=j+l
        avu=avu+us[li,lj,li,lj]
        avj=avj+(us[li,lj,li,lj]-us[li,lj,lj,li])
    avu=avu/(2*l+1)/(2*l+1)
    avj=avj/(2*l+1)/(2*l)
    avj=avu-avj
    print*,'u-av:',avu
    print*,'j-av:',avj
    #-------------- Calculate P-H symmetry energy (reference energy for Go)
    eph=0.0
    for m1 in range(nlm):
      for m2 in range(nlm):
        eph=eph+us[m1,m2,m1,m2]-us[m1,m2,m2,m1]/2.0
    #---------- Average over m:
    eph=eph/float(nlm)
    print 'Eph=',eph
    ctormt(yor,yoi,l)
    u = np.zeros((nlm,nlm,nlm,nlm))
    if cub:
      for ms1 in range(nlm):
        for ms2 in range(nlm):
          for ms3 in range(nlm):
            for ms4 in range(nlm):
              for ms5 in range(nlm):
                am1 = yor[ms1,ms5]-yoi[ms1,ms5]*1.0j
                if am1==amz: continue
                for ms6 in range(nlm):
                  am2 = yor[ms2,ms6]-yoi[ms2,ms6]*1.0j
                  if am2 == amz: continue
                  for ms7 in range(nlm):
                    am3 = yor[ms3,ms7]+yoi[ms3,ms7]*1.0j
                    if am3 == amz: continue
                    for ms8 in range(nlm):
                      am4 = yor[ms4,ms8]+yoi[ms4,ms8]*1.0j
                      if am4 == amz: continue
                      u[ms1,ms2,ms3,ms4]= u[ms1,ms2,ms3,ms4] + am1*am2*am3*am4*us[ms5,ms6,ms7,ms8]
    else:
      for ms1 in range(nlm):
        for ms2 in range(nlm):
          for ms3 in range(nlm):
            for ms4 in range(nlm)
              u[ms1,ms2,ms3,ms4] = us[ms1,ms2,ms3,ms4]


def rcl_init(l,uvalue,jvalue,rcl):
    #local variables
    ev2ry = 0.0
    uv = 0.0
    jv = 0.0
    i = 0
    print 'RCL_INIT: '
    #ev2ry = 1.0
    ev2ry = 13.6058
    # uv = uvalue / ev2ry
    #jv = jvalue / ev2ry
    uv = uvalue
    jv = jvalue
    rcl[0] = uv
    #p-shell
    if l == 1:
      rcl[1] = jv *5.0
      #d-shell
    elif l == 2:
      rcl[1] = jv * 14.0 / (1.0 + 0.63)
      rcl[2] = 0.63 * rcl[1]
      #f-shell
    elif l==3 :
      rcl[1] = 6435.0 * jv / (286.0 + 195.0 * 451.0 / 675.0 + 250.0 * 1001.0 / 2025.0)
      rcl[2] = 451.0 * rcl[1] / 675.0
      rcl[3] = 1001.0 * rcl[1] / 2025.0
    print '  F0   F1   F2   F3'
    print (rcl[i] for i in range(1, l+2))
    print 'RCL_INIT: exit '

#------------------------------------------------------
#************************************************************
#*   Calculation of the Gaunt coefficients C(L2M2,L1M1,LM)  *
#*							                                *
#*    l2m2                    /         *		            *
#*   C       =C(l2m2,l1m1,lm)=\dr*Y(r)*Y(r)*Y(r)		    *
#*    lm,l1m1		     /    lm   l1m1 l2m2	            *
#*							                                *
#*    and C.ne.0 when l2=/l1-l/,/l1-l/+2,...,l1+l,m2=m1-m   *
#*    Y(lm) is a complex spherical garmonic with a phase    *
#*    after Condon and Shortley				                *
#* Written by S.Yu.Savrasov (P.N.Lebedev Physical Institute)*
#************************************************************
def SGAUNT(C,LMM):
    LMAXW=4
    LMMAX=(2*LMAXW+2)**2
    C = np.array((2*LMAXW+2,LMMAX,LMMAX))
    PI = np.pi
    for L1 in range(0,LMM+1):
      for M1 in range(-L1,L1+1):
        for L  in range(0,LMM+1):
          for M  in range(-L,L+1):
            L1M1=L1*(L1+1)+M1+1
            LM=L*(L+1)+M+1
            for L2 in range(ABS(L1-L),L1+L+1,2):
              LL2=L2/2
              M2=M1-M
              if np.abs(M2)<=L2: #!!! selection rule
                AJ=L
                BJ=L2
                AM=M
                BM=M2
                CJ=L1
                CM=M1
                A1=CLEBSCH(AJ,BJ,0.0,0.0,CJ,0.0)  # Clebsch-Gordan coefficients
                A2=CLEBSCH(AJ,BJ,AM,BM,CJ,CM)        # Clebsch-Gordan coefficients
                DL1=2*L +1
                DL2=2*L2+1
                DL3=2*L1+1
                C[LL2,L1M1,LM]=A1*A2*np.sqrt(DL1*DL2/DL3/4.00/PI)
              elif np.abs(M2)>L2:
                C[LL2,L1M1,LM]=0.0


#******************************************************************
#* Program calculates Clebsch-Gordan coefficients                 *
#* See: Landau and Lifshitz, Vol.3                                *
#*      cj,cm                                                     *
#*     C                                                          *
#*      aj,am,bj,bm                                               *
#* Written by A.Soldatov (IAE)                                    *
#******************************************************************
def CLEBSCH(AJ,BJ,AM,BM,CJ,CM):
  F = np.array((100))
  N = 100
  K = 0
  X = 2.0
  F[0]=0.
  F[1]=0.
  if K <= 0:
   K=1
   for I in range(2,N):
     F[I]=F[I-1]+np.log(X)
     X=X+1.0
  I=AM+BM-CM+.10
  if I == 0:
    I1=AJ+BJ-CJ+1.10
  else:
    return 0.0
  if I1>0:
    I2=AJ-BJ+CJ+1.10
  else:
    return 0.0
  if I2>0.0:
    I3=BJ+CJ-AJ+1.10
  else:
    return 0.0
  if I3>0:
    X=AJ+BJ+CJ+2.10
  else:
    return 0.0
  I4=X
  I=X+.60
  I=I4-I
  if I==0:
    X=AJ+AM+1.10
  else:
    return 0.0
  I5=X
  if I5>0:
    I=X+.60
  else:
    return 0.0
  I=I-I5
  if I==0:
    I6=AJ-AM+1.10
  else:
    return 0.0
  if I6>0:
    X=BJ+BM+1.10
  else:
    return 0.0
  I7=X
  if I7>0:
    I=X+.60
  I=I-I7
  if I==0:
    I8=BJ-BM+1.10
  else:
    return 0.0
  if I8>0:
    X=CJ+CM+1.10
  I9=X
  if I9>0:
    I=X+.60
  I=I-I9
  if I == 0:
    I10=CJ-CM+1.10
  if I10>0:
    X=F(I1)+F(I2)+F(I3)-F(I4)
  I=I5-I6
  if I==0:
    I=I7-I8
    IF(I)18,200,18
    X=X+F(I5)+F(I6)+F(I7)+F(I8)+F(I9)+F(I10)
  else:
    X=X+F(I5)+F(I6)+F(I7)+F(I8)+F(I9)+F(I10)
X=X*.5D0
I10=MIN0(I1,I6,I7)
I2=I1-I5
I3=I1-I8
I9=MAX0(0,I2,I3)+1
I1=I1+1
I6=I6+1
I7=I7+1
I8=I9/2
E=1.D0
I5=I9*.5D0+.6D0
I8=I8-I5
IF(I8)20,21,20
21 E=-1.D0
20 S=0.D0
DO 22 I=I9,I10
C=X-F(I)-F(I1-I)-F(I6-I)-F(I7-I)-
*F(I-I2)-F(I-I3)
S=S+E*DEXP(C)
22 E=1.D0-E-1.D0
CLEBSCH=DSQRT(CJ+CJ+1.0)*S
RETURN
200 I=I4/2
I5=I4*.5D0+.6D0
I=I-I5
IF(I)100,201,100
201 I6=I5-I6+1
I7=I5-I8+1
I8=I5-I10+1
S=X*0.5D0+F(I5)-F(I6)-F(I7)-F(I8)
S=DEXP(S)
I5=I8/2
I6=I8*.5D0+.6D0
I5=I5-I6
IF(I5)202,203,202
203 S=1.D0-S-1.D0
202 CLEBSCH=S*DSQRT(CJ+CJ+1.D0)
RETURN
100 CLEBSCH=0.D0
RETURN
END

def ctormt(yor,yoi,l):
  #.................................................................ctormt
  #
  #---->    transformation from (ms) to real harmonics basis set
  #
  yor = np.zeros((l,l))
  yoi = np.zeros((l,l))
  sqtwo=1.0/np.sqrt(2.0)
  if l == 0 :
      yor[0,0]=1.0
  elif l==1:
      yoi[0,0]= sqtwo
      yoi[0,2]= sqtwo
      yor[1,1]=1.0
      yor[2,0]= sqtwo
      yor[2,2]=-sqtwo
  elif l== 2:
      yoi[0,0]= sqtwo
      yoi[0,4]=-sqtwo
      yoi[1,1]= sqtwo
      yoi[1,3]= sqtwo
      yor[2,2]=1.0
      yor[3,1]= sqtwo
      yor[3,3]=-sqtwo
      yor[4,0]= sqtwo
      yor[4,4]= sqtwo
  elif l == 3:
      yoi[0,0]= sqtwo
      yoi[0,6]= sqtwo
      yoi[1,1]= sqtwo
      yoi[1,5]=-sqtwo
      yoi[2,2]= sqtwo
      yoi[2,4]= sqtwo
      yor[3,3]=1.0
      yor[4,2]= sqtwo
      yor[4,4]=-sqtwo
      yor[5,1]= sqtwo
      yor[5,5]= sqtwo
      yor[6,0]= sqtwo
      yor[6,6]=-sqtwo



U = np.array([2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0])
xmu = np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
t = np.array([[0.0,1.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0],
              [1.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0],
              [0.0,1.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0],
              [1.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0],
              [1.0,0.0,0.0,0.0,0.0,1.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
              [0.0,1.0,0.0,0.0,1.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0],
              [0.0,0.0,1.0,0.0,0.0,1.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0],
              [0.0,0.0,0.0,1.0,1.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0],
              [0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,1.0,1.0,0.0,0.0,0.0],
              [0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0],
              [0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,1.0,0.0,0.0,1.0,0.0],
              [0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0],
              [1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,1.0],
              [0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,1.0,0.0],
              [0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,1.0],
              [0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,1.0,0.0]])

Ns = len(U)
sectors = np.array([[Ns/2,Ns/2],])

data = h5py.File("input.h5", "w");

beta = data.create_dataset("BETA", shape=(), dtype='f', data=10.0)

hop_g = data.create_group("sectors")
hop_g.create_dataset("values", data=sectors)


hop_g = data.create_group("hopping")
hop_g.create_dataset("values", data=t)

int_g = data.create_group("interaction")
int_ds = int_g.create_dataset("values", shape=(Ns,), data=U)

int_g = data.create_group("chemical_potential")
int_ds = int_g.create_dataset("values", shape=(Ns,), data=xmu)

