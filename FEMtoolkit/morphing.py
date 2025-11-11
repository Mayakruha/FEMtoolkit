import numpy as np

#==========================================================================
# Changes radius by value of Mov inside range [R0, R1]. Radius axis: x-axis
#--------------------------------------------------------------------------
def radius_change(coord, R0=1.0, R1=2.1, Mov=0.01):
    Angle=atan(coord[i]/coord[2])
    R_init=(coord[1]**2+coord[2]**2)**0.5
    if R_init>R1:
        R=R_init+Mov
    else:
        R=R_init+Mov*(R_init-R0)/(R1-R0)
    return [coord[0], R*sin(Angle), R*cos(Angle)]
#--------------------------------------------------------------------------
# Moves along the vetor Mov between point0 and point1.
# Linear decrease of the move outside region between point0 and point1
# inside disctance 2*Mov from point0 and point1
#--------------------------------------------------------------------------
def move(coord, point0=(6.2,0.0,0.0), point1=(6.9,0.0,0.0), Mov=(-0.02,0.0,0.0)):
    NormVec=np.zeros(3)
    coord_new=coord.copy()
    for i in range(3): NormVec[i]=Mov[i]
    Disp=np.linalg.norm(NormVec)
    NormVec=NormVec/Disp
    S0=0
    S1=0
    for i in range(3):
        S0+=NormVec[i]*(coord[i]-point0[i])
        S1+=NormVec[i]*(coord[i]-point1[i])
    absS0=abs(S0)
    absS1=abs(S1)
    if S0*S1<0:
        for i in range(3): coord_new[i]+=Mov[i]
    elif (absS0<2*Disp)or(absS1<2*Disp):
        if absS0<absS1:
            for i in range(3): coord_new[i]+=(1-absS0/2/Disp)*Mov[i]
        else:
            for i in range(3): coord_new[i]+=(1-absS1/2/Disp)*Mov[i]
    return coord_new
