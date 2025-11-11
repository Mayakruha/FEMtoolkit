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
#===================================================================
#
#         Shift a hole
#
# Variables:
# D0 - shift at the start of the hole
# D1 - shift at the end of the hole
# Axis - unit vector of the hole axis
# SiftDir - unit vector for shift
# NodeSet - set of nodes on hole surface
#===================================================================
def ShiftSurf(mesh,D0,D1,Axis,ShiftDir,NodeSet):
    for i in range(len(mesh.point_sets[NodeSet])):
        Coords=mesh.points[mesh.point_sets[NodeSet][i]]
        AxisCoord=np.dot(Axis,Coords)
        if i==0:
            CoordMin=AxisCoord
            CoordMax=AxisCoord
        elif AxisCoord<CoordMin:
            CoordMin=AxisCoord
        elif AxisCoord>CoordMax:
            CoordMax=AxisCoord
    setNode=set(mesh.point_sets[NodeSet])
    MiddleNodes={}
    for Node in mesh.point_sets[NodeSet]: MiddleNodes[Node]=0
    for El in range(1,mesh.MaxElemNum+1):
        if mesh.cells[El]!=1:
            for i in range(len(FacesNodes[mesh.Eltype[El]])):
                Flag=True
                for Node_i in FacesNodes[mesh.Eltype[El]][i]:
                    if not mesh.cells[El][Node_i] in setNode: Flag=False
                if Flag:
                    for j in range(len(FacesNodes[mesh.Eltype[El]][i])):
                        if MiddleNodes[mesh.cells[El][FacesNodes[mesh.Eltype[El]][i][j]]]==0:                    
                            if mesh.Eltype[El]==6:
                                if j==0:
                                    if i==0 and not mesh.cells[El][7] in setNode: MiddleNodes[mesh.cells[El][0]]=mesh.cells[El][7]
                                    if i==1 and not mesh.cells[El][6] in setNode: MiddleNodes[mesh.cells[El][0]]=mesh.cells[El][6]
                                    if i==3 and not mesh.cells[El][4] in setNode: MiddleNodes[mesh.cells[El][0]]=mesh.cells[El][4]
                                if j==1:
                                    if i==0 and not mesh.cells[El][8] in setNode: MiddleNodes[mesh.cells[El][1]]=mesh.cells[El][8]
                                    if i==1 and not mesh.cells[El][5] in setNode: MiddleNodes[mesh.cells[El][1]]=mesh.cells[El][5]
                                    if i==2 and not mesh.cells[El][4] in setNode: MiddleNodes[mesh.cells[El][1]]=mesh.cells[El][4]
                                if j==2:
                                    if i==0 and not mesh.cells[El][9] in setNode: MiddleNodes[mesh.cells[El][2]]=mesh.cells[El][9]
                                    if i==2 and not mesh.cells[El][6] in setNode: MiddleNodes[mesh.cells[El][2]]=mesh.cells[El][6]
                                    if i==3 and not mesh.cells[El][5] in setNode: MiddleNodes[mesh.cells[El][2]]=mesh.cells[El][5]
                                if j==3:
                                    if i==1 and not mesh.cells[El][9] in setNode: MiddleNodes[mesh.cells[El][3]]=mesh.cells[El][9]
                                    if i==2 and not mesh.cells[El][7] in setNode: MiddleNodes[mesh.cells[El][3]]=mesh.cells[El][7]
                                    if i==3 and not mesh.cells[El][8] in setNode: MiddleNodes[mesh.cells[El][3]]=mesh.cells[El][8]
                            if mesh.Eltype[El]==7:                    
                                if j==0:
                                    if i==0 and not mesh.cells[El][12] in setNode: MiddleNodes[mesh.cells[El][0]]=mesh.cells[El][12]
                                    if i==2 and not mesh.cells[El][8] in setNode: MiddleNodes[mesh.cells[El][0]]=mesh.cells[El][8]
                                    if i==4 and not mesh.cells[El][6] in setNode: MiddleNodes[mesh.cells[El][0]]=mesh.cells[El][6]
                                if j==1:
                                    if i==0 and not mesh.cells[El][13] in setNode: MiddleNodes[mesh.cells[El][1]]=mesh.cells[El][13]
                                    if i==2 and not mesh.cells[El][7] in setNode: MiddleNodes[mesh.cells[El][1]]=mesh.cells[El][7]
                                    if i==3 and not mesh.cells[El][6] in setNode: MiddleNodes[mesh.cells[El][1]]=mesh.cells[El][6]
                                if j==2:
                                    if i==0 and not mesh.cells[El][14] in setNode: MiddleNodes[mesh.cells[El][2]]=mesh.cells[El][14]
                                    if i==3 and not mesh.cells[El][8] in setNode: MiddleNodes[mesh.cells[El][2]]=mesh.cells[El][8]
                                    if i==4 and not mesh.cells[El][7] in setNode: MiddleNodes[mesh.cells[El][2]]=mesh.cells[El][7]
                                if j==3:
                                    if i==1 and not mesh.cells[El][12] in setNode: MiddleNodes[mesh.cells[El][3]]=mesh.cells[El][12]
                                    if i==2 and not mesh.cells[El][11] in setNode: MiddleNodes[mesh.cells[El][3]]=mesh.cells[El][11]
                                    if i==4 and not mesh.cells[El][9] in setNode: MiddleNodes[mesh.cells[El][3]]=mesh.cells[El][9]
                                if j==4:
                                    if i==1 and not mesh.cells[El][13] in setNode: MiddleNodes[mesh.cells[El][4]]=mesh.cells[El][13]
                                    if i==2 and not mesh.cells[El][10] in setNode: MiddleNodes[mesh.cells[El][4]]=mesh.cells[El][10]
                                    if i==3 and not mesh.cells[El][9] in setNode: MiddleNodes[mesh.cells[El][4]]=mesh.cells[El][9]
                                if j==5:
                                    if i==1 and not mesh.cells[El][14] in setNode: MiddleNodes[mesh.cells[El][5]]=mesh.cells[El][14]
                                    if i==3 and not mesh.cells[El][11] in setNode: MiddleNodes[mesh.cells[El][5]]=mesh.cells[El][11]
                                    if i==4 and not mesh.cells[El][10] in setNode: MiddleNodes[mesh.cells[El][5]]=mesh.cells[El][10]
                            if mesh.Eltype[El]==8 or mesh.Eltype[El]==9:                    
                                if j==0:
                                    if i==0 and not mesh.cells[El][16] in setNode: MiddleNodes[mesh.cells[El][0]]=mesh.cells[El][16]
                                    if i==2 and not mesh.cells[El][11] in setNode: MiddleNodes[mesh.cells[El][0]]=mesh.cells[El][11]
                                    if i==5 and not mesh.cells[El][8] in setNode: MiddleNodes[mesh.cells[El][0]]=mesh.cells[El][8]
                                if j==1:
                                    if i==0 and not mesh.cells[El][17] in setNode: MiddleNodes[mesh.cells[El][1]]=mesh.cells[El][17]
                                    if i==2 and not mesh.cells[El][9] in setNode: MiddleNodes[mesh.cells[El][1]]=mesh.cells[El][9]
                                    if i==3 and not mesh.cells[El][8] in setNode: MiddleNodes[mesh.cells[El][1]]=mesh.cells[El][8]
                                if j==2:
                                    if i==0 and not mesh.cells[El][18] in setNode: MiddleNodes[mesh.cells[El][2]]=mesh.cells[El][18]
                                    if i==3 and not mesh.cells[El][10] in setNode: MiddleNodes[mesh.cells[El][2]]=mesh.cells[El][10]
                                    if i==4 and not mesh.cells[El][9] in setNode: MiddleNodes[mesh.cells[El][2]]=mesh.cells[El][9]
                                if j==3:
                                    if i==0 and not mesh.cells[El][19] in setNode: MiddleNodes[mesh.cells[El][3]]=mesh.cells[El][19]
                                    if i==4 and not mesh.cells[El][11] in setNode: MiddleNodes[mesh.cells[El][3]]=mesh.cells[El][11]
                                    if i==5 and not mesh.cells[El][10] in setNode: MiddleNodes[mesh.cells[El][3]]=mesh.cells[El][10]
                                if j==4:
                                    if i==1 and not mesh.cells[El][16] in setNode: MiddleNodes[mesh.cells[El][4]]=mesh.cells[El][16]
                                    if i==2 and not mesh.cells[El][15] in setNode: MiddleNodes[mesh.cells[El][4]]=mesh.cells[El][15]
                                    if i==5 and not mesh.cells[El][12] in setNode: MiddleNodes[mesh.cells[El][4]]=mesh.cells[El][12]
                                if j==5:
                                    if i==1 and not mesh.cells[El][17] in setNode: MiddleNodes[mesh.cells[El][5]]=mesh.cells[El][17]
                                    if i==2 and not mesh.cells[El][13] in setNode: MiddleNodes[mesh.cells[El][5]]=mesh.cells[El][13]
                                    if i==3 and not mesh.cells[El][12] in setNode: MiddleNodes[mesh.cells[El][5]]=mesh.cells[El][12]
                                if j==6:
                                    if i==1 and not mesh.cells[El][18] in setNode: MiddleNodes[mesh.cells[El][6]]=mesh.cells[El][18]
                                    if i==3 and not mesh.cells[El][14] in setNode: MiddleNodes[mesh.cells[El][6]]=mesh.cells[El][14]
                                    if i==4 and not mesh.cells[El][13] in setNode: MiddleNodes[mesh.cells[El][6]]=mesh.cells[El][13]
                                if j==7:
                                    if i==1 and not mesh.cells[El][19] in setNode: MiddleNodes[mesh.cells[El][7]]=mesh.cells[El][19]
                                    if i==4 and not mesh.cells[El][15] in setNode: MiddleNodes[mesh.cells[El][7]]=mesh.cells[El][15]
                                    if i==5 and not mesh.cells[El][14] in setNode: MiddleNodes[mesh.cells[El][7]]=mesh.cells[El][14]
    for Node in mesh.point_sets[NodeSet]:
        Coords=np.array(mesh.points[Node])
        AxisCoord=np.dot(Axis,Coords)
        Vect=ShiftDir*(D0+(D1-D0)*(AxisCoord-CoordMin)/(CoordMax-CoordMin))
        mesh.points[Node][:]+=Vect[:]
        if MiddleNodes[Node]!=0: mesh.points[MiddleNodes[Node]][:]+=0.5*Vect[:]
