from FEMtoolkit import EstFaces, Normals
import numpy as np
from numpy import atan, sin, cos, pi
#===================================================================
#
#         Morphing (for Abaqus / Calculix)
#
# Variables:
# mesh        - mesh
# NodeSet     - set of nodes on hole surface
# func(coord) - function that returns new coordinates if Normal=False or displacement value if Normal=True 
# FreeNodeSet - set of nodes without any boundary conditions
#===================================================================
def morph(mesh, NodeSet, func, FreeNodeSet='', Normal=False, NLGEOM=True):
    MinNum=0
    MaxNum=0
    for cellblock in mesh.cell_data['Element_Ids']:
        for ElLabel in cellblock:
            if MinNum==0 or MinNum>ElLabel: MinNum=ElLabel
            if MaxNum<ElLabel: MaxNum=ElLabel
    Faces=EstFaces(mesh)
    if Normal: norms=Normals(mesh, NodeSet)
    print('WRITING A FILE FOR MORPHING...')
    OuterNodes=np.zeros(mesh.points.shape[0], dtype=np.int8)
    for ElType in Faces:
        for NodeMin in Faces[ElType]:
            for NodeMax in Faces[ElType][NodeMin]:
                for face in Faces[ElType][NodeMin][NodeMax]:
                    if face[0]==1:
                        for Num in list(face[1]):
                            OuterNodes[Num]=1
    MovedNodes={}
    if Normal:
        for i in range(len(mesh.point_sets[NodeSet])):
            Num=mesh.point_sets[NodeSet][i]
            Vect=norms[i]*func(mesh.points[Num])
            if np.linalg.norm(Vect)>0:
                MovedNodes[mesh.point_data['Node_Ids'][Num]]=Vect
                OuterNodes[Num]=0
    else:
        for Num in mesh.point_sets[NodeSet]:
            Vect=func(mesh.points[Num])-mesh.points[Num]
            if np.linalg.norm(Vect)>0:
                MovedNodes[mesh.point_data['Node_Ids'][Num]]=Vect
                OuterNodes[Num]=0
    if FreeNodeSet:
        for Num in mesh.point_sets[FreeNodeSet]:
            OuterNodes[Num]=0
    f=open('Run_Morphing.inp','w')
    f.write('** Mesh file **\n')
    f.write('*INCLUDE, INPUT=......\n')
    f.write('*MATERIAL, NAME=AUXETIC_MAT\n')
    f.write('*ELASTIC, TYPE=ISOTROPIC\n')
    f.write('1,-0.99\n')
    f.write('*ELSET, ELSET=EALL, generate\n')
    f.write(str(MinNum)+', '+str(MaxNum)+', 1\n')
    f.write('*SOLID SECTION, ELSET=EALL, MATERIAL=AUXETIC_MAT\n')
    f.write('*NSET, NSET=FIXED')
    NumsInLine=0
    for Num in range(mesh.points.shape[0]):
        if OuterNodes[Num]:
            if NumsInLine==0: f.write('\n')
            f.write(str(mesh.point_data['Node_Ids'][Num])+',')
            if NumsInLine==15:
                NumsInLine=0
            else:
                NumsInLine+=1
    f.write('\n')
    f.write('*BOUNDARY\n')
    f.write('FIXED, 1, 3\n')
    f.write('*STEP, NAME=MORPHING, NLGEOM=')
    if NLGEOM: f.write('YES\n*STATIC\n0.1,')
    else: f.write('NO\n*STATIC\n1.0,')
    f.write(' 1.0, 1e-05, 1\n')
    f.write('*BOUNDARY\n')
    for NodeLabel in MovedNodes:
        for i in range(3):
            f.write(str(NodeLabel)+', '+str(i+1)+', '+str(i+1)+', '+str(MovedNodes[NodeLabel][i])+'\n')
    f.write('*OUTPUT, FIELD, FREQUENCY=999\n')
    f.write('*NODE OUTPUT\n')
    f.write('U\n')
    f.write('*END STEP\n')
    f.close()
#==========================================================================
# Changes radius by value of Mov around axis X inside radius range [R0, R1]
# in the X range [X0, X1] with smooth transition dX, dR
#--------------------------------------------------------------------------
def radius_change(coord, R=[1.0, 2.0], X=[0.0, 1.0], dR=0.001, dX=0.002, Mov=0.01):
    Angle=atan(coord[1]/coord[2])
    R_init=(coord[1]**2+coord[2]**2)**0.5
    ScaleX=0
    if coord[0]>=X[0] and coord[0]<=X[1]:
        ScaleX=1
    else:
        S=min(abs(X[0]-coord[0]),abs(X[1]-coord[0]))/dX
        if S<1: ScaleX=1-S
    if R_init>=R[0] and R_init<=R[1]:
        R=R_init+ScaleX*Mov
    else:
        S=min(abs(R[0]-R_init),abs(R[0]-R_init))/dR
        if S<1:
            R=R_init+ScaleX*(1-S)*Mov
        else:
            R=R_init
    return [coord[0], R*sin(Angle), R*cos(Angle)]
#==========================================================================
# Changes cylidrical fillet radius by value of dR. Fillet axis position X0, R0
#--------------------------------------------------------------------------
def torus_rad_change(coord, X0=6.816, R0=1.184, dR=-0.002):
    Angle=atan(coord[1]/coord[2])
    R_init(coord[1]**2+coord[2]**2)**0.5
    if coord[0]==X0:
        if R_init-R0>0: Torus_ang=pi/2
        else: Torus_ang=-pi/2
    elif coord[0]-X0>0: Torus_ang=atan((R_init-R0)/(coord[0]-X0))
    else: Torus_ang=atan((R_init-R0)/(coord[0]-X0))+pi
    R_tor=((R0-R_init)**2+(X0-coord[0])**2)**0.5
    R_new=R0+dR+(R_tor+dR)*sin(Torus_ang)
    return [X0+dR+(R_tor+dR)*cos(Torus_ang), R_new*sin(Angle), R_new*cos(Angle)]
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
