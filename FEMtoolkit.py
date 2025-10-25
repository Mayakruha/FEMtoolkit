import numpy as np
from meshio import read, Mesh, CellBlock
import vtk
FacesNodes={'triangle':((0,1),(1,2),(2,0)),'quad':((0,1),(1,2),(2,3),(3,0)),'triangle6':((0,1,3),(1,2,4),(2,0,5)),\
'quad8':((0,1,4),(1,2,5),(2,3,6),(3,0,7)),'tetra':((0,1,2),(0,3,1),(1,3,2),(2,3,0)),'wedge':((0,1,2),(3,5,4),(0,3,4,1),(1,4,5,2),(2,5,3,0)),\
'hexahedron':((0,1,2,3),(4,7,6,5),(0,4,5,1),(1,5,6,2),(2,6,7,3),(3,7,4,0)),'tetra10':((0,1,2,4,5,6),(0,3,1,7,8,4),(1,3,2,8,9,5),(2,3,0,9,7,6)),\
'wedge15':((0,1,2,6,7,8),(3,5,4,9,10,11),(0,3,4,1,12,9,13,6),(1,4,5,2,13,10,14,7),(2,5,3,0,14,11,12,8)),\
'hexahedron20':((0,1,2,3,8,9,10,11),(4,7,6,5,12,13,14,15),(0,4,5,1,16,12,17,8),(1,5,6,2,17,13,18,9),(2,6,7,3,18,14,19,10),(3,7,4,0,19,15,16,11))}
#---------------
NodesForNormals={'tetra10':(((6,4),(4,5),(5,6),(6,5),(4,6),(5,4)),((4,7),(7,8),(8,4),(4,8),(7,4),(8,7)),((5,8),(8,9),(9,5),(5,9),(8,5),(9,8)),\
((6,9),(9,7),(7,6),(6,7),(9,6),(7,9)))}
#------GENERAL FUNCTIONS------------------
def NormToTri(Nodes):
    X=(Nodes[1][1]-Nodes[1][0])*(Nodes[2][2]-Nodes[2][0])-(Nodes[1][2]-Nodes[1][0])*(Nodes[2][1]-Nodes[2][0])
    Y=(Nodes[0][2]-Nodes[0][0])*(Nodes[2][1]-Nodes[2][0])-(Nodes[0][1]-Nodes[0][0])*(Nodes[2][2]-Nodes[2][0])
    Z=(Nodes[0][1]-Nodes[0][0])*(Nodes[1][2]-Nodes[1][0])-(Nodes[0][2]-Nodes[0][0])*(Nodes[1][1]-Nodes[1][0])
    L=np.linalg.norm((X,Y,Z))
    return (X/L, Y/L, Z/L)
def normalize_v3(arr):
	lens = np.sqrt( arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 )     
	arr[:,0] /= lens     
	arr[:,1] /= lens     
	arr[:,2] /= lens
	return arr
#===================================================================
# Collect info about Faces
# Faces={} # the firts key - element type; the second key - min mumber of nodes; the second key - max mumber of nodes
#          # list (number of faces, set of numbers of nodes)
#===================================================================
def EstFaces(mesh):
    print('ANALYSING FACES...')
    Faces={}
    for cell_block in mesh.cells:
        if cell_block.type in FacesNodes:
            if not cell_block.type in Faces: Faces[cell_block.type]={}
            for El in cell_block.data:
                for Indx in FacesNodes[cell_block.type]:
                    Flag=True
                    Nodes=set()
                    for i in Indx: Nodes.add(El[i])
                    minNode=min(Nodes)
                    maxNode=max(Nodes)
                    if minNode in Faces[cell_block.type]:
                        if maxNode in Faces[cell_block.type][minNode]:
                            for Fc in Faces[cell_block.type][minNode][maxNode]:
                                if Fc[1]==Nodes:
                                    Fc[0]+=1
                                    Flag=False
                        else: Faces[cell_block.type][minNode][maxNode]=[]
                    else:
                        Faces[cell_block.type][minNode]={}
                        Faces[cell_block.type][minNode][maxNode]=[]
                    if Flag: Faces[cell_block.type][minNode][maxNode].append([1,Nodes])
    return Faces
#===================================================================
#         Return normals in nodes
# NSet - name of a node set 
#===================================================================
def Normals(mesh, NSet):
    print('ASSESING NORMALS...')
    N=len(mesh.point_sets[NSet])
    norm=np.zeros((N,3))
    List=np.full(mesh.points.shape[0],-1,dtype=np.int32)
    for i in range(N):
        Node=mesh.point_sets[NSet][i]
        List[Node]=i
    for block in mesh.cells:
        for Elem in block.data:
            for i in range(len(FacesNodes[block.type])):
                Flag=True
                for NdIndx in FacesNodes[block.type][i]:
                    if List[Elem[NdIndx]]==-1:
                        Flag=False
                        break
                if Flag:
                    for j in range(len(FacesNodes[block.type][i])):
                        Node=Elem[FacesNodes[block.type][i][j]]
                        Vc1=mesh.points[Elem[NodesForNormals[block.type][i][j][0]]]-mesh.points[Node]
                        Vc2=mesh.points[Elem[NodesForNormals[block.type][i][j][1]]]-mesh.points[Node]
                        norm[List[Node]]+=np.cross(Vc1,Vc2)
    for i range(N):
        length=np.linalg.norm(norm[i])
        if length!=0:
            norm[i][:]=norm[i][:]/length
    return norm
#===================================================================
#         import Point/Node data
#===================================================================
def import_ndload(mesh,FileName,LoadName):
    NodeValue={}
    f=open(FileName,'r')
    txt=f.readline()
    while txt:
        Values=txt.split(',')
        Node=int(Values[0])
        NodeValue[Node]=float(Values[1])
        txt=f.readline()      
    f.close()
    size=len(mesh.points)
    mesh.point_data[LoadName]=np.zeros(size)
    if 'Node_Num' in mesh.point_data:
        for i in range(size):
            mesh.point_data[LoadName][i]=NodeValue[mesh.point_data['Node_Num'][i]]
    else:
        Values=list(NodeValue.values())
        for i in range(min(size,len(Values))):
            mesh.point_data[LoadName][i]=Values[i]        
#===================================================================
#         import Point/Node data from 2ndFlow
#===================================================================
def import_ndload2ndFlow(mesh,FileName,LoadName):
    NodeValue={}
    f=open(FileName,'r')
    txt=f.readline()
    while txt:
        if 'Time' in 'txt':
            LoadNm=LoadName+'_'+txt[txt.find(':')+3:-1].zfill(5)
            NodeValue[LoadNm]={}
            txt=f.readline()
            txt=f.readline()
            txt=f.readline()
        Values=txt.split(',')
        if len(Values)==5:
            Node=int(Values[2])
            NodeValue[LoadNm][Node]=float(Values[4])
        txt=f.readline()
    f.close()
    size=len(mesh.points)
    for LoadNm in NodeValue:
        mesh.point_data[LoadNm]=np.zeros(size)
        if 'Node_Num' in mesh.point_data:
            for i in range(size):
                mesh.point_data[LoadNm][i]=NodeValue[LoadNm][mesh.point_data['Node_Num'][i]]
        else:
            Values=list(NodeValue[LoadNm].values())
            for i in range(min(size,len(Values))):
                mesh.point_data[LoadNm][i]=Values[i]               
#===================================================================
#         export Point/Node data
#===================================================================
def export_ndload(mesh,FileName,LoadName,separator=','):
    f=open(FileName,'w')
    List=[]
    for Name in mesh.point_data:
        if LoadName in Name:
            List.append(Name)
    Flag=len(List)>1
    size=len(mesh.points)
    for Name in List:
        if Flag:
            f.write('Time='+Name[Name.rfind('_')+1:]+'\n')
        if 'Node_Num' in mesh.point_data:
            for i in range(size):
                f.write(str(mesh.point_data['Node_Num'][i])+separator+str(mesh.point_data[Name][i])+'\n')
        else:
            for i in range(size):
                f.write(str(i+1)+separator+str(mesh.point_data[Name][i])+'\n')
    f.close()
#===================================================================
#         change Point/Node data
#===================================================================
def change_ndload(mesh, LoadName, dValue, Raff, Point, Vect, ChangeType):
    size=len(mesh.points)
    if ChangeType=='PLANE':            
        VecLen=(Vect[0]**2+Vect[1]**2+Vect[2]**2)**0.5
        for j in range(size):
            Sum=0
            for i in range(3):
                Sum+=(mesh.points[j][i]-Point[i])*Vect[i]/VecLen
            R=abs(Sum)
            if R<Raff:
                mesh.point_data[LoadName][j]+=dValue*(1-(R/Raff)**2)
    elif ChangeType=='CYL':
        Raff=Raff**2
        for j in range(size):
            Sum=0
            for i in range(3):
                Sum+=Vect[i]*(mesh.points[j][i]-Point[i])
            Sum/=Vect[0]+Vect[1]+Vect[2]
            R=0
            for i in range(3):
                R+=(mesh.points[j][i]-Point[i]-Vect[i]*Sum)**2
            if R<Raff:
                mesh.point_data[LoadName][j]+=dValue*(1-R/Raff)
#===================================================================
#         Replace quadratic elements by linear elements
#===================================================================
def Make3DLinearMesh(mesh):
    cells_triang=[]
    cells_tetr=[]
    tri_oldnums=[]
    tetr_oldnums=[]
    Elems={}
    for CellBlock in mesh.cells:
        Elems[CellBlock.type]={}
        for i in range(len(CellBlock.data)):
            Nodelist=CellBlock.data[i]
            if CellBlock.type=='triangle':
                cells_triang.append(Nodelist)
                tri_oldnums.append((CellBlock.type,i))
                Elems_tri[CellBlock.type][i]=[]
            elif CellBlock.type=='quad':
                cells_triang.append([Nodelist[0],Nodelist[1],Nodelist[2]])
                cells_triang.append([Nodelist[2],Nodelist[3],Nodelist[0]])
                for j in range(2):tri_oldnums.append((CellBlock.type,i))
                Elems_tri[CellBlock.type][i]=[]
            elif CellBlock.type=='triangle6':
                cells_triang.append([Nodelist[0],Nodelist[3],Nodelist[5]])
                cells_triang.append([Nodelist[3],Nodelist[1],Nodelist[4]])
                cells_triang.append([Nodelist[4],Nodelist[2],Nodelist[5]])
                cells_triang.append([Nodelist[3],Nodelist[4],Nodelist[5]])
                for j in range(4):tri_oldnums.append((CellBlock.type,i))
                Elems_tri[CellBlock.type][i]=[]
            elif CellBlock.type=='quad8':
                cells_triang.append([Nodelist[7],Nodelist[0],Nodelist[4]])
                cells_triang.append([Nodelist[4],Nodelist[1],Nodelist[5]])
                cells_triang.append([Nodelist[5],Nodelist[2],Nodelist[6]])
                cells_triang.append([Nodelist[6],Nodelist[3],Nodelist[7]])
                cells_triang.append([Nodelist[4],Nodelist[6],Nodelist[7]])
                cells_triang.append([Nodelist[4],Nodelist[5],Nodelist[6]])
                for j in range(6):tri_oldnums.append((CellBlock.type,i))
                Elems_tri[CellBlock.type][i]=[]
            elif CellBlock.type=='tetra':
                cells_tetr.append(Nodelist)
                tetr_oldnums.append((CellBlock.type,i))
                Elems_tet[CellBlock.type][i]=[]
            elif CellBlock.type=='wedge':
                cells_tetr.append([Nodelist[0],Nodelist[1],Nodelist[3],Nodelist[2]])
                cells_tetr.append([Nodelist[1],Nodelist[4],Nodelist[3],Nodelist[2]])
                cells_tetr.append([Nodelist[3],Nodelist[2],Nodelist[4],Nodelist[5]])
                for j in range(3):tetr_oldnums.append((CellBlock.type,i))
                Elems_tet[CellBlock.type][i]=[]
            elif CellBlock.type=='hexahedron':
                cells_tetr.append([Nodelist[0],Nodelist[1],Nodelist[3],Nodelist[4]])
                cells_tetr.append([Nodelist[1],Nodelist[2],Nodelist[3],Nodelist[4]])
                cells_tetr.append([Nodelist[3],Nodelist[4],Nodelist[2],Nodelist[7]])
                cells_tetr.append([Nodelist[5],Nodelist[4],Nodelist[6],Nodelist[1]])
                cells_tetr.append([Nodelist[4],Nodelist[7],Nodelist[6],Nodelist[1]])
                cells_tetr.append([Nodelist[6],Nodelist[1],Nodelist[7],Nodelist[2]])
                for j in range(6):tetr_oldnums.append((CellBlock.type,i))
                Elems_tet[CellBlock.type][i]=[]
            elif CellBlock.type=='tetra10':
                cells_tetr.append([Nodelist[0],Nodelist[4],Nodelist[6],Nodelist[7]])
                cells_tetr.append([Nodelist[4],Nodelist[1],Nodelist[5],Nodelist[8]])
                cells_tetr.append([Nodelist[5],Nodelist[2],Nodelist[6],Nodelist[9]))
                cells_tetr.append([Nodelist[7],Nodelist[8],Nodelist[9],Nodelist[3]))                
                cells_tetr.append([Nodelist[6],Nodelist[4],Nodelist[5],Nodelist[7]))
                cells_tetr.append([Nodelist[4],Nodelist[8],Nodelist[5],Nodelist[7]))
                cells_tetr.append([Nodelist[5],Nodelist[8],Nodelist[9],Nodelist[7]))
                cells_tetr.append([Nodelist[5],Nodelist[9],Nodelist[6],Nodelist[7]))
                for j in range(8):tetr_oldnums.append((CellBlock.type,i))
                Elems_tet[CellBlock.type][i]=[]
            elif CellBlock.type=='wedge15':
                cells_tetr.append([Nodelist[0],Nodelist[6],Nodelist[12],Nodelist[8]])
                cells_tetr.append([Nodelist[1],Nodelist[13],Nodelist[6],Nodelist[7]])
                cells_tetr.append([Nodelist[3],Nodelist[12],Nodelist[9],Nodelist[11]])
                cells_tetr.append([Nodelist[4],Nodelist[9],Nodelist[13],Nodelist[10]])                
                cells_tetr.append([Nodelist[11],Nodelist[12],Nodelist[9],Nodelist[5]])
                cells_tetr.append([Nodelist[10],Nodelist[9],Nodelist[13],Nodelist[5]])
                cells_tetr.append([Nodelist[12],Nodelist[8],Nodelist[6],Nodelist[2]])
                cells_tetr.append([Nodelist[6],Nodelist[7],Nodelist[13],Nodelist[2]])                
                cells_tetr.append([Nodelist[9],Nodelist[12],Nodelist[13],Nodelist[5]])
                cells_tetr.append([Nodelist[12],Nodelist[6],Nodelist[13],Nodelist[2]])
                cells_tetr.append([Nodelist[5],Nodelist[2],Nodelist[12],Nodelist[13]])
                for j in range(11):tetr_oldnums.append((CellBlock.type,i))
                Elems_tet[CellBlock.type][i]=[]
            elif CellBlock.type=='hexahedron20':
                cells_tetr.append([Nodelist[0],Nodelist[8],Nodelist[11],Nodelist[16]])
                cells_tetr.append([Nodelist[1],Nodelist[9],Nodelist[8],Nodelist[17]])
                cells_tetr.append([Nodelist[2],Nodelist[10],Nodelist[9],Nodelist[18]])
                cells_tetr.append([Nodelist[3],Nodelist[11],Nodelist[10],Nodelist[19]])               
                cells_tetr.append([Nodelist[4],Nodelist[15],Nodelist[12],Nodelist[16]])
                cells_tetr.append([Nodelist[5],Nodelist[12],Nodelist[13],Nodelist[17]])
                cells_tetr.append([Nodelist[6],Nodelist[13],Nodelist[14],Nodelist[18]])
                cells_tetr.append([Nodelist[7],Nodelist[15],Nodelist[14],Nodelist[19]])
                cells_tetr.append([Nodelist[11],Nodelist[12],Nodelist[15],Nodelist[16]])
                cells_tetr.append([Nodelist[11],Nodelist[8],Nodelist[12],Nodelist[16]])
                cells_tetr.append([Nodelist[12],Nodelist[8],Nodelist[9],Nodelist[17]])
                cells_tetr.append([Nodelist[9],Nodelist[13],Nodelist[12],Nodelist[17]])                
                cells_tetr.append([Nodelist[9],Nodelist[10],Nodelist[13],Nodelist[18]])
                cells_tetr.append([Nodelist[10],Nodelist[14],Nodelist[13],Nodelist[18]])
                cells_tetr.append([Nodelist[10],Nodelist[14],Nodelist[11],Nodelist[19]])
                cells_tetr.append([Nodelist[14],Nodelist[15],Nodelist[11],Nodelist[19]])               
                cells_tetr.append([Nodelist[11],Nodelist[12],Nodelist[14],Nodelist[15]])
                cells_tetr.append([Nodelist[14],Nodelist[12],Nodelist[9],Nodelist[13]])
                cells_tetr.append([Nodelist[11],Nodelist[8],Nodelist[9],Nodelist[12]])
                cells_tetr.append([Nodelist[11],Nodelist[9],Nodelist[10],Nodelist[14]])
                cells_tetr.append([Nodelist[14],Nodelist[12],Nodelist[11],Nodelist[9]])
                for j in range(21):tetr_oldnums.append((CellBlock.type,i))
                Elems_tet[CellBlock.type][i]=[]
    #------CELLS
    cells=[]
    TriNum=len(cells_triang)
    TetNum=len(cells_tetr)
    if TriNum>0:
        cells.append(CellBlock('triangle',np.array(cells_triang)))
    if TetNum>0:
        cells.append(CellBlock('tetra',np.array(cells_tetr)))
    #------CELLS_DATA
    cell_data={}
    for Name in [Name for Name in mesh.cell_data.keys() if Name != 'Elem_Num']:
        cell_data[Name]=[[],[]]
        for OldEl in tri_oldnums: cell_data[Name][0].append(mesh.cell_data_dict[Name][OldEl[0]][OldEl[1]])
        for OldEl in tetr_oldnums: cell_data[Name][1].append(mesh.cell_data_dict[Name][OldEl[0]][OldEl[1]])
    #------CELL_SETS
    for i in range(TriNum):
        Elems_tri[tri_oldnums[i][0]][tri_oldnums[i][1]].append(i)
    for i in range(TetNum):
        Elems_tet[tetr_oldnums[i][0]][tetr_oldnums[i][1]].append(i)
    cell_sets={}
    for Name in mesh.cell_sets_dict:
        cell_sets[Name]=[[],[]]
        for ElType in mesh.cell_sets_dict[Name]:
            for Num in mesh.cell_sets_dict[Name][ElType]:
                cell_sets[Name][0]+=Elems_tri[ElType][Num]
                cell_sets[Name][1]+=Elems_tet[ElType][Num]
    return Mesh(mesh.points.copy(), cells, point_data=mesh.point_data.copy(), cell_data=cell_data, point_sets=mesh.point_sets.copy(), cell_sets=cell_sets)
#===================================================================
#         import Face load
# LoadType: 'P' - Pressure; 'S' - Heat Flux; 'F' - HTC
#===================================================================
def import_fcload(mesh,FileName,LoadType):
    face_data[LoadType]={}
    f=open(FileName,'r')
    txt=f.readline()
    while txt:
        Values=txt.split(',')
        El=int(Values[0])
        face_data[LoadType][int(Values[0])]=[]
        Val=[int(Values[1][-1:]),float(Values[2])]
        if len(Values)>3:Val.append(float(Values[3]))
        face_data[LoadType][int(Values[0])].append(Val)
        txt=f.readline()        
    f.close()
    if not LoadType in mesh.face_data: mesh.face_data[LoadType]={}
    if 'Elem_Num' in mesh.cells_data_dict:
        for ElType in mesh.cells_data_dict['Elem_Num']:
            size=len(mesh.cells_data_dict['Elem_Num'][ElType])
            if not ElType in mesh.face_data[LoadType]: mesh.face_data[LoadType][ElType]=np.zeros((size,len(FacesNodes[ElType],2)))
            for i in range(0,size):
                if mesh.cells_data_dict['Elem_Num'][ElType][i] in face_data[LoadType]:
                    for data in face_data[LoadType][mesh.cells_data_dict['Elem_Num'][ElType][i]]:
                        for j in range(len(data)-1):
                            mesh.face_data[LoadType][ElType][i][data[0]-1][j]=data[j+1]
    else:
        i0=0
        for ElType in mesh.cells_dict:
            size=len(mesh.cells_dict[ElType])
            if not ElType in mesh.face_data[LoadType]: mesh.face_data[LoadType][ElType]=np.zeros((size,len(FacesNodes[ElType],2)))
            for i in range(0,size):
                if i+i0+1 in face_data[LoadType]:
                    for data in face_data[LoadType][i+i0+1]:
                        for j in range(len(data)-1):
                            mesh.face_data[LoadType][ElType][i][data[0]-1][j]=data[j+1]
            i0+=size
#===================================================================
#         export Face load
#===================================================================
def export_fcload(mesh,FileName,LoadName):
    f=open(FileName,'w')
    if 'Elem_Num' in mesh.cells_data_dict:
        for ElType in mesh.cells_data_dict['Elem_Num']:
            for i in range(0,len(mesh.cells_data_dict['Elem_Num'][ElType])):
                for fc_num in range(len(FacesNodes[ElType])):
                    if mesh.face_data[LoadName][ElType][i][fc_num][0]!=0:
                        f.write(str(mesh.cells_data_dict['Elem_Num'][ElType][i])+', '+LoadName+str(fc_num+1))
                        for Val in mesh.face_data[LoadName][ElType][i][fc_num]:
                            if Val!=0:f.write(', '+str(Val))
                        f.write('\n')
    else:
        i0=0
        for ElType in mesh.cells_dict:
            size=len(mesh.cells_dict[ElType])
            for i in range(0,size):
                for fc_num in range(len(FacesNodes[ElType])):
                    if mesh.face_data[LoadName][ElType][i][fc_num][0]!=0:
                        f.write(str(i0+i+1)+', '+LoadName+str(fc_num+1))
                        for Val in mesh.face_data[LoadName][ElType][i][fc_num]:
                            if Val!=0:f.write(', '+str(Val))
                        f.write('\n')
            i0+=size
    f.close()
#===================================================================
# Creates mesh for Face load (for mapping)
# Variables:
# mesh - Original mesh
# Surf - Name of a surface
#===================================================================
def MeshFromFaceLoad(mesh, Surf):
    size=len(mesh.points)
    Nums=np.full(size,size,dtype=np.int32)
    NumsEl={}
    Points=[]
    Cells=[]
    Node_count=0
    Elm_count=0
    for Face in mesh.surfaces[Surf]:
        for ElType in mesh.cell_sets_dict[Face[0]]:
            NumsEl[ElType]={}
            for ElemNum in mesh.cell_sets_dict[Face[0]][ElType]:            
                Flag=False
                for LoadName in mesh.face_data:
                    if not ElType in mesh.face_data[LoadName]: break 
                    if mesh.face_data[LoadName][ElType][ElemNum][Face[1]]!=0:
                        Flag=True
                if Flag:
                    if not Face[1] in NumsEl[ElType]: NumsEl[ElType][Face[1]]={}
                    Nodelist=[]
                    for i in FacesNodes[ElType][Face[1]]:
                        Node=mesh.cells_dict[ElType][ElemNum][i]
                        if Nums[Node]==size:
                            Nums[Node]=Node_count
                            Points.append(mesh.points[Node])
                            Node_count+=1
                        Nodelist.append(Nums[Node])
                    if ElType=='tetra':
                        Cells.append(Nodelist)
                        NumsEl[ElType][Face[1]][ElemNum]=[Elm_count,]
                        Elm_count+=1
                    elif ElType=='tetra10':
                        Cells.append([Nodelist[0],Nodelist[3],Nodelist[5]])
                        Cells.append([Nodelist[1],Nodelist[4],Nodelist[3]])
                        Cells.append([Nodelist[3],Nodelist[4],Nodelist[5]])
                        Cells.append([Nodelist[2],Nodelist[5],Nodelist[4]])
                        NumsEl[ElType][Face[1]][ElemNum]=[Elm_count,Elm_count+1,Elm_count+2,Elm_count+3]
                        Elm_count+=4
    #-----------Field---------------------------------------
    cell_data={}
    for LoadName in mesh.face_data:
        cell_data[LoadName]=[np.zeros(Elm_count),]
        for ElType in mesh.face_data[LoadName]:
            for ElemNum in range(len(mesh.face_data[LoadName][ElType])):
                for Face in range(len(mesh.face_data[LoadName][ElType][ElemNum])):
                    for j in NumsEl[ElType][Face][ElemNum]:
                        cell_data[LoadName][0][j]=mesh.face_data[LoadName][ElType][ElemNum][0]
    return Mesh(Points(), [CellBlock('triangle',np.array(Cells)),], cell_data=cell_data)
#===================================================================
#         Node set -> Surface
#===================================================================
def NodeIntoSurf(mesh,NSet):
    self.surfaces[NSet]=[]
    for ElType in mesh.cells_dict:
        for i in range(len(mesh.cells_dict[ElType]):
            for FaceIndx in range(len(FacesNodes[ElType])):
                Flag=True
                for NdIndx in FacesNodes[ElType][FaceIndx]:
                    if not mesh.cells_dict[ElType][i][NdIndx] in mesh.points_sets[NSet]: Flag=False
                if Flag:
                    SetFaceName=NSet+'_S'+str(FaceIndx+1)
                    if not (SetFaceName,FaceIndx) in mesh.surfaces[NSet]: mesh.surfaces[NSet].append((SetFaceName,FaceIndx))
                    if not SetFaceName in mesh.cell_sets_dict: mesh.cell_sets_dict[SetFaceName]={}
                    if not ElType in mesh.cell_sets_dict[SetFaceName]:mesh.cell_sets_dict[SetFaceName][ElType]=[]
                    mesh.cell_sets_dict[SetFaceName][ElType].append(i)
#===================================================================
#
#         Extract thickness of coating that is simulated by linear triangular prism
#
# Variables:
# NsetName  - Name of a set of nodes on internal surface (surface betwen base material and coating)
# EsetNames - List of Names of a set of elements for layers of coating
#===================================================================
def ExtractCoating(mesh, NsetName, EsetNames):
    points=[]
    cells=[]
    ElemRef=[]
    StackElem={}
    for i in range(1,len(EsetNames)):
        StackElem[EsetNames[i]]=[]
        for jj in range(len(mesh.cell_sets[EsetNames[i]])):
            StackElem[EsetNames[i]].append(set(mesh.cell_sets[EsetNames[i]][jj]))
    i=0
    for jj in range(len(mesh.cell_sets[EsetNames[0]])):
        ElemRef.append({})
        for ElemNum in mesh.cell_sets[EsetNames[0]][jj]:
            ElemExist=False
            if mesh.cells[jj].data[ElemNum][0] in mesh.point_sets[NsetName] and mesh.cells[jj].data[ElemNum][1] in mesh.point_sets[NsetName] and mesh.cells[jj].data[ElemNum][2] in mesh.point_sets[NsetName]:
                for j in range(3): points.append(mesh.points[mesh.cells[jj].data[ElemNum][j]])
                ElemRef[jj][ElemNum]=(i*3,i*3+1,i*3+2,0)
                Nodes0=(mesh.cells[jj].data[ElemNum][3], mesh.cells[jj].data[ElemNum]4], mesh.cells[jj].data[ElemNum][5])
                ElemExist=True
            elif mesh.cells[jj].data[ElemNum][3] in mesh.point_sets[NsetName] and mesh.cells[jj].data[ElemNum][4] in mesh.point_sets[NsetName] and mesh.cells[jj].data[ElemNum][5] in mesh.point_sets[NsetName]:
                for j in range(3): points.append(mesh.points[mesh.cells[jj].data[ElemNum][j+3]])
                ElemRef[jj][ElemNum]=(i*3,i*3+1,i*3+2,1)
                Nodes0=(mesh.cells[jj].data[ElemNum][0],mesh.cells[jj].data[ElemNum][1],mesh.cells[jj].data[ElemNum][2])
                ElemExist=True
            else:
                print('The nodes haven\'t been found in '+NsetName+' for element '+str(ElemNum))
            for EName in StackElem:
                for ElNum in list(StackElem[EName][jj]):
                    Nodes=set(mesh.cells[jj].data[ElNum])
                    if len(Nodes.intersection(set(Nodes0)))==3:
                        if mesh.cells[jj].data[ElNum][0]==Nodes0[0]:
                            ElemRef[jj][ElNum]=(i*3,i*3+1,i*3+2,0)
                            Nodes0=(mesh.cells[jj].data[ElemNum][3],mesh.cells[jj].data[ElemNum][4],mesh.cells[jj].data[ElemNum][5])
                        elif mesh.cells[jj].data[ElNum][3]==Nodes0[0]:
                            ElemRef[jj][ElNum]=(i*3,i*3+1,i*3+2,1)
                            Nodes0=(mesh.cells[jj].data[ElemNum][0],mesh.cells[jj].data[ElemNum][1],mesh.cells[jj].data[ElemNum][2])
                        elif mesh.cells[jj].data[ElNum][1]==Nodes0[0]:
                            ElemRef[jj][ElNum]=(i*3+1,i*3+2,i*3,0)
                            Nodes0=(mesh.cells[jj].data[ElemNum][4],mesh.cells[jj].data[ElemNum][5],mesh.cells[jj].data[ElemNum][3])
                        elif mesh.cells[jj].data[ElNum][4]==Nodes0[0]:
                            ElemRef[jj][ElNum]=(i*3+1,i*3+2,i*3,1)
                            Nodes0=(mesh.cells[jj].data[ElemNum][1],mesh.cells[jj].data[ElemNum][2],mesh.cells[jj].data[ElemNum][0])
                        elif mesh.cells[jj].data[ElNum][2]==Nodes0[0]:
                            ElemRef[jj][ElNum]=(i*3+2,i*3,i*3+1,0)
                            Nodes0=(mesh.cells[jj].data[ElemNum][5],mesh.cells[jj].data[ElemNum][3],mesh.cells[jj].data[ElemNum][2])
                        elif mesh.cells[jj].data[ElNum][5]==Nodes0[0]:
                            ElemRef[jj][ElNum]=(i*3+2,i*3,i*3+1,1)
                            Nodes0=(mesh.cells[jj].data[ElemNum][2],mesh.cells[jj].data[ElemNum][0],mesh.cells[jj].data[ElemNum][1])
                        StackElem[EName][jj].remove(ElNum)
                        break
            if ElemExist:
                cells.append((i*3,i*3+1,i*3+2))                
                i+=1
    #============= thickness analysis
    Thick={}
    for EName in EsetNames:
        Thick[EName]=np.zeros(3*i)
            for jj in range(len(mesh.cell_sets[EName])):
                for ENum in mesh.cell_sets[EName][jj]:
                    if ElemRef[jj][ENum]!=None:
                        Nodes=mesh.cells[jj].data[ENum]
                        Vb1=np.array((mesh.points[Nodes[1]][0]-mesh.points[Nodes[0]][0],mesh.points[Nodes[1]][1]-mesh.points[Nodes[0]][1],mesh.points[Nodes[1]][2]-mesh.points[Nodes[0]][2]))
                        Vb2=np.array((mesh.points[Nodes[2]][0]-mesh.points[Nodes[0]][0],mesh.points[Nodes[2]][1]-mesh.points[Nodes[0]][1],mesh.points[Nodes[2]][2]-mesh.points[Nodes[0]][2]))
                        NormB=np.cross(Vb1, Vb2)
                        Vt1=np.array((mesh.points[Nodes[4]][0]-mesh.points[Nodes[3]][0],mesh.points[Nodes[4]][1]-mesh.points[Nodes[3]][1],mesh.points[Nodes[4]][2]-mesh.points[Nodes[3]][2]))
                        Vt2=np.array((mesh.points[Nodes[5]][0]-mesh.points[Nodes[3]][0],mesh.points[Nodes[5]][1]-mesh.points[Nodes[3]][1],mesh.points[Nodes[5]][2]-mesh.points[Nodes[3]][2]))
                        NormT=np.cross(Vt1, Vt2)
                        if ElemRef[jj][ENum][3]==0:
                            NormB=NormB/np.linalg.norm(NormB)
                            for j in range(0,3):
                                Vec=np.array((mesh.points[Nodes[3]][0]-mesh.points[Nodes[j]][0],mesh.points[Nodes[3]][1]-mesh.points[Nodes[j]][1],mesh.points[Nodes[3]][2]-mesh.points[Nodes[j]][2]))
                                Dist=abs(np.dot(NormT,Vec)/np.dot(NormB,NormT))
                                Thick[EName][ElemRef[jj][ENum][j]]=Dist
                        elif ElemRef[jj][ENum][3]==1:
                            NormT=NormT/np.linalg.norm(NormT)  
                            for j in range(3,6):
                                Vec=np.array((mesh.points[Nodes[0]][0]-mesh.points[Nodes[j]][0],mesh.points[Nodes[0]][1]-mesh.points[Nodes[j]][1],mesh.points[Nodes[0]][2]-mesh.points[Nodes[j]][2]))
                                Dist=abs(np.dot(NormB,Vec)/np.dot(NormB,NormT))
                                Thick[EName][ElemRef[jj][ENum][j-3]]=Dist
    return Mesh(points, [CellBlock('triangle',np.array(cells))], point_data=Thick)
#===================================================================
#
#         Create coating
#
# Variables:
# ThickNames - Names of NodeValue sets with thickness
# NSet       - Set of nodes for coating
# Inside     - If False the function creates coating
#===================================================================
def CreateLayers(mesh, ThickNames, NSet, Inside=False):
    N=len(mesh.point_sets[NSet])
    LayerNum=len(ThickNames)
    thickness=np.zeros((LayerNum,N))
    vertices = np.zeros((LayerNum+1,N,3))
    #-------------
    List=np.full(mesh.points.shape[0],-1,dtype=np.int32)
    for i in range(N):
        Node=mesh.point_sets[NSet][i]
        List[Node]=i
        vertices[0][i][0]=mesh.points[Node][0]
        vertices[0][i][1]=mesh.points[Node][1]
        vertices[0][i][2]=mesh.points[Node][2]
        for j in range(LayerNum):thickness[j][i]=mesh.point_data[ThickNames[j]][Node]
#-----Reading faces
    Faces_dyn=[]
    for block in mesh.cells:
        for Elem in block.data:
            Node0=Elem[0]
            Node1=Elem[1]
            Node2=Elem[2]
            Node3=Elem[3]
            if List[Node0]>-1:
                if List[Node3]>-1:
                    if List[Node1]>-1:
                        Faces_dyn.append((List[Node0],List[Node1],List[Node3]))
                    if List[Node2]>-1:
                        Faces_dyn.append((List[Node0],List[Node3],List[Node2]))
                if List[Node1]>-1 and List[Node2]>-1:
                    Faces_dyn.append((List[Node0],List[Node2],List[Node1]))
            if List[Node1]>-1 and List[Node2]>-1 and List[Node3]>-1:
                Faces_dyn.append((List[Node1],List[Node2],List[Node3]))
    faces=np.array(Faces_dyn)
    FacesNum=len(faces)
    del Faces_dyn
#-----------------------------------
    norm = np.zeros((N,3))
    tris = vertices[0][faces]
    n = np.cross( tris[::,1 ] - tris[::,0]  , tris[::,2 ] - tris[::,0] )
    normalize_v3(n)
    norm[ faces[:,0] ] += n 
    norm[ faces[:,1] ] += n 
    norm[ faces[:,2] ] += n 
    normalize_v3(norm)
#-----------------------------------
    if Inside:
        for j in range(LayerNum):
            vertices[0,:,0]-=norm[:,0] * thickness[j]
            vertices[0,:,1]-=norm[:,1] * thickness[j]
            vertices[0,:,2]-=norm[:,2] * thickness[j]
        for i in range(N):
            Node=mesh.point_sets[NSet][i]
            mesh.points[Node][0]=vertices[0][i][0]
            mesh.points[Node][1]=vertices[0][i][1]
            mesh.points[Node][2]=vertices[0][i][2]
    for j in range(LayerNum):
        vertices[j+1,:,0]=norm[:,0] * thickness[j] + vertices[j,:,0]
        vertices[j+1,:,1]=norm[:,1] * thickness[j] + vertices[j,:,1]
        vertices[j+1,:,2]=norm[:,2] * thickness[j] + vertices[j,:,2]
    #------------output-----------------------
    NodesNum=np.zeros((LayerNum+1,FacesNum,3),dtype=np.int64)
    for i in range(FacesNum):
        NodesNum[0][i]=[mesh.point_sets[NSet][faces[i,0]],mesh.point_sets[NSet][faces[i,1]],mesh.point_sets[NSet][faces[i,2]]]
    points=list(mesh.points)
    NdNum=len(points)
    ElmNum=0
    NewNodeNums={}
    Elems=mesh.cells.copy()
    BlockNum=len(Elems)
    cell_sets=mesh.cell_sets.copy()
    for cellset in cell_sets:
        cell_sets[cell_set].append([])
    for j in range(LayerNum):
        for i in range(N):    
            if thickness[j,i]>0:
                points.append(np.array((vertices[j+1][i][0],vertices[j+1][i][1],vertices[j+1][i][2])))
                NewNodeNums[j*N+i]=NdNum
                NdNum+=1
        cell_sets[ThickNames[j]]=[]
        for i in range(BlockNum+1):
            cell_sets[ThickNames[j]].append([])
        NodesNum[j+1]=self.MaxNodeNum+1+j*N+faces
        for i in range(FacesNum):
            cell_sets[ThickNames[j]][BlockNum].append(ElmNum)
            ElmNum+=1
            self.Elems[self.MaxElemNum+1+FacesNum*j+i]=list()
            for Node in NodesNum[j][i]: self.Elems[self.MaxElemNum+1+FacesNum*j+i]+=(Node,)
            for Node in NodesNum[j+1][i]: self.Elems[self.MaxElemNum+1+FacesNum*j+i]+=(Node,)
    self.MaxNodeNum+=N*LayerNum
    self.MaxElemNum+=FacesNum*LayerNum
    return Mesh(np.array(points), Elems, cell_sets=cell_sets)
#===================================================================
#
#         Create a file with equations for Abaqus for nodes on two symmetrical cuts around X-axis
#
# Variables:
# FileName - File Name for output
# NSet1 - Name of a set of nodes on the first cut
# NSet2 - Name of a set of nodes on the second cut
# tolerance - Distance error over R and over axis (pair of nodes is chosen inside the tolerance)
# method - 'Tolerance' is the fastest method, but points can be omitted;
#          'Nearest' - equations for all points in NSet1
#===================================================================
def SymmetryEquations(mesh,FileName,NSet1,NSet2,method='Tolerance',tolerance=(0.00001,0.00001)):
    if not NSet1 in mesh.point_sets:
        print(NSet1+' hasnt been found')
        return
    else: print(NSet1+' contains '+str(len(mesh.point_sets[NSet1])))
    if not NSet2 in mesh.point_sets:
        print(NSet2+' hasnt been found')
        return        
    else: print(NSet2+' contains '+str(len(mesh.point_sets[NSet2])))
    NodLabels={}
    if 'Node_Num' in mesh.point_data:
        for Node in mesh.point_sets[NSet1]: NodeLabels[Node]=mesh.point_data['Node_Num'][Node]
        for Node in mesh.point_sets[NSet2]: NodeLabels[Node]=mesh.point_data['Node_Num'][Node]
    else:
        for Node in mesh.point_sets[NSet1]: NodeLabels[Node]=Node+1
        for Node in mesh.point_sets[NSet2]: NodeLabels[Node]=Node+1
    Sym1=[]
    Sym2=list(mesh.point_sets[NSet2])
    #------------------------------------------
    #--------------NSET------------------------
    f=open(FileName,'w')
    f.write('*Nset, nset=SYMNODES\n')
    Count=0        
    for i in range(len(mesh.point_sets[NSet1])):
        f.write(str(NodeLabels[mesh.point_sets[NSet1][i]]))
        if Count<15:
            f.write(',')
            Count+=1
        else:
            f.write('\n')
            Count=0
    Num=len(mesh.point_sets[NSet2])
    for i in range(len(mesh.point_sets[NSet2])):
        f.write(str(NodeLabels[mesh.point_sets[NSet2][i]]))
        if Count<15 and i<Num-1:
            f.write(',')
            Count+=1
        else:
            f.write('\n')
            Count=0
    #------------CSYS-------------------------
    f.write('*TRANSFORM, nset=SYMNODES, TYPE=C\n')
    f.write('0,0,0,1,0,0\n')
    #------------EQUATIONS--------------------
    f.write('*equation\n')
    #-------------METHOD: Tolerance--------------------
    if method=='Tolerance':
        for Node1 in mesh.point_sets[NSet1]:
            R1=(mesh.points[Node1][1]**2+mesh.points[Node1][2]**2)**0.5
            Flag=True
            i=0
            while Flag:
                Node2=Sym2[i]
                R2=(mesh.points[Node2][1]**2+mesh.points[Node2][2]**2)**0.5
                if abs(R1-R2)<tolerance[0] and abs(mesh.points[Node1][0]-mesh.points[Node2][0])<tolerance[1]:
                    f.write('2\n')
                    f.write(str(NodeLabels[Node1])+',1,-1,'+str(NodeLabels[Node2])+',1,1\n')
                    f.write('2\n')
                    f.write(str(NodeLabels[Node1])+',2,-1,'+str(NodeLabels[Node2])+',2,1\n')
                    f.write('2\n')
                    f.write(str(NodeLabels[Node1])+',3,-1,'+str(NodeLabels[Node2])+',3,1\n')
                    Sym2.remove(Node2)
                    Flag=False
                i+=1
                if i==len(Sym2):                    
                    if Flag:
                        Sym1.append(Node1)
                        Flag=False                            
    #-------------METHOD: Nearest--------------------
    if method=='Nearest':
        DistMin=0
        DistMax=0
        for Node1 in mesh.point_sets[NSet1]:
            R1=(mesh.points[Node1][1]**2+mesh.points[Node1][2]**2)**0.5
            for i in range(len(Sym2)):
                Node2=Sym2[i]
                R2=(mesh.points[Node2][1]**2+mesh.points[Node2][2]**2)**0.5
                Dist=((R1-R2)**2+(mesh.points[Node1][0]-mesh.points[Node2][0])**2)**0.5
                if i==0 or DistMin>Dist:
                    DistMin=Dist
                    Node2min=Node2
            f.write('2\n')
            f.write(str(NodeLabels[Node1])+',1,-1,'+str(NodeLabels[Node2min])+',1,1\n')
            f.write('2\n')
            f.write(str(NodeLabels[Node1])+',2,-1,'+str(NodeLabels[Node2min])+',2,1\n')
            f.write('2\n')
            f.write(str(NodeLabels[Node1])+',3,-1,'+str(NodeLabels[Node2min])+',3,1\n')
            Sym2.remove(Node2min)
            if DistMax<DistMin:DistMax=DistMin
        print('Maximum distance: '+str(DistMax))
    f.close()
    #-------------Statistics--------------------
    print(str(len(Sym1))+' nodes where a pair has not been found for NSet1:')
    print(str(len(Sym2))+' nodes where a pair has not been found for NSet2:')
    print('Use '+FileName+'_err for error details')        
    f=open(FileName+'_err','w')
    Num=len(Sym1)
    if Num>0:            
        Count=0
        f.write('*Nset, nset=SYM1_err\n')
        for i in range(len(Sym1)):
            f.write(str(NodeLabels[Sym1[i]]))
            if Count<15 and i<Num-1:
                f.write(',')
                Count+=1
            else:
                f.write('\n')
                Count=0            
    Num=len(Sym2)
    if Num>0:
        Count=0
        f.write('*Nset, nset=SYM2_err\n')
        for i in range(len(Sym2)):
           f.write(str(NodeLabels[Sym2[i]]))
           if Count<15 and i<Num-1:
                f.write(',')
                Count+=1
           else:
                f.write('\n')
                Count=0
     f.close()
#===================================================================
#
#         Volume mapping by means of a vtu-file
#
# Variables:
# FileName    - Name of a vtu-file with field (vtkXMLUnstructuredGridReader with vtkFloatArrays)
# NodeSet     - Name of a set of nodes for mapping 
# Tlrnc       - Tolerance. If relative distance is more than Tlrnc, the minimum distance method is used
#===================================================================
def mapping(mesh_targ,FileName,NodeSet,Tlrnc=0.005):
    Reader=vtk.vtkXMLUnstructuredGridReader()
    Reader.SetFileName(FileName)
    Reader.Update()
    vtkData=Reader.GetOutput()
    FieldNum=vtkData.GetPointData().GetNumberOfArrays()
    for i in range(FieldNum):
        mesh.points_data[vtkData.GetPointData().GetArrayName(i)]=np.zeros(mesh.points.shape[0])
    Cell_Num=vtkData.GetNumberOfCells()
    Mtrxs=np.zeros((Cell_Num,3,3))
    M=np.zeros((3,3))
    V=np.zeros((FieldNum,4))
    Value=np.zeros(FieldNum)
    DX=0
    DY=0
    DZ=0
    Xmin=vtkData.GetCell(0).GetPoints().GetPoint(0)[0]
    Xmax=vtkData.GetCell(0).GetPoints().GetPoint(0)[0]
    Ymin=vtkData.GetCell(0).GetPoints().GetPoint(0)[1]
    Ymax=vtkData.GetCell(0).GetPoints().GetPoint(0)[1]
    Zmin=vtkData.GetCell(0).GetPoints().GetPoint(0)[2]
    Zmax=vtkData.GetCell(0).GetPoints().GetPoint(0)[2]
    for i in range(Cell_Num):
        Points=vtkData.GetCell(i).GetPoints()
        XElmin=Points.GetPoint(0)[0]
        XElmax=Points.GetPoint(0)[0]
        YElmin=Points.GetPoint(0)[1]
        YElmax=Points.GetPoint(0)[1]
        ZElmin=Points.GetPoint(0)[2]
        ZElmax=Points.GetPoint(0)[2]    
        for k in range(3):
            for j in range(3):
                M[j][k]=Points.GetPoint(1+k)[j]-Points.GetPoint(0)[j]
            if XElmin>Points.GetPoint(1+k)[0]:XElmin=Points.GetPoint(1+k)[0]
            if XElmax<Points.GetPoint(1+k)[0]:XElmax=Points.GetPoint(1+k)[0]
            if YElmin>Points.GetPoint(1+k)[1]:YElmin=Points.GetPoint(1+k)[1]
            if YElmax<Points.GetPoint(1+k)[1]:YElmax=Points.GetPoint(1+k)[1]
            if ZElmin>Points.GetPoint(1+k)[2]:ZElmin=Points.GetPoint(1+k)[2]
            if ZElmax<Points.GetPoint(1+k)[2]:ZElmax=Points.GetPoint(1+k)[2]
        Mtrxs[i]=np.linalg.inv(M)
        DX+=XElmax-XElmin
        DY+=YElmax-YElmin
        DZ+=ZElmax-ZElmin
        if Xmin>XElmin:Xmin=XElmin
        if Xmax<XElmax:Xmax=XElmax
        if Ymin>YElmin:Ymin=YElmin
        if Ymax<YElmax:Ymax=YElmax
        if Zmin>ZElmin:Zmin=ZElmin
        if Zmax<ZElmax:Zmax= ZElmax
    for Node in self.NSets[NodeSet]:
        if Xmin>mesh.points[Node][0]:Xmin=mesh.points[Node][0]
        if Xmax<mesh.points[Node][0]:Xmax=mesh.points[Node][0]
        if Ymin>mesh.points[Node][1]:Ymin=mesh.points[Node][1]
        if Ymax<mesh.points[Node][1]:Ymax=mesh.points[Node][1]
        if Zmin>mesh.points[Node][2]:Zmin=mesh.points[Node][2]
        if Zmax<mesh.points[Node][2]:Zmax=mesh.points[Node][2]
    DX/=Cell_Num
    DY/=Cell_Num
    DZ/=Cell_Num
    Xmax+=DX*Tlrnc
    Ymax+=DY*Tlrnc
    Zmax+=DZ*Tlrnc
    Nx=int((Xmax-Xmin)/DX)
    Ny=int((Ymax-Ymin)/DY)
    Nz=int((Zmax-Zmin)/DZ)
    DX=(Xmax-Xmin)/Nx
    DY=(Ymax-Ymin)/Ny
    DZ=(Zmax-Zmin)/Nz
    CellDistr=[]
    for i in range(Nx):
        CellDistr.append([])
        for j in range(Ny):
            CellDistr[i].append([])
            for k in range(Nz):
                CellDistr[i][j].append([])
    for Cell_i in range(Cell_Num):
        Points=vtkData.GetCell(Cell_i).GetPoints()
        XElmin=Points.GetPoint(0)[0]
        XElmax=Points.GetPoint(0)[0]
        YElmin=Points.GetPoint(0)[1]
        YElmax=Points.GetPoint(0)[1]
        ZElmin=Points.GetPoint(0)[2]
        ZElmax=Points.GetPoint(0)[2]
        for k in range(3):
            if XElmin>Points.GetPoint(1+k)[0]:XElmin=Points.GetPoint(1+k)[0]
            if XElmax<Points.GetPoint(1+k)[0]:XElmax=Points.GetPoint(1+k)[0]
            if YElmin>Points.GetPoint(1+k)[1]:YElmin=Points.GetPoint(1+k)[1]
            if YElmax<Points.GetPoint(1+k)[1]:YElmax=Points.GetPoint(1+k)[1]
            if ZElmin>Points.GetPoint(1+k)[2]:ZElmin=Points.GetPoint(1+k)[2]
            if ZElmax<Points.GetPoint(1+k)[2]:ZElmax=Points.GetPoint(1+k)[2]        
        for i in range(int((XElmin-Xmin)/DX),int((XElmax-Xmin)/DX)+1):
            for j in range(int((YElmin-Ymin)/DY),int((YElmax-Ymin)/DY)+1):
                for k in range(int((ZElmin-Zmin)/DZ),int((ZElmax-Zmin)/DZ)+1):
                    CellDistr[i][j][k].append(Cell_i)
    #======================================
    MinDistNodes=[]
    NodeWOEl={}
    for Node in mesh.point_sets[NodeSet]:
        GlPoint=mesh.points[Node]
        ip=int((GlPoint[0]-Xmin)/DX)
        jp=int((GlPoint[1]-Ymin)/DY)
        kp=int((GlPoint[2]-Zmin)/DZ)
        if len(CellDistr[ip][jp][kp])==0:
            if not ip in NodeWOEl: NodeWOEl[ip]={}
            if not jp in NodeWOEl[ip]: NodeWOEl[ip][jp]={}
            if not kp in NodeWOEl[ip][jp]: NodeWOEl[ip][jp][kp]=[]
            NodeWOEl[ip][jp][kp].append(Node)
        MinDist=0
        MinCell=0
        MinNode=0
        Flag=True
        for i in CellDistr[ip][jp][kp]:
            Points=vtkData.GetCell(i).GetPoints()
            for j in range(4):
                PntCoord=Points.GetPoint(j)
                Dist=((GlPoint[0]-PntCoord[0])**2+(GlPoint[1]-PntCoord[1])**2+(GlPoint[2]-PntCoord[2])**2)**0.5            
                if (i==CellDistr[ip][jp][kp][0] and j==0) or Dist<MinDist:
                    MinDist=Dist
                    MinCell=i
                    MinNode=j
            LcPoint=np.dot(Mtrxs[i],GlPoint-np.array(Points.GetPoint(0)))
            if LcPoint[0]>=-Tlrnc and LcPoint[1]>=-Tlrnc and LcPoint[2]>=-Tlrnc and (LcPoint[0]+LcPoint[1]+LcPoint[2])<=1+Tlrnc:
                for j in range(4):
                    CellNode=vtkData.GetCell(i).GetPointIds().GetId(j)
                    for FN in range(FieldNum):
                        V[FN][j]=vtkData.GetPointData().GetArray(FN).GetValue(CellNode)
                for FN in range(FieldNum):
                    Value[FN]=V[FN][0]+(V[FN][1]-V[FN][0])*LcPoint[0]+(V[FN][2]-V[FN][0])*LcPoint[1]+(V[FN][3]-V[FN][0])*LcPoint[2]            
                Flag=False
                break
        if Flag:
            CellNode=vtkData.GetCell(MinCell).GetPointIds().GetId(MinNode)
            for FN in range(FieldNum):
                Value[FN]=vtkData.GetPointData().GetArray(FN).GetValue(CellNode)
            MinDistNodes.append(Node)
        for FN in range(FieldNum):
            mesh.point_data[vtkData.GetPointData().GetArrayName(FN)][Node]=Value[FN]
#-----Nodes in the cells of the grid without field elements-----
    for ip in NodeWOEl:
        for jp in NodeWOEL[ip]:
            for kp in NodeWOEL[ip][jp]:
                Box=[]
                shift=1
                icrit=max(ip,jp,kp,Nx-ip-1,Ny-jp-1,Nz-kp-1)
                while len(Box)==0 and shift<=icrit:
                    ip0=max(0,ip-shift+1)
                    ip1=min(Nx,ip+shift-1)
                    jp0=max(0,jp-shift)
                    jp1=min(Ny,jp+shift)
                    kp0=max(0,kp-shift)
                    kp1=min(Nz,kp+shift)
                    if ip-shift>=0:
                        for j in range(jp0,jp1):
                            for k in range(kp0,kp1):
                                Box+=CellDistr[ip-shift][j][k]
                    if ip+shift<Nx:
                        for j in range(jp0,jp1):
                            for k in range(kp0,kp1):
                               Box+=CellDistr[ip+shift][j][k]
                    if jp-shift>=0:
                        for i in range(ip0,ip1):
                            for k in range(kp0,kp1):
                               Box+=CellDistr[i][jp-shift][k]
                    if jp+shift<Ny:
                        for i in range(ip0,ip1):
                            for k in range(kp0,kp1):
                               Box+=CellDistr[i][jp+shift][k]
                    jp0=max(0,jp-shift+1)
                    jp1=min(Ny,jp+shift-1)
                    if kp-shift>=0:
                        for i in range(ip0,ip1):
                            for j in range(jp0,jp1):
                                Box+=CellDistr[i][j][kp-shift]
                    if kp+shift<Nz:
                        for i in range(ip0,ip1):
                           for j in range(jp0,jp1):
                               Box+=CellDistr[i][j][kp+shift]
                    shift+=1
                for Node in NodeWOEl[ip][jp][kp]:
                    GlPoint=mesh.points[Node]
                    MinDist=0
                    MinCell=0
                    MinNode=0
                    for i in Box:
                        Points=vtkData.GetCell(i).GetPoints()
                        for j in range(4):
                            PntCoord=Points.GetPoint(j)
                            Dist=((GlPoint[0]-PntCoord[0])**2+(GlPoint[1]-PntCoord[1])**2+(GlPoint[2]-PntCoord[2])**2)**0.5
                            if (i==Box[0] and j==0) or (Dist<MinDist):
                                MinDist=Dist
                                MinCell=i
                                MinNode=j
                    CellNode=vtkData.GetCell(MinCell).GetPointIds().GetId(MinNode)
                    for FN in range(FieldNum):
                        Value[FN]=vtkData.GetPointData().GetArray(FN).GetValue(CellNode)
                    for FN in range(FieldNum):
                        mesh.point_data[vtkData.GetPointData().GetArrayName(FN)][Node]=Value[FN]
#----------------------------------------------
    if len(MinDistNodes)>0:
        print('The nearest method was applied to '+str(len(MinDistNodes))+' nodes')
        print('See MinDistNodes list')
        mesh.data_sets['MinDistNodes']=MinDistNodes
#===================================================================
#
#         Mapping from surface data on nodes
#
# Variables:
# FileName - Name of a vtu-file (vtkXMLUnstructuredGridReader with vtkFloatArrays)
# NodeSet - Name of a set of nodes for mapping
# DistError - Distance error (just for messaging)
#===================================================================
def map_surf(mesh,FileName,NodeSet,DistError=0.0001):
    Reader=vtk.vtkXMLUnstructuredGridReader()
    Reader.SetFileName(FileName)
    Reader.Update()
    vtkSurfdData=Reader.GetOutput()
    Cell_Num=vtkSurfdData.GetNumberOfCells()
    Mtrxs=np.zeros((Cell_Num,3,3))
    M=np.zeros((3,3))
    V1=np.zeros(3)
    V2=np.zeros(3)
    V3=np.zeros(3)
    for i in range(Cell_Num):
        Points=vtkSurfdData.GetCell(i).GetPoints()
        for j in range(3):
            V1[j]=Points.GetPoint(1)[j]-Points.GetPoint(0)[j]
            V2[j]=Points.GetPoint(2)[j]-Points.GetPoint(0)[j]
            M[j][0]=V1[j]
            M[j][1]=V2[j]
        Norm=np.cross(V1,V2)
        Norm=Norm/np.linalg.norm(Norm)
        for j in range(3): M[j][2]=Norm[j]
        Mtrxs[i]=np.linalg.inv(M)
    #======================================
    for j in range(vtkSurfdData.GetPointData().GetNumberOfArrays()):
        mesh.point_data[vtkSurfdData.GetPointData().GetArray(j).GetName()]={}
    mesh.point_sets['NodesOutOfTolerance']=[]
    for Node in self.NSets[NodeSet]:
        GlPoint=mesh.points[Node]
        Flag=True
        i=0
        MinDist=0
        TolFlag=True
        while Flag:
            Points=vtkSurfdData.GetCell(i).GetPoints()
            for j in range(3):
                V1[j]=Points.GetPoint(0)[j]-GlPoint[j]
                V2[j]=Points.GetPoint(1)[j]-GlPoint[j]
                V3[j]=Points.GetPoint(2)[j]-GlPoint[j]
            Dist=min((np.linalg.norm(V1),np.linalg.norm(V2),np.linalg.norm(V3)))
            LcPoint=np.dot(Mtrxs[i],GlPoint-np.array(Points.GetPoint(0)))
            if LcPoint[0]>=0 and LcPoint[1]>=0 and (LcPoint[0]+LcPoint[1])<=1 and abs(LcPoint[2])<=DistError:
                Flag=False
                i_Cell=i
                Ksi=LcPoint[0]
                Nu=LcPoint[1]
                TolFlag=False
            elif LcPoint[0]>=0 and LcPoint[1]>=0 and (LcPoint[0]+LcPoint[1])<=1 and abs(LcPoint[2])>DistError:
                if i==0 or MinDist>LcPoint[2]:
                    i_Cell=i
                    Ksi=LcPoint[0]
                    Nu=LcPoint[1]
                    MinDist=LcPoint[2]
            else:
                if i==0 or MinDist>Dist:
                    i_Cell=i
                    Ksi=LcPoint[0]
                    Nu=LcPoint[1]
                    MinDist=Dist
            i+=1
            if i==Cell_Num: Flag=False
        if TolFlag: mesh.point_sets['NodesOutOfTolerance'].append(Node)
        for j in range(vtkSurfdData.GetPointData().GetNumberOfArrays()):
            for k in range(3):
                CellNode=vtkSurfdData.GetCell(i_Cell).GetPointIds().GetId(k)
                V1[k]=vtkSurfdData.GetPointData().GetArray(j).GetValue(CellNode)
            Value=V1[0]+(V1[1]-V1[0])*Ksi+(V1[2]-V1[0])*Nu
            mesh.point_data[vtkSurfdData.GetPointData().GetArray(j).GetName()][Node]=Value
    if len(mesh.point_data['NodesOutOfTolerance'])>0:
        print(str(len(mesh.point_sets['NodesOutOfTolerance']))+' nodes are out of tolerance')
        print('Look at "NodesOutOfTolerance" node set')
#===================================================================
#
#         Mapping from surface data on faces
#
# Variables:
# FileName - Name of a vtu-file (vtkXMLUnstructuredGridReader with vtkFloatArrays)
# SetName - Name of a set of nodes/surfaces for mapping
# DistError - Distance error (just for messaging)
# method - 'NODE', 'FACE'
#===================================================================
def map_surf(mesh,FileName,SetName,DistError=0.0001,method='FACE'):
    Reader=vtk.vtkXMLUnstructuredGridReader()
    Reader.SetFileName(FileName)
    Reader.Update()
    vtkSurfdData=Reader.GetOutput()
    Cell_Num=vtkSurfdData.GetNumberOfCells()
    Mtrxs=np.zeros((Cell_Num,3,3))
    M=np.zeros((3,3))
    V1=np.zeros(3)
    V2=np.zeros(3)
    V3=np.zeros(3)
    GlPoint=np.zeros(3)
	DX=0
	DY=0
	DZ=0
    Xmin=vtkData.GetCell(0).GetPoints().GetPoint(0)[0]
    Xmax=vtkData.GetCell(0).GetPoints().GetPoint(0)[0]
    Ymin=vtkData.GetCell(0).GetPoints().GetPoint(0)[1]
    Ymax=vtkData.GetCell(0).GetPoints().GetPoint(0)[1]
    Zmin=vtkData.GetCell(0).GetPoints().GetPoint(0)[2]
    Zmax=vtkData.GetCell(0).GetPoints().GetPoint(0)[2]
    for i in range(Cell_Num):
        Points=vtkData.GetCell(i).GetPoints()
        XElmin=Points.GetPoint(0)[0]
        XElmax=Points.GetPoint(0)[0]
        YElmin=Points.GetPoint(0)[1]
        YElmax=Points.GetPoint(0)[1]
        ZElmin=Points.GetPoint(0)[2]
        ZElmax=Points.GetPoint(0)[2]
        for k in range(2):
            if XElmin>Points.GetPoint(1+k)[0]:XElmin=Points.GetPoint(1+k)[0]
            if XElmax<Points.GetPoint(1+k)[0]:XElmax=Points.GetPoint(1+k)[0]
            if YElmin>Points.GetPoint(1+k)[1]:YElmin=Points.GetPoint(1+k)[1]
            if YElmax<Points.GetPoint(1+k)[1]:YElmax=Points.GetPoint(1+k)[1]
            if ZElmin>Points.GetPoint(1+k)[2]:ZElmin=Points.GetPoint(1+k)[2]
            if ZElmax<Points.GetPoint(1+k)[2]:ZElmax=Points.GetPoint(1+k)[2]
        for j in range(3):
            V1[j]=Points.GetPoint(1)[j]-Points.GetPoint(0)[j]
            V2[j]=Points.GetPoint(2)[j]-Points.GetPoint(0)[j]
            M[j][0]=V1[j]
            M[j][1]=V2[j]
        Norm=np.cross(V1,V2)
        Norm=Norm/np.linalg.norm(Norm)
        for j in range(3): M[j][2]=Norm[j]
        Mtrxs[i]=np.linalg.inv(M)
        DX+=XElmax-XElmin
        DY+=YElmax-YElmin
        DZ+=ZElmax-ZElmin
        if Xmin>XElmin:Xmin=XElmin
        if Xmax<XElmax:Xmax=XElmax
        if Ymin>YElmin:Ymin=YElmin
        if Ymax<YElmax:Ymax=YElmax
        if Zmin>ZElmin:Zmin=ZElmin
        if Zmax<ZElmax:Zmax= ZElmax
    DX/=Cell_Num
    DY/=Cell_Num
    DZ/=Cell_Num
	if method=='NODE':
	    FieldNum=vtkSurfdData.GetPointData().GetNumberOfArrays()
        for j in range(FieldNum):
	        if not vtkSurfdData.GetCellData().GetArray(j).GetName() in mesh.point_data:
                mesh.point_data[vtkSurfdData.GetPointData().GetArray(j).GetName()]={}
        for Node in mesh.point_sets[NodeSet]:
            if Xmin>mesh.points[Node][0]:Xmin=mesh.points[Node][0]
            if Xmax<mesh.points[Node][0]:Xmax=mesh.points[Node][0]
            if Ymin>mesh.points[Node][1]:Ymin=mesh.points[Node][1]
            if Ymax<mesh.points[Node][1]:Ymax=mesh.points[Node][1]
            if Zmin>mesh.points[Node][2]:Zmin=mesh.points[Node][2]
            if Zmax<mesh.points[Node][2]:Zmax=mesh.points[Node][2]
    if method=='FACE':
        FieldNum=vtkSurfdData.GetCellData().GetNumberOfArrays()
        for j in range(FieldNum):
	        if not vtkSurfdData.GetCellData().GetArray(j).GetName() in mesh.FaceLoad:
                mesh.FaceLoad[vtkSurfdData.GetCellData().GetArray(j).GetName()]={}
        for Face in mesh.Surfs[SetName]:
            for ElemNum in mesh.cell_sets[Face[0]]:
                for i in FacesNodes[mesh.Eltype[ElemNum]][Face[1]]:
                    Node=mesh.cells[ElemNum][i]
                    if Xmin>mesh.points[Node][0]:Xmin=mesh.points[Node][0]
                    if Xmax<mesh.points[Node][0]:Xmax=mesh.points[Node][0]
                    if Ymin>mesh.points[Node][1]:Ymin=mesh.points[Node][1]
                    if Ymax<mesh.points[Node][1]:Ymax=mesh.points[Node][1]
                    if Zmin>mesh.points[Node][2]:Zmin=mesh.points[Node][2]
                    if Zmax<mesh.points[Node][2]:Zmax=mesh.points[Node][2]
    Xmax+=DX*DistError
    Ymax+=DY*DistError
    Zmax+=DZ*DistError
    Nx=int((Xmax-Xmin)/DX)
    Ny=int((Ymax-Ymin)/DY)
    Nz=int((Zmax-Zmin)/DZ)
    DX=(Xmax-Xmin)/Nx
    DY=(Ymax-Ymin)/Ny
    DZ=(Zmax-Zmin)/Nz
    CellDistr=[]
    for i in range(Nx):
        CellDistr.append([])
        for j in range(Ny):
            CellDistr[i].append([])
            for k in range(Nz):
                CellDistr[i][j].append([])
    for Cell_i in range(Cell_Num):
        Points=vtkData.GetCell(Cell_i).GetPoints()
        XElmin=Points.GetPoint(0)[0]
        XElmax=Points.GetPoint(0)[0]
        YElmin=Points.GetPoint(0)[1]
        YElmax=Points.GetPoint(0)[1]
        ZElmin=Points.GetPoint(0)[2]
        ZElmax=Points.GetPoint(0)[2]
        for k in range(3):
            if XElmin>Points.GetPoint(1+k)[0]:XElmin=Points.GetPoint(1+k)[0]
            if XElmax<Points.GetPoint(1+k)[0]:XElmax=Points.GetPoint(1+k)[0]
            if YElmin>Points.GetPoint(1+k)[1]:YElmin=Points.GetPoint(1+k)[1]
            if YElmax<Points.GetPoint(1+k)[1]:YElmax=Points.GetPoint(1+k)[1]
            if ZElmin>Points.GetPoint(1+k)[2]:ZElmin=Points.GetPoint(1+k)[2]
            if ZElmax<Points.GetPoint(1+k)[2]:ZElmax=Points.GetPoint(1+k)[2]        
        for i in range(int((XElmin-Xmin)/DX),int((XElmax-Xmin)/DX)+1):
            for j in range(int((YElmin-Ymin)/DY),int((YElmax-Ymin)/DY)+1):
                for k in range(int((ZElmin-Zmin)/DZ),int((ZElmax-Zmin)/DZ)+1):
                    CellDistr[i][j][k].append(Cell_i)
    #======================================
    mesh.point_sets['NodesOutOfTolerance']=[]
    mesh.Surfs['FacesOutOfTolerance']=[]
    MinDist=[]
    NodeWOEl={}
    Data=[]
    if method=='NODE':
        Data.append([mesh.point_sets[SetName],0])
    if method=='FACE':
        for Face in mesh.Surfs[SetName]:
            Data.append([mesh.cell_sets[Face[0]],Face[1]])
    for Face in Data:
        for Indx in Face[0]:
            GlPoint.fill(0)
            if method=='NODE':
                GlPoint+=mesh.points[Indx]
            if method=='FACE':
                for i in FacesNodes[mesh.Eltype[Indx]][Face[1]]:
                    GlPoint+=mesh.points[mesh.cells[Indx][i]]
                GlPoint/=len(FacesNodes[mesh.Eltype[Indx]][Face[1]])
            ip=int((GlPoint[0]-Xmin)/DX)
            jp=int((GlPoint[1]-Ymin)/DY)
            kp=int((GlPoint[2]-Zmin)/DZ)
            if len(CellDistr[ip][jp][kp])==0:
                if not ip in NodeWOEl: NodeWOEl[ip]={}
                if not jp in NodeWOEl[ip]: NodeWOEl[ip][jp]={}
                if not kp in NodeWOEl[ip][jp]: NodeWOEl[ip][jp][kp]=[]
                NodeWOEl[ip][jp][kp].append([Indx,Face[1]])
                if method=='NODE': mesh.point_sets['NodesOutOfTolerance'].append(Indx)
                if method=='FACE':
                    if not 'FacesOutOfTolerance_S'+str(Face[1]) in mesh.cell_sets:
                        mesh.cell_sets['FacesOutOfTolerance_S'+str(Face[1])]=[]
                        mesh.Surfs['FacesOutOfTolerance'].append(['FacesOutOfTolerance_S'+str(Face[1]),Face[1]])
                    mesh.cell_sets['FacesOutOfTolerance_S'+str(Face[1])].append(Indx)                  
            MinDist=0
            Flag=True
            for i in range(len(CellDistr[ip][jp][kp])):
                Points=vtkSurfdData.GetCell(CellDistr[ip][jp][kp][i]).GetPoints()
                for j in range(3):
                    V1[j]=Points.GetPoint(0)[j]-GlPoint[j]
                    V2[j]=Points.GetPoint(1)[j]-GlPoint[j]
                    V3[j]=Points.GetPoint(2)[j]-GlPoint[j]
                Dist=min((np.linalg.norm(V1),np.linalg.norm(V2),np.linalg.norm(V3)))
                LcPoint=np.dot(Mtrxs[CellDistr[ip][jp][kp][i]],GlPoint-np.array(Points.GetPoint(0)))
                if LcPoint[0]>=-DistError and LcPoint[1]>=-DistError and (LcPoint[0]+LcPoint[1])<=1+DistError and abs(LcPoint[2])<=DistError:
                    i_Cell=CellDistr[ip][jp][kp][i]
                    Ksi=LcPoint[0]
                    Nu=LcPoint[1]
                    Flag=False
                    break
                elif LcPoint[0]>=0 and LcPoint[1]>=0 and (LcPoint[0]+LcPoint[1])<=1 and abs(LcPoint[2])>DistError:
                    if i==0 or MinDist>LcPoint[2]:
                        i_Cell=CellDistr[ip][jp][kp][i]
                        Ksi=LcPoint[0]
                        Nu=LcPoint[1]
                        MinDist=LcPoint[2]
                else:
                    if i==0 or MinDist>Dist:
                        i_Cell=CellDistr[ip][jp][kp][i]
                        Ksi=LcPoint[0]
                        Nu=LcPoint[1]
                        MinDist=Dist
            if method=='NODE':
                if Flag: mesh.point_sets['NodesOutOfTolerance'].append(Indx)
                for j in range(FieldNum):
                    for k in range(3):
                        CellNode=vtkSurfdData.GetCell(i_Cell).GetPointIds().GetId(k)
                        V1[k]=vtkSurfdData.GetPointData().GetArray(j).GetValue(CellNode)
                    Value=V1[0]+(V1[1]-V1[0])*Ksi+(V1[2]-V1[0])*Nu
                    mesh.point_data[vtkSurfdData.GetPointData().GetArray(j).GetName()][Indx]=Value
            if method=='FACE':
               if Flag:
                    if not 'FacesOutOfTolerance_S'+str(Face[1]) in mesh.cell_sets:
                        mesh.cell_sets['FacesOutOfTolerance_S'+str(Face[1])]=[]
                        mesh.Surfs['FacesOutOfTolerance'].append(['FacesOutOfTolerance_S'+str(Face[1]),Face[1]])
                    mesh.cell_sets['FacesOutOfTolerance_S'+str(Face[1])].append(Indx)
                for j in range(FieldNum):
                    Value=vtkSurfdData.GetCellData().GetArray(j).GetValue(i_Cell)
                    if not Indx in mesh.FaceLoad[vtkSurfdData.GetCellData().GetArray(j).GetName()]:
                        mesh.FaceLoad[vtkSurfdData.GetCellData().GetArray(j).GetName()][Indx]=[]
                    mesh.FaceLoad[vtkSurfdData.GetCellData().GetArray(j).GetName()][Indx].append([Face[1],Value])
#-----Nodes/Faces in the cells of the grid without field elements-----
    for ip in NodeWOEl:
        for jp in NodeWOEL[ip]:
            for kp in NodeWOEL[ip][jp]:
                Box=[]
                shift=1
                icrit=max(ip,jp,kp,Nx-ip-1,Ny-jp-1,Nz-kp-1)
                while len(Box)==0 and shift<=icrit:
                    ip0=max(0,ip-shift+1)
                    ip1=min(Nx,ip+shift-1)
                    jp0=max(0,jp-shift)
                    jp1=min(Ny,jp+shift)
                    kp0=max(0,kp-shift)
                    kp1=min(Nz,kp+shift)
                    if ip-shift>=0:
                        for j in range(jp0,jp1):
                            for k in range(kp0,kp1):
                                Box+=CellDistr[ip-shift][j][k]
                    if ip+shift<Nx:
                        for j in range(jp0,jp1):
                             for k in range(kp0,kp1):
                                Box+=CellDistr[ip+shift][j][k]
                    if jp-shift>=0:
                        for i in range(ip0,ip1):
                            for k in range(kp0,kp1):
                                Box+=CellDistr[i][jp-shift][k]
                    if jp+shift<Ny:
                        for i in range(ip0,ip1):
                            for k in range(kp0,kp1):
                                Box+=CellDistr[i][jp+shift][k]
                    jp0=max(0,jp-shift+1)
                    jp1=min(Ny,jp+shift-1)
                    if kp-shift>=0:
                        for i in range(ip0,ip1):
                            for j in range(jp0,jp1):
                                Box+=CellDistr[i][j][kp-shift]
                    if kp+shift<Nz:
                        for i in range(ip0,ip1):
                            for j in range(jp0,jp1):
                                Box+=CellDistr[i][j][kp+shift]
                    shift+=1
                for Face in NodeWOEl[ip][jp][kp]:
                    Indx=Face[0]
                    GlPoint.fill(0)
                    if method=='NODE':
                        GlPoint+=mesh.points[Indx]
                    if method=='FACE':
                        for i in FacesNodes[mesh.Eltype[Indx]][Face[1]]:
                            GlPoint+=mesh.points[mesh.cells[Indx][i]]
                        GlPoint/=len(FacesNodes[mesh.Eltype[Indx]][Face[1]])
                    MinDist=0
                    MinCell=0
                    MinNode=0
                    for i in Box:
                        Points=vtkSurfdData.GetCell(i).GetPoints()
                        for j in range(3):
                            PntCoord=Points.GetPoint(j)
                            Dist=((GlPoint[0]-PntCoord[0])**2+(GlPoint[1]-PntCoord[1])**2+(GlPoint[2]-PntCoord[2])**2)**0.5
                            if (i==Box[0] and j==0) or (Dist<MinDist):
                                MinDist=Dist
                                MinCell=i
                                MinNode=j
                    if method=='NODE':
                        CellNode=vtkSurfdData.GetCell(MinCell).GetPointIds().GetId(MinNode)
                        for j in range(FieldNum):
                            Value=vtkSurfdData.GetPointData().GetArray(j).GetValue(CellNode)
                            mesh.point_data[vtkSurfdData.GetPointData().GetArrayName(j)][Indx]=Value
                    if method=='FACE':
                        for j in range(FieldNum):
                            Value=vtkSurfdData.GetCellData().GetArray(j).GetValue(MinCell)
                            if not Indx in mesh.FaceLoad[vtkSurfdData.GetCellData().GetArray(j).GetName()]:
                                mesh.FaceLoad[vtkSurfdData.GetCellData().GetArray(j).GetName()][Indx]=[]
                            mesh.FaceLoad[vtkSurfdData.GetCellData().GetArray(j).GetName()][Indx].append([Face[1],Value])
#--------------------------------------------------
        if len(mesh.point_sets['NodesOutOfTolerance'])>0:
            print(str(len(mesh.point_sets['NodesOutOfTolerance']))+' nodes are out of tolerance')
            print('Look at "NodesOutOfTolerance" node set')    
        if len(mesh.Surfs['FacesOutOfTolerance'])>0:
            Count=0
            for Face in mesh.Surfs['FacesOutOfTolerance']:
                Count+=len(mesh.cell_sets['FacesOutOfTolerance_S'+str(Face[1])])
            print(str(Count)+' faces are out of tolerance')
            print('Look at "FacesOutOfTolerance" set')  
#===================================================================
#
#         Mapping for 2D tasks
# Create FaceLoad from csv-file
#
# Variables:
# FileName  - Name of a csv-file (coordinate;Value 1; Value 2)
# SurfName  - Name of a surface for mapping
# i_coord   - Index of coordinate (0-x; 1-y; 2-z)
# LoadType  - 'P' - Pressure; 'S' - Heat Flux; 'F' - HTC
# separator - Separator in the csv-file
#===================================================================
def map_edge(mesh,FileName,SurfName,i_coord,LoadType,separator=';'):
    if not LoadType in mesh.FaceLoad: mesh.FaceLoad[LoadType]={}
    f=open(FileName,'r')
    Values0=list(map(float,f.readline().split(separator)))
    ValN=len(Values0)
    txt=f.readline()
    while txt:
        Values=list(map(float,txt.split(separator)))
        for SSet in mesh.Surfs[SurfName]:
            for El in mesh.cell_sets[SSet[0]]:
                Coord=0
                for Nd_i in FacesNodes[mesh.Eltype[El]][SSet[1]]:
                    Coord+=mesh.points[mesh.cells[El][Nd_i]][i_coord]
                Coord/=len(FacesNodes[mesh.Eltype[El]][SSet[1]])
                if Coord>=Values0[0] and Coord<Values[0]:
                    if not El in mesh.FaceLoad[LoadType]: mesh.FaceLoad[LoadType][El]=[]
                    Val=[SSet[1]+1,]
                    for i in range(1,ValN):
                        Val.append(Values0[i]+(Values[i]-Values0[i])/(Values[0]-Values0[0])*(Coord-Values0[0]))
                    mesh.FaceLoad[LoadType][El].append(Val)
        txt=f.readline()
        Values0=Values.copy()
    f.close()
#===================================================================
#
#         Divide mesh to simulate crack
#
# Variables:
# NSet - Name of sets of divided nodes
# ESet - Name of elements set on one side of crack
#===================================================================
def crack(mesh,NSet,ESet):
    NodeList={}
    for Node in mesh.point_sets[NSet]:
        NodeList[Node]=0
    for i in range(mesh.cells.shape[0]):
        if mesh.Elems[i]!=1:
            if not i in mesh.cell_sets[ESet]:
                for j in range(len(mesh.cells[i])):
                    if mesh.cells[i][j] in mesh.point_sets[NSet]:
                        if NodeList[mesh.cells[i][j]]==0:
                            mesh.MaxNodeNum+=1
                            NodeList[mesh.cells[i][j]]=mesh.MaxNodeNum                                
                        mesh.cells[i][j]=NodeList[mesh.cells[i][j]]
    mesh.points=np.resize(mesh.points,mesh.MaxNodeNum+1)
    for Node in mesh.point_sets[NSet]:
        if NodeList[Node]!=0:
            for j in range(3): mesh.points[NodeList[Node]]=mesh.points[Node].copy()
            for LoadName in mesh.point_data: mesh.point_data[LoadName][NodeList[Node]]=mesh.point_data[LoadName][Node]
#===================================================================
#
#         Write coordinates of a node set
#
# Variables:
# FileName - Name of file
# NSet     - Name of sets of divided nodes
#===================================================================
def NodesCoord(mesh,FileName,NSet):
    f=open(FileName,'w')
    for Node in mesh.point_sets[NSet]:
        f.write(str(mesh.points[Node][0])+', '+str(mesh.points[Node][1])+', '+str(mesh.points[Node][2])+'\n')
    f.close()
#===================================================================
#
#         Project nodes on the surface defined by elements and nodes
#
# Variables:
# PrjctNodes - Name of set for projected nodes
# SurfElems  - Name of set for elements under the surface
# SurfNodes  - Name of set for nodes on the surface
#===================================================================
def ProjectNodesToSurf(mesh,PrjctNodes,SurfElems,SurfNodes):
        FaceList=[]
        Mtrx=np.zeros((3,3))
        Pnt0=np.zeros(3)
        for i in mesh.cell_sets[SurfElems]:
            for Indx in FacesNodes[mesh.Eltype[i]]:
                Flag=True
                for j in Indx:
                    if not mesh.cells[i][j] in mesh.point_sets[SurfNodes]: Flag=False
                if Flag:
                    Pnt0[:]=mesh.points[mesh.cells[i][Indx[0]]][:]
                    Mtrx[0]=np.cross(mesh.points[mesh.cells[i][Indx[1]]][:]-Pnt0[:],mesh.points[mesh.cells[i][Indx[2]]][:]-Pnt0[:])
                    Mtrx[1][:]=mesh.points[mesh.cells[i][Indx[1]]][:]-Pnt0[:]
                    Mtrx[2][:]=mesh.points[mesh.cells[i][Indx[2]]][:]-Pnt0[:]
                    FaceList.append(((i,Indx),Mtrx,np.dot(Mtrx[0],Pnt0)))
        #-----------------------------
        for Node_pr in mesh.point_sets[PrjctNodes]:
            minDist2=0
            for ElemFace in FaceList:
                minFaceDist2=np.linalg.norm(mesh.points[mesh.cells[ElemFace[0][0]][ElemFace[0][1][0]]][:]-mesh.points[Node_pr][:])
                for k in range(1, len(ElemFace[0][1])):
                    Dist2=np.linalg.norm(mesh.points[mesh.cells[ElemFace[0][0]][ElemFace[0][1][k]]][:]-mesh.points[Node_pr][:])
                    if Dist2<minFaceDist2:minFaceDist2=Dist2
                if (minDist2==0)or minDist2>minFaceDist2:
                    minDist2=minFaceDist2
                    minFace=ElemFace
            PrjctCrd=np.linalg.solve(minFace[1],np.array((minFace[2],np.dot(minFace[1][1],mesh.points[Node_pr]),np.dot(minFace[1][2],mesh.points[Node_pr]))))
            mesh.points[Node_pr]=tuple(PrjctCrd)
#===================================================================
#
#         Create submodel
#
# Variables:
# CentralNodes - List of Nodes nummbers
# Radius       - Radius around nodes to catch elements
#===================================================================
def CreateSubModel(mesh,CentralNodes,Radius):
    if len(mesh.Faces)==0: mesh.EstFaces()
    mesh.point_sets['NAll']=[]
    mesh.cell_sets['EAll']=[]
    mesh.point_sets['SubmodelNodes']=[]
    NodeFlag=np.full(mesh.MaxNodeNum+1,False)
    ElFlag=np.full(mesh.MaxElemNum+1,False)
    for Node in range(1,mesh.MaxNodeNum+1):
        if type(mesh.points[Node])==np.ndarray:
            for CNode in CentralNodes:
                if np.linalg.norm(np.array(mesh.points[Node])[:]-np.array(mesh.points[CNode])[:])<=Radius:
                    NodeFlag[Node]=True                
    for El in range(1,mesh.MaxElemNum+1):
        if mesh.cells[El]!=1:
            for Node in mesh.cells[El]:
                if NodeFlag[Node]:
                    ElFlag[El]=True
                    mesh.cell_sets['EAll'].append(El)
                    break
    for El in mesh.cell_sets['EAll']:
        for Node in mesh.cells[El]:
            if not NodeFlag[Node]:NodeFlag[Node]=True            
    #cleaning
    for Node in range(1,mesh.MaxNodeNum+1):
        if not NodeFlag[Node]: mesh.points[Node]=None
    for El in range(1,mesh.MaxElemNum+1):
        if not ElFlag[El]:
            mesh.cells[El]=1
            mesh.Eltype[El]=0
    for NSet in mesh.point_sets:
        Dict=[]
        for Node in mesh.point_sets[NSet]:
            if type(mesh.points[Node])==np.ndarray: Dict.append(Node)
        mesh.point_sets[NSet]=Dict
    for ESet in mesh.cell_sets.keys():
        Dict=[]
        for El in mesh.cell_sets[ESet]:
            if mesh.cells[El]!=1: Dict.append(El)
        mesh.cell_sets[ESet]=Dict
     #------------------
    Faces0=mesh.Faces.copy()
    mesh.EstFaces()
    for El in mesh.cell_sets['EAll']:
        for Indx in FacesNodes[mesh.Eltype[El]]:
            Nodes=set()
            for i in Indx: Nodes.add(mesh.cells[El][i])
            minNode=min(Nodes)
            maxNode=max(Nodes)
            if minNode in mesh.Faces:
                if maxNode in mesh.Faces[minNode]:
                    for Fc in mesh.Faces[minNode][maxNode]:
                        if Fc[0]==1:
                            for Fc0 in Faces0[minNode][maxNode]:
                                if Fc[1]==Fc0[1] and Fc0[0]==2:
                                    mesh.point_sets['SubmodelNodes']+=list(Fc[1])
    print('Submodel has been prepared')
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
def morph(mesh, NodeSet, func, FreeNodeSet='', Normal=False):
    MinNum=0
    MaxNum=0
    for cellblock in mesh.cell_data['Element_Ids']:
        for ElLabel in cellblock:
            if MinNum==0 or MinNum>ElLabel: MinNum=ElLabel
            if MaxNum<ElLabel: MaxNum=ElLabl
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
                MovedNodes[mesh.point_data['Node_Ids']]=Vect
                OuterNodes[Num]=0
    else:
        for Num in mesh.point_sets[NodeSet]:
            Vect=func(mesh.points[Num])-mesh.points[Num]
            if np.linalg.norm(Vect)>0:
                MovedNodes[mesh.point_data['Node_Ids']]=Vect
                OuterNodes[Num]=0
    if FreeNodeSet:
        for Num in mesh.point_sets[FreeNodeSet]:
            OuterNodes[Num]=0
    f=open('Run_Morphing.inp','w')
    f.write('*INCLUDE, INPUT='+filename+'\n')
    f.write('*MATERIAL, NAME=AUXETIC_MAT\n')
    f.write('*ELASTIC, TYPE=ISOTROPIC\n')
    f.write('1,-0.99\n')
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
    f.write('*STEP, NAME=MORPHING, NLGEOM=YES\n')
    f.write('*STATIC\n')
    f.write('0.1, 1.0, 1e-05, 1\n')
    f.write('*BOUNDARY\n')
    for NobelLabel in MovedNodes:
        for i in range(3):
            f.write(str(NodeLabel)+', '+str(i+1)+', '+str(i+1)+', '+str(MovedNodes[NodeLabel][i])+'\n')
    f.write('*OUTPUT, FIELD, FREQUENCY=999\n')
    f.write('*NODE OUTPUT\n')
    f.write('U\n')
    f.write('*END STEP\n')
    f.close()
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
#===================================================================
#
#         Merge nodes
#
# MasterNds - Nodes with constant coordinates
# SlaveNds  - Movable nodes
#===================================================================
def MergeNodes(mesh, MasterNds, SlaveNds):
    Cnct=np.zeros(mesh.MaxNodeNum+1,dtype=np.int32)
    for SlaveNode in mesh.point_sets[SlaveNds]:
        Cnct[SlaveNode]=mesh.point_sets[MasterNds][0]
        R2min=(mesh.points[SlaveNode][0]-mesh.points[Cnct[SlaveNode]][0])**2+\
        (mesh.points[SlaveNode][1]-mesh.points[Cnct[SlaveNode]][1])**2+\
        (mesh.points[SlaveNode][2]-mesh.points[Cnct[SlaveNode]][2])**2
        for MasterNode in mesh.point_sets[MasterNds]:
            R2=(mesh.points[MasterNode][0]-mesh.points[SlaveNode][0])**2+\
            (mesh.points[MasterNode][1]-mesh.points[SlaveNode][1])**2+\
            (mesh.points[MasterNode][2]-mesh.points[SlaveNode][2])**2
            if R2min>R2:
                R2min=R2
                Cnct[SlaveNode]=MasterNode
    for SlaveNode in mesh.point_sets[SlaveNds]:
        mesh.points[SlaveNode]=None
    for NSet in mesh.point_sets:
        for i in range(len(mesh.point_sets[NSet])):
            Node=mesh.point_sets[NSet][i]
            if Cnct[Node]>0: mesh.point_sets[NSet][i]=Cnct[Node]
    for El in range(1,mesh.MaxElemNum+1):
        if mesh.cells[El]!=1:
            for i in range(len(mesh.cells[El])):
                Node=mesh.cells[El][i]
                if Cnct[Node]>0: mesh.cells[El][i]=Cnct[Node]
