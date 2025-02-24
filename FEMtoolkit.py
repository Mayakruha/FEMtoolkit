import numpy as np
import vtk
#------CONSTANTS------------------
AbaqElemTypes={'MASS':1,'S3':2, 'CAX3':2, 'SPS3':2,'CAX4':3,'CPS4':3, 'CAX6':4, 'CPS6':4, 'CAX8':5, 'SPC8':5, 'S8R':6,\
'C3D4':7,'C3D6':8,'C3D8':9,'C3D10':10,'C3D15':11,'C3D20R':12,'C3D20':12}
CalculiXElemTypes={'9':3}
FacesNodes=(None,None,((0,1),(1,2),(2,0)),((0,1),(1,2),(2,3),(3,0)),((0,1,3),(1,2,4),(2,0,5)),\
((0,1,2),(0,3,1),(1,3,2),(2,3,0)),((0,1,4),(1,2,5),(2,3,6),(3,0,7)),\
((0,1,4),(1,2,5),(2,3,6),(3,0,7)),((0,1,2),(3,5,4),(0,3,4,1),(1,4,5,2),(2,5,3,0)),\
((0,1,2,3),(4,7,6,5),(0,4,5,1),(1,5,6,2),(2,6,7,3),(3,7,4,0)),\
((0,1,2,4,5,6),(0,3,1,7,8,4),(1,3,2,8,9,5),(2,3,0,9,7,6)),\
((0,1,2,6,7,8),(3,5,4,9,10,11),(0,3,4,1,12,9,13,6),(1,4,5,2,13,10,14,7),(2,5,3,0,14,11,12,8)),\
((0,1,2,3,8,9,10,11),(4,7,6,5,12,13,14,15),(0,4,5,1,16,12,17,8),(1,5,6,2,17,13,18,9),(2,6,7,3,18,14,19,10),(3,7,4,0,19,15,16,11)))
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
#------CLASS------------------
class FEMtoolkit:
    def __init__(self):
        self.NSets={}
        self.ESets={}
        self.Surfs={}
        self.MaxNodeNum=0
        self.MaxElemNum=0
        self.Coord=np.full(1,None)
        self.Elems=np.ones(1,dtype=tuple)
        self.Eltype=np.zeros(1,dtype=np.int8)
        self.TypeList={}
        self.NodeValue={} # key: Name of load; key: Node (int); Value
        self.FaceLoad={}
        self.Faces={} # the first key - min mumber of nodes; the second key - max mumber of nodes
                      # list (number of faces, set of numbers of nodes)
#===================================================================
#         Collect info about Faces
#===================================================================
    def EstFaces(self):
        self.Faces={}
        for El in range(1,self.MaxElemNum+1):
            if self.Elems[El]!=1 and self.Eltype[El]>1:
                for Indx in FacesNodes[self.Eltype[El]]:
                    Flag=True
                    Nodes=set()
                    for i in Indx: Nodes.add(self.Elems[El][i])
                    minNode=min(Nodes)
                    maxNode=max(Nodes)
                    if minNode in self.Faces:
                        if maxNode in self.Faces[minNode]:
                            for Fc in self.Faces[minNode][maxNode]:
                                if Fc[1]==Nodes:
                                    Fc[0]+=1
                                    Flag=False
                        else: self.Faces[minNode][maxNode]=[]
                    else:
                        self.Faces[minNode]={}
                        self.Faces[minNode][maxNode]=[]
                    if Flag: self.Faces[minNode][maxNode].append([1,Nodes])
#===================================================================
#         export Abaqus inp-file
#===================================================================
    def export_abq(self,FileName):
        f=open(FileName,'w')
        f.write('*Node\n')
        for i in range(self.MaxNodeNum+1):
            if type(self.Coord[i])==np.ndarray: f.write(str(i)+', '+str(self.Coord[i][0])+', '+str(self.Coord[i][1])+', '+str(self.Coord[i][2])+'\n')
        Num_prev=0
        for i in range(self.MaxElemNum+1):
            if self.Elems[i]!=1:
                Num=len(self.Elems[i])
                if Num_prev!=Num:
                    f.write('*Element, type='+self.TypeList[self.Eltype[i]]+'\n')
                    Num_prev=Num
                f.write(str(i))
                if Num<=15:
                    for Node in self.Elems[i]: f.write(', '+str(Node))
                else:
                    for j in range(15): f.write(', '+str(self.Elems[i][j]))
                    f.write(',\n')
                    f.write(str(self.Elems[i][15]))
                    for j in range(16,Num): f.write(', '+str(self.Elems[i][j]))
                f.write('\n')
        for SetName in self.NSets:
            if len(self.NSets[SetName])!=0:
                f.write('*Nset, nset='+SetName+'\n')
                Count=0
                Num=len(self.NSets[SetName])
                List=self.NSets[SetName].copy()
                for j in range(Num):
                    if Count==0: f.write(str(List.pop()))
                    else: f.write(', '+str(List.pop()))
                    Count+=1
                    if Count==16 and j!=Num-1:
                        f.write('\n')
                        Count=0
                f.write('\n')
        for SetName in self.ESets:
            if len(self.ESets[SetName])!=0:
                f.write('*Elset, elset='+SetName+'\n')
                Count=0
                Num=len(self.ESets[SetName])
                List=self.ESets[SetName].copy()
                for j in range(Num):
                    if Count==0: f.write(str(List.pop()))
                    else: f.write(', '+str(List.pop()))
                    Count+=1
                    if Count==16 and j!=Num-1:
                        f.write('\n')
                        Count=0
                f.write('\n')
        for SetName in self.Surfs:
            f.write('*Surface, type=ELEMENT, name='+SetName+'\n')
            for Face in self.Surfs[SetName]: f.write(Face[0]+', S'+str(Face[1]+1)+'\n')
        f.close()
#===================================================================
#         import Node load
#===================================================================
    def import_ndload(self,FileName,LoadName):
        f=open(FileName,'r')
        self.NodeValue[LoadName]={}
        txt=f.readline()
        while txt:
            Values=txt.split(',')
            Node=int(Values[0])
            if Node<=self.MaxNodeNum: self.NodeValue[LoadName][Node]=float(Values[1])
            txt=f.readline()
        f.close()
#===================================================================
#         export Node load
#===================================================================
    def export_ndload(self,FileName,LoadName,separator=','):
        f=open(FileName,'w')
        for Node in range(1,self.MaxNodeNum+1):
            if (Node in self.NodeValue[LoadName])and(type(self.Coord[Node])==np.ndarray):
                f.write(str(Node)+separator+str(self.NodeValue[LoadName][Node])+'\n')
        f.close()
#===================================================================
#         import rpt
#===================================================================
    def read_rpt(self, file_name):
        Clmns = {}
        Flag = False
        RowLen = 0
        Alg = 0
        rpt_f = open(file_name,'r')
        rawline = rpt_f.readline()
        while rawline and Alg == 0:
            if 'EXTRAPOLATE_AVERAGE_COMPUTE' in rawline: Alg=1
            rawline=rpt_f.readline()
        #-----------AVERAGE NODAL-------------------------------
        if Alg==1:
            while rawline:
                txt_len=len(rawline)
                if Flag and txt_len!=RowLen: Flag=False
                if 'Node' in rawline:
                    count=0
                    while count!=-1:
                        label_st=txt_len-len(rawline[count:].lstrip())
                        label_en=rawline.find(' ', label_st)
                        if rawline[label_st:label_en]=='Node':
                            node_clmn=(count,label_en)
                        else:
                            Clmns[rawline[label_st:label_en]]=(count,label_en)
                            self.NodeValue[rawline[label_st:label_en]]={}
                        count=label_en
                    rpt_f.readline()
                    rpt_f.readline()
                    Flag=True
                    RowLen=txt_len
                elif Flag:
                    Node=int(rawline[node_clmn[0]:node_clmn[1]])
                    for Param in Clmns:
                        self.NodeValue[Param][Node]=float(rawline[Clmns[Param][0]:Clmns[Param][1]])
                rawline=rpt_f.readline()
        rpt_f.close()
        print('The rpt-file has been imported')
#===================================================================
#         import Face load
# LoadType: 'P' - Pressure; 'S' - Heat Flux; 'F' - HTC
#===================================================================
    def import_fcload(self,FileName,LoadType):
        f=open(FileName,'r')
        if not LoadType in self.FaceLoad: self.FaceLoad[LoadType]={}
        txt=f.readline()
        while txt:
            Values=txt.split(',')
            El=int(Values[0])
            if not El in self.FaceLoad[LoadType]:self.FaceLoad[LoadType][int(Values[0])]=[]
            Val=[int(Values[1][-1:]),float(Values[2])]
            if len(Values)>3:Val.append(float(Values[3]))
            self.FaceLoad[LoadType][int(Values[0])].append(Val)
            txt=f.readline()
        f.close()
#===================================================================
#         export Face load
#===================================================================
    def export_fcload(self,FileName,LoadName):
        f=open(FileName,'w')
        for El in range(1,self.MaxElemNum+1):
            if (El in self.FaceLoad[LoadName])and(self.Elems[El]!=1):
                for Fc in self.FaceLoad[LoadName][El]:
                    f.write(str(El)+', '+LoadName+str(Fc[0]))
                    for i in range(1,len(Fc)):f.write(', '+str(Fc[i]))
                    f.write('\n')
        f.close()
#===================================================================
#         Node set -> Surface
#===================================================================
    def NodeIntoSurf(self,NSet):
        self.Surfs[NSet]=[]
        for i in range(1, self.MaxElemNum+1):
            if self.Elems[i]!=1:
                for FaceIndx in range(len(FacesNodes[self.Eltype[i]])):
                    Flag=True
                    for NdIndx in FacesNodes[self.Eltype[i]][FaceIndx]:
                        if not self.Elems[i][NdIndx] in self.NSets[NSet]: Flag=False
                    if Flag:
                        SetFaceName=NSet+'_S'+str(FaceIndx+1)
                        if not (SetFaceName,FaceIndx) in self.Surfs[NSet]: self.Surfs[NSet].append((SetFaceName,FaceIndx))
                        if not SetFaceName in self.ESets: self.ESets[SetFaceName]=[]
                        self.ESets[SetFaceName].append(i)
#===================================================================
#         Changes nodal load locally
#
# ValueName - Name of Nodal load for the changes
# dValue    - Maximum value of the changes
# AreaType  - Type of distribution ('Plane',)
# Size      - Size of area of the changes
# Point     - Coordinates of a point
# Vector    - Coordinates of a vector
#===================================================================
    def ChangeField(self, ValueName, dValue, AreaType, Size, Point, Vector):
        VecLen=(Vector[0]**2+Vector[1]**2+Vector[2]**2)**0.5
        #===================================================================
        # AreaType  - 'Plane'
        # Size      - Maximal distance of affect from the plane
        # Point     - Coordinates of any point on the plane
        # Vector    - a normal ti the plane
        #===================================================================
        if AreaType=='Plane':
            for Node in self.NodeValue[ValueName]:
                Sum=0
                for i in range(3):
                    Sum+=(self.Coord[Node][i]-Point[i])*Vector[i]/VecLen
                R=abs(Sum)
                if R<Size:
                    self.NodeValue[ValueName][Node]+=dValue*(1-(R/Size)**2)
#===================================================================
#         caling of coordinates
#===================================================================
    def Scale(self,Scale):
        for i in range(1,self.MaxNodeNum+1):
            if type(self.Coord[i])==np.ndarray:
                for j in range(3):
                    self.Coord[i][j]*=Scale
#===================================================================
#
#         write command file with edges to apply in Abaqus
#
# Variables:
# FileName - File Name for a Python file to be applied in Abaqus
# NSet - Name of a set of nodes applied to the edges
# ESet - Name of a set of elements applied to the edges
# Scale - scale for coordinates to be applied for output
# PrjctMtrx - Matrix that transfers a global coordinate vector to a local coordinate vector of the sketch
#===================================================================
    def WriteEdges(self,FileName,NSet,ESet, Scale=1, PrjctMtrx=[[1,0,0],[0,1,0]]):
        Edges={}
        f=open(FileName,'w')
        f.write('#Select the following: \n')
        f.write('#-Partition face: Sketch \n')
        f.write('#-face, "Specify" for Sketch Origin \n')
        f.write('#-0,0,0\n')
        f.write('####### Edit Model Name and copy commands from the file to the command line')
        f.write('ModelName="Model-1"\n')
        for ElemNum in self.ESets[ESet]:
            Nodelist=self.Elems[ElemNum].copy()
            Nodelist.sort()
            Num=len(Nodelist)
            for i in range(Num-1):
                for j in Nodelist[i+1:Num]:
                    print(str(Nodelist[i])+' '+str(j))
                    if Nodelist[i] in self.NSets[NSet] and j in self.NSets[NSet]:
                        if not Nodelist[i] in Edges: Edges[Nodelist[i]]={}
                        if not j in Edges[Nodelist[i]]: Edges[Nodelist[i]][j]=0
                        Edges[Nodelist[i]][j]+=1
        for Node1 in Edges:
            for Node2 in Edges[Node1]:
                if Edges[Node1][Node2]>1:
                    Point1=np.dot(PrjctMtrx,np.array(self.Coord[Node1]))
                    Point2=np.dot(PrjctMtrx,np.array(self.Coord[Node2]))
                    f.write('mdb.models[ModelName].sketches["__profile__"].Line(point1=('+str(Point1[0]*Scale)+','+str(Point1[1]*Scale)+'),\
point2=('+str(Point2[0]*Scale)+','+str(Point2[1]*Scale)+'))\n')
        f.close()
#===================================================================
#
#         Extract thickness of coating that is simulated by linear triangular prism
#
# Variables:
# FileName - File Name of a vtu-file for output
#           vtkXMLUnstructuredGridReader should be used for reading
#           Reader=vtk.vtkXMLUnstructuredGridReader()
#           Reader.SetFileName('File.vtu')
#           Reader.Update()
#           mesh=Reader.GetOutput()
# NsetName - Name of a set of nodes on internal surface (surface betwen base material and coating)
# EsetNames - Array of Names of a set of elements for layers of coating
#===================================================================
    def ExtractCoating(self,FileName,NsetName,EsetNames):
        ElemRef=np.full(self.MaxElemNum+1,None,dtype=tuple)
        StackElem={}
        for i in range(1,len(EsetNames)):
            StackElem[EsetNames[i]]=self.ESets[EsetNames[i]].copy()
        Points=vtk.vtkPoints()
        mesh=vtk.vtkUnstructuredGrid()
        ElemLen=len(self.ESets[EsetNames[0]])
        i=0
        for ElemNum in self.ESets[EsetNames[0]]: 
            ElemExist=False
            if self.Elems[ElemNum][0] in self.NSets[NsetName] and self.Elems[ElemNum][1] in self.NSets[NsetName] and self.Elems[ElemNum][2] in self.NSets[NsetName]:
                for j in range(3):Points.InsertNextPoint(self.Coord[self.Elems[ElemNum][j]][0],self.Coord[self.Elems[ElemNum][j]][1],self.Coord[self.Elems[ElemNum][j]][2])
                ElemRef[ElemNum]=(i*3,i*3+1,i*3+2,0)
                Nodes0=(self.Elems[ElemNum][3],self.Elems[ElemNum][4],self.Elems[ElemNum][5])
                ElemExist=True
            elif self.Elems[ElemNum][3] in self.NSets[NsetName] and self.Elems[ElemNum][4] in self.NSets[NsetName] and self.Elems[ElemNum][5] in self.NSets[NsetName]:
                for j in range(3):Points.InsertNextPoint(self.Coord[self.Elems[ElemNum][j+3]][0],self.Coord[self.Elems[ElemNum][j+3]][1],self.Coord[self.Elems[ElemNum][j+3]][2])
                ElemRef[ElemNum]=(i*3,i*3+1,i*3+2,1)
                Nodes0=(self.Elems[ElemNum][0],self.Elems[ElemNum][1],self.Elems[ElemNum][2])
                ElemExist=True
            else:
                print('The nodes havent been found in '+NsetName+' for element '+str(ElemNum))
            for j in range(len(StackElem)):
                EName=EsetNames[j+1]
                Flag=len(StackElem[EName])>0
                i_el=0
                while Flag:
                    Nodes=set(self.Elems[StackElem[EName][i_el]])
                    if len(Nodes.intersection(set(Nodes0)))==3:
                        if self.Elems[StackElem[EName][i_el]][0]==Nodes0[0]:
                            ElemRef[StackElem[EName][i_el]]=(i*3,i*3+1,i*3+2,0)
                            Nodes0=(self.Elems[ElemNum][3],self.Elems[ElemNum][4],self.Elems[ElemNum][5])
                        elif self.Elems[StackElem[EName][i_el]][3]==Nodes0[0]:
                            ElemRef[StackElem[EName][i_el]]=(i*3,i*3+1,i*3+2,1)
                            Nodes0=(self.Elems[ElemNum][0],self.Elems[ElemNum][1],self.Elems[ElemNum][2])
                        elif self.Elems[StackElem[EName][i_el]][1]==Nodes0[0]:
                            ElemRef[StackElem[EName][i_el]]=(i*3+1,i*3+2,i*3,0)
                            Nodes0=(self.Elems[ElemNum][4],self.Elems[ElemNum][5],self.Elems[ElemNum][3])
                        elif self.Elems[StackElem[EName][i_el]][4]==Nodes0[0]:
                            ElemRef[StackElem[EName][i_el]]=(i*3+1,i*3+2,i*3,1)
                            Nodes0=(self.Elems[ElemNum][1],self.Elems[ElemNum][2],self.Elems[ElemNum][0])
                        elif self.Elems[StackElem[EName][i_el]][2]==Nodes0[0]:
                            ElemRef[StackElem[EName][i_el]]=(i*3+2,i*3,i*3+1,0)
                            Nodes0=(self.Elems[ElemNum][5],self.Elems[ElemNum][3],self.Elems[ElemNum][2])
                        elif self.Elems[StackElem[EName][i_el]][5]==Nodes0[0]:
                            ElemRef[StackElem[EName][i_el]]=(i*3+2,i*3,i*3+1,1)
                            Nodes0=(self.Elems[ElemNum][2],self.Elems[ElemNum][0],self.Elems[ElemNum][1])
                        StackElem[EName].remove(StackElem[EName][i_el])
                        Flag=False
                    i_el+=1
                    if i_el==len(StackElem[EName]):Flag=False
            if ElemExist:i+=1
        ElemLen=i
        mesh.Allocate(ElemLen)
        for i in range(ElemLen): mesh.InsertNextCell(vtk.VTK_TRIANGLE,3,(i*3,i*3+1,i*3+2))
        mesh.SetPoints(Points)
        #============= thickness analysis
        Thick=[]
        for i in range(len(EsetNames)):
            EName=EsetNames[i]
            Thick.append(vtk.vtkFloatArray())
            Thick[i].SetName(EName)
            Thick[i].SetNumberOfValues(ElemLen*3)
            Thick[i].Fill(0)
            for ENum in self.ESets[EName]:
                if ElemRef[ENum]!=None:
                    Vb1=np.array((self.Coord[self.Elems[ENum][1]][0]-self.Coord[self.Elems[ENum][0]][0],self.Coord[self.Elems[ENum][1]][1]-self.Coord[self.Elems[ENum][0]][1],self.Coord[self.Elems[ENum][1]][2]-self.Coord[self.Elems[ENum][0]][2]))
                    Vb2=np.array((self.Coord[self.Elems[ENum][2]][0]-self.Coord[self.Elems[ENum][0]][0],self.Coord[self.Elems[ENum][2]][1]-self.Coord[self.Elems[ENum][0]][1],self.Coord[self.Elems[ENum][2]][2]-self.Coord[self.Elems[ENum][0]][2]))
                    NormB=np.cross(Vb1, Vb2)
                    Vt1=np.array((self.Coord[self.Elems[ENum][4]][0]-self.Coord[self.Elems[ENum][3]][0],self.Coord[self.Elems[ENum][4]][1]-self.Coord[self.Elems[ENum][3]][1],self.Coord[self.Elems[ENum][4]][2]-self.Coord[self.Elems[ENum][3]][2]))
                    Vt2=np.array((self.Coord[self.Elems[ENum][5]][0]-self.Coord[self.Elems[ENum][3]][0],self.Coord[self.Elems[ENum][5]][1]-self.Coord[self.Elems[ENum][3]][1],self.Coord[self.Elems[ENum][5]][2]-self.Coord[self.Elems[ENum][3]][2]))
                    NormT=np.cross(Vt1, Vt2)
                    if ElemRef[ENum][3]==0:
                        NormB=NormB/np.linalg.norm(NormB)
                        for j in range(0,3):
                            Vec=np.array((self.Coord[self.Elems[ENum][3]][0]-self.Coord[self.Elems[ENum][j]][0],self.Coord[self.Elems[ENum][3]][1]-self.Coord[self.Elems[ENum][j]][1],self.Coord[self.Elems[ENum][3]][2]-self.Coord[self.Elems[ENum][j]][2]))
                            Dist=abs(np.dot(NormT,Vec)/np.dot(NormB,NormT))
                            Thick[i].SetValue(ElemRef[ENum][j],Dist)
                    elif ElemRef[ENum][3]==1:
                        NormT=NormT/np.linalg.norm(NormT)  
                        for j in range(3,6):
                            Vec=np.array((self.Coord[self.Elems[ENum][0]][0]-self.Coord[self.Elems[ENum][j]][0],self.Coord[self.Elems[ENum][0]][1]-self.Coord[self.Elems[ENum][j]][1],self.Coord[self.Elems[ENum][0]][2]-self.Coord[self.Elems[ENum][j]][2]))
                            Dist=abs(np.dot(NormB,Vec)/np.dot(NormB,NormT))
                            Thick[i].SetValue(ElemRef[ENum][j-3],Dist)
        for i in range(len(EsetNames)): mesh.GetPointData().AddArray(Thick[i])
        #============= output
        output=vtk.vtkXMLUnstructuredGridWriter()
        output.SetInputData(mesh)
        output.SetFileName(FileName)
        output.Write()
        return mesh
#===================================================================
#
#         Create coating
#
# Variables:
# ThickNames - Names of NodeValue sets with thickness
# NSet       - Set of nodes for coating
# Inside     - If False the function creates coating
#===================================================================
    def CreateLayers(self,ThickNames,NSet,Inside=False):
        self.TypeList[5]='C3D6'
        N=len(self.NSets[NSet])
        LayerNum=len(ThickNames)
        thickness=np.zeros((LayerNum,N))
        vertices = np.zeros((LayerNum+1,N,3))
        #-------------
        List=np.full(self.MaxNodeNum+1,-1,dtype=np.int32)
        for i in range(N):
            Node=self.NSets[NSet][i]
            List[Node]=i
            vertices[0][i][0]=self.Coord[Node][0]
            vertices[0][i][1]=self.Coord[Node][1]
            vertices[0][i][2]=self.Coord[Node][2]
            for j in range(LayerNum):thickness[j][i]=self.NodeValue[ThickNames[j]][Node]
#-----Reading faces
        Faces_dyn=[]
        for i in range(self.MaxElemNum+1):
            if self.Elems[i]!=1:
                Node0=self.Elems[i][0]
                Node1=self.Elems[i][1]
                Node2=self.Elems[i][2]
                Node3=self.Elems[i][3]
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
                Node=self.NSets[NSet][i]
                self.Coord[Node][0]=vertices[0][i][0]
                self.Coord[Node][1]=vertices[0][i][1]
                self.Coord[Node][2]=vertices[0][i][2]
        for j in range(LayerNum):
            vertices[j+1,:,0]=norm[:,0] * thickness[j] + vertices[j,:,0]
            vertices[j+1,:,1]=norm[:,1] * thickness[j] + vertices[j,:,1]
            vertices[j+1,:,2]=norm[:,2] * thickness[j] + vertices[j,:,2]
        #------------output-----------------------
        NodesNum=np.zeros((LayerNum+1,FacesNum,3),dtype=np.int64)
        for i in range(FacesNum):
            NodesNum[0][i]=[self.NSets[NSet][faces[i,0]],self.NSets[NSet][faces[i,1]],self.NSets[NSet][faces[i,2]]]
        self.Coord=np.resize(self.Coord,self.MaxNodeNum+N*LayerNum+1)
        self.Elems=np.resize(self.Elems,self.MaxElemNum+FacesNum*LayerNum+1)
        self.Eltype=np.resize(self.Eltype,self.MaxElemNum+FacesNum*LayerNum+1)
        for j in range(LayerNum):
            for i in range(N):    
                if thickness[j,i]>0:
                    self.Coord[self.MaxNodeNum+1+j*N+i]=np.array((vertices[j+1][i][0],vertices[j+1][i][1],vertices[j+1][i][2]))
            self.ESets[ThickNames[j]]=[]
            NodesNum[j+1]=self.MaxNodeNum+1+j*N+faces
            for i in range(FacesNum):
                self.ESets[ThickNames[j]].append(self.MaxElemNum+1+FacesNum*j+i)
                self.Elems[self.MaxElemNum+1+FacesNum*j+i]=list()
                self.Eltype[self.MaxElemNum+1+FacesNum*j+i]=5
                for Node in NodesNum[j][i]: self.Elems[self.MaxElemNum+1+FacesNum*j+i]+=(Node,)
                for Node in NodesNum[j+1][i]: self.Elems[self.MaxElemNum+1+FacesNum*j+i]+=(Node,)
        self.MaxNodeNum+=N*LayerNum
        self.MaxElemNum+=FacesNum*LayerNum
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
    def SymmetryEquations(self,FileName,NSet1,NSet2,method='Tolerance',tolerance=(0.00001,0.00001)):
        if not NSet1 in self.NSets:
            print(NSet1+' hasnt been found')
            return
        else: print(NSet1+' contains '+str(len(self.NSets[NSet1])))
        if not NSet2 in self.NSets:
            print(NSet2+' hasnt been found')
            return        
        else: print(NSet2+' contains '+str(len(self.NSets[NSet2])))
        Sym1=[]
        Sym2=self.NSets[NSet2].copy()
        #------------------------------------------
        #--------------NSET------------------------
        f=open(FileName,'w')
        f.write('*Nset, nset=SYMNODES\n')
        Count=0        
        for i in range(len(self.NSets[NSet1])):
            f.write(str(self.NSets[NSet1][i]))
            if Count<15:
                f.write(',')
                Count+=1
            else:
                f.write('\n')
                Count=0
        Num=len(self.NSets[NSet2])
        for i in range(len(self.NSets[NSet2])):
            f.write(str(self.NSets[NSet2][i]))
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
            for Node1 in self.NSets[NSet1]:
                R1=(self.Coord[Node1][1]**2+self.Coord[Node1][2]**2)**0.5
                Flag=True
                i=0
                while Flag:
                    Node2=Sym2[i]
                    R2=(self.Coord[Node2][1]**2+self.Coord[Node2][2]**2)**0.5
                    if abs(R1-R2)<tolerance[0] and abs(self.Coord[Node1][0]-self.Coord[Node2][0])<tolerance[1]:
                        f.write('2\n')
                        f.write(str(Node1)+',1,-1,'+str(Node2)+',1,1\n')
                        f.write('2\n')
                        f.write(str(Node1)+',2,-1,'+str(Node2)+',2,1\n')
                        f.write('2\n')
                        f.write(str(Node1)+',3,-1,'+str(Node2)+',3,1\n')
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
            for Node1 in self.NSets[NSet1]:
                R1=(self.Coord[Node1][1]**2+self.Coord[Node1][2]**2)**0.5
                for i in range(len(Sym2)):
                    Node2=Sym2[i]
                    R2=(self.Coord[Node2][1]**2+self.Coord[Node2][2]**2)**0.5
                    Dist=((R1-R2)**2+(self.Coord[Node1][0]-self.Coord[Node2][0])**2)**0.5
                    if i==0 or DistMin>Dist:
                        DistMin=Dist
                        Node2min=Node2
                f.write('2\n')
                f.write(str(Node1)+',1,-1,'+str(Node2min)+',1,1\n')
                f.write('2\n')
                f.write(str(Node1)+',2,-1,'+str(Node2min)+',2,1\n')
                f.write('2\n')
                f.write(str(Node1)+',3,-1,'+str(Node2min)+',3,1\n')
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
                f.write(str(Sym1[i]))
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
                f.write(str(Sym2[i]))
                if Count<15 and i<Num-1:
                    f.write(',')
                    Count+=1
                else:
                    f.write('\n')
                    Count=0
        f.close()
#===================================================================
#
#         Create a vtu-file with values for nodes defined in NodeValue
#
# Variables:
# FileName   - Name of a vtu-file (vtkXMLUnstructuredGridReader with vtkFloatArrays)
# Eset       - Name of a set of elements 
#===================================================================
    def LinearVTU(self,FileName,Eset):
        Nums=np.full(len(self.Coord),self.MaxNodeNum+1,dtype=np.int32)
        Points=vtk.vtkPoints()
        mesh=vtk.vtkUnstructuredGrid()
        Node_count=0
        for ElemNum in self.ESets[Eset]:
            Nodelist=[]
            for Node in self.Elems[ElemNum]:
                if Nums[Node]==self.MaxNodeNum+1:
                    Nums[Node]=Node_count
                    Points.InsertNextPoint(self.Coord[Node][0],self.Coord[Node][1],self.Coord[Node][2])
                    Node_count+=1
                Nodelist.append(Nums[Node])
            if self.Eltype[ElemNum]==3:#linear quad
                mesh.InsertNextCell(vtk.VTK_TRIANGLE,3,(Nodelist[0],Nodelist[1],Nodelist[2]))
                mesh.InsertNextCell(vtk.VTK_TRIANGLE,3,(Nodelist[2],Nodelist[3],Nodelist[0]))
            elif self.Eltype[ElemNum]==4: mesh.InsertNextCell(vtk.VTK_TETRA,4,Nodelist)
            elif self.Eltype[ElemNum]==5:#linear wedge
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[0],Nodelist[1],Nodelist[3],Nodelist[2]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[1],Nodelist[4],Nodelist[3],Nodelist[2]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[3],Nodelist[2],Nodelist[4],Nodelist[5]))
            elif self.Eltype[ElemNum]==6:#linear hexahedron
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[0],Nodelist[1],Nodelist[3],Nodelist[4]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[1],Nodelist[2],Nodelist[3],Nodelist[4]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[3],Nodelist[4],Nodelist[2],Nodelist[7]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[5],Nodelist[4],Nodelist[6],Nodelist[1]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[4],Nodelist[7],Nodelist[6],Nodelist[1]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[6],Nodelist[1],Nodelist[7],Nodelist[2]))
            elif self.Eltype[ElemNum]==7:#quadratic tetra
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[0],Nodelist[4],Nodelist[6],Nodelist[7]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[4],Nodelist[1],Nodelist[5],Nodelist[8]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[5],Nodelist[2],Nodelist[6],Nodelist[9]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[7],Nodelist[8],Nodelist[9],Nodelist[3]))                
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[6],Nodelist[4],Nodelist[5],Nodelist[7]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[4],Nodelist[8],Nodelist[5],Nodelist[7]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[5],Nodelist[8],Nodelist[9],Nodelist[7]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[5],Nodelist[9],Nodelist[6],Nodelist[7]))
            elif self.Eltype[ElemNum]==8:#quadratic wedge
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[0],Nodelist[6],Nodelist[12],Nodelist[8]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[1],Nodelist[13],Nodelist[6],Nodelist[7]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[3],Nodelist[12],Nodelist[9],Nodelist[11]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[4],Nodelist[9],Nodelist[13],Nodelist[10]))                
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[11],Nodelist[12],Nodelist[9],Nodelist[5]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[10],Nodelist[9],Nodelist[13],Nodelist[5]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[12],Nodelist[8],Nodelist[6],Nodelist[2]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[6],Nodelist[7],Nodelist[13],Nodelist[2]))                
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[9],Nodelist[12],Nodelist[13],Nodelist[5]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[12],Nodelist[6],Nodelist[13],Nodelist[2]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[5],Nodelist[2],Nodelist[12],Nodelist[13]))               
            elif self.Eltype[ElemNum]==9:#quadratic hexahedron
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[0],Nodelist[8],Nodelist[11],Nodelist[16]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[1],Nodelist[9],Nodelist[8],Nodelist[17]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[2],Nodelist[10],Nodelist[9],Nodelist[18]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[3],Nodelist[11],Nodelist[10],Nodelist[19]))               
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[4],Nodelist[15],Nodelist[12],Nodelist[16]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[5],Nodelist[12],Nodelist[13],Nodelist[17]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[6],Nodelist[13],Nodelist[14],Nodelist[18]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[7],Nodelist[15],Nodelist[14],Nodelist[19]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[11],Nodelist[12],Nodelist[15],Nodelist[16]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[11],Nodelist[8],Nodelist[12],Nodelist[16]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[12],Nodelist[8],Nodelist[9],Nodelist[17]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[9],Nodelist[13],Nodelist[12],Nodelist[17]))                
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[9],Nodelist[10],Nodelist[13],Nodelist[18]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[10],Nodelist[14],Nodelist[13],Nodelist[18]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[10],Nodelist[14],Nodelist[11],Nodelist[19]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[14],Nodelist[15],Nodelist[11],Nodelist[19]))               
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[11],Nodelist[12],Nodelist[14],Nodelist[15]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[14],Nodelist[12],Nodelist[9],Nodelist[13]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[11],Nodelist[8],Nodelist[9],Nodelist[12]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[11],Nodelist[9],Nodelist[10],Nodelist[14]))
                mesh.InsertNextCell(vtk.VTK_TETRA,4,(Nodelist[14],Nodelist[12],Nodelist[11],Nodelist[9]))
        mesh.SetPoints(Points)
        #------------Field-------------------------------
        for i in range(len(self.NodeValue)):
            ValueName=list(self.NodeValue.keys())[i]            
            Values=vtk.vtkFloatArray()
            Values.SetName(ValueName)
            Values.SetNumberOfValues(Node_count)
            for NodeGl in self.NodeValue[ValueName]:
                Node=Nums[NodeGl]
                if Node<self.MaxNodeNum+1: Values.SetValue(Node,self.NodeValue[ValueName][NodeGl])
            mesh.GetPointData().AddArray(Values)
        output=vtk.vtkXMLUnstructuredGridWriter()
        output.SetInputData(mesh)
        output.SetFileName(FileName)
        output.Write()
#===================================================================
#
#         Create a vtu-file with values for nodes defined in NodeValue
#
# Variables:
# FileName - Name of a vtu-file (vtkXMLUnstructuredGridReader with vtkFloatArrays)
# Eset     - Name of a set of elements
# Id       - If true then Id numbers are in output
# Scalars  - List of Names of NodeValue for Output
# Vectors  - dictionary of Vectors for Output {Name: Vector - list of three Names of NodeValue}
#===================================================================
    def export_VTU(self,FileName,Eset,Scalars=[],Vectors={},Id=False):
        Nums=np.full(len(self.Coord),self.MaxNodeNum+1,dtype=np.int32)
        Points=vtk.vtkPoints()
        mesh=vtk.vtkUnstructuredGrid()
        Node_count=0
        for ElemNum in self.ESets[Eset]:
            Nodelist=[]
            for Node in self.Elems[ElemNum]:
                if Nums[Node]==self.MaxNodeNum+1:
                    Nums[Node]=Node_count
                    Points.InsertNextPoint(self.Coord[Node][0],self.Coord[Node][1],self.Coord[Node][2])
                    Node_count+=1
                Nodelist.append(Nums[Node])
            if self.Eltype[ElemNum]==3: mesh.InsertNextCell(vtk.VTK_QUAD,4,Nodelist)
            elif self.Eltype[ElemNum]==4: mesh.InsertNextCell(vtk.VTK_TETRA,4,Nodelist)
            elif self.Eltype[ElemNum]==5: mesh.InsertNextCell(vtk.VTK_WEDGE,6,Nodelist)
            elif self.Eltype[ElemNum]==6: mesh.InsertNextCell(vtk.VTK_HEXAHEDRON,8,Nodelist)
            elif self.Eltype[ElemNum]==7: mesh.InsertNextCell(vtk.VTK_QUADRATIC_TETRA,10,Nodelist)
            elif self.Eltype[ElemNum]==8: mesh.InsertNextCell(vtk.VTK_QUADRATIC_WEDGE,15,Nodelist)
            elif self.Eltype[ElemNum]==9: mesh.InsertNextCell(vtk.VTK_QUADRATIC_HEXAHEDRON,20,Nodelist)
        mesh.SetPoints(Points)
        #------------Scalars-------------------------------
        for i in range(len(Scalars)):
            ValueName=Scalars[i]            
            Values=vtk.vtkFloatArray()
            Values.SetName(ValueName)
            Values.SetNumberOfValues(Node_count)
            for NodeGl in self.NodeValue[ValueName]:
                Node=Nums[NodeGl]
                if Node<self.MaxNodeNum+1: Values.SetValue(Node,self.NodeValue[ValueName][NodeGl])
            mesh.GetPointData().AddArray(Values)
        #------------Vectors-------------------------------
        for VecName in Vectors:
            Values=vtk.vtkFloatArray()
            Values.SetName(VecName)
            Values.SetNumberOfComponents(3)
            Values.SetNumberOfTuples(Node_count)
            for NodeGl in self.NodeValue[ValueName]:
                Node=Nums[NodeGl]
                if Node<self.MaxNodeNum+1:
                    ValX=self.NodeValue[Vectors[VecName][0]][NodeGl]
                    ValY=self.NodeValue[Vectors[VecName][1]][NodeGl]
                    ValZ=self.NodeValue[Vectors[VecName][2]][NodeGl]
                    Values.SetTuple3(Node,ValX,ValY,ValZ)
            mesh.GetPointData().AddArray(Values)
        if Id:
            #----Node------
            Values=vtk.vtkFloatArray()
            Values.SetName('Node Id')
            Values.SetNumberOfValues(Node_count)
            for NodeGl in range(1,self.MaxNodeNum+1):
                Node=Nums[NodeGl]
                if Node<self.MaxNodeNum+1: Values.SetValue(Node,NodeGl)
            mesh.GetPointData().AddArray(Values)
            #----Element-----
            Values=vtk.vtkFloatArray()
            Values.SetName('Elem Id')
            Values.SetNumberOfValues(mesh.GetNumberOfCells())
            for i in range(len(self.ESets[Eset])):
                ElemNum=self.ESets[Eset][i]
                if self.Eltype[ElemNum]>1: Values.SetValue(i,ElemNum)
            mesh.GetCellData().SetScalars(Values)
        output=vtk.vtkXMLUnstructuredGridWriter()
        output.SetInputData(mesh)
        output.SetFileName(FileName)
        output.Write()
#===================================================================
#
#         Volume mapping by means of a vtu-file
#
# Variables:
# FileName    - Name of a vtu-file with field (vtkXMLUnstructuredGridReader with vtkFloatArrays)
# NodeSet     - Name of a set of nodes for 
# Tlrnc       - Tolerance. If relative distance is more than Tlrnc, the minimum distance method is used
#===================================================================
    def mapping(self,FileName,NodeSet,Tlrnc=0.005):
        Reader=vtk.vtkXMLUnstructuredGridReader()
        Reader.SetFileName(FileName)
        Reader.Update()
        vtkData=Reader.GetOutput()
        FieldNum=vtkData.GetPointData().GetNumberOfArrays()
        for i in range(FieldNum):
            self.NodeValue[vtkData.GetPointData().GetArrayName(i)]={}
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
            if Xmin>self.Coord[Node][0]:Xmin=self.Coord[Node][0]
            if Xmax<self.Coord[Node][0]:Xmax=self.Coord[Node][0]
            if Ymin>self.Coord[Node][1]:Ymin=self.Coord[Node][1]
            if Ymax<self.Coord[Node][1]:Ymax=self.Coord[Node][1]
            if Zmin>self.Coord[Node][2]:Zmin=self.Coord[Node][2]
            if Zmax<self.Coord[Node][2]:Zmax=self.Coord[Node][2]
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
        for Node in self.NSets[NodeSet]:
            GlPoint=np.array((self.Coord[Node][0],self.Coord[Node][1],self.Coord[Node][2]))
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
                self.NodeValue[vtkData.GetPointData().GetArrayName(FN)][Node]=Value[FN]
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
                        GlPoint=np.array((self.Coord[Node][0],self.Coord[Node][1],self.Coord[Node][2]))
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
                            self.NodeValue[vtkData.GetPointData().GetArrayName(FN)][Node]=Value[FN]
#----------------------------------------------
        if len(MinDistNodes)>0:
            print('The nearest method was applied to '+str(len(MinDistNodes))+' nodes')
            print('See MinDistNodes list')
            self.NSets['MinDistNodes']=MinDistNodes
#===================================================================
#
#         Mapping from surface data
#
# Variables:
# FileName - Name of a vtu-file (vtkXMLUnstructuredGridReader with vtkFloatArrays)
# NodeSet - Name of a set of nodes for mapping
# DistError - Distance error (just for messaging)
#===================================================================
    def map_surf(self,FileName,NodeSet,DistError=0.0001):
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
            self.NodeValue[vtkSurfdData.GetPointData().GetArray(j).GetName()]={}
        self.NSets['NodesOutOfTolerance']=[]
        for Node in self.NSets[NodeSet]:
            GlPoint=np.array((self.Coord[Node][0],self.Coord[Node][1],self.Coord[Node][2]))
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
            if TolFlag: self.NSets['NodesOutOfTolerance'].append(Node)
            for j in range(vtkSurfdData.GetPointData().GetNumberOfArrays()):
                for k in range(3):
                    CellNode=vtkSurfdData.GetCell(i_Cell).GetPointIds().GetId(k)
                    V1[k]=vtkSurfdData.GetPointData().GetArray(j).GetValue(CellNode)
                Value=V1[0]+(V1[1]-V1[0])*Ksi+(V1[2]-V1[0])*Nu
                self.NodeValue[vtkSurfdData.GetPointData().GetArray(j).GetName()][Node]=Value
        if len(self.NSets['NodesOutOfTolerance'])>0:
            print(str(len(self.NSets['NodesOutOfTolerance']))+' nodes are out of tolerance')
            print('Look at "NodesOutOfTolerance" node set')
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
    def map_edge(self,FileName,SurfName,i_coord,LoadType,separator=';'):
        if not LoadType in self.FaceLoad: self.FaceLoad[LoadType]={}
        f=open(FileName,'r')
        Values0=list(map(float,f.readline().split(separator)))
        ValN=len(Values0)
        txt=f.readline()
        while txt:
            Values=list(map(float,txt.split(separator)))
            for SSet in self.Surfs[SurfName]:
                for El in self.ESets[SSet[0]]:
                    Coord=0
                    for Nd_i in FacesNodes[self.Eltype[El]][SSet[1]]:
                        Coord+=self.Coord[self.Elems[El][Nd_i]][i_coord]
                    Coord/=len(FacesNodes[self.Eltype[El]][SSet[1]])
                    if Coord>=Values0[0] and Coord<Values[0]:
                        if not El in self.FaceLoad[LoadType]: self.FaceLoad[LoadType][El]=[]
                        Val=[SSet[1]+1,]
                        for i in range(1,ValN):
                            Val.append(Values0[i]+(Values[i]-Values0[i])/(Values[0]-Values0[0])*(Coord-Values0[0]))
                        self.FaceLoad[LoadType][El].append(Val)
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
    def crack(self,NSet,ESet):
        NodeList={}
        for Node in self.NSets[NSet]:
            NodeList[Node]=0
        for i in range(self.MaxElemNum+1):
            if self.Elems[i]!=1:
                if not i in self.ESets[ESet]:
                    for j in range(len(self.Elems[i])):
                        if self.Elems[i][j] in self.NSets[NSet]:
                            if NodeList[self.Elems[i][j]]==0:
                                self.MaxNodeNum+=1
                                NodeList[self.Elems[i][j]]=self.MaxNodeNum                                
                            self.Elems[i][j]=NodeList[self.Elems[i][j]]
        self.Coord=np.resize(self.Coord,self.MaxNodeNum+1)
        for Node in self.NSets[NSet]:
            if NodeList[Node]!=0:
                for j in range(3): self.Coord[NodeList[Node]]=np.array((self.Coord[Node][0],self.Coord[Node][1],self.Coord[Node][2]))
                for LoadName in self.NodeValue: self.NodeValue[LoadName][NodeList[Node]]=self.NodeValue[LoadName][Node]
#===================================================================
#
#         Write coordinates of a node set
#
# Variables:
# FileName - Name of file
# NSet     - Name of sets of divided nodes
#===================================================================
    def NodesCoord(self,FileName,NSet):
        f=open(FileName,'w')
        for Node in self.NSets[NSet]:
            f.write(str(self.Coord[Node][0])+', '+str(self.Coord[Node][1])+', '+str(self.Coord[Node][2])+'\n')
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
    def ProjectNodesToSurf(self,PrjctNodes,SurfElems,SurfNodes):
            FaceList=[]
            Mtrx=np.zeros((3,3))
            Pnt0=np.zeros(3)
            for i in self.ESets[SurfElems]:
                for Indx in FacesNodes[self.Eltype[i]]:
                    Flag=True
                    for j in Indx:
                        if not self.Elems[i][j] in self.NSets[SurfNodes]: Flag=False
                    if Flag:
                        Pnt0[:]=self.Coord[self.Elems[i][Indx[0]]][:]
                        Mtrx[0]=np.cross(self.Coord[self.Elems[i][Indx[1]]][:]-Pnt0[:],self.Coord[self.Elems[i][Indx[2]]][:]-Pnt0[:])
                        Mtrx[1][:]=self.Coord[self.Elems[i][Indx[1]]][:]-Pnt0[:]
                        Mtrx[2][:]=self.Coord[self.Elems[i][Indx[2]]][:]-Pnt0[:]
                        FaceList.append(((i,Indx),Mtrx,np.dot(Mtrx[0],Pnt0)))
            #-----------------------------
            for Node_pr in self.NSets[PrjctNodes]:
                minDist2=0
                for ElemFace in FaceList:
                    minFaceDist2=np.linalg.norm(self.Coord[self.Elems[ElemFace[0][0]][ElemFace[0][1][0]]][:]-self.Coord[Node_pr][:])
                    for k in range(1, len(ElemFace[0][1])):
                        Dist2=np.linalg.norm(self.Coord[self.Elems[ElemFace[0][0]][ElemFace[0][1][k]]][:]-self.Coord[Node_pr][:])
                        if Dist2<minFaceDist2:minFaceDist2=Dist2
                    if (minDist2==0)or minDist2>minFaceDist2:
                        minDist2=minFaceDist2
                        minFace=ElemFace
                PrjctCrd=np.linalg.solve(minFace[1],np.array((minFace[2],np.dot(minFace[1][1],self.Coord[Node_pr]),np.dot(minFace[1][2],self.Coord[Node_pr]))))
                self.Coord[Node_pr]=tuple(PrjctCrd)
#===================================================================
#
#         Delete elements
#
# Variables:
# EsetNames - List of Names for element sets
#===================================================================
    def DeleteElements(self,EsetNames):
        NodeFlag=np.full(self.MaxNodeNum+1,False)
        for ESetName in EsetNames:
            for El in self.ESets[ESetName]:
                for Node in self.Elems[El]:
                    NodeFlag[Node]=True
                self.Elems[El]=1
            self.ESets.pop(ESetName)
        for El in range(1,self.MaxElemNum+1):
            if self.Elems[El]!=1:
                for Node in self.Elems[El]:
                    if NodeFlag[Node]: NodeFlag[Node]=False 
        for i in range(1,self.MaxNodeNum+1):
            if NodeFlag[i]: self.Coord[i]=None
        for NSetName in list(self.NSets.keys()):
            NumList=self.NSets[NSetName].copy()
            self.NSets[NSetName]=[]
            for Node in NumList:
                if not NodeFlag[Node]:
                    self.NSets[NSetName].append(Node)
            if len(self.NSets[NSetName])==0:self.NSets.pop(NSetName)
        for ESetName in list(self.ESets.keys()):
            NumList=self.ESets[ESetName].copy()
            self.ESets[ESetName]=[]
            for El in NumList:
                if self.Elems[El]!=1: self.ESets[ESetName].append(El)
            if len(self.ESets[ESetName])==0:self.ESets.pop(ESetName)
        print('Elements have been deleted')
#===================================================================
#
#         Create submodel
#
# Variables:
# CentralNodes - List of Nodes nummbers
# Radius       - Radius around nodes to catch elements
#===================================================================
    def CreateSubModel(self,CentralNodes,Radius):
        if len(self.Faces)==0: self.EstFaces()
        self.NSets['NAll']=[]
        self.ESets['EAll']=[]
        self.NSets['SubmodelNodes']=[]
        NodeFlag=np.full(self.MaxNodeNum+1,False)
        ElFlag=np.full(self.MaxElemNum+1,False)
        for Node in range(1,self.MaxNodeNum+1):
            if type(self.Coord[Node])==np.ndarray:
                for CNode in CentralNodes:
                    if np.linalg.norm(np.array(self.Coord[Node])[:]-np.array(self.Coord[CNode])[:])<=Radius:
                        NodeFlag[Node]=True                
        for El in range(1,self.MaxElemNum+1):
            if self.Elems[El]!=1:
                for Node in self.Elems[El]:
                    if NodeFlag[Node]:
                        ElFlag[El]=True
                        self.ESets['EAll'].append(El)
                        break
        for El in self.ESets['EAll']:
            for Node in self.Elems[El]:
                if not NodeFlag[Node]:NodeFlag[Node]=True            
        #cleaning
        for Node in range(1,self.MaxNodeNum+1):
            if not NodeFlag[Node]: self.Coord[Node]=None
        for El in range(1,self.MaxElemNum+1):
            if not ElFlag[El]:
                self.Elems[El]=1
                self.Eltype[El]=0
        for NSet in self.NSets:
            Dict=[]
            for Node in self.NSets[NSet]:
                if type(self.Coord[Node])==np.ndarray: Dict.append(Node)
            self.NSets[NSet]=Dict
        for ESet in self.ESets.keys():
            Dict=[]
            for El in self.ESets[ESet]:
                if self.Elems[El]!=1: Dict.append(El)
            self.ESets[ESet]=Dict
        #------------------
        Faces0=self.Faces.copy()
        self.EstFaces()
        for El in self.ESets['EAll']:
            for Indx in FacesNodes[self.Eltype[El]]:
                Nodes=set()
                for i in Indx: Nodes.add(self.Elems[El][i])
                minNode=min(Nodes)
                maxNode=max(Nodes)
                if minNode in self.Faces:
                    if maxNode in self.Faces[minNode]:
                        for Fc in self.Faces[minNode][maxNode]:
                            if Fc[0]==1:
                                for Fc0 in Faces0[minNode][maxNode]:
                                    if Fc[1]==Fc0[1] and Fc0[0]==2:
                                        self.NSets['SubmodelNodes']+=list(Fc[1])
        print('Submodel has been prepared')
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
    def ShiftSurf(self,D0,D1,Axis,ShiftDir,NodeSet):
        for i in range(len(self.NSets[NodeSet])):
            Coords=self.Coord[self.NSets[NodeSet][i]]
            AxisCoord=np.dot(Axis,Coords)
            if i==0:
                CoordMin=AxisCoord
                CoordMax=AxisCoord
            elif AxisCoord<CoordMin:
                CoordMin=AxisCoord
            elif AxisCoord>CoordMax:
                CoordMax=AxisCoord
        setNode=set(self.NSets[NodeSet])
        MiddleNodes={}
        for Node in self.NSets[NodeSet]: MiddleNodes[Node]=0
        for El in range(1,self.MaxElemNum+1):
            if self.Elems[El]!=1:
                for i in range(len(FacesNodes[self.Eltype[El]])):
                    Flag=True
                    for Node_i in FacesNodes[self.Eltype[El]][i]:
                        if not self.Elems[El][Node_i] in setNode: Flag=False
                    if Flag:
                        for j in range(len(FacesNodes[self.Eltype[El]][i])):
                            if MiddleNodes[self.Elems[El][FacesNodes[self.Eltype[El]][i][j]]]==0:                    
                                if self.Eltype[El]==6:
                                    if j==0:
                                        if i==0 and not self.Elems[El][7] in setNode: MiddleNodes[self.Elems[El][0]]=self.Elems[El][7]
                                        if i==1 and not self.Elems[El][6] in setNode: MiddleNodes[self.Elems[El][0]]=self.Elems[El][6]
                                        if i==3 and not self.Elems[El][4] in setNode: MiddleNodes[self.Elems[El][0]]=self.Elems[El][4]
                                    if j==1:
                                        if i==0 and not self.Elems[El][8] in setNode: MiddleNodes[self.Elems[El][1]]=self.Elems[El][8]
                                        if i==1 and not self.Elems[El][5] in setNode: MiddleNodes[self.Elems[El][1]]=self.Elems[El][5]
                                        if i==2 and not self.Elems[El][4] in setNode: MiddleNodes[self.Elems[El][1]]=self.Elems[El][4]
                                    if j==2:
                                        if i==0 and not self.Elems[El][9] in setNode: MiddleNodes[self.Elems[El][2]]=self.Elems[El][9]
                                        if i==2 and not self.Elems[El][6] in setNode: MiddleNodes[self.Elems[El][2]]=self.Elems[El][6]
                                        if i==3 and not self.Elems[El][5] in setNode: MiddleNodes[self.Elems[El][2]]=self.Elems[El][5]
                                    if j==3:
                                        if i==1 and not self.Elems[El][9] in setNode: MiddleNodes[self.Elems[El][3]]=self.Elems[El][9]
                                        if i==2 and not self.Elems[El][7] in setNode: MiddleNodes[self.Elems[El][3]]=self.Elems[El][7]
                                        if i==3 and not self.Elems[El][8] in setNode: MiddleNodes[self.Elems[El][3]]=self.Elems[El][8]
                                if self.Eltype[El]==7:                    
                                    if j==0:
                                        if i==0 and not self.Elems[El][12] in setNode: MiddleNodes[self.Elems[El][0]]=self.Elems[El][12]
                                        if i==2 and not self.Elems[El][8] in setNode: MiddleNodes[self.Elems[El][0]]=self.Elems[El][8]
                                        if i==4 and not self.Elems[El][6] in setNode: MiddleNodes[self.Elems[El][0]]=self.Elems[El][6]
                                    if j==1:
                                        if i==0 and not self.Elems[El][13] in setNode: MiddleNodes[self.Elems[El][1]]=self.Elems[El][13]
                                        if i==2 and not self.Elems[El][7] in setNode: MiddleNodes[self.Elems[El][1]]=self.Elems[El][7]
                                        if i==3 and not self.Elems[El][6] in setNode: MiddleNodes[self.Elems[El][1]]=self.Elems[El][6]
                                    if j==2:
                                        if i==0 and not self.Elems[El][14] in setNode: MiddleNodes[self.Elems[El][2]]=self.Elems[El][14]
                                        if i==3 and not self.Elems[El][8] in setNode: MiddleNodes[self.Elems[El][2]]=self.Elems[El][8]
                                        if i==4 and not self.Elems[El][7] in setNode: MiddleNodes[self.Elems[El][2]]=self.Elems[El][7]
                                    if j==3:
                                        if i==1 and not self.Elems[El][12] in setNode: MiddleNodes[self.Elems[El][3]]=self.Elems[El][12]
                                        if i==2 and not self.Elems[El][11] in setNode: MiddleNodes[self.Elems[El][3]]=self.Elems[El][11]
                                        if i==4 and not self.Elems[El][9] in setNode: MiddleNodes[self.Elems[El][3]]=self.Elems[El][9]
                                    if j==4:
                                        if i==1 and not self.Elems[El][13] in setNode: MiddleNodes[self.Elems[El][4]]=self.Elems[El][13]
                                        if i==2 and not self.Elems[El][10] in setNode: MiddleNodes[self.Elems[El][4]]=self.Elems[El][10]
                                        if i==3 and not self.Elems[El][9] in setNode: MiddleNodes[self.Elems[El][4]]=self.Elems[El][9]
                                    if j==5:
                                        if i==1 and not self.Elems[El][14] in setNode: MiddleNodes[self.Elems[El][5]]=self.Elems[El][14]
                                        if i==3 and not self.Elems[El][11] in setNode: MiddleNodes[self.Elems[El][5]]=self.Elems[El][11]
                                        if i==4 and not self.Elems[El][10] in setNode: MiddleNodes[self.Elems[El][5]]=self.Elems[El][10]
                                if self.Eltype[El]==8 or self.Eltype[El]==9:                    
                                    if j==0:
                                        if i==0 and not self.Elems[El][16] in setNode: MiddleNodes[self.Elems[El][0]]=self.Elems[El][16]
                                        if i==2 and not self.Elems[El][11] in setNode: MiddleNodes[self.Elems[El][0]]=self.Elems[El][11]
                                        if i==5 and not self.Elems[El][8] in setNode: MiddleNodes[self.Elems[El][0]]=self.Elems[El][8]
                                    if j==1:
                                        if i==0 and not self.Elems[El][17] in setNode: MiddleNodes[self.Elems[El][1]]=self.Elems[El][17]
                                        if i==2 and not self.Elems[El][9] in setNode: MiddleNodes[self.Elems[El][1]]=self.Elems[El][9]
                                        if i==3 and not self.Elems[El][8] in setNode: MiddleNodes[self.Elems[El][1]]=self.Elems[El][8]
                                    if j==2:
                                        if i==0 and not self.Elems[El][18] in setNode: MiddleNodes[self.Elems[El][2]]=self.Elems[El][18]
                                        if i==3 and not self.Elems[El][10] in setNode: MiddleNodes[self.Elems[El][2]]=self.Elems[El][10]
                                        if i==4 and not self.Elems[El][9] in setNode: MiddleNodes[self.Elems[El][2]]=self.Elems[El][9]
                                    if j==3:
                                        if i==0 and not self.Elems[El][19] in setNode: MiddleNodes[self.Elems[El][3]]=self.Elems[El][19]
                                        if i==4 and not self.Elems[El][11] in setNode: MiddleNodes[self.Elems[El][3]]=self.Elems[El][11]
                                        if i==5 and not self.Elems[El][10] in setNode: MiddleNodes[self.Elems[El][3]]=self.Elems[El][10]
                                    if j==4:
                                        if i==1 and not self.Elems[El][16] in setNode: MiddleNodes[self.Elems[El][4]]=self.Elems[El][16]
                                        if i==2 and not self.Elems[El][15] in setNode: MiddleNodes[self.Elems[El][4]]=self.Elems[El][15]
                                        if i==5 and not self.Elems[El][12] in setNode: MiddleNodes[self.Elems[El][4]]=self.Elems[El][12]
                                    if j==5:
                                        if i==1 and not self.Elems[El][17] in setNode: MiddleNodes[self.Elems[El][5]]=self.Elems[El][17]
                                        if i==2 and not self.Elems[El][13] in setNode: MiddleNodes[self.Elems[El][5]]=self.Elems[El][13]
                                        if i==3 and not self.Elems[El][12] in setNode: MiddleNodes[self.Elems[El][5]]=self.Elems[El][12]
                                    if j==6:
                                        if i==1 and not self.Elems[El][18] in setNode: MiddleNodes[self.Elems[El][6]]=self.Elems[El][18]
                                        if i==3 and not self.Elems[El][14] in setNode: MiddleNodes[self.Elems[El][6]]=self.Elems[El][14]
                                        if i==4 and not self.Elems[El][13] in setNode: MiddleNodes[self.Elems[El][6]]=self.Elems[El][13]
                                    if j==7:
                                        if i==1 and not self.Elems[El][19] in setNode: MiddleNodes[self.Elems[El][7]]=self.Elems[El][19]
                                        if i==4 and not self.Elems[El][15] in setNode: MiddleNodes[self.Elems[El][7]]=self.Elems[El][15]
                                        if i==5 and not self.Elems[El][14] in setNode: MiddleNodes[self.Elems[El][7]]=self.Elems[El][14]
        for Node in self.NSets[NodeSet]:
            Coords=np.array(self.Coord[Node])
            AxisCoord=np.dot(Axis,Coords)
            Vect=ShiftDir*(D0+(D1-D0)*(AxisCoord-CoordMin)/(CoordMax-CoordMin))
            self.Coord[Node][:]+=Vect[:]
            if MiddleNodes[Node]!=0:self.Coord[MiddleNodes[Node]][:]+=0.5*Vect[:]
#===================================================================
#
#         Merge nodes
#
# MasterNds - Nodes with constant coordinates
# SlaveNds  - Movable nodes
#===================================================================
    def MergeNodes(self,MasterNds,SlaveNds):
        Cnct=np.zeros(self.MaxNodeNum+1,dtype=np.int32)
        for SlaveNode in self.NSets[SlaveNds]:
            Cnct[SlaveNode]=self.NSets[MasterNds][0]
            R2min=(self.Coord[SlaveNode][0]-self.Coord[Cnct[SlaveNode]][0])**2+\
            (self.Coord[SlaveNode][1]-self.Coord[Cnct[SlaveNode]][1])**2+\
            (self.Coord[SlaveNode][2]-self.Coord[Cnct[SlaveNode]][2])**2
            for MasterNode in self.NSets[MasterNds]:
                R2=(self.Coord[MasterNode][0]-self.Coord[SlaveNode][0])**2+\
                (self.Coord[MasterNode][1]-self.Coord[SlaveNode][1])**2+\
                (self.Coord[MasterNode][2]-self.Coord[SlaveNode][2])**2
                if R2min>R2:
                    R2min=R2
                    Cnct[SlaveNode]=MasterNode
        for SlaveNode in self.NSets[SlaveNds]:
            self.Coord[SlaveNode]=None
        for NSet in self.NSets:
            for i in range(len(self.NSets[NSet])):
                Node=self.NSets[NSet][i]
                if Cnct[Node]>0: self.NSets[NSet][i]=Cnct[Node]
        for El in range(1,self.MaxElemNum+1):
            if self.Elems[El]!=1:
                for i in range(len(self.Elems[El])):
                    Node=self.Elems[El][i]
                    if Cnct[Node]>0: self.Elems[El][i]=Cnct[Node]
#===================================================================
#
#         Renumbering
#
#===================================================================
    def ReNumb(self):
        mesh=FEMtoolkit()
        #-----node
        NodeRef=np.zeros(self.MaxNodeNum+1,dtype=np.int32)
        Count=0
        for i in range(1,self.MaxNodeNum+1):
            if type(self.Coord[i])==np.ndarray:
                Count+=1
                NodeRef[i]=Count
        mesh.MaxNodeNum=Count
        mesh.Coord=np.full(mesh.MaxNodeNum+1,None)
        for i in range(1,self.MaxNodeNum+1):
            if type(self.Coord[i])==np.ndarray: mesh.Coord[NodeRef[i]]=self.Coord[i].copy()
        for SetName in self.NSets.keys():
            mesh.NSets[SetName]=[]
            for Node in self.NSets[SetName]: mesh.NSets[SetName].append(NodeRef[Node])
        #-----elements
        ElRef=np.zeros(self.MaxElemNum+1,dtype=np.int32)
        Count=0
        for i in range(self.MaxElemNum+1):
            if self.Elems[i]!=1:
                Count+=1
                ElRef[i]=Count
        mesh.MaxElemNum=Count
        mesh.Elems=np.ones(mesh.MaxElemNum+1,dtype=tuple)
        mesh.Eltype=np.zeros(mesh.MaxElemNum+1,dtype=np.int8)
        for i in range(self.MaxElemNum+1):
            if self.Elems[i]!=1:
                Nodes=[]
                for Node in self.Elems[i]: Nodes.append(NodeRef[Node])
                mesh.Elems[ElRef[i]]=list(Nodes)
                mesh.Eltype[ElRef[i]]=self.Eltype[i].copy()
        for AbaqElemType in self.TypeList.keys():
            mesh.TypeList[AbaqElemType]=self.TypeList[AbaqElemType]
        for SetName in self.ESets.keys():
            mesh.ESets[SetName]=[]
            for i in self.ESets[SetName]:mesh.ESets[SetName].append(ElRef[i])
        for SetName in self.Surfs.keys():
            mesh.Surfs[SetName]=self.Surfs[SetName].copy()
        return mesh
#===================================================================
#            READERS:
#===================================================================
#         import Abaqus inp-file
#===================================================================
def import_abq(FileName):
    mesh=FEMtoolkit()
    mesh.MaxNodeNum=0
    mesh.MaxElemNum=0
    Section=''
    NSet=''
    ESet=''
    ElemNodeNum=0
    ElementType=''
    f=open(FileName,'r')
    txt=f.readline()[:-1]
    while txt:
        if Section=='node':
            while txt and not '*' in txt:
                ValueTxt=txt.split(',')
                NodeNum=int(ValueTxt[0])
                if NodeNum>mesh.MaxNodeNum:mesh.MaxNodeNum=NodeNum
                txt=f.readline()[:-1]
                while '**' in txt: txt=f.readline()[:-1]
            Section=''
        if Section=='element':
            while txt and not '*' in txt:
                ValueTxt=txt.split(',')
                ElemNum=int(ValueTxt[0])
                if ElemNum>mesh.MaxElemNum:mesh.MaxElemNum=ElemNum
                Num=len(ValueTxt)
                if ValueTxt[Num-1]=='':Num-=1                    
                while Num<ElemNodeNum:
                    txt=f.readline()[:-1]
                    ValueTxt=txt.split(',')
                    for Val in ValueTxt:
                        if Val!='':
                            Num+=1                    
                txt=f.readline()[:-1]
                while '**' in txt: txt=f.readline()[:-1]
            Section=''
        if '*node' in txt.lower() and not '*node output' in txt.lower(): Section='node'
        if '*element' in txt.lower() and not '*element output' in txt.lower():
            Section='element'
            SetNamePos=txt.lower().find('type')+5
            if ',' in txt[SetNamePos:]: ElementType=txt[SetNamePos:txt.find(',',SetNamePos)].replace(' ','')
            else: ElementType=txt[SetNamePos:].replace(' ','')
            mesh.TypeList[AbaqElemTypes[ElementType]]=ElementType
            if FacesNodes[AbaqElemTypes[ElementType]]!=None:
                ElemNodeNum=0
                for Nodes in FacesNodes[AbaqElemTypes[ElementType]]:
                    maxNode=max(Nodes)
                    if ElemNodeNum<maxNode:ElemNodeNum=maxNode
                ElemNodeNum+=1
            else:ElemNodeNum=1
        txt=f.readline()[:-1]
    f.close()
    mesh.Coord=np.full(mesh.MaxNodeNum+1,None)
    mesh.Elems=np.ones(mesh.MaxElemNum+1,dtype=tuple)
    mesh.Eltype=np.zeros(mesh.MaxElemNum+1,dtype=np.int8)
    f=open(FileName,'r')
    Section=''
    NSet=''
    ESet=''
    Surf=''        
    txt=f.readline()[:-1]
    while txt:
        if Section=='node':
            while txt and not '*' in txt:
                ValueTxt=txt.split(',')
                NodeNum=int(ValueTxt[0])
                if NSet!='': mesh.NSets[NSet].append(NodeNum)
                mesh.Coord[NodeNum]=np.array(list(map(float,ValueTxt[1:])))
                txt=f.readline()[:-1]
                while '**' in txt: txt=f.readline()[:-1]
            Section=''
            NSet=''
        elif NSet!='':
            while txt and not '*' in txt:
                for Val in txt.replace(' ','').split(','):
                    if Val in mesh.NSets:
                        for NodeNum in mesh.NSets[Val]: mesh.NSets[NSet].append(NodeNum)
                    elif Val!='': mesh.NSets[NSet].append(int(Val))
                txt=f.readline()[:-1]
                while '**' in txt: txt=f.readline()[:-1]
            NSet=''           
        if Section=='element':
            while txt and not '*' in txt:
                ValueTxt=txt.split(',')
                ElemNum=int(ValueTxt[0])
                mesh.Eltype[ElemNum]=AbaqElemTypes[ElementType]
                if ESet!='': mesh.ESets[ESet].append(ElemNum)
                Num=len(ValueTxt)
                if ValueTxt[Num-1]=='':Num-=1
                mesh.Elems[ElemNum]=list(map(int,ValueTxt[1:Num]))
                while Num<ElemNodeNum:
                    txt=f.readline()[:-1]
                    ValueTxt=txt.split(',')
                    for Val in ValueTxt:
                        if Val!='':
                            mesh.Elems[ElemNum].append(int(Val))
                            Num+=1
                txt=f.readline()[:-1]
                while '**' in txt: txt=f.readline()[:-1]
            Section=''
            ESet=''
            ElementType=''
        elif ESet!='':
            while txt and not '*' in txt:
                for Val in txt.replace(' ','').split(','):
                    if Val in mesh.ESets:
                        for ElemNum in mesh.ESets[Val]: mesh.ESets[ESet].append(ElemNum)
                    elif Val!='': mesh.ESets[ESet].append(int(Val))                   
                txt=f.readline()[:-1]
                while '**' in txt: txt=f.readline()[:-1]
            ESet=''
        if Surf!='':
            while txt and not '*' in txt:
                ValueTxt=txt.replace(' ','').split(',')
                if ValueTxt[1][0]=='S':ValueTxt[1]=ValueTxt[1][1:]
                mesh.Surfs[Surf].append((ValueTxt[0],int(float(ValueTxt[1]))-1))               
                txt=f.readline()[:-1]
                while '**' in txt: txt=f.readline()[:-1]
            Surf=''            
        if '*node' in txt.lower() and not '*node output' in txt.lower():
            Section='node'
            txt.replace(' ','')
            SetNamePos=txt.lower().find('nset')
            if SetNamePos>4:
                SetNamePos=txt.find('=',SetNamePos)
                if ',' in txt[SetNamePos:]: NSet=txt[SetNamePos+1:txt.find(',',SetNamePos)]
                else: NSet=txt[SetNamePos+1:]                    
                if not NSet in mesh.NSets: mesh.NSets[NSet]=[]
        if '*element' in txt.lower() and not '*element output' in txt.lower():
            Section='element'
            txt.replace(' ','')
            SetNamePos=txt.lower().find('elset')
            if SetNamePos>4:
                SetNamePos=txt.find('=',SetNamePos)
                if ',' in txt[SetNamePos:]: ESet=txt[SetNamePos+1:txt.find(',',SetNamePos)]
                else: ESet=txt[SetNamePos+1:]
                if not ESet in mesh.ESets: mesh.ESets[ESet]=[]
            SetNamePos=txt.lower().find('type')+5
            if ',' in txt[SetNamePos:]: ElementType=txt[SetNamePos:txt.find(',',SetNamePos)].replace(' ','')
            else: ElementType=txt[SetNamePos:].replace(' ','')
            if FacesNodes[AbaqElemTypes[ElementType]]!=None:
                ElemNodeNum=0
                for Nodes in FacesNodes[AbaqElemTypes[ElementType]]:
                    maxNode=max(Nodes)
                    if ElemNodeNum<maxNode:ElemNodeNum=maxNode
                ElemNodeNum+=1
            else:ElemNodeNum=1
        if '*nset' in txt.lower():
            txt.replace(' ','')
            SetNamePos=txt.lower().find('nset',3)
            SetNamePos=txt.find('=',SetNamePos)
            if ',' in txt[SetNamePos:]: NSet=txt[SetNamePos+1:txt.find(',',SetNamePos)]
            else: NSet=txt[SetNamePos+1:]
            if not NSet in mesh.NSets: mesh.NSets[NSet]=[]
            if 'generate' in txt.lower():
                txt=f.readline()[:-1]
                ValueTxt=txt.split(',')
                for NodeNum in range(int(ValueTxt[0]),int(ValueTxt[1])+int(ValueTxt[2]),int(ValueTxt[2])):
                    mesh.NSets[NSet].append(NodeNum)
                NSet=''
        if '*elset' in txt.lower():
            txt.replace(' ','')
            SetNamePos=txt.lower().find('elset',4)
            SetNamePos=txt.find('=',SetNamePos)
            if ',' in txt[SetNamePos:]: ESet=txt[SetNamePos+1:txt.find(',',SetNamePos)]
            else: ESet=txt[SetNamePos+1:]
            if not ESet in mesh.ESets: mesh.ESets[ESet]=[]
            if 'generate' in txt.lower():
                txt=f.readline()[:-1]
                ValueTxt=txt.split(',')
                for ElemNum in range(int(ValueTxt[0]),int(ValueTxt[1])+int(ValueTxt[2]),int(ValueTxt[2])):
                    mesh.ESets[ESet].append(ElemNum)
                ESet=''
        if '*surface' in txt.lower() and not '*surface interaction' in txt.lower():
            txt.replace(' ','')
            SetNamePos=txt.lower().find('name',7)
            if ',' in txt[SetNamePos:]: Surf=txt[SetNamePos+5:txt.find(',',SetNamePos)]
            else: Surf=txt[SetNamePos+5:]
            if not Surf in mesh.Surfs: mesh.Surfs[Surf]=[]
        txt=f.readline()[:-1]
    f.close()
    print('Model has been imported')
    return mesh
#===================================================================
#         import CalculiX frd-file (ASCII)
#===================================================================
def import_CalcX(FileName):
    mesh=FEMtoolkit()
    mesh.MaxNodeNum=0
    mesh.MaxElemNum=0
    f=open(FileName,'r')    
    mesh.TypeList[2]='CAX4'
    mesh.ESets['EAll']=[]
    txt=f.readline()
    while txt:
        if txt[4:6]=='2C':
            RecNum=int(txt[6:36])
            for i in range(RecNum):
                txt=f.readline()
                NodeNum=int(txt[3:13])
                if NodeNum>mesh.MaxNodeNum:mesh.MaxNodeNum=NodeNum
        elif txt[4:6]=='3C':
            RecNum=int(txt[6:36])
            for i in range(RecNum):
                txt=f.readline()
                ElemNum=int(txt[3:13])
                mesh.ESets['EAll'].append(ElemNum)
                if ElemNum>mesh.MaxElemNum:mesh.MaxElemNum=ElemNum
                txt=f.readline()
        txt=f.readline()
    f.close()
    mesh.Coord=np.full(mesh.MaxNodeNum+1,None)
    mesh.Elems=np.ones(mesh.MaxElemNum+1,dtype=tuple)
    mesh.Eltype=np.zeros(mesh.MaxElemNum+1,dtype=np.int8)
    f=open(FileName,'r')
    Section=''
    txt=f.readline()
    while txt:
        if txt[4:6]=='2C':            
            RecNum=int(txt[24:36])
            for i in range(RecNum):
                txt=f.readline()
                NodeNum=int(txt[3:13])
                mesh.Coord[NodeNum]=np.array((float(txt[13:25]),float(txt[25:37]),float(txt[37:])))
        if txt[4:6]=='3C':
            RecNum=int(txt[24:36])
            for i in range(RecNum):
                txt=f.readline()
                ElemNum=int(txt[3:13])
                mesh.Eltype[ElemNum]=CalculiXElemTypes[txt[13:18].replace(' ','')]
                txt=f.readline()
                NodeList=[]
                j0=3
                while j0<len(txt)-1:
                    NodeList.append(int(txt[j0:j0+10]))
                    j0+=10
                mesh.Elems[ElemNum]=NodeList
        if txt[2:7]=='100CL':
            RecNum=int(txt[24:36])
            txt=f.readline()
            ColNum=int(txt[11:18])
            ResNames=[]
            for i in range(ColNum):
                txt=f.readline()
                ResNames.append(txt[3:11].replace(' ',''))
                mesh.NodeValue[ResNames[i]]={}
            for i in range(RecNum):
                txt=f.readline()
                NodeNum=int(txt[3:13])
                j0=0                
                while j0*12+13<len(txt)-1:
                    mesh.NodeValue[ResNames[j0]][NodeNum]=float(txt[j0*12+13:j0*12+25])
                    j0+=1 
        txt=f.readline()
    f.close()
    print('Model has been imported')
    return mesh
