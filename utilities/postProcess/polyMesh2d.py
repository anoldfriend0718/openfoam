import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import linecache
import os
import vtk
from vtk.util import numpy_support as VN

# read internal data file generated by foamToVTK utility
def readInternal(caseDir, m2mm=1, meshGenerator='blockMesh'):
    caseDir=os.path.abspath(caseDir)
    vtkfile = '%s/VTK/%s_0.vtk'%(caseDir,caseDir.split('/')[-1])
    if(not os.path.exists(vtkfile)):
        cmd='foamToVTK -useTimeName -case %s'%(caseDir)
        os.system(cmd)
        if(not os.path.exists(vtkfile)):
            print("The vtk file doesn't exist, please run `foamToVTK` for the case.\n%s"%(vtkfile))
            print('Tried to run command for you: %s'%(cmd))
            print("but the problem is still not solved, do it by hand please")
            return
    data=meshio.read(vtkfile)
    cells=data.cells[0].data
    points=data.points
    cells_poly=[]
    if(meshGenerator=='blockMesh'):
        cells_poly=np.zeros((cells.shape[0],4),dtype=int)
        for i in range(0,cells.shape[0]):
            points_cell=points[cells[i]]
            ind_front=(points_cell[:,2]==points_cell[:,2].max())
            cells_poly[i,:]=cells[i,ind_front]
        cells_poly=cells_poly[:,[0,1,3,2]] # correct node connection order to construct rectangule
        cells_poly=cells_poly
    elif(meshGenerator=='gmsh'):
        cells_tri=np.zeros((cells.shape[0],3),dtype=int)
        for i in range(0,cells.shape[0]):
            points_cell=points[cells[i]]
            ind_front=(points_cell[:,2]==points_cell[:,2].max())
            cells_tri[i,:]=cells[i,ind_front]
        cells_poly=cells_tri
    x,y,z=points[:,0],points[:,1],points[:,2]
    print("nPoints: %d, nPoints_2D: %d, nCells: %d"%(len(points), len(points)/2, len(cells)))
    return x*m2mm, y*m2mm,z*m2mm, cells_poly
# return python dict, {'x':[],'y':[],'z':[]}
def readPoints(caseDir,timeName=None):
    caseDir=os.path.abspath(caseDir)
    file_points='%s/constant/polyMesh/points'%(caseDir)
    if(timeName==None):
        file_points='%s/constant/polyMesh/points'%(caseDir)
    else:
        file_points='%s/%s/polyMesh/points'%(caseDir,timeName)
    if(not os.path.exists(file_points)):
        print("File does not exist: %s"%(file_points))
        return
    alldata=linecache.getlines(file_points)
    linecache.clearcache()
    nPoints, start, end=0,0,0
    points={'x':[],'y':[],'z':[]}
    for i in range(0,len(alldata)):
        alldata[i]=alldata[i].replace('\n','')
        if(alldata[i]=='('):
            nPoints=int(alldata[i-1])
            start=i+1
        elif(alldata[i]==')'):
            end=i
    # read coordinates
    for i in range(start, end):
        xyz=alldata[i].replace('(','').replace(')','').split()
        for key, ind in zip(points.keys(),[0,1,2]):
            points[key].append(xyz[ind])
    for key in points.keys():
        points[key]=np.array(points[key],dtype=float)
    return points
# read OpenFOAM poly Mesh
# return a python dict {'nNodes':[], 'index':[]}
def readFaces(caseDir):
    caseDir=os.path.abspath(caseDir)
    file_faces='%s/constant/polyMesh/faces'%(caseDir)
    if(not os.path.exists(file_faces)):
        print("File does not exist: %s"%(file_faces))
        return
    alldata=linecache.getlines(file_faces)
    linecache.clearcache()
    nFaces,face_start,face_end=0,0,0
    for i in range(0,len(alldata)):
        alldata[i]=alldata[i].replace('\n','')
        if(alldata[i]=='('):
            nFaces=int(alldata[i-1])
            face_start=i+1
        elif(alldata[i]==')'):
            face_end=i-1
    if(not nFaces==(face_end-face_start+1)):
        print('The faces file of case %s maybe not correct, because the face number is not consistant'%(caseDir))
    # read point index of each face
    faces={'nNodes':[],'index':[]}
    for i in range(face_start, face_end+1):
        n=int(alldata[i].split('(')[0])
        index=alldata[i].split('(')[1].split(')')[0].split(' ')
        faces['nNodes'].append(n)
        faces['index'].append(np.array(index,dtype=int))
    return faces
# read cellZones
def readCellZones(caseDir):
    caseDir=os.path.abspath(caseDir)
    file_faces='%s/constant/polyMesh/cellZones'%(caseDir)
    if(not os.path.exists(file_faces)):
        # print("File does not exist: %s"%(file_faces))
        print('There is no cellZones in the case of %s'%(caseDir.split('/')[-1]))
        return []
    alldata=linecache.getlines(file_faces)
    linecache.clearcache()
    nCellZones,startCellZones=0,0
    for i in range(0,len(alldata)):
        if(alldata[i][0]=='('):
            alldata[i]=alldata[i].replace('\n','')
            nCellZones=int(alldata[i-1])
            startCellZones=i+1
            break
    names,starts,ends=[],[],[]
    for i in range(startCellZones,len(alldata)):
        if(alldata[i][0]=='{'):
            name_cellZone=alldata[i-1].replace('\n','')
            names.append(name_cellZone)
            starts.append(i+1)
        elif(alldata[i][0]=='}'):
            ends.append(i-1)
    if(not ((len(names)==len(starts)) & (len(names)==len(ends)))):
        print('The cellZones of %s case is not correct, because it has %d cell zone names, %d starts({), %d ends(}), they not equal'%(caseDir,len(names),len(starts),len(ends)))
        exit(0)
    cellZones={}
    for name,start,end in zip(names,starts,ends):
        nCells,cell_start,cell_end=0,0,0
        for i in range(start,end):
            alldata[i]=alldata[i].replace('\n','')
            if(alldata[i][0]=='('):
                nCells=int(alldata[i-1])
                cell_start=i+1
            elif(alldata[i][0]==')'):
                cell_end=i-1
        # print(nCells,cell_start,cell_end)
        cellZones[name]=np.array(alldata[cell_start:cell_end+1], dtype=int)
    return cellZones
# read owner 
# return a int array, the length is equal to number of all faces
def readOwner(caseDir):
    caseDir=os.path.abspath(caseDir)
    file_faces='%s/constant/polyMesh/owner'%(caseDir)
    if(not os.path.exists(file_faces)):
        print("File does not exist: %s"%(file_faces))
        return
    alldata=linecache.getlines(file_faces)
    linecache.clearcache()
    nFaces,face_start,face_end=0,0,0
    for i in range(0,len(alldata)):
        alldata[i]=alldata[i].replace('\n','')
        if(alldata[i]=='('):
            nFaces=int(alldata[i-1])
            face_start=i+1
        elif(alldata[i]==')'):
            face_end=i-1
    if(not nFaces==(face_end-face_start+1)):
        print('The faces file of case %s maybe not correct, because the face number is not consistant'%(caseDir))
    # read point index of each face
    owners=np.array(alldata[face_start:face_end+1], dtype=int)
    return owners
# read neighbour 
# return a int array, the length is equal to number of internal faces
def readNeighbour(caseDir):
    caseDir=os.path.abspath(caseDir)
    file_faces='%s/constant/polyMesh/neighbour'%(caseDir)
    if(not os.path.exists(file_faces)):
        print("File does not exist: %s"%(file_faces))
        return
    alldata=linecache.getlines(file_faces)
    linecache.clearcache()
    nFaces,face_start,face_end=0,0,0
    for i in range(0,len(alldata)):
        alldata[i]=alldata[i].replace('\n','')
        if(alldata[i]=='('):
            nFaces=int(alldata[i-1])
            face_start=i+1
        elif(alldata[i]==')'):
            face_end=i-1
    if(not nFaces==(face_end-face_start+1)):
        print('The faces file of case %s maybe not correct, because the face number is not consistant'%(caseDir))
    # read point index of each face
    neighbours=np.array(alldata[face_start:face_end+1], dtype=int)
    return neighbours
# read boundary patches information 
# return a python dict {'name':[], 'nFaces':[], 'startFace':[], 'type':[], 'index':[]}
def readBoundary_(caseDir, nAllFaces=None):
    caseDir=os.path.abspath(caseDir)
    file_boundary='%s/constant/polyMesh/boundary'%(caseDir)
    if(not os.path.exists(file_boundary)):
        print("File does not exist: %s"%(file_boundary))
        return
    alldata=linecache.getlines(file_boundary)
    linecache.clearcache()
    nBoundaries, start, end=0,0,0
    boundaries={'name':[], 'nFaces':[], 'startFace':[], 'type':[],'index':[]}
    for i in range(0,len(alldata)):
        alldata[i]=alldata[i].replace('\n','')
        if(alldata[i]=='('):
            nBoundaries=int(alldata[i-1])
            start=i+1
        elif(alldata[i]==')'):
            end=i
    # print(start, end)
    # get start line and end line of each patch
    start_patch,end_patch=[],[]
    for i in range(start,end):
        if(len(alldata[i])>0):
            if(alldata[i][-1]=='{'): # a patch start
                start_patch.append(i)
            if(alldata[i][-1]=='}'):
                end_patch.append(i)
    if((not (nBoundaries==len(start_patch))) & (not (nBoundaries==len(start_patch)))):
        print('boundary file parse failure, because boundary number are not consistant: %f'%(caseDir))
    for start, end in zip(start_patch, end_patch):
        boundaries['name'].append(alldata[start-1].split()[0])
        for i in range(start, end):
            if('type' in alldata[i]):
                boundaries['type'].append(alldata[i].split()[1].split(';')[0])
            if('nFaces' in alldata[i]):
                boundaries['nFaces'].append(int(alldata[i].split()[1].split(';')[0]))
            if('startFace' in alldata[i]):
                boundaries['startFace'].append(int(alldata[i].split()[1].split(';')[0]))
    index_allBoundaries=[]
    for n, start in zip(boundaries['nFaces'],boundaries['startFace']):
        index = np.arange(start, n+start).tolist()
        boundaries['index'].append(index)
        index_allBoundaries = index_allBoundaries +index
    # print('nBoundaries: %d, '%(nBoundaries), boundaries['name'])
    # calculate all internal faces index
    index_internalFaces=[]
    name_faces=[]
    if(not nAllFaces==None):
        name_faces=['internal']*nAllFaces
        # get all index of all internal faces
        inds_faces=np.array([True]*nAllFaces)
        inds_faces[index_allBoundaries]=False
        inds_internalFaces=(inds_faces==True)
        index_faces=np.arange(0,nAllFaces)
        index_internalFaces=index_faces[inds_internalFaces]
        for name, index in zip(boundaries['name'],boundaries['index']):
            for ind in index:
                name_faces[ind]=name
        print('nInternalFaces: %d'%(len(index_internalFaces)))
    return boundaries,index_internalFaces,name_faces
def readBoundary(caseDir):
    caseDir=os.path.abspath(caseDir)
    file_boundary='%s/constant/polyMesh/boundary'%(caseDir)
    if(not os.path.exists(file_boundary)):
        print("File does not exist: %s"%(file_boundary))
        exit(0)
    alldata=linecache.getlines(file_boundary)
    linecache.clearcache()
    nBoundaries, start, end=0,0,0
    boundaries={}
    for i in range(0,len(alldata)):
        alldata[i]=alldata[i].replace('\n','')
        if(alldata[i]=='('):
            nBoundaries=int(alldata[i-1])
            start=i+1
        elif(alldata[i]==')'):
            end=i
    # print(start, end)
    # get start line and end line of each patch
    start_patch,end_patch=[],[]
    for i in range(start,end):
        if(len(alldata[i])>0):
            if(alldata[i][-1]=='{'): # a patch start
                start_patch.append(i)
            if(alldata[i][-1]=='}'):
                end_patch.append(i)
    if((not (nBoundaries==len(start_patch))) & (not (nBoundaries==len(start_patch)))):
        print('boundary file parse failure, because boundary number are not consistant: %f'%(caseDir))
    for start, end in zip(start_patch, end_patch):
        name_bd=alldata[start-1].split()[0]
        nFaces,startFaces,type_bd=0,0,''
        for i in range(start, end):
            if('type' in alldata[i]):
                type_bd=alldata[i].split()[1].split(';')[0]
            if('nFaces' in alldata[i]):
                nFaces=int(alldata[i].split()[1].split(';')[0])
            if('startFace' in alldata[i]):
                startFaces=int(alldata[i].split()[1].split(';')[0])
        boundaries[name_bd]={'nFaces':nFaces, 'startFace':startFaces, 'type':type_bd}
    return boundaries
# return python 2D list, [faces], faces=[face_1, face_2, ..., face_n]
def getCells_(owners, neighbours):
    nCells=np.max([owners.max(),neighbours.max()]) + 1
    cells=[]
    for i in range(0,nCells):
        cells.append([])
    for i in range(0, len(owners)):
        cells[owners[i]].append(i)
    for i in range(0, len(neighbours)):
        cells[neighbours[i]].append(i)
    return cells
def getCells(caseDir):
    cells=getCells_(readOwner(caseDir), readNeighbour(caseDir))
    return cells
def read(caseDir):
    caseDir=os.path.abspath(caseDir)
    points=readPoints(caseDir)
    x,y,z=points['x'],points['y'],points['z']
    faces = readFaces(caseDir)
    owners = readOwner(caseDir)
    neighbours = readNeighbour(caseDir)
    faceIndex_cells=getCells_(owners, neighbours)
    boundaries, index_internalFaces, name_faces = readBoundary_(caseDir,len(faces['nNodes']))
    faces['name']=name_faces
    # 1. get empty patch name
    name_emptyPatch=[]
    for name_patch, type_patch in zip(boundaries['name'],boundaries['type']):
        if(type_patch=='empty'):
            name_emptyPatch.append(name_patch)
    # 2. neighbour cells of a cell
    cells={'faces':faceIndex_cells,'neighbour':[],'owner':[]}
    for i in range(0,len(faceIndex_cells)):
        cells['owner'].append([])
        cells['neighbour'].append([])
    # 2.1 faces own to cell
    for ind_cell in range(0,len(cells['faces'])):
        faces_own_to_cell = np.where(owners==ind_cell)[0]
        # print('cell %d is the owner of faces: '%(ind_cell),faces_own_to_cell)
        for face in faces_own_to_cell:
            if(faces['name'][face] in name_emptyPatch):
                continue
            if(faces['name'][face]=='internal'):
                cells['neighbour'][ind_cell].append(neighbours[face])
        # 2.2 faces neighbour to cell
        faces_neighbour_to_cell=np.where(neighbours==ind_cell)[0]
        for face in faces_neighbour_to_cell:
            if(faces['name'][face] in name_emptyPatch):
                continue
            cells['owner'][ind_cell].append(owners[face])
    mesh={'points':points, 'cells':cells, 'faces':faces, 'owners': owners, 'neighbours':neighbours}
    return mesh
def readField_(caseDir,timeName,fieldName,regionName='internalField'):
    caseDir=os.path.abspath(caseDir)
    file_points='%s/%s/%s'%(caseDir,timeName,fieldName)
    if(not os.path.exists(file_points)):
        print("File does not exist: %s"%(file_points))
        return
    alldata=linecache.getlines(file_points)
    linecache.clearcache()
    # find start of region, e.g. internalField, or boundary , e.g. seafloor, ...
    start_region=-1
    for i in range(0,len(alldata)):
        if(regionName in alldata[i].replace('\n','')):
            start_region=i 
            break
    if(start_region==-1):
        print('Did not find %s in the %s field file'%(regionName,fieldName))
        exit(0)
    # print('Start row of %s is %d'%(regionName,start_region))
    nPoints, start, end=0,0,0
    field={}
    for i in range(0,9):
        field['%s_%d'%(fieldName,i)]=[]
    for i in range(start_region,len(alldata)):
        alldata[i]=alldata[i].replace('\n','')
        if(alldata[i]=='('):
            if(start!=0):
                continue
            nPoints=int(alldata[i-1])
            start=i+1
        elif(alldata[i]==')'):
            if(end!=0):
                continue
            end=i
    # read coordinates
    #TODO: cannot process uniform vector (potential bug)
    if "internalField   uniform" in  alldata[start_region]:
        # print(alldata[start_region])
        values=alldata[start_region].replace(" ","").replace("internalField","") \
                                .replace("uniform","").replace(";","") \
                                .replace('(','').replace(')','').split()
        # print(values)
        for i in range(0,len(values)):
            # print(values[i])
            field['%s_%.0f'%(fieldName,i)].append(values[i])
    else:
        for i in range(start, end):
            values=alldata[i].replace('(','').replace(')','').split()
            for i in range(0,len(values)):
                field['%s_%.0f'%(fieldName,i)].append(values[i])
    
    keys=list(field.keys())
    for i in range(0,len(keys)):
        key=keys[i]
        if(field[key]==[]):
            del field[key]
            continue
        field[key]=np.array(field[key],dtype=float)
    # if it is a scalar field, use field name as key 
    keys=list(field.keys())
    if(len(keys)==1):
        field[fieldName]=field[keys[0]]
        del field[keys[0]]
    return field

def readField(caseDir,timeName,fieldName,regionName='internalField'):
    field=readField_(caseDir,timeName,fieldName,regionName)
    # convert field to array
    keys=list(field.keys())
    nArray=len(field[keys[0]])
    fieldArray=np.zeros((nArray,len(keys)))
    for i in range(0,len(keys)):
        fieldArray[:,i]=field[keys[i]]
    return fieldArray
def getMesh(caseDir, patchName):
    # 1. read mesh as triangles
    points=readPoints(caseDir)
    x,y,z=points['x'],points['y'],points['z']
    cells=getCells(caseDir)
    faces=readFaces(caseDir)
    boundaries,index_internalFaces,name_faces=readBoundary_(caseDir)
    ind_bd=0
    try:
        ind_bd=(boundaries['name'].index(patchName))
    except:
        print('%sError: %s%s is not a valid patch name: '%(C_RED,C_DEFAULT,patchName),boundaries['name'])
        exit(0)
    startFace=boundaries['startFace'][ind_bd]
    nFaces=boundaries['nFaces'][ind_bd]
    endFaces=startFace+nFaces
    polygons=[]
    validCells_index=[False]*len(cells)
    for i in range(0,len(cells)):
        cell=np.array(cells[i])
        # for face in cell:
        validFaces=np.array(np.where((cell>=startFace)  & (cell<endFaces) ))[0]
        if(len(validFaces)>0):
            polygons.append(faces['index'][cell[validFaces[0]]])
            validCells_index[i]=True
    polygons=np.array(polygons)
    validCells_index=np.array(validCells_index)
    triangles=[]
    polygon2triangle=[]
    for i in range(0,len(polygons)):
        polygon=polygons[i]
        for j in range(0,len(polygon)-2):
            triangles.append([polygon[0],polygon[1+j],polygon[2+j]])
        polygon2triangle.append(len(polygon)-2)
    triangles=np.array(triangles)
    polygon2triangle=np.array(polygon2triangle)
    # 1.3 construct unstructured mesh grid
    vtk_points = vtk.vtkPoints()
    vtk_triangle = vtk.vtkTriangle()
    vtk_cells = vtk.vtkCellArray()
    # vtk points
    x=x[0:int(len(x)/2)] # only retain one side , either front or back
    y=y[0:int(len(y)/2)]
    for i in range(len(x)):
        vtk_points.InsertNextPoint(x[i],y[i],z[0])
    # vtk cells
    for cell in triangles:
        for i in range(0,3):
            vtk_triangle.GetPointIds().SetId(i,cell[i])
        vtk_cells.InsertNextCell(vtk_triangle)
    # construct unstructedgrid
    VTU = vtk.vtkUnstructuredGrid()
    VTU.SetPoints(vtk_points)
    VTU.SetCells(vtk.VTK_TRIANGLE, vtk_cells)
    MeshData={'x':x,'y':y,'z':z,'triangles':triangles,'polygons':polygons,'validCells_index':validCells_index,'poly2tri':polygon2triangle,'FieldData_usg':VTU}
    return MeshData #x,y,triangles,polygons,validCells_index,polygon2triangle,VTU
def readCellData_to_pointData(caseDir, timename, fieldNames,MeshData):
    validCells_index=MeshData['validCells_index']
    polygon2triangle=MeshData['poly2tri']
    VTU=MeshData['FieldData_usg']
    # 1.2 read cell data
    dataNames=[]
    cellData={}
    for fieldName in fieldNames:
        field=readField_(caseDir,timename,fieldName)
        
        for i,dataName in zip(range(0,len(field)),field.keys()):
            dataNames.append(dataName)
            celldata = vtk.vtkFloatArray()
            celldata.SetNumberOfValues(MeshData['triangles'].shape[0])
            celldata.SetName(dataName)
            if field[dataName].shape[0]==1:
                num=validCells_index.shape[0]
                cellData[dataName]=np.repeat(field[dataName],num)
            else:
                cellData[dataName]=field[dataName][validCells_index] # only extract field values associated with polygons
            # print(cellData.keys())
            index_cell=0
            for i in range(0,len(cellData[dataName])):
                for j in range(0,polygon2triangle[i]):
                    celldata.SetValue( index_cell, cellData[dataName][i] )
                    index_cell=index_cell+1
            VTU.GetCellData().AddArray(celldata)
    c2p = vtk.vtkCellDataToPointData()
    c2p.SetInputData(VTU)
    c2p.Update()
    pointdata={}
    for dataName in dataNames:
        pointdata[dataName]=VN.vtk_to_numpy(c2p.GetOutput().GetPointData().GetArray(dataName))
    return {'pointData':pointdata,'cellData':cellData,"x":MeshData["x"],"y":MeshData["y"]}
def getTimes(caseDir):
    times=os.listdir(caseDir)
    timeDirs=[]
    for t in times:
        if(os.path.isdir(caseDir+'/'+t)):
            if(t.replace('.','',1).isdigit()):
                timeDirs.append(t)
    #         else:
    #             print(t,'is not a directory')
    timeDirs=np.array(timeDirs)
    times=np.array(timeDirs,dtype=float)
    ind=np.argsort(times)
    return timeDirs[ind],times[ind]
#unit normal vector of plane defined by points a, b, and c
def unit_normal(a, b, c):
    x = np.linalg.det([[1,a[1],a[2]],
         [1,b[1],b[2]],
         [1,c[1],c[2]]])
    y = np.linalg.det([[a[0],1,a[2]],
         [b[0],1,b[2]],
         [c[0],1,c[2]]])
    z = np.linalg.det([[a[0],a[1],1],
         [b[0],b[1],1],
         [c[0],c[1],1]])
    magnitude = (x**2 + y**2 + z**2)**.5
    return (x/magnitude, y/magnitude, z/magnitude)

#area of polygon poly
def poly_area_normal(poly):
    if len(poly) < 3: # not a plane - no area
        return 0
    total = [0, 0, 0]
    N = len(poly)
    for i in range(N):
        vi1 = poly[i]
        vi2 = poly[(i+1) % N]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    #unit normal
    nf=unit_normal(poly[0], poly[1], poly[2])
    area = np.dot(total, nf)/2
    return area,np.array(nf)
# convert polygons to triangles, and/or convert polygon cell data to triangle cell data and point data
def polyData2TriData(polygons,x=[],y=[],z=[],cellFields=[],pointFields=[]):
    triangles,cellData_tri,pointData_tri=[],{},{}
    # tri mesh
    for polygon in polygons:
        for j in range(0,len(polygon)-2):
            triangles.append([polygon[0],polygon[1+j],polygon[2+j]])
    triangles=np.array(triangles,dtype=int)
    # cell data
    numPolygons=len(polygons)
    if(cellFields!=[]):
        for fieldName in cellFields.keys():
            cellData_tri[fieldName]=[]
            numCellData=len(cellFields[fieldName])
            if(numCellData!=numPolygons):
                print('Error: polygon cell field %s has %d value, but there are all %d polygons'%(fieldName,numCellData,numPolygons))
                exit(0)
            for cellField, polygon in zip(cellFields[fieldName],polygons):
                for j in range(0,len(polygon)-2):
                    cellData_tri[fieldName].append(cellField)
            cellData_tri[fieldName]=np.array(cellData_tri[fieldName])
    # cell data to point data
    if((len(x)>0) & (len(y)>0) & (len(z)>0) & (len(cellData_tri)>0)):
        # construct unstructured mesh grid
        vtk_points = vtk.vtkPoints()
        vtk_triangle = vtk.vtkTriangle()
        vtk_cells = vtk.vtkCellArray()
        # vtk points
        for i in range(len(x)):
            vtk_points.InsertNextPoint(x[i],y[i],z[i])
        # vtk cells
        for cell in triangles:
            for i in range(0,3):
                vtk_triangle.GetPointIds().SetId(i,cell[i])
            vtk_cells.InsertNextCell(vtk_triangle)
        # construct unstructedgrid
        VTU = vtk.vtkUnstructuredGrid()
        VTU.SetPoints(vtk_points)
        VTU.SetCells(vtk.VTK_TRIANGLE, vtk_cells)
        # add cell data
        for dataName in cellData_tri.keys():
            celldata = vtk.vtkFloatArray()
            numComponents=len(cellData_tri[dataName][0])
            celldata.SetNumberOfComponents(numComponents)
            celldata.SetNumberOfTuples(triangles.shape[0])
            celldata.SetName(dataName)
            if(numComponents==1):
                for i in range(0,len(cellData_tri[dataName])):
                    celldata.SetValue( i, cellData_tri[dataName][i] )
                VTU.GetCellData().AddArray(celldata)
            elif(numComponents==3):
                for i in range(0,len(cellData_tri[dataName])):
                    celldata.SetTuple3( i, cellData_tri[dataName][i][0], cellData_tri[dataName][i][1],cellData_tri[dataName][i][2])
                VTU.GetCellData().AddArray(celldata)
            else:
                print('field %s has %d components, the polyData2TriData function can only deal with scalar and vector field sofar'%(dataName,numComponents))
                exit(0)
        # celldata to point data
        c2p = vtk.vtkCellDataToPointData()
        c2p.SetInputData(VTU)
        c2p.Update()
        pointdata=c2p.GetOutput().GetPointData()
        for i in range(0,pointdata.GetNumberOfArrays()):
            arrayname=pointdata.GetArrayName(i)
            pointData_tri[arrayname]=VN.vtk_to_numpy(pointdata.GetArray(i))
    return {'triangles':triangles,'pointData':pointData_tri,'cellData':cellData_tri}
# read cell zone index and fields on a patch
def getDataOnPatch_cellZones(caseDir,timeName,patchName,cellZoneNames,fieldNames,calSf=False,triMesh=True):
    # variables need to return
    x,y,z,Sf,fields,ind_cellZones,polygons_patch=[],[],[],[],{},{},[]
    # read native OF polyMesh
    # when calculate Sf or return trimesh, read x,y,z
    if(calSf | triMesh):
        points=readPoints(caseDir)
        x,y,z=points['x'],points['y'],points['z']
    faces=readFaces(caseDir)
    boundaries=readBoundary(caseDir)
    owners=readOwner(caseDir)
    cellZones=readCellZones(caseDir)
    # 1. read fields
    for fieldName in fieldNames:
        fields[fieldName]=readField(caseDir,timeName,fieldName,regionName=patchName)
    # 2. cell zones
    nCells=owners.max()+100 #BUG: fix it later
    cellNames=np.array([' '*50]*nCells)
    for name in cellZones:
        cellNames[cellZones[name]]=name
    db_patch=boundaries[patchName]
    faces_patch=np.array(range(db_patch['startFace'],db_patch['startFace']+db_patch['nFaces']))
    for cellZoneName in cellZoneNames:
        ownerCells_patch=owners[faces_patch]
        cellNames_patch=cellNames[ownerCells_patch]
        ind_cellZones[cellZoneName]=(cellNames_patch == cellZoneName)
    # # polygons
    faces['index']=np.array(faces['index'])
    polygons_patch=faces['index'][faces_patch]
    # # 3. calculate face normal vector
    if(calSf):
        for face in polygons_patch:
            x_face,y_face,z_face=x[face],y[face],z[face]
            xyz_poly=np.hstack((x_face.reshape(-1,1),y_face.reshape(-1,1)))
            xyz_poly=np.hstack((xyz_poly,z_face.reshape(-1,1)))
            area,nf=poly_area_normal(xyz_poly)
            Sf.append(area*nf)
    # polydata to tri data
    triData=polyData2TriData(polygons_patch,x=x,y=y,z=z,cellFields=fields,pointFields=[])
    return {'x':x,'y':y,'z':z,'triData':triData,'polygons':polygons_patch,'fields':fields,'cellZones':ind_cellZones,'Sf':Sf}

def plotMeshTopology(ax,caseDir, ind_cell=None,index_intFace=None, meshGenerator='blockMesh',**kwargs):
    caseDir=os.path.abspath(caseDir)
    points=readPoints(caseDir)
    x,y,z=points['x'],points['y'],points['z']
    faces = readFaces(caseDir)
    owners = readOwner(caseDir)
    neighbours = readNeighbour(caseDir)
    cells=getCells_(owners, neighbours)
    boundaries, index_internalFaces, name_faces = readBoundary_(caseDir,len(faces['nNodes']))
    faces['name']=name_faces
    # get empty patch name
    name_emptyPatch=[]
    for name_patch, type_patch in zip(boundaries['name'],boundaries['type']):
        if(type_patch=='empty'):
            name_emptyPatch.append(name_patch)
    # 1. plot front face (rectangle) of each cell and cell index in the rect center
    cells_poly=[0]*len(cells)
    label='Cell index'
    for i in range(0,len(cells)):
        x_cell,y_cell=[],[]
        for face in cells[i]:
            if(name_faces[face] in name_emptyPatch):
                cells_poly[i]=faces['index'][face]
                x_cell,y_cell=x[faces['index'][face]],y[faces['index'][face]]
                ax.plot(x_cell.mean(), y_cell.mean(),'o',mfc='lightskyblue',mec='k',ms=15,label=label,**kwargs)
                ax.text(x_cell.mean(), y_cell.mean(),str('%d'%(i)), va='center',ha='center')
                label=None
                break
    # 2. plot face of the front patch, startFace and nFaces of a patch can be found in constant/polyMesh/boundary file
    # index_face_front=115+11
    # ax.fill(x[faces['index'][index_face_front]],y[faces['index'][index_face_front]],fc='lightgray', alpha=0.5, label='The %d$_{th}$ face on front patch'%(index_face_front))
    # 3. plot all internal face 
    for i in range(0,len(index_internalFaces)):
        index_face_internal = index_internalFaces[i]
        if(i==0):
            label='Internal face: %d'%(len(index_internalFaces))
        else:
            label=None
        index_points_internalFace=faces['index'][index_face_internal]
        x_face,y_face,z_face=x[index_points_internalFace], y[index_points_internalFace],z[index_points_internalFace]
        ax.plot(x_face,y_face,'k',label=label,lw=1)
        # norm=np.cross([x_face[1] - x_face[0], y_face[1] - y_face[0], z_face[1] - z_face[0]],
        #              [x_face[2] - x_face[1], y_face[2] - y_face[1], z_face[2] - z_face[1]])
        # norm=norm[0:2]/np.sqrt(np.sum(norm**2))
        # theta=90-np.arccos(norm[1])/np.pi*180
        ax.text(x_face.mean(),y_face.mean(),'%d'%(index_face_internal),va='center',ha='center',color='k',bbox={'color':'lightgray'})
    # 4. plot all boundary patches
    for name, patchType, patchIndex,lc in zip(boundaries['name'],boundaries['type'],boundaries['index'],plt.rcParams['axes.prop_cycle'].by_key()['color']):

        if(name in name_emptyPatch): # skip front and back patches, this is a 2D case!!!
            continue
        for i in range(0,len(patchIndex)):
            index_face_patch = patchIndex[i]
            if(i==0):
                label='%s(%s): %d'%(name,patchType,len(patchIndex))
            else:
                label=None
            index_points_patchFace=faces['index'][index_face_patch]
            x_tmp,y_tmp=x[index_points_patchFace], y[index_points_patchFace]
            ax.plot(x_tmp,y_tmp,color=lc,label=label,**kwargs)
            rot= 90 if(x_tmp.min()==x_tmp.max()) else 0
            ax.text(x_tmp.mean(),y_tmp.mean(),'%d'%(index_face_patch),va='center',ha='center',rotation=rot, color=lc, bbox={'color':'lightgray'}, alpha=0.5)
    # 5. plot a internal face and marker its owner and neighbour cell
    index_intFace = int(len(index_internalFaces)/2) if (index_intFace==None) else index_intFace
    ax.plot(x[faces['index'][index_intFace]], y[faces['index'][index_intFace]],'r', label='The %d$_{th}$ internal face'%(index_intFace),**kwargs)
    # print(owners[index_intFace],neighbours[index_intFace])
    ax.fill(x[cells_poly[owners[index_intFace]]],  y[cells_poly[owners[index_intFace]]], fc='dodgerblue',label='Owner cell of face %d'%(index_intFace))
    ax.fill(x[cells_poly[neighbours[index_intFace]]],  y[cells_poly[neighbours[index_intFace]]], fc='purple',label='Neighbour cell of face %d'%(index_intFace))

    # 6. plot a cell and its neighbour cells and faces
    ind_cell= int(len(cells_poly)/2) if(ind_cell==None) else ind_cell
    x_cell,y_cell=x[cells_poly[ind_cell]], y[cells_poly[ind_cell]]
    ax.fill(x_cell,y_cell, label='The %d$_{th}$ cell'%(ind_cell),fc='limegreen')
    len_diag=np.sqrt(np.sum(np.array([x_cell.max() - x_cell.min(), y_cell.max()-y_cell.min()])**2))
    ax.text(x_cell.mean()-len_diag/8,y_cell.mean()-len_diag/8, '$C$', color='w', fontsize=16, fontweight='bold', ha='center', va='center')
    # 6.1 faces own to cell
    index_local_cells, index_local_faces=1,1
    faces_own_to_cell = np.where(owners==ind_cell)[0]
    for face in faces_own_to_cell:
        lc='r'
        x_face,y_face,z_face=x[faces['index'][face]], y[faces['index'][face]], z[faces['index'][face]]
        # calculate normal vector of the face
        norm=np.cross([x_face[1] - x_face[0], y_face[1] - y_face[0], z_face[1] - z_face[0]],
                     [x_face[2] - x_face[1], y_face[2] - y_face[1], z_face[2] - z_face[1]])
        norm=norm[0:2]/np.sqrt(np.sum(norm**2))
        if(faces['name'][face] in name_emptyPatch):
            continue
        if(faces['name'][face]=='internal'):
            x_cell_w,y_cell_w=x[cells_poly[neighbours[face]]], y[cells_poly[neighbours[face]]]
            ax.fill(x_cell_w,y_cell_w,fc='gray',alpha=0.8)
            dist_cells = np.sqrt(np.sum((np.array([x_cell_w.mean(), y_cell_w.mean()]) - np.array([x_cell.mean(), y_cell.mean()]))**2))
            ax.text(x_cell_w.mean()-norm[0]*dist_cells*0.2, y_cell_w.mean()-norm[1]*dist_cells*0.2,
                    '$F_{\mathregular{%d}}$'%(index_local_cells),va='center',ha='center',color='w', fontsize=16)
            index_local_cells = index_local_cells+1
        if(faces['name'][face] in boundaries['name']): # point out if face is a boundary face
            lc='cyan'
            # norm=-norm
        lf,=ax.plot(x_face,y_face, color=lc,**kwargs)
        # plot the normal vector of the face
        xy_cf=np.array([x_face.mean(), y_face.mean()])
        len_arrow=np.sqrt(np.sum(np.array([x_face.max()-x_face.min(),y_face.max()-y_face.min()])**2))/3
        ax.annotate("", xy=xy_cf+norm*len_arrow/2, xytext=xy_cf-norm*len_arrow/2,arrowprops=dict(arrowstyle="->",color=lf.get_color()))
        # plot local face index
        norm_orth=[-norm[1],norm[0]] # one orthogonal vector of the norm vector
        ax.text(xy_cf[0]+norm_orth[0]*len_arrow/2, xy_cf[1]+norm_orth[1]*len_arrow/2, '$f_{%d}$'%(index_local_faces), fontsize=16, ha='center', va='center',color='w')
        index_local_faces = index_local_faces+1
    # 6.2 faces neighbour to cell
    faces_neighbour_to_cell=np.where(neighbours==ind_cell)[0]
    for face in faces_neighbour_to_cell:
        if(faces['name'][face] in name_emptyPatch):
            continue
        x_face,y_face,z_face=x[faces['index'][face]], y[faces['index'][face]], z[faces['index'][face]]
        norm=np.cross([x_face[1] - x_face[0], y_face[1] - y_face[0], z_face[1] - z_face[0]],
                     [x_face[2] - x_face[1], y_face[2] - y_face[1], z_face[2] - z_face[1]])
        norm=norm[0:2]/np.sqrt(np.sum(norm**2))
        x_cell_n,y_cell_n=x[cells_poly[owners[face]]], y[cells_poly[owners[face]]]
        ax.fill(x_cell_n, y_cell_n,fc='gray',alpha=0.8)
        dist_cells = np.sqrt(np.sum((np.array([x_cell_n.mean(), y_cell_n.mean()]) - np.array([x_cell.mean(), y_cell.mean()]))**2))
        ax.text(x_cell_n.mean()+norm[0]*dist_cells*0.2, y_cell_n.mean()+norm[1]*dist_cells*0.2,
                '$F_{\mathregular{%d}}$'%(index_local_cells),va='center',ha='center',color='w', fontsize=16)
        index_local_cells = index_local_cells+1
        lf,=ax.plot(x_face,y_face, color='b',**kwargs)
        # plot the normal vector of the face
        xy_cf=np.array([x_face.mean(), y_face.mean()])
        len_arrow=np.sqrt(np.sum(np.array([x_face.max()-x_face.min(),y_face.max()-y_face.min()])**2))/3
        ax.annotate("", xy=xy_cf+norm*len_arrow/2, xytext=xy_cf-norm*len_arrow/2,arrowprops=dict(arrowstyle="->",color=lf.get_color()))
        norm_orth=[-norm[1],norm[0]] # one orthogonal vector of the norm vector
        ax.text(xy_cf[0]+norm_orth[0]*len_arrow/2, xy_cf[1]+norm_orth[1]*len_arrow/2, '$f_{%d}$'%(index_local_faces), fontsize=16, ha='center', va='center',color='w')
        index_local_faces = index_local_faces+1
    
    return x,y,z,cells_poly,faces,boundaries,owners,neighbours