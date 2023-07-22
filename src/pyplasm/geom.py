import sys,functools,math,copy
import numpy

TET_ORIENTED_TRIANGLES=[[0,1,3],[0,3,2],[0,2,1],[1,2,3]]

# /////////////////////////////////////////////////////////////////////////
class Point3d:
	
	def __init__(self,x:float=0.0, y:float=0.0, z:float=0.0):
		self.x = float(x)
		self.y = float(y)
		self.z = float(z)

	# toList
	def toList(self):
		return [self.x,self.y,self.z]
		
	# __repr__
	def __repr__(self):
		return f"Point3d({self.x}, {self.y}, {self.z})"
	
	# valid
	def valid(self): 
		return True

	def __getitem__(self, key):
		return {0:self.x,1:self.y,2:self.z}[key]
			
	def __setitem__(self, key, value):
		if key==0: 
			self.x=value
		elif key==1:
			self.y=value
		elif key==2:
			self.z=value
		else:
			raise Exception("iternal error")
	
	def __add__(self, o):
		return Point3d((self.x + o.x), (self.y + o.y), (self.z + o.z))
	
	def __sub__(self, o):
		return Point3d((self.x - o.x), (self.y - o.y), (self.z - o.z))
		
	def __mul__(self, vs):
		return Point3d(self.x*vs,self.y*vs,self.z*vs)		

	def __rmul__(self, vs):
		return Point3d(self.x*vs,self.y*vs,self.z*vs)
	
	def __iadd__(self, o):
		self.x += o.x
		self.y += o.y
		self.z += o.z
		return self
	
	def __isub__(self, o):
		self.x -= o.x
		self.y -= o.y
		self.z -= o.z
		return self
	
	def __neg__(self):
		return Point3d(-self.x, -self.y, -self.z)
		
	def length(self):
		return math.sqrt((self.x*self.x) + (self.y*self.y) + (self.z*self.z))

	def normalized(self):
		len=self.length()
		if len==0: len=1.0
		return Point3d((self.x / len), (self.y / len), (self.z / len))

	@staticmethod
	def dotProduct(a, b):
		return (a.x * b.x) + (a.y * b.y) + (a.z * b.z)

	# dot
	def dot(self,b):
		return Point3d.dotProduct(self,b)

	@staticmethod
	def crossProduct(a,b): 
		return Point3d(
			a.y * b.z - b.y * a.z, 
			a.z * b.x - b.z * a.x, 
			a.x * b.y - b.x * a.y)
	
	# cross
	def cross(self,b):
		return Point3d.crossProduct(self,b)
	
# /////////////////////////////////////////////////////////////////////////
class Point4d:
	
	def __init__(self,x:float=0.0, y:float=0.0, z:float=0.0, w:float=0.0):
		self.x = float(x)
		self.y = float(y)
		self.z = float(z)
		self.w = float(w)	

	# toList
	def toList(self):
		return [self.x,self.y,self.z,self.w]

	# __repr__
	def __repr__(self):
		return f"Point4d({self.x}, {self.y}, {self.z}, {self.w})"
	
	def __getitem__(self, key):
		return {0:self.x, 1:self.y, 2:self.z, 3:self.w}[key]
			
	def __setitem__(self, key, value):
		if key==0: 
			self.x=value
		elif key==1:
			self.y=value
		elif key==2:
			self.z=value
		elif key==3:
			self.w=value
		else:
			raise Exception("iternal error")	
	
	# withoutHomo
	def withoutHomo(self):
		w=self.w
		if w==0: w=1.0
		return Point3d(self.x/w,self.y/w,self.z/w)

# /////////////////////////////////////////////////////////////////////
class Quaternion:

	def __init__(self,fW:float =1.0,fX:float=0.0,fY:float=0.0,fZ:float=0.0):
		self.w = fW
		self.x = fX
		self.y = fY
		self.z = fZ

	# toList
	def toList(self):
		return [self.w,self.x,self.y,self.z]		

	# __repr__
	def __repr__(self):
		return f"Quaternion({self.w}, {self.x}, {self.y}, {self.z})" 

	# getAxis
	def getAxis(self)  :
		return Point3d(self.x, self.y, self.z).normalized() if self.x!=0.0 or self.y!=0.0 or self.z!=0.0 else Point3d(0,0,1)

	# getAngle
	def getAngle(self) :
		w=self.w
		if w<-1.0: w=1.0
		if w>+1.0: w=1.0
		return 2.0*math.acos(w)

	@staticmethod
	def FromAxisAndAngle(axis,angle):
		axis=axis.normalized()
		halfangle = 0.5*angle
		fSin = math.sin(halfangle)
		return Quaternion(math.cos(halfangle),fSin*axis.x,fSin*axis.y,fSin*axis.z)
		
	def toMatrix4d(self):
		fTx  = 2.0*self.x
		fTy  = 2.0*self.y
		fTz  = 2.0*self.z
		fTwx = fTx*self.w
		fTwy = fTy*self.w
		fTwz = fTz*self.w
		fTxx = fTx*self.x
		fTxy = fTy*self.x
		fTxz = fTz*self.x
		fTyy = fTy*self.y
		fTyz = fTz*self.y
		fTzz = fTz*self.z
		return Matrix4d([
			1.0-(fTyy+fTzz),    (fTxy-fTwz),    (fTxz+fTwy),0.0,
					(fTxy+fTwz),1.0-(fTxx+fTzz),    (fTyz-fTwx),0.0,
					(fTxz-fTwy),    (fTyz+fTwx),1.0-(fTxx+fTyy),0.0,
			0.0,0.0,0.0,1.0])

	# product
	def __mul__(self, b):
		a=self
		return Quaternion(
			a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
			a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
			a.w * b.y + a.y * b.w + a.z * b.x - a.x * b.z,
			a.w * b.z + a.z * b.w + a.x * b.y - a.y * b.x)		
	
	
		
# /////////////////////////////////////////////////////////////////////
class Box3d:

	# constructor
	def __init__(self,p1 : Point3d =Point3d(),p2 : Point3d=Point3d()):
		self.p1=Point3d(*p1.toList())
		self.p2=Point3d(*p2.toList())
		
	@staticmethod
	def invalid():
		m,M=sys.float_info.min,sys.float_info.max
		return Box3d(
			Point3d(M,M,M),
			Point3d(m,m,m))	

	# valid
	def valid(self):
		for i in range(3):
			if (self.p1[i]>self.p2[i]): return False
		return True

	# __repr__
	def __repr__(self):
		return f"Box3d({repr(self.p1)}, {repr(self.p2)})"

	# addPoint
	def addPoint(self,point):
		for i in range(3):
			self.p1[i]=min(self.p1[i],point[i])
			self.p2[i]=max(self.p2[i],point[i])
		return self

	# addBox
	def addBox(self,box):
		self.addPoint(box.p1)
		self.addPoint(box.p2)
		return self

	# addPoints
	def addPoints(self,points):
		for point in points:
				self.addPoint(point)
		return self
		


# //////////////////////////////////////////////////////////////////////
class BoxNd:
	
	# constructor
	def __init__(self,*arg):

		# dimension
		if len(arg)==1 and isinstance(arg[0],int):
			self.p1=[sys.float_info.max] * arg[0]
			self.p2=[sys.float_info.min] * arg[0]
			
		# two lists
		elif len(arg)==2 and isinstance(arg[0],(list, tuple)) and isinstance(arg[1],(list, tuple)):
			assert(len(arg[0])==len(arg[1]))
			self.p1=copy.copy(arg[0])
			self.p2=copy.copy(arg[1])
		 
		else:
			raise Exception("Invalid constructor call")
		
	# toList
	def toList(self):
		return [copy.copy(self.p1),copy.copy(self.p2)]

	# valid
	def valid(self):
		for i in range(self.dim()):
			if (self.p1[i]>self.p2[i]): return False
		return True
		
	# __eq__
	def __eq__(self, other):
		return (isinstance(other, self.__class__) and self.__dict__==other.__dict__)
	
	# __ne__
	def __ne__(self, other):
		return not self.__eq__(other)  
		
	#fuzzyEqual
	def fuzzyEqual(self,other,Epsilon=1e-4):
		if not isinstance(other, self.__class__) or other.dim()!=self.dim(): 
			return False
		p1=[abs(a-b)<=Epsilon for a,b in zip(self.p1,other.p1)]
		p2=[abs(a-b)<=Epsilon for a,b in zip(self.p2,other.p2)]
		return (not False in p1) and (not False in p2)
		
	# repr
	def __repr__(self):
			return f"BoxNd({repr(self.p1)}, {repr(self.p2)})"
		
	# dim
	def dim(self):
		return len(self.p1)      
		
	# size
	def size(self):
		return [(To-From) for From,To in zip(self.p1,self.p2)]
		
	# center
	def center(self):
		return [0.5*(From+To) for From,To in zip(self.p1,self.p2)]
	
	# addPoint
	def addPoint(self,point):
		for i in range(self.dim()):
			self.p1[i]=min(self.p1[i],point[i])
			self.p2[i]=max(self.p2[i],point[i])
		return self
		
	# addPoint
	def addPoints(self,points):
		for point in points: self.addPoint(point)
		return self    
		
	# addBox
	def addBox(self,other):
		return self.addPoint(other.p1).addPoint(other.p2)
	
	
# /////////////////////////////////////////////////////////////////////
class Matrix4d:
	def __init__(self, mat=[1.0, 0.0, 0.0, 0.0,  0.0, 1.0, 0.0, 0.0,  0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 1.0]):
		self.mat = copy.copy(mat)
		
	# __repr__
	def __repr__(self):
		s=",".join([str(it) for it in self.mat])
		return f"Matrix4d([{s}])"
	
	# toList
	def toList(self):
		return copy.copy(self.mat)

	# get
	def get(self,row,col):
			return self.mat[row*4+col]

	# toQuaternion
	def toQuaternion(self):
		
		kRot=[
				[self.mat[ 0],self.mat[ 1],self.mat[ 2]],
				[self.mat[ 4],self.mat[ 5],self.mat[ 6]],
				[self.mat[ 8],self.mat[ 9],self.mat[10]]]
				
		fTrace = kRot[0][0]+kRot[1][1]+kRot[2][2]
		
		if fTrace>0.0:
			fRoot = math.sqrt(fTrace + 1.0)
			w = 0.5*fRoot
			fRoot = 0.5/fRoot
			x = (kRot[2][1]-kRot[1][2])*fRoot
			y = (kRot[0][2]-kRot[2][0])*fRoot
			z = (kRot[1][0]-kRot[0][1])*fRoot
			return Quaternion(w,x,y,z)
		else:
			s_iNext = [1, 2, 0]
			i = 0
			if ( kRot[1][1] > kRot[0][0] ) : i = 1
			if ( kRot[2][2] > kRot[i][i] ) : i = 2
			j = s_iNext[i]
			k = s_iNext[j]
			fRoot = math.sqrt(kRot[i][i]-kRot[j][j]-kRot[k][k] + 1.0)
			Q = [0,0,0]
			Q[i] = 0.5 * fRoot
			fRoot = 0.5/fRoot
			w = (kRot[k][j]-kRot[j][k])*fRoot
			Q[j] = (kRot[j][i]+kRot[i][j])*fRoot
			Q[k] = (kRot[k][i]+kRot[i][k])*fRoot
			return Quaternion(w,Q[0],Q[1],Q[2])
		

	# transposed
	def transposed(self):
		return Matrix4d([
			self.mat[ 0],self.mat[ 4],self.mat[ 8],self.mat[12],
			self.mat[ 1],self.mat[ 5],self.mat[ 9],self.mat[13],
			self.mat[ 2],self.mat[ 6],self.mat[10],self.mat[14],
			self.mat[ 3],self.mat[ 7],self.mat[11],self.mat[15],
		])

	# transformPoint
	def transformPoint(self,p):

		if isinstance(p,list) or isinstance(p,tuple):
			x,y,z=p
		elif isinstance(p,Point3d):
			x,y,z=p.x,p.y,p.z
		else:
			raise Exception(f"wrong argument {type(p)}")

		X=(self.mat[ 0]*(x)+self.mat[ 1]*(y)+self.mat[ 2]*(z)+ self.mat[ 3]*(1.0))
		Y=(self.mat[ 4]*(x)+self.mat[ 5]*(y)+self.mat[ 6]*(z)+ self.mat[ 7]*(1.0))
		Z=(self.mat[ 8]*(x)+self.mat[ 9]*(y)+self.mat[10]*(z)+ self.mat[11]*(1.0))
		W=(self.mat[12]*(x)+self.mat[13]*(y)+self.mat[14]*(z)+ self.mat[15]*(1.0))
		return Point4d(X,Y,Z,W)	
	
	# __mul__
	def __mul__(self, other):
		A=self.mat
		B=other.mat
		return Matrix4d([
			A[ 0]*B[ 0]+A[ 1]*B[ 4]+A[ 2]*B[ 8]+A[ 3]*B[12],
			A[ 0]*B[ 1]+A[ 1]*B[ 5]+A[ 2]*B[ 9]+A[ 3]*B[13],
			A[ 0]*B[ 2]+A[ 1]*B[ 6]+A[ 2]*B[10]+A[ 3]*B[14],
			A[ 0]*B[ 3]+A[ 1]*B[ 7]+A[ 2]*B[11]+A[ 3]*B[15], 
			A[ 4]*B[ 0]+A[ 5]*B[ 4]+A[ 6]*B[ 8]+A[ 7]*B[12], 
			A[ 4]*B[ 1]+A[ 5]*B[ 5]+A[ 6]*B[ 9]+A[ 7]*B[13], 
			A[ 4]*B[ 2]+A[ 5]*B[ 6]+A[ 6]*B[10]+A[ 7]*B[14], 
			A[ 4]*B[ 3]+A[ 5]*B[ 7]+A[ 6]*B[11]+A[ 7]*B[15], 
			A[ 8]*B[ 0]+A[ 9]*B[ 4]+A[10]*B[ 8]+A[11]*B[12], 
			A[ 8]*B[ 1]+A[ 9]*B[ 5]+A[10]*B[ 9]+A[11]*B[13], 
			A[ 8]*B[ 2]+A[ 9]*B[ 6]+A[10]*B[10]+A[11]*B[14], 
			A[ 8]*B[ 3]+A[ 9]*B[ 7]+A[10]*B[11]+A[11]*B[15], 
			A[12]*B[ 0]+A[13]*B[ 4]+A[14]*B[ 8]+A[15]*B[12], 
			A[12]*B[ 1]+A[13]*B[ 5]+A[14]*B[ 9]+A[15]*B[13], 
			A[12]*B[ 2]+A[13]*B[ 6]+A[14]*B[10]+A[15]*B[14], 
			A[12]*B[ 3]+A[13]*B[ 7]+A[14]*B[11]+A[15]*B[15]])

	# inverted
	def inverted(self):
		trans=self.transposed()
		m=trans.mat
		inv=[0.0]*16
		inv[ 0] =  m[5]*m[10]*m[15] - m[5]*m[11]*m[14] - m[9]*m[6]*m[15]+ m[9]*m[7]*m[14] + m[13]*m[6]*m[11] - m[13]*m[7]*m[10]
		inv[ 4] = -m[4]*m[10]*m[15] + m[4]*m[11]*m[14] + m[8]*m[6]*m[15]- m[8]*m[7]*m[14] - m[12]*m[6]*m[11] + m[12]*m[7]*m[10]
		inv[ 8] =  m[4]*m[ 9]*m[15] - m[4]*m[11]*m[13] - m[8]*m[5]*m[15]+ m[8]*m[7]*m[13] + m[12]*m[5]*m[11] - m[12]*m[7]*m[ 9]
		inv[12] = -m[4]*m[ 9]*m[14] + m[4]*m[10]*m[13] + m[8]*m[5]*m[14]- m[8]*m[6]*m[13] - m[12]*m[5]*m[10] + m[12]*m[6]*m[ 9]
		inv[ 1] = -m[1]*m[10]*m[15] + m[1]*m[11]*m[14] + m[9]*m[2]*m[15]- m[9]*m[3]*m[14] - m[13]*m[2]*m[11] + m[13]*m[3]*m[10]
		inv[ 5] =  m[0]*m[10]*m[15] - m[0]*m[11]*m[14] - m[8]*m[2]*m[15]+ m[8]*m[3]*m[14] + m[12]*m[2]*m[11] - m[12]*m[3]*m[10]
		inv[ 9] = -m[0]*m[ 9]*m[15] + m[0]*m[11]*m[13] + m[8]*m[1]*m[15]- m[8]*m[3]*m[13] - m[12]*m[1]*m[11] + m[12]*m[3]*m[ 9]
		inv[13] =  m[0]*m[ 9]*m[14] - m[0]*m[10]*m[13] - m[8]*m[1]*m[14]+ m[8]*m[2]*m[13] + m[12]*m[1]*m[10] - m[12]*m[2]*m[ 9]
		inv[ 2] =  m[1]*m[ 6]*m[15] - m[1]*m[ 7]*m[14] - m[5]*m[2]*m[15]+ m[5]*m[3]*m[14] + m[13]*m[2]*m[ 7] - m[13]*m[3]*m[ 6]
		inv[ 6] = -m[0]*m[ 6]*m[15] + m[0]*m[ 7]*m[14] + m[4]*m[2]*m[15]- m[4]*m[3]*m[14] - m[12]*m[2]*m[ 7] + m[12]*m[3]*m[ 6]
		inv[10] =  m[0]*m[ 5]*m[15] - m[0]*m[ 7]*m[13] - m[4]*m[1]*m[15]+ m[4]*m[3]*m[13] + m[12]*m[1]*m[ 7] - m[12]*m[3]*m[ 5]
		inv[14] = -m[0]*m[ 5]*m[14] + m[0]*m[ 6]*m[13] + m[4]*m[1]*m[14]- m[4]*m[2]*m[13] - m[12]*m[1]*m[ 6] + m[12]*m[2]*m[ 5]
		inv[ 3] = -m[1]*m[ 6]*m[11] + m[1]*m[ 7]*m[10] + m[5]*m[2]*m[11]- m[5]*m[3]*m[10] - m[ 9]*m[2]*m[ 7] + m[ 9]*m[3]*m[ 6]
		inv[ 7] =  m[0]*m[ 6]*m[11] - m[0]*m[ 7]*m[10] - m[4]*m[2]*m[11]+ m[4]*m[3]*m[10] + m[ 8]*m[2]*m[ 7] - m[ 8]*m[3]*m[ 6]
		inv[11] = -m[0]*m[ 5]*m[11] + m[0]*m[ 7]*m[ 9] + m[4]*m[1]*m[11]- m[4]*m[3]*m[ 9] - m[ 8]*m[1]*m[ 7] + m[ 8]*m[3]*m[ 5]
		inv[15] =  m[0]*m[ 5]*m[10] - m[0]*m[ 6]*m[ 9] - m[4]*m[1]*m[10]+ m[4]*m[2]*m[ 9] + m[ 8]*m[1]*m[ 6] - m[ 8]*m[2]*m[ 5]
		det = m[0]*inv[0] + m[1]*inv[4] + m[2]*inv[8] + m[3]*inv[12]
		if det==0: return Matrix4d()
		return Matrix4d([it * (1.0/det) for it in inv]).transposed() 
	
	# dropW
	def dropW(self):
		return [
			self.mat[0],self.mat[1],self.mat[ 2],
			self.mat[4],self.mat[5],self.mat[ 6],
			self.mat[8],self.mat[9],self.mat[10]]
		
	@ staticmethod
	def Translate(vt):
		return Matrix4d([
				1.0, 0.0, 0.0, vt[0],
				0.0, 1.0, 0.0, vt[1],
				0.0, 0.0, 1.0, vt[2],
				0.0, 0.0, 0.0, 1.0])	
	
	# translate
	def translate(self,x,y,z):
		tmp=self* Matrix4d.Translate([x,y,z])
		self.mat=tmp.mat

	@ staticmethod
	def Scale(vs):
		return Matrix4d([
				vs[0], 0.0, 0.0, 0.0,
				0.0, vs[1], 0.0, 0.0,
				0.0, 0.0, vs[2], 0.0,
				0.0, 0.0, 0.0, 1.0])
	
	# scale
	def scale(self,x,y,z):
		tmp=self* Matrix4d.Scale([x,y,z])
		self.mat=tmp.mat
	
	@ staticmethod
	def LookAt(eye, center, up):
		forward=(center-eye).normalized()
		side   = Point3d.crossProduct(forward,up).normalized()
		up     = Point3d.crossProduct(side,forward)
		m = Matrix4d([
			side[0],up[0],-forward[0], 0.0,
			side[1],up[1],-forward[1], 0.0,
			side[2],up[2],-forward[2], 0.0,
			    0.0,  0.0,        0.0, 1.0
		])
		ret= m.transposed() * Matrix4d.Translate(-1*eye)
		return ret
	

	# getLookAt
	def getLookAt(self):
		vmat=self.inverted().mat
		pos=Point3d(  vmat[3], vmat[7], vmat[11])
		vup=Point3d(  vmat[1], vmat[5], vmat[ 9]).normalized()
		dir=Point3d( -vmat[2],-vmat[6],-vmat[10]).normalized()
		return [pos,dir,vup]

	@staticmethod
	def Rotate(q):
		return q.toMatrix4d()
	
	@staticmethod
	def RotateAroundCenter(center, q):
		if q.getAngle()==0: 
			return Matrix4d()
		else:
			return Matrix4d.Translate(center) * Matrix4d.Rotate(q) * Matrix4d.Translate(-center)
		
	# rotate
	def rotate(self,q):
		tmp=self* Matrix4d.Rotate(q)
		self.mat=tmp.mat	

	@ staticmethod
	def Perspective(fovy, aspect, zNear, zFar):
		radians =  math.radians(fovy/2.0)
		cotangent = math.cos(radians) / math.sin(radians)
		m=Matrix4d()
		m.mat[ 0] = cotangent / aspect
		m.mat[ 5] = cotangent
		m.mat[10] = -(zFar + zNear) / (zFar - zNear)
		m.mat[11] = -1
		m.mat[14] = -2 * zNear * zFar / (zFar - zNear)
		m.mat[15] = 0
		ret=m.transposed()
		return ret
	

# ////////////////////////////////////////////////////////////
# MatrixNd class (first row/column are in omogeneous coordinates!)
class MatrixNd:
	
	def __init__(self,arg=None):
		if (arg is None):  arg=0
		self.T=numpy.matrix(numpy.identity(arg,dtype=float)) if isinstance(arg,int) else numpy.matrix(arg,dtype=float)
		
	def __getitem__(self, *args):
		return self.T.__getitem__(*args)
	
	def __setitem__(self, *args):
		return self.T.__setitem__(*args)
		
	# repr
	def __repr__(self):
		
		if (self.isIdentity()):
			return 'MatrixNd('+str(self.dim())+')'  
		else:
			return 'MatrixNd('+repr(self.toList())+')'  
			
	# isIdentity
	def isIdentity(self):
		return numpy.array_equal(self.T,numpy.identity(self.dim()))
		
	# toList
	def toList(self):
		return self.T.tolist()
		
	# transpose
	def transpose(self):
		return MatrixNd(self.T.transpose())
		
	# invert
	def invert(self):
		return MatrixNd(numpy.linalg.inv(self.T))
		
	# dim
	def dim(self):
		assert(self.T.shape[0]==self.T.shape[1])
		return self.T.shape[0]     
			 
	# embed
	def embed(self,dim):
		if (self.dim()>=dim): return self
		ret=MatrixNd(dim)
		ret.T[0:self.dim(),0:self.dim()]=self.T
		return ret
		
	# adjoin
	def adjoin(self,other):
		M,N=self.dim(),other.dim()
		ret=MatrixNd(M+N-1)
	 
		ret[1:M,1:M]=self [1:,1:]
		for I in range(1,M): 
			ret[I,0]=self [I,0]
			ret[0,I]=self [0,I]
		 
		ret[M: ,M: ]=other[1:,1:] 
		for I in range(1,N): 
			ret[I+M-1,0]=other[I,0]
			ret[0,I+M-1]=other[0,I]
		
		return ret
		
	# __mul__
	def __mul__(self, other):
		assert (isinstance(other,MatrixNd))
		return MatrixNd(self.T * other.T) 
		
	# transformPoint
	def transformPoint(self,point):
		assert(isinstance(point,(list, tuple)))
		point=self.T * numpy.matrix([[1.0] + list(point) + [0.0]*(self.dim()-len(point)-1)]).transpose()
		return [point[i,0]/point[0,0] for i in range(1,self.dim())]
		
	# translate (example in 2D: vt([1.0,2.0]))
	@staticmethod
	def translate(vt):
		T=MatrixNd(len(vt)+1)
		for I in range(1,T.dim()): T[I,0]=vt[I-1]
		return T
		
	# scale (example in 2D: vs([1.0,2.0]))
	@staticmethod
	def scale(vs):
		T=MatrixNd(len(vs)+1)
		for I in range(1,T.dim()): T[I,I]=vs[I-1]
		return T

	# rotate (example in 2D i=1 j=2)
	@staticmethod
	def rotate(i,j,angle):
		T=MatrixNd(max([i,j])+1) 
		T[i,i]=+math.cos(angle) ; T[i,j]=-math.sin(angle)
		T[j,i]=+math.sin(angle) ; T[j,j]=+math.cos(angle)
		return T 
	
# /////////////////////////////////////////////////////////////////
def ComputeNormal(p0,p1,p2):
  p0=list(p0) + [0]*(3-len(p0))
  p1=list(p1) + [0]*(3-len(p1))
  p2=list(p2) + [0]*(3-len(p2))
  n=Point3d.crossProduct(Point3d(*p1)-Point3d(*p0),Point3d(*p2)-Point3d(*p0)).normalized()
  return [n.x,n.y,n.z]    
  
# /////////////////////////////////////////////////////////////////
# computeTetOrientation
#  see http://math.stackexchange.com/questions/183030/given-a-tetrahedron-how-to-find-the-outward-surface-normals-for-each-side)
#  see http://www.geuz.org/pipermail/gmsh/2012/007251.html

def GoodTetOrientation(v0,v1,v2,v3):
  v0=list(v0) + [0]*(3-len(v0));v0=Point3d(*v0)
  v1=list(v1) + [0]*(3-len(v1));v1=Point3d(*v1)
  v2=list(v2) + [0]*(3-len(v2));v2=Point3d(*v2)
  v3=list(v3) + [0]*(3-len(v3));v3=Point3d(*v3)
  n=Point3d.crossProduct(v3-v1,v2-v1)
  return Point3d.dotProduct(n,v0-v1)>0
  