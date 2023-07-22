import math,sys,types
import numpy

from pyplasm.geom import *

# viewer dependencies
from OpenGL.GL			import *
from OpenGL.GLU		  import *
import glfw

# ////////////////////////////////////////////
def Flatten(l):
	return functools.reduce(lambda x, y: x + y, l)

# /////////////////////////////////////////////////////////////////////
class GLVertexBuffer:

	destroy_list = []	
	
	# constructor
	def __init__(self, data):
		data=Flatten(data)
		self.id=glGenBuffers(1)
		glBindBuffer(GL_ARRAY_BUFFER, self.id)
		glBufferData (GL_ARRAY_BUFFER, len(data)*4, numpy.array(data, dtype="float32"), GL_STATIC_DRAW)
		glBindBuffer(GL_ARRAY_BUFFER, 0)			
		
	# destructor
	def __del__(self):
		GLVertexBuffer.destroy_list.append(self.id)

	# bind
	def bind(self, location, num_components):
		glBindBuffer(GL_ARRAY_BUFFER, self.id)
		glVertexAttribPointer(location,num_components,GL_FLOAT,GL_FALSE,0,None)
		glEnableVertexAttribArray(location)	
		glBindBuffer(GL_ARRAY_BUFFER, 0)	
	

# /////////////////////////////////////////////////////////////////////
class GLBatch:

	POINTS    = 0x0000
	LINES     = 0x0001
	TRIANGLES  = 0x0004

	# constructor
	def __init__(self,primitive):
		self.primitive=primitive
		self.T=Matrix4d() 
		self.vertex_array=None
		self.vertices=None; self.gpu_vertices=None
		self.normals =None; self.gpu_normals=None
		self.colors  =None; self.gpu_colors=None

		# used only for direct rendering
		self.ambient  =(0.2,0.2,0.2,1.0) 
		self.diffuse  =(0.8,0.8,0.8,1.0) 
		self.specular =(0.1,0.1,0.1,1.0)
		self.emission =(0.0,0.0,0.0,0.0)
		self.shininess=100.0

		self.properties={}

	# __str__
	def __str__(self):
		return str({"type":"GLBatch", "primitive":self.primitive,"vertices":self.vertices,"normals":self.normals})       

	# dim  
	def dim(self):
		return 3 	

	# getBoundingBox
	def getBoundingBox(self):
		return Box3d().addPoints([self.T.transformPoint(p) for p in self.vertices])
	
	@staticmethod
	def GetBoundingBox(batches):
		ret=Box3d()
		for batch in batches:
			ret.addBox(batch.getBoundingBox())
		return ret if ret.valid() else Box3d(Point3d([-1,-1,-1]),Point3d([+1,+1,+1]))
			
	# prependTransformation
	def prependTransformation(self,T):
		self.T=T * self.T
		
	# writeProperties
	def writeProperties(self,properties):
		for key,value in properties.items():
			self.properties[key]=value

	# glRenderBatchShaderMode
	def glRenderBatchShaderMode(self,PROJECTION,MODELVIEW, lightpos):

			modelview=MODELVIEW * self.T
			normal_matrix = modelview.inverted().transposed().dropW()
				
			for bRenderLines in [False,True]:

				shader=GLShaderProgram.CreatePhongShader(
					lighting_enabled=True if self.normals else False,
					color_attribute_enabled=True if self.colors else False)

				if self.primitive==GLBatch.TRIANGLES:
					glPolygonOffset(-1.0,-1.0)
					glEnable(GL_POLYGON_OFFSET_LINE)
					if bRenderLines:
						glPolygonMode(GL_FRONT_AND_BACK,GL_LINE)	
						shader=GLShaderProgram.CreatePhongShader(lighting_enabled=False,color_attribute_enabled=False)
				else:
					if bRenderLines: 
						continue
					
				a_position          = glGetAttribLocation(shader.program_id, "a_position")
				a_normal            = glGetAttribLocation(shader.program_id, "a_normal")
				a_color             = glGetAttribLocation(shader.program_id, "a_color")
				
				u_modelview_matrix  = glGetUniformLocation(shader.program_id, "u_modelview_matrix" )
				u_projection_matrix = glGetUniformLocation(shader.program_id, "u_projection_matrix")	
				u_normal_matrix     = glGetUniformLocation(shader.program_id, "u_normal_matrix" )	
				u_color             = glGetUniformLocation(shader.program_id, "u_color")
				u_light_position    = glGetUniformLocation(shader.program_id, "u_light_position")
				
				glUseProgram(shader.program_id)
				
				glUniformMatrix4fv(u_modelview_matrix ,1, GL_TRUE, numpy.array(modelview .mat   , dtype="float32"))
				glUniformMatrix4fv(u_projection_matrix,1, GL_TRUE, numpy.array(PROJECTION.mat   , dtype="float32"))
				glUniformMatrix3fv(u_normal_matrix    ,1, GL_TRUE, numpy.array(normal_matrix,     dtype="float32"))
				
				if u_light_position>=0:
					glUniform3f(u_light_position,lightpos[0]/lightpos[3],lightpos[1]/lightpos[3],lightpos[2]/lightpos[3])				
				
				if self.vertex_array is None: self.vertex_array=glGenVertexArrays(1)
				glBindVertexArray(self.vertex_array)
				
				# positions
				if self.vertices and a_position>=0: 
					if self.gpu_vertices is None: self.gpu_vertices=GLVertexBuffer(self.vertices)
					self.gpu_vertices.bind(a_position,3)
				
				# normals
				if self.normals and a_normal>=0:
					if self.gpu_normals is None: self.gpu_normals=GLVertexBuffer(self.normals)
					self.gpu_normals.bind(a_normal, 3)
				
				# colors
				if self.colors and a_color>=0:
					if self.gpu_colors is None: self.gpu_colors=GLVertexBuffer(self.colors)
					self.gpu_colors.bind(a_color, 4)
					
				# render lines
				if u_color>=0:
					if bRenderLines:
						glUniform4f(u_color,0.0,0.0,0.0,1.0)	
					else:
						glUniform4f(u_color,0.5,0.5,0.5,1.0)	
			 
				num_vertices=int(len(self.vertices))
				glDrawArrays(self.primitive, 0, num_vertices)
				
				if self.vertices and a_position>=0:
					glDisableVertexAttribArray(a_position)
				
				if self.normals and a_normal>=0:
					glDisableVertexAttribArray(a_normal)
					
				if self.colors and a_color>=0:	
					glDisableVertexAttribArray(a_color)
					
				glBindVertexArray(0)
				glUseProgram(0)	
				
				glDepthMask(True)
				glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)
				glDisable(GL_POLYGON_OFFSET_LINE)		


	# glRenderBatchImmediateMode
	def glRenderBatchImmediateMode(self):

		if (self.colors is None):
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT  , self.ambient)
			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE  , self.diffuse)
			glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR , self.specular)
			glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, self.shininess )
			if (self.diffuse[3]<1):
				glEnable(GL_BLEND)
				glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
			 
		if (not self.vertices is None):
			if self.gpu_vertices is None: self.gpu_vertices=GLVertexBuffer(self.vertices)
			glBindBuffer(GL_ARRAY_BUFFER, self.gpu_vertices.id)
			glVertexPointer(3,GL_FLOAT,0,None)
			glBindBuffer(GL_ARRAY_BUFFER,0)  
			glEnableClientState(GL_VERTEX_ARRAY)
	
		if (not  self.normals  is None):
			if self.gpu_normals is None: self.gpu_normals=GLVertexBuffer(self.normals)
			glBindBuffer(GL_ARRAY_BUFFER, self.gpu_normals.id)
			glNormalPointer(GL_FLOAT,0,None)
			glBindBuffer(GL_ARRAY_BUFFER,0)  
			glEnableClientState(GL_NORMAL_ARRAY)      
	
		if (not  self.colors  is None):
			if self.gpu_colors is None: self.gpu_colors=GLVertexBuffer(self.colors)
			glBindBuffer(GL_ARRAY_BUFFER, self.gpu_colors.id)
			glColorPointer(3,GL_FLOAT,0,None)
			glBindBuffer(GL_ARRAY_BUFFER,0)  
			glEnableClientState(GL_COLOR_ARRAY)
			glEnable(GL_COLOR_MATERIAL)
	
		glPushMatrix()

		matrix=self.T.transposed().toList()

		glMultMatrixf(matrix)
		num_vertices=len(self.vertices)
		glDrawArrays(self.primitive, 0, num_vertices)
		glPopMatrix()
	
		if (self.diffuse[3]<1) :
			glDisable(GL_BLEND)
	
		if (self.colors):
			glDisableClientState(GL_COLOR_ARRAY)
			glDisable(GL_COLOR_MATERIAL)
	
		if (self.normals) :
			glDisableClientState(GL_NORMAL_ARRAY)
	
		if (self.vertices):
			glDisableClientState(GL_VERTEX_ARRAY)
					

# //////////////////////////////////////////////////////////////////
def GLCuboid(box):
	
	points=[
		Point3d(box.p1.x,box.p1.y,box.p1.z),
		Point3d(box.p2.x,box.p1.y,box.p1.z),
		Point3d(box.p2.x,box.p2.y,box.p1.z),
		Point3d(box.p1.x,box.p2.y,box.p1.z),
		Point3d(box.p1.x,box.p1.y,box.p2.z),
		Point3d(box.p2.x,box.p1.y,box.p2.z),
		Point3d(box.p2.x,box.p2.y,box.p2.z),
		Point3d(box.p1.x,box.p2.y,box.p2.z)
	]
	
	faces=[[0, 1, 2, 3],[3, 2, 6, 7],[7, 6, 5, 4],[4, 5, 1, 0],[5, 6, 2, 1],[7, 4, 0, 3]]
	
	vertices=[]
	for face in faces:
		p3,p2,p1,p0 = points[face[0]],points[face[1]],points[face[2]],points[face[3]] # reverse order
		vertices.append([p0.x,p0.y,p0.z])
		vertices.append([p1.x,p1.y,p1.z])
		vertices.append([p2.x,p2.y,p2.z])
		vertices.append([p0.x,p0.y,p0.z])
		vertices.append([p2.x,p2.y,p2.z])
		vertices.append([p3.x,p3.y,p3.z])
		
	ret=GLBatch(GLBatch.TRIANGLES)
	ret.vertices = vertices
	ret.normals  = []

	for I in range(0,len(vertices),3):
		p0=Point3d(*vertices[I+0])
		p1=Point3d(*vertices[I+1])
		p2=Point3d(*vertices[I+2])
		n=Point3d.crossProduct(p1-p0,p2-p0).normalized()
		ret.normals.append([n.x,n.y,n.z])	
		ret.normals.append([n.x,n.y,n.z])
		ret.normals.append([n.x,n.y,n.z])

	return ret

# //////////////////////////////////////////////////////////////////
def GLAxis(box):
	ret=GLBatch(GLBatch.LINES)
	p0=[box.p1[0],box.p1[1],box.p1[2]]
	p1=[box.p2[0],box.p2[1],box.p2[2]]
	p0=[p0[0],p0[1],p0[2]]
	px=[p1[0],p0[1],p0[2]]
	py=[p0[0],p1[1],p0[2]]
	pz=[p0[0],p0[1],p1[2]]
	r=[1.0,0.0,0.0,1.0]
	g=[0.0,1.0,0.0,1.0]
	b=[0.0,0.0,1.0,1.0]
	ret.vertices=[p0,px,p0,py,p0,pz]
	ret.colors=[r,r,g,g,b,b]
	return ret


# /////////////////////////////////////////////////////////////
class FrustumMap:

	# constructor
	def __init__(self, viewport,projection,modelview):

		self.viewport=Matrix4d(
			viewport[2] / 2.0 ,                   0,       0, viewport[0] + viewport[2] / 2.0,
												0,  viewport[3] / 2.0,       0, viewport[1] + viewport[3] / 2.0,
												0,                  0, 1 / 2.0, 1 / 2.0)	
		self.projection=projection
		self.modelview=modelview
		
		self.inv_viewport   = self.viewport.inverted()
		self.inv_projection = self.projection.inverted()
		self.inv_modelview  = self.modelview.inverted()
		
	# projectPoint
	def projectPoint(self,p3):
		p4=self.viewport.T.transformPoint(self.projection.T.transformPoint(self.modelview .T.transformPoint(Point4d(p3[0],p3[1],p3[2],1.0))))
		return Point3d(p4[0]/p4[3],p4[1]/p4[3],p4[2]/p4[3])

	# unprojectPoint
	def unprojectPoint(self,p2, Z = 0.0):
		p4 = self.inv_modelview.transformPoint(self.inv_projection.transformPoint(self.inv_viewport.transformPoint(Point4d(p2[0], p2[1], Z, 1.0))))
		if p4.w==0: p4.w=1
		return Point3d(p4[0]/p4[3],p4[1]/p4[3],p4[2]/p4[3])	



# /////////////////////////////////////////////////////////////////////x
cached_shaders={}

class GLShaderProgram:
	def __init__(self,vertex, fragment):
		self.program_id = glCreateProgram()
		glAttachShader(self.program_id, GLShaderProgram.CreateShader(GL_VERTEX_SHADER  ,vertex  ))
		glAttachShader(self.program_id, GLShaderProgram.CreateShader(GL_FRAGMENT_SHADER,fragment))
		glLinkProgram(self.program_id)
		if glGetProgramiv(self.program_id, GL_LINK_STATUS) != GL_TRUE:
			raise RuntimeError('Error linking program: %s' % (glGetProgramInfoLog(self.program_id),))

	@staticmethod
	def CreateShader(type,source):
		id = glCreateShader(type)
		glShaderSource(id, source)
		glCompileShader(id)
		if glGetShaderiv(id, GL_COMPILE_STATUS) != GL_TRUE:
			raise RuntimeError('Shader compilation failed: %s' % (glGetShaderInfoLog(id),))
		return id
	
	@staticmethod
	def CreatePhongShader(lighting_enabled,color_attribute_enabled):
		global cached_shaders

		key=f"phong/{lighting_enabled}/{color_attribute_enabled}"
		if key in cached_shaders:
			return cached_shaders[key]
		
		defines=""
		if lighting_enabled:
			defines+="#define LIGHTING_ENABLED 1\n"
			
		if color_attribute_enabled:
			defines+="#define COLOR_ATTRIBUTE_ENABLED 1\n"
			
		v=defines + """
		uniform mat4 u_modelview_matrix;
		uniform mat4 u_projection_matrix;
		uniform vec4 u_color;

		attribute  vec4 a_position;

		#if LIGHTING_ENABLED
		attribute  vec3 a_normal;
		#endif

		#if COLOR_ATTRIBUTE_ENABLED
		attribute vec4 a_color;
		#endif

		#if LIGHTING_ENABLED
		uniform mat3 u_normal_matrix;
		uniform vec3 u_light_position;
		varying vec3 v_normal;
		varying vec3 v_light_dir;
		varying vec3 v_eye_vec;
		#endif

		#if COLOR_ATTRIBUTE_ENABLED
			varying vec4 v_color;
		#endif

		void main() 
		{
			vec4 eye_pos= u_modelview_matrix * a_position;
			
		#if LIGHTING_ENABLED	
			v_normal = u_normal_matrix * a_normal;
			vec3 vVertex = vec3(u_modelview_matrix * a_position);
			v_light_dir  = normalize(u_light_position - vVertex);
			v_eye_vec    = normalize(-vVertex);
		#endif	

		#if COLOR_ATTRIBUTE_ENABLED
			v_color=a_color;
		#endif
			
			gl_Position = u_projection_matrix * eye_pos;
		}
		"""

		f=defines + """
		uniform vec4 u_color;

		#if LIGHTING_ENABLED
		varying vec3 v_normal;
		varying vec3 v_light_dir;
		varying vec3 v_eye_vec;
		#endif

		#if COLOR_ATTRIBUTE_ENABLED
			varying vec4 v_color;
		#endif

		void main() 
		{
			vec4 frag_color=u_color; 
			
			#if LIGHTING_ENABLED
			vec3 N = normalize(v_normal   );
			vec3 L = normalize(v_light_dir);
			vec3 E = normalize(v_eye_vec  );

			vec4  u_material_ambient  = vec4(0.2,0.2,0.2,1.0);
			vec4  u_material_diffuse  = vec4(0.8,0.8,0.8,1.0);
			vec4  u_material_specular = vec4(0.1,0.1,0.1,1.0);
			float u_material_shininess=100.0;	
			
			if(gl_FrontFacing)
			{
				frag_color = u_material_ambient;
				float NdotL = abs(dot(N,L));
				if (NdotL>0.0)
					{
					vec3 R = reflect(-L, N);
					float NdotHV = abs(dot(R, E));
					frag_color += u_material_diffuse * NdotL;
					frag_color += u_material_specular * pow(NdotHV,u_material_shininess);
				}
			}
			else
			{
				frag_color = u_material_ambient;
				float NdotL = abs(dot(-N,L));
				if (NdotL>0.0);
				{
					vec3 R = reflect(-L, -N);
					float NdotHV=abs(dot(R, E));
					frag_color += u_material_diffuse * NdotL;
					frag_color += u_material_specular * pow(NdotHV,u_material_shininess);
				}
			}
		#endif

		#if COLOR_ATTRIBUTE_ENABLED
			frag_color =v_color;
		#endif

			gl_FragColor = frag_color;
		}
		"""
		cached_shaders[key]=GLShaderProgram(v,f)	
		return cached_shaders[key]		
		


# /////////////////////////////////////////////////////////////////////
class Viewer:

	# constructor
	def __init__(self,batches,title=""):
		self.button=0
		self.batches=batches
		self.W=1024
		self.H=768
		self.fov=60.0
		self.button=0
		self.shader=None
		self.shader_mode=False
		self.title=title
		
		self.guessBestPosition()
		
		glfw.init()
		
		#glfw.window_hint(glfw.CONTEXT_VERSION_MAJOR, 3)
		#glfw.window_hint(glfw.CONTEXT_VERSION_MINOR, 2)
		#glfw.window_hint(glfw.OPENGL_FORWARD_COMPAT, GL_TRUE)
		#glfw.window_hint(glfw.OPENGL_PROFILE, glfw.OPENGL_CORE_PROFILE)
		
		self.win = glfw.create_window(self.W, self.H, "pyplasm" if not title else title, None, None)	
		glfw.make_context_current(self.win)			
		
		glfw.set_window_size_callback  (self.win, lambda win, w,h : self.onResize(w,h))  
		glfw.set_key_callback          (self.win, lambda win, key, scancode, action, mods: self.onKeyPress(key, scancode, action, mods))
		glfw.set_cursor_pos_callback   (self.win, lambda win, x,y : self.onMouseMove(x,y))  
		glfw.set_mouse_button_callback (self.win, lambda win, button,action,mods: self.onMouseButton(button,action,mods))	
		glfw.set_scroll_callback       (self.win, lambda win, dx,dy: self.onMouseWheel(dx,dy))


	# guessBestPosition
	def guessBestPosition(self):

		self.bounds=GLBatch.GetBoundingBox(self.batches)
		size=self.bounds.p2-self.bounds.p1
		
		self.center = (self.bounds.p1+self.bounds.p2)*0.5
		self.pos = self.center+2.1*size
		self.vup = Point3d(0,0,1)

		self.center_of_rotation = self.center
		self.rotation=Quaternion()

		dir=(self.center-self.pos).normalized()

		print(f"{self.title} bounds={self.bounds}")

		camera_plane=Point4d(dir[0],dir[1],dir[2],-1.0*Point3d.dotProduct(dir,self.pos))

		fPoint = Point3d(
			self.bounds.p2[0] if (dir[0] >= 0) else self.bounds.p1[0], 
			self.bounds.p2[1] if (dir[1] >= 0) else self.bounds.p1[1], 
			self.bounds.p2[2] if (dir[2] >= 0) else self.bounds.p1[2])
		
		nPoint = Point3d(
			self.bounds.p1[0] if (dir[0] >= 0) else self.bounds.p2[0], 
			self.bounds.p1[1] if (dir[1] >= 0) else self.bounds.p2[1], 
			self.bounds.p1[2] if (dir[2] >= 0) else self.bounds.p2[2])
		
		def GetDistance(p):
			return p[0]*camera_plane[0]+p[1]*camera_plane[1]+p[2]*camera_plane[2]+1.0*camera_plane[3]
		
		nDist = GetDistance(nPoint)
		fDist = GetDistance(fPoint)
		self.zNear = 0.1 if nDist <= 0 else nDist
		self.zFar  = 2 * fDist - nDist
		self.walk_speed=(self.zFar-self.zNear)*0.004

		self.redisplay()

	# run
	def run(self):
		while not glfw.window_should_close(self.win) :
			self.glRenderAll()
			glfw.swap_buffers(self.win)
			glfw.poll_events()
		glfw.terminate()		
				
		
	# redisplay
	def redisplay(self):
		pass			

	# onResize
	def onResize(self, width, height):
		if (width==0 or height==0): return
		self.W = width
		self.H = height
		self.redisplay()
		
	# onMouseButton
	def onMouseButton(self,button,action,mods):
		button={0:1,1:3,2:2}[button]
		
		if action == glfw.PRESS and self.button==0:
			self.button=button
			self.redisplay()		
			return
			
		if action==glfw.RELEASE and button==self.button:
			self.button=0
			self.redisplay()		
			return
		
	# onMouseMove
	def onMouseMove(self,x,y):
		
		button=self.button
		
		if (button==0):
			self.mouse_beginx = x
			self.mouse_beginy = y
			return
			
		modelview=self.getModelview()
		
		screen1=Point3d(self.mouse_beginx,self.H-self.mouse_beginy,0)
		screen2=Point3d(x,self.H-y,0)

		if (button==1):

			 # /note: going 'inside' rotation
			Ti = self.getViewportDirectTransformation() * \
				self.getProjection() * \
				Matrix4d.LookAt(self.pos, self.center, self.vup) * \
				Matrix4d.Translate(self.center_of_rotation)
			Ti = Ti.inverted()
			a = (Ti.transformPoint(screen1).withoutHomo()).normalized()
			b = (Ti.transformPoint(screen2).withoutHomo()).normalized()
			axis = a.cross(b).normalized()
			angle = math.acos(a.dot(b))

			if axis.length() and axis.valid():
				self.rotation=Quaternion.FromAxisAndAngle(axis, angle * 1.0) * self.rotation

		elif (button==3):

			# NOTE: excluding the rotation part;
			Ti = self.getViewportDirectTransformation() * \
				self.getProjection() * \
				Matrix4d.LookAt(self.pos, self.center, self.vup) 
			Ti = Ti.inverted()

			vt = (Ti.transformPoint(screen2).withoutHomo() - Ti.transformPoint(screen1).withoutHomo()) * 1.0;

			if vt.valid():
				self.pos=self.pos-vt
				self.center=self.center-vt

		self.mouse_beginx = x
		self.mouse_beginy = y
		self.redisplay()
		
	# onMouseWheel
	def onMouseWheel(self, dx,dy):
		K=(+1 if dy>0 else -1)*self.walk_speed
		dir=(self.center-self.pos).normalized()
		self.pos   =self.pos   +dir*K
		self.center=self.center+dir*K
		self.redisplay()

	# onKeyPress
	def onKeyPress(self,key, scancode, action, mods):
		
		if action != glfw.PRESS:
			return	
			
		x,y=0.5*self.W,0.5*self.H	
		
		if key == glfw.KEY_ESCAPE: 
			glfw.set_window_should_close(self.win, True);

			return		
		
		if (key=="+" or key=="="):
			self.walk_speed*=0.95
			return 

		if (key=='-' or key=="_"):
			self.walk_speed*=(1.0/0.95)
			return 

		if (key=='w'):
				dir=self.unproject(x,y)
				self.pos=self.pos+dir*self.walk_speed
				self.redisplay()		
				return 

		if (key=='s'):
				dir=self.unproject(x,y)
				self.pos=self.pos-dir*self.walk_speed
				self.redisplay()		
				return 
				
		if (key=='a'	or key==glfw.KEY_LEFT):
			right=Point3d.crossProduct(self.dir,self.vup).normalized()
			self.pos=self.pos-right*self.walk_speed
			self.redisplay()		
			return 

		if (key=='d' or key==glfw.KEY_RIGHT):
			right=Point3d.crossProduct(self.dir,self.vup).normalized()
			self.pos=self.pos+right*self.walk_speed
			self.redisplay()	
			return				 

		if (key==glfw.KEY_UP):
			self.pos=self.pos+self.vup*self.walk_speed
			self.redisplay()	
			return 

		if (key==glfw.KEY_DOWN):
			self.pos=self.pos-self.vup*self.walk_speed
			self.redisplay()	
			return 		
			
	# getModelview
	def getModelview(self):
		return \
			Matrix4d.LookAt(self.pos, self.center, self.vup) * \
			Matrix4d.RotateAroundCenter(self.center_of_rotation, self.rotation)
				
	# getProjection
	def getProjection(self):
		return Matrix4d.Perspective(self.fov,self.W/float(self.H),self.zNear,self.zFar)

	# getViewport
	def getViewport(self):
		return [ 0, 0, self.W, self.H]

	# getViewportDirectTransformation
	def getViewportDirectTransformation(self):
		viewport=self.getViewport()
		sx = viewport[2] / 2.0; ox = viewport[0] + viewport[2]  / 2.0
		sy = viewport[3] / 2.0; oy = viewport[1] + viewport[3] / 2.0
		sz = 1 / 2.0; oz = 1 / 2.0

		return Matrix4d([
			sx, 0, 0, ox,
			0, sy, 0, oy,
			0, 0, sz, oz,
			0, 0, 0, 1])


	# project
	def project(self,pos):
		viewport,projection,modelview = self.getViewport(), self.getProjection(), self.getModelview()
		map=FrustumMap(viewport,projection,modelview)
		return map.projectPoint(pos.x,pos.y,pos.z)
		
	# unproject
	def unproject(self,x,y):
		viewport,projection,modelview = self.getViewport(), self.getProjection(), self.getModelview()
		map=FrustumMap(viewport,projection,modelview)
		P1=map.unprojectPoint(x,self.H-y,-1)
		P2=map.unprojectPoint(x,self.H-y,+1)
		return (P2-P1).normalized() 
		
	# glRenderAll
	def glRenderAll(self):	

		# destruction of old vertex buffers
		for it in GLVertexBuffer.destroy_list:
			glDeleteBuffers(1, GLuint(it))   
		GLVertexBuffer.destroy_list=[]

		glEnable(GL_DEPTH_TEST)
		glDepthFunc(GL_LEQUAL)
		glDisable(GL_CULL_FACE)
		glClearDepth(1.0)
		glClearColor(0.3,0.4,0.5, 0.00)

		glViewport(0,0,self.W,self.H)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		
		projection = self.getProjection()
		modelview  = self.getModelview ()

		if not self.shader_mode:
			glEnable(GL_NORMALIZE)
			glShadeModel(GL_SMOOTH)
			glDisable(GL_COLOR_MATERIAL)
			glEnable(GL_LIGHTING)

			glMatrixMode(GL_PROJECTION)
			glLoadMatrixf(projection.transposed().toList())
			
			glMatrixMode(GL_MODELVIEW)  
			glLoadMatrixf(modelview.transposed().toList())

			glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,0)
			glLightfv(GL_LIGHT0,GL_AMBIENT  , [+1.00,+1.00,+1.00,+1.00])
			glLightfv(GL_LIGHT0,GL_DIFFUSE  , [+1.00,+1.00,+1.00,+1.00])
			glLightfv(GL_LIGHT0,GL_SPECULAR , [+1.00,+1.00,+1.00,+1.00])

			# set default light
			if (True):
				glEnable(GL_LIGHTING)
				glEnable(GL_LIGHT0)
				glLightfv(GL_LIGHT0, GL_POSITION, [self.pos.x,self.pos.y,self.pos.z,1.0])
				glLightfv(GL_LIGHT0, GL_AMBIENT,  [1, 1, 1, 1])
				glLightfv(GL_LIGHT0, GL_DIFFUSE,  [1, 1, 1, 1])
				glLightfv(GL_LIGHT0, GL_SPECULAR, [1, 1, 1, 1])			

	
		if self.shader_mode:
			lightpos=modelview.transformPoint(self.pos)

			for batch in self.batches:
				batch.glRenderBatchShaderMode(projection,modelview,lightpos)
		else:

			# draw axis
			if (True):
				glDisable(GL_LIGHTING)
				glLineWidth(3)
				glBegin(GL_LINES)
				glColor3f(1,0,0);glVertex3i(0,0,0);glVertex3i(1,0,0)
				glColor3f(0,1,0);glVertex3i(0,0,0);glVertex3i(0,1,0)
				glColor3f(0,0,1);glVertex3i(0,0,0);glVertex3i(0,0,1)
				glEnd()
				glLineWidth(1)
				glEnable(GL_LIGHTING)
				
			for batch in self.batches:
				batch.glRenderBatchImmediateMode()
					
			 # draw lines
			if True:
				glDisable(GL_LIGHTING)
				glDepthMask(False) 
				glLineWidth(2)
				glColor3f(0.5,0.5,0.5)
				glPolygonMode(GL_FRONT_AND_BACK,GL_LINE)
				for batch in self.batches:
					if batch.dim()>=2: batch.glRenderBatchImmediateMode()
				glDepthMask(True)
				glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)
				glLineWidth(1) 


# ///////////////////////////////////////////////////////////////////  
def GLView(batches,title=""):
	viewer=Viewer(batches,title=title)
	viewer.run()

# //////////////////////////////////////////////////////////////////////////////////////
def Main():

	box=Box3d(Point3d(0,0,0),Point3d(1,1,1))
	box=Box3d(Point3d(-1.000000000000001, -11.999999999999998, -2.999999999999999), Point3d(32.25, 29.999999999999996, 11.390625))

	batches=[
		GLCuboid(box),
		GLAxis(Box3d(Point3d(-1,-1,-1),Point3d(2,2,2)))
	]
	viewer=Viewer(batches)
	viewer.run()
	 
# //////////////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
	Main()