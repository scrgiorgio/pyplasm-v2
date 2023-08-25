using LinearAlgebra
using Test
import Base.:(==)
import Base.:*
import Base.size
import Base.transpose

using PyCall
spatial = pyimport_conda("scipy.spatial", "scipy") # the second argument is the conda package name


# /////////////////////////////////////////////////////////////
function ComputeTriangleNormal(p0, p1, p2)
	 p0 = vcat(p0, zeros(3 - length(p0)))
	 p1 = vcat(p1, zeros(3 - length(p1)))
	 p2 = vcat(p2, zeros(3 - length(p2)))
	 a = [p1[i] - p0[i] for i in 1:3]
	 b = [p2[i] - p0[i] for i in 1:3]
	 ret = cross(a, b)
	 N = norm(ret)
	 if N == 0.0
		  N = 1.0
	 end
	 return [ret[i] / N for i in 1:3]
end

# /////////////////////////////////////////////////////////////
function GoodTetOrientation(v0, v1, v2, v3)
	 v0 = vcat(v0, zeros(3 - length(v0)))
	 v1 = vcat(v1, zeros(3 - length(v1)))
	 v2 = vcat(v2, zeros(3 - length(v2)))
	 v3 = vcat(v3, zeros(3 - length(v3)))
	 a = [v3[i] - v1[i] for i in 1:3]
	 b = [v2[i] - v1[i] for i in 1:3]
	 c = [v0[i] - v1[i] for i in 1:3]
	 n = cross(a, b)
	 return dot(n, c) > 0
end

# /////////////////////////////////////////////////////////////
mutable struct BoxNd
	 p1::Vector{Float64}
	 p2::Vector{Float64}

	function BoxNd(dim::Int)
		p1 = [+floatmax(Float64) for I in 1:dim]
		p2 = [-floatmax(Float64) for I in 1:dim]	
		new(p1, p2)
	end
	
	function BoxNd(p1::Vector{Float64},p2::Vector{Float64})
		@assert length(p1) == length(p2)
		p1 = copy(p1)
		p2 = copy(p2)
		new(p1, p2)
	end

end


toList(self::BoxNd) = [copy(self.p1), copy(self.p2)]

function valid(self::BoxNd)
	 for i in 1:dim(self)
		  if self.p1[i] > self.p2[i]
				return false
		  end
	 end
	 return true
end

==(box1::BoxNd, box2::BoxNd) = isa(box1, typeof(box2)) && box1.p1 == box2.p1 && box1.p2 == box2.p2

function fuzzyEqual(box1::BoxNd, box2::BoxNd, Epsilon=1e-4)
	 if !(isa(box2, typeof(box1)) && dim(box1) == dim(box2))
		  return false
	 end
	 p1 = [abs(a - b) <= Epsilon for (a, b) in zip(box1.p1, box2.p1)]
	 p2 = [abs(a - b) <= Epsilon for (a, b) in zip(box1.p2, box2.p2)]
	 return !(false in p1) && !(false in p2)
end

Base.show(io::IO, self::BoxNd) = print(io, "BoxNd(", repr(self.p1), ", ", repr(self.p2), ")")

dim(self::BoxNd) = length(self.p1)

function size(self::BoxNd)
	 return [To - From for (From, To) in zip(self.p1, self.p2)]
end

function center(self::BoxNd)
	 return [0.5 * (From + To) for (From, To) in zip(self.p1, self.p2)]
end

function addPoint(self::BoxNd, point::Vector{Float64})
	 for i in 1:dim(self)
		self.p1[i] = min(self.p1[i], point[i])
		self.p2[i] = max(self.p2[i], point[i])
	 end
	 return self
end

function addPoints(self::BoxNd, points)
	 for point in points
		  addPoint(self, point)
	 end
	 return self
end

function addBox(box1::BoxNd, box2::BoxNd)
	 return addPoint(box1, box2.p1).addPoint(box2.p2)
end

# /////////////////////////////////////////////////////////////
mutable struct MatrixNd
	T::Matrix{Float64}

	function MatrixNd(dim::Int=0)
		T = Matrix{Float64}(I, dim, dim)
		new(T)
	end

	function MatrixNd(other::MatrixNd)
		new(copy(other.T))
	end	
	
	function MatrixNd(T::Matrix{Float64})
		new(copy(T))
	end

	function MatrixNd(arg::Vector{Vector{Float64}})
		T = reduce(vcat, arg')
		new(T)
	end

	function MatrixNd(arg::Vector{Vector{Int64}})
		T = reduce(vcat, arg')
		new(T)
	end	

end


Base.getindex(self::MatrixNd, args...) = getindex(self.T, args...)
Base.setindex!(self::MatrixNd, args...) = setindex!(self.T, args...)

==(matrix1::MatrixNd, matrix2::MatrixNd) = isa(matrix1, typeof(matrix2)) && matrix1.T == matrix2.T

Base.show(io::IO, self::MatrixNd) = print(io, "MatrixNd(", repr(toList(self)), ")")

function isIdentity(self::MatrixNd)
	 return self.T == I
end

toList(self::MatrixNd) = [self.T[R,:] for R in 1:size(self.T,1)]

function transpose(self::MatrixNd)
	 return MatrixNd(Matrix{Float64}(transpose(self.T)))
end

function invert(self::MatrixNd)
	 return MatrixNd(inv(self.T))
end

dim(self::MatrixNd) = size(self.T, 1)

function embed(self::MatrixNd, target_dim)
	current_dim=dim(self)
	 if target_dim <= current_dim
		  return self
	 end
	 ret = MatrixNd(target_dim)
	 ret.T[1:current_dim, 1:current_dim] = self.T
	 return ret
end

function adjoin(matrix1::MatrixNd, matrix2::MatrixNd)
	M, N = dim(matrix1), dim(matrix2)
	T=M + N - 1
	ret = MatrixNd(T)
	ret[2:M, 2:M] = matrix1[2:M, 2:M]
	for I in 2:M
		ret[I, 1] = matrix1[I, 1]
		ret[1, I] = matrix1[1, I]
	end
	ret[M+1:T, M+1:T] = matrix2[2:N, 2:N]
	for I in 2:N
		ret[I+M-1, 1] = matrix2[I, 1]
		ret[1, I+M-1] = matrix2[1, I]
	end
	return ret
end

function *(matrix1::MatrixNd, matrix2::MatrixNd)
	 return MatrixNd(matrix1.T * matrix2.T)
end

function transformPoint(self::MatrixNd, point::Vector{Float64})
	 point = self.T * [1.0; point; zeros(dim(self) - length(point) - 1)]
	 return [point[i, 1] / point[1, 1] for i in 2:dim(self)]
end

function translate(vt)
	 T = MatrixNd(length(vt) + 1)
	 for I in 2:dim(T)
		  T[I, 1] = vt[I-1]
	 end
	 return T
end

function scale(vs)
	 T = MatrixNd(length(vs) + 1)
	 for I in 2:dim(T)
		  T[I, I] = vs[I-1]
	 end
	 return T
end

function rotate(i, j, angle)
	i+=1
	j+=1
	T = MatrixNd(max(i, j))
	T[i, i] = +cos(angle)
	T[i, j] = -sin(angle)
	T[j, i] = +sin(angle)
	T[j, j] = +cos(angle)
	return T
end


# /////////////////////////////////////////////////////////////
mutable struct BuildMkPol
	 db::Dict{Vector{Float64}, Int}
	 points::Vector{Vector{Float64}}
	 hulls::Vector{Vector{Int}}

	# constructor
	function BuildMkPol()
		self=new(Dict{Vector{Float64}, Int}(),Vector{Vector{Float64}}(),Vector{Vector{Int}}())
		return self
	end

	# constructor
	function BuildMkPol(points::Vector{Vector{Float64}}, hulls::Vector{Vector{Int}}=Vector{Vector{Int}}())

		self=BuildMkPol()

		if isempty(points)
			 return self
		end
		
		if isempty(hulls)
			 hulls = [collect(1:length(points))]
		end

		pdim = length(points[1])	

		for hull in hulls

			 if isempty(hull)
				  continue
			 end

			 if pdim == 1
				  @assert length(hull) >= 2
				  box = BoxNd(1)
				  for idx in hull
						box = addPoint(box, points[idx])
				  end
				  addHull(self, Vector{Vector{Float64}}([box.p1, box.p2]))
			 else
				  hull_points = [points[idx] for idx in hull]

				  try
						h = spatial.ConvexHull([points[idx] for idx in hull])
						hull_points = [Vector{Float64}(h.points[idx+1,:]) for idx in h.vertices]
				  catch
				  end
				  addHull(self, hull_points)
			 end
		end

		return self

  end	

end


function addPoint(self::BuildMkPol, p::Vector{Float64})
	 idx = get(self.db, p, 0)
	 if idx>=1
		  return idx
	 else
		  idx = length(self.db) + 1
		  self.db[p] = idx
		  push!(self.points, p)
		  return idx
	 end
end

function addHull(self::BuildMkPol, points::Vector{Vector{Float64}})
	push!(self.hulls, [addPoint(self, p) for p in points])
end

dim(self::BuildMkPol) = isempty(self.points) ? 0 : length(self.points[1])

function box(self::BuildMkPol)
	 ret = BoxNd(dim(self))
	 if !isempty(self.points)
		  addPoints(ret, self.points)
	 end
	 return ret
end

Base.show(io::IO, self::BuildMkPol) = print(io, "BuildMkPol(points=", repr(self.points), ", hulls=", repr(self.hulls), ")")

function toSimplicialForm(self::BuildMkPol)

	 if isempty(self.points) || isempty(self.hulls)
		  return self
	 end
	 pdim = dim(self)
	 
	 if pdim <= 1
		  return self
	 end
	 
	 ret = BuildMkPol()
	 for hull in self.hulls
		  if length(hull) <= pdim + 1
				addHull(ret, [self.points[idx] for idx in hull])
		  else
				try
					 d = spatial.Delaunay([self.points[idx] for idx in hull])
					 for simplex in [ d.simplices[R,:] for R in 1:size(d.simplices,1)]
							simplex_points = [Vector{Float64}(d.points[idx+1,:]) for idx in simplex]
						  addHull(ret, simplex_points)
					 end
				catch 
				end
		  end
	 end
	 fixOrientation!(ret)
	 return ret
end

function fixOrientation!(self::BuildMkPol)
	 pdim = dim(self)
	 if pdim != 2 && pdim != 3
		  return
	 end
	 
	 if pdim == 2
		  fixed = Vector{Vector{Int}}()
		  for simplex in self.hulls
				if length(simplex) == 3
					 p0 = self.points[simplex[1]]
					 p1 = self.points[simplex[2]]
					 p2 = self.points[simplex[3]]
					 n = ComputeTriangleNormal(p0, p1, p2)
					 if n[3] < 0
						  simplex = [simplex[3], simplex[2], simplex[1]]
					 end
				end
				push!(fixed, simplex)
		  end
		  self.hulls = fixed
	 end
	 
	 if pdim == 3
		  fixed = Vector{Vector{Int}}()
		  for simplex in self.hulls
				if length(simplex) == 4
					 p0 = self.points[simplex[1]]
					 p1 = self.points[simplex[2]]
					 p2 = self.points[simplex[3]]
					 p3 = self.points[simplex[4]]
					 if !GoodTetOrientation(p0, p1, p2, p3)
						  simplex = [simplex[3], simplex[2], simplex[1], simplex[4]]
					 end
				end
				push!(fixed, simplex)
		  end
		  self.hulls = fixed
	 end
	 return self
end

function getBatches(self::BuildMkPol)
	 dim = dim(self)
	 simplicial_form = toSimplicialForm(self)
	 
	 points = GLBatch(GLBatch.POINTS)
	 lines = GLBatch(GLBatch.LINES)
	 triangles = GLBatch(GLBatch.TRIANGLES)
	 
	 for hull in simplicial_form.hulls
		  hull_dim = length(hull)
		  if hull_dim == 1
				p0 = self.points[hull[1]]
				push!(points.vertices, p0)
		  elseif hull_dim == 2
				p0 = self.points[hull[1]]
				p1 = self.points[hull[2]]
				push!(lines.vertices, p0)
				push!(lines.vertices, p1)
		  elseif hull_dim == 3
				p0 = self.points[hull[1]]
				p1 = self.points[hull[2]]
				p2 = self.points[hull[3]]
				n = ComputeTriangleNormal(p0, p1, p2)
				push!(triangles.vertices, p0)
				push!(triangles.normals, n)
				push!(triangles.vertices, p1)
				push!(triangles.normals, n)
				push!(triangles.vertices, p2)
				push!(triangles.normals, n)
		  elseif hull_dim == 4
				for T in [[1, 2, 4], [1, 4, 3], [1, 3, 2], [2, 3, 4]]
					 p0 = self.points[hull[T[1]]]
					 p1 = self.points[hull[T[2]]]
					 p2 = self.points[hull[T[3]]]
					 n = ComputeTriangleNormal(p0, p1, p2)
					 push!(triangles.vertices, p0)
					 push!(triangles.normals, n)
					 push!(triangles.vertices, p1)
					 push!(triangles.normals, n)
					 push!(triangles.vertices, p2)
					 push!(triangles.normals, n)
				end
		  else
				throw("Cannot handle geometry with dim > 3")
		  end
	 end
	 
	 ret = GLBatch[]
	 if !isempty(points.vertices)
		  push!(ret, points)
	 end
	 if !isempty(lines.vertices)
		  push!(ret, lines)
	 end
	 if !isempty(triangles.vertices)
		  push!(ret, triangles)
	 end
	 return ret
end


# /////////////////////////////////////////////////////////////
mutable struct Hpc
	 T::MatrixNd
	 childs::Union{Vector{Hpc}, Vector{BuildMkPol}}
	 properties::Dict{Any, Any}
	 
	 # constructor
	 function Hpc(T::MatrixNd=MatrixNd(0), childs:: Union{Vector{Hpc}, Vector{BuildMkPol}}=[], properties=Dict())

		  self = new()
		  self.childs = childs
		  self.properties = properties
		  if length(childs) > 0
				Tdim = maximum([dim(child) for child in childs]) + 1
				self.T = embed(T, Tdim)
		  else
				self.T = T
		  end
		  return self
	 end
end
	 
function MkPol(points::Vector{Vector{Float64}}, hulls::Vector{Vector{Int}}=Vector{Vector{Int}}())
	obj=BuildMkPol(points, hulls)
	return Hpc(MatrixNd(), [obj])
end

function Struct(pols::Vector{Hpc})
	return Hpc(MatrixNd(), pols)
end

function Cube(dim::Int, From::Float64=0.0, To::Float64=1.0)
	@assert dim>=1
	points = [[From],[To]]
	for I in 2:dim
		a=[ [p;From] for p in points]
		b=[ [p;To  ] for p in points]
		points =  [a;b]
	end
	return MkPol(points)
end

function Simplex(dim::Int)
	@assert dim>=1
	points = [ [0.0 for _ in 1:dim] ]
	for I in 1:dim
		point=[0.0 for _ in 1:dim]
		point[I]=1.0
		push!(points,point)
	end
	return MkPol(points)
end

function Join(pols::Vector{Hpc})
	points = []
	for (T, properties, obj) in toList(Hpc(MatrixNd(), pols))
		points = [transformPoint(T,p) for p in obj.points]
	end
	return MkPol(points)
end

function Quote(sequence::Vector{Float64})
	pos = 0.0
	points = [[pos]]
	hulls = Vector{Vector{Int}}()
	for value in sequence
		next = pos + abs(value)
		push!(points, [next])
		if value >= 0
				push!(hulls, [length(points)-1, length(points)])
		end
		pos = next
	end
	return MkPol(points, hulls)
end

function Base.show(io::IO, self::Hpc)
	print(io, "Hpc(", self.T, ", ", self.childs, ", ", self.properties, ")")
end

function dim(self::Hpc)
	return dim(self.T) - 1
end

function toList(self::Hpc)
	ret = []
	Tdim = dim(self) + 1
	stack = [[MatrixNd(Tdim), Dict(), self]]
	while !isempty(stack)
		T, properties, node = pop!(stack)
		if isa(node, Hpc)
				T = T * embed(node.T, Tdim)
				if !isempty(node.properties)
					properties = copy(properties)
					for (key, value) in node.properties
						properties[key] = value
					end
				end
				for child in node.childs
					push!(stack, [T, properties, child])
				end
		else
				push!(ret, [T, properties, node])
		end
	end
	return ret
end

function box(self::Hpc)
	box = BoxNd(dim(self))
	for (T, properties, obj) in toList(self)
		addPoints(box, [transformPoint(T,p) for p in obj.points])
	end
	return box
end

function Transform(self::Hpc, T::MatrixNd)
	return Hpc(T, [self])
end

function Translate(self::Hpc, vt::Vector{Float64})
	return Hpc(translate(vt), [self])
end

function Scale(self::Hpc, vs::Vector{Float64})
	return Hpc(scale(vs), [self])
end

function Rotate(self::Hpc, i::Int, j::Int, angle::Float64)
	return Hpc(rotate(i, j, angle), [self])
end

function Power(a::Hpc, b::Hpc)
	childs = Vector{Hpc}()
	for (T2, properties2, obj2) in toList(b)
		for (T1, properties1, obj1) in toList(a)
				points = Vector{Vector{Float64}}()
				for py in obj2.points
					for px in obj1.points
						push!(points, [px;py])
					end
				end
				hulls = Vector{Vector{Int}}()
				nx, ny = length(obj1.points), length(obj2.points)
				for hy in obj2.hulls
					for hx in obj1.hulls
						hull=Vector{Int}()
						for A in hx
							for B in hy
								push!(hull,1+((B-1)*nx + (A-1)))
							end
						end
						push!(hulls,hull)
					end
				end
				T = adjoin(T1, T2)
				push!(childs, Hpc(T, [BuildMkPol(points, hulls)]))
		end
	end
	return Hpc(MatrixNd(), childs)
end

function UkPol(self::Hpc)
	points, hulls = [], []
	for (T, properties, obj) in toList(self)
		offset = length(points)
		push!(points, [transformPoint(T,p) for p in obj.points]...)
		push!(hulls, [[offset + idx for idx in hull] for hull in obj.hulls]...)
	end
	return [points, hulls]
end

function View(self::Hpc, title::String="")
	batches = []
	for (T, properties, obj) in toList(self)
		T = embed(T, 4)
		T3D = [
				T[1, 1], T[1, 2], T[1, 3], T[1, 0],
				T[2, 1], T[2, 2], T[2, 3], T[2, 0],
				T[3, 1], T[3, 2], T[3, 3], T[3, 0],
				T[0, 1], T[0, 2], T[0, 3], T[0, 0]
		]
		for batch in obj.getBatches()
				prependTransformation(batch, T3D)
				writeProperties(batch, properties)
				push!(batches, batch)
		end
	end
	GLView(batches, title=title)
end

function MapFn(self::Hpc, fn)
	childs = Vector{Hpc}()
	for (T, properties, obj) in toList(self)
		obj = toSimplicialForm(obj)
		points = [fn(transformPoint(T,p)) for p in obj.points]
		hulls = obj.hulls
		push!(childs, Hpc(MatrixNd(), [BuildMkPol(points, hulls)], properties))
	end
	ret = Hpc(MatrixNd(), childs)
	return ret
end

function ToBoundaryForm(self::Hpc)
	POINTDB = Dict()
	faces = []
	for (T, properties, obj) in toList(self)
		obj = toSimplicialForm(obj)
		points, hulls = [transformPoint(T,p) for p in obj.points], obj.hulls
		dim = length(points[1])
		mapped = Dict()
		for P in 1:length(points)
				point = tuple(points[P]...)
				if !(point in POINTDB)
					POINTDB[point] = length(POINTDB)
				end
				mapped[P] = POINTDB[point]
		end
		for hull in hulls
				bfaces = []
				if length(hull) < (dim + 1)
					bfaces = [collect(1:length(hull))]
				elseif length(hull) == (dim + 1)
					if dim == 0
						bfaces = [[1]]
					elseif dim == 1
						bfaces = [[1], [2]]
					elseif dim == 2
						bfaces = [[1, 2], [2, 3], [3, 1]]
					elseif dim == 3
						bfaces = [[1, 2, 4], [1, 4, 3], [1, 3, 2], [2, 3, 4]]
					else
						error("not supported")
					end
				else
					error("internal error")
				end
				for face in bfaces
					push!(faces, [mapped[hull[it]] for it in face])
				end
		end
	end
	num_occurrence = Dict()
	for face in faces
		key = tuple(sort(face))
		if !(key in num_occurrence)
				num_occurrence[key] = 0
		end
		num_occurrence[key] += 1
	end
	faces = [face for face in faces if num_occurrence[tuple(sort(face))] == 1]
	points = [nothing for _ in 1:length(POINTDB)]
	for (point, num) in POINTDB
		points[num] = point
	end
	ret = MkPol(points, faces)
	return ret
end


function TestComputeNormal()
	 @test ComputeTriangleNormal([0, 0, 0], [1, 0, 0], [0, 1, 0]) == [0.0, 0.0, 1.0]
	 @test ComputeTriangleNormal([0, 0, 0], [0, 1, 0], [0, 0, 1]) == [1.0, 0.0, 0.0]
	 @test ComputeTriangleNormal([0, 0, 0], [1, 0, 1], [1, 0, 0]) == [0.0, 1.0, 0.0]
end

function TestGoodTet()
	 @test GoodTetOrientation([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]) == true
	 @test GoodTetOrientation([0, 0, 0], [0, 1, 0], [1, 0, 0], [0, 0, 1]) == false
end

function TestBox()
	 x1, y1, z1 = 1, 2, 3
	 x2, y2, z2 = 4, 5, 6
	 b = BoxNd([1.0, x1, y1, z1], [1.0, x2, y2, z2])
	 @test dim(b) == 4
	 @test b.p1 == [1.0, x1, y1, z1]
	 @test b.p2 == [1.0, x2, y2, z2]
	 @test valid(b) == true
	 @test b == BoxNd(toList(b)...)
	 @test size(b) == [0.0, x2 - x1, y2 - y1, z2 - z1]
	 @test center(b) == [1.0, (x1 + x2) / 2, (y1 + y2) / 2, (z1 + z2) / 2]
	 x1, y1, z1 = -1, -2, -3
	 x2, y2, z2 = 10, 20, 30
	 addPoint(b, [1.0, x1, y1, z1])
	 addPoint(b, [1.0, x2, y2, z2])
	 @test b == BoxNd([1.0, x1, y1, z1], [1.0, x2, y2, z2])
end

function TestMat()
	T = MatrixNd(4)
	v = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
	T[1, 1] += 0.0
	@test T == MatrixNd(v)
	@test dim(T) == 4
	@test T[1, 1] == 1.0
	@test isIdentity(T) == true
	@test toList(T) == v
	@test transpose(T) == T
	@test invert(T) == T
	@test embed(T, 5) == MatrixNd(5)
	@test T * T == T
	
	# transform point does not want homo coordinate
	p = [2.0, 3.0, 4.0]
	@test transformPoint(T, p) == p
	@test translate([10, 20, 30]) == MatrixNd([[1, 0, 0, 0], [10, 1, 0, 0], [20, 0, 1, 0], [30, 0, 0, 1]])
	@test scale([10, 20, 30]) == MatrixNd([[1, 0, 0, 0], [0, 10, 0, 0], [0, 0, 20, 0], [0, 0, 0, 30]])
	angle = π / 2
	@test rotate(1, 2, angle) == MatrixNd([[1, 0, 0], [0.0, cos(angle), -sin(angle)], [0.0, sin(angle), cos(angle)]])
end

function TestMkPol()

	# 1D
	points=[[0.0],[1.0],[2.0],   [8.0],[9.0],[10.0]]
	hulls=[[1,2,3],[4,5,6]]
	obj=BuildMkPol(points,hulls)
	@test dim(obj)==1
	@test box(obj)==BoxNd([0.0],[10.0])
	obj=toSimplicialForm(obj)
	@test length(obj.points)==4
	@test length(obj.hulls)==2
	@test box(obj)==BoxNd([0.0],[10.0])

	# 2D
	points = [[0.0, 0.0], [0.2, 0.2], [1.0, 0.0], [0.3, 0.3], [1.0, 1.0], [0.4, 0.4], [0.0, 1.0], [0.5, 0.5], [0.2, 0.8]]
	hulls = [collect(1:length(points))]
	obj = BuildMkPol(points, hulls)
	@test dim(obj) == 2
	@test box(obj) == BoxNd([0.0,0.0], [1.0,1.0])
	obj = toSimplicialForm(obj)
	@test length(obj.points) == 4
	@test length(obj.hulls) == 2
	@test box(obj) == BoxNd([0.0,0.0], [1.0,1.0])

	# 3D
	points = [
		[0.0, 0.0, 0.0],
		[1.0, 0.0, 0.0],
		[1.0, 1.0, 0.0],
		[0.0, 1.0, 0.0],
		[0.0, 0.0, 1.0],
		[1.0, 0.0, 1.0],
		[1.0, 1.0, 1.0],
		[0.0, 1.0, 1.0],
		[0.1, 0.1, 0.1],
		[0.2, 0.2, 0.2],
		[0.3, 0.3, 0.3]
	]
	hulls = [collect(1:length(points))]
	obj = BuildMkPol(points, hulls)
	@test dim(obj) == 3
	@test box(obj) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
	obj = toSimplicialForm(obj)
	@test dim(obj) == 3
	@test box(obj) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
	@test length(obj.points) == 8
	@test length(obj.hulls) == 6
end

function TestHpc()
	points = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.2, 0.2], [0.3, 0.3], [0.4, 0.4], [0.5, 0.5], [0.2, 0.8]]
	hulls = [collect(1:length(points))]
	obj = MkPol(points, hulls)
	@test dim(obj) == 2
	@test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])

	# struct
	obj = Struct([obj])
	@test dim(obj) == 2
	@test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])
	
	# cube
	obj = Cube(3, 0.0, 1.0)
	@test dim(obj) == 3
	@test box(obj) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
	@test box(obj) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
	v = toList(obj)
	@test length(v) == 1
	(T, properties, obj) = v[1]
	obj = toSimplicialForm(obj)
	@test dim(obj) == 3
	@test box(obj) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
	@test length(obj.points) == 8
	@test length(obj.hulls) == 6
	
	# simplex
	obj = Simplex(3)
	@test dim(obj) == 3
	@test box(obj) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
	v = toList(obj)
	@test length(v) == 1
	(T, properties, obj) = v[1]
	obj = toSimplicialForm(obj)
	@test dim(obj) == 3
	@test box(obj) == BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
	@test length(obj.points) == 4
	@test length(obj.hulls) == 1
	
	# quote
	obj = Quote([1.0, -2.0, 1.0])
	@test box(obj) == BoxNd([0.0], [4.0])
	(T, properties, obj) = toList(obj)[1]
	@test length(obj.points) == 4
	@test obj.hulls == [[1, 2], [3, 4]]
	
	# join
	obj = Join([Cube(2, 0.0, 1.0), Cube(2, 0.0, 1.0)])
	(T, properties, obj) = toList(obj)[1]
	obj = toSimplicialForm(obj)
	@test length(obj.points) == 4
	@test length(obj.hulls) == 2
	
	# Transform
	obj = Transform(Cube(2, 0.0, 1.0),MatrixNd(3))
	(T, properties, obj) = toList(obj)[1]
	@test dim(T) == 3
	@test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])

	# power
	obj = Power(Cube(2, 0.0, 1.0), Cube(1, 10.0, 20.0))
	l = toList(obj)
	@test length(l) == 1
	(T, properties, obj) = l[1]
	@test length(obj.points) == 8
	@test box(obj) == BoxNd([0.0, 0.0, 10.0], [1.0, 1.0, 20.0])

	# Scale
	obj = Scale(Cube(2, 0.0, 1.0),[1.0, 1.0])
	@test dim(obj) == 2
	@test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])

	# Rotate
	obj = Rotate(Cube(2, 0.0, 1.0), 1,2, 0.0)
	@test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])

	# Translate
	obj = Translate(Cube(2, 0.0, 1.0),[0.0, 0.0])
	@test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])

	# MapFn
	function fn(p)
		return [p[I]+0.0 for I in 1:length(p)]
	end

	pbj=MapFn(Cube(2, 0.0, 1.0),fn)
	@test box(obj) == BoxNd([0.0, 0.0], [1.0, 1.0])

	# TODO 
	# UkPol
	# ToBoundaryForm
	 
end

if abspath(PROGRAM_FILE) == @__FILE__
	TestComputeNormal()
	TestGoodTet()
	TestBox()
	TestMat()
	TestMkPol()
	TestHpc()
	println("all test ok")
end