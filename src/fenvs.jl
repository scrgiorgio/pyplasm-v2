PI = pi
COS = cos
LEN = length
AND = all
OR = any

using Combinatorics

include("viewer.jl")
include("hpc.jl")

# /////////////////////////////////////////////////////////////////
function ToFloat64(value)
	if isa(value,Vector)
		return [ToFloat64(it) for it in value]
	else
		return Float64(value)
	end
end

# /////////////////////////////////////////////////////////////////
function C(fun)
	function C1(arg1)
		function C2(arg2)
			return fun([arg1,arg2])
		end
		return C2
	end
	return C1
end

# /////////////////////////////////////////////////////////////////
ATAN2(l) = atan(l[2],l[1])

# /////////////////////////////////////////////////////////////////
MOD(l) = float(l[1] % l[2])

# /////////////////////////////////////////////////////////////////
function CAT(args)
	return reduce(vcat, args)
end

function ISMAT(v::Vector{Vector{Float64}})
	return true
end

function ISMAT(v::MatrixNd)
	return true
end

function ISMAT(v::Any)
	return false
end

# /////////////////////////////////////////////////////////////////
function INV(T::MatrixNd)
	return invert(T)
end

function INV(T::Vector{Vector{Float64}})
	return invert(MatrixNd(T))
end

# /////////////////////////////////////////////////////////////////
function EVERY(predicate, iterable)
	for x in iterable
		if !predicate(x)
			return false
		end
	end
	return true
end

# /////////////////////////////////////////////////////////////////
#function CURRY(fn, cargs..., ckwargs...)
#	function call_fn(fargs..., fkwargs...)
#		d = Dict(ckwargs...)
#		d = merge(d, fkwargs)
#		return fn(cargs..., fargs...; d...)
#	end
#	return call_fn
#end

# /////////////////////////////////////////////////////////////////
#function AND(list)
#	for i in list
#		if !i
#			return false
#		end
#	end
##	return true
#end

# /////////////////////////////////////////////////////////////////
ID(anyValue) = anyValue

# /////////////////////////////////////////////////////////////////
function K(AnyValue)
	function K0(obj)
		return AnyValue
	end
	return K0
end
TT = K(true)

# /////////////////////////////////////////////////////////////////
function DISTL(args)
	Element, List = args
	return [[Element, e] for e in List]
end

# /////////////////////////////////////////////////////////////////
function DISTR(args)
	List, Element = args
	return [[e,Element] for e in List]
end

# /////////////////////////////////////////////////////////////////
function COMP(Funs)
	function compose(f,g)
		function h(x)
			return f(g(x))
		end
		return h
	end
	return reduce(compose,Funs)
end

# /////////////////////////////////////////////////////////////////
function AA(f)
	function AA0(args)
		return map(f, args)
	end
	return AA0
end

Eq(x,y) = x == y

# /////////////////////////////////////////////////////////////////
function EQ(List)
	for i in List
		if !(i == List[1])
			return false
		end
	end
	return true
end

# /////////////////////////////////////////////////////////////////
function NEQ(List)
	return !EQ(List)
end

LT(a) = function(b)
	return b < a
end

LE(a) = function(b)
	return b <= a
end

GT(a) = function(b)
	return b > a
end

GE(a) = function(b)
	return b >= a
end

# /////////////////////////////////////////////////////////////////
function ISGT(args)
A , B = args
return GT(A)(B)
end

function ISLT(args)
A , B = args
return LT(A)(B)
end

function ISGE(args)
A , B = args
return GE(A)(B)
end

function ISLE(args)
	A , B = args
	return LE(A)(B)
end

function BIGGER(args)
	A , B = args
	return A >= B ? A : B
end

function SMALLER(args)
	A , B = args
	return A <= B ? A : B
end

# /////////////////////////////////////////////////////////////////
function FILTER(predicate)
	function FILTER0(sequence)
		ret=[]
		for item in sequence
			if predicate(item)
					push!(ret, item)
			end
		end
		return ret
	end
	return FILTER0
end

# /////////////////////////////////////////////////////////////////
function APPLY(args)
	f,x = args
	return f(x)
end

# /////////////////////////////////////////////////////////////////
function INSR(f)
	function INSR0(seq)
		length=length(seq)
		res = seq[end]
		for i in length-1:-1:1
			res = f([seq[i], res])
		end
		return res
	end
	return INSR0
end

# /////////////////////////////////////////////////////////////////
function INSL(f)
	function INSL0(seq)
		res = seq[1]
		for item in seq[2:end]
			res = f([res,item])
		end
		return res
	end
	return INSL0
end

# /////////////////////////////////////////////////////////////////
function CONS(Funs)
	return function(x)
		return [f(x) for f in Funs]
	end
end

# /////////////////////////////////////////////////////////////////
function IF(funs)
	function IF1(arg)
		f1, f2, f3 = funs
		return f1(arg) ? f2(arg) : f3(arg)
	end
	return IF1
end

# /////////////////////////////////////////////////////////////////
function LIFT(f)
	return function(funs)
		return COMP([f, CONS(funs)])
	end
end

# /////////////////////////////////////////////////////////////////
function RAISE(f)
	function RAISE0(args)
		return IF([ISSEQOF(ISFUN), LIFT(f), f])(args)
	end
	return RAISE0
end

# /////////////////////////////////////////////////////////////////
ISNUM(x) =  isa(x, Int) || isa(x, Float64) || isa(x, Complex)

ISFUN(x) = isa(x, Function)
ISREAL(x) = isa(x, Float64)
ISSEQ(x) = isa(x, Array)

# /////////////////////////////////////////////////////////////////
function ISSEQOF(type_checker)
	function ISSEQOF0(arg)
		if !isa(arg, Array)
			return false
		end
		for item in arg
			if !type_checker(item)
					return false
			end
		end
		return true
	end
	return ISSEQOF0
end

# /////////////////////////////////////////////////////////////////
function VECTSUM(vects)
	 return [sum(x) for x in zip(vects...)]
end

function VECTDIFF(vects)
	 return [l[1] - sum(l[2:end]) for l in zip(vects...)]
end

# /////////////////////////////////////////////////////////////////
function MEANPOINT(points)
	coeff=1.0/length(points)
	return [coeff*x for x in VECTSUM(points)]
end

# /////////////////////////////////////////////////////////////////
function SUM(args)
	if isa(args, Array) && ISPOL(args[1])
		return UNION(args)
	end
	if isa(args, Array) && ISNUM(args[1])
		return sum(args)
	end
	if isa(args, Array) && isa(args[1], Array)
		if isa(args[1][1], Array)
			return AA(VECTSUM)(zip(args...))
		else
			return VECTSUM(args)
		end
	end
	error("'+\' function has been applied to $args!")
end
ADD = SUM

# /////////////////////////////////////////////////////////////////
function DIFF(args)
	if isa(args, Array) && ISPOL(args[1])
		return DIFFERENCE(args)
	end
	if ISNUM(args)
		return -1 * args
	end
	if isa(args, Array) && ISNUM(args[1])
		return reduce((x,y) -> x - y, args)
	end
	if isa(args, Array) && isa(args[1], Array)
		if isa(args[1][1], Array)
			return AA(VECTDIFF)(zip(args...))
		else
			return VECTDIFF(args)
		end
	end
	error("\'-\' function has been applied to $args!")
end

# /////////////////////////////////////////////////////////////////
function PROD(args)
	if isa(args, Array) && ISPOL(args[1])
		return POWER(args)
	end
	if isa(args, Array) && ISSEQOF(ISNUM)(args)
		return reduce((x,y) -> x * y, args)
	end
	if isa(args, Array) && length(args) == 2 && ISSEQOF(ISNUM)(args[1]) && ISSEQOF(ISNUM)(args[2])
		return sum([a*b for (a,b) in zip(args[1],args[2])])
	end
	error("PROD function has been applied to $args!")
end

# /////////////////////////////////////////////////////////////////
SQR = RAISE(RAISE(PROD))([ID,ID])

function DIV(args)
	 return reduce((x,y) -> x/Float64(y), args)
end

REVERSE(List) = reverse(List)

# /////////////////////////////////////////////////////////////////
function TRANS(List)
	if isa(List, MatrixNd)
		return List'
	elseif isa(List, Tuple) || isa(List, Array)
		ret = zip(List...)
		cast_to_float = true
		for it in ret
			for jt in it
					if !(isa(jt, Int) || isa(jt, Float64))
						cast_to_float = false
					end
			end
		end
		if cast_to_float
			return [[Float64(jt) for jt in it] for it in ret]
		end
		return ret
	else
		error("invalid argument")
	end
end

# /////////////////////////////////////////////////////////////////
FIRST(List) = List[1]
LAST(List) = List[end]
TAIL(List) = List[2:end]
RTAIL(List) = List[1:end-1]

# /////////////////////////////////////////////////////////////////
function AR(args)
	v,value=args
	 return [v; value]
end

function AL(args)
	value,v=args
	 return [value;v]
end

LIST(x) = [x]

function Not(x)
	 return !x
end

NOT = AA(Not)

# /////////////////////////////////////////////////////////////////
function PROGRESSIVESUM(arg)
	ret,acc=[],0
	for value in arg
		acc+=value
		push!(ret, acc)
	end
	return ret
end

# /////////////////////////////////////////////////////////////////
function INTSTO(n)
	return collect(1:n)
end

function FROMTO(args)
	return collect(args[1]:args[2])
end

# /////////////////////////////////////////////////////////////////
SEL(n) = function(v)
	return v[n]
end

# /////////////////////////////////////////////////////////////////
S1 = SEL(1)
S2 = SEL(2)
S3 = SEL(3)
S4 = SEL(4)
S5 = SEL(5)
S6 = SEL(6)
S7 = SEL(7)
S8 = SEL(8)
S9 = SEL(9)
S10 = SEL(10)

# /////////////////////////////////////////////////////////////////
function N(n)
	return function(List)
		return [List for _ in 1:n]
	end
end

# /////////////////////////////////////////////////////////////////
function DIESIS(n)
	return function(List)
		return [List for _ in 1:n]
	end
end

# /////////////////////////////////////////////////////////////////
function NN(n)
	return function(v)
		return repeat(v,n)
	end
end

# /////////////////////////////////////////////////////////////////
function DOUBLE_DIESIS(n)
	return function(v)
		return repeat(v,n)
	end
end

# /////////////////////////////////////////////////////////////////
function AS(fun)
	return function(args)
		return COMP([CONS, AA(fun)])(args)
	end
end

# /////////////////////////////////////////////////////////////////
function AC(fun)
	return function(args)
		return COMP(AA(fun)(args))
	end
end

# /////////////////////////////////////////////////////////////////
function CHARSEQ(String)
	return [String[i] for i in 1:length(String)]
end

STRING(Charseq) = join(Charseq)

# /////////////////////////////////////////////////////////////////
function RANGE(Pair)
	if Pair[end] - Pair[1] >= 0
		return collect(Pair[1]:Pair[end])
	end
	return collect(Pair[1]:-1:Pair[end])
end

# /////////////////////////////////////////////////////////////////
SIGN(Number) = Number >= 0 ? 1 : -1

# /////////////////////////////////////////////////////////////////
function PRINT(AnyValue)
	 println(AnyValue)
	 return AnyValue
end

function PRINTPOL(PolValue)
	 println(PolValue)
	 flush(stdout)
	 return PolValue
end

# /////////////////////////////////////////////////////////////////
function TREE(f)
	function TREE1(v)
		N = length(v)
		if N == 1
			return v[1]
		end
		middle = trunc(Int,N/2)
		return f([TREE1(v[1:middle]);TREE1(v[middle+1:N])])
	end
	return TREE1
end

# /////////////////////////////////////////////////////////////////
function MERGE(f)
	function MERGE1(v)
		a, b = v
		if length(a) == 0
			return b
		end
		if length(b) == 0
			return a
		end
		res = f(a[1], b[1])
		if !res
			return [a[1];MERGE1([a[2:end],b])]
		else
			return [b[1];MERGE1([a,b[2:end]])]
		end
	end
	return MERGE1
end

# /////////////////////////////////////////////////////////////////
function CASE(ListPredFuns)
	function CASE_NO_CURRIED(ListPredFuns, x)
		for p in ListPredFuns
			if p[1](x)
					return p[2](x)
			end
		end
	end
	return function(arg)
		return CASE_NO_CURRIED(ListPredFuns, arg)
	end
end

# /////////////////////////////////////////////////////////////////
function PERMUTATIONS(SEQ)
	if length(SEQ) <= 1
		return [SEQ]
	end
	ret=[]
	for i in 1:length(SEQ)
		element =SEQ[i]
		rest  =PERMUTATIONS([ SEQ[1:i-1] ; SEQ[i+1:end] ])
		for r in rest
			push!(ret, [element; r])
		end
	end
	return ret
end

# /////////////////////////////////////////////////////////////////
function IDNT(N)
	return MatrixNd(N)
end

# /////////////////////////////////////////////////////////////////
function SPLIT_2PI(N)
	delta=2*PI/N
	return [i*delta for i in 0:N-1]
end

# /////////////////////////////////////////////////////////////////
function VECTPROD(args)
	a,b=args
	ax,ay,az=[Float64(it) for it in a]
	bx,by,bz=[Float64(it) for it in b]
	return [
		ay * bz - by * az,
		az * bx - bz * ax,
		ax * by - bx * ay
	]
end

# /////////////////////////////////////////////////////////////////
function VECTNORM(u)
	return sqrt(sum(x*x for x in u))
end

# /////////////////////////////////////////////////////////////////
function INNERPROD(args)
	return COMP([COMP([RAISE(SUM), AA(RAISE(PROD))]), TRANS])(args)
end

# /////////////////////////////////////////////////////////////////
function SCALARVECTPROD(args)
	s,l=args
	if !isa(l, Array)
		s,l=l,s
	end
	return [s*l[i] for i in 1:length(l)]
end

# /////////////////////////////////////////////////////////////////
function MIXEDPROD(args)
	A , B , C = args
	return INNERPROD([VECTPROD([A,B]),C])
end

# /////////////////////////////////////////////////////////////////
function UNITVECT(V)
	norm=VECTNORM(V)
	return [coord/norm for coord in V]
end

# /////////////////////////////////////////////////////////////////
function DIRPROJECT(E)
	function DIRPROJECT0(V)
		return SCALARVECTPROD([(INNERPROD([E,V])),E])
	end
	return DIRPROJECT0
end

# /////////////////////////////////////////////////////////////////
function ORTHOPROJECT(E)
	function ORTHOPROJECT0(V)
		return VECTDIFF([V,DIRPROJECT((E))(V)])
	end
	return ORTHOPROJECT0
end


# /////////////////////////////////////////////////////////////////
function FACT(N)
	if N>0
		return PROD(INTSTO(N)) 
	else 
		return 1
	end
end

ISREALVECT = ISSEQOF(ISREAL)
ISFUNVECT = ISSEQOF(ISFUN)
ISVECT = COMP([OR, CONS([ISREALVECT, ISFUNVECT])])
ISPOINT = ISVECT

# /////////////////////////////////////////////////////////////////
function CHOOSE(args)
	N , K = args
	return FACT(N)/float(FACT(K)*FACT(N-K))
end

# /////////////////////////////////////////////////////////////////
function TRACE(MATRIX)
	acc=0.0
	N=dim(MATRIX)
	for I in 1:N
		acc+=MATRIX[I,I]
	end
	return acc
end

# /////////////////////////////////////////////////////////////////
function MATHOM(T)
	N=dim(T)
	ret=MatrixNd(N+1)
	ret[2:N+1,2:N+1]=T[1:N,1:N]
	return ret
end

# /////////////////////////////////////////////////////////////////
function SCALARMATPROD(args)
	scalar,mat=float(args[1]),args[2]
	return [[scalar*coord for coord in row] for row in mat]
end

# /////////////////////////////////////////////////////////////////
function MATDOTPROD(args)
	return COMP([INNERPROD, AA(CAT)])(args)
end

# /////////////////////////////////////////////////////////////////
function ORTHO(matrix)
	return SCALARMATPROD([0.5,SUM([matrix,TRANS(matrix)])])
end

# /////////////////////////////////////////////////////////////////
function SKEW(matrix)
	return SCALARMATPROD([0.5,DIFF([matrix,TRANS(matrix)])])
end

# /////////////////////////////////////////////////////////////////
function SUBSEQ(I_J)
	function SUBSEQ0(SEQ)
		return SEQ[I_J[1]:I_J[2]]
	end
	return SUBSEQ0
end

# /////////////////////////////////////////////////////////////////
function VECT2MAT(v)
	n=length(v)
	return [[r != c ? 0 : v[r] for c in 1:n] for r in 1:n]
end

# /////////////////////////////////////////////////////////////////
function VECT2DTOANGLE(v)
	v=UNITVECT(v)
	return acos(v[1])*(v[2]>=0 ? 1 : -1)
end

# /////////////////////////////////////////////////////////////////
function CART(l)
	return vec(collect(Iterators.product(l...)))
end

# /////////////////////////////////////////////////////////////////
function POWERSET(l)
	return collect(powerset([1,2,3]))
end

# /////////////////////////////////////////////////////////////////
function ARC(args)
	degrees , cents = args
	return PI*(degrees+cents)/(100.0*180.0)
end

# /////////////////////////////////////////////////////////////////
function ISPOL(obj)
	return isa(obj, Hpc)
end

# /////////////////////////////////////////////////////////////////
function PRINTPOL(obj)
	println(obj)
	return obj
end

# /////////////////////////////////////////////////////////////////
function PRINT(obj)
	println(obj)
	return obj
end

# /////////////////////////////////////////////////////////////////
function VIEW(obj,title::String="")
	View(obj, title)
end

# /////////////////////////////////////////////////////////////////
function GRID(sequence)
	cursor=0
	points=[[0.0]]
	hulls=Vector{Vector{Int}}()
	N=1
	for value in sequence
		cursor+=abs(value)
		push!(points,[cursor])
		N+=1
		if value>=0
			push!(hulls,[N-1,N])
		end
	end
	return  MkPol(points, hulls)
end
QUOTE = GRID

# Q = COMP([QUOTE, IF([ISSEQ, ID, CONS([ID])])])

# /////////////////////////////////////////////////////////////////
function INTERVALS(A)
	A=Float64(A)
	function INTERVALS0(N::Int64)
		return QUOTE([A/N for I in 1:N])
	end
	return INTERVALS0
end

# /////////////////////////////////////////////////////////////////
function CUBOID(vs)
	return Scale(Cube(length(vs)),[Float64(it) for it in vs])
end

function CUBE(size)
	return CUBOID([Float64(size), Float64(size), Float64(size)])
end


# /////////////////////////////////////////////////////////////////
function HEXAHEDRON()
	return Cube(3, -1.0/sqrt(3.0), +1.0/sqrt(3.0))
end

# /////////////////////////////////////////////////////////////////
function SIMPLEX(dim)
	return Simplex(dim)
end

RN(pol) = dim(pol)

DIM(pol) = dim(pol)

# /////////////////////////////////////////////////////////////////
function ISPOLDIM(dims)
	function ISPOLDIM1(pol)
		d = dims[1]
		n = dims[2]
		return (d == DIM(pol)) && (n == RN(pol))
	end
	return ISPOLDIM1
end

# /////////////////////////////////////////////////////////////////
function MKPOL(points::Vector{Vector{Float64}}, hulls::Vector{Vector{Int}},__pols=Nothing)
	return MkPol(points, hulls)
end

# deprecated
function MKPOL(args)
	return MKPOL(args[1], args[2])
end

# /////////////////////////////////////////////////////////////////
MK = COMP([MKPOL, CONS([LIST, K([[1]]), K([[1]])])])

# /////////////////////////////////////////////////////////////////
function CONVEXHULL(points)
	return MKPOL(points, [collect(1:length(points))], [[1]])
end

# /////////////////////////////////////////////////////////////////
function UKPOL(pol)
	points, hulls = UkPol(pol)
	ret=Vector{Any}()
	push!(ret,points)
	push!(ret,hulls)
	push!(ret,[[1]])
	return ret
end

UK = COMP([COMP([S1, S1]), UKPOL])

# /////////////////////////////////////////////////////////////////
OPTIMIZE(pol) = pol

# /////////////////////////////////////////////////////////////////
function TRANSLATE(axis)
	function TRANSLATE1(values)
		function TRANSLATE2(pol)
			axis = ISNUM(axis) ? [axis] : axis
			values = ISNUM(values) ? [values] : values
			vt = [0.0 for I in 1:max(axis...)]
			for (a, t) in zip(axis, values)
					vt[a] = t
			end
			return Translate(pol,vt)
		end
		return TRANSLATE2
	end
	return TRANSLATE1
end
T = TRANSLATE

# /////////////////////////////////////////////////////////////////
function SCALE(axis)
	function SCALE1(values)
		function SCALE2(pol)
			axis = ISNUM(axis) ? [axis] : axis
			values = ISNUM(values) ? [values] : values
			vs = [1.0 for I in 1:max(axis...)] 
			for (a, t) in zip(axis, values)
				vs[a] = t
			end
			return Scale(pol,vs)
		end
		return SCALE2
	end
	return SCALE1
end
S = SCALE

# /////////////////////////////////////////////////////////////////
function ROTATE(plane_indexes)
	function ROTATE1(angle)
		function ROTATE2(pol)
			return Rotate(pol, plane_indexes[1], plane_indexes[2], angle)
		end
		return ROTATE2
	end
	return ROTATE1
end
R = ROTATE

# /////////////////////////////////////////////////////////////////
function SHEARING(i)
	function SHEARING1(shearing_vector_list)
		function SHEARING2(pol)
			error("shearing not implemented!")
		end
		return SHEARING2
	end
	return SHEARING1
end
H = SHEARING

# /////////////////////////////////////////////////////////////////
function MAT(matrix)
	function MAT0(pol)
		return Transform(pol,MatrixNd(matrix))
	end
	return MAT0
end

# /////////////////////////////////////////////////////////////////
function EMBED(up_dim)
	function EMBED1(pol)
		return Hpc(MatrixNd(dim(pol) + up_dim+1), [pol])
	end
	return EMBED1
end

# /////////////////////////////////////////////////////////////////
function STRUCT(seq, nrec=0)
	if !isa(seq, Vector)
		error("STRUCT must be applied to a list")
	end
	if isempty(seq)
		error("STRUCT must be applied to a not empty list")
	end
	if nrec == 0
		seq = copy(seq)
	end

	# collect all geometry wihout transformations
	pols = Vector{Hpc}()
	while !isempty(seq) && ISPOL(seq[1])
		push!(pols, seq[1])
		seq = seq[2:end]
	end

	# callect all transformation on the right
	transformations = []
	while !isempty(seq) && ISFUN(seq[1])
		push!(transformations, seq[1])
		seq = seq[2:end]
	end

	if !isempty(seq) && !ISPOL(seq[1]) && !ISFUN(seq[1])
		error("STRUCT arguments not valid, not all elements are pols or transformations")
	end
	
	# recursive on the right, apply transformations
	if !isempty(seq)
		@assert(ISPOL(seq[1]))
		child = STRUCT(seq, nrec+1)
		@assert(ISPOL(child))
		if !isempty(transformations)
			child = COMP(transformations)(child)
		end
		push!(pols, child)
	end

	if isempty(pols)
		error("Cannot find geometry in STRUCT, found only transformations")
	end
	return Struct(pols)
end

# /////////////////////////////////////////////////////////////////////////////
struct BspFace
	vertices::Vector{Vector{Float64}}
	plane::Vector{Float64}
	
	function BspFace(vertices, plane=nothing)

		# automatic computation of plane
		if plane==nothing
			dim = length(vertices[1])
			if dim == 2
				@assert length(vertices) == 2
				x1, y1 = vertices[1][1], vertices[1][2]
				x2, y2 = vertices[2][1], vertices[2][2]
				normal = [(y2-y1), -(x2-x1)]
			elseif dim == 3
				normal = ComputeTriangleNormal(vertices[1], vertices[2], vertices[3])
			else
				error("todo")
			end

			w = sum([normal[i]*vertices[1][i] for i in 1:dim])
			plane = [normal; -w]
		 end
		 new(vertices, plane)
	end
	
end


function splitEdge(plane::Vector{Float64}, vi::Vector{Float64}, vj::Vector{Float64})
	dim = length(vi)
	valuev1 = abs(plane[end] + sum([plane[i]*vi[i] for i in 1:dim]))
	valuev2 = abs(plane[end] + sum([plane[i]*vj[i] for i in 1:dim]))
	alpha = 1.0 / (valuev1 + valuev2)
	beta = valuev1 * alpha
	alpha = alpha * valuev2
	return [alpha*vi[i] + beta*vj[i] for i in 1:dim]
end

function splitFace(self::BspFace, plane::Vector{Float64}, EPSILON::Float64=1e-5)
	dim = length(plane) - 1
	COPLANAR, ABOVE, BELOW, SPANNING = 0, 1, 2, 3
	ptype, types = COPLANAR, []
	for v in self.vertices
		 @assert length(v) == dim
		 t = plane[end] + sum([plane[i]*v[i] for i in 1:dim])
		 if t < -EPSILON
			  type = BELOW
		 elseif t > EPSILON
			  type = ABOVE
		 else
			  type = COPLANAR
		 end
		 ptype |= type
		 push!(types, type)
	end
	if ptype == BELOW
		 return [self, nothing, nothing, nothing]
	end
	if ptype == ABOVE
		 return [nothing, nothing, nothing, self]
	end
	if ptype == COPLANAR
		 if sum([plane[i]*self.plane[i] for i in 1:dim+1]) > 0
			  return [nothing, nothing, self, nothing]
		 else
			  return [nothing, self, nothing, nothing]
		 end
	end
	@assert ptype == SPANNING
	b, f = [], []
	if dim == 2
		 @assert length(self.vertices) == 2
		 ti, tj = types[1], types[2]
		 vi, vj = self.vertices[1], self.vertices[2]
		 if ti != BELOW
			  push!(f, vi)
		 end
		 if ti != ABOVE
			  push!(b, vi)
		 end
		 if tj != BELOW
			  push!(f, vj)
		 end
		 if tj != ABOVE
			  push!(b, vj)
		 end
		 if (ti | tj) == SPANNING
			  v = splitEdge(plane, vi, vj)
			  push!(b, v)
			  push!(f, v)
		 end
	elseif dim == 3
		 for i in 1:length(self.vertices)
			  j = (i + 1) % length(self.vertices)
			  ti, tj = types[i], types[j]
			  vi, vj = self.vertices[i], self.vertices[j]
			  if ti != BELOW
					push!(f, vi)
			  end
			  if ti != ABOVE
					push!(b, vi)
			  end
			  if (ti | tj) == SPANNING
					v = splitEdge(plane, vi, vj)
					push!(b, v)
					push!(f, v)
			  end
		 end
	else
		 error("not supported")
	end
	@assert length(b) >= dim && length(f) >= dim
	return [BspFace(b, self.plane), nothing, nothing, BspFace(f, self.plane)]
end

# //////////////////////////////////////////////////////////////////////////////////////
struct Bsp
	 
	 plane::Vector{Float64}
	 faces::Vector{BspFace}
	 below::Bsp
	 above::Bsp

	 function Bsp()
		new()
	 end
	 
end

function getFaces(self::Bsp)
	if self==nothing
		return []
	end
	return [self.faces;getFaces(self.below);getFaces(self.above)]
end

function insertFaces(self::Bsp, faces::Vector{BspFace})
	if length(faces)==0
		 return self
	end

	if self.plane == nothing
		 @assert self.below === nothing && self.above === nothing
		 self.plane = faces[1].plane
	end

	below, above = [], []
	for face in faces
		 b, cb, ca, a = splitFace(face, self.plane)
		 if b != nothing
			  push!(below, b)
		 end
		 if cb != nothing
			  push!(self.faces, cb)
		 end
		 if ca != nothing
			  push!(self.faces, ca)
		 end
		 if a != nothing
			  push!(above, a)
		 end
	end

	if !isempty(above)
		 if self.above === nothing
			  self.above = Bsp()
		 end
		 insertFaces(self.above, above)
	end

	if !isempty(below)
		 if self.below === nothing
			  self.below = Bsp()
		 end
		 insertFaces(self.below, below)
	end
end

function clipFaces(self::Bsp, faces)
	if self.plane === nothing
		 return faces
	end
	below, above = [], []
	for p in faces
		 b, cb, ca, a = splitFace(p, self.plane)
		 if b !== nothing
			  push!(below, b)
		 end
		 if cb !== nothing
			  push!(below, cb)
		 end
		 if ca !== nothing
			  push!(above, ca)
		 end
		 if a !== nothing
			  push!(above, a)
		 end
	end
	below = self.below !== nothing ? clipFaces(self.below, below) : []
	above = self.above !== nothing ? clipFaces(self.above, above) : []
	return above + below
end

function clipTo(self::Bsp, other)
	self.faces = clipFaces(other, self.faces)
	if self.below !== nothing
		 clipTo(self.below, other)
	end
	if self.above !== nothing
		 clipTo(self.above, other)
	end
end

function Complement(bsp::Bsp)
	if bsp === nothing
		 return nothing
	end
	ret = Bsp()
	if bsp.plane !== nothing
		 ret.plane = [-1*it for it in bsp.plane]
	end
	for p in bsp.faces
		 new_p = BspFace(reverse(p.vertices), [-1*c for c in p.plane])
		 push!(ret.faces, new_p)
	end
	ret.below = Complement(bsp.above)
	ret.above = Complement(bsp.below)
	return ret
end

function Union(a::Bsp, b::Bsp)
	clipTo(a, b)
	clipTo(b, a)
	b = Bsp.Complement(b)
	clipTo(b, a)
	b = Bsp.Complement(b)
	insertFaces(a, getFaces(b))
	return a
end

function Intersection(a::Bsp, b::Bsp)
	return Bsp.Complement(Bsp.Union(Bsp.Complement(a), Bsp.Complement(b)))
end

function Difference(a::Bsp, b::Bsp)
	return Bsp.Complement(Bsp.Union(Bsp.Complement(a), b))
end

function Xor(a::Bsp, b::Bsp)
	return Bsp.Union(Bsp.Intersection(a, Bsp.Complement(b)), Bsp.Intersection(Bsp.Complement(a), b))
end

function fromHpc(hpc::Hpc)
	ret = Bsp()
	faces = []
	for (T, properties, obj) in toList(ToBoundaryForm(hpc))
		for hull in obj.hulls
			points=[transformPoint(T,obj.points[I]) for I in hull]
			push!(faces, BspFace(points))
		end
	end
	insertFaces(ret, faces)
	return ret
end

function toHpc(self::Bsp)
	batches, faces = [], getFaces(self)
	dim = self.plane !== nothing ? length(self.plane) - 1 : 0
	if dim == 0
		 return Hpc()
	end
	@assert dim == 1 || dim == 2 || dim == 3
	points, hulls = [], []
	for face in faces
		 if dim == 1
			  @assert length(face.vertices) == 1
			  N = length(points)
			  points += face.vertices
			  push!(hulls, collect(N:N+length(face.vertices)-1))
		 elseif dim == 2
			  @assert length(face.vertices) == 2
			  N = length(points)
			  points += face.vertices
			  push!(hulls, collect(N:N+length(face.vertices)-1))
		 elseif dim == 3
			  @assert length(face.vertices) >= 3
			  for I in 2:length(face.vertices)-1
					N = length(points)
					points += [face.vertices[1], face.vertices[I], face.vertices[I+1]]
					push!(hulls, collect(N:N+2))
			  end
		 else
			  error("not supported")
		 end
	end
	return MkPol(points, hulls)
end

function UNION(objs::Vector{Hpc})
	objs = [fromHpc(obj) for obj in objs]
	res = objs[1]
	for I in 2:length(objs)
		res = Union(res, objs[I])
	end
	return toHpc(res)
end

function INTERSECTION(objs::Vector{Hpc})
	objs = [fromHpc(obj) for obj in objs]
	res = objs[1]
	for I in 2:length(objs)
		res = Intersection(res, objs[I])
	end
	return toHpc(res)
end

function DIFFERENCE(objs::Vector{Hpc})
	objs = [fromHpc(obj) for obj in objs]
	res = objs[1]
	for I in 2:length(objs)
		res = Difference(res, objs[I])
	end
	return toHpc(res)
end

function XOR(objs::Vector{Hpc})
	objs = [fromHpc(obj) for obj in objs]
	res = objs[1]
	for I in 2:length(objs)
		res = Xor(res, objs[I])
	end
	return toHpc(res)
end

# ///////////////////////////////////////////////////////////
function JOIN(pol_list)
	 if ISPOL(pol_list)
		  pol_list = [pol_list]
	 end
	 return Join(pol_list)
end

# ///////////////////////////////////////////////////////////
function POWER(objs_list)
	 if !isa(objs_list, Vector) || length(objs_list) != 2
		  error("POWER can only be applied to a list of 2 arguments")
	 end
	 if ISNUM(objs_list[1]) && ISNUM(objs_list[2])
		  return objs_list[1] ^ objs_list[2]
	 end
	 return Power(objs_list[1], objs_list[2])
end

# ///////////////////////////////////////////////////////////
function SKELETON(ord)
	error("not implemented")
	SKEL_0 = SKELETON(0)
	SKEL_1 = SKELETON(1)
	SKEL_2 = SKELETON(2)
	SKEL_3 = SKELETON(3)
	SKEL_4 = SKELETON(4)
	SKEL_5 = SKELETON(5)
	SKEL_6 = SKELETON(6)
	SKEL_7 = SKELETON(7)
	SKEL_8 = SKELETON(8)
	SKEL_9 = SKELETON(9)
end

# ///////////////////////////////////////////////////////////
function SIZE(sel)
	function SIZE1(pol)
		S = size(box(pol))
		return isa(sel, Vector) ? [S[i] for i in sel] : S[sel]
	end
	return SIZE1
end

# ///////////////////////////////////////////////////////////
function MIN(sel)
	function MIN1(pol)
		b = box(pol)
		return isa(sel, Vector) ? [b.p1[i] for i in sel] : b.p1[sel]
	end
	return MIN1
end

# ///////////////////////////////////////////////////////////
function MAX(sel)
	function MAX1(pol)
		b = box(pol)
		return isa(sel, Vector) ? [b.p2[i] for i in sel] : b.p2[sel]
	end
	return MAX1
end

# ///////////////////////////////////////////////////////////
function MED(sel)
	function MED1(pol)
		c = center(box(pol))
		return isa(sel, Vector) ? [c[i] for i in sel] : c[sel]
	end
	return MED1
end

# ///////////////////////////////////////////////////////////
function ALIGN(args)
	function ALIGN0(args, pols)
		pol1, pol2 = pols
		box1, box2 = box(pol1), box(pol2)
		if isa(args, Vector) && !isempty(args) && ISNUM(args[1])
			args = [args]
		end
		max_index = max([index for (index, pos1, po2) in args]...)
		vt = zeros(max_index)
		for (index, pos1, pos2) in args
			p1 = ifelse(pos1 == MIN, box1.p1, ifelse(pos1 == MAX, box1.p2, center(box1)))
			p1 = ifelse(index <= length(p1), p1[index], 0.0)
			p2 = ifelse(pos2 == MIN, box2.p1, ifelse(pos2 == MAX, box2.p2, center(box2)))
			p2 = ifelse(index <= length(p2), p2[index], 0.0)
			vt[index] -= (p2 - p1)
		end
		return Struct([pol1, Translate(pol2,vt)])
	end
	return pols -> ALIGN0(args, pols)
end
TOP = ALIGN([[3, MAX, MIN], [1, MED, MED], [2, MED, MED]])
BOTTOM = ALIGN([[3, MIN, MAX], [1, MED, MED], [2, MED, MED]])
LEFT = ALIGN([[1, MIN, MAX], [3, MIN, MIN]])
RIGHT = ALIGN([[1, MAX, MIN], [3, MIN, MIN]])
UP = ALIGN([[2, MAX, MIN], [3, MIN, MIN]])
DOWN = ALIGN([[2, MIN, MAX], [3, MIN, MIN]])

# ///////////////////////////////////////////////////////////
function BOX(sel)
	function BOX0(pol)
		if !isa(sel, Vector)
			sel = [sel]
		end
		dim = length(sel)
		b = box(pol)
		vt = [  b.p1[i] for i in sel]
		vs = [size(b)[i] for i in sel]
		return Translate(Scale(Cube(dim),vs),vt)
	end
	return BOX0
end

# ///////////////////////////////////////////////////////////
function MAP(fn)
	function MAP0(pol::Hpc)
		if isa(fn, Tuple) || isa(fn, Vector)
			return MapFn(pol, p -> [f(p) for f in fn])
		else
			return MapFn(pol,fn)
		end
	end
	return MAP0
end


# /////////////////////////////////////////////////////////////////
function CIRCLE_POINTS(R, N)
	return [[R*cos(i*2*pi/N), R*sin(i*2*pi/N)] for i in 0:N-1]
end

# ///////////////////////////////////////////////////////////
function CIRCUMFERENCE(R)
	function CIRCUMFERENCE1(N)
		domain=INTERVALS(2*pi)(N)
		fn=p -> [R*cos(p[1]), R*sin(p[1])]
		return MAP(fn)(domain)
	end
	return CIRCUMFERENCE1
end

# ///////////////////////////////////////////////////////////
function NGON(N)
	return CIRCUMFERENCE(1)(N)
end

# ///////////////////////////////////////////////////////////
function RING(radius::Vector{Float64})
	R1, R2 = radius
	function RING0(subds)
		N, M = subds
		domain = Translate(POWER([INTERVALS(2*pi)(N), INTERVALS(R2-R1)(M)]),[0.0, R1])
		fun = p -> [p[2]*cos(p[1]), p[2]*sin(p[1])]
		return MAP(fun)(domain)
	end
	return RING0
end

# ///////////////////////////////////////////////////////////
function TUBE(args::Vector{Float64})
	r1, r2, height = args
	function TUBE0(N)
		return Power(RING([r1, r2])([N, 1]), QUOTE([height]))
	end
	return TUBE0
end

# ///////////////////////////////////////////////////////////
function CIRCLE(R::Float64)
	function CIRCLE0(subs)
		N, M = subs
		domain = POWER([INTERVALS(2*pi)(N), INTERVALS(R)(M)])
		fun = p -> [p[2]*cos(p[1]), p[2]*sin(p[1])]
		return MAP(fun)(domain)
	end
	return CIRCLE0
end

# ///////////////////////////////////////////////////////////
function MY_CYLINDER(args::Vector{Float64})
	R, H = args
	function MY_CYLINDER0(N)
		points = CIRCLE_POINTS(R, N)
		circle = MkPol(points, [collect(1:N)])
		return Power(circle, MkPol([[0], [H]], [[1, 2]]))
	end
	return MY_CYLINDER0
end
CYLINDER = MY_CYLINDER

# /////////////////////////////////////////////////////////////
function SPHERE(radius::Float64)
	function SPHERE0(subds)
		N, M = subds
		domain = Translate(Power(INTERVALS(pi)(N), INTERVALS(2*pi)(M)),[-pi/2, 0])
		fx = p -> radius * cos(p[1]) * sin(p[2])
		fy = p -> radius * cos(p[1]) * cos(p[2])
		fz = p -> radius * sin(p[1])
		return MAP([fx, fy, fz])(domain)
	end
	return SPHERE0
end

# /////////////////////////////////////////////////////////////
function TORUS(radius)
	r1, r2 = radius
	function TORUS0(subds)
		N, M = subds
		a = 0.5*(r2-r1)
		c = 0.5*(r1+r2)
		domain = Power(INTERVALS(2*pi)(N), INTERVALS(2*pi)(M))
		fx = p -> (c+a*cos(p[2])) * cos(p[1])
		fy = p -> (c+a*cos(p[2])) * sin(p[1])
		fz = p -> a*sin(p[2])
		return MAP([fx, fy, fz])(domain)
	end
	return TORUS0
end

# /////////////////////////////////////////////////////////////
function CONE(args)
	radius, height = args
	function CONE0(N)
		basis = CIRCLE(radius)([N, 1])
		apex = T(3)(height)(SIMPLEX(0))
		return JOIN([basis, apex])
	end
	return CONE0
end

# /////////////////////////////////////////////////////////////
function TRUNCONE(args)
	R1, R2, H = args
	function TRUNCONE0(N)
		domain = Power(QUOTE([2*pi/N for i in 1:N]), QUOTE([1]))
		fn = p -> [
			(R1+p[2]*(R2-R1))*cos(p[1]),
			(R1+p[2]*(R2-R1))*sin(p[1]),
			H*p[2]
		]
		return MAP(fn)(domain)
	end
	return TRUNCONE0
end

# /////////////////////////////////////////////////////////////
function DODECAHEDRON()
	a = 1.0/(sqrt(3.0))
	g = 0.5*(sqrt(5.0)-1)

	top = MKPOL(
		[[1-g, 1, 0-g], [1+g, 1, 0-g]], 
		[[1, 2]], 
		[[1]]
	)
	basis = EMBED(1)(CUBOID([2.0, 2.0]))
	roof = T([1, 2, 3])([-1, -1, -1])(JOIN([basis, top]))
	roofpair = STRUCT([roof, R([2, 3])(pi), roof])
	return S([1, 2, 3])([a, a, a])(STRUCT([
		Cube(3, -1, +1),
		roofpair,
		R([1, 3])(pi/2), R([1, 2])(pi/2),
		roofpair,
		R([1, 2])(pi/2), R([2, 3])(pi/2),
		roofpair
	]))
end

# /////////////////////////////////////////////////////////////
function ICOSAHEDRON()
	g = 0.5*(sqrt(5)-1)
	b = 2.0/(sqrt(5*sqrt(5)))
	rectx = T([1, 2])([-g, -1])(CUBOID([2*g, 2]))
	recty = R([1, 3])(pi/2)(R([1, 2])(pi/2)(rectx))
	rectz = R([2, 3])(pi/2)(R([1, 2])(pi/2)(rectx))
	return S([1, 2, 3])([b, b, b])(JOIN([rectx, recty, rectz]))
end

function TETRAHEDRON()
	 return JOIN([T(3)(-1.0/3.0)(NGON(3)), MK([0, 0, 1])])
end

function POLYPOINT(points)
	return MkPol(points, [[i] for i in 1:length(points)])
end

function POLYLINE(points)
	return MkPol(points, [[i, i+1] for i in 1:length(points)-1])
end

function TRIANGLESTRIPE(points)
	cells = [i%2==0 ? [i, i+1, i+2] : [i+1, i, i+2] for i in 1:length(points)-2]
	return MkPol(points, cells)
end

function TRIANGLEFAN(points)
	cells = [[1, i-1, i] for i in 2:length(points)]
	return MkPol(points, cells)
end

function MIRROR(D)
	function MIRROR0(pol)
		return STRUCT([S(D)(-1)(pol), pol])
	end
	return MIRROR0
end

# /////////////////////////////////////////////////////////////
function POLYMARKER(type, MARKERSIZE=0.1)
	A, B = MARKERSIZE, -MARKERSIZE
	marker0 = MkPol([[A], [0], [0], [A], [B], [0], [0], [B]], [[0, 1], [1, 2], [2, 3], [3, 0]])
	marker1 = MkPol([[A], [A], [B], [A], [B], [B], [A], [B]], [[0, 2], [1, 3]])
	marker2 = MkPol([[A], [A], [B], [A], [B], [B], [A], [B]], [[0, 1], [1, 2], [2, 3], [3, 0]])
	marker3 = STRUCT([marker0, marker1])
	marker4 = STRUCT([marker0, marker2])
	marker5 = STRUCT([marker1, marker2])
	marker = [marker0, marker1, marker2, marker3, marker4, marker5][mod(type, 6)+1]
	function POLYMARKER_POINTS(points)
		dim = length(points[1])
		axis = collect(1:dim)
		return Struct([T(axis, point)(marker) for point in points])
	end
	return POLYMARKER_POINTS
end

# /////////////////////////////////////////////////////////////
function BEZIER(U)
	function BEZIER0(controldata_fn)
		N = length(controldata_fn)-1
		function map_fn(point)
			t = U(point)
			controldata = [ isa(fun, Function) ? fun(point) : fun for fun in controldata_fn]
			ret = [0.0 for i in 1:length(controldata[1])]
			for I in 0:N
					weight = CHOOSE([N, I])*((1-t)^(N-I))*(t^I)
					for K in 1:length(ret)
						ret[K] += weight*controldata[I+1][K]
					end
			end
			return ret
		end
		return map_fn
	end
	return BEZIER0
end

# /////////////////////////////////////////////////////////////
function BEZIERCURVE(controlpoints)
	return BEZIER(S1)(controlpoints)
end

# /////////////////////////////////////////////////////////////
function COONSPATCH(args)
	su0_fn, su1_fn, s0v_fn, s1v_fn = args
	function map_fn(point)
		u, v = point
		su0 = isa(su0_fn, Function) ? su0_fn(point) : su0_fn
		su1 = isa(su1_fn, Function) ? su1_fn(point) : su1_fn
		s0v = isa(s0v_fn, Function) ? s0v_fn(point) : s0v_fn
		s1v = isa(s1v_fn, Function) ? s1v_fn(point) : s1v_fn
		ret = [0.0 for i in 1:length(su0)]
		for K in 1:length(ret)
			ret[K] = (1-u)*s0v[K] + u*s1v[K] + (1-v)*su0[K] + v*su1[K] + (1-u)*(1-v)*s0v[K] + (1-u)*v*s0v[K] + u*(1-v)*s1v[K] + u*v*s1v[K]
		end
		return ret
	end
	return map_fn
end

# /////////////////////////////////////////////////////////////
function RULEDSURFACE(args)
	alpha_fn, beta_fn = args
	function map_fn(point)
		u, v = point
		alpha, beta = alpha_fn(point), beta_fn(point)
		ret = [0.0 for i in 1:length(alpha)]
		for K in 1:length(ret)
			ret[K] = alpha[K] + v*beta[K]
		end
		return ret
	end
	return map_fn
end

# /////////////////////////////////////////////////////////////
function PROFILEPRODSURFACE(args)
	profile_fn, section_fn = args
	function map_fun(point)
		u, v = point
		profile, section = profile_fn(point), section_fn(point)
		ret = [profile[1]*section[1], profile[1]*section[2], profile[3]]
		return ret
	end
	return map_fun
end

# /////////////////////////////////////////////////////////////
function ROTATIONALSURFACE(args)
profile = args
function map_fn(point)
	u, v = point
	f, h, g = profile(point)
	ret = [f*cos(v), f*sin(v), g]
	return ret
end
return map_fn
end

# /////////////////////////////////////////////////////////////
function CYLINDRICALSURFACE(args)
	alpha_fun = args[1]
	beta_fun = map(K, args[2])
	return RULEDSURFACE([alpha_fun, beta_fun])
end

# /////////////////////////////////////////////////////////////
function CONICALSURFACE(args)
	apex = args[1]
	alpha_fn = point -> apex
	beta_fn = point -> [args[2](point)[i]-apex[i] for i in 1:length(apex)]
	return RULEDSURFACE([alpha_fn, beta_fn])
end

# /////////////////////////////////////////////////////////////
function CUBICHERMITE(U)
	function CUBICHERMITE0(args)
		p1_fn, p2_fn, s1_fn, s2_fn = args
		function map_fn(point)
			u = U(point)
			u2 = u*u
			u3 = u2*u
			p1, p2, s1, s2 = [ isa(f, Function) ? f(point) : f for f in [p1_fn, p2_fn, s1_fn, s2_fn]]
			ret = [0.0 for i in 1:length(p1)]
			for i in 1:length(ret)
					ret[i] = (2*u3-3*u2+1)*p1[i] + (-2*u3+3*u2)*p2[i] + (u3-2*u2+u)*s1[i] + (u3-u2)*s2[i]
			end
			return ret
		end
		return map_fn
	end
	return CUBICHERMITE0
end

# /////////////////////////////////////////////////////////////
function HERMITE(args)
	P1, P2, T1, T2 = args
	return CUBICHERMITE(S1)([P1, P2, T1, T2])
end

# /////////////////////////////////////////////////////////////
function Q(H)
	return MkPol([[0], [H]], [[1, 2]])
end

# /////////////////////////////////////////////////////////////
function EXTRUDE(args)
	__N, Pol, H = args
	return Power(Pol, Q(H))
end

# /////////////////////////////////////////////////////////////
function MULTEXTRUDE(P)
	 function MULTEXTRUDE0(H)
		  return Power(P, Q(H))
	 end
	 return MULTEXTRUDE0
end

# /////////////////////////////////////////////////////////////
function PROJECT(M)
	function PROJECT0(POL)
		vertices, cells, pols = UKPOL(POL)
		vertices = [vert[1:end-M] for vert in vertices]
		return MKPOL(vertices, cells, pols)
	end
	return PROJECT0
end

# /////////////////////////////////////////////////////////////
function SPLITCELLS(scene)
	vertices, cells, pols = UKPOL(scene)
	ret = []
	for c in cells
		push!(ret, MKPOL(vertices, [c], [[1]]))
	end
	return ret
end

# /////////////////////////////////////////////////////////////
function EXTRACT_WIRES(scene)
	return SPLITCELLS(SKEL_1(scene))
end

SPLITPOLS = SPLITCELLS


# /////////////////////////////////////////////////////////////
function PERMUTAHEDRON(d)
	vertices = ToFloat64(PERMUTATIONS(collect(1:d+1)))
	center = MEANPOINT(vertices)
	cells = [collect(1:length(vertices))]
	object = MKPOL(vertices, cells, [[1]])
	object = Translate(object, [-coord for coord in center])
	for i in 1:d
		object = R([i, d+1])(pi/4)(object)
	end
	object = PROJECT(1)(object)
	return object
end

# /////////////////////////////////////////////////////////////
function STAR(N)
	function CIRCLEPOINTS(STARTANGLE)
		function CIRCLEPOINTS1(R)
			function CIRCLEPOINTS0(N)
					return [COMP([CONS([RAISE(PROD)([K(R), COS]), RAISE(PROD)([K(R), SIN])]), RAISE(SUM)([ID, K(STARTANGLE)])]) for STARTANGLE in COMP([COMP([AA(RAISE(PROD)), TRANS]), CONS([K(collect(1:N)), DIESIS(N)])])((2*pi/N))]
			end
			return CIRCLEPOINTS0
		end
		return CIRCLEPOINTS1
	end
	return COMP([COMP([TRIANGLEFAN, CAT]), TRANS])([CIRCLEPOINTS(0)(1)(N), CIRCLEPOINTS((pi/N))(2.5)(N)])
end


# /////////////////////////////////////////////////////////////
function SCHLEGEL2D(D)
	function map_fn(point)
		return [D*point[1]/point[3], D*point[2]/point[3]]
	end
	return MAP(map_fn)
end

# /////////////////////////////////////////////////////////////
function SCHLEGEL3D(D)
	function map_fn(point)
		return [D*point[1]/point[4], D*point[2]/point[4], D*point[3]/point[4]]
	end
	return MAP(map_fn)
end

# /////////////////////////////////////////////////////////////
function FINITECONE(pol)
	 point = [0 for i in 1:RN(pol)]
	 return JOIN([pol, MK(point)])
end

# /////////////////////////////////////////////////////////////
function PRISM(HEIGHT)
	function PRISM0(BASIS)
		return Power(BASIS, QUOTE([HEIGHT]))
	end
	return PRISM0
end

# /////////////////////////////////////////////////////////////
function CROSSPOLYTOPE(D)
	points = []
	for i in 1:D
		point_pos = [0 for x in 1:D]; point_pos[i] = +1
		point_neg = [0 for x in 1:D]; point_neg[i] = -1
		push!(points, point_pos, point_neg)
	end
	cells = [collect(1:D*2)]
	pols = [[1]]
	return MKPOL(points, cells, pols)
end

# /////////////////////////////////////////////////////////////
function OCTAHEDRON()
	return CROSSPOLYTOPE(2)
end

# /////////////////////////////////////////////////////////////
function ROTN(args)
	alpha, N = args
	N = UNITVECT(N)
	QX = UNITVECT(VECTPROD([[0, 0, 1], N]))
	QZ = UNITVECT(N)
	QY = VECTPROD([QZ, QX])
	Q = MATHOM([QX, QY, QZ])
	ISUP = COMP([AND, CONS([COMP([C(EQ)(0), S1]), COMP([C(EQ)(0), S2]), COMP([COMP([NOT, C(EQ)(0)]), S3])])])
	if N[1]==0 && N[2]==0
		return R([1, 2])(alpha)
	else
		return COMP([MAT(TRANS(Q)), R([1, 2])(alpha), MAT(Q)])
	end
end

# /////////////////////////////////////////////////////////////
function  MKVERSORK()
	return TOP([CYLINDER([1.0/100.0, 7.0/8.0])(6), CONE([1.0/16.0, 1.0/8])(8)])
end

function MKVECTOR(P1)
	function MKVECTOR0(P2)
		TR = T([1, 2, 3])(P1)
		U = VECTDIFF([P2, P1])
		ALPHA = acos(INNERPROD([[0, 0, 1], UNITVECT(U)]))
		B = VECTNORM(U)
		SC = S([1, 2, 3])([B, B, B])
		N = VECTPROD([[0, 0, 1], U])
		ROT = ROTN([ALPHA, N])
		return COMP([COMP([TR, ROT]), SC])(MKVERSORK())
	end
	return MKVECTOR0
end



# /////////////////////////////////////////////////////////////
function CUBICUBSPLINE(domain)
	function CUBICUBSPLINE0(args)
		q1_fn, q2_fn, q3_fn, q4_fn = args
		function map_fn(point)
			u = S1(point)
			u2 = u*u
			u3 = u2*u
			q1, q2, q3, q4 = [ isa(f, Function)  ? f(point) : f for f in [q1_fn, q2_fn, q3_fn, q4_fn]]
			ret = [0 for x in 1:length(q1)]
			for i in 1:length(ret)
					ret[i] = (1.0/6.0)*((-u3+3*u2-3*u+1)*q1[i] + (3*u3-6*u2+4)*q2[i] + (-3*u3+3*u2+3*u+1)*q3[i] + (u3)*q4[i])
			end
			return ret
		end
		return MAP(map_fn)(domain)
	end
	return CUBICUBSPLINE0
end

# //////////////////////////////////////////////////////////////
function CUBICCARDINAL(domain, h=1)
	function CUBICCARDINAL0(args)
		q1_fn, q2_fn, q3_fn, q4_fn = args
		function map_fn(point)
			u = S1(point)
			u2 = u*u
			u3 = u2*u
			q1, q2, q3, q4 = [ isa(f, Function) ? f(point) : f for f in [q1_fn, q2_fn, q3_fn, q4_fn]]
			ret = [0.0 for i in 1:length(q1)]
			for i in 1:length(ret)
					ret[i] = (-h*u3+2*h*u2-h*u)*q1[i] + ((2-h)*u3+(h-3)*u2+1)*q2[i] + ((h-2)*u3+(3-2*h)*u2+h*u)*q3[i] + (h*u3-h*u2)*q4[i]
			end
			return ret
		end
		return MAP(map_fn)(domain)
	end
	return CUBICCARDINAL0
end

# //////////////////////////////////////////////////////////////
function SPLINE(curve)
	function SPLINE0(points)
		ret = Vector{Hpc}()
		for i in 1:length(points)-4+1
			P = points[i:i+4-1]
			push!(ret, curve(P))
		end
		return Struct(ret)
	end
	return SPLINE0
end

# //////////////////////////////////////////////////////////////
function JOINTS(curve)
	knotzero = MK([0])
	function JOINTS0(points)
		points, cells, pols = UKPOL(SPLINE(curve(knotzero)))
		return POLYMARKER(2)(points)
	end
	return JOINTS0
end

# //////////////////////////////////////////////////////////////
function BERNSTEINBASIS(U)
	function BERNSTEIN0(N)
		function BERNSTEIN1(I)
			function map_fn(point)
					t = U(point)
					ret = CHOOSE([N, I])*((1-t)^(N-I))*(t^I)
					return ret
			end
			return map_fn
		end
		return [BERNSTEIN1(I) for I in 0:N]
	end
	return BERNSTEIN0
end

# //////////////////////////////////////////////////////////////
function TENSORPRODSURFACE(args)
	ubasis, vbasis = args
	function TENSORPRODSURFACE0(controlpoints_fn)
		controlpoints_fn=ToFloat64(controlpoints_fn)
		function map_fn(point)
			u, v = point
			U = [f([u]) for f in ubasis]
			V = [f([v]) for f in vbasis]
			controlpoints = [ isa(f, Function) ? f(point) : f for f in controlpoints_fn]
			target_dim = length(controlpoints[1][1])
			ret = [0.0 for x in 1:target_dim]
			for i in 1:length(ubasis)
					for j in 1:length(vbasis)
						for M in 1:target_dim
							ret[M] += U[i]*V[j] * controlpoints[i][j][M]
						end
					end
			end
			return ret
		end
		return map_fn
	end
	return TENSORPRODSURFACE0
end

# //////////////////////////////////////////////////////////////
function BILINEARSURFACE(controlpoints)
	return TENSORPRODSURFACE([BERNSTEINBASIS(S1)(1), BERNSTEINBASIS(S1)(1)])(controlpoints)
end

# //////////////////////////////////////////////////////////////
function BIQUADRATICSURFACE(controlpoints)

	function u0(point) 
		u = S1(point)
		return 2*u*u-u
	end
	function u1(point) 
		u = S1(point)
		return 4*u-4*u*u
	end
	function u2(point) 
		u = S1(point)
		return 2*u*u-3*u+1
	end
	basis = [u0, u1, u2]
	return TENSORPRODSURFACE([basis, basis])(controlpoints)
end

# //////////////////////////////////////////////////////////////
function HERMITESURFACE(controlpoints)
	function H0(point) 
		u = S1(point)
		u2 = u*u
		u3 = u2*u
		return u3-u2
	end
	function H1(point) 
		u = S1(point)
		u2 = u*u
		u3 = u2*u
		return u3-2*u2+u
	end
	function H2(point) 
		u = S1(point)
		u2 = u*u
		u3 = u2*u
		return 3*u2-2*u3
	end
	function H3(point) 
		u = S1(point)
		u2 = u*u
		u3 = u2*u
		return 2*u3-3*u2+1
	end

	basis = [H3, H2, H1, H0]
	return TENSORPRODSURFACE([basis, basis])(controlpoints)
end

# //////////////////////////////////////////////////////////////
function BEZIERSURFACE(controlpoints)
	M = length(controlpoints) - 1
	N = length(controlpoints[1]) - 1
	return TENSORPRODSURFACE([BERNSTEINBASIS(S1)(M), BERNSTEINBASIS(S1)(N)])(controlpoints)
end

# //////////////////////////////////////////////////////////////
function TENSORPRODSOLID(args)
	ubasis, vbasis, wbasis = args
	function TENSORPRODSOLID0(controlpoints_fn)
		controlpoints_fn=ToFloat64(controlpoints_fn)
		function map_fn(point)
			u, v, w = point
			U = [f([u]) for f in ubasis]
			V = [f([v]) for f in vbasis]
			W = [f([w]) for f in wbasis]
			controlpoints = [ isa(f, Function) ? f(point) : f for f in controlpoints_fn]
			target_dim = length(controlpoints[1][1][1])
			ret = [0 for x in 1:target_dim]
			for i in 1:length(ubasis)
					for j in 1:length(vbasis)
						for k in 1:length(wbasis)
							for M in 1:target_dim
									ret[M] += U[i]*V[j]*W[k] * controlpoints[M][i][j][k]
							end
						end
					end
			end
			return ret
		end
		return map_fn
	end
	return TENSORPRODSOLID0
end

# //////////////////////////////////////////////////////////////
function BEZIERMANIFOLD(degrees)
	basis = [BERNSTEINBASIS(S1)(d) for d in degrees]
	return TENSORPRODSOLID(basis)
end

function LOCATE(args)
	pol, a, distances = args
	ret = []
	for d in distances
		push!(ret, T(a)(d), pol)
	end
	return STRUCT(ret)
end

function RIF(size)
	thin = 0.01*size
	x = COLOR(RED)(CUBOID([size, thin, thin]))
	y = COLOR(GREEN)(CUBOID([thin, size, thin]))
	z = COLOR(BLUE)(CUBOID([thin, thin, size]))
	return STRUCT([x, y, z])
end

# //////////////////////////////////////////////////////////////
function FRACTALSIMPLEX(D)
	function FRACTALSIMPLEX0(N)
		mkpols = COMP([COMP([COMP([COMP([STRUCT, AA(MKPOL)]), AA(AL)]), DISTR]), CONS([ID, K([[FROMTO([1,D+1])], [[1]]])])])
		function COMPONENT(args)
			i, seq = args
			firstseq = seq[1:i-1]
			pivot  = seq[i-1]
			lastseq = seq[i:length(seq)]
			firstpart = AA(MEANPOINT)(DISTR([firstseq, pivot]))
			lastpart  = AA(MEANPOINT)(DISTR([lastseq , pivot]))
			return CAT([firstpart, [pivot], lastpart])
		end
		expand = COMP([COMP([AA(COMPONENT), DISTR]), CONS([COMP([INTSTO, LEN]), ID])])
		splitting = COMP([COMP, DIESIS(N)])(COMP([CAT, AA(expand)]))
		return COMP([COMP([COMP([COMP([mkpols, splitting]), CONS([S1])])])])(UKPOL(SIMPLEX(D)))
	end
	return FRACTALSIMPLEX0
end

# //////////////////////////////////////////////////////////////
function PYRAMID(H)
	function PYRAMID0(pol)
		barycenter = MEANPOINT(UKPOL(pol)[1])
		return JOIN([MK(barycenter+[H]), pol])
	end
	return PYRAMID0
end

function MESH(seq)
	return INSL(RAISE(PROD))([QUOTE(i) for i in seq])
end

function NU_GRID(data)
	polylines = [POLYLINE(i) for i in data]
	return INSL(RAISE(PROD))(polylines)
end


function SEGMENT(sx)
	function SEGMENT0(args)
		N = length(args[1])
		A, B = args
		P0 = A
		P1 = [A[i]+(B[i]-A[i])*sx for i in 1:N]
		return POLYLINE([P0,P1])
	end
	return SEGMENT0
end

# //////////////////////////////////////////////////////////////
function SOLIDIFY(pol)
	box = box(pol)
	min = box.p1[1]
	max = box.p2[1]
	siz = max-min
	far_point = max+siz*100
	function InftyProject(pol)
		verts, cells, pols = UKPOL(pol)
		verts = [[far_point] + v[2:end] for v in verts]
		return MKPOL(verts, cells, pols)
	end
	function IsFull(pol)
		return DIM(pol) == RN(pol)
	end
	ret = SPLITCELLS(pol)
	ret = [JOIN([pol, InftyProject(pol)]) for pol in ret]
	return XOR(FILTER(IsFull)(ret))
end

# //////////////////////////////////////////////////////////////
function EXTRUSION(angle)
	function EXTRUSION1(height)
		function EXTRUSION0(pol)
			dim = DIM(pol)
			cells = SPLITCELLS(SKELETON(dim)(pol))
			slice = [EMBED(1)(c) for c in cells]
			tensor = COMP([T(dim+1)(1.0/height), R([dim-1,dim])(angle/height)])
			layer = STRUCT([JOIN([p,tensor(p)]) for p in slice])
			return COMP([COMP([STRUCT, CAT]), DIESIS(height)])([layer, tensor])
		end
		return EXTRUSION0
	end
	return EXTRUSION1
end

# //////////////////////////////////////////////////////////////
function EX(args)
	x1, x2 = args
	function EX0(pol)
		dim = DIM(pol)
		return T(dim+1)(x1)(S(dim+1)(x2-x1)(EXTRUSION(0.0)(1.0)(pol)))
	end
	return EX0
end

# //////////////////////////////////////////////////////////////
function LEX(args)
	x1, x2 = args
	function LEX0(pol)
		function SHEARTENSOR(A)
			function SHEARTENSOR0(POL)
					dim = DIM(POL)
					newrow = K((AR([CAT([[0, 1],DIESIS((dim-2))(0)]),A])))
					update = (COMP([CONS, CAT]))([[S1, newrow],AA(SEL)((FROMTO([3,dim+1])))])
					matrix = update(IDNT(dim+1))        
					return(MAT(matrix))(POL)
			end
			return SHEARTENSOR0
		end
		ret = EXTRUSION(0)(1)(pol)
		ret = SHEARTENSOR(x2-x1)(ret)
		ret = S(DIM(pol)+1)(x2-x1)(ret)
		ret = T(DIM(pol)+1)(x1)(ret)
		return ret
	end
	return LEX0
end

# //////////////////////////////////////////////////////////////
function SEX(args)
	x1, x2 = args
	function SEX1(height)
		function SEX0(pol)
			dim = DIM(pol)
			ret = EXTRUSION(x2-x1)(height)(pol)
			ret = S(dim+1)(x2-x1)(ret)
			ret = R([dim,dim-1])(x1)(ret)
			return ret
		end
		return SEX0
	end
	return SEX1
end

# //////////////////////////////////////////////////////////////
function UKPOLF(pol)
	error("not implemented")
end

# //////////////////////////////////////////////////////////////
function POLAR(pol, precision=1e-6)
	 faces, cells, pols = UKPOLF(pol)
	 for i in 1:length(faces)
		  mod = -1*faces[i][1]
		  if abs(mod) < precision
				mod = 1
		  end
		  faces[i] = [value/mod for value in faces[i][2:end]]
	 end
	 return MKPOL(faces, cells, pols)
end

# //////////////////////////////////////////////////////////////
function SWEEP(v)
	function SWEEP0(pol)
		ret = Power(pol, QUOTE([1]))
		mat = IDNT(length(v)+2)
		for i in 1:length(v)
			mat[i+1][length(v)+1] = v[i]
		end
		ret = MAT(mat)(ret)
		return PROJECT(1)(ret)
	end
	return SWEEP0
end

# //////////////////////////////////////////////////////////////
function MINKOWSKI(vects)
	function MINKOWSKI0(pol)
		ret = pol
		for i in range(length(vects)-1,-1,-1)
			ret = SWEEP(vects[i])(ret)
		end
		return ret
	end
	return MINKOWSKI0
end

# //////////////////////////////////////////////////////////////
function OFFSET(v)
	function OFFSET0(pol)
		ret = pol
		for i in 1:length(v)
			shear = [j!=i ? 0 : v[i] for j in 1:length(v)] + [0 for j in 1:i]
			mat = IDNT(length(shear)+2)
			for i in 1:length(shear)
				mat[i+1][length(shear)+1] = shear[i]
			end
			ret = MAT(mat)((Power(ret, QUOTE([1]))))
		end
		return PROJECT(length(v))(ret)
	end
	return OFFSET0
end

# //////////////////////////////////////////////////////////////
function THINSOLID(surface, delta=1e-4)
	function map_fn(point)
		u, v, w = point
		P0 = surface([u, v])
		PX = surface([u+delta, v])
		PY = surface([u, v+delta])
		GX = [PX[i]-P0[i] for i in 1:3]
		GY = [PY[i]-P0[i] for i in 1:3]
		normal = UNITVECT(VECTPROD([GX, GY]))
		ret = [P0[i]+w*normal[i] for i in 1:3]
		return ret
	end
	return map_fn
end

# //////////////////////////////////////////////////////////////
function PLANE(args)
	p0, p1, p2 = args
	v1 = VECTDIFF([p1, p0])
	v2 = VECTDIFF([p2, p0])
	side1 = VECTNORM(v1)
	side2 = VECTNORM(v2)
	normal = UNITVECT(VECTPROD([v1, v2]))
	axis = VECTPROD([[0, 0, 1], normal])
	angle = acos(INNERPROD([[0, 0, 1], normal]))
	geometry = T([1,2,3])(p0)(ROTN([angle, axis])(T([1,2])([-1*side1,-1*side2])(CUBOID([2*side1, 2*side2]))))
	return  [normal, p0, geometry]
end

# //////////////////////////////////////////////////////////////
function RATIONALBEZIER(controlpoints_fn)
	controlpoints_fn=ToFloat64(controlpoints_fn)
	degree = length(controlpoints_fn)-1
	basis = BERNSTEINBASIS(S1)(degree)
	function map_fn(point)
		controlpoints = [ isa(f, Function) ? f(point) : f for f in controlpoints_fn]
		target_dim = length(controlpoints[1])
		ret = [0 for i in 1:target_dim]
		for i in 1:length(basis)
			coeff = basis[i](point)
			for M in 1:target_dim
				ret[M] += coeff * controlpoints[i][M] 
			end
		end
		last = ret[end]
		if last != 0
			ret = [value/last for value in ret]
		end
		ret = ret[1:end-1]
		return  ret
	end
	return map_fn
end

# //////////////////////////////////////////////////////////////
function ELLIPSE(args)
	A, B = args
	function ELLIPSE0(N)
		C = 0.5*sqrt(2)
		mapping = RATIONALBEZIER([[A, 0, 1], [A*C, B*C, C], [0, B, 1]])
		quarter = MAP(mapping)((INTERVALS(1.0)(N)))
		half = STRUCT([quarter, S(2)(-1)(quarter)])
		return STRUCT([half, S(1)(-1)(half)])
	end
	return ELLIPSE0
end

# //////////////////////////////////////////////////////////////
function CURVE_NORMAL(curve)
	function map_fn(point)
		xu, yu = curve(point)
		mod2 = xu*xu+yu*yu
		den = mod2 > 0 ? sqrt(mod2) : 0
		return [-yu/den, xu/den]
	end
	return map_fn
end

# //////////////////////////////////////////////////////////////
function DERBEZIER(controlpoints_fn)
	controlpoints_fn=ToFloat64(controlpoints_fn)
	degree = length(controlpoints_fn)-1
	
	function DERBERNSTEIN(N)
		function DERBERNSTEIN0(I)
			function map_fn(point)
				t = S1(point)
				return CHOOSE([N,I]) * t^(I-1) * (1-t)^(N-I-1) * (I-N*t)
			end
			return  map_fn
		end
		return DERBERNSTEIN0
	end

	basis = [DERBERNSTEIN(degree)(i) for i in 1:degree+1]
	function map_fn(point)
		controlpoints = [isa(f, Function)  ? f(point) : f for f in controlpoints_fn]
		target_dim = length(controlpoints[1])
		ret = [0 for i in 1:target_dim]
		for i in 1:length(basis)
			coeff = basis[i](point)
			for M in 1:target_dim
				ret[M] += coeff * controlpoints[i][M] 
			end
		end
		return ret
	end
	return map_fn
end

# //////////////////////////////////////////////////////////////
function BEZIERSTRIPE(args)
	controlpoints, width, n = args
	bezier = BEZIERCURVE(controlpoints)
	normal = CURVE_NORMAL(DERBEZIER(controlpoints))
	function map_fn(point)
		u, v = point
		bx, by = bezier(point)
		nx, ny = normal(point)
		ret = [bx+v*nx, by+v*ny]
		return ret
	end
	domain = S(2)(width)(T(1)(0.00001)(Power(INTERVALS(1.0)(n),INTERVALS(1.0)(1))))
	return MAP(map_fn)(domain)
end

# //////////////////////////////////////////////////////////////
function BSPLINE(degree)
	function BSPLINE0(knots)
		function BSPLINE1(points_fn)
			n = length(points_fn)-1
			m = length(knots)-1
			k = degree+1
			T = knots
			tmin, tmax = T[k], T[n+1]
			
			if length(knots) != (n+k+1)
					error("Invalid point/knots/degree for bspline!")
			end
			
			function N(i, k, t)
					
					if k == 1
						if (t >= T[i] && t < T[i+1]) || (t == tmax && t >= T[i] && t <= T[i+1])
							return 1
						else
							return 0
						end
					end
					
					ret = 0
					num1, div1 = t-T[i], T[i+k-1]-T[i]  
					if div1 != 0
						ret += (num1/div1) * N(i, k-1, t)
					end
					num2, div2 = T[i+k]-t, T[i+k]-T[i+1]
					if div2 != 0
						ret += (num2/div2) * N(i+1, k-1, t)
					end
					return ret
			end
			
			function map_fn(point)
					t = point[1]
					
					points = [ isa(f, Function) ? f(point) : f for f in points_fn]
					target_dim = length(points[1])
					ret = [0 for i in 1:target_dim]
					for i in 1:n+1
						coeff = N(i, k, t) 
						for M in 1:target_dim
							ret[M] += points[i][M] * coeff
						end
					end
					return ret
			end
			return map_fn
		end
		return BSPLINE1
	end
	return BSPLINE0
end

# //////////////////////////////////////////////////////////////
function NUBSPLINE(degree, totpoints=80)
	function NUBSPLINE1(knots)
		function NUBSPLINE2(points)
			m = length(knots)
			tmin = min(knots...)
			tmax = max(knots...)
			tsiz = tmax-tmin
			v = [tsiz/float(totpoints-1) for i in 1:totpoints-1]
			@assert length(v)+1 == totpoints
			v = [-tmin] + v
			domain = QUOTE(v)
			return MAP(BSPLINE(degree)(knots)(points))(domain)
		end
		return NUBSPLINE2
	end
	return NUBSPLINE1
end

# //////////////////////////////////////////////////////////////
function DISPLAYNUBSPLINE(args, marker_size=0.1)
	degree, knots, points = args
	spline_view_knots = POLYMARKER(2, marker_size)(UKPOL(NUBSPLINE(degree, length(knots))(knots)(points))[1])
	return  STRUCT([
		degree > 0 ? 
			NUBSPLINE(degree)(knots)(points) : 
			POLYMARKER(3, marker_size)(points)
		,spline_view_knots
		,POLYLINE(points)
		,POLYMARKER(1, marker_size)(points)
	])
end

# //////////////////////////////////////////////////////////////
function RATIONALBSPLINE(degree)
	function RATIONALBSPLINE0(knots)
		function RATIONALBSPLINE1(points)
			bspline = BSPLINE(degree)(knots)(points)
			function map_fn(point)
					ret = bspline(point)
					last = ret[end]
					if last != 0
						ret = [value/last for value in ret]
					end
					ret = ret[1:end-1]
					return ret
			end
			return map_fn
		end
		return RATIONALBSPLINE1
	end
	return RATIONALBSPLINE0
end

# //////////////////////////////////////////////////////////////
function NURBSPLINE(degree, totpoints=80)
	function NURBSPLINE1(knots)
		function NURBSPLINE2(points)
			m = length(knots)
			tmin = min(knots...)
			tmax = max(knots...)
			tsiz = tmax-tmin
			v = [tsiz/float(totpoints-1) for i in 1:totpoints-1]
			@assert length(v)+1 == totpoints
			v = [-tmin] + v
			domain = QUOTE(v)
			return MAP(RATIONALBSPLINE(degree)(knots)(points))(domain)
		end
		return NURBSPLINE2
	end
	return NURBSPLINE1
end

# //////////////////////////////////////////////////////////////
function DISPLAYNURBSPLINE(args, marker_size=0.1)
	degree, knots, points = args
	spline_view_knots = POLYMARKER(2, marker_size)(UKPOL(NURBSPLINE(degree, length(knots))(knots)(points))[1])
	return  STRUCT([
		degree > 0 ? 
			NURBSPLINE(degree)(knots)(points) :
			POLYMARKER(3, marker_size)(points)
		,spline_view_knots
		,POLYLINE(points)
		,POLYMARKER(1, marker_size)(points)
	])
end

# /////////////////////////////////////////////////////////////
function TestSphere()
	VIEW(SPHERE(1.0)([16,16]), "TestSphere")
end

function TestTorus()
	VIEW(TORUS([1.0,2.0])([20,20]), "TestTorus")
end

function TestBezier()
	VIEW(MAP(BEZIER(S1)([[-0,0],[1,0],[1,1],[2,1],[3,1]]))(INTERVALS(1.0)(32)), "TestBezier-1")
	C0 = BEZIER(S1)([[0,0,0],[10,0,0]])
	C1 = BEZIER(S1)([[0,2,0],[8,3,0],[9,2,0]])
	C2 = BEZIER(S1)([[0,4,1],[7,5,-1],[8,5,1],[12,4,0]])
	C3 = BEZIER(S1)([[0,6,0],[9,6,3],[10,6,-1]])
	VIEW(MAP(BEZIER(S2)([C0,C1,C2,C3]))(Power(INTERVALS(1.0)(10),INTERVALS(1.0)(10))), "TestBezier-2")
end

function TestCoonsPatch()
	Su0 = BEZIER(S1)([[0,0,0],[10,0,0]])
	Su1 = BEZIER(S1)([[0,10,0],[2.5,10,3],[5,10,-3],[7.5,10,3],[10,10,0]])
	Sv0 = BEZIER(S2)([[0,0,0],[0,0,3],[0,10,3],[0,10,0]])
	Sv1 = BEZIER(S2)([[10,0,0],[10,5,3],[10,10,0]])
	VIEW(MAP(COONSPATCH([Su0,Su1,Sv0,Sv1]))(Power(INTERVALS(1.0)(10),INTERVALS(1.0)(10))), "TestCoonsPatch")
end

function TestRuledSurface()
	alpha = point -> [point[1],point[1],       0 ]
	beta = point -> [      -1,      +1,point[1] ]
	domain = T([1,2])([-1,-1])(Power(INTERVALS(2.0)(10),INTERVALS(2.0)(10)))
	VIEW(MAP(RULEDSURFACE([alpha,beta]))(domain), "TestRuledSurface")
end

function TestProfileProdSurface()
	alpha = BEZIER(S1)([[0.1,0,0],[2,0,0],[0,0,4],[1,0,5]])
	beta = BEZIER(S2)([[0,0,0],[3,-0.5,0],[3,3.5,0],[0,3,0]])
	domain = Power(INTERVALS(1.0)(20),INTERVALS(1.0)(20))
	VIEW(Struct([MAP(alpha)(domain),MAP(beta )(domain),MAP(PROFILEPRODSURFACE([alpha,beta]))(domain)]), "TestProfileProdSurface")
end

function TestRotationalSurface()
	profile = BEZIER(S1)([[0,0,0],[2,0,1],[3,0,4]]) 
	domain = Power(INTERVALS(1.0)(10),INTERVALS(2*PI)(30)) 
	VIEW(MAP(ROTATIONALSURFACE(profile))(domain), "TestRotationalSurface")
end

function TestCylindricalSurface()
	alpha = BEZIER(S1)([[1,1,0],[-1,1,0],[1,-1,0],[-1,-1,0]])
	Udomain = INTERVALS(1.0)(20)
	Vdomain = INTERVALS(1.0)(6)
	domain = Power(Udomain,Vdomain)
	fn = CYLINDRICALSURFACE([alpha,[0.0,0.0,1.0]])
	VIEW(MAP(fn)(domain), "TestCylindricalSurface")
end

function TestConicalSurface()
	domain = Power(INTERVALS(1.0)(20),INTERVALS(1.0)(6))
	beta = BEZIER(S1)([ [1,1,0],[-1,1,0],[1,-1,0],[-1,-1,0] ])
	VIEW(MAP(CONICALSURFACE([[0,0,1],beta]))(domain), "TestConicalSurface")
end

function TestCubicHermite()
	domain = INTERVALS(1.0)(20)
	VIEW(Struct([
		MAP(CUBICHERMITE(S1)([[1,0],[1,1],[ -1, 1],[ 1,0]]))(domain),
		MAP(CUBICHERMITE(S1)([[1,0],[1,1],[ -2, 2],[ 2,0]]))(domain),
		MAP(CUBICHERMITE(S1)([[1,0],[1,1],[ -4, 4],[ 4,0]]))(domain),
		MAP(CUBICHERMITE(S1)([[1,0],[1,1],[-10,10],[10,0]]))(domain)])
		, "TestCubicHermite-1")
	c1 = CUBICHERMITE(S1)([[1  ,0,0],[0  ,1,0],[0,3,0],[-3,0,0]])
	c2 = CUBICHERMITE(S1)([[0.5,0,0],[0,0.5,0],[0,1,0],[-1,0,0]])
	sur3 = CUBICHERMITE(S2)([c1,c2,[1,1,1],[-1,-1,-1]])
	domain = Power(INTERVALS(1.0)(14),INTERVALS(1.0)(14))
	VIEW(MAP(sur3)(domain), "TestCubicHermite-2")
end

function TestPermutahedron()
	VIEW(Struct([PERMUTAHEDRON(2),(PERMUTAHEDRON(2))]), "TestPermutahedron-1")
	VIEW(Struct([PERMUTAHEDRON(3),(PERMUTAHEDRON(3))]), "TestPermutahedron-2")
end

function TestSchegel3d()
	VIEW(SCHLEGEL3D(0.2)((T([1,2,3,4])([-1.0/3.0,-1.0/3.0,-1,+1])(SIMPLEX(4)))), "TestSchegel3d-1")
	VIEW(SCHLEGEL3D(0.2)((T([1,2,3,4])([-1,-1,-1,1])(CUBOID([2,2,2,2])))), "TestSchegel3d-2")
	VIEW(SCHLEGEL3D(0.2)((T([1,2,3,4])([-1.0/3.0,-1.0/3.0,-1,+1])(Power(SIMPLEX(2),SIMPLEX(2))))), "TestSchegel3d-3")
end

function TestCubicSpline()
	domain = INTERVALS(1.0)(20)
	points = [
		[-3.0,6.0],
		[-4.0,2.0],
		[-3.0,-1.0],
		[-1.0,1.0],
		[1.5,1.5],
		[3.0,4.0],
		[5.0,5.0],
		[7.0,2.0],
		[6.0,-2.0],
		[2.0,-3.0]
		]
	VIEW(SPLINE(CUBICCARDINAL(domain))(points), "TestCubicSpline-1")
	VIEW(SPLINE(CUBICUBSPLINE(domain))(points), "TestCubicSpline-2")
end

function TestBilinarSurface()
	controlpoints = [
		[[0.0,0.0,0.0],[2.0,-4.0,2.0]],
		[[0.0,3.0,1.0],[4.0,0.0,0.0]]
	]
	domain = Power(INTERVALS(1.0)(10),INTERVALS(1.0)(10))
	mapping = BILINEARSURFACE(controlpoints)
	VIEW(MAP(mapping)(domain), "TestBilinarSurface")
end

function TestBiquadraticSurface()
	controlpoints = [[[0,0,0],[2,0,1],[3,1,1]],[[1,3,-1],[3,2,0],[4,2,0]],[[0,9,0],[2,5,1],[3,3,2]]]
	domain = Power(INTERVALS(1.0)(10),INTERVALS(1.0)(10))
	mapping = BIQUADRATICSURFACE(controlpoints)
	VIEW(MAP(mapping)(domain), "TestBiquadraticSurface")
end



function TestHermiteSurface()
	controlpoints = ToFloat64([[[0,0,0 ],[2,0,1],[3,1,1],[4,1,1]],[[1,3,-1],[3,2,0],[4,2,0],[4,2,0]],[[0,4,0 ],[2,4,1],[3,3,2],[5,3,2]],[[0,6,0 ],[2,5,1],[3,4,1],[4,4,0]]])
	println("!!!",controlpoints)
	domain = Power(INTERVALS(1.0)(10),INTERVALS(1.0)(10))
	mapping = HERMITESURFACE(controlpoints)
	VIEW(MAP(mapping)(domain), "TestHermiteSurface")
end

function TestBezierSurface()
	controlpoints = [[[ 0,0,0],[0 ,3  ,4],[0,6,3],[0,10,0]],[[ 3,0,2],[2 ,2.5,5],[3,6,5],[4,8,2]],[[ 6,0,2],[8 ,3 , 5],[7,6,4.5],[6,10,2.5]],[[10,0,0],[11,3  ,4],[11,6,3],[10,9,0]]]
	domain = Power(INTERVALS(1.0)(10),INTERVALS(1.0)(10))
	mapping = BEZIERSURFACE(controlpoints)
	VIEW(MAP(mapping)(domain), "TestBezierSurface")
end

function TestBezierManifold()
	grid1D = INTERVALS(1.0)(5)
	domain3D = Power(Power(grid1D,grid1D),grid1D)
	degrees = [2,2,2]
	Xtensor =  [[[0,1,2],[-1,0,1],[0,1,2]],[[0,1,2],[-1,0,1],[0,1,2]],[[0,1,2],[-1,0,1],[0,1,2]]]
	Ytensor =  [[[0,0,0.8],[1,1,1],[2,3,2]],[[0,0,0.8],[1,1,1],[2,3,2]],[[0,0,0.8],[1,1,1],[2,3,2]]]
	Ztensor =  [[[0,0,0],[0,0,0],[0,0,0]],[[1,1,1],[1,1,1],[1,1,1]],[[2,2,1],[2,2,1],[2,2,1]]] 
	mapping = BEZIERMANIFOLD(degrees)([Xtensor,Ytensor,Ztensor])
	VIEW(MAP(mapping)(domain3D), "TestBezierManifold")
end

function TestOffset()
	verts = [[0,0,0],[3,0,0],[3,2,0],[0,2,0],[0,0,1.5],[3,0,1.5],[3,2,1.5],[0,2,1.5],[0,1,2.2],[3,1,2.2]]
	cells = [[1,2],[2,3],[3,4],[4,1],[5,6],[6,7],[7,8],[8,5],[1,5],[2,6],[3,7],[4,8],[5,9],[8,9],[6,10],[7,10], [9,10]]
	pols = [[1]]
	House = MKPOL(verts,cells,pols)
	VIEW(Struct([OFFSET([0.1,0.2,0.1])(House), T(1)(1.2*SIZE(1)(House))(House)]), "TestOffset")
end

function TestThinSolid()
	Su0 = COMP([BEZIERCURVE([[0,0,0],[10,0,0]]),CONS([S1])])
	Su1 = COMP([BEZIERCURVE([[0,10,0],[2.5,10,3],[5,10,-3],[7.5,10,3],[10,10,0]]),CONS([S1]) ])
	S0v = COMP([BEZIERCURVE([[0,0,0],[0,0,3],[0,10,3],[0,10,0]]) , CONS([S2]) ]) 
	S1v = COMP([BEZIERCURVE([[10,0,0],[10,5,3],[10,10,0]]) ,CONS([S2])   ])
	surface = COONSPATCH([Su0,Su1,S0v,S1v])
	VIEW(MAP(surface)(Power(INTERVALS(1.0)(10),INTERVALS(1.0)(10))), "TestThinSolid-1")
	solidMapping = THINSOLID(surface)
	Domain3D = Power(Power(INTERVALS(1.0)(5),INTERVALS(1.0)(5)),INTERVALS(0.5)(5))
	VIEW(MAP(solidMapping)(Domain3D), "TestThinSolid-2")
end

function TestEllipse()
	VIEW(ELLIPSE([1,2])(8), "TestEllipse")
end

function TestBezierStripe()
	vertices = [[0,0],[1.5,0],[-1,2],[2,2],[2,0]]
	VIEW(Struct([POLYLINE(vertices),Power(BEZIERSTRIPE([vertices,0.25,22]),QUOTE([0.9]))]), "TestBezierStripe")
end

function TestDisplayNubSpline()
	ControlPoints = [[0,0],[-1,2],[1,4],[2,3],[1,1],[1,2],[2.5,1], [2.5,3], [4,4],[5,0]]
	VIEW(DISPLAYNUBSPLINE([3,[0,0,0,0, 1,2,3,4,5, 6    ,7,7,7,7], ControlPoints]), "TestDisplayNubSpline")
end

function TestDisplayNurbsSpline()
	knots = [0,0,0,1,1,2,2,3,3,4,4,4]
	_p=sqrt(2)/2.0
	controlpoints = [[-1,0,1], [-_p,_p,_p], [0,1,1], [_p,_p,_p],[1,0,1], [_p,-_p,_p], [0,-1,1], [-_p,-_p,_p], [-1,0,1]]
	VIEW(DISPLAYNURBSPLINE([2, knots, controlpoints]), "TestDisplayNurbsSpline")
end

function TestMinkowski()
	p = MKPOL(
		[[0.0,0.0]],
		[[1]],
		[[1]])
	B = MINKOWSKI([  [-1.0/2.0,-1*sqrt(3.0/2.0)] , [-1.0/2.0,sqrt(3.0/2.0)] , [1,0] ])(p)
	vertices = [[0,0],[1,0],[1,0.5],[0.5,0.5],[0.5,1],[0,1]]
	pol1D = MKPOL(vertices,[[1,2],[2,3],[3,4],[4,5],[5,6],[6,1]],[[1],[2],[3],[4],[5],[6]])
	pol2D = MKPOL(vertices,[[1,2,3,4],[4,5,6,1]],[[1,2]])
	Min0 = STRUCT([T([1,2])(v)(S([1,2])([0.1,0.1])(B)) for v in vertices ])
	Min1 = MINKOWSKI([[0.1*-1.0/2.0,0.1*-1*sqrt(3.0/2.0)],[0.1*-1.0/2.0,0.1*sqrt(3.0/2.0)] , [0.1*1,0.1*0] ])(pol1D)
	Min2 = MINKOWSKI([[0.1*-1.0/2.0,0.1*-1*sqrt(3.0/2.0)],[0.1*-1.0/2.0,0.1*sqrt(3.0/2.0)] , [0.1*1,0.1*0] ])(pol2D)
	VIEW(Struct([Min0,Min1,Min2]), "TestMinkowski")
end

function TestThinsolid()
	Su0 = COMP([BEZIERCURVE([[0,0,0],[10,0,0]]),CONS([S1])])
	Su1 = COMP([BEZIERCURVE([[0,10,0],[2.5,10,3],[5,10,-3],[7.5,10,3],[10,10,0]]),CONS([S1]) ])
	S0v = COMP([BEZIERCURVE([[0,0,0],[0,0,3],[0,10,3],[0,10,0]]) , CONS([S2]) ]) 
	S1v = COMP([BEZIERCURVE([[10,0,0],[10,5,3],[10,10,0]]) ,CONS([S2])   ])
	surface = COONSPATCH([Su0,Su1,S0v,S1v])
	VIEW(MAP(surface)(Power(INTERVALS(1.0)(10),INTERVALS(1.0)(10))), "TestThinsolid-1")
	solidMapping = THINSOLID(surface)
	Domain3D = Power(Power(INTERVALS(1.0)(5),INTERVALS(1.0)(5)),INTERVALS(0.5)(5))
	VIEW(MAP(solidMapping)(Domain3D), "TestThinsolid-2")
end

function TestEllipse()
	VIEW(ELLIPSE([1,2])(8), "TestEllipse")
end

function TestBezierStripe()
	vertices = [[0,0],[1.5,0],[-1,2],[2,2],[2,0]]
	VIEW(Struct([POLYLINE(vertices),Power(BEZIERSTRIPE([vertices,0.25,22]),QUOTE([0.9]))]), "TestBezierStripe")
end

function TestDisplayNubSpline()
	ControlPoints = [[0,0],[-1,2],[1,4],[2,3],[1,1],[1,2],[2.5,1], [2.5,3], [4,4],[5,0]]
	VIEW(DISPLAYNUBSPLINE([3,[0,0,0,0, 1,2,3,4,5, 6,7,7,7,7], ControlPoints]), "TestDisplayNubSpline")
end

function TestDisplayNurbsSpline()
	knots = [0,0,0,1,1,2,2,3,3,4,4,4]
	_p = sqrt(2)/2.0
	controlpoints = [[-1,0,1], [-_p,_p,_p], [0,1,1], [_p,_p,_p],[1,0,1], [_p,-_p,_p], [0,-1,1], [-_p,-_p,_p], [-1,0,1]]
	VIEW(DISPLAYNURBSPLINE([2, knots, controlpoints]), "TestDisplayNurbsSpline")
end

function TestMinkowski()
	p = MKPOL([[0.0,0.0]],[[1]],[[1]])
	B = MINKOWSKI([[-1.0/2.0,-sqrt(3.0/2.0)], [-1.0/2.0,sqrt(3.0/2.0)], [1,0]])(p)
	vertices = [[0.0,0.0],[1.0,0.0],[1.0,0.5],[0.5,0.5],[0.5,1.0],[0.0,1.0]]
	pol1D = MKPOL(vertices,[[1,2],[2,3],[3,4],[4,5],[5,6],[6,1]],[[1],[2],[3],[4],[5],[6]])
	pol2D = MKPOL(vertices,[[1,2,3,4],[4,5,6,1]],[[1,2]])
	Min0 = STRUCT([T([1,2])(v)(S([1,2])([0.1,0.1])(B)) for v in vertices])
	Min1 = MINKOWSKI([[0.1*-1.0/2.0,0.1*-sqrt(3.0/2.0)],[0.1*-1.0/2.0,0.1*sqrt(3.0/2.0)],[0.1*1,0.1*0]])(pol1D)
	Min2 = MINKOWSKI([[0.1*-1.0/2.0,0.1*-sqrt(3.0/2.0)],[0.1*-1.0/2.0,0.1*sqrt(3.0/2.0)],[0.1*1,0.1*0]])(pol2D)
	A = Power(Min2,Q(0.05))
	B = Power(Min0,Q(0.70))
	C = Power(Min1,Q(0.05))
	VIEW(TOP([TOP([A,B]),C]), "TestMINKOWSKI")
end

function TestPolar()
	VIEW(POLAR(CUBOID([1,1,1])), "TestPolar")
end

function TestSolidify()
VIEW(SOLIDIFY(STRUCT(AA(POLYLINE)([
		[[0,0],[4,2],[2.5,3],[4,5],[2,5],[0,3],[-3,3],[0,0]],
		[[0,3],[0,1],[2,2],[2,4],[0,3]],
		[[2,2],[1,3],[1,2],[2,2]]]))), "TestSolidify")
end

function TestDiff()
	mypol1 = T([1,2])([-5,-5])(CUBOID([10,10]))
	mypol2 = S([1,2])([0.9,0.9])(mypol1)
	mypol3 = DIFF([mypol1,mypol2])
	VIEW(STRUCT([
			EX([0,10])(mypol3), T(1)(12),
			LEX([0,10])(mypol3), T(1)(25),
			S(3)(3)(SEX([0,PI])(16)(mypol3))
	]),"TestDiff")
end

function TestCube()
	VIEW(Cube(3))
end

function TestMapSphere()
	N,M = 16,16
	domain = Power(
		Translate(Quote([pi/N for I in 1:N]),[-pi/2]), 
		Quote([2*pi/M for I in 1:M]))
	fn=p -> [
		cos(p[1])*sin(p[2]), 
		cos(p[1])*cos(p[2]), 
		sin(p[1])
	]
	obj = MAP(fn)(domain)
	VIEW(obj)
end

function TestMkPol()
	out = MkPol(
		[[0.0],[1.0],[2.0],[3.0],[4.0],[5.0]],
		[[6,4],[1,2]])
	VIEW(out)
end



# ////////////////////////////////////////////////
if abspath(PROGRAM_FILE) == @__FILE__

	if false

		@assert C( v -> sum(v))(1)(2)==3
		@assert CAT([[1,2],[3,4]])==[1, 2, 3, 4]
		@assert INV([[1.0,0.0],[0.0,1.0]])==MatrixNd([[1.0, 0.0], [0.0, 1.0]])
		@assert EVERY(x -> x>=0,[1,2,3,4])==true
		@assert EVERY(x -> x>0,[1,-2,3,4])==false
		@assert AND([true,true])==true
		@assert AND([true,false])==false
		@assert ID(10)==10
		@assert DISTL([1,[2,3,4]])==[[1, 2], [1, 3], [1, 4]]
		@assert DISTR([[1,2,3],0])==[[1, 0], [2, 0], [3, 0]]
		@assert COMP([v -> [v;3], v -> [v;2], v -> [v;1]])([0])==[0,1,2,3]
		@assert AA(x -> x*2)([1,2,3])==[2, 4, 6]
		@assert EQ([1,1,1])==true
		@assert EQ([1,1,2])==false
		@assert NEQ([1,1,2])==true
		@assert NEQ([1,1,1])==false
		@assert FILTER(LE(0))([-1,0,1,2,3,4])==[-1, 0]
		@assert FILTER(GE(0))([-1,0,1,2,3,4])==[0, 1, 2, 3, 4]
		@assert APPLY([ x -> x*2, 2])==4
		@assert INSL(x -> x[1]-x[2])([1,2,3])==-4 # (1-2)-4
		@assert CONS([x -> x+1,x -> x+2])(0)==[1, 2]
		@assert IF([x -> x>0, K(10),K(20)])(+1)==10
		@assert IF([x -> x>0, K(10),K(20)])(-1)==20
		@assert LIFT(ADD)([cos,sin])(pi/2.0)==1.0
		@assert RAISE(ADD)([1,2])==3
		@assert RAISE(ADD)([cos,sin])(pi/2)==1.0
		@assert ISNUM(0.0)==true
		@assert ISFUN(x -> x) and ISFUN(abs)
		@assert !ISFUN(3) 
		@assert ISSEQOF(x -> isa(x,Int))([1,2,3  ])==true
		@assert ISSEQOF(x -> isa(x,Int))([1,2,3.0])==false
		@assert ISMAT([[1.0,2.0],[3.0,4.0]]) && ISMAT(MatrixNd(3))
		@assert ISMAT([1,2,3,4])==false
		@assert VECTSUM([[10,11,12],[0,1,2],[1,1,1]])==[11, 13, 15]
		@assert VECTDIFF([[10,11,12],[0,1,2],[1,1,1]])==[9, 9, 9]
		@assert MEANPOINT([[0,0,0],[1,1,1],[2,2,2]])==[1.0, 1.0, 1.0]
		@assert SUM([1,2,3])==6
		@assert SUM([[1,2,3],[4,5,6]])==[5, 7, 9]
		@assert SUM([[[1,2],[3,4]],[[10,20],[30,40]],[[100,200],[300,400]] ])==[[111, 222], [333, 444]]
		@assert DIFF(2)==-2
		@assert DIFF([1,2,3])==-4
		@assert DIFF([[1,2,3],[1,2,3]])==[0,0,0]
		@assert PROD([1,2,3,4])==24
		@assert PROD([[1,2,3],[4,5,6]])==32
		@assert DIV([10.0,2.0,5.0])==1
		@assert REVERSE([1,2,3])==[3, 2, 1]
		@assert REVERSE([1])==[1]
		@assert TRANS([[1.0,2.0],[3.0,4.0]])==[[1.0, 3.0], [2.0, 4.0]]
		@assert AR([[1,2,3],0])==[1, 2, 3, 0]
		@assert AL([0,[1,2,3]])==[0, 1, 2, 3]
		@assert PROGRESSIVESUM([1,2,3,4])==[1, 3, 6, 10]
		@assert INTSTO(5)==[1, 2, 3, 4, 5]
		@assert FROMTO([1,4])==[1, 2, 3, 4]
		@assert S1([1,2,3])==1
		@assert S2([1,2,3])==2
		@assert N(3)(10)==[10, 10, 10]
		@assert DIESIS(3)(10)==[10, 10, 10]
		@assert NN(3)([10])==[10, 10, 10]
		@assert DOUBLE_DIESIS(3)([10])==[10, 10, 10]
		@assert AS(SEL)([1,2,3])([10,11,12])==[10, 11, 12]
		@assert AC(SEL)([1,2,3])([10,11,[12,[13]]])==13
		#@assert CHARSEQ('hello')==['h', 'e', 'l', 'l', 'o']
		#@assert STRING(CHARSEQ('hello'))=='hello'
		@assert RANGE([1,3])==[1, 2, 3]
		@assert RANGE([3,1])==[3, 2, 1]
		@assert SIGN(10)==1
		@assert SIGN(-10)==-1
		@assert TREE(x -> x[1]>=x[end] ? x[1] : x[end])([1,2,3,4,3,2,1])==4
		@assert TREE(x -> x[1]>=x[end] ? x[1] : x[end])([1,2,3,4,3,2])==4
		@assert MERGE((x,y) -> x>y)([[1,3,4,5],[2,4,8]])==[1, 2, 3, 4, 4, 5, 8]
		@assert CASE([[LT(0),K(-1)],[C(EQ)(0),K(0)],[GT(0),K(+1)]])(-10)==-1
		@assert CASE([[LT(0),K(-1)],[C(EQ)(0),K(0)],[GT(0),K(+1)]])(0)==0
		@assert CASE([[LT(0),K(-1)],[C(EQ)(0),K(0)],[GT(0),K(+1)]])(10)==1
		@assert length(PERMUTATIONS([1,2,3]))==6
		@assert toList(IDNT(0))==[]
		@assert toList(IDNT(2))==[[1.0, 0.0], [0.0, 1.0]]
		@assert abs(SPLIT_2PI(4)[3]-PI)<1e-4
		@assert VECTPROD([[1,0,0],[0,1,0]])==[0.0, 0.0, 1.0]
		@assert VECTPROD([[0,1,0],[0,0,1]])==[1.0, 0.0, 0.0]
		@assert VECTPROD([[0,0,1],[1,0,0]])==[0.0, 1.0, 0.0]
		@assert VECTNORM([1,0,0])==1.0
		@assert INNERPROD([[1,2,3],[4,5,6]])==32.0
		@assert SCALARVECTPROD([2,[0,1,2]])==[0, 2, 4]
		@assert SCALARVECTPROD([[0,1,2],2])==[0, 2, 4]
		@assert MIXEDPROD([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])==1.0
		@assert UNITVECT([2,0,0])==[1.0, 0.0, 0.0]
		@assert UNITVECT([1,1,1])==UNITVECT([2,2,2])
		@assert DIRPROJECT([1.0,0.0,0.0])([2.0,0.0,0.0])==[2.0, 0.0, 0.0]
		@assert DIRPROJECT([1.0,0.0,0.0])([0.0,1.0,0.0])==[0.0, 0.0, 0.0]
		@assert ORTHOPROJECT([1.0,0.0,0.0])([1.0,1.0,0.0])==[0.0, 1.0, 0.0]
		@assert FACT(4)==24
		@assert FACT(0)==	1
		@assert CHOOSE([7,3])==35.0
		@assert TRACE(MatrixNd([[5.0,0.0],[0.0,10.0]]))==15.0
		@assert toList(MATHOM(MatrixNd( [[1,2],[3,4]])))==[[1.0, 0.0, 0.0], [0.0, 1.0, 2.0], [0.0, 3.0, 4.0]]
		@assert SCALARMATPROD([10.0,[[1,2],[3,4]]])==[[10.0, 20.0], [30.0, 40.0]]
		@assert ORTHO([[1,0],[0,1]])==[[1.0, 0.0], [0.0, 1.0]]
		@assert SKEW([[1,0],[0,1]])==[[0.0, 0.0], [0.0, 0.0]]
		@assert length(CART([ [1, 2, 3], ['a', 'b'],[10,11] ]))==12
		@assert length(POWERSET([1, 2, 3]))==8
		@assert ISPOL(Cube(2))==true
		@assert box(GRID([1,-1,1])) == BoxNd([0.0], [3.0])
		@assert box(GRID([-1,1,-1,1]))==BoxNd([1.0], [4.0])
		@assert box(INTERVALS(10.0)(8))==BoxNd([0.0], [10.0])
		@assert box(CUBOID([1,2,3]))==BoxNd([0.0, 0.0, 0.0], [1.0, 2.0, 3.0])
		@assert box(SIMPLEX(3))==BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
		@assert RN(Cube(2))==2
		@assert DIM(Cube(2))==2
		@assert box(MKPOL([[0.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,1.0]], [[1,2,3,4]] ))==BoxNd([0.0, 0.0], [1.0, 1.0])
		@assert box((TRANSLATE(3)(2)(Cube(2))))== BoxNd([0.0, 0.0, 2.0], [1.0, 1.0, 2.0])
		@assert box((TRANSLATE([1,3])([1,2])(Cube(2))))== BoxNd([1.0, 0.0, 2.0], [2.0, 1.0, 2.0])
		@assert box((SCALE(3)(2)(Cube(3))))==BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 2.0])
		@assert box((SCALE([3,1])([4,2])(Cube(3))))==BoxNd([0.0, 0.0, 0.0], [2.0, 1.0, 4.0])
		@assert fuzzyEqual(box((ROTATE([1,2])(PI/2)(Cube(2)))),BoxNd([-1.0,0.0],[0.0,+1.0]))
		@assert box((MAT([[1,0,0],[1,1,0],[2,0,1]])(Cube(2))))==BoxNd([1.0, 2.0], [2.0, 3.0])

		@assert box((STRUCT([
			Cube(2),
			T([1,2])([1.0,1.0]),
			T([1,2])([1.0,1.0]),
			Cube(2),
			Cube(2,1.0,2.0)])))==BoxNd([0.0, 0.0], [4.0, 4.0])

		@assert box((STRUCT([
			T([1,2])([1,1]),
			T([1,2])([1,1]),
			Cube(2),
			T([1,2])([1,1]),
			T([1,2])([1,1]),
			Cube(2),
			Cube(2,1.0,2.0)])))==BoxNd([2.0, 2.0], [6.0, 6.0])

		@assert box(JOIN([Cube(2,0.0,1.0)]))==BoxNd([0.0, 0.0], [1.0, 1.0])
		@assert POWER([2,2])==4.0
		@assert box((POWER([Cube(2),Cube(1)])))==BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
		@assert SIZE(1)(Cube(2))==1.0
		@assert SIZE([1,3])(SCALE([1,2,3])([1,2,3])(Cube(3)))==[1.0, 3.0]
		@assert MIN(1)(Cube(2))==0.0
		@assert MIN([1,3])(TRANSLATE([1,2,3])([10.0,20.0,30.0])(Cube(3)))==[10.0, 30.0]
		@assert MAX(1)(Cube(2))==1.0
		@assert MAX([1,3])(TRANSLATE([1,2,3])([10.0,20.0,30.0])(Cube(3)))==[11.0, 31.0]
		@assert MED(1)(Cube(2))==0.5
		@assert MED([1,3])(Cube(3))==[0.5, 0.5]
		@assert box((ALIGN([3,MAX,MIN])([Cube(3),Cube(3)])))==BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 2.0])
		@assert box((TOP([Cube(3),Cube(3)])))==BoxNd([0.0, 0.0, 0.0], [1.0, 1.0, 2.0])
		@assert box((BOTTOM([Cube(3),Cube(3)])))==BoxNd([0.0, 0.0, -1.0], [1.0, 1.0, 1.0])
		@assert box((LEFT([Cube(3),Cube(3)])))==BoxNd([-1.0, 0.0, 0.0], [1.0, 1.0, 1.0])
		@assert box((RIGHT([Cube(3),Cube(3)])))==BoxNd([0.0, 0.0, 0.0], [2.0, 1.0, 1.0])
		@assert box((UP([Cube(3,0.0,1.0),Cube(3,5.0,6.0)])))==BoxNd([0.0, 0.0, 0.0], [6.0, 2.0, 1.0])
		@assert box((DOWN([Cube(3,0.0,1.0),Cube(3,5.0,6.0)])))==BoxNd([0.0, -1.0, 0.0], [6.0, 1.0, 1.0])

		obj=Translate(Cube(3),[1.0,2.0,3.0])
		@assert box(BOX([1,3])(obj))==BoxNd([1.0, 3.0], [2.0, 4.0])

		obj=Translate(Cube(3),[1.0,2.0,3.0])
		@assert box(BOX(3)(obj))==BoxNd([3.0], [4.0])

		@assert box((MAP([S1,S2])(Cube(2))))==BoxNd([0.0, 0.0], [1.0, 1.0])
		@assert box((MAP(ID)(Cube(2))))==BoxNd([0.0, 0.0], [1.0, 1.0])
		@assert box(CIRCUMFERENCE(1)(8))==BoxNd([-1.0, -1.0], [1.0, 1.0])
		@assert box((RING([0.5,1])([8,8])))==BoxNd([-1.0, -1.0], [1.0, 1.0])
		@assert box( CIRCLE(1.0)([8,8]))==BoxNd([-1.0, -1.0], [1.0, 1.0])
		@assert fuzzyEqual(box(CYLINDER([1.0,2.0])(8)),BoxNd([-1.0,-1.0,0.0],[+1.0,+1.0,2.0]))
		@assert box((SPHERE(1.0)([8,8])))==BoxNd([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0])
		@assert fuzzyEqual(box((TORUS([1.0,2.0])([8,8]))),BoxNd([-2.0,-2.0,-0.5],[+2.0,+2.0,+0.5]))
		@assert fuzzyEqual(box((CONE([1.0,3.0])(16))),BoxNd([-1.0,-1.0,0.0],[+1.0,+1.0,3.0]))

		# BROKEN
		if false
			@assert fuzzyEqual(box(UNION([
				Cube(2,0.0,1.0),
				Cube(2,0.5,1.5)])),BoxNd([0.0,0.0],[1.5,1.5]))
			@assert fuzzyEqual(box(INTERSECTION([
				Cube(2,0.0,1.0),
				Cube(2,0.5,1.5)])),BoxNd([0.5,0.5],[1.0,1.0]))
			@assert fuzzyEqual(box(DIFFERENCE([
				Cube(2,0.0,1.0),
				Cube(2,0.5,1.5)])),BoxNd([0.0,0.0],[1.0,1.0]))
			@assert fuzzyEqual(box(XOR([
				Cube(2,0,1),
				Cube(2,0.5,1.5)])),BoxNd([0.0,0.0],[1.5,1.5]))
		end

	else

		# ok
		if false
			#TestCube()
			#TestMkPol()
			#TestSphere()
			#TestMapSphere()
			#TestTorus()
			#TestBezier()
			TestCoonsPatch()
			TestRuledSurface()
			TestProfileProdSurface()
			TestRotationalSurface()
			TestConicalSurface()
			TestCubicHermite()
			TestSchegel3d()
			TestHermiteSurface()
			TestPermutahedron() 
			TestCubicSpline() 
			TestBilinarSurface()
			TestBiquadraticSurface()
			TestBezierSurface()
		end

		# TestOffset()
		#TestThinSolid()
		#TestEllipse()
		#TestBezierStripe()
		#TestDisplayNubSpline()
		#TestDisplayNurbsSpline()
		#TestMinkowski()

		# BROKEN
		# TestCylindricalSurface() 
		# TestBezierManifold()

	end
end

