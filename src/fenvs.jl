PI = pi
COS = cos
LEN = length
AND = all
OR = any

include("viewer.jl")
include("hpc.jl")

# /////////////////////////////////////////////////////////////////
"""
```jldoctest
julia> C(x,y -> x+y)(1)(2)
3
```
"""
function C(fun)
	return function(arg1)
		return function(arg2)
			return fun([arg1,arg2])
		end
	end
end

# /////////////////////////////////////////////////////////////////
ATAN2(l) = atan2(l[2],l[1])

# /////////////////////////////////////////////////////////////////
MOD(l) = float(l[1] % l[2])

# /////////////////////////////////////////////////////////////////
function CAT(args)
	return reduce(vcat, args)
end

# /////////////////////////////////////////////////////////////////
function INV(mat)
	return inv(mat)
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
	return apply(f,[x])
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
	 return reduce((x,y) -> x/float(y), args)
end

REVERSE(List) = reverse(List)

# /////////////////////////////////////////////////////////////////
function TRANS(List)
	if isa(List, Matrix)
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
			return [[float(jt) for jt in it] for it in ret]
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
	 return [args[1] + abs(args[2])]
end

function AL(args)
	 return [args[1]] + args[2]
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
	return function(List)
		return List^n
	end
end

# /////////////////////////////////////////////////////////////////
function DOUBLE_DIESIS(n)
	return function(List)
		return List^n
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
	function TREE_NO_CURRIED(fun, List)
		length = length(List)
		if length == 1
			return List[1]
		end
		k = div(length, 2)
		return f([TREE_NO_CURRIED(f, List[1:k])] + [TREE_NO_CURRIED(f, List[k+1:end])])
	end
	return function(x)
		return TREE_NO_CURRIED(f,x)
	end
end

# /////////////////////////////////////////////////////////////////
function MERGE(f)
	function MERGE_NO_CURRIED(f, List)
		list_a, list_b = List
		if length(list_a) == 0
			return list_b
		end
		if length(list_b) == 0
			return list_a
		end
		res = f(list_a[1], list_b[1])
		if !res
			return [list_a[1]] + MERGE_NO_CURRIED(f,[list_a[2:end], list_b])
		else
			return [list_b[1]] + MERGE_NO_CURRIED(f,[list_a, list_b[2:end]])
		end
	end
	return function(x)
		return MERGE_NO_CURRIED(f,x)
	end
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
		rest  =PERMUTATIONS(SEQ[1:i-1] + SEQ[i+1:end])
		for r in rest
			push!(ret, [element; r])
		end
	end
	return ret
end

# /////////////////////////////////////////////////////////////////
function IDNT(N)
	return Matrix(N)
end

# /////////////////////////////////////////////////////////////////
function SPLIT_2PI(N)
	delta=2*PI/N
	return [i*delta for i in 0:N-1]
end

# /////////////////////////////////////////////////////////////////
function VECTPROD(args)
	a,b=args
	ax,ay,az=[float(it) for it in a]
	bx,by,bz=[float(it) for it in b]
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
	dim=length(MATRIX)
	for i in 1:dim
		acc+=MATRIX[i,i]
	end
	return acc
end

# /////////////////////////////////////////////////////////////////
function PASCALTRIANGLE(N)
	if N==0
		return [[1]]
	end
	if N==1
		return [[1],[1,1]]
	end
	prev=PASCALTRIANGLE(N-1)
	last_row=prev[end]
	cur=[1]+[last_row[i-1]+last_row[i]  for i in 2:length(last_row)] + [1]
	return prev+[cur]
end

# /////////////////////////////////////////////////////////////////
function MATHOM(M)
	dim=length(M)
	ret=Matrix(dim+1)
	for R in 1:dim
		for C in 1:dim
			ret[R+1,C+1]=M[R,C]
		end
	end
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
	CART2 = COMP([COMP([CAT, AA(DISTL)]), DISTR])
	F1 = AA((AA(CONS([ID]))))
	return TREE(COMP([AA(CAT), CART2]))(F1(l))
end

# /////////////////////////////////////////////////////////////////
#function POWERSET(l)
#	return COMP([COMP([AA(CAT), CART]), AA((CONS([CONS([ID]), K([]))]))])(l)
#end

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
function VIEW(obj,title="")
	 obj.view(title=title)
end

# /////////////////////////////////////////////////////////////////
function GRID(sequence)
	cursor,points,hulls=(0,[[0]],[])
	for value in sequence
		points = points + [[cursor + abs(value)]]
		if value>=0
			hulls += [[length(points)-2,length(points)-1]]
		end
		cursor = cursor + abs(value)
	end
	return  MkPol(points, hulls)
end
QUOTE = GRID

# Q = COMP([QUOTE, IF([ISSEQ, ID, CONS([ID])])])

# /////////////////////////////////////////////////////////////////
function INTERVALS(A)
	function INTERVALS0(N)
		return QUOTE([float(A)/float(N) for i in 1:N])
	end
	return INTERVALS0
end

# /////////////////////////////////////////////////////////////////
function CUBOID(vs)
	return Scale(Cube(length(vs)),vs)
end

function CUBE(size)
	return CUBOID([size, size, size])
end


# /////////////////////////////////////////////////////////////////
function HEXAHEDRON()
	return Cube(3, -1.0/sqrt(3.0), +1.0/sqrt(3.0))
end

# /////////////////////////////////////////////////////////////////
function SIMPLEX(dim)
	return Simplex(dim)
end

RN(pol) = pol.dim()

DIM(pol) = pol.dim()

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
function MKPOL(args_list)
	points,hulls,__pols=args_list
	return MKPOL(points, hulls, __pols)
end

# /////////////////////////////////////////////////////////////////
MK = COMP([MKPOL, CONS([LIST, K([[1]]), K([[1]])])])

# /////////////////////////////////////////////////////////////////
function CONVEXHULL(points)
	return MKPOL(points, [collect(1:length(points))], [[1]])
end

# /////////////////////////////////////////////////////////////////
function UKPOL(pol)
	points, hulls = pol.ukpol()
	hulls = [[idx+1 for idx in hull] for hull in hulls]
	return [points, hulls, [[1]]]
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
			vt = [0.0] * maximum(axis)
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
	function SCALE1(axis, values)
		function SCALE2(axis, values, pol)
			axis = ISNUM(axis) ? [axis] : axis
			values = ISNUM(values) ? [values] : values
			vs = [1.0] * maximum(axis)
			for (a, t) in zip(axis, values)
					vs[a] = t
			end
			return pol.scale(vs)
		end
		return SCALE2
	end
	return SCALE1(axis)
end
S = SCALE

# /////////////////////////////////////////////////////////////////
function ROTATE(plane_indexes)
	function ROTATE1(angle)
		function ROTATE2(pol)
			dim = maximum(plane_indexes)
			return pol.rotate(plane_indexes[1], plane_indexes[2], angle)
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
		return pol.transform(MatrixNd(matrix))
	end
	return MAT0
end

# /////////////////////////////////////////////////////////////////
function EMBED(up_dim)
	function EMBED1(pol)
		new_dim_pol = pol.dim() + up_dim
		return Hpc(MatrixNd(new_dim_pol+1), [pol])
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
	pols = []
	while !isempty(seq) && ISPOL(seq[1])
		push!(pols, seq[1])
		seq = seq[2:end]
	end
	transformations = []
	while !isempty(seq) && ISFUN(seq[1])
		push!(transformations, seq[1])
		seq = seq[2:end]
	end
	if !isempty(seq) && !ISPOL(seq[1]) && !ISFUN(seq[1])
		error("STRUCT arguments not valid, not all elements are pols or transformations")
	end
	if !isempty(seq)
		assert(ISPOL(seq[1]))
		child = STRUCT(seq, nrec+1)
		assert(ISPOL(child))
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
	
	function BspFace(vertices, plane)
		 if isempty(plane)
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
			  plane = [normal..., -w]
		 end
		 new(vertices, plane)
	end
	
	function splitEdge(plane, vi, vj)
		 dim = length(vi)
		 valuev1 = abs(plane[end] + sum([plane[i]*vi[i] for i in 1:dim]))
		 valuev2 = abs(plane[end] + sum([plane[i]*vj[i] for i in 1:dim]))
		 alpha = 1.0 / (valuev1 + valuev2)
		 beta = valuev1 * alpha
		 alpha = alpha * valuev2
		 return [alpha*vi[i] + beta*vj[i] for i in 1:dim]
	end
	
	function split(self, plane, EPSILON=1e-5)
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
end


# //////////////////////////////////////////////////////////////////////////////////////
struct Bsp
	 
	 
	 plane::Union{Nothing, Vector{Float64}}
	 faces::Vector{BspFace}
	 below::Union{Nothing, Bsp}
	 above::Union{Nothing, Bsp}
	 
end


function getFaces(self::Bsp)
	return self.faces
		 + (self.below !== nothing ? getFaces(self.below) : [])
		 + (self.above !== nothing ? getFaces(self.above) : [])
end

function insertFaces(self::Bsp, faces)
	if isempty(faces)
		 return self
	end
	if self.plane === nothing
		 @assert self.below === nothing && self.above === nothing
		 self.plane = faces[1].plane
	end
	below, above = [], []
	for p in faces
		 b, cb, ca, a = split(p, self.plane)
		 if b !== nothing
			  push!(below, b)
		 end
		 if cb !== nothing
			  push!(self.faces, cb)
		 end
		 if ca !== nothing
			  push!(self.faces, ca)
		 end
		 if a !== nothing
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
		 b, cb, ca, a = split(p, self.plane)
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
	for (T, properties, obj) in toBoundaryForm(hpc).toList()
		 faces += [BspFace([T.transformPoint(obj.points[I]) for I in hull]) for hull in obj.hulls]
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

UNION(objs) = begin
	 objs = [Bsp.fromHpc(obj) for obj in objs]
	 res = objs[1]
	 for I in 2:length(objs)
		  res = Bsp.Union(res, objs[I])
	 end
	 return toHpc(res)
end

INTERSECTION(objs) = begin
	 objs = [Bsp.fromHpc(obj) for obj in objs]
	 res = objs[1]
	 for I in 2:length(objs)
		  res = Bsp.Intersection(res, objs[I])
	 end
	 return toHpc(res)
end

DIFFERENCE(objs) = begin
	 objs = [Bsp.fromHpc(obj) for obj in objs]
	 res = objs[1]
	 for I in 2:length(objs)
		  res = Bsp.Difference(res, objs[I])
	 end
	 return toHpc(res)
end

XOR(objs) = begin
	 objs = [Bsp.fromHpc(obj) for obj in objs]
	 res = objs[1]
	 for I in 2:length(objs)
		  res = Bsp.Xor(res, objs[I])
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
function SIZE(List)
	 function SIZE1(pol)
		  size = box(pol).size()
		  return isa(List, Vector) ? [size[i] for i in List] : size[List]
	 end
	 return SIZE1
end

# ///////////////////////////////////////////////////////////
function MIN(List)
	 function MIN1(pol)
		  box = box(pol)
		  return isa(List, Vector) ? [box.p1[i] for i in List] : box.p1[List]
	 end
	 return MIN1
end

# ///////////////////////////////////////////////////////////
function MAX(List)
	 function MAX1(pol)
		  box = box(pol)
		  return isa(List, Vector) ? [box.p2[i] for i in List] : box.p2[List]
	 end
	 return MAX1
end

# ///////////////////////////////////////////////////////////
function MED(List)
	 function MED1(pol)
		  center = box(pol).center()
		  return isa(List, Vector) ? [center[i] for i in List] : center[List]
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
		  max_index = maximum([index for (index, pos1, po2) in args])
		  vt = zeros(max_index)
		  for (index, pos1, pos2) in args
				p1 = ifelse(pos1 == MIN, box1.p1, ifelse(pos1 == MAX, box1.p2, box1.center()))
				p1 = ifelse(index <= length(p1), p1[index], 0.0)
				p2 = ifelse(pos2 == MIN, box2.p1, ifelse(pos2 == MAX, box2.p2, box2.center()))
				p2 = ifelse(index <= length(p2), p2[index], 0.0)
				vt[index] -= (p2 - p1)
		  end
		  return Struct([pol1, pol2.translate(vt)])
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
function BOX(List)
	 function BOX0(List, pol)
		  if !isa(List, Vector)
				List = [List]
		  end
		  dim = length(List)
		  box = box(pol)
		  vt = [box.p1[i] for i in List]
		  vs = [box.size()[i] for i in List]
		  return Cube(dim).scale(vs).translate(vt)
	 end
	 return pol -> BOX0(List, pol)
end

# ///////////////////////////////////////////////////////////
function MAP(fn)
	function MAP0(fn, pol)
		if isa(fn, Tuple) || isa(fn, Vector)
			map_fn = point -> [f(point) for f in fn]
		else
			map_fn = fn
		end
		return pol.mapFn(map_fn)
	end
	return pol -> MAP0(fn, pol)
end



# /////////////////////////////////////////////////////////////////
function CIRCLE_POINTS(R, N)
	return [[R*cos(i*2*pi/N), R*sin(i*2*pi/N)] for i in 0:N-1]
end

# ///////////////////////////////////////////////////////////
function CIRCUMFERENCE(R)
	return N -> begin
		MAP(p -> [R*cos(p[1]), R*sin(p[1])], INTERVALS(2*pi)(N))
	end
end

# ///////////////////////////////////////////////////////////
function NGON(N)
	return CIRCUMFERENCE(1)(N)
end

# ///////////////////////////////////////////////////////////
function RING(radius)
	R1, R2 = radius
	function RING0(subds)
		N, M = subds
		domain = POWER([INTERVALS(2*pi)(N), INTERVALS(R2-R1)(M)]).translate([0.0, R1])
		fun = p -> [p[2]*cos(p[1]), p[2]*sin(p[1])]
		return MAP(fun, domain)
	end
	return RING0
end

# ///////////////////////////////////////////////////////////
function TUBE(args)
	r1, r2, height = args
	function TUBE0(N)
		return Power(RING([r1, r2])([N, 1]), QUOTE([height]))
	end
	return TUBE0
end

# ///////////////////////////////////////////////////////////
function CIRCLE(R)
	function CIRCLE0(subs)
		N, M = subs
		domain = POWER([INTERVALS(2*pi)(N), INTERVALS(R)(M)])
		fun = p -> [p[2]*cos(p[1]), p[2]*sin(p[1])]
		return MAP(fun, domain)
	end
	return CIRCLE0
end

# ///////////////////////////////////////////////////////////
function MY_CYLINDER(args)
	R, H = args
	function MY_CYLINDER0(N)
		points = CIRCLE_POINTS(R, N)
		circle = Mkpol(points, [collect(1:N)])
		return Power(circle, MkPol([[0], [H]], [[0, 1]]))
	end
	return MY_CYLINDER0
end
CYLINDER = MY_CYLINDER

# /////////////////////////////////////////////////////////////
function SPHERE(radius)
	function SPHERE0(subds)
		N, M = subds
		domain = Power(INTERVALS(pi)(N), INTERVALS(2*pi)(M)).translate([-pi/2, 0])
		fx = p -> radius * cos(p[1]) * sin(p[2])
		fy = p -> radius * cos(p[1]) * cos(p[2])
		fz = p -> radius * sin(p[1])
		ret = MAP([fx, fy, fz], domain)
		return ret
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
		return MAP([fx, fy, fz], domain)
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
		return MAP(fn, domain)
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
					weight = CHOOSE([N, I])*pow(1-t, N-I)*pow(t, I)
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
	vertices = PERMUTATIONS(collect(1:d+1))
	center = MEANPOINT(vertices)
	cells = [collect(1:length(vertices))]
	object = MKPOL(vertices, cells, [[1]])
	object = object.translate([-coord for coord in center])
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
		ret = []
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
					ret = CHOOSE([N, I])*pow(1-t, N-I)*pow(t, I)
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
		function map_fn(point)
			u, v = point
			U = [f([u]) for f in ubasis]
			V = [f([v]) for f in vbasis]
			controlpoints = [ isa(f, Function) ? f(point) : f for f in controlpoints_fn]
			target_dim = length(controlpoints[1][1])
			ret = [0 for x in 1:target_dim]
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
	u0(point) = u = S1(point); 2*u*u-u
	u1(point) = u = S1(point); 4*u-4*u*u
	u2(point) = u = S1(point); 2*u*u-3*u+1
	basis = [u0, u1, u2]
	return TENSORPRODSURFACE([basis, basis])(controlpoints)
end

# //////////////////////////////////////////////////////////////
function HERMITESURFACE(controlpoints)
	H0(point) = u = S1(point); u2 = u*u; u3 = u2*u; u3-u2
	H1(point) = u = S1(point); u2 = u*u; u3 = u2*u; u3-2*u2+u
	H2(point) = u = S1(point); u2 = u*u; u3 = u2*u; 3*u2-2*u3
	H3(point) = u = S1(point); u2 = u*u; u3 = u2*u; 2*u3-3*u2+1
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

function CURVE2MAPVECT(CURVE)
	D = length((CURVE([0])))
	return [ COMP([SEL(i),CURVE]) for i in FROMTO([1,D]) ]
end

function SEGMENT(sx)
	function SEGMENT0(args)
		N = length(args[1])
		A, B = args
		P0 = A
		P1 = [A[i]+(B[i]-A[i])*sx for i in range(N)]
		return POLYLINE([P0,P1])
	end
	return SEGMENT0
end

# //////////////////////////////////////////////////////////////
function SOLIDIFY(pol)
	box = pol.box()
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
	 for i in range(length(faces))
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
		for i in range(length(v))
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
		for i in range(length(v))
			shear = [j!=i ? 0 : v[i] for j in range(length(v))] + [0 for j in range(i)]
			mat = IDNT(length(shear)+2)
			for i in range(length(shear))
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
		GX = [PX[i]-P0[i] for i in range(3)]
		GY = [PY[i]-P0[i] for i in range(3)]
		normal = UNITVECT(VECTPROD([GX, GY]))
		ret = [P0[i]+w*normal[i] for i in range(3)]
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
	degree = length(controlpoints_fn)-1
	basis = BERNSTEINBASIS(S1)(degree)
	function map_fn(point)
		controlpoints = [ isa(f, Function) ? f(point) : f for f in controlpoints_fn]
		target_dim = length(controlpoints[1])
		ret = [0 for i in range(target_dim)]
		for i in range(length(basis))
			coeff = basis[i](point)
			for M in range(target_dim)
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

	basis = [DERBERNSTEIN(degree)(i) for i in range(degree+1)]
	function map_fn(point)
		controlpoints = [isa(f, Function)  ? f(point) : f for f in controlpoints_fn]
		target_dim = length(controlpoints[1])
		ret = [0 for i in range(target_dim)]
		for i in range(length(basis))
			coeff = basis[i](point)
			for M in range(target_dim)
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
	domain = S(2)(width)(T(1)(0.00001)(Power(INTERVALS(1)(n),INTERVALS(1)(1))))
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
					ret = [0 for i in range(target_dim)]
					for i in range(n+1)
						coeff = N(i, k, t) 
						for M in range(target_dim)
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
			tmin = minimum(knots)
			tmax = maximum(knots)
			tsiz = tmax-tmin
			v = [tsiz/float(totpoints-1) for i in range(totpoints-1)]
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
			tmin = minimum(knots)
			tmax = maximum(knots)
			tsiz = tmax-tmin
			v = [tsiz/float(totpoints-1) for i in range(totpoints-1)]
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
	VIEW(SPHERE(1)([16,16]), title="TestSphere")
end

function TestTorus()
	VIEW(TORUS([1,2])([20,20]), title="TestTorus")
end

function TestBezier()
	VIEW(MAP(BEZIER(S1)([[-0,0],[1,0],[1,1],[2,1],[3,1]]))(INTERVALS(1)(32)), title="TestBezier-1")
	C0 = BEZIER(S1)([[0,0,0],[10,0,0]])
	C1 = BEZIER(S1)([[0,2,0],[8,3,0],[9,2,0]])
	C2 = BEZIER(S1)([[0,4,1],[7,5,-1],[8,5,1],[12,4,0]])
	C3 = BEZIER(S1)([[0,6,0],[9,6,3],[10,6,-1]])
	VIEW(MAP(BEZIER(S2)([C0,C1,C2,C3]))(Power(INTERVALS(1)(10),INTERVALS(1)(10))), title="TestBezier-2")
end

function TestCoonsPatch()
	Su0 = BEZIER(S1)([[0,0,0],[10,0,0]])
	Su1 = BEZIER(S1)([[0,10,0],[2.5,10,3],[5,10,-3],[7.5,10,3],[10,10,0]])
	Sv0 = BEZIER(S2)([[0,0,0],[0,0,3],[0,10,3],[0,10,0]])
	Sv1 = BEZIER(S2)([[10,0,0],[10,5,3],[10,10,0]])
	VIEW(MAP(COONSPATCH([Su0,Su1,Sv0,Sv1]))(Power(INTERVALS(1)(10),INTERVALS(1)(10))), title="TestCoonsPatch")
end

function TestRuledSurface()
	alpha = point -> [point[0],point[0],       0 ]
	beta = point -> [      -1,      +1,point[0] ]
	domain = T([1,2])([-1,-1])(Power(INTERVALS(2)(10),INTERVALS(2)(10)))
	VIEW(MAP(RULEDSURFACE([alpha,beta]))(domain), title="TestRuledSurface")
end

function TestProfileProdSurface()
	alpha = BEZIER(S1)([[0.1,0,0],[2,0,0],[0,0,4],[1,0,5]])
	beta = BEZIER(S2)([[0,0,0],[3,-0.5,0],[3,3.5,0],[0,3,0]])
	domain = Power(INTERVALS(1)(20),INTERVALS(1)(20))
	VIEW(Struct([MAP(alpha)(domain),MAP(beta )(domain),MAP(PROFILEPRODSURFACE([alpha,beta]))(domain)]), title="TestProfileProdSurface")
end

function TestRotationalSurface()
	profile = BEZIER(S1)([[0,0,0],[2,0,1],[3,0,4]]) 
	domain = Power(INTERVALS(1)(10),INTERVALS(2*PI)(30)) 
	VIEW(MAP(ROTATIONALSURFACE(profile))(domain), title="TestRotationalSurface")
end

function TestCylindricalSurface()
	alpha = BEZIER(S1)([[1,1,0],[-1,1,0],[1,-1,0],[-1,-1,0]])
	Udomain = INTERVALS(1)(20)
	Vdomain = INTERVALS(1)(6)
	domain = Power(Udomain,Vdomain)
	fn = CYLINDRICALSURFACE([alpha,[0,0,1]])
	VIEW(MAP(fn)(domain), title="TestCylindricalSurface")
end

function TestConicalSurface()
	domain = Power(INTERVALS(1)(20),INTERVALS(1)(6))
	beta = BEZIER(S1)([ [1,1,0],[-1,1,0],[1,-1,0],[-1,-1,0] ])
	VIEW(MAP(CONICALSURFACE([[0,0,1],beta]))(domain), title="TestConicalSurface")
end

function TestCubicHermite()
	domain = INTERVALS(1)(20)
	VIEW(Struct([
		MAP(CUBICHERMITE(S1)([[1,0],[1,1],[ -1, 1],[ 1,0]]))(domain),
		MAP(CUBICHERMITE(S1)([[1,0],[1,1],[ -2, 2],[ 2,0]]))(domain),
		MAP(CUBICHERMITE(S1)([[1,0],[1,1],[ -4, 4],[ 4,0]]))(domain),
		MAP(CUBICHERMITE(S1)([[1,0],[1,1],[-10,10],[10,0]]))(domain)])
		, title="TestCubicHermite-1")
	c1 = CUBICHERMITE(S1)([[1  ,0,0],[0  ,1,0],[0,3,0],[-3,0,0]])
	c2 = CUBICHERMITE(S1)([[0.5,0,0],[0,0.5,0],[0,1,0],[-1,0,0]])
	sur3 = CUBICHERMITE(S2)([c1,c2,[1,1,1],[-1,-1,-1]])
	domain = Power(INTERVALS(1)(14),INTERVALS(1)(14))
	VIEW(MAP(sur3)(domain), title="TestCubicHermite-2")
end

function TestPermutahedron()
	VIEW(Struct([PERMUTAHEDRON(2),(PERMUTAHEDRON(2))]), title="TestPermutahedron-1")
	VIEW(Struct([PERMUTAHEDRON(3),(PERMUTAHEDRON(3))]), title="TestPermutahedron-2")
end

function TestSchegel3d()
	VIEW(SCHLEGEL3D(0.2)((T([1,2,3,4])([-1.0/3.0,-1.0/3.0,-1,+1])(SIMPLEX(4)))), title="TestSchegel3d-1")
	VIEW(SCHLEGEL3D(0.2)((T([1,2,3,4])([-1,-1,-1,1])(CUBOID([2,2,2,2])))), title="TestSchegel3d-2")
	VIEW(SCHLEGEL3D(0.2)((T([1,2,3,4])([-1.0/3.0,-1.0/3.0,-1,+1])(Power(SIMPLEX(2),SIMPLEX(2))))), title="TestSchegel3d-3")
end

function TestCubicSpline()
	domain = INTERVALS(1)(20)
	points = [[-3,6],[-4,2],[-3,-1],[-1,1],[1.5,1.5],[3,4],[5,5],[7,2],[6,-2],[2,-3]]
	VIEW(SPLINE(CUBICCARDINAL(domain))(points), title="TestCubicSpline-1")
	VIEW(SPLINE(CUBICUBSPLINE(domain))(points), title="TestCubicSpline-2")
end

function TestBilinarSurface()
	controlpoints = [[[0,0,0],[2,-4,2]],[[0,3,1],[4,0,0]]]
	domain = Power(INTERVALS(1)(10),INTERVALS(1)(10))
	mapping = BILINEARSURFACE(controlpoints)
	VIEW(MAP(mapping)(domain), title="TestBilinarSurface")
end

function TestBiquadraticSurface()
	controlpoints = [[[0,0,0],[2,0,1],[3,1,1]],[[1,3,-1],[3,2,0],[4,2,0]],[[0,9,0],[2,5,1],[3,3,2]]]
	domain = Power(INTERVALS(1)(10),INTERVALS(1)(10))
	mapping = BIQUADRATICSURFACE(controlpoints)
	VIEW(MAP(mapping)(domain), title="TestBiquadraticSurface")
end

function TestHermiteSurface()
	controlpoints = [[[0,0,0 ],[2,0,1],[3,1,1],[4,1,1]],[[1,3,-1],[3,2,0],[4,2,0],[4,2,0]],[[0,4,0 ],[2,4,1],[3,3,2],[5,3,2]],[[0,6,0 ],[2,5,1],[3,4,1],[4,4,0]]]
	domain = Power(INTERVALS(1)(10),INTERVALS(1)(10))
	mapping = HERMITESURFACE(controlpoints)
	VIEW(MAP(mapping)(domain), title="TestHermiteSurface")
end

function TestBezierSurface()
	controlpoints = [[[ 0,0,0],[0 ,3  ,4],[0,6,3],[0,10,0]],[[ 3,0,2],[2 ,2.5,5],[3,6,5],[4,8,2]],[[ 6,0,2],[8 ,3 , 5],[7,6,4.5],[6,10,2.5]],[[10,0,0],[11,3  ,4],[11,6,3],[10,9,0]]]
	domain = Power(INTERVALS(1)(10),INTERVALS(1)(10))
	mapping = BEZIERSURFACE(controlpoints)
	VIEW(MAP(mapping)(domain), title="TestBezierSurface")
end

function TestBezierManifold()
	grid1D = INTERVALS(1)(5)
	domain3D = Power(Power(grid1D,grid1D),grid1D)
	degrees = [2,2,2]
	Xtensor =  [[[0,1,2],[-1,0,1],[0,1,2]],[[0,1,2],[-1,0,1],[0,1,2]],[[0,1,2],[-1,0,1],[0,1,2]]]
	Ytensor =  [[[0,0,0.8],[1,1,1],[2,3,2]],[[0,0,0.8],[1,1,1],[2,3,2]],[[0,0,0.8],[1,1,1],[2,3,2]]]
	Ztensor =  [[[0,0,0],[0,0,0],[0,0,0]],[[1,1,1],[1,1,1],[1,1,1]],[[2,2,1],[2,2,1],[2,2,1]]] 
	mapping = BEZIERMANIFOLD(degrees)([Xtensor,Ytensor,Ztensor])
	VIEW(MAP(mapping)(domain3D), title="TestBezierManifold")
end

function TestOffset()
	verts = [[0,0,0],[3,0,0],[3,2,0],[0,2,0],[0,0,1.5],[3,0,1.5],[3,2,1.5],[0,2,1.5],[0,1,2.2],[3,1,2.2]]
	cells = [[1,2],[2,3],[3,4],[4,1],[5,6],[6,7],[7,8],[8,5],[1,5],[2,6],[3,7],[4,8],[5,9],[8,9],[6,10],[7,10], [9,10]]
	pols = [[1]]
	House = MKPOL(verts,cells,pols)
	VIEW(Struct([OFFSET([0.1,0.2,0.1])(House), T(1)(1.2*SIZE(1)(House))(House)]), title="TestOffset")
end

function TestThinSolid()
	Su0 = COMP([BEZIERCURVE([[0,0,0],[10,0,0]]),CONS([S1])])
	Su1 = COMP([BEZIERCURVE([[0,10,0],[2.5,10,3],[5,10,-3],[7.5,10,3],[10,10,0]]),CONS([S1]) ])
	S0v = COMP([BEZIERCURVE([[0,0,0],[0,0,3],[0,10,3],[0,10,0]]) , CONS([S2]) ]) 
	S1v = COMP([BEZIERCURVE([[10,0,0],[10,5,3],[10,10,0]]) ,CONS([S2])   ])
	surface = COONSPATCH([Su0,Su1,S0v,S1v])
	VIEW(MAP(surface)(Power(INTERVALS(1)(10),INTERVALS(1)(10))), title="TestThinSolid-1")
	solidMapping = THINSOLID(surface)
	Domain3D = Power(Power(INTERVALS(1)(5),INTERVALS(1)(5)),INTERVALS(0.5)(5))
	VIEW(MAP(solidMapping)(Domain3D), title="TestThinSolid-2")
end

function TestEllipse()
	VIEW(ELLIPSE([1,2])(8), title="TestEllipse")
end

function TestBezierStripe()
	vertices = [[0,0],[1.5,0],[-1,2],[2,2],[2,0]]
	VIEW(Struct([POLYLINE(vertices),Power(BEZIERSTRIPE([vertices,0.25,22]),QUOTE([0.9]))]), title="TestBezierStripe")
end

function TestDisplayNubSpline()
	ControlPoints = [[0,0],[-1,2],[1,4],[2,3],[1,1],[1,2],[2.5,1], [2.5,3], [4,4],[5,0]]
	VIEW(DISPLAYNUBSPLINE([3,[0,0,0,0, 1,2,3,4,5, 6    ,7,7,7,7], ControlPoints]), title="TestDisplayNubSpline")
end

function TestDisplayNurbsSpline()
	knots = [0,0,0,1,1,2,2,3,3,4,4,4]
	_p=sqrt(2)/2.0
	controlpoints = [[-1,0,1], [-_p,_p,_p], [0,1,1], [_p,_p,_p],[1,0,1], [_p,-_p,_p], [0,-1,1], [-_p,-_p,_p], [-1,0,1]]
	VIEW(DISPLAYNURBSPLINE([2, knots, controlpoints]), title="TestDisplayNurbsSpline")
end

function TestMinkowski()
	p = MKPOL(
		[[0,0]],
		[[1]],
		[[1]])
	B = MINKOWSKI([  [-1.0/2.0,-1*sqrt(3.0/2.0)] , [-1.0/2.0,sqrt(3.0/2.0)] , [1,0] ])(p)
	vertices = [[0,0],[1,0],[1,0.5],[0.5,0.5],[0.5,1],[0,1]]
	pol1D = MKPOL(vertices,[[1,2],[2,3],[3,4],[4,5],[5,6],[6,1]],[[1],[2],[3],[4],[5],[6]])
	pol2D = MKPOL(vertices,[[1,2,3,4],[4,5,6,1]],[[1,2]])
	Min0 = STRUCT([T([1,2])(v)(S([1,2])([0.1,0.1])(B)) for v in vertices ])
	Min1 = MINKOWSKI([[0.1*-1.0/2.0,0.1*-1*sqrt(3.0/2.0)],[0.1*-1.0/2.0,0.1*sqrt(3.0/2.0)] , [0.1*1,0.1*0] ])(pol1D)
	Min2 = MINKOWSKI([[0.1*-1.0/2.0,0.1*-1*sqrt(3.0/2.0)],[0.1*-1.0/2.0,0.1*sqrt(3.0/2.0)] , [0.1*1,0.1*0] ])(pol2D)
	VIEW(Struct([Min0,Min1,Min2]), title="TestMinkowski")
end

function TestThinsolid()
	Su0 = COMP([BEZIERCURVE([[0,0,0],[10,0,0]]),CONS([S1])])
	Su1 = COMP([BEZIERCURVE([[0,10,0],[2.5,10,3],[5,10,-3],[7.5,10,3],[10,10,0]]),CONS([S1]) ])
	S0v = COMP([BEZIERCURVE([[0,0,0],[0,0,3],[0,10,3],[0,10,0]]) , CONS([S2]) ]) 
	S1v = COMP([BEZIERCURVE([[10,0,0],[10,5,3],[10,10,0]]) ,CONS([S2])   ])
	surface = COONSPATCH([Su0,Su1,S0v,S1v])
	VIEW(MAP(surface)(Power(INTERVALS(1)(10),INTERVALS(1)(10))), title="TestThinsolid-1")
	solidMapping = THINSOLID(surface)
	Domain3D = Power(Power(INTERVALS(1)(5),INTERVALS(1)(5)),INTERVALS(0.5)(5))
	VIEW(MAP(solidMapping)(Domain3D), title="TestThinsolid-2")
end

function TestEllipse()
	VIEW(ELLIPSE([1,2])(8), title="TestEllipse")
end

function TestBezierStripe()
	vertices = [[0,0],[1.5,0],[-1,2],[2,2],[2,0]]
	VIEW(Struct([POLYLINE(vertices),Power(BEZIERSTRIPE([vertices,0.25,22]),QUOTE([0.9]))]), title="TestBezierStripe")
end

function TestDisplayNubSpline()
ControlPoints = [[0,0],[-1,2],[1,4],[2,3],[1,1],[1,2],[2.5,1], [2.5,3], [4,4],[5,0]]
VIEW(DISPLAYNUBSPLINE([3,[0,0,0,0, 1,2,3,4,5, 6,7,7,7,7], ControlPoints]), title="TestDisplayNubSpline")
end

function TestDisplayNurbsSpline()
knots = [0,0,0,1,1,2,2,3,3,4,4,4]
_p = sqrt(2)/2.0
controlpoints = [[-1,0,1], [-_p,_p,_p], [0,1,1], [_p,_p,_p],[1,0,1], [_p,-_p,_p], [0,-1,1], [-_p,-_p,_p], [-1,0,1]]
VIEW(DISPLAYNURBSPLINE([2, knots, controlpoints]), title="TestDisplayNurbsSpline")
end

function TestMinkowski()
p = MKPOL([[0,0]],[[1]],[[1]])
B = MINKOWSKI([[-1.0/2.0,-sqrt(3.0/2.0)], [-1.0/2.0,sqrt(3.0/2.0)], [1,0]])(p)
vertices = [[0,0],[1,0],[1,0.5],[0.5,0.5],[0.5,1],[0,1]]
pol1D = MKPOL(vertices,[[1,2],[2,3],[3,4],[4,5],[5,6],[6,1]],[[1],[2],[3],[4],[5],[6]])
pol2D = MKPOL(vertices,[[1,2,3,4],[4,5,6,1]],[[1,2]])
Min0 = STRUCT([T([1,2])(v)(S([1,2])([0.1,0.1])(B)) for v in vertices])
Min1 = MINKOWSKI([[0.1*-1.0/2.0,0.1*-sqrt(3.0/2.0)],[0.1*-1.0/2.0,0.1*sqrt(3.0/2.0)],[0.1*1,0.1*0]])(pol1D)
Min2 = MINKOWSKI([[0.1*-1.0/2.0,0.1*-sqrt(3.0/2.0)],[0.1*-1.0/2.0,0.1*sqrt(3.0/2.0)],[0.1*1,0.1*0]])(pol2D)
A = Power(Min2,Q(0.05))
B = Power(Min0,Q(0.70))
C = Power(Min1,Q(0.05))
VIEW(TOP([TOP([A,B]),C]), title="TestMINKOWSKI")
end

function TestPolar()
	VIEW(POLAR(CUBOID([1,1,1])), title="TestPolar")
end

function TestSolidify()
VIEW(SOLIDIFY(STRUCT(AA(POLYLINE)([
		[[0,0],[4,2],[2.5,3],[4,5],[2,5],[0,3],[-3,3],[0,0]],
		[[0,3],[0,1],[2,2],[2,4],[0,3]],
		[[2,2],[1,3],[1,2],[2,2]]]))), title="TestSolidify")
end

function TestDiff()
	mypol1 = T([1,2])([-5,-5])(CUBOID([10,10]))
	mypol2 = S([1,2])([0.9,0.9])(mypol1)
	mypol3 = DIFF([mypol1,mypol2])
	VIEW(STRUCT([
			EX([0,10])(mypol3), T(1)(12),
			LEX([0,10])(mypol3), T(1)(25),
			S(3)(3)(SEX([0,PI])(16)(mypol3))
	]), title="TestDiff")
end

function TestCube()
	VIEW(Cube(3))
end

function TestMapSphere()
	N,M = 16,16
	domain = Power(Quote([pi/N]*N).translate([-pi/2]), Quote([2*pi/M]*M))
	obj = domain.mapFn(p -> [cos(p[1])*sin(p[2]), cos(p[1])*cos(p[2]), sin(p[1])])
	VIEW(obj)
end

function TestMkPol()
	out = MkPol([[0],[1],[2],[3],[4],[5]],[[5,3],[0,1]])
	VIEW(out)
end


