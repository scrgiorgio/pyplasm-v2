using Plasm
using Test

@testset "Plasm.jl" begin

	using Test, Documenter, Plasm
	doctest(Plasm)
	
	if false

	# viewer
	GLView([
		GLCuboid(Box3d(Point3d(0,0,0),Point3d(1,1,1)))
		GLAxis(Point3d(0,0,0),Point3d(+1.1,+1.1,+1.1))
		])	

	# hpc
	TestComputeNormal()
	TestGoodTet()
	TestBox()
	TestMat()
	TestMkPol()
	TestHpc()
	println("all test ok")

	# fenvs


	tests = MyTests()
	tests.TestCube()
	tests.TestMapSphere()
	tests.TestMkPol()
	tests.TestSphere()
	tests.TestTorus()
	tests.TestBezier()
	tests.TestCoonsPatch()
	tests.TestRuledSurface()
	tests.TestProfileProdSurface()
	tests.TestRotationalSurface()
	tests.TestCylindricalSurface()
	tests.TestConicalSurface()
	tests.TestCubicHermite()
	tests.TestPermutahedron()
	tests.TestSchegel3d()
	tests.TestCubicSpline()
	tests.TestBilinarSurface()
	tests.TestBiquadraticSurface()
	tests.TestHermiteSurface()
	tests.TestBezierSurface()
	tests.TestBezierManifold()
	tests.TestOffset()
	tests.TestThinSolid()
	tests.TestEllipse()
	tests.TestBezierStripe()
	tests.TestDisplayNubSpline()
	tests.TestDisplayNurbsSpline()
	tests.TestMinkowski()

	end


	
end
