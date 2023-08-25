struct Test
	x::Int64

	function Test(x::Int64)
		println("constructor")
		new(x)
	end

	function maybeInside()
		println("maybeInside")
	end

end

function testFunction1()
	println("otherMethod")
end

function testFunction2(self::Test)
	println("otherMethod")
end


t=Test(10)
testFunction1()
testFunction2(t)
Test.maybeInside()