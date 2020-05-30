#=	program8.jl

Chebyshev插值法 =#

"""
	function getChebyshevPointsX(a, b, n)

给定闭区间的端点及插值点数，返回Chebyshev插值点. 

...
# Arguments
- 'a': 区间左端点.
- 'b': 区间右端点.
- 'n': 区间内的插值点数量（注意不包含闭区间的两个端点，但返回值会添加）.
...
# Examples
```julia-repl
julia> x0=[0,2,3];
julia> y0=[1,2,4];
julia> c=newTdd(x0,y0,3)
3-element Array{Float64,1}:
 1.0
 0.5
 0.5
```
"""
function getChebyshevPointsX(a, b, n_, addEndPoints=true)
	if a<b
		factor1, factor2 = (b-a)/2.0, (b+a)/2.0; 
		Chebyshev_f(n)=[factor1*cos((2i-1)*pi/(2*n)) + factor2 for i =1:n];
		PointsX=Chebyshev_f.(n_)
		sort!(PointsX);
		if addEndPoints	
			insert!(PointsX, 1, a);
			insert!(PointsX, n_+2, b);
			return PointsX;
		else
			return PointsX;
		end
	else
		error("Pls input a<b");
	end
end