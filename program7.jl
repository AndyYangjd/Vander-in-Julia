#=	program7.jl

牛顿差商插值法 =#

"""
	newTdd(x,y, n)

给定两个向量，求出使用牛顿差商插值多项式的系数. 

...
# Arguments
- 'x': 横轴x向量.
- 'y': 竖轴y向量.
- 'n': 向量的长度.
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
function newTdd(x, y, n)
	v=zeros(Float64, n, n);
	for i=1:n
		v[i, 1]=y[i];
	end

	for i=2:n
		for j=1:n+1-i
			v[j,i]=(v[j+1, i-1]-v[j, i-1])/(x[j+i-1]-x[j]);
		end
	end

	return v[1,:];
end