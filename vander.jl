"""
	Vander(x, increasing::Bool=true)

Create a Vander Matrix depends on the input Vector. 

...
# Arguments
- 'x::Union{AbstractVector, Array{::Real,2}': one dimension Vector or Matrix.
- 'increasing::Bool=true':
...
# Examples
```julia-repl
julia> Vander([1 2 3])
[ 1 1 1
  1 2 4
  1 3 9 ]
julia> Vander([1, 2, 3], false)
[ 1 1 1
  1 2 3
  1 4 9 ]
```
"""
function Vander(x::Union{AbstractVector, Array{<:Real,2}}, increasing::Bool=true)
	x_len=length(x);

	if x_len==0
		error("Pls input a vector which length > 0");
	end

	if typeof(x)<:Array{<:Real,2}
		if size(x)[1] != 1
			error("Plse input 1 dimension Vector");
		end
	end

	T=typeof(x[1]);
	if T <: Number
		vander_matrix=Matrix{T}(undef, x_len, x_len);
		for i in 1:x_len
			vander_matrix[:, i]=@. x^(i-1);
		end
	
		if increasing
			return vander_matrix;
		else
			return transpose(vander_matrix);
		end
	else
		error("Pls input Numerical Vector");
	end
end