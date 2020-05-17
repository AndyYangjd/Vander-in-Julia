"""
	Vander(x, order::Int=0, increasing::Bool=true)

Create a Vander Matrix depends on the input Vector. 

...
# Arguments
- 'x::Union{AbstractVector, Array{::Real,2}': one dimension Vector or Matrix.
- 'order': the row or column dimension of Matrix depends on increasing.
- 'increasing::Bool=true': decide the orientation of Matrix.
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
function Vander(x::Union{AbstractVector, Array{<:Real,2}}, order::Int=0, increasing::Bool=true)
	row_len=length(x);

	if row_len==0
		error("Pls input a vector which length > 0");
	end

	if typeof(x)<:Array{<:Real,2}
		if size(x)[1] != 1
			error("Plse input 1 dimension Vector");
		end
	end

	if order<0
		error("Pls input a order>0")
	end
	col_len= order>0 ? order : row_len

	T=typeof(x[1]);
	if T <: Number
		vander_matrix=Matrix{T}(undef, row_len, col_len);

		for i in 1:col_len
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