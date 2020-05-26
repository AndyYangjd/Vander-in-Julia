#=	program6.jl

手写二次插值法，参考：http://mathonline.wikidot.com/deleted:quadratic-polynomial-interpolation =#

using DataInterpolations, Plots;gr()

t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3];
u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22];

n=length(t);
if n>=3
	data=[];
	t_list=[];
	for i in 1:(n-2)
		L0(x)=(x-t[i+1])*(x-t[i+2])/((t[i]-t[i+1])*(t[i]-t[i+2]));
		L1(x)=(x-t[i])*(x-t[i+2])/((t[i+1]-t[i])*(t[i+1]-t[i+2]));
		L2(x)=(x-t[i])*(x-t[i+1])/((t[i+2]-t[i])*(t[i+2]-t[i+1]));
		P(x)=u[i]*L0(x)+u[i+1]*L1(x)+u[i+2]*L2(x);
		if i != (n-2)
			piece_t=t[i]:0.01:t[i+1];
		else
			piece_t=t[i]:0.01:t[i+2];
		end
		push!(t_list, piece_t);
		piece_data=P.(piece_t);
		push!(data, piece_data);
	end
end

flat_data =  collect(Iterators.flatten(data));
flat_t = collect(Iterators.flatten(t_list));

plot(flat_t, flat_data, label="Fitting")
scatter!(t, u, label="input data")