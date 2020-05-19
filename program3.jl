#=	program3.jl

使用阻尼牛顿法求函数 f(x)=x1^2+25x2^2的极小点。 =#

using LinearAlgebra;

# 将 f(x)转换为矩阵形式，也可以不用转，但求梯度方便。
G = [2 0;0 50];
gradient_fx(x) = G * x;
f(x) = 0.5 * x' * G * x;
f2(x1, x2) = x1^2 + 25x2^2;

x_next = [2;2];
x_list = copy(x_next);
y_list = [f(x_next)];
grad_list = [0.0;0.0];
alpha_list = [0.0];

TOL = 0.0001;
while norm(x_next) > TOL
	global x_next, x_list, y_list, grad_list, alpha_list;
	# 计算当前点梯度
	gradient_k = gradient_fx(x_next);
	# 设置搜索方向为牛顿方向
	d_k = -inv(G) * gradient_k;
	# 计算步长为最佳步长
	alpha = (x_next' * G * [1;1]) / (gradient_k' * [1;1]);
	# 计算下一点
	x_next = x_next + alpha * d_k;
	# 将运算过程中的值保存
	x_list = hcat(x_list, x_next);
	y_list = push!(y_list, f(x_next));
	grad_list = hcat(grad_list, gradient_k);
	alpha_list = push!(alpha_list, alpha);
end

data = hcat(hcat(hcat(x_list', y_list), alpha_list), grad_list');
data = @. round(data, digits = 4);

println("k \t x1 \t x2 \t f(x) \t alpha \t Gradient_x1 \t Gradient_x2");
row_len, col_len = size(data);
for row in 1:row_len
	print(row, "\t")
	for col in 1:col_len
		print(data[row, col], "\t");
	end
	println("");
end