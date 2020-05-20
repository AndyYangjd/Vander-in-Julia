#=	program5.jl

使用共轭梯度法求函数 f(X)=0.5*X^T*G*X+b^TX的极小点,
其中G=[2 -2; -2 4], b=[-4;0]。 =#

using LinearAlgebra;

# 定义函数
G = [2 -2; -2 4];
b = [-4;0];
f(X) = 0.5 * X' * G * X + b' * X;
# 求梯度
g_fx(X) = G * X + b;

# 设置初始参数
x_k = [1.0;1];
# 设置误差
TOL = 0.0001;
# 设置最大循环次数
N_MAX = 30;
# 设置布尔量控制循环
gate = true
# 设置数组来保存运算过程中的结果
x_list = [];
y_list = [];
g_list = [];
alpha_list = [];
conj_list = [];

while gate
	global x_k, conj_k, x_list, y_list, g_list, alpha_list, conj_list, gate
	# 保存当前的 x 向量
	push!(x_list, x_k);
	# 保存当前的 y 函数值
	push!(y_list, f(x_k));
	# 计算当前点梯度
	g_k = g_fx(x_k);
	push!(g_list, g_k);
	# 如果梯度为0，则退出
	if norm(g_k) < TOL
		gate = false;
	end
	# 求共轭方向
	if length(conj_list) > 0
		beta_k = norm(g_k)^2 / norm(g_list[end - 1])^2;
		conj_k = -g_k + beta_k * conj_k;
		push!(conj_list, conj_k);
	else
		conj_k = -g_k;
		push!(conj_list, conj_k);
	end
	# 设置搜索方向
	d_k = conj_k;
	# 计算步长为最佳步长
	alpha = -(d_k' * (G * x_k + b)) / (d_k' * G * d_k);
	push!(alpha_list, alpha);
	# 计算下一点
	x_k = x_k + alpha * d_k;
	# 判断是否达到精度要求
	if norm(x_k - x_list[end]) < TOL || length(x_list) > N_MAX
		gate = false;
	end
end

# 将数据整合在一起
row_len = length(x_list);
col_len = sum(map(length, [x_list[1], g_list[1], conj_list[1]])) + 2;
data = zeros(Float64, (row_len, col_len));
for i in 1:length(x_list)
	row_x, row_g, row_H = map(vec, [x_list[i], g_list[i], conj_list[i]]);
	row = hcat(hcat(hcat(hcat(row_x', y_list[i]), alpha_list[i]), row_g'), row_H');
	data[i, :] = row;
end
# 为了显示美观，将数据保留4位小数
data = @. round(data, digits = 4);

println("k \t x1 \t x2 \t f(x) \t alpha \t Gradient_x1 \t Gradient_x2 \t Conj-Vector(2)");
for row in 1:row_len
	print(row, "\t")
	for col in 1:col_len
		print(data[row, col], "\t");
	end
	println("");
end