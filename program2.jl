#=	program2.jl

使用梯度方向+Goldstein-Armijo准则求函数 f(x)=x1^2+25x2^2的极小点。
迭代起始点为 x0=[2;2]，Armijo准则参数 u1=0.2, u2=0.8, p1=0.5, p2=1.5。
显然已知最小点为(0,0)，按前向误差|x*-x|的二范数最小设置终止条件,
也可按照|x[k+1]-x[k]|< TOL来设定终止准则。 =#

using LinearAlgebra

# 将 f(x)转换为矩阵形式，也可以不用转，但求梯度方便。
G = [2 0;0 50];
gradient_fx(x) = G * x;
f(x) = 0.5 * x' * G * x;

x_next = [2;2];
x_list = copy(x_next);
y_list = [f(x_next)];
grad_list = [0.0;0.0];
alpha_list = [0.0];

# 设置Armijo准则参数
u1 = 0.2;
u2 = 0.8;
p1 = 0.5;
p2 = 1.5;

TOL = 0.0001;
while norm(x_next) > TOL
	global x_next, x_list, y_list, grad_list, alpha_list;
	# 计算当前点梯度
	G_k = gradient_fx(x_next);
	y_k = f(x_next);
	# 设置搜索方向
	dk = -G_k;
	# 计算步长
	alpha = 1.;
	# 设置一个布尔量控制循环跳出
	gate_alpha = true;
	while gate_alpha
		if (-alpha * u1 * G_k' * dk) > (y_k - f(x_next + alpha * dk))
			alpha = p1 * alpha;
			continue;
		end
    	
		if (y_k - f(x_next + alpha * dk)) > (-alpha * u2 * G_k' * dk)
			alpha = p2 * alpha;
			continue;
		else
			gate_alpha=false;
		end 
	end

	# 计算下一点
	x_next = x_next - alpha * G_k;
	# 将运算过程中的值保存
	x_list = hcat(x_list, x_next);
	y_list = push!(y_list, y_k);
	grad_list = hcat(grad_list, G_k);
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