#=	image_of_f.jl

绘制等高线图解释Armijo-Goldstein准则，函数f(x,y)=x^2+25y^2 =#

using  LinearAlgebra, Plots;gr()

f(x,y) = x^2 + 25y^2
f1(x, z) = sqrt(z - x^2) / 5;
f2(x, z) = -sqrt(z - x^2) / 5;
# 绘制(2,2)点处的等高线
z = f(2, 2);
x = -3:0.001:3;
plot(x, hcat(f1.(x, z), f2.(x, z)),c = [:blue :blue],legend = false, framestyle = :zerolines)

# 绘制(2,2)到下一迭代点(1.875, -1.125)的长度
plot!([2;1.875], [2;-1.125], c = [:red], marker = :dot, arrow = 2)

# 绘制(1.875, -1.125)处的等高线
z = f(1.875, -1.125);
plot!(x, hcat(f1.(x, z), f2.(x, z)),c = [:yellow :yellow],legend = false, framestyle = :zerolines)

# 绘制函数值下降最多处(1.9199, -0.0031)点的等高线
z = f(1.9199, -0.0031);
x_max = floor(sqrt(z), digits = 6);
x2 = -x_max:0.01:x_max;
plot!(x2, hcat(f1.(x2, z), f2.(x2, z)),c = [:black :black],legend = false, framestyle = :zerolines)

# 绘制第二个点的一阶泰勒展开的等高线
alpha = 0.0312;
z = f(1.875, -1.125) - alpha * dot([4;100], [4;100]);

gui()