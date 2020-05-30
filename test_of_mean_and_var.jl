#=  test_of_mean_and_var.jl

    从一个指数分布（lambda=1/4.5）中进行10^6次实验。使用Monte Carlo方法，
每次抽取10个数据组成一组，求取其均值和方差，然后将所有组的均值和方差直方图画出并分析。
注：指数分布是一种连续分布，其均值为1/lambda,方差为 (1/lambda)^2。 
程序参考并修改了StatisticsWithJulia Listing 5.1
=#

using Random, Distributions, Plots;gr()
Random.seed!(0);

lambda = 1 / 4.5;
expDist = Exponential(1 / lambda);
n, N = 10, 10^6;

means = Array{Float64}(undef, N);
variances = Array{Float64}(undef, N);
data_list = [];

for i in 1:N
    data = rand(expDist, n);
    push!(data_list, data);
    means[i] = mean(data);
    variances[i] = var(data);
end

flat_data = collect(Iterators.flatten(data_list));
flat_mean = mean(flat_data);
flat_var = var(flat_data, corrected = false, mean = flat_mean);

println("实际均值: ", mean(expDist), "\t实际方差: ", var(expDist))
println("所有数据的均值:", flat_mean, "\t方差: ", flat_var)
println("分组采样均值的均值: ", mean(means))
println("分组采样方差的均值: ", mean(variances))

stephist(means, bins = 200, c = :blue, normed = true, label = "Histogram of Sample Means");
stephist!(variances, bins = 600, c = :red, normed = true, label = "Histogram of Sample variances");
plot!(xlims = (0, 40),ylims = (0, 0.4),xlabel = "Statistic Value", ylabel = "Density");
png("test_mean_and_var");