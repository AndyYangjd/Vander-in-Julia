#= 
参考：https://blog.csdn.net/wordwarwordwar/article/details/62078726 
胡广书的《数字信号处理－理论、算法与实现（第二版）》第三章第八节《关于正弦信号抽样的讨论》，
得出了关于正弦信号抽样的六个结论，最后总结了一个总的原则：抽样频率应为信号频率的整数倍，抽样点数应包含整周期。
DFT就是DFS，只不过DFT先将有限长信号进行周期延拓，然后求DFS，再截取一个周期。 =#

using Plots,FFTW;pyplot()

# 待绘制频谱信号的函数表达式，由一个直流信号和两个余弦信号构成
# Ts1=0.02s ,Ts2=0.01s; fs1=50, fs2=100;
s(t) = 0.2 + 0.7 * cos(2 * pi * 50t + pi / 9) + 0.2 * cos(2 * pi * 100t + 70 / 180 * pi);

# 采样率是s(t)最大频率的整数倍，保证了采样之后的序列是周期序列，当然也要满足采样定理要求
Fs = 1000;
T = 1 / Fs;
# s(t)的周期是0.02s，以上述采样率采样一个周期需要20个采样长度，因此1000个采样长度保证至少采样一个周期，而且总采样长度是完整周期
L = 1000;

# 求信号的采样点处的函数值
t = [i * T for i in 0:L - 1];
value_s = s.(t);

# 绘制原始信号频谱
p1 = plot(t[1:60], value_s[1:60], title = "Original Signal s(t)", xlabel = "t/s", ylabel = "s", label = "No Sampling");
plot!(p1, t[1:60], value_s[1:60], line = :stem, marker = :circle, label = "Sampling Points");

# 对采样的信号值进行FFT
sk = fft(value_s); 
# 对FFT结果归一化
normed_sk = sk * T;

# 进行IFFT并绘制还原后的信号
is = ifft(sk);
p2 = plot(t[1:60], real(is)[1:60], title = "IFFT Signal", xlabel = "t/s", ylabel = "ifft_s", label = "No Sampling");
plot!(p2, t[1:60], real(is)[1:60], line = :stem, marker = :circle, label = "Sampling Points");

# fftfreq实现了频率的对称移位功能，不用再使用fftshift函数把频谱移动到以0为对称轴了。
# 也就是说，fftfreq和fftshift两个函数只能二选一使用，强烈建议使用前者，省去很多代码量。
freqs = fftfreq(L, Fs);

# 绘制信号的幅频特性
mag = abs.(sk);
p3 = plot(freqs, mag, line = :stem, marker = :circle, label = "No normalization", title = "Mag", xlabel = "f/Hz");

# 绘制归一化后原始信号的频谱
normed_mag = abs.(normed_sk);
p4 = plot(freqs, normed_mag, line = :stem, marker = :circle, label = "Normalization", title="Normalized Mag", ylim = (0, 0.4));
annotate!([
	(freqs[1], normed_mag[1], Plots.text("$(round(normed_mag[1];sigdigits = 2))", :red, :left)), 
	(freqs[51], normed_mag[51], Plots.text("$(round(normed_mag[51];sigdigits = 2))", :red, :left)), 
	(freqs[101], normed_mag[101], Plots.text("$(round(normed_mag[101];sigdigits = 2))", :red, :left))
]);

# 求信号的相位
normed_phase = angle.(normed_sk);
p5 = plot(freqs, normed_phase, line = :stem, marker = :circle, title = "Phase", xlim = (0, maximum(freqs) + 1), xlabel = "f/Hz", ylabel = "radian");
annotate!([
	(freqs[1], normed_phase[1], Plots.text("$(round(normed_phase[1];sigdigits = 2))", :red, :left)), 
	(freqs[51], normed_phase[51], Plots.text("$(round(normed_phase[51];sigdigits = 2))", :red, :left)), 
	(freqs[101], normed_phase[101], Plots.text("$(round(normed_phase[101];sigdigits = 2))", :red, :left))
]);

# 打印指定频率处的值
mag_bool = isapprox(normed_mag[[1,51, 101]], [0.2, 0.35, 0.1]);
println("FFT得出的幅频值与原信号是否一致： $mag_bool");
phase_bool = isapprox(normed_phase[[1,51, 101]], [0, pi / 9, 7 * pi / 18]);
println("FFT得出的相角值与原信号是否一致： $phase_bool");

display(plot(p1, p2, p4, p5, layout = (2, 2)));