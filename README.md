## 介绍
本仓库使用Julia编写练习使用的代码，包括FFT、优化算法等。

## 文件内容简介
- fft_practice.jl: 使用FFT绘制一个周期信号的频谱
- vander.jl: 仅包含一个函数，通过输入向量得出其对应的范德蒙矩阵，对LinearALagebra库中特殊函数的补充
- program1.jl：采用梯度下降+最大下降步长的方法求一个二次函数极小值
- program2.jl：采用梯度下降+Armijo-Goldstein准则算法求一个二次函数极小值
- program3.jl：采用阻尼牛顿法求一个二次函数极小值，结果证明最佳，一步到位。
- program4.jl：采用拟牛顿法（基于DFP算法的变尺度法）。