# ppkdmeans3
to be continue...

10/21
有个很严肃的问题是
我现在没法把参数设置成我想要的样子
比如控制one ciphertext one number

以及涉及乘法深度和噪声的参数设计（麻了，我好菜！

❗虽然原论文的array_min函数可以用
但是比较鸡肋，只是一个简单的demo，最好还是自己去写

10/22

kd tree构造部分

除了加密以外所有东西的框架都构造好了

这里发现了一个问题也就是

本来类个方差嘛

我不是简化成了
$s = n\sum{x_i^2} + (\sum{x_i})^2$

但是发现中间结果很容易就超过了加密的范围

eg，现在范围是65535，然后算出来400000+

这不是**n**已知嘛

我们再简化一下
$s = \sum{x_i^2} + (\sum{x_i})^2/2$

再看看