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
$s = \sum{x_i^2} + (\sum{x_i})^2/n$

再看看

10/23

就是说，方差加密超过范围的问题

经过上面的简化，还是会超过嘛

试了方法一：数组中所有数据都减去最小值（无法保证

方法二：所有数都除以100->>缩小100倍，暂时能用，但是没啥技术含量，而且要根据数据的范围来变化，以后肯定会出问题，💡那到时候再说吧！

❗❗大问题

也就是说，方差比较的时候要用到，全为0和1的密文

但是不能直接加密1，因为这样子得到的编码结果就是
[1]
[]
[]
[]
[1]
[]
[]
[]
[1]
[]
[]
[]
[1]
[]
[]
[]
[1]
[]
[]
[]

然后俺再和比如说，200的密文乘一下
[0 2]
[0 3]
[]
[]
[0 2]
[0 3]
[]
[]
[0 2]
[0 3]
[]
[]
[0 2]
[0 3]
[]
[]
[0 2]
[0 3]
[]
[]
整一个就是直接GG

为什么呢，你看结果
[0 2]
[]
[]
[]
[0 2]
[]
[]
[]
[0 2]
[]
[]
[]
[0 2]
[]
[]
[]
[0 2]
[]
[]
[]
每个数字第二行无了

现在想了一个很简陋的办法，也就是，咱不加密1，咱们加密$0*16^0 + 1*16^1+1*16^2+1*16^3 = 4369$，然后结果就是

[1]
[1]
[1]
[1]
[1]
[1]
[1]
[1]
[1]
[1]
[1]
[1]
[1]
[1]
[1]
[1]
[1]
[1]
[1]
[1]

再去和200的密文做乘法，欸，成了

[0 2]
[0 3]
[]
[]
[0 2]
[0 3]
[]
[]
[0 2]
[0 3]
[]
[]
[0 2]
[0 3]
[]
[]
[0 2]
[0 3]
[]
[]

❗❗❗这个方法不知道以后会不有啥大问题，暂且先用着，以及这样耦合程度就增加了，一旦俺后面改变了编码规则，比如说每个数只能编码到$a*16^0+b*16^1+c*16^2$空间里，那就要随时改这个值，到时候再看吧！

## 11/03
$\sum_{i=1}^{d}\left(\alpha c_{i}-\alpha_{j} k_{j i}\right)^{2}=\alpha^{2} \sum c_{i}^{2}-2 \alpha \alpha_{j} \sum c_{i} k_{j i}+\alpha_{j}^{2} \sum k_{j i}^{2}$
关于这一部分乘法，用不同的方法，可以优化乘法的交互次数

如果按照公式左边，纯纯的计算距离

乘法的次数是((1+1)+1)*d = O(3d)次交互，没办法并行

按照右边的化简公式：一次交互同步计算 $c_i^2, c_ik_{ji}, k_{ji}^2,\alpha^2, \alpha_j^2,\alpha\alpha_j$，然后计算

