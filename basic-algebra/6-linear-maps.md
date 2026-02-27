# 线性映射

本章节我们讨论线性映射。首先，线性映射是 $F^n$ 到 $F^k$ 的函数集，其次这些函数是 $\mathbb{F}$ 域中两个线性空间的函数集。这些讨论会围绕着这些函数的矩阵意义构造计算。

我们观察任意的 k-by-n 的矩阵，定义函数 $L: \mathbb{F}^n\, to\, \mathbb{F}^k$ 为 $L(v)=Av$，并且这个函数对所有 $\mathbb{F}^n$ 中的 u 和 v 以及所有的标量 c 都满足

$$
\begin{aligned}
L(u+v)&=L(u)+L(v) \\
L(cv)&=cL(v)
\end{aligned}
$$

如一个函数 $\mathbb{F}^n \mapsto \mathbb{F}^k$满足这两个条件都可以成为线性的：要么强调 $\mathbb{F}$是线性的，要么强调标量。传统上这样的函数称为线性映射或者线性变换。因此线性映射对应矩阵，反之亦然。

## 命题 2.11

如果 $L:\mathbb{F}^n \mapsto \mathbb{F}^k$ 是一个线性映射，那么存在一个唯一的 k-by-n 的矩阵，对所有 $\mathbb{F}$ 中的$L(v)=Av$ 满足 $L(v)=Av$。

### 强调

这个证明会展示如何得到矩阵 A。

### 证明

对于 $1 \le j \le n$，令 $e_j$ 为 $\mathbb{F}^n$ 的第j个标准基向量，它的 $j^{th}$ 项是1，其它项为0。令 A 的 $j^{th}$ 列为 k 维列向量 $L(e_j)$。如果 v 是列向量 $(c_1,c_2,\cdots,c_n)$，那么

$$
\begin{aligned}
L(v)&=L(\Sigma^n_{j=1})=\Sigma^n_{j=1}L(c_je_j)\\
&=\Sigma^n_{j=1}c_jL(e_j)=\Sigma^n_{j=1}c_j(j^{th}\, columns\, of \, A).
\end{aligned}
$$

令 $L(v)_i$ 为列向量 $L(v)$ 的$i^{th}$项，这个方程就成为

$$
L(v)_i = \Sigma^n_{j=1}c_jA_{ij}
$$

右边是 $Av$ 的 $i^{th}$ 项，因此 $L(v) = Av$。存在性得证。至于唯一性，我们观察 A 的 $j^{th}$ 列，对于每个 j 都有形如 $L(e_j)=Ae_j$ 成立，唯一性得证。

### 例1

令 $L: \mathbb{R}^2 \mapsto \mathbb{R}^2$ 为将复数向量逆时针旋转 $\theta$ 的计算。考虑 L 的几何意义，由向量加法的平行四边形法则可以得到，L是线性的。计算过程显示 $L(1\,0)=\left(\begin{matrix} cos\, \theta\\ sin\, \theta \end{matrix}\right)$。应用命题 2.11和命题证明中构造矩阵的方法，我们可以得到对于 $\mathbb{R}^2$ 中所有的 v 都有 $L(v)=\left(\begin{matrix}cos\,\theta\,& -sin \, \theta \\ sin\,\theta &\,\cos\,\theta\end{matrix}\right)$。

我们可以通过将对应位置的系数相加来将两个向量 $L: \mathbb{F}^n \mapsto \mathbb{F}^k $ 和 $M: \mathbb{F}^n+\mathbb{F}^k$ 相加： $(L+M)(v)=L(v)+M(v)$。我们也可以通过将每个位置的系数都与系数相乘来实现一个向量与标量的相乘。那么 $L+M$ 和 $cL$ 是线性的。并且它们作为 $\mathbb{F}^n $ 到 $\mathbb{F}^k$ 的所有函数构成的向量集的向量子空间。因此它本身也是一个向量空间。传统上这个向量空间的符号是 $Hom_{\mathbb{F}}(\mathbb{F}^n, \mathbb{F}^k)$，符号 Hom 表示规则 $L(u+v)=L(u)+L(v)$，而子式 $\mathbb{F}$ 表示对于所有 $\mathbb{F}$ 中的 c 都满足规则 $L(cv)=cL(v)$。

加入 L 表达为矩阵 A，M 表达为矩阵 B，那么 $L + M$对应 $A+B$，$cL$ 对应 $cA$。下一个命题展示了线性映射的组合对应为矩阵乘法的过程。

## 命题 2.12

令 $L: \mathbb{F}^n \to \mathbb{F}^m$ 为 m 行 n 列矩阵 A 对应的线性映射，且 k 行 m 列矩阵 B 对应的线性映射为 $ M: \mathbb{F}^m \to \mathbb{F}^k $。那么组合函数 $ M \circ L : \mathbb{F}^n \to \mathbb{F}^k$ 是线性的，且它对应 k 行 n 列的矩阵 BA 。

### 证明

函数 $ M \circ L$ 满足 $ (M \circ L)(u+v) = M(L(u+v)) = M(Lu + Lv)= M(Lu)+M(Lv)=(L \circ L)(u)+(M \circ L)(v)$ 并且它满足 $(M \circ L)(cv)=c(M \circ L)(v)$。因此它是线性的。线性映射对应矩阵的乘法结合律由 $ (M \circ L)(v) = M(L(v)) = (B)(Lv)=B(Av) = (BA)v $ 给出，因此 $ M \circ L $ 对应于 BA 。

现在我们扩大讨论范围，对任意 $\mathbb{F}$ 的向量空间之间的 $ L : U \to V$。对于 U 中所有的 u 和 v，以及所有的标量 c ，如果满足下列条件，我们称 $ L : U \to V $ 是线性的，或者 $ \mathbb{F}$ 是线性的：

$$
\begin{aligned}
L(u+v) &= L(u)+L(v)\\
L(cv) &= cL(v)
\end{aligned}
$$

对于特定的 $U=\mathbb{F}^n$和 $V=\mathbb{F}^k$，线性函数被称为线性映射或者线性变换。所有线性映射的集合 $ L:V \to V$ 是一个 $\mathbb{F}$域的线性空间，由 $Hom_\mathbb{F}(U, V)$ 定义。以下结果是处理线性映射的基础。

## 命题 2.13

令 U 和 V 是 $\mathbb{F}$ 的向量空间，$\Gamma$ 是 U 的基。那么每个函数 $\ell : \Gamma \to V$ 对应一个且仅有一个线性映射 $ L : U \to V$ 严格限定于 $L\big |_\Gamma = \ell$。

### 强调

L 是 $\ell$ 的线性扩展。