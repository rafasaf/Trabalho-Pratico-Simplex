# Trabalho-Pratico-Simplex

O objetivo deste trabalho é resolver PLs gerais, a serem fornecidas e cujo formato será especificado abaixo. Em outras palavras, vamos fazer uma aplicação do método simplex.

Resolva a programação linear definida por

$$
\max \ \mathbf{c}^T \mathbf{x} \\
\text{sujeito a} \ A\mathbf{x} = \mathbf{b} \\
\mathbf{x} \geq 0
$$

e encontre o certificado que comprove seu resultado.

$$
A = \begin{pmatrix}
a_{1,1} & a_{1,2} & \cdots & a_{1,m} \\
a_{2,1} & a_{2,2} & \cdots & a_{2,m} \\
\vdots & \vdots & \ddots & \vdots \\
a_{n,1} & a_{n,2} & \cdots & a_{n,m}
\end{pmatrix} \quad
\mathbf{b} = \begin{pmatrix}
b_1 \\
b_2 \\
\vdots \\
b_n
\end{pmatrix} \quad
\mathbf{c} = \begin{pmatrix}
c_1 \\
c_2 \\
\vdots \\
c_m
\end{pmatrix}
$$
