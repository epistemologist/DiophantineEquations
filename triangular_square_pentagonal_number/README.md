We show that there does not exist a [triangular square pentagonal number](https://mathworld.wolfram.com/PentagonalSquareTriangularNumber.html) following the approach given in [this paper](https://www-math.st.tokushima-u.ac.jp/journal/2020/2020-1.pdf) as well as extending to other cases.

We define the $n^{th}$ polygonal number of $s$ sides to be $$P_s(n) = \frac{1}{2}\left( (s-2) n^2 - (s-4)n\right) $$

### Pentagonal Triangular Square Number
Assume that there existed a pentagonal triangular square number; i.e. there existed a solution to $P_3(a) = P_4(b) = P_5(c) \iff a^2+a = 2b^2 = 3c^2 - c$. We have:

$$\begin{cases}X = 2a+1 \\\\ Y = 2b  \\\\ Z = 6c-1\end{cases}\implies \begin{cases}X^2 = 4a^2+4a+1 \\\\ Y^2 = 4b^2 \\\\ Z^2 = 36c^2 -12c+1\end{cases}$$

Therefore, we have:
```math
\begin{aligned}
X^2 = 4a^2+4a+1 &=4(a^2+a)+1 \\ 
&= 4(2b^2) + 1 \\ 
&= 8b^2+1 \\ 
&=2Y^2 + 1
\end{aligned}
```
and so $X^2-2Y^2=1$.

Similarly, 
```math
\begin{aligned}Z^2 = 36c^2-12c+1&=12(3c^2-c)+1 \\
&=  12(2b^2)+1 \\
&= 24b^2+1 \\
&= 6Y^2+1\end{aligned}
``` 
and so $Z^2-6Y^2=1$.

We therefore have the system of Pell equations $$X^2-2Y^2=1, Z^2-6Y^2=1$$ 

Notice that:
```math
\begin{aligned}
X^2 Z^2 &= (1+2Y^2)(1+6Y^2) \\
&\implies (XYZ)^2 = Y^2(1+2Y^2)(1+6Y^2) \\
&\implies 144(XYZ)^2 = (12XYZ)^2 = 12Y^2 (6+12Y^2)(2 + 12 Y^2)
\end{aligned}
```

Letting $u = 12XYZ, v = 12Y^2$, we have the equation $u^2 = v(v+2)(v+6)$ which is an elliptic curve.

We sanity check in Sage:
```python
sage: R.<a_, b_, c_> = QQ[]
....: I = Ideal( a_^2+a_-2*b_^2, 3*c_^2-c_-2*b_^2 )
....: Q = R.quotient(I); a,b,c = Q.gens()
....: X = 2*a+1; Y = 2*b; Z = 6*c-1
....: u = 12*X*Y*Z; v = 12*Y^2
sage: u^2 - v*(v+2)*(v+6)
0
```
We therefore have that if $a,b,c$ are integers satisfying the above equation that $u=24 \cdot (6 c - 1) \cdot b \cdot (2 a + 1)=288 a b c - 48 a b + 144 b c - 24 b$ and $v = 48b^2$ satisfies the equation $u^2 = v(v+2)(v+6)$. We check all possible integral points on this elliptic curve with Sage:

```python
sage: R.<u,v> = QQ[]
....: a,b,c = var("a b c", domain="integer")
....: sols = []
....: E = EllipticCurve(u^2 - v*(v+2)*(v+6))
....: for p in E.integral_points(both_signs = True):
....:     v_sol, u_sol = p.xy()
....:     # v = 48b^2 => b = sqrt(v/48)
....:     b = sqrt(v_sol/48)
....:     if b.is_integer():
....:         s = solve([24*(6*c-1)*b*(2*a+1) - u_sol], a,c)
....:         try:
....:             for a_,c_ in s:
....:                 assert 24*(6*c_-1)*b*(2*a_+1) == u_sol
....:                 sols.append( (int(a_),b,int(c_)) )
....:                 print(sols[-1])
....:         except:
....:             continue
....: 
(-2, 1, 1)
(7, 1, 0)
(-8, 1, 0)
(1, 1, 1)
```

Note for the case of `b=u_sol=0`, we have that `24*(6*c-1)*b*(2*a+1) - u_sol` reduces to 0, so we have the additional solutiion `(0,0,0)`. We therefore have the only pentagonal triangular square number corresponds to the solution $a=b=c=1$ (i.e 1 is the only solution). Extending to negative numbers, we have:
```python
sage: sols.append( ( 0, 0, 0 ) )
....: for a,b,c in sols:
....:     print(f"a={a}, b={b}, c={c}, {a^2+a} = {2*b^2} = {3*c^2-c}")
....: 
a=-2, b=1, c=1, 2 = 2 = 2
a=7, b=1, c=0, 56 = 2 = 0
a=-8, b=1, c=0, 56 = 2 = 0
a=1, b=1, c=1, 2 = 2 = 2
a=0, b=0, c=0, 0 = 0 = 0
```

### Polygonal Trianglar Square Number
We attempt to extend the above argument to the equation $P_3(a) = P_4(b) = P_s(c) \iff a^2+a = 2b^2 = (s-2)c^2 - (s-4)c$ where $s > 6$ is an integer. (Note that if $s=6$, since all hexagonal numbers are triangular numbers, this reduces to finding all [hexagonal square numbers](https://mathworld.wolfram.com/HexagonalSquareNumber.html) )

For $t_1, \cdots ,t_5 \in \mathbb{Z}$, we let $$\begin{cases} X = t_1 a + t_2 \\ Y = t_3b \\ Z = t_4 c + t_5 \end{cases} \implies \begin{cases} X^2 = t_1^2a^2 + 2t_1t_2a + t_2^2 \\ Y^2 = t_3^2 b^2  \\ Z^2= t_4^2 c^2 + 2 t_4t_5c + t_5^2 \end{cases}$$

Using the relation $a^2+a=2b^2 \implies b^2 = (a^2+a)/2$, we have:
$$\begin{aligned} X^2 - kY^2 &=  t_1^2a^2 + 2t_1t_2a + t_2^2 -k t_3^2b^2 \\&= t_1^2a^2 + 2t_1t_2a + t_2^2 - k t_3^2(a^2+a)/2  \\ &= \underbrace{ (t_1^2 - kt_3^2/2)}_{=0}a^2 + \underbrace{ (2t_1t_2-kt_3^2/2) }_{=0}a + t_3^2 \end{aligned}$$

and similarly, using the relation $(s-2)c^2 - (s-4)c = 2b^2 \implies b^2 = ((s-2)c^2 - (s-4)c)/2$, we have:
$$\begin{aligned} Z^2 - l Y^2 &= t_4^2 c^2 + 2 t_4t_5c + t_5^2 -lt_3^2b^2 \\ &= t_4^2 c^2 + 2 t_4 t_5 c + t_5^2 - lt_3^2\left(\frac{s-2}{2} c^2 -  \frac{s-4}{2} c\right) \\&= \underbrace{ \left( t_4^2 - lt_3^2(s-2)/2\right)}_{=0} c^2 + \underbrace{ \left(2t_4t_5 -lt_3^2(s-4)/2\right)}_{=0} c + t_5^2 \end{aligned}$$

We therefore have the polynomial system $\begin{cases} 2 t_1^2 - k t_3^2 = 0 \\ 4 t_1 t_2 - k t_3^2 = 0 \\ 2t_4^2 - l t_3^2 (s-2) = 0 \\ 4 t_4 t_5 - l t_3^2 (s-4) = 0 \end{cases}$


