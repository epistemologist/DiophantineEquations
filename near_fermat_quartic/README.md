We attempt to find solutions to the equation $x^4 + y^4 + 1 = z^2$ (see [here](https://mathoverflow.net/questions/61794/the-diophantine-eq-x4-y4-1-z2/85754#85754) for initial discussion

Source gives approach which we take from [here](https://mathoverflow.net/a/85754), but briefly, mod 20 analysis yields that $10|x, 10|y$; we therefore have the following procedure:
 - let $y = 10t$, then loop over all divisors $a | y^4+1$
 - let $b = (y^4+1)/a$, we set $a = z-x^2, b = z+x^2 \implies z=(b+a)/2, x^2 = (b-a)/2$
 - if $x^2$ calculated above is actually a square, we have found a solution

We use [this library](github.com/hurchalla/factoring) to calculate prime factorization/divisors for 64-128 bit integers and `openmp` for parallelization - compile the program with `g++ -std="c++17" -O2 -Wall -Wextra -DNDEBUG -I./factoring_headers/include -fopenmp` (including headers to aforementioned library, making sure to not set `DEBUG` for performance reasons)

On my machine, I was unable to find any solutions with $y < 10^9$ in 117 core hours (~16 hours real time). 
