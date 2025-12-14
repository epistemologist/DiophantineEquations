from sage.all import *
from collections import namedtuple
from itertools import product
from tqdm import tqdm


Solution = namedtuple("Solution", ["a", "b", "c", "m", "n", "l", "common_val"])

def get_sols(a, b, c):
	# We search for solutions to P_a(m) = P_b(n) = P_c(l)
	R = PolynomialRing( QQ, names = ('X', 'Y') )
	(X, Y) = R.gens()
	S = PolynomialRing( ZZ, names = ("m", "n", "l"))
	(m,n,l) = S.gens()
	# Following Thm 2.4 from https://tokushima-u.repo.nii.ac.jp/record/2009587/files/jm_55_1.pdf
	A = a-2; B = b-2; C = c-2;
	E = EllipticCurve(
		Y**2 - X * (X - C*(A-B)*(A*B-4)) * (X - B*(A-C)*(A*C-4))
	)
	possible_m_vals = set()
	for pt in tqdm( E.integral_points(both_signs = True) ):
		X_sol, Y_sol = pt.xy()
		poly1 = -X_sol + B*C*(2*A*m - (A-2))**2
		poly2 = -Y_sol + A*B*C*(2*A*m-(A-2))*(2*B*n-(B-2))*(2*C*l-(C-2))
		# poly1 is univariate quadratic in m; get roots 
		m_vals = poly1.univariate_polynomial().roots(multiplicities=False)
		possible_m_vals |= set(m_vals)
	# iterate over possible m values to get possible n,l values
	P = lambda s, n: ( (s-2)*n**2 - (s-4)*n ) / 2
	solutions = []
	for m_val in possible_m_vals:
		common_val = P(a, m_val)
		# Solve for b,n
		n_vals = ( P(b, n) - common_val ).univariate_polynomial().roots(multiplicities = False)
		# Solve for c,l
		l_vals = ( P(c, l) - common_val ).univariate_polynomial().roots(multiplicities = False)
		solutions.extend([ 
			Solution(a=a,b=b,c=c,m=m_,n=n_,l=l_,common_val=common_val)
			for m_, n_, l_ in product([m_val], n_vals, l_vals) ]
		)
	return [ s for s in solutions if all(i.is_integer() for i in tuple(s))]

out = []
for i in range(7, 200):
	tmp = get_sols(3, 4, i)
	print(tmp)
	out.extend(tmp)
