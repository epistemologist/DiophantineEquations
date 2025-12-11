#include<iostream>
#include<cstdint>
#include<vector>
#include<map>
#include<set>

#include "hurchalla/factoring/factorize.h"
#include<omp.h>

using namespace std;

typedef uint32_t u32;
typedef uint64_t u64;
typedef unsigned __int128 u128;
typedef __int128 i128; 

typedef map<u128, size_t> Factorization;

// Utility functions to deal with u128 integers
#define U128(LO, HI) ( (((u128) HI ) << 64L ) | LO )
#define U128_HI(N) ( (u64) (N>>64L) )
#define U128_LO(N) ( (u64) (N & 0xffffffffffffffffL) )

// Debug to print u128
void print(u128 x) {
    if (x > 9) print(x / 10);
    putchar(x % 10 + '0');
}

/*
We find solutions to the equation
x^4 + y^4 + 1 = z^2

# mod 20 analysis
In [7]: [(x,y,z) for x,y,z in product(range(20), repeat=3) if (x**4 + y**4 - z**2 + 1) % 20 == 0]
Out[7]:
[(0, 0, 1),
 (0, 0, 9),
 (0, 0, 11),
 (0, 0, 19),
 (0, 10, 1),
 (0, 10, 9),
 (0, 10, 11),
 (0, 10, 19),
 (10, 0, 1),
 (10, 0, 9),
 (10, 0, 11),
 (10, 0, 19),
 (10, 10, 1),
 (10, 10, 9),
 (10, 10, 11),
 (10, 10, 19)]

 => We must have that 10|x, 10|y

See https://mathoverflow.net/a/85754 for approach 
In more detail:
 - let y = 10t, loop over all divisors a | y^4+1
 - let b = (y^4+1)/a, set a = z - x^2, b = z + x^2
 - 		=> z = (b+a)/2, x^2 = (b-a)/2
 - note: we need a <= b = (y^4+1)/a => a^2 <= y^4+1
 - check if x^2 is actually a square
*/
// Search for solutions with y <= Y_MAX
#define Y_MAX 1'000'000'000L
// Chunk size (size of interval to search for on each core)
#define CHUNK_SIZE 200000

set<u128> gen_divisors(u128 N) {
	// Generate factors and put in Factorization object
	unsigned int num_factors;
	auto factors = hurchalla::factorize(N, num_factors);
	Factorization f;
	for (auto p: factors) {
		if (p > 0) f[p] += 1;
	}
	// Generate divisors
	set<u128> curr_divisors = {1};
	set<u128> tmp;
	for (auto pair: f) {
		auto p = pair.first; auto e = pair.second;
		u128 p_to_e = 1;
		for (int i = 0; i <= e; i++) {
			for (auto d: curr_divisors) tmp.insert( d * p_to_e );
			p_to_e *= p;
		}
		curr_divisors.insert( tmp.begin(), tmp.end() );
	}
	return curr_divisors;
}

bool is_square(u128 N) {
	// Only works for N < 2^126ish
	u128 N_orig = N, res = 0;
	int bit_length = U128_HI(N) == 0 ?
	(64 - __builtin_clzll(U128_LO(N))) :
	(128 - __builtin_clzll(U128_HI(N)) );
	u128 bit = ( (u128) 1 ) << (bit_length/2*2);
	while (bit > 0) {
		if (N >= res + bit) {
			N -= res + bit;
			res = (res >> 1) + bit;
		} else {
			res >>= 1;
		}
		bit >>= 2;
	}
	return (res * res == N_orig);
}

bool check_solution(u128 y) {
	// Given a y, checks if there exists a solution to
	// x^4 + y^4 + 1 = z^2
	// Discard the case of trivial solution
	if (y == 0) return false ;
	u128 Y = y*y*y*y+1;
	for (auto a: gen_divisors(Y)) {
		if (a * a > Y) break;
		u128 b = Y / a;
		if (a % 2 != b % 2) continue;
		u128 z = (b+a)/2 ;
		u128 x_squared = (b-a)/2;
		if (is_square(x_squared)) { 
			print(x_squared); cout << endl;		
			return true;
		}
	}
	return false;
}

u128 check_interval(u128 lo, u128 hi) {
	// Checks for solutions to
	// x^4 + y^4 + 1 = z^2
	// where lo <= y < hi
	while (lo % 10 != 0) lo++;
	u128 out = 1; // Return 1 in case of no solution
	#pragma omp parallel for schedule(static)
	for (u128 i = lo; i < hi; i += 10) {
		if ( check_solution( i ) ) out = i; 
	}
	return out; // Failure
}

int main() {
	int max_i = Y_MAX / CHUNK_SIZE;
	for (int i = max_i; i > 0; i--) {
		int lo = (i-1) * CHUNK_SIZE;
		int hi = (i == max_i) ? Y_MAX : i*CHUNK_SIZE;
		// Construct interval (lo , hi)
		fprintf(stderr, "[+] Searching interval [ %d, %d ) for solutions\n", lo, hi);
		u128 sol = check_interval(lo, hi);
		if (sol != 1) {
			cout << "Solution found!!!: " << U128_HI(sol) << " " << U128_LO(sol) << endl;
			exit(0);
		}
	}
	cout << "No solutions found!" << endl;
}
