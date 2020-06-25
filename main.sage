# lol haven't tested on non-squarefree N
# (edit: someone has tested on non-squarefree N, and it seems to work.
#  yay me?)
# if you see any bugs let me know or something i guess

# implementation of algorithm detailed at section 4.2 of
# https://www.math.leidenuniv.nl/scripties/Broker.pdf,
# for redpwnctf 2019

# note:
# when i got a solution for d = 3, so the root of the hilbert
# class polynomial was trivial, stuff didn't work... not sure if
# this is an issue that arises with the particular thingo that
# the author uses to simplify computation of roots. maybe
# a closer (read: any) reading of this paper reveals why this
# happens, and i just missed something.
# anyways, it's resolved by just iterating higher :^)

# this impl is way shorter if you don't worry about non squarefree
# things, which i didn't originally, but since i'm doing it again
# might as well do things right yknow

def four_point_two(N):
	f = factor(N)

	# for the square divisors of N
	Ns = 1
	for p, e in f:
		Ns *= p^(e//2)
	S = Ns.divisors()

	d = 0
	while True:
		# fifth step: wait, five is larger than two?
		d += 1

		# second step
		if not is_squarefree(d):
			continue

		# preprocessing stuffs
		K.<a> = QuadraticField(-d)
		Zw = K.ring_of_integers()
		w = a
		if (-d) % 4 == 1:
			w = (-1 + a) / 2

		# the units in Z[w]
		U = [Zw(1), Zw(-1)]
		if d == 1:
			U += [w, -w]
		if d == 3:
			U += [w, -w, w*w, -w*w]

		# third step: splitting behaviour of primes
		# k1 and k2 are always 1 for squarefree N,
		# which is probably what you have, but whatever
		k1 = 1
		D = K.discriminant()
		fail = False
		for p, e in f:
			l = kronecker(D, p)
			if l == -1:
				# inert
				if e % 2 == 1:
					fail = True
					break
				k1 *= p^(e//2)
			if l == 0:
				# ramifies
				k1 *= p^(e//2)

		if fail:
			continue

		N1 = N // k1^2
		ZN1 = ZZ.quo(ZZ*N1)

		# fourth step: stuff
		# r is a root of f, found via quad formula
		for r in ZN1(-d).sqrt(all=True):
			if (-d) % 4 == 1:
				r = (-1 + r) / 2

			# square divisors (trivial if N squarefree)
			for k2 in S:
				k = k1 * k2
				N0 = N1 // k2^2

				# sage inbuilts ;)
				I = Ideal(N0, w - Zw(r))
				G = I.gens_reduced()
				if len(G) > 1:
					# couldn't find a single generator
					continue

				# iterate over "all generators"
				# i.e. chosen generator times units
				for g in G:
					for u in U:
						alpha0 = g * u
						alpha = k*alpha0
						# casting to integer is important lol
						n = ZZ((1-alpha).norm())
						if is_prime(n):
							yield d, alpha

# just my value, put whatevs here
N = 648000029051172520969507317226754839528660854785798909825348737
F = four_point_two(N)

# i'm doing this because sometimes d = 3 dies, like i mentioned
while True:
	d, alpha = next(F)
	print("trying", d, alpha)

	p = ZZ((1-alpha).norm())
	K.<a> = QuadraticField(-d)

	# CM algorithm, also detailed in the paper at blah blah dot blah
	# somewhere in 3 or something

	R = hilbert_class_polynomial(K.discriminant()).roots(GF(p))
	if not R:
		print("failure :(")
		continue
	j = R[0][0]

	if j == 0:
		E = EllipticCurve(GF(p), [0, 1])
	elif j == 1728:
		E = EllipticCurve(GF(p), [1, 0])
	else:
		a = 27*j / (4 * (1728 - j))
		E = EllipticCurve(GF(p), [a, -a])
	
	if E.order() == N:
		print("success!", E)
		break
	E = E.quadratic_twist()
	if E.order() == N:
		print("success!", E)
		break

	# rip, what a shame
	# this isn't really meant to happen as far as i can tell
	# but life is like that sometimes
	print("failure :(")
