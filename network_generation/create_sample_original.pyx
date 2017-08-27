	def create_sample(self):
		Kappa = self.kappa
		N = self.n
		d = self.d
		shortcut_prob = self.shortcut_prob
		boundary = self.boundary

		#Define the variables used:
		cdef double dom_size = 1
		cdef int i
		cdef int j
		cdef double dense
		cdef double r_c
		cdef double[:,:] positions
		cdef int p = 2
		cdef double dij
		cdef double dist

		#Make sure the mean degree is a float so that the next computation works:
		Kappa = float(Kappa)
		# inverted analytic function for r = f(K)
		r_c = (1.0/((3.141592)**0.5) )*(( ((Kappa)/N)*scipy.special.gamma( (d +2.0)/2.0  )   )**(1.0/d ) )

		#Define arrays to store the integer labels of the connected pairs:
		source= [ ]
		target= [ ]

		# N-list of d-lists
		# Randomly generate positions in a d-dimensional standard Normal distribution
		if boundary = 'g':
			positions = np.random.multivariate_normal(np.zeros(d), np.eye(d), size=N)
		#Randomly generate the node positions in the hypercube:
		else:
			positions = np.random.uniform(0, 1.0, (N, d))

		# number of nodes
		cdef int n = positions.shape[0]
		# number of dimensions
		cdef int q = positions.shape[1]

		for i in range(n) :
			for j in range(i+1,n) :

				dij = 0

				#Loop over number of dimensions
				for k in range(q):
					# Compute the absolute distance
					dist = abs( positions[i, k] - positions[j, k] )

					#Extra condition for periodic BCs:
					if boundary == 'p' :
						if dist>0.5*dom_size :
							dist = dom_size - dist
					# Add to the total distance
					dij = dij + dist**2

				dij = dij**0.5

				#Compute the connection probability:
				# returns 1 if within radius, shortcut_prob otherwise
				probability = Top_Hat(dij , r_c , shortcut_prob)

				u = np.random.uniform(0,1.0)

				# for zero shortcut probability, always true in the radius
				# for nonzero baseline, might be connecting long-range shortcuts
				if u < probability :
					# 1/2 probability for direction
					if random() >= 0.5: # i --> j
						source.append(i)
						target.append(j)
					else: # j --> i
						source.append(j)
						target.append(i)

		return RGGSample(
			    source=source,
				target=target,
				positions=np.asarray(positions),
				kappa=Kappa,
				n=N,
				d=d,
				boundary=boundary,
				shortcut_prob=shortcut_prob,
				)