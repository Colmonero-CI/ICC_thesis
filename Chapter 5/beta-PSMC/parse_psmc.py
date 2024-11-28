import sys
import matplotlib.pyplot as mpl

# file as first input
file = sys.argv[1]

# if we want, take name (bootstrap or main) and iteratio and sample as second and third and fourth input
if len(sys.argv) > 2:
	method = sys.argv[2]
	name = sys.argv[3]
	iteration = sys.argv[4]
	sample = sys.argv[5]
else:
	name = "na"
	iteration = "na"
	sample = "na"

# open file
f = open(file, "r")

# dictionary for storing the results
results = []
r = None
# read file line by line
for line in f.readlines():
	l = line.strip().split()
	if l[0] == 'LK':
		#print(l[0] + ": " + l[1])
		if r:
			results.append(r)
		r = {'log_lik': float(l[1]), 'results': {'k': [], 't_k': [], 'lambda_k': [], 'pi_k': [], 'sum_l_not=k_A_kl': [], 'A_kk': [], 'method': method, 'name': name, 'iteration': iteration, 'theta_0': None, 'rho_0': None}}
		# start a new iteration of results for appending to the main dict
	if l[0] == 'TR':
		theta_0 = float(l[1])
		rho_0 = float(l[2])
	if l[0] == 'RS':
		r['results']['k'].append(l[1])
		r['results']['t_k'].append(l[2])
		r['results']['lambda_k'].append(l[3])
		r['results']['pi_k'].append(l[4])
		r['results']['sum_l_not=k_A_kl'].append(l[5])
		r['results']['A_kk'].append(l[6])
		r['results']['method'] = method
		r['results']['name'] = name
		r['results']['iteration'] = iteration
		r['results']['sample'] = sample
		r['results']['theta_0'] = theta_0
		r['results']['rho_0'] = rho_0
	
# append the last iteration
results.append(r)

# check if the last iteration has the highest log likelihood
if results[-1]['log_lik'] != max([r['log_lik'] for r in results if r['log_lik'] != 0]):
	print("Warning: the last iteration does not have the highest log likelihood!")

# write results to file
f_out = open(file + ".parsed.txt", "w")

# write header
f_out.write("k\tt_k\tlambda_k\tpi_k\tsum_l_not=k_A_kl\tA_kk\ttheta_0\trho_0\tmethod\tname\titeration\tsample\n")

# write results
r = results[-1]['results']
for i in range(len(r['k'])):
	f_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(r['k'][i], r['t_k'][i], r['lambda_k'][i], r['pi_k'][i], r['sum_l_not=k_A_kl'][i], r['A_kk'][i], r['theta_0'], r['rho_0'], r['method'], r['name'], r['iteration'], r['sample']))
