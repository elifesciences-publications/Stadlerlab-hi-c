import math

size=50
outfilename = "test.txt"

E = {}

E[2] = {}
E[2][0] = 0
E[2][1] = 1
E[2][2] = 0

for t in range(4,size,2):
    E[t] = {}
    E[t][0] = 0
    E[t][t] = 0
    E[t][1] = 1
    E[t][t-1] = 1
    
    for a in range(2,t-1):
        b = t - a
        p = math.factorial(t) / (math.factorial(a) * math.factorial(b))
        
        E[t][a] = (math.factorial(t-2) / math.factorial(t)) *    (a * (a-1) * E[t-2][a-2] + 2 * a * b * (1 + E[t-2][a-1]) +  b * (b-1) * E[t-2][a])
        print(str(E[t][a]))

outfile = open(outfilename, 'w')
for i in range(4,size,2):
	outfile.write('row' + str(i))
	for j in range(2,t-1):
		#i_new = min(i,j)
		#j_new = max(i,j)
		outfile.write('\t' + str(E[i][j]))
	outfile.write('\n')
outfile.close()
