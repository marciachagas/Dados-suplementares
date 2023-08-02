import sys, re
import numpy as np
import matplotlib.pyplot as plt

# check stationary genes in boolean evolution for last n time steps:
nsteps=1000
if len(sys.argv) < 2:
  print("Call this script parsing the file containing the trajectories calculated with ASSA-PBN")

file = sys.argv[1]

evolution = []
data = open(file, 'r')
for d in data:
    d = d.strip()
    state = re.split('\t', d)
    evolution.append(state)
data.close()
evolution.pop()
evolution.reverse()

if nsteps >= len(evolution):
    sys.exit("Number of steps is not valid, it is larger than the total time!")

genes = evolution.pop()
genes.pop()
ngenes = len(genes)

#Counter initialization
counter =  []
for j in range(ngenes): 
    counter.append(0)

# Counter incrementation
for t in range(nsteps): 
    for g in range(ngenes):
        if evolution[t][g] == evolution[t + 1][g]:
            counter[g] += 1
        else:
            counter[g] = 0 
    #print("#############################")
    #print(evolution[t], "   ", counter)
    #print(evolution[t+1])

# Normalization of the counter variable:
counter = [ c / nsteps for c in counter]

#removing conserved genes from file conserved_genes.csv
cons_file = open("conserved_genes.csv", "r")
conserved = []
for l in cons_file:
    conserved.append(l.strip())

dynamic_genes={}

for i in range(len(genes)):
    if genes[i] not in conserved:
            dynamic_genes[genes[i]]  = counter[i]
    

out = open("reached_basin.csv", "w")
for g in dynamic_genes.keys():        
    print(g, "\t", dynamic_genes[g], file = out)
out.close()



#generate a bar plot
y_pos = np.arange(len(dynamic_genes.keys()))

# Create bars
plt.bar(y_pos, dynamic_genes.values() )

# Create names on the x-axis
#plt.xticks(y_pos, dynamic_genes)
text = "Conserved genes in last " + str(nsteps) + " steps from total = " + str(len(evolution) -1)  
plt.text(0.2, 1.1, text, horizontalalignment='left', verticalalignment='top')
# Show graphic
plt.show()

exit()

fig = plt.figure()
ax = fig.add_axes([0,0,1,1])
#langs = ['C', 'C++', 'Java', 'Python', 'PHP']
#students = [23,17,35,29,12]
#ax.bar(langs,students)
ax.bar(dynamic_genes.keys(), dynamic_genes.values() )

#n= range(len(dynamic_genes))
#print(n)
plt.show()
exit()
