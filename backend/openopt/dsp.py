from openopt import DSP
import sys
import networkx as nx

G = nx.Graph()

network_file = open(sys.argv[1])
solver = sys.argv[2]
output = open(sys.argv[1].replace(".txt", "_dsp_" + solver + ".txt"), "w")


for line in network_file:
    split_line = line.split()
    p1 = split_line[0]
    p2 = split_line[1]
    if (p1 == "FROM" and p2 == "TO"):
        continue
    G.add_edge(p1,p2)

p = DSP(G)
r = p.solve(solver,iprint=0)


for i in r.solution:
    output.write(i + "\n")

output.close()
network_file.close()



