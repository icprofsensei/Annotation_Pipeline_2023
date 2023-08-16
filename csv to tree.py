
list = [['a', 'b', 'c', 'd', 'z'],['e', 'b', 'f', 'd', 'z'], ['g', 'i', 'f', 'd', 'z'], ['l', 'm', 'n', 'o', 'z'], ['p', 'b', 'f', 'd', 'z'], ['q', 'r', 's', 't', 'z'], ['v', 'w', 's', 't', 'z'], ['x', 'y', 'n', 'o', 'z'],['i', 'o', 'z'],['f', 'd', 'z'], ['k', 't', 'z']]
longest = len(max(list, key = len))
node = max(list, key = len)[len(max(list, key = len)) - 1]

for i in list: 
    nodeset = 'fixed'
    i = i[::-1]
    if len(i) == longest:
        if i[0] != node[0]:   
            print(i[0])
            nodeset == 'unfixed'
if nodeset == 'fixed':
    transitions = []
    reversed_list = [i[::-1] for i in list if len(i) == longest]
    for i in range(longest):
        if i <= longest - 2:
            transls = [rl[i] + ' --> ' + rl[i + 1] for rl in reversed_list]
            transls = set(transls)
            transitions.append(transls)
    print(transitions)

 
    height = len(transitions[0])
    items = len(transitions)
    allnodes = []
    for j in range(0, items + 1):   
        if j <= items -1: 
            nodes = []
            for i in transitions[j]:
                node = i.split(" ")[0]
                nodes.append(node)
        elif j == items:
            nodes = []
            for i in transitions[j - 1]:
                node = i.split(" ")[2]
                nodes.append(node)
        allnodes.append(set(nodes))
    print(allnodes)

