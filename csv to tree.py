
list = [['a', 'b', 'c', 'd', 'z'],['e', 'b', 'f', 'd', 'z'], ['g', 'i', 'f', 'd', 'z'], ['l', 'm', 'n', 'o', 'z'], ['p', 'b', 'f', 'd', 'z'], ['q', 'r', 's', 't', 'z'], ['v', 'w', 's', 't', 'z'], ['x', 'y', 'n', 'o', 'z'],['i', 'o', 'z'],['f', 'd', 'z'], ['k', 't', 'z']]
longest = len(max(list, key = len))
node = max(list, key = len)[len(max(list, key = len)) - 1]
for i in list: 
    node = 'fixed'
    i = i[::-1]
    if len(i) == longest:
        if i[0] != 0:   
            node == 'unfixed'
if node == 'fixed':
    transitions = []
    reversed_list = [i[::-1] for i in list if len(i) == longest]
    for i in range(longest):
        if i <= longest - 2:
            transls = [rl[i] + ' --> ' + rl[i + 1] for rl in reversed_list]
            transls = set(transls)
            transitions.append(transls)
    print(transitions)
    for t in transitions[0]:
        nextnode = t.split(" ")[2]
        for tr in transitions[1]:
            if tr.split(" ")[0] == nextnode:
                nextnextnode = tr.split(" ")[2]
                for tra in transitions[2]:
                    if tra.split(" ")[0] == nextnextnode:
                        nextnextnextnode = tra.split(" ")[2]
                        for tran in transitions[3]:
                            if tran.split(" ")[0] == nextnextnextnode:
                                print(t, tr, tra, tran)