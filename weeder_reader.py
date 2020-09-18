
def weederpwm(weederoutput):
    motifs = []
    nt = ['A', 'G', "C", "T"]
    with open(weederoutput) as wo:
        while True:
            line = wo.readline()
            if line == '': break
            while True:
                if line == '': break
                line = wo.readline()
                if line.startswith('>MAT'):
                    lineA = wo.readline()
                    lineC = wo.readline()
                    lineG = wo.readline()
                    lineT = wo.readline()
                    m = []
                    lA = lineA.split()
                    lC = lineC.split()
                    lG = lineG.split()
                    lT = lineT.split()
                    size = len(lA) - 1
                    for i in range(size):
                        m.append({})
                    for i in range(size):
                        m[i]["A"] = float(lA[i + 1])
                        m[i]["C"] = float(lC[i + 1])
                        m[i]["G"] = float(lG[i + 1])
                        m[i]["T"] = float(lT[i + 1])

                    motifs.append(m)
    return motifs




