import numpy as np
import matplotlib.pyplot as plt

data = []
counter = 0

f = open("water.out", "r")

for i in f:
    if counter == 0:
        pass
    else:
        tokens = i.strip().split()
        tempt = float(tokens[7])
        if counter % 100 == 0:
            step = int(tokens[0])
            data.append([step, tempt])
    counter += 1

data = np.array(data)
data = data.transpose()

plt.plot(data[0], data[1])
plt.show()

np.savetxt("tempt", data)
           
        

