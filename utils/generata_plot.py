#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

x = np.random.randn(10, 3)

fig, ax = plt.subplots()
colors = ['red', 'tan', 'lime']
ax.hist([[2,5,3], [2,12,15], [2, 7, 15]], bins=25,  histtype='bar', color=colors )
plt.xlabel('Impostor scores')
plt.ylabel('Genuine scores')
plt.title('YD_classes')
plt.legend = ax.legend(loc='lower left', shadow=True)

plt.show()

fig.savefig('myimage.svg', format='svg')