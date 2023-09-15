import svgutils.compose as sc
from IPython.display import SVG # /!\ note the 'SVG' function also in svgutils.compose
import numpy as np
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, figsize=(4,4))
ax.plot(np.sin(np.linspace(0,2.*np.pi)), np.cos(np.linspace(0,2.*np.pi)), 'k--', lw=2.)
ax.plot(np.random.randn(20)*.3, np.random.randn(20)*.3, 'ro', label='random sampling')
ax.legend()
ax2 = plt.axes([.2, .2, .2, .2])
ax2.bar([0,1], [70,30])
plt.xticks([0.5,1.5], ['water  ', ' ground'])
plt.yticks([0,50])
plt.title('ratio (%)')
fig.savefig('cover.svg', transparent=True)
# here starts the assembling using svgutils 
sc.Figure("8cm", "8cm", 
    SVG("results2/trees/img_tree.svg")
    ).save("compose.svg")
SVG('compose.svg')
