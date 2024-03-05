from .mtool import LPModes
from .viewer import ModeShow
import matplotlib.pyplot as plt
import numpy as np

lps = LPModes(v = 2, d=10)
# Sort the uset items by value in accending order
sorted_uset = sorted(lps.uset.items(), key=lambda x: x[1])
for key, value in sorted_uset:
    print(f"{key}: {value:.4f}")
print(sorted_uset)

ms = ModeShow(lps)
lp01 = lps.LP(0, 1, rot=90, jones=(1/np.sqrt(2), 1j/np.sqrt(2)))

fig, axe = plt.subplots()
axe.set_aspect("equal")
ms.plot2d(fig, axe, lp01, core='on', tickoff=True, vector='on')
plt.show()
