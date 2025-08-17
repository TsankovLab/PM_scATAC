import numpy as np

# Monkey-patch for deepdish compatibility with NumPy >=1.24
if not hasattr(np, "object"):
    np.object = object

import chrombpnet.CHROMBPNET as chrombpnet

chrombpnet.main()
