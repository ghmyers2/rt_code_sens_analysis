from SALib.analyze import sobol
import numpy as np
from SALib.sample import saltelli
from SALib.analyze import sobol
import numpy as np


if __name__ == "__main__":
    problem = {
        'num_vars': 11,
        'names': ['ar_1', 'd_1', 'f_1', 'Reff_1', 'ar_2', 'd_2', 'f_2', 'Reff_2', 'thick_top', 'dens_top', 'dens_bottom'],
        'bounds': [[0.037, 26.738],
                   [0.05, 0.7],
                   [0.0, 0.9999],
                   [50.0, 2560.0],
                   [0.037, 26.738],
                   [0.05, 0.7],
                   [0.0, 0.9999],
                   [50.0, 2560.0],
                   [0.00001, 0.05],
                   [0.1, 0.5],
                   [0.1, 0.5]]
    }
    
    # , 'soot_1', 'soot_2', 'aod'
    # [0.01, 100.00], [0.01, 100.00], [0.0, 0.4]
    
    num_samples = 2048 # 2^(11)
    param_values = saltelli.sample(problem, num_samples)
    dir = '/home/accurt/Shared/SOBOL/'
    
    np.savetxt(dir + "param_values_SOBOL.txt", param_values)
