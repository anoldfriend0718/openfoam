import numpy as np

class AMDomainReader(object):
    def __init__(self) -> None:
        super().__init__()
        self._skiprows=15
    
    def read(self,structure_file,nx,ny,nz):
        domain=np.loadtxt(structure_file,skiprows=self._skiprows)
        domain=domain.reshape(nz,ny,nx)
        domain=domain.T
        return domain
