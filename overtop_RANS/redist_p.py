from proteus.default_p import *
from proteus.mprans import RDLS
from proteus import Context
from proteus import WaveTools as wt

ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions


genMesh = mesh.genMesh
movingDomain = ct.movingDomain
T = ct.T
mwl = ct.waterLevel #water_level
x_dry = ct.x_dry

# generate a solitary wave
waves = wt.SolitaryWave(waveHeight=ct.height,
			mwl=ct.mwl,
                    	depth=ct.depth,
                   	g=np.array([0., -9.81, 0.]),
                   	waveDir=ct.direction,
                        #trans = np.zeros(3,"d"),
			trans = np.array([ct.x0, 0., 0.]),
                    	fast = ct.fast)

"""
The redistancing equation in the sloshbox test problem.
"""

LevelModelType = RDLS.LevelModel

coefficients = RDLS.Coefficients(applyRedistancing=ct.applyRedistancing,
                                 epsFact=ct.epsFact_redistance,
                                 nModelId=int(ct.movingDomain)+2,
                                 rdModelId=int(ct.movingDomain)+3,
                                 useMetrics=ct.useMetrics,
                                 backgroundDiffusionFactor=ct.backgroundDiffusionFactor)

def getDBC_rd(x, flag):
    pass

dirichletConditions     = {0: getDBC_rd}
weakDirichletConditions = {0: RDLS.setZeroLSweakDirichletBCsSimple}

advectiveFluxBoundaryConditions = {}
diffusiveFluxBoundaryConditions = {0: {}}


class PHI_IC:
    def uOfXT(self, x, t):
        phi_x = x[nd-2] - x_dry
        if phi_x < 0.0:
            phi_y = x[nd-1] - mwl - waves.eta(x,0)
            if phi_y < 0.0:
                return max(phi_x, phi_y)
            else:
                return phi_y
        else:
            phi_y = x[nd-1] - mwl
            if phi_y < 0.0:
                return phi_x
            else:
                return (phi_x ** 2 + phi_y ** 2)**0.5

#class PHI_IC:
#    def uOfXT(self, x, t):
#        if x[nd-2]<=x_dry:
#            return x[nd-1] - ct.waterLevel - waves.eta(x,0)
#        else:
#            return x[nd-1]

#class PHI_IC:
#    def uOfXT(self, x, t):
#        return x[nd-1] - ct.waterLevel - waves.eta(x,0)

#class PHI_IC:
#    def uOfXT(self, x, t):
#        return x[nd-1] - ct.waterLevel

initialConditions  = {0: PHI_IC()}
