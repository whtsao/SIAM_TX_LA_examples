from proteus.default_p import *
from proteus.mprans import NCLS
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


LevelModelType = NCLS.LevelModel

coefficients = NCLS.Coefficients(V_model=int(ct.movingDomain)+0,
                                 RD_model=int(ct.movingDomain)+3,
                                 ME_model=int(ct.movingDomain)+2,
                                 checkMass=False,
                                 useMetrics=ct.useMetrics,
                                 epsFact=ct.epsFact_consrv_heaviside,
                                 sc_uref=ct.ls_sc_uref,
                                 sc_beta=ct.ls_sc_beta,
                                 movingDomain=ct.movingDomain)

dirichletConditions = {0: lambda x, flag: None}

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

initialConditions = {0: PHI_IC()}
