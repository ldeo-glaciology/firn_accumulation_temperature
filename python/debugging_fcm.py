 # %% 
# upwind_difference_matrix
def upwind_difference_matrix(z0,zL,n,v):
    from scipy.sparse import spdiags
    import numpy as np
    dz=(zL-z0)/(n-1)
    r1fdz=1/dz

    #     (1)  finite difference approximation for positive v     
    if v > 0:
    #
    #             sparse discretization matrix      
    #
    #             interior points      
        e = np.ones(n)
        D = spdiags(np.stack([-e, e]), np.array([-1, 0]), n, n).toarray()
    # 
    #             boundary point      
        D[0,0:2] = np.array([-1, +1])
            

    #     (2)  finite difference approximation for negative v
    if v < 0:
    #
    #            sparse discretization matrix      
    #
    #             interior points      
        e = np.ones(n)
        D = spdiags(np.stack([-e, e]), np.array([-1, 0]), n, n).toarray()
    #
    #             boundary point      
        D[n-1,(n-2):(n)] = np.array([-1, +1])
                            
    
    return r1fdz*D


# %
# model equations      


def eqns(t, y):
    #print('eqns_test_1')
    from scipy import integrate
    Ly = len(y)
    
    N = p['N']   
    z_h = p['z_h']
    D1 = p['D1']
    lambda_g  = p['lambda_g']
    lambda_c  = p['lambda_c']
    delta = p['delta']
    Ar = p['ArthenNumber']
    w_s = p['w_s']

    ### Collect the simulation variables.
    phi = y[:-1:4]
    r2 = y[1:-1:4]
    T = y[2:-1:4]
    A = y[3:-1:4]
    H = y[-1]

    ### create dydt vector
    dydt = np.empty_like(y)

    ### create views of dydt
    dphidt = dydt[:-1:4]
    dr2dt = dydt[1:-1:4]
    dTdt = dydt[2:-1:4]
    dAdt = dydt[3:-1:4]
    dHdt = dydt[Ly-1:Ly]

    ### Compute gradients.
    dphidz = np.matmul(D1, phi) 
    dr2dz  = np.matmul(D1, r2)  
    dAdz   = np.matmul(D1, A)  

    ### Compute stress.
    s_int = H*(1 - phi) 
    sigma = integrate.cumulative_trapezoid(s_int,z_h, initial=0) 

    ### Compute the ice velocity.
    v_int = -(H/Ar)*sigma**n*phi**m*np.exp(lambda_c*T)/r2 
    w = integrate.cumulative_trapezoid(v_int, z_h, initial=0) + w_s 

    ### Column height.
    dHdt[:] = w[-1] - beta/(1 - phi[-1])    

    ### Change in porosity.
    dphidt[:] = (1/H)*(D1@((1 - phi)*w) + dHdt*z_h*dphidz) 

    ### Change in square of the grain size.
    if p['sim_r']:   # (only if sim_r ==1) 
        dr2dt[:] = (1/H)*(dHdt*z_h - w)*dr2dz + (1 - delta*r2)*np.exp(lambda_g*T) 
    else:
        dr2dt[:] = 0*r2 


    ### Change in temperature.
    dTdt[:] = 0 

    ### Change in age.
    dAdt[:] = 1 + (1/H)*(dHdt*z_h - w)*dAdz 
    #    dAdt = zeros(N,1) 

    ### Upper surface boundary conditions.
    dphidt[0] = 0 
    dr2dt[0] = 0 
    dAdt[0] = 0 

    #### Collected ODE vector.    
    #dVarsdt = np.concatenate(dphidt,  dr2dt,  dAdt,  dTdt,  dHdt)        
            
    return dydt




# %%
# setup + run
#class fs:   # short for firn simulation
# """This is firn model class"""

# def __init__(self,**kwargs):
#     #import numpy as np
#     #self.sim_label = sim_label
#     p = {}

# def fullRun(self)
#    self.setup()
#    self.integrate()

#def setup(self,**kwargs):
print('starting setup')

import numpy as np

p = {}
# 1. define default paramaters   
p['sim_label'] = 'default_label'
p['b0_mpy'] = 0.1
p['beta'] = 1
p['nu'] = 1
p['T_s_dim'] = 253.15
p['plotting'] = 1
p['saving_xt'] = 1
#        p['plotting_period'] = 3000
#        p['sampling_period'] = 20
p['dz'] = 0.01
#        p['t_total'] = 20
p['save_dir'] = 'results_scratch'
p['save'] = 1
p['sim_T'] = True
p['sim_r'] = True
p['PauseGrainEvolution_t'] = np.nan
p['z0'] = 100
p['r_s_dim'] = 2.5000e-07   # grain size at the surface (0.5 mm)**2 
p['phi_s'] = 0.5   
p['n'] = 1
p['simDuration'] = 10
p['scaleDuration'] = False

# 2. Replace any parameters that are defined by the user (in kwargs) 
# for key, values in kwargs.items():
#     p[key] = values     
    
# 3. Define (or extract from p) the dimensional parameters of the system.
b0_mpy = p['b0_mpy']        # ice equivalent accumulation rate [m / yr]
T_s_dim = p['T_s_dim']      # upper surface temperature [K]
z_0 = p['z0']               # initial column height [m]
dz_dim = p['dz']*z_0        # dimensional numerical grid spacing [m]
r2_s_dim = p['r_s_dim']     # upper surface grain size [m**2] (0.5 mm)**2 
r2_f = 0.01**2                        # maximum grain size [m**2] (1 cm)**2 
phi_s = p['phi_s']          # upper surface porosity
c_i = 2009                            # heat capacity [J / (kg K)]
E_c = 60e3                            # compaction activation energy [J]
E_g = 42e3                            # grain growth activation energy [J]
k_c = 9.2e-9                          # flow parameter [kg?1 m3 s] from Arthern kc
k_g = 1.3e-7/r2_f                     # grain growth rate constant [m**2 / s]
m = 1                                 # porosity exponent
n = p['n']                       # stress exponent
kappa_0 = 2.1                         # thermal conductivity [W / (m K)], # arthern and wingham 1998 pg. 18
G = 0.05                              # geothermal heat flux  [W/m**2]

## 4. Constants, scales and parameters
### 4.1 constants
g = 9.81                           # acceleration due to gravity [m s**-2]
spy = 24*365*3600                  # seconds per year
R = 8.31                           # ideal gas constant
rho_i = 918
p['rho_i'] = rho_i                        # ice density [kg / m**3]

### 4.2 Scaling parameters.
b_0 = b0_mpy/spy 
h_0 = z_0 
r2_0 = (h_0*k_g*r2_f/b_0)*np.exp(-E_g/(R*T_s_dim)) 
t_0 = h_0/b_0 
T_0 = G*z_0/kappa_0 
sigma_0 = g*rho_i*h_0 
w_0 = b_0                 # scale of the vertical velocity is the accumulation rate

# save scales for use in plotting later
p['b_0'] = b_0
p['r2_0'] = r2_0
p['t_0'] = t_0
p['sigma_0'] = sigma_0
p['w_0'] = w_0
p['h_0'] = h_0
p['T_0'] = T_0

### 4.3 Non-dimensional parameters.
lambda_c = E_c*T_0/(R*T_s_dim**2.0) 
lambda_g = E_g*T_0/(R*T_s_dim**2.0) 
gamma = (sigma_0/4.0)**(1-n)    # factor to account for non-linear rheology, under the default value of n = 1, gamma = 1 and it has no effect on the value of Ar
Ar = r2_0/(k_c*t_0*sigma_0*np.exp(-E_c/(R*T_s_dim)))/gamma  # Arthern number
Fl = h_0*G/(kappa_0*T_0) 
Pe = rho_i*c_i*b_0*h_0/kappa_0 
beta = p['beta'] 
delta = r2_0/r2_f 
nu = p['nu']       # accumulation multiplier

## 5. Compute vertical velocity upper boundary conditon
w_s = nu*beta/(1 - phi_s)  

## 6. Set up real space grid in height coordinates.
z_init = np.arange(z_0,-dz_dim,-dz_dim)
N = z_init.size

### Normalized depth coordinates.
z_h = np.flip(z_init)/z_0 
#print(len(z_h))


## 6.1. time span
simDuration = p['simDuration']
beta = p['beta']

if p['scaleDuration']:
    p['t_span'] = np.array([0, simDuration/beta])
else:
    p['t_span']  = np.array([0, simDuration])


## 7. Initial conditions
Ly0 = len(z_h)*4+1

y0 = np.zeros(Ly0,dtype=np.float64)
phi_init = y0[:-1:4]
r2_hat_init = y0[1:-1:4]
T_hat_init = y0[2:-1:4]
A_hat_init = y0[3:-1:4]
H_init = y0[Ly0-1:Ly0] 

#print(len(phi_init))
### Initial porosity
phi_init[:] = (1-z_h)*phi_s 

### Dimensional  grain size squared.
if p['sim_r']:
    r2_hat_init[:] = r2_s_dim/r2_0 + z_h 
else:
    r2_hat_init[:] = r2_s_dim/r2_0 + 0*z_h 

### Dimensionless temperature (mostly not used)
T_hat_init[:] = np.zeros(N) 

### Dimensionless firn age.
A_hat_init[:] = z_h 

### domain height inital condition     
H_init[:] = z_0/h_0

### Initial conditions vector.
#Vars_init = [phi_init  r2_hat_init  A_hat_init  T_hat_init  z_0/h_0] 

## 8. Define gradient operator
### Finite difference gradient operator using three point upwind scheme.
D1 = upwind_difference_matrix(z_h[0], z_h[-1], N, 1) 

## 9. Save outputs for use at the solution stage

# 10. Add the newly created values to p
p['N'] = N 
p['D1'] = D1
p['spy'] = spy 
p['lambda_c'] = lambda_c 
p['lambda_g'] = lambda_g 
p['delta'] = delta 
p['PecletNumber'] = Pe 
p['FluxNumber'] = Fl 
p['ArthenNumber'] = Ar 
p['y0'] = y0
p['z_h'] = z_h
p['w_s'] = w_s


print('finished setup!')


# 

# integrateun      


#def integrate(self):
from scipy import integrate
from scipy import interpolate

### Solve the system

#options = odeset('Events', @SteadyState, 'RelTol',1e-8,'AbsTol', 1e-8);
#[Time, Vars] = ode15s(@GovEq, time, Vars_init, options);
# %%
print('starting integration')
sol = integrate.solve_ivp(eqns, p['t_span'], p['y0']);    # need to change back to self.eqns
print('finished integration!')

# 
### Collect the variables.
phi = sol.y[:-1:4]
r2 = sol.y[1:-1:4]
T = sol.y[2:-1:4]
A = sol.y[3:-1:4]
Height = sol.y[-1]
t = sol.t
### compute density from porosity
rho = p['rho_i']*(1 - phi)        
            
### Un-normalized depth coordinates
import numpy.matlib

Hg, zg = np.meshgrid(Height,p['z_h'])
Depth = p['h_0']*Hg*zg

### Compute the firn thickness as a function of time

#z830 = interp1(Phi(:,end),Depth(:,end)/p.h_0,1-830/p.rho_i);  # normalized
# %%
### Compute normalized ice velocity
W = np.empty_like(phi)
M = np.empty_like(Height)
z830 = np.empty_like(Height)
FAC = np.empty_like(Height)
#Mass = nan(numel(Time),1);
for i in range(len(t)):
### Compute stress    
    S_int = Height[i]*(1 - phi[:,i])
    Sigma = integrate.cumulative_trapezoid(S_int,z_h,initial=0)
### Compute velocity    
    W_int = -(Height[i]/Ar)*Sigma**n*phi[:,i]*m*np.exp(lambda_c*T[:,i])/r2[:,i]
    W[:,i] = integrate.cumulative_trapezoid(W_int, z_h, initial=0) + w_s
### Compute total mass in the column.
    M[i] = integrate.trapezoid(1 - phi[:,i], Depth[:,i])
### Compute firn air conten
    FAC[i] = integrate.trapezoid(phi[:,i], Depth[:,i])
### Compute the firn thickness (depth to 830 kg/m^3).
    f = interpolate.interp1d(phi[:,i],Depth[:,i]/p["h_0"])
    z830[i] = f(1-830/p["rho_i"])
# %%

### Create output xarray
import xarray as xr
xr.set_options(display_expand_attrs=False)
# %%
out = xr.Dataset(
    data_vars=dict(
        phi=(["z_h", "t"], phi),
        r2=(["z_h", "t"], r2),
        rho=(["z_h", "t"], rho),
        A=(["z_h", "t"], A),
        T=(["z_h", "t"], T),
        w=(["z_h", "t"], W),
        H=(["t"], Height),
        M=(["t"], M),
        FAC=(["t"], FAC),
        z830=(["t"], z830),
    ),
    coords=dict(
        z_h=(["z_h"], p['z_h']),
        depth=(["z_h", "t"],Depth),
        t=t,
    ),
    attrs=dict(simulation_parameters=p),
)

# add attributes to variables
out.phi.attrs = dict(name="porosity", 
    long_name="nondimensional porosity",
    scale=1,
    scale_units="none")

out.rho.attrs = dict(name="density", 
    long_name="density",
    units="m^3/s")

out.r2.attrs = dict(name="grainsize", 
    long_name="nondimensional squared grain radius",
    scale=p["r2_0"],
    scale_units="m^2")

out.A.attrs = dict(name="age", 
    long_name="nondimensional age",
    scale=p["t_0"],
    scale_units="s")

out.T.attrs = dict(name="temperature", 
    long_name="nondimensional temperature",
    scale="to be added")

out.H.attrs = dict(name="height", 
    long_name="nondimensional domain height",
    scale=p["z0"],
    scale_units="m")

out.w.attrs = dict(name="velocity", 
    long_name="nondimensional velocity",
    scale=p['w_0'],
    scale_units="m/s")

out.M.attrs = dict(name="mass", 
    long_name="total ice-equivelent mass",
    scale_units="m ice-equivelent"
    notes="this is the total mass in the domain expressed as meters ice-equivelent. It is computed as the depth integral of the 1-phi at each time step.")

out.z830.attrs = dict(name="firn air content", 
    long_name="firn air content",
    units="m"
    notes="computed as the depth integral of the porosity. Note that if the firn is not completed compacted by the bottom of the domain (i.e. phi(z_h=1,t) > 0) then this FAC value does not accurately represent the total firn air content, as some air nominally exists beneath the bottom of the domain.")

out.z830.attrs = dict(name="firn thickness", 
    long_name="nondimensional firn thickenss",
    scale=p['z_0'],
    scale_units="m")

# add attributes to coordinates
out.z_h.attrs = dict(name="non_dim_depth",
    long_name="scaled depth",
    scale="H(t) * z_0",
    scale_units="m",
    notes="the scaling of this variable changes with time through the simulation. So z_h = 1 means a different dimensional depth at each time step. Specifically, z_h = 1 corresponds to H(t)*z_0. H(t) is one of the other variables, and z_0 is the dimensional scale for depth [m].")
    
out.t.attrs = dict(name="time", 
    long_name="nondimensional time",
    scale=p["t_0"],
    scale_units="s")

out.depth.attrs = dict(name="depth", 
    long_name="dimensional depth",
    units="m")

# plot
# out.phi.isel(t=-1).plot() # profile of last time slice
# out.phi.hvplot()   # static z-t plot of phi
# out.phi.hvplot(groupby='t') # profile plot with a time-slider
# out.H.hvplot()

### Create output stucture.
#Out = struct('Time', Time, 'Rho', Rho, 'Temp', Temp, 'GrainSize', GrainSize,...
#    'Sigma', Sigma, 'Phi', Phi, 'Age', Age, 'Depth', Depth, 'W', W, 'p', p,...
#    'H', Height, 'Mass', Mass, 'zeta830final',zeta830final);

# return out

   

# %%
