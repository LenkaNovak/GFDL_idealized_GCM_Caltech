# exec(open("uf_analytical.py").read())

#modules
import numpy as np
import matplotlib.pyplot as plt

# GFDL FORMULATION


zeta = np.arange(-10,50,1)

def phi_m(zeta, r_crit = 2.0):
    b_stab = 1. / r_crit
    if zeta <= 0:
        phi_m = (1.-16. * zeta)**(-1/2)
    else:
        phi_m = 1. + zeta * (5. + b_stab * zeta)/(1. + zeta)
    return phi_m

def phi_h(zeta, r_crit = 2.0):
    b_stab = 1. / r_crit
    if zeta <= 0:
        phi_h = (1. - 16. * zeta)**(-1/2)
    else:
        phi_h = 1. + zeta * (5. + b_stab * zeta)/(1. + zeta)
    return phi_h


phi_m_c = [phi_m(z) for z in zeta]
r_crit = 8.
phi_m_8 = [phi_m(z, r_crit) for z in zeta]
r_crit = 40.
phi_m_40 = [phi_m(z, r_crit) for z in zeta]

phi_h_c = [phi_h(z) for z in zeta]
r_crit = 8.
phi_h_8 = [phi_h(z, r_crit) for z in zeta]
r_crit = 40.
phi_h_40 = [phi_h(z, r_crit) for z in zeta]

plt.figure()
ax = plt.subplot(121)
plt.plot(zeta, phi_m_8, label = 'phi_m_8')
plt.plot(zeta, phi_m_40, label = 'phi_m_40')
plt.plot(zeta, phi_m_c, label = 'phi_m_c')
plt.ylabel('phi_m')
plt.xlabel('zeta')
plt.legend()
plt.subplot(122)
plt.plot(zeta, phi_h_8, label = 'phi_h_8')
plt.plot(zeta, phi_h_40, label = 'phi_h_40')
plt.plot(zeta, phi_h_c, label = 'phi_h_c')
plt.legend()
plt.savefig('phi_m_phi_h.png')


plt.figure()
ax = plt.subplot(121)
plt.plot(zeta / phi_m_8, phi_m_8, label = 'phi_m_8')
plt.plot(zeta / phi_m_40, phi_m_40, label = 'phi_m_40')
plt.plot(zeta / phi_m_c, phi_m_c, label = 'phi_m_c')
plt.ylabel('phi_m')
plt.xlabel('ri_b')
plt.legend()
plt.xlim(-10, 10)
plt.savefig('phi_m_phi_h_ri.png')