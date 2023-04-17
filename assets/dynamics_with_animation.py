# necham zatizeni ale zvetsim kostrukci na 5 pater

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)}, threshold=np.inf)
# matplotlib.use("Agg")

pocet_pater = 5
E = 3.43 * 10 ** 7  # kN/m2
I = 6.75 * 10 ** -4  # m4
l = 3.1  # m
m = 60  # kNs^2/m = 60 tun vaha kazdeho patra

# zavedeni utlumu
pom_utlum = 0.05
omega1 = 4.75  # s-1
alfa = pom_utlum * omega1
# beta = pom_utlum / omega1
# betu dam rovnou nule, abych to mel snazsi, ale prijde mi to divny, z matice tuhosti by tam sly o rad vetsi cisla

k = 2 * 12 * E * I / l ** 3


# lokalni matice tuhost
K_lok = np.array([[k, -k],
                  [-k, k]])
# globalni matice tuhosti
K_glob = np.zeros([pocet_pater+1, pocet_pater+1])
for i in range(pocet_pater):
    K_glob[i, i] += K_lok[0, 0]
    K_glob[i+1, i] += K_lok[0+1, 0]
    K_glob[i, i+1] += K_lok[0, 0+1]
    K_glob[i+1, i+1] += K_lok[0+1, 0+1]

# smazu spodni patro, tam posuny znam
K_glob = np.delete(K_glob, -1, 0)
K_glob = np.delete(K_glob, -1, 1)
print('Globalni matice tuhosti')
print(K_glob)

# Matice hmotnosti, bude diagonalni matice (hmotnost soustredim na vrcholu prvku)
M_glob = np.zeros([pocet_pater, pocet_pater])
for i in range(pocet_pater):
    M_glob[i, i] = m
print('Globalni matice hmotnosti')
print(M_glob)

# Matice utlumu zavedu jako C = alfa*M
C = alfa * M_glob

# casovy krok
dt = 0.02  # s
# dt = 0.059125  # dela hezkej obrazek

# nahrani akcelerogramu pro el centro
el_centro = np.loadtxt("ElCentro.txt")
# pridam tam 10 s ticho za to
# tady musim zajistit aby to interpolovalo pro pripad ze mam desli krok nez je casovy krok akcelerogramu
if dt != 0.02:
    kroku_org = len(el_centro)
    ted_kroku = int(kroku_org * 0.02/dt)
    print(ted_kroku)
    el_centro_mod = np.zeros([ted_kroku, 2])
    for i in range(ted_kroku):
        cas = dt * i
        el_centro_mod[i, 0] = cas
        el_centro_mod[i, 1] = np.interp(cas, xp=el_centro[:, 0], fp=el_centro[:, 1])
    print(el_centro_mod)
    el_centro = el_centro_mod

pridavek = np.zeros([int(10/dt), 2])
for i in range(len(pridavek)):
    pridavek[i, 0] = el_centro[-1, 0] + dt * (i+1)

el_centro = np.append(el_centro, pridavek, axis=0)

vsechny_kroky = el_centro[:, 0]
a_zakl = el_centro[:, 1]
# celkovy pocet kroku dt
kroku_celkem = len(a_zakl)
print('celkem casovych kroku =', kroku_celkem)

# celkovy cas
t_celkem = kroku_celkem * dt
print('celkovy sledovany cas =',t_celkem, 's')

# zatizeni rozdelim do bodu, jako prislusnou hmotu v matici hmotnosti * zrychleni v danem okamziku
f = np.zeros([kroku_celkem, pocet_pater])


# tady nemam uplne 100% proc to takhle funguje
for i in range(kroku_celkem):
    for j in range(pocet_pater):
        f[i, j] = - M_glob[j, j] * a_zakl[i]
# print('f\n', f)

####################################################################
# vypocet
####################################################################

RHS = 1/(dt**2) * M_glob + 1/(2 * dt) * C

RHS_inv = np.linalg.inv(RHS)

u = np.zeros([kroku_celkem, pocet_pater])
print('RHS inv')
print(RHS_inv)
print(RHS_inv.shape)
for i in range(2, kroku_celkem):
    LHS = f[i-1] - (K_glob - 2/(dt**2) * M_glob) @ u[i-1] - (1/(dt**2) * M_glob - 1/(2*dt) * C) @ u[i-2]
    u[i] = RHS_inv @ LHS

print(u)

pridavek_u = np.zeros([1, len(u)])
print(pridavek_u)

u_celkem = np.zeros([len(u), pocet_pater+1])


for i in range(len(u)):
    for j in range(pocet_pater):
        u_celkem[i,j] = u[i,j]

print(len(u_celkem))
print('maximalni posuny v patrech')
for i in range(pocet_pater):
    print('maximalni posun patra',pocet_pater-i, '=',np.amax(u[:,i]), 'm')
#


# y = np.linspace(15.5, 0, num=6)
y = np.linspace(pocet_pater*3.1, 3.1, num=pocet_pater)
#animace
fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, figsize=(9, 6))

def animate(i):
    # print('frame =', i, '/', len(u_celkem))
    ax1.clear()
    ax1.set_title('relative deflection of floors')
    ax1.set_ylabel('floors (height) [m]')
    ax1.set_xlabel('deflection [mm]')
    ax1.set_xlim(-9, 9)
    # ax1.plot(u_celkem[i], y, 'o')
    ax1.plot(np.array([0, 0]), np.array([0, 15.5]), color='grey')
    ax1.plot(u[i] * 100, y, 'o')

    ax2.clear()
    ax2.set_title('El Centro')
    ax2.set_ylabel('acceleration [m/s2]')
    ax2.set_xlabel('time [s]')
    ax2.plot(vsechny_kroky, a_zakl, color='grey')
    # ax.set_title("Zrychleni v zakladove spare")
    ax2.set_ylim(-0.35, 0.35)
    ukazatel_y = np.array([-0.3, 0.3])
    ukazatel_x = np.array([i*dt, i*dt])

    ax2.plot(ukazatel_x,ukazatel_y, color='red')





ani = FuncAnimation(fig, func=animate, frames=np.arange(0, len(u_celkem)), interval=20, repeat_delay=2000, save_count=2060)
# plt.show()

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=50, metadata=dict(artist='Me'), bitrate=1800)

print('cekej...')
ani.save('finaloutput.mp4', writer=writer)




