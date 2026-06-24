"""
SIMULACION TPI - Control de Temperatura en Incubadora Neonatal
Teoria de Control - UTN FRBA - Civale
Modelo de 3 nodos termicos + PID sintonizado por Ziegler-Nichols (curva de reaccion)
Genera: curva de reaccion, extraccion K/tau/L, y respuesta de lazo cerrado.
"""
import numpy as np
import matplotlib.pyplot as plt

# ====================================================================
# 1) PARAMETROS FISICOS DEL CASO NOMINAL
# ====================================================================
Cair  = 2289.0     # capacidad termica aire interno efectivo [J/°C]
Cstr  = 8000.0     # estructura + colchon + cupula [J/°C]
Cbaby = 5250.0     # neonato 1.5 kg x 3500 J/(kg°C) [J/°C]
Ras   = 0.020      # R aire-estructura [°C/W] (forzador, acople fuerte)
Rba   = 0.090      # R aire-bebe [°C/W]
Renv  = 0.180      # R aire-sala [°C/W] (perdida principal)
Rstr  = 0.400      # R estructura-sala [°C/W]
Qmet  = 4.0        # calor metabolico neonato [W]
Pmax  = 250.0      # potencia maxima calefactor [W]
eta   = 1.0
Tamb_nom = 24.0    # sala nominal [°C]
L_d   = 8.0        # retardo medicion+proc+distribucion [s]
DT    = 0.5        # paso de integracion [s]

def derivadas(Ta,Ts,Tb,u,Tamb,Qp,Qevap):
    Qh = eta*Pmax*u
    dTa=(Qh+(Ts-Ta)/Ras+(Tb-Ta)/Rba-(Ta-Tamb)/Renv-Qp)/Cair
    dTs=((Ta-Ts)/Ras-(Ts-Tamb)/Rstr)/Cstr
    dTb=(Qmet-(Tb-Ta)/Rba-Qevap)/Cbaby
    return dTa,dTs,dTb

# ====================================================================
# 2) CURVA DE REACCION (lazo abierto) -> K, tau, L
# ====================================================================
def curva_reaccion(u_step=0.5, T_END=4000.0):
    n=int(T_END/DT)
    Ta=Ts=Tb=Tamb_nom
    nd=int(L_d/DT); buf=[0.0]*nd
    t=[]; y=[]
    for i in range(n):
        buf.append(u_step); ua=buf.pop(0)
        dTa,dTs,dTb=derivadas(Ta,Ts,Tb,ua,Tamb_nom,0,0.8)
        Ta+=dTa*DT;Ts+=dTs*DT;Tb+=dTb*DT
        t.append(i*DT); y.append(Ta)
    return np.array(t),np.array(y),u_step

t_cr,y_cr,u_step=curva_reaccion()
T0,Tinf=y_cr[0],y_cr[-1]
K=(Tinf-T0)/(u_step*100)              # °C por % de potencia
dydt=np.gradient(y_cr,t_cr); im=np.argmax(dydt); slope=dydt[im]
L=t_cr[im]-(y_cr[im]-T0)/slope
tau=(t_cr[im]+(Tinf-y_cr[im])/slope)-L

# ====================================================================
# 3) ZIEGLER-NICHOLS (curva de reaccion) -> Kp, Ki, Kd
# ====================================================================
Kp_zn=1.2*tau/(K*L); Ti=2.0*L; Td=0.5*L
Ki_zn=Kp_zn/Ti; Kd_zn=Kp_zn*Td

print("="*60)
print("IDENTIFICACION DE LA PLANTA (curva de reaccion)")
print("="*60)
print(f"K   = {K:.4f} °C/%   (ganancia estatica)")
print(f"tau = {tau:.1f} s     (constante de tiempo)")
print(f"L   = {L:.1f} s        (retardo)")
print(f"L/tau = {L/tau:.3f}  (<0.3 -> facilmente controlable)")
print("\nZIEGLER-NICHOLS (PID):")
print(f"Kp = {Kp_zn:.2f}")
print(f"Ki = {Ki_zn:.3f}   (Ti={Ti:.2f}s)")
print(f"Kd = {Kd_zn:.2f}   (Td={Td:.2f}s)")

# ====================================================================
# 4) LAZO CERRADO con escenarios del informe
# ====================================================================
def lazo_cerrado(Kp,Ki,Kd,T_END=4000.0):
    n=int(T_END/DT)
    Ta=Ts=Tb=Tamb_nom
    integ=0.0; pe=0.0
    nd=int(L_d/DT); buf=[0.0]*nd
    R={k:[] for k in['t','Ta','Tref','u','e','Tb','Ts','perturb']}
    for i in range(n):
        t=i*DT
        Tref=35.0 if t<200 else 36.5
        Qp=12.0 if 1500<=t<1560 else 0.0
        Tamb=22.0 if 2300<=t<3000 else Tamb_nom
        Qevap=2.0 if t>=3300 else 0.8
        pert=0
        if Qp>0:pert=1
        elif Tamb<Tamb_nom:pert=2
        elif Qevap>1:pert=3
        e=Tref-Ta; d=(e-pe)/DT; pe=e
        it=integ+e*DT
        ur=Kp*e+Ki*it+Kd*d
        us=np.clip(ur/100.0,0,1)
        sat=(ur/100.0>1)or(ur/100.0<0)
        if not sat or (e*ur<0): integ=it
        u=us
        buf.append(u); ua=buf.pop(0)
        dTa,dTs,dTb=derivadas(Ta,Ts,Tb,ua,Tamb,Qp,Qevap)
        Ta+=dTa*DT;Ts+=dTs*DT;Tb+=dTb*DT
        for k,v in[('t',t),('Ta',Ta),('Tref',Tref),('u',u*100),('e',e),('Tb',Tb),('Ts',Ts),('perturb',pert)]:
            R[k].append(v)
    return {k:np.array(v) for k,v in R.items()}

r=lazo_cerrado(Kp_zn,Ki_zn,Kd_zn)

# metricas tras cambio de setpoint
m=(r['t']>=200)&(r['t']<1500)
Ta_m=r['Ta'][m]; t_m=r['t'][m]
pico=Ta_m.max(); over=pico-36.5; ess=abs(Ta_m[-1]-36.5)
dentro=np.abs(Ta_m-36.5)<=0.5; ts=None
for i in range(len(dentro)):
    if dentro[i:].all(): ts=t_m[i]-200; break
print("\nDESEMPENO (cambio de setpoint 35 -> 36.5 °C):")
print(f"  Sobreimpulso     = {over:.2f} °C")
print(f"  Error estable    = {ess:.3f} °C")
print(f"  Tiempo establec. = {ts:.0f} s ({ts/60:.1f} min)")

# ====================================================================
# 5) GRAFICOS
# ====================================================================
plt.rcParams.update({'font.size':9,'figure.facecolor':'white'})

# --- Figura A: curva de reaccion con tangente ---
fig1,ax=plt.subplots(figsize=(8,4.5))
ax.plot(t_cr,y_cr,color='#1f6feb',lw=2,label='Respuesta del aire $T_a$ (lazo abierto)')
tang_t=np.array([L,L+tau])
tang_y=np.array([T0,Tinf])
ax.plot(tang_t,tang_y,'r--',lw=1.5,label='Tangente de máxima pendiente')
ax.axhline(T0,color='gray',ls=':',lw=1); ax.axhline(Tinf,color='gray',ls=':',lw=1)
ax.axvline(L,color='green',ls=':',lw=1)
ax.annotate(f'L = {L:.1f}s',(L,T0),textcoords="offset points",xytext=(8,8),color='green')
ax.annotate(f'τ = {tau:.0f}s',(L+tau/2,T0+1),color='red')
ax.set_xlabel('Tiempo [s]'); ax.set_ylabel('Temperatura del aire [°C]')
ax.set_title('Curva de Reacción del Proceso (escalón 50% potencia)\nMétodo de Ziegler-Nichols')
ax.legend(); ax.grid(alpha=0.3); ax.set_xlim(0,2500)
fig1.tight_layout(); fig1.savefig('fig1_curva_reaccion.png',dpi=130)

# --- Figura B: respuesta de lazo cerrado completa ---
fig2,axes=plt.subplots(3,1,figsize=(9,8),sharex=True)
# panel 1: temperaturas
ax=axes[0]
ax.plot(r['t'],r['Tref'],'g--',lw=1.5,label='Referencia $T_{ref}$')
ax.plot(r['t'],r['Ta'],color='#1f6feb',lw=1.8,label='Aire $T_a$ (controlada)')
ax.plot(r['t'],r['Tb'],color='#f0883e',lw=1.2,alpha=0.8,label='Neonato $T_b$')
ax.plot(r['t'],r['Ts'],color='#8957e5',lw=1,alpha=0.6,label='Estructura $T_s$')
ax.fill_between(r['t'],36.0,37.0,color='green',alpha=0.08,label='Banda ±0.5°C')
ax.set_ylabel('Temperatura [°C]'); ax.set_ylim(21,38)
ax.set_title('Respuesta del Lazo Cerrado con PID (Ziegler-Nichols + anti-windup)')
ax.legend(loc='lower right',ncol=2,fontsize=8); ax.grid(alpha=0.3)
# panel 2: potencia
ax=axes[1]
ax.plot(r['t'],r['u'],color='#cf222e',lw=1.2)
ax.fill_between(r['t'],0,r['u'],color='#cf222e',alpha=0.12)
ax.set_ylabel('Potencia u [%]'); ax.set_ylim(0,105); ax.grid(alpha=0.3)
ax.set_title('Señal de control (potencia del calefactor)',fontsize=9)
# panel 3: error
ax=axes[2]
ax.plot(r['t'],r['e'],color='#8250df',lw=1.2)
ax.axhline(0,color='gray',lw=0.8); ax.axhline(0.5,color='green',ls=':',lw=1)
ax.axhline(-0.5,color='green',ls=':',lw=1)
ax.set_ylabel('Error [°C]'); ax.set_xlabel('Tiempo [s]'); ax.grid(alpha=0.3)
ax.set_ylim(-2,2); ax.set_title('Error e(t) = $T_{ref}$ - $T_a$',fontsize=9)
# marcar perturbaciones
for axx in axes:
    axx.axvspan(1500,1560,color='orange',alpha=0.10)
    axx.axvspan(2300,3000,color='cyan',alpha=0.08)
    axx.axvspan(3300,4000,color='magenta',alpha=0.06)
axes[0].annotate('Apertura\nportillo',(1530,33.5),fontsize=7,ha='center',color='darkorange',fontweight='bold')
axes[0].annotate('Sala\n22°C',(2650,25.5),fontsize=7,ha='center',color='teal',fontweight='bold')
axes[0].annotate('+ Evapo-\nración',(3650,25.5),fontsize=7,ha='center',color='purple',fontweight='bold')
fig2.tight_layout(); fig2.savefig('fig2_lazo_cerrado.png',dpi=130)

print("\nGraficos guardados: fig1_curva_reaccion.png, fig2_lazo_cerrado.png")
