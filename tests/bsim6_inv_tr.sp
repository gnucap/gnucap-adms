*Sample netlist for BSIM6
.load lang_adms.so
*Inverter Transient
*.option abstol=1e-4 reltol=1e-3 post ingold
*.option trtol=7
.option nobypass noincmode
.verilog

attach ./BSIM6.1.1.so
`include "modelcard.nmos"
`include "modelcard.pmos"
spice
* --- Voltage Sources ---
vdd   supply  0 dc=1.0
vin  vi  0 dc=0.5 sin (0.5 0.5 1MEG)

* --- Inverter Subcircuit ---
.subckt inverter vin vout vdd gnd
    Xp1 vout vin vdd gnd pmos W=10u L=10u 
    Xn1 vout vin gnd gnd nmos W=10u L=10u 
.ends

* --- Inverter ---
Xinv1  vi vo supply 0 inverter

* --- Transient Analysis ---

.print tran v(vi) v(vo)
+hidden(0) iter(0)
.tran 0 5u 5u basic trace=a

.status notime
.end
