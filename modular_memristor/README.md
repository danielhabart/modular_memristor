
We present a novel modular memristor model that effectively integrates a voltage-controlled memristive core, volatile memory, synaptic-like plasticity and saturation effects into one computationally efficient and modular framework.

This repo contains two minimal working examples for simulating the response for pulse-train stimulation and spike timing-dependent plasticity-like effects. 



################

In the case of pulse-train stimulation, data input is considered in the form of pulses, that have a stimulation pulse of amplitude Va and read pulse of amplitude Vr. 

```
Amplitude
   ^         ↑    ┌──────┐
   |         │    │      │
   |     Va  │    │      │
   |         │    │      │
   |         │    │      │
  0|──────── ↓ ───┘      └────────────────┐        ↑     ┌───────────┘
   |                                      │     Vr │     │ 
   |                                      └─────── ↓ ────┘
   +------------------------------------------------------------------------------------> t
   
                  <- Ta -><---- Tgap ----><----- Tr ----->

                  <------------------------- T ---------------------->   

```
Va      =   stimulation pulse amplitude
Ta      =   stimulation pulse duration

Tgap    =   time between stimulation and read pulses

Vr      =   read pulse amplitude
Tr      =   read pulse duration

T       =   period of stimulation /read cycle.



