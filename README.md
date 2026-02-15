# modular_memristor

This package contains a memristor model that effectively integrates a voltage-controlled memristive core, volatile memory, synaptic-like plasticity and saturation effects into one computationally efficient and modular framework.

Data input 

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

Va      =   stimulation pulse amplitude
Ta      =   stimulation pulse duration

Tgap    =   time between stimulation and read pulses

Vr      =   read pulse amplitude
Tr      =   read pulse duration

T       =   period of stimulation /read cycle.

this is handled by the package io.py

%%%%%%%%%%%%%%%%%%%%%%%

memmodel.py 


