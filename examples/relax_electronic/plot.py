#!/bin/env python3

import numpy
import math

def eq_E( T ):
    boltz = 1.38e-23
    E_state = [0.0, 72231.4, 85782.7, 98057.1133]
    degen = [9, 27, 54, 18]
    E = 0
    part_denom = 0
    for i in range(len(E_state)):
         E += E_state[i]*degen[i]*math.exp(-E_state[i]/T)
         part_denom += degen[i]*math.exp(-E_state[i]/T)
    E = boltz*E/part_denom
    return E


with open('log.sparta','r') as f:
    header = True
    footer = False
    data = {}
    for row in f:
        if header:
            if 'Step' in row:
                header = False
                col_headers = row.split()
                for col in col_headers:
                    data[col] = []
            continue

        if not footer:
            if 'Loop' in row:
                footer = True
                continue
        else:
            continue

        values = row.split()
        for i in range(len(col_headers)):
            data[col_headers[i]].append(float(values[i]))
            
import matplotlib.pyplot


calced_eq_E = []
for T in data['c_Ttrans']:
   calced_eq_E.append(eq_E(T)) 

fig, ax = matplotlib.pyplot.subplots()
ax.plot(data['Step'], calced_eq_E, label='Elec(Eq)')
ax.plot(data['Step'], data['c_Eelec'], label='Elec')
ax.set_xlabel( 'Timestep (ns)')
ax.set_ylabel( 'T (K)')
ax.legend(loc='upper right')
fig.savefig('temps.png',dpi=800,bbox_inches='tight')
