set term qt

set xlabel "Timestep"
set ylabel "mflow"

plot [*:*][0:6e-10] "log.sparta" u 1:12, 4e-10 lw 5 t "Target: 4e-10"


