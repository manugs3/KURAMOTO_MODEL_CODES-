reset
set terminal pngcairo size 1600,800 enhanced font 'Arial,12'
set output 'phase_frequency_evolution.png'

set multiplot layout 1,2 title "N=5000, k=4, r= 0.946"  font ",14"

# --- FASES ---
set title ''
set xlabel 't[t.u.]'
set ylabel '{/Symbol q}'
set key outside

plot for [i=2:51] 'thetas-t.dat' using 1:i with lines notitle 

# --- FRECUENCIAS ---
set title ''
set xlabel 't[t.u.]'
set ylabel 'd{/Symbol q}/dt'
set key outside

plot for [i=2:51] 'omegas-t.dat' using 1:i with lines notitle

unset multiplot