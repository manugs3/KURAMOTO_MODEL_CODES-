set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'comparacion_R_k_Gaus.png'

set xlabel 'K'
set ylabel 'r'
set key top left
set grid

set logscale x
set format x "10^{%L}"

# Marcador en x = 1.596 con etiqueta k_c
set label "K_c=1.6" at 1.596, graph -0.04 center textcolor rgb "black"
set arrow 5 from 1.596, graph -0.02 to 1.596, graph 0 nohead lw 1 lc rgb "black"

# Estilos
set style line 1 lt 1 lw 0.5 pt 7 ps 1 lc rgb "green"     # n=1000
set style line 2 lt 1 lw 0.5 pt 5 ps 1 lc rgb "red"    # n=5000


# Gráfico
plot \
    	'R-k(n=1000).dat' using 1:2 with points ls 1 title 'N = 1000', \
    	'R-k(n=5000).dat' using 1:2 with points ls 2 title 'N = 5000'
   