set terminal pngcairo size 1200,1000 enhanced font 'Helvetica,12'
set output '4-tiempos.png'

# Paleta de color para fase
set palette defined (-3.2 "black", -1 "blue", 0 "purple", 1 "orange", 3.2 "yellow")
set cbrange [-pi:pi]
set cblabel "Oscillator phase"

# Formato general
unset key
unset xtics
unset ytics
set xlabel "Oscillator index i"
set ylabel "Oscillator index j"
set size ratio 1

# Multiplot: layout de 2x2 con títulos por encima de cada gráfico
set multiplot layout 2,2 rowsfirst title "Oscillator phases for different times" font ",14"

# Subplot 1
set label 1 "t = 0" at screen 0.20,0.98 center font ",12"
plot "4-tiempos.dat" index 0 matrix with image
unset label 1

# Subplot 2
set label 1 "t = 5" at screen 0.70,0.98 center font ",12"
plot "4-tiempos.dat" index 1 matrix with image
unset label 1

# Subplot 3
set label 1 "t = 200" at screen 0.20,0.48 center font ",12"
plot "4-tiempos.dat" index 2 matrix with image
unset label 1

# Subplot 4
set label 1 "t = 5000" at screen 0.70,0.48 center font ",12"
plot "4-tiempos.dat" index 3 matrix with image
unset label 1

unset multiplot