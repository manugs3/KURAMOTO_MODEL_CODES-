set term gif animate optimize delay 10 size 600,600
set output 'osciladores_red.gif'

# Rango de colores y paleta
set cbrange [-pi:pi]
set palette defined (-pi "black", -1 "blue", 0 "purple", 1 "orange", pi "yellow")

# Axes
set size ratio 1
unset key
unset xtics
unset ytics

do for[ii=0:2000]{plot "thetas-t-2D.dat" i ii matrix with image notitle}
