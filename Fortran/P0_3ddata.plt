set xlabel "x"
set ylabel "y"
set dgrid3d 29,19
set hidden3d
set terminal x11 0
set nokey
set grid
set mouse
set title "My 3D"
splot "3ddata.dat" u 1:2:3 with lines;pause mouse
keypress
