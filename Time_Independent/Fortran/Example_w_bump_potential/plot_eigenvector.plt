set xrange [-10:10]
set tics font "Helvetica,15"
set key font ",15"
set title "Eigenvectors of simple SHO" font "Helvetica,16"

plot   'Eigenvectors.dat' u 1:2 w l lw 5 title 'state-1'
replot 'Eigenvectors.dat' u 1:3 w l lw 5 title 'state-2' 
replot 'Eigenvectors.dat' u 1:4 w l lw 5 dt 2 title 'state-3'
replot 'Eigenvectors.dat' u 1:5 w l lw 5 dt 4 title 'state-4'
replot 'Eigenvectors.dat' u 1:6 w l lw 5 dt 6 title 'state-5'
