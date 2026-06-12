set terminal pdfcairo enhanced color size 10cm,25cm font "Arial,12"
set output "nilsson_diagram_tracked.pdf"

unset key
set xlabel "{/Symbol d}"
set ylabel "Energy"
set xrange [-0.30:0.30]
set grid

plot "tracks.dat" using 1:2 with lines lw 1