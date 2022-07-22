# Gnuplot script file for plotting data
set encoding utf8
set autoscale                                                       # scale axes automatically
unset logscale                                                      # remove any log-scaling
unset xlabel                                                        # remove any previous xlabel
unset ylabel                                                        # remove any previous ylabel
set xtics auto                                                      # set xtics automatically
set ytics auto                                                      # set ytics automatically
unset key
unset size

load 'palette.pal'

cd sprintf('%s/%.3f/%.3f', path, temp, dens)
set terminal postscript eps enhanced color
set size 0.5,0.5
set key ins top left

set xlabel 'n'
set ylabel 'MSD(n)'
set output 'mc_msd.eps'
plot 'mc_msd' using 1:($2+$3):($2-$3) notitle with filledcurve fc rgb "gray", \
'mc_msd' title "MSD(n)" with lines lw 2 lc rgb "black", \
'mc_msdlin' title "lin" w lines lw 2 dt 2 lc rgb "black", \
'mc_msdlinzero' title "linzero" w lines lw 2 dt 3 lc rgb "black"

set xlabel 't^*'
set ylabel 'MSD(t^*)'
set output 'md_msd.eps'
plot 'md_msd' using 1:($2+$3):($2-$3) notitle with filledcurve fc rgb "gray", \
'md_msd' title "MSD(t^*)" with lines lw 2 lc rgb "black", \
'md_msdlin' title "lin" w lines lw 2 dt 2 lc rgb "black", \
'md_msdlinzero' title "linzero" w lines lw 2 dt 3 lc rgb "black"

unset key

set xlabel 'n'
set ylabel 'MSD(n)'
set output 'mc_msd_sampl.eps'
plot for [i=1:20] sprintf('mc_msd%03d', i) with lines ls i
set output 'mc_msdlin_sampl.eps'
plot for [i=1:20] sprintf('mc_msdlin%03d', i) with lines ls i
set output 'mc_msdlinzero_sampl.eps'
plot for [i=1:20] sprintf('mc_msdlinzero%03d', i) with lines ls i

set xlabel 't^*'
set ylabel 'MSD(t^*)'
set output 'md_msd_sampl.eps'
plot for [i=1:20] sprintf('md_msd%03d', i) with lines ls i
set output 'md_msdlin_sampl.eps'
plot for [i=1:20] sprintf('md_msdlin%03d', i) with lines ls i
set output 'md_msdlinzero_sampl.eps'
plot for [i=1:20] sprintf('md_msdlinzero%03d', i) with lines ls i

cd path