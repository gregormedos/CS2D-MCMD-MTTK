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

cd sprintf('%s/%.3f/%.3f', path, temp, dens)
set terminal postscript eps enhanced color
set size 0.5,0.5
set xlabel 'i^{th} move'
set ylabel 'U^*'
set output 'mc_min_pot.eps'
plot 'mc_min' using 1:2 with lines lc rgb "black"
set ylabel 'F^*r^*'
set output 'mc_min_vir.eps'
plot 'mc_min' using 1:3 with lines lc rgb "black"

cd path