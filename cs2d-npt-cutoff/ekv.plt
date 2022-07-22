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

cd sprintf('%s/%.3f/%.3f', path, pres, temp)
set terminal postscript eps enhanced color
set size 0.5,0.5
set xlabel 'i^{th} move'
set ylabel 'U^*+p^*V^*-NT^*ln(V^*)'
set output 'mc_ekv_pot.eps'
plot 'mc_ekv' using 1:2 with lines lc rgb "black"
set ylabel 'H^*'
set output 'mc_ekv_ent.eps'
plot 'mc_ekv' using 1:3 with lines lc rgb "black"
set ylabel 'F^*r^*'
set output 'mc_ekv_vir.eps'
plot 'mc_ekv' using 1:4 with lines lc rgb "black"
set ylabel 'L^*'
set output 'mc_ekv_vol.eps'
plot 'mc_ekv' using 1:5 with lines lc rgb "black"
set xlabel 't^*'
set ylabel 'U^*+p^*V^*-NT^*ln(V^*)'
set output 'md_ekv_pot.eps'
plot 'md_ekv' using 1:2 with lines lc rgb "black"
set ylabel 'H^*'
set output 'md_ekv_ent.eps'
plot 'md_ekv' using 1:3 with lines lc rgb "black"
set ylabel 'E_{kin}^*'
set output 'md_ekv_kin.eps'
plot 'md_ekv' using 1:4 with lines lc rgb "black"
set ylabel 'E_{tot}^*'
set output 'md_ekv_tot.eps'
plot 'md_ekv' using 1:5 with lines lc rgb "black"
set ylabel 'F^*r^*'
set output 'md_ekv_vir.eps'
plot 'md_ekv' using 1:6 with lines lc rgb "black"
set ylabel 'L^*'
set output 'md_ekv_vol.eps'
plot 'md_ekv' using 1:7 with lines lc rgb "black"
set ylabel 'T^*'
set output 'md_ekv_tem.eps'
plot 'md_ekv' using 1:8 with lines lc rgb "black"

cd path