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
do for [i=1:nseries] {
  set ylabel 'U^*'
  set output sprintf('mc_sam%03d_pot.eps', i)
  plot sprintf('mc_sam%03d', i) using 1:2 with lines lc rgb "black"
  set ylabel 'F^*r^*'
  set output sprintf('mc_sam%03d_vir.eps', i)
  plot sprintf('mc_sam%03d', i) using 1:3 with lines lc rgb "black"
}
set xlabel 't^*'
do for [i=1:nseries] {
  set ylabel 'U^*'
  set output sprintf('md_sam%03d_pot.eps', i)
  plot sprintf('md_sam%03d', i) using 1:2 with lines lc rgb "black"
  set ylabel 'E_{kin}^*'
  set output sprintf('md_sam%03d_kin.eps', i)
  plot sprintf('md_sam%03d', i) using 1:3 with lines lc rgb "black"
  set ylabel 'E_{tot}^*'
  set output sprintf('md_sam%03d_tot.eps', i)
  plot sprintf('md_sam%03d', i) using 1:4 with lines lc rgb "black"
  set ylabel 'F^*r^*'
  set output sprintf('md_sam%03d_vir.eps', i)
  plot sprintf('md_sam%03d', i) using 1:5 with lines lc rgb "black"
  set ylabel 'T^*'
  set output sprintf('md_sam%03d_tem.eps', i)
  plot sprintf('md_sam%03d', i) using 1:6 with lines lc rgb "black"
}

cd path