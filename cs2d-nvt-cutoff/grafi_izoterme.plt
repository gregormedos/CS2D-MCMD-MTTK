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

n=3
array izoterme[n] = [0.120, 0.500, 1.000]

set terminal postscript eps enhanced color
set size 0.5,0.5
set key ins vert left top

set xlabel '{/Symbol r}^*'

set ylabel 'U^*/N'
set output '!mc_energija_izoterme.eps'
plot for [j=1:n] '!mc_data_izoterma'.j u 2:3 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output '!md_energija_izoterme.eps'
plot for [j=1:n] '!md_data_izoterma'.j u 2:3 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output 'mc_energija_izoterme.eps'
plot for [j=1:n] 'mc_data_izoterma'.j u 2:3:4 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j
set output 'md_energija_izoterme.eps'
plot for [j=1:n] 'md_data_izoterma'.j u 2:3:4 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j

set ylabel 'p^*'
set output '!mc_tlak_izoterme.eps'
plot for [j=1:n] '!mc_data_izoterma'.j u 2:4 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output '!md_tlak_izoterme.eps'
plot for [j=1:n] '!md_data_izoterma'.j u 2:4 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output 'mc_tlak_izoterme.eps'
plot for [j=1:n] 'mc_data_izoterma'.j u 2:5:6 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j
set output 'md_tlak_izoterme.eps'
plot for [j=1:n] 'md_data_izoterma'.j u 2:5:6 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j

set ylabel 'C_V^*/N'
set output '!mc_kapaciteta_izoterme.eps'
plot for [j=1:n] '!mc_data_izoterma'.j u 2:5 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output '!md_kapaciteta_izoterme.eps'
plot for [j=1:n] '!md_data_izoterma'.j u 2:5 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output 'mc_kapaciteta_izoterme.eps'
plot for [j=1:n] 'mc_data_izoterma'.j u 2:7:8 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j
set output 'md_kapaciteta_izoterme.eps'
plot for [j=1:n] 'md_data_izoterma'.j u 2:7:8 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j

set ylabel '{/Symbol t}^*'
set output '!mc_translacijski-parameter-urejenosti_izoterme.eps'
plot for [j=1:n] '!mc_data_izoterma'.j u 2:6 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output '!md_translacijski-parameter-urejenosti_izoterme.eps'
plot for [j=1:n] '!md_data_izoterma'.j u 2:6 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output 'mc_translacijski-parameter-urejenosti_izoterme.eps'
plot for [j=1:n] 'mc_data_izoterma'.j u 2:9:10 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j
set output 'md_translacijski-parameter-urejenosti_izoterme.eps'
plot for [j=1:n] 'md_data_izoterma'.j u 2:9:10 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j

set ylabel 'D^*'
set output '!mc_difuzijski-koeficient_izoterme.eps'
plot for [j=1:n] '!mc_data_izoterma'.j u 2:7 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output '!md_difuzijski-koeficient_izoterme.eps'
plot for [j=1:n] '!md_data_izoterma'.j u 2:7 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output 'mc_difuzijski-koeficient_izoterme.eps'
plot for [j=1:n] 'mc_data_izoterma'.j u 2:11:12 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j
set output 'md_difuzijski-koeficient_izoterme.eps'
plot for [j=1:n] 'md_data_izoterma'.j u 2:11:12 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j

set ylabel 'D^*'
set output '!mc_difuzijski-koeficient-nicla_izoterme.eps'
plot for [j=1:n] '!mc_data_izoterma'.j u 2:8 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output '!md_difuzijski-koeficient-nicla_izoterme.eps'
plot for [j=1:n] '!md_data_izoterma'.j u 2:8 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output 'mc_difuzijski-koeficient-nicla_izoterme.eps'
plot for [j=1:n] 'mc_data_izoterma'.j u 2:13:14 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j
set output 'md_difuzijski-koeficient-nicla_izoterme.eps'
plot for [j=1:n] 'md_data_izoterma'.j u 2:13:14 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j
