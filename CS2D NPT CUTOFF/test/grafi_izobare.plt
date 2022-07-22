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

n=1
array izoterme[n] = [1.000]

set terminal postscript eps enhanced color
set size 0.5,0.5
set key ins vert left top


set output '!mc_gostota_izoterme.eps'
set xlabel 'p^*'
set ylabel '{/Symbol r}*'
plot for [j=1:n] '!mc_data_izoterma'.j u 2:3 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output '!md_gostota_izoterme.eps'
set xlabel 'p^*'
set ylabel '{/Symbol r}*'
plot for [j=1:n] '!md_data_izoterma'.j u 2:3 title sprintf("%.3f", izoterme[j]) with linespoints ls j

set output '!mc_volumen_izoterme.eps'
set xlabel 'p^*'
set ylabel 'V*/N'
plot for [j=1:n] '!mc_data_izoterma'.j u 2:4 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output '!md_volumen_izoterme.eps'
set xlabel 'p^*'
set ylabel 'V*/N'
plot for [j=1:n] '!md_data_izoterma'.j u 2:4 title sprintf("%.3f", izoterme[j]) with linespoints ls j

set output '!mc_potencialna_izoterme.eps'
set xlabel 'p^*'
set ylabel 'H^*/N'
plot for [j=1:n] '!mc_data_izoterma'.j u 2:5 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output '!md_potencialna_izoterme.eps'
set xlabel 'p^*'
set ylabel 'H^*/N'
plot for [j=1:n] '!md_data_izoterma'.j u 2:5 title sprintf("%.3f", izoterme[j]) with linespoints ls j

set output '!mc_kapaciteta_izoterme.eps'
set xlabel 'p^*'
set ylabel 'C_P^*/N'
plot for [j=1:n] '!mc_data_izoterma'.j u 2:6 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output '!md_kapaciteta_izoterme.eps'
set xlabel 'p^*'
set ylabel 'C_P^*/N'
plot for [j=1:n] '!md_data_izoterma'.j u 2:6 title sprintf("%.3f", izoterme[j]) with linespoints ls j

set output '!mc_kappa_izoterme.eps'
set xlabel 'p^*'
set ylabel '{/Symbol k}^*'
plot for [j=1:n] '!mc_data_izoterma'.j u 2:7 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output '!md_kappa_izoterme.eps'
set xlabel 'p^*'
set ylabel '{/Symbol k}^*'
plot for [j=1:n] '!md_data_izoterma'.j u 2:7 title sprintf("%.3f", izoterme[j]) with linespoints ls j

set output '!mc_alpha_izoterme.eps'
set xlabel 'p^*'
set ylabel '{/Symbol a}^*'
plot for [j=1:n] '!mc_data_izoterma'.j u 2:8 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output '!md_alpha_izoterme.eps'
set xlabel 'p^*'
set ylabel '{/Symbol a}^*'
plot for [j=1:n] '!md_data_izoterma'.j u 2:8 title sprintf("%.3f", izoterme[j]) with linespoints ls j

set output '!mc_kemijski_izoterme.eps'
set xlabel 'p^*'
set ylabel '{/Symbol m}^*'
plot for [j=1:n] '!mc_data_izoterma'.j u 2:9 title sprintf("%.3f", izoterme[j]) with linespoints ls j
set output '!md_kemijski_izoterme.eps'
set xlabel 'p^*'
set ylabel '{/Symbol m}^*'
plot for [j=1:n] '!md_data_izoterma'.j u 2:9 title sprintf("%.3f", izoterme[j]) with linespoints ls j


set output 'mc_gostota_izoterme.eps'
set xlabel 'p^*'
set ylabel '{/Symbol r}*'
plot for [j=1:n] 'mc_data_izoterma'.j u 2:3:4 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j
set output 'md_gostota_izoterme.eps'
set xlabel 'p^*'
set ylabel '{/Symbol r}*'
plot for [j=1:n] 'md_data_izoterma'.j u 2:3:4 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j

set output 'mc_volumen_izoterme.eps'
set xlabel 'p^*'
set ylabel 'V*/N'
plot for [j=1:n] 'mc_data_izoterma'.j u 2:5:6 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j
set output 'md_volumen_izoterme.eps'
set xlabel 'p^*'
set ylabel 'V*/N'
plot for [j=1:n] 'md_data_izoterma'.j u 2:5:6 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j

set output 'mc_potencialna_izoterme.eps'
set xlabel 'p^*'
set ylabel 'H^*/N'
plot for [j=1:n] 'mc_data_izoterma'.j u 2:7:8 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j
set output 'md_potencialna_izoterme.eps'
set xlabel 'p^*'
set ylabel 'H^*/N'
plot for [j=1:n] 'md_data_izoterma'.j u 2:7:8 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j

set output 'mc_kapaciteta_izoterme.eps'
set xlabel 'p^*'
set ylabel 'C_P^*/N'
plot for [j=1:n] 'mc_data_izoterma'.j u 2:9:10 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j
set output 'md_kapaciteta_izoterme.eps'
set xlabel 'p^*'
set ylabel 'C_P^*/N'
plot for [j=1:n] 'md_data_izoterma'.j u 2:9:10 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j

set output 'mc_kappa_izoterme.eps'
set xlabel 'p^*'
set ylabel '{/Symbol k}^*'
plot for [j=1:n] 'mc_data_izoterma'.j u 2:11:12 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j
set output 'md_kappa_izoterme.eps'
set xlabel 'p^*'
set ylabel '{/Symbol k}^*'
plot for [j=1:n] 'md_data_izoterma'.j u 2:11:12 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j

set output 'mc_alpha_izoterme.eps'
set xlabel 'p^*'
set ylabel '{/Symbol a}^*'
plot for [j=1:n] 'mc_data_izoterma'.j u 2:13:14 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j
set output 'md_alpha_izoterme.eps'
set xlabel 'p^*'
set ylabel '{/Symbol a}^*'
plot for [j=1:n] 'md_data_izoterma'.j u 2:13:14 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j

set output 'mc_kemijski_izoterme.eps'
set xlabel 'p^*'
set ylabel '{/Symbol m}^*'
plot for [j=1:n] 'mc_data_izoterma'.j u 2:15:16 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j
set output 'md_kemijski_izoterme.eps'
set xlabel 'p^*'
set ylabel '{/Symbol m}^*'
plot for [j=1:n] 'md_data_izoterma'.j u 2:15:16 title sprintf("%.3f", izoterme[j]) with yerrorbars ls j