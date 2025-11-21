#
# gnuplot command file for plotting flames
#
#
#
#set xrange [0:1.0E-2]
#set yrange [-1.0:1.0]
set data style lines
#
plot "fort.3"
pause 5
#
plot "dens000000.res"
pause 5
#
plot "uvel000000.res"
pause 5
#
plot "ener000000.res"
pause 5
#
plot "pres000000.res"
pause 5
#
plot "temp000000.res"
pause 5
#
plot "yf01000000.res"
pause 5
#
plot "yf02000000.res"
pause 5
#
plot "yf03000000.res"
pause 5
#
plot "yf04000000.res"
pause 5
#
plot "yf05000000.res"
pause 5
#
plot "yf06000000.res"
pause 5
#
plot "yf07000000.res"
pause 5
#
plot "yf08000000.res"
pause 5
#
#plot "yf09000000.res"
#pause 5
##
#plot "yf10000000.res"
#pause 5
##
#plot "yf11000000.res"
#pause 5
##
#plot "yf12000000.res"
#pause 5
##
#plot "yf13000000.res"
#pause 5
##
#plot "yf14000000.res"
#pause 5
##
#plot "yf15000000.res"
#pause 5
##
#plot "yf16000000.res"
#pause 5
##
#plot "yf17000000.res"
#pause 5
##
#plot "yf18000000.res"
#pause 5
#
exit
