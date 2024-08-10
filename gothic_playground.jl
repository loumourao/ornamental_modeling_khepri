using KhepriAutoCAD

include("./gothic_windows.jl")
using .GothicWindows

delete_all_shapes()

#with(current_cs, cs_from_o_vz(u0(), vx())) do
#    arch(xy(-10, -16), xy(10, 16), 1, 0.75, 0.75, 1, 3)
#end

gothic_window(xy(-10, -16), xy(10, 16), 1, 2, 3, 1, 1)

#a = u0()
#b = xy(1, -3)
#c = xy(2, -1)
#d = xy(3, -4.3)
#e = xy(5.2, -6)
#f = xy(5.2, -7)
#g = xy(2.3, -7)
#h = y(-11)
#profile = line(a, b, c, d, e, f, g, h, a)
#
#sweep(circle(u0(), 6), profile)
