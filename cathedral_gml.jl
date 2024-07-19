using KhepriAutoCAD

delete_all_shapes()

# Pillars are stored in an array and oriented N-S (row-wise) and W-E (column-wise)
# dir is important because they are usually not symmetric as they will be attached to different typse of walls
# We will use the following convention in the xy-plane (which can later be switched to other planes through appropriate trasnformations)
#      N (y^)
#      |
# W -- O -- E (x>)
#      |   
#      S
# where 'O' is the origin
struct pillar
    type
    coord
    pos
    top
    bot
    dir
    wallN
    wallS
    wallE
    wallW
end

# Start pillars
pillars = Array{pillar}(nothing, 10, 17)

# Instantiate green section 
# TODO: Instantiate row 1-2, columns 8-11
# TODO: Instantiate row 9-10, columns 8-11

# Instantiate orange section
# TODO: Instantiate row 3-8, columns 1-3

# Instantiate blue section (main nave)
# TODO: Instantiate row 3-8, columns 4-13

# Instantiate yellow section
# TODO: Instantiate row 3-8, columns 15-17