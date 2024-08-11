using KhepriAutoCAD

delete_all_shapes()

# Pillars are stored in an array and oriented N-S (row-wise) and W-E (column-wise)
# We will use the following convention in the xy-plane (which can later be switched to other planes through appropriate trasnformations)
#      N (y^)
#      |
# W -- O -- E (x>)
#      |   
#      S
# where 'O' is the origin

const EQUIDISTANT_SECTIONS_LENGTH = 7.53
const ORANGE_SECTION_JUMP_LENGTH = EQUIDISTANT_SECTIONS_LENGTH + EQUIDISTANT_SECTIONS_LENGTH / 2

const START_INDEX = 1
const END_INDEX = -1

const ORANGE_SECTION_ROW_RANGE = range(3, 8)
const ORANGE_SECTION_COLUMN_RANGE = range(1, 3)

const BLUE_SECTION_ROW_RANGE = range(3, 8)
const BLUE_SECTION_COLUMN_RANGE = range(4, 13)

const FIRST_GREEN_SECTION_ROW_RANGE = range(1, 2)
const SECOND_GREEN_SECTION_ROW_RANGE = range(9, 10)
const GREEN_SECTION_COLUMN_RANGE = range(8, 11)

const YELLOW_SECTION_ROW_RANGE = range(3, 8)
const YELLOW_SECTION_COLUMN_RANGE = range(14, 17)

# the remaining attributes are filled once the pillar geometry is created
struct pillar
    type #mandatory
    coord #mandatory
    pos #mandatory
    top
    bot
    dir
    wallN
    wallS
    wallE
    wallW
end

struct wall
    type
    colR
    colL
end

function buttress_pillar(center, height)
    half_height = height / 2
    bottom_left_corner = center - vxy(1, 1)
    upper_right_corner = center + vxy(1, 1)

    pillar_base_diagonal_vector = upper_right_corner - bottom_left_corner
    pillar_base_width = abs(pillar_base_diagonal_vector.x)
    pillar_base_height = abs(pillar_base_diagonal_vector.y)

    pillar_shape = surface_rectangle(bottom_left_corner, upper_right_corner)

    pillar_half_position = center + vz(half_height)
    pillar_top_position = center + vz(height)

    pillar_bottom_path = line(center, pillar_half_position)
    pillar_top_path = line(pillar_half_position, pillar_top_position)

    pillar_bottom = sweep(pillar_bottom_path, pillar_shape)
    pillar_top = sweep(pillar_top_path, pillar_shape)
    
    buttress_pillar_shape_bottom_left_corner = center - vx(pillar_base_width * 1.5)
    buttress_pillar_shape_upper_right_corner = center + vxy(pillar_base_width * 1.5, pillar_base_height * 5)

    buttress_pillar_mid_shape_bottom_left_corner = buttress_pillar_shape_bottom_left_corner - vx(pillar_base_width / 2) + vz(half_height)
    buttress_pillar_mid_shape_upper_right_corner = buttress_pillar_shape_upper_right_corner + vxyz(pillar_base_width / 2, pillar_base_height, half_height)

    buttress_pillar_base = surface_rectangle(buttress_pillar_shape_bottom_left_corner, buttress_pillar_shape_upper_right_corner)
    buttress_pillar_pre_mid = surface_rectangle(buttress_pillar_shape_bottom_left_corner + vz(half_height - pillar_base_width / 8), 
                                                    buttress_pillar_shape_upper_right_corner + vz(half_height - pillar_base_width / 8))
    buttress_pillar_mid = surface_rectangle(buttress_pillar_mid_shape_bottom_left_corner, buttress_pillar_mid_shape_upper_right_corner)
    buttress_pillar_post_mid = surface_rectangle(buttress_pillar_shape_bottom_left_corner + vz(half_height + pillar_base_width / 2), 
                                                    buttress_pillar_shape_upper_right_corner + vz(half_height + pillar_base_width / 2))
    buttress_pillar_top = surface_rectangle(buttress_pillar_shape_bottom_left_corner + vz(height), buttress_pillar_shape_upper_right_corner + vz(height))

    buttress_pillar = loft_ruled([buttress_pillar_base, buttress_pillar_pre_mid, 
                                    buttress_pillar_mid, buttress_pillar_post_mid, buttress_pillar_top])

    return union(pillar_bottom, pillar_top, buttress_pillar)
end

#pillars = Array{pillar}(nothing, 10, 17)
#
#function instantiate_orange_section(pillars)
#end
#
#function instantiate_blue_section(pillars)
#    for row = BLUE_SECTION_ROW_RANGE
#        for column = BLUE_SECTION_COLUMN_RANGE
#        end
#    end
#    return pillars
#end
#
#function instantiate_green_section(pillars)
#end
#
#function instantiate_yellow_section(pillars)
#end
#
#function instantiate_pillars(pillars)
#    pillars = instantiate_orange_section(pillars)
#    pillars = instantiate_blue_section(pillars)
#    pillars = instantiate_green_section(pillars)
#    pillars = instantiate_yellow_section(pillars)
#
#    return pillars
#end

buttress_pillar(u0(), 100)
# Start pillars

# Instantiate green section
# TODO: Instantiate row 1-2, columns 8-11
# TODO: Instantiate row 9-10, columns 8-11

# Instantiate orange section
# TODO: Instantiate row 3-8, columns 1-3

# Instantiate blue section (main nave)
# TODO: Instantiate row 3-8, columns 4-13

# Instantiate yellow section
# TODO: Instantiate row 3-8, columns 15-17
