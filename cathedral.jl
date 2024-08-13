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

# == INDEXES == #
ROW_INDEX = 1
COLUMN_INDEX = 2
ORIENTATION_INDEX = 3
# == INDEXES == #

# == ORIENTATIONS == #
NORTH = vy(1)
SOUTH = vy(-1)
WEST = vx(-1)
EAST = vx(1)
# == ORIENTATIONS == #

# == BASIC PILLAR INFO == #
BASIC_PILLAR_WIDTH = 1
BASIC_PILLAR_DEPTH = 1
BASIC_PILLAR_HEIGHT = 105
# == BASIC PILLAR INFO == #

# == MEASUREMENTS == #
EQUIDISTANT_SECTIONS_LENGTH = 7.53
WEST_END_JUMP_LENGTH = EQUIDISTANT_SECTIONS_LENGTH + EQUIDISTANT_SECTIONS_LENGTH / 2
# == MEASUREMENTS == #

# == COMPLEX PILLARS INFO == #
UPPER_RIGHT_OUTER_CROSSING_BUTTRESS_INFO = (3, 12, EAST)
UPPER_LEFT_OUTER_CROSSING_BUTTRESS_INFO = (3, 8, NORTH)
LOWER_LEFT_OUTER_CROSSING_BUTTRESS_INFO = (9, 8, WEST)
LOWER_RIGHT_OUTER_CROSSING_BUTTRESS_INFO = (9, 12, SOUTH)
OUTER_CROSSING_BUTTRESSES_INFO = [UPPER_RIGHT_OUTER_CROSSING_BUTTRESS_INFO, UPPER_LEFT_OUTER_CROSSING_BUTTRESS_INFO,
                                    LOWER_LEFT_OUTER_CROSSING_BUTTRESS_INFO, LOWER_RIGHT_OUTER_CROSSING_BUTTRESS_INFO]

UPPER_RIGHT_CROSSING_PILLAR_INFO = (5, 11, EAST)
UPPER_LEFT_CROSSING_PILLAR_INFO = (5, 9, WEST)
LOWER_LEFT_CROSSING_PILLAR_INFO = (7, 9, WEST)
LOWER_RIGHT_CROSSING_PILLAR_INFO = (7, 11, EAST)
CROSSING_PILLARS_INFO = [UPPER_RIGHT_CROSSING_PILLAR_INFO, UPPER_LEFT_CROSSING_PILLAR_INFO, 
                            LOWER_LEFT_CROSSING_PILLAR_INFO, LOWER_RIGHT_CROSSING_PILLAR_INFO]
# == COMPLEX PILLARS INFO == #

struct Pillar
    type
    row
    column
    center
    orientation
    bottom_anchor
    top_anchor
    north_wall
    south_wall
    west_wall
    east_wall
end

function Pillar(type, row, column, center, orientation;
                    bottom_anchor = nothing, top_anchor = nothing, 
                        north_wall = nothing, south_wall = nothing, west_wall = nothing, east_wall = nothing)
    return Pillar(type, row, column, center, orientation, 
                    bottom_anchor, top_anchor, 
                        north_wall, south_wall, west_wall, east_wall)
end

struct Wall
    type
    colR
    colL
end

# == PILLARS == #
function get_orientation_polar_angle(orientation)
    return atan(orientation.y, orientation.x)
end

# == MODELERS == #
function basic_pillar(center; 
                        height = BASIC_PILLAR_HEIGHT)
    half_width = BASIC_PILLAR_WIDTH / 2
    half_depth = BASIC_PILLAR_DEPTH / 2

    bottom_left_corner = u0() - vxy(half_width, half_depth)
    upper_right_corner = u0() + vxy(half_width, half_depth)

    pillar_shape = surface_rectangle(bottom_left_corner, upper_right_corner)
    pillar_path = line(center, center + vz(height))

    return sweep(pillar_path, pillar_shape)
end

function crossing_pillar(center; 
                            height = BASIC_PILLAR_HEIGHT * 1.15)
    three_quarters_width = BASIC_PILLAR_WIDTH * 0.75
    three_quarters_depth = BASIC_PILLAR_DEPTH * 0.75

    bottom_left_corner = u0() - vxy(three_quarters_width, three_quarters_depth)
    upper_right_corner = u0() + vxy(three_quarters_width, three_quarters_depth)

    pillar_shape = surface_rectangle(bottom_left_corner, upper_right_corner)
    pillar_path = line(center, center + vz(height))

    return sweep(pillar_path, pillar_shape)
end

function half_outer_crossing_buttress(center, height)
    half_height = height / 2

    half_width = BASIC_PILLAR_WIDTH / 2
    eighth_width = BASIC_PILLAR_WIDTH / 8
    quadruple_width = BASIC_PILLAR_WIDTH * 4

    half_depth = BASIC_PILLAR_DEPTH / 2

    base_bottom_left_corner = center - vxy(half_width, half_depth)
    base_upper_right_corner = base_bottom_left_corner + vxy(quadruple_width, BASIC_PILLAR_DEPTH)

    pre_mid_bottom_left_corner = base_bottom_left_corner + vz(half_height - eighth_width)
    pre_mid_upper_right_corner = base_upper_right_corner + vz(half_height - eighth_width)

    mid_bottom_left_corner = base_bottom_left_corner + vz(half_height)
    mid_upper_right_corner = base_upper_right_corner + vxyz(half_width, BASIC_PILLAR_DEPTH, half_height)

    post_mid_bottom_left_corner = base_bottom_left_corner + vz(half_height + half_width)
    post_mid_upper_right_corner = base_upper_right_corner + vz(half_height + half_width)

    top_bottom_left_corner = base_bottom_left_corner + vz(height)
    top_upper_right_corner = base_upper_right_corner + vz(height)

    pillar_base = surface_rectangle(base_bottom_left_corner, base_upper_right_corner)
    pillar_pre_mid = surface_rectangle(pre_mid_bottom_left_corner, pre_mid_upper_right_corner)
    pillar_mid = surface_rectangle(mid_bottom_left_corner, mid_upper_right_corner)
    pillar_post_mid = surface_rectangle(post_mid_bottom_left_corner, post_mid_upper_right_corner)
    pillar_top = surface_rectangle(top_bottom_left_corner, top_upper_right_corner)

    return loft_ruled([pillar_base, pillar_pre_mid, pillar_mid, pillar_post_mid, pillar_top])
end

function outer_crossing_buttress(center, orientation; 
                                    height = BASIC_PILLAR_HEIGHT * 0.95)
    orientation = get_orientation_polar_angle(orientation)

    quarter_width = BASIC_PILLAR_WIDTH / 4
    three_quarters_width = quarter_width * 3

    quarter_depth = BASIC_PILLAR_DEPTH / 4
    three_quarters_depth = quarter_depth * 3

    pillar = basic_pillar(center; height = height)

    transformation_point = center + vxy(quarter_width, quarter_depth)
    right_buttress = move(half_outer_crossing_buttress(center, height), vxy(three_quarters_width, three_quarters_depth))
    left_buttress = mirror(rotate(deepcopy(right_buttress), π/2, transformation_point), transformation_point + vz())

    outer_crossing_buttress = union(pillar, right_buttress, left_buttress)

    return rotate(outer_crossing_buttress, orientation, center)
end

function aisle_inner_pillar(center, orientation; 
                                height = BASIC_PILLAR_HEIGHT * 1.10)
    orientation = get_orientation_polar_angle(orientation)

    half_width = BASIC_PILLAR_WIDTH / 2
    three_quarters_depth = BASIC_PILLAR_DEPTH * 0.75

    bottom_left_corner = u0() - vxy(half_width, three_quarters_depth)
    upper_right_corner = u0() + vxy(half_width, three_quarters_depth)

    pillar_shape = surface_rectangle(bottom_left_corner, upper_right_corner)
    pillar_path = line(center, center + vz(height))

    aisle_inner_pillar = sweep(pillar_path, pillar_shape)
    return rotate(aisle_inner_pillar, orientation, center)
end

function aisle_buttress(center, orientation; 
                            height = BASIC_PILLAR_HEIGHT * 0.95)
    orientation = get_orientation_polar_angle(orientation)
    half_height = height / 2

    eighth_width = BASIC_PILLAR_WIDTH / 8
    half_width = BASIC_PILLAR_WIDTH / 2
    quintuple_depth = BASIC_PILLAR_DEPTH * 5

    pillar = basic_pillar(center; height = height)

    base_bottom_left_corner = center - vx(BASIC_PILLAR_WIDTH)
    base_upper_right_corner = center + vxy(BASIC_PILLAR_WIDTH, quintuple_depth)

    pre_mid_bottom_left_corner = base_bottom_left_corner + vz(half_height - eighth_width)
    pre_mid_upper_right_corner = base_upper_right_corner + vz(half_height - eighth_width)
    
    mid_bottom_left_corner = base_bottom_left_corner - vx(half_width) + vz(half_height)
    mid_upper_right_corner = base_upper_right_corner + vxyz(half_width, BASIC_PILLAR_DEPTH, half_height)

    post_mid_bottom_left_corner = base_bottom_left_corner + vz(half_height + half_width)
    post_mid_upper_right_corner = base_upper_right_corner + vz(half_height + half_width)

    top_bottom_left_corner = base_bottom_left_corner + vz(height)
    top_upper_right_corner = base_upper_right_corner + vz(height)

    pillar_base = surface_rectangle(base_bottom_left_corner, base_upper_right_corner)
    pillar_pre_mid = surface_rectangle(pre_mid_bottom_left_corner, pre_mid_upper_right_corner)
    pillar_mid = surface_rectangle(mid_bottom_left_corner, mid_upper_right_corner)
    pillar_post_mid = surface_rectangle(post_mid_bottom_left_corner, post_mid_upper_right_corner)
    pillar_top = surface_rectangle(top_bottom_left_corner, top_upper_right_corner)

    aisle_buttress = rotate(union(pillar, loft_ruled([pillar_base, pillar_pre_mid, pillar_mid, pillar_post_mid, pillar_top])), -π/2, center)

    return rotate(aisle_buttress, orientation, center)
end 

function ambulatory_inner_pillar(center, orientation; 
                                    height = BASIC_PILLAR_HEIGHT * 1.10)
    orientation = get_orientation_polar_angle(orientation)

    third_width = BASIC_PILLAR_WIDTH / 3
    eighth_width = BASIC_PILLAR_WIDTH / 8

    half_depth = BASIC_PILLAR_DEPTH / 2

    upper_left_corner = u0() - vx(third_width) + vy(half_depth)
    upper_right_corner = upper_left_corner + vx(third_width * 2)
    bottom_left_corner = u0() - vxy(eighth_width, half_depth)
    bottom_right_corner = bottom_left_corner + vx(eighth_width * 2)

    pillar_shape = surface_polygon(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner)
    pillar_path = line(center, center + vz(height))
    
    ambulatory_inner_pillar = rotate(sweep(pillar_path, pillar_shape), -π, center)

    return rotate(ambulatory_inner_pillar, orientation, center)
end

function ambulatory_outer_pillar(center, orientation; 
                                    height = BASIC_PILLAR_HEIGHT)
    orientation = get_orientation_polar_angle(orientation)

    half_width = BASIC_PILLAR_WIDTH / 2
    third_width = BASIC_PILLAR_WIDTH / 3

    half_depth = BASIC_PILLAR_DEPTH / 2

    upper_left_corner = u0() - vx(BASIC_PILLAR_WIDTH - third_width) + vy(half_depth)
    upper_right_corner = upper_left_corner + vx((BASIC_PILLAR_WIDTH - third_width) * 2)
    bottom_left_corner = u0() - vxy(half_width, half_depth)
    bottom_right_corner = bottom_left_corner + vx(BASIC_PILLAR_WIDTH)

    pillar_shape = surface_polygon(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner)
    pillar_path = line(center, center + vz(height))
    
    ambulatory_inner_pillar = rotate(sweep(pillar_path, pillar_shape), -π, center)

    return rotate(ambulatory_inner_pillar, orientation, center)
end

function ambulatory_buttress(center, orientation; 
                                height = BASIC_PILLAR_HEIGHT * 0.64)
    orientation = get_orientation_polar_angle(orientation)

    half_width = BASIC_PILLAR_WIDTH / 2

    third_depth = BASIC_PILLAR_DEPTH / 3

    rectangle_bottom_left_corner = u0() - vx(BASIC_PILLAR_WIDTH) + vy(third_depth)
    rectangle_upper_right_corner = u0() + vxy(BASIC_PILLAR_WIDTH, BASIC_PILLAR_DEPTH + third_depth)
    rectangle_bottom_right_corner = xy(rectangle_upper_right_corner.x, rectangle_bottom_left_corner.y)
    rectangle_upper_left_corner = xy(rectangle_bottom_left_corner.x, rectangle_upper_right_corner.y)
    trapezoid_bottom_left_corner = rectangle_bottom_left_corner - vxy(half_width, third_depth)
    trapezoid_bottom_right_corner = rectangle_bottom_right_corner + vx(half_width) - vy(third_depth)

    pillar_shape = surface_polygon(trapezoid_bottom_left_corner, trapezoid_bottom_right_corner, 
                                        rectangle_bottom_right_corner, rectangle_upper_right_corner, 
                                            rectangle_upper_left_corner, rectangle_bottom_left_corner)
    pillar_path = line(center, center + vz(height))

    ambulatory_buttress = rotate(sweep(pillar_path, pillar_shape), -π, center)
    
    return rotate(ambulatory_buttress, orientation, center)
end
# == MODELERS == #

# == INSTANTIATORS == #
# Pillar instantiation for display purposes
#crossing_pillar(u0() - vxy(EQUIDISTANT_SECTIONS_LENGTH, EQUIDISTANT_SECTIONS_LENGTH * 2))
#outer_crossing_buttress(u0(), EAST)
#aisle_buttress(u0() + vx(EQUIDISTANT_SECTIONS_LENGTH * 2), NORTH)
#basic_pillar(u0() + vx(EQUIDISTANT_SECTIONS_LENGTH * 2) - vy(EQUIDISTANT_SECTIONS_LENGTH))
#aisle_inner_pillar(u0() + vx(EQUIDISTANT_SECTIONS_LENGTH * 2) - vy(EQUIDISTANT_SECTIONS_LENGTH * 2), NORTH)
#ambulatory_buttress(u0() + vx(EQUIDISTANT_SECTIONS_LENGTH * 4), NORTH)
#ambulatory_outer_pillar(u0() + vx(EQUIDISTANT_SECTIONS_LENGTH * 4) - vy(EQUIDISTANT_SECTIONS_LENGTH), NORTH)
#ambulatory_inner_pillar(u0() + vx(EQUIDISTANT_SECTIONS_LENGTH * 4) - vy(EQUIDISTANT_SECTIONS_LENGTH * 2), NORTH)

pillars = Array{Union{Pillar, Nothing}}(nothing, 11, 18)

function get_pillar_coordinates(row, column)

    return xy(EQUIDISTANT_SECTIONS_LENGTH * (column - 1), 
                EQUIDISTANT_SECTIONS_LENGTH * (-row + 1))
end

function instantiate_crossing_pillars(pillars)
    for pillar_info in CROSSING_PILLARS_INFO
        row = pillar_info[ROW_INDEX]
        column = pillar_info[COLUMN_INDEX]
        center = get_pillar_coordinates(row, column)
        orientation = pillar_info[ORIENTATION_INDEX]
        type = crossing_pillar(center)
        
        pillars[row, column] = Pillar(type, row, column, center, orientation)
    end
end

function instantiate_outer_crossing_buttresses(pillars)
    for pillar_info in OUTER_CROSSING_BUTTRESSES_INFO
        row = pillar_info[ROW_INDEX]
        column = pillar_info[COLUMN_INDEX]
        center = get_pillar_coordinates(row, column)
        orientation = pillar_info[ORIENTATION_INDEX]
        type = outer_crossing_buttress(center, orientation)
        
        pillars[row, column] = Pillar(type, row, column, center, orientation)
    end
end

function instantiate_aisle_buttresses(pillars)
end

instantiate_crossing_pillars(pillars)
instantiate_outer_crossing_buttresses(pillars)

#function get_pillar_orientation(row, column, section)
#    if section == CROSSING_ENDS_ID
#        return column in range(CROSSING_ENDS_COLUMN_RANGE[START_INDEX], CROSSING_ENDS_COLUMN_RANGE[])
#    end
#end
#
#function instantiate_crossing_ends(pillars)
#    for row = LOWER_CROSSING_END_ROW_RANGE
#        for column = CROSSING_ENDS_COLUMN_RANGE
#            center = get_pillar_coordinates(row, center)
#            aisle_buttress(center, 100, )
#        end
#    end
#end
#
#function instantiate_west_end(pillars)
#end
#
#function instantiate_aisle(pillars)
#    for row = aisle_ROW_RANGE
#        for column = aisle_COLUMN_RANGE
#        end
#    end
#    return pillars
#end
#
#function instantiate_ambulatory(pillars)
#end
# == INSTANTIATORS == #
# == PILLARS == #
