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

# == MEASUREMENTS == #
EQUIDISTANT_SECTIONS_LENGTH = 7.53
WEST_END_JUMP_COLUMN_INDEX = 3
AMBULATORY_SECTION_START_COLUMN_INDEX = 16
AMBULATORY_CENTER_X = EQUIDISTANT_SECTIONS_LENGTH * (AMBULATORY_SECTION_START_COLUMN_INDEX - 1)
AMBULATORY_CENTER_Y = EQUIDISTANT_SECTIONS_LENGTH * (-6 + 1)
AMBULATORY_CENTER = xy(AMBULATORY_CENTER_X, AMBULATORY_CENTER_Y)
AMBULATORY_START_ANGLE = π/2
AMBULATORY_ANGLE_INCREMENT = AMBULATORY_START_ANGLE / 4
# == MEASUREMENTS == #

# == PILLARS INFO == #
# == BASIC PILLAR == #
BASIC_PILLAR_WIDTH = 1
BASIC_PILLAR_DEPTH = 1
BASIC_PILLAR_HEIGHT = 105
# == BASIC PILLAR == #

# == WEST END PILLARS == #
N_WEST_END_PILLARS = (range(3, 5), range(1, 3), WEST)
S_WEST_END_PILLARS = (range(7, 9), range(1, 3), WEST)
# == WEST END PILLARS == #

# == AISLE BUTTRESSES == #
WN_AISLE_BUTTRESSES_INFO = (range(1, 2), 8, WEST)
WS_AISLE_BUTTRESSES_INFO = (range(10, 11), 8, WEST)
EN_AISLE_BUTTRESSES_INFO = (range(1, 2), 12, EAST)
ES_AISLE_BUTTRESSES_INFO = (range(10, 11), 12, EAST)
AISLE_BUTTRESSES_DIFFERENT_ROWS_INFO = [WN_AISLE_BUTTRESSES_INFO, WS_AISLE_BUTTRESSES_INFO, 
                                            EN_AISLE_BUTTRESSES_INFO, ES_AISLE_BUTTRESSES_INFO]
NW_AISLE_BUTTRESSES_INFO = (3, range(4, 7), NORTH)
NE_AISLE_BUTTRESSES_INFO = (3, range(13, 15), NORTH)
SW_AISLE_BUTTRESSES_INFO = (9, range(4, 7), SOUTH)
SE_AISLE_BUTTRESSES_INFO = (9, range(13, 15), SOUTH)
AISLE_BUTTRESSES_DIFFERENT_COLUMNS_INFO = [NW_AISLE_BUTTRESSES_INFO, NE_AISLE_BUTTRESSES_INFO, 
                                            SW_AISLE_BUTTRESSES_INFO, SE_AISLE_BUTTRESSES_INFO]
# == AISLE BUTTRESSES == #

# == AISLE OUTER PILLARS == #
NW_AISLE_OUTER_PILLARS = (4, range(4, 8), NORTH)
NE_AISLE_OUTER_PILLARS = (4, range(12, 15), NORTH)
SW_AISLE_OUTER_PILLARS = (8, range(4, 8), SOUTH)
SE_AISLE_OUTER_PILLARS = (8, range(12, 15), SOUTH)
AISLE_OUTER_PILLARS_INFO = [NW_AISLE_OUTER_PILLARS, NE_AISLE_OUTER_PILLARS, 
                                SW_AISLE_OUTER_PILLARS, SE_AISLE_OUTER_PILLARS]
# == AISLE OUTER PILLARS == #

# == AISLE INNER PILLARS == #
WN_AISLE_INNER_PILLARS = (range(1, 4), 9, WEST)
EN_AISLE_INNER_PILLARS = (range(1, 4), 11, EAST)
WS_AISLE_INNER_PILLARS = (range(8, 11), 9, WEST)
ES_AISLE_INNER_PILLARS = (range(8, 11), 11, EAST)
AISLE_INNER_PILLARS_DIFFERENT_ROWS_INFO = [WN_AISLE_INNER_PILLARS, EN_AISLE_INNER_PILLARS, 
                                            WS_AISLE_INNER_PILLARS, ES_AISLE_INNER_PILLARS]
NW_AISLE_INNER_PILLARS = (5, range(4, 8), NORTH)
NE_AISLE_INNER_PILLARS = (5, range(12, 15), NORTH)
SW_AISLE_INNER_PILLARS = (7, range(4, 8), SOUTH)
SE_AISLE_INNER_PILLARS = (7, range(12, 15), SOUTH)
AISLE_INNER_PILLARS_DIFFERENT_COLUMNS_INFO = [NW_AISLE_INNER_PILLARS, NE_AISLE_INNER_PILLARS, 
                                                SW_AISLE_INNER_PILLARS, SE_AISLE_INNER_PILLARS]
# == AISLE INNER PILLARS == #

# == OUTER CROSSING BUTTRESSES == #
UPPER_RIGHT_OUTER_CROSSING_BUTTRESS_INFO = (3, 12, EAST)
UPPER_LEFT_OUTER_CROSSING_BUTTRESS_INFO = (3, 8, NORTH)
LOWER_LEFT_OUTER_CROSSING_BUTTRESS_INFO = (9, 8, WEST)
LOWER_RIGHT_OUTER_CROSSING_BUTTRESS_INFO = (9, 12, SOUTH)
OUTER_CROSSING_BUTTRESSES_INFO = [UPPER_RIGHT_OUTER_CROSSING_BUTTRESS_INFO, UPPER_LEFT_OUTER_CROSSING_BUTTRESS_INFO,
                                    LOWER_LEFT_OUTER_CROSSING_BUTTRESS_INFO, LOWER_RIGHT_OUTER_CROSSING_BUTTRESS_INFO]
# == OUTER CROSSING BUTTRESSES == #

# == CROSSING PILLARS == #
UPPER_RIGHT_CROSSING_PILLAR_INFO = (5, 11, EAST)
UPPER_LEFT_CROSSING_PILLAR_INFO = (5, 9, WEST)
LOWER_LEFT_CROSSING_PILLAR_INFO = (7, 9, WEST)
LOWER_RIGHT_CROSSING_PILLAR_INFO = (7, 11, EAST)
CROSSING_PILLARS_INFO = [UPPER_RIGHT_CROSSING_PILLAR_INFO, UPPER_LEFT_CROSSING_PILLAR_INFO, 
                            LOWER_LEFT_CROSSING_PILLAR_INFO, LOWER_RIGHT_CROSSING_PILLAR_INFO]
# == CROSSING PILLARS == #

# == AMBULATORY BUTTRESSES == #
N_AMBULATORY_BUTTRESSES = (3, range(16, 18), nothing)
S_AMBULATORY_BUTTRESSES = (9, range(16, 18), nothing)
AMBULATORY_BUTTRESSES_INFO = [N_AMBULATORY_BUTTRESSES, S_AMBULATORY_BUTTRESSES]
# == AMBULATORY BUTTRESSES == #

# == AMBULATORY OUTER PILLARS == #
N_AMBULATORY_OUTER_PILLARS = (4, range(16, 18), nothing)
S_AMBULATORY_OUTER_PILLARS = (8, range(16, 18), nothing)
AMBULATORY_OUTER_PILLARS_INFO = [N_AMBULATORY_OUTER_PILLARS, S_AMBULATORY_OUTER_PILLARS]
# == AMBULATORY OUTER PILLARS == #

# == AMBULATORY INNER PILLARS == #
N_AMBULATORY_INNER_PILLARS = (5, range(16, 18), nothing)
S_AMBULATORY_INNER_PILLARS = (7, range(16, 18), nothing)
AMBULATORY_INNER_PILLARS_INFO = [N_AMBULATORY_INNER_PILLARS, S_AMBULATORY_INNER_PILLARS]
# == AMBULATORY INNER PILLARS == #
# == PILLARS INFO == #

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
function get_orientation_polar_angle(orientation::Union{VX, VXY})
    return atan(orientation.y, orientation.x)
end

function get_orientation_polar_angle(center::Union{X, XY})
    ambulatory_center_to_center_vector = center - AMBULATORY_CENTER
    return atan(ambulatory_center_to_center_vector.y, ambulatory_center_to_center_vector.x)
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

function west_end_pillar(center, orientation;
                            height = BASIC_PILLAR_HEIGHT * 0.475)
    orientation = get_orientation_polar_angle(orientation)
    west_end_pillar = basic_pillar(center; height)
    return rotate(west_end_pillar, orientation, center)
end

function aisle_buttress(center, orientation; 
                            height = BASIC_PILLAR_HEIGHT * 0.95)
    orientation = get_orientation_polar_angle(orientation)
    half_height = height / 2

    eighth_width = BASIC_PILLAR_WIDTH / 8
    half_width = BASIC_PILLAR_WIDTH / 2
    quarter_width = BASIC_PILLAR_WIDTH / 4

    half_depth = BASIC_PILLAR_DEPTH / 2
    quadruple_depth = BASIC_PILLAR_DEPTH * 4

    pillar = basic_pillar(center; height = height)

    base_bottom_left_corner = center - vx(BASIC_PILLAR_WIDTH)
    base_upper_right_corner = center + vxy(BASIC_PILLAR_WIDTH, quadruple_depth)

    pre_mid_bottom_left_corner = base_bottom_left_corner + vz(half_height - eighth_width)
    pre_mid_upper_right_corner = base_upper_right_corner + vz(half_height - eighth_width)
    
    mid_bottom_left_corner = base_bottom_left_corner - vx(quarter_width) + vz(half_height)
    mid_upper_right_corner = base_upper_right_corner + vxyz(quarter_width, half_depth, half_height)

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

function aisle_outer_pillar(center, orientation;
                                height = BASIC_PILLAR_HEIGHT)
    orientation = get_orientation_polar_angle(orientation)
    aisle_outer_pillar = basic_pillar(center; height = height)
    return rotate(aisle_outer_pillar, orientation, center)
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

function half_outer_crossing_buttress(center, height)
    half_height = height / 2

    half_width = BASIC_PILLAR_WIDTH / 2
    quarter_width = BASIC_PILLAR_WIDTH / 4
    eighth_width = BASIC_PILLAR_WIDTH / 8
    quadruple_width = BASIC_PILLAR_WIDTH * 4

    half_depth = BASIC_PILLAR_DEPTH / 2

    base_bottom_left_corner = center - vxy(half_width, half_depth)
    base_upper_right_corner = base_bottom_left_corner + vxy(quadruple_width, BASIC_PILLAR_DEPTH)

    pre_mid_bottom_left_corner = base_bottom_left_corner + vz(half_height - eighth_width)
    pre_mid_upper_right_corner = base_upper_right_corner + vz(half_height - eighth_width)

    mid_bottom_left_corner = base_bottom_left_corner + vz(half_height)
    mid_upper_right_corner = base_upper_right_corner + vxyz(quarter_width, half_depth, half_height)

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

function crossing_pillar(center, orientation; 
                            height = BASIC_PILLAR_HEIGHT * 1.15)
    orientation = get_orientation_polar_angle(orientation)
    three_quarters_width = BASIC_PILLAR_WIDTH * 0.75
    three_quarters_depth = BASIC_PILLAR_DEPTH * 0.75

    bottom_left_corner = u0() - vxy(three_quarters_width, three_quarters_depth)
    upper_right_corner = u0() + vxy(three_quarters_width, three_quarters_depth)

    pillar_shape = surface_rectangle(bottom_left_corner, upper_right_corner)
    pillar_path = line(center, center + vz(height))

    crossing_pillar = sweep(pillar_path, pillar_shape)

    return rotate(crossing_pillar, orientation, center)
end

function ambulatory_buttress(center; 
                                height = BASIC_PILLAR_HEIGHT * 0.64)
    orientation = get_orientation_polar_angle(center)

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

function ambulatory_outer_pillar(center; 
                                    height = BASIC_PILLAR_HEIGHT)
    orientation = get_orientation_polar_angle(center)

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

function ambulatory_inner_pillar(center; 
                                    height = BASIC_PILLAR_HEIGHT * 1.10)
    orientation = get_orientation_polar_angle(center)

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
# == MODELERS == #

# == INSTANTIATORS == #
# For demonstration purposes only
#crossing_pillar(u0() - vxy(EQUIDISTANT_SECTIONS_LENGTH, EQUIDISTANT_SECTIONS_LENGTH * 2))
#outer_crossing_buttress(u0(), EAST)
#aisle_buttress(u0() + vx(EQUIDISTANT_SECTIONS_LENGTH * 2), NORTH)
#basic_pillar(u0() + vx(EQUIDISTANT_SECTIONS_LENGTH * 2) - vy(EQUIDISTANT_SECTIONS_LENGTH))
#aisle_inner_pillar(u0() + vx(EQUIDISTANT_SECTIONS_LENGTH * 2) - vy(EQUIDISTANT_SECTIONS_LENGTH * 2), NORTH)
#ambulatory_buttress(u0() + vx(EQUIDISTANT_SECTIONS_LENGTH * 4), NORTH)
#ambulatory_outer_pillar(u0() + vx(EQUIDISTANT_SECTIONS_LENGTH * 4) - vy(EQUIDISTANT_SECTIONS_LENGTH), NORTH)
#ambulatory_inner_pillar(u0() + vx(EQUIDISTANT_SECTIONS_LENGTH * 4) - vy(EQUIDISTANT_SECTIONS_LENGTH * 2), NORTH)

function get_ambulatory_pillar_coordinates(row, column)
    if row < 6
        polar_angle = AMBULATORY_START_ANGLE - AMBULATORY_ANGLE_INCREMENT * (column - 15)
    else
        polar_angle = -AMBULATORY_START_ANGLE + AMBULATORY_ANGLE_INCREMENT * (column - 15)
    end

    if row == 3 || row == 9
        distance = EQUIDISTANT_SECTIONS_LENGTH * 2 + EQUIDISTANT_SECTIONS_LENGTH / 2
    elseif row == 4 || row == 8
        distance = EQUIDISTANT_SECTIONS_LENGTH * 2
    elseif row == 5 || row == 7
        distance = EQUIDISTANT_SECTIONS_LENGTH
    end

    return AMBULATORY_CENTER + vpol(distance, polar_angle)
end

function get_pillar_coordinates(row, column)
    x = EQUIDISTANT_SECTIONS_LENGTH * (column - 1)
    y = EQUIDISTANT_SECTIONS_LENGTH * (-row + 1)

    if column >= WEST_END_JUMP_COLUMN_INDEX && column < AMBULATORY_SECTION_START_COLUMN_INDEX
        x += EQUIDISTANT_SECTIONS_LENGTH / 2
    elseif column >= AMBULATORY_SECTION_START_COLUMN_INDEX
        return get_ambulatory_pillar_coordinates(row, column)
    end

    return xy(x, y)
end

function rangeless_instantiator(pillars, pillar_info_collection, pillar_type_instantiator)
    for pillar_info in pillar_info_collection
        row = pillar_info[ROW_INDEX]
        column = pillar_info[COLUMN_INDEX]
        center = get_pillar_coordinates(row, column)
        orientation = pillar_info[ORIENTATION_INDEX]
        type = pillar_type_instantiator(center, orientation)
        pillars[row, column] = Pillar(type, row, column, center, orientation)
    end  
end

function row_range_instantiator(pillars, pillar_info_collection, pillar_type_instantiator)
    for pillar_info in pillar_info_collection
        row_range = pillar_info[ROW_INDEX]
        column = pillar_info[COLUMN_INDEX]

        for row in row_range
            center = get_pillar_coordinates(row, column)
            orientation = pillar_info[ORIENTATION_INDEX]
            type = pillar_type_instantiator(center, orientation)
            pillars[row, column] = Pillar(type, row, column, center, orientation)
        end
    end
end

function column_range_instantiator(pillars, pillar_info_collection, pillar_type_instantiator)
    for pillar_info in pillar_info_collection
        row = pillar_info[ROW_INDEX]
        column_range = pillar_info[COLUMN_INDEX]

        for column in column_range
            center = get_pillar_coordinates(row, column)
            orientation = pillar_info[ORIENTATION_INDEX]

            if isnothing(orientation)
                type = pillar_type_instantiator(center)
            else
                type = pillar_type_instantiator(center, orientation)
            end

            pillars[row, column] = Pillar(type, row, column, center, orientation)
        end
    end
end

function instantiate_west_end_pillars(pillars)
    row_range = N_WEST_END_PILLARS[ROW_INDEX]
    column_range = N_WEST_END_PILLARS[COLUMN_INDEX]
    orientation = N_WEST_END_PILLARS[ORIENTATION_INDEX]

    for row in row_range
        for column in column_range
            center = get_pillar_coordinates(row, column)
            type = west_end_pillar(center, orientation)
            pillars[row, column] = Pillar(type, row, column, center, orientation)
        end
    end

    row_range = S_WEST_END_PILLARS[ROW_INDEX]
    column_range = S_WEST_END_PILLARS[COLUMN_INDEX]
    orientation = S_WEST_END_PILLARS[ORIENTATION_INDEX]

    for row in row_range
        for column in column_range
            center = get_pillar_coordinates(row, column)
            type = west_end_pillar(center, orientation)
            pillars[row, column] = Pillar(type, row, column, center, orientation)
        end
    end
end

function instantiate_aisle_buttresses(pillars)
    row_range_instantiator(pillars, AISLE_BUTTRESSES_DIFFERENT_ROWS_INFO, aisle_buttress)
    column_range_instantiator(pillars, AISLE_BUTTRESSES_DIFFERENT_COLUMNS_INFO, aisle_buttress)
end

function instantiate_aisle_outer_pillars(pillars)
    column_range_instantiator(pillars, AISLE_OUTER_PILLARS_INFO, aisle_outer_pillar)
end

function instantiate_aisle_inner_pillars(pillars)
    row_range_instantiator(pillars, AISLE_INNER_PILLARS_DIFFERENT_ROWS_INFO, aisle_inner_pillar)
    column_range_instantiator(pillars, AISLE_INNER_PILLARS_DIFFERENT_COLUMNS_INFO, aisle_inner_pillar)
end

function instantiate_outer_crossing_buttresses(pillars)
    rangeless_instantiator(pillars, OUTER_CROSSING_BUTTRESSES_INFO, outer_crossing_buttress)
end

function instantiate_crossing_pillars(pillars)
    rangeless_instantiator(pillars, CROSSING_PILLARS_INFO, crossing_pillar)
end

function instantiate_ambulatory_buttresses(pillars)
    column_range_instantiator(pillars, AMBULATORY_BUTTRESSES_INFO, ambulatory_buttress)
end

function instantiate_ambulatory_outer_pillars(pillars)
    column_range_instantiator(pillars, AMBULATORY_OUTER_PILLARS_INFO, ambulatory_outer_pillar)
end

function instantiate_ambulatory_inner_pillars(pillars)
    column_range_instantiator(pillars, AMBULATORY_INNER_PILLARS_INFO, ambulatory_inner_pillar)
end

function instantiate_all_pillars(pillars)
    instantiate_west_end_pillars(pillars)
    instantiate_aisle_buttresses(pillars)
    instantiate_aisle_outer_pillars(pillars)
    instantiate_aisle_inner_pillars(pillars)
    instantiate_outer_crossing_buttresses(pillars)
    instantiate_crossing_pillars(pillars)
    instantiate_ambulatory_buttresses(pillars)    
    instantiate_ambulatory_outer_pillars(pillars)
    instantiate_ambulatory_inner_pillars(pillars)
end
# == INSTANTIATORS == #
# == PILLARS == #

# == PLAYGROUND == #
pillars = Array{Union{Pillar, Nothing}}(nothing, 11, 18)

instantiate_all_pillars(pillars)
# == PLAYGROUND == #
