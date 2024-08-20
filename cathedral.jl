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

# == DIRECTIONS == #
NORTH = vy(1)
SOUTH = -NORTH
WEST = vx(-1)
EAST = -WEST

GROWING_HEIGHT_DIRECTION = vz(1)
DECREASING_HEIGHT_DIRECTION = -GROWING_HEIGHT_DIRECTION
# == DIRECTIONS == #

# == MEASUREMENTS & DELIMIETERS == #
EQUIDISTANT_SECTIONS_LENGTH = 7.53 * 2

WEST_END_JUMP_COLUMN_INDEX = 3
MIDDLE_HALLWAY_JUMP_ROW_INDEX = 6
MIDDLE_HALLWAY_JUMP_COLUMN_INDEX = 10
AMBULATORY_SECTION_START_COLUMN_INDEX = 15

AMBULATORY_CENTER_X = EQUIDISTANT_SECTIONS_LENGTH * AMBULATORY_SECTION_START_COLUMN_INDEX
AMBULATORY_CENTER_Y = EQUIDISTANT_SECTIONS_LENGTH * -5
AMBULATORY_CENTER = xy(AMBULATORY_CENTER_X, AMBULATORY_CENTER_Y)

AMBULATORY_START_ANGLE = π/2
AMBULATORY_ANGLE_INCREMENT = AMBULATORY_START_ANGLE / 4
# == MEASUREMENTS & DELIMITERS == #

# == PILLARS INFO == #
# == BASIC PILLAR == #
BASIC_PILLAR_WIDTH = 1
BASIC_PILLAR_DEPTH = 1
BASIC_PILLAR_HEIGHT = 105
# == BASIC PILLAR == #

# == WEST END PILLARS == #
WEST_END_PILLAR_WIDTH = BASIC_PILLAR_WIDTH * 1.5
WEST_END_PILLAR_DEPTH = BASIC_PILLAR_DEPTH * 1.5
WEST_END_PILLAR_HEIGHT = BASIC_PILLAR_HEIGHT * 0.475

N_WEST_END_PILLARS = (range(3, 5), range(1, 3), WEST)
S_WEST_END_PILLARS = (range(6, 8), range(1, 3), WEST)
# == WEST END PILLARS == #

# == AISLE BUTTRESSES == #
AISLE_BUTTRESS_WIDTH = BASIC_PILLAR_WIDTH
AISLE_BUTTRESS_DEPTH = BASIC_PILLAR_DEPTH
AISLE_BUTTRESS_HEIGHT = BASIC_PILLAR_HEIGHT * 0.95

WN_AISLE_BUTTRESSES_INFO = (range(1, 2), 8, WEST)
EN_AISLE_BUTTRESSES_INFO = (range(1, 2), 11, EAST)
WS_AISLE_BUTTRESSES_INFO = (range(9, 10), 8, WEST)
ES_AISLE_BUTTRESSES_INFO = (range(9, 10), 11, EAST)
NW_AISLE_BUTTRESSES_INFO = (3, range(4, 7), NORTH)
NE_AISLE_BUTTRESSES_INFO = (3, range(12, 14), NORTH)
SW_AISLE_BUTTRESSES_INFO = (8, range(4, 7), SOUTH)
SE_AISLE_BUTTRESSES_INFO = (8, range(12, 14), SOUTH)

AISLE_BUTTRESSES_DIFFERENT_ROWS_INFO = [WN_AISLE_BUTTRESSES_INFO, WS_AISLE_BUTTRESSES_INFO, 
                                            EN_AISLE_BUTTRESSES_INFO, ES_AISLE_BUTTRESSES_INFO]
AISLE_BUTTRESSES_DIFFERENT_COLUMNS_INFO = [NW_AISLE_BUTTRESSES_INFO, NE_AISLE_BUTTRESSES_INFO, 
                                            SW_AISLE_BUTTRESSES_INFO, SE_AISLE_BUTTRESSES_INFO]
# == AISLE BUTTRESSES == #

# == AISLE OUTER PILLARS == #
AISLE_OUTER_PILLAR_WIDTH = AISLE_BUTTRESS_WIDTH
AISLE_OUTER_PILLAR_DEPTH = AISLE_BUTTRESS_DEPTH
AISLE_OUTER_PILLAR_HEIGHT = BASIC_PILLAR_HEIGHT * 1.025

NW_AISLE_OUTER_PILLARS = (4, range(4, 8), NORTH)
NE_AISLE_OUTER_PILLARS = (4, range(11, 14), NORTH)
SW_AISLE_OUTER_PILLARS = (7, range(4, 8), SOUTH)
SE_AISLE_OUTER_PILLARS = (7, range(11, 14), SOUTH)

AISLE_OUTER_PILLARS_INFO = [NW_AISLE_OUTER_PILLARS, NE_AISLE_OUTER_PILLARS, 
                                SW_AISLE_OUTER_PILLARS, SE_AISLE_OUTER_PILLARS]
# == AISLE OUTER PILLARS == #

# == AISLE INNER PILLARS == #
AISLE_INNER_PILLAR_WIDTH = AISLE_OUTER_PILLAR_WIDTH * 1.5
AISLE_INNER_PILLAR_DEPTH = AISLE_OUTER_PILLAR_DEPTH
AISLE_INNER_PILLAR_HEIGHT = BASIC_PILLAR_HEIGHT * 1.10

WN_AISLE_INNER_PILLARS = (range(1, 4), 9, WEST)
EN_AISLE_INNER_PILLARS = (range(1, 4), 10, EAST)
WS_AISLE_INNER_PILLARS = (range(7, 10), 9, WEST)
ES_AISLE_INNER_PILLARS = (range(7, 10), 10, EAST)
NW_AISLE_INNER_PILLARS = (5, range(4, 8), NORTH)
NE_AISLE_INNER_PILLARS = (5, range(11, 14), NORTH)
SW_AISLE_INNER_PILLARS = (6, range(4, 8), SOUTH)
SE_AISLE_INNER_PILLARS = (6, range(11, 14), SOUTH)

AISLE_INNER_PILLARS_DIFFERENT_ROWS_INFO = [WN_AISLE_INNER_PILLARS, EN_AISLE_INNER_PILLARS, 
                                            WS_AISLE_INNER_PILLARS, ES_AISLE_INNER_PILLARS]
AISLE_INNER_PILLARS_DIFFERENT_COLUMNS_INFO = [NW_AISLE_INNER_PILLARS, NE_AISLE_INNER_PILLARS, 
                                                SW_AISLE_INNER_PILLARS, SE_AISLE_INNER_PILLARS]
# == AISLE INNER PILLARS == #

# == OUTER CROSSING BUTTRESSES == #
OUTER_CROSSING_BUTTRESS_WIDTH = BASIC_PILLAR_WIDTH
OUTER_CROSSING_BUTTRESS_DEPTH = BASIC_PILLAR_DEPTH
OUTER_CROSSING_BUTTRESS_HEIGHT = AISLE_BUTTRESS_HEIGHT

UPPER_RIGHT_OUTER_CROSSING_BUTTRESS_INFO = (3, 11, EAST)
UPPER_LEFT_OUTER_CROSSING_BUTTRESS_INFO = (3, 8, NORTH)
LOWER_LEFT_OUTER_CROSSING_BUTTRESS_INFO = (8, 8, WEST)
LOWER_RIGHT_OUTER_CROSSING_BUTTRESS_INFO = (8, 11, SOUTH)

OUTER_CROSSING_BUTTRESSES_INFO = [UPPER_RIGHT_OUTER_CROSSING_BUTTRESS_INFO, UPPER_LEFT_OUTER_CROSSING_BUTTRESS_INFO,
                                    LOWER_LEFT_OUTER_CROSSING_BUTTRESS_INFO, LOWER_RIGHT_OUTER_CROSSING_BUTTRESS_INFO]
# == OUTER CROSSING BUTTRESSES == #

# == CROSSING PILLARS == #
CROSSING_PILLAR_WIDTH = AISLE_INNER_PILLAR_WIDTH
CROSSING_PILLAR_DEPTH = AISLE_INNER_PILLAR_WIDTH
CROSSING_PILLAR_HEIGHT = AISLE_INNER_PILLAR_HEIGHT

UPPER_RIGHT_CROSSING_PILLAR_INFO = (5, 10, EAST)
UPPER_LEFT_CROSSING_PILLAR_INFO = (5, 9, WEST)
LOWER_LEFT_CROSSING_PILLAR_INFO = (6, 9, WEST)
LOWER_RIGHT_CROSSING_PILLAR_INFO = (6, 10, EAST)

CROSSING_PILLARS_INFO = [UPPER_RIGHT_CROSSING_PILLAR_INFO, UPPER_LEFT_CROSSING_PILLAR_INFO, 
                            LOWER_LEFT_CROSSING_PILLAR_INFO, LOWER_RIGHT_CROSSING_PILLAR_INFO]
# == CROSSING PILLARS == #

# == AMBULATORY BUTTRESSES == #
AMBULATORY_BUTTRESS_WIDTH = 0
AMBULATORY_BUTTRESS_DEPTH = 0
AMBULATORY_BUTTRESS_HEIGHT = BASIC_PILLAR_HEIGHT * 0.64

N_AMBULATORY_BUTTRESSES = (3, range(15, 17), nothing)
S_AMBULATORY_BUTTRESSES = (8, range(15, 17), nothing)

AMBULATORY_BUTTRESSES_INFO = [N_AMBULATORY_BUTTRESSES, S_AMBULATORY_BUTTRESSES]
# == AMBULATORY BUTTRESSES == #

# == AMBULATORY OUTER PILLARS == #
AMBULATORY_OUTER_PILLAR_WIDTH = AISLE_OUTER_PILLAR_WIDTH
AMBULATORY_OUTER_PILLAR_DEPTH = AISLE_OUTER_PILLAR_DEPTH
AMBULATORY_OUTER_PILLAR_HEIGHT = AISLE_OUTER_PILLAR_HEIGHT

N_AMBULATORY_OUTER_PILLARS = (4, range(15, 17), nothing)
S_AMBULATORY_OUTER_PILLARS = (7, range(15, 17), nothing)

AMBULATORY_OUTER_PILLARS_INFO = [N_AMBULATORY_OUTER_PILLARS, S_AMBULATORY_OUTER_PILLARS]
# == AMBULATORY OUTER PILLARS == #

# == AMBULATORY INNER PILLARS == #
AMBULATORY_INNER_PILLAR_WIDTH = 0
AMBULATORY_INNER_PILLAR_DEPTH = 0
AMBULATORY_INNER_PILLAR_HEIGHT = AISLE_INNER_PILLAR_HEIGHT

N_AMBULATORY_INNER_PILLARS = (5, range(15, 17), nothing)
S_AMBULATORY_INNER_PILLARS = (6, range(15, 17), nothing)

AMBULATORY_INNER_PILLARS_INFO = [N_AMBULATORY_INNER_PILLARS, S_AMBULATORY_INNER_PILLARS]
# == AMBULATORY INNER PILLARS == #
# == PILLARS INFO == #

# == WALLS INFO == #
# == FLYING BUTTRESSES == #
NW_FLYING_BUTTRESS_LEFT_PILLARS = (3, range(4, 8))
NW_FLYING_BUTTRESS_RIGHT_PILLARS = (4, range(4, 8))
NE_FLYING_BUTTRESS_LEFT_PILLARS = (3, range(11, 14))
NE_FLYING_BUTTRESS_RIGHT_PILLARS = (4, range(11, 14))
N_FLYING_BUTTRESSES_LEFT_ANCHOR_HEIGHT = AISLE_BUTTRESS_HEIGHT / 2 + AISLE_BUTTRESS_HEIGHT / 6
N_FLYING_BUTTRESSES_RIGHT_ANCHOR_HEIGHT = AISLE_OUTER_PILLAR_HEIGHT / 2 + AISLE_OUTER_PILLAR_HEIGHT / 4

SW_FLYING_BUTTRESS_LEFT_PILLARS = (7, range(4, 8))
SW_FLYING_BUTTRESS_RIGHT_PILLARS = (8, range(4, 8))
SE_FLYING_BUTTRESS_LEFT_PILLARS = (7, range(11, 14))
SE_FLYING_BUTTRESS_RIGHT_PILLARS = (8, range(11, 14))
S_FLYING_BUTTRESSES_LEFT_ANCHOR_HEIGHT = N_FLYING_BUTTRESSES_RIGHT_ANCHOR_HEIGHT
S_FLYING_BUTTRESSES_RIGHT_ANCHOR_HEIGHT = N_FLYING_BUTTRESSES_LEFT_ANCHOR_HEIGHT
# == FLYING BUTTRESSES == #

# == MAIN HALLWAY WALLS == #
MAIN_HALLWAYS_WALL_HEIGHT = AISLE_INNER_PILLAR_HEIGHT

HORIZONTAL_HALLWAY_WALL_LEFT_PILLARS = (5, range(4, 14))
HORIZONTAL_HALLWAY_WALL_RIGHT_PILLARS = (6, range(4, 14))

VERTICAL_HALLWAY_WALL_LEFT_PILLARS = (range(1, 10), 9)
VERTICAL_HALLWAY_WALL_RIGHT_PILLARS = (range(1, 10), 10)
# == MAIN HALLWAY WALLS == #

# == AISLE OUTER WALLS == #
AISLE_OUTER_WALL_HEIGHT = AISLE_BUTTRESS_HEIGHT / 2

NW_AISLE_OUTER_WALL_LEFT_PILLARS = (3, range(4, 8))
NW_AISLE_OUTER_WALL_RIGHT_PILLARS = (4, range(4, 8))
NE_AISLE_OUTER_WALL_LEFT_PILLARS = (3, range(11, 14))
NE_AISLE_OUTER_WALL_RIGHT_PILLARS = (4, range(11, 14))

SW_AISLE_OUTER_WALL_LEFT_PILLARS = (7, range(4, 8))
SW_AISLE_OUTER_WALL_RIGHT_PILLARS = (8, range(4, 8))
SE_AISLE_OUTER_WALL_LEFT_PILLARS = (7, range(11, 14))
SE_AISLE_OUTER_WALL_RIGHT_PILLARS = (8, range(11, 14))
# == AISLE OUTER WALLS == #

# == AISLE MIDDLE WALLS == #
AISLE_MIDDLE_WALL_HEIGHT = AISLE_OUTER_WALL_HEIGHT

NW_AISLE_MIDDLE_WALL_LEFT_PILLARS = (4, range(4, 8))
NW_AISLE_MIDDLE_WALL_RIGHT_PILLARS = (4, range(5, 9))
NE_AISLE_MIDDLE_WALL_LEFT_PILLARS = (4, range(10, 13))
NE_AISLE_MIDDLE_WALL_RIGHT_PILLARS = (4, range(11, 14))

SW_AISLE_MIDDLE_WALL_LEFT_PILLARS = (7, range(4, 8))
SW_AISLE_MIDDLE_WALL_RIGHT_PILLARS = (7, range(5, 9))
SE_AISLE_MIDDLE_WALL_LEFT_PILLARS = (7, range(10, 13))
SE_AISLE_MIDDLE_WALL_RIGHT_PILLARS = (7, range(11, 14))
# == AISLE MIDDLE WALLS == #

# == AISLE INNER WALLS == #
AISLE_INNER_WALL_HEIGHT = AISLE_OUTER_WALL_HEIGHT

NW_AISLE_INNER_WALL_LEFT_PILLARS = (4, range(4, 8))
NW_AISLE_INNER_WALL_RIGHT_PILLARS = (5, range(4, 8))
NE_AISLE_INNER_WALL_LEFT_PILLARS = (4, range(11, 14))
NE_AISLE_INNER_WALL_RIGHT_PILLARS = (5, range(11, 14))

SW_AISLE_INNER_WALL_LEFT_PILLARS = (6, range(4, 8))
SW_AISLE_INNER_WALL_RIGHT_PILLARS = (7, range(4, 8))
SE_AISLE_INNER_WALL_LEFT_PILLARS = (6, range(11, 14))
SE_AISLE_INNER_WALL_RIGHT_PILLARS = (7, range(11, 14))
# == AISLE INNER WALLS == #

# == TRANSEPT OUTER WALLS == #
TRANSEPT_OUTER_WALL_HEIGHT = AISLE_OUTER_WALL_HEIGHT

WN_TRANSEPT_OUTER_WALL_LEFT_PILLARS = (range(1, 3), 8)
WN_TRANSEPT_OUTER_WALL_RIGHT_PILLARS = (range(1, 3), 9)
EN_TRANSEPT_OUTER_WALL_LEFT_PILLARS = (range(1, 3), 10)
EN_TRANSEPT_OUTER_WALL_RIGHT_PILLARS = (range(1, 3), 11)

WS_TRANSEPT_OUTER_WALL_LEFT_PILLARS = (range(8, 10), 8)
WS_TRANSEPT_OUTER_WALL_RIGHT_PILLARS = (range(8, 10), 9)
ES_TRANSEPT_OUTER_WALL_LEFT_PILLARS = (range(8, 10), 10)
ES_TRANSEPT_OUTER_WALL_RIGHT_PILLARS = (range(8, 10), 11)
# == TRANSEPT OUTER WALLS == #

# == AMBULATORY OUTER WALLS == #
AMBULATORY_OUTER_WALL_HEIGHT = AISLE_OUTER_WALL_HEIGHT

TOP_AMBULATORY_OUTER_WALL_LEFT_PILLARS = (4, range(14, 17))
TOP_AMBULATORY_OUTER_WALL_RIGHT_PILLARS = (4, range(15, 18))

BOTTOM_AMBULATORY_OUTER_WALL_LEFT_PILLARS = (7, range(14, 16))
BOTTOM_AMBULATORY_OUTER_WALL_RIGHT_PILLARS = (7, range(15, 17))
# == AMBULATORY OUTER WALLS == #
# == WALLS INFO == #

# == CATHEDRAL ASSETS DATA STRUCTURES == #
mutable struct Pillar
    type
    row
    column
    center
    orientation
    width
    depth
    height
    north_wall
    south_wall
    west_wall
    east_wall
    ambulatory_oriented_wall
end

function Pillar(type, row, column, center, orientation, width, depth, height)
    Pillar(type, row, column, center, orientation, width, depth, height, nothing, nothing, nothing, nothing, nothing)
end

function get_pillar_center(pillar)
    return pillar.center
end

function get_pillar_orientation(pillar)
    return pillar.orientation
end

function get_pillar_width(pillar)
    return pillar.width
end

function get_pillar_depth(pillar)
    return pillar.depth
end

function get_pillar_height(pillar)
    return pillar.height
end

function get_pillar_north_wall(pillar)
    return pillar.north_wall
end

function get_pillar_south_wall(pillar)
    return pillar.south_wall
end

function get_pillar_west_wall(pillar)
    return pillar.west_wall
end

function get_pillar_east_wall(pillar)
    return pillar.east_wall
end

function set_pillar_north_wall(pillar, wall)
    pillar.north_wall = wall
end

function set_pillar_south_wall(pillar, wall)
    pillar.south_wall = wall
end

function set_pillar_west_wall(pillar, wall)
    pillar.west_wall = wall
end

function set_pillar_east_wall(pillar, wall)
    pillar.east_wall = wall
end

mutable struct Wall
    type
    left_pillar
    right_pillar
end

function get_wall_model(wall)
    return wall.type
end

function get_wall_left_pillar(wall)
    return wall.left_pillar
end

function get_wall_right_pillar(wall)
    return wall.right_pillar
end
# == CATHEDRAL ASSETS DATA STRUCTURES == #

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
                        width = BASIC_PILLAR_WIDTH, 
                        depth = BASIC_PILLAR_DEPTH, 
                        height = BASIC_PILLAR_HEIGHT)
    half_width = width / 2
    half_depth = depth / 2

    bottom_left_corner = u0() - vxy(half_width, half_depth)
    upper_right_corner = u0() + vxy(half_width, half_depth)

    pillar_shape = surface_rectangle(bottom_left_corner, upper_right_corner)
    pillar_path = line(center, center + vz(height))
    basic_pillar = rotate(sweep(pillar_path, pillar_shape), π/2, center)

    return basic_pillar
end

function west_end_pillar(center, orientation;
                            width = WEST_END_PILLAR_WIDTH,
                            depth = WEST_END_PILLAR_DEPTH,
                            height = WEST_END_PILLAR_HEIGHT)
    orientation = get_orientation_polar_angle(orientation)
    west_end_pillar = rotate(basic_pillar(center; width = width, depth = depth, height = height), orientation, center)
    return (model = west_end_pillar, width = width, depth = depth, height = height)
end

function aisle_buttress(center, orientation;
                            width = AISLE_BUTTRESS_WIDTH,
                            depth = AISLE_BUTTRESS_DEPTH, 
                            height = AISLE_BUTTRESS_HEIGHT)
    orientation = get_orientation_polar_angle(orientation)
    half_height = height / 2

    eighth_width = BASIC_PILLAR_WIDTH / 8
    half_width = BASIC_PILLAR_WIDTH / 2
    quarter_width = BASIC_PILLAR_WIDTH / 4

    half_depth = BASIC_PILLAR_DEPTH / 2
    quadruple_depth = BASIC_PILLAR_DEPTH * 4

    pillar = basic_pillar(center; width = width, depth = depth, height = height)

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

    aisle_buttress = rotate(union(pillar, loft_ruled([pillar_base, pillar_pre_mid, pillar_mid, pillar_post_mid, pillar_top])), -π/2 + orientation, center)

    return (model = aisle_buttress, width = width, depth = depth, height = height)
end

function aisle_outer_pillar(center, orientation;
                                width = AISLE_OUTER_PILLAR_WIDTH, 
                                depth = AISLE_OUTER_PILLAR_DEPTH, 
                                height = AISLE_OUTER_PILLAR_HEIGHT)
    orientation = get_orientation_polar_angle(orientation)
    aisle_outer_pillar = rotate(basic_pillar(center; width = width, depth = depth, height = height), orientation, center)
    return (model = aisle_outer_pillar, width = width, depth = depth, height = height)
end

function aisle_inner_pillar(center, orientation; 
                                width = AISLE_INNER_PILLAR_WIDTH,
                                depth = AISLE_INNER_PILLAR_DEPTH,
                                height = AISLE_INNER_PILLAR_HEIGHT)
    orientation = get_orientation_polar_angle(orientation)
    aisle_inner_pillar = rotate(basic_pillar(center; width = width, depth = depth, height = height), orientation, center)
    return (model = aisle_inner_pillar, width = width, depth = depth, height = height)
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
                                    width = OUTER_CROSSING_BUTTRESS_WIDTH,
                                    depth = OUTER_CROSSING_BUTTRESS_DEPTH,
                                    height = OUTER_CROSSING_BUTTRESS_HEIGHT)
    orientation = get_orientation_polar_angle(orientation)

    quarter_width = BASIC_PILLAR_WIDTH / 4
    three_quarters_width = quarter_width * 3

    quarter_depth = BASIC_PILLAR_DEPTH / 4
    three_quarters_depth = quarter_depth * 3

    pillar = basic_pillar(center; width = width, depth = depth, height = height)

    transformation_point = center + vxy(quarter_width, quarter_depth)
    right_buttress = move(half_outer_crossing_buttress(center, height), vxy(three_quarters_width, three_quarters_depth))
    left_buttress = mirror(rotate(deepcopy(right_buttress), π/2, transformation_point), transformation_point + vz())

    outer_crossing_buttress = rotate(union(pillar, right_buttress, left_buttress), orientation, center)

    return (model = outer_crossing_buttress, width = width, depth = depth, height = height)
end

function crossing_pillar(center, orientation; 
                            width = CROSSING_PILLAR_WIDTH,
                            depth = CROSSING_PILLAR_DEPTH,
                            height = CROSSING_PILLAR_HEIGHT)
    orientation = get_orientation_polar_angle(orientation)
    crossing_pillar = rotate(basic_pillar(center; width = width, depth = depth, height = height), orientation, center)
    return (model = crossing_pillar, width = width, depth = depth, height = height)
end

function ambulatory_buttress(center, orientation;
                                width = AISLE_BUTTRESS_WIDTH,
                                depth = AISLE_BUTTRESS_DEPTH, 
                                height = AISLE_BUTTRESS_HEIGHT)
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

    ambulatory_buttress = rotate(sweep(pillar_path, pillar_shape), -π + orientation, center)
    
    return (model = ambulatory_buttress, width = width, depth = depth, height = height)
end

function ambulatory_outer_pillar(center, orientation;
                                    width = AMBULATORY_OUTER_PILLAR_WIDTH,
                                    depth = AMBULATORY_OUTER_PILLAR_DEPTH,
                                    height = AMBULATORY_OUTER_PILLAR_HEIGHT)
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
    
    ambulatory_outer_pillar = rotate(sweep(pillar_path, pillar_shape), -π + orientation, center)

    return (model = ambulatory_outer_pillar, width = width, depth = depth, height = height)
end

function ambulatory_inner_pillar(center, orientation;
                                    width = AMBULATORY_INNER_PILLAR_WIDTH,
                                    depth = AMBULATORY_INNER_PILLAR_DEPTH,
                                    height = AMBULATORY_INNER_PILLAR_HEIGHT)
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
    
    ambulatory_inner_pillar = rotate(sweep(pillar_path, pillar_shape), -π + orientation, center)

    return (model = ambulatory_inner_pillar, width = width, depth = depth, height = height)
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
    if row < MIDDLE_HALLWAY_JUMP_ROW_INDEX
        polar_angle = AMBULATORY_START_ANGLE - AMBULATORY_ANGLE_INCREMENT * (column - 14)
    else
        polar_angle = -AMBULATORY_START_ANGLE + AMBULATORY_ANGLE_INCREMENT * (column - 14)
    end

    if row == 3 || row == 8
        distance = EQUIDISTANT_SECTIONS_LENGTH * 2 + EQUIDISTANT_SECTIONS_LENGTH / 2
    elseif row == 4 || row == 7
        distance = EQUIDISTANT_SECTIONS_LENGTH * 2
    elseif row == 5 || row == 6
        distance = EQUIDISTANT_SECTIONS_LENGTH
    end

    return AMBULATORY_CENTER + vpol(distance, polar_angle)
end

function get_pillar_coordinates(row, column)
    x = EQUIDISTANT_SECTIONS_LENGTH * (column - 1)
    y = EQUIDISTANT_SECTIONS_LENGTH * (-row + 1)

    if row >= MIDDLE_HALLWAY_JUMP_ROW_INDEX
        y -= EQUIDISTANT_SECTIONS_LENGTH
    end

    if column >= WEST_END_JUMP_COLUMN_INDEX && column < MIDDLE_HALLWAY_JUMP_COLUMN_INDEX
        x += EQUIDISTANT_SECTIONS_LENGTH / 2
    elseif column >= MIDDLE_HALLWAY_JUMP_COLUMN_INDEX && column < AMBULATORY_SECTION_START_COLUMN_INDEX
        x += EQUIDISTANT_SECTIONS_LENGTH / 2 + EQUIDISTANT_SECTIONS_LENGTH
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
        model = pillar_type_instantiator(center, orientation)
        type = model.model
        width = model.width
        depth = model.depth
        height = model.height
        pillars[row, column] = Pillar(type, row, column, center, orientation, width, depth, height)
    end  
end

function row_range_instantiator(pillars, pillar_info_collection, pillar_type_instantiator)
    for pillar_info in pillar_info_collection
        row_range = pillar_info[ROW_INDEX]
        column = pillar_info[COLUMN_INDEX]

        for row in row_range
            center = get_pillar_coordinates(row, column)
            orientation = pillar_info[ORIENTATION_INDEX]
            model = pillar_type_instantiator(center, orientation)
            type = model.model
            width = model.width
            depth = model.depth
            height = model.height
            pillars[row, column] = Pillar(type, row, column, center, orientation, width, depth, height)
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
            model = pillar_type_instantiator(center, orientation)
            type = model.model
            width = model.width
            depth = model.depth
            height = model.height
            pillars[row, column] = Pillar(type, row, column, center, orientation, width, depth, height)
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
            model = west_end_pillar(center, orientation)
            type = model.model
            width = model.width
            depth = model.depth
            height = model.height
            pillars[row, column] = Pillar(type, row, column, center, orientation, width, depth, height)
        end
    end

    row_range = S_WEST_END_PILLARS[ROW_INDEX]
    column_range = S_WEST_END_PILLARS[COLUMN_INDEX]
    orientation = S_WEST_END_PILLARS[ORIENTATION_INDEX]

    for row in row_range
        for column in column_range
            center = get_pillar_coordinates(row, column)
            model = west_end_pillar(center, orientation)
            type = model.model
            width = model.width
            depth = model.depth
            height = model.height
            pillars[row, column] = Pillar(type, row, column, center, orientation, width, depth, height)
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

# == WALLS == #
# == UTILITIES == #
function normalize_vector(v)
    return v / norm(v)
end

function displace_point_by_vector(point, vector, magnitude)
    unit_vector = vector / norm(vector)

    return point + unit_vector * magnitude
end

function arc(center, start_point, end_point)
    a = start_point - center
    b = end_point - center

    ccw_start_angle = atan(a.y, a.x)
    ccw_end_angle = atan(b.y, b.x)
    radius = norm(a)

    if ccw_end_angle > ccw_start_angle
        amplitude = ccw_end_angle - ccw_start_angle
    else
        amplitude = 2π - ccw_start_angle + ccw_end_angle
    end

    return KhepriAutoCAD.arc(center, radius, ccw_start_angle, amplitude)
end

function surface_arc(arc)
    center = arc_center(arc)
    radius = arc_radius(arc)
    start_angle = arc_start_angle(arc)
    amplitude = arc_amplitude(arc)

    return KhepriAutoCAD.surface_arc(center, radius, start_angle, amplitude)
end

function get_opposite_orientation(orientation)
    if orientation == NORTH
        return SOUTH
    elseif orientation == SOUTH
        return NORTH
    elseif orientation == WEST
        return EAST
    elseif orientation == EAST
        return WEST
    end

    return nothing
end

function get_ccw_perpendicular_orientation(orientation)
    if orientation == NORTH
        return WEST
    elseif orientation == SOUTH
        return EAST
    elseif orientation == WEST
        return SOUTH
    elseif orientation == EAST
        return NORTH
    end

    return nothing
end

function get_bidimensional_perpendicular_vector(v)
    return vxy(-v.y, v.x)
end

function get_perpendicular_to_vectors(v, w)
    return cross(v, w)
end

function get_wall_anchors(left_pillar, right_pillar)
    left_pillar_center = get_pillar_center(left_pillar)
    left_pillar_width = get_pillar_width(left_pillar)
    right_pillar_center = get_pillar_center(right_pillar)
    right_pillar_width = get_pillar_width(right_pillar)

    left_to_right_pillar_vector = right_pillar_center - left_pillar_center

    left_anchor = displace_point_by_vector(left_pillar_center, left_to_right_pillar_vector, left_pillar_width / 2)
    right_anchor = displace_point_by_vector(right_pillar_center, left_to_right_pillar_vector, -right_pillar_width / 2)

    return (left_anchor = left_anchor, right_anchor = right_anchor)
end

function get_arch_top_height(left_point, right_point, excess)
    arch_midpoint = intermediate_loc(left_point, right_point)

    if isapprox(excess, 0.5)
        return distance(left_point, arch_midpoint)
    else
        excess_displacement_vector = right_point - left_point
        arcs_radius = norm(excess_displacement_vector) * excess
        left_arc_center = displace_point_by_vector(left_point, excess_displacement_vector, arcs_radius)
        height = sqrt(arcs_radius^2 - distance(left_arc_center, arch_midpoint)^2)

        return height
    end

    return nothing
end
# == UTILITIES == #

# == MODELERS == #
function flying_buttress(left_pillar, right_pillar, left_anchor_height, right_anchor_height;
                            depth = get_pillar_depth(left_pillar))
    half_depth = depth / 2
    flying_buttress_anchors = get_wall_anchors(left_pillar, right_pillar)
    left_anchor = flying_buttress_anchors.left_anchor + GROWING_HEIGHT_DIRECTION * left_anchor_height
    right_anchor = flying_buttress_anchors.right_anchor + GROWING_HEIGHT_DIRECTION * right_anchor_height
    flying_buttress_medial_line = line(left_anchor, right_anchor)

    surface_top_half = extrusion(flying_buttress_medial_line, GROWING_HEIGHT_DIRECTION * depth * 4)
    surface_bottom_half = extrusion(flying_buttress_medial_line, DECREASING_HEIGHT_DIRECTION * depth * 4)
    surface_flying_buttress = union(surface_top_half, surface_bottom_half)

    front_depth_direction = EAST
    back_depth_direction = WEST
    surface_depth_front_half = extrusion(surface_flying_buttress, front_depth_direction * half_depth)
    surface_depth_back_half = extrusion(surface_flying_buttress, back_depth_direction * half_depth)
    first_flying_buttress = union(surface_depth_front_half, surface_depth_back_half)
    second_flying_buttress = move(deepcopy(first_flying_buttress), GROWING_HEIGHT_DIRECTION * depth * 16)
end

function standing_lancet_arch_top_block(left_point, right_point, depth, excess)
    arch_depth_vector = GROWING_HEIGHT_DIRECTION * depth / 2
    arch_midpoint = intermediate_loc(left_point, right_point)
    rotation_axis = right_point - arch_midpoint

    if isapprox(excess, 0.5)
        right_arc_center = left_arc_center = arch_midpoint
        arcs = arc(arch_midpoint, right_point, left_point)
        delete_shape(arcs)

        lancet_arch_surface = surface_arc(arcs)
        lancet_arch_depth_first_half = extrusion(lancet_arch_surface, arch_depth_vector)
        lancet_arch_depth_second_half = extrusion(lancet_arch_surface, -arch_depth_vector)

        standing_lancet_arch = union(lancet_arch_depth_first_half, lancet_arch_depth_second_half)

        return rotate(standing_lancet_arch, π/2, arch_midpoint, rotation_axis)
    else
        excess_displacement_vector = right_point - left_point
        arcs_radius = norm(excess_displacement_vector) * excess

        left_arc_center = displace_point_by_vector(left_point, excess_displacement_vector, arcs_radius)
        right_arc_center = displace_point_by_vector(right_point, excess_displacement_vector, -arcs_radius)

        arc_intersection_displacement_vector = get_perpendicular_to_vectors(GROWING_HEIGHT_DIRECTION, excess_displacement_vector)
        arc_intersection_displacement_magnitude = sqrt(arcs_radius^2 - distance(left_arc_center, arch_midpoint)^2)
        arc_intersection = displace_point_by_vector(arch_midpoint, arc_intersection_displacement_vector, arc_intersection_displacement_magnitude)

        right_arc = arc(right_arc_center, right_point, arc_intersection)
        left_arc = arc(left_arc_center, arc_intersection, left_point)
        delete_shape(right_arc)
        delete_shape(left_arc)

        right_arc = surface_arc(right_arc)
        left_arc = surface_arc(left_arc)
        triangle_between_arcs = surface_polygon(right_point, arc_intersection, left_point)

        lancet_arch_surface = union(right_arc, left_arc, triangle_between_arcs)
        lancet_arch_depth_first_half = extrusion(lancet_arch_surface, arch_depth_vector)
        lancet_arch_depth_second_half = extrusion(lancet_arch_surface, -arch_depth_vector)

        standing_lancet_arch = union(lancet_arch_depth_first_half, lancet_arch_depth_second_half)

        return rotate(standing_lancet_arch, π/2, arch_midpoint, rotation_axis)
    end

    return nothing
end

function standing_wall_block(left_anchor, right_anchor, depth, height)
    wall_base_medial_line = line(left_anchor, right_anchor)

    left_to_right_anchor_vector = right_anchor - left_anchor
    wall_depth_vector = normalize_vector(get_perpendicular_to_vectors(left_to_right_anchor_vector, GROWING_HEIGHT_DIRECTION)) * depth / 2

    wall_base_first_half = extrusion(wall_base_medial_line, wall_depth_vector)
    wall_base_second_half = extrusion(wall_base_medial_line, -wall_depth_vector)
    wall_base = union(wall_base_first_half, wall_base_second_half)

    wall = extrusion(wall_base, GROWING_HEIGHT_DIRECTION * height)

    return wall
end

function standing_hollow_wall_block(left_anchor, right_anchor, depth, height, hole_height, offset)
    wall = standing_wall_block(left_anchor, right_anchor, depth, height)

    left_to_right_anchor_vector = right_anchor - left_anchor
    left_anchor_offset = displace_point_by_vector(left_anchor, left_to_right_anchor_vector, offset)
    right_anchor_offset = displace_point_by_vector(right_anchor, left_to_right_anchor_vector, -offset)

    hole = standing_wall_block(left_anchor_offset, right_anchor_offset, depth, hole_height)

    return subtraction(wall, hole)
end

function standing_arch_wall_block(left_pillar, right_pillar, depth, height, excess, offset)
    wall_anchors = get_wall_anchors(left_pillar, right_pillar)
    left_anchor = wall_anchors.left_anchor
    right_anchor = wall_anchors.right_anchor

    left_to_right_anchor_vector = right_anchor - left_anchor
    left_anchor_offset = displace_point_by_vector(left_anchor, left_to_right_anchor_vector, offset)
    right_anchor_offset = displace_point_by_vector(right_anchor, left_to_right_anchor_vector, -offset)

    arch_top_height = get_arch_top_height(left_anchor_offset, right_anchor_offset, excess)
    wall_hole_height = height - arch_top_height - distance(get_pillar_center(left_pillar), left_anchor) * 2

    hollow_wall = standing_hollow_wall_block(left_anchor, right_anchor, depth, height, wall_hole_height, offset)
    arch_top = standing_lancet_arch_top_block(left_anchor_offset + GROWING_HEIGHT_DIRECTION * wall_hole_height, 
                                                right_anchor_offset + GROWING_HEIGHT_DIRECTION * wall_hole_height, 
                                                    depth, excess)
    arch_wall = subtraction(hollow_wall, arch_top)

    return arch_wall
end

function standing_window_wall_block(left_pillar, right_pillar, orientation, height, depth, offset)
    wall_anchors = get_wall_anchors(left_pillar, right_pillar)
    left_anchor = wall_anchors.left_anchor
    right_anchor = wall_anchors.right_anchor
    wall_center = intermediate_loc(left_anchor, right_anchor)
    width = distance(left_anchor, right_anchor)

    wall = standing_wall_block(center, width, depth, height, orientation)

    return wall
end

function main_hallway_arch(left_pillar, right_pillar;
                                depth = get_pillar_depth(left_pillar),
                                height = MAIN_HALLWAYS_WALL_HEIGHT,
                                excess = 1,
                                offset = 0)
    standing_arch_wall_block(left_pillar, right_pillar, depth, height, excess, offset)
end

function aisle_outer_arch(left_pillar, right_pillar;
                                depth = get_pillar_depth(left_pillar),
                                height = AISLE_OUTER_WALL_HEIGHT,
                                excess = 1,
                                offset = 0)
    standing_arch_wall_block(left_pillar, right_pillar, depth, height, excess, offset)
end

function aisle_middle_arch(left_pillar, right_pillar;
                                depth = get_pillar_depth(left_pillar),
                                height = AISLE_MIDDLE_WALL_HEIGHT,
                                excess = 1,
                                offset = 0)
    standing_arch_wall_block(left_pillar, right_pillar, depth, height, excess, offset)
end

function aisle_inner_arch(left_pillar, right_pillar;
                                depth = get_pillar_depth(left_pillar),
                                height = AISLE_INNER_WALL_HEIGHT,
                                excess = 1,
                                offset = 0)
    standing_arch_wall_block(left_pillar, right_pillar, depth, height, excess, offset)
end

function transept_outer_arch(left_pillar, right_pillar;
                                depth = get_pillar_depth(left_pillar),
                                height = TRANSEPT_OUTER_WALL_HEIGHT,
                                excess = 1,
                                offset = 0)
    standing_arch_wall_block(left_pillar, right_pillar, depth, height, excess, offset)
end

function arch_and_window_wall(left_pillar, right_pillar, orientation;
                                height = ARCH_AND_WINDOW_WALL_HEIGHT,
                                depth = get_pillar_depth(left_pillar),
                                offset = 0)
    standing_window_wall_block(left_pillar, right_pillar, orientation, height, depth, offset)
end

#function ambulatory_outer_arch(left_pillar, right_pillar;
#                                height = AMBULATORY_OUTER_WALL_HEIGHT,
#                                depth = AISLE_OUTER_PILLAR_DEPTH,
#                                hole_height = height * 0.85,
#                                offset = 0)
#    # Calculate orientation based upon the right pillar's orientation and decrease by half the value of angle increment centered on the ambulatory center
#    # with the same radius value as that of the current row
#    # better alternativeto solve wall orientation issues is to find the perpendicular vector to the line connecting both wall anchors - its 
#    standing_arch_wall_block(left_pillar, right_pillar, orientation, height, depth, hole_height, offset)
#end
# Wall types:
# bottom arch
# complete arch DONE
# the 2 above can be combined into 1 function (height, depth, excess)
# intermediate block
# complete block
# the 2 above can be combined into 1 function (height and depth)
# bottom window
# top window
# the 2 above can be combined into 1 function
# trapezoid block (outer to ambulatory buttress connection)
# == MODELERS == #

# == INSTANTIATORS == #
function set_pillars_wall_attributes(left_pillar, right_pillar, wall)
    left_pillar_center = get_pillar_center(left_pillar)
    right_pillar_center = get_pillar_center(right_pillar)

    direction = right_pillar_center - left_pillar_center
    normalized_direction = direction / norm(direction)

    if isapprox(normalized_direction.x, SOUTH.x) && isapprox(normalized_direction.y, SOUTH.y)
        set_pillar_south_wall(left_pillar, wall)
        set_pillar_north_wall(right_pillar, wall)
    elseif isapprox(normalized_direction.x, EAST.x) && isapprox(normalized_direction.y, EAST.y)
        set_pillar_east_wall(left_pillar, wall)
        set_pillar_west_wall(right_pillar, wall)
    end
end

function column_range_instantiator(pillars, wall_type_instantiator, column_range, left_pillars_row_index, right_pillars_row_index)
    for column in column_range
        left_pillar = pillars[left_pillars_row_index, column]
        right_pillar = pillars[right_pillars_row_index, column]
        model = wall_type_instantiator(left_pillar, right_pillar)
        wall_type = Wall(model, left_pillar, right_pillar)
        set_pillars_wall_attributes(left_pillar, right_pillar, wall_type)
    end
end

function right_pillar_increment_column_range_instantiator(pillars, wall_type_instantiator, column_range, left_pillars_row_index, right_pillars_row_index)
    for column in column_range
        left_pillar = pillars[left_pillars_row_index, column]
        right_pillar = pillars[right_pillars_row_index, column + 1]
        model = wall_type_instantiator(left_pillar, right_pillar)
        wall_type = Wall(model, left_pillar, right_pillar)
        set_pillars_wall_attributes(left_pillar, right_pillar, wall_type)
    end
end

function flying_buttress_column_range_instantiator(pillars, column_range, left_pillars_row_index, right_pillars_row_index, left_anchor_height, right_anchor_height)
    for column in column_range
        left_pillar = pillars[left_pillars_row_index, column]
        right_pillar = pillars[right_pillars_row_index, column]
        flying_buttress(left_pillar, right_pillar, left_anchor_height, right_anchor_height)
    end
end

function row_range_instantiator(pillars, wall_type_instantiator, row_range, left_pillars_column_index, right_pillars_column_index)
    for row in row_range
        left_pillar = pillars[row, left_pillars_column_index]
        right_pillar = pillars[row, right_pillars_column_index]
        model = wall_type_instantiator(left_pillar, right_pillar)
        wall_type = Wall(model, left_pillar, right_pillar)
        set_pillars_wall_attributes(left_pillar, right_pillar, wall_type)
    end
end

function instantiate_flying_buttresses(pillars)
    flying_buttress_column_range_instantiator(pillars, NW_FLYING_BUTTRESS_LEFT_PILLARS[COLUMN_INDEX],
                                                NW_FLYING_BUTTRESS_LEFT_PILLARS[ROW_INDEX], NW_FLYING_BUTTRESS_RIGHT_PILLARS[ROW_INDEX],
                                                N_FLYING_BUTTRESSES_LEFT_ANCHOR_HEIGHT, N_FLYING_BUTTRESSES_RIGHT_ANCHOR_HEIGHT)
    flying_buttress_column_range_instantiator(pillars, NE_FLYING_BUTTRESS_LEFT_PILLARS[COLUMN_INDEX],
                                                NE_FLYING_BUTTRESS_LEFT_PILLARS[ROW_INDEX], NE_FLYING_BUTTRESS_RIGHT_PILLARS[ROW_INDEX],
                                                N_FLYING_BUTTRESSES_LEFT_ANCHOR_HEIGHT, N_FLYING_BUTTRESSES_RIGHT_ANCHOR_HEIGHT)

    flying_buttress_column_range_instantiator(pillars, SW_FLYING_BUTTRESS_LEFT_PILLARS[COLUMN_INDEX],
                                                SW_FLYING_BUTTRESS_LEFT_PILLARS[ROW_INDEX], SW_FLYING_BUTTRESS_RIGHT_PILLARS[ROW_INDEX],
                                                S_FLYING_BUTTRESSES_LEFT_ANCHOR_HEIGHT, S_FLYING_BUTTRESSES_RIGHT_ANCHOR_HEIGHT)
    flying_buttress_column_range_instantiator(pillars, SE_FLYING_BUTTRESS_LEFT_PILLARS[COLUMN_INDEX],
                                                SE_FLYING_BUTTRESS_LEFT_PILLARS[ROW_INDEX], SE_FLYING_BUTTRESS_RIGHT_PILLARS[ROW_INDEX],
                                                S_FLYING_BUTTRESSES_LEFT_ANCHOR_HEIGHT, S_FLYING_BUTTRESSES_RIGHT_ANCHOR_HEIGHT)
end

function instantiate_main_hallways_walls(pillars)
    column_range_instantiator(pillars, main_hallway_arch, HORIZONTAL_HALLWAY_WALL_LEFT_PILLARS[COLUMN_INDEX], 
                                HORIZONTAL_HALLWAY_WALL_LEFT_PILLARS[ROW_INDEX], HORIZONTAL_HALLWAY_WALL_RIGHT_PILLARS[ROW_INDEX])
    row_range_instantiator(pillars, main_hallway_arch, VERTICAL_HALLWAY_WALL_LEFT_PILLARS[ROW_INDEX], 
                                VERTICAL_HALLWAY_WALL_LEFT_PILLARS[COLUMN_INDEX], VERTICAL_HALLWAY_WALL_RIGHT_PILLARS[COLUMN_INDEX])
end

function instantiate_aisles_outer_walls(pillars)
    column_range_instantiator(pillars, aisle_outer_arch, NW_AISLE_OUTER_WALL_LEFT_PILLARS[COLUMN_INDEX],
                                NW_AISLE_OUTER_WALL_LEFT_PILLARS[ROW_INDEX], NW_AISLE_OUTER_WALL_RIGHT_PILLARS[ROW_INDEX])
    column_range_instantiator(pillars, aisle_outer_arch, NE_AISLE_OUTER_WALL_LEFT_PILLARS[COLUMN_INDEX],
                                NE_AISLE_OUTER_WALL_LEFT_PILLARS[ROW_INDEX], NE_AISLE_OUTER_WALL_RIGHT_PILLARS[ROW_INDEX])

    column_range_instantiator(pillars, aisle_outer_arch, SW_AISLE_OUTER_WALL_LEFT_PILLARS[COLUMN_INDEX],
                                SW_AISLE_OUTER_WALL_LEFT_PILLARS[ROW_INDEX], SW_AISLE_OUTER_WALL_RIGHT_PILLARS[ROW_INDEX])
    column_range_instantiator(pillars, aisle_outer_arch, SE_AISLE_OUTER_WALL_LEFT_PILLARS[COLUMN_INDEX],
                                SE_AISLE_OUTER_WALL_LEFT_PILLARS[ROW_INDEX], SE_AISLE_OUTER_WALL_RIGHT_PILLARS[ROW_INDEX])
end

function instantiate_aisles_middle_walls(pillars)
    right_pillar_increment_column_range_instantiator(pillars, aisle_middle_arch, NW_AISLE_MIDDLE_WALL_LEFT_PILLARS[COLUMN_INDEX],
                                                        NW_AISLE_MIDDLE_WALL_LEFT_PILLARS[ROW_INDEX], NW_AISLE_MIDDLE_WALL_RIGHT_PILLARS[ROW_INDEX])
    right_pillar_increment_column_range_instantiator(pillars, aisle_middle_arch, NE_AISLE_MIDDLE_WALL_LEFT_PILLARS[COLUMN_INDEX],
                                                        NE_AISLE_MIDDLE_WALL_LEFT_PILLARS[ROW_INDEX], NE_AISLE_MIDDLE_WALL_RIGHT_PILLARS[ROW_INDEX])

    right_pillar_increment_column_range_instantiator(pillars, aisle_middle_arch, SW_AISLE_MIDDLE_WALL_LEFT_PILLARS[COLUMN_INDEX],
                                                        SW_AISLE_MIDDLE_WALL_LEFT_PILLARS[ROW_INDEX], SW_AISLE_MIDDLE_WALL_RIGHT_PILLARS[ROW_INDEX])
    right_pillar_increment_column_range_instantiator(pillars, aisle_middle_arch, SE_AISLE_MIDDLE_WALL_LEFT_PILLARS[COLUMN_INDEX],
                                                        SE_AISLE_MIDDLE_WALL_LEFT_PILLARS[ROW_INDEX], SE_AISLE_MIDDLE_WALL_RIGHT_PILLARS[ROW_INDEX])  
end

function instantiate_aisles_inner_walls(pillars)
    column_range_instantiator(pillars, aisle_inner_arch, NW_AISLE_INNER_WALL_LEFT_PILLARS[COLUMN_INDEX],
                                NW_AISLE_INNER_WALL_LEFT_PILLARS[ROW_INDEX], NW_AISLE_INNER_WALL_RIGHT_PILLARS[ROW_INDEX])
    column_range_instantiator(pillars, aisle_inner_arch, NE_AISLE_INNER_WALL_LEFT_PILLARS[COLUMN_INDEX],
                                NE_AISLE_INNER_WALL_LEFT_PILLARS[ROW_INDEX], NE_AISLE_INNER_WALL_RIGHT_PILLARS[ROW_INDEX])

    column_range_instantiator(pillars, aisle_inner_arch, SW_AISLE_INNER_WALL_LEFT_PILLARS[COLUMN_INDEX],
                                SW_AISLE_INNER_WALL_LEFT_PILLARS[ROW_INDEX], SW_AISLE_INNER_WALL_RIGHT_PILLARS[ROW_INDEX])
    column_range_instantiator(pillars, aisle_inner_arch, SE_AISLE_INNER_WALL_LEFT_PILLARS[COLUMN_INDEX],
                                SE_AISLE_INNER_WALL_LEFT_PILLARS[ROW_INDEX], SE_AISLE_INNER_WALL_RIGHT_PILLARS[ROW_INDEX])
end

function instantiate_transept_outer_walls(pillars)
    row_range_instantiator(pillars, transept_outer_arch, WN_TRANSEPT_OUTER_WALL_LEFT_PILLARS[ROW_INDEX],
                                WN_TRANSEPT_OUTER_WALL_LEFT_PILLARS[COLUMN_INDEX], WN_TRANSEPT_OUTER_WALL_RIGHT_PILLARS[COLUMN_INDEX])
    row_range_instantiator(pillars, transept_outer_arch, EN_TRANSEPT_OUTER_WALL_LEFT_PILLARS[ROW_INDEX],
                                EN_TRANSEPT_OUTER_WALL_LEFT_PILLARS[COLUMN_INDEX], EN_TRANSEPT_OUTER_WALL_RIGHT_PILLARS[COLUMN_INDEX])

    row_range_instantiator(pillars, transept_outer_arch, WS_TRANSEPT_OUTER_WALL_LEFT_PILLARS[ROW_INDEX],
                                WS_TRANSEPT_OUTER_WALL_LEFT_PILLARS[COLUMN_INDEX], WS_TRANSEPT_OUTER_WALL_RIGHT_PILLARS[COLUMN_INDEX])
    row_range_instantiator(pillars, transept_outer_arch, ES_TRANSEPT_OUTER_WALL_LEFT_PILLARS[ROW_INDEX],
                                ES_TRANSEPT_OUTER_WALL_LEFT_PILLARS[COLUMN_INDEX], ES_TRANSEPT_OUTER_WALL_RIGHT_PILLARS[COLUMN_INDEX])
end

function instantiate_all_walls(pillars)
    instantiate_flying_buttresses(pillars)
    instantiate_main_hallways_walls(pillars)
    instantiate_aisles_outer_walls(pillars)
    instantiate_aisles_middle_walls(pillars)
    instantiate_aisles_inner_walls(pillars)
    instantiate_transept_outer_walls(pillars)
end
# == INSTANTIATORS == #
# == WALLS == #

# == PLAYGROUND == #
pillars = Array{Union{Pillar, Nothing}}(nothing, 10, 17)

instantiate_all_pillars(pillars)
instantiate_all_walls(pillars)

#pillar_test = pillars[1, 9]
#pillar_test2 = pillars[5, 8]
#pillar_test3 = pillars[5, 9]
#
#delete_shape(get_wall_model(get_pillar_east_wall(pillar_test)))
#delete_shape(get_wall_model(get_pillar_south_wall(pillar_test2)))
#delete_shape(get_wall_model(get_pillar_south_wall(pillar_test3)))
#delete_shape(get_wall_model(get_pillar_east_wall(pillar_test3)))
# == PLAYGROUND == #
