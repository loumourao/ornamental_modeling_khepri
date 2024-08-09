module GothicStructuralGeometry
# == ARCHES == #
function get_offset_excess(left_point, right_point, previous_excess, offset_value)
    width_vector = right_point - left_point
    width = norm(width_vector)
    left_arc_center = displace_point_by_vector(left_point, width_vector, width * previous_excess)
    
    offset_left_point = left_point + vx(offset_value)
    offset_right_point = right_point - vx(offset_value)
    offset_width_vector = offset_right_point - offset_left_point
    offset_width = norm(offset_width_vector)
    offset_radius = distance(left_point, left_arc_center) - offset_value

    return offset_radius / offset_width
end

function get_left_sub_arch_body(bottom_left_corner, upper_right_corner, outer_offset, inner_offset, vertical_distance_to_sub_arches)
    bottom_right_corner = xy(upper_right_corner.x, bottom_left_corner.y)
    arch_bottom_midpoint = intermediate_loc(bottom_left_corner, bottom_right_corner)

    bottom_left_corner = bottom_left_corner + vxy(outer_offset, outer_offset)
    bottom_right_corner = arch_bottom_midpoint - vx(inner_offset / 2) + vy(outer_offset)
    upper_left_corner = xy(bottom_left_corner.x, upper_right_corner.y - vertical_distance_to_sub_arches)
    upper_right_corner = xy(bottom_right_corner.x, upper_left_corner.y)

    return (bottom_left_corner = bottom_left_corner, bottom_right_corner = bottom_right_corner, 
                upper_left_corner = upper_left_corner, upper_right_corner = upper_right_corner)
end

function get_right_sub_arch_body(bottom_left_corner, upper_right_corner, outer_offset, inner_offset, vertical_distance_to_sub_arches)
    upper_left_corner = xy(bottom_left_corner.x, upper_right_corner.y)
    arch_upper_midpoint = intermediate_loc(upper_left_corner, upper_right_corner)

    upper_right_corner = upper_right_corner - vxy(outer_offset, vertical_distance_to_sub_arches)
    upper_left_corner = arch_upper_midpoint + vx(inner_offset / 2) - vy(vertical_distance_to_sub_arches)
    bottom_right_corner = xy(upper_right_corner.x, bottom_left_corner.y + outer_offset)
    bottom_left_corner = xy(upper_left_corner.x, bottom_right_corner.y)

    return (bottom_left_corner = bottom_left_corner, bottom_right_corner = bottom_right_corner, 
                upper_left_corner = upper_left_corner, upper_right_corner = upper_right_corner)
end

export get_offset_excess, get_left_sub_arch_body, get_right_sub_arch_body
# == ARCH TOPS == #
function lancet_arch_top(left_point, right_point, excess)
    arch_midpoint = intermediate_loc(left_point, right_point)
    
    if isapprox(excess, 0.5)
        right_arc_center = left_arc_center = arch_midpoint
        arc(arch_midpoint, right_point, left_point)
    else
        excess_displacement_vector = right_point - left_point
        arcs_radius = norm(excess_displacement_vector) * excess

        right_arc_center = displace_point_by_vector(right_point, excess_displacement_vector, -arcs_radius)
        left_arc_center = displace_point_by_vector(left_point, excess_displacement_vector, arcs_radius)
        
        arc_intersection_x = arch_midpoint.x
        arc_intersection_y = sqrt(arcs_radius^2 - distance(left_arc_center, arch_midpoint)^2) + arch_midpoint.y
        arc_intersection = xy(arc_intersection_x, arc_intersection_y)

        right_arc = arc(right_arc_center, right_point, arc_intersection)
        left_arc = arc(left_arc_center, arc_intersection, left_point)
    end

    return (right_arc = right_arc, left_arc = left_arc)
end

export lancet_arch_top
# == ARCH TOPS == #
# == ARCHES == #

# == ROSETTE == #
function compute_rosette(bottom_left_corner, upper_right_corner, 
                            main_arch_right_arc, left_sub_arch_right_arc, 
                                outer_offset, inner_offset)
    main_arch_right_arc_center = arc_center(main_arch_right_arc)
    main_arch_right_arc_radius = arc_radius(main_arch_right_arc)
    left_sub_arch_right_arc_center = arc_center(left_sub_arch_right_arc)
    left_sub_arch_right_arc_radius = arc_radius(left_sub_arch_right_arc)

    ellipse_center = intermediate_loc(left_sub_arch_right_arc_center, main_arch_right_arc_center)
    ellipse_radius = (main_arch_right_arc_radius - outer_offset + left_sub_arch_right_arc_radius) / 2

    vertical_axis = vy(1)
    vertical_axis_point = intermediate_loc(bottom_left_corner, upper_right_corner)

    rosette_center = intersect_line_ellipse(vertical_axis_point, vertical_axis, ellipse_center, ellipse_radius).second_intersection_point
    rosette_radius = distance(rosette_center, left_sub_arch_right_arc_center) - left_sub_arch_right_arc_radius - inner_offset

    return (rosette_center = rosette_center, rosette_radius = rosette_radius)
end

export compute_rosette
# == ROSETTE == #
end
