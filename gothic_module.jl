using KhepriAutoCAD
const K = KhepriAutoCAD

delete_all_shapes()

# == UTILITY FUNCTIONS == #
# == MATH == #
function solve_quadratic(a, b, c)
    discriminant = b^2 - 4 * a * c
    
    if discriminant < 0
        return nothing
    elseif discriminant == 0
        return [-b / (2 * a)]
    else
        sqrt_discriminant = sqrt(discriminant)
        return [(-b - sqrt_discriminant) / (2 * a), (-b + sqrt_discriminant) / (2 * a)]
    end
end
# == MATH == #

# == VECTORS == #
function displace_point_by_vector(point, vector, magnitude)
    unit_vector = vector / norm(vector)

    return point + unit_vector * magnitude
end
# == VECTORS == #

# == CIRCLES == #
function intersect_circles(m0, r0, m1, r1)
    d = distance(m0, m1)
    
    # d == 0 prevents the set of interesection points when both circles overlap
    if d > r0 + r1 || d < abs(r0 - r1) || d == 0
        return nothing
    end

    a = (r0^2 - r1^2 + d^2) / (2 * d)
    p = m0 + a * (m1 - m0) / d
    c = sqrt(r0^2 - a^2)
    c_vector = m1 - m0
    c_vector = vxy(-c_vector.y, c_vector.x)
    c_vector = c_vector / norm(c_vector) * c

    first_intersection_point = p + c_vector
    second_intersection_point = p - c_vector

    return first_intersection_point.y >= second_intersection_point.y ? 
            [first_intersection_point, second_intersection_point] : [second_intersection_point, first_intersection_point]
end

function intersect_line_ellipse(p0, v, m, r)
    a = norm(v)^2
    b = 2 * dot((p0 - m), v)
    c = norm(p0 - m)^2 - r^2

    t_values = solve_quadratic(a, b, c)
    
    if length(t_values) == 0
        return nothing
    end

    t_values = sort(t_values)
    t0 = t_values[1]
    t1 = t_values[2]

    q0 = p0 + v * t0
    q1 = p0 + v * t1

    return [q0, q1, t0, t1]
end

# == ARCS == #
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

    return K.arc(center, radius, ccw_start_angle, amplitude)
end

function arc_start_point(arc)
    center = arc_center(arc)
    radius = arc_radius(arc)
    start_angle = arc_start_angle(arc)

    return center + vpol(radius, start_angle)
end

function arc_end_point(arc)
    center = arc_center(arc)
    radius = arc_radius(arc)
    start_angle = arc_start_angle(arc)
    amplitude = arc_amplitude(arc)

    return center + vpol(radius, start_angle + amplitude)
end

function intersect_arcs(arc0, arc1)
    m0 = arc_center(arc0)
    r0 = arc_radius(arc0)
    arc0_start_angle = arc_start_angle(arc0)
    arc0_amplitude = arc_amplitude(arc0)

    m1 = arc_center(arc1)
    r1 = arc_radius(arc1)
    arc1_start_angle = arc_start_angle(arc1)
    arc1_amplitude = arc_amplitude(arc1)

    circle_intersection_points = intersect_circles(m0, r0, m1, r1)

    if isnothing(circle_intersection_points)
        return nothing
    end

    valid_intersection_points = []

    for intersection_point in circle_intersection_points
        a = intersection_point - m0
        a_θ = atan(a.y, a.x)
        a_θ_amplitude = arc0_start_angle > a_θ ? 2π - arc0_start_angle + a_θ : a_θ - arc0_start_angle

        b = intersection_point - m1
        b_θ = atan(b.y, b.x)
        b_θ_amplitude = arc1_start_angle > b_θ ? 2π - arc1_start_angle + b_θ : b_θ - arc1_start_angle

        if a_θ_amplitude <= arc0_amplitude && b_θ_amplitude <= arc1_amplitude
            push!(valid_intersection_points, intersection_point)
        end
    end

    if isempty(valid_intersection_points)
        return nothing
    end

    return length(valid_intersection_points) > 1 ? valid_intersection_points : valid_intersection_points[1]
end

function offset_arc(previous_arc, offset_value)
    center = arc_center(previous_arc)
    radius = arc_radius(previous_arc)
    start_point = arc_start_point(previous_arc)
    end_point = arc_end_point(previous_arc)

    ca = start_point - center
    cb = end_point - center
    new_radius = radius - offset_value
    new_start_point = displace_point_by_vector(center, ca, new_radius)
    new_end_point = displace_point_by_vector(center, cb, new_radius)

    return arc(center, new_start_point, new_end_point)
end

function arc_bidirectionally_extended_uniform_offset(arc, offset_value)
    delete_shape(arc)
    center = arc_center(arc)
    radius = arc_radius(arc)
    start_angle = arc_start_angle(arc)
    amplitude = arc_amplitude(arc)

    scaling_factor = (radius - offset_value) / radius
    scaled_amplitude = amplitude * scaling_factor
    bidirectional_extension = amplitude - scaled_amplitude

    new_radius = radius - offset_value
    new_amplitude = amplitude + bidirectional_extension
    new_start_angle = start_angle - bidirectional_extension / 2

    return K.arc(center, new_radius, new_start_angle, new_amplitude)
end

function arc_amplitude_extended_offset(arc, offset_value)
    delete_shape(arc)
    center = arc_center(arc)
    radius = arc_radius(arc)
    start_angle = arc_start_angle(arc)
    amplitude = arc_amplitude(arc)

    scaling_factor = (radius - offset_value) / radius
    scaled_amplitude = amplitude * scaling_factor
    amplitude_extension = (amplitude - scaled_amplitude) / 2

    new_radius = radius - offset_value
    new_start_angle = start_angle + amplitude_extension

    return K.arc(center, new_radius, new_start_angle, amplitude)
end

function arc_start_angle_extended_offset(arc, offset_value)
    delete_shape(arc)
    center = arc_center(arc)
    radius = arc_radius(arc)
    start_angle = arc_start_angle(arc)
    amplitude = arc_amplitude(arc)

    scaling_factor = (radius - offset_value) / radius
    scaled_amplitude = amplitude * scaling_factor
    amplitude_extension = (amplitude - scaled_amplitude) / 2

    new_radius = radius - offset_value
    new_start_angle = start_angle - amplitude_extension

    return K.arc(center, new_radius, new_start_angle, amplitude)
end
# == ARCS == #
# == CIRCLES == #

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
# == ARCHES == #
# == UTILITY FUNCTIONS == #

# == ORNAMENTATION FUNCTIONS == #
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
# == ARCH TOPS == #

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

    rosette_center = intersect_line_ellipse(vertical_axis_point, vertical_axis, ellipse_center, ellipse_radius)[2]
    rosette_radius = distance(rosette_center, left_sub_arch_right_arc_center) - left_sub_arch_right_arc_radius - inner_offset

    return (rosette_center = rosette_center, rosette_radius = rosette_radius)
end

function rosette_rounded_foils(center, radius, n_foils, orientation, inner_offset)
    # Check if n_foils >= 1, otherwise raise an error
    Δα = 2π / n_foils
    
    foil_radius = (radius * sin(Δα/2)) / (1 + sin(Δα/2))
    center_to_foil_center_length = radius - foil_radius
    center_to_foil_end_points = center_to_foil_center_length * cos(Δα/2)

    foil_center = center + vpol(center_to_foil_center_length, orientation)
    foil_start_point = center + vpol(center_to_foil_end_points, orientation - (Δα/2))
    foil_end_point = center + vpol(center_to_foil_end_points, orientation + (Δα/2))

    
    fillet_points = get_rosette_rounded_foils_fillet_points(center, radius, 
                                                                foil_center, foil_radius, 
                                                                    Δα, inner_offset)
    fillet_right_point = fillet_points.fillet_right_point
    fillet_left_point = fillet_points.fillet_left_point
    fillet_bottom_point = fillet_points.fillet_bottom_point
    fillet_right_arc_center = fillet_points.fillet_right_arc_center
    fillet_left_arc_center = fillet_points.fillet_left_arc_center
    fillet_upper_arc_center = fillet_points.fillet_upper_arc_center
    
    current_rotation_angle = 0

    while n_foils > 0
        # Fillets
        fillet_arc_position = current_rotation_angle + orientation
        rotate(arc(fillet_right_arc_center, fillet_right_point, fillet_bottom_point), fillet_arc_position, center)
        rotate(arc(fillet_upper_arc_center, fillet_right_point, fillet_left_point), fillet_arc_position, center)
        rotate(arc(fillet_left_arc_center, fillet_bottom_point, fillet_left_point), fillet_arc_position, center)

        # Rounded foils
        foil = arc_bidirectionally_extended_uniform_offset(arc(foil_center, foil_start_point, foil_end_point), inner_offset)
        rotate(foil, current_rotation_angle, center)

        # Rounded foils connections
        connection_start_point = arc_end_point(foil)
        center_to_foil_start_point = arc_start_point(foil) - center
        center_to_foil_start_point_polar_angle = atan(center_to_foil_start_point.y, center_to_foil_start_point.x)
        connection_end_point = center + vpol(norm(center_to_foil_start_point), center_to_foil_start_point_polar_angle + Δα)
        connection = line(connection_start_point, connection_end_point)
        rotate(connection, current_rotation_angle, center)

        # 3D rosette rounded foils
        three_dimensionalized_foil_and_connection = three_dimensionalize_rosette_rounded_foils(foil, connection, -inner_offset)
        rotate(three_dimensionalized_foil_and_connection.foil, current_rotation_angle, center)
        #rotate(three_dimensionalized_foil_and_connection.connection, current_rotation_angle, center)

        current_rotation_angle += Δα
        n_foils -= 1
    end
end

function rosette_pointed_foils(center, radius, n_foils, displacement_ratio, orientation, inner_offset)
    # Check if n_foils >= 1 and if displacement_ratio > 1, otherwise raise an error
    Δα = 2π / n_foils

    rounded_foil_radius = (radius * sin(Δα/2)) / (1 + sin(Δα/2))
    center_to_foil_center_length = radius - rounded_foil_radius
    center_to_foil_end_points = center_to_foil_center_length * cos(Δα/2)

    rounded_foil_center = center + vpol(center_to_foil_center_length, orientation)
    foil_start_point = center + vpol(center_to_foil_end_points, orientation - (Δα/2))
    foil_end_point = center + vpol(center_to_foil_end_points, orientation + (Δα/2))

    foil_right_arc_center_displacement_vector = rounded_foil_center - foil_start_point
    foil_left_arc_center_displacement_vector = rounded_foil_center - foil_end_point
    displacement_vectors_magnitude = norm(foil_right_arc_center_displacement_vector) * displacement_ratio

    foil_right_arc_center = displace_point_by_vector(foil_start_point, foil_right_arc_center_displacement_vector, displacement_vectors_magnitude)
    foil_left_arc_center = displace_point_by_vector(foil_end_point, foil_left_arc_center_displacement_vector, displacement_vectors_magnitude)

    pointed_foil_radius = distance(foil_right_arc_center, foil_start_point)
    arc_intersection_axis_point = intermediate_loc(foil_right_arc_center, foil_left_arc_center)
    foil_right_arc_center_to_arc_intersection_axis_point_length = distance(foil_right_arc_center, arc_intersection_axis_point)

    rounded_foil_center_to_arc_intersection_vector = arc_intersection_axis_point - rounded_foil_center
    rounded_foil_center_to_arc_intersection_vector_magnitude = sqrt(pointed_foil_radius^2 - foil_right_arc_center_to_arc_intersection_axis_point_length^2)
    arc_intersection = displace_point_by_vector(arc_intersection_axis_point, rounded_foil_center_to_arc_intersection_vector, rounded_foil_center_to_arc_intersection_vector_magnitude)

    scaling_factor = radius / distance(center, arc_intersection)

    fillet_points = get_rosette_pointed_foils_fillet_points(center, radius, 
                                                                rounded_foil_center, foil_right_arc_center, pointed_foil_radius, 
                                                                    Δα, scaling_factor, inner_offset)
    fillet_right_point = fillet_points.fillet_right_point
    fillet_left_point = fillet_points.fillet_left_point
    fillet_bottom_point = fillet_points.fillet_bottom_point
    fillet_right_arc_center = fillet_points.fillet_right_arc_center
    fillet_left_arc_center = fillet_points.fillet_left_arc_center
    fillet_upper_arc_center = fillet_points.fillet_upper_arc_center

    current_rotation_angle = 0

    while n_foils > 0
        # Fillets
        fillet_arc_position = current_rotation_angle + orientation
        scale(rotate(arc(fillet_right_arc_center, fillet_right_point, fillet_bottom_point), fillet_arc_position, center), scaling_factor, center)
        scale(rotate(arc(fillet_upper_arc_center, fillet_right_point, fillet_left_point), fillet_arc_position, center), scaling_factor, center)
        scale(rotate(arc(fillet_left_arc_center, fillet_bottom_point, fillet_left_point), fillet_arc_position, center), scaling_factor, center)

        # Pointed foils
        foil_right_arc = arc_start_angle_extended_offset(arc(foil_right_arc_center, foil_start_point, arc_intersection), inner_offset)
        foil_left_arc = arc_amplitude_extended_offset(arc(foil_left_arc_center, arc_intersection, foil_end_point), inner_offset)
        scale(rotate(foil_left_arc, current_rotation_angle, center), scaling_factor, center)
        scale(rotate(foil_right_arc, current_rotation_angle, center), scaling_factor, center)

        # Pointed foils connections
        connection_start_point = arc_end_point(foil_left_arc)
        center_to_foil_start_point = arc_start_point(foil_right_arc) - center
        center_to_foil_start_point_polar_angle = atan(center_to_foil_start_point.y, center_to_foil_start_point.x)
        connection_end_point = center + vpol(norm(center_to_foil_start_point), center_to_foil_start_point_polar_angle + Δα)
        connection = line(connection_start_point, connection_end_point)
        scale(rotate(connection, current_rotation_angle, center), scaling_factor, center)

        # 3D rosette pointed foils
        three_dimensionalized_foil_and_connection = three_dimensionalize_rosette_pointed_foils(foil_right_arc, foil_left_arc, connection, -inner_offset)
        scale(rotate(three_dimensionalized_foil_and_connection.foil_right_arc, current_rotation_angle, center), scaling_factor, center)
        scale(rotate(three_dimensionalized_foil_and_connection.foil_left_arc, current_rotation_angle, center), scaling_factor, center)
        #scale(rotate(three_dimensionalized_foil_and_connection.connection, current_rotation_angle, center), scaling_factor, center)

        n_foils -= 1
        current_rotation_angle += Δα
    end
end
# == ROSETTE == #

# == FILLETS == #
function circular_rosette_fillets(right_arc_center, left_arc_center, arcs_radius,
                                    rosette_center, rosette_radius, 
                                        right_sub_arch_right_arc_center, right_sub_arch_left_arc_center, 
                                            left_sub_arch_left_arc_center, left_sub_arch_right_arc_center, sub_arcs_radius)
    # Right fillet points calculations
    right_fillet_bottom_point = intersect_circles(right_arc_center, arcs_radius, right_sub_arch_right_arc_center, sub_arcs_radius)[1]
    right_fillet_top_point = intersect_circles(right_arc_center, arcs_radius, rosette_center, rosette_radius)[2]
    right_fillet_left_point = intersect_circles(right_sub_arch_right_arc_center, sub_arcs_radius, rosette_center, rosette_radius)[1]

    # Right fillet modeling
    right_fillet_right_arc = arc(right_arc_center, right_fillet_bottom_point, right_fillet_top_point)
    right_fillet_left_arc = arc(rosette_center, right_fillet_left_point, right_fillet_top_point)
    right_fillet_bottom_arc = arc(right_sub_arch_right_arc_center, right_fillet_bottom_point, right_fillet_left_point)

    # Top fillet points calculations
    top_fillet_right_point = intersect_circles(right_arc_center, arcs_radius, rosette_center, rosette_radius)[1]
    top_fillet_top_point = intersect_circles(left_arc_center, arcs_radius, right_arc_center, arcs_radius)[1]
    top_fillet_left_point = intersect_circles(left_arc_center, arcs_radius, rosette_center, rosette_radius)[1]

    # Top fillet modeling
    arc(right_arc_center, top_fillet_right_point, top_fillet_top_point)
    arc(left_arc_center, top_fillet_top_point, top_fillet_left_point)
    arc(rosette_center, top_fillet_right_point, top_fillet_left_point)

    # Left fillet modeling
    left_fillet_right_arc = deepcopy(right_fillet_left_arc)
    left_fillet_left_arc = deepcopy(right_fillet_right_arc)
    left_fillet_bottom_arc = deepcopy(right_fillet_bottom_arc)

    mirror(left_fillet_right_arc, rosette_center)
    mirror(left_fillet_left_arc, rosette_center)
    mirror(left_fillet_bottom_arc, rosette_center)

    # I left this here in the event of code refactoring due to simplifying measures linked to symmetric analysis of the sub-arches
    # Left fillet points calculations
    #left_fillet_right_point = intersect_circles(left_sub_arch_left_arc_center, sub_arcs_radius, rosette_center, rosette_radius)[1]
    #left_fillet_top_point = intersect_circles(left_arc_center, arcs_radius, rosette_center, rosette_radius)[2]
    #left_fillet_bottom_point = intersect_circles(left_arc_center, arcs_radius, left_sub_arch_left_arc_center, sub_arcs_radius)[1]

    ## Left fillet modeling
    #arc(rosette_center, left_fillet_top_point, left_fillet_right_point)
    #arc(left_arc_center, left_fillet_top_point, left_fillet_bottom_point)
    #arc(left_sub_arch_left_arc_center, left_fillet_right_point, left_fillet_bottom_point)

    # Bottom fillet points calculations
    bottom_fillet_right_point = intersect_circles(right_sub_arch_left_arc_center, sub_arcs_radius, rosette_center, rosette_radius)[2]
    bottom_fillet_left_point = intersect_circles(left_sub_arch_right_arc_center, sub_arcs_radius, rosette_center, rosette_radius)[2]
    bottom_fillet_bottom_point = intersect_circles(left_sub_arch_right_arc_center, sub_arcs_radius, right_sub_arch_left_arc_center, sub_arcs_radius)[1]

    # Bottom fillet modeling
    arc(right_sub_arch_left_arc_center, bottom_fillet_right_point, bottom_fillet_bottom_point)
    arc(rosette_center, bottom_fillet_left_point, bottom_fillet_right_point)
    arc(left_sub_arch_right_arc_center, bottom_fillet_bottom_point, bottom_fillet_left_point)
end

function get_rosette_rounded_foils_fillet_points(rosette_center, rosette_radius, 
                                                    foil_center, foil_radius, 
                                                        Δα, inner_offset)
    rosette_radius = rosette_radius - inner_offset
    rosette_center_to_foil_center_length = distance(rosette_center, foil_center)

    right_foil_center = rosette_center + vpol(rosette_center_to_foil_center_length, 0)
    left_foil_center = rosette_center + vpol(rosette_center_to_foil_center_length, Δα)
    displacement_vector = left_foil_center - right_foil_center

    fillet_right_point = intersect_circles(rosette_center, rosette_radius, right_foil_center, foil_radius)[1]
    fillet_left_point = intersect_circles(rosette_center, rosette_radius, left_foil_center, foil_radius)[2]
    fillet_bottom_point = displace_point_by_vector(right_foil_center, displacement_vector, norm(displacement_vector) / 2)

    return (fillet_right_point = fillet_right_point, fillet_left_point = fillet_left_point, fillet_bottom_point = fillet_bottom_point, 
                fillet_right_arc_center = right_foil_center, fillet_left_arc_center = left_foil_center, fillet_upper_arc_center = rosette_center)
end

function get_rosette_pointed_foils_fillet_points(rosette_center, rosette_radius, 
                                                    foil_center, right_arc_center, foil_radius, 
                                                        Δα, scaling_factor, inner_offset)
    scaling_factor = 1 / scaling_factor
    rosette_radius = rosette_radius * scaling_factor - inner_offset

    center_to_foil_center_vector = foil_center - rosette_center
    center_to_foil_arc_center_vector = right_arc_center - rosette_center
    center_to_foil_center_length = norm(center_to_foil_center_vector)
    center_to_foil_arc_center_length = norm(center_to_foil_arc_center_vector)
    Δβ = angle_between(center_to_foil_center_vector, center_to_foil_arc_center_vector)

    right_foil_center = rosette_center + vpol(center_to_foil_center_length, 0)
    left_foil_center = rosette_center + vpol(center_to_foil_center_length, Δα)
    right_foil_left_arc_center = rosette_center + vpol(center_to_foil_arc_center_length, -Δβ)
    left_foil_right_arc_center = rosette_center + vpol(center_to_foil_arc_center_length, Δα + Δβ)
    displacement_vector = left_foil_center - right_foil_center

    fillet_right_point = intersect_circles(rosette_center, rosette_radius, right_foil_left_arc_center, foil_radius)[1]
    fillet_left_point = intersect_circles(rosette_center, rosette_radius, left_foil_right_arc_center, foil_radius)[2]
    fillet_bottom_point = displace_point_by_vector(right_foil_center, displacement_vector, norm(displacement_vector) / 2)

    return (fillet_right_point = fillet_right_point, fillet_left_point = fillet_left_point, fillet_bottom_point = fillet_bottom_point, 
                fillet_right_arc_center = right_foil_left_arc_center, fillet_left_arc_center = left_foil_right_arc_center, fillet_upper_arc_center = rosette_center)
end
# == FILLETS == #
# == ORNAMENTATION FUNCTIONS == #

function three_dimensionalize_arch_top(left_arc, right_arc, offset_value, profile = surface_circle(u0(), offset_value / 2))
    left_arc = offset_arc(left_arc, offset_value / 2)
    right_arc = offset_arc(right_arc, offset_value / 2)

    sweep(left_arc, profile)
    sweep(right_arc, profile)
end

function three_dimensionalize_arch_body(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner, 
                                                offset_value, profile = surface_circle(x(offset_value / 2), offset_value / 2))
    arch_body = line(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner)
    sweep(arch_body, profile)
end

function three_dimensionalize_arch_middle(bottom_left_corner, upper_right_corner, vertical_distance_to_sub_arch, 
                                                    offset_value, profile = surface_circle(u0(), offset_value / 2))
    bottom_right_corner = xy(upper_right_corner.x, bottom_left_corner.y)
    bottom_midpoint = intermediate_loc(bottom_left_corner, bottom_right_corner)

    height = (upper_right_corner - bottom_right_corner).y
    vertical_displacement = height - vertical_distance_to_sub_arch

    upper_midpoint = bottom_midpoint + vy(vertical_displacement)
    arch_middle = line(upper_midpoint, bottom_midpoint)

    sweep(arch_middle, profile)
end

function three_dimensionalize_rosette(rosette_center, rosette_radius, 
                                        outer_offset, inner_offset, 
                                            profile = surface_circle(u0(), (outer_offset + inner_offset) / 2))
    rosette = circle(rosette_center, rosette_radius - outer_offset + ((outer_offset + inner_offset) / 2))
    sweep(rosette, profile)
end

function three_dimensionalize_rosette_rounded_foils(foil, connection, offset_value, profile = surface_circle(u0(), abs(offset_value) / 2))
    foil = arc_bidirectionally_extended_uniform_offset(foil, offset_value / 2)
    connection = nothing
    
    foil = sweep(foil, profile)

    return (foil = foil, connection = connection)
end

function three_dimensionalize_rosette_pointed_foils(foil_right_arc, foil_left_arc, connection, offset_value, profile = surface_circle(u0(), abs(offset_value) / 2))
    foil_right_arc = arc_start_angle_extended_offset(foil_right_arc, offset_value / 2)
    foil_left_arc = arc_amplitude_extended_offset(foil_left_arc, offset_value / 2)
    connection = nothing

    foil_right_arc = sweep(foil_right_arc, profile)
    foil_left_arc = sweep(foil_left_arc, profile)

    return (foil_right_arc = foil_right_arc, foil_left_arc = foil_left_arc, connection = connection)
end
# == MAIN ARCH == #
function three_dimensional_arch(bottom_left_corner, upper_right_corner, excess, 
                                    recursion_level, vertical_distance_to_sub_arch, 
                                        outer_offset, inner_offset)
    # Arch body auxiliary coordinates
    upper_left_corner = xy(bottom_left_corner.x, upper_right_corner.y)
    bottom_right_corner = xy(upper_right_corner.x, bottom_left_corner.y)

    # Arch body
    #arch_body = line(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner)

    # Arch top portion
    arcs = lancet_arch_top(upper_left_corner, upper_right_corner, excess)
    right_arc_center = arc_center(arcs.right_arc)
    left_arc_center = arc_center(arcs.left_arc)
    arcs_radius = arc_radius(arcs.right_arc)

    # 3D Arch
    three_dimensionalize_arch_top(arcs.left_arc, arcs.right_arc, outer_offset)
    three_dimensionalize_arch_body(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner, outer_offset)
    
    # Sub-Arches
    if recursion_level > 0
        # 3D Arch continuation
        three_dimensionalize_arch_middle(bottom_left_corner, upper_right_corner, vertical_distance_to_sub_arch, inner_offset)

        # Auxiliary calculations
        arch_width = distance(bottom_left_corner, bottom_right_corner)
        outer_offset_ratio = outer_offset / arch_width
        inner_offset_ratio = inner_offset / arch_width
        
        # Sub-Arches auxiliary calculations
        sub_arch_excess = get_offset_excess(upper_left_corner, upper_right_corner, excess, outer_offset)
        left_sub_arch_body = get_left_sub_arch_body(bottom_left_corner, upper_right_corner, outer_offset, inner_offset, vertical_distance_to_sub_arch)
        right_sub_arch_body = get_right_sub_arch_body(bottom_left_corner, upper_right_corner, outer_offset, inner_offset, vertical_distance_to_sub_arch)

        # Sub-arches offset values
        sub_arch_width = distance(left_sub_arch_body.bottom_left_corner, left_sub_arch_body.bottom_right_corner)
        sub_arch_outer_offset = sub_arch_width * outer_offset_ratio
        sub_arch_inner_offset = sub_arch_width * inner_offset_ratio

        # Outer sub-arches calculationss
        outer_sub_arch_excess = get_offset_excess(left_sub_arch_body.upper_left_corner, left_sub_arch_body.upper_right_corner, 
                                                    sub_arch_excess, -inner_offset)

        left_outer_sub_arch_top = lancet_arch_top(left_sub_arch_body.upper_left_corner - vx(inner_offset), left_sub_arch_body.upper_right_corner + vx(inner_offset), 
                                                    outer_sub_arch_excess)

        right_outer_sub_arch_top = lancet_arch_top(right_sub_arch_body.upper_left_corner - vx(inner_offset), right_sub_arch_body.upper_right_corner + vx(inner_offset),
                                                        outer_sub_arch_excess)
        
        # 3D outer arches
        three_dimensionalize_arch_top(left_outer_sub_arch_top.left_arc, left_outer_sub_arch_top.right_arc, inner_offset)
        three_dimensionalize_arch_top(right_outer_sub_arch_top.left_arc, right_outer_sub_arch_top.right_arc, inner_offset)

        # Sub-arches
        left_sub_arch = three_dimensional_arch(left_sub_arch_body.bottom_left_corner, left_sub_arch_body.upper_right_corner, sub_arch_excess, 
                                                    recursion_level - 1, vertical_distance_to_sub_arch, 
                                                        sub_arch_outer_offset, sub_arch_inner_offset)
        right_sub_arch = three_dimensional_arch(right_sub_arch_body.bottom_left_corner, right_sub_arch_body.upper_right_corner, sub_arch_excess, 
                                                    recursion_level - 1, vertical_distance_to_sub_arch, 
                                                        sub_arch_outer_offset, sub_arch_inner_offset)
        left_sub_arch_left_arc_center = arc_center(left_sub_arch.left_arc)
        left_sub_arch_right_arc_center = arc_center(left_sub_arch.right_arc)

        right_sub_arch_left_arc_center = arc_center(right_sub_arch.left_arc)
        right_sub_arch_right_arc_center = arc_center(right_sub_arch.right_arc)

        sub_arcs_radius = arc_radius(left_sub_arch.left_arc)

        # Rosette
        rosette = compute_rosette(bottom_left_corner, upper_right_corner, 
                                    arcs.right_arc, left_sub_arch.right_arc, 
                                        outer_offset, inner_offset)
        rosette_center = rosette.rosette_center
        rosette_radius = rosette.rosette_radius
        three_dimensionalize_rosette(rosette_center, rosette_radius, rosette_radius * outer_offset_ratio, inner_offset)
        #rosette_rounded_foils(rosette_center, rosette_radius, 9, π/2, rosette_radius * inner_offset_ratio)
        rosette_pointed_foils(rosette_center, rosette_radius, 9, 2, π/2, rosette_radius * inner_offset_ratio)

        # Fillets
        circular_rosette_fillets(right_arc_center, left_arc_center, arcs_radius - outer_offset, 
                            rosette_center, rosette_radius + inner_offset, 
                                right_sub_arch_right_arc_center, right_sub_arch_left_arc_center,
                                    left_sub_arch_left_arc_center, left_sub_arch_right_arc_center, sub_arcs_radius + inner_offset)
    end

    return (left_arc = arcs.left_arc, right_arc = arcs.right_arc)
end
# == MAIN ARCH == #

#with(current_cs, cs_from_o_vz(u0(), vx())) do
#    arch(xy(-10, -16), xy(10, 16), 1, 0.75, 0.75, 1, 3)
#end

three_dimensional_arch(xy(-10, -16), xy(10, 16), 1, 2, 3, 1, 1)
