using KhepriAutoCAD

delete_all_shapes()

DEFAULT_PROFILE_RADIUS = 1

# == GEOMETRY UTILITY FUNCTIONS == #
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

    return nothing
end

function normalize_vector(v)
    return v / norm(v)
end

function displace_point_by_vector(point, vector, magnitude)
    unit_vector = normalize_vector(vector)

    return point + unit_vector * magnitude
end

function get_bidimensional_perpendicular_vector(v)
    return vxy(-v.y, v.x)
end

function get_perpendicular_to_vectors(v, w)
    return cross(v, w)
end

# == CIRCLES == #
function intersect_circles_ccw_intersection_order(m0, r0, m1, r1)
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

    m0_to_first_intersection_point = first_intersection_point - m0
    m0_to_second_intersection_point = second_intersection_point - m0

    θ0 = atan(m0_to_first_intersection_point.y, m0_to_first_intersection_point.x)
    θ1 = atan(m0_to_second_intersection_point.y, m0_to_second_intersection_point.x)

    if θ0 <= θ1
        ccw_first_intersection_point = first_intersection_point
        ccw_second_intersection_point = second_intersection_point
    else
        ccw_first_intersection_point = second_intersection_point
        ccw_second_intersection_point = first_intersection_point
    end

    return (ccw_first_intersection_point = ccw_first_intersection_point, ccw_second_intersection_point = ccw_second_intersection_point)
end

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

    if first_intersection_point.y >= second_intersection_point.y
        greater_y_intersection_point = first_intersection_point
        lower_y_intersection_point = second_intersection_point
    else
        greater_y_intersection_point = second_intersection_point
        lower_y_intersection_point = first_intersection_point
    end

    return (greater_y_intersection_point = greater_y_intersection_point, lower_y_intersection_point = lower_y_intersection_point)
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

    return (first_intersection_point = q0, second_intersection_point = q1, 
                first_intersection_point_parameter = t0, second_intersection_point_parameter = t1)
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

    return KhepriAutoCAD.arc(center, radius, ccw_start_angle, amplitude)
end

function surface_arc(arc)
    center = arc_center(arc)
    radius = arc_radius(arc)
    start_angle = arc_start_angle(arc)
    amplitude = arc_amplitude(arc)

    return KhepriAutoCAD.surface_arc(center, radius, start_angle, amplitude)
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

function get_angular_adjustment_from_offset(arc, offset_value)
    amplitude = arc_amplitude(arc)
    radius = arc_radius(arc)

    Δθ = amplitude / radius
    tangential_offset = offset_value * sin(Δθ / 2)
    angular_adjustment = tangential_offset / radius

    return angular_adjustment
end

function arc_bidirectionally_extended_uniform_offset(arc, offset_value)
    delete_shape(arc)
    center = arc_center(arc)
    radius = arc_radius(arc)
    start_angle = arc_start_angle(arc)
    amplitude = arc_amplitude(arc)

    angular_adjustment = get_angular_adjustment_from_offset(arc, offset_value)

    new_radius = radius - offset_value
    new_start_angle = start_angle - angular_adjustment
    new_amplitude = amplitude + angular_adjustment * 2

    return KhepriAutoCAD.arc(center, new_radius, new_start_angle, new_amplitude)
end

function arc_amplitude_extended_offset(arc, offset_value)
    delete_shape(arc)
    center = arc_center(arc)
    radius = arc_radius(arc)
    start_angle = arc_start_angle(arc)
    amplitude = arc_amplitude(arc)

    angular_adjustment = get_angular_adjustment_from_offset(arc, offset_value) * 2

    new_radius = radius - offset_value
    new_start_angle = start_angle + angular_adjustment

    return KhepriAutoCAD.arc(center, new_radius, new_start_angle, amplitude)
end

function arc_start_angle_extended_offset(arc, offset_value)
    delete_shape(arc)
    center = arc_center(arc)
    radius = arc_radius(arc)
    start_angle = arc_start_angle(arc)
    amplitude = arc_amplitude(arc)

    angular_adjustment = get_angular_adjustment_from_offset(arc, offset_value) * 2

    new_radius = radius - offset_value
    new_start_angle = start_angle - angular_adjustment

    return KhepriAutoCAD.arc(center, new_radius, new_start_angle, amplitude)
end
# == ARCS == #
# == CIRCLES == #
# == GEOMETRY UTILITY FUNCTIONS == #

# == GOTHIC STRUCTURAL GEOMETRY == #
# == ARCHES == #
# == BODIES == #
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
# == BODIES == #

# == TOPS == #
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

function lancet_arch_top(left_point, right_point, excess)
    arch_midpoint = intermediate_loc(left_point, right_point)
    
    if isapprox(excess, 0.5)
        right_arc_center = left_arc_center = arch_midpoint
        arc_intersection = xy(arch_midpoint.x, arch_midpoint + vy(distance(left_point, arch_midpoint)))
    else
        excess_displacement_vector = right_point - left_point
        arcs_radius = norm(excess_displacement_vector) * excess

        right_arc_center = displace_point_by_vector(right_point, excess_displacement_vector, -arcs_radius)
        left_arc_center = displace_point_by_vector(left_point, excess_displacement_vector, arcs_radius)
        
        arc_intersection_x = arch_midpoint.x
        arc_intersection_y = sqrt(arcs_radius^2 - distance(left_arc_center, arch_midpoint)^2) + arch_midpoint.y
        arc_intersection = xy(arc_intersection_x, arc_intersection_y)
    end

    right_arc = arc(right_arc_center, right_point, arc_intersection)
    left_arc = arc(left_arc_center, arc_intersection, left_point)
    delete_shapes([right_arc, left_arc])

    return (right_arc = right_arc, left_arc = left_arc)
end
# == TOPS == #
# == ARCHES == #

# == ROSETTES == #
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
# == ROSETTES == #

# == FILLETS == #
function circular_rosette_fillets(right_arc_center, left_arc_center, arcs_radius,
                                    rosette_center, rosette_radius, 
                                        right_sub_arch_right_arc_center, right_sub_arch_left_arc_center, 
                                            left_sub_arch_left_arc_center, left_sub_arch_right_arc_center, sub_arcs_radius)
    # Circles intersections
    left_arc_and_right_arc_intersection = intersect_circles(left_arc_center, arcs_radius, right_arc_center, arcs_radius)
    top_fillet_top_point = left_arc_and_right_arc_intersection.greater_y_intersection_point

    left_arc_and_rosette_intersection = intersect_circles(left_arc_center, arcs_radius, rosette_center, rosette_radius)
    left_fillet_top_point = left_arc_and_rosette_intersection.lower_y_intersection_point
    top_fillet_left_point = left_arc_and_rosette_intersection.greater_y_intersection_point

    right_arc_and_rosette_intersection = intersect_circles(right_arc_center, arcs_radius, rosette_center, rosette_radius)
    right_fillet_top_point = right_arc_and_rosette_intersection.lower_y_intersection_point
    top_fillet_right_point = right_arc_and_rosette_intersection.greater_y_intersection_point

    left_arc_and_left_sub_arch_left_arc_intersection = intersect_circles(left_arc_center, arcs_radius, left_sub_arch_left_arc_center, sub_arcs_radius)
    left_fillet_bottom_point = left_arc_and_left_sub_arch_left_arc_intersection.greater_y_intersection_point

    right_arc_and_right_sub_arch_right_arc_intersection = intersect_circles(right_arc_center, arcs_radius, right_sub_arch_right_arc_center, sub_arcs_radius)
    right_fillet_bottom_point = right_arc_and_right_sub_arch_right_arc_intersection.greater_y_intersection_point

    left_sub_arch_left_arc_and_rosette_intersection = intersect_circles_ccw_intersection_order(left_sub_arch_left_arc_center, sub_arcs_radius,
                                                                                                rosette_center, rosette_radius)
    left_fillet_right_point = left_sub_arch_left_arc_and_rosette_intersection.ccw_second_intersection_point

    right_sub_arch_right_arc_and_rosette_intersection = intersect_circles_ccw_intersection_order(right_sub_arch_right_arc_center, sub_arcs_radius,
                                                                                                    rosette_center, rosette_radius)
    right_fillet_left_point = right_sub_arch_right_arc_and_rosette_intersection.ccw_first_intersection_point

    left_sub_arch_right_arc_and_rosette_intersection = intersect_circles(left_sub_arch_right_arc_center, sub_arcs_radius,
                                                                            rosette_center, rosette_radius)
    bottom_fillet_left_point = left_sub_arch_right_arc_and_rosette_intersection.lower_y_intersection_point

    right_sub_arch_left_arc_and_rosette_intersection = intersect_circles(right_sub_arch_left_arc_center, sub_arcs_radius,
                                                                            rosette_center, rosette_radius)
    bottom_fillet_right_point = right_sub_arch_left_arc_and_rosette_intersection.lower_y_intersection_point

    left_sub_arch_right_arc_and_right_sub_arch_left_arc_intersection = intersect_circles(left_sub_arch_right_arc_center, sub_arcs_radius,
                                                                                            right_sub_arch_left_arc_center, sub_arcs_radius)
    bottom_fillet_bottom_point = left_sub_arch_right_arc_and_right_sub_arch_left_arc_intersection.greater_y_intersection_point

    # Top fillet modeling
    arc(right_arc_center, top_fillet_right_point, top_fillet_top_point)
    arc(left_arc_center, top_fillet_top_point, top_fillet_left_point)
    arc(rosette_center, top_fillet_right_point, top_fillet_left_point)

    # Left fillet modeling
    arc(rosette_center, left_fillet_top_point, left_fillet_right_point)
    arc(left_arc_center, left_fillet_top_point, left_fillet_bottom_point)
    arc(left_sub_arch_left_arc_center, left_fillet_right_point, left_fillet_bottom_point)

    # Right fillet modeling
    arc(right_arc_center, right_fillet_bottom_point, right_fillet_top_point)
    arc(rosette_center, right_fillet_left_point, right_fillet_top_point)
    arc(right_sub_arch_right_arc_center, right_fillet_bottom_point, right_fillet_left_point)

    # Bottom fillet modeling
    arc(right_sub_arch_left_arc_center, bottom_fillet_right_point, bottom_fillet_bottom_point)
    arc(rosette_center, bottom_fillet_left_point, bottom_fillet_right_point)
    arc(left_sub_arch_right_arc_center, bottom_fillet_bottom_point, bottom_fillet_left_point)
end
# == FILLETS == #
# == GOTHIC STRUCTURAL GEOMETRY == #

# == GOTHIC ORNAMENTAL GEOMETRY == #
# == SOLID ORNAMENTATIONS == #
# == ARCHES == #
function three_dimensionalize_arch_top(left_arc, right_arc, offset_value, profile)
    half_offset = offset_value / 2
    profile = surface(scale(profile, abs(half_offset) / DEFAULT_PROFILE_RADIUS, u0()))
    left_arc = offset_arc(left_arc, half_offset)
    right_arc = offset_arc(right_arc, half_offset)
    delete_shapes([left_arc, right_arc])

    return union(sweep(left_arc, profile), sweep(right_arc, profile))
end

function three_dimensionalize_arch_body(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner, offset_value, profile)
    half_offset = offset_value / 2
    profile = surface(move(scale(profile, abs(half_offset) / DEFAULT_PROFILE_RADIUS, u0()), vx(half_offset)))
    arch_body = line(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner)
    delete_shape(arch_body)

    return sweep(arch_body, profile)
end

function three_dimensionalize_arch_middle(bottom_left_corner, upper_right_corner, vertical_distance_to_sub_arch, offset_value, profile)
    half_offset = offset_value / 2
    profile = surface(scale(profile, abs(half_offset) / DEFAULT_PROFILE_RADIUS, u0()))
    bottom_right_corner = xy(upper_right_corner.x, bottom_left_corner.y)
    bottom_midpoint = intermediate_loc(bottom_left_corner, bottom_right_corner)

    height = distance(upper_right_corner, bottom_right_corner)
    vertical_displacement = height - vertical_distance_to_sub_arch

    upper_midpoint = bottom_midpoint + vy(vertical_displacement)
    arch_middle = line(upper_midpoint, bottom_midpoint)
    delete_shape(arch_middle)

    return sweep(arch_middle, profile)
end
# == ARCHES == #

# == ROSETTES == #
function three_dimensionalize_rosette(rosette_center, rosette_radius, offset_value, profile)
    half_offset = offset_value / 2
    profile = surface(scale(profile, half_offset / DEFAULT_PROFILE_RADIUS, u0()))
    rosette = circle(rosette_center, rosette_radius + half_offset)
    delete_shape(rosette)

    return sweep(rosette, profile)
end

function three_dimensionalize_rosette_rounded_foils(foil, connection, offset_value, profile)
    half_offset = offset_value / 2
    profile = surface(scale(profile, abs(half_offset) / DEFAULT_PROFILE_RADIUS, u0()))
    foil = arc_bidirectionally_extended_uniform_offset(foil, half_offset)
    connection = nothing
    
    foil = sweep(foil, profile)

    return (foil = foil, connection = connection)
end

function three_dimensionalize_rosette_pointed_foils(foil_right_arc, foil_left_arc, connection, offset_value, profile)
    half_offset = offset_value / 2
    profile = surface(scale(profile, abs(half_offset) / DEFAULT_PROFILE_RADIUS, u0()))
    foil_right_arc = arc_bidirectionally_extended_uniform_offset(foil_right_arc, half_offset)
    foil_left_arc = arc_bidirectionally_extended_uniform_offset(foil_left_arc, half_offset)
    connection = nothing

    foil_right_arc = sweep(foil_right_arc, profile)
    foil_left_arc = sweep(foil_left_arc, profile)

    return (foil_right_arc = foil_right_arc, foil_left_arc = foil_left_arc, connection = connection)
end
# == ROSETTES == #
# == SOLID ORNAMENTATIONS == #

# == ROSETTES == #
function get_rounded_foil(rosette_center, rosette_radius, Δα, orientation)
    foil_radius = (rosette_radius * sin(Δα/2)) / (1 + sin(Δα/2))
    rosette_center_to_foil_center_length = rosette_radius - foil_radius
    rosette_center_to_foil_end_points = rosette_center_to_foil_center_length * cos(Δα/2)

    foil_center = rosette_center + vpol(rosette_center_to_foil_center_length, orientation)
    foil_start_point = rosette_center + vpol(rosette_center_to_foil_end_points, orientation - (Δα/2))
    foil_end_point = rosette_center + vpol(rosette_center_to_foil_end_points, orientation + (Δα/2))

    return arc(foil_center, foil_start_point, foil_end_point)
end

function get_pointed_foil(rosette_center, rosette_radius, Δα, displacement_ratio, orientation)
    rounded_foil = get_rounded_foil(rosette_center, rosette_radius, Δα, orientation)
    delete_shape(rounded_foil)
    rounded_foil_center = arc_center(rounded_foil)
    rounded_foil_start_point = arc_start_point(rounded_foil)
    rounded_foil_end_point = arc_end_point(rounded_foil)

    pointed_foil_right_arc_center_displacement_vector = rounded_foil_center - rounded_foil_start_point
    pointed_foil_left_arc_center_displacement_vector = rounded_foil_center - rounded_foil_end_point
    displacement_vectors_magnitude = norm(pointed_foil_right_arc_center_displacement_vector) * displacement_ratio

    pointed_foil_right_arc_center = displace_point_by_vector(rounded_foil_start_point, pointed_foil_right_arc_center_displacement_vector, displacement_vectors_magnitude)
    pointed_foil_left_arc_center = displace_point_by_vector(rounded_foil_end_point, pointed_foil_left_arc_center_displacement_vector, displacement_vectors_magnitude)
    pointed_foil_radius = distance(pointed_foil_right_arc_center, rounded_foil_start_point)

    arc_intersection_axis_point = intermediate_loc(pointed_foil_right_arc_center, pointed_foil_left_arc_center)
    pointed_foil_right_arc_center_to_arc_intersection_axis_point_length = distance(pointed_foil_right_arc_center, arc_intersection_axis_point)
    rounded_foil_center_to_arc_intersection_vector = arc_intersection_axis_point - rounded_foil_center
    rounded_foil_center_to_arc_intersection_vector_magnitude = sqrt(pointed_foil_radius^2 - pointed_foil_right_arc_center_to_arc_intersection_axis_point_length^2)
    arc_intersection = displace_point_by_vector(arc_intersection_axis_point, rounded_foil_center_to_arc_intersection_vector, rounded_foil_center_to_arc_intersection_vector_magnitude)

    pointed_foil_right_arc = arc(pointed_foil_right_arc_center, rounded_foil_start_point, arc_intersection) 
    pointed_foil_left_arc = arc(pointed_foil_left_arc_center, arc_intersection, rounded_foil_end_point)

    return (right_arc = pointed_foil_right_arc, left_arc = pointed_foil_left_arc, rounded_foil_center = rounded_foil_center)
end

# == FILLETS == #
function get_rosette_rounded_foils_fillet(rosette_center, rosette_radius, 
                                                    foil_center, foil_radius, 
                                                        Δα, inner_offset)
    rosette_radius = rosette_radius - inner_offset
    rosette_center_to_foil_center_length = distance(rosette_center, foil_center)

    fillet_right_arc_center = rosette_center + vpol(rosette_center_to_foil_center_length, 0)
    fillet_upper_arc_center = rosette_center
    fillet_left_arc_center = rosette_center + vpol(rosette_center_to_foil_center_length, Δα)

    displacement_vector = fillet_left_arc_center - fillet_right_arc_center

    fillet_right_point = intersect_circles(rosette_center, rosette_radius, fillet_right_arc_center, foil_radius).greater_y_intersection_point
    fillet_left_point = intersect_circles(rosette_center, rosette_radius, fillet_left_arc_center, foil_radius).lower_y_intersection_point
    fillet_bottom_point = displace_point_by_vector(fillet_right_arc_center, displacement_vector, norm(displacement_vector) / 2)

    fillet_right_arc = arc(fillet_right_arc_center, fillet_right_point, fillet_bottom_point)
    fillet_upper_arc = arc(fillet_upper_arc_center, fillet_right_point, fillet_left_point)
    fillet_left_arc = arc(fillet_left_arc_center, fillet_bottom_point, fillet_left_point)

    return (right_arc = fillet_right_arc, upper_arc = fillet_upper_arc, left_arc = fillet_left_arc)
end

function get_rosette_pointed_foils_fillet(rosette_center, rosette_radius, 
                                                    foil_center, right_arc_center, foil_radius, 
                                                        Δα, scaling_factor, inner_offset)
    scaling_factor = 1 / scaling_factor
    rosette_radius = rosette_radius * scaling_factor - inner_offset

    rosette_center_to_foil_center_vector = foil_center - rosette_center
    rosette_center_to_foil_arc_center_vector = right_arc_center - rosette_center
    rosette_center_to_foil_center_length = norm(rosette_center_to_foil_center_vector)
    rosette_center_to_foil_arc_center_length = norm(rosette_center_to_foil_arc_center_vector)
    Δβ = angle_between(rosette_center_to_foil_center_vector, rosette_center_to_foil_arc_center_vector)

    right_foil_center = rosette_center + vpol(rosette_center_to_foil_center_length, 0)
    left_foil_center = rosette_center + vpol(rosette_center_to_foil_center_length, Δα)

    fillet_right_arc_center = rosette_center + vpol(rosette_center_to_foil_arc_center_length, -Δβ)
    fillet_upper_arc_center = rosette_center
    fillet_left_arc_center = rosette_center + vpol(rosette_center_to_foil_arc_center_length, Δα + Δβ)

    displacement_vector = left_foil_center - right_foil_center

    fillet_right_point = intersect_circles(rosette_center, rosette_radius, fillet_right_arc_center, foil_radius).greater_y_intersection_point
    fillet_left_point = intersect_circles(rosette_center, rosette_radius, fillet_left_arc_center, foil_radius).lower_y_intersection_point
    fillet_bottom_point = displace_point_by_vector(right_foil_center, displacement_vector, norm(displacement_vector) / 2)

    fillet_right_arc = arc(fillet_right_arc_center, fillet_right_point, fillet_bottom_point)
    fillet_upper_arc = arc(fillet_upper_arc_center, fillet_right_point, fillet_left_point)
    fillet_left_arc = arc(fillet_left_arc_center, fillet_bottom_point, fillet_left_point)

    return (right_arc = fillet_right_arc, upper_arc = fillet_upper_arc, left_arc = fillet_left_arc)
end
# == FILLETS == #

function rosette_rounded_foils(rosette_center, rosette_radius, n_foils, orientation, inner_offset, profile)
    # Check if n_foils >= 1, otherwise raise an error
    Δα = 2π / n_foils

    foil = get_rounded_foil(rosette_center, rosette_radius, Δα, orientation)
    #fillet = get_rosette_rounded_foils_fillet(rosette_center, rosette_radius, 
    #                                            arc_center(foil), arc_radius(foil), 
    #                                                Δα, inner_offset)
    foil = arc_bidirectionally_extended_uniform_offset(foil, inner_offset)
    current_rotation_angle = 0

    while n_foils > 0
        # Fillets
        #fillet_arc_position = current_rotation_angle + orientation
        #rotate(fillet.right_arc, fillet_arc_position, rosette_center)
        #rotate(fillet.upper_arc, fillet_arc_position, rosette_center)
        #rotate(fillet.left_arc, fillet_arc_position, rosette_center)

        # Rounded foils
        rotate(foil, current_rotation_angle, rosette_center)

        # Rounded foils connections
        #connection_start_point = arc_end_point(foil)
        #rosette_center_to_foil_start_point = arc_start_point(foil) - rosette_center
        #rosette_center_to_foil_start_point_polar_angle = atan(rosette_center_to_foil_start_point.y, rosette_center_to_foil_start_point.x)
        #connection_end_point = rosette_center + vpol(norm(rosette_center_to_foil_start_point), rosette_center_to_foil_start_point_polar_angle + Δα)
        #connection = line(connection_start_point, connection_end_point)
        #rotate(connection, current_rotation_angle, rosette_center)

        # 3D rosette rounded foils
        three_dimensionalized_foil_and_connection = three_dimensionalize_rosette_rounded_foils(foil, connection, -inner_offset, profile)
        rotate(three_dimensionalized_foil_and_connection.foil, current_rotation_angle, rosette_center)
        #rotate(three_dimensionalized_foil_and_connection.connection, current_rotation_angle, rosette_center)

        current_rotation_angle += Δα
        n_foils -= 1
    end
end

function rosette_pointed_foils(rosette_center, rosette_radius, n_foils, displacement_ratio, orientation, inner_offset, profile)
    # Check if n_foils >= 1 and if displacement_ratio > 1, otherwise raise an error
    Δα = 2π / n_foils

    foil = get_pointed_foil(rosette_center, rosette_radius, Δα, displacement_ratio, orientation)
    rounded_foil_center = foil.rounded_foil_center
    outer_foil_right_arc = foil.right_arc 
    outer_foil_left_arc = foil.left_arc
    foil_radius = arc_radius(outer_foil_right_arc)
    foil_right_arc = arc_start_angle_extended_offset(outer_foil_right_arc, inner_offset)
    foil_left_arc = arc_amplitude_extended_offset(outer_foil_left_arc, inner_offset)
    
    scaling_factor = rosette_radius / distance(rosette_center, arc_end_point(foil_right_arc))
    #fillet = get_rosette_pointed_foils_fillet(rosette_center, rosette_radius, 
    #                                            rounded_foil_center, arc_center(outer_foil_right_arc), foil_radius, 
    #                                                Δα, scaling_factor, inner_offset)

    current_rotation_angle = 0

    while n_foils > 0
        # Fillets
        #fillet_arc_position = current_rotation_angle + orientation
        #scale(rotate(fillet.right_arc, fillet_arc_position, rosette_center), scaling_factor, rosette_center)
        #scale(rotate(fillet.upper_arc, fillet_arc_position, rosette_center), scaling_factor, rosette_center)
        #scale(rotate(fillet.left_arc, fillet_arc_position, rosette_center), scaling_factor, rosette_center)

        # Pointed foils
        scale(rotate(foil_right_arc, current_rotation_angle, rosette_center), scaling_factor, rosette_center)
        scale(rotate(foil_left_arc, current_rotation_angle, rosette_center), scaling_factor, rosette_center)

        # Pointed foils connections
        #connection_start_point = arc_end_point(foil_left_arc)
        #rosette_center_to_foil_start_point = arc_start_point(foil_right_arc) - rosette_center
        #rosette_center_to_foil_start_point_polar_angle = atan(rosette_center_to_foil_start_point.y, rosette_center_to_foil_start_point.x)
        #connection_end_point = rosette_center + vpol(norm(rosette_center_to_foil_start_point), rosette_center_to_foil_start_point_polar_angle + Δα)
        #connection = line(connection_start_point, connection_end_point)
        #scale(rotate(connection, current_rotation_angle, rosette_center), scaling_factor, rosette_center)

        # 3D rosette pointed foils
        three_dimensionalized_foil_and_connection = three_dimensionalize_rosette_pointed_foils(outer_foil_right_arc, outer_foil_left_arc, 
                                                                                                    connection, inner_offset, profile)
        scale(rotate(three_dimensionalized_foil_and_connection.foil_right_arc, current_rotation_angle, rosette_center), scaling_factor, rosette_center)
        scale(rotate(three_dimensionalized_foil_and_connection.foil_left_arc, current_rotation_angle, rosette_center), scaling_factor, rosette_center)
        #scale(rotate(three_dimensionalized_foil_and_connection.connection, current_rotation_angle, rosette_center), scaling_factor, rosette_center)

        n_foils -= 1
        current_rotation_angle += Δα
    end
end
# == ROSETTES == #
# == GOTHIC ORNAMENTAL GEOMETRY == #

# == GOTHIC WINDOWS == #
function gothic_window(window)
    bottom_left_corner = get_window_bottom_left_corner(window)
    upper_right_corner = get_window_upper_right_corner(window)
    excess = get_window_excess(window)
    vertical_distance_to_sub_arch = get_window_vertical_distance_to_sub_arch(window)
    outer_offset = get_window_outer_offset(window)
    inner_offset = get_window_inner_offset(window)
    profile = get_window_profile(window)
    rosette_style = get_window_rosette_style(window)
    left_sub_arch_style = deepcopy(get_window_left_sub_arch_style(window))
    right_sub_arch_style = deepcopy(get_window_right_sub_arch_style(window))

    # Arch body auxiliary coordinates
    upper_left_corner = xy(bottom_left_corner.x, upper_right_corner.y)
    bottom_right_corner = xy(upper_right_corner.x, bottom_left_corner.y)

    # Arch top portion
    arcs = lancet_arch_top(upper_left_corner, upper_right_corner, excess)
    right_arc_center = arc_center(arcs.right_arc)
    left_arc_center = arc_center(arcs.left_arc)
    arcs_radius = arc_radius(arcs.right_arc)

    # 3D Arch
    three_dimensional_window = three_dimensionalize_arch_top(arcs.left_arc, arcs.right_arc, outer_offset, profile)
    three_dimensional_window = union(three_dimensional_window,
                                        three_dimensionalize_arch_body(upper_left_corner, bottom_left_corner, 
                                                                        bottom_right_corner, upper_right_corner, 
                                                                            outer_offset, profile))

    # Sub-Arches
    if !isnothing(left_sub_arch_style) && !isnothing(right_sub_arch_style)
        # 3D Arch continuation
        three_dimensional_window = union(three_dimensional_window, 
                                            three_dimensionalize_arch_middle(bottom_left_corner, upper_right_corner, 
                                                                                vertical_distance_to_sub_arch, inner_offset, profile))

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
        three_dimensional_window = union(three_dimensional_window, 
                                            three_dimensionalize_arch_top(left_outer_sub_arch_top.left_arc, 
                                                                            left_outer_sub_arch_top.right_arc, inner_offset, profile))
        three_dimensional_window = union(three_dimensional_window, 
                                            three_dimensionalize_arch_top(right_outer_sub_arch_top.left_arc, 
                                                                            right_outer_sub_arch_top.right_arc, inner_offset, profile))

        # Left sub-arch style additions
        set_window_bottom_left_corner(left_sub_arch_style, left_sub_arch_body.bottom_left_corner)
        set_window_upper_right_corner(left_sub_arch_style, left_sub_arch_body.upper_right_corner)
        set_window_excess(left_sub_arch_style, sub_arch_excess)
        set_window_vertical_distance_to_sub_arch(left_sub_arch_style, vertical_distance_to_sub_arch)
        set_window_outer_offset(left_sub_arch_style, sub_arch_outer_offset)
        set_window_inner_offset(left_sub_arch_style, sub_arch_inner_offset)

        # Right sub-arch style additions
        set_window_bottom_left_corner(right_sub_arch_style, right_sub_arch_body.bottom_left_corner)
        set_window_upper_right_corner(right_sub_arch_style, right_sub_arch_body.upper_right_corner)
        set_window_excess(right_sub_arch_style, sub_arch_excess)
        set_window_vertical_distance_to_sub_arch(right_sub_arch_style, vertical_distance_to_sub_arch)
        set_window_outer_offset(right_sub_arch_style, sub_arch_outer_offset)
        set_window_inner_offset(right_sub_arch_style, sub_arch_inner_offset)

        # Sub-arches
        left_sub_arch = gothic_window(left_sub_arch_style)
        right_sub_arch = gothic_window(right_sub_arch_style)

        three_dimensional_window = union(three_dimensional_window, left_sub_arch.three_dimensional_window)
        three_dimensional_window = union(three_dimensional_window, right_sub_arch.three_dimensional_window)

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
        rosette_profile = get_rosette_profile(rosette_style)
        rosette_foil_instantiator = get_rosette_foil_instantiator(rosette_style)

        three_dimensional_window = union(three_dimensional_window, three_dimensionalize_rosette(rosette_center, rosette_radius, 
                                                                                                    inner_offset, rosette_profile))

        if rosette_foil_instantiator === rosette_rounded_foils
            rosette_n_foils = get_rosette_n_foils(rosette_style)
            rosette_starting_foil_orientation = get_rosette_starting_foil_orientation(rosette_style)
            rosette_foil_instantiator(rosette_center, rosette_radius, rosette_n_foils, 
                                        rosette_starting_foil_orientation, rosette_radius * inner_offset_ratio, rosette_profile)
        elseif rosette_foil_instantiator === rosette_pointed_foils
            rosette_n_foils = get_rosette_n_foils(rosette_style)
            rosette_starting_foil_orientation = get_rosette_starting_foil_orientation(rosette_style)
            rosette_displacement_ratio = get_rosette_displacement_ratio(rosette_style)
            rosette_foil_instantiator(rosette_center, rosette_radius, rosette_n_foils,
                                        rosette_displacement_ratio, rosette_starting_foil_orientation, 
                                            rosette_radius * inner_offset_ratio, rosette_profile)
        end

        # Fillets
        #circular_rosette_fillets(right_arc_center, left_arc_center, arcs_radius - outer_offset, 
        #                            rosette_center, rosette_radius + inner_offset, 
        #                                right_sub_arch_right_arc_center, right_sub_arch_left_arc_center,
        #                                    left_sub_arch_left_arc_center, left_sub_arch_right_arc_center, sub_arcs_radius + inner_offset)
    end

    return (left_arc = arcs.left_arc, right_arc = arcs.right_arc, three_dimensional_window = three_dimensional_window)
end
# == GOTHIC WINDOWS == #

# == GOTHIC PLAYGROUND == #
mutable struct GothicWindow
    bottom_left_corner
    upper_right_corner
    excess
    vertical_distance_to_sub_arch
    outer_offset
    inner_offset
    window_profile
    rosette_style
    left_sub_arch_style
    right_sub_arch_style
end

function get_window_bottom_left_corner(window)
    return window.bottom_left_corner
end

function get_window_upper_right_corner(window)
    return window.upper_right_corner
end

function get_window_excess(window)
    return window.excess
end

function get_window_vertical_distance_to_sub_arch(window)
    return window.vertical_distance_to_sub_arch
end

function get_window_outer_offset(window)
    return window.outer_offset
end

function get_window_inner_offset(window)
    return window.inner_offset
end

function get_window_profile(window)
    return window.window_profile
end

function get_window_rosette_style(window)
    return window.rosette_style
end

function get_window_left_sub_arch_style(window)
    return window.left_sub_arch_style
end

function get_window_right_sub_arch_style(window)
    return window.right_sub_arch_style
end

function set_window_bottom_left_corner(window, bottom_left_corner)
    window.bottom_left_corner = bottom_left_corner
end

function set_window_upper_right_corner(window, upper_right_corner)
    window.upper_right_corner = upper_right_corner
end

function set_window_excess(window, excess)
    window.excess = excess
end

function set_window_vertical_distance_to_sub_arch(window, vertical_distance_to_sub_arch)
    window.vertical_distance_to_sub_arch = vertical_distance_to_sub_arch
end

function set_window_outer_offset(window, outer_offset)
    window.outer_offset = outer_offset
end

function set_window_inner_offset(window, inner_offset)
    window.inner_offset = inner_offset
end

function Gothic_Window(bottom_left_corner, upper_right_corner, 
                        excess, vertical_distance_to_sub_arch, 
                            outer_offset, inner_offset, window_profile, 
                                rosette_style, left_sub_arch_style, right_sub_arch_style)
    main_arch = GothicWindow(bottom_left_corner, upper_right_corner,
                                excess, vertical_distance_to_sub_arch,
                                    outer_offset, inner_offset, window_profile, 
                                        rosette_style, left_sub_arch_style, right_sub_arch_style)
    return gothic_window(main_arch)
end

function Sub_Gothic_Window(window_profile, rosette_style, left_sub_arch_style, right_sub_arch_style)
    return GothicWindow(nothing, nothing, nothing, nothing, nothing, nothing, 
                            window_profile, rosette_style, left_sub_arch_style, right_sub_arch_style)
end

mutable struct Rosette
    foils_instantiator
    n_foils
    starting_foil_orientation
    displacement_ratio
    rosette_profile
end

function get_rosette_profile(rosette)
    return rosette.rosette_profile
end

function get_rosette_foil_instantiator(rosette)
    return rosette.foils_instantiator
end

function get_rosette_n_foils(rosette)
    return rosette.n_foils
end

function get_rosette_starting_foil_orientation(rosette)
    return rosette.starting_foil_orientation
end

function get_rosette_displacement_ratio(rosette)
    return rosette.displacement_ratio
end

function Empty_Rosette_Style(rosette_profile)
    return Rosette(nothing, nothing, nothing, nothing, rosette_profile)
end

function Rosette_Pointed_Style(n_foils, starting_foil_orientation, displacement_ratio, rosette_profile)
    return Rosette(rosette_pointed_foils, n_foils, starting_foil_orientation, displacement_ratio, rosette_profile)
end

function Rosette_Rounded_Style(n_foils, starting_foil_orientation, rosette_profile)
    return Rosette(rosette_rounded_foils, n_foils, starting_foil_orientation, nothing, rosette_profile)
end

#struct Fillets
#    right_fillet_profile
#    upper_fillet_profile
#    left_fillet_profile
#    bottom_fillet_profile
#end

# IMPORTANT - PROFILES MUST BE CENTERED IN THE ORIGIN AND SHOULD EITHER BE A CIRCLE OR A SHAPE INSCRIBED IN A CIRCLES
# THE CIRCLE'S RADIUS SHOULD ALWAYS BE EQUAL TO 1 SINCE IT WILL BE SCALED PROVIDED THE OUTER/INNER OFFSET AS IT PROPAGATES THROUGH THE MODELING PROCESS
# TO EXEMPLIFY WE WILL RESORT TO 3 DIFFERENT PROFILES - A CIRCLE AND TWO DIFFERENT KINDS OF STARS
# A WINDOW'S STYLE IS DIVIDED INTO TWO FACTORS - ITS 3D SHAPE (I.E., THE MODEL'S GEOMETRIC OUTLINE) AND ITS SUB-FIELDS
# THE MODEL'S 3D SHAPE IS MATERIALIZED VIA A 2D PROFILE THAT IS SWEPT ALONG THE MEDIAL LINE BETWEEN THE WINDOW'S BORDER, FILLETS AND SUB-ARCHES
# SUCH IS DETERMINED BY ITS OUTER AND INNER OFFSETS
# WHILST THE SUB-FIELDS, COMPOSED OF SUB-ARCHES, FILLETS, AND ROSETTES, CAN INCORPORATE TWO DIFFERENT STYLISTIC APPROACHES IN A RECURSIVE FASHION,
# (I.E., THE SAME AS THE MAIN ARCH, IN THE SENSE THAT ITS 3D SHAPE CAN ALSO BE PROVIDED AS WELL AS ITS SHAPES, WITH THE LATTER RELYING ON PARAMETRIC
# VARIATIONS, COMBINATIONS, AND HYBRIDIZATIONS)
# FOLLOWING THIS RATIONALE, ONE CAM DERIVE A MULTITUDE OF DIFFERENT STYLES PROVIDED ENDLESS COMBINATIONS OF PARAMETRIC AND 3D SWEEPING NATURE
# TO EXEMPLIFY THE EXPRESSIVENESS THAT ONE CAN ACCOMPLISH WHEN MODELING GOTHIC WINDOWS WITH OUR STRATEGY, WE WILL RESORT TO:
# 1ST STYLE - A DEFAULT WINDOW THAT SHARES THE EXACT SAME VALUES ACCROSS ITS RECURSIVE STAGES (CIRCLE_PROFILE, EMPTY ROSETTE)
# 2ND STYLE - A MAIN ARCH (CIRCLE_PROFILE, ROUNDED 6-FOIL ROSETTE), SUB-ARCHES (QUAD_STAR_PROFILE, POINTED 3-FOIL ROSETTE)
# 3RD STYLE - A MAIN ARCH (PENTAGON_STAR_PROFILE, POINTED 3-FOIL ROSETTE WITH CIRCLE_PROFILE, PI/2 ORIENTED, AND 2 DISPLACEMENT RATIO), 
# LEFT SUB-ARCH (CIRCLE_PROFILE, ROUNDED 6-FOIL ROSETTE WITH QUAD_START_PROFILE 0 ORIENTED),
# RIGHT SUB-ARCH (QUAD_STAR_PROFULE, POINTED 9-FOIL ROSETTE WITH CIRCLE_PROFILE, PI/4 ORIENTED, AND 3 DISPLACEMENT RATIO)
# HOWEVER, THE EXPRESSIVENESS DOESN'T END HERE. WE CAN REFER TO DIFFERENT PARAMETERS WITH RESPECT TO THE WINDOWS' EXCESS, VERTICAL DISTANCE TO SUB ARCHES
# AND OUTER AND INNER OFFSET! TO EACH OF THESE STYLES WE WILL APPLY DIFFERENT VALUES PERTAINING TO THESE FACTORS

# == 1ST STYLE == #
function Gothic_Window_First_Style(bottom_left_corner, upper_right_corner, excess, vertical_distance_to_sub_arch, outer_offset, inner_offset)
    window_profile = surface_circle(u0(), DEFAULT_PROFILE_RADIUS)

    rosette_style = Empty_Rosette_Style(window_profile)

    sub_arches_style = Sub_Gothic_Window(window_profile, rosette_style, nothing, nothing)

    main_arch = Gothic_Window(bottom_left_corner, upper_right_corner,
                                excess, vertical_distance_to_sub_arch,
                                    outer_offset, inner_offset, window_profile,
                                        rosette_style, sub_arches_style, sub_arches_style)

    return main_arch
end

first_style_window = Gothic_Window_First_Style(xy(-10, -16), xy(10, 16), 1, 3, 1, 1)
#delete_shape(first_style_window.three_dimensional_window)
# == 1ST STYLE == #

# == 2ND STYLE == #
#function Gothic_Window_Second_Style(bottom_left_corner, upper_right_corner, excess, vertical_distance_to_sub_arch, outer_offset, inner_offset)
#    main_arch_profile = surface_circle(u0(), DEFAULT_PROFILE_RADIUS)
#    sub_arches_profile = union(surface(regular_polygon(4, u0(), DEFAULT_PROFILE_RADIUS, 0)), 
#                                surface(regular_polygon(4, u0(), DEFAULT_PROFILE_RADIUS, π/4)))
#
#    main_arch_rosette_style = Rosette_Rounded_Style(6, 0, main_arch_profile)
#    sub_arches_rosette_style = Rosette_Pointed_Style(3, 0, 2, sub_arches_profile)
#
#    sub_sub_arches_style = Sub_Gothic_Window(sub_arches_profile, nothing, nothing, nothing)
#
#    sub_arches_style = Sub_Gothic_Window(sub_arches_profile, sub_arches_rosette_style, sub_sub_arches_style, sub_sub_arches_style)
#
#    main_arch = Gothic_Window(bottom_left_corner, upper_right_corner,
#                                excess, vertical_distance_to_sub_arch,
#                                    outer_offset, inner_offset, main_arch_profile,
#                                        main_arch_rosette_style, sub_arches_style, sub_arches_style)
#
#    return main_arch
#end
#
#second_style_window = Gothic_Window_Second_Style(xy(-10, -16), xy(10, 16), 1, 3, 1, 1)
# == 2ND STYLE == #

# == 3RD STYLE == #
#function Gothic_Window_Third_Style(bottom_left_corner, upper_right_corner,
#                                        excess, vertical_distance_to_sub_arch,
#                                            outer_offset, inner_offset)
#    main_arch_profile = union(surface(regular_polygon(4, u0(), DEFAULT_PROFILE_RADIUS, 0)), 
#                                surface(regular_polygon(4, u0(), DEFAULT_PROFILE_RADIUS, π/4)))
#
#    left_sub_arch_profile = surface_circle(u0(), DEFAULT_PROFILE_RADIUS)
#    right_sub_arch_profile = surface(regular_polygon(4, u0(), 1))
#
#    main_rosette_profile = left_sub_arch_profile
#
#    sub_left_rosette_profile = right_sub_arch_profile
#    sub_right_rosette_profile = main_arch_profile
#
#    main_arch_rosette_style = Rosette_Pointed_Style(3, π/2, 2, main_rosette_profile)
#    sub_left_rosette_style = Rosette_Rounded_Style(6, 0, sub_left_rosette_profile)
#    sub_right_rosette_style = Rosette_Pointed_Style(9, π/4, 5, sub_right_rosette_profile)
#
#    sub_sub_arches_style = Sub_Gothic_Window(main_arch_profile, nothing, nothing, nothing)
#
#    left_sub_arch_style = Sub_Gothic_Window(left_sub_arch_profile, sub_left_rosette_style, sub_sub_arches_style, sub_sub_arches_style)
#    right_sub_arch_style = Sub_Gothic_Window(right_sub_arch_profile, sub_right_rosette_style, sub_sub_arches_style, sub_sub_arches_style)
#
#    main_arch = Gothic_Window(bottom_left_corner, upper_right_corner, 
#                                excess, vertical_distance_to_sub_arch, 
#                                    outer_offset, inner_offset, main_arch_profile, 
#                                        main_arch_rosette_style, left_sub_arch_style, right_sub_arch_style)
#
#    return main_arch
#end
#
#third_style_window = Gothic_Window_Third_Style(xy(-10, -16), xy(10, 16), 1, 3, 1, 1)
# == 3RD STYLE == #

#gothic_window(xy(-10, -16), xy(10, 16), 1, 2, 3, 1, 1; three_dimensionality_enabled = false)
#gothic_window(xy(-10, -16), xy(10, 16), 2, 3, 3, 1, 1; three_dimensionality_enabled = false)
#gothic_window(xy(-10, -16), xy(10, 16), 0.8, 3, 3, 1, 1; three_dimensionality_enabled = false)

# == GOTHIC PLAYGROUND == #

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
AMBULATORY_INNER_PILLAR_WIDTH = AISLE_INNER_PILLAR_DEPTH
AMBULATORY_INNER_PILLAR_DEPTH = AISLE_INNER_PILLAR_WIDTH
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
AISLE_INNER_WALL_HEIGHT = AISLE_MIDDLE_WALL_HEIGHT

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
TRANSEPT_OUTER_WALL_HEIGHT = AISLE_INNER_WALL_HEIGHT

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
AMBULATORY_OUTER_WALL_HEIGHT = TRANSEPT_OUTER_WALL_HEIGHT

TOP_AMBULATORY_OUTER_WALL_LEFT_PILLARS = (4, range(14, 16))
TOP_AMBULATORY_OUTER_WALL_RIGHT_PILLARS = (4, range(15, 17))

BOTTOM_AMBULATORY_OUTER_WALL_LEFT_PILLARS = (7, range(14, 16))
BOTTOM_AMBULATORY_OUTER_WALL_RIGHT_PILLARS = (7, range(15, 17))
# == AMBULATORY OUTER WALLS == #

# == AMBULATORY MIDDLE WALLS == #
AMBULATORY_MIDDLE_WALL_HEIGHT = AMBULATORY_OUTER_WALL_HEIGHT

TOP_AMBULATORY_MIDDLE_WALL_LEFT_PILLARS = (4, range(14, 17))
TOP_AMBULATORY_MIDDLE_WALL_RIGHT_PILLARS = (5, range(14, 17))

BOTTOM_AMBULATORY_MIDDLE_WALL_LEFT_PILLARS = (6, range(14, 17))
BOTTOM_AMBULATORY_MIDDLE_WALL_RIGHT_PILLARS = (7, range(14, 17))
# == AMBULATORY MIDDLE WALLS == #
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
end

function Pillar(type, row, column, center, orientation, width, depth, height)
    Pillar(type, row, column, center, orientation, width, depth, height, nothing, nothing, nothing, nothing)
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

function standing_arch_wall_block(left_pillar, right_pillar, depth, height, arch_excess, offset)
    wall_anchors = get_wall_anchors(left_pillar, right_pillar)
    left_anchor = wall_anchors.left_anchor
    right_anchor = wall_anchors.right_anchor

    left_to_right_anchor_vector = right_anchor - left_anchor
    left_anchor_offset = displace_point_by_vector(left_anchor, left_to_right_anchor_vector, offset)
    right_anchor_offset = displace_point_by_vector(right_anchor, left_to_right_anchor_vector, -offset)

    arch_top_height = get_arch_top_height(left_anchor_offset, right_anchor_offset, arch_excess)
    wall_hole_height = height - arch_top_height - distance(get_pillar_center(left_pillar), left_anchor) * 2

    hollow_wall = standing_hollow_wall_block(left_anchor, right_anchor, depth, height, wall_hole_height, offset)
    arch_top = standing_lancet_arch_top_block(left_anchor_offset + GROWING_HEIGHT_DIRECTION * wall_hole_height, 
                                                right_anchor_offset + GROWING_HEIGHT_DIRECTION * wall_hole_height, 
                                                    depth, arch_excess)
    arch_wall = subtraction(hollow_wall, arch_top)

    return arch_wall
end

function standing_window_wall_block(left_pillar, right_pillar, depth, height, 
                                        arch_excess, vertical_distance_to_sub_arch, window_style_instantiator, offset)
    wall_anchors = get_wall_anchors(left_pillar, right_pillar)
    left_anchor = wall_anchors.left_anchor
    right_anchor = wall_anchors.right_anchor

    left_to_right_anchor_vector = right_anchor - left_anchor
    left_anchor_offset = displace_point_by_vector(left_anchor, left_to_right_anchor_vector, offset)
    right_anchor_offset = displace_point_by_vector(right_anchor, left_to_right_anchor_vector, -offset)
    window_to_wall_general_offset = distance(left_anchor_offset, right_anchor_offset) / 8

    bottom_left_corner = displace_point_by_vector(left_anchor_offset, left_to_right_anchor_vector, 
                            window_to_wall_general_offset) + GROWING_HEIGHT_DIRECTION * window_to_wall_general_offset
    bottom_right_corner = displace_point_by_vector(right_anchor_offset, left_to_right_anchor_vector,
                            -window_to_wall_general_offset) + GROWING_HEIGHT_DIRECTION * window_to_wall_general_offset

    arch_top_height = get_arch_top_height(bottom_left_corner, bottom_right_corner, arch_excess)

    upper_left_corner = bottom_left_corner + GROWING_HEIGHT_DIRECTION * (height - arch_top_height - window_to_wall_general_offset * 2)
    upper_right_corner = bottom_right_corner + GROWING_HEIGHT_DIRECTION * (height - arch_top_height - window_to_wall_general_offset * 2)
    
    window_body_hole = standing_wall_block(bottom_left_corner, bottom_right_corner, depth, height - window_to_wall_general_offset * 2 - arch_top_height)
    window_arch_hole = standing_lancet_arch_top_block(upper_left_corner, upper_right_corner, depth, arch_excess)
    window_hole = union(window_body_hole, window_arch_hole)

    wall = standing_wall_block(left_anchor_offset, right_anchor_offset, depth, height)
    window_wall = subtraction(wall, window_hole)
    
    #window = window_style_instantiator(bottom_left_corner, upper_right_corner, 
    #                                    excess, vertical_distance_to_sub_arch, depth, depth)
#
    #return nothing
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

function ambulatory_outer_arch(left_pillar, right_pillar;
                                depth = get_pillar_depth(left_pillar),
                                height = AMBULATORY_OUTER_WALL_HEIGHT,
                                excess = 1,
                                offset = 0)
    standing_arch_wall_block(left_pillar, right_pillar, depth, height, excess, offset)
end

function ambulatory_middle_arch(left_pillar, right_pillar;
                                    depth = get_pillar_depth(left_pillar),
                                    height = AMBULATORY_MIDDLE_WALL_HEIGHT,
                                    excess = 1,
                                    offset = 0)
    standing_arch_wall_block(left_pillar, right_pillar, depth, height, excess, offset)
end

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
function set_pillars_wall_attributes(left_pillar, right_pillar, wall, wall_type_instantiator)
    left_pillar_center = get_pillar_center(left_pillar)
    right_pillar_center = get_pillar_center(right_pillar)

    direction = right_pillar_center - left_pillar_center
    normalized_direction = direction / norm(direction)

    if wall_type_instantiator === ambulatory_outer_arch
        set_pillar_east_wall(left_pillar, wall)
        set_pillar_west_wall(right_pillar, wall)
    elseif wall_type_instantiator === ambulatory_middle_arch
        set_pillar_north_wall(left_pillar, wall)
        set_pillar_south_wall(right_pillar, wall)
    elseif isapprox(normalized_direction.x, SOUTH.x) && isapprox(normalized_direction.y, SOUTH.y)
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
        set_pillars_wall_attributes(left_pillar, right_pillar, wall_type, wall_type_instantiator)
    end
end

function right_pillar_increment_column_range_instantiator(pillars, wall_type_instantiator, column_range, left_pillars_row_index, right_pillars_row_index)
    for column in column_range
        left_pillar = pillars[left_pillars_row_index, column]
        right_pillar = pillars[right_pillars_row_index, column + 1]
        model = wall_type_instantiator(left_pillar, right_pillar)
        wall_type = Wall(model, left_pillar, right_pillar)
        set_pillars_wall_attributes(left_pillar, right_pillar, wall_type, wall_type_instantiator)
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
        set_pillars_wall_attributes(left_pillar, right_pillar, wall_type, wall_type_instantiator)
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

function instantiate_ambulatory_outer_walls(pillars)
    right_pillar_increment_column_range_instantiator(pillars, ambulatory_outer_arch, TOP_AMBULATORY_OUTER_WALL_LEFT_PILLARS[COLUMN_INDEX],
                                                        TOP_AMBULATORY_OUTER_WALL_LEFT_PILLARS[ROW_INDEX], TOP_AMBULATORY_OUTER_WALL_RIGHT_PILLARS[ROW_INDEX])
    right_pillar_increment_column_range_instantiator(pillars, ambulatory_outer_arch, BOTTOM_AMBULATORY_OUTER_WALL_LEFT_PILLARS[COLUMN_INDEX],
                                                        BOTTOM_AMBULATORY_OUTER_WALL_LEFT_PILLARS[ROW_INDEX], BOTTOM_AMBULATORY_OUTER_WALL_RIGHT_PILLARS[ROW_INDEX])
end

function instantiate_ambulatory_middle_walls(pillars)
    column_range_instantiator(pillars, ambulatory_middle_arch, TOP_AMBULATORY_MIDDLE_WALL_LEFT_PILLARS[COLUMN_INDEX],
                                TOP_AMBULATORY_MIDDLE_WALL_LEFT_PILLARS[ROW_INDEX], TOP_AMBULATORY_MIDDLE_WALL_RIGHT_PILLARS[ROW_INDEX])
    column_range_instantiator(pillars, ambulatory_middle_arch, BOTTOM_AMBULATORY_MIDDLE_WALL_LEFT_PILLARS[COLUMN_INDEX],
                                BOTTOM_AMBULATORY_MIDDLE_WALL_LEFT_PILLARS[ROW_INDEX], BOTTOM_AMBULATORY_MIDDLE_WALL_RIGHT_PILLARS[ROW_INDEX])
end

function instantiate_all_walls(pillars)
    instantiate_flying_buttresses(pillars)
    instantiate_main_hallways_walls(pillars)
    instantiate_aisles_outer_walls(pillars)
    instantiate_aisles_middle_walls(pillars)
    instantiate_aisles_inner_walls(pillars)
    instantiate_transept_outer_walls(pillars)
    instantiate_ambulatory_outer_walls(pillars)
    instantiate_ambulatory_middle_walls(pillars)
end
# == INSTANTIATORS == #
# == WALLS == #

# == PLAYGROUND == #
#pillars = Array{Union{Pillar, Nothing}}(nothing, 10, 17)
#
#instantiate_all_pillars(pillars)
#instantiate_all_walls(pillars)
# == PLAYGROUND == #

#circle = surface_circle(u0(), 2)
#pathA = line(x(5), xz(5, 10))
#pathB = line(x(15), xz(15, 10))
#sweep(pathA, circle)
#sweep(pathB, circle)
#pillarA = Pillar(nothing, nothing, nothing, x(5), NORTH, 4, 4, 10)
#pillarB = Pillar(nothing, nothing, nothing, x(15), NORTH, 4, 4, 10)
#standing_window_wall_block(pillarA, pillarB, get_pillar_depth(pillarA), get_pillar_height(pillarA), 1, 3, nothing, 0)
