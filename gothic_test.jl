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

function displace_point_by_vector(point, vector, magnitude)
    unit_vector = vector / norm(vector)

    return point + unit_vector * magnitude
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

    sweep(left_arc, profile)
    sweep(right_arc, profile)
end

function three_dimensionalize_arch_body(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner, offset_value, profile)
    half_offset = offset_value / 2
    profile = surface(move(scale(profile, abs(half_offset) / DEFAULT_PROFILE_RADIUS, u0()), vx(half_offset)))
    arch_body = line(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner)
    sweep(arch_body, profile)
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

    sweep(arch_middle, profile)
end
# == ARCHES == #

# == ROSETTES == #
function three_dimensionalize_rosette(rosette_center, rosette_radius, offset_value, profile)
    half_offset = offset_value / 2
    profile = surface(scale(profile, half_offset / DEFAULT_PROFILE_RADIUS, u0()))
    rosette = circle(rosette_center, rosette_radius + half_offset)
    sweep(rosette, profile)
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
function gothic_window(bottom_left_corner, upper_right_corner, 
                        excess, recursion_level, vertical_distance_to_sub_arch, 
                            outer_offset, inner_offset;
                                three_dimensionality_enabled = true)
    # Arch body auxiliary coordinates
    upper_left_corner = xy(bottom_left_corner.x, upper_right_corner.y)
    bottom_right_corner = xy(upper_right_corner.x, bottom_left_corner.y)

    # Arch body
    line(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner)

    # Arch top portion
    arcs = lancet_arch_top(upper_left_corner, upper_right_corner, excess)
    right_arc_center = arc_center(arcs.right_arc)
    left_arc_center = arc_center(arcs.left_arc)
    arcs_radius = arc_radius(arcs.right_arc)

    # 3D Arch
    if three_dimensionality_enabled
        three_dimensionalize_arch_top(arcs.left_arc, arcs.right_arc, outer_offset)
        three_dimensionalize_arch_body(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner, outer_offset)
    end
        
    # Sub-Arches
    if recursion_level > 0
        # 3D Arch continuation
        if three_dimensionality_enabled
            three_dimensionalize_arch_middle(bottom_left_corner, upper_right_corner, vertical_distance_to_sub_arch, inner_offset)
        end

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

        if three_dimensionality_enabled
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
        end

        # Sub-arches
        left_sub_arch = gothic_window(left_sub_arch_body.bottom_left_corner, left_sub_arch_body.upper_right_corner, sub_arch_excess, 
                                            recursion_level - 1, vertical_distance_to_sub_arch, 
                                                sub_arch_outer_offset, sub_arch_inner_offset; 
                                                    three_dimensionality_enabled = three_dimensionality_enabled)
        right_sub_arch = gothic_window(right_sub_arch_body.bottom_left_corner, right_sub_arch_body.upper_right_corner, sub_arch_excess, 
                                            recursion_level - 1, vertical_distance_to_sub_arch, 
                                                sub_arch_outer_offset, sub_arch_inner_offset;
                                                    three_dimensionality_enabled = three_dimensionality_enabled)

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

        if three_dimensionality_enabled
            three_dimensionalize_rosette(rosette_center, rosette_radius, inner_offset)
        end
        rosette_rounded_foils(rosette_center, rosette_radius, 9, π/2, rosette_radius * inner_offset_ratio)
        #rosette_pointed_foils(rosette_center, rosette_radius, 9, 2, π/2, rosette_radius * inner_offset_ratio)

        # Fillets
        circular_rosette_fillets(right_arc_center, left_arc_center, arcs_radius - outer_offset, 
                                    rosette_center, rosette_radius + inner_offset, 
                                        right_sub_arch_right_arc_center, right_sub_arch_left_arc_center,
                                            left_sub_arch_left_arc_center, left_sub_arch_right_arc_center, sub_arcs_radius + inner_offset)
    end

    return (left_arc = arcs.left_arc, right_arc = arcs.right_arc)
end

function gothic_window(window; three_dimensionality_enabled = true)
    bottom_left_corner = get_window_bottom_left_corner(window)
    upper_right_corner = get_window_upper_right_corner(window)
    excess = get_window_excess(window)
    recursion_level = get_window_recursion_level(window)
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

    # Arch body
    line(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner)

    # Arch top portion
    arcs = lancet_arch_top(upper_left_corner, upper_right_corner, excess)
    right_arc_center = arc_center(arcs.right_arc)
    left_arc_center = arc_center(arcs.left_arc)
    arcs_radius = arc_radius(arcs.right_arc)

    # 3D Arch
    if three_dimensionality_enabled
        three_dimensionalize_arch_top(arcs.left_arc, arcs.right_arc, outer_offset, profile)
        three_dimensionalize_arch_body(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner, outer_offset, profile)
    end

    # Sub-Arches
    if recursion_level > 0
        # 3D Arch continuation
        if three_dimensionality_enabled
            three_dimensionalize_arch_middle(bottom_left_corner, upper_right_corner, vertical_distance_to_sub_arch, inner_offset, profile)
        end

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

        if three_dimensionality_enabled
            # Outer sub-arches calculationss
            outer_sub_arch_excess = get_offset_excess(left_sub_arch_body.upper_left_corner, left_sub_arch_body.upper_right_corner, 
                                                        sub_arch_excess, -inner_offset)

            left_outer_sub_arch_top = lancet_arch_top(left_sub_arch_body.upper_left_corner - vx(inner_offset), left_sub_arch_body.upper_right_corner + vx(inner_offset), 
                                                        outer_sub_arch_excess)

            right_outer_sub_arch_top = lancet_arch_top(right_sub_arch_body.upper_left_corner - vx(inner_offset), right_sub_arch_body.upper_right_corner + vx(inner_offset),
                                                            outer_sub_arch_excess)

            # 3D outer arches
            three_dimensionalize_arch_top(left_outer_sub_arch_top.left_arc, left_outer_sub_arch_top.right_arc, inner_offset, profile)
            three_dimensionalize_arch_top(right_outer_sub_arch_top.left_arc, right_outer_sub_arch_top.right_arc, inner_offset, profile)
        end

        # Left sub-arch style additions
        set_window_bottom_left_corner(left_sub_arch_style, left_sub_arch_body.bottom_left_corner)
        set_window_upper_right_corner(left_sub_arch_style, left_sub_arch_body.upper_right_corner)
        set_window_excess(left_sub_arch_style, sub_arch_excess)
        set_window_recursion_level(left_sub_arch_style, recursion_level - 1)
        set_window_vertical_distance_to_sub_arch(left_sub_arch_style, vertical_distance_to_sub_arch)
        set_window_outer_offset(left_sub_arch_style, sub_arch_outer_offset)
        set_window_inner_offset(left_sub_arch_style, sub_arch_inner_offset)

        # Right sub-arch style additions
        set_window_bottom_left_corner(right_sub_arch_style, right_sub_arch_body.bottom_left_corner)
        set_window_upper_right_corner(right_sub_arch_style, right_sub_arch_body.upper_right_corner)
        set_window_excess(right_sub_arch_style, sub_arch_excess)
        set_window_recursion_level(right_sub_arch_style, recursion_level - 1)
        set_window_vertical_distance_to_sub_arch(right_sub_arch_style, vertical_distance_to_sub_arch)
        set_window_outer_offset(right_sub_arch_style, sub_arch_outer_offset)
        set_window_inner_offset(right_sub_arch_style, sub_arch_inner_offset)

        # Sub-arches
        left_sub_arch = gothic_window(left_sub_arch_style; three_dimensionality_enabled = three_dimensionality_enabled)
        right_sub_arch = gothic_window(right_sub_arch_style; three_dimensionality_enabled = three_dimensionality_enabled)

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

        if three_dimensionality_enabled
            three_dimensionalize_rosette(rosette_center, rosette_radius, inner_offset, rosette_profile)
        end

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
        circular_rosette_fillets(right_arc_center, left_arc_center, arcs_radius - outer_offset, 
                                    rosette_center, rosette_radius + inner_offset, 
                                        right_sub_arch_right_arc_center, right_sub_arch_left_arc_center,
                                            left_sub_arch_left_arc_center, left_sub_arch_right_arc_center, sub_arcs_radius + inner_offset)
    end

    return (left_arc = arcs.left_arc, right_arc = arcs.right_arc)
end

function get_extruded_gothic_window_outline(window)
    bottom_left_corner = get_window_bottom_left_corner(window)
    upper_right_corner = get_window_upper_right_corner(window)
    upper_left_corner = xy(bottom_left_corner.x, upper_right_corner.y)
    bottom_right_corner = xy(upper_right_corner.x, bottom_left_corner.y)
    depth = get_window_outer_offset(window)
    excess = get_window_excess(window)

    arch_top = standing_lancet_arch_top_block(upper_left_corner, upper_right_corner, depth, excess)
    arch_body_first_half = extrude(surface(polygon(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner)), GROWING_HEIGHT_DIRECTION * depth)
    arch_body_second_half = extrude(surface(polygon(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner)), DECREASING_HEIGHT_DIRECTION * depth)
    arch = union(arch_top, arch_body_first_half, arch_body_second_half)

    return arch
end
# == GOTHIC WINDOWS == #

# == GOTHIC PLAYGROUND == #
mutable struct GothicWindow
    bottom_left_corner
    upper_right_corner
    excess
    recursion_level
    vertical_distance_to_sub_arch
    outer_offset
    inner_offset
    profile
    rosette_style
    #fillets_style
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

function get_window_recursion_level(window)
    return window.recursion_level
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
    return window.profile
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

function set_window_recursion_level(window, recursion_level)
    window.recursion_level = recursion_level
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
                        excess, recursion_level, vertical_distance_to_sub_arch, 
                            outer_offset, inner_offset, profile, rosette_style,
                                left_sub_arch_style, right_sub_arch_style)
    main_arch = GothicWindow(bottom_left_corner, upper_right_corner,
                                excess, recursion_level, vertical_distance_to_sub_arch,
                                    outer_offset, inner_offset, profile, rosette_style,
                                        left_sub_arch_style, right_sub_arch_style)
    gothic_window(main_arch)
end

function Sub_Gothic_Window(profile, rosette_style, left_sub_arch_style, right_sub_arch_style)
    return GothicWindow(nothing, nothing, nothing, nothing, nothing, nothing, nothing, profile, rosette_style, left_sub_arch_style, right_sub_arch_style)
end

mutable struct Rosette
    profile
    foils_instantiator
    n_foils
    starting_foil_orientation
    displacement_ratio
end

function get_rosette_profile(rosette)
    return rosette.profile
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

function Empty_Rosette_Style(profile)
    return Rosette(profile, nothing, nothing, nothing, nothing)
end

function Rosette_Pointed_Style(profile, n_foils, starting_foil_orientation, displacement_ratio)
    return Rosette(profile, rosette_pointed_foils, n_foils, starting_foil_orientation, displacement_ratio)
end

function Rosette_Rounded_Style(profile, n_foils, starting_foil_orientation)
    return Rosette(profile, rosette_rounded_foils, n_foils, starting_foil_orientation, nothing)
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
#function Gothic_Window_First_Style(profile, excess, recursion_level, vertical_distance_to_sub_arch, outer_offset, inner_offset)
#    rosette_style = Empty_Rosette_Style(profile)
#
#    sub_sub_arches_style = Sub_Gothic_Window(profile, nothing, nothing, nothing)
#
#    sub_arches_style = Sub_Gothic_Window(profile, rosette_style, sub_sub_arches_style, sub_sub_arches_style)
#
#    main_arch = Gothic_Window(xy(-10, -16), xy(10, 16), excess, recursion_level, vertical_distance_to_sub_arch, 
#                                outer_offset, inner_offset, profile, 
#                                    rosette_style, sub_arches_style, sub_arches_style)
#
#    return main_arch
#end
#
#circle_profile = surface_circle(u0(), DEFAULT_PROFILE_RADIUS)
#first_style_window = Gothic_Window_First_Style(circle_profile, 1, 2, 3, 1, 1)
# == 1ST STYLE == #

# == 2ND STYLE == #
#function Gothic_Window_Second_Style(main_arch_profile, sub_arches_profile, excess, recursion_level, vertical_distance_to_sub_arch, outer_offset, inner_offset)
#    main_arch_rosette_style = Rosette_Rounded_Style(main_arch_profile, 6, 0)
#    sub_arches_rosette_style = Rosette_Pointed_Style(sub_arches_profile, 3, 0, 2)
#
#    sub_sub_arches_style = Sub_Gothic_Window(sub_arches_profile, nothing, nothing, nothing)
#
#    left_sub_arch_style = Sub_Gothic_Window(sub_arches_profile, sub_arches_rosette_style, sub_sub_arches_style, sub_sub_arches_style)
#    right_sub_arch_style = Sub_Gothic_Window(sub_arches_profile, sub_arches_rosette_style, sub_sub_arches_style, sub_sub_arches_style)
#
#    main_arch = Gothic_Window(xy(-10, -16), xy(10, 16), excess, recursion_level, vertical_distance_to_sub_arch, 
#                                outer_offset, inner_offset, main_arch_profile, main_arch_rosette_style, 
#                                    left_sub_arch_style, right_sub_arch_style)
#
#    return main_arch
#end
#
#circle_profile = surface_circle(u0(), DEFAULT_PROFILE_RADIUS)
#quad_star_profile = union(surface(regular_polygon(4, u0(), DEFAULT_PROFILE_RADIUS, 0)), surface(regular_polygon(4, u0(), DEFAULT_PROFILE_RADIUS, π/4)))
#Gothic_Window_Second_Style(circle_profile, quad_star_profile, 1, 2, 3, 1, 1)
# == 2ND STYLE == #

# 3RD STYLE - A MAIN ARCH (QUAD_STAR_PROFILE, POINTED 3-FOIL ROSETTE WITH CIRCLE_PROFILE, PI/2 ORIENTED, AND 2 DISPLACEMENT RATIO), 
# LEFT SUB-ARCH (CIRCLE_PROFILE, ROUNDED 6-FOIL ROSETTE WITH SQUARE_PROFILE 0 ORIENTED),
# RIGHT SUB-ARCH (SQUARE_PROFILE, POINTED 9-FOIL ROSETTE WITH CIRCLE_PROFILE, PI/4 ORIENTED, AND 5 DISPLACEMENT RATIO)
# == 3RD STYLE == #
#function Gothic_Window_Third_Style(main_arch_profile, left_sub_arch_profile, right_sub_arch_profile,
#                                    main_rosette_profile, sub_left_rosette_profile, sub_right_rosette_profile,
#                                        excess, recursion_level, vertical_distance_to_sub_arch,
#                                            outer_offset, inner_offset)
#    main_arch_rosette_style = Rosette_Pointed_Style(main_rosette_profile, 3, π/2, 2)
#    sub_left_rosette_style = Rosette_Rounded_Style(sub_left_rosette_profile, 6, 0)
#    sub_right_rosette_style = Rosette_Pointed_Style(sub_right_rosette_profile, 9, π/4, 5)
#
#    sub_sub_arches_style = Sub_Gothic_Window(main_arch_profile, nothing, nothing, nothing)
#
#    left_sub_arch_style = Sub_Gothic_Window(left_sub_arch_profile, sub_left_rosette_style, sub_sub_arches_style, sub_sub_arches_style)
#    right_sub_arch_style = Sub_Gothic_Window(right_sub_arch_profile, sub_right_rosette_style, sub_sub_arches_style, sub_sub_arches_style)
#
#    main_arch = Gothic_Window(xy(-10, -16), xy(10, 16), excess, recursion_level, vertical_distance_to_sub_arch, 
#                                outer_offset, inner_offset, main_arch_profile, main_arch_rosette_style, 
#                                    left_sub_arch_style, right_sub_arch_style)
#
#    return main_arch
#end
#
#circle_profile = surface_circle(u0(), DEFAULT_PROFILE_RADIUS)
#quad_star_profile = union(surface(regular_polygon(4, u0(), DEFAULT_PROFILE_RADIUS, 0)), surface(regular_polygon(4, u0(), DEFAULT_PROFILE_RADIUS, π/4)))
#square_profile = surface(regular_polygon(4, u0(), 1))
#Gothic_Window_Third_Style(quad_star_profile, 
#                            circle_profile, square_profile, 
#                                circle_profile, square_profile, quad_star_profile,
#                                    1, 2, 3, 1, 1)
# == 3RD STYLE == #

#gothic_window(xy(-10, -16), xy(10, 16), 1, 2, 3, 1, 1; three_dimensionality_enabled = false)
#gothic_window(xy(-10, -16), xy(10, 16), 2, 3, 3, 1, 1; three_dimensionality_enabled = false)
#gothic_window(xy(-10, -16), xy(10, 16), 0.8, 3, 3, 1, 1; three_dimensionality_enabled = false)

# == GOTHIC PLAYGROUND == #
