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
function angle_between_ccw_oriented(v, w)
    ccw_start_angle = atan(v.y, v.x)
    ccw_end_angle = atan(w.y, w.x)

    if ccw_end_angle > ccw_start_angle
        amplitude = ccw_end_angle - ccw_start_angle
    else
        amplitude = 2π - ccw_start_angle + ccw_end_angle
    end

    return amplitude
end

function arc(center, start_point, end_point)
    a = start_point - center
    b = end_point - center
    radius = norm(a)
    ccw_start_angle = atan(a.y, a.x)
    amplitude = angle_between_ccw_oriented(a, b)

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

function three_dimensionalize_rosette_rounded_foils(foil, offset_value, profile)
    half_offset = offset_value / 2
    profile = surface(scale(profile, abs(half_offset) / DEFAULT_PROFILE_RADIUS, u0()))
    foil = offset_arc(foil, half_offset)
    delete_shape(foil)

    foil = sweep(foil, profile)

    return foil
end

function three_dimensionalize_rosette_pointed_foils(foil_right_arc, foil_left_arc, offset_value, profile)
    half_offset = offset_value / 2
    profile = surface(scale(profile, abs(half_offset) / DEFAULT_PROFILE_RADIUS, u0()))
    foil_right_arc = offset_arc(foil_right_arc, half_offset)
    foil_left_arc = offset_arc(foil_left_arc, half_offset)
    delete_shapes([foil_right_arc, foil_left_arc])

    foil_right_arc = sweep(foil_right_arc, profile)
    foil_left_arc = sweep(foil_left_arc, profile)

    return (foil_right_arc = foil_right_arc, foil_left_arc = foil_left_arc)
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
    delete_shapes([pointed_foil_right_arc, pointed_foil_left_arc])

    return (right_arc = pointed_foil_right_arc, left_arc = pointed_foil_left_arc, rounded_foil_center = rounded_foil_center)
end

function rosette_rounded_foils(rosette_center, rosette_radius, n_foils, orientation, inner_offset, profile)
    three_dimensional_rounded_foils = []

    # Check if n_foils >= 1, otherwise raise an error
    Δα = 2π / n_foils

    while n_foils > 0
        # Rounded foils
        foil = get_rounded_foil(rosette_center, rosette_radius, Δα, orientation)
        delete_shape(foil)
        foil = offset_arc(foil, inner_offset)
        delete_shape(foil)

        # 3D rosette rounded foils
        push!(three_dimensional_rounded_foils, three_dimensionalize_rosette_rounded_foils(foil, -inner_offset, profile))

        orientation += Δα
        n_foils -= 1
    end

    return three_dimensional_rounded_foils
end

function rosette_pointed_foils(rosette_center, rosette_radius, n_foils, displacement_ratio, orientation, inner_offset, profile)
    three_dimensional_pointed_foils = []

    # Check if n_foils >= 1 and if displacement_ratio > 1, otherwise raise an error
    Δα = 2π / n_foils

    while n_foils > 0
        # Pointed foils
        foil = get_pointed_foil(rosette_center, rosette_radius, Δα, displacement_ratio, orientation)
        outer_foil_right_arc = foil.right_arc 
        outer_foil_left_arc = foil.left_arc

        scaling_factor = rosette_radius / distance(rosette_center, arc_end_point(outer_foil_right_arc))

        # 3D rosette pointed foils
        three_dimensionalized_foil = three_dimensionalize_rosette_pointed_foils(outer_foil_right_arc, outer_foil_left_arc, 
                                                                                    inner_offset, profile)
        
        push!(three_dimensional_pointed_foils, scale(three_dimensionalized_foil.foil_right_arc, scaling_factor, rosette_center))
        push!(three_dimensional_pointed_foils, scale(three_dimensionalized_foil.foil_left_arc, scaling_factor, rosette_center))

        n_foils -= 1
        orientation += Δα
    end

    return three_dimensional_pointed_foils
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
    three_dimensional_window = []

    # Arch body auxiliary coordinates
    upper_left_corner = xy(bottom_left_corner.x, upper_right_corner.y)
    bottom_right_corner = xy(upper_right_corner.x, bottom_left_corner.y)

    # Arch top portion
    arcs = lancet_arch_top(upper_left_corner, upper_right_corner, excess)

    # 3D Arch
    push!(three_dimensional_window, three_dimensionalize_arch_top(arcs.left_arc, arcs.right_arc, outer_offset, profile))
    push!(three_dimensional_window, three_dimensionalize_arch_body(upper_left_corner, bottom_left_corner, 
                                                                    bottom_right_corner, upper_right_corner, 
                                                                        outer_offset, profile))

    # Sub-Arches
    if !isnothing(left_sub_arch_style) && !isnothing(right_sub_arch_style)
        # 3D Arch continuation
        push!(three_dimensional_window, three_dimensionalize_arch_middle(bottom_left_corner, upper_right_corner, 
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
        push!(three_dimensional_window, three_dimensionalize_arch_top(left_outer_sub_arch_top.left_arc, left_outer_sub_arch_top.right_arc, 
                                                                        inner_offset, profile))
        push!(three_dimensional_window, three_dimensionalize_arch_top(right_outer_sub_arch_top.left_arc, right_outer_sub_arch_top.right_arc, 
                                                                        inner_offset, profile))

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

        push!(three_dimensional_window, left_sub_arch.three_dimensional_window)
        push!(three_dimensional_window, right_sub_arch.three_dimensional_window)

        # Rosette
        rosette = compute_rosette(bottom_left_corner, upper_right_corner, 
                                    arcs.right_arc, left_sub_arch.right_arc, 
                                        outer_offset, inner_offset)
        rosette_center = rosette.rosette_center
        rosette_radius = rosette.rosette_radius
        rosette_profile = get_rosette_profile(rosette_style)
        rosette_foil_instantiator = get_rosette_foil_instantiator(rosette_style)

        push!(three_dimensional_window, three_dimensionalize_rosette(rosette_center, rosette_radius, 
                                                                        inner_offset, rosette_profile))

        if rosette_foil_instantiator === rosette_rounded_foils
            rosette_n_foils = get_rosette_n_foils(rosette_style)
            rosette_starting_foil_orientation = get_rosette_starting_foil_orientation(rosette_style)

            push!(three_dimensional_window, rosette_foil_instantiator(rosette_center, rosette_radius, rosette_n_foils, 
                                                                        rosette_starting_foil_orientation, 
                                                                            rosette_radius * inner_offset_ratio, rosette_profile)...)
        elseif rosette_foil_instantiator === rosette_pointed_foils
            rosette_n_foils = get_rosette_n_foils(rosette_style)
            rosette_starting_foil_orientation = get_rosette_starting_foil_orientation(rosette_style)
            rosette_displacement_ratio = get_rosette_displacement_ratio(rosette_style)

            push!(three_dimensional_window, rosette_foil_instantiator(rosette_center, rosette_radius, rosette_n_foils,
                                                                        rosette_displacement_ratio, rosette_starting_foil_orientation, 
                                                                            rosette_radius * inner_offset_ratio, rosette_profile)...)
        end
    end

    return (left_arc = arcs.left_arc, right_arc = arcs.right_arc, three_dimensional_window = union(three_dimensional_window...))
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

    sub_sub_arches_style = Sub_Gothic_Window(window_profile, nothing, nothing, nothing)

    sub_arches_style = Sub_Gothic_Window(window_profile, rosette_style, sub_sub_arches_style, sub_sub_arches_style)

    main_arch = Gothic_Window(bottom_left_corner, upper_right_corner,
                                excess, vertical_distance_to_sub_arch,
                                    outer_offset, inner_offset, window_profile,
                                        rosette_style, sub_arches_style, sub_arches_style)

    return main_arch
end

#first_style_window = Gothic_Window_First_Style(xy(-10, -16), xy(10, 16), 1, 3, 1, 1)
#delete_shape(first_style_window.three_dimensional_window)
# == 1ST STYLE == #

# == 2ND STYLE == #
function Gothic_Window_Second_Style(bottom_left_corner, upper_right_corner, excess, vertical_distance_to_sub_arch, outer_offset, inner_offset)
    main_arch_profile = surface_circle(u0(), DEFAULT_PROFILE_RADIUS)
    sub_arches_profile = union(surface(regular_polygon(4, u0(), DEFAULT_PROFILE_RADIUS, 0)), 
                                surface(regular_polygon(4, u0(), DEFAULT_PROFILE_RADIUS, π/4)))
    delete_shape(sub_arches_profile)
    
    main_arch_rosette_style = Rosette_Rounded_Style(6, 0, main_arch_profile)
    sub_arches_rosette_style = Rosette_Pointed_Style(3, π/2, 2, sub_arches_profile)

    sub_sub_arches_style = Sub_Gothic_Window(sub_arches_profile, nothing, nothing, nothing)

    sub_arches_style = Sub_Gothic_Window(sub_arches_profile, sub_arches_rosette_style, sub_sub_arches_style, sub_sub_arches_style)

    main_arch = Gothic_Window(bottom_left_corner, upper_right_corner,
                                excess, vertical_distance_to_sub_arch,
                                    outer_offset, inner_offset, main_arch_profile,
                                        main_arch_rosette_style, sub_arches_style, sub_arches_style)

    return main_arch
end
#
#second_style_window = Gothic_Window_Second_Style(xy(-10, -16), xy(10, 16), 1, 3, 1, 1)
#delete_shape(second_style_window.three_dimensional_window)
# == 2ND STYLE == #

# == 3RD STYLE == #
function Gothic_Window_Third_Style(bottom_left_corner, upper_right_corner,
                                        excess, vertical_distance_to_sub_arch,
                                            outer_offset, inner_offset)
    main_arch_profile = union(surface(regular_polygon(4, u0(), DEFAULT_PROFILE_RADIUS, 0)), 
                                surface(regular_polygon(4, u0(), DEFAULT_PROFILE_RADIUS, π/4)))

    left_sub_arch_profile = surface_circle(u0(), DEFAULT_PROFILE_RADIUS)
    right_sub_arch_profile = surface(regular_polygon(4, u0(), 1))
    delete_shape(right_sub_arch_profile)

    main_rosette_profile = left_sub_arch_profile

    sub_left_rosette_profile = right_sub_arch_profile
    sub_right_rosette_profile = main_arch_profile

    main_arch_rosette_style = Rosette_Pointed_Style(3, π/2, 2, main_rosette_profile)
    sub_left_rosette_style = Rosette_Rounded_Style(6, 0, sub_left_rosette_profile)
    sub_right_rosette_style = Rosette_Pointed_Style(9, π/4, 5, sub_right_rosette_profile)

    sub_sub_arches_style = Sub_Gothic_Window(main_arch_profile, nothing, nothing, nothing)

    left_sub_arch_style = Sub_Gothic_Window(left_sub_arch_profile, sub_left_rosette_style, sub_sub_arches_style, sub_sub_arches_style)
    right_sub_arch_style = Sub_Gothic_Window(right_sub_arch_profile, sub_right_rosette_style, sub_sub_arches_style, sub_sub_arches_style)

    main_arch = Gothic_Window(bottom_left_corner, upper_right_corner, 
                                excess, vertical_distance_to_sub_arch, 
                                    outer_offset, inner_offset, main_arch_profile, 
                                        main_arch_rosette_style, left_sub_arch_style, right_sub_arch_style)

    return main_arch
end

#third_style_window = Gothic_Window_Third_Style(xy(-10, -16), xy(10, 16), 1, 3, 1, 1)
#delete_shape(third_style_window.three_dimensional_window)
# == 3RD STYLE == #

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
HORIZONTAL_WALL = :HORIZONTAL_WALL
VERTICAL_WALL = :VERTICAL_WALL
# == INDEXES == #

# == DIRECTIONS == #
NORTH = vy(1)
SOUTH = -NORTH
WEST = vx(-1)
EAST = -WEST
GROWING_HEIGHT_DIRECTION = vz(1)
DECREASING_HEIGHT_DIRECTION = -GROWING_HEIGHT_DIRECTION
# == DIRECTIONS == #

# == MEASUREMENTS AND DELIMITERS == #

# == MEASUREMENTS AND DELIMITERS DATA STRUCTURES == #
MEASUREMENTS_AND_DELIMITERS = Dict{Symbol, Any}()

function assign_measurements_and_delimiters_key_value(key, value)
    MEASUREMENTS_AND_DELIMITERS[key] =  value
end

function instantiate_measurements_and_delimiters(default_pillars_distance,
                                                    w_aisle_row_range, w_aisle_col_range,
                                                    e_aisle_row_range, e_aisle_col_range,
                                                    s_transept_row_range, s_transept_col_range,
                                                    n_transept_row_range, n_transept_col_range,
                                                    ambulatory_row_range, ambulatory_col_range)

    assign_measurements_and_delimiters_key_value(:DEFAULT_PILLARS_DISTANCE, default_pillars_distance)

    assign_measurements_and_delimiters_key_value(:W_AISLE_ROW_RANGE, w_aisle_row_range)
    assign_measurements_and_delimiters_key_value(:W_AISLE_COL_RANGE, w_aisle_col_range)
    assign_measurements_and_delimiters_key_value(:E_AISLE_ROW_RANGE, e_aisle_row_range)
    assign_measurements_and_delimiters_key_value(:E_AISLE_COL_RANGE, e_aisle_col_range)

    assign_measurements_and_delimiters_key_value(:S_TRANSEPT_ROW_RANGE, s_transept_row_range)
    assign_measurements_and_delimiters_key_value(:S_TRANSEPT_COL_RANGE, s_transept_col_range)
    assign_measurements_and_delimiters_key_value(:N_TRANSEPT_ROW_RANGE, n_transept_row_range)
    assign_measurements_and_delimiters_key_value(:N_TRANSEPT_COL_RANGE, n_transept_col_range)

    assign_measurements_and_delimiters_key_value(:AMBULATORY_ROW_RANGE, ambulatory_row_range)
    assign_measurements_and_delimiters_key_value(:AMBULATORY_COL_RANGE, ambulatory_col_range)

    ambulatory_center_x = default_pillars_distance * last(e_aisle_col_range) + default_pillars_distance
    ambulatory_center_y = default_pillars_distance * last(s_transept_row_range) + default_pillars_distance
    ambulatory_start_angle = atan(EAST.y, EAST.x)
    ambulatory_angle_increment = atan(NORTH.y, NORTH.x) / (length(ambulatory_col_range) + 1)
    assign_measurements_and_delimiters_key_value(:AMBULATORY_CENTER, xy(ambulatory_center_x, ambulatory_center_y))
    assign_measurements_and_delimiters_key_value(:AMBULATORY_START_ANGLE, ambulatory_start_angle)
    assign_measurements_and_delimiters_key_value(:AMBULATORY_ANGLE_INCREMENT, ambulatory_angle_increment)
end
# == MEASUREMENTS AND DELIMITERS DATA STRUCTURES == #

# == MEASUREMENTS AND DELIMITERS GETTERS == #
function get_measurements_and_delimiters_value(key)
    return MEASUREMENTS_AND_DELIMITERS[key]
end

function get_default_pillars_distance()
    return get_measurements_and_delimiters_value(:DEFAULT_PILLARS_DISTANCE)
end

function get_w_aisle_row_range()
    return get_measurements_and_delimiters_value(:W_AISLE_ROW_RANGE)
end

function get_w_aisle_col_range()
    return get_measurements_and_delimiters_value(:W_AISLE_COL_RANGE)
end

function get_e_aisle_row_range()
    return get_measurements_and_delimiters_value(:E_AISLE_ROW_RANGE)
end

function get_e_aisle_col_range()
    return get_measurements_and_delimiters_value(:E_AISLE_COL_RANGE)
end

function get_n_transept_row_range()
    return get_measurements_and_delimiters_value(:N_TRANSEPT_ROW_RANGE)
end

function get_n_transept_col_range()
    return get_measurements_and_delimiters_value(:N_TRANSEPT_COL_RANGE)
end

function get_s_transept_row_range()
    return get_measurements_and_delimiters_value(:S_TRANSEPT_ROW_RANGE)
end

function get_s_transept_col_range()
    return get_measurements_and_delimiters_value(:S_TRANSEPT_COL_RANGE)
end

function get_ambulatory_row_range()
    return get_measurements_and_delimiters_value(:AMBULATORY_ROW_RANGE)
end

function get_ambulatory_col_range()
    return get_measurements_and_delimiters_value(:AMBULATORY_COL_RANGE)
end

function get_ambulatory_center()
    return get_measurements_and_delimiters_value(:AMBULATORY_CENTER)
end

function get_ambulatory_start_angle()
    return get_measurements_and_delimiters_value(:AMBULATORY_START_ANGLE)
end

function get_ambulatory_angle_increment()
    return get_measurements_and_delimiters_value(:AMBULATORY_ANGLE_INCREMENT)
end

function get_horizontal_hallway_jump_prev_row()
    return last(get_s_transept_row_range())
end

function get_vertical_hallway_jump_prev_col()
    s_transept_col_range = get_s_transept_col_range()
    return s_transept_col_range[div(length(s_transept_col_range), 2)]
end
# == MEASUREMENTS AND DELIMITERS GETTERS == #
# == MEASUREMENTS & DELIMITERS == #

# == PILLARS == #
PILLARS_INFO = Dict{Symbol, Any}()

mutable struct Pillar
    model
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

function Pillar(model, row, column, center, orientation, width, depth, height)
    Pillar(model, row, column, center, orientation, width, depth, height, nothing, nothing, nothing, nothing)
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

function assign_pillars_info_key_value(key, value)
    PILLARS_INFO[key] = value
end

function instantiate_pillars_info(basic_pillar_width, basic_pillar_depth, basic_pillar_height)
    # == RANGES == #
    w_aisle_row_range = get_w_aisle_row_range()
    w_aisle_col_range = get_w_aisle_col_range()
    e_aisle_row_range = get_e_aisle_row_range()
    e_aisle_col_range = get_e_aisle_col_range()

    n_transept_row_range = get_n_transept_row_range()
    n_transept_col_range = get_n_transept_col_range()
    s_transept_row_range = get_s_transept_row_range()
    s_transept_col_range = get_s_transept_col_range()

    ambulatory_row_range = get_ambulatory_row_range()
    ambulatory_col_range = get_ambulatory_col_range()
    # == RANGES == #

    # == BASIC PILLARS == #
    assign_pillars_info_key_value(:BASIC_PILLAR_WIDTH, basic_pillar_width)
    assign_pillars_info_key_value(:BASIC_PILLAR_DEPTH, basic_pillar_depth)
    assign_pillars_info_key_value(:BASIC_PILLAR_HEIGHT, basic_pillar_height)
    # == BASIC PILLARS == #

    # == AISLE BUTTRESSES == #
    aisle_buttress_width = basic_pillar_width
    aisle_buttress_depth = basic_pillar_depth
    aisle_buttress_height = basic_pillar_height * 0.95

    ws_aisle_buttresses = (first(w_aisle_row_range), w_aisle_col_range[1:end-1], SOUTH)
    wn_aisle_buttresses = (last(w_aisle_row_range), w_aisle_col_range[1:end-1], NORTH)
    es_aisle_buttresses = (first(e_aisle_row_range), e_aisle_col_range[2:end], SOUTH)
    en_aisle_buttresses = (last(e_aisle_row_range), e_aisle_col_range[2:end], NORTH)
    aisles_buttresses = [ws_aisle_buttresses, wn_aisle_buttresses, 
                            es_aisle_buttresses, en_aisle_buttresses]

    assign_pillars_info_key_value(:AISLE_BUTTRESS_WIDTH, aisle_buttress_width)
    assign_pillars_info_key_value(:AISLE_BUTTRESS_DEPTH, aisle_buttress_depth)
    assign_pillars_info_key_value(:AISLE_BUTTRESS_HEIGHT, aisle_buttress_height)
    assign_pillars_info_key_value(:AISLES_BUTTRESSES, aisles_buttresses)
    # == AISLE BUTTRESSES == #

    # == OUTER CROSSING BUTTRESSES == #
    outer_crossing_buttress_width = aisle_buttress_width
    outer_crossing_buttress_depth = aisle_buttress_depth
    outer_crossing_buttress_height = aisle_buttress_height

    ws_outer_crossing_buttress = (first(w_aisle_row_range), last(w_aisle_col_range), WEST)
    wn_outer_crossing_buttress = (last(w_aisle_row_range), last(w_aisle_col_range), NORTH)
    es_outer_crossing_buttress = (first(e_aisle_row_range), first(e_aisle_col_range), SOUTH)
    en_outer_crossing_buttress = (last(e_aisle_row_range), first(e_aisle_col_range), EAST)
    outer_crossing_buttresses = [ws_outer_crossing_buttress, wn_outer_crossing_buttress, 
                                    es_outer_crossing_buttress, en_outer_crossing_buttress]

    assign_pillars_info_key_value(:OUTER_CROSSING_BUTTRESS_WIDTH, outer_crossing_buttress_width)
    assign_pillars_info_key_value(:OUTER_CROSSING_BUTTRESS_DEPTH, outer_crossing_buttress_depth)
    assign_pillars_info_key_value(:OUTER_CROSSING_BUTTRESS_HEIGHT, outer_crossing_buttress_height)
    assign_pillars_info_key_value(:OUTER_CROSSING_BUTTRESSES, outer_crossing_buttresses)
    # == OUTER CROSSING BUTTRESSES == #

    # == TRANSEPT BUTTRESSES == #
    transept_buttress_width = outer_crossing_buttress_width
    transept_buttress_depth = outer_crossing_buttress_depth
    transept_buttress_height = outer_crossing_buttress_height
    
    sw_transept_buttresses = (range(first(s_transept_row_range), first(w_aisle_row_range) - 1), last(w_aisle_col_range), WEST)
    se_transept_buttresses = (range(first(s_transept_row_range), first(e_aisle_row_range) - 1), first(e_aisle_col_range), EAST)
    nw_transept_buttresses = (range(last(w_aisle_row_range) + 1, last(n_transept_row_range)), last(w_aisle_col_range), WEST)
    ne_transept_buttresses = (range(last(e_aisle_row_range) + 1, last(n_transept_row_range)), first(e_aisle_col_range), EAST)
    transept_buttresses = [sw_transept_buttresses, se_transept_buttresses,
                            nw_transept_buttresses, ne_transept_buttresses]

    assign_pillars_info_key_value(:TRANSEPT_BUTTRESS_WIDTH, transept_buttress_width)
    assign_pillars_info_key_value(:TRANSEPT_BUTTRESS_DEPTH, transept_buttress_depth)
    assign_pillars_info_key_value(:TRANSEPT_BUTTRESS_HEIGHT, transept_buttress_height)
    assign_pillars_info_key_value(:TRANSEPT_BUTTRESSES, transept_buttresses)
    # == TRANSEPT BUTTRESSES == #

    # == AISLE OUTER PILLARS == #
    aisle_outer_pillar_width = aisle_buttress_width
    aisle_outer_pillar_depth = aisle_buttress_depth
    aisle_outer_pillar_height = basic_pillar_height * 1.025

    ws_aisle_outer_pillars = (w_aisle_row_range[2:div(length(w_aisle_row_range), 2) - 1], w_aisle_col_range, SOUTH)
    wn_aisle_outer_pillars = (w_aisle_row_range[div(length(w_aisle_row_range), 2) + 2:end-1], w_aisle_col_range, NORTH)
    es_aisle_outer_pillars = (e_aisle_row_range[2:div(length(e_aisle_row_range), 2)-1], e_aisle_col_range, SOUTH)
    en_aisle_outer_pillars = (e_aisle_row_range[div(length(e_aisle_row_range), 2) + 2:end-1], e_aisle_col_range, NORTH)
    aisle_outer_pillars = [ws_aisle_outer_pillars, wn_aisle_outer_pillars,
                            es_aisle_outer_pillars, en_aisle_outer_pillars]

    assign_pillars_info_key_value(:AISLE_OUTER_PILLAR_WIDTH, aisle_outer_pillar_width)
    assign_pillars_info_key_value(:AISLE_OUTER_PILLAR_DEPTH, aisle_outer_pillar_depth)
    assign_pillars_info_key_value(:AISLE_OUTER_PILLAR_HEIGHT, aisle_outer_pillar_height)
    assign_pillars_info_key_value(:AISLE_OUTER_PILLARS, aisle_outer_pillars)
    # == AISLE OUTER PILLARS == #

    # == AISLE INNER PILLARS == #
    aisle_inner_pillar_width = aisle_outer_pillar_width * 1.5
    aisle_inner_pillar_depth = aisle_outer_pillar_depth
    aisle_inner_pillar_height = basic_pillar_height * 1.10

    ws_aisle_inner_pillars = (w_aisle_row_range[div(length(w_aisle_row_range), 2)], w_aisle_col_range, SOUTH)
    wn_aisle_inner_pillars = (w_aisle_row_range[div(length(w_aisle_row_range), 2) + 1], w_aisle_col_range, NORTH)
    es_aisle_inner_pillars = (e_aisle_row_range[div(length(e_aisle_row_range), 2)], e_aisle_col_range, SOUTH)
    en_aisle_inner_pillars = (e_aisle_row_range[div(length(e_aisle_row_range), 2) + 1], e_aisle_col_range, NORTH)
    aisle_inner_pillars = [ws_aisle_inner_pillars, wn_aisle_inner_pillars,
                            es_aisle_inner_pillars, en_aisle_inner_pillars]

    assign_pillars_info_key_value(:AISLE_INNER_PILLAR_WIDTH, aisle_inner_pillar_width)
    assign_pillars_info_key_value(:AISLE_INNER_PILLAR_DEPTH, aisle_inner_pillar_depth)
    assign_pillars_info_key_value(:AISLE_INNER_PILLAR_HEIGHT, aisle_inner_pillar_height)
    assign_pillars_info_key_value(:AISLE_INNER_PILLARS, aisle_inner_pillars)
    # == AISLE INNER PILLARS == #

    # == TRANSEPT PILLARS == #
    transept_pillar_width = aisle_inner_pillar_width
    transept_pillar_depth = aisle_inner_pillar_depth
    transept_pillar_height = aisle_inner_pillar_height

    sw_transept_pillars = (s_transept_row_range[1:end-1], s_transept_col_range[1:div(length(s_transept_col_range), 2)], WEST)
    se_transept_pillars = (s_transept_row_range[1:end-1], s_transept_col_range[div(length(s_transept_col_range), 2) + 1:end], EAST)
    nw_transept_pillars = (n_transept_row_range[2:end], n_transept_col_range[1:div(length(n_transept_col_range), 2)], WEST)
    ne_transept_pillars = (n_transept_row_range[2:end], n_transept_col_range[div(length(n_transept_col_range), 2) + 1:end], EAST)
    transept_pillars = [sw_transept_pillars, se_transept_pillars,
                            nw_transept_pillars, ne_transept_pillars]

    assign_pillars_info_key_value(:TRANSEPT_PILLAR_WIDTH, transept_pillar_width)
    assign_pillars_info_key_value(:TRANSEPT_PILLAR_DEPTH, transept_pillar_depth)
    assign_pillars_info_key_value(:TRANSEPT_PILLAR_HEIGHT, transept_pillar_height)
    assign_pillars_info_key_value(:TRANSEPT_PILLARS, transept_pillars)
    # == TRANSEPT PILLARS == #

    # == CROSSING PILLARS == #
    crossing_pillar_width = transept_pillar_width
    crossing_pillar_depth = aisle_inner_pillar_width
    crossing_pillar_height = aisle_inner_pillar_height
    
    sw_crossing_pillar = (last(s_transept_row_range), s_transept_col_range[div(length(s_transept_col_range), 2)], WEST)
    se_crossing_pillar = (last(s_transept_row_range), s_transept_col_range[div(length(s_transept_col_range), 2) + 1], EAST)
    nw_crossing_pillar = (first(n_transept_row_range), n_transept_col_range[div(length(n_transept_col_range), 2)], WEST)
    ne_crossing_pillar = (first(n_transept_row_range), n_transept_col_range[div(length(n_transept_col_range), 2) + 1], EAST)
    crossing_pillars = [sw_crossing_pillar, se_crossing_pillar, 
                            nw_crossing_pillar, ne_crossing_pillar]

    assign_pillars_info_key_value(:CROSSING_PILLAR_WIDTH, crossing_pillar_width)
    assign_pillars_info_key_value(:CROSSING_PILLAR_DEPTH, crossing_pillar_depth)
    assign_pillars_info_key_value(:CROSSING_PILLAR_HEIGHT, crossing_pillar_height)
    assign_pillars_info_key_value(:CROSSING_PILLARS, crossing_pillars)
    # == CROSSING PILLARS == #

    # == AMBULATORY BUTTRESSES == #
    ambulatory_buttress_width = aisle_buttress_width
    ambulatory_buttress_depth = aisle_buttress_depth
    ambulatory_buttress_height = basic_pillar_height * 0.64
    
    n_ambulatory_buttresses = (last(ambulatory_row_range), ambulatory_col_range, nothing)
    s_ambulatory_buttresses = (first(ambulatory_row_range), ambulatory_col_range, nothing)
    ambulatory_buttresses = [n_ambulatory_buttresses, s_ambulatory_buttresses]

    assign_pillars_info_key_value(:AMBULATORY_BUTTRESS_WIDTH, ambulatory_buttress_width)
    assign_pillars_info_key_value(:AMBULATORY_BUTTRESS_DEPTH, ambulatory_buttress_depth)
    assign_pillars_info_key_value(:AMBULATORY_BUTTRESS_HEIGHT, ambulatory_buttress_height)
    assign_pillars_info_key_value(:AMBULATORY_BUTTRESSES, ambulatory_buttresses)
    # == AMBULATORY BUTTRESSES == #
    
    # == AMBULATORY OUTER PILLARS == #
    ambulatory_outer_pillar_width = aisle_outer_pillar_width
    ambulatory_outer_pillar_depth = aisle_outer_pillar_depth
    ambulatory_outer_pillar_height = aisle_outer_pillar_height
    
    s_ambulatory_outer_pillars = (ambulatory_row_range[2:div(length(ambulatory_row_range), 2) - 1], ambulatory_col_range, nothing)
    n_ambulatory_outer_pillars = (ambulatory_row_range[div(length(ambulatory_row_range), 2) + 2:end-1], ambulatory_col_range, nothing)
    ambulatory_outer_pillars = [s_ambulatory_outer_pillars, n_ambulatory_outer_pillars]

    assign_pillars_info_key_value(:AMBULATORY_OUTER_PILLAR_WIDTH, ambulatory_outer_pillar_width)
    assign_pillars_info_key_value(:AMBULATORY_OUTER_PILLAR_DEPTH, ambulatory_outer_pillar_depth)
    assign_pillars_info_key_value(:AMBULATORY_OUTER_PILLAR_HEIGHT, ambulatory_outer_pillar_height)
    assign_pillars_info_key_value(:AMBULATORY_OUTER_PILLARS, ambulatory_outer_pillars)
    # == AMBULATORY OUTER PILLARS == #
    
    # == AMBULATORY INNER PILLARS == #
    ambulatory_inner_pillar_width = aisle_inner_pillar_depth
    ambulatory_inner_pillar_depth = aisle_inner_pillar_width
    ambulatory_inner_pillar_height = aisle_inner_pillar_height
    
    s_ambulatory_inner_pillars = (ambulatory_row_range[div(length(ambulatory_row_range), 2)], ambulatory_col_range, nothing)
    n_ambulatory_inner_pillars = (ambulatory_row_range[div(length(ambulatory_row_range), 2) + 1], ambulatory_col_range, nothing)
    ambulatory_inner_pillars = [s_ambulatory_inner_pillars, n_ambulatory_inner_pillars]
    
    assign_pillars_info_key_value(:AMBULATORY_INNER_PILLAR_WIDTH, ambulatory_inner_pillar_width)
    assign_pillars_info_key_value(:AMBULATORY_INNER_PILLAR_DEPTH, ambulatory_inner_pillar_depth)
    assign_pillars_info_key_value(:AMBULATORY_INNER_PILLAR_HEIGHT, ambulatory_inner_pillar_height)
    assign_pillars_info_key_value(:AMBULATORY_INNER_PILLARS, ambulatory_inner_pillars)
    # == AMBULATORY INNER PILLARS == #
end

function get_pillars_info_value(key)
    return PILLARS_INFO[key]
end

function get_pillars_row_info(pillars)
    return pillars[ROW_INDEX]
end

function get_pillars_col_info(pillars)
    return pillars[COLUMN_INDEX]
end

function get_pillars_orientation_info(pillars)
    return pillars[ORIENTATION_INDEX]
end

# == PILLAR INFO GETTERS == #
function get_basic_pillar_width()
    return get_pillars_info_value(:BASIC_PILLAR_WIDTH)
end

function get_basic_pillar_depth()
    return get_pillars_info_value(:BASIC_PILLAR_DEPTH)
end

function get_basic_pillar_height()
    return get_pillars_info_value(:BASIC_PILLAR_HEIGHT)
end

function get_aisle_buttress_width()
    return get_pillars_info_value(:AISLE_BUTTRESS_WIDTH)
end

function get_aisle_buttress_depth()
    return get_pillars_info_value(:AISLE_BUTTRESS_DEPTH)
end

function get_aisle_buttress_height()
    return get_pillars_info_value(:AISLE_BUTTRESS_HEIGHT)
end

function get_aisles_buttresses()
    return get_pillars_info_value(:AISLES_BUTTRESSES)
end

function get_outer_crossing_buttress_width()
    return get_pillars_info_value(:OUTER_CROSSING_BUTTRESS_WIDTH)
end

function get_outer_crossing_buttress_depth()
    return get_pillars_info_value(:OUTER_CROSSING_BUTTRESS_DEPTH)
end

function get_outer_crossing_buttress_height()
    return get_pillars_info_value(:OUTER_CROSSING_BUTTRESS_HEIGHT)
end

function get_outer_crossing_buttresses()
    return get_pillars_info_value(:OUTER_CROSSING_BUTTRESSES)
end

function get_transept_buttress_width()
    return get_pillars_info_value(:TRANSEPT_BUTTRESS_WIDTH)
end

function get_transept_buttress_depth()
    return get_pillars_info_value(:TRANSEPT_BUTTRESS_DEPTH)
end

function get_transept_buttress_height()
    return get_pillars_info_value(:TRANSEPT_BUTTRESS_HEIGHT)
end

function get_transept_buttresses()
    return get_pillars_info_value(:TRANSEPT_BUTTRESSES)
end

function get_aisle_outer_pillar_width()
    return get_pillars_info_value(:AISLE_OUTER_PILLAR_WIDTH)
end

function get_aisle_outer_pillar_depth()
    return get_pillars_info_value(:AISLE_OUTER_PILLAR_DEPTH)
end

function get_aisle_outer_pillar_height()
    return get_pillars_info_value(:AISLE_OUTER_PILLAR_HEIGHT)
end

function get_aisle_outer_pillars()
    return get_pillars_info_value(:AISLE_OUTER_PILLARS)
end

function get_aisle_inner_pillar_width()
    return get_pillars_info_value(:AISLE_INNER_PILLAR_WIDTH)
end

function get_aisle_inner_pillar_depth()
    return get_pillars_info_value(:AISLE_INNER_PILLAR_DEPTH)
end

function get_aisle_inner_pillar_height()
    return get_pillars_info_value(:AISLE_INNER_PILLAR_HEIGHT)
end

function get_aisle_inner_pillars()
    return get_pillars_info_value(:AISLE_INNER_PILLARS)
end

function get_transept_pillar_width()
    return get_pillars_info_value(:TRANSEPT_PILLAR_WIDTH)
end

function get_transept_pillar_depth()
    return get_pillars_info_value(:TRANSEPT_PILLAR_DEPTH)
end

function get_transept_pillar_height()
    return get_pillars_info_value(:TRANSEPT_PILLAR_HEIGHT)
end

function get_transept_pillars()
    return get_pillars_info_value(:TRANSEPT_PILLARS)
end

function get_crossing_pillar_width()
    return get_pillars_info_value(:CROSSING_PILLAR_WIDTH)
end

function get_crossing_pillar_depth()
    return get_pillars_info_value(:CROSSING_PILLAR_DEPTH)
end

function get_crossing_pillar_height()
    return get_pillars_info_value(:CROSSING_PILLAR_HEIGHT)
end

function get_crossing_pillars()
    return get_pillars_info_value(:CROSSING_PILLARS)
end

function get_ambulatory_buttress_width()
    return get_pillars_info_value(:AMBULATORY_BUTTRESS_WIDTH)
end

function get_ambulatory_buttress_depth()
    return get_pillars_info_value(:AMBULATORY_BUTTRESS_DEPTH)
end

function get_ambulatory_buttress_height()
    return get_pillars_info_value(:AMBULATORY_BUTTRESS_HEIGHT)
end

function get_ambulatory_buttresses()
    return get_pillars_info_value(:AMBULATORY_BUTTRESSES)
end

function get_ambulatory_outer_pillar_width()
    return get_pillars_info_value(:AMBULATORY_OUTER_PILLAR_WIDTH)
end

function get_ambulatory_outer_pillar_depth()
    return get_pillars_info_value(:AMBULATORY_OUTER_PILLAR_DEPTH)
end

function get_ambulatory_outer_pillar_height()
    return get_pillars_info_value(:AMBULATORY_OUTER_PILLAR_HEIGHT)
end

function get_ambulatory_outer_pillars()
    return get_pillars_info_value(:AMBULATORY_OUTER_PILLARS)
end

function get_ambulatory_inner_pillar_width()
    return get_pillars_info_value(:AMBULATORY_INNER_PILLAR_WIDTH)
end

function get_ambulatory_inner_pillar_depth()
    return get_pillars_info_value(:AMBULATORY_INNER_PILLAR_DEPTH)
end

function get_ambulatory_inner_pillar_height()
    return get_pillars_info_value(:AMBULATORY_INNER_PILLAR_HEIGHT)
end

function get_ambulatory_inner_pillars()
    return get_pillars_info_value(:AMBULATORY_INNER_PILLARS)
end
# == PILLAR INFO GETTERS == #

function get_orientation_positional_angle(orientation::Union{VX, VXY})
    return atan(orientation.y, orientation.x)
end

function get_orientation_positional_angle(center::Union{X, XY})
    ambulatory_center_to_pillar_center_vector = center - get_ambulatory_center()
    return atan(ambulatory_center_to_pillar_center_vector.y, ambulatory_center_to_pillar_center_vector.x)
end

# == MODELERS == #
function basic_pillar(center; 
                        width = get_basic_pillar_width(), 
                        depth = get_basic_pillar_depth(), 
                        height = get_basic_pillar_height())
    half_width = width / 2
    half_depth = depth / 2

    bottom_left_corner = u0() - vxy(half_width, half_depth)
    upper_right_corner = u0() + vxy(half_width, half_depth)

    pillar_shape = surface_rectangle(bottom_left_corner, upper_right_corner)
    pillar_path = line(center, center + vz(height))
    basic_pillar = rotate(sweep(pillar_path, pillar_shape), π/2, center)

    return basic_pillar
end

function buttress(center, orientation;
                    width = get_aisle_buttress_width(),
                    depth = get_aisle_buttress_depth(), 
                    height = get_aisle_buttress_height())
    orientation = get_orientation_positional_angle(orientation)
    basic_pillar_width = get_basic_pillar_width()
    basic_pillar_depth = get_basic_pillar_depth()
    half_width = basic_pillar_width / 2
    quarter_width = basic_pillar_width / 4
    half_height = height / 2

    pillar = basic_pillar(center; width = width, depth = depth, height = height)

    base_bottom_left_corner = center - vx(basic_pillar_width)
    base_upper_right_corner = center + vxy(basic_pillar_width, basic_pillar_depth * 4)
    
    mid_bottom_left_corner = base_bottom_left_corner - vx(quarter_width) + vz(half_height)
    mid_upper_right_corner = base_upper_right_corner + vxyz(quarter_width, basic_pillar_depth / 2, half_height)

    post_mid_bottom_left_corner = base_bottom_left_corner + vz(half_height + half_width)
    post_mid_upper_right_corner = base_upper_right_corner + vz(half_height + half_width)

    top_bottom_left_corner = base_bottom_left_corner + vz(height)
    top_upper_right_corner = base_upper_right_corner + vz(height)

    buttress_base = extrusion(surface_rectangle(base_bottom_left_corner, base_upper_right_corner), vz(half_height))
    buttress_mid = surface_rectangle(mid_bottom_left_corner, mid_upper_right_corner)
    buttress_post_mid = surface_rectangle(post_mid_bottom_left_corner, post_mid_upper_right_corner)
    buttress_top = surface_rectangle(top_bottom_left_corner, top_upper_right_corner)

    buttress = rotate(union(pillar, buttress_base, loft_ruled([buttress_mid, buttress_post_mid, buttress_top])), -π/2 + orientation, center)

    return (model = buttress, width = width, depth = depth, height = height)
end

function aisle_outer_pillar(center, orientation;
                                width = get_aisle_outer_pillar_width(), 
                                depth = get_aisle_outer_pillar_depth(), 
                                height = get_aisle_outer_pillar_height())
    orientation = get_orientation_positional_angle(orientation)
    aisle_outer_pillar = rotate(basic_pillar(center; width = width, depth = depth, height = height), orientation, center)
    return (model = aisle_outer_pillar, width = width, depth = depth, height = height)
end

function aisle_inner_pillar(center, orientation; 
                                width = get_aisle_inner_pillar_width(),
                                depth = get_aisle_inner_pillar_depth(),
                                height = get_aisle_inner_pillar_height())
    orientation = get_orientation_positional_angle(orientation)
    aisle_inner_pillar = rotate(basic_pillar(center; width = width, depth = depth, height = height), orientation, center)
    return (model = aisle_inner_pillar, width = width, depth = depth, height = height)
end

function half_outer_crossing_buttress(center, height)
    basic_pillar_width = get_basic_pillar_width()
    basic_pillar_depth = get_basic_pillar_depth()
    half_width = basic_pillar_width / 2
    half_depth = basic_pillar_depth / 2
    half_height = height / 2

    base_bottom_left_corner = center - vxy(half_width, half_depth)
    base_upper_right_corner = base_bottom_left_corner + vxy(basic_pillar_width * 4, basic_pillar_depth)

    mid_bottom_left_corner = base_bottom_left_corner + vz(half_height)
    mid_upper_right_corner = base_upper_right_corner + vxyz(basic_pillar_width / 4, half_depth, half_height)

    post_mid_bottom_left_corner = base_bottom_left_corner + vz(half_height + half_width)
    post_mid_upper_right_corner = base_upper_right_corner + vz(half_height + half_width)

    top_bottom_left_corner = base_bottom_left_corner + vz(height)
    top_upper_right_corner = base_upper_right_corner + vz(height)

    half_outer_crossing_buttress_base = extrusion(surface_rectangle(base_bottom_left_corner, base_upper_right_corner), vz(half_height))
    half_outer_crossing_buttress_mid = surface_rectangle(mid_bottom_left_corner, mid_upper_right_corner)
    half_outer_crossing_buttress_post_mid = surface_rectangle(post_mid_bottom_left_corner, post_mid_upper_right_corner)
    half_outer_crossing_buttress_top = surface_rectangle(top_bottom_left_corner, top_upper_right_corner)

    return union(half_outer_crossing_buttress_base, loft_ruled([half_outer_crossing_buttress_mid, 
                                                                    half_outer_crossing_buttress_post_mid, half_outer_crossing_buttress_top]))
end

function outer_crossing_buttress(center, orientation;
                                    width = get_outer_crossing_buttress_width(),
                                    depth = get_outer_crossing_buttress_depth(),
                                    height = get_outer_crossing_buttress_height())
    orientation = get_orientation_positional_angle(orientation)
    basic_pillar_width = get_basic_pillar_width()
    basic_pillar_depth = get_basic_pillar_depth()

    pillar = basic_pillar(center; width = width, depth = depth, height = height)

    transformation_point = center + vxy(basic_pillar_width / 4, basic_pillar_depth / 4)
    right_buttress = move(half_outer_crossing_buttress(center, height), vxy(basic_pillar_width * 0.75, basic_pillar_depth * 0.75))
    left_buttress = mirror(rotate(deepcopy(right_buttress), π/2, transformation_point), transformation_point + vz())

    outer_crossing_buttress = rotate(union(pillar, right_buttress, left_buttress), orientation, center)

    return (model = outer_crossing_buttress, width = width, depth = depth, height = height)
end

function transept_pillar(center, orientation; 
                            width = get_transept_pillar_width(),
                            depth = get_transept_pillar_depth(),
                            height = get_transept_pillar_height())
    aisle_inner_pillar(center, orientation; width = width, depth = depth, height = height)
end

function crossing_pillar(center, orientation; 
                            width = get_crossing_pillar_width(),
                            depth = get_crossing_pillar_depth(),
                            height = get_crossing_pillar_height())
    orientation = get_orientation_positional_angle(orientation)
    crossing_pillar = rotate(basic_pillar(center; width = width, depth = depth, height = height), orientation, center)
    return (model = crossing_pillar, width = width, depth = depth, height = height)
end

function ambulatory_buttress(center, orientation;
                                width = get_ambulatory_buttress_width(),
                                depth = get_ambulatory_buttress_depth(), 
                                height = get_ambulatory_buttress_height())
    orientation = get_orientation_positional_angle(center)
    basic_pillar_width = get_basic_pillar_width()
    half_width = basic_pillar_width / 2
    basic_pillar_depth = get_basic_pillar_depth()
    third_depth = basic_pillar_depth / 3

    rectangle_bottom_left_corner = u0() - vx(basic_pillar_width) + vy(third_depth)
    rectangle_upper_right_corner = u0() + vxy(basic_pillar_width, basic_pillar_depth + third_depth)
    rectangle_bottom_right_corner = xy(rectangle_upper_right_corner.x, rectangle_bottom_left_corner.y)
    rectangle_upper_left_corner = xy(rectangle_bottom_left_corner.x, rectangle_upper_right_corner.y)

    trapezoid_bottom_left_corner = rectangle_bottom_left_corner - vxy(half_width, third_depth)
    trapezoid_bottom_right_corner = rectangle_bottom_right_corner + vx(half_width) - vy(third_depth)

    ambulatory_buttress_shape = surface_polygon(trapezoid_bottom_left_corner, trapezoid_bottom_right_corner, 
                                                    rectangle_bottom_right_corner, rectangle_upper_right_corner, 
                                                        rectangle_upper_left_corner, rectangle_bottom_left_corner)
    ambulatory_buttress_path = line(center, center + vz(height))

    ambulatory_buttress = rotate(sweep(ambulatory_buttress_path, ambulatory_buttress_shape), -π + orientation, center)
    
    return (model = ambulatory_buttress, width = width, depth = depth, height = height)
end

function ambulatory_outer_pillar(center, orientation;
                                    width = get_ambulatory_outer_pillar_width(),
                                    depth = get_ambulatory_outer_pillar_depth(),
                                    height = get_ambulatory_outer_pillar_height())
    orientation = get_orientation_positional_angle(center)
    basic_pillar_width = get_basic_pillar_width()
    third_width = basic_pillar_width / 3
    basic_pillar_depth = get_basic_pillar_depth()
    half_depth = basic_pillar_depth / 2

    upper_left_corner = u0() - vx(basic_pillar_width - third_width) + vy(half_depth)
    upper_right_corner = upper_left_corner + vx((basic_pillar_width - third_width) * 2)
    bottom_left_corner = u0() - vxy(basic_pillar_width / 2, half_depth)
    bottom_right_corner = bottom_left_corner + vx(basic_pillar_width)

    ambulatory_outer_pillar_shape = surface_polygon(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner)
    ambulatory_outer_pillar_path = line(center, center + vz(height))
    
    ambulatory_outer_pillar = rotate(sweep(ambulatory_outer_pillar_path, ambulatory_outer_pillar_shape), -π + orientation, center)

    return (model = ambulatory_outer_pillar, width = width, depth = depth, height = height)
end

function ambulatory_inner_pillar(center, orientation;
                                    width = get_ambulatory_inner_pillar_width(),
                                    depth = get_ambulatory_inner_pillar_depth(),
                                    height = get_ambulatory_inner_pillar_height())
    orientation = get_orientation_positional_angle(center)
    basic_pillar_width = get_basic_pillar_width()
    third_width = basic_pillar_width / 3
    eighth_width = basic_pillar_width / 8
    basic_pillar_depth = get_basic_pillar_depth()
    half_depth = basic_pillar_depth / 2

    upper_left_corner = u0() - vx(third_width) + vy(half_depth)
    upper_right_corner = upper_left_corner + vx(third_width * 2)
    bottom_left_corner = u0() - vxy(eighth_width, half_depth)
    bottom_right_corner = bottom_left_corner + vx(eighth_width * 2)

    ambulatory_inner_pillar_shape = surface_polygon(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner)
    ambulatory_inner_pillar_path = line(center, center + vz(height))
    
    ambulatory_inner_pillar = rotate(sweep(ambulatory_inner_pillar_path, ambulatory_inner_pillar_shape), -π + orientation, center)

    return (model = ambulatory_inner_pillar, width = width, depth = depth, height = height)
end
# == MODELERS == #

# == INSTANTIATORS == #
# For demonstration purposes only
#crossing_pillar(u0() - vxy(DEFAULT_PILLARS_DISTANCE, DEFAULT_PILLARS_DISTANCE * 2))
#outer_crossing_buttress(u0(), EAST)
#aisle_buttress(u0() + vx(DEFAULT_PILLARS_DISTANCE * 2), NORTH)
#basic_pillar(u0() + vx(DEFAULT_PILLARS_DISTANCE * 2) - vy(DEFAULT_PILLARS_DISTANCE))
#aisle_inner_pillar(u0() + vx(DEFAULT_PILLARS_DISTANCE * 2) - vy(DEFAULT_PILLARS_DISTANCE * 2), NORTH)
#ambulatory_buttress(u0() + vx(DEFAULT_PILLARS_DISTANCE * 4), NORTH)
#ambulatory_outer_pillar(u0() + vx(DEFAULT_PILLARS_DISTANCE * 4) - vy(DEFAULT_PILLARS_DISTANCE), NORTH)
#ambulatory_inner_pillar(u0() + vx(DEFAULT_PILLARS_DISTANCE * 4) - vy(DEFAULT_PILLARS_DISTANCE * 2), NORTH)

function get_upper_ambulatory_positional_angle(column)
    return get_ambulatory_start_angle() + get_ambulatory_angle_increment() * (last(get_ambulatory_col_range()) - column + 1)
end

function get_lower_ambulatory_positional_angle(column)
    return -get_upper_ambulatory_positional_angle(column)
end

function get_ambulatory_positional_angle(row, column)
    horizontal_hallway_jump_prev_row = get_horizontal_hallway_jump_prev_row()

    if row > horizontal_hallway_jump_prev_row
        return get_upper_ambulatory_positional_angle(column)
    else
        return get_lower_ambulatory_positional_angle(column)
    end

    return nothing
end

function get_upper_ambulatory_positional_radial_length(row)
    default_pillars_distance = get_default_pillars_distance()
    horizontal_hallway_jump_prev_row = get_horizontal_hallway_jump_prev_row()
    ambulatory_last_row = last(get_ambulatory_row_range())
    positional_radial_length = (row - horizontal_hallway_jump_prev_row) * default_pillars_distance

    if row == ambulatory_last_row
        positional_radial_length -= default_pillars_distance / 2
    end

    return positional_radial_length
end

function get_lower_ambulatory_positional_radial_length(row)
    default_pillars_distance = get_default_pillars_distance()
    horizontal_hallway_jump_prev_row = get_horizontal_hallway_jump_prev_row()
    ambulatory_first_row = first(get_ambulatory_row_range())
    positional_radial_length = (horizontal_hallway_jump_prev_row - row + 1) * default_pillars_distance

    if row == ambulatory_first_row
        positional_radial_length -= default_pillars_distance / 2
    end

    return positional_radial_length
end

function get_ambulatory_positional_radial_length(row)
    horizontal_hallway_jump_prev_row = get_horizontal_hallway_jump_prev_row()

    if row > horizontal_hallway_jump_prev_row
        return get_upper_ambulatory_positional_radial_length(row)
    else
        return get_lower_ambulatory_positional_radial_length(row)
    end

    return nothing
end

function get_ambulatory_pillar_coordinates(row, column)
    return get_ambulatory_center() + vpol(get_ambulatory_positional_radial_length(row), 
                                            get_ambulatory_positional_angle(row, column))
end

function get_pillar_coordinates(row, column)
    default_pillars_distance = get_default_pillars_distance()
    horizontal_hallway_jump_prev_row = get_horizontal_hallway_jump_prev_row()
    vertical_hallway_jump_prev_col = get_vertical_hallway_jump_prev_col()
    ambulatory_first_column = first(get_ambulatory_col_range())

    x = default_pillars_distance * column
    y = default_pillars_distance * row

    if row > horizontal_hallway_jump_prev_row
        y += default_pillars_distance
    end

    if column > vertical_hallway_jump_prev_col && column < ambulatory_first_column
        x += default_pillars_distance
    elseif column >= ambulatory_first_column
        return get_ambulatory_pillar_coordinates(row, column)
    end

    return xy(x, y)
end

function static_row_col_range_pillar_instantiator(pillars, pillar_info_collection, pillar_type_instantiator)
    for pillar_info in pillar_info_collection
        row = get_pillars_row_info(pillar_info)
        col_range = get_pillars_col_info(pillar_info)
        orientation = get_pillars_orientation_info(pillar_info)

        for col in col_range
            center = get_pillar_coordinates(row, col)
            pillar = pillar_type_instantiator(center, orientation)
            pillar_model = pillar.model
            pillar_width = pillar.width
            pillar_depth = pillar.depth
            pillar_height = pillar.height
            pillars[row, col] = Pillar(pillar_model, row, col, center, orientation, pillar_width, pillar_depth, pillar_height)
        end
    end
end

function row_range_static_col_pillar_instantiator(pillars, pillar_info_collection, pillar_type_instantiator)
    for pillar_info in pillar_info_collection
        row_range = get_pillars_row_info(pillar_info)
        col = get_pillars_col_info(pillar_info)
        orientation = get_pillars_orientation_info(pillar_info)

        for row in row_range
            center = get_pillar_coordinates(row, col)
            pillar = pillar_type_instantiator(center, orientation)
            pillar_model = pillar.model
            pillar_width = pillar.width
            pillar_depth = pillar.depth
            pillar_height = pillar.height
            pillars[row, col] = Pillar(pillar_model, row, col, center, orientation, pillar_width, pillar_depth, pillar_height)
        end
    end
end

function static_row_static_col_pillar_instantiator(pillars, pillar_info_collection, pillar_type_instantiator)
    for pillar_info in pillar_info_collection
        row = get_pillars_row_info(pillar_info)
        col = get_pillars_col_info(pillar_info)
        orientation = get_pillars_orientation_info(pillar_info)
        center = get_pillar_coordinates(row, col)
        pillar = pillar_type_instantiator(center, orientation)
        pillar_model = pillar.model
        pillar_width = pillar.width
        pillar_depth = pillar.depth
        pillar_height = pillar.height
        pillars[row, col] = Pillar(pillar_model, row, col, center, orientation, pillar_width, pillar_depth, pillar_height)
    end
end

function row_range_col_range_pillar_instantiator(pillars, pillar_info_collection, pillar_type_instantiator)
    for pillar_info in pillar_info_collection
        row_range = get_pillars_row_info(pillar_info)
        col_range = get_pillars_col_info(pillar_info)
        orientation = get_pillars_orientation_info(pillar_info)

        for row in row_range
            for col in col_range
                center = get_pillar_coordinates(row, col)
                pillar = pillar_type_instantiator(center, orientation)
                pillar_model = pillar.model
                pillar_width = pillar.width
                pillar_depth = pillar.depth
                pillar_height = pillar.height
                pillars[row, col] = Pillar(pillar_model, row, col, center, orientation, pillar_width, pillar_depth, pillar_height)
            end
        end
    end
end

function instantiate_aisles_buttresses(pillars)
    static_row_col_range_pillar_instantiator(pillars, get_aisles_buttresses(), buttress)
end

function instantiate_outer_crossing_buttresses(pillars)
    static_row_static_col_pillar_instantiator(pillars, get_outer_crossing_buttresses(), outer_crossing_buttress)
end

function instantiate_transept_buttresses(pillars)
    row_range_static_col_pillar_instantiator(pillars, get_transept_buttresses(), buttress)
end

function instantiate_aisle_outer_pillars(pillars)
    row_range_col_range_pillar_instantiator(pillars, get_aisle_outer_pillars(), aisle_outer_pillar)
end

function instantiate_aisle_inner_pillars(pillars)
    static_row_col_range_pillar_instantiator(pillars, get_aisle_inner_pillars(), aisle_inner_pillar)
end

function instantiate_transept_pillars(pillars)
    row_range_col_range_pillar_instantiator(pillars, get_transept_pillars(), transept_pillar)
end

function instantiate_crossing_pillars(pillars)
    static_row_static_col_pillar_instantiator(pillars, get_crossing_pillars(), crossing_pillar)
end

function instantiate_ambulatory_buttresses(pillars)
    static_row_col_range_pillar_instantiator(pillars, get_ambulatory_buttresses(), ambulatory_buttress)
end

function instantiate_ambulatory_outer_pillars(pillars)
    row_range_col_range_pillar_instantiator(pillars, get_ambulatory_outer_pillars(), ambulatory_outer_pillar)
end

function instantiate_ambulatory_inner_pillars(pillars)
    static_row_col_range_pillar_instantiator(pillars, get_ambulatory_inner_pillars(), ambulatory_inner_pillar)
end

function instantiate_all_pillars(pillars)
    instantiate_aisles_buttresses(pillars)
    instantiate_outer_crossing_buttresses(pillars)
    instantiate_transept_buttresses(pillars)
    instantiate_aisle_outer_pillars(pillars)
    instantiate_aisle_inner_pillars(pillars)
    instantiate_transept_pillars(pillars)
    instantiate_crossing_pillars(pillars)
    instantiate_ambulatory_buttresses(pillars)
    instantiate_ambulatory_outer_pillars(pillars)
    instantiate_ambulatory_inner_pillars(pillars)
end
# == INSTANTIATORS == #
# == PILLARS == #

# == WALLS INFO == #
## == FLYING BUTTRESSES == #
#NW_FLYING_BUTTRESS_LEFT_PILLARS = (3, range(4, 8))
#NW_FLYING_BUTTRESS_RIGHT_PILLARS = (4, range(4, 8))
#NE_FLYING_BUTTRESS_LEFT_PILLARS = (3, range(11, 14))
#NE_FLYING_BUTTRESS_RIGHT_PILLARS = (4, range(11, 14))
#N_FLYING_BUTTRESSES_LEFT_ANCHOR_HEIGHT = AISLE_BUTTRESS_HEIGHT / 2 + AISLE_BUTTRESS_HEIGHT / 6
#N_FLYING_BUTTRESSES_RIGHT_ANCHOR_HEIGHT = AISLE_OUTER_PILLAR_HEIGHT / 2 + AISLE_OUTER_PILLAR_HEIGHT / 4
#
#SW_FLYING_BUTTRESS_LEFT_PILLARS = (7, range(4, 8))
#SW_FLYING_BUTTRESS_RIGHT_PILLARS = (8, range(4, 8))
#SE_FLYING_BUTTRESS_LEFT_PILLARS = (7, range(11, 14))
#SE_FLYING_BUTTRESS_RIGHT_PILLARS = (8, range(11, 14))
#S_FLYING_BUTTRESSES_LEFT_ANCHOR_HEIGHT = N_FLYING_BUTTRESSES_RIGHT_ANCHOR_HEIGHT
#S_FLYING_BUTTRESSES_RIGHT_ANCHOR_HEIGHT = N_FLYING_BUTTRESSES_LEFT_ANCHOR_HEIGHT
## == FLYING BUTTRESSES == #

# == WALLS == #
WALLS_INFO = Dict{Symbol, Any}()

mutable struct Wall
    model
    left_pillar
    right_pillar
end

function get_wall_left_pillar(wall)
    return wall.left_pillar
end

function get_wall_right_pillar(wall)
    return wall.right_pillar
end

function set_wall_left_pillar(wall, pillar)
    wall.left_pillar = pillar
end

function set_wall_right_pillar(wall, pillar)
    wall.right_pillar = pillar
end

function assign_walls_info_key_value(key, value)
    WALLS_INFO[key] = value
end

function instantiate_walls_info()
    # == RANGES == #
    w_aisle_row_range = get_w_aisle_row_range()
    w_aisle_col_range = get_w_aisle_col_range()
    e_aisle_row_range = get_e_aisle_row_range()
    e_aisle_col_range = get_e_aisle_col_range()

    n_transept_row_range = get_n_transept_row_range()
    n_transept_col_range = get_n_transept_col_range()
    s_transept_row_range = get_s_transept_row_range()
    s_transept_col_range = get_s_transept_col_range()

    ambulatory_row_range = get_ambulatory_row_range()
    ambulatory_col_range = get_ambulatory_col_range()
    # == RANGES == #

    # == AISLE BUTTRESS WALLS == #
    aisle_buttress_wall_height = get_aisle_buttress_height() / 2

    ws_aisle_buttress_walls = (first(w_aisle_row_range), w_aisle_col_range[1:end-1])
    wn_aisle_buttress_walls = (last(w_aisle_row_range), w_aisle_col_range[1:end-1])
    es_aisle_buttress_walls = (first(e_aisle_row_range), e_aisle_col_range[1:end-1])
    en_aisle_buttress_walls = (last(e_aisle_row_range), e_aisle_col_range[1:end-1])

    aisle_buttress_walls = [ws_aisle_buttress_walls, wn_aisle_buttress_walls,
                                es_aisle_buttress_walls, en_aisle_buttress_walls]

    assign_walls_info_key_value(:AISLE_BUTTRESS_WALL_HEIGHT, aisle_buttress_wall_height)
    assign_walls_info_key_value(:AISLE_BUTTRESS_WALLS, aisle_buttress_walls)
    # == AISLE BUTTRESS WALLS == #

    # == TRANSEPT BUTTRESS WALLS == #
    transept_buttress_wall_height = get_transept_buttress_height() / 2

    sw_transept_buttress_walls = (range(first(s_transept_row_range), first(w_aisle_row_range) - 1), last(w_aisle_col_range))
    se_transept_buttress_walls = (range(first(s_transept_row_range), first(e_aisle_row_range) - 1), first(e_aisle_col_range))
    nw_transept_buttress_walls = (range(last(w_aisle_row_range), last(n_transept_row_range) - 1), last(w_aisle_col_range))
    ne_transept_buttress_walls = (range(last(e_aisle_row_range), last(n_transept_row_range) - 1), first(e_aisle_col_range))
    transept_buttress_walls = [sw_transept_buttress_walls, se_transept_buttress_walls,
                                nw_transept_buttress_walls, ne_transept_buttress_walls]

    assign_walls_info_key_value(:TRANSEPT_BUTTRESS_WALL_HEIGHT, transept_buttress_wall_height)
    assign_walls_info_key_value(:TRANSEPT_BUTTRESS_WALLS, transept_buttress_walls)
    # == TRANSEPT BUTTRESS WALLS == #

    # == AMBULATORY BUTTRESS WALLS == #
    ambulatory_buttress_wall_height = aisle_buttress_wall_height

    s_ambulatory_buttress_walls = (first(ambulatory_row_range), ambulatory_col_range)
    n_ambulatory_buttress_walls = (last(ambulatory_row_range), ambulatory_col_range)
    ambulatory_buttress_walls = [s_ambulatory_buttress_walls, n_ambulatory_buttress_walls]

    assign_walls_info_key_value(:AMBULATORY_BUTTRESS_WALL_HEIGHT, ambulatory_buttress_wall_height)
    assign_walls_info_key_value(:AMBULATORY_BUTTRESS_WALLS, ambulatory_buttress_walls)
    # == AMBULATORY BUTTRESS WALLS == #
    
    # == AISLE OUTER PILLAR HORIZONTAL WALLS == #
    aisle_outer_pillar_horizontal_wall_height = aisle_buttress_wall_height

    ws_aisle_outer_pillar_horizontal_walls = (range(w_aisle_row_range[2], get_horizontal_hallway_jump_prev_row() - 1), w_aisle_col_range[1:end-1])
    wn_aisle_outer_pillar_horizontal_walls = (range(get_horizontal_hallway_jump_prev_row() + 2, w_aisle_row_range[end-1]), w_aisle_col_range[1:end-1])
    es_aisle_outer_pillar_horizontal_walls = (range(e_aisle_row_range[2], get_horizontal_hallway_jump_prev_row() - 1), e_aisle_col_range[1:end-1])
    en_aisle_outer_pillar_horizontal_walls = (range(get_horizontal_hallway_jump_prev_row() + 2, e_aisle_row_range[end-1]), e_aisle_col_range[1:end-1])
    aisle_outer_pillar_horizontal_walls = [ws_aisle_outer_pillar_horizontal_walls, wn_aisle_outer_pillar_horizontal_walls,
                                            es_aisle_outer_pillar_horizontal_walls, en_aisle_outer_pillar_horizontal_walls]

    assign_walls_info_key_value(:AISLE_OUTER_PILLAR_HORIZONTAL_WALL_HEIGHT, aisle_outer_pillar_horizontal_wall_height)
    assign_walls_info_key_value(:AISLE_OUTER_PILLAR_HORIZONTAL_WALLS, aisle_outer_pillar_horizontal_walls)
    # == AISLE OUTER PILLAR HORIZONTAL WALLS == #

    # == AISLE OUTER PILLAR VERTICAL WALLS == #
    aisle_outer_pillar_vertical_wall_height = aisle_outer_pillar_horizontal_wall_height

    ws_aisle_outer_pillar_vertical_walls = (range(first(w_aisle_row_range), get_horizontal_hallway_jump_prev_row() - 1), w_aisle_col_range)
    wn_aisle_outer_pillar_vertical_walls = (range(get_horizontal_hallway_jump_prev_row() + 1, w_aisle_row_range[end-1]), w_aisle_col_range)
    es_aisle_outer_pillar_vertical_walls = (range(first(e_aisle_row_range), get_horizontal_hallway_jump_prev_row() - 1), e_aisle_col_range)
    en_aisle_outer_pillar_vertical_walls = (range(get_horizontal_hallway_jump_prev_row() + 1, e_aisle_row_range[end-1]), e_aisle_col_range)
    aisle_outer_pillar_vertical_walls = [ws_aisle_outer_pillar_vertical_walls, wn_aisle_outer_pillar_vertical_walls,
                                            es_aisle_outer_pillar_vertical_walls, en_aisle_outer_pillar_vertical_walls]

    assign_walls_info_key_value(:AISLE_OUTER_PILLAR_VERTICAL_WALL_HEIGHT, aisle_outer_pillar_vertical_wall_height)
    assign_walls_info_key_value(:AISLE_OUTER_PILLAR_VERTICAL_WALLS, aisle_outer_pillar_vertical_walls)
    # == AISLE OUTER PILLAR VERTICAL WALLS == #

    # == AISLE INNER PILLAR HORIZONTAL WALLS == #
    aisle_inner_pillar_horizontal_wall_height = get_aisle_inner_pillar_height()

    ws_aisle_inner_pillar_horizontal_walls = (get_horizontal_hallway_jump_prev_row(), w_aisle_col_range)
    wn_aisle_inner_pillar_horizontal_walls = (get_horizontal_hallway_jump_prev_row() + 1, w_aisle_col_range)
    es_aisle_inner_pillar_horizontal_walls = (get_horizontal_hallway_jump_prev_row(), range(get_vertical_hallway_jump_prev_col() + 1, e_aisle_col_range[end-1]))
    en_aisle_inner_pillar_horizontal_walls = (get_horizontal_hallway_jump_prev_row() + 1, range(get_vertical_hallway_jump_prev_col() + 1, e_aisle_col_range[end-1]))
    aisle_inner_pillar_horizontal_walls = [ws_aisle_inner_pillar_horizontal_walls, wn_aisle_inner_pillar_horizontal_walls,
                                            es_aisle_inner_pillar_horizontal_walls, en_aisle_inner_pillar_horizontal_walls]
    
    assign_walls_info_key_value(:AISLE_INNER_PILLAR_HORIZONTAL_WALL_HEIGHT, aisle_inner_pillar_horizontal_wall_height)
    assign_walls_info_key_value(:AISLE_INNER_PILLAR_HORIZONTAL_WALLS, aisle_inner_pillar_horizontal_walls)
    # == AISLE INNER PILLAR HORIZONTAL WALLS == #

    # == AISLE INNER PILLAR VERTICAL WALLS == #
    aisle_inner_pillar_vertical_wall_height = aisle_inner_pillar_horizontal_wall_height
    
    w_aisle_inner_pillar_vertical_walls = (get_horizontal_hallway_jump_prev_row(), range(first(w_aisle_col_range), get_vertical_hallway_jump_prev_col()))
    e_aisle_inner_pillar_vertical_walls = (get_horizontal_hallway_jump_prev_row(), range(get_vertical_hallway_jump_prev_col() + 1, last(e_aisle_col_range)))
    aisle_inner_pillar_vertical_walls = [w_aisle_inner_pillar_vertical_walls, e_aisle_inner_pillar_vertical_walls]

    assign_walls_info_key_value(:AISLE_INNER_PILLAR_VERTICAL_WALL_HEIGHT, aisle_inner_pillar_vertical_wall_height)
    assign_walls_info_key_value(:AISLE_INNER_PILLAR_VERTICAL_WALLS, aisle_inner_pillar_vertical_walls)
    # == AISLE INNER PILLAR VERTICAL WALLS == #

    # == TRANSEPT OUTER PILLAR HORIZONTAL WALLS == #
    transept_outer_pillar_horizontal_wall_height = aisle_buttress_wall_height

    sw_transept_outer_pillar_horizontal_walls = (s_transept_row_range[1:end-1], range(last(w_aisle_col_range), get_vertical_hallway_jump_prev_col() - 1))
    se_transept_outer_pillar_horizontal_walls = (s_transept_row_range[1:end-1], range(get_vertical_hallway_jump_prev_col() + 1, last(s_transept_col_range)))
    nw_transept_outer_pillar_horizontal_walls = (n_transept_row_range[2:end], range(last(w_aisle_col_range), get_vertical_hallway_jump_prev_col() - 1))
    ne_transept_outer_pillar_horizontal_walls = (n_transept_row_range[2:end], range(get_vertical_hallway_jump_prev_col() + 1, last(s_transept_col_range)))
    transept_outer_pillar_horizontal_walls = [sw_transept_outer_pillar_horizontal_walls, se_transept_outer_pillar_horizontal_walls,
                                                nw_transept_outer_pillar_horizontal_walls, ne_transept_outer_pillar_horizontal_walls]

    assign_walls_info_key_value(:TRANSEPT_OUTER_PILLAR_HORIZONTAL_WALL_HEIGHT, transept_outer_pillar_horizontal_wall_height)
    assign_walls_info_key_value(:TRANSEPT_OUTER_PILLAR_HORIZONTAL_WALLS, transept_outer_pillar_horizontal_walls)
    # == TRANSEPT OUTER PILLAR HORIZONTAL WALLS == #

    # == TRANSEPT OUTER PILLAR VERTICAL WALLS == #
    transept_outer_pillar_vertical_wall_height = transept_outer_pillar_horizontal_wall_height

    sw_transept_outer_pillar_vertical_walls = (s_transept_row_range[1:end-1], range(first(s_transept_col_range), get_vertical_hallway_jump_prev_col() - 1))
    se_transept_outer_pillar_vertical_walls = (s_transept_row_range[1:end-1], range(get_vertical_hallway_jump_prev_col() + 2, last(s_transept_col_range)))
    nw_transept_outer_pillar_vertical_walls = (n_transept_row_range[1:end-1], range(first(n_transept_col_range), get_vertical_hallway_jump_prev_col() - 1))
    ne_transept_outer_pillar_vertical_walls = (n_transept_row_range[1:end-1], range(get_vertical_hallway_jump_prev_col() + 2, last(n_transept_col_range)))
    transept_outer_pillar_vertical_walls = [sw_transept_outer_pillar_vertical_walls, se_transept_outer_pillar_vertical_walls,
                                                nw_transept_outer_pillar_vertical_walls, ne_transept_outer_pillar_vertical_walls]

    assign_walls_info_key_value(:TRANSEPT_OUTER_PILLAR_VERTICAL_WALL_HEIGHT, transept_outer_pillar_vertical_wall_height)
    assign_walls_info_key_value(:TRANSEPT_OUTER_PILLAR_VERTICAL_WALLS, transept_outer_pillar_vertical_walls)
    # == TRANSEPT OUTER PILLAR VERTICAL WALLS == #

    # == TRANSEPT INNER PILLAR HORIZONTAL WALLS == #
    transept_inner_pillar_horizontal_wall_height = aisle_inner_pillar_horizontal_wall_height

    s_transept_inner_pillar_horizontal_walls = (s_transept_row_range, get_vertical_hallway_jump_prev_col())
    n_transept_inner_pillar_horizontal_walls = (n_transept_row_range, get_vertical_hallway_jump_prev_col())
    transept_inner_pillar_horizontal_walls = [s_transept_inner_pillar_horizontal_walls, n_transept_inner_pillar_horizontal_walls]

    assign_walls_info_key_value(:TRANSEPT_INNER_PILLAR_HORIZONTAL_WALL_HEIGHT, transept_inner_pillar_horizontal_wall_height)
    assign_walls_info_key_value(:TRANSEPT_INNER_PILLAR_HORIZONTAL_WALLS, transept_inner_pillar_horizontal_walls)
    # == TRANSEPT INNER PILLAR HORIZONTAL WALLS == #

    # == TRANSEPT INNER PILLAR VERTICAL WALLS == #
    transept_inner_pillar_vertical_wall_height = aisle_inner_pillar_vertical_wall_height

    sw_transept_inner_pillar_vertical_walls = (s_transept_row_range[1:end-1], get_vertical_hallway_jump_prev_col())
    se_transept_inner_pillar_vertical_walls = (s_transept_row_range[1:end-1], get_vertical_hallway_jump_prev_col() + 1)
    nw_transept_inner_pillar_vertical_walls = (n_transept_row_range[1:end-1], get_vertical_hallway_jump_prev_col())
    ne_transept_inner_pillar_vertical_walls = (n_transept_row_range[1:end-1], get_vertical_hallway_jump_prev_col() + 1)
    transept_inner_pillar_vertical_walls = [sw_transept_inner_pillar_vertical_walls, se_transept_inner_pillar_vertical_walls,
                                                nw_transept_inner_pillar_vertical_walls, ne_transept_inner_pillar_vertical_walls]

    assign_walls_info_key_value(:TRANSEPT_INNER_PILLAR_VERTICAL_WALL_HEIGHT, transept_inner_pillar_vertical_wall_height)
    assign_walls_info_key_value(:TRANSEPT_INNER_PILLAR_VERTICAL_WALLS, transept_inner_pillar_vertical_walls)
    # == TRANSEPT INNER PILLAR VERTICAL WALLS == #

    # == AMBULATORY OUTER PILLAR HORIZONTAL WALLS == #
    ambulatory_outer_pillar_horizontal_wall_height = aisle_buttress_wall_height

    s_ambulatory_outer_pillar_horizontal_walls = (range(ambulatory_row_range[2], get_horizontal_hallway_jump_prev_row() - 1), 
                                                    range(last(e_aisle_col_range), ambulatory_col_range[end]))
    n_ambulatory_outer_pillar_horizontal_walls = (range(get_horizontal_hallway_jump_prev_row() + 2, ambulatory_row_range[end-1]), 
                                                    range(last(e_aisle_col_range), ambulatory_col_range[end-1]))
    ambulatory_outer_pillar_horizontal_walls = [s_ambulatory_outer_pillar_horizontal_walls, n_ambulatory_outer_pillar_horizontal_walls]

    assign_walls_info_key_value(:AMBULATORY_OUTER_PILLAR_HORIZONTAL_WALL_HEIGHT, ambulatory_outer_pillar_horizontal_wall_height)
    assign_walls_info_key_value(:AMBULATORY_OUTER_PILLAR_HORIZONTAL_WALLS, ambulatory_outer_pillar_horizontal_walls)
    # == AMBULATORY OUTER PILLAR HORIZONTAL WALLS == #

    # == AMBULATORY OUTER PILLAR VERTICAL WALLS == #
    ambulatory_outer_pillar_vertical_wall_height = ambulatory_outer_pillar_horizontal_wall_height

    s_ambulatory_outer_pillar_vertical_walls = (range(ambulatory_row_range[2], get_horizontal_hallway_jump_prev_row() - 1), 
                                                    range(last(e_aisle_col_range), ambulatory_col_range[end]))
    n_ambulatory_outer_pillar_vertical_walls = (range(get_horizontal_hallway_jump_prev_row() + 1, ambulatory_row_range[end-2]),
                                                    range(last(e_aisle_col_range), ambulatory_col_range[end]))
    ambulatory_outer_pillar_vertical_walls = [s_ambulatory_outer_pillar_vertical_walls, n_ambulatory_outer_pillar_vertical_walls]

    assign_walls_info_key_value(:AMBULATORY_OUTER_PILLAR_VERTICAL_WALL_HEIGHT, ambulatory_outer_pillar_vertical_wall_height)
    assign_walls_info_key_value(:AMBULATORY_OUTER_PILLAR_VERTICAL_WALLS, ambulatory_outer_pillar_vertical_walls)
    # == AMBULATORY OUTER PILLAR VERTICAL WALLS == #

    # == AMBULATORY INNER PILLAR WALLS == #
    ambulatory_inner_pillar_wall_height = aisle_inner_pillar_horizontal_wall_height

    s_ambulatory_inner_pillar_walls = (get_horizontal_hallway_jump_prev_row(), range(last(e_aisle_col_range), ambulatory_col_range[end]))
    n_ambulatory_inner_pillar_walls = (get_horizontal_hallway_jump_prev_row() + 1, range(last(e_aisle_col_range), ambulatory_col_range[end-1]))
    ambulatory_inner_pillar_walls = [s_ambulatory_inner_pillar_walls, n_ambulatory_inner_pillar_walls]

    assign_walls_info_key_value(:AMBULATORY_INNER_PILLAR_WALL_HEIGHT, ambulatory_inner_pillar_wall_height)
    assign_walls_info_key_value(:AMBULATORY_INNER_PILLAR_WALLS, ambulatory_inner_pillar_walls)
    # == AMBULATORY INNER PILLAR WALLS == #
end
# == WALLS == #

# == WALL GETTERS == #
function get_walls_info_value(key)
    return WALLS_INFO[key]
end

function get_aisle_buttress_wall_height()
    return get_walls_info_value(:AISLE_BUTTRESS_WALL_HEIGHT)
end

function get_aisle_buttress_walls()
    return get_walls_info_value(:AISLE_BUTTRESS_WALLS)
end

function get_transept_buttress_wall_height()
    return get_walls_info_value(:TRANSEPT_BUTTRESS_WALL_HEIGHT)
end

function get_transept_buttress_walls()
    return get_walls_info_value(:TRANSEPT_BUTTRESS_WALLS)
end

function get_aisle_outer_pillar_horizontal_wall_height()
    return get_walls_info_value(:AISLE_OUTER_PILLAR_HORIZONTAL_WALL_HEIGHT)
end

function get_aisle_outer_pillar_horizontal_walls()
    return get_walls_info_value(:AISLE_OUTER_PILLAR_HORIZONTAL_WALLS)
end

function get_aisle_outer_pillar_vertical_wall_height()
    return get_walls_info_value(:AISLE_OUTER_PILLAR_VERTICAL_WALL_HEIGHT)
end

function get_aisle_outer_pillar_vertical_walls()
    return get_walls_info_value(:AISLE_OUTER_PILLAR_VERTICAL_WALLS)
end

function get_aisle_inner_pillar_horizontal_wall_height()
    return get_walls_info_value(:AISLE_INNER_PILLAR_HORIZONTAL_WALL_HEIGHT)
end

function get_aisle_inner_pillar_horizontal_walls()
    return get_walls_info_value(:AISLE_INNER_PILLAR_HORIZONTAL_WALLS)
end

function get_aisle_inner_pillar_vertical_wall_height()
    return get_walls_info_value(:AISLE_INNER_PILLAR_VERTICAL_WALL_HEIGHT)
end

function get_aisle_inner_pillar_vertical_walls()
    return get_walls_info_value(:AISLE_INNER_PILLAR_VERTICAL_WALLS)
end

function get_transept_outer_pillar_horizontal_wall_height()
    return get_walls_info_value(:TRANSEPT_OUTER_PILLAR_HORIZONTAL_WALL_HEIGHT)
end

function get_transept_outer_pillar_horizontal_walls()
    return get_walls_info_value(:TRANSEPT_OUTER_PILLAR_HORIZONTAL_WALLS)
end

function get_transept_outer_pillar_vertical_wall_height()
    return get_walls_info_value(:TRANSEPT_OUTER_PILLAR_VERTICAL_WALL_HEIGHT)
end

function get_transept_outer_pillar_vertical_walls()
    return get_walls_info_value(:TRANSEPT_OUTER_PILLAR_VERTICAL_WALLS)
end

function get_transept_inner_pillar_horizontal_wall_height()
    return get_walls_info_value(:TRANSEPT_INNER_PILLAR_HORIZONTAL_WALL_HEIGHT)
end

function get_transept_inner_pillar_horizontal_walls()
    return get_walls_info_value(:TRANSEPT_INNER_PILLAR_HORIZONTAL_WALLS)
end

function get_transept_inner_pillar_vertical_wall_height()
    return get_walls_info_value(:TRANSEPT_INNER_PILLAR_VERTICAL_WALL_HEIGHT)
end

function get_transept_inner_pillar_vertical_walls()
    return get_walls_info_value(:TRANSEPT_INNER_PILLAR_VERTICAL_WALLS)
end

function get_ambulatory_outer_pillar_horizontal_wall_height()
    return get_walls_info_value(:AMBULATORY_OUTER_PILLAR_HORIZONTAL_WALL_HEIGHT)
end

function get_ambulatory_outer_pillar_horizontal_walls()
    return get_walls_info_value(:AMBULATORY_OUTER_PILLAR_HORIZONTAL_WALLS)
end

function get_ambulatory_outer_pillar_vertical_wall_height()
    return get_walls_info_value(:AMBULATORY_OUTER_PILLAR_VERTICAL_WALL_HEIGHT)
end

function get_ambulatory_outer_pillar_vertical_walls()
    return get_walls_info_value(:AMBULATORY_OUTER_PILLAR_VERTICAL_WALLS)
end

function get_ambulatory_inner_pillar_wall_height()
    return get_walls_info_value(:AMBULATORY_INNER_PILLAR_WALL_HEIGHT)
end

function get_ambulatory_inner_pillar_walls()
    return get_walls_info_value(:AMBULATORY_INNER_PILLAR_WALLS)
end
# == WALL GETTERS == #
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

function standing_arch_block(left_pillar, right_pillar, depth, height, arch_excess, offset)
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
    
    # WINDOW STUFF TO BE REVISED OR SIMPLIFIED OR MODULARIZED - WHATEVER'S BEST UNDER TIME CONSTRAINTS
    planar_gothic_window_width = distance(bottom_left_corner, bottom_right_corner)
    planar_gothic_window_height_without_arch = distance(upper_left_corner, bottom_left_corner)

    planar_upper_left_corner = u0() - vx(planar_gothic_window_width / 2)
    planar_upper_right_corner = planar_upper_left_corner + vx(planar_gothic_window_width)
    planar_bottom_left_corner = planar_upper_left_corner - vy(planar_gothic_window_height_without_arch)

    planar_gothic_window = window_style_instantiator(planar_bottom_left_corner, planar_upper_right_corner, 
                                                        arch_excess, vertical_distance_to_sub_arch, 
                                                            depth / 2, depth / 2).three_dimensional_window

    # Such that it is now standing
    planar_gothic_window = rotate(planar_gothic_window, π/2, u0(), EAST)

    # Such that it now shares the same orientation as that of the wall to which it belongs
    wall_base_vector = right_anchor_offset - left_anchor_offset
    window_base_vector = planar_upper_right_corner - planar_upper_left_corner
    planar_gothic_window = rotate(planar_gothic_window, angle_between_ccw_oriented(wall_base_vector, window_base_vector), u0(), DECREASING_HEIGHT_DIRECTION)
    planar_gothic_window = move(planar_gothic_window, intermediate_loc(upper_left_corner, upper_right_corner) - u0())

    return union(window_wall, planar_gothic_window)
end

function standing_arch_wall_window_block(left_pillar, right_pillar, depth, 
                                            arch_block_height, blank_block_height, window_block_height, 
                                                arch_excess, vertical_distance_to_sub_arch, window_style_instantiator, offset)
    wall_anchors = get_wall_anchors(left_pillar, right_pillar)
    left_anchor = wall_anchors.left_anchor
    right_anchor = wall_anchors.right_anchor

    left_to_right_anchor_vector = right_anchor - left_anchor
    left_anchor_offset = displace_point_by_vector(left_anchor, left_to_right_anchor_vector, offset)
    right_anchor_offset = displace_point_by_vector(right_anchor, left_to_right_anchor_vector, -offset)

    arch_block = standing_arch_block(left_pillar, right_pillar, depth, arch_block_height, arch_excess, offset)
    blank_block = standing_wall_block(left_anchor_offset, right_anchor_offset, depth, blank_block_height)
    window_block = standing_window_wall_block(left_pillar, right_pillar, depth, window_block_height, 
                                            arch_excess, vertical_distance_to_sub_arch, window_style_instantiator, offset)
    blank_block = move(blank_block, GROWING_HEIGHT_DIRECTION * arch_block_height)
    window_block = move(window_block, GROWING_HEIGHT_DIRECTION * (arch_block_height + blank_block_height))

    return union(arch_block, blank_block, window_block)
end

function aisle_buttress_wall(left_pillar, right_pillar;
                                depth = get_pillar_depth(left_pillar),
                                height = get_aisle_buttress_wall_height(),
                                excess = 1,
                                vertical_distance_to_sub_arch = 1,
                                window_style_instantiator = Gothic_Window_First_Style,
                                offset = 0)
    standing_window_wall_block(left_pillar, right_pillar, depth, height,
                                excess, vertical_distance_to_sub_arch, window_style_instantiator, offset)
end

function transept_buttress_wall(left_pillar, right_pillar;
                                    depth = get_pillar_depth(left_pillar),
                                    height = get_transept_buttress_wall_height(),
                                    excess = 1,
                                    vertical_distance_to_sub_arch = 1,
                                    window_style_instantiator = Gothic_Window_First_Style,
                                    offset = 0)
    standing_window_wall_block(left_pillar, right_pillar, depth, height,
                                excess, vertical_distance_to_sub_arch, window_style_instantiator, offset)
end

function aisle_outer_pillar_horizontal_wall(left_pillar, right_pillar;
                                                depth = get_pillar_width(left_pillar),
                                                height = get_aisle_outer_pillar_horizontal_wall_height(),
                                                excess = 1,
                                                offset = 0)
    standing_arch_block(left_pillar, right_pillar, depth, height, excess, offset)
end

function aisle_outer_pillar_vertical_wall(left_pillar, right_pillar;
                                            depth = get_pillar_depth(left_pillar),
                                            height = get_aisle_outer_pillar_vertical_wall_height(),
                                            excess = 1,
                                            offset = 0)
    standing_arch_block(left_pillar, right_pillar, depth, height, excess, offset)
end

function aisle_inner_pillar_horizontal_wall(left_pillar, right_pillar;
                                                depth = get_pillar_width(left_pillar),
                                                height = get_aisle_inner_pillar_horizontal_wall_height(),
                                                excess = 1,
                                                vertical_distance_to_sub_arch = 1,
                                                window_style_instantiator = Gothic_Window_First_Style,
                                                offset = 0)
    arch_block_height = height * 0.4318
    blank_block_height = (height - arch_block_height) / 3
    window_block_height = height - arch_block_height - blank_block_height
    standing_arch_wall_window_block(left_pillar, right_pillar, depth, 
                                        arch_block_height, blank_block_height, window_block_height,
                                            excess, vertical_distance_to_sub_arch, window_style_instantiator, offset)
end

function aisle_inner_pillar_vertical_wall(left_pillar, right_pillar;
                                            depth = get_pillar_depth(left_pillar),
                                            height = get_aisle_inner_pillar_vertical_wall_height(),
                                            excess = 1,
                                            offset = 0)
    standing_arch_block(left_pillar, right_pillar, depth, height, excess, offset)
end

function transept_outer_pillar_horizontal_wall(left_pillar, right_pillar;
                                                depth = get_pillar_depth(left_pillar),
                                                height = get_transept_outer_pillar_horizontal_wall_height(),
                                                excess = 1,
                                                offset = 0)
    standing_arch_block(left_pillar, right_pillar, depth, height, excess, offset)
end

function transept_outer_pillar_vertical_wall(left_pillar, right_pillar;
                                                depth = get_pillar_depth(left_pillar),
                                                height = get_transept_outer_pillar_vertical_wall_height(),
                                                excess = 1,
                                                offset = 0)
    standing_arch_block(left_pillar, right_pillar, depth, height, excess, offset)
end

function transept_inner_pillar_horizontal_wall(left_pillar, right_pillar;
                                                depth = get_pillar_depth(left_pillar),
                                                height = get_aisle_inner_pillar_vertical_wall_height(),
                                                excess = 1,
                                                offset = 0)
    standing_arch_block(left_pillar, right_pillar, depth, height, excess, offset)
end

function transept_inner_pillar_vertical_wall(left_pillar, right_pillar;
                                                depth = get_pillar_width(left_pillar),
                                                height = get_aisle_inner_pillar_horizontal_wall_height(),
                                                excess = 1,
                                                vertical_distance_to_sub_arch = 1,
                                                window_style_instantiator = Gothic_Window_First_Style,
                                                offset = 0)
    arch_block_height = height * 0.4318
    blank_block_height = (height - arch_block_height) / 3
    window_block_height = height - arch_block_height - blank_block_height
    standing_arch_wall_window_block(left_pillar, right_pillar, depth, 
                                        arch_block_height, blank_block_height, window_block_height,
                                            excess, vertical_distance_to_sub_arch, window_style_instantiator, offset)
end

function ambulatory_outer_pillar_horizontal_wall(left_pillar, right_pillar;
                                                    depth = get_pillar_depth(left_pillar),
                                                    height = get_ambulatory_outer_pillar_horizontal_wall_height(),
                                                    excess = 1,
                                                    offset = 0)
    standing_arch_block(left_pillar, right_pillar, depth, height, excess, offset)
end

function ambulatory_outer_pillar_vertical_wall(left_pillar, right_pillar;
                                                    depth = get_pillar_depth(left_pillar),
                                                    height = get_ambulatory_outer_pillar_vertical_wall_height(),
                                                    excess = 1,
                                                    offset = 0)
    standing_arch_block(left_pillar, right_pillar, depth, height, excess, offset)
end

function ambulatory_inner_pillar_wall(left_pillar, right_pillar;
                                        depth = get_pillar_width(left_pillar),
                                        height = get_ambulatory_inner_pillar_wall_height(),
                                        excess = 1,
                                        vertical_distance_to_sub_arch = 1,
                                        window_style_instantiator = Gothic_Window_First_Style,
                                        offset = 0)
    arch_block_height = height * 0.4318
    blank_block_height = (height - arch_block_height) / 3
    window_block_height = height - arch_block_height - blank_block_height
    standing_arch_wall_window_block(left_pillar, right_pillar, depth, 
                                        arch_block_height, blank_block_height, window_block_height,
                                            excess, vertical_distance_to_sub_arch, window_style_instantiator, offset)
end
# == MODELERS == #

# == INSTANTIATORS == #
function ambulatory_outer_arch()
    return nothing
end

function ambulatory_inner_window()
    return nothing
end

function ambulatory_middle_arch()
    return nothing
end

function set_pillars_wall_attributes(left_pillar, right_pillar, wall, wall_type_instantiator)
    left_pillar_center = get_pillar_center(left_pillar)
    right_pillar_center = get_pillar_center(right_pillar)

    direction = right_pillar_center - left_pillar_center
    normalized_direction = normalize_vector(direction)

    if wall_type_instantiator === ambulatory_outer_arch || wall_type_instantiator === ambulatory_inner_window
        set_pillar_east_wall(left_pillar, wall)
        set_pillar_west_wall(right_pillar, wall)
    elseif wall_type_instantiator === ambulatory_middle_arch
        set_pillar_north_wall(left_pillar, wall)
        set_pillar_south_wall(right_pillar, wall)
    elseif isapprox(normalized_direction.x, NORTH.x) && isapprox(normalized_direction.y, NORTH.y)
        set_pillar_north_wall(left_pillar, wall)
        set_pillar_south_wall(right_pillar, wall)
    elseif isapprox(normalized_direction.x, EAST.x) && isapprox(normalized_direction.y, EAST.y)
        set_pillar_east_wall(left_pillar, wall)
        set_pillar_west_wall(right_pillar, wall)
    end
end

function static_row_col_range_wall_instantiator(pillars, wall_info_collection, wall_type_instantiator, wall_orientation)
    for wall_info in wall_info_collection
        row = wall_info[ROW_INDEX]
        col_range = wall_info[COLUMN_INDEX]
        for col in col_range
            left_pillar = pillars[row, col]

            if wall_orientation == HORIZONTAL_WALL
                if col + 1 <= last(size(pillars))
                    right_pillar = pillars[row, col+1]
                else
                    row_delta = get_horizontal_hallway_jump_prev_row() - row
                    right_pillar = pillars[get_horizontal_hallway_jump_prev_row() + row_delta + 1, col]
                end
            else
                right_pillar = pillars[row+1, col]
            end

            wall_model = wall_type_instantiator(left_pillar, right_pillar)
            wall = Wall(wall_model, left_pillar, right_pillar)
            set_pillars_wall_attributes(left_pillar, right_pillar, wall, wall_type_instantiator)
        end
    end
end

function row_range_col_range_wall_instantiator(pillars, wall_info_collection, wall_type_instantiator, wall_orientation)
    for wall_info in wall_info_collection
        row_range = wall_info[ROW_INDEX]
        col_range = wall_info[COLUMN_INDEX]
        for row in row_range
            for col in col_range
                left_pillar = pillars[row, col]

                if wall_orientation == HORIZONTAL_WALL
                    if col + 1 <= last(size(pillars))
                        right_pillar = pillars[row, col+1]
                    else
                        row_delta = get_horizontal_hallway_jump_prev_row() - row
                        right_pillar = pillars[get_horizontal_hallway_jump_prev_row() + row_delta + 1, col]
                    end
                else
                    right_pillar = pillars[row+1, col]
                end

                wall_model = wall_type_instantiator(left_pillar, right_pillar)
                wall = Wall(wall_model, left_pillar, right_pillar)
                set_pillars_wall_attributes(left_pillar, right_pillar, wall, wall_type_instantiator)
            end
        end
    end
end

function row_range_static_col_wall_instantiator(pillars, wall_info_collection, wall_type_instantiator, wall_orientation)
    for wall_info in wall_info_collection
        row_range = wall_info[ROW_INDEX]
        col = wall_info[COLUMN_INDEX]
        for row in row_range
            left_pillar = pillars[row, col]
            right_pillar = wall_orientation == HORIZONTAL_WALL ? pillars[row, col+1] : pillars[row+1, col]
            wall_model = wall_type_instantiator(left_pillar, right_pillar)
            wall = Wall(wall_model, left_pillar, right_pillar)
            set_pillars_wall_attributes(left_pillar, right_pillar, wall, wall_type_instantiator)
        end
    end
end

function instantiate_aisle_buttress_walls(pillars)
    static_row_col_range_wall_instantiator(pillars, get_aisle_buttress_walls(), aisle_buttress_wall, HORIZONTAL_WALL)
end

function instantiate_transept_buttress_walls(pillars)
    row_range_static_col_wall_instantiator(pillars, get_transept_buttress_walls(), transept_buttress_wall, VERTICAL_WALL)
end

function instantiate_aisle_outer_pillar_horizontal_walls(pillars)
    row_range_col_range_wall_instantiator(pillars, get_aisle_outer_pillar_horizontal_walls(), aisle_outer_pillar_horizontal_wall, HORIZONTAL_WALL)
end

function instantiate_aisle_outer_pillar_vertical_walls(pillars)
    row_range_col_range_wall_instantiator(pillars, get_aisle_outer_pillar_vertical_walls(), aisle_outer_pillar_vertical_wall, VERTICAL_WALL)
end

function instantiate_aisle_inner_pillar_horizontal_walls(pillars)
    static_row_col_range_wall_instantiator(pillars, get_aisle_inner_pillar_horizontal_walls(), aisle_inner_pillar_horizontal_wall, HORIZONTAL_WALL)
end

function instantiate_aisle_inner_pillar_vertical_walls(pillars)
    static_row_col_range_wall_instantiator(pillars, get_aisle_inner_pillar_vertical_walls(), aisle_inner_pillar_vertical_wall, VERTICAL_WALL)
end

function instantiate_transept_outer_pillar_horizontal_walls(pillars)
    row_range_col_range_wall_instantiator(pillars, get_transept_outer_pillar_horizontal_walls(), transept_outer_pillar_horizontal_wall, HORIZONTAL_WALL)
end

function instantiate_transept_outer_pillar_vertical_walls(pillars)
    row_range_col_range_wall_instantiator(pillars, get_transept_outer_pillar_vertical_walls(), transept_outer_pillar_vertical_wall, VERTICAL_WALL)
end

function instantiate_transept_inner_pillar_horizontal_walls(pillars)
    row_range_static_col_wall_instantiator(pillars, get_transept_inner_pillar_horizontal_walls(), transept_inner_pillar_horizontal_wall, HORIZONTAL_WALL)
end

function instantiate_transept_inner_pillar_vertical_walls(pillars)
    row_range_static_col_wall_instantiator(pillars, get_transept_inner_pillar_vertical_walls(), transept_inner_pillar_vertical_wall, VERTICAL_WALL)
end

function instantiate_ambulatory_outer_pillar_horizontal_walls(pillars)
    row_range_col_range_wall_instantiator(pillars, get_ambulatory_outer_pillar_horizontal_walls(), ambulatory_outer_pillar_horizontal_wall, HORIZONTAL_WALL)
end

function instantiate_ambulatory_outer_pillar_vertical_walls(pillars)
    row_range_col_range_wall_instantiator(pillars, get_ambulatory_outer_pillar_vertical_walls(), ambulatory_outer_pillar_vertical_wall, VERTICAL_WALL)
end

function instantiate_ambulatory_inner_pillar_walls(pillars)
    static_row_col_range_wall_instantiator(pillars, get_ambulatory_inner_pillar_walls(), ambulatory_inner_pillar_wall, HORIZONTAL_WALL)
end

function instantiate_all_walls(pillars)
    #instantiate_aisle_buttress_walls(pillars)
    #instantiate_transept_buttress_walls(pillars)
    #instantiate_aisle_outer_pillar_horizontal_walls(pillars)
    #instantiate_aisle_outer_pillar_vertical_walls(pillars)
    #instantiate_aisle_inner_pillar_horizontal_walls(pillars)
    #instantiate_aisle_inner_pillar_vertical_walls(pillars)
    #instantiate_transept_outer_pillar_horizontal_walls(pillars)
    #instantiate_transept_outer_pillar_vertical_walls(pillars)
    #instantiate_transept_inner_pillar_horizontal_walls(pillars)
    #instantiate_transept_inner_pillar_vertical_walls(pillars)
    instantiate_ambulatory_outer_pillar_horizontal_walls(pillars)
    instantiate_ambulatory_outer_pillar_vertical_walls(pillars)
    instantiate_ambulatory_inner_pillar_walls(pillars)
end
# == INSTANTIATORS == #
# == WALLS == #

# == PLAYGROUND == #
# == GML CATHEDRAL == #
instantiate_measurements_and_delimiters(7.53 * 2,
                                            3:8, 1:8,
                                            3:8, 11:14,
                                            1:5, 9:10,
                                            6:10, 9:10,
                                            3:8, 15:17)
instantiate_pillars_info(1, 1, 105)
pillars = Array{Union{Pillar, Nothing}}(nothing, 10, 17)
# == GML CATHEDRAL == #

# == AMIENS, CATHEDRALE NOTRE-DAME == #
#instantiate_measurements_and_delimiters(7.53 * 2,
#                                            3:6, 1:7,
#                                            2:7, 10:13,
#                                            1:4, 8:9,
#                                            5:8, 8:9,
#                                            2:7, 14:16)
#instantiate_pillars_info(1, 1, 105)
#pillars = Array{Union{Pillar, Nothing}}(nothing, 8, 16)
# == AMIENS, CATHEDRALE NOTRE-DAME == #

# == REIMS CATHEDRALE NOTRE-DAME == #
#instantiate_measurements_and_delimiters(7.53 * 2,
#                                            2:5, 1:9,
#                                            1:6, 12:14,
#                                            1:3, 10:11,
#                                            4:6, 10:11,
#                                            1:6, 15:16)
#instantiate_pillars_info(1, 1, 105)
#pillars = Array{Union{Pillar, Nothing}}(nothing, 6, 16)
# == REIMS CATHEDRALE NOTRE-DAME == #

# == PARIS CATHEDRALE NOTRE-DAME == #
#instantiate_measurements_and_delimiters(7.53 * 2,
#                                            1:6, 1:8,
#                                            1:6, 11:15,
#                                            1:3, 9:10,
#                                            4:6, 9:10,
#                                            1:6, 16:17)
#instantiate_pillars_info(1, 1, 105)
#pillars = Array{Union{Pillar, Nothing}}(nothing, 6, 17)
# == PARIS CATHEDRALE NOTRE-DAME == #

function measure_time(f, args...)
    @time f(args...)
end

measure_time(instantiate_all_pillars, pillars)

instantiate_walls_info()
measure_time(instantiate_all_walls, pillars)
# == PLAYGROUND == #
