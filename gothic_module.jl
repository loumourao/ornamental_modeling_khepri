using KhepriAutoCAD
const K = KhepriAutoCAD

delete_all_shapes()

# == Utility Functions == #
function displace_point_by_vector(point, vector, magnitude)
    unit_vector = vector / norm(vector)
    return point + unit_vector * magnitude
end

function arc(center, starting_point, ending_point)
    ca = starting_point - center
    cb = ending_point - center

    ccw_starting_angle = atan(ca.y, ca.x)
    ccw_ending_angle = atan(cb.y, cb.x)
    radius = norm(ca)

    if ccw_ending_angle >= ccw_starting_angle
        amplitude = ccw_ending_angle - ccw_starting_angle
    else
        amplitude = 2π - (ccw_starting_angle - ccw_ending_angle)
    end

    return K.arc(center, radius, ccw_starting_angle, amplitude)
end

# The functions below might very well (and should be) removed or, at least, replaced by better suited alternatives
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

function intersect_circles(m0, r0, m1, r1)
    d = distance(m0, m1)
    
    # IMPORTANT: d == 0 prevents the set of interesection points when both circles overlap
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

    return first_intersection_point.y >= second_intersection_point.y ? [first_intersection_point, second_intersection_point] : [second_intersection_point, first_intersection_point]
end

function intersect_circles(c1, c2)
    m0 = circle_center(c1)
    r0 = circle_radius(c1)
    m1 = circle_center(c2)
    r1 = circle_center(c2)

    d = distance(m0, m1)

    # IMPORTANT: d == 0 prevents the set of interesection points when both circles overlap
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

    return first_intersection_point.y >= second_intersection_point.y ? [first_intersection_point, second_intersection_point] : [second_intersection_point, first_intersection_point]
end

function uniformly_extended_inner_offset(offset_value, center, starting_point, ending_point)
    ca = starting_point - center
    cb = ending_point - center
    radius = norm(ca)
    starting_angle = atan(ca.y, ca.x)
    ending_angle = atan(cb.y, cb.x)

    if ending_angle >= starting_angle
        amplitude = ending_angle - starting_angle
    else
        amplitude = 2π - (starting_angle - ending_angle)
    end

    scaling_factor = (radius - offset_value) / radius
    scaled_amplitude = amplitude * scaling_factor
    
    new_radius = radius - offset_value
    new_amplitude = amplitude + (amplitude - scaled_amplitude)
    new_starting_angle = starting_angle - ((amplitude - scaled_amplitude) / 2)

    return K.arc(center, new_radius, new_starting_angle, new_amplitude)
end
# == Utility Functions == #

# Check w/ Prof. what he thinks about this - reserve for final implementation stage tweaks only
# == Predicate Functions for Exception Handling == #
#function is_point(a)
#    return typeof(a) == K::XYZ
#end
#
#function is_number(a)
#    return typeof(a) == Core::Real
#end
#
#function is_point_to_the_left(a, b)
#    return a.x < b.x 
#end
#
#function has_equal_y(a, b)
#    return a.y == b.y
#end
#
#function is_arc_construction_possible(left_point, right_point, excess)
#    return is_point(left_point) && is_point(right_point) && is_number(excess) && excess >= 0.5 &&
#            is_point_to_the_left(left_point, right_point) && has_equal_y(left_point, right_point)
#end
# == Predicate Functions for Exception Handling == #

# == Ornamentation Functions == #
# == Arch Tops == #
function lancet_arch_top(left_point, right_point, excess)
    # Intersection should be performed by CGAL
    # As an alternative, one can extend Khepri with a new algorithm that captures circle intersections
    # The best case for this algorithm, though, is to find a way to automatically get the y of the intersection point (x is the midpoint of ab)
    # this way we would get automatic information about the final height of the window in relation to its excess + bottom "body" height
    # Search for alternative in the event that excess == 0.5

    #if !is_arc_construction_possible(left_point, right_point, excess)
    #    throw(ArgumentError("Arc construction is not possible with the arguments provided"))
    #end

    arch_midpoint = intermediate_loc(left_point, right_point)
    excess_displacement_vector = right_point - left_point
    arcs_radius = norm(excess_displacement_vector) * excess

    right_arc_center = displace_point_by_vector(right_point, excess_displacement_vector, -arcs_radius)
    left_arc_center = displace_point_by_vector(left_point, excess_displacement_vector, arcs_radius)
    
    # This way we can avoid computationally complex operations such as circle_intersection
    # But we should be able to further abstract this concept by resorting to CGAL to be as close as possible to constructive methods
    arc_intersection_x = arch_midpoint.x
    arc_intersection_y = sqrt(arcs_radius^2 - distance(left_arc_center, arch_midpoint)^2) + arch_midpoint.y
    arc_intersection = xy(arc_intersection_x, arc_intersection_y)

    arc(right_arc_center, right_point, arc_intersection)
    arc(left_arc_center, arc_intersection, left_point)

    return (right_arc_center = right_arc_center, left_arc_center = left_arc_center, arcs_radius = arcs_radius)
end

#function pointed_trefoil_arch_top(left_point, right_point, excess)
#    
#end
# == Arch Tops == #

# == Rosette == #
# The code below begs for changes that conform to those mentioned in the corresponding labeled utility functions
function compute_rosette(main_arch_right_arc_center, main_arch_right_arc_radius, 
                    left_sub_arch_right_arc_center, left_sub_arch_right_arc_radius, 
                        outer_offset, inner_offset, vertical_axis, vertical_axis_point)
    ellipse_center = intermediate_loc(left_sub_arch_right_arc_center, main_arch_right_arc_center)
    ellipse_radius = (main_arch_right_arc_radius - outer_offset + left_sub_arch_right_arc_radius) / 2

    rosette_center = intersect_line_ellipse(vertical_axis_point, vertical_axis, ellipse_center, ellipse_radius)[2]
    rosette_radius = distance(rosette_center, left_sub_arch_right_arc_center) - left_sub_arch_right_arc_radius - inner_offset

    circle(rosette_center, rosette_radius)

    return (rosette_center = rosette_center, rosette_radius = rosette_radius)
end

function rosette_rounded_foils(center, radius, n_foils, orientation, outer_offset, inner_offset)
    # Check if n_foils >= 1, otherwise raise an error
    Δα = 2π / n_foils
    foil_radius = (radius * sin(Δα/2)) / (1 + sin(Δα/2))
    center_to_foil_center_length = radius - foil_radius

    foil_center = center + vpol(center_to_foil_center_length, orientation)
    starting_foil_point = center + vpol(center_to_foil_center_length * cos(Δα/2), orientation - (Δα/2))
    ending_foil_point = center + vpol(center_to_foil_center_length * cos(Δα/2), orientation + (Δα/2))

    current_rotation_angle = 0

    fillet_right_point = intersect_circles(center, radius - outer_offset, center + vpol(center_to_foil_center_length, 0), foil_radius)[1]
    fillet_left_point = intersect_circles(center, radius - outer_offset, center + vpol(center_to_foil_center_length, Δα), foil_radius)[2]
    displacement_vector = (center + vpol(center_to_foil_center_length, Δα)) - (center + vpol(center_to_foil_center_length, 0))
    fillet_bottom_point = displace_point_by_vector(center + vpol(center_to_foil_center_length, 0), displacement_vector, norm(displacement_vector) / 2)
    circle(center, radius)
    while n_foils > 0
        rotate(arc(center + vpol(center_to_foil_center_length, 0), fillet_right_point, fillet_bottom_point), current_rotation_angle + orientation, center)
        rotate(arc(center, fillet_right_point, fillet_left_point), current_rotation_angle + orientation, center)
        rotate(arc(center + vpol(center_to_foil_center_length, Δα), fillet_bottom_point, fillet_left_point), current_rotation_angle + orientation, center)
        rotate(uniformly_extended_inner_offset(inner_offset, foil_center, starting_foil_point, ending_foil_point), current_rotation_angle, center)

        current_rotation_angle += Δα
        n_foils -= 1
    end
end

function rosette_pointed_foils(center, radius, n_foils, displacement_ratio, orientation)
    # Check if n_foils >= 1 and if displacement_ratio > 1, otherwise raise an error
    Δα = 2π / n_foils
    rounded_foil_radius = (radius * sin(Δα/2)) / (1 + sin(Δα/2))
    center_to_foil_center_length = radius - rounded_foil_radius

    rounded_foil_center = center + vpol(center_to_foil_center_length, orientation)
    starting_foil_point = center + vpol(center_to_foil_center_length * cos(Δα/2), orientation - (Δα/2))
    ending_foil_point = center + vpol(center_to_foil_center_length * cos(Δα/2), orientation + (Δα/2))

    right_arc_center_displacement_vector = rounded_foil_center - starting_foil_point
    left_arc_center_displacement_vector = rounded_foil_center - ending_foil_point
    displacement_vectors_magnitude = norm(right_arc_center_displacement_vector) * displacement_ratio

    right_arc_center = displace_point_by_vector(starting_foil_point, right_arc_center_displacement_vector, displacement_vectors_magnitude)
    left_arc_center = displace_point_by_vector(ending_foil_point, left_arc_center_displacement_vector, displacement_vectors_magnitude)

    pointed_foil_radius = distance(right_arc_center, starting_foil_point)
    B = distance(right_arc_center, left_arc_center) / 2
    Z = intermediate_loc(right_arc_center, left_arc_center)

    rounded_foil_center_to_arc_intersection_vector = Z - rounded_foil_center
    rounded_foil_center_to_arc_intersection_vector_magnitude = sqrt(pointed_foil_radius^2 - B^2)
    arc_intersection = displace_point_by_vector(Z, rounded_foil_center_to_arc_intersection_vector, rounded_foil_center_to_arc_intersection_vector_magnitude)
    
    scaling_factor = radius / distance(center, arc_intersection)
    current_rotation_angle = 0

    while n_foils > 0
        scale(rotate(arc(right_arc_center, starting_foil_point, arc_intersection), current_rotation_angle, center), scaling_factor, center)
        scale(rotate(arc(left_arc_center, arc_intersection, ending_foil_point), current_rotation_angle, center), scaling_factor, center)

        n_foils -= 1
        current_rotation_angle += Δα
    end
end
# == Rosette == #

# == Fillets == #
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
#
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
# == Fillets == #
# == Ornamentation Functions == #

# == Main Arch == #
function arch(bottom_left_corner, upper_right_corner, excess, outer_offset, inner_offset, recursion_level, vertical_distance_to_sub_arches)
    # Arch body coordinates
    upper_left_corner = xy(bottom_left_corner.x, upper_right_corner.y)
    bottom_right_corner = xy(upper_right_corner.x, bottom_left_corner.y)

    # Arch top portion
    arcs = lancet_arch_top(upper_left_corner, upper_right_corner, excess)
    right_arc_center = arcs.right_arc_center
    left_arc_center = arcs.left_arc_center
    arcs_radius = arcs.arcs_radius

    # Arch body
    line(upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner)

    # 3D profile sweeps
    #slice(swept_surface(arcs[1], circular_path(u0(), 1)), u0(), vx(-1))
    #slice(swept_surface(arcs[2], circular_path(u0(), 1)), u0(), vx(1))
    #swept_surface(polygon, circular_path(u0(), 1))
    
    # Sub-arch calculations
    if recursion_level > 0
        # Offset values
        inner_offset_displacement = inner_offset / 2

        # Sub-arches coordinate auxiliary variables
        left_sub_arch_bottom_left_corner = bottom_left_corner + vxy(outer_offset, outer_offset)
        right_sub_arch_upper_right_corner = upper_right_corner - vxy(outer_offset, vertical_distance_to_sub_arches)
        offset_width_midpoint = intermediate_loc(left_sub_arch_bottom_left_corner, right_sub_arch_upper_right_corner)
        left_sub_arch_upper_right_corner = xy(offset_width_midpoint.x - inner_offset_displacement, right_sub_arch_upper_right_corner.y)
        right_sub_arch_bottom_left_corner = xy(offset_width_midpoint.x + inner_offset_displacement, left_sub_arch_bottom_left_corner.y)

        # TODO: Find a clever way to calculate the new excess
        # For now, the excess will be the same
        sub_arches_excess = excess

        # Sub-arches
        left_sub_arch = arch(left_sub_arch_bottom_left_corner, left_sub_arch_upper_right_corner, sub_arches_excess, outer_offset, inner_offset, recursion_level - 1, vertical_distance_to_sub_arches)
        right_sub_arch = arch(right_sub_arch_bottom_left_corner, right_sub_arch_upper_right_corner, sub_arches_excess, outer_offset, inner_offset, recursion_level - 1, vertical_distance_to_sub_arches)

        # Auxiliary parameters for rosette
        vertical_axis_point = intermediate_loc(left_sub_arch_bottom_left_corner, right_sub_arch_upper_right_corner)
        vertical_axis = intermediate_loc(bottom_left_corner, upper_right_corner) - vertical_axis_point

        # Rosette
        rosette = compute_rosette(right_arc_center, arcs_radius, left_sub_arch.right_arc_center, left_sub_arch.arcs_radius, 
                                    outer_offset, inner_offset, vertical_axis, vertical_axis_point)
        rosette_center = rosette.rosette_center
        rosette_radius = rosette.rosette_radius

        rosette_rounded_foils(rosette_center, rosette_radius, 6, π/2, 0.3, 0.3)
    end

    if recursion_level >= 1
        # Fillets
        circular_rosette_fillets(right_arc_center, left_arc_center, arcs_radius - outer_offset, 
                            rosette_center, rosette_radius + inner_offset, 
                                right_sub_arch.right_arc_center, right_sub_arch.left_arc_center,
                                    left_sub_arch.left_arc_center, left_sub_arch.right_arc_center, left_sub_arch.arcs_radius + inner_offset)
    end

    return (right_arc_center = right_arc_center, left_arc_center = left_arc_center, arcs_radius = arcs_radius)
end
# == Main Arch == #

#with(current_cs, cs_from_o_vz(u0(), vx())) do
#    arch(xy(-10, -16), xy(10, 16), 1, 0.75, 0.75, 1, 3)
#end

arch(xy(-10, -16), xy(10, 16), 1, 0.75, 0.75, 2, 3)

# Instead of resorting to simple arrays to store information about arcs, circles, polygons (as GML does) we resort to Khepri's objects and types
# in such a way that allows us to perform type-checking operations as well as to better aid in the understanding of the program on the user's end
