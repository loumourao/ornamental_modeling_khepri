# TODO: Solve corner joining issues (look into PostScript's setlinejoin function - learn how it functions)
# and look at Sven Havemann's way of handling this issue (profile plane is orthogonal to the tangent of the curve)
# at corner points, the sweep is continued onto the bisector plane - only in cases where the curve is locally symmetric to the bisector plane
# in more general cases on should follow along the medial axis of the two parts of the curve

using KhepriAutoCAD
const K = KhepriAutoCAD
delete_all_shapes()

# == Utility Functions == #
function displace_point_by_vector(point, vector, magnitude)
    unit_vector = vector / K.norm(vector)
    return point + unit_vector * magnitude
end

function midpoint(a, b)
    x = (a.x + b.x) / 2
    y = (a.y + b.y) / 2
    z = (a.z + b.z) / 2

    return K.xyz(x, y, z)
end

function arc(center, starting_point, ending_point)
    ca = starting_point - center
    cb = ending_point - center

    ccw_starting_angle = K.atan(ca.y, ca.x)
    ccw_ending_angle = K.atan(cb.y, cb.x)

    if ccw_ending_angle >= ccw_starting_angle
        amplitude = ccw_ending_angle - ccw_starting_angle
    else
        amplitude = 2π - (ccw_starting_angle - ccw_ending_angle)
    end

    radius = K.norm(ca)

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
    a = K.norm(v)^2
    b = 2 * K.dot((p0 - m), v)
    c = K.norm(p0 - m)^2 - r^2

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

# == Rosette == #
function rosette_rounded_foils(center, radius, n_foils, orientation)
    # Check if n_foils >= 1, otherwise raise an error
    Δα = 2π / n_foils
    foil_radius = (radius * sin(Δα/2)) / (1 + sin(Δα/2))
    center_to_foil_center_length = radius - foil_radius

    foil_center = center + K.vpol(center_to_foil_center_length, orientation)
    starting_foil_point = center + K.vpol(center_to_foil_center_length * cos(Δα/2), orientation - (Δα/2))
    ending_foil_point = center + K.vpol(center_to_foil_center_length * cos(Δα/2), orientation + (Δα/2))

    current_rotation_angle = 0
    
    circle(center, radius)

    while n_foils > 0
        K.rotate(arc(foil_center, starting_foil_point, ending_foil_point), current_rotation_angle, center)

        current_rotation_angle += Δα
        n_foils -= 1
    end
end

function rosette_pointed_foils(center, radius, n_foils, displacement_ratio, orientation)
    # Check if n_foils >= 1 and if displacement_ratio > 1, otherwise raise an error
    Δα = 2π / n_foils
    rounded_foil_radius = (radius * sin(Δα/2)) / (1 + sin(Δα/2))
    center_to_foil_center_length = radius - rounded_foil_radius

    rounded_foil_center = center + K.vpol(center_to_foil_center_length, orientation)
    starting_foil_point = center + K.vpol(center_to_foil_center_length * cos(Δα/2), orientation - (Δα/2))
    ending_foil_point = center + K.vpol(center_to_foil_center_length * cos(Δα/2), orientation + (Δα/2))

    right_arc_center_displacement_vector = rounded_foil_center - starting_foil_point
    left_arc_center_displacement_vector = rounded_foil_center - ending_foil_point
    displacement_vectors_magnitude = K.norm(right_arc_center_displacement_vector) * displacement_ratio

    right_arc_center = displace_point_by_vector(starting_foil_point, right_arc_center_displacement_vector, displacement_vectors_magnitude)
    left_arc_center = displace_point_by_vector(ending_foil_point, left_arc_center_displacement_vector, displacement_vectors_magnitude)

    pointed_foil_radius = K.distance(right_arc_center, starting_foil_point)
    B = K.distance(right_arc_center, left_arc_center) / 2
    Z = midpoint(right_arc_center, left_arc_center)

    rounded_foil_center_to_arc_intersection_vector = Z - rounded_foil_center
    rounded_foil_center_to_arc_intersection_vector_magnitude = sqrt(pointed_foil_radius^2 - B^2)
    arc_intersection = displace_point_by_vector(Z, rounded_foil_center_to_arc_intersection_vector, rounded_foil_center_to_arc_intersection_vector_magnitude)
    
    scaling_factor = radius / K.distance(center, arc_intersection)
    current_rotation_angle = 0
    
    circle(center, radius)

    while n_foils > 0
        K.scale(K.rotate(arc(right_arc_center, starting_foil_point, arc_intersection), current_rotation_angle, center), scaling_factor)
        K.scale(K.rotate(arc(left_arc_center, arc_intersection, ending_foil_point), current_rotation_angle, center), scaling_factor)

        n_foils -= 1
        current_rotation_angle += Δα
    end
end
# == Rosette == #

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

    lancet_arch_midpoint = midpoint(left_point, right_point)
    excess_displacement_vector = right_point - left_point
    arcs_radius = K.norm(excess_displacement_vector) * excess

    left_arc_center = displace_point_by_vector(right_point, excess_displacement_vector, -arcs_radius)
    right_arc_center = displace_point_by_vector(left_point, excess_displacement_vector, arcs_radius)
    
    # This way we can avoid computationally complex operations such as circle_intersection
    # But we should be able to further abstract this concept by resorting to CGAL to be as close as possible to constructive methods
    arc_intersection_x = lancet_arch_midpoint.x
    arc_intersection_y = sqrt(arcs_radius^2 - K.distance(left_arc_center, lancet_arch_midpoint)^2) + lancet_arch_midpoint.y
    arc_intersection = K.xy(arc_intersection_x, arc_intersection_y)

    right_arc = arc(left_arc_center, right_point, arc_intersection)
    left_arc = arc(right_arc_center, arc_intersection, left_point)

    return [right_arc, left_arc, left_arc_center, arcs_radius]
end

# Or just the trefoil for simplification purposes
#function pointed_n_foil_arch_top(left_point, right_point)
#end
#
#function ogee_arch_top(left_point, right_point)
#end
# == Arch Tops == #

# == Main Arch == #
# Ask Prof. if I should specify the height of the sub_arches itself or if the ratio is a better practice
# Maybe resort to the multiple dispatch capabilites of Julia (?)
function arch(bottom_left_corner, upper_right_corner, excess, outer_offset, inner_offset, recursion_level, vertical_distance_to_sub_arches)
    # Arch body coordinates
    upper_left_corner = K.xy(bottom_left_corner.x, upper_right_corner.y)
    bottom_right_corner = K.xy(upper_right_corner.x, bottom_left_corner.y)

    # Arch top portion
    arcs = lancet_arch_top(upper_left_corner, upper_right_corner, excess)

    # Main arch
    main_arch = [arcs[1], arcs[2], upper_left_corner, bottom_left_corner, bottom_right_corner, upper_right_corner]
    polygon = K.line(main_arch[3:end])

    # 3D profile sweeps
    #profile = polygon()
    #sweep(arcs[1], profile)
    #sweep(arcs[2], profile)
    #sweep(polygon, circle(xy(0.2, 0), 0.2))
    
    if recursion_level > 0
        offset_bottom_left_corner = bottom_left_corner + K.vxy(outer_offset, outer_offset)
        offset_upper_right_corner = upper_right_corner - K.vxy(outer_offset, vertical_distance_to_sub_arches)
        offset_width_midpoint = midpoint(offset_bottom_left_corner, offset_upper_right_corner)
        inner_offset_displacement = inner_offset/2

        # TODO: Find a clever way to calculate the new excess
        # For now, the excess will be the same
        sub_arches_excess = excess

        left_sub_arch_upper_right_corner = K.xy(offset_width_midpoint.x - inner_offset_displacement, offset_upper_right_corner.y)
        right_sub_arch_bottom_left_corner = K.xy(offset_width_midpoint.x + inner_offset_displacement, offset_bottom_left_corner.y)

        left_sub_arch = arch(offset_bottom_left_corner, left_sub_arch_upper_right_corner, sub_arches_excess, outer_offset, inner_offset, recursion_level - 1, vertical_distance_to_sub_arches)
        right_sub_arch = arch(right_sub_arch_bottom_left_corner, offset_upper_right_corner, sub_arches_excess, outer_offset, inner_offset, recursion_level - 1, vertical_distance_to_sub_arches)

        # The code below begs for changes that conform to those mentioned in the corresponding labeled utility functions
        ellipse_center = midpoint(left_sub_arch[1], arcs[3])
        ellipse_rad = (arcs[4] + left_sub_arch[2]) / 2
        rosetteMid = intersect_line_ellipse(midpoint(K.xy(offset_bottom_left_corner.x, offset_upper_right_corner.y), offset_upper_right_corner), vy(1), ellipse_center, ellipse_rad)[2]
        rosetteRad = distance(rosetteMid, left_sub_arch[1]) - left_sub_arch[2] - inner_offset
        rosette_rounded_foils(rosetteMid, rosetteRad, 6, 0)
    end


    return [arcs[3], arcs[4]]
end
# == Main Arch == #

arch(K.xy(-10, -16), K.xy(10, 16), 1, 0.75, 0.75, 2, 3)

# Rosette decorations:
# Instead of solely developing standing and laying orientations for the foils, set these as two default is_arc_construction_possible
# and develop a new, general, construction that aims at allowing for different rotational configurations based upon their default positions
# Another thing that GML doesn't resort to as a simplification measure is polar coordinates! Something we can make use of to simplify these constructions
# as opposed to the GML approach of Cartesian coordinates

# Profiles (i.e., 3D)
# 2D planes swept along the region border curves (plane being orthogonal to the line or the tangent of the curve, in the case that the latter is a curve)
# at corner points, where the tangent is discontinuous, the sweep is continued onto the bisector plane)

# Instead of resorting to simple arrays to store information about arcs, circles, polygons (as GML does) we resort to Khepri's objects and types
# in such a way that allows us to perform type-checking operations as well as to better aid in the understanding of the program on the user's end

# Perguntar PROF se devo, ou não, recorrer à normal de um determinado plano de forma a especificar o arco em questão (de notar que, para a modelacao
# das janelas em questao, recorro sempre ao plano XY para depois fazer a transformacao desejada
