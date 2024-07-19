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

    amplitude = ccw_starting_angle < ccw_ending_angle ? ccw_ending_angle - ccw_starting_angle : 2*pi - K.angle_between(ca, cb)
    radius = K.norm(ca)

    return K.arc(center, radius, ccw_starting_angle, amplitude)
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
    arcs_radius = K.distance(left_point, right_point) * excess

    left_arc_center = displace_point_by_vector(right_point, excess_displacement_vector, -arcs_radius)
    right_arc_center = displace_point_by_vector(left_point, excess_displacement_vector, arcs_radius)
    
    # This way we can avoid computationally complex operations such as circle_intersection
    # But we should be able to further abstract this concept by resorting to CGAL to be as close as possible to constructive methods
    arc_intersection_x = lancet_arch_midpoint.x
    arc_intersection_y = sqrt(arcs_radius^2 - K.distance(left_arc_center, lancet_arch_midpoint)^2) + lancet_arch_midpoint.y
    arc_intersection = K.xy(arc_intersection_x, arc_intersection_y)

    right_arc = arc(left_arc_center, right_point, arc_intersection)
    left_arc = arc(right_arc_center, arc_intersection, left_point)

    return [right_arc, left_arc]
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
    sweep(arcs[1], circle(u0(), 0.4))
    sweep(arcs[2], circle(u0(), 0.4))
    sweep(polygon, circle(u0(), 0.4))
    
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
    end
end
# == Main Arch == #

arch(K.xy(-4, -6), K.xy(4, 6), 1, 0.4, 0.4, 1, 3)

# For sub-arches, the circle midpoints of the original arch are kept the same
# the only thing that changes is the offset of the upper corners
# TODO: Try to automatically calculate the new excess provided this piece of information

# Rosette decorations:
# Develop an n-foil generative function
# Develop sub-types for n-foil (pointed and rounded n-foils)
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
