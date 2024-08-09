module GothicOrnamentalGeometry
# == ROSETTE ORNAMENTATION == #
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
    rounded_foil_center = arc_center(rounded_foil)
    rounded_foil_start_point = arc_start_point(rounded_foil)
    rounded_foil_end_point = arc_end_point(rounded_foil)

    pointed_foil_right_arc_center_displacement_vector = rounded_foil_center - rounded_foil_start_point
    pointed_foil_left_arc_center_displacement_vector = rounded_foil_center - rounded_foil_end_point
    displacement_vectors_magnitude = norm(foil_right_arc_center_displacement_vector) * displacement_ratio

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

    fillet_right_point = intersect_circles(rosette_center, rosette_radius, fillet_right_arc_center, foil_radius)[1]
    fillet_left_point = intersect_circles(rosette_center, rosette_radius, fillet_left_arc_center, foil_radius)[2]
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
    fillet_left_arc_center = rosette_center + vpol(rosette_center_to_foil_arc_center_length, Δα + Δβ)

    displacement_vector = left_foil_center - right_foil_center

    fillet_right_point = intersect_circles(rosette_center, rosette_radius, fillet_right_arc_center, foil_radius)[1]
    fillet_left_point = intersect_circles(rosette_center, rosette_radius, fillet_left_arc_center, foil_radius)[2]
    fillet_bottom_point = displace_point_by_vector(right_foil_center, displacement_vector, norm(displacement_vector) / 2)

    fillet_right_arc = arc(fillet_right_arc_center, fillet_right_point, fillet_bottom_point)
    fillet_upper_arc = arc(fillet_upper_arc_center, fillet_right_point, fillet_left_point)
    fillet_left_arc = arc(fillet_left_arc_center, fillet_bottom_point, fillet_left_point)

    return (right_arc = fillet_right_arc, upper_arc = fillet_upper_arc, left_arc = fillet_left_arc)
end
# == FILLETS == #

function rosette_rounded_foils(rosette_center, rosette_radius, n_foils, orientation, inner_offset)
    # Check if n_foils >= 1, otherwise raise an error
    Δα = 2π / n_foils

    foil = get_rounded_foil(rosette_center, rosette_radius, Δα, orientation)
    fillet = get_rosette_rounded_foils_fillet(rosette_center, rosette_radius, 
                                                arc_center(foil), arc_radius(foil), 
                                                    Δα, inner_offset)

    current_rotation_angle = 0

    while n_foils > 0
        # Fillets
        fillet_arc_position = current_rotation_angle + orientation
        rotate(fillet.right_arc, fillet_arc_position, rosette_center)
        rotate(fillet.upper_arc, fillet_arc_position, rosette_center)
        rotate(fillet.left_arc, fillet_arc_position, rosette_center)

        # Rounded foils
        foil = arc_bidirectionally_extended_uniform_offset(foil, inner_offset)
        rotate(foil, current_rotation_angle, rosette_center)

        # Rounded foils connections
        connection_start_point = arc_end_point(foil)
        rosette_center_to_foil_start_point = arc_start_point(foil) - rosette_center
        rosette_center_to_foil_start_point_polar_angle = atan(rosette_center_to_foil_start_point.y, rosette_center_to_foil_start_point.x)
        connection_end_point = rosette_center + vpol(norm(rosette_center_to_foil_start_point), rosette_center_to_foil_start_point_polar_angle + Δα)
        connection = line(connection_start_point, connection_end_point)
        rotate(connection, current_rotation_angle, rosette_center)

        # 3D rosette rounded foils
        three_dimensionalized_foil_and_connection = three_dimensionalize_rosette_rounded_foils(foil, connection, -inner_offset)
        rotate(three_dimensionalized_foil_and_connection.foil, current_rotation_angle, rosette_center)
        #rotate(three_dimensionalized_foil_and_connection.connection, current_rotation_angle, rosette_center)

        current_rotation_angle += Δα
        n_foils -= 1
    end
end

function rosette_pointed_foils(rosette_center, rosette_radius, n_foils, displacement_ratio, orientation, inner_offset)
    # Check if n_foils >= 1 and if displacement_ratio > 1, otherwise raise an error
    Δα = 2π / n_foils

    scaling_factor = rosette_radius / distance(rosette_center, arc_intersection)
    
    foil = get_pointed_foil(rosette_center, rosette_radius, Δα, displacement_ratio, orientation)
    rounded_foil_center = foil.rounded_foil_center
    outer_foil_right_arc = foil.right_arc 
    outer_foil_left_arc = foil.left_arc
    foil_radius = arc_radius(outer_foil_right_arc)

    fillet = get_rosette_pointed_foils_fillet_points(rosette_center, rosette_radius, 
                                                        rounded_foil_center, arc_center(outer_foil_right_arc), foil_radius, 
                                                            Δα, scaling_factor, inner_offset)

    current_rotation_angle = 0

    while n_foils > 0
        # Fillets
        fillet_arc_position = current_rotation_angle + orientation
        scale(rotate(fillet.right_arc, fillet_arc_position, rosette_center), scaling_factor, rosette_center)
        scale(rotate(fillet.upper_arc, fillet_arc_position, rosette_center), scaling_factor, rosette_center)
        scale(rotate(fillet.left_arc, fillet_bottom_point, fillet_left_point), scaling_factor, rosette_center)

        # Pointed foils
        foil_right_arc = arc_start_angle_extended_offset(outer_foil_right_arc, inner_offset)
        foil_left_arc = arc_amplitude_extended_offset(outer_foil_left_arc, inner_offset)
        scale(rotate(foil_right_arc, current_rotation_angle, rosette_center), scaling_factor, rosette_center)
        scale(rotate(foil_left_arc, current_rotation_angle, rosette_center), scaling_factor, rosette_center)

        # Pointed foils connections
        connection_start_point = arc_end_point(foil_left_arc)
        rosette_center_to_foil_start_point = arc_start_point(foil_right_arc) - rosette_center
        rosette_center_to_foil_start_point_polar_angle = atan(rosette_center_to_foil_start_point.y, rosette_center_to_foil_start_point.x)
        connection_end_point = rosette_center + vpol(norm(rosette_center_to_foil_start_point), rosette_center_to_foil_start_point_polar_angle + Δα)
        connection = line(connection_start_point, connection_end_point)
        scale(rotate(connection, current_rotation_angle, rosette_center), scaling_factor, rosette_center)

        # 3D rosette pointed foils
        three_dimensionalized_foil_and_connection = three_dimensionalize_rosette_pointed_foils(outer_foil_right_arc, outer_foil_left_arc, 
                                                                                                    connection, inner_offset)
        scale(rotate(three_dimensionalized_foil_and_connection.foil_right_arc, current_rotation_angle, rosette_center), scaling_factor, rosette_center)
        scale(rotate(three_dimensionalized_foil_and_connection.foil_left_arc, current_rotation_angle, rosette_center), scaling_factor, rosette_center)
        #scale(rotate(three_dimensionalized_foil_and_connection.connection, current_rotation_angle, rosette_center), scaling_factor, rosette_center)

        n_foils -= 1
        current_rotation_angle += Δα
    end
end
# == ROSETTE ORNAMENTATION == #
end
