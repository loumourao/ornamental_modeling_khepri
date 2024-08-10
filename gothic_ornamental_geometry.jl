module GothicOrnamentalGeometry
    module SolidOrnamentations
        module Arches
            using KhepriAutoCAD
            include("./geometry_utility_functions.jl")
            using .GeometryUtilityFunctions

            export three_dimensionalize_arch_top, three_dimensionalize_arch_body, three_dimensionalize_arch_middle

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
        end
        using .Arches
        export three_dimensionalize_arch_top, three_dimensionalize_arch_body, three_dimensionalize_arch_middle

        module Rosettes
            using KhepriAutoCAD
            include("./geometry_utility_functions.jl")
            using .GeometryUtilityFunctions

            export three_dimensionalize_rosette, three_dimensionalize_rosette_rounded_foils, three_dimensionalize_rosette_pointed_foils

            function three_dimensionalize_rosette(rosette_center, rosette_radius, 
                                                    outer_offset, inner_offset, 
                                                        profile = surface_circle(u0(), (outer_offset + inner_offset) / 2))
                rosette = circle(rosette_center, rosette_radius - outer_offset + ((outer_offset + inner_offset) / 2))
                sweep(rosette, profile)
            end

            function three_dimensionalize_rosette_rounded_foils(foil, connection, 
                                                                    offset_value, profile = surface_circle(u0(), abs(offset_value) / 2))
                foil = arc_bidirectionally_extended_uniform_offset(foil, offset_value / 2)
                connection = nothing
                
                foil = sweep(foil, profile)

                return (foil = foil, connection = connection)
            end

            function three_dimensionalize_rosette_pointed_foils(foil_right_arc, foil_left_arc, connection, 
                                                                    offset_value, profile = surface_circle(u0(), abs(offset_value) / 2))
                foil_right_arc = arc_bidirectionally_extended_uniform_offset(foil_right_arc, offset_value / 2)
                foil_left_arc = arc_bidirectionally_extended_uniform_offset(foil_left_arc, offset_value / 2)
                connection = nothing

                foil_right_arc = sweep(foil_right_arc, profile)
                foil_left_arc = sweep(foil_left_arc, profile)

                return (foil_right_arc = foil_right_arc, foil_left_arc = foil_left_arc, connection = connection)
            end
        end
        using .Rosettes
        export three_dimensionalize_rosette, three_dimensionalize_rosette_rounded_foils, three_dimensionalize_rosette_pointed_foils
    end
    using .SolidOrnamentations
    export three_dimensionalize_arch_top, three_dimensionalize_arch_body, three_dimensionalize_arch_middle
    export three_dimensionalize_rosette, three_dimensionalize_rosette_rounded_foils, three_dimensionalize_rosette_pointed_foils

    module Rosettes
        using KhepriAutoCAD
        include("./geometry_utility_functions.jl")
        using .GeometryUtilityFunctions
        using ..SolidOrnamentations

        export get_rounded_foil, get_pointed_foil, rosette_rounded_foils, rosette_pointed_foils

        function get_rounded_foil(rosette_center, rosette_radius, Δα, orientation)
            foil_radius = (rosette_radius * sin(Δα/2)) / (1 + sin(Δα/2))
            rosette_center_to_foil_center_length = rosette_radius - foil_radius
            rosette_center_to_foil_end_points = rosette_center_to_foil_center_length * cos(Δα/2)
        
            foil_center = rosette_center + vpol(rosette_center_to_foil_center_length, orientation)
            foil_start_point = rosette_center + vpol(rosette_center_to_foil_end_points, orientation - (Δα/2))
            foil_end_point = rosette_center + vpol(rosette_center_to_foil_end_points, orientation + (Δα/2))
        
            return GeometryUtilityFunctions.Circles.Arcs.arc(foil_center, foil_start_point, foil_end_point)
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
        
            pointed_foil_right_arc = GeometryUtilityFunctions.Circles.Arcs.arc(pointed_foil_right_arc_center, rounded_foil_start_point, arc_intersection) 
            pointed_foil_left_arc = GeometryUtilityFunctions.Circles.Arcs.arc(pointed_foil_left_arc_center, arc_intersection, rounded_foil_end_point)
        
            return (right_arc = pointed_foil_right_arc, left_arc = pointed_foil_left_arc, rounded_foil_center = rounded_foil_center)
        end

        module Fillets
            using KhepriAutoCAD
            include("./geometry_utility_functions.jl")
            using .GeometryUtilityFunctions
            export get_rosette_rounded_foils_fillet, get_rosette_pointed_foils_fillet

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
            
                fillet_right_arc = GeometryUtilityFunctions.Circles.Arcs.arc(fillet_right_arc_center, fillet_right_point, fillet_bottom_point)
                fillet_upper_arc = GeometryUtilityFunctions.Circles.Arcs.arc(fillet_upper_arc_center, fillet_right_point, fillet_left_point)
                fillet_left_arc = GeometryUtilityFunctions.Circles.Arcs.arc(fillet_left_arc_center, fillet_bottom_point, fillet_left_point)
            
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

                fillet_right_point = intersect_circles(rosette_center, rosette_radius, fillet_right_arc_center, foil_radius).greater_y_intersection_point
                fillet_left_point = intersect_circles(rosette_center, rosette_radius, fillet_left_arc_center, foil_radius).lower_y_intersection_point
                fillet_bottom_point = displace_point_by_vector(right_foil_center, displacement_vector, norm(displacement_vector) / 2)

                fillet_right_arc = GeometryUtilityFunctions.Circles.Arcs.arc(fillet_right_arc_center, fillet_right_point, fillet_bottom_point)
                fillet_upper_arc = GeometryUtilityFunctions.Circles.Arcs.arc(fillet_upper_arc_center, fillet_right_point, fillet_left_point)
                fillet_left_arc = GeometryUtilityFunctions.Circles.Arcs.arc(fillet_left_arc_center, fillet_bottom_point, fillet_left_point)

                return (right_arc = fillet_right_arc, upper_arc = fillet_upper_arc, left_arc = fillet_left_arc)
            end
        end
        using .Fillets
        export get_rosette_rounded_foils_fillet, get_rosette_pointed_foils_fillet

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
    end
    using .Rosettes
    export get_rounded_foil, get_pointed_foil, rosette_rounded_foils, rosette_pointed_foils
    export get_rosette_rounded_foils_fillet, get_rosette_pointed_foils_fillet
end
