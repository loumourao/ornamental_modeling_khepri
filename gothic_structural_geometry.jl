module GothicStructuralGeometry
    module Arches
        using KhepriAutoCAD

        module Bodies
            using KhepriAutoCAD
            include("./geometry_utility_functions.jl")
            using .GeometryUtilityFunctions

            export get_left_sub_arch_body, get_right_sub_arch_body

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
        end
        using .Bodies
        export get_left_sub_arch_body, get_right_sub_arch_body

        module Tops
            using KhepriAutoCAD
            include("./geometry_utility_functions.jl")
            using .GeometryUtilityFunctions

            export get_offset_excess, lancet_arch_top

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
                    GeometryUtilityFunctions.Circles.Arcs.arc(arch_midpoint, right_point, left_point)
                else
                    excess_displacement_vector = right_point - left_point
                    arcs_radius = norm(excess_displacement_vector) * excess
            
                    right_arc_center = displace_point_by_vector(right_point, excess_displacement_vector, -arcs_radius)
                    left_arc_center = displace_point_by_vector(left_point, excess_displacement_vector, arcs_radius)
                    
                    arc_intersection_x = arch_midpoint.x
                    arc_intersection_y = sqrt(arcs_radius^2 - distance(left_arc_center, arch_midpoint)^2) + arch_midpoint.y
                    arc_intersection = xy(arc_intersection_x, arc_intersection_y)
            
                    right_arc = GeometryUtilityFunctions.Circles.Arcs.arc(right_arc_center, right_point, arc_intersection)
                    left_arc = GeometryUtilityFunctions.Circles.Arcs.arc(left_arc_center, arc_intersection, left_point)
                end
            
                return (right_arc = right_arc, left_arc = left_arc)
            end
        end
        using .Tops
        export get_offset_excess, lancet_arch_top
    end
    using .Arches
    export get_left_sub_arch_body, get_right_sub_arch_body
    export get_offset_excess, lancet_arch_top

    module Rosettes
        using KhepriAutoCAD
        include("./geometry_utility_functions.jl")
        using .GeometryUtilityFunctions

        export compute_rosette

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
    end
    using .Rosettes
    export compute_rosette

    module Fillets
        using KhepriAutoCAD
        include("./geometry_utility_functions.jl")
        using .GeometryUtilityFunctions

        export circular_rosette_fillets

        function circular_rosette_fillets(right_arc_center, left_arc_center, arcs_radius,
                                            rosette_center, rosette_radius, 
                                                right_sub_arch_right_arc_center, right_sub_arch_left_arc_center, 
                                                    left_sub_arch_left_arc_center, left_sub_arch_right_arc_center, sub_arcs_radius)
            # Right fillet points calculations
            right_fillet_bottom_point = intersect_circles(right_arc_center, arcs_radius, right_sub_arch_right_arc_center, sub_arcs_radius).greater_y_intersection_point
            right_fillet_top_point = intersect_circles(right_arc_center, arcs_radius, rosette_center, rosette_radius).lower_y_intersection_point
            right_fillet_left_point = intersect_circles(right_sub_arch_right_arc_center, sub_arcs_radius, rosette_center, rosette_radius).greater_y_intersection_point

            # Right fillet modeling
            right_fillet_right_arc = GeometryUtilityFunctions.Circles.Arcs.arc(right_arc_center, right_fillet_bottom_point, right_fillet_top_point)
            right_fillet_left_arc = GeometryUtilityFunctions.Circles.Arcs.arc(rosette_center, right_fillet_left_point, right_fillet_top_point)
            right_fillet_bottom_arc = GeometryUtilityFunctions.Circles.Arcs.arc(right_sub_arch_right_arc_center, right_fillet_bottom_point, right_fillet_left_point)

            # Top fillet points calculations
            top_fillet_right_point = intersect_circles(right_arc_center, arcs_radius, rosette_center, rosette_radius).greater_y_intersection_point
            top_fillet_top_point = intersect_circles(left_arc_center, arcs_radius, right_arc_center, arcs_radius).greater_y_intersection_point
            top_fillet_left_point = intersect_circles(left_arc_center, arcs_radius, rosette_center, rosette_radius).greater_y_intersection_point

            # Top fillet modeling
            GeometryUtilityFunctions.Circles.Arcs.arc(right_arc_center, top_fillet_right_point, top_fillet_top_point)
            GeometryUtilityFunctions.Circles.Arcs.arc(left_arc_center, top_fillet_top_point, top_fillet_left_point)
            GeometryUtilityFunctions.Circles.Arcs.arc(rosette_center, top_fillet_right_point, top_fillet_left_point)

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
            #GeometryUtilityFunctions.Circles.Arcs.arc(rosette_center, left_fillet_top_point, left_fillet_right_point)
            #GeometryUtilityFunctions.Circles.Arcs.arc(left_arc_center, left_fillet_top_point, left_fillet_bottom_point)
            #GeometryUtilityFunctions.Circles.Arcs.arc(left_sub_arch_left_arc_center, left_fillet_right_point, left_fillet_bottom_point)

            # Bottom fillet points calculations
            bottom_fillet_right_point = intersect_circles(right_sub_arch_left_arc_center, sub_arcs_radius, rosette_center, rosette_radius).lower_y_intersection_point
            bottom_fillet_left_point = intersect_circles(left_sub_arch_right_arc_center, sub_arcs_radius, rosette_center, rosette_radius).lower_y_intersection_point
            bottom_fillet_bottom_point = intersect_circles(left_sub_arch_right_arc_center, sub_arcs_radius, right_sub_arch_left_arc_center, sub_arcs_radius).greater_y_intersection_point

            # Bottom fillet modeling
            GeometryUtilityFunctions.Circles.Arcs.arc(right_sub_arch_left_arc_center, bottom_fillet_right_point, bottom_fillet_bottom_point)
            GeometryUtilityFunctions.Circles.Arcs.arc(rosette_center, bottom_fillet_left_point, bottom_fillet_right_point)
            GeometryUtilityFunctions.Circles.Arcs.arc(left_sub_arch_right_arc_center, bottom_fillet_bottom_point, bottom_fillet_left_point)
        end
    end
    using .Fillets
    export circular_rosette_fillets
end
