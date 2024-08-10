module GothicWindows
    using KhepriAutoCAD
    include("./geometry_utility_functions.jl")
    using .GeometryUtilityFunctions
    include("./gothic_structural_geometry.jl")
    using .GothicStructuralGeometry
    include("./gothic_ornamental_geometry.jl")
    using .GothicOrnamentalGeometry

    export gothic_window

    function gothic_window(bottom_left_corner, upper_right_corner, excess, 
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
            left_sub_arch = gothic_window(left_sub_arch_body.bottom_left_corner, left_sub_arch_body.upper_right_corner, sub_arch_excess, 
                                                        recursion_level - 1, vertical_distance_to_sub_arch, 
                                                            sub_arch_outer_offset, sub_arch_inner_offset)
            right_sub_arch = gothic_window(right_sub_arch_body.bottom_left_corner, right_sub_arch_body.upper_right_corner, sub_arch_excess, 
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
end
