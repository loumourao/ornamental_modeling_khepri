using Khepri
backend(autocad)

delete_all_shapes()

function gothic_window(center, width, height, arch_excess, recursion_level, distance_to_sub_arch)
    #=
    center - the normal of a given window
    width - the width of the main arch
    height - the height from bottom to arch init (not the whole height - should discover a different way to include the whole arch later on - most likely related to excess)
    =#
    function main_arch(center, width, height, arch_excess)
        # body calculations
        body_bottom_left_corner = center - vxy(width/2, height/2)
        body = surface(rectangle(body_bottom_left_corner, width, height))

        # arch calculations
        body_upper_left_corner = body_bottom_left_corner + vy(height)
        body_upper_right_corner = body_upper_left_corner + vx(width)
        circle_radius = distance(body_upper_left_corner, body_upper_right_corner) * arch_excess
        left_circle_midpoint = body_upper_left_corner + vx(circle_radius)
        right_circle_midpoint = body_upper_right_corner - vx(circle_radius)
        left_circle = surface(circle(left_circle_midpoint, circle_radius))
        right_circle = surface(circle(right_circle_midpoint, circle_radius))
        arch = intersection(left_circle, right_circle)

        # arch
        union(body, arch)
    end

    function get_displaced_sub_archs_centers(previous_center, new_width, distance_to_parent_arch)
        displaced_center_x = new_width/2
        displaced_center_y = -distance_to_parent_arch/2
        left_sub_arch_center = previous_center + vxy(-displaced_center_x, displaced_center_y)
        right_sub_arch_center = previous_center + vxy(displaced_center_x, displaced_center_y)
        return (left_sub_arch_center, right_sub_arch_center)
    end

    main = main_arch(center, width, height, arch_excess)

    if recursion_level > 0
        sub_arch_width = width/2
        sub_arch_height = height - distance_to_sub_arch[1]
        sub_archs_centers = get_displaced_sub_archs_centers(center, sub_arch_width, distance_to_sub_arch[1])

        left_sub_arch = gothic_window(sub_archs_centers[1], sub_arch_width, sub_arch_height, arch_excess, recursion_level - 1, distance_to_sub_arch[2:end])
        right_sub_arch = gothic_window(sub_archs_centers[2], sub_arch_width, sub_arch_height, arch_excess, recursion_level - 1, distance_to_sub_arch[2:end])
    end
end

gothic_window(xy(0,0), 10, 10, 1, 2, [3,3])