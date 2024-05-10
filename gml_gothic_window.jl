using Khepri
backend(autocad)

delete_all_shapes()

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

    # main arch calculations
    union(body, arch)
end

function sub_arch(center, width, height, arch_excess)

end

main_arch(xy(0,0), 10, 10, 1)