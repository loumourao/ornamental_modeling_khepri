using Khepri
backend(autocad)

window_dict = Dict(
                    "excess" => 1.25,
                    "arcDown" => 2.0,
                    "bdInner" => 0.4,
                    "bdOuter" => 0.5,
                    "wallSetback" => 0.1,
                    "heightBott" => 2.0,
                    "kseg" => 6,
                    "Style" => 0
                  )

# Define 5.32 below

# Define 5.33 below

# Define 5.36 below

# Define line_2pt below
function line_2pt(p0, p1, t)
    return (1 - t) * p0 + t * p1
end

# Define intersect_circles below
function intersect_circles(m0, r0, m1, r1, nrml)
    circle_center_distance = distance(m0, m1)

    if circle_center_distance > r0 + r1 || circle_center_distance < abs(r0 - r1)
        return nothing
    
    m0_to_intersection_axis_midpoint = (r0^2 - r1^2 + circle_center_distance^2) / (2 * circle_center_distance)
    intersection_axis_midpoint = m0 + m0_to_intersection_axis_midpoint * (m1 - m0) / d
    intersection_axis_midpoint_to_intersection_point = sqrt(r^2 - m0_to_intersection_axis_midpoint^2)
    intersection_axis_midpoint_to_intersection_point = intersection_axis_midpoint_to_intersection_point * [m1[1] - m0[1], m1[0] - m0[0]] / circle_center_distance

    first_intersection_point = intersection_axis_midpoint + intersection_axis_midpoint_to_intersection_point
    second_intersection_point = intersection_axis_midpoint - intersection_axis_midpoint_to_intersection_point

    return [first_intersection_point, second_intersection_point]
end

# Pointed Arch
function gw_pointed_arch(pL, pR, excess, offset, nrml)
    mR = line_2pt(pR, pL, excess)
    mL = line_2pt(pL, pR, excess)
    rad = (distance(pL, pR) * excess) - offset
    qT = intersect_circles(mL, rad, mR, rad, nrml) # pop the value that doesn't interest us
end

# Define 5.38 below

# Define 5.39 below

# Define circle_seg below

# Define move_2pt below
