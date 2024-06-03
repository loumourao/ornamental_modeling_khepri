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

# Utility functions begin
# IMPORTANT: For now, the normal is always pointing towards the viewer
function line_2pt(p0, p1, t)
    return p0 + (p1 - p0) * t
end

function intersect_circles(m0, r0, m1, r1, nrml)
    d = distance(m0, m1)

    # IMPORTANT: d == 0 prevents the set of interesection points when both circles overlap
    if d > r0 + r1 || d < abs(r0 - r1) || d == 0
        return nothing
    
    a = (r0^2 - r1^2 + d^2) / (2 * d)
    p = m0 + a * (m1 - m0) / d
    c = sqrt(r0^2 - a^2)
    c_vector = cross(m1 - m0, nrml)
    c_vector = c_vector / norm(c_vector) * c
    
    # TODO: Refactor the line below to execute the dot product correctly
    if cross(nrml, c_vector) * (m1 - m0) > 0
        first_intersection_point = p + c_vector
        second_intersection_point = p - c_vector
    else
        first_intersection_point = p - c_vector
        second_intersection_point = p + c_vector
    end

    return [first_intersection_point, second_intersection_point]
end

function move_2pt(p, q, t)
    v = q - p
    return p + v / norm(v) * t
end

function circle_seg_aux(a, m, b, n_points)
    ma = a - m
    mb = b - m
    # TODO: Refactor the line below to execute the dot product correctly
    angle = acos(ma * mb / (norm(ma) * norm(mb))) / (n_points + 1)
    rad = norm(ma)
    circle_seg = [a]

    for i = 1:n_points
        push!(circle_seg, m + vpol(rad, angle * i))
    end

    push!(circle_seg, b)

    return circle_seg
end

function circle_seg(point_array, nrml, n, mode)
    a = point_array[1]
    m = point_array[2]
    b = point_array[3]

    ma = a - m
    mb = b - m

    # TODO: Refactor the line below to execute the dot product correctly
    if cross(ma, mb) * nrml < 0
        prev_b = b
        b = a
        a = prev_b
    end

    if mode == 0
        circle_seg = circle_seg_aux(a, m, b, n)
    elseif mode == 1
        circle_seg = circle_seg_aux(a, m, b, n + 1)
    elseif mode == 2
        # TODO: Ask for help here
    end

    return circle_seg
end
# Utility functions end

# Pointed Arch
function gw_pointed_arch(pL, pR, excess, offset, nrml)
    mR = line_2pt(pR, pL, excess)
    mL = line_2pt(pL, pR, excess)
    rad = distance(pL, pR) * excess - offset
    qT = intersect_circles(mL, rad, mR, rad, nrml)[1]
    
    left_circle_segment = [qT, mL, move_2pt(pL, pR, offset)]
    right_circle_segment = [move_2pt(pR, pL, offset), mR, qT]
    
    return [left_circle_segment, right_circle_segment, rad]
end

function gw_polygon_2arcs_height(kseg, nrml, hb, arcR, arcL)
    bR = arcR[1]
    bL = arcL[3]
    
    bR = xyz(bR.x, bR.y, hb)
    bL = xyz(bL.x, bL.y, hb)

    polygon = [bR]
    push!(polygon, circle_seg(arcR, nrml, kseg, 1)[1:end-1])
    push!(polygon, circle_seg(arcL, nrml, kseg, 1))
    push!(polygon, bL)
    
    return collect(Iterators.flatten(polygon))
end
