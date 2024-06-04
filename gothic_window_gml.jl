using Khepri
backend(autocad)

windowdict = Dict(
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

function midpoint_2pt(a, b)
    x = (a.x + b.x) / 2
    y = (a.y + b.y) / 2
    z = (a.z + b.z) / 2

    return xyz(x, y, z)
end

function setlength_vec(vec, length)
    return vec / norm(vec) * length
end

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

function intersect_line_ellipse(p0, p1, m0, m1, r)
    d = p1 - p0
    
    a = d.x^2 + d.y^2 + d.z^2
    b = 2 * (d.x * (p0.x - m0.x) + d.y * (p0.y - m0.y) + d.z * (p0.z - m0.z) + d.x * (p0.x - m1.x) + d.y * (p0.y - m1.y) + d.z * (p0.z - m1.z))
    c = (p0.x - m0.x)^2 + (p0.y - m0.y)^2 + (p0.z - m0.z)^2 + (p0.x - m1.x)^2 + (p0.y - m1.y)^2 + (p0.z - m1.z)^2 - r^2

    t_values = solve_quadratic(a, b, c)

    if length(t_values) == 0
        return nothing
    
    t_values = sort(t_values)
    t0 = t_values[1]
    t1 = t_values[2]

    q0 = p0 + d * t0
    q1 = p1 + d * t1

    return [q0, q1, t0, t1]
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

function gw_polygon_2arcs_height(arcL, arcR, hb, nrml, kseg)
    bR = arcR[1]
    bL = arcL[3]
    
    bR = xyz(bR.x, bR.y, hb)
    bL = xyz(bL.x, bL.y, hb)

    polygon = [bR]
    push!(polygon, circle_seg(arcR, nrml, kseg, 1)[1:end-1]...)
    push!(polygon, circle_seg(arcL, nrml, kseg, 1)...)
    push!(polygon, bL)
    
    return polygon
end

function gw_compute_arcs_rosette(pL, pR, nrml, excess, bdOuter, arcDown, bdInner)
    arch = gw_pointed_arch(pL, pR, excess, bdOuter, nrml)
    rad = pop!(arch)
    arcR = pop!(arch)
    arcL = pop!(arch)

    pLL = arcL[3] + vxyz(0, 0, -1) * arcDown
    pRR = arcR[1] + vyz(0, 0, -1) * arcDown
    pM = midpoint_2pt(pLL, pRR)
    dpM = setlength_vec(pRR - pLL, bdInner * 0.5)

    left_arch = gw_pointed_arch(pLL, pM - dpM, excess, 0.0, nrml)
    radL = pop!(left_arch)
    arcLR = pop!(left_arch)
    arcLL = pop!(left_arch)

    right_arch = gw_pointed_arch(pM + dpM, pRR, excess, 0.0, nrml)
    radR = pop!(right_arch)
    arcRR = pop!(right_arch)
    arcRL = pop!(right_arch)

    rosetteMid = intersect_line_ellipse(pM + vxyz(0, 0, 1), pM - vxyz(0, 0, 1), arcLR[2], arcR[2], rad + radL + bdInner)[1]
    rosetteRad = distance(rosetteMid, arcLR[2]) - radL - bdInner

    return [arcLL, arcLR, arcRL, rosetteMid, rosetteRad]
end

function gw_gothic_window(edgeWall, edgeBack, pBaseL, pBaseR, windowdict)
    nrml = facenormal(edgeWall)

    # DECORATE MAIN ARCH
    arcs = gw_pointed_arch(pBaseL, pBaseR, get(windowdict, "excess"), 0.0, nrml)
    main_arch = gw_polygon_2arcs_height(arcs[1], arcs[2], get(windowdict, "heightBott"), nrml, get(windowdict, "kseg"))
    edgeArch = style_main_arch(main_arch, get(windowdict, "bdOuter"), get(windowdict, "wallSetback"), edgeWall, edgeBack)

    # COMPUTE SUB-ARCS AND ROSETTE
    setBack = nrml * -get(windowdict, "wallSetback")
    pL = pBaseL + setBack
    pR = pBaseR + setBack
    sub_arcs_and_rosette = gw_compute_arcs_rosette(pL, pR, nrml, get(windowdict, "excess"), get(windowdict, "bdOuter"), get(windowdict, "arcDown"), get(windowdict, "bdInner"))
    rosetteRad = pop!(sub_arcs_and_rosette)
    rosetteMid = pop!(sub_arcs_and_rosette)
    arcRR = pop!(sub_arcs_and_rosette)
    arcRL = pop!(sub_arcs_and_rosette)
    arcLR = pop!(sub_arcs_and_rosette)
    arcLL = pop!(sub_arcs_and_rosette)

    # COMPUTE AND DECORATE THE FOUR FILLETS

    # DECORATE THE ROSETTE

    # DECORATE THE TWO SUB-ARCHES
end