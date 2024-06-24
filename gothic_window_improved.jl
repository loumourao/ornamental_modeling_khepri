using KhepriAutoCAD

delete_all_shapes()

# Utility Functions #
function line_2pt(p0, p1, t)
    return p0 + (p1 - p0) * t
end

function intersect_circles(m0, r0, m1, r1, nrml)
    d = distance(m0, m1)

    # IMPORTANT: d == 0 prevents the set of interesection points when both circles overlap
    if d > r0 + r1 || d < abs(r0 - r1) || d == 0
        return nothing
    end

    a = (r0^2 - r1^2 + d^2) / (2 * d)
    p = m0 + a * (m1 - m0) / d
    c = sqrt(r0^2 - a^2)
    c_vector = cross(m1 - m0, nrml)
    c_vector = c_vector / norm(c_vector) * c
    
    if dot(cross(nrml, c_vector), (m1 - m0)) > 0
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

function circle_seg_aux(a, m, b, rad, n_points)
    ma = a - m
    mb = b - m
    m = loc_from_o_vx_vy(m, ma, mb)
    angle_increment = angle_between(ma, mb) / n_points
    current_angle = 0

    circle_seg = [a]

    while n_points > 1
        current_angle += angle_increment
        push!(circle_seg, m + vpol(rad, current_angle, m.cs))
        n_points -= 1
    end

    push!(circle_seg, b)

    return circle_seg
end

function circle_seg(point_array, nrml, n, mode)
    a = point_array[1]
    m = point_array[2]
    b = point_array[3]
    rad = distance(m, a)

    if mode == 0
        circle_seg = circle_seg_aux(a, m, b, rad, n)
    elseif mode == 1
        circle_seg = circle_seg_aux(a, m, b, rad, n + 1)
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

function intersect_line_ellipse(p0, v, m, r)
    a = norm(v)^2
    b = 2 * dot((p0 - m), v)
    c = norm(p0 - m)^2 - r^2

    t_values = solve_quadratic(a, b, c)
    
    if length(t_values) == 0
        return nothing
    end

    t_values = sort(t_values)
    t0 = t_values[1]
    t1 = t_values[2]

    q0 = p0 + v * t0
    q1 = p0 + v * t1

    return [q0, q1, t0, t1]
end

function circle(center, nrml, rad, n_points)
    cx = center + vx(rad) - center
    cy = center + vy(rad) - center
    m = loc_from_o_vx_vy(center, cx, cy)
    angle_increment = 2 * pi / n_points
    current_angle = 0

    circle = [center + cx]

    while n_points > 1
        current_angle += angle_increment
        push!(circle, m + vpol(rad, current_angle, m.cs))
        n_points -= 1
    end

    return circle
end
# Utility Functions #

# Style Functions #
function style_main_arch(poly, bdOuter, wallSetback)
    polygon(poly)
end

function style_fillet(poly)
    polygon(poly)
end

function style_rosette(poly, mid, rad)
    polygon(poly)
end

function style_sub_arch(poly, arcL, bh)
    polygon(poly)
end
# Style Functions #

# Structure Functions #
struct gothic_window
    excess
    arcDown
    bdInner
    bdOuter
    wallSetback
    heightBott
    kseg
    Style
end

function gw_pointed_arch(pL, pR, excess, offset, nrml)
    mR = line_2pt(pR, pL, excess)
    mL = line_2pt(pL, pR, excess)
    rad = distance(pL, pR) * excess - offset
    qT = intersect_circles(mL, rad, mR, rad, nrml)[1]

    arcL = [qT, mL, move_2pt(pL, pR, offset)]
    arcR = [move_2pt(pR, pL, offset), mR, qT]

    return [arcL, arcR, rad]
end

function gw_polygon_2arcs_height(arcL, arcR, hb, nrml, kseg)
    bR = xy(arcR[1].x, -hb)
    bL = xy(arcL[3].x, -hb)

    polygon = [bR]
    push!(polygon, circle_seg(arcR, nrml, kseg, 1)[1:end-1]...)
    push!(polygon, circle_seg(arcL, nrml, kseg, 1)...)
    push!(polygon, bL)
    push!(polygon, bR)

    return polygon
end

function gw_compute_arcs_rosette(pL, pR, nrml, windowdict)
    sub_arch = gw_pointed_arch(pL, pR, windowdict.excess, windowdict.bdOuter, nrml)
    arcL = sub_arch[1]
    arcR = sub_arch[2]
    rad = sub_arch[3]

    pLL = arcL[3] + vy(-1) * windowdict.arcDown
    pRR = arcR[1] + vy(-1) * windowdict.arcDown
    pM = midpoint_2pt(pLL, pRR)
    dpM = setlength_vec(pRR - pLL, windowdict.bdInner * 0.5)

    sub_archL = gw_pointed_arch(pLL, pM - dpM, windowdict.excess, 0.0, nrml)
    arcLL = sub_archL[1]
    arcLR = sub_archL[2]
    radL = sub_archL[3]

    sub_archR = gw_pointed_arch(pM + dpM, pRR, windowdict.excess, 0.0, nrml)
    arcRL = sub_archR[1]
    arcRR = sub_archR[2]
    radR = sub_archR[3]

    ellipse_center = midpoint_2pt(arcLR[2], arcR[2])
    ellipse_rad = (rad + radL) / 2
    rosetteMid = intersect_line_ellipse(pM, vy(1), ellipse_center, ellipse_rad)[2]
    rosetteRad = distance(rosetteMid, arcLR[2]) - radL - windowdict.bdInner

    return [arcLL, arcLR, arcRL, arcRR, rosetteMid, rosetteRad]
end

function gw_compute_fillets(arcLL, arcLR, arcRL, arcRR, rosetteMid, rosetteRad, pL, pR, nrml, windowdict)
    arcs = gw_pointed_arch(pL, pR, windowdict.excess, windowdict.bdOuter, nrml)
    arcL = arcs[1]
    arcR = arcs[2]
    
    circleL_mid = arcL[2]
    circleR_mid = arcR[2]
    circleLL_mid = arcLL[2]
    circleLR_mid = arcLR[2]
    circleRL_mid = arcRL[2]
    circleRR_mid = arcRR[2]

    circleL_rad = distance(arcL[1], circleL_mid)
    circleR_rad = distance(arcR[1], circleR_mid)
    circleLL_rad = distance(arcLL[1], circleLL_mid) + windowdict.bdOuter
    circleLR_rad = distance(arcLR[1], circleLR_mid) + windowdict.bdInner
    circleRL_rad = distance(arcRL[1], circleRL_mid) + windowdict.bdInner
    circleRR_rad = distance(arcRR[1], circleRR_mid) + windowdict.bdOuter

    rosetteRad += windowdict.bdInner
    
    ros_arcL = intersect_circles(rosetteMid, rosetteRad, circleL_mid, circleL_rad, nrml)
    ros_arcR = intersect_circles(rosetteMid, rosetteRad, circleR_mid, circleR_rad, nrml)
    arcL_arcR = intersect_circles(circleL_mid, circleL_rad, circleR_mid, circleR_rad, nrml)
    arcL_arcLL = intersect_circles(circleL_mid, circleL_rad, circleLL_mid, circleLL_rad, nrml)
    arcR_arcRR = intersect_circles(circleR_mid, circleR_rad, circleRR_mid, circleRR_rad, nrml)
    arcLR_arcRL = intersect_circles(circleLR_mid, circleLR_rad, circleRL_mid, circleRL_rad, nrml)
    arcLL_ros = intersect_circles(circleLL_mid, circleLL_rad, rosetteMid, rosetteRad, nrml)
    arcRR_ros = intersect_circles(circleRR_mid, circleRR_rad, rosetteMid, rosetteRad, nrml)
    arcLR_ros = intersect_circles(circleLR_mid, circleLR_rad, rosetteMid, rosetteRad, nrml)
    arcRL_ros = intersect_circles(circleRL_mid, circleRL_rad, rosetteMid, rosetteRad, nrml)

    right_fillet_points = [arcR_arcRR[2], ros_arcR[2], arcRR_ros[1]]
    upper_fillet_points = [ros_arcR[1], arcL_arcR[1], ros_arcL[2]]
    left_fillet_points = [arcLL_ros[2], ros_arcL[1], arcL_arcLL[1]]
    lower_fillet_points = [arcRL_ros[2], arcLR_ros[1], arcLR_arcRL[2]]

    right_fillet = [[right_fillet_points[1], circleR_mid, right_fillet_points[2]], 
                    [right_fillet_points[2], rosetteMid, right_fillet_points[3]], 
                    [right_fillet_points[3], circleRR_mid, right_fillet_points[1]]]

    upper_fillet = [[upper_fillet_points[1], circleR_mid, upper_fillet_points[2]],
                    [upper_fillet_points[2], circleL_mid, upper_fillet_points[3]],
                    [upper_fillet_points[3], rosetteMid, upper_fillet_points[1]]]

    left_fillet = [[left_fillet_points[1], rosetteMid, left_fillet_points[2]],
                   [left_fillet_points[2], circleL_mid, left_fillet_points[3]],
                   [left_fillet_points[3], circleLL_mid, left_fillet_points[1]]]

    lower_fillet = [[lower_fillet_points[1], rosetteMid, lower_fillet_points[2]],
                    [lower_fillet_points[2], circleLR_mid, lower_fillet_points[3]], 
                    [lower_fillet_points[3], circleRL_mid, lower_fillet_points[1]]]

    return [right_fillet, upper_fillet, left_fillet, lower_fillet]
end

function gw_polygon_fillets(right_fillet, upper_fillet, left_fillet, lower_fillet, nrml, windowdict)
    print(upper_fillet)
    new_right_fillet = vcat([circle_seg(fillet, nrml, windowdict.kseg, 1) for fillet in right_fillet]...)
    new_upper_fillet = vcat([circle_seg(fillet, nrml, windowdict.kseg, 1) for fillet in upper_fillet]...)
    new_left_fillet = vcat([circle_seg(fillet, nrml, windowdict.kseg, 1) for fillet in left_fillet]...)
    new_lower_fillet = vcat([circle_seg(fillet, nrml, windowdict.kseg, 1) for fillet in lower_fillet]...)

    return [new_right_fillet, new_upper_fillet, new_left_fillet, new_lower_fillet]
end

function gw_gothic_window(windowdict, pBaseL, pBaseR, nrml)
    # DECORATE MAIN ARCH
    arcs = gw_pointed_arch(pBaseL, pBaseR, windowdict.excess, 0.0, nrml)
    arcL = arcs[1]
    arcR = arcs[2]
    main_arch = gw_polygon_2arcs_height(arcL, arcR, windowdict.heightBott, nrml, windowdict.kseg)
    style_main_arch(main_arch, windowdict.bdOuter, windowdict.wallSetback)

    # COMPUTE SUB-ARCS AND ROSETTE
    window_setback = nrml * -windowdict.wallSetback
    pL = pBaseL + window_setback
    pR = pBaseR + window_setback
    arcs_and_rosette = gw_compute_arcs_rosette(pL, pR, nrml, windowdict)
    arcLL = arcs_and_rosette[1]
    arcLR = arcs_and_rosette[2]
    arcRL = arcs_and_rosette[3]
    arcRR = arcs_and_rosette[4]
    rosetteMid = arcs_and_rosette[5]
    rosetteRad = arcs_and_rosette[6]

    # COMPUTE AND DECORATE THE FOUR FILLETS
    fillets = gw_compute_fillets(arcLL, arcLR, arcRL, arcRR, rosetteMid, rosetteRad, pL, pR, nrml, windowdict)
    polygon_fillets = gw_polygon_fillets(fillets..., nrml, windowdict)

    for fillet in polygon_fillets
        style_fillet(fillet)
    end

    # DECORATE THE ROSETTE
    rosette_circle = circle(rosetteMid, nrml, rosetteRad, windowdict.kseg * 4)
    style_rosette(rosette_circle, rosetteMid, rosetteRad)

    # DECORATE THE TWO SUB-ARCHES
    heightInnerArc = windowdict.heightBott - windowdict.bdOuter

    archR = gw_polygon_2arcs_height(arcRL, arcRR, heightInnerArc, nrml, windowdict.kseg)
    style_sub_arch(archR, arcRL, nothing)

    archL = gw_polygon_2arcs_height(arcLL, arcLR, heightInnerArc, nrml, windowdict.kseg)
    style_sub_arch(archL, arcLL, nothing)
end
# Structure Functions #

# Modeling Test #
gw = gothic_window(1.25, 2.0, 0.4, 0.5, 0.1, 10, 6, nothing)
nrml = vz(1)
pL = xy(-5, 0)
pR = xy(5, 0)
gw_gothic_window(gw, pL, pR, nrml)
# Modeling Test #
