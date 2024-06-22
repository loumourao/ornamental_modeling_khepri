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

function intersect_line_ellipse(p0, p1, m0, m1, r)
    d = p1 - p0
    
    a = d.x^2 + d.y^2 + d.z^2
    b = 2 * (d.x * (p0.x - m0.x) + d.y * (p0.y - m0.y) + d.z * (p0.z - m0.z) + d.x * (p0.x - m1.x) + d.y * (p0.y - m1.y) + d.z * (p0.z - m1.z))
    c = (p0.x - m0.x)^2 + (p0.y - m0.y)^2 + (p0.z - m0.z)^2 + (p0.x - m1.x)^2 + (p0.y - m1.y)^2 + (p0.z - m1.z)^2 - r^2

    t_values = solve_quadratic(a, b, c)

    if length(t_values) == 0
        return nothing
    end

    t_values = sort(t_values)
    t0 = t_values[1]
    t1 = t_values[2]

    q0 = p0 + d * t0
    q1 = p1 + d * t1

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

    rosetteMid = intersect_line_ellipse(pM + vy(1), pM - vy(1), arcLR[2], arcR[2], rad + radL + windowdict.bdInner)[1]
    rosetteRad = distance(rosetteMid, arcLR[2]) - radL - windowdict.bdInner

    return [arcLL, arcLR, arcRL, arcRR, rosetteMid, rosetteRad]
end

function gw_compute_fillets(arcLL, arcLR, arcRL, arcRR, rosetteMid, rosetteRad, pL, pR, nrml)
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
