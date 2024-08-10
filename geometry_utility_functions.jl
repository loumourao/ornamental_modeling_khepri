module GeometryUtilityFunctions
    using KhepriAutoCAD

    export solve_quadratic, displace_point_by_vector

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

        return nothing
    end

    function displace_point_by_vector(point, vector, magnitude)
        unit_vector = vector / norm(vector)

        return point + unit_vector * magnitude
    end

    module Circles
        using KhepriAutoCAD
        using ..GeometryUtilityFunctions

        export intersect_circles, intersect_line_ellipse

        function intersect_circles(m0, r0, m1, r1)
            d = distance(m0, m1)

            # d == 0 prevents the set of interesection points when both circles overlap
            if d > r0 + r1 || d < abs(r0 - r1) || d == 0
                return nothing
            end

            a = (r0^2 - r1^2 + d^2) / (2 * d)
            p = m0 + a * (m1 - m0) / d
            c = sqrt(r0^2 - a^2)
            c_vector = m1 - m0
            c_vector = vxy(-c_vector.y, c_vector.x)
            c_vector = c_vector / norm(c_vector) * c
        
            first_intersection_point = p + c_vector
            second_intersection_point = p - c_vector
        
            if first_intersection_point.y >= second_intersection_point.y
                greater_y_intersection_point = first_intersection_point
                lower_y_intersection_point = second_intersection_point
            else
                greater_y_intersection_point = second_intersection_point
                lower_y_intersection_point = first_intersection_point
            end
        
            return (greater_y_intersection_point = greater_y_intersection_point, lower_y_intersection_point = lower_y_intersection_point)
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
        
            return (first_intersection_point = q0, second_intersection_point = q1, 
                        first_intersection_point_parameter = t0, second_intersection_point_parameter = t1)
        end

        module Arcs
            using KhepriAutoCAD
            using ..GeometryUtilityFunctions

            export arc, arc_start_point, arc_end_point, intersect_arcs, offset_arc, 
                    get_angular_adjustment_from_offset, arc_bidirectionally_extended_uniform_offset, 
                        arc_amplitude_extended_offset, arc_start_angle_extended_offset

            function arc(center, start_point, end_point)
                a = start_point - center
                b = end_point - center
            
                ccw_start_angle = atan(a.y, a.x)
                ccw_end_angle = atan(b.y, b.x)
                radius = norm(a)
            
                if ccw_end_angle > ccw_start_angle
                    amplitude = ccw_end_angle - ccw_start_angle
                else
                    amplitude = 2π - ccw_start_angle + ccw_end_angle
                end
            
                return KhepriAutoCAD.arc(center, radius, ccw_start_angle, amplitude)
            end

            function arc_start_point(arc)
                center = arc_center(arc)
                radius = arc_radius(arc)
                start_angle = arc_start_angle(arc)
            
                return center + vpol(radius, start_angle)
            end

            function arc_end_point(arc)
                center = arc_center(arc)
                radius = arc_radius(arc)
                start_angle = arc_start_angle(arc)
                amplitude = arc_amplitude(arc)
            
                return center + vpol(radius, start_angle + amplitude)
            end

            function intersect_arcs(arc0, arc1)
                m0 = arc_center(arc0)
                r0 = arc_radius(arc0)
                arc0_start_angle = arc_start_angle(arc0)
                arc0_amplitude = arc_amplitude(arc0)
            
                m1 = arc_center(arc1)
                r1 = arc_radius(arc1)
                arc1_start_angle = arc_start_angle(arc1)
                arc1_amplitude = arc_amplitude(arc1)
            
                circle_intersection_points = intersect_circles(m0, r0, m1, r1)
            
                if isnothing(circle_intersection_points)
                    return nothing
                end
            
                valid_intersection_points = []
            
                for intersection_point in circle_intersection_points
                    a = intersection_point - m0
                    a_θ = atan(a.y, a.x)
                    a_θ_amplitude = arc0_start_angle > a_θ ? 2π - arc0_start_angle + a_θ : a_θ - arc0_start_angle
            
                    b = intersection_point - m1
                    b_θ = atan(b.y, b.x)
                    b_θ_amplitude = arc1_start_angle > b_θ ? 2π - arc1_start_angle + b_θ : b_θ - arc1_start_angle
            
                    if a_θ_amplitude <= arc0_amplitude && b_θ_amplitude <= arc1_amplitude
                        push!(valid_intersection_points, intersection_point)
                    end
                end
            
                if isempty(valid_intersection_points)
                    return nothing
                end
            
                return length(valid_intersection_points) > 1 ? valid_intersection_points : valid_intersection_points[1]
            end

            function offset_arc(arc, offset_value)
                center = arc_center(arc)
                radius = arc_radius(arc)
                start_point = arc_start_point(arc)
                end_point = arc_end_point(arc)
            
                ca = start_point - center
                cb = end_point - center
                new_radius = radius - offset_value
                new_start_point = displace_point_by_vector(center, ca, new_radius)
                new_end_point = displace_point_by_vector(center, cb, new_radius)
            
                return GeometryUtilityFunctions.Circles.Arcs.arc(center, new_start_point, new_end_point)
            end

            function get_angular_adjustment_from_offset(arc, offset_value)
                amplitude = arc_amplitude(arc)
                radius = arc_radius(arc)
            
                Δθ = amplitude / radius
                tangential_offset = offset_value * sin(Δθ / 2)
                angular_adjustment = tangential_offset / radius
            
                return angular_adjustment
            end

            function arc_bidirectionally_extended_uniform_offset(arc, offset_value)
                delete_shape(arc)
                center = arc_center(arc)
                radius = arc_radius(arc)
                start_angle = arc_start_angle(arc)
                amplitude = arc_amplitude(arc)
            
                angular_adjustment = get_angular_adjustment_from_offset(arc, offset_value)
            
                new_radius = radius - offset_value
                new_start_angle = start_angle - angular_adjustment
                new_amplitude = amplitude + angular_adjustment * 2
            
                return KhepriAutoCAD.arc(center, new_radius, new_start_angle, new_amplitude)
            end

            function arc_amplitude_extended_offset(arc, offset_value)
                delete_shape(arc)
                center = arc_center(arc)
                radius = arc_radius(arc)
                start_angle = arc_start_angle(arc)
                amplitude = arc_amplitude(arc)
            
                angular_adjustment = get_angular_adjustment_from_offset(arc, offset_value) * 2
            
                new_radius = radius - offset_value
                new_start_angle = start_angle + angular_adjustment
            
                return KhepriAutoCAD.arc(center, new_radius, new_start_angle, amplitude)
            end

            function arc_start_angle_extended_offset(arc, offset_value)
                delete_shape(arc)
                center = arc_center(arc)
                radius = arc_radius(arc)
                start_angle = arc_start_angle(arc)
                amplitude = arc_amplitude(arc)
            
                angular_adjustment = get_angular_adjustment_from_offset(arc, offset_value) * 2
            
                new_radius = radius - offset_value
                new_start_angle = start_angle - angular_adjustment
            
                return KhepriAutoCAD.arc(center, new_radius, new_start_angle, amplitude)
            end
        end
        using .Arcs
        export arc, arc_start_point, arc_end_point, intersect_arcs, offset_arc, 
                get_angular_adjustment_from_offset, arc_bidirectionally_extended_uniform_offset, 
                    arc_amplitude_extended_offset, arc_start_angle_extended_offset
    end
    using .Circles
    export intersect_circles, intersect_line_ellipse
    export arc, arc_start_point, arc_end_point, intersect_arcs, offset_arc, 
            get_angular_adjustment_from_offset, arc_bidirectionally_extended_uniform_offset, 
                arc_amplitude_extended_offset, arc_start_angle_extended_offset
end
