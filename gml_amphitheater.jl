using Khepri
backend(autocad)

function n_gon_solid(n_sides, center, radius, height, angle, inscribed)
    n_gon = surface(regular_polygon(n_sides, center, radius, angle, inscribed))
    extrusion(n_gon, height)
end

function stair_base(bottom_left_corner, width, length, height)
    extrusion(surface(rectangle(bottom_left_corner, width, length)), height)
end

function stair_step(bottom_left_corner, width, length, height)
    extrusion(surface(rectangle(bottom_left_corner, width, length)), height)
end

let 
    n_gon_radius = 2
    number_of_sides = 7
    n_gon_center = xy(0, 0)
    n_gon_height = 2
    n_gon_starting_angle = (2 * pi) / (number_of_sides * 2)
    
    amphitheater = n_gon_solid(number_of_sides, n_gon_center, n_gon_radius, n_gon_height, n_gon_starting_angle, true)
    
    stair_base_bottom_left_corner = pol(n_gon_radius, -n_gon_starting_angle)
    stair_base_width = 0.25
    stair_base_length = distance(stair_base_bottom_left_corner, pol(n_gon_radius, n_gon_starting_angle))
    stair_base_height = n_gon_height

    stair_step_bottom_left_corner = stair_base_bottom_left_corner + vz(stair_base_height)
    stair_step_width = stair_base_width
    stair_step_length = stair_base_length
    stair_step_height = 0.25

    for i = 1:7
        base = stair_base(stair_base_bottom_left_corner, stair_base_width, stair_base_length, stair_base_height)
        step = stair_step(stair_step_bottom_left_corner, stair_step_width, stair_step_length, stair_step_height)
        amphitheater = union(amphitheater, base, step)

        for i = 1:10
            if i % 2 != 0
                stair_step_bottom_left_corner = stair_step_bottom_left_corner + vx(stair_step_width)
                step = stair_step(stair_step_bottom_left_corner, stair_step_width, stair_step_length, stair_step_height)
            else
                stair_step_bottom_left_corner = stair_step_bottom_left_corner + vz(stair_step_height)
                step = stair_step(stair_step_bottom_left_corner, stair_step_width, stair_step_length, stair_step_height)
            end
            amphitheater = union(amphitheater, step)
        end

        stair_step_bottom_left_corner = stair_base_bottom_left_corner + vz(stair_base_height)
        amphitheater = rotate(amphitheater, (2 * pi) / (number_of_sides))
    end
end
