"""
    Builder module

    This module provides functions for constructing initial configuration with the geometrical elements used in FEM scheme.

"""

function circle_points(R::Float64, N::Int, h::Float64)
    theta = acos(h/R)
    dtheta = theta / (N - 1)
    return [R * sin(i * dtheta) for i in 0:N-1], [R * cos(i * dtheta) for i in 0:N-1]
end

function circle_points(R::Float64, N::Int, h_start::Float64, h_end::Float64, shift::Float64=0.0)
    theta_start = acos(h_start/R)
    theta_end = acos(h_end/R)
    dtheta = (theta_end - theta_start) / (N - 1)
    return [R * sin(theta_start + i * dtheta) for i in 0:N-1], [R * cos(theta_start + i * dtheta) + shift for i in 0:N-1]
end

function linear(r0::Float64, z0::Float64, r1::Float64, z1::Float64, N::Int)
    rs = collect(range(r0, r1, N))
    zs = collect(range(z0, z1, N))
    return rs, zs
end

# simple vesicle
function vesicle_simple(Vupper_extra, Vwater, is_open = true,; r_tot=50.0)
    global counter_coordinate, counter_volume_element, counter_volume_patch

    counter_coordinate = 1
    counter_volume_element = 1
    counter_volume_patch = 1

    # elastic parameters
    κ = 40.0
    Jₛ = 0.0
    kₜ = 40.0
    Kₐ = 120.0
    h₀ = 1.5
    # resolution control
    ds_mid = 2.0

    v₀_membrane_upper = 4*π*(r_tot+0.5*h₀)*(r_tot+0.5*h₀)*h₀*(1.0+Vupper_extra)
    v₀_membrane_lower = 4*π*(r_tot-0.5*h₀)*(r_tot-0.5*h₀)*h₀

    v₀_water = 4.0*π*(r_tot-h₀)*(r_tot-h₀)*(r_tot-h₀)/3.0*Vwater

    if is_open
        rm, zm = circle_points(r_tot, 100, -r_tot)
        ru, zu = circle_points(r_tot+h₀, 100, -r_tot-h₀)
        rl, zl = circle_points(r_tot-h₀, 100, -r_tot+h₀)
    else
        r_initial_bud = 10.0
        y1, y2 = 0.0, 0.0
        rm, zm = circle_points(r_initial_bud, 7, r_initial_bud, - r_initial_bud+0.5, r_tot + r_initial_bud)
        rm1, zm1 = circle_points(r_tot, 7, r_tot - 0.5, -r_tot, 0.0)
        rm = vcat(rm, rm1)
        zm = vcat(zm, zm1)
        ru, zu = circle_points(r_initial_bud + h₀, 7, r_initial_bud + h₀, -r_initial_bud - h₀ + 0.5, r_tot + r_initial_bud + h₀)
        ru1, zu1 = circle_points(r_tot + h₀, 7, r_tot + h₀ - 0.5, -(r_tot + h₀), 0.0)
        ru = vcat(ru, ru1)
        zu = vcat(zu, zu1)
        rl, zl = circle_points(r_initial_bud - h₀, 7, r_initial_bud - h₀, -r_initial_bud + h₀ + 0.5, r_tot+ r_initial_bud - h₀)
        rl1, zl1 = circle_points(r_tot - h₀, 7, r_tot - h₀ - 0.5, -(r_tot - h₀), 0.0)
        rl = vcat(rl, rl1)
        zl = vcat(zl, zl1)
    end

    ds = [sqrt((rm[i+1] - rm[i]) * (rm[i+1] - rm[i]) + (zm[i+1] - zm[i]) * (zm[i+1] - zm[i])) for i in 1:length(rm)-1]
    pushfirst!(ds, 0.0)
    
    ds = cumsum(ds)

    spl_rus = Spline1D(ds, ru)
    spl_zus = Spline1D(ds, zu)
    spl_rls = Spline1D(ds, rl)
    spl_zls = Spline1D(ds, zl)
    spl_rms = Spline1D(ds, rm)
    spl_zms = Spline1D(ds, zm)

    ds = LinRange(0.0, ds[end], 100)

    ru, zu = spl_rus(ds), spl_zus(ds)
    rm, zm = spl_rms(ds), spl_zms(ds)
    rl, zl = spl_rls(ds), spl_zls(ds)

    # now let's build Coordinates out of them
    coords_u = [SurfaceCoordinate(ru[i], zu[i]) for i in 2:99]
    push!(coords_u, SurfaceCoordinate(ru[end], zu[end], fixed_r=true))
    pushfirst!(coords_u, SurfaceCoordinate(ru[1], zu[1], fixed_r=true))
    coords_m = [MidplaneCoordinate(rm[i], zm[i]) for i in 2:99]
    push!(coords_m, MidplaneCoordinate(rm[end], zm[end], fixed_r=true))
    pushfirst!(coords_m, MidplaneCoordinate(rm[1], zm[1], fixed_r=true))
    coords_l = [SurfaceCoordinate(rl[i], zl[i]) for i in 2:99]
    push!(coords_l, SurfaceCoordinate(rl[end], zl[end], fixed_r=true, fixed_z=true))
    pushfirst!(coords_l, SurfaceCoordinate(rl[1], zl[1], fixed_r=true))

    # create the patches
    patch_u = LeafletPatch(κ, Jₛ, 0.0, kₜ, Kₐ, h₀, constrain_v = true, v₀=v₀_membrane_upper)
    patch_l = LeafletPatch(κ, Jₛ, 0.0, kₜ, Kₐ, h₀, constrain_v = true, v₀=v₀_membrane_lower)
    patch_water = WaterPatch(v₀_water)

    # make the remeshing patches
    rem_patch = RemeshingPatch(coords_u, coords_m, coords_l, ds_mid, patch_u, patch_l, nothing, patch_water, fixed_midplane_z=false)

    return remesh([rem_patch])
end

# manuscript vesicle
function vesicle_manuscript(;Vupper_extra=0.0, Vwater=1.0, is_open = true, is_binary = false, Φ₀=0.0, Jₛᵘ=0.0, Jₛ_upper=0.0, Jₛ_lower=0.0, r_tot=50.0)
    global counter_coordinate, counter_volume_element, counter_volume_patch

    counter_coordinate = 1
    counter_volume_element = 1
    counter_volume_patch = 1

    # elastic parameters
    κ = 59.6
    kₜ = 39.1
    Kₐ = 123.0
    h₀ = 1.51
    # resolution control
    ds_mid = 1.0

    v₀_membrane_upper = 4*π*(r_tot+0.5*h₀)*(r_tot+0.5*h₀)*h₀*(1.0+Vupper_extra)
    v₀_membrane_lower = 4*π*(r_tot-0.5*h₀)*(r_tot-0.5*h₀)*h₀

    v₀_water = 4.0*π*(r_tot-h₀)*(r_tot-h₀)*(r_tot-h₀)/3.0*Vwater

    if is_open
        rm, zm = circle_points(r_tot, 100, -r_tot)
        ru, zu = circle_points(r_tot+h₀, 100, -r_tot-h₀)
        rl, zl = circle_points(r_tot-h₀, 100, -r_tot+h₀)
    else
        r_initial_bud = 10.0
        y1, y2 = 0.0, 0.0
        rm, zm = circle_points(r_initial_bud, 7, r_initial_bud, - r_initial_bud+0.5, r_tot + r_initial_bud)
        rm1, zm1 = circle_points(r_tot, 7, r_tot - 0.5, -r_tot, 0.0)
        rm = vcat(rm, rm1)
        zm = vcat(zm, zm1)
        ru, zu = circle_points(r_initial_bud + h₀, 7, r_initial_bud + h₀, -r_initial_bud - h₀ + 0.5, r_tot + r_initial_bud + h₀)
        ru1, zu1 = circle_points(r_tot + h₀, 7, r_tot + h₀ - 0.5, -(r_tot + h₀), 0.0)
        ru = vcat(ru, ru1)
        zu = vcat(zu, zu1)
        rl, zl = circle_points(r_initial_bud - h₀, 7, r_initial_bud - h₀, -r_initial_bud + h₀ + 0.5, r_tot+ r_initial_bud - h₀)
        rl1, zl1 = circle_points(r_tot - h₀, 7, r_tot - h₀ - 0.5, -(r_tot - h₀), 0.0)
        rl = vcat(rl, rl1)
        zl = vcat(zl, zl1)
    end

    ds = [sqrt((rm[i+1] - rm[i]) * (rm[i+1] - rm[i]) + (zm[i+1] - zm[i]) * (zm[i+1] - zm[i])) for i in 1:length(rm)-1]
    pushfirst!(ds, 0.0)
    
    ds = cumsum(ds)

    spl_rus = Spline1D(ds, ru)
    spl_zus = Spline1D(ds, zu)
    spl_rls = Spline1D(ds, rl)
    spl_zls = Spline1D(ds, zl)
    spl_rms = Spline1D(ds, rm)
    spl_zms = Spline1D(ds, zm)

    ds = LinRange(0.0, ds[end], 100)

    ru, zu = spl_rus(ds), spl_zus(ds)
    rm, zm = spl_rms(ds), spl_zms(ds)
    rl, zl = spl_rls(ds), spl_zls(ds)

    # now let's build Coordinates out of them
    if is_binary 
        coords_u = [SurfaceCoordinate(ru[i], zu[i], Φ₀) for i in 2:99]
        push!(coords_u, SurfaceCoordinate(ru[end], zu[end], Φ₀, fixed_r=true))
        pushfirst!(coords_u, SurfaceCoordinate(ru[1], zu[1], Φ₀, fixed_r=true))
    else
        coords_u = [SurfaceCoordinate(ru[i], zu[i]) for i in 2:99]
        push!(coords_u, SurfaceCoordinate(ru[end], zu[end], fixed_r=true))
        pushfirst!(coords_u, SurfaceCoordinate(ru[1], zu[1], fixed_r=true))
    end
    coords_m = [MidplaneCoordinate(rm[i], zm[i]) for i in 2:99]
    push!(coords_m, MidplaneCoordinate(rm[end], zm[end], fixed_r=true))
    pushfirst!(coords_m, MidplaneCoordinate(rm[1], zm[1], fixed_r=true))
    coords_l = [SurfaceCoordinate(rl[i], zl[i]) for i in 2:99]
    push!(coords_l, SurfaceCoordinate(rl[end], zl[end], fixed_r=true, fixed_z=true))
    pushfirst!(coords_l, SurfaceCoordinate(rl[1], zl[1], fixed_r=true))

    # create the patches
    if is_binary
        patch_u = LeafletPatch(κ, 0.0, Jₛ_upper, kₜ, Kₐ, h₀, κ, 0.0, Jₛᵘ, kₜ, Kₐ, h₀, constrain_v = true, v₀=v₀_membrane_upper, constrain_Φ=true, Φ₀=Φ₀)
    else
        patch_u = LeafletPatch(κ, 0.0, Jₛ_upper, kₜ, Kₐ, h₀, constrain_v = true, v₀=v₀_membrane_upper)
    end
    patch_l = LeafletPatch(κ, 0.0, Jₛ_lower, kₜ, Kₐ, h₀, constrain_v = true, v₀=v₀_membrane_lower)
    patch_water = WaterPatch(v₀_water)

    # make the remeshing patches
    rem_patch = RemeshingPatch(coords_u, coords_m, coords_l, ds_mid, patch_u, patch_l, nothing, patch_water, fixed_midplane_z=false)

    return remesh([rem_patch])
end

function vesicle_manuscript_ADE(;ν=1.0, m₀_4π=1.0, is_open = true, Rₐ=50.0, h₀=1.5)
    global counter_coordinate, counter_volume_element, counter_volume_patch

    counter_coordinate = 1
    counter_volume_element = 1
    counter_volume_patch = 1

    # elastic parameters
    κ = 10 * 4.114
    kₜ = 10 * 4.114
    D = 2*h₀
    Kₐ = κ * 16 * π / D / D # α = 4
    # resolution control
    ds_mid = 1.0

    ΔA₀ = (2*D*Rₐ)*m₀_4π*4*π
    A₀ = 4*π*Rₐ^2

    v₀_membrane_upper = h₀*(A₀ + ΔA₀/4)
    v₀_membrane_lower = h₀*(A₀ - ΔA₀/4)

    v₀_water = 4.0*π*(Rₐ-h₀)*(Rₐ-h₀)*(Rₐ-h₀)/3.0*ν

    if is_open
        rm, zm = circle_points(Rₐ, 100, -Rₐ)
        ru, zu = circle_points(Rₐ+h₀, 100, -Rₐ-h₀)
        rl, zl = circle_points(Rₐ-h₀, 100, -Rₐ+h₀)
    else
        r_initial_bud = 10.0
        y1, y2 = 0.0, 0.0
        rm, zm = circle_points(r_initial_bud, 7, r_initial_bud, - r_initial_bud+0.5, Rₐ + r_initial_bud)
        rm1, zm1 = circle_points(Rₐ, 7, Rₐ - 0.5, -Rₐ, 0.0)
        rm = vcat(rm, rm1)
        zm = vcat(zm, zm1)
        ru, zu = circle_points(r_initial_bud + h₀, 7, r_initial_bud + h₀, -r_initial_bud - h₀ + 0.5, Rₐ + r_initial_bud + h₀)
        ru1, zu1 = circle_points(Rₐ + h₀, 7, Rₐ + h₀ - 0.5, -(Rₐ + h₀), 0.0)
        ru = vcat(ru, ru1)
        zu = vcat(zu, zu1)
        rl, zl = circle_points(r_initial_bud - h₀, 7, r_initial_bud - h₀, -r_initial_bud + h₀ + 0.5, Rₐ+ r_initial_bud - h₀)
        rl1, zl1 = circle_points(Rₐ - h₀, 7, Rₐ - h₀ - 0.5, -(Rₐ - h₀), 0.0)
        rl = vcat(rl, rl1)
        zl = vcat(zl, zl1)
    end

    ds = [sqrt((rm[i+1] - rm[i]) * (rm[i+1] - rm[i]) + (zm[i+1] - zm[i]) * (zm[i+1] - zm[i])) for i in 1:length(rm)-1]
    pushfirst!(ds, 0.0)
    
    ds = cumsum(ds)

    spl_rus = Spline1D(ds, ru)
    spl_zus = Spline1D(ds, zu)
    spl_rls = Spline1D(ds, rl)
    spl_zls = Spline1D(ds, zl)
    spl_rms = Spline1D(ds, rm)
    spl_zms = Spline1D(ds, zm)

    ds = LinRange(0.0, ds[end], 100)

    ru, zu = spl_rus(ds), spl_zus(ds)
    rm, zm = spl_rms(ds), spl_zms(ds)
    rl, zl = spl_rls(ds), spl_zls(ds)

    # now let's build Coordinates out of them
    coords_u = [SurfaceCoordinate(ru[i], zu[i]) for i in 2:99]
    push!(coords_u, SurfaceCoordinate(ru[end], zu[end], fixed_r=true))
    pushfirst!(coords_u, SurfaceCoordinate(ru[1], zu[1], fixed_r=true))
    coords_m = [MidplaneCoordinate(rm[i], zm[i]) for i in 2:99]
    push!(coords_m, MidplaneCoordinate(rm[end], zm[end], fixed_r=true))
    pushfirst!(coords_m, MidplaneCoordinate(rm[1], zm[1], fixed_r=true))
    coords_l = [SurfaceCoordinate(rl[i], zl[i]) for i in 2:99]
    push!(coords_l, SurfaceCoordinate(rl[end], zl[end], fixed_r=true, fixed_z=false))
    pushfirst!(coords_l, SurfaceCoordinate(rl[1], zl[1], fixed_r=true))

    # create the patches
    patch_u = LeafletPatch(κ, 0.0, 0.0, kₜ, Kₐ, h₀, constrain_v = true, v₀=v₀_membrane_upper)
    patch_l = LeafletPatch(κ, 0.0, 0.0, kₜ, Kₐ, h₀, constrain_v = true, v₀=v₀_membrane_lower)
    patch_water = WaterPatch(v₀_water)

    # make the remeshing patches
    rem_patch = RemeshingPatch(coords_u, coords_m, coords_l, ds_mid, patch_u, patch_l, nothing, patch_water, fixed_midplane_z=false)

    return remesh([rem_patch])
end