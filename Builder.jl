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

function teardrop()
    global counter_coordinate, counter_volume_element, counter_volume_patch
    counter_coordinate = 1
    counter_volume_element = 1
    counter_volume_patch = 1

    # elastic parameters
    κ = 10 * 4.114
    Jₛu = -0.1
    Jₛl = -0.1
    kₜ = 10 * 4.114
    Kₐ = 30 * 4.114
    h₀ = 1.5

    # resolution control
    r_low_res = 15.0
    ds_mid_high_res = 0.5
    ds_mid_low_res = 1.0

    rp = 3.0
    rb = 1.5
    shift = 5.0
    r_tot = 40.0

    rim0, zim0 = linear(rp + h₀, 0.0, rp + shift, rb + h₀, 50)
    riu, ziu = linear(rp, 0.0, rp + shift, rb + 2 * h₀, 50)
    rii0, zii0 = linear(rp + 2 * h₀, 0.0, rp + shift, rb, 50)

    rom_high_res = range(rim0[end], rim0[end] + r_low_res, length=50)[2:end]
    rom_low_res = range(rim0[end] + r_low_res + 1.0, r_tot, length=50)
    zom_high_res = fill(zim0[end], 50 - 1)
    zom_low_res = fill(zim0[end], 50)
    rou_high_res = range(riu[end], rim0[end] + r_low_res, length=50)[2:end]
    rou_low_res = range(rim0[end] + r_low_res + 1.0, r_tot, length=50)
    zou_high_res = fill(ziu[end], 50 - 1)
    zou_low_res = fill(ziu[end], 50)
    rol_high_res = range(rii0[end], rim0[end] + r_low_res, length=50)[2:end]
    rol_low_res = range(rim0[end] + r_low_res + 1.0, r_tot, length=50)
    zol_high_res = fill(zii0[end], 50 - 1)
    zol_low_res = fill(zii0[end], 50)

    # now let's build Coordinates out of them
    coords_iu = [SurfaceCoordinate(riu[1], ziu[1], fixed_r = true, fixed_z = true)]
    coords_iu = vcat(coords_iu, [SurfaceCoordinate(riu[i], ziu[i]) for i in 2:50])
    coords_ou_high_res = [SurfaceCoordinate(rou_high_res[i], zou_high_res[i]) for i in 1:50-1]
    coords_ou_low_res = [SurfaceCoordinate(rou_low_res[i], zou_low_res[i]) for i in 1:49]
    push!(coords_ou_low_res, SurfaceCoordinate(rou_low_res[50], zou_low_res[50], fixed_r = true, fixed_z = false))
    pushfirst!(coords_ou_high_res, coords_iu[end])
    pushfirst!(coords_ou_low_res, coords_ou_high_res[end])

    coords_im0 = [MidplaneCoordinate(rim0[1], zim0[1], fixed_r = false, fixed_z = true)]
    coords_im0 = vcat(coords_im0, [MidplaneCoordinate(rim0[i], zim0[i]) for i in 2:50-1])
    push!(coords_im0, MidplaneCoordinate(rim0[end], zim0[end], fixed_r = false, fixed_z = false))
    coords_om_high_res = [MidplaneCoordinate(rom_high_res[i], zom_high_res[i]) for i in 1:50-1]
    coords_om_low_res = [MidplaneCoordinate(rom_low_res[i], zom_low_res[i]) for i in 1:49]
    push!(coords_om_low_res, MidplaneCoordinate(rom_low_res[50], zom_low_res[50], fixed_r = true, fixed_z = false))
    pushfirst!(coords_om_high_res, coords_im0[end])
    pushfirst!(coords_om_low_res, coords_om_high_res[end])

    coords_ii0 = [SurfaceCoordinate(rii0[i], zii0[i]) for i in 2:50]
    pushfirst!(coords_ii0, SurfaceCoordinate(rii0[1], zii0[1], fixed_r = false, fixed_z = true))
    coords_ol_high_res = [SurfaceCoordinate(rol_high_res[i], zol_high_res[i]) for i in 1:50-1]
    coords_ol_low_res = [SurfaceCoordinate(rol_low_res[i], zol_low_res[i]) for i in 1:49]
    push!(coords_ol_low_res, SurfaceCoordinate(rol_low_res[50], zol_low_res[50], fixed_r = true, fixed_z = true))
    pushfirst!(coords_ol_high_res, coords_ii0[end])
    pushfirst!(coords_ol_low_res, coords_ol_high_res[end])

    # create the patches
    patch_u = LeafletPatch(κ, 0.0, Jₛu, kₜ, Kₐ, h₀)
    patch_l = LeafletPatch(κ, 0.0, Jₛl, kₜ, Kₐ, h₀)

    # make the remeshing patches
    rem_patch_iu = RemeshingPatch(coords_iu, coords_im0, coords_ii0, ds_mid_high_res, patch_u, patch_l)
    rem_patch_o_high_res = RemeshingPatch(coords_ou_high_res, coords_om_high_res, coords_ol_high_res, ds_mid_high_res, patch_u, patch_l)
    rem_patch_o_low_res = RemeshingPatch(coords_ou_low_res, coords_om_low_res, coords_ol_low_res, ds_mid_low_res, patch_u, patch_l)

    return remesh([rem_patch_iu, rem_patch_o_high_res, rem_patch_o_low_res])
end

function monolayer_stretch(a=1.0)
    global counter_coordinate, counter_volume_element, counter_volume_patch
    counter_coordinate = 1
    counter_volume_element = 1
    counter_volume_patch = 1

    R = 50.0

    # elastic parameters
    κ = 10 * 4.114
    Jₛ = 2.0/R
    kₜ = 10 * 4.114
    Kₐ = 1.0
    h₀ = 1.5

    # resolution control
    ds = 4.0

    V_tot = π * R * R * h₀ / a

    rm, zm = range(0.0, R, 50), fill(0.0, 50)
    ru, zu = range(0.0, R, 50), fill(h₀, 50)
 
    
    # now let's build Coordinates out of them
    coords_u = [SurfaceCoordinate(ru[1], zu[1], fixed_r = true, fixed_z = false)]
    coords_u = vcat(coords_u, [SurfaceCoordinate(ru[i], zu[i]) for i in 2:49])
    coords_u = push!(coords_u, SurfaceCoordinate(ru[50], zu[50], fixed_r = true, fixed_z = false))
    
    coords_m = [MidplaneCoordinate(rm[1], zm[1], fixed_r = true, fixed_z = false)]
    coords_m = vcat(coords_m, [MidplaneCoordinate(rm[i], zm[i]) for i in 2:49])
    coords_m = push!(coords_m, MidplaneCoordinate(rm[50], zm[50], fixed_r = true, fixed_z = false))

    # create the patches
    patch_u = LeafletPatch(κ, 0.0, Jₛ, kₜ, Kₐ, h₀, constrain_v = true, v₀=V_tot)

    # make the remeshing patches
    rem_patch = RemeshingPatch(coords_u, coords_m, nothing, ds, patch_u, nothing)

    return remesh([rem_patch])
end

function bilayer()
    global counter_coordinate, counter_volume_element, counter_volume_patch
    counter_coordinate = 1
    counter_volume_element = 1
    counter_volume_patch = 1

    R = 20.0

    # elastic parameters
    κ = 10 * 4.114
    kₜ = 10 * 4.114
    Kₐ = 30 * 4.114
    h₀ = 1.5
    Jₛ1 = 2.0/(R + h₀)
    Jₛ2 = -2.0/(R - h₀)

    # resolution control
    ds = 0.5

    V_u = 0.5 * (4 * π * (R + 0.5 * h₀) * (R + 0.5 * h₀) * h₀ + π*h₀*h₀*h₀/3)
    V_l = 0.5 * (4 * π * (R - 0.5 * h₀) * (R - 0.5 * h₀) * h₀ + π*h₀*h₀*h₀/3)

    r_tot = R - h₀

    rm, zm = range(0.0, r_tot, 50), fill(0.0, 50)
    ru, zu = range(0.0, r_tot, 50), fill(h₀, 50)
    ri, zi = range(0.0, r_tot, 50), fill(-h₀, 50)

    
    # now let's build Coordinates out of them
    coords_u = [SurfaceCoordinate(ru[1], zu[1], fixed_r = true, fixed_z = false)]
    coords_u = vcat(coords_u, [SurfaceCoordinate(ru[i], zu[i]) for i in 2:50])
    
    coords_m = [MidplaneCoordinate(rm[1], zm[1], fixed_r = true, fixed_z = false)]
    coords_m = vcat(coords_m, [MidplaneCoordinate(rm[i], zm[i]) for i in 2:50])

    coords_l = [SurfaceCoordinate(ri[1], zi[1], fixed_r = true, fixed_z = false)]
    coords_l = vcat(coords_l, [SurfaceCoordinate(ri[i], zi[i]) for i in 2:50])

    # create the patches
    patch_u = LeafletPatch(κ, 0.0, Jₛ1, kₜ, Kₐ, h₀, constrain_v = true, v₀=V_u)
    patch_l = LeafletPatch(κ, 0.0, Jₛ2, kₜ, Kₐ, h₀, constrain_v = true, v₀=V_l)

    # make the remeshing patches
    rem_patch = RemeshingPatch(coords_u, coords_m, coords_l, ds, patch_u, patch_l)

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

function vesicle_manuscript_ADE(;ν=1.0, m₀_4π=1.0, is_open = true, Rₐ=50.0, h₀=1.5, κ=10.0, α=4.0)
    global counter_coordinate, counter_volume_element, counter_volume_patch

    counter_coordinate = 1
    counter_volume_element = 1
    counter_volume_patch = 1

    # elastic parameters
    κ *= 4.114
    kₜ = 10 * 4.114
    D = 2*h₀
    Kₐ = κ * 4 * α * π / D / D # α = 4
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