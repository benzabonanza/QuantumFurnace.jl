#!/usr/bin/env julia
# analyze_qf_e4z_41.jl
#
# Post-sweep analyzer for qf-e4z.41 GNS smooth-Metro sweep.
# Produces a side-by-side summary against the qf-e4z.34 CKG arm at the
# matching (β_phys=0.5, seed=46) cells.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BSON
using QuantumFurnace
using Printf

const GNS_DIR = joinpath(@__DIR__, "output", "sweep_qf_e4z_41_gns_smooth_metro_betaphys0_5")
const CKG_DIR = joinpath(@__DIR__, "output", "sweep_qf_e4z_34_ckg_vs_dll_plot", "ckg")

function _load(path)
    isfile(path) || return nothing
    return BSON.load(path, QuantumFurnace)[:result]
end

function _gns_path(n, sigma_factor, seed)
    s = @sprintf("%.6f", float(sigma_factor)); s = rstrip(s, '0'); s = rstrip(s, '.')
    s = isempty(s) ? "0" : s
    return joinpath(GNS_DIR, "sweep_n$(n)_sigma$(s)_seed$(seed)_L_GNS_Energy.bson")
end

function _ckg_path(n, beta_phys, seed)
    return joinpath(CKG_DIR, "sweep_n$(n)_betaphys$(beta_phys)_seed$(seed)_L_KMS_Energy.bson")
end

function main()
    @printf("\n%s\n", "="^96)
    @printf("qf-e4z.41 GNS vs qf-e4z.34 CKG side-by-side (β_phys=0.5, seed=46)\n")
    @printf("%s\n", "="^96)

    @printf("\n[1] n=3 σ-ramp tail (GNS approx-DB floor scaling)\n")
    @printf("%s\n", "-"^96)
    @printf("  σ_f      σ_alg     r_D  gap_phys  τ_mix_phys     floor       src\n")
    for σf in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125)
        r = _load(_gns_path(3, σf, 46))
        r === nothing && continue
        tmp = isfinite(r[:mixing_time_phys]) ? @sprintf("%9.4g", r[:mixing_time_phys]) : "     Inf "
        @printf("  %-7g  %.5f  %2d   %.4f    %s     %.3e   %s\n",
                σf, r[:sigma_alg], r[:r_D], r[:gap_phys], tmp,
                r[:floor_distance], string(r[:mixing_time_source]))
    end

    @printf("\n[2] σ_f=1/16 comparison cells: GNS (this sweep) vs CKG (qf-e4z.34, σ_f=1)\n")
    @printf("%s\n", "-"^96)
    @printf("  n  | GNS @ σ_f=1/16            | CKG @ σ_f=1               | ratios (GNS/CKG)\n")
    @printf("     | gap_phys  τ_mix  floor    | gap_phys  τ_mix           | gap     τ_mix\n")
    for n in 3:8
        rg = _load(_gns_path(n, 0.0625, 46))
        rc = _load(_ckg_path(n, "0.5", 46))
        if rg === nothing || rc === nothing
            @printf("  %d  | %s\n", n, rg === nothing ? "GNS missing" : "CKG missing")
            continue
        end
        gns_tmp = isfinite(rg[:mixing_time_phys]) ? @sprintf("%6.3f", rg[:mixing_time_phys]) : "  Inf "
        ckg_tmp = isfinite(rc[:mixing_time_phys]) ? @sprintf("%6.3f", rc[:mixing_time_phys]) : "  Inf "
        gap_ratio = rg[:gap_phys] / rc[:gap_phys]
        tau_ratio = isfinite(rg[:mixing_time_phys]) && isfinite(rc[:mixing_time_phys]) ?
                    rg[:mixing_time_phys] / rc[:mixing_time_phys] : NaN
        @printf("  %d  | %7.4f   %s  %.2e  | %7.4f   %s          | %.3f   %.3f\n",
                n, rg[:gap_phys], gns_tmp, rg[:floor_distance],
                rc[:gap_phys], ckg_tmp, gap_ratio, tau_ratio)
    end

    @printf("\n[3] σ_f=1/16 cells: convergence status\n")
    @printf("%s\n", "-"^96)
    @printf("  n   r_D  kdim  src           floor       all_converged  matvecs  wall\n")
    for n in 3:8
        r = _load(_gns_path(n, 0.0625, 46))
        r === nothing && (continue)
        kdim = haskey(r, :krylovdim_p1) ? r[:krylovdim_p1] :
               haskey(r, :krylovdim)    ? r[:krylovdim]    : -1
        @printf("  %d   %d   %3d   %-12s  %.3e   %-6s         %5d    %.1fs\n",
                n, r[:r_D], kdim, string(r[:mixing_time_source]),
                r[:floor_distance], string(r[:all_converged]),
                r[:total_matvecs], r[:wall_time])
    end

    # σ_f=1/32 bail-out cells, if any
    @printf("\n[4] σ_f=1/32 bail-out cells (if any)\n")
    @printf("%s\n", "-"^96)
    any_bail = false
    for n in 3:8
        r = _load(_gns_path(n, 0.03125, 46))
        r === nothing && continue
        any_bail = true
        @printf("  n=%d r_D=%d kdim=%d  src=%-12s  floor=%.3e  gap_phys=%.4f  τ_mix_phys=%g\n",
                n, r[:r_D], r[:krylovdim_p1], string(r[:mixing_time_source]),
                r[:floor_distance], r[:gap_phys], r[:mixing_time_phys])
    end
    any_bail || @printf("  (none — σ_f=1/16 was sufficient for all n)\n")

    @printf("\n%s\n", "="^96)
end

main()
