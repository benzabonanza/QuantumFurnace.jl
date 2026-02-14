struct KrausScratch{T}
    jump_oft::Matrix{T}
    LdagL::Matrix{T}
    R::Matrix{T}
    rho_jump::Matrix{T}
    K0::Matrix{T}
    tmp1::Matrix{T}
    tmp2::Matrix{T}
    rho_next::Matrix{T}
end

function KrausScratch(T::Type{ComplexF64}, dim::Int)
    Zm() = zeros(T, dim, dim)
    return KrausScratch(Zm(), Zm(), Zm(), Zm(), Zm(), Zm(), Zm(), Zm())
end