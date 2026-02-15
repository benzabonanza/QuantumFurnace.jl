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

function KrausScratch(::Type{CT}, dim::Int) where {CT<:Complex}
    Zm() = zeros(CT, dim, dim)
    return KrausScratch{CT}(Zm(), Zm(), Zm(), Zm(), Zm(), Zm(), Zm(), Zm())
end