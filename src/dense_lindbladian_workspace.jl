"""
    LiouvillianScratch{T<:Complex}

Reusable matrix buffers for dense Liouvillian construction.
"""
struct LiouvillianScratch{T<:Complex}
    jump_tmp::Matrix{T}
    jump_conj::Matrix{T}
    jump_dag_jump::Matrix{T}
    jump2_jump1::Matrix{T}
end


function LiouvillianScratch(::Type{CT}, dim::Int) where {CT<:Complex}
    Zm() = zeros(CT, dim, dim)
    return LiouvillianScratch{CT}(Zm(), Zm(), Zm(), Zm())
end


"""
    DenseLindbladianWorkspace{T<:Complex}

Identity and scratch buffers for the independent dense constructor.
"""
struct DenseLindbladianWorkspace{T<:Complex}
    Id::Matrix{T}
    scratch::LiouvillianScratch{T}
end

function DenseLindbladianWorkspace(::Type{CT}, dim::Int) where {CT<:Complex}
    return DenseLindbladianWorkspace{CT}(
        Matrix{CT}(I, dim, dim),
        LiouvillianScratch(CT, dim),
    )
end
