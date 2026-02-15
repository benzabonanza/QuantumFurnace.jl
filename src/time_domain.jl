function _truncate_time_labels_for_oft(time_labels::AbstractVector{<:Real}, sigma::Real; tolerance::Real = 1e-12)
    
    time_oft_prefactor = sqrt(sigma * sqrt(2 / pi) / (2 * pi))
    cutoff = sqrt(log(time_oft_prefactor / tolerance)) / sigma

    truncated_labels = filter(t -> abs(t) <= cutoff, time_labels)
    t_max = maximum(time_labels)
    if !isempty(truncated_labels) && (maximum(truncated_labels) == t_max)
        t_max = maximum(time_labels)
        residual = exp(-t_max^2 * sigma^2) * time_oft_prefactor

        @warn """
        OFT Integration Warning: Time array was not truncated.
        Gaussian at the ends should be small but it is: $(residual)
        """
    end

    return truncated_labels
end
