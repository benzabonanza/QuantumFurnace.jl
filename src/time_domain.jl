function _truncate_time_labels_for_oft(
    time_labels::AbstractVector{<:Real},
    sigma::Real;
    tolerance::Real = 1e-12,
    filter::AbstractFilter = GaussianFilter(sigma),
    )

    cutoff = filter_time_cutoff(filter, tolerance)

    # Qualify `Base.filter` because `filter` names the kernel argument.
    truncated_labels = Base.filter(t -> abs(t) <= cutoff, time_labels)
    t_max = maximum(time_labels)
    if !isempty(truncated_labels) && (maximum(truncated_labels) == t_max)
        t_max = maximum(time_labels)
        residual = abs(time_kernel(filter, t_max))

        @warn """
        OFT Integration Warning: Time array was not truncated.
        Filter kernel at the ends should be small but it is: $(residual)
        """
    end

    return truncated_labels
end
