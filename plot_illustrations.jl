using CairoMakie


#Mean free path

x = 0:0.02:4
fig, ax, lin = lines(x, exp.(-x),
    axis = (xlabel = "Number of mean free paths [1]",
        ylabel = "PDF [1]"
        ),)
save("figures/mean_free_path_distribution.png", fig, px_per_unit = 3.3)


lines(x, cumsum(exp.(-x))*0.02,
    axis = (xlabel = "Number of mean free paths [1]",
        ylabel = "CDF [1]"
        ),)
