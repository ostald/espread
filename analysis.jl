



dir = readdir("results/")
folders = filter(x-> contains(x, "14:29"), dir)
for f in folders
    open(joinpath("results", f, "results.txt"), "r") do file
        read(file)
    end
end

