
function save_commit_hash(res_dir)
    # saves current commit hash for reproducability
    current_commit = read(`git rev-parse HEAD`, String)
    open(joinpath(res_dir, "git_commit.txt"),"a") do io
        println(io,"current_commit=",current_commit)
    end
end

function auto_commit(resdir, commit_message)
    run(`git add -u`)
    current_commit = read(`git commit -m "auto_commit: " + commit_message`)